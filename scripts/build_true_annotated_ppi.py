"""Build target-only PPI contacts with reversible source residue identities.

Experimental assembly contacts, not independently curated interaction hotspots.
The old PPI dataset is read only. Exported query numbering is chain A, 1..N;
source chain/number/insertion code are retained in every label and mapping.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil
import warnings

import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
from scipy.spatial.distance import cdist

from build_investigation_db import read_block
from build_reviewed_enzymes import ROOT, save_json, table
from build_antibody_targets import contact_tables, digest, key, load

LABEL_HEADER = ["instance", "chain", "res_id", "ins_code", "res_name", "is_annotated",
                "label_type", "min_partner_distance_angstrom", "partner_chains",
                "source_chain", "source_res_id", "source_ins_code"]
MAP_HEADER = ["instance", "chain", "res_id", "ins_code", "res_name", "source_chain",
              "source_res_id", "source_ins_code", "source_entity_id"]
CONTACT_HEADER = ["instance", "target_chain", "target_res_id", "target_ins_code", "target_res_name",
                  "partner_chain", "partner_res_id", "partner_ins_code", "partner_res_name", "min_distance_angstrom"]
ANTIBODY = re.compile(r"antibody|nanobody|\bvhh\b|\bfab\b|\bfv\b|immunoglobulin (?:heavy|light)|(?:heavy|light) chain.*(?:variable|immunoglobulin)", re.I)


def residue_keys(array):
    return list(dict.fromkeys(key(array, i) for i in range(len(array))))


def protein_arrays(path):
    cif, block = read_block(path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        arr = pdbx.get_structure(cif, model=1, extra_fields=["label_entity_id", "b_factor"])
    arr = arr[struc.filter_amino_acids(arr) & ~np.isin(arr.element, ["H", "D"])]
    names = []
    if "entity" in block and "pdbx_description" in block["entity"]:
        names = block["entity"]["pdbx_description"].as_array(str).tolist()
    if any(ANTIBODY.search(n) for n in names):
        raise ValueError("antibody_entity_description")
    for chain in set(arr.chain_id):
        part = arr[arr.chain_id == chain]
        if len(set(part.label_entity_id)) != 1:
            raise ValueError("ambiguous_chain_entity")
        identities = [(key(part, i), str(part.atom_name[i])) for i in range(len(part))]
        if len(set(identities)) != len(identities):
            raise ValueError("duplicate_residue_atoms")
        residue_names = {}
        for k in residue_keys(part):
            if k[:3] in residue_names and residue_names[k[:3]] != k[3]:
                raise ValueError("ambiguous_source_residue")
            residue_names[k[:3]] = k[3]
    return arr


def select_target(arr):
    # At least 80 resolved residues on both sides excludes peptide-only partners;
    # this is an explicit operational size restriction, not proof of folding.
    parts = {str(c): arr[arr.chain_id == c] for c in set(arr.chain_id)}
    eligible = {c: a for c, a in parts.items() if len(residue_keys(a)) >= 80}
    if len(eligible) < 2:
        raise ValueError("fewer_than_two_80_residue_protein_chains")
    for chain in sorted(eligible, key=lambda c: (-len(residue_keys(eligible[c])), c)):
        target = eligible[chain]
        if len(residue_keys(target)) > 1500:
            continue
        others = [c for c in sorted(eligible) if c != chain]
        partners = arr[np.isin(arr.chain_id, others)]
        _, contacts = contact_tables(target, partners, "candidate", 5.)
        if contacts:
            touching = sorted({r[5] for r in contacts})
            return target, arr[np.isin(arr.chain_id, touching)], touching
    raise ValueError("no_eligible_resolved_interface")


def remap(target, instance):
    keys = residue_keys(target)
    lookup = {k: i+1 for i, k in enumerate(keys)}
    exported = target.copy()
    exported.chain_id[:] = "A"
    exported.res_id = np.array([lookup[key(target, i)] for i in range(len(target))])
    exported.ins_code[:] = ""
    entities = {key(target, i): str(target.label_entity_id[i]) for i in range(len(target))}
    mapping = [[instance, "A", lookup[k], "", k[3], k[0], k[1], k[2], entities[k]] for k in keys]
    return exported, mapping, lookup


def make_tables(target, partners, instance):
    exported, mapping, lookup = remap(target, instance)
    rs, cs = contact_tables(target, partners, instance, 5.)
    labels = []
    for r in rs:
        source_key = (r[1], r[2], r[3], r[4])
        labels.append([instance, "A", lookup[source_key], "", r[4], r[5],
                       "observed_ppi_contact" if r[5] else "no_observed_contact_to_selected_partners",
                       r[7], r[8], r[1], r[2], r[3]])
    contacts = []
    for r in cs:
        k = (r[1], r[2], r[3], r[4])
        contacts.append([instance, "A", lookup[k], "", r[4], *r[5:]])
    return exported, mapping, sorted(labels, key=lambda r: r[2]), sorted(contacts, key=lambda r: (r[2], r[5], r[6], r[7]))


def canonical_sequence(target):
    # Includes modified residue identity; never collapses insertion-coded residues.
    return "-".join(k[3] for k in residue_keys(target))


def verify_tables(out, name, header, rows):
    with (out / name).open() as fh:
        stored = list(csv.reader(fh, delimiter="\t"))
    if stored != [header] + [[str(v) for v in r] for r in rows]:
        raise ValueError("table mismatch: " + name)


def validate(out):
    manifest = json.loads((out / "build_manifest.json").read_text())
    records = json.loads((out / "provenance.json").read_text())
    if not manifest["complete"] or len(records) != manifest["requested"]:
        raise ValueError("incomplete build")
    mappings = []; labels = []; contacts = []; pdbs = set(); sequences = set(); instances = set(); mapped_insertions = 0
    for record in records:
        inst = record["instance"]; source = out / record["source_complex"]
        dest = out / "structures" / (inst + ".cif")
        if digest(source) != record["source_sha256"] or digest(dest) != record["target_sha256"]:
            raise ValueError("checksum mismatch")
        arr = protein_arrays(source)
        target = arr[arr.chain_id == record["source_target_chain"]]
        partners = arr[np.isin(arr.chain_id, record["source_partner_chains"])]
        if len(residue_keys(target)) < 80 or not len(partners):
            raise ValueError("invalid chain sizes")
        if set(target.chain_id) & set(partners.chain_id):
            raise ValueError("target self contact")
        for c in set(partners.chain_id):
            if len(residue_keys(partners[partners.chain_id == c])) < 80:
                raise ValueError("peptide-only partner")
        expected, mapping, rs, cs = make_tables(target, partners, inst)
        exported = load(dest)
        if not np.array_equal(expected.coord, exported.coord):
            raise ValueError("target coordinates changed")
        for name in ["chain_id", "res_id", "ins_code", "res_name", "atom_name", "element", "hetero"]:
            if not np.array_equal(getattr(expected, name), getattr(exported, name)):
                raise ValueError("export identity mismatch")
        # Independent dense-distance calculation, not the KD-tree label method.
        atom_min = []
        for start in range(0, len(target), 256):
            atom_min.extend(cdist(target.coord[start:start+256].astype(float), partners.coord.astype(float)).min(axis=1))
        positive = {key(target, i) for i, d in enumerate(atom_min) if d <= 5.}
        labelled = {(r[9], r[10], r[11], r[4]) for r in rs if r[5]}
        if positive != labelled or not positive:
            raise ValueError("independent contact audit failed")
        seq = canonical_sequence(target)
        if inst in instances or record["pdb_id"] in pdbs or seq in sequences:
            raise ValueError("duplicate target")
        pdbs.add(record["pdb_id"]); sequences.add(seq); instances.add(inst)
        mapped_insertions += sum(bool(r[7]) for r in mapping)
        mappings.extend(mapping); labels.extend(rs); contacts.extend(cs)
    verify_tables(out, "residue_mapping.tsv", MAP_HEADER, mappings)
    verify_tables(out, "residue_labels.tsv", LABEL_HEADER, labels)
    verify_tables(out, "annotated_residues.tsv", LABEL_HEADER, [r for r in labels if r[5]])
    verify_tables(out, "contacts.tsv", CONTACT_HEADER, contacts)
    for name in ["panel.tsv", "structures.tsv"]:
        with (out / name).open() as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        if len(rows) != len(instances) or {r["instance"] for r in rows} != instances:
            raise ValueError("membership mismatch")
    if {p.stem for p in (out / "structures").glob("*.cif")} != instances:
        raise ValueError("export file mismatch")
    result = {"valid": True, "target_proteins": len(instances), "unique_pdb_entries": len(pdbs),
              "resolved_target_residues": len(labels), "contact_residues": sum(r[5] for r in labels),
              "residue_contact_pairs": len(contacts), "insertion_coded_source_residues_preserved": mapped_insertions,
              "partner_atoms_exported": 0, "independent_dense_distance_audit": True,
              "interpretation": "experimental assembly protein contacts, not independent curated biological validation"}
    save_json(out / "validation.json", result)
    (out / "REPORT.md").write_text(
        "# true_annotated_ppi\n\n" +
        f"Validated {len(instances)} target-only proteins from {len(pdbs)} distinct source PDB entries.\n" +
        f"{result['contact_residues']} contact residues; {len(labels)} resolved target residues.\n\n" +
        "Labels are observed 5.0 A heavy-atom contacts to specified protein partners in stored biological assembly 1.\n" +
        "Target and partner chains each have >=80 resolved amino acids; targets have <=1500.\n" +
        "Antibody entries excluded using local SAbDab2 membership and entity-name screening.\n" +
        "Target-only exports use chain A and sequential residue numbering; residue_mapping.tsv retains all source insertion codes.\n" +
        "No partner atoms are exported. Source complexes remain separate.\n\n" +
        "Not curated energetic hotspots or independent validation of every biological interface. No reviewed UniProt requirement.\n" +
        "No family clustering; homomeric contacts can be present. Unresolved/other-partner contacts remain unassessed.\n")
    print(json.dumps(result, indent=2), flush=True)
    return result


def build(args):
    out = args.out.resolve(); source = args.source.resolve()
    if (out / "build_manifest.json").exists():
        raise ValueError("finished output exists; validate it or choose another output")
    for name in ["structures", "source_complexes"]:
        (out / name).mkdir(parents=True, exist_ok=True)
    metadata = source / "structures.tsv"
    with metadata.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    with args.sabdab.open() as fh:
        antibody_pdbs = {r["PDB_ID"][-4:].upper() for r in csv.DictReader(fh)}
    rows.sort(key=lambda r: hashlib.sha256(("true-ppi-150:" + r["pdb_id"]).encode()).hexdigest())
    records = []; labels = []; contacts = []; mappings = []; structures = []; panel = []; sequences = set(); excluded = Counter()
    for row in rows:
        pdb = row["pdb_id"]
        if pdb.upper() in antibody_pdbs:
            excluded["antibody_pdb_in_sabdab"] += 1; continue
        if row["structure_source"] != "biological assembly 1":
            excluded["not_biological_assembly_1"] += 1; continue
        src = source / "structures" / row["structure_file"]
        try:
            arr = protein_arrays(src)
            target, partners, partner_chains = select_target(arr)
            sequence = canonical_sequence(target)
            if sequence in sequences:
                excluded["exact_target_sequence_duplicate"] += 1; continue
            inst = pdb + "_target"
            exported, mapping, rs, cs = make_tables(target, partners, inst)
        except (ValueError, KeyError, OSError) as exc:
            excluded[str(exc)] += 1; continue
        dest = out / "structures" / (inst + ".cif")
        cif = pdbx.CIFFile(); pdbx.set_structure(cif, exported); cif.write(dest)
        copy = out / "source_complexes" / src.name; shutil.copyfile(src, copy)
        source_chain = str(target.chain_id[0])
        records.append({"instance": inst, "pdb_id": pdb, "source_metadata": row,
                        "source_target_chain": source_chain, "source_partner_chains": partner_chains,
                        "source_target_entity": str(target.label_entity_id[0]),
                        "source_complex": str(copy.relative_to(out)), "original_source": str(src),
                        "source_sha256": digest(copy), "target_sha256": digest(dest)})
        labels.extend(rs); contacts.extend(cs); mappings.extend(mapping); sequences.add(sequence)
        structures.append([inst, pdb, str(dest), source_chain, "|".join(partner_chains), len(rs), sum(r[5] for r in rs), "biological assembly 1"])
        panel.append([inst, str(dest), "true_annotated_ppi"])
        print("selected", len(records), "of", args.count, inst, flush=True)
        if len(records) == args.count:
            break
    table(out / "panel.tsv", ["instance", "structure_file", "category"], panel)
    table(out / "structures.tsv", ["instance", "pdb_id", "structure_file", "source_target_chain", "source_partner_chains",
                                   "n_resolved_target_residues", "n_annotated", "structure_source"], structures)
    table(out / "residue_mapping.tsv", MAP_HEADER, mappings)
    table(out / "residue_labels.tsv", LABEL_HEADER, labels)
    table(out / "annotated_residues.tsv", LABEL_HEADER, [r for r in labels if r[5]])
    table(out / "contacts.tsv", CONTACT_HEADER, contacts)
    save_json(out / "provenance.json", records)
    save_json(out / "build_manifest.json", {
        "built_utc": datetime.now(timezone.utc).isoformat(), "requested": args.count, "selected": len(records),
        "complete": len(records) == args.count, "source_metadata": str(metadata), "source_metadata_sha256": digest(metadata),
        "source_build_manifest": json.loads((source / "build_manifest.json").read_text()),
        "antibody_exclusion_metadata": str(args.sabdab.resolve()), "antibody_exclusion_sha256": digest(args.sabdab),
        "contact_distance_angstrom": 5.0, "excluded": dict(excluded),
        "selection": "ID hash true-ppi-150 order; one target per entry; longest eligible contacting chain, ties source chain ID; exact sequence dedup",
        "size_policy": "target 80..1500 resolved amino acids; each retained partner >=80; all touching eligible partners retained",
        "antibody_policy": "exclude local SAbDab2 PDB members and antibody entity descriptions",
        "mapping": "target chain A, sequential 1..N, blank insertion codes; reversible source identity mapping",
        "limitations": ["assembly contacts not independent biological interaction curation", "source experimental quality inherited from old dataset",
                        "source manifest asymmetric-unit caveat is stale: actual files are assembly 1",
                        "no family-level deduplication", "antibody exclusion not guaranteed exhaustive beyond metadata",
                        "short proteins/peptide partners excluded by explicit >=80 residue rule", "unresolved residues unknown",
                        "source complexes include full assemblies; only structures/ is target-only"]})
    if len(records) != args.count:
        raise RuntimeError("quota shortfall; see manifest")
    validate(out)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, default=ROOT / "investigation_databases/protein_protein_interactions")
    p.add_argument("--out", type=Path, default=ROOT / "investigation_databases/true_annotated_ppi")
    p.add_argument("--sabdab", type=Path, default=Path.home() / "documents_local_only/master_thesis_data_files/splits_final/abag_split.csv")
    p.add_argument("--count", type=int, default=150)
    p.add_argument("--validate-only", action="store_true")
    args = p.parse_args()
    if args.count < 1:
        p.error("count must be positive")
    if args.validate_only:
        validate(args.out.resolve())
    else:
        build(args)


if __name__ == "__main__":
    main()