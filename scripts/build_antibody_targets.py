"""Target-only experimental antibody interfaces from the local SAbDab2 release.

Labels are resolved direct protein contacts, not curated energetic hotspots or
all possible epitopes. Source complexes are retained outside structures/.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import csv
import hashlib
import json
from pathlib import Path
import shutil
import warnings

import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
from scipy.spatial import cKDTree

from build_reviewed_enzymes import ROOT, save_json, table

RESIDUE_HEADER = ["instance", "chain", "res_id", "ins_code", "res_name", "is_annotated",
                  "label_type", "min_antibody_distance_angstrom", "contacting_antibody_chains"]
CONTACT_HEADER = ["instance", "target_chain", "target_res_id", "target_ins_code", "target_res_name",
                  "antibody_chain", "antibody_res_id", "antibody_ins_code", "antibody_res_name",
                  "min_distance_angstrom"]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        return pdbx.get_structure(pdbx.CIFFile.read(path), model=1)


def key(array, i):
    return str(array.chain_id[i]), int(array.res_id[i]), str(array.ins_code[i]).strip(), str(array.res_name[i])


def chain_selection(row):
    chains = row.get("agchains", "").split("/")
    types = row.get("agtypes", "").split("/")
    clusters = row.get("agclusters", "").split("/")
    if len(chains) != len(types) or len(chains) != len(clusters):
        raise ValueError("ambiguous antigen chain/type/cluster assignment")
    target = [c for c, t in zip(chains, types) if t == "PROTEIN"]
    cluster = {v for v, t in zip(clusters, types) if t == "PROTEIN"}
    antibody = [row.get(k, "").strip() for k in ["Hchain", "Lchain"]]
    antibody = [c for c in antibody if c and c.lower() not in {"nan", "none", "na"}]
    if not target or not antibody or not all(cluster) or set(target) & set(antibody):
        raise ValueError("invalid target/antibody assignments")
    if len(set(target)) != len(target) or len(set(antibody)) != len(antibody):
        raise ValueError("duplicate chain assignments")
    return target, antibody, cluster


def extract(array, row):
    target, antibody, _ = chain_selection(row)
    # Explicit amino-acid filter: no water, glycans, nucleic acid or other ligands.
    heavy = ~np.isin(array.element, ["H", "D"])
    protein = struc.filter_amino_acids(array)
    ag = array[heavy & protein & np.isin(array.chain_id, target)]
    ab = array[heavy & protein & np.isin(array.chain_id, antibody)]
    if set(ag.chain_id) != set(target) or set(ab.chain_id) != set(antibody):
        raise ValueError("missing selected protein chain atoms")
    if any(len(c) != 1 for c in target):
        raise ValueError("target chain not FoldDisco expressible")
    for arr in [ag, ab]:
        ids = [(key(arr, i), str(arr.atom_name[i])) for i in range(len(arr))]
        if len(set(ids)) != len(ids):
            raise ValueError("duplicate atoms/alternate conformers after parsing")
    return ag, ab


def contact_tables(ag, ab, instance, cutoff):
    if cutoff <= 0 or not np.isfinite(cutoff):
        raise ValueError("positive finite contact cutoff required")
    tree = cKDTree(ab.coord)
    distances = tree.query(ag.coord)[0]
    minima = {}; pairs = {}
    for i, partners in enumerate(tree.query_ball_point(ag.coord, cutoff)):
        target = key(ag, i)
        minima[target] = min(minima.get(target, float("inf")), float(distances[i]))
        for j in partners:
            pair = target, key(ab, j)
            distance = float(np.linalg.norm(ag.coord[i].astype(float) - ab.coord[j].astype(float)))
            pairs[pair] = min(pairs.get(pair, float("inf")), distance)
    partners = {}
    for target, antibody in pairs:
        partners.setdefault(target, set()).add(antibody[0])
    residues = []
    for target in sorted(minima):
        positive = target in partners
        residues.append([instance, *target, int(positive),
                         "observed_antibody_contact" if positive else "no_observed_contact_in_selected_complex",
                         f"{minima[target]:.6f}", "|".join(sorted(partners.get(target, set())))])
    contacts = [[instance, *target, *antibody, f"{distance:.6f}"]
                for (target, antibody), distance in sorted(pairs.items())]
    return residues, contacts


def fingerprint(ag):
    # Invariant to chain names/numbering/order; conservative exact sequence dedup.
    sequences = []
    for chain in sorted(set(ag.chain_id)):
        arr = ag[ag.chain_id == chain]; seen = set(); names = []
        for i in range(len(arr)):
            k = key(arr, i)
            if k not in seen:
                names.append(k[3]); seen.add(k)
        sequences.append("-".join(names))
    return hashlib.sha256("|".join(sorted(sequences)).encode()).hexdigest()


def validate(out):
    manifest = json.loads((out / "build_manifest.json").read_text())
    records = json.loads((out / "provenance.json").read_text())
    if not manifest["complete"] or len(records) != manifest["requested"]:
        raise ValueError("incomplete dataset")
    expected_residues = []; expected_contacts = []; clusters = set(); pdbs = set(); sequences = set(); instances = set()
    for record in records:
        row = record["metadata"]; inst = row["INSTANCE"]
        if inst in instances:
            raise ValueError("duplicate instance")
        instances.add(inst)
        source = out / record["source_file"]; dest = out / "structures" / (inst + ".cif")
        if digest(source) != record["source_sha256"] or digest(dest) != record["target_sha256"]:
            raise ValueError("checksum mismatch")
        ag, ab = extract(load(source), row); exported = load(dest)
        if len(ag) != len(exported) or not np.array_equal(ag.coord, exported.coord):
            raise ValueError("target coordinates changed")
        for annotation in ["chain_id", "res_id", "ins_code", "res_name", "atom_name", "element", "hetero"]:
            if not np.array_equal(getattr(ag, annotation), getattr(exported, annotation)):
                raise ValueError("target identities changed")
        target, antibody, group = chain_selection(row)
        if set(exported.chain_id) & set(antibody) or set(exported.chain_id) != set(target):
            raise ValueError("antibody/unselected chain in export")
        fp = fingerprint(ag)
        if clusters & group or row["PDB_ID"] in pdbs or fp in sequences:
            raise ValueError("target redundancy")
        clusters.update(group); pdbs.add(row["PDB_ID"]); sequences.add(fp)
        rs, cs = contact_tables(ag, ab, inst, manifest["contact_distance_angstrom"])
        if not cs:
            raise ValueError("no observed interface")
        expected_residues.extend(rs); expected_contacts.extend(cs)
    for name, header, values in [
        ("residue_labels.tsv", RESIDUE_HEADER, expected_residues),
        ("annotated_residues.tsv", RESIDUE_HEADER, [r for r in expected_residues if r[5]]),
        ("contacts.tsv", CONTACT_HEADER, expected_contacts),
    ]:
        with (out / name).open() as fh:
            actual = list(csv.reader(fh, delimiter="\t"))
        if actual != [header] + [[str(v) for v in r] for r in values]:
            raise ValueError("labels do not reproduce: " + name)
    for name in ["panel.tsv", "structures.tsv"]:
        with (out / name).open() as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        if len(rows) != len(instances) or {r["instance"] for r in rows} != instances:
            raise ValueError("membership mismatch: " + name)
    if {p.stem for p in (out / "structures").glob("*.cif")} != instances:
        raise ValueError("export directory membership mismatch")
    result = {"valid": True, "target_structures": len(records), "resolved_target_residues": len(expected_residues),
              "contact_residues": sum(r[5] for r in expected_residues), "residue_contact_pairs": len(expected_contacts),
              "disjoint_antigen_clusters": len(clusters), "unique_pdb_entries": len(pdbs),
              "antibody_atoms_exported": 0, "labels_reproduced": True}
    save_json(out / "validation.json", result)
    (out / "REPORT.md").write_text(
        "# true_annotated_antibody_targets\n\n" +
        f"Validated {len(records)} experimental target-only structures; {result['contact_residues']} resolved contact residues.\n\n" +
        "Labels: any target amino-acid heavy atom within 5.0 A of the selected antibody variable-region heavy atoms.\n" +
        "Structures contain target amino-acid heavy atoms only; source complexes are stored separately.\n" +
        "Annotations are reproducible structural contacts, NOT curated energetic hotspots or every possible epitope.\n" +
        "Unresolved residues, water-mediated and glycan contacts are not assessed. Source constructs may be truncated.\n" +
        "One selected antibody instance per PDB/antigen cluster set; labels are specific to that partner, not all antibodies in the original deposition.\n")
    print(json.dumps(result, indent=2), flush=True)
    return result


def build(args):
    out = args.out.resolve(); source = args.source.resolve()
    if (out / "build_manifest.json").exists():
        raise ValueError("finished output exists; use --validate-only or a new output directory")
    for name in ["structures", "source_complexes"]:
        (out / name).mkdir(parents=True, exist_ok=True)
    csv_path = source / "abag_split.csv"
    with csv_path.open() as fh:
        rows = list(csv.DictReader(fh))
    candidates = []; excluded = Counter()
    for row in rows:
        try:
            resolution = float(row["resolution"])
            if row["method"] != "XRAY" or not 0 < resolution <= 3.0 or row["holo"].lower() != "true":
                continue
            chain_selection(row)
            candidates.append(row)
        except (ValueError, KeyError):
            excluded["invalid_metadata_or_nonprotein_target"] += 1
    candidates.sort(key=lambda r: (float(r["resolution"]), r["INSTANCE"]))
    records = []; residue_rows = []; contact_rows = []; panel = []; structures = []
    clusters = set(); pdbs = set(); sequences = set()
    for row in candidates:
        target, antibody, group = chain_selection(row)
        if clusters & group or row["PDB_ID"] in pdbs:
            excluded["pdb_or_cluster_redundancy"] += 1; continue
        inst = row["INSTANCE"]; src = source / (inst + ".cif")
        if not src.exists():
            excluded["missing_source"] += 1; continue
        try:
            ag, ab = extract(load(src), row)
            fp = fingerprint(ag)
            if fp in sequences:
                excluded["exact_target_sequence_duplicate"] += 1; continue
            rs, cs = contact_tables(ag, ab, inst, args.contact)
            if not cs:
                excluded["no_resolved_interface"] += 1; continue
        except (ValueError, KeyError) as exc:
            excluded[str(exc)] += 1; continue
        dest = out / "structures" / (inst + ".cif")
        cif = pdbx.CIFFile(); pdbx.set_structure(cif, ag); cif.write(dest)
        src_copy = out / "source_complexes" / src.name; shutil.copyfile(src, src_copy)
        record = {"metadata": row, "source_file": str(src_copy.relative_to(out)),
                  "original_source_file": str(src), "source_sha256": digest(src_copy),
                  "target_sha256": digest(dest), "target_sequence_fingerprint": fp,
                  "interpretation": "contacts to this SAbDab2 selected antibody instance; variable-region cropped complex"}
        records.append(record); residue_rows.extend(rs); contact_rows.extend(cs)
        panel.append([inst, str(dest), "true_annotated_antibody_targets"])
        structures.append([inst, row["PDB_ID"], str(dest), "/".join(target), "/".join(antibody),
                           row["resolution"], len(rs), sum(r[5] for r in rs), "/".join(sorted(group)),
                           row.get("ab_ag_split", ""), row.get("agresolvedseqs", ""), row.get("agexpectedseqs", "")])
        clusters.update(group); pdbs.add(row["PDB_ID"]); sequences.add(fp)
        print("selected", len(records), "of", args.count, inst, flush=True)
        if len(records) == args.count:
            break
    table(out / "panel.tsv", ["instance", "structure_file", "category"], panel)
    table(out / "structures.tsv", ["instance", "pdb_id", "structure_file", "target_chains", "antibody_chains",
                                  "resolution_angstrom", "n_resolved_target_residues", "n_annotated", "antigen_clusters",
                                  "source_split", "source_resolved_sequences", "source_expected_sequences"], structures)
    table(out / "residue_labels.tsv", RESIDUE_HEADER, residue_rows)
    table(out / "annotated_residues.tsv", RESIDUE_HEADER, [r for r in residue_rows if r[5]])
    table(out / "contacts.tsv", CONTACT_HEADER, contact_rows)
    save_json(out / "provenance.json", records)
    save_json(out / "build_manifest.json", {
        "built_utc": datetime.now(timezone.utc).isoformat(), "requested": args.count, "selected": len(records),
        "complete": len(records) == args.count, "source": "local SAbDab2 experimental cropped complexes",
        "source_metadata": str(csv_path), "source_metadata_sha256": digest(csv_path),
        "contact_distance_angstrom": args.contact, "source_rows": len(rows), "candidate_rows": len(candidates),
        "selection": "XRAY <=3.0 A; holo; resolution ascending then instance; no shared antigen cluster, PDB entry or exact target sequence",
        "excluded": dict(excluded), "source_complex_directory": "source_complexes", "target_only_directory": "structures",
        "label_definition": "resolved target amino-acid heavy atom within cutoff of selected antibody variable-region amino-acid heavy atom",
        "limitations": ["structural contact labels, not curated energetic hotspots", "no glycan or water-mediated labels",
                        "unresolved residues unknown; source resolved/expected sequences retained",
                        "source complexes contain selected variable regions, not necessarily all antibodies in the PDB",
                        "cluster IDs inherited from local SAbDab2 release; no new family clustering",
                        "source coordinate numbering preserved, not assumed original PDB/UniProt numbering",
                        "experimental bound conformation; target export is not an apo measurement"]})
    if len(records) != args.count:
        raise RuntimeError("quota shortfall; see manifest")
    validate(out)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, default=Path.home() / "documents_local_only/master_thesis_data_files/splits_final")
    p.add_argument("--out", type=Path, default=ROOT / "investigation_databases/true_annotated_antibody_targets")
    p.add_argument("--count", type=int, default=70)
    p.add_argument("--contact", type=float, default=5.0)
    p.add_argument("--validate-only", action="store_true")
    args = p.parse_args()
    if args.count < 1 or args.contact != 5.0:
        p.error("positive count and agreed 5.0 A contact definition required")
    if args.validate_only:
        validate(args.out.resolve())
    else:
        build(args)


if __name__ == "__main__":
    main()