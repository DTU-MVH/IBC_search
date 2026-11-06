"""Reviewed zinc fingers with explicit, domain-contained UniProt Zn ligands.

No proximity labels or whole-domain positives. 'True annotated' means traceable
reviewed-source annotations, not experimental confirmation of every residue.
"""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import csv
import hashlib
import json
from pathlib import Path
import re
import urllib.parse

from build_reviewed_enzymes import (
    AA3, ROOT, fetch, model_for, position, save_json, table, validate_model,
)

QUERY = r"reviewed:true AND ft_zn_fing:* AND ft_binding:CHEBI\:29105"


def exact_domain(feature, length):
    loc = feature.get("location", {})
    a = loc.get("start", {}); b = loc.get("end", {})
    start = a.get("value"); end = b.get("value")
    if (a.get("modifier") != "EXACT" or b.get("modifier") != "EXACT"
            or not isinstance(start, int) or not isinstance(end, int)
            or not 1 <= start <= end <= length):
        return None
    return start, end


def select_labels(entry):
    if entry.get("entryType") != "UniProtKB reviewed (Swiss-Prot)":
        return [], "not_reviewed"
    sequence = entry["sequence"]["value"]
    domains = []
    for f in entry.get("features", []):
        if f["type"] == "Zinc finger":
            bounds = exact_domain(f, len(sequence))
            if bounds:
                domains.append((bounds, f))
    if not domains:
        return [], "no_exact_zinc_finger_domain"
    labels = []; seen = set()
    for feature in entry.get("features", []):
        if (feature["type"] != "Binding site"
                or feature.get("ligand", {}).get("id") != "ChEBI:CHEBI:29105"):
            continue
        p = position(feature, sequence)
        if p is None:
            continue
        overlapping = [(bounds, f) for bounds, f in domains if bounds[0] <= p <= bounds[1]]
        if not overlapping:
            continue
        if len(overlapping) != 1:
            return [], "ambiguous_overlapping_domains"
        if p in seen:
            return [], "duplicate_binding_residue"
        if sequence[p-1] not in AA3:
            return [], "nonstandard_binding_residue"
        seen.add(p)
        bounds, domain = overlapping[0]
        labels.append({"position": p, "residue": sequence[p-1],
                       "label_type": "zinc_finger_zinc_ligand", "feature_type": "Binding site",
                       "description": feature.get("description", ""),
                       "ligand": feature["ligand"], "evidence": feature.get("evidences", []),
                       "domain_start": bounds[0], "domain_end": bounds[1],
                       "domain_description": domain.get("description", ""),
                       "domain_evidence": domain.get("evidences", []),
                       "source_feature": feature, "source_domain": domain})
    return sorted(labels, key=lambda x: x["position"]), "eligible" if labels else "no_explicit_zinc_ligands_inside_exact_domains"


ANNOTATION_HEADER = ["instance", "uniprot_accession", "chain", "res_id", "res_name", "is_annotated",
                     "label_type", "site_ligand", "domain_start", "domain_end", "domain_description",
                     "evidence_json", "domain_evidence_json"]


def annotation_rows(record, chain):
    ac = record["accession"]
    return [[ac, ac, chain, l["position"], AA3[l["residue"]], 1, l["label_type"], "ZN",
             l["domain_start"], l["domain_end"], l["domain_description"],
             json.dumps(l["evidence"], separators=(",", ":")),
             json.dumps(l["domain_evidence"], separators=(",", ":"))] for l in record["labels"]]


def validate(out):
    manifest = json.loads((out / "build_manifest.json").read_text())
    records = json.loads((out / "provenance.json").read_text())
    if not manifest["complete"] or len(records) != manifest["requested"]:
        raise ValueError("incomplete build")
    seen = set(); expected = []; evidence = Counter(); domains = Counter()
    for record in records:
        ac = record["accession"]
        if ac in seen:
            raise ValueError("duplicate accession")
        seen.add(ac)
        labels, reason = select_labels(record["uniprot_entry"])
        if reason != "eligible" or labels != record["labels"]:
            raise ValueError("source labels do not reproduce")
        raw = (out / "structures" / (ac + ".pdb")).read_bytes()
        if hashlib.sha256(raw).hexdigest() != record["structure_sha256"]:
            raise ValueError("structure checksum mismatch")
        chain = validate_model(raw.decode(), record["uniprot_entry"]["sequence"]["value"], labels)
        expected.extend(annotation_rows(record, chain))
        for l in labels:
            if not l["evidence"]:
                evidence["no_feature_evidence"] += 1
            for e in l["evidence"]:
                evidence[e["evidenceCode"]] += 1
        for _, _, description in {(l["domain_start"], l["domain_end"], l["domain_description"]) for l in labels}:
            domains[description] += 1
    with (out / "annotated_residues.tsv").open() as fh:
        actual = list(csv.reader(fh, delimiter="\t"))
    if actual != [ANNOTATION_HEADER] + [[str(x) for x in row] for row in expected]:
        raise ValueError("annotation table differs from source")
    for name in ["panel.tsv", "structures.tsv"]:
        with (out / name).open() as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        if len(rows) != len(seen) or {r["instance"] for r in rows} != seen:
            raise ValueError("panel/structure membership mismatch")
    result = {"valid": True, "reviewed_proteins": len(seen), "annotated_residue_rows": len(expected),
              "binding_evidence_codes": dict(evidence), "labelled_domains": dict(domains),
              "interpretation": "explicit reviewed-source domain-contained Zn binding annotations; not experimental-only and not guaranteed complete coordination shells"}
    save_json(out / "validation.json", result)
    (out / "REPORT.md").write_text(
        "# true_annotated_zincfingers\n\n"
        f"Validated {len(seen)} unique reviewed UniProt proteins and {len(expected)} explicit Zn-binding residue labels.\n\n"
        "Only exact single-residue Zn binding features inside uniquely assigned, exact zinc-finger domains are labelled.\n"
        "No proximity-derived labels, generic ligand pockets, or whole-domain positives.\n\n"
        "Structures are AlphaFold DB predicted monomers, sequence/numbering checked against UniProt.\n"
        "Reviewed does not mean experimental-only. Domain completeness, metal geometry and family independence are not established.\n"
        "Unlabelled residues are unknown, not proven non-binding. See validation.json and provenance.json for evidence.\n")
    print(json.dumps(result, indent=2), flush=True)
    return result


def build(args):
    out = args.out.resolve(); out.mkdir(parents=True, exist_ok=True)
    if (out / "build_manifest.json").exists():
        raise ValueError("finished output exists; use --validate-only or choose a new directory")
    (out / "structures").mkdir(exist_ok=True)
    url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode({
        "query": QUERY, "format": "json", "size": 500, "sort": "accession asc"})
    selected = []; excluded = Counter(); pages = []; page = 0
    while url and len(selected) < args.count:
        cache = out / "cache/uniprot" / f"page_{page}.json"
        raw, headers = fetch(url, cache)
        hp = cache.with_suffix(".headers.json")
        if headers:
            save_json(hp, headers)
        else:
            headers = json.loads(hp.read_text())
        pages.append({"url": url, "sha256": hashlib.sha256(raw).hexdigest(),
                      "release": headers.get("x-uniprot-release"), "total_results": headers.get("x-total-results")})
        candidates = []
        for entry in json.loads(raw)["results"]:
            labels, reason = select_labels(entry)
            if labels:
                candidates.append((entry, labels, out))
            else:
                excluded[reason] += 1
        print("page", page, "eligible", len(candidates), flush=True)
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            for start in range(0, len(candidates), args.jobs):
                for result, reason in pool.map(model_for, candidates[start:start+args.jobs]):
                    if result is None:
                        excluded[reason] += 1
                    elif len(selected) < args.count:
                        selected.append(result)
                save_json(out / "build_progress.json", {"selected": len(selected), "requested": args.count,
                                                         "excluded": dict(excluded), "pages": pages})
                print("selected", len(selected), "of", args.count, flush=True)
                if len(selected) >= args.count:
                    break
        nxt = re.search(r'<([^>]+)>; rel="next"', headers.get("link", ""))
        url = nxt.group(1) if nxt else None
        page += 1
    records = []; annotations = []; structures = []; panel = []; fasta = []
    for item in selected:
        entry = item["entry"]; ac = entry["primaryAccession"]
        dest = out / "structures" / (ac + ".pdb")
        dest.write_bytes(Path(item["cached_structure"]).read_bytes())
        record = {"accession": ac, "uniprot_entry": entry, "labels": item["labels"],
                  "alphafold_metadata": item["model"], "structure_sha256": item["sha256"]}
        records.append(record); annotations.extend(annotation_rows(record, item["chain"]))
        panel.append([ac, str(dest), "true_annotated_zincfingers"])
        structures.append([ac, ac, str(dest), len(entry["sequence"]["value"]), len(item["labels"]),
                           entry["organism"]["scientificName"], "AlphaFold DB predicted monomer",
                           item["model"]["pdbUrl"], item["sha256"]])
        fasta.append(f">{ac}\n{entry['sequence']['value']}\n")
    table(out / "panel.tsv", ["instance", "structure_file", "category"], panel)
    table(out / "structures.tsv", ["instance", "uniprot_accession", "structure_file", "sequence_length",
                                   "n_annotated", "organism", "structure_source", "structure_url", "sha256"], structures)
    table(out / "annotated_residues.tsv", ANNOTATION_HEADER, annotations)
    (out / "sequences.fasta").write_text("".join(fasta))
    save_json(out / "provenance.json", records)
    save_json(out / "build_manifest.json", {
        "built_utc": datetime.now(timezone.utc).isoformat(), "requested": args.count,
        "selected": len(selected), "complete": len(selected) == args.count, "query": QUERY,
        "pages": pages, "excluded": dict(excluded), "reviewed_only": True,
        "selection": "accession ascending; first eligible exact full-length AFDB models; no sequence clustering",
        "label_rule": "explicit exact Zn Binding site uniquely inside exact Zinc finger domain; no proximity labels",
        "evidence_policy": "preserve UniProt residue and domain evidence; not experimental-only",
        "structure_policy": "predicted full-length canonical monomers; entire sequence and residue numbering validated",
        "limitations": ["No complete coordination-shell requirement", "No pLDDT cutoff", "Not family independent",
                        "Domain containment is the operational assignment, not an independent structural validation",
                        "Unlabelled residues remain unknown"]})
    if len(selected) != args.count:
        raise RuntimeError("quota shortfall; no relaxed labels or unreviewed padding used")
    validate(out)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, default=ROOT / "investigation_databases/true_annotated_zincfingers")
    p.add_argument("--count", type=int, default=100)
    p.add_argument("--jobs", type=int, default=2)
    p.add_argument("--validate-only", action="store_true")
    args = p.parse_args()
    if args.count < 1 or args.jobs < 1:
        p.error("positive count and jobs required")
    if args.validate_only:
        validate(args.out.resolve())
    else:
        build(args)


if __name__ == "__main__":
    main()