"""Re-validate a finished reviewed-enzyme build against cached source features."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path

from build_reviewed_enzymes import QUERIES, select_labels, validate_model, save_json


def validate(out):
    manifest = json.loads((out / "build_manifest.json").read_text())
    if not manifest["complete"]:
        raise ValueError("build is incomplete")
    summary = {}; all_accessions = set()
    for category in QUERIES:
        directory = out / category
        provenance = json.loads((directory / "provenance.json").read_text())
        if len(provenance) != manifest["requested_per_class"]:
            raise ValueError("incorrect class count")
        with (directory / "annotated_residues.tsv").open() as fh:
            annotations = list(csv.DictReader(fh, delimiter="\t"))
        expected = Counter(); evidence = Counter()
        for record in provenance:
            ac = record["accession"]
            if ac in all_accessions:
                raise ValueError("duplicate accession across/within classes")
            all_accessions.add(ac)
            labels, reason = select_labels(record["uniprot_entry"], category)
            if reason != "eligible" or labels != record["labels"]:
                raise ValueError("labels do not reproduce from reviewed source features")
            raw = (directory / "structures" / (ac + ".pdb")).read_bytes()
            if hashlib.sha256(raw).hexdigest() != record["structure_sha256"]:
                raise ValueError("structure checksum changed")
            chain = validate_model(raw.decode(), record["uniprot_entry"]["sequence"]["value"], labels)
            for label in labels:
                expected[(ac, chain, str(label["position"]), label["label_type"])] += 1
                if not label["evidence"]:
                    evidence["no_feature_evidence"] += 1
                for ev in label["evidence"]:
                    evidence[ev["evidenceCode"]] += 1
        observed = Counter((r["instance"], r["chain"], r["res_id"], r["label_type"]) for r in annotations)
        if observed != expected:
            raise ValueError("annotation table differs from source features")
        summary[category] = {"proteins": len(provenance), "annotations": len(annotations),
                             "label_types": dict(Counter(r["label_type"] for r in annotations)),
                             "evidence_codes": dict(evidence)}
    result = {"valid": True, "unique_reviewed_accessions": len(all_accessions), "classes": summary}
    save_json(out / "validation.json", result)
    print(json.dumps(result, indent=2))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    validate(parser.parse_args().out.resolve())