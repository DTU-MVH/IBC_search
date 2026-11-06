"""Offline tests for explicit, reviewed, domain-contained zinc labels."""
from copy import deepcopy
from pathlib import Path
import sys
import hashlib
import json

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from build_true_annotated_zincfingers import (
    select_labels, exact_domain, validate, table, annotation_rows, ANNOTATION_HEADER,
)


def feature(kind, start, end=None, **kwargs):
    return {"type": kind, "location": {"start": {"value": start, "modifier": "EXACT"},
            "end": {"value": start if end is None else end, "modifier": "EXACT"}}, **kwargs}


def example():
    return {"entryType": "UniProtKB reviewed (Swiss-Prot)", "sequence": {"value": "CACHHCAC"},
            "features": [feature("Zinc finger", 1, 5, description="C2H2-type"),
                         feature("Binding site", 1, ligand={"id": "ChEBI:CHEBI:29105", "name": "Zn(2+)"},
                                 evidences=[{"evidenceCode": "ECO:0000250"}]),
                         feature("Binding site", 8, ligand={"id": "ChEBI:CHEBI:29105"}),
                         feature("Binding site", 4, ligand={"id": "ChEBI:CHEBI:29108", "name": "Ca(2+)"})]}


def test_only_explicit_zinc_inside_domain():
    labels, reason = select_labels(example())
    assert reason == "eligible"
    assert [x["position"] for x in labels] == [1]
    assert labels[0]["domain_description"] == "C2H2-type"
    assert labels[0]["evidence"] == [{"evidenceCode": "ECO:0000250"}]
    assert labels[0]["source_feature"] == example()["features"][1]


def test_reviewed_and_exact_domain_required():
    e = example(); e["entryType"] = "UniProtKB unreviewed (TrEMBL)"
    assert select_labels(e)[1] == "not_reviewed"
    e = example(); e["features"][0]["location"]["start"]["modifier"] = "UNCERTAIN"
    assert select_labels(e)[1] == "no_exact_zinc_finger_domain"
    assert exact_domain(feature("Zinc finger", 0, 5), 8) is None


def test_range_and_uncertain_binding_not_used():
    e = example(); e["features"][1]["location"]["end"]["value"] = 3
    assert not select_labels(e)[0]
    e = example(); e["features"][1]["location"]["start"]["modifier"] = "UNCERTAIN"
    assert not select_labels(e)[0]


def test_overlapping_domains_rejected():
    e = example(); e["features"].append(feature("Zinc finger", 1, 8))
    assert select_labels(e)[1] == "ambiguous_overlapping_domains"


def test_duplicates_rejected():
    e = example(); e["features"].append(deepcopy(e["features"][1]))
    assert select_labels(e)[1] == "duplicate_binding_residue"


def test_domain_without_binding_is_not_annotated():
    e = example(); e["features"] = e["features"][:1]
    assert not select_labels(e)[0]


def test_finished_dataset_audit_detects_label_tampering(tmp_path):
    e = example(); e["sequence"]["value"] = "CHS"
    e["features"] = [feature("Zinc finger", 1, 3, description="test"),
                     feature("Binding site", 1, ligand={"id": "ChEBI:CHEBI:29105"})]
    labels, _ = select_labels(e)
    pdb = ("ATOM      1  CA  CYS A   1       0.000   0.000   0.000  1.00 90.00           C  \n"
           "ATOM      2  CA  HIS A   2       3.800   0.000   0.000  1.00 90.00           C  \n"
           "ATOM      3  CA  SER A   3       7.600   0.000   0.000  1.00 90.00           C  \nEND\n")
    (tmp_path / "structures").mkdir()
    (tmp_path / "structures/TEST.pdb").write_text(pdb)
    record = {"accession": "TEST", "uniprot_entry": e, "labels": labels,
              "structure_sha256": hashlib.sha256(pdb.encode()).hexdigest()}
    (tmp_path / "provenance.json").write_text(json.dumps([record]))
    (tmp_path / "build_manifest.json").write_text(json.dumps({"complete": True, "requested": 1}))
    for name in ["panel.tsv", "structures.tsv"]:
        table(tmp_path / name, ["instance"], [["TEST"]])
    table(tmp_path / "annotated_residues.tsv", ANNOTATION_HEADER, annotation_rows(record, "A"))
    assert validate(tmp_path)["reviewed_proteins"] == 1
    table(tmp_path / "annotated_residues.tsv", ANNOTATION_HEADER, [])
    with pytest.raises(ValueError, match="annotation table differs"):
        validate(tmp_path)