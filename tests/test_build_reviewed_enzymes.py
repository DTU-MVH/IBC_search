"""Offline tests for the restricted enzyme annotation policy."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from build_reviewed_enzymes import select_labels, validate_model


def feature(n, kind="Active site", description="", ligand=None, end=None):
    return {"type": kind, "description": description,
            "location": {"start": {"value": n, "modifier": "EXACT"},
                         "end": {"value": n if end is None else end, "modifier": "EXACT"}},
            "ligand": ligand or {}, "evidences": [{"evidenceCode": "ECO:0000250"}]}


def entry(features, sequence="HDS"):
    return {"entryType": "UniProtKB reviewed (Swiss-Prot)",
            "sequence": {"value": sequence}, "features": features}


def test_serine_requires_one_complete_explicit_triad():
    fs = [feature(n, description="Charge relay system") for n in (1, 2, 3)]
    labels, reason = select_labels(entry(fs), "serine_proteases")
    assert reason == "eligible"
    assert [x["position"] for x in labels] == [1, 2, 3]
    assert all(x["label_type"] == "catalytic_triad" for x in labels)
    assert not select_labels(entry(fs[:2]), "serine_proteases")[0]
    assert not select_labels(entry(fs, "AAA"), "serine_proteases")[0]
    assert not select_labels(entry(fs + [feature(1)]), "serine_proteases")[0]


def test_nitrogenase_ignores_atp_and_separates_cofactors():
    fs = [feature(1, "Binding site", ligand={"name": "ATP"}),
          feature(2, "Binding site", ligand={"name": "[7Fe-Mo-9S-C-homocitryl] cluster"}),
          feature(3, "Binding site", ligand={"name": "[4Fe-4S] cluster"})]
    labels, _ = select_labels(entry(fs), "nitrogenases")
    assert [x["position"] for x in labels] == [2, 3]
    assert [x["label_type"] for x in labels] == ["catalytic_cofactor_ligand", "fe_protein_electron_transfer_ligand"]


def test_metallo_excludes_structural_metal_and_general_binding():
    fs = [feature(1), feature(2, "Binding site", ligand={"name": "Zn(2+)", "note": "catalytic"}),
          feature(3, "Binding site", ligand={"name": "Zn(2+)", "note": "structural"})]
    labels, _ = select_labels(entry(fs), "metalloproteases")
    assert [x["position"] for x in labels] == [1, 2]
    assert not select_labels(entry(fs[1:]), "metalloproteases")[0]


def test_rejects_unreviewed_and_ranged_features():
    e = entry([feature(1)])
    e["entryType"] = "UniProtKB unreviewed (TrEMBL)"
    assert not select_labels(e, "metalloproteases")[0]
    assert not select_labels(entry([feature(1, end=2)]), "metalloproteases")[0]


def test_model_mapping_requires_exact_sequence():
    pdb = ("ATOM      1  CA  HIS A   1       0.000   0.000   0.000  1.00 90.00           C  \n"
           "ATOM      2  CA  ASP A   2       3.800   0.000   0.000  1.00 90.00           C  \n"
           "ATOM      3  CA  SER A   3       7.600   0.000   0.000  1.00 90.00           C  \nEND\n")
    assert validate_model(pdb, "HDS", [{"position": 3, "residue": "S"}]) == "A"
    with pytest.raises(ValueError, match="sequence mismatch"):
        validate_model(pdb, "HAS", [])
    with pytest.raises(ValueError, match="coverage mismatch"):
        validate_model(pdb, "HD", [])