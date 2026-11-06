"""Offline checks for target-site identity, unknowns and matched controls."""

import json
from pathlib import Path
import sys

import biotite.structure as struc
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from validate_target_sites import (  # noqa: E402
    assembly_interfaces, load_mcsa, overlap, qualified_hit,
    residue_records, standardized_difference,
)


def hit(**kw):
    return {"query_residues": "A1,A2,A3,A4,A5,A6",
            "matching_residues": "B1,B2,B3,B4,B5,_", "rmsd": "1", **kw}


def test_quality_filter_requires_injective_and_alignable_mapping():
    assert qualified_hit(hit())[1] == "qualified"
    assert qualified_hit(hit(matching_residues="B1,B2,B3,B4,B5"))[1] == "unalignable"
    assert qualified_hit(hit(matching_residues="B1,B1,B3,B4,B5,_"))[1] == "noninjective"
    assert qualified_hit(hit(rmsd="nan"))[1] == "rmsd"
    assert qualified_hit(hit(rmsd="-1"))[1] == "rmsd"
    assert qualified_hit(hit(matching_residues="B1,B2,B3,B4,_,_"))[1] == "coverage"


def test_unknown_is_not_zero_and_control_preserves_chain_and_residue_type():
    ann = {"residues": {
        "A1": {"name": "HIS", "zinc": 1},
        "A2": {"name": "HIS", "zinc": 0},
        "A3": {"name": "ALA", "zinc": 0},
        "B1": {"name": "HIS", "zinc": 1},
        "A4": {"name": "HIS", "zinc": None},
    }}
    assert overlap(["A1"], ann, "zinc") == (1, 0.5)
    assert overlap(["A1", "A3"], ann, "zinc") == (0.5, 0.25)
    assert overlap(["A4"], ann, "zinc") is None
    assert overlap(["A99"], ann, "zinc") is None


def atoms():
    a = struc.AtomArray(3)
    a.chain_id = np.array(["A", "A", "A"])
    a.res_id = np.array([1, 1, 2]); a.res_name = np.array(["ALA"] * 3)
    a.atom_name = np.array(["CA"] * 3); a.element = np.array(["C"] * 3)
    a.coord = np.array([[0, 0, 0], [2, 0, 0], [20, 0, 0]], dtype=float)
    a.set_annotation("sym_id", np.array([0, 1, 0]))
    return a


def test_homomeric_symmetry_copies_count_as_different_protein_chains():
    labels, present = assembly_interfaces(atoms())
    assert labels == {"A1"}
    assert present == {"A1", "A2"}


def test_insertion_code_collision_cannot_map_to_folddisco_token():
    a = atoms(); a.ins_code = np.array(["", "A", ""])
    _, valid, ambiguous = residue_records(a)
    assert "A1" not in valid
    assert ambiguous == ["A1"]
    assert "A2" in valid


def test_mcsa_uses_author_numbers_and_only_reference_residues(tmp_path):
    p = tmp_path / "mcsa.json"
    p.write_text(json.dumps([{"mcsa_id": 1, "residue_chains": [
        {"pdb_id": "1ABC", "chain_name": "A", "resid": 5,
         "auth_resid": 105, "code": "His", "is_reference": True},
        {"pdb_id": "2ABC", "chain_name": "A", "auth_resid": 5,
         "code": "His", "is_reference": False},
    ]}]))
    result = load_mcsa(p)
    assert list(result) == ["1abc"]
    assert result["1abc"]["A105"]["name"] == "HIS"


def test_cluster_bootstrap_and_size_standardization_are_reproducible():
    a = [(f"a{i//10}", 6, 1.0) for i in range(30)]
    b = [(f"b{i//10}", 6, 0.0) for i in range(30)]
    result = standardized_difference(a, b, bootstrap=30)
    assert result["difference"] == 1
    assert result["difference_ci95"] == [1, 1]
    assert result["query_structures"] == [3, 3]
    assert result == standardized_difference(a, b, bootstrap=30)
    assert standardized_difference(a[:1], b)["status"] == "insufficient_common_size_strata"