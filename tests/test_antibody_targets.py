"""Offline interface geometry and target-only extraction tests."""
import sys
from pathlib import Path

import numpy as np
import pytest
import biotite.structure as struc

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from build_antibody_targets import chain_selection, contact_tables, extract, fingerprint


def metadata():
    return {"agchains": "A", "agtypes": "PROTEIN", "agclusters": "P1", "Hchain": "H", "Lchain": ""}


def atoms():
    arr = struc.AtomArray(6)
    arr.chain_id = ["A", "A", "A", "H", "H", "A"]
    arr.res_id = [1, 2, 3, 10, 11, 4]
    arr.res_name = ["ALA", "CYS", "GLY", "TYR", "SER", "HOH"]
    arr.atom_name = ["CA", "CA", "CA", "CA", "CA", "O"]
    arr.element = ["C", "C", "C", "C", "C", "O"]
    arr.hetero = [False, False, False, False, False, True]
    arr.ins_code = ["", "A", "", "", "", ""]
    arr.coord = np.array([[0,0,0], [9,0,0], [30,0,0], [4,0,0], [13,0,0], [4,1,0]])
    return arr


def test_target_only_no_water_and_all_partner_pairs():
    ag, ab = extract(atoms(), metadata())
    assert set(ag.chain_id) == {"A"}
    assert len(ag) == 3
    rows, pairs = contact_tables(ag, ab, "test", 5.)
    assert [r[5] for r in rows] == [1, 1, 0]
    assert rows[1][3] == "A"  # insertion code retained
    assert len(pairs) == 3  # residue 2 contacts two antibody residues
    assert rows[2][7] == "17.000000"


def test_cutoff_boundary_and_no_self_labels():
    ag, ab = extract(atoms(), metadata())
    rows, pairs = contact_tables(ag, ab, "test", 4.)
    assert len(pairs) == 2
    assert all(r[1] != r[5] for r in pairs)
    with pytest.raises(ValueError):
        contact_tables(ag, ab, "test", 0)


def test_chain_overlap_and_ambiguous_metadata_rejected():
    row = metadata(); row["Hchain"] = "A"
    with pytest.raises(ValueError):
        chain_selection(row)
    row = metadata(); row["agclusters"] = "P1/P2"
    with pytest.raises(ValueError):
        chain_selection(row)


def test_exact_sequence_fingerprint_ignores_numbering():
    ag, _ = extract(atoms(), metadata())
    other = ag.copy(); other.res_id += 20; other.chain_id[:] = "B"
    assert fingerprint(ag) == fingerprint(other)


def test_missing_antibody_rejected():
    row = metadata(); row["Lchain"] = "L"
    with pytest.raises(ValueError, match="missing selected"):
        extract(atoms(), row)


def test_contacts_match_independent_brute_force():
    ag, ab = extract(atoms(), metadata())
    rs, cs = contact_tables(ag, ab, "test", 5.)
    matrix = np.linalg.norm(ag.coord[:, None, :].astype(float) - ab.coord[None, :, :].astype(float), axis=2)
    assert [r[5] for r in rs] == list((matrix <= 5.).any(axis=1).astype(int))
    assert [float(r[7]) for r in rs] == list(matrix.min(axis=1))
    assert len(cs) == int((matrix <= 5.).sum())