"""Insertion-aware PPI identities, partner selection and output validation."""
from pathlib import Path
import sys

import numpy as np
import pytest
import biotite.structure as struc

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from build_true_annotated_ppi import (
    ANTIBODY, canonical_sequence, make_tables, remap, select_target,
)


def small_arrays():
    a = struc.AtomArray(3)
    a.chain_id = ["LONG", "LONG", "LONG"]
    a.res_id = [52, 52, 53]
    a.ins_code = ["", "A", ""]
    a.res_name = ["ALA", "CYS", "GLY"]
    a.atom_name = ["CA", "CA", "CA"]
    a.element = ["C", "C", "C"]
    a.coord = np.array([[0, 0, 0], [20, 0, 0], [40, 0, 0]], dtype=float)
    a.set_annotation("label_entity_id", np.array(["1", "1", "1"]))
    b = a[:1].copy(); b.chain_id = np.array(["PARTNER"]); b.res_id[:] = 1; b.coord = np.array([[4, 0, 0]], dtype=float)
    return a, b


def test_insertion_codes_not_merged_and_coordinates_preserved():
    target, partner = small_arrays()
    export, mapping, labels, contacts = make_tables(target, partner, "test")
    assert list(export.chain_id) == ["A"] * 3
    assert list(export.res_id) == [1, 2, 3]
    assert np.array_equal(export.coord, target.coord)
    assert [r[5] for r in labels] == [1, 0, 0]
    assert mapping[0][6:8] == [52, ""]
    assert mapping[1][6:8] == [52, "A"]
    assert len(contacts) == 1
    assert contacts[0][5] == "PARTNER"
    assert canonical_sequence(target) == "ALA-CYS-GLY"


def test_mapping_is_reversible():
    target, _ = small_arrays()
    _, mapping, lookup = remap(target, "test")
    source_keys = [(r[5], r[6], r[7], r[4]) for r in mapping]
    assert len(set(source_keys)) == 3
    assert [lookup[k] for k in source_keys] == [1, 2, 3]


def long_array(counts):
    a = struc.AtomArray(sum(counts))
    chains = []; numbers = []; coords = []
    for i, n in enumerate(counts):
        chains += [chr(65+i)] * n
        numbers += list(range(1, n+1))
        coords += [[float(j)*3.8, float(i)*4, 0.] for j in range(n)]
    a.chain_id = chains; a.res_id = numbers; a.res_name[:] = "ALA"
    a.atom_name[:] = "CA"; a.element[:] = "C"; a.coord = np.array(coords, dtype=float)
    a.set_annotation("label_entity_id", np.array(chains))
    return a


def test_peptide_partner_excluded():
    with pytest.raises(ValueError, match="fewer_than_two"):
        select_target(long_array([100, 20]))
    target, partners, chains = select_target(long_array([100, 90, 20]))
    assert set(target.chain_id) == {"A"}
    assert set(partners.chain_id) == {"B"}
    assert chains == ["B"]


def test_antibody_names_detected():
    assert ANTIBODY.search("Fab heavy chain")
    assert ANTIBODY.search("Nanobody")
    assert not ANTIBODY.search("DNA-binding protein")