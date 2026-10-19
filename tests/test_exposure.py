"""SASA / RSA computation and the "missing is not zero" rule."""

from __future__ import annotations

import pytest

from ibc_search.errors import Missing
from ibc_search.exposure import (
    TIEN2013_THEORETICAL_MAX_ACC,
    compute_exposure,
)
from ibc_search.structure_io import parse_structure, resolve_input

from .conftest import Atom, Model, Residue, backbone_atoms, write_pdb


def _exposure(path, config):
    parsed = parse_structure(resolve_input(str(path), config), config)
    return parsed, compute_exposure(parsed, config)


def test_rsa_reference_table_is_complete_and_named():
    # All 20 canonical residue types, from the named reference.
    assert len(TIEN2013_THEORETICAL_MAX_ACC) == 20
    assert TIEN2013_THEORETICAL_MAX_ACC["GLY"] == 104.0
    assert TIEN2013_THEORETICAL_MAX_ACC["TRP"] == 285.0


def test_exposed_residues_have_positive_sasa_and_rsa(linear_chain_pdb, base_config):
    parsed, exposure = _exposure(linear_chain_pdb, base_config)
    assert exposure.probe_radius == 1.4
    assert exposure.rsa_reference == "tien2013_theoretical"
    assert exposure.n_context_atoms == len(parsed.coords)

    for i in range(len(parsed.residues)):
        sasa = exposure.sasa[i]
        rsa = exposure.rsa[i]
        assert not isinstance(sasa, Missing)
        assert not isinstance(rsa, Missing)
        # An isolated extended chain is almost entirely solvent-exposed.
        assert sasa > 0.0
        assert rsa == pytest.approx(sasa / TIEN2013_THEORETICAL_MAX_ACC["ALA"])


def test_probe_radius_is_configurable_and_changes_sasa(linear_chain_pdb, base_config):
    _, small = _exposure(linear_chain_pdb, base_config)
    base_config["exposure"]["probe_radius_angstrom"] = 2.5
    _, large = _exposure(linear_chain_pdb, base_config)

    total_small = sum(v for v in small.sasa if not isinstance(v, Missing))
    total_large = sum(v for v in large.sasa if not isinstance(v, Missing))
    assert large.probe_radius == 2.5
    # A bigger probe rolls over crevices, changing the measured area.
    assert total_small != pytest.approx(total_large)


def test_nonstandard_residue_rsa_is_missing_with_a_reason_not_zero(
    base_config, tmp_path
):
    """A residue type with no reference maximum must yield MISSING, not 0."""
    residues = [
        Residue("A", 1, "ALA", backbone_atoms([0.0, 0.0, 0.0])),
        # MSE (selenomethionine) is a genuine polymer residue with no Tien
        # et al. maximum.
        Residue(
            "A", 2, "MSE",
            backbone_atoms([4.0, 0.0, 0.0]) + [Atom("CG", (4.0, -2.9, 0.0))],
            hetero=True,
        ),
        Residue("A", 3, "ALA", backbone_atoms([8.0, 0.0, 0.0])),
    ]
    path = write_pdb(tmp_path / "mse.pdb", [Model(residues)])
    parsed, exposure = _exposure(path, base_config)

    mse_index = next(
        i for i, r in enumerate(parsed.residues) if r.rid.res_name == "MSE"
    )
    # MSE is classified as a polymer residue despite being a HETATM record.
    assert parsed.residues[mse_index].residue_class == "polymer_nonstandard_aa"

    # Its SASA is a real number, but its RSA is unavailable with a reason.
    assert not isinstance(exposure.sasa[mse_index], Missing)
    rsa = exposure.rsa[mse_index]
    assert isinstance(rsa, Missing)
    assert "Tien" in rsa.reason
    assert str(rsa) == ""  # renders as an empty cell, not "0"
    assert exposure.rsa_value_or_none(mse_index) is None
    # A MISSING RSA must never read as exposed at any threshold.
    assert exposure.is_exposed(mse_index, 0.0) is False


def test_residue_outside_the_sasa_context_is_missing_with_a_reason(
    two_chain_pdb, base_config
):
    base_config["input"]["chains"] = ["A"]
    base_config["exposure"]["context"] = "scanned_chains"
    parsed, exposure = _exposure(two_chain_pdb, base_config)

    chain_b = [
        i for i, r in enumerate(parsed.residues) if r.rid.auth_chain_id == "B"
    ]
    assert chain_b
    for i in chain_b:
        sasa = exposure.sasa[i]
        assert isinstance(sasa, Missing)
        assert "not in the SASA context" in sasa.reason


def test_waters_are_never_sasa_occluders(waters_and_ligand_pdb, base_config):
    parsed, exposure = _exposure(waters_and_ligand_pdb, base_config)
    water_indices = [
        i for i, r in enumerate(parsed.residues) if r.rid.res_name == "HOH"
    ]
    assert water_indices
    for i in water_indices:
        assert isinstance(exposure.sasa[i], Missing)


def test_hetero_occluder_setting_changes_exposure(waters_and_ligand_pdb, base_config):
    base_config["input"]["chains"] = ["A"]
    base_config["exposure"]["include_hetero_as_occluder"] = False
    _, without = _exposure(waters_and_ligand_pdb, base_config)
    assert any("excluded from the SASA context" in w for w in without.warnings)

    base_config["exposure"]["include_hetero_as_occluder"] = True
    _, with_hetero = _exposure(waters_and_ligand_pdb, base_config)
    # The ion and sugar now occlude, so the context grows.
    assert with_hetero.n_context_atoms > without.n_context_atoms


def test_exposure_context_choice_is_recorded(two_chain_pdb, base_config):
    base_config["input"]["chains"] = ["A"]
    _, scanned = _exposure(two_chain_pdb, base_config)
    assert "context=scanned_chains" in scanned.context_description

    base_config["exposure"]["context"] = "all_polymer_chains"
    _, everything = _exposure(two_chain_pdb, base_config)
    assert "context=all_polymer_chains" in everything.context_description
    assert everything.n_context_atoms > scanned.n_context_atoms


def test_buried_residue_has_low_but_real_sasa(base_config, tmp_path):
    """A buried residue must report a small number, distinct from MISSING."""
    residues = []
    n = 1
    for x in range(3):
        for y in range(3):
            for z in range(3):
                residues.append(
                    Residue("A", n, "ALA", backbone_atoms([4.5 * x, 4.5 * y, 4.5 * z]))
                )
                n += 1
    path = write_pdb(tmp_path / "lattice.pdb", [Model(residues)])
    parsed, exposure = _exposure(path, base_config)

    values = [exposure.rsa_value_or_none(i) for i in range(len(parsed.residues))]
    assert all(v is not None for v in values)
    assert min(values) < 0.05  # the enclosed centre
    assert max(values) > 0.20  # the shell
