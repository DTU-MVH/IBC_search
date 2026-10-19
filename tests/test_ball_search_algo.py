"""``ball_search_algo`` behaviour: centres, neighbours, discontinuity,
connectivity, deduplication, invariance and determinism."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from ibc_search.ball_search_algo import (
    CENTER_KIND_CA_FALLBACK,
    CENTER_KIND_CB,
    CENTER_KIND_SIDECHAIN_CENTROID,
    CENTER_KIND_VIRTUAL_CB_GLY,
    CENTER_KIND_VIRTUAL_CB_MISSING,
    compute_ball_center,
    run_ball_search,
    virtual_cb,
)
from ibc_search.exposure import compute_exposure
from ibc_search.structure_io import parse_structure, resolve_input

from .conftest import write_rotated_pdb


def pipeline(path: Path, config: dict):
    parsed = parse_structure(resolve_input(str(path), config), config)
    exposure = compute_exposure(parsed, config)
    return parsed, exposure, run_ball_search(parsed, exposure, config)


def labels_of(parsed, patch) -> set[str]:
    return {parsed.residues[m.residue_index].rid.label() for m in patch.members}


# --- C. ball centres --------------------------------------------------------

def test_virtual_cb_formula_reproduces_real_cb():
    """Validate the ideal-geometry CB construction against observed CB atoms.

    Uses a real structure, because the only honest justification for the
    constants is that they reproduce measured C-beta positions.
    """
    import biotite.structure.io.pdb as pdb_io

    real = Path(
        "/Users/mvh_mac/Documents/master_thesis_2026/folddisco_test/pdb_files/4CHA.pdb"
    )
    if not real.is_file():
        pytest.skip("reference structure 4CHA.pdb not available")

    array = pdb_io.PDBFile.read(str(real)).get_structure(model=1, altloc="first")
    deviations = []
    n = len(array)
    i = 0
    while i < n:
        chain, res_id = array.chain_id[i], array.res_id[i]
        atoms: dict[str, np.ndarray] = {}
        j = i
        while j < n and array.chain_id[j] == chain and array.res_id[j] == res_id:
            atoms[str(array.atom_name[j])] = np.asarray(array.coord[j], dtype=float)
            j += 1
        if {"N", "CA", "C", "CB"} <= atoms.keys():
            built = virtual_cb(atoms["N"], atoms["CA"], atoms["C"])
            deviations.append(float(np.linalg.norm(built - atoms["CB"])))
        i = j

    assert len(deviations) > 100
    deviations = np.array(deviations)
    # Ideal geometry cannot match every real side chain exactly, but it must be
    # well inside a bond length, otherwise the constants are wrong.
    assert deviations.mean() < 0.25, f"mean deviation {deviations.mean():.3f} A"
    assert np.percentile(deviations, 95) < 0.5


def test_glycine_uses_flagged_virtual_cb(glycine_pdb, base_config):
    parsed, _, _ = pipeline(glycine_pdb, base_config)
    by_label = {r.rid.label(): i for i, r in enumerate(parsed.residues)}

    gly_center = compute_ball_center(parsed, by_label["A3"], "cb")
    assert gly_center.kind == CENTER_KIND_VIRTUAL_CB_GLY
    assert gly_center.is_virtual is True
    assert "VIRTUAL" in gly_center.note

    ala_center = compute_ball_center(parsed, by_label["A1"], "cb")
    assert ala_center.kind == CENTER_KIND_CB
    assert ala_center.is_virtual is False


def test_missing_cb_in_non_glycine_is_distinguished_from_glycine(
    glycine_pdb, base_config
):
    parsed, _, _ = pipeline(glycine_pdb, base_config)
    by_label = {r.rid.label(): i for i, r in enumerate(parsed.residues)}

    leu = compute_ball_center(parsed, by_label["A5"], "cb")
    # A LEU without a CB is a defect; a GLY without one is not. The kinds must
    # differ so the two are never conflated downstream.
    assert leu.kind == CENTER_KIND_VIRTUAL_CB_MISSING
    assert leu.is_virtual is True
    assert leu.kind != CENTER_KIND_VIRTUAL_CB_GLY


def test_ca_only_residue_falls_back_to_ca(glycine_pdb, base_config):
    parsed, _, _ = pipeline(glycine_pdb, base_config)
    by_label = {r.rid.label(): i for i, r in enumerate(parsed.residues)}
    center = compute_ball_center(parsed, by_label["A7"], "cb")
    assert center.kind == CENTER_KIND_CA_FALLBACK
    assert center.is_virtual is False


def test_sidechain_centroid_mode_uses_sidechain_atoms(altloc_pdb, base_config):
    parsed, _, _ = pipeline(altloc_pdb, base_config)
    center = compute_ball_center(parsed, 0, "sidechain_centroid")
    assert center.kind == CENTER_KIND_SIDECHAIN_CENTROID
    assert center.n_atoms_used == 2  # CB and OG of the chosen altloc

    # The centroid must sit between CB and OG, not at the CA.
    atoms = {
        str(parsed.atom_names[i]): parsed.coords[i]
        for i in range(parsed.residues[0].atom_start, parsed.residues[0].atom_stop)
    }
    expected = (atoms["CB"] + atoms["OG"]) / 2
    assert center.position == pytest.approx(expected, abs=1e-6)


def test_sidechain_centroid_falls_back_for_glycine(glycine_pdb, base_config):
    parsed, _, _ = pipeline(glycine_pdb, base_config)
    by_label = {r.rid.label(): i for i, r in enumerate(parsed.residues)}
    center = compute_ball_center(parsed, by_label["A3"], "sidechain_centroid")
    # Glycine has no side chain, so the documented fallback must engage.
    assert center.kind == CENTER_KIND_VIRTUAL_CB_GLY
    assert center.is_virtual is True


# --- D. neighbour extraction ------------------------------------------------

def test_selecting_atom_and_distance_recorded_within_radius(
    linear_chain_pdb, base_config
):
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["ball_search"]["min_patch_size"] = 2
    parsed, _, result = pipeline(linear_chain_pdb, base_config)
    assert result.patches

    for patch in result.patches:
        for member in patch.members:
            assert member.selecting_atom_name
            # Membership is on minimum heavy-atom distance to the centre, so
            # every recorded distance must be inside the ball.
            assert member.selecting_atom_distance <= patch.radius + 1e-9
            # And the recorded atom must really belong to the member residue.
            rec = parsed.residues[member.residue_index]
            names = {
                str(parsed.atom_names[i])
                for i in range(rec.atom_start, rec.atom_stop)
            }
            assert member.selecting_atom_name in names


def test_larger_radius_is_a_superset(linear_chain_pdb, base_config):
    base_config["ball_search"]["radii_angstrom"] = [5.0, 10.0]
    base_config["ball_search"]["min_patch_size"] = 2
    parsed, _, result = pipeline(linear_chain_pdb, base_config)

    def by_seed(radius):
        out = {}
        for patch in result.patches:
            if patch.radius != radius:
                continue
            for seed in patch.seed_residue_indices:
                out[seed] = set(patch.residue_indices())
        return out

    small, large = by_seed(5.0), by_seed(10.0)
    common = set(small) & set(large)
    assert common
    for seed in common:
        # A ball is monotone in its radius; connectivity pruning cannot remove
        # a member that a smaller ball kept, because the smaller patch is a
        # subset of the larger one and the seed component only grows.
        assert small[seed] <= large[seed], f"seed {seed} not monotone"


def test_all_atoms_of_selected_residues_remain_available(
    linear_chain_pdb, base_config
):
    """Extracting a motif must not discard atoms of its member residues."""
    from ibc_search.ball_search_algo import _patch_atom_indices

    base_config["ball_search"]["radii_angstrom"] = [8.0]
    parsed, _, result = pipeline(linear_chain_pdb, base_config)
    patch = result.patches[0]
    atom_idx = _patch_atom_indices(parsed, list(patch.residue_indices()))

    expected = sum(
        parsed.residues[i].atom_stop - parsed.residues[i].atom_start
        for i in patch.residue_indices()
    )
    assert len(atom_idx) == expected
    # Backbone atoms are still there, not just CB.
    names = {str(parsed.atom_names[i]) for i in atom_idx}
    assert {"N", "CA", "C", "O"} <= names


# --- E. discontinuity -------------------------------------------------------

def test_hairpin_yields_sequence_discontinuous_patches(hairpin_pdb, base_config):
    """The extractor must be a 3D search, not a sequence window."""
    base_config["ball_search"]["radii_angstrom"] = [8.0]
    parsed, _, result = pipeline(hairpin_pdb, base_config)

    discontinuous = [p for p in result.patches if p.is_discontinuous]
    assert discontinuous, "a hairpin must produce discontinuous patches"

    patch = discontinuous[0]
    assert patch.n_segments >= 2
    assert len(patch.segment_lengths) == patch.n_segments
    assert sum(patch.segment_lengths) == patch.size
    assert len(patch.author_numbering_gaps) >= 1
    assert all(g > 0 for g in patch.author_numbering_gaps)

    # At least one patch must pair residues from both arms of the hairpin.
    assert any(
        any(int(lbl[1:]) <= 6 for lbl in labels_of(parsed, p))
        and any(int(lbl[1:]) >= 7 for lbl in labels_of(parsed, p))
        for p in result.patches
    )


def test_contiguous_run_reports_one_segment(linear_chain_pdb, base_config):
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["ball_search"]["min_patch_size"] = 3
    parsed, _, result = pipeline(linear_chain_pdb, base_config)
    middles = [
        p for p in result.patches
        if p.n_segments == 1 and p.size == 3
    ]
    assert middles
    for patch in middles:
        assert patch.segment_lengths == (3,)
        assert patch.author_numbering_gaps == ()
        assert patch.is_discontinuous is False


# --- F. connectivity --------------------------------------------------------

def test_disconnected_cluster_is_dropped_by_default(bridged_pdb, base_config):
    """A ball spanning two non-contacting clusters must not merge them."""
    base_config["ball_search"]["radii_angstrom"] = [10.0]
    base_config["ball_search"]["min_patch_size"] = 2
    base_config["connectivity"]["disconnected_policy"] = "drop"
    base_config["connectivity"]["contact_cutoff_angstrom"] = 5.0
    parsed, _, result = pipeline(bridged_pdb, base_config)

    near = {"A1", "A2"}
    far = {"A10", "A11"}
    assert result.patches
    for patch in result.patches:
        members = labels_of(parsed, patch)
        # Every patch must lie wholly within one cluster.
        assert members <= near or members <= far, members
    assert any(p.n_dropped_disconnected > 0 for p in result.patches)
    assert any(p.n_components_before_connectivity > 1 for p in result.patches)


def test_disconnected_cluster_can_be_flagged_instead(bridged_pdb, base_config):
    base_config["ball_search"]["radii_angstrom"] = [10.0]
    base_config["ball_search"]["min_patch_size"] = 2
    base_config["connectivity"]["disconnected_policy"] = "flag"
    parsed, _, result = pipeline(bridged_pdb, base_config)

    merged = [
        p for p in result.patches
        if any(m.disconnected_flagged for m in p.members)
    ]
    assert merged, "flag policy must retain the far cluster, marked"
    patch = merged[0]
    assert len(labels_of(parsed, patch)) == 4
    assert patch.n_dropped_disconnected == 0


def test_keep_policy_ignores_connectivity(bridged_pdb, base_config):
    base_config["ball_search"]["radii_angstrom"] = [10.0]
    base_config["ball_search"]["min_patch_size"] = 2
    base_config["connectivity"]["disconnected_policy"] = "keep"
    parsed, _, result = pipeline(bridged_pdb, base_config)
    assert any(len(labels_of(parsed, p)) == 4 for p in result.patches)


def test_buried_members_are_excluded_or_flagged(base_config, tmp_path):
    """A residue below the member RSA threshold must be dropped or flagged."""
    from .conftest import Model, Residue, backbone_atoms, write_pdb

    # A compact 3x3x3 lattice: the centre residue is enclosed and so buried,
    # while the shell residues are exposed.
    residues = []
    n = 1
    for x in range(3):
        for y in range(3):
            for z in range(3):
                residues.append(
                    Residue(
                        "A", n, "ALA",
                        backbone_atoms(np.array([4.5 * x, 4.5 * y, 4.5 * z])),
                    )
                )
                n += 1
    path = write_pdb(tmp_path / "lattice.pdb", [Model(residues)])

    base_config["ball_search"]["radii_angstrom"] = [8.0]
    base_config["ball_search"]["min_patch_size"] = 2
    base_config["ball_search"]["buried_member_policy"] = "drop"
    parsed, exposure, result = pipeline(path, base_config)

    threshold = base_config["ball_search"]["member_rsa_threshold"]
    buried = {
        i for i in range(len(parsed.residues))
        if (exposure.rsa_value_or_none(i) or 0.0) < threshold
    }
    assert buried, "the lattice must contain at least one buried residue"
    for patch in result.patches:
        for member in patch.members:
            if member.residue_index in buried and not member.is_seed:
                raise AssertionError("buried member survived the drop policy")
    assert sum(p.n_dropped_buried for p in result.patches) > 0

    base_config["ball_search"]["buried_member_policy"] = "flag"
    _, _, flagged_result = pipeline(path, base_config)
    assert any(
        m.buried_flagged for p in flagged_result.patches for m in p.members
    )


# --- G. deduplication -------------------------------------------------------

def test_identical_member_sets_are_deduplicated_retaining_seeds(
    linear_chain_pdb, base_config
):
    # A radius larger than the whole fixture guarantees that every seed
    # selects the identical member set, which is precisely the case
    # deduplication must collapse without losing any seed.
    base_config["ball_search"]["radii_angstrom"] = [50.0]
    base_config["ball_search"]["min_patch_size"] = 2
    parsed, _, result = pipeline(linear_chain_pdb, base_config)

    # Member sets must be unique within a radius.
    sets = [frozenset(p.residue_indices()) for p in result.patches]
    assert len(sets) == len(set(sets)), "duplicate member sets were emitted"

    # All ten residues, reachable from every seed, must collapse to one patch
    # that retains every contributing seed.
    assert len(result.patches) == 1
    patch = result.patches[0]
    assert patch.size == 10
    assert len(patch.seed_residue_indices) == len(result.seed_residue_indices) > 1
    for patch in result.patches:
        assert len(set(patch.seed_residue_indices)) == len(patch.seed_residue_indices)
        # Every retained seed must itself be a member, marked as a seed.
        seed_members = {m.residue_index for m in patch.members if m.is_seed}
        assert set(patch.seed_residue_indices) <= seed_members


def test_jaccard_clustering_groups_overlapping_patches(linear_chain_pdb, base_config):
    base_config["ball_search"]["radii_angstrom"] = [10.0]
    base_config["dedup"]["cluster_overlapping"] = True
    base_config["dedup"]["jaccard_threshold"] = 0.6
    _, _, result = pipeline(linear_chain_pdb, base_config)

    assert all(p.cluster_id for p in result.patches)
    clusters = {p.cluster_id for p in result.patches}
    # Overlapping balls on a short chain must collapse to fewer clusters than
    # patches, but not to a single cluster at this threshold.
    assert 1 <= len(clusters) < len(result.patches)

    base_config["dedup"]["jaccard_threshold"] = 1.0
    _, _, strict = pipeline(linear_chain_pdb, base_config)
    # At threshold 1.0 only identical sets cluster, and those are already
    # deduplicated, so every patch stands alone.
    assert len({p.cluster_id for p in strict.patches}) == len(strict.patches)


# --- H. scales --------------------------------------------------------------

def test_multiple_scales_are_recorded_and_linked(linear_chain_pdb, base_config):
    # scale_labels is set explicitly here rather than relying on whatever the
    # shipped defaults happen to map a given radius to -- the production
    # default set has changed twice already (A-3.7/A-3.7b) as radii were
    # found to saturate on minibinder-sized inputs, and this test is about
    # the *lookup mechanism*, not about any particular production radius.
    base_config["ball_search"]["radii_angstrom"] = [5.0, 12.0]
    base_config["ball_search"]["scale_labels"] = {
        "5.0": "submotif", "12.0": "epitope_scale",
    }
    base_config["ball_search"]["min_patch_size"] = 2
    _, _, result = pipeline(linear_chain_pdb, base_config)

    radii = {p.radius for p in result.patches}
    assert radii == {5.0, 12.0}
    scales = {p.radius: p.scale for p in result.patches}
    assert scales[5.0] == "submotif"
    assert scales[12.0] == "epitope_scale"

    # The member_set_hash is keyed on author identity, so the same residue set
    # found at two radii is detectably the same set.
    small = {p.member_set_hash for p in result.patches if p.radius == 5.0}
    assert all(isinstance(h, str) and len(h) == 16 for h in small)


def test_patch_size_bounds_are_enforced(linear_chain_pdb, base_config):
    base_config["ball_search"]["radii_angstrom"] = [15.0]
    base_config["ball_search"]["min_patch_size"] = 2
    base_config["ball_search"]["max_patch_size"] = 4
    _, _, result = pipeline(linear_chain_pdb, base_config)

    oversized = [p for p in result.patches if p.size > 4]
    assert oversized, "the fixture should produce patches above the cap"
    for patch in oversized:
        # Oversized patches stay in the table but emit no query, and the
        # reason must point at the configuration, not at FoldDisco.
        assert patch.emit_query is False
        assert "max_patch_size" in patch.emit_block_reason
        assert patch.folddisco_expressible is True
    assert all(p.size >= 2 for p in result.patches)


# --- invariance and determinism --------------------------------------------

def test_patches_are_invariant_under_rigid_motion(hairpin_pdb, base_config, tmp_path):
    """A rigid rotation + translation must not change what is extracted.

    Tolerance note: the rotated fixture is written back out as a PDB file,
    whose coordinates are quantised to 3 decimal places. That round-trip alone
    perturbs a distance by up to a few times 1e-3 A, so distances are compared
    at 1e-2 A. The substantive invariance claims -- identical member sets,
    identical selecting atoms, identical segment structure -- are compared
    exactly.
    """
    coordinate_tolerance = 1e-2
    base_config["ball_search"]["radii_angstrom"] = [8.0]
    rotated = write_rotated_pdb(hairpin_pdb, tmp_path / "hairpin_rot.pdb")

    parsed_a, _, result_a = pipeline(hairpin_pdb, base_config)
    parsed_b, _, result_b = pipeline(rotated, base_config)

    assert len(result_a.patches) == len(result_b.patches)
    for patch_a, patch_b in zip(result_a.patches, result_b.patches):
        assert labels_of(parsed_a, patch_a) == labels_of(parsed_b, patch_b)
        assert patch_a.n_segments == patch_b.n_segments
        assert patch_a.segment_lengths == patch_b.segment_lengths
        # Geometry metrics are rotation/translation invariant quantities.
        assert patch_a.spatial_diameter == pytest.approx(
            patch_b.spatial_diameter, abs=coordinate_tolerance
        )
        assert patch_a.radius_of_gyration == pytest.approx(
            patch_b.radius_of_gyration, abs=coordinate_tolerance
        )
        for member_a, member_b in zip(patch_a.members, patch_b.members):
            assert member_a.selecting_atom_name == member_b.selecting_atom_name
            assert member_a.selecting_atom_distance == pytest.approx(
                member_b.selecting_atom_distance, abs=coordinate_tolerance
            )


def test_output_is_deterministic_across_repeated_runs(hairpin_pdb, base_config):
    base_config["ball_search"]["radii_angstrom"] = [8.0]
    _, _, first = pipeline(hairpin_pdb, base_config)
    _, _, second = pipeline(hairpin_pdb, base_config)

    assert [p.patch_id for p in first.patches] == [p.patch_id for p in second.patches]
    assert [p.query_string for p in first.patches] == [
        p.query_string for p in second.patches
    ]
    assert [p.cluster_id for p in first.patches] == [
        p.cluster_id for p in second.patches
    ]


@pytest.mark.parametrize("workers", [1, 2, 4])
def test_output_is_identical_across_worker_counts(hairpin_pdb, base_config, workers):
    base_config["ball_search"]["radii_angstrom"] = [8.0, 12.0]
    base_config["performance"]["workers"] = 1
    _, _, reference = pipeline(hairpin_pdb, base_config)

    base_config["performance"]["workers"] = workers
    parsed, _, actual = pipeline(hairpin_pdb, base_config)

    assert [p.patch_id for p in reference.patches] == [
        p.patch_id for p in actual.patches
    ]
    for patch_a, patch_b in zip(reference.patches, actual.patches):
        assert patch_a.residue_indices() == patch_b.residue_indices()
        assert patch_a.seed_residue_indices == patch_b.seed_residue_indices
        assert patch_a.query_string == patch_b.query_string
        assert patch_a.cluster_id == patch_b.cluster_id


# --- seeds ------------------------------------------------------------------

def test_no_seeds_is_handled_gracefully(linear_chain_pdb, base_config):
    base_config["ball_search"]["seed_rsa_threshold"] = 1.5  # unreachable
    _, _, result = pipeline(linear_chain_pdb, base_config)
    assert result.patches == []
    assert result.stats["n_seeds"] == 0
    assert any("no residue reached the seed RSA threshold" in w for w in result.warnings)


def test_seed_is_always_a_member_of_its_own_patch(hairpin_pdb, base_config):
    base_config["ball_search"]["radii_angstrom"] = [5.0, 8.0]
    _, _, result = pipeline(hairpin_pdb, base_config)
    assert result.patches
    for patch in result.patches:
        members = set(patch.residue_indices())
        for seed in patch.seed_residue_indices:
            assert seed in members


def test_seed_normal_filter_runs_and_can_only_remove_members(
    hairpin_pdb, base_config
):
    base_config["ball_search"]["radii_angstrom"] = [10.0]
    base_config["connectivity"]["use_seed_normal_filter"] = False
    _, _, without = pipeline(hairpin_pdb, base_config)

    base_config["connectivity"]["use_seed_normal_filter"] = True
    _, _, with_filter = pipeline(hairpin_pdb, base_config)

    def sets_by_seed(result):
        out: dict[int, set[int]] = {}
        for patch in result.patches:
            for seed in patch.seed_residue_indices:
                out[seed] = set(patch.residue_indices())
        return out

    a, b = sets_by_seed(without), sets_by_seed(with_filter)
    for seed in set(a) & set(b):
        assert b[seed] <= a[seed], "the normal filter must only remove members"
