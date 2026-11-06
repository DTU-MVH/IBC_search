"""Where the background balls are allowed to sit.

The `surface` mode exists to remove a confound that was measured, not
imagined: with the default `protein` mode the two classes landed at opposite
ends of the exposure scale wherever the annotation was extreme. So these
tests pin the properties that make `surface` a control rather than a label --
every chosen ball is genuinely exposed, an uncomputable RSA is never assumed
exposed, and a protein that cannot supply enough surface says so instead of
quietly reaching into the interior.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from ibc_search.background_selection import (
    MODE_PROTEIN,
    MODE_SURFACE,
    STRATUM_BURIED,
    STRATUM_EXPOSED,
)
from ibc_search.exposure import compute_exposure
from ibc_search.few_balls_search import (
    SEED_CLASS_ANNOTATED,
    SEED_CLASS_BACKGROUND,
    run_few_balls_search,
)
from ibc_search.structure_io import parse_structure, resolve_input


def seeds_file(tmp_path: Path, labels, name="seeds.tsv") -> Path:
    path = tmp_path / name
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["chain", "res_id", "is_annotated"])
        for chain, res in labels:
            w.writerow([chain, res, 1])
    return path


def run(path: Path, config: dict, seeds: Path):
    config["ball_search"]["seed_residues_path"] = str(seeds)
    parsed = parse_structure(resolve_input(str(path), config), config)
    exposure = compute_exposure(parsed, config)
    return parsed, exposure, run_few_balls_search(parsed, exposure, config)


def background_indices(result) -> list[int]:
    return [p.seed_residue_indices[0] for p in result.patches
            if p.seed_class == SEED_CLASS_BACKGROUND]


def test_surface_mode_picks_only_exposed_residues(buried_core_pdb, base_config):
    """Every background ball must sit above the exposure threshold."""
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["few_balls_search"]["background_mode"] = MODE_SURFACE
    base_config["few_balls_search"]["min_background_balls"] = 6
    base_config["few_balls_search"]["background_rsa_threshold"] = 0.20
    seeds = seeds_file(buried_core_pdb.parent, [("A", 1)])
    _, exposure, result = run(buried_core_pdb, base_config, seeds)

    chosen = background_indices(result)
    assert chosen
    for i in chosen:
        rsa = exposure.rsa_value_or_none(i)
        assert rsa is not None and rsa >= 0.20


def test_protein_mode_still_reaches_the_buried_core(buried_core_pdb,
                                                    base_config):
    """The default is unchanged: it samples the whole protein."""
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["few_balls_search"]["background_mode"] = MODE_PROTEIN
    base_config["few_balls_search"]["min_background_balls"] = 6
    seeds = seeds_file(buried_core_pdb.parent, [("A", 1)])
    _, _, result = run(buried_core_pdb, base_config, seeds)

    strata = result.stats["per_stratum"]
    assert strata[STRATUM_BURIED]["n_selected"] >= 1
    assert strata[STRATUM_EXPOSED]["n_selected"] >= 1


def test_surface_mode_does_not_stratify_and_says_so(buried_core_pdb,
                                                    base_config):
    """A pool that is exposed by construction has nothing to stratify."""
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["few_balls_search"]["background_mode"] = MODE_SURFACE
    base_config["few_balls_search"]["min_background_balls"] = 4
    seeds = seeds_file(buried_core_pdb.parent, [("A", 1)])
    _, _, result = run(buried_core_pdb, base_config, seeds)

    assert result.stats["background_mode"] == MODE_SURFACE
    # Not left as a misleading `true` in the manifest.
    assert result.stats["background_stratified_by_burial"] is False


def test_surface_mode_shrinks_the_candidate_pool(buried_core_pdb, base_config):
    """The filter must be visible in the audit, not just in its effect."""
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["few_balls_search"]["background_mode"] = MODE_SURFACE
    base_config["few_balls_search"]["min_background_balls"] = 4
    seeds = seeds_file(buried_core_pdb.parent, [("A", 1)])
    _, _, result = run(buried_core_pdb, base_config, seeds)

    before = result.stats["n_candidates_before_surface_filter"]
    after = result.stats["n_background_candidates"]
    assert after <= before


def test_surface_mode_under_delivers_rather_than_going_inside(base_config,
                                                              buried_core_pdb):
    """Asking for more surface than exists is reported, never back-filled."""
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["few_balls_search"]["background_mode"] = MODE_SURFACE
    base_config["few_balls_search"]["min_background_balls"] = 10_000
    seeds = seeds_file(buried_core_pdb.parent, [("A", 1)])
    _, exposure, result = run(buried_core_pdb, base_config, seeds)

    assert result.stats["n_background_shortfall"] > 0
    assert any("background balls and placed" in w for w in result.warnings)
    # Whatever WAS placed is still all surface -- the shortfall did not
    # license reaching into the interior.
    for i in background_indices(result):
        assert (exposure.rsa_value_or_none(i) or -1.0) >= 0.20


def test_both_modes_record_the_realised_exposure_of_both_classes(
    buried_core_pdb, base_config
):
    """`surface` equalises the pool, not the distribution: show the residual."""
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["few_balls_search"]["min_background_balls"] = 4
    seeds = seeds_file(buried_core_pdb.parent, [("A", 1)])
    for mode in (MODE_PROTEIN, MODE_SURFACE):
        base_config["few_balls_search"]["background_mode"] = mode
        _, _, result = run(buried_core_pdb, base_config, seeds)
        assert result.stats["mean_rsa_annotated"] is not None
        assert result.stats["mean_rsa_background"] is not None


def test_surface_background_is_more_exposed_than_protein_background(
    buried_core_pdb, base_config
):
    """The mode must actually change the control, not just label it."""
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["few_balls_search"]["min_background_balls"] = 6
    seeds = seeds_file(buried_core_pdb.parent, [("A", 1)])

    base_config["few_balls_search"]["background_mode"] = MODE_PROTEIN
    _, _, a = run(buried_core_pdb, base_config, seeds)
    base_config["few_balls_search"]["background_mode"] = MODE_SURFACE
    _, _, b = run(buried_core_pdb, base_config, seeds)

    assert b.stats["mean_rsa_background"] >= a.stats["mean_rsa_background"]


def test_an_unknown_mode_is_rejected(buried_core_pdb, base_config):
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["few_balls_search"]["background_mode"] = "whatever"
    seeds = seeds_file(buried_core_pdb.parent, [("A", 1)])
    with pytest.raises(ValueError, match="background_mode"):
        run(buried_core_pdb, base_config, seeds)


def test_classes_are_still_labelled_and_matched_in_surface_mode(
    buried_core_pdb, base_config
):
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["few_balls_search"]["background_mode"] = MODE_SURFACE
    base_config["few_balls_search"]["min_background_balls"] = 2
    seeds = seeds_file(buried_core_pdb.parent, [("A", 2), ("A", 3), ("A", 4)])
    _, _, result = run(buried_core_pdb, base_config, seeds)

    ann = sum(1 for p in result.patches
              if p.seed_class == SEED_CLASS_ANNOTATED)
    bg = sum(1 for p in result.patches
             if p.seed_class == SEED_CLASS_BACKGROUND)
    assert ann == 3
    assert bg == 3
