"""``few_balls_search``: annotated balls plus a matched, dispersed background.

The method's value is entirely in the properties of the background set, so
that is what these tests pin: it is the same size as the annotated set, it
honours the floor, it never lands on an annotated residue, it spans buried and
exposed positions rather than only the surface, it is spread out, and every
motif says which class it belongs to. A background set that quietly became a
surface sample, or that silently shrank, would still produce a plausible table
-- so each of those is asserted rather than assumed.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

from ibc_search.exposure import compute_exposure
from ibc_search.few_balls_search import (
    SEED_CLASS_ANNOTATED,
    SEED_CLASS_BACKGROUND,
    STRATUM_BURIED,
    STRATUM_EXPOSED,
    _allocate,
    _farthest_point_select,
    run_few_balls_search,
)
from ibc_search.structure_io import parse_structure, resolve_input


def seeds_file(tmp_path: Path, labels: list[tuple[str, int]],
               name: str = "seeds.tsv") -> Path:
    path = tmp_path / name
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["chain", "res_id", "is_annotated"])
        for chain, res in labels:
            w.writerow([chain, res, 1])
    return path


def pipeline(path: Path, config: dict, seeds: Path | None):
    config["ball_search"]["seed_residues_path"] = str(seeds) if seeds else ""
    parsed = parse_structure(resolve_input(str(path), config), config)
    exposure = compute_exposure(parsed, config)
    return parsed, run_few_balls_search(parsed, exposure, config)


def classes(result) -> dict[str, int]:
    out = {SEED_CLASS_ANNOTATED: 0, SEED_CLASS_BACKGROUND: 0}
    for p in result.patches:
        out[p.seed_class] += 1
    return out


def seed_labels(parsed, result, seed_class: str) -> set[str]:
    return {
        parsed.residues[p.seed_residue_indices[0]].rid.label()
        for p in result.patches if p.seed_class == seed_class
    }


# --- the matched count ------------------------------------------------------

def test_background_matches_the_annotated_count(linear_chain_pdb, base_config):
    """Four annotated residues, above the floor: four background balls.

    Four of ten, so six non-annotated residues remain and the match is
    actually available -- annotating more than half the chain would make the
    count unreachable and test the shortfall path instead.
    """
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["few_balls_search"]["min_background_balls"] = 2
    seeds = seeds_file(linear_chain_pdb.parent,
                       [("A", i) for i in range(1, 5)])
    _, result = pipeline(linear_chain_pdb, base_config, seeds)

    assert result.stats["n_annotated_balls"] == 4
    assert result.stats["n_background_balls"] == 4
    assert result.stats["n_background_shortfall"] == 0
    assert result.stats["classes_equal_size"] is True
    assert classes(result) == {SEED_CLASS_ANNOTATED: 4,
                               SEED_CLASS_BACKGROUND: 4}


def test_the_floor_applies_when_there_are_few_annotations(
    linear_chain_pdb, base_config
):
    """Two annotated residues still get the minimum background set."""
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["few_balls_search"]["min_background_balls"] = 5
    seeds = seeds_file(linear_chain_pdb.parent, [("A", 1), ("A", 2)])
    _, result = pipeline(linear_chain_pdb, base_config, seeds)

    assert result.stats["n_annotated_balls"] == 2
    assert result.stats["n_background_balls"] == 5
    assert result.stats["classes_equal_size"] is False
    # The inequality is deliberate, so it must be said out loud rather than
    # left for the reader to notice in the counts.
    assert any("min_background_balls" in w for w in result.warnings)


def test_a_shortfall_is_recorded_not_hidden(linear_chain_pdb, base_config):
    """Asking for more background than the protein can supply is reported."""
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["few_balls_search"]["min_background_balls"] = 50
    seeds = seeds_file(linear_chain_pdb.parent, [("A", 1)])
    _, result = pipeline(linear_chain_pdb, base_config, seeds)

    # Ten residues, one annotated: nine candidates at most.
    assert result.stats["n_background_balls"] == 9
    assert result.stats["n_background_shortfall"] == 41
    assert any("background balls and placed" in w for w in result.warnings)


# --- what the background may and may not be ---------------------------------

def test_background_never_lands_on_an_annotated_residue(
    linear_chain_pdb, base_config
):
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["few_balls_search"]["min_background_balls"] = 4
    annotated = [("A", i) for i in range(1, 5)]
    seeds = seeds_file(linear_chain_pdb.parent, annotated)
    parsed, result = pipeline(linear_chain_pdb, base_config, seeds)

    ann = seed_labels(parsed, result, SEED_CLASS_ANNOTATED)
    bg = seed_labels(parsed, result, SEED_CLASS_BACKGROUND)
    assert ann == {"A1", "A2", "A3", "A4"}
    assert not (ann & bg)


def test_every_motif_carries_its_class(linear_chain_pdb, base_config):
    """No motif may be unlabelled: the comparison depends on this column."""
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["few_balls_search"]["min_background_balls"] = 3
    seeds = seeds_file(linear_chain_pdb.parent, [("A", 1), ("A", 2), ("A", 3)])
    _, result = pipeline(linear_chain_pdb, base_config, seeds)

    assert result.patches
    assert all(p.seed_class in (SEED_CLASS_ANNOTATED, SEED_CLASS_BACKGROUND)
               for p in result.patches)


def test_both_radii_label_both_classes(linear_chain_pdb, base_config):
    base_config["ball_search"]["radii_angstrom"] = [5.0, 7.0]
    base_config["few_balls_search"]["min_background_balls"] = 3
    seeds = seeds_file(linear_chain_pdb.parent, [("A", 1), ("A", 2), ("A", 3)])
    _, result = pipeline(linear_chain_pdb, base_config, seeds)

    for radius in (5.0, 7.0):
        at_r = [p for p in result.patches if p.radius == radius]
        assert sum(1 for p in at_r if p.seed_class == SEED_CLASS_ANNOTATED) == 3
        assert sum(1 for p in at_r
                   if p.seed_class == SEED_CLASS_BACKGROUND) == 3


# --- the method is defined by an annotated set ------------------------------

def test_no_seed_list_is_an_error_not_a_whole_protein_scan(
    linear_chain_pdb, base_config
):
    """Falling back to every residue would answer a different question."""
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    with pytest.raises(ValueError, match="needs --seed-residues"):
        pipeline(linear_chain_pdb, base_config, None)


def test_annotations_absent_from_the_structure_are_reported(
    linear_chain_pdb, base_config
):
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["few_balls_search"]["min_background_balls"] = 5
    seeds = seeds_file(linear_chain_pdb.parent,
                       [("A", 1), ("A", 2), ("Z", 99)])
    _, result = pipeline(linear_chain_pdb, base_config, seeds)

    assert result.stats["n_annotated_requested"] == 3
    assert result.stats["n_annotated_not_in_structure"] == 1
    assert result.stats["n_annotated_balls"] == 2
    assert any("not\npresent" in w.replace(" ", "\n") or
               "not present in this structure" in w for w in result.warnings)


def test_a_structure_with_no_matching_annotation_says_it_has_no_control(
    linear_chain_pdb, base_config
):
    """Background balls with no annotated class control for nothing."""
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["few_balls_search"]["min_background_balls"] = 3
    seeds = seeds_file(linear_chain_pdb.parent, [("Z", 1)])
    _, result = pipeline(linear_chain_pdb, base_config, seeds)

    assert result.stats["n_annotated_balls"] == 0
    assert result.stats["n_background_balls"] == 3
    assert any("NO annotated class" in w for w in result.warnings)


# --- dispersion and stratification ------------------------------------------

def test_the_background_is_spread_out_not_clustered(
    linear_chain_pdb, base_config
):
    """Three background balls on a 10-residue chain should not be adjacent.

    Residues sit 4 A apart along x, so a clustered choice would give a
    minimum separation of 4 A. Farthest-point sampling must beat that.
    """
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["few_balls_search"]["min_background_balls"] = 3
    seeds = seeds_file(linear_chain_pdb.parent, [("A", 5)])
    _, result = pipeline(linear_chain_pdb, base_config, seeds)

    sep = result.stats["min_separation_between_background_centres_angstrom"]
    assert sep is not None
    assert sep > 4.0


def test_background_is_pushed_away_from_the_annotated_sites(
    linear_chain_pdb, base_config
):
    """The annotated centres are in the chosen set before the first pick."""
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["few_balls_search"]["min_background_balls"] = 1
    seeds = seeds_file(linear_chain_pdb.parent, [("A", 1)])
    parsed, result = pipeline(linear_chain_pdb, base_config, seeds)

    # A1 is at one end, so the single background ball belongs at the other.
    assert seed_labels(parsed, result, SEED_CLASS_BACKGROUND) == {"A10"}
    assert result.stats[
        "min_distance_background_to_annotated_centre_angstrom"] > 30.0


def test_stratification_reaches_buried_positions(buried_core_pdb, base_config):
    """A background set must not be only surface residues.

    Farthest-point sampling alone favours the hull, so without the burial
    strata this is the test that would fail.
    """
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["few_balls_search"]["min_background_balls"] = 4
    base_config["few_balls_search"]["background_stratify_by_burial"] = True
    seeds = seeds_file(buried_core_pdb.parent, [("A", 1)])
    _, result = pipeline(buried_core_pdb, base_config, seeds)

    strata = result.stats["per_stratum"]
    assert strata[STRATUM_BURIED]["n_candidates"] > 0, "fixture has no core"
    assert strata[STRATUM_BURIED]["n_selected"] >= 1
    assert strata[STRATUM_EXPOSED]["n_selected"] >= 1


def test_stratification_can_be_turned_off(buried_core_pdb, base_config):
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["few_balls_search"]["min_background_balls"] = 4
    base_config["few_balls_search"]["background_stratify_by_burial"] = False
    seeds = seeds_file(buried_core_pdb.parent, [("A", 1)])
    _, result = pipeline(buried_core_pdb, base_config, seeds)

    assert result.stats["background_stratified_by_burial"] is False
    assert result.stats["n_background_balls"] == 4


# --- insertion codes ----------------------------------------------------

def test_insertion_code_siblings_do_not_enter_the_walk_unclassified(
    icode_pdb, base_config
):
    """H10, H10A and H10B share an author number but are three residues.

    Selecting seeds by (chain, res_id) would admit all three whenever any one
    was chosen, putting a residue in the walk that belongs to neither class.
    That aborted a real scan on true_annotated_antibody_targets, so it is
    pinned here: every placed ball must carry a class, and the counts must be
    the ones the method chose.
    """
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["few_balls_search"]["min_background_balls"] = 2
    seeds = seeds_file(icode_pdb.parent, [("H", 10)])
    _, result = pipeline(icode_pdb, base_config, seeds)

    assert result.patches
    assert all(p.seed_class in (SEED_CLASS_ANNOTATED, SEED_CLASS_BACKGROUND)
               for p in result.patches)
    counts = classes(result)
    assert counts[SEED_CLASS_ANNOTATED] == result.stats["n_annotated_balls"]
    assert counts[SEED_CLASS_BACKGROUND] == result.stats["n_background_balls"]
    # A seed and its insertion-code siblings are distinct residues: a sibling
    # may be chosen as background, but never silently as a third class.
    assert counts[SEED_CLASS_ANNOTATED] + counts[SEED_CLASS_BACKGROUND] \
        == len(result.patches)


def test_place_centers_refuses_both_seed_selectors():
    from ibc_search.ball_placement import place_centers

    with pytest.raises(ValueError, match="not both"):
        place_centers(None, {}, seed_keys=set(), seed_indices=set())


# --- reproducibility --------------------------------------------------------

def test_the_same_input_gives_the_same_background(
    linear_chain_pdb, base_config
):
    """Selection must not depend on dict order, worker count or chance."""
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["few_balls_search"]["min_background_balls"] = 4
    seeds = seeds_file(linear_chain_pdb.parent, [("A", 1), ("A", 2)])

    parsed_a, first = pipeline(linear_chain_pdb, dict(base_config), seeds)
    parsed_b, second = pipeline(linear_chain_pdb, dict(base_config), seeds)

    assert (seed_labels(parsed_a, first, SEED_CLASS_BACKGROUND)
            == seed_labels(parsed_b, second, SEED_CLASS_BACKGROUND))
    assert ([p.patch_id for p in first.patches]
            == [p.patch_id for p in second.patches])


def test_patch_ids_mark_the_method(linear_chain_pdb, base_config):
    """`f` keeps a few-balls table from colliding with a dumb-search one."""
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["few_balls_search"]["min_background_balls"] = 2
    seeds = seeds_file(linear_chain_pdb.parent, [("A", 1), ("A", 2)])
    _, result = pipeline(linear_chain_pdb, base_config, seeds)

    assert all("_f" in p.patch_id for p in result.patches)


# --- membership is untouched ------------------------------------------------

def test_background_balls_collect_annotated_residues_too(
    linear_chain_pdb, base_config
):
    """No member filter: a ball keeps what is inside it, class regardless."""
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["few_balls_search"]["min_background_balls"] = 9
    seeds = seeds_file(linear_chain_pdb.parent, [("A", 5)])
    parsed, result = pipeline(linear_chain_pdb, base_config, seeds)

    bg = [p for p in result.patches if p.seed_class == SEED_CLASS_BACKGROUND]
    neighbours = {
        parsed.residues[m.residue_index].rid.label()
        for p in bg for m in p.members
    }
    assert "A5" in neighbours, "an annotated residue must still be a member"


def test_overlap_with_annotated_balls_is_measured(
    linear_chain_pdb, base_config
):
    """Overlap is reported rather than prevented, so the choice is visible."""
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["few_balls_search"]["min_background_balls"] = 9
    seeds = seeds_file(linear_chain_pdb.parent, [("A", 5)])
    _, result = pipeline(linear_chain_pdb, base_config, seeds)

    # A4 and A6 are inside A5's 5 A ball and are eligible background seeds.
    assert result.stats["n_background_seeds_inside_an_annotated_ball"] >= 1


# --- the two helpers, directly ----------------------------------------------

def test_allocate_sums_to_the_total_and_is_proportional():
    strata = {"a": list(range(90)), "b": list(range(10)), "c": []}
    alloc = _allocate(strata, 10)
    assert sum(alloc.values()) == 10
    assert alloc["c"] == 0
    assert alloc["a"] > alloc["b"] > 0


def test_allocate_gives_a_small_stratum_at_least_one():
    """A 96%-buried protein must still get an exposed background ball."""
    strata = {"buried": list(range(960)), "exposed": list(range(40))}
    alloc = _allocate(strata, 5)
    assert sum(alloc.values()) == 5
    assert alloc["exposed"] >= 1


def test_allocate_cannot_exceed_a_stratum_or_the_pool():
    strata = {"a": [1, 2], "b": [3]}
    alloc = _allocate(strata, 10)
    assert alloc == {"a": 2, "b": 1}
    assert _allocate({"a": []}, 5) == {"a": 0}
    assert _allocate({"a": [1]}, 0) == {"a": 0}


def test_farthest_point_picks_the_extremes():
    """Four points on a line, staying away from an occupied origin.

    x = 0, 1, 2, 10 with the origin already taken. The first pick is x = 10,
    the furthest from the origin. The second maximises the distance to the
    NEAREST chosen point, which is x = 2 (2 A from the origin, 8 A from
    x = 10) rather than x = 0 or x = 1, both of which sit on top of the
    origin. Picking an end here would mean the objective was distance to the
    last pick instead of to the chosen set.
    """
    pts = np.array([[0.0, 0, 0], [1.0, 0, 0], [2.0, 0, 0], [10.0, 0, 0]])
    got = _farthest_point_select(pts, [0, 1, 2, 3], 2, np.array([[0.0, 0, 0]]))
    assert got == [3, 2]


def test_farthest_point_is_deterministic_under_ties():
    """Symmetric input must resolve to the lower index, every time."""
    pts = np.array([[-1.0, 0, 0], [1.0, 0, 0]])
    chosen = np.array([[0.0, 0.0, 0.0]])
    runs = {tuple(_farthest_point_select(pts, [0, 1], 1, chosen))
            for _ in range(5)}
    assert runs == {(0,)}


def test_farthest_point_respects_an_empty_request():
    pts = np.array([[0.0, 0, 0]])
    assert _farthest_point_select(pts, [0], 0, None) == []
    assert _farthest_point_select(pts, [], 3, None) == []
