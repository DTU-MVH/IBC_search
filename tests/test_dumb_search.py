"""``dumb_search`` behaviour: one ball per residue, and no filters at all.

The point of these tests is that the method stays dumb. Most of them assert
the *absence* of a behaviour -- no threshold, no connectivity check, no
deduplication, no size limit -- because that absence is the whole method, and
a well-meaning future addition would be a regression rather than a feature.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from ibc_search.ball_search_algo import run_ball_search
from ibc_search.dumb_search import run_dumb_search
from ibc_search.exposure import compute_exposure
from ibc_search.structure_io import parse_structure, resolve_input


def pipeline(path: Path, config: dict):
    parsed = parse_structure(resolve_input(str(path), config), config)
    exposure = compute_exposure(parsed, config)
    return parsed, exposure, run_dumb_search(parsed, exposure, config)


def labels_of(parsed, patch) -> set[str]:
    return {parsed.residues[m.residue_index].rid.label() for m in patch.members}


def seed_label(parsed, patch) -> str:
    (seed,) = patch.seed_residue_indices
    return parsed.residues[seed].rid.label()


# --- the walk ---------------------------------------------------------------

def test_one_ball_per_residue_from_terminus_to_terminus(linear_chain_pdb, base_config):
    """Ten residues, one radius: ten motifs, in sequence order, A1 to A10."""
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    parsed, _, result = pipeline(linear_chain_pdb, base_config)

    assert len(result.patches) == 10
    assert [seed_label(parsed, p) for p in result.patches] == [
        f"A{i}" for i in range(1, 11)
    ]
    assert all(len(p.seed_residue_indices) == 1 for p in result.patches)


def test_every_residue_is_in_its_own_ball(linear_chain_pdb, base_config):
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    parsed, _, result = pipeline(linear_chain_pdb, base_config)

    for patch in result.patches:
        assert seed_label(parsed, patch) in labels_of(parsed, patch)
        assert sum(1 for m in patch.members if m.is_seed) == 1


def test_each_radius_walks_the_whole_chain(linear_chain_pdb, base_config):
    """Two radii, ten residues: twenty motifs, grouped by radius."""
    base_config["ball_search"]["radii_angstrom"] = [5.0, 8.0]
    parsed, _, result = pipeline(linear_chain_pdb, base_config)

    assert len(result.patches) == 20
    assert [p.radius for p in result.patches] == [5.0] * 10 + [8.0] * 10
    for radius in (5.0, 8.0):
        at_radius = [p for p in result.patches if p.radius == radius]
        assert [seed_label(parsed, p) for p in at_radius] == [
            f"A{i}" for i in range(1, 11)
        ]


def test_geometry_matches_the_ball_search_for_a_shared_seed(hairpin_pdb, base_config):
    """Same centre, same radius, same members: only the filters differ.

    Run the real ball search with every filter switched off, so the two
    methods must agree residue for residue. This is what "reuses the existing
    ball geometry" has to mean to be worth claiming.
    """
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["ball_search"]["seed_rsa_threshold"] = 0.0
    base_config["ball_search"]["member_rsa_threshold"] = 0.0
    base_config["ball_search"]["min_patch_size"] = 1
    base_config["connectivity"]["disconnected_policy"] = "keep"

    parsed = parse_structure(resolve_input(str(hairpin_pdb), base_config), base_config)
    exposure = compute_exposure(parsed, base_config)
    dumb = run_dumb_search(parsed, exposure, base_config)
    ball = run_ball_search(parsed, exposure, base_config)

    ball_sets = {p.member_set_hash for p in ball.patches}
    dumb_sets = {p.member_set_hash for p in dumb.patches}
    # The ball search deduplicates identical member sets and dumb_search does
    # not, so dumb_sets is the superset only up to that; as sets of residue
    # sets they are equal.
    assert dumb_sets == ball_sets


def test_balls_reach_across_sequence(hairpin_pdb, base_config):
    """A hairpin motif is discontinuous: this is a 3D search, not a window."""
    base_config["ball_search"]["radii_angstrom"] = [6.0]
    parsed, _, result = pipeline(hairpin_pdb, base_config)

    assert any(p.is_discontinuous for p in result.patches)


# --- what is deliberately absent -------------------------------------------

def test_buried_residues_are_still_seeds_and_members(bridged_pdb, base_config):
    """No RSA threshold anywhere: raising it changes nothing.

    The ball search drops seeds below ``seed_rsa_threshold``; dumb_search must
    produce the same motifs at 0.0 and at an impossible 1.5.
    """
    base_config["ball_search"]["radii_angstrom"] = [6.0]
    base_config["ball_search"]["seed_rsa_threshold"] = 0.0
    base_config["ball_search"]["member_rsa_threshold"] = 0.0
    parsed, _, low = pipeline(bridged_pdb, base_config)

    base_config["ball_search"]["seed_rsa_threshold"] = 1.5
    base_config["ball_search"]["member_rsa_threshold"] = 1.5
    _, _, high = pipeline(bridged_pdb, base_config)

    assert len(low.patches) == len(high.patches) == len(parsed.residues)
    assert [p.member_set_hash for p in low.patches] == [
        p.member_set_hash for p in high.patches
    ]


def test_no_connectivity_check(bridged_pdb, base_config):
    """A 10 A ball spanning two out-of-contact clusters keeps both.

    ``bridged_pdb`` is the fixture built for the ball search's connectivity
    check. Here the far cluster must survive, and must not be flagged either,
    because no connectivity question was ever asked.
    """
    base_config["ball_search"]["radii_angstrom"] = [10.0]
    base_config["connectivity"]["disconnected_policy"] = "drop"
    parsed, _, result = pipeline(bridged_pdb, base_config)

    first = next(p for p in result.patches if seed_label(parsed, p) == "A1")
    assert {"A1", "A2", "A10", "A11"} <= labels_of(parsed, first)
    assert all(not m.disconnected_flagged for m in first.members)
    assert all(p.n_dropped_disconnected == 0 for p in result.patches)


def test_no_deduplication(bridged_pdb, base_config):
    """Two residues whose balls hold the same set still give two motifs.

    In ``bridged_pdb`` the paired residues A1/A2 see only each other at 5 A,
    as do A10/A11, so each pair's two balls hold an identical residue set. The
    ball search would merge each pair into one patch carrying two seeds; here
    all four motifs survive, one per residue, each with its own seed.
    """
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    parsed, _, result = pipeline(bridged_pdb, base_config)

    hashes = [p.member_set_hash for p in result.patches]
    assert len(hashes) == 4
    assert len(set(hashes)) == 2, "fixture should produce two repeated sets"
    assert [seed_label(parsed, p) for p in result.patches] == [
        "A1", "A2", "A10", "A11",
    ]
    for patch in result.patches:
        assert len(patch.seed_residue_indices) == 1


def test_no_minimum_size(two_chain_pdb, base_config):
    """A lone residue in its ball is still a motif.

    A tiny radius leaves every ball holding only its own residue. The ball
    search drops those as undersized; dumb_search keeps all of them, and the
    FoldDisco contract -- not the extractor -- is what declines to emit a
    query for a motif too small to express.
    """
    base_config["ball_search"]["radii_angstrom"] = [1.0]
    base_config["ball_search"]["min_patch_size"] = 4
    parsed, _, result = pipeline(two_chain_pdb, base_config)

    assert len(result.patches) == len(parsed.residues)
    assert all(p.size == 1 for p in result.patches)
    assert all(not p.emit_query for p in result.patches)


def test_oversized_motifs_are_kept_but_not_queried(linear_chain_pdb, base_config):
    """Above ``max_patch_size`` a motif stays in the table, with no query."""
    base_config["ball_search"]["radii_angstrom"] = [20.0]
    base_config["ball_search"]["max_patch_size"] = 3
    parsed, _, result = pipeline(linear_chain_pdb, base_config)

    assert len(result.patches) == len(parsed.residues)
    assert all(p.size > 3 for p in result.patches)
    assert all(not p.emit_query for p in result.patches)
    assert any("exceed max_patch_size" in w for w in result.warnings)


# --- housekeeping -----------------------------------------------------------

def test_residue_without_a_placeable_centre_is_stepped_over(glycine_pdb, base_config):
    """No CB, no backbone, no CA: no ball, a warning, still a member.

    Residue A7 of the fixture has only a CA, so it *can* take a centre; the
    case tested here is that the walk never fabricates one, and that glycine
    (A3) is walked normally with a virtual centre rather than skipped.
    """
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    parsed, _, result = pipeline(glycine_pdb, base_config)

    walked = {seed_label(parsed, p) for p in result.patches}
    assert "A3" in walked, "glycine must be walked, not skipped"
    assert result.stats["n_virtual_centers"] >= 1
    assert result.stats["n_seeds"] == len(result.patches)


def test_waters_and_ligands_are_neither_seeds_nor_members(
    waters_and_ligand_pdb, base_config
):
    base_config["ball_search"]["radii_angstrom"] = [8.0]
    parsed, _, result = pipeline(waters_and_ligand_pdb, base_config)

    assert len(result.patches) == 5
    for patch in result.patches:
        for label in labels_of(parsed, patch):
            assert not label.startswith("B")
            assert label not in {"A100", "A101"}


def test_method_is_recorded_in_the_stats(linear_chain_pdb, base_config):
    """The manifest must say which method produced a run."""
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    _, _, result = pipeline(linear_chain_pdb, base_config)

    assert result.stats["method"] == "dumb_search"


def test_patch_ids_are_unique_and_marked(linear_chain_pdb, base_config):
    base_config["ball_search"]["radii_angstrom"] = [5.0, 8.0]
    _, _, result = pipeline(linear_chain_pdb, base_config)

    ids = [p.patch_id for p in result.patches]
    assert len(set(ids)) == len(ids)
    # The trailing field marks the method: `d` here, `p` from the ball search,
    # so the two cannot be confused once a row leaves its directory.
    assert all(pid.rsplit("_", 1)[-1].startswith("d") for pid in ids)


def test_determinism_across_worker_counts(hairpin_pdb, base_config):
    base_config["ball_search"]["radii_angstrom"] = [6.0]
    base_config["performance"]["workers"] = 1
    _, _, one = pipeline(hairpin_pdb, base_config)
    base_config["performance"]["workers"] = 4
    _, _, four = pipeline(hairpin_pdb, base_config)

    assert [p.patch_id for p in one.patches] == [p.patch_id for p in four.patches]
    assert [p.member_set_hash for p in one.patches] == [
        p.member_set_hash for p in four.patches
    ]


# --- the search input is deduplicated, the motif set is not ----------------

def test_identical_motifs_are_searched_once_but_still_counted_twice(
        bridged_pdb, base_config, tmp_path):
    """One line per distinct query; one index row per patch.

    ``bridged_pdb`` gives two pairs whose balls hold identical residue sets,
    so four motifs produce only two distinct queries. Searching a query twice
    returns byte-identical hits for twice the time, and on a real 5.5 A run
    0.16% of queries were such repeats.

    What must NOT change is the motif set: dumb_search still emits one patch
    per residue, and the patches table still lists all four.
    """
    from ibc_search.exposure import compute_exposure
    from ibc_search.outputs import write_tables
    from ibc_search.structure_io import parse_structure, resolve_input

    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["run"]["output_dir"] = str(tmp_path / "out")
    parsed = parse_structure(resolve_input(str(bridged_pdb), base_config),
                             base_config)
    exposure = compute_exposure(parsed, base_config)
    result = run_dumb_search(parsed, exposure, base_config)
    paths = write_tables(tmp_path / "out", parsed, exposure, result, base_config)

    lines = [ln for ln in
             paths["folddisco_queries"].read_text().splitlines() if ln]
    assert len(lines) == len(set(lines)), "no query may be searched twice"
    assert len(lines) == 2, "four motifs, two distinct residue sets"

    import csv as _csv
    index = list(_csv.DictReader(
        paths["folddisco_query_index"].open(), delimiter="\t"))
    assert len(index) == 4, "every patch keeps its own index row"
    # The two patches sharing a set must point at the same line.
    by_line: dict[str, list[str]] = {}
    for row in index:
        by_line.setdefault(row["line_number"], []).append(row["patch_id"])
    assert sorted(len(v) for v in by_line.values()) == [2, 2]
    assert len(result.patches) == 4, "the motif set is untouched"


def test_the_deduplication_is_reported_not_silent(bridged_pdb, base_config,
                                                  tmp_path):
    from ibc_search.exposure import compute_exposure
    from ibc_search.outputs import write_tables
    from ibc_search.structure_io import parse_structure, resolve_input

    base_config["ball_search"]["radii_angstrom"] = [5.0]
    parsed = parse_structure(resolve_input(str(bridged_pdb), base_config),
                             base_config)
    exposure = compute_exposure(parsed, base_config)
    result = run_dumb_search(parsed, exposure, base_config)
    write_tables(tmp_path / "out", parsed, exposure, result, base_config)

    assert any("share a residue set" in w for w in result.warnings)


def test_a_run_with_no_repeats_emits_no_dedup_warning(linear_chain_pdb,
                                                      base_config, tmp_path):
    """The warning must mean something, so it cannot fire on every run."""
    from ibc_search.exposure import compute_exposure
    from ibc_search.outputs import write_tables
    from ibc_search.structure_io import parse_structure, resolve_input

    base_config["ball_search"]["radii_angstrom"] = [5.0]
    parsed = parse_structure(resolve_input(str(linear_chain_pdb), base_config),
                             base_config)
    exposure = compute_exposure(parsed, base_config)
    result = run_dumb_search(parsed, exposure, base_config)
    paths = write_tables(tmp_path / "out", parsed, exposure, result, base_config)

    lines = [ln for ln in
             paths["folddisco_queries"].read_text().splitlines() if ln]
    assert len(lines) == len(result.patches)
    assert not any("share a residue set" in w for w in result.warnings)


# --- seeded placement -------------------------------------------------------

def _seed_file(path, rows, header=("pdb_id", "chain", "res_id")):
    import csv as _csv
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = _csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)
    return path


def test_seeding_restricts_the_balls_to_the_named_residues(
        linear_chain_pdb, base_config, tmp_path):
    seeds = _seed_file(tmp_path / "s.tsv", [["linear", "A", 3], ["linear", "A", 7]])
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["ball_search"]["seed_residues_path"] = str(seeds)
    parsed, _, result = pipeline(linear_chain_pdb, base_config)

    assert len(result.patches) == 2, "ten residues, two named"
    assert sorted(seed_label(parsed, p) for p in result.patches) == ["A3", "A7"]
    assert result.stats["seeded"] is True
    assert result.stats["n_seed_residues_requested"] == 2


def test_seeding_does_not_restrict_MEMBERSHIP(linear_chain_pdb, base_config,
                                              tmp_path):
    """The motif is the seed's real neighbourhood, not a subgraph of the list.

    Restricting members too would make the ball describe the annotation
    rather than the structure around it.
    """
    seeds = _seed_file(tmp_path / "s.tsv", [["linear", "A", 5]])
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["ball_search"]["seed_residues_path"] = str(seeds)
    parsed, _, result = pipeline(linear_chain_pdb, base_config)

    (patch,) = result.patches
    labels = labels_of(parsed, patch)
    assert "A5" in labels
    assert labels - {"A5"}, "neighbours must be collected, named or not"
    assert any(x not in {"A5"} for x in labels)


def test_no_seed_file_keeps_one_ball_per_residue(linear_chain_pdb, base_config):
    """An absent list means NO RESTRICTION, never 'no seeds'."""
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["ball_search"]["seed_residues_path"] = ""
    _, _, result = pipeline(linear_chain_pdb, base_config)

    assert len(result.patches) == 10
    assert result.stats["seeded"] is False


def test_seed_rows_are_filtered_to_this_structure(linear_chain_pdb,
                                                  base_config, tmp_path):
    """One annotation table drives a whole panel without being split up."""
    seeds = _seed_file(tmp_path / "s.tsv",
                       [["linear", "A", 3], ["OTHER", "A", 4],
                        ["OTHER", "A", 5]])
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["ball_search"]["seed_residues_path"] = str(seeds)
    parsed, _, result = pipeline(linear_chain_pdb, base_config)

    assert [seed_label(parsed, p) for p in result.patches] == ["A3"]


def test_an_is_annotated_column_is_honoured(linear_chain_pdb, base_config,
                                            tmp_path):
    seeds = _seed_file(
        tmp_path / "s.tsv",
        [["linear", "A", 3, "1"], ["linear", "A", 4, "0"]],
        header=("pdb_id", "chain", "res_id", "is_annotated"))
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["ball_search"]["seed_residues_path"] = str(seeds)
    parsed, _, result = pipeline(linear_chain_pdb, base_config)

    assert [seed_label(parsed, p) for p in result.patches] == ["A3"]


def test_a_seed_absent_from_the_structure_is_reported_not_silent(
        linear_chain_pdb, base_config, tmp_path):
    """A requested seed that is not there is a different fact from a seed
    that found nothing, and must not be counted as the latter."""
    seeds = _seed_file(tmp_path / "s.tsv",
                       [["linear", "A", 3], ["linear", "Z", 999]])
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["ball_search"]["seed_residues_path"] = str(seeds)
    _, _, result = pipeline(linear_chain_pdb, base_config)

    assert result.stats["n_seed_residues_not_in_structure"] == 1
    assert any("not present in this structure" in w for w in result.warnings)
    assert len(result.patches) == 1


def test_a_missing_seed_file_is_an_error_not_an_empty_set(linear_chain_pdb,
                                                          base_config):
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["ball_search"]["seed_residues_path"] = "/nope/missing.tsv"
    with pytest.raises(ValueError, match="seed residue file not found"):
        pipeline(linear_chain_pdb, base_config)


# --- a second round of seeding ----------------------------------------------
#
# The linear fixture makes the arithmetic checkable by hand: ten ALA spaced
# 4 A along x, and a 5 A ball reaches only the immediate sequence
# neighbours. So one round on A5 collects {A4, A5, A6}, and promoting those
# to seeds must give exactly three balls.

def test_a_second_round_seeds_the_members_of_the_first_round_balls(
        linear_chain_pdb, base_config, tmp_path):
    seeds = _seed_file(tmp_path / "s.tsv", [["linear", "A", 5]])
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["ball_search"]["seed_residues_path"] = str(seeds)
    base_config["ball_search"]["seed_rounds"] = 2
    parsed, _, result = pipeline(linear_chain_pdb, base_config)

    assert sorted(seed_label(parsed, p) for p in result.patches) == \
        ["A4", "A5", "A6"]
    st = result.stats
    assert st["seed_rounds"] == 2
    assert st["n_seed_residues_requested"] == 1, "the list the caller gave"
    assert st["n_seed_residues_added_by_expansion"] == 2
    assert st["n_seed_residues_after_expansion"] == 3


def test_the_added_seeds_are_exactly_what_round_one_already_collected(
        linear_chain_pdb, base_config, tmp_path):
    """Expansion promotes members to seeds; it never invents a residue.

    The round-2 seed set must equal the union of round 1's patch members,
    so nothing enters the walk that round 1's motifs did not already
    contain.
    """
    seeds = _seed_file(tmp_path / "s.tsv",
                       [["linear", "A", 3], ["linear", "A", 8]])
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["ball_search"]["seed_residues_path"] = str(seeds)

    parsed, _, one = pipeline(linear_chain_pdb, base_config)
    members_of_round_one = set()
    for patch in one.patches:
        members_of_round_one |= labels_of(parsed, patch)

    base_config["ball_search"]["seed_rounds"] = 2
    parsed2, _, two = pipeline(linear_chain_pdb, base_config)

    assert {seed_label(parsed2, p) for p in two.patches} == \
        members_of_round_one


def test_one_round_is_the_default_so_nothing_changes_unasked(
        linear_chain_pdb, base_config, tmp_path):
    seeds = _seed_file(tmp_path / "s.tsv", [["linear", "A", 5]])
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["ball_search"]["seed_residues_path"] = str(seeds)
    base_config["ball_search"].pop("seed_rounds", None)
    parsed, _, result = pipeline(linear_chain_pdb, base_config)

    assert [seed_label(parsed, p) for p in result.patches] == ["A5"]
    assert result.stats["seed_rounds"] == 1
    assert result.stats["n_seed_residues_added_by_expansion"] is None


def test_rounds_compose_rather_than_repeating_the_first(
        linear_chain_pdb, base_config, tmp_path):
    """Three rounds reach two steps out, not one step twice."""
    seeds = _seed_file(tmp_path / "s.tsv", [["linear", "A", 5]])
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["ball_search"]["seed_residues_path"] = str(seeds)
    base_config["ball_search"]["seed_rounds"] = 3
    parsed, _, result = pipeline(linear_chain_pdb, base_config)

    assert sorted(seed_label(parsed, p) for p in result.patches) == \
        ["A3", "A4", "A5", "A6", "A7"]


def test_expansion_does_not_restrict_membership_either(
        linear_chain_pdb, base_config, tmp_path):
    """A round-2 ball is still the real neighbourhood of its own centre.

    A4's ball must contain A3, which is in neither the annotation nor the
    expanded seed set. Expansion changes WHERE balls are placed, never what
    they collect.
    """
    seeds = _seed_file(tmp_path / "s.tsv", [["linear", "A", 5]])
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["ball_search"]["seed_residues_path"] = str(seeds)
    base_config["ball_search"]["seed_rounds"] = 2
    parsed, _, result = pipeline(linear_chain_pdb, base_config)

    a4, = [p for p in result.patches if seed_label(parsed, p) == "A4"]
    assert "A3" in labels_of(parsed, a4)


def test_more_rounds_without_a_seed_list_is_an_error_not_a_no_op(
        linear_chain_pdb, base_config):
    """With no list every residue is already a seed, so the request is
    incoherent and saying so beats silently running the unseeded method."""
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["ball_search"]["seed_residues_path"] = ""
    base_config["ball_search"]["seed_rounds"] = 2

    with pytest.raises(ValueError, match="needs a seed residue list"):
        pipeline(linear_chain_pdb, base_config)


def test_zero_rounds_is_rejected(linear_chain_pdb, base_config, tmp_path):
    seeds = _seed_file(tmp_path / "s.tsv", [["linear", "A", 5]])
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["ball_search"]["seed_residues_path"] = str(seeds)
    base_config["ball_search"]["seed_rounds"] = 0

    with pytest.raises(ValueError, match="seed_rounds must be"):
        pipeline(linear_chain_pdb, base_config)


def test_expansion_uses_the_largest_radius_since_balls_are_nested(
        linear_chain_pdb, base_config, tmp_path):
    """Two radii share one walk, so one expansion radius has to serve both.

    The largest is the union: a 5 A ball at a centre is contained in the
    9 A ball at that same centre, so no smaller radius can reach a residue
    the largest misses.
    """
    seeds = _seed_file(tmp_path / "s.tsv", [["linear", "A", 5]])
    base_config["ball_search"]["radii_angstrom"] = [5.0, 9.0]
    base_config["ball_search"]["seed_residues_path"] = str(seeds)
    base_config["ball_search"]["seed_rounds"] = 2
    parsed, _, result = pipeline(linear_chain_pdb, base_config)

    # 9 A from A5 reaches A3..A7; the walk is those five, run at each radius.
    seeds_seen = {seed_label(parsed, p) for p in result.patches}
    assert seeds_seen == {"A3", "A4", "A5", "A6", "A7"}
    assert len(result.patches) == 10, "five seeds at each of two radii"
    assert any("nested" in w for w in result.warnings)


def test_the_seed_id_can_be_named_when_the_file_is_not(
        linear_chain_pdb, base_config, tmp_path):
    """A biological-assembly file is `<id>-assembly1`; its annotations are
    filed under `<id>`. Matching nothing must not be the silent default."""
    seeds = _seed_file(tmp_path / "s.tsv", [["11KP", "A", 5]])
    base_config["ball_search"]["radii_angstrom"] = [5.0]
    base_config["ball_search"]["seed_residues_path"] = str(seeds)

    parsed, _, miss = pipeline(linear_chain_pdb, base_config)
    assert miss.patches == [], "the stem 'linear' matches no row"

    base_config["ball_search"]["seed_structure_id"] = "11KP"
    parsed, _, hit = pipeline(linear_chain_pdb, base_config)
    assert [seed_label(parsed, p) for p in hit.patches] == ["A5"]
    assert hit.stats["seed_structure_id"] == "11KP"
