"""The within-structure annotated-vs-background verdict.

The verdict is a single bit per annotated residue, so almost everything that
could go wrong with it is silent: comparing against the wrong structure's
background, counting a tie as a win, turning an unanswerable comparison into
a 0, or mixing radii that are not comparable. Each of those is pinned here.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from annotated_vs_background import (  # noqa: E402
    SEED_CLASS_ANNOTATED,
    SEED_CLASS_BACKGROUND,
    annotated_beats_background,
    statistics,
)


def rec(inst, patch_id, seed_class, n, radius=5.5, seed="A1"):
    return {"instance": inst, "radius": radius, "patch_id": patch_id,
            "seed_residue": seed, "seed_class": seed_class,
            "n_distinct_sites": n}


def ann(inst, pid, n, **kw):
    return rec(inst, pid, SEED_CLASS_ANNOTATED, n, **kw)


def bg(inst, pid, n, **kw):
    return rec(inst, pid, SEED_CLASS_BACKGROUND, n, **kw)


# --- the comparison is within one structure ---------------------------------

def test_the_background_compared_against_is_the_same_structures():
    """A strong background elsewhere must not decide another structure."""
    records = [
        ann("AAA", "p1", 50), bg("AAA", "p2", 10), bg("AAA", "p3", 20),
        ann("BBB", "p4", 50), bg("BBB", "p5", 900), bg("BBB", "p6", 800),
    ]
    out = {r["instance"]: r["beats_background"]
           for r in annotated_beats_background(records)}
    assert out == {"AAA": 1, "BBB": 0}


def test_each_radius_is_judged_separately():
    """Counts at different radii are different scales and never pooled."""
    records = [
        ann("AAA", "p1", 50, radius=5.0), bg("AAA", "p2", 10, radius=5.0),
        ann("AAA", "p3", 50, radius=8.0), bg("AAA", "p4", 900, radius=8.0),
    ]
    out = {r["radius_angstrom"]: r["beats_background"]
           for r in annotated_beats_background(records)}
    assert out == {5.0: 1, 8.0: 0}


# --- the rule itself --------------------------------------------------------

def test_the_test_is_strictly_greater_so_a_tie_is_not_a_win():
    records = [ann("A", "p1", 10), bg("A", "p2", 10), bg("A", "p3", 10)]
    (row,) = annotated_beats_background(records)
    assert row["background_median"] == 10
    assert row["beats_background"] == 0


def test_the_reference_is_the_median_not_the_mean_by_default():
    """One huge background ball must not decide the whole structure."""
    records = [ann("A", "p1", 50),
               bg("A", "p2", 10), bg("A", "p3", 10), bg("A", "p4", 10_000)]
    (by_median,) = annotated_beats_background(records)
    (by_mean,) = annotated_beats_background(records, reference="mean")
    assert by_median["beats_background"] == 1
    assert by_mean["beats_background"] == 0


def test_mean_reference_can_be_asked_for_explicitly():
    records = [ann("A", "p1", 50), bg("A", "p2", 10), bg("A", "p3", 20)]
    (row,) = annotated_beats_background(records, reference="mean")
    assert row["background_mean"] == 15
    assert row["beats_background"] == 1


def test_an_unknown_reference_is_rejected():
    with pytest.raises(ValueError, match="median.*mean"):
        annotated_beats_background([], reference="mode")


# --- unanswerable is not the same as failed ---------------------------------

def test_no_background_gives_no_verdict_rather_than_a_zero():
    """A structure with no control answers nothing; 0 would be a claim."""
    records = [ann("A", "p1", 50), ann("A", "p2", 60)]
    rows = annotated_beats_background(records)
    assert len(rows) == 2
    assert all(r["beats_background"] is None for r in rows)
    assert all(r["n_background_balls"] == 0 for r in rows)
    assert all(r["percentile_within_background"] is None for r in rows)


def test_rows_without_a_verdict_are_excluded_from_the_rate():
    records = [ann("A", "p1", 50),                       # no background
               ann("B", "p2", 50), bg("B", "p3", 10)]    # beats
    verdicts = annotated_beats_background(records)
    stats = statistics(records, verdicts)
    assert stats["pooled_annotated_residues"] == 2
    assert stats["pooled_annotated_with_a_verdict"] == 1
    assert stats["pooled_fraction_beating"] == 1.0


# --- the percentile tells a near miss from a clear one ----------------------

def test_the_percentile_distinguishes_a_near_miss():
    records = [ann("A", "p1", 9),
               bg("A", "p2", 1), bg("A", "p3", 10), bg("A", "p4", 11)]
    (row,) = annotated_beats_background(records)
    assert row["beats_background"] == 0
    # Beat one of three background balls: not a win, but not bottom either.
    assert row["percentile_within_background"] == pytest.approx(1 / 3, abs=1e-4)


def test_the_percentile_is_one_when_it_beats_everything():
    records = [ann("A", "p1", 99), bg("A", "p2", 1), bg("A", "p3", 2)]
    (row,) = annotated_beats_background(records)
    assert row["percentile_within_background"] == 1.0
    assert row["beats_background"] == 1


# --- statistics -------------------------------------------------------------

def test_statistics_report_mean_variation_and_distinct_counts():
    records = [ann("A", "p1", 10), ann("A", "p2", 20),
               bg("A", "p3", 1), bg("A", "p4", 1)]
    stats = statistics(records, annotated_beats_background(records))
    a = stats["by_seed_class"][SEED_CLASS_ANNOTATED]
    b = stats["by_seed_class"][SEED_CLASS_BACKGROUND]
    assert a["n"] == 2 and a["mean"] == 15.0
    assert a["sd"] == pytest.approx(7.0711, abs=1e-3)
    assert a["n_distinct_values"] == 2
    # Two balls with the SAME count are one distinct value, not two.
    assert b["n_distinct_values"] == 1


def test_per_structure_fractions_are_reported_beside_the_pooled_one():
    """One structure with many annotations must not set the headline alone."""
    records = [
        # 4 annotated, all winning
        *[ann("BIG", f"a{i}", 90) for i in range(4)], bg("BIG", "b", 1),
        # 1 annotated, losing
        ann("SMALL", "a9", 1), bg("SMALL", "b9", 90),
    ]
    stats = statistics(records, annotated_beats_background(records))
    assert stats["pooled_fraction_beating"] == pytest.approx(4 / 5)
    # Per structure the two count equally: one at 1.0, one at 0.0.
    assert stats["mean_of_per_structure_fractions"] == pytest.approx(0.5)
    assert stats["n_structures_with_a_verdict"] == 2


def test_annotated_are_split_by_verdict_in_the_statistics():
    records = [ann("A", "p1", 90), ann("A", "p2", 1), bg("A", "p3", 50)]
    stats = statistics(records, annotated_beats_background(records))
    assert stats["annotated_by_verdict"]["beats_background"]["n"] == 1
    assert stats["annotated_by_verdict"]["does_not_beat"]["n"] == 1


# --- reading a run ----------------------------------------------------------

def test_unlabelled_patches_are_refused(tmp_path):
    """A run whose balls carry no class cannot be analysed this way."""
    from annotated_vs_background import motif_scores

    d = tmp_path / "XXXX"
    d.mkdir()
    with (d / "patches.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["patch_id", "radius_angstrom", "seed_class",
                    "seed_residue_labels", "emit_query"])
        w.writerow(["XXXX_m1_r5.5_d00000", "5.5", "", "A1", "true"])
    with pytest.raises(ValueError, match="few_balls_search"):
        motif_scores(tmp_path, "XXXX")


def test_balls_with_no_query_are_absent_not_zero(tmp_path):
    """Never assessed is a different fact from assessed and empty."""
    from annotated_vs_background import motif_scores

    d = tmp_path / "XXXX"
    d.mkdir()
    with (d / "patches.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["patch_id", "radius_angstrom", "seed_class",
                    "seed_residue_labels", "emit_query"])
        w.writerow(["p1", "5.5", SEED_CLASS_ANNOTATED, "A1", "true"])
        w.writerow(["p2", "5.5", SEED_CLASS_ANNOTATED, "A2", "false"])
    got = motif_scores(tmp_path, "XXXX")
    assert [g["patch_id"] for g in got] == ["p1"]
    # Searched and matched nothing IS zero, because it was assessed.
    assert got[0]["n_distinct_sites"] == 0
