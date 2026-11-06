"""Tests for counting distinct hits per residue.

The behaviours worth pinning are the ones that fail silently and still
produce a plausible-looking number:

1. a site reached through several of a residue's overlapping motifs counts
   **once** for that residue, so the count measures the database rather than
   how many motifs happened to contain the residue;
2. matched tokens are compared as a set, so token order cannot split one site
   into two, and ``_`` placeholders never become part of a site's identity;
3. a searched motif that matched nothing scores **0**, while a residue that
   was never queried is **absent** -- the two mean different things;
4. nothing anywhere is z-scored.
"""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import pytest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "count_scoring.py"
sys.path.insert(0, str(SCRIPT.parent))

from count_scoring import (  # noqa: E402
    auc,
    residue_hit_counts,
    residue_motif_counts,
    site_key,
)


def _write(path: Path, header: list[str], rows: list[list]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def _run(root: Path, name: str, hits: list[list], *,
         members: list[list] | None = None,
         queried: list[str] | None = None) -> Path:
    """A results directory with two overlapping motifs sharing residue A2.

    p1 = {A1, A2}, p2 = {A2, A3}. Residue A2 is the interesting one: whatever
    both motifs find, it must count once for A2.
    """
    d = root / name
    _write(d / "patch_members.tsv",
           ["patch_id", "auth_chain_id", "auth_seq_id"],
           members or [["p1", "A", 1], ["p1", "A", 2],
                       ["p2", "A", 2], ["p2", "A", 3]])
    _write(d / "folddisco_queries.index.tsv",
           ["patch_id", "query_string"],
           [[p, "A1,A2"] for p in (queried if queried is not None else ["p1", "p2"])])
    _write(d / "folddisco_hits.tsv",
           ["patch_id", "tid", "idf", "rmsd", "query_residues",
            "matching_residues"],
           hits)
    return d


# --- deduplication ----------------------------------------------------------

def test_the_same_site_through_two_motifs_counts_once(tmp_path):
    """The central rule: one site, reached twice, is one hit for the residue."""
    d = _run(tmp_path, "m", [
        ["p1", "pdb/aa/1abc.ent.gz", "100", "0.5", "A1,A2", "B10,B11"],
        ["p2", "pdb/aa/1abc.ent.gz", "100", "0.5", "A2,A3", "B10,B11"],
    ])
    counts = residue_hit_counts(d)

    assert counts[("A", 2)] == 1, "A2 is in both motifs; the site is still one"
    assert counts[("A", 1)] == 1
    assert counts[("A", 3)] == 1


def test_different_sites_in_the_same_structure_count_separately(tmp_path):
    d = _run(tmp_path, "m", [
        ["p1", "pdb/aa/1abc.ent.gz", "100", "0.5", "A1,A2", "B10,B11"],
        ["p2", "pdb/aa/1abc.ent.gz", "100", "0.5", "A2,A3", "B40,B41"],
    ])
    counts = residue_hit_counts(d)

    assert counts[("A", 2)] == 2, "same structure, two places, two hits"
    assert counts[("A", 1)] == 1
    assert counts[("A", 3)] == 1


def test_the_same_site_in_different_structures_counts_separately(tmp_path):
    d = _run(tmp_path, "m", [
        ["p1", "pdb/aa/1abc.ent.gz", "100", "0.5", "A1,A2", "B10,B11"],
        ["p1", "pdb/bb/2xyz.ent.gz", "100", "0.5", "A1,A2", "B10,B11"],
    ])
    assert residue_hit_counts(d)[("A", 1)] == 2


def test_deduplication_is_per_residue_not_global(tmp_path):
    """One site shared by two residues is one piece of evidence about each."""
    d = _run(tmp_path, "m", [
        ["p1", "pdb/aa/1abc.ent.gz", "100", "0.5", "A1,A2", "B10,B11"],
        ["p2", "pdb/aa/1abc.ent.gz", "100", "0.5", "A2,A3", "B10,B11"],
    ])
    counts = residue_hit_counts(d)

    assert counts[("A", 1)] == counts[("A", 2)] == counts[("A", 3)] == 1


def test_a_residue_is_not_scored_by_how_many_motifs_contain_it(tmp_path):
    """Ten motifs finding one site give 1, not 10.

    This is the failure the deduplication exists to prevent: without it the
    count would track motif overlap, which is a property of the radius grid.
    """
    members = [["p1", "A", 1]]
    for n in range(2, 12):
        members.append([f"p{n}", "A", 1])
    hits = [[f"p{n}", "pdb/aa/1abc.ent.gz", "100", "0.5", "A1", "B10"]
            for n in range(1, 12)]
    d = _run(tmp_path, "m", hits, members=members,
             queried=[f"p{n}" for n in range(1, 12)])

    assert residue_hit_counts(d)[("A", 1)] == 1
    assert residue_motif_counts(d)[("A", 1)] == 11, "coverage is still reported"


# --- what identifies a site -------------------------------------------------

def test_token_order_does_not_split_one_site_in_two(tmp_path):
    d = _run(tmp_path, "m", [
        ["p1", "pdb/aa/1abc.ent.gz", "100", "0.5", "A1,A2", "B10,B11"],
        ["p2", "pdb/aa/1abc.ent.gz", "100", "0.5", "A2,A3", "B11,B10"],
    ])
    assert residue_hit_counts(d)[("A", 2)] == 1


def test_unmatched_placeholders_are_not_part_of_a_site(tmp_path):
    """``B10,B11,_`` and ``B10,B11`` matched the same two residues."""
    assert site_key("t", "B10,B11,_") == site_key("t", "B10,B11")
    assert site_key("t", "B10,_,B11")[1] == ("B10", "B11")


def test_a_row_that_matched_nothing_is_not_a_site(tmp_path):
    d = _run(tmp_path, "m", [
        ["p1", "pdb/aa/1abc.ent.gz", "100", "0.5", "A1,A2", "_,_"],
    ])
    assert residue_hit_counts(d)[("A", 1)] == 0


# --- zero versus absent -----------------------------------------------------

def test_a_searched_motif_that_matched_nothing_scores_zero(tmp_path):
    d = _run(tmp_path, "m", [])
    counts = residue_hit_counts(d)

    assert counts[("A", 1)] == 0, "assessed and found nothing is 0, not absent"
    assert set(counts) == {("A", 1), ("A", 2), ("A", 3)}


def test_a_residue_in_no_emitted_query_is_absent_rather_than_zero(tmp_path):
    """p2 was never queried, so A3 was never assessed at all."""
    d = _run(tmp_path, "m", [
        ["p1", "pdb/aa/1abc.ent.gz", "100", "0.5", "A1,A2", "B10,B11"],
    ], queried=["p1"])
    counts = residue_hit_counts(d)

    assert ("A", 3) not in counts
    assert counts[("A", 2)] == 1


# --- the summary ------------------------------------------------------------

def test_auc_counts_ties_as_half():
    """Ties are the normal case for counts, not an edge case."""
    assert auc([1.0], [1.0]) == 0.5
    assert auc([2.0], [1.0]) == 1.0
    assert auc([], [1.0]) is None


def test_an_unlabelled_run_reports_counts_and_no_auc(tmp_path):
    _run(tmp_path, "design1", [
        ["p1", "pdb/aa/1abc.ent.gz", "100", "0.5", "A1,A2", "B10,B11"],
        ["p2", "pdb/aa/1abc.ent.gz", "100", "0.5", "A2,A3", "B40,B41"],
    ])
    _write(tmp_path / "panel.tsv", ["instance"], [["design1"]])

    proc = subprocess.run(
        [sys.executable, str(SCRIPT), str(tmp_path), "--threshold", "2"],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr

    summary = json.loads(
        (tmp_path / "residue_hit_count_summary.json").read_text())["summary"]
    assert summary["ground_truth_epitopes"] is False
    assert summary["attribution"] == "matched"
    assert "auc_mean" not in summary
    assert summary["residues_scored"] == 3
    # Only A2 reaches two distinct sites.
    assert summary["residues_called"] == 1

    rows = list(csv.DictReader(
        (tmp_path / "residue_hit_counts.tsv").open(), delimiter="\t"))
    by_res = {r["res_id"]: r for r in rows}
    assert by_res["2"]["n_distinct_hits"] == "2"
    assert by_res["2"]["n_motifs"] == "2"
    assert {r["is_epitope"] for r in rows} == {""}, \
        "an unlabelled residue must be blank, never 0 -- 0 means 'not an epitope'"


def test_nothing_is_z_scored(tmp_path):
    """The output must not acquire a z column or a z-like summary key."""
    _run(tmp_path, "design1", [
        ["p1", "pdb/aa/1abc.ent.gz", "100", "0.5", "A1,A2", "B10,B11"],
    ])
    _write(tmp_path / "panel.tsv", ["instance"], [["design1"]])

    proc = subprocess.run(
        [sys.executable, str(SCRIPT), str(tmp_path)],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr

    banned = ("_z", "zscore", "z_score", "sigma", "stdev")
    header = (tmp_path / "residue_hit_counts.tsv").read_text().splitlines()[0]
    assert not any(b in header for b in banned)

    payload = json.loads(
        (tmp_path / "residue_hit_count_summary.json").read_text())
    summary = payload["summary"]
    assert summary["z_scored"] is False
    # `z_scored` is the flag saying there is no z-score; every *other* key
    # must stay clear of the vocabulary. (`residues_with_zero_hits` is a
    # count of zeros, which is why the check is on whole key names.)
    assert not any(
        b in k for k in summary if k != "z_scored" for b in ("zscore", "z_score")
    )
    assert "mean_z" not in summary and "max_z" not in summary

    counts = residue_hit_counts(tmp_path / "design1")
    assert all(isinstance(v, int) for v in counts.values()), \
        "a count must stay an integer count"


def test_the_raw_count_stays_a_whole_count(tmp_path):
    """``residue_hit_counts`` is not divided by coverage.

    The normalised view lives in the separate ``hits_per_motif`` column; the
    numerator must stay an integer count of distinct sites, so that the two
    can be read against each other.
    """
    d = _run(tmp_path, "m", [
        ["p1", "pdb/aa/1abc.ent.gz", "100", "0.5", "A1,A2", "B10,B11"],
        ["p2", "pdb/aa/1abc.ent.gz", "100", "0.5", "A2,A3", "B40,B41"],
    ])
    counts = residue_hit_counts(d)

    assert counts == {("A", 1): 1, ("A", 2): 2, ("A", 3): 1}


def test_only_matched_query_residues_receive_credit(tmp_path):
    d = _run(tmp_path, "m", [
        ["p1", "t", "1", "0", "A1,A2", "B10,_"],
        ["p2", "t", "1", "0", "A2,A3", "_,B11"],
    ])
    assert residue_hit_counts(d) == {("A", 1): 1, ("A", 2): 0, ("A", 3): 1}
    assert residue_hit_counts(d, attribution="motif")[("A", 2)] == 2


def test_unmatched_occurrence_does_not_gain_credit_from_shared_site(tmp_path):
    d = _run(tmp_path, "m", [
        ["p1", "t", "1", "0", "A1,A2", "B10,_"],
        ["p2", "t", "1", "0", "A2,A3", "B10,_"],
    ])
    assert residue_hit_counts(d) == {("A", 1): 1, ("A", 2): 1, ("A", 3): 0}


def test_query_order_not_patch_member_order_controls_attribution(tmp_path):
    d = _run(tmp_path, "m", [["p1", "t", "1", "0", "A2,A1", "_,B10"]])
    assert residue_hit_counts(d)[("A", 1)] == 1
    assert residue_hit_counts(d)[("A", 2)] == 0


@pytest.mark.parametrize("query,target", [
    ("A1,A9", "B10,B11"), ("A1,A1", "B10,B11"),
    ("", "B10"), ("A1,A2", "B10,invalid"),
])
def test_invalid_mapping_is_rejected(tmp_path, query, target):
    """Still fatal: a row that is WRONG, as opposed to merely unalignable.

    A query naming residues that are not the motif's members, a repeated
    token, an empty query or an unparseable target token all mean the row
    cannot be trusted at all. Those stay errors.
    """
    d = _run(tmp_path, "m", [["p1", "t", "1", "0", query, target]])
    with pytest.raises(ValueError):
        residue_hit_counts(d)


@pytest.mark.parametrize("query,target", [
    ("A1,A2", "B10"), ("A1,A2", "B10,_,B11"),
])
def test_an_unalignable_row_is_skipped_and_counted(tmp_path, query, target):
    """Length mismatch means no positional mapping, not a corrupt row.

    Aborting the structure over it discarded 43 of 187 real structures, most
    of them over well under 1% of their rows, so the row is dropped alone.
    The drop is counted, never silent.
    """
    d = _run(tmp_path, "m", [["p1", "t", "1", "0", query, target]])
    report = {}
    counts = residue_hit_counts(d, report=report)

    assert report["rows_ambiguous_skipped"] == 1
    assert report["rows_total"] == 1
    # The residues were assessed and nothing creditable was found, so 0 --
    # not absent, which would mean never searched.
    assert counts[("A", 1)] == 0 and counts[("A", 2)] == 0


def test_good_rows_still_count_when_a_sibling_row_is_skipped(tmp_path):
    """The whole point: one bad row must not cost the structure its data."""
    d = _run(tmp_path, "m", [
        ["p1", "t", "1", "0", "A1,A2", "B10"],          # unalignable
        ["p1", "u", "1", "0", "A1,A2", "B20,B21"],      # fine
    ])
    report = {}
    counts = residue_hit_counts(d, report=report)

    assert report["rows_ambiguous_skipped"] == 1
    assert counts[("A", 1)] == 1, "the good row must still be credited"
    assert counts[("A", 2)] == 1


def test_missing_search_output_is_not_zero(tmp_path):
    d = _run(tmp_path, "m", [])
    (d / "folddisco_hits.tsv").unlink()
    with pytest.raises(ValueError, match="missing search output"):
        residue_hit_counts(d)


def test_cli_reports_skipped_rows_rather_than_failing(tmp_path):
    """An ambiguous structure no longer costs the panel its output.

    What must not happen is the skip going unmentioned: the count reaches
    stderr and the summary, so a reader of the numbers sees what they cost.
    """
    _run(tmp_path, "ok", [["p1", "t", "1", "0", "A1,A2", "B10,B11"]])
    _run(tmp_path, "ambiguous", [["p1", "t", "1", "0", "A1,A2", "B10"]])
    _write(tmp_path / "panel.tsv", ["instance"], [["ok"], ["ambiguous"]])
    proc = subprocess.run([sys.executable, str(SCRIPT), str(tmp_path)],
                          capture_output=True, text=True)

    assert proc.returncode == 0, proc.stderr
    assert "skipped 1" in proc.stderr and "no positional" in proc.stderr
    summary = json.loads(
        (tmp_path / "residue_hit_count_summary.json").read_text())["summary"]
    assert summary["rows_ambiguous_skipped"] == 1
    assert summary["structures_with_skipped_rows"] == 1
    # Both structures are present; the ambiguous one was not dropped.
    rows = list(csv.DictReader(
        (tmp_path / "residue_hit_counts.tsv").open(), delimiter="\t"))
    assert {r["instance"] for r in rows} == {"ok", "ambiguous"}

    proc = subprocess.run([sys.executable, str(SCRIPT), str(tmp_path),
                           "--attribution", "motif"], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    summary = json.loads((tmp_path / "residue_hit_count_summary.json").read_text())
    assert summary["summary"]["attribution"] == "motif"


def test_labelled_partial_match_scores_only_matching_residue(tmp_path):
    _run(tmp_path, "m", [["p1", "t", "1", "0", "A1,A2", "B10,_"]])
    _write(tmp_path / "panel.tsv", ["instance"], [["m"]])
    _write(tmp_path / "epitopes_test.tsv",
           ["instance", "chain", "res_id", "is_epitope"],
           [["m", "A", 1, 1], ["m", "A", 2, 0], ["m", "A", 3, 0]])
    proc = subprocess.run([sys.executable, str(SCRIPT), str(tmp_path)],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    summary = json.loads((tmp_path / "residue_hit_count_summary.json").read_text())["summary"]
    assert summary["attribution"] == "matched"
    assert summary["epitope_mean_hits"] == 1
    assert summary["other_mean_hits"] == 0
    assert summary["auc_mean"] == 1


# --- pipe-joined patch ids --------------------------------------------------

def test_split_patch_ids_handles_the_joined_form():
    from count_scoring import split_patch_ids

    assert split_patch_ids("p1") == ["p1"]
    assert split_patch_ids("p1|p2") == ["p1", "p2"]
    assert split_patch_ids("p1|p2|p3") == ["p1", "p2", "p3"]
    assert split_patch_ids("") == []
    assert split_patch_ids("p1||p2") == ["p1", "p2"], "empty segments dropped"


def _joined_run(root: Path, name: str = "m") -> Path:
    """Two patches sharing one residue set, so the hit row names both.

    This is what `folddisco_search.sh` writes when several patches carry the
    same query string: it cannot tell them apart, so it reports `id1|id2`.
    Under dumb_search that happens whenever two balls enclose the same
    members.
    """
    d = root / name
    _write(d / "patch_members.tsv",
           ["patch_id", "auth_chain_id", "auth_seq_id"],
           [["p1", "A", 1], ["p1", "A", 2],
            ["p2", "A", 1], ["p2", "A", 2]])
    _write(d / "folddisco_queries.index.tsv",
           ["patch_id", "query_string"], [["p1", "A1,A2"], ["p2", "A1,A2"]])
    _write(d / "folddisco_hits.tsv",
           ["patch_id", "tid", "idf", "rmsd", "query_residues",
            "matching_residues"],
           [["p1|p2", "pdb/aa/1abc.ent.gz", "9", "0.3", "A1,A2", "B10,B11"]])
    return d


def test_a_joined_row_is_credited_not_dropped_motif(tmp_path):
    """The regression: 0.43% of real hit rows were silently lost."""
    d = _joined_run(tmp_path)
    counts = residue_hit_counts(d, attribution="motif")

    assert counts[("A", 1)] == 1, "the hit must reach the residue, not vanish"
    assert counts[("A", 2)] == 1


def test_a_joined_row_is_credited_not_dropped_matched(tmp_path):
    """`matched` raised on these rows, aborting 113 of 187 structures."""
    d = _joined_run(tmp_path)
    counts = residue_hit_counts(d, attribution="matched")

    assert counts[("A", 1)] == 1
    assert counts[("A", 2)] == 1


def test_a_joined_row_counts_once_per_residue_not_once_per_patch(tmp_path):
    """Both patches are the SAME residue set, so the site is one site.

    Crediting per patch without the per-residue set union would double this
    residue's count purely because the motif was extracted twice.
    """
    d = _joined_run(tmp_path)
    for mode in ("motif", "matched"):
        counts = residue_hit_counts(d, attribution=mode)
        assert counts[("A", 1)] == 1, mode


def test_an_unqueried_patch_inside_a_joined_row_still_raises(tmp_path):
    """Splitting must not become a way for an unknown patch to slip through."""
    d = _joined_run(tmp_path)
    _write(d / "folddisco_queries.index.tsv",
           ["patch_id", "query_string"], [["p1", "A1,A2"]])

    with pytest.raises(ValueError, match="unknown or unqueried"):
        residue_hit_counts(d, attribution="matched")


def test_joined_patches_must_agree_on_their_member_set(tmp_path):
    """If they disagree they could not have shared a query string."""
    d = _joined_run(tmp_path)
    _write(d / "patch_members.tsv",
           ["patch_id", "auth_chain_id", "auth_seq_id"],
           [["p1", "A", 1], ["p1", "A", 2],
            ["p2", "A", 1], ["p2", "A", 3]])

    with pytest.raises(ValueError, match="uniquely"):
        residue_hit_counts(d, attribution="matched")


# --- normalising by search effort ------------------------------------------

def test_hits_per_motif_divides_the_count_by_the_coverage(tmp_path):
    """The reported rate is exactly n_distinct_hits / n_motifs."""
    _run(tmp_path, "design1", [
        ["p1", "pdb/aa/1abc.ent.gz", "100", "0.5", "A1,A2", "B10,B11"],
        ["p2", "pdb/aa/1abc.ent.gz", "100", "0.5", "A2,A3", "B40,B41"],
    ])
    _write(tmp_path / "panel.tsv", ["instance"], [["design1"]])

    proc = subprocess.run([sys.executable, str(SCRIPT), str(tmp_path)],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr

    rows = {r["res_id"]: r for r in csv.DictReader(
        (tmp_path / "residue_hit_counts.tsv").open(), delimiter="\t")}
    for r in rows.values():
        assert float(r["hits_per_motif"]) == pytest.approx(
            int(r["n_distinct_hits"]) / int(r["n_motifs"]), abs=1e-4)
    # A2 sits in both motifs and reached both sites; A1 in one and reached one.
    assert float(rows["2"]["hits_per_motif"]) == pytest.approx(1.0)
    assert float(rows["1"]["hits_per_motif"]) == pytest.approx(1.0)


def test_normalisation_cancels_pure_search_effort(tmp_path):
    """A residue in eleven motifs that all find one site does not look rich.

    The raw count already deduplicates to 1; the point here is that the rate
    falls as the coverage rises, which is the whole purpose of the column.
    Eleven searches that turned up one site is a rate of 1/11, not of 1.
    """
    members = [[f"p{n}", "A", 1] for n in range(1, 12)]
    hits = [[f"p{n}", "pdb/aa/1abc.ent.gz", "100", "0.5", "A1", "B10"]
            for n in range(1, 12)]
    _run(tmp_path, "design1", hits, members=members,
         queried=[f"p{n}" for n in range(1, 12)])
    _write(tmp_path / "panel.tsv", ["instance"], [["design1"]])

    proc = subprocess.run([sys.executable, str(SCRIPT), str(tmp_path)],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr

    row, = list(csv.DictReader(
        (tmp_path / "residue_hit_counts.tsv").open(), delimiter="\t"))
    assert row["n_distinct_hits"] == "1"
    assert row["n_motifs"] == "11"
    assert float(row["hits_per_motif"]) == pytest.approx(1 / 11, abs=1e-4)


def test_a_rate_is_never_blank_when_the_residue_was_assessed(tmp_path):
    """Every scored residue has a denominator, so every row carries a rate.

    An assessed residue is by definition a member of at least one emitted
    query, so ``n_motifs`` cannot be zero. The blank branch is defensive
    only -- this test is what says so.
    """
    _run(tmp_path, "design1", [
        ["p1", "pdb/aa/1abc.ent.gz", "100", "0.5", "A1,A2", "B10,B11"],
    ])
    _write(tmp_path / "panel.tsv", ["instance"], [["design1"]])

    proc = subprocess.run([sys.executable, str(SCRIPT), str(tmp_path)],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr

    rows = list(csv.DictReader(
        (tmp_path / "residue_hit_counts.tsv").open(), delimiter="\t"))
    assert rows
    assert all(r["hits_per_motif"] != "" for r in rows)
    assert all(int(r["n_motifs"]) >= 1 for r in rows)


def test_the_summary_reports_the_pooled_rate_and_says_what_it_is(tmp_path):
    _run(tmp_path, "design1", [
        ["p1", "pdb/aa/1abc.ent.gz", "100", "0.5", "A1,A2", "B10,B11"],
        ["p2", "pdb/aa/1abc.ent.gz", "100", "0.5", "A2,A3", "B40,B41"],
    ])
    _write(tmp_path / "panel.tsv", ["instance"], [["design1"]])

    proc = subprocess.run([sys.executable, str(SCRIPT), str(tmp_path)],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr

    payload = json.loads(
        (tmp_path / "residue_hit_count_summary.json").read_text())
    summary = payload["summary"]
    assert summary["residues_with_a_rate"] == 3
    assert summary["pooled_mean_hits_per_motif"] == pytest.approx(1.0)
    assert summary["pooled_median_hits_per_motif"] == pytest.approx(1.0)
    # The caveat travels with the number: a reader must not mistake the
    # ratio for the mean hit count of the residue's motifs.
    assert "deduplicated" in summary["rate_note"]
    one, = payload["per_structure"]
    assert one["mean_hits_per_motif"] == pytest.approx(1.0)
    assert one["max_hits_per_motif"] == pytest.approx(1.0)
