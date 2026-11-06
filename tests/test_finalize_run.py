"""The standard run layout: what every FoldDisco run must contain.

These pin the contract rather than the look of any one figure. The failure
they guard against is drift -- two runs finalised at different times coming
out with different figures, different statistics and different thresholds,
which is what made runs incomparable before there was an orchestrator.
"""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "finalize_run.py"
sys.path.insert(0, str(SCRIPT.parent))

from finalize_run import (  # noqa: E402
    labelled_run,
    read_scores,
    write_score_statistics,
)


def _write(path: Path, header: list[str], rows: list[list]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def test_labelled_run_is_read_from_the_data_not_a_flag(tmp_path):
    unlabelled = tmp_path / "a.tsv"
    _write(unlabelled, ["instance", "is_epitope"], [["x", "0"], ["x", ""]])
    assert not labelled_run(unlabelled)

    labelled = tmp_path / "b.tsv"
    _write(labelled, ["instance", "is_epitope"], [["x", "0"], ["x", "1"]])
    assert labelled_run(labelled)


def test_a_missing_scores_file_is_not_a_labelled_run(tmp_path):
    assert not labelled_run(tmp_path / "nope.tsv")


def test_read_scores_skips_rows_it_cannot_parse(tmp_path):
    path = tmp_path / "s.tsv"
    _write(path, ["instance", "res_id", "n_distinct_hits"],
           [["a", "1", "5"], ["a", "2", ""], ["a", "x", "7"], ["b", "3", "9"]])
    got = read_scores(path, "n_distinct_hits")

    assert got["a"] == [(1, 5.0)]
    assert got["b"] == [(3, 9.0)]


def test_descriptive_statistics_are_written_without_labels(tmp_path):
    """An unlabelled run must still populate 01_statistics/."""
    scores = {"s1": [(i, float(i)) for i in range(10)],
              "s2": [(i, 4.0) for i in range(5)]}
    summary = write_score_statistics(
        scores, tmp_path, stem="distinct_hits",
        value_name="distinct matched sites", threshold=1.0)

    assert summary["structures"] == 2
    assert summary["residues_scored"] == 15
    assert summary["pooled_max"] == 9.0
    assert summary["residues_scoring_zero"] == 1
    assert (tmp_path / "distinct_hits_per_structure.tsv").is_file()


def test_descriptive_statistics_report_no_performance_number(tmp_path):
    """Without a positive class, an AUC would be measuring nothing."""
    summary = write_score_statistics(
        {"s": [(1, 2.0), (2, 3.0)]}, tmp_path, stem="x",
        value_name="v", threshold=1.0)

    assert not any(k in summary for k in ("auc", "precision", "recall", "f1"))
    assert "no positive class" in summary["note"]


def _count_run(root: Path) -> Path:
    """A minimal finished run: one structure, hits, no epitope labels."""
    _write(root / "panel.tsv", ["instance", "pdb_id"], [["s1", "1ABC"]])
    d = root / "s1"
    _write(d / "patch_members.tsv",
           ["patch_id", "auth_chain_id", "auth_seq_id"],
           [["p1", "A", 1], ["p1", "A", 2], ["p2", "A", 2], ["p2", "A", 3]])
    _write(d / "folddisco_queries.index.tsv",
           ["patch_id", "query_string"], [["p1", "A1,A2"], ["p2", "A2,A3"]])
    _write(d / "folddisco_hits.tsv",
           ["patch_id", "tid", "idf", "rmsd", "query_residues",
            "matching_residues"],
           [["p1", "t/1abc.gz", "9", "0.3", "A1,A2", "B10,B11"],
            ["p2", "t/2xyz.gz", "9", "0.3", "A2,A3", "B40,B41"]])
    return root


def test_a_run_gets_the_standard_layout(tmp_path):
    _count_run(tmp_path)
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), str(tmp_path), "--no-pymol"],
        capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr

    assert (tmp_path / "00_figures").is_dir()
    assert (tmp_path / "01_statistics").is_dir()
    # The binned dot map and the count heatmap are the stated minimum, in
    # both themes.
    for name in ("count_bins.png", "count_bins_dark.png",
                 "count_heatmap.png", "count_heatmap_dark.png"):
        assert (tmp_path / "00_figures" / name).is_file(), name
    assert (tmp_path / "01_statistics" / "run_summary.json").is_file()
    assert (tmp_path / "01_statistics" / "distributions.json").is_file()


def test_an_unlabelled_run_reports_no_ground_truth_and_no_auc(tmp_path):
    _count_run(tmp_path)
    subprocess.run([sys.executable, str(SCRIPT), str(tmp_path), "--no-pymol"],
                   capture_output=True, text=True, check=True)

    summary = json.loads(
        (tmp_path / "01_statistics" / "run_summary.json").read_text())
    assert summary["has_epitope_labels"] is False
    assert "no_ground_truth_note" in summary
    assert not (tmp_path / "01_statistics" / "summary.json").is_file(), \
        "the epitope battery must not run without a positive class"


def test_the_summary_lists_what_it_actually_produced(tmp_path):
    _count_run(tmp_path)
    subprocess.run([sys.executable, str(SCRIPT), str(tmp_path), "--no-pymol"],
                   capture_output=True, text=True, check=True)

    summary = json.loads(
        (tmp_path / "01_statistics" / "run_summary.json").read_text())
    on_disk = {p.name for p in (tmp_path / "00_figures").glob("*.png")}
    assert set(summary["produced"]["figures"]) == on_disk
    assert "run_summary.json" in summary["produced"]["statistics"]


def test_failed_count_scoring_stops_finalization_before_stale_figures(tmp_path):
    """Figures must never be drawn over counts that did not recompute.

    The trigger is a row that is genuinely CORRUPT -- an unparseable target
    token -- not merely one without a positional mapping. An unalignable row
    is skipped and counted now, so it no longer stops anything; see
    test_an_ambiguous_row_no_longer_stops_finalization.
    """
    _count_run(tmp_path)
    path = tmp_path / "s1" / "folddisco_hits.tsv"
    path.write_text(path.read_text().replace("B10,B11", "B10,!!bad"))
    _write(tmp_path / "residue_hit_counts.tsv",
           ["instance", "res_id", "n_distinct_hits"], [["s1", 1, 999]])
    proc = subprocess.run([sys.executable, str(SCRIPT), str(tmp_path), "--no-pymol"],
                          capture_output=True, text=True)
    assert proc.returncode == 2
    assert "refusing to visualize stale" in proc.stderr
    assert not list((tmp_path / "00_figures").glob("*.png"))
    proc = subprocess.run([sys.executable, str(SCRIPT), str(tmp_path), "--no-pymol",
                           "--count-attribution", "motif"], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    assert (tmp_path / "00_figures" / "count_bins.png").is_file()


def test_an_ambiguous_row_no_longer_stops_finalization(tmp_path):
    """A row with no positional mapping costs that row, not the run."""
    _count_run(tmp_path)
    path = tmp_path / "s1" / "folddisco_hits.tsv"
    path.write_text(path.read_text().replace("B10,B11", "B10"))

    proc = subprocess.run([sys.executable, str(SCRIPT), str(tmp_path), "--no-pymol"],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    assert (tmp_path / "00_figures" / "count_bins.png").is_file()

    summary = json.loads(
        (tmp_path / "residue_hit_count_summary.json").read_text())["summary"]
    assert summary["rows_ambiguous_skipped"] >= 1, "the cost must be recorded"


# --- concurrency ------------------------------------------------------------

def test_a_second_finalize_on_the_same_directory_refuses(tmp_path):
    """Two concurrent runs interleave and mix provenance invisibly.

    This really happened: a `motif` run and a `matched` run overlapped on
    one directory and left run_summary.json saying one thing while the
    figures came from the other, with 87 of 187 renders built from the
    superseded counts. Every file was well-formed; only mtimes showed it.
    """
    import json as _json
    import os as _os

    from finalize_run import RunLock

    _count_run(tmp_path)
    # A lock held by a process that really is alive: this interpreter.
    (tmp_path / ".finalize.lock").write_text(
        _json.dumps({"pid": _os.getpid(), "started": "now"}))

    proc = subprocess.run(
        [sys.executable, str(SCRIPT), str(tmp_path), "--no-pymol"],
        capture_output=True, text=True)

    assert proc.returncode != 0
    assert "already being finalised" in (proc.stdout + proc.stderr)
    assert not list((tmp_path / "00_figures").glob("*.png")), \
        "it must refuse before writing anything"


def test_a_stale_lock_is_taken_over(tmp_path):
    """A killed run leaves a lock behind; that must not need manual cleanup."""
    import json as _json

    _count_run(tmp_path)
    # pid 2**22 is above the macOS/Linux maximum, so it cannot be running.
    (tmp_path / ".finalize.lock").write_text(
        _json.dumps({"pid": 2 ** 22, "started": "long ago"}))

    proc = subprocess.run(
        [sys.executable, str(SCRIPT), str(tmp_path), "--no-pymol"],
        capture_output=True, text=True)

    assert proc.returncode == 0, proc.stderr
    assert "stale lock" in proc.stdout
    assert (tmp_path / "00_figures" / "count_bins.png").is_file()


def test_the_lock_is_released_on_success(tmp_path):
    _count_run(tmp_path)
    subprocess.run([sys.executable, str(SCRIPT), str(tmp_path), "--no-pymol"],
                   capture_output=True, text=True, check=True)

    assert not (tmp_path / ".finalize.lock").exists()


# --- whole-set statistics ---------------------------------------------------

def test_red_residues_matches_the_figure_rule(tmp_path):
    """The statistics and the 3-D figure must agree on what "red" is.

    render_pymol computes the cutoff independently so a figure never depends
    on the statistics module to draw. That makes agreement a thing to test,
    not a thing to assume.
    """
    from finalize_run import red_residues

    counts = [((("A"), i), float(v)) for i, v in enumerate([1, 2, 3, 4, 5,
                                                            6, 7, 8, 9, 10])]
    red = red_residues(counts, 0.20)
    assert len(red) == 2, "top 20% of ten residues"
    assert {k[1] for k in red} == {8, 9}, "the two highest"


def test_red_includes_ties_so_it_can_exceed_the_fraction(tmp_path):
    """Integer counts tie often; excluding ties would make the marked set
    depend on dict order rather than on the data."""
    from finalize_run import red_residues

    counts = [(("A", i), 5.0) for i in range(10)]
    assert len(red_residues(counts, 0.20)) == 10, "all tied at the cutoff"


def test_set_statistics_reports_rsa_for_red_and_other(tmp_path):
    from finalize_run import write_set_statistics

    _count_run(tmp_path)
    _write(tmp_path / "s1" / "residues.tsv",
           ["auth_chain_id", "auth_seq_id", "rsa", "rsa_missing_reason"],
           [["A", 1, "0.80", ""], ["A", 2, "0.10", ""], ["A", 3, "0.50", ""]])
    keyed = {"s1": [(("A", 1), 10.0), (("A", 2), 1.0), (("A", 3), 2.0)]}

    out = write_set_statistics(tmp_path, ["s1"], keyed,
                               tmp_path, value_name="v", top_fraction=0.34)

    assert out["red_definition"]["n_red"] == 1
    assert out["rsa_of_red_residues"]["mean"] == 0.80, "the top residue"
    assert out["rsa_of_other_residues"]["n"] == 2
    assert (tmp_path / "set_statistics.json").is_file()


def test_missing_rsa_is_excluded_not_counted_as_zero(tmp_path):
    """Zero is a real RSA -- a buried residue -- so it must not stand in
    for 'could not be measured'."""
    from finalize_run import write_set_statistics

    _count_run(tmp_path)
    _write(tmp_path / "s1" / "residues.tsv",
           ["auth_chain_id", "auth_seq_id", "rsa", "rsa_missing_reason"],
           [["A", 1, "0.80", ""], ["A", 2, "", "no reference maximum"]])
    keyed = {"s1": [(("A", 1), 10.0), (("A", 2), 1.0)]}

    out = write_set_statistics(tmp_path, ["s1"], keyed,
                               tmp_path, value_name="v", top_fraction=0.5)

    assert out["rsa_of_other_residues"]["n"] == 0, "the missing one is absent"
    assert out["rsa_missing"]["not_red"] == 1, "and it is counted as missing"
    assert out["rsa_of_red_residues"]["mean"] == 0.80


def test_set_statistics_separates_residues_from_structures(tmp_path):
    """Pooled-over-residues and per-structure-medians answer different
    questions and must not be conflated."""
    from finalize_run import write_set_statistics

    _count_run(tmp_path)
    (tmp_path / "s2").mkdir(exist_ok=True)
    keyed = {"s1": [(("A", i), 1.0) for i in range(100)],
             "s2": [(("A", 1), 500.0)]}

    out = write_set_statistics(tmp_path, ["s1", "s2"], keyed,
                               tmp_path, value_name="v", top_fraction=0.2)

    # 101 residues pooled, but only 2 structures.
    assert out["pooled_over_residues"]["n"] == 101
    assert out["per_structure_medians"]["n"] == 2
    assert out["pooled_over_residues"]["mean"] != \
        out["per_structure_medians"]["mean"]


def test_set_statistics_also_writes_a_flat_tsv(tmp_path):
    """The JSON nests for the caveats; the TSV is for reading the four
    distributions side by side, which is the actual comparison."""
    import csv as _csv

    from finalize_run import write_set_statistics

    _count_run(tmp_path)
    _write(tmp_path / "s1" / "residues.tsv",
           ["auth_chain_id", "auth_seq_id", "rsa", "rsa_missing_reason"],
           [["A", 1, "0.80", ""], ["A", 2, "0.10", ""], ["A", 3, "0.50", ""]])
    keyed = {"s1": [(("A", 1), 10.0), (("A", 2), 1.0), (("A", 3), 2.0)]}

    write_set_statistics(tmp_path, ["s1"], keyed, tmp_path,
                         value_name="v", top_fraction=0.34)

    rows = list(_csv.DictReader((tmp_path / "set_statistics.tsv").open(),
                                delimiter="\t"))
    assert [r["group"] for r in rows] == [
        "score_pooled_over_residues", "score_per_structure_medians",
        "rsa_red_residues", "rsa_other_residues",
    ]
    assert float(rows[2]["mean"]) == 0.80, "red RSA row"
    assert int(rows[3]["n"]) == 2, "other RSA row"


# --- the effort-normalised score -------------------------------------------

def test_the_normalised_score_gets_its_own_figures_and_statistics(tmp_path):
    """Normalising adds files; it never overwrites the raw-count ones.

    Both views must survive side by side: the 3-D figures colour on the raw
    count, so a statistics file that silently became the normalised one
    would describe a different red set from the one on screen.
    """
    _count_run(tmp_path)
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), str(tmp_path), "--no-pymol"],
        capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr

    figs = tmp_path / "00_figures"
    stats = tmp_path / "01_statistics"
    for name in ("count_bins.png", "count_heatmap.png",
                 "rate_bins.png", "rate_bins_dark.png",
                 "rate_heatmap.png", "rate_heatmap_dark.png"):
        assert (figs / name).is_file(), name
    assert (stats / "set_statistics.json").is_file()
    assert (stats / "set_statistics_hits_per_motif.json").is_file()
    assert (stats / "distinct_hits_per_structure.tsv").is_file()
    assert (stats / "hits_per_motif_per_structure.tsv").is_file()

    descriptive = json.loads((stats / "distributions.json").read_text())
    assert "distinct_hits" in descriptive and "hits_per_motif" in descriptive


def test_the_normalised_red_set_warns_that_it_is_not_the_figures(tmp_path):
    _count_run(tmp_path)
    subprocess.run([sys.executable, str(SCRIPT), str(tmp_path), "--no-pymol"],
                   capture_output=True, text=True, check=True)

    out = json.loads((tmp_path / "01_statistics"
                      / "set_statistics_hits_per_motif.json").read_text())
    assert "n_distinct_hits" in out["red_definition"]["warning"]


def test_a_blank_rate_is_absent_rather_than_zero(tmp_path):
    """read_scores must drop an unparseable value, not read it as 0.0.

    A residue with no known denominator has no rate. Entering it as zero
    would put it at the bottom of every distribution and of the red ranking.
    """
    from finalize_run import read_scores

    path = tmp_path / "counts.tsv"
    _write(path, ["instance", "res_id", "hits_per_motif"],
           [["s1", 1, "2.5"], ["s1", 2, ""], ["s1", 3, "0.5"]])

    assert read_scores(path, "hits_per_motif") == {"s1": [(1, 2.5), (3, 0.5)]}
