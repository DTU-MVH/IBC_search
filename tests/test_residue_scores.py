"""Tests for projecting motif scores onto residues.

Two behaviours are worth pinning, because both fail silently and both would
produce a number that looks like performance:

1. a residue takes the MAXIMUM z of the motifs containing it, so a strong
   motif is not averaged away by weak ones that happen to overlap it;
2. a run with no epitope labels reports call counts and **no AUC**, rather
   than an AUC computed against an empty positive class.
"""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "residue_scores.py"
sys.path.insert(0, str(SCRIPT.parent))

from residue_scores import auc, residue_scores  # noqa: E402


def _write(path: Path, header: list[str], rows: list[list]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def _member_run(root: Path, name: str) -> Path:
    """A member directory whose residue 2 sits in a weak and a strong motif."""
    d = root / name
    _write(d / "patch_members.tsv",
           ["patch_id", "auth_chain_id", "auth_seq_id"],
           [["p1", "A", 1], ["p1", "A", 2],
            ["p2", "A", 2], ["p2", "A", 3]])
    _write(d / "folddisco_queries.index.tsv",
           ["patch_id", "query_string"],
           [["p1", "A1,A2"], ["p2", "A2,A3"]])
    _write(d / "idf_scored.tsv",
           ["query_residues", "idf_z"],
           [["A1,A2", "-0.5000"], ["A2,A3", "+2.4000"]])
    return d


def test_a_residue_takes_the_highest_scoring_motif_it_belongs_to(tmp_path):
    scores = residue_scores(_member_run(tmp_path, "m"))
    assert scores[("A", 2)] == 2.4, "the strong motif must win the overlap"
    assert scores[("A", 1)] == -0.5
    assert scores[("A", 3)] == 2.4


def test_residues_in_no_motif_are_absent_rather_than_zero(tmp_path):
    scores = residue_scores(_member_run(tmp_path, "m"))
    assert ("A", 99) not in scores


def test_auc_counts_ties_as_half():
    assert auc([1.0], [1.0]) == 0.5
    assert auc([2.0], [1.0]) == 1.0
    assert auc([], [1.0]) is None


def test_an_unlabelled_run_reports_calls_and_no_auc(tmp_path):
    """No antibody in the structure means no positive class to score against."""
    _member_run(tmp_path, "design1")
    _write(tmp_path / "panel.tsv", ["instance"], [["design1"]])

    proc = subprocess.run(
        [sys.executable, str(SCRIPT), str(tmp_path), "--threshold", "1.7"],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr

    summary = json.loads(
        (tmp_path / "residue_score_summary.json").read_text())["summary"]
    assert summary["ground_truth_epitopes"] is False
    assert "auc_mean" not in summary
    assert "precision_at_threshold" not in summary
    # A1 scores -0.5; A2 and A3 both take +2.4 from the strong motif.
    assert summary["residues_scored"] == 3
    assert summary["residues_called"] == 2
    assert summary["structures_with_a_call"] == 1

    rows = list(csv.DictReader(
        (tmp_path / "residue_scores.tsv").open(), delimiter="\t"))
    assert {r["is_epitope"] for r in rows} == {""}, \
        "an unlabelled residue must be blank, never 0 -- 0 means 'not an epitope'"
