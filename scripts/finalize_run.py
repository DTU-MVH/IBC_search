"""Produce the standard output layout for a finished FoldDisco run.

Run this after ``run_panel.py`` (or let ``run_panel.py --finalize`` do it) and
every run ends up with the same shape, whatever it was searched against:

::

    <run>/00_figures/            run-level figures
    <run>/01_statistics/         run-level statistics
    <run>/<instance>/00_figures/ per-structure figures, PyMOL session, PDB

Why an orchestrator rather than a list of commands in a README
---------------------------------------------------------------
The pieces already existed -- ``residue_scores.py``, ``count_scoring.py``,
``epitope_statistics.py``, the plotters, ``render_all_structures.py`` -- and
were run by hand in whatever order and with whatever flags the moment
called for. Two runs then came out with different figures, different
statistics and different thresholds, and comparing them meant first working
out what each had actually been given. This file is the answer to "what does
a run always contain", so that question has one answer.

The minimum, produced for every run
-----------------------------------
``00_figures/``
    ``score_bins.png`` (and ``_dark``) -- the **binned dot map**: one bin per
    structure, every residue a point, the **median a black line** with a
    vertical black deviation bar and horizontal caps. On a count run the same
    figure is drawn over hit counts.

    ``count_heatmap.png`` (and ``_dark``) -- structures down, residue
    position across, green to blue by distinct-hit count. Drawn whenever
    count scores exist.

``01_statistics/``
    ``summary.json`` always. When the run has epitope labels, the full
    ``epitope_statistics.py`` battery that ``results/sabdab40_full_r6``
    established as the baseline -- pooled AUC with a bootstrap interval
    resampled over structures, the permutation null, the threshold sweep,
    the per-structure intervals and the coverage artifact -- plus their
    figures in ``00_figures/``.

``<instance>/00_figures/``
    The PyMOL still with its see-through surface, the ``.pse`` session, the
    PyMOL script and a PDB carrying the score in its B-factor column so any
    viewer can reproduce the colouring. CA spheres mark called or epitope
    residues where there are any.

Anything beyond that minimum is opt-in, by flag.

What this file does NOT decide
------------------------------
It runs the existing tools and does not reimplement any of their statistics,
colour choices or interpretation rules. If a number here looks wrong, it is
wrong in the script that produced it, which is where the reasoning for it
lives. Nothing here turns a hit into evidence of antigenicity.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from datetime import datetime
from collections import defaultdict
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[1]
PY = str(REPO / ".venv/bin/python")

sys.path.insert(0, str(REPO / "src"))

from ibc_search.visualization import (  # noqa: E402
    binned_point_plot,
    count_heatmap,
)

FIGURES = "00_figures"
STATISTICS = "01_statistics"

Z_SCORES = "residue_scores.tsv"
COUNTS = "residue_hit_counts.tsv"


def run(cmd: list[str], label: str) -> bool:
    """Run one step. A failure is reported and does not abort the rest.

    One missing optional input should not cost a run its figures, so each
    step is independent and the summary records which ones ran.
    """
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=str(REPO))
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout).strip().splitlines()
        print(f"  {label}: FAILED {tail[-1] if tail else ''}"[:200])
        return False
    print(f"  {label}: ok")
    return True


def instances(run_dir: Path) -> list[str]:
    panel = run_dir / "panel.tsv"
    if not panel.is_file():
        return []
    return [r["instance"] for r in csv.DictReader(panel.open(), delimiter="\t")
            if r.get("instance")]


def has_any(run_dir: Path, insts: list[str], name: str) -> bool:
    return any((run_dir / i / name).is_file() for i in insts)


def read_scores(path: Path, column: str) -> dict[str, list[tuple[int, float]]]:
    """``instance -> [(res_id, value)]``, in file order.

    Residue id is kept alongside the value because the heatmap needs
    sequence order and the dot map does not; reading once serves both.
    """
    out: dict[str, list[tuple[int, float]]] = defaultdict(list)
    if not path.is_file():
        return out
    for r in csv.DictReader(path.open(), delimiter="\t"):
        try:
            out[r["instance"]].append((int(r["res_id"]), float(r[column])))
        except (KeyError, ValueError):
            continue
    return out


def labelled_run(path: Path) -> bool:
    """True when any residue carries a positive epitope label.

    Read from the data rather than from a flag, so a run cannot be told it
    has ground truth that it does not have.
    """
    if not path.is_file():
        return False
    for r in csv.DictReader(path.open(), delimiter="\t"):
        if r.get("is_epitope") == "1":
            return True
    return False


# ---------------------------------------------------------------------------
# The minimum figures
# ---------------------------------------------------------------------------

def draw_binned_dot_map(scores: dict[str, list[tuple[int, float]]],
                        fig_dir: Path, *, stem: str, y_label: str,
                        title: str, subtitle: str, reference: float | None,
                        footnote: str) -> list[str]:
    """One bin per structure, every residue a point, median in black.

    ``summary="median"`` with ``summary_ink="variance"`` is what makes the
    bar the median rather than the top-decile marker the exploratory figures
    use, and draws it in the theme's ink -- black on the light surface. The
    vertical bar spans the IQR and carries horizontal caps at both ends.
    """
    bins = [(inst, [v for _, v in vals])
            for inst, vals in sorted(scores.items()) if vals]
    if not bins:
        return []
    written = []
    for theme in ("light", "dark"):
        out = fig_dir / (f"{stem}.png" if theme == "light" else f"{stem}_dark.png")
        binned_point_plot(
            [("all structures", bins)], out,
            y_label=y_label, title=title, subtitle=subtitle,
            reference=reference,
            reference_label="zero" if reference is not None else "",
            summary="median", summary_ink="variance",
            theme=theme, footnote=footnote,
        )
        written.append(out.name)
    return written


def draw_count_heatmap(counts: dict[str, list[tuple[int, float]]],
                       fig_dir: Path, *, stem: str = "count_heatmap",
                       value_label: str = "distinct hits",
                       title: str = "Distinct matched sites per residue",
                       subtitle: str | None = None,
                       footnote: str | None = None) -> list[str]:
    rows = [(inst, [v for _, v in sorted(vals)])
            for inst, vals in sorted(counts.items()) if vals]
    if not rows:
        return []
    written = []
    for theme in ("light", "dark"):
        out = fig_dir / (f"{stem}.png" if theme == "light"
                         else f"{stem}_dark.png")
        count_heatmap(
            rows, out,
            value_label=value_label,
            title=title,
            subtitle=(subtitle if subtitle is not None else
                      ("one row per structure, residues in sequence order; "
                       "a site is one (reference structure, matched residue "
                       "set) pair, deduplicated per residue")),
            footnote=(footnote if footnote is not None else
                      ("A COUNT, not a z-score: a high value means the "
                       "motifs around that residue are COMMON in the "
                       "searched database. Grey = no residue at that "
                       "position in a shorter chain, which is not the same "
                       "as a count of zero. Nothing here is evidence of "
                       "antigenicity.")),
            theme=theme,
        )
        written.append(out.name)
    return written


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

class RunLock:
    """Refuse to finalise a directory another finalise is already writing.

    Two concurrent runs over one directory interleave their output and
    produce a set that is internally inconsistent in a way nothing in the
    files reveals: this happened once, leaving a ``run_summary.json`` saying
    ``motif`` beside figures built from ``matched`` counts, plus 87 of 187
    per-structure renders coloured from the superseded table. Every file was
    individually well-formed, so only the mtimes gave it away.

    A lock whose pid is gone is stale -- a killed run leaves one behind --
    and is taken over with a message, rather than requiring a manual delete
    that would only teach the reader to delete locks reflexively.
    """

    def __init__(self, run_dir: Path):
        self.path = run_dir / ".finalize.lock"

    def __enter__(self) -> "RunLock":
        if self.path.exists():
            try:
                held = json.loads(self.path.read_text())
            except (json.JSONDecodeError, OSError):
                held = {}
            pid = held.get("pid")
            if pid and _pid_alive(int(pid)):
                raise SystemExit(
                    f"error: {self.path.parent} is already being finalised "
                    f"by pid {pid} (started {held.get('started', 'unknown')})."
                    "\n  Two runs over one directory interleave their output "
                    "and silently mix provenance."
                    "\n  Wait for it, or stop it and re-run."
                )
            print(f"  (taking over a stale lock from pid {pid})")
        self.path.write_text(json.dumps({
            "pid": os.getpid(),
            "started": datetime.now().isoformat(timespec="seconds"),
        }))
        return self

    def __exit__(self, *exc: object) -> None:
        # Only ever remove our OWN lock, so a takeover race cannot delete the
        # lock belonging to the run that won it.
        try:
            held = json.loads(self.path.read_text())
            if held.get("pid") == os.getpid():
                self.path.unlink()
        except (json.JSONDecodeError, OSError):
            pass


def _pid_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def read_scores_with_chain(path: Path, column: str
                           ) -> dict[str, list[tuple[tuple[str, int], float]]]:
    """``instance -> [((chain, res_id), value)]``.

    :func:`read_scores` drops the chain, which is enough for the figures but
    not for joining to ``residues.tsv``: a multi-chain structure reuses
    residue numbers across chains, so ``res_id`` alone is ambiguous there.
    """
    out: dict[str, list[tuple[tuple[str, int], float]]] = defaultdict(list)
    if not path.is_file():
        return out
    for r in csv.DictReader(path.open(), delimiter="\t"):
        try:
            out[r["instance"]].append(
                ((r["chain"], int(r["res_id"])), float(r[column])))
        except (KeyError, ValueError):
            continue
    return out


def load_rsa(results_dir: Path) -> dict[tuple[str, int], float]:
    """``(chain, res_id) -> RSA`` for one structure, missing values omitted.

    A residue whose RSA could not be computed is ABSENT from the mapping
    rather than present as zero. Zero is a real RSA -- a fully buried
    residue -- so substituting it for "not measurable" would pull every mean
    that touches it downward and make buriedness look commoner than it is.
    """
    out: dict[tuple[str, int], float] = {}
    path = results_dir / "residues.tsv"
    if not path.is_file():
        return out
    for r in csv.DictReader(path.open(), delimiter="\t"):
        if r.get("rsa_missing_reason"):
            continue
        try:
            out[(r["auth_chain_id"], int(r["auth_seq_id"]))] = float(r["rsa"])
        except (KeyError, ValueError):
            continue
    return out


def _describe(values: list[float]) -> dict:
    """Mean, sd, quartiles and shape for one pooled sample."""
    import statistics as st

    if not values:
        return {"n": 0}
    v = sorted(values)
    n = len(v)
    mean = st.fmean(v)
    sd = st.pstdev(v) if n > 1 else 0.0

    def q(f: float) -> float:
        return round(v[min(int(f * n), n - 1)], 4)

    out = {
        "n": n,
        "mean": round(mean, 4),
        "sd": round(sd, 4),
        "min": round(v[0], 4),
        "q05": q(0.05), "q25": q(0.25),
        "median": round(st.median(v), 4),
        "q75": q(0.75), "q95": q(0.95),
        "max": round(v[-1], 4),
        "iqr": round(q(0.75) - q(0.25), 4),
    }
    # Reported because the count distribution is visibly right-tailed, and a
    # mean plus sd alone would imply a symmetry it does not have.
    if sd > 0:
        out["skew"] = round(sum(((x - mean) / sd) ** 3 for x in v) / n, 4)
        out["cv"] = round(sd / mean, 4) if mean else None
    return out


def red_residues(counts: list[tuple[Any, float]], top_fraction: float) -> set:
    """The residues a 3-D figure would colour red.

    The same rule :mod:`render_pymol` applies: the top fraction of THIS
    structure's own scored residues, with ties at the boundary INCLUDED.
    Implemented separately there so a figure never depends on the statistics
    module to draw -- which makes their agreement something to test rather
    than assume.
    """
    if not counts:
        return set()
    ranked = sorted((v for _, v in counts), reverse=True)
    cut_index = max(int(round(len(ranked) * top_fraction)) - 1, 0)
    cutoff = ranked[cut_index]
    return {key for key, v in counts if v >= cutoff}


def write_score_statistics(scores: dict[str, list[tuple[int, float]]],
                           stat_dir: Path, *, stem: str, value_name: str,
                           threshold: float) -> dict:
    """Per-structure and pooled distribution statistics, always written.

    The sabdab40 baseline's inferential battery -- pooled AUC, its
    permutation null, the threshold sweep -- needs a positive class and is
    produced only for a labelled run. These are the statistics that are
    defined without one, so that an unlabelled run still has a populated
    ``01_statistics/`` instead of a directory holding one summary file.

    They are DESCRIPTIVE. Nothing here is a performance number, because
    with no labels there is nothing to perform against.
    """
    import statistics as st

    rows = []
    pooled: list[float] = []
    for inst, vals in sorted(scores.items()):
        v = sorted(x for _, x in vals)
        if not v:
            continue
        pooled += v
        rows.append({
            "instance": inst,
            "n_scored": len(v),
            "min": round(v[0], 4),
            "q25": round(v[len(v) // 4], 4),
            "median": round(st.median(v), 4),
            "q75": round(v[(3 * len(v)) // 4], 4),
            "max": round(v[-1], 4),
            "mean": round(st.fmean(v), 4),
            "sd": round(st.pstdev(v), 4) if len(v) > 1 else 0.0,
            "n_at_or_above_threshold": sum(1 for x in v if x >= threshold),
            "n_zero": sum(1 for x in v if x == 0),
        })
    if not rows:
        return {}

    out_tsv = stat_dir / f"{stem}_per_structure.tsv"
    with out_tsv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    pooled.sort()
    return {
        "value": value_name,
        "structures": len(rows),
        "residues_scored": len(pooled),
        "threshold": threshold,
        "residues_at_or_above_threshold":
            sum(1 for x in pooled if x >= threshold),
        "pooled_min": round(pooled[0], 4),
        "pooled_median": round(st.median(pooled), 4),
        "pooled_mean": round(st.fmean(pooled), 4),
        "pooled_sd": round(st.pstdev(pooled), 4) if len(pooled) > 1 else 0.0,
        "pooled_max": round(pooled[-1], 4),
        "residues_scoring_zero": sum(1 for x in pooled if x == 0),
        "per_structure_file": out_tsv.name,
        "note": ("Descriptive only. With no epitope labels there is no "
                 "positive class, so no AUC, precision or recall is defined "
                 "and none is reported here."),
    }


def write_set_statistics(run_dir: Path, insts: list[str],
                         keyed: dict[str, list[tuple[tuple[str, int], float]]],
                         stat_dir: Path, *, value_name: str,
                         top_fraction: float, stem: str = "set_statistics",
                         ) -> dict:
    """Whole-set statistics, including the RSA of the residues drawn red.

    Two things are reported that a per-structure table cannot show:

    * the **pooled** distribution of the score over every assessed residue in
      the set, and separately the distribution of per-structure medians.
      Those answer different questions -- the first describes residues, the
      second describes structures -- and on this data they differ, because
      structures vary about threefold in their typical count.
    * the **solvent accessibility of the residues a figure colours red**,
      against the residues it does not. "Red" is per structure (the top
      fraction of its own scored residues), so the comparison is pooled over
      a consistently-defined set rather than over an absolute cutoff that
      would mean different things in different folds.

    The RSA comparison is DESCRIPTIVE. A difference between the two groups
    says the high-count residues sit at a different mean exposure; it does
    not say exposure causes the count, nor the reverse, and nothing here
    establishes antigenicity.
    """
    import statistics as st

    pooled: list[float] = []
    medians: list[float] = []
    red_rsa: list[float] = []
    other_rsa: list[float] = []
    n_red = n_other = 0
    n_red_no_rsa = n_other_no_rsa = 0

    for inst in insts:
        vals = keyed.get(inst) or []
        if not vals:
            continue
        pooled += [v for _, v in vals]
        medians.append(st.median([v for _, v in vals]))

        rsa = load_rsa(run_dir / inst)
        red = red_residues(vals, top_fraction)
        for key, _ in vals:
            is_red = key in red
            if is_red:
                n_red += 1
            else:
                n_other += 1
            value = rsa.get(key)
            if value is None:
                # Counted, not quietly dropped: an RSA that could not be
                # computed is a different fact from a low one.
                if is_red:
                    n_red_no_rsa += 1
                else:
                    n_other_no_rsa += 1
                continue
            (red_rsa if is_red else other_rsa).append(value)

    out = {
        "value": value_name,
        "structures": len([i for i in insts if keyed.get(i)]),
        "pooled_over_residues": _describe(pooled),
        "per_structure_medians": _describe(medians),
        "red_definition": {
            "rule": ("top fraction of EACH structure's own scored residues, "
                     "ties at the boundary included"),
            "top_fraction": top_fraction,
            "n_red": n_red,
            "n_not_red": n_other,
            "actual_red_fraction": (round(n_red / (n_red + n_other), 4)
                                    if n_red + n_other else None),
        },
        "rsa_of_red_residues": _describe(red_rsa),
        "rsa_of_other_residues": _describe(other_rsa),
        "rsa_missing": {
            "red": n_red_no_rsa,
            "not_red": n_other_no_rsa,
            "note": ("residues whose RSA could not be computed, or that are "
                     "absent from residues.tsv, are EXCLUDED from the RSA "
                     "means rather than entered as zero: zero is a real RSA "
                     "and would bias the means downward"),
        },
        "notes": [
            ("pooled_over_residues describes RESIDUES; per_structure_medians "
             "describes STRUCTURES. They answer different questions and on "
             "this data they differ, because structures vary severalfold in "
             "their typical count."),
            ("The RSA comparison is descriptive. It does not establish that "
             "exposure drives the count, or the reverse."),
            ("The actual red fraction can exceed top_fraction: ties at the "
             "boundary are included, and integer counts tie often."),
        ],
    }
    (stat_dir / f"{stem}.json").write_text(json.dumps(out, indent=2))

    # A flat companion to the JSON. The JSON nests, which is right for the
    # caveats but wrong for reading four distributions side by side, and
    # side by side is the whole question: whether the red residues differ
    # from the rest. One row per distribution, one column per statistic.
    groups = [
        ("score_pooled_over_residues", out["pooled_over_residues"]),
        ("score_per_structure_medians", out["per_structure_medians"]),
        ("rsa_red_residues", out["rsa_of_red_residues"]),
        ("rsa_other_residues", out["rsa_of_other_residues"]),
    ]
    columns = ["group", "n", "mean", "sd", "min", "q05", "q25", "median",
               "q75", "q95", "max", "iqr", "skew", "cv"]
    tsv = stat_dir / f"{stem}.tsv"
    with tsv.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(columns)
        for name, d in groups:
            w.writerow([name] + [d.get(c, "") for c in columns[1:]])
    out["set_statistics_file"] = tsv.name
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", type=Path, help="a directory holding panel.tsv")
    ap.add_argument("--threshold", type=float, default=1.6,
                    help="z at or above which a residue is called [1.6, the "
                         "sabdab40 baseline's headline threshold]")
    ap.add_argument("--top-fraction", type=float, default=0.20,
                    help="fraction of EACH structure's own scored residues "
                         "marked in the 3-D figures. Per-structure by "
                         "design: a hit count has no absolute scale that "
                         "transfers between structures [0.20]")
    ap.add_argument("--count-threshold", type=int, default=1,
                    help="distinct hits at or above which a residue is "
                         "called [1]")
    ap.add_argument("--count-attribution", choices=("matched", "motif"), default="matched",
                    help="count only matched query residues [matched]; motif selects legacy attribution")
    ap.add_argument("-j", "--jobs", type=int, default=6)
    ap.add_argument("--no-pymol", action="store_true",
                    help="skip the per-structure renders (they need PyMOL)")
    ap.add_argument("--bootstrap", type=int, default=2000)
    ap.add_argument("--permutations", type=int, default=2000)
    args = ap.parse_args(argv)

    run_dir = args.run_dir
    insts = instances(run_dir)
    if not insts:
        print(f"error: no panel.tsv with instances at {run_dir}",
              file=sys.stderr)
        return 2

    fig_dir = run_dir / FIGURES
    stat_dir = run_dir / STATISTICS
    fig_dir.mkdir(parents=True, exist_ok=True)
    stat_dir.mkdir(parents=True, exist_ok=True)

    with RunLock(run_dir):
        return _finalize(args, run_dir, insts, fig_dir, stat_dir)


def _finalize(args, run_dir: Path, insts: list[str], fig_dir: Path,
              stat_dir: Path) -> int:
    """The body, held under the directory lock for its whole duration."""

    produced: dict[str, list[str]] = {"figures": [], "statistics": [],
                                      "per_structure": []}
    steps: dict[str, bool] = {}

    # --- scoring tables, built only where their inputs exist -------------
    print("scoring")
    if has_any(run_dir, insts, "idf_scored.tsv"):
        steps["residue_scores"] = run(
            [PY, str(REPO / "scripts/residue_scores.py"), str(run_dir),
             "--threshold", str(args.threshold)], "residue_scores")
    else:
        print("  residue_scores: skipped (no idf_scored.tsv; the run was not "
              "z-scored, which is expected against a database with no "
              "fitted IDF parameters)")
    if has_any(run_dir, insts, "folddisco_hits.tsv"):
        steps["count_scoring"] = run(
            [PY, str(REPO / "scripts/count_scoring.py"), str(run_dir),
             "--threshold", str(args.count_threshold),
             "--attribution", args.count_attribution], "count_scoring")
        if not steps["count_scoring"]:
            print("error: count scoring failed; refusing to visualize stale or incomplete counts",
                  file=sys.stderr)
            return 2
    else:
        print("  count_scoring: skipped (no folddisco_hits.tsv)")

    z_path, count_path = run_dir / Z_SCORES, run_dir / COUNTS
    z_scores = read_scores(z_path, "max_motif_z")
    counts = read_scores(count_path, "n_distinct_hits")
    # Effort-normalised counts: distinct sites per motif the residue was
    # searched in. Residues whose denominator is unknown carry a blank and
    # are absent here rather than present as zero.
    rates = read_scores(count_path, "hits_per_motif")
    labelled = labelled_run(z_path) or labelled_run(count_path)

    # --- the minimum figures ---------------------------------------------
    print("figures")
    if z_scores:
        produced["figures"] += draw_binned_dot_map(
            z_scores, fig_dir, stem="score_bins",
            y_label="max motif z-score",
            title="Motif z-score per residue, by structure",
            subtitle="every assessed residue is a point; one bin per structure",
            reference=0.0,
            footnote=("A z-score is a position in a fitted null, not a "
                      "probability. Residues in no emitted motif are absent, "
                      "not zero."))
        print(f"  score_bins: ok ({len(z_scores)} structures)")
    if counts:
        produced["figures"] += draw_binned_dot_map(
            counts, fig_dir, stem="count_bins",
            y_label="distinct matched sites",
            title="Distinct hits per residue, by structure",
            subtitle="every assessed residue is a point; one bin per structure",
            reference=None,
            footnote=("A COUNT, deduplicated per residue, NOT z-scored. "
                      "Zero means assessed and unmatched; a residue in no "
                      "emitted query is absent instead."))
        produced["figures"] += draw_count_heatmap(counts, fig_dir)
        print(f"  count_bins + count_heatmap: ok ({len(counts)} structures)")
    if rates:
        produced["figures"] += draw_binned_dot_map(
            rates, fig_dir, stem="rate_bins",
            y_label="distinct sites per motif searched",
            title="Distinct hits per motif searched, by structure",
            subtitle="every assessed residue is a point; one bin per structure",
            reference=None,
            footnote=("The count DIVIDED by the number of emitted queries "
                      "containing the residue, so a residue is not scored "
                      "higher merely for sitting in more balls. Not a "
                      "per-motif mean: the numerator is deduplicated across "
                      "those motifs. Bounded above by the search --top cap."))
        produced["figures"] += draw_count_heatmap(
            rates, fig_dir, stem="rate_heatmap",
            value_label="sites per motif",
            title="Distinct matched sites per motif searched",
            subtitle=("one row per structure, residues in sequence order; "
                      "the count divided by the number of emitted queries "
                      "containing the residue"),
            footnote=("Normalised for search effort, so a mid-chain residue "
                      "is not scored higher merely for falling inside more "
                      "balls. Not a per-motif mean: the numerator is "
                      "deduplicated across those motifs. Bounded above by "
                      "the search --top cap. Grey = no residue at that "
                      "position, or a residue with no known denominator -- "
                      "neither is a rate of zero."))
        print(f"  rate_bins + rate_heatmap: ok ({len(rates)} structures)")
    if not z_scores and not counts:
        print("  nothing scored, so no run-level figures", file=sys.stderr)

    # --- statistics -------------------------------------------------------
    print("statistics")
    if labelled:
        # The sabdab40 baseline battery. Only defined with a positive class,
        # which is why it is gated on the labels being present in the data.
        steps["epitope_statistics"] = run(
            [PY, str(REPO / "scripts/epitope_statistics.py"), str(run_dir),
             "--out-dir", str(stat_dir), "--threshold", str(args.threshold),
             "--bootstrap", str(args.bootstrap),
             "--permutations", str(args.permutations)], "epitope_statistics")
        steps["plot_epitope_statistics"] = run(
            [PY, str(REPO / "scripts/plot_epitope_statistics.py"), str(run_dir),
             "--stats-dir", str(stat_dir), "--out-dir", str(fig_dir)],
            "plot_epitope_statistics")
        steps["plot_epitope_scores"] = run(
            [PY, str(REPO / "scripts/plot_epitope_scores.py"), str(run_dir),
             "--out-dir", str(fig_dir), "--threshold", str(args.threshold)],
            "plot_epitope_scores")
    else:
        print("  epitope battery: skipped -- this run has NO epitope labels, "
              "so there is no positive class and an AUC, a precision or a "
              "recall computed here would measure nothing")

    descriptive: dict = {}
    if z_scores:
        descriptive["motif_z"] = write_score_statistics(
            z_scores, stat_dir, stem="motif_z",
            value_name="max motif z-score", threshold=args.threshold)
    if counts:
        descriptive["distinct_hits"] = write_score_statistics(
            counts, stat_dir, stem="distinct_hits",
            value_name="distinct matched sites",
            threshold=float(args.count_threshold))
    if rates:
        descriptive["hits_per_motif"] = write_score_statistics(
            rates, stat_dir, stem="hits_per_motif",
            value_name="distinct matched sites per motif searched",
            threshold=1.0)
    if descriptive:
        (stat_dir / "distributions.json").write_text(
            json.dumps(descriptive, indent=2))
        print(f"  distributions: ok ({', '.join(descriptive)})")

    # Whole-set statistics, produced for EVERY run: the pooled distribution,
    # the spread between structures, and the solvent accessibility of the
    # residues a figure colours red against those it does not.
    if counts or z_scores:
        path = count_path if counts else z_path
        column = "n_distinct_hits" if counts else "max_motif_z"
        keyed = read_scores_with_chain(path, column)
        setstats = write_set_statistics(
            run_dir, insts, keyed, stat_dir,
            value_name=("distinct matched sites" if counts
                        else "max motif z-score"),
            top_fraction=args.top_fraction)
        red = setstats["rsa_of_red_residues"]
        other = setstats["rsa_of_other_residues"]
        print(f"  set_statistics: ok (pooled n="
              f"{setstats['pooled_over_residues'].get('n', 0):,}; "
              f"mean RSA red {red.get('mean')} vs other {other.get('mean')})")

    # The same whole-set battery on the effort-normalised score, written
    # separately rather than in place of the raw one. "Red" is re-derived
    # from the ratio here, so these red residues are NOT the ones the 3-D
    # figures colour -- those follow the raw count. Two files, so neither
    # table can be mistaken for the other.
    if rates:
        keyed_rates = read_scores_with_chain(count_path, "hits_per_motif")
        ratestats = write_set_statistics(
            run_dir, insts, keyed_rates, stat_dir,
            value_name="distinct matched sites per motif searched",
            top_fraction=args.top_fraction,
            stem="set_statistics_hits_per_motif")
        ratestats["red_definition"]["warning"] = (
            "ranked on hits_per_motif, so this red set differs from the one "
            "the 3-D figures colour, which ranks on n_distinct_hits")
        (stat_dir / "set_statistics_hits_per_motif.json").write_text(
            json.dumps(ratestats, indent=2))
        rred = ratestats["rsa_of_red_residues"]
        rother = ratestats["rsa_of_other_residues"]
        print(f"  set_statistics_hits_per_motif: ok (pooled n="
              f"{ratestats['pooled_over_residues'].get('n', 0):,}; "
              f"mean RSA red {rred.get('mean')} vs other {rother.get('mean')})")

    # --- per-structure figures -------------------------------------------
    print("per-structure figures")
    if args.no_pymol:
        print("  skipped (--no-pymol)")
    elif z_scores or counts:
        scores_for_render = z_path if z_scores else count_path
        column = "max_motif_z" if z_scores else "n_distinct_hits"

        # The 3-D figures keep the red/grey scheme whichever scoring drove
        # them: one convention in 3-D, so a reader never has to work out
        # which scale a structure was coloured on. The green-to-blue ramp is
        # the heatmap's alone.
        steps["render_all_structures"] = run(
            [PY, str(REPO / "scripts/render_all_structures.py"), str(run_dir),
             "--scores", str(scores_for_render),
             "--score-column", column,
             "--top-fraction", str(args.top_fraction),
             "--threshold", str(args.threshold if z_scores
                                else args.count_threshold),
             "--figures-dir", FIGURES, "-j", str(args.jobs)],
            "render_all_structures")
        produced["per_structure"] = sorted(
            p.name for i in insts
            for p in (run_dir / i / FIGURES).glob("*") if p.is_file()
        )[:8]
        print(f"  coloured by {column}")
    else:
        print("  skipped (nothing scored)")

    produced["figures"] = sorted({p.name for p in fig_dir.glob("*.png")})

    summary = {
        "run": str(run_dir),
        "structures": len(insts),
        "scored_z": len(z_scores),
        "scored_counts": len(counts),
        "scored_rates": len(rates),
        "has_epitope_labels": labelled,
        "threshold_z": args.threshold,
        "threshold_counts": args.count_threshold,
        "figure_top_fraction": args.top_fraction,
        "count_attribution": args.count_attribution,
        "steps": steps,
        "produced": produced,
        "layout": {
            "00_figures": "run-level figures",
            "01_statistics": "run-level statistics",
            "<instance>/00_figures": ("per-structure PyMOL still, .pse "
                                      "session, script and B-factor PDB"),
        },
        "notes": [
            ("The binned dot map draws the MEDIAN as a black line with a "
             "vertical black IQR bar and horizontal caps. On the dark theme "
             "that ink is white: black would be invisible."),
            ("Hit counts are NOT z-scored. A high count means the motif is "
             "COMMON in the searched database, close to the opposite of what "
             "a high IDF z-score means."),
            ("3-D figures use the red/grey scheme for both scorings; the "
             "green-to-blue ramp is used only in the 2-D count heatmap."),
            ("Residues are marked red in the 3-D figures when they fall in "
             "the top fraction of THAT structure's own scored residues, not "
             "against an absolute cutoff: a hit count depends on how much of "
             "the database resembles that particular fold, so a count that "
             "is high for one structure is unremarkable in another. Ties at "
             "the boundary are included, so the marked set can exceed the "
             "nominal fraction."),
        ],
        "interpretation_notice": (
            "Nothing in this run establishes antigenicity, antibody binding "
            "or anti-drug-antibody reactivity. A hit is a geometric motif "
            "match against a structure database."),
    }
    if not labelled:
        summary["no_ground_truth_note"] = (
            "No epitope labels: no AUC, precision or recall is defined, and "
            "none was computed. Every call is a prediction with nothing to "
            "score it against.")
    summary["produced"]["statistics"] = sorted(
        {p.name for p in stat_dir.glob("*")} | {"run_summary.json"})
    (stat_dir / "run_summary.json").write_text(json.dumps(summary, indent=2))
    produced["statistics"] = summary["produced"]["statistics"]

    print(f"\n{run_dir}")
    print(f"  {FIGURES}/     {len(produced['figures'])} figures")
    print(f"  {STATISTICS}/  {len(produced['statistics'])} files")
    failed = [k for k, v in steps.items() if not v]
    if failed:
        print(f"  FAILED steps: {', '.join(failed)}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
