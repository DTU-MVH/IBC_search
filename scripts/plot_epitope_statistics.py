"""Figures for the statistics in ``01_statistics/``, written to ``00_figures/``.

``plot_epitope_scores.py`` draws the score distributions. This draws the
*inferential* views that go with ``epitope_statistics.py``: the pooled ROC
with its bootstrap band, the threshold sweep, the per-structure intervals,
and the coverage artifact.

Design choices carried over from ``ibc_search.visualization``
-------------------------------------------------------------
The themes are imported rather than redefined, so every figure in the
repository shares one surface, one ink and one validated hue. Light and dark
are each selected against their own surface and are not inversions.

The sweep is **faceted, one metric per panel**, rather than four lines in
four colours on one axis. Four categorical hues compared all-pairs do not
clear the palette's separation floor, and the documented remedy is to facet:
each panel carries a single series, so identity comes from the panel title
and the figure survives greyscale and colour-vision deficiency. It also
avoids putting precision, recall, F1 and MCC -- which do not share a range,
MCC being signed -- on one scale.

Error bars, not bare points
---------------------------
The per-structure panel draws each AUC with its bootstrap interval. The
point of that figure is that the intervals are wide enough to overlap almost
everything, so a reader does not come away ranking the structures. Drawing
the points alone would invite exactly the reading the interval forbids.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import textwrap
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ibc_search.visualization import THEMES  # noqa: E402


def _style(ax, t, *, x_label="", y_label="", title=None):
    ax.set_facecolor(t["surface"])
    ax.tick_params(colors=t["text_secondary"], labelsize=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(t["grid"])
    ax.grid(True, color=t["grid"], linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    if x_label:
        ax.set_xlabel(x_label, color=t["text_secondary"], fontsize=9)
    if y_label:
        ax.set_ylabel(y_label, color=t["text_secondary"], fontsize=9)
    if title:
        ax.set_title(title, color=t["text_primary"], fontsize=10, loc="left")


def _footnote(fig, text: str, t, width_in: float) -> float:
    """Wrap a footnote to the figure width and return the rect bottom it needs.

    Matplotlib does not wrap `fig.text`, so a long note runs off the canvas
    and is silently cropped in the saved PNG. Wrapping to a character count
    derived from the figure width keeps it inside, and the caller reserves
    the matching strip with `tight_layout(rect=...)`.
    """
    chars = max(40, int(width_in * 20))
    lines = textwrap.wrap(text, chars)
    fig.text(0.015, 0.012, "\n".join(lines), color=t["text_secondary"],
             fontsize=7, va="bottom")
    return 0.022 * len(lines) + 0.02


def _roc(pos: list[float], neg: list[float]) -> tuple[list[float], list[float]]:
    pts = sorted([(z, 1) for z in pos] + [(z, 0) for z in neg], reverse=True)
    tp = fp = 0
    xs, ys = [0.0], [0.0]
    for _, g in pts:
        tp += g
        fp += 1 - g
        xs.append(fp / len(neg))
        ys.append(tp / len(pos))
    return xs, ys


def plot_roc(pos, neg, summary, out_path: Path, theme: str) -> Path:
    t = THEMES[theme]
    xs, ys = _roc(pos, neg)
    fig, ax = plt.subplots(figsize=(5.2, 5.0), facecolor=t["surface"])
    ax.plot([0, 1], [0, 1], color=t["grid"], linewidth=1.5,
            linestyle=(0, (4, 3)), zorder=1)
    ax.plot(xs, ys, color=t["series"], linewidth=2.0, zorder=3)
    a = summary["pooled_auc"]
    ci = summary.get("pooled_auc_ci95")
    null = summary.get("pooled_auc_null_mean")
    sub = f"pooled AUC {a:.3f}"
    if ci:
        sub += f", 95% CI [{ci[0]:.3f}, {ci[1]:.3f}]"
    _style(ax, t, x_label="false positive rate", y_label="true positive rate")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    fig.suptitle("Pooled ROC, epitope vs other surface residue",
                 color=t["text_primary"], fontsize=11, x=0.02, ha="left")
    ax.text(0.97, 0.10, sub, transform=ax.transAxes, ha="right",
            color=t["text_primary"], fontsize=9)
    if null is not None:
        ax.text(0.97, 0.04,
                f"within-structure null {null:.3f}  (p = "
                f"{summary['pooled_auc_permutation_p']:.2f})",
                transform=ax.transAxes, ha="right",
                color=t["summary"], fontsize=9)
    # The diagonal is the 0.5 reference, which is NOT this statistic's null:
    # pooling residues across structures of differing composition buys AUC on
    # its own. Saying so on the figure stops the diagonal being read as the
    # bar the curve has to clear.
    note = ("Interval resampled over STRUCTURES, the unit that was sampled; "
            "residue-level resampling would be anticonservative. The dashed "
            "diagonal marks 0.5, which is not the null here -- pooling across "
            "structures of differing score and epitope composition yields "
            f"{null:.3f} on permuted labels, so that is the bar."
            if null is not None else
            "Interval resampled over STRUCTURES, the unit that was sampled.")
    bottom = _footnote(fig, note, t, 5.2)
    fig.tight_layout(rect=(0, bottom, 1, 0.96))
    fig.savefig(out_path, dpi=200, facecolor=t["surface"])
    plt.close(fig)
    return out_path


def plot_sweep(sweep, base, out_path: Path, theme: str) -> Path:
    t = THEMES[theme]
    metrics = [("precision", "precision"), ("recall", "recall"),
               ("f1", "F1"), ("mcc", "MCC")]
    fig, axes = plt.subplots(2, 2, figsize=(8.0, 6.0), facecolor=t["surface"])
    xs = [d["threshold"] for d in sweep]
    for ax, (key, label) in zip(axes.flat, metrics):
        ys = [d[key] if d[key] is not None else float("nan") for d in sweep]
        ax.plot(xs, ys, color=t["series"], linewidth=2.0,
                marker="o", markersize=4.5, zorder=3)
        if key == "precision":
            ax.axhline(base, color=t["summary"], linewidth=1.5,
                       linestyle=(0, (4, 3)), zorder=2)
            ax.text(0.99, base, f"base rate {base:.3f}",
                    transform=ax.get_yaxis_transform(),
                    va="bottom", ha="right",
                    color=t["summary"], fontsize=7.5)
        if key == "mcc":
            ax.axhline(0.0, color=t["grid"], linewidth=1.2, zorder=2)
        _style(ax, t, x_label="call threshold (z)", y_label=label, title=label)
    fig.suptitle("Threshold sweep: every metric against where the cut is placed",
                 color=t["text_primary"], fontsize=11, x=0.02, ha="left")
    bottom = _footnote(fig,
                       "One panel per metric rather than four hues on one "
                       "axis: the metrics do not share a range and four "
                       "all-pairs hues do not clear the palette's separation "
                       "floor.", t, 8.0)
    fig.tight_layout(rect=(0, bottom, 1, 0.955))
    fig.savefig(out_path, dpi=200, facecolor=t["surface"])
    plt.close(fig)
    return out_path


def plot_per_structure(per_structure, out_path: Path, theme: str) -> Path:
    t = THEMES[theme]
    rows = [d for d in per_structure if d["ci_lo"] is not None]
    rows.sort(key=lambda d: d["auc"])
    n = len(rows)
    fig, ax = plt.subplots(figsize=(6.4, max(4.0, 0.22 * n + 1.4)),
                           facecolor=t["surface"])
    ys = range(n)
    for y, d in zip(ys, rows):
        ax.plot([d["ci_lo"], d["ci_hi"]], [y, y],
                color=t["variance"], linewidth=1.4, alpha=0.55, zorder=2)
        ax.plot([d["auc"]], [y], marker="o", markersize=5,
                color=t["series"], zorder=3,
                markeredgecolor=t["surface"], markeredgewidth=1.2)
    ax.axvline(0.5, color=t["summary"], linewidth=1.5,
               linestyle=(0, (4, 3)), zorder=1)
    ax.set_yticks(list(ys))
    ax.set_yticklabels([f"{d['instance'].replace('pdb_0000', '')} "
                        f"(n+={d['n_epitope']})" for d in rows], fontsize=6.5)
    _style(ax, t, x_label="AUC, epitope vs other", y_label="")
    ax.set_xlim(0, 1)
    fig.suptitle("Per-structure AUC with its bootstrap interval",
                 color=t["text_primary"], fontsize=11, x=0.02, ha="left")
    bottom = _footnote(fig,
                       "The intervals overlap almost everything. No ranking "
                       "from this panel is supported; only the pooled figure "
                       "has usable precision.", t, 6.4)
    fig.tight_layout(rect=(0, bottom, 1, 0.965))
    fig.savefig(out_path, dpi=200, facecolor=t["surface"])
    plt.close(fig)
    return out_path


def plot_coverage(cov_rows, scores_by_k, out_path: Path, theme: str) -> Path:
    """Observed mean score against the order-statistic null, per coverage k."""
    t = THEMES[theme]
    ks = [r["k"] for r in cov_rows if r["k"] in scores_by_k and scores_by_k[r["k"]]]
    if not ks:
        return out_path
    obs = [sum(scores_by_k[k]) / len(scores_by_k[k]) for k in ks]
    null = [next(r["null"] for r in cov_rows if r["k"] == k) for k in ks]
    fig, ax = plt.subplots(figsize=(6.4, 4.6), facecolor=t["surface"])
    ax.plot(ks, null, color=t["summary"], linewidth=2.0, marker="s",
            markersize=5, zorder=2, label="E[max | k], order-statistic null")
    ax.plot(ks, obs, color=t["series"], linewidth=2.0, marker="o",
            markersize=5.5, zorder=3, label="observed mean score")
    for k, o, nl in zip(ks, obs, null):
        if k in (ks[0], ks[-1]):
            ax.annotate(f"{o:+.2f}", (k, o), textcoords="offset points",
                        xytext=(0, 8), ha="center", fontsize=7.5,
                        color=t["text_primary"])
    leg = ax.legend(frameon=False, fontsize=8, loc="lower right")
    for txt in leg.get_texts():
        txt.set_color(t["text_primary"])
    _style(ax, t, x_label="motifs containing the residue (k)",
           y_label="max motif z")
    fig.suptitle("The coverage artifact: how much of the score is counting",
                 color=t["text_primary"], fontsize=11, x=0.02, ha="left")
    bottom = _footnote(fig,
                       "A residue in k motifs takes the max of k draws, so "
                       "the score rises with k even when nothing differs. "
                       "Where observed tracks the null, the rise carries no "
                       "information beyond the counting.", t, 6.4)
    fig.tight_layout(rect=(0, bottom, 1, 0.955))
    fig.savefig(out_path, dpi=200, facecolor=t["surface"])
    plt.close(fig)
    return out_path


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--stats-dir", type=Path, default=None,
                    help="[<run_dir>/01_statistics]")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="[<run_dir>/00_figures]")
    ap.add_argument("--theme", default="light", choices=("light", "dark"))
    args = ap.parse_args(argv)

    stats_dir = args.stats_dir or (args.run_dir / "01_statistics")
    out_dir = args.out_dir or (args.run_dir / "00_figures")
    out_dir.mkdir(parents=True, exist_ok=True)

    sfile = stats_dir / "summary.json"
    if not sfile.is_file():
        print(f"error: {sfile} missing; run epitope_statistics.py first",
              file=sys.stderr)
        return 2
    blob = json.loads(sfile.read_text())
    summary, sweep = blob["summary"], blob["threshold_sweep"]
    per_structure = blob["per_structure"]

    pos, neg = [], []
    for r in csv.DictReader((args.run_dir / "residue_scores.tsv").open(),
                            delimiter="\t"):
        if r["is_epitope"] == "":
            continue
        (pos if r["is_epitope"] == "1" else neg).append(float(r["max_motif_z"]))

    suffix = "" if args.theme == "light" else "_dark"
    written = [
        plot_roc(pos, neg, summary, out_dir / f"roc_pooled{suffix}.png", args.theme),
        plot_sweep(sweep, summary["epitope_base_rate"],
                   out_dir / f"threshold_sweep{suffix}.png", args.theme),
        plot_per_structure(per_structure,
                           out_dir / f"per_structure_auc_ci{suffix}.png",
                           args.theme),
    ]

    cov_path = stats_dir / "coverage_artifact.tsv"
    if cov_path.is_file():
        cov_rows = [{"k": int(r["k_motifs"]), "null": float(r["expected_max_null"])}
                    for r in csv.DictReader(cov_path.open(), delimiter="\t")]
        # Re-derive observed score by k from the residue table joined to the
        # per-k counts the statistics step wrote.
        by_k: dict[int, list[float]] = defaultdict(list)
        counts = {r["k"]: r for r in cov_rows}
        srows = [r for r in csv.DictReader(
            (args.run_dir / "residue_scores.tsv").open(), delimiter="\t")
            if r["is_epitope"] != ""]
        kmap = _coverage_from_run(args.run_dir, {r["instance"] for r in srows})
        for r in srows:
            k = kmap.get((r["instance"], r["chain"], int(r["res_id"])))
            if k and k in counts:
                by_k[k].append(float(r["max_motif_z"]))
        written.append(plot_coverage(
            cov_rows, by_k, out_dir / f"coverage_artifact{suffix}.png",
            args.theme))

    for p in written:
        print(f"wrote {p}")
    return 0


def _coverage_from_run(run_dir: Path, instances) -> dict[tuple[str, str, int], int]:
    """Motifs-per-residue, counting only motifs that carry a z."""
    cover: dict[tuple[str, str, int], int] = defaultdict(int)
    for inst in instances:
        d = run_dir / inst
        scored, index = d / "idf_scored.tsv", d / "folddisco_queries.index.tsv"
        members = d / "patch_members.tsv"
        if not (scored.is_file() and index.is_file() and members.is_file()):
            continue
        by_query = {}
        for r in csv.DictReader(scored.open(), delimiter="\t"):
            try:
                by_query[r["query_residues"]] = float(r["idf_z"])
            except (ValueError, KeyError):
                continue
        keep = {r["patch_id"] for r in csv.DictReader(index.open(), delimiter="\t")
                if r["query_string"] in by_query}
        for r in csv.DictReader(members.open(), delimiter="\t"):
            if r["patch_id"] in keep:
                cover[(inst, r["auth_chain_id"], int(r["auth_seq_id"]))] += 1
    return dict(cover)


if __name__ == "__main__":
    sys.exit(main())
