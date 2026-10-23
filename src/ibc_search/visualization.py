"""Plots for motif scores. Figures only -- no scoring logic lives here.

The one form implemented so far is a **binned point plot**: one bin per
structure along x, the per-motif score on y, every motif drawn as a point
with the median marked.

Why it is faceted by category rather than coloured by it
--------------------------------------------------------
Four categories on one panel would need four categorical hues compared
against each other in any pairing. The validated palette clears its gates for
*three* slots under all-pairs comparison; the fourth puts yellow beside
orange, which fails both the normal-vision and colour-vision-deficiency
floors. The documented remedy is to facet instead of adding a hue, so each
category gets its own panel and a single series -- identity comes from the
panel title, not from colour, and the figure is readable in greyscale and
under every form of colour blindness.

Density
-------
A structure can contribute well over a thousand motifs, which as opaque dots
is a solid block. Points are drawn small, semi-transparent and jittered, with
the median and interquartile range overlaid so the summary survives the
overplotting. The point cloud shows the spread; the bars carry the reading.
"""

from __future__ import annotations

import argparse
import csv
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Sequence

# Light and dark are each selected against their own surface, not an
# automatic inversion of one another.
THEMES = {
    "light": {
        "surface": "#fcfcfb",
        "text_primary": "#0b0b0b",
        "text_secondary": "#52514e",
        "grid": "#d9d8d4",
        "series": "#2a78d6",
        "summary": "#e34948",
    },
    "dark": {
        "surface": "#1a1a19",
        "text_primary": "#f0efec",
        "text_secondary": "#c3c2b7",
        "grid": "#3a3a38",
        "series": "#3987e5",
        "summary": "#e66767",
    },
}


def binned_point_plot(
    groups: Sequence[tuple[str, Sequence[tuple[str, Sequence[float]]]]],
    out_path: Path,
    *,
    y_label: str = "score",
    title: str | None = None,
    subtitle: str | None = None,
    reference: float | None = 0.0,
    reference_label: str = "reference mean",
    summary: str = "top_decile",
    top_fraction: float = 0.10,
    theme: str = "light",
    max_points: int = 400,
    seed: int = 0,
) -> Path:
    """Draw one panel per group, one bin per structure within it.

    ``groups`` is ``[(group_name, [(bin_name, values), ...]), ...]``. Panels
    share a y-axis so bins are comparable across groups, which is the whole
    point of the figure.

    ``max_points`` subsamples the cloud per bin for legibility; the median and
    IQR are always computed on the **full** set, so the summary never depends
    on the subsample.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    palette = THEMES[theme]
    rng = np.random.default_rng(seed)

    widths = [max(len(bins), 1) for _, bins in groups]
    fig, axes = plt.subplots(
        1, len(groups), figsize=(2.4 + 0.62 * sum(widths), 5.0),
        sharey=True, gridspec_kw={"width_ratios": widths},
    )
    if len(groups) == 1:
        axes = [axes]
    fig.patch.set_facecolor(palette["surface"])

    for ax, (group_name, bins) in zip(axes, groups):
        ax.set_facecolor(palette["surface"])
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color(palette["grid"])
        ax.tick_params(colors=palette["text_secondary"], length=0)
        ax.yaxis.grid(True, color=palette["grid"], linewidth=0.6, alpha=0.7)
        ax.set_axisbelow(True)

        if reference is not None:
            ax.axhline(reference, color=palette["text_secondary"],
                       linewidth=1.0, linestyle=(0, (4, 3)), alpha=0.8, zorder=1)

        for i, (bin_name, values) in enumerate(bins):
            vals = np.asarray([v for v in values], dtype=float)
            if vals.size == 0:
                continue
            shown = vals if vals.size <= max_points else rng.choice(
                vals, max_points, replace=False
            )
            # Jitter is cosmetic spread only; it carries no information.
            x = i + rng.uniform(-0.17, 0.17, size=shown.size)
            ax.scatter(x, shown, s=5, color=palette["series"],
                       alpha=0.30, linewidths=0, zorder=2, rasterized=True)

            # The bulk of the distribution is not the quantity of interest
            # when only a few motifs are expected to be real candidates, so
            # the default summary is the top decile: the bar is its median
            # (about the 95th percentile, which unlike the maximum does not
            # drift upward simply because a structure has more motifs) and
            # the whisker spans the decile itself.
            if summary == "top_decile":
                cut = float(np.percentile(vals, 100 * (1 - top_fraction)))
                tail = vals[vals >= cut]
                bar = float(np.median(tail)) if tail.size else float(np.max(vals))
                lo, hi = cut, float(np.max(vals))
            else:
                bar = float(np.median(vals))
                lo = float(np.percentile(vals, 25))
                hi = float(np.percentile(vals, 75))
            ax.hlines(bar, i - 0.32, i + 0.26, color=palette["summary"],
                      linewidth=2.4, zorder=4)
            # The value sits beside the line it belongs to, so the reading
            # does not depend on tracing a mark back to the y-axis.
            ax.annotate(
                f"{bar:+.2f}", xy=(i + 0.28, bar), xytext=(1, 0),
                textcoords="offset points", va="center", ha="left",
                fontsize=8.2, color=palette["text_primary"], zorder=5,
                annotation_clip=False,
            )

        ax.set_xticks(range(len(bins)))
        ax.set_xticklabels(
            [b for b, _ in bins], rotation=45, ha="right", fontsize=9,
            color=palette["text_primary"],
        )
        ax.set_xlim(-0.6, len(bins) - 0.12)
        ax.set_title(group_name, fontsize=10.5, color=palette["text_primary"],
                     pad=9)
        counts = sum(len(v) for _, v in bins)
        ax.set_xlabel(f"n = {counts:,}", fontsize=8.5,
                      color=palette["text_secondary"], labelpad=7)

    axes[0].set_ylabel(y_label, fontsize=10, color=palette["text_primary"])

    if title:
        fig.text(0.012, 0.962, title, fontsize=13,
                 color=palette["text_primary"], ha="left", va="top")
    if subtitle:
        fig.text(0.012, 0.902, subtitle, fontsize=9.5,
                 color=palette["text_secondary"], ha="left", va="top")
    if reference is not None:
        note = (f"red line: median of the top {top_fraction:.0%} of motifs"
                if summary == "top_decile" else "red line: median")
        fig.text(0.988, 0.902, f"dashed line: {reference_label}   {note}",
                 fontsize=8.5, color=palette["text_secondary"], ha="right",
                 va="top")

    fig.tight_layout(rect=(0, 0, 1, 0.86 if (title or subtitle) else 1))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, facecolor=palette["surface"])
    plt.close(fig)
    return out_path


# ---------------------------------------------------------------------------
# Loading scored motifs
# ---------------------------------------------------------------------------


def load_scored(
    panel_tsv: Path, results_root: Path, score_file: str, score_column: str
) -> list[tuple[str, list[tuple[str, list[float]]]]]:
    """Read one scored file per structure, grouped by the panel's category."""
    by_group: dict[str, list[tuple[str, list[float]]]] = defaultdict(list)
    order: list[str] = []
    missing: list[str] = []

    for row in csv.DictReader(panel_tsv.open(), delimiter="\t"):
        pdb_id, category = row["pdb_id"], row["category"]
        path = results_root / pdb_id / score_file
        if not path.is_file():
            missing.append(pdb_id)
            continue
        values = []
        for r in csv.DictReader(path.open(), delimiter="\t"):
            try:
                values.append(float(r[score_column]))
            except (ValueError, KeyError):
                continue
        if category not in order:
            order.append(category)
        by_group[category].append((pdb_id, values))

    if missing:
        print(f"warning: no {score_file} for {', '.join(missing)}", file=sys.stderr)
    return [(g, by_group[g]) for g in order]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("panel", type=Path, help="TSV with pdb_id and category columns")
    ap.add_argument("results_root", type=Path, help="directory of per-structure results")
    ap.add_argument("-o", "--out", type=Path, required=True)
    ap.add_argument("--score-file", default="idf_scored.tsv")
    ap.add_argument("--score-column", default="idf_z")
    ap.add_argument("--y-label", default="IDF z-score (size-normalised)")
    ap.add_argument("--title", default=None)
    ap.add_argument("--subtitle", default=None)
    ap.add_argument("--theme", choices=tuple(THEMES), default="light")
    ap.add_argument("--reference", type=float, default=0.0)
    ap.add_argument("--max-points", type=int, default=400)
    ap.add_argument("--summary", choices=("top_decile", "median"),
                    default="top_decile",
                    help="what the bar marks [top_decile]")
    args = ap.parse_args(argv)

    groups = load_scored(args.panel, args.results_root,
                         args.score_file, args.score_column)
    if not groups:
        print("error: nothing to plot", file=sys.stderr)
        return 2
    out = binned_point_plot(
        groups, args.out, y_label=args.y_label, title=args.title,
        subtitle=args.subtitle, theme=args.theme, reference=args.reference,
        max_points=args.max_points, summary=args.summary,
    )
    total = sum(len(v) for _, bins in groups for _, v in bins)
    print(f"wrote {out}  ({len(groups)} panels, "
          f"{sum(len(b) for _, b in groups)} bins, {total:,} motifs)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
