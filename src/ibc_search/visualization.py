"""Plots for motif scores. Figures only -- no scoring logic lives here.

Three forms, all point plots, because the spread of a structure's motif
scores is the thing worth seeing and a box or a bar hides it:

``binned_point_plot``
    one bin per structure, every motif a point, median and IQR overlaid.
``paired_point_plot``
    two labelled classes offset within each bin, nothing subsampled.
``value_strip_plot``
    one value per structure against a reference line.

Rendering a *structure* is not done here -- that is
``scripts/render_pymol.py``, which drives PyMOL. A matplotlib approximation
of a molecular surface used to live in this module as a fallback for machines
without PyMOL; it was removed once PyMOL became the renderer, because a
space-filling projection of heavy atoms reads like a surface without being
one, and two renderers disagreeing about the same structure is worse than
one dependency.

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
import sys
from collections import defaultdict
from pathlib import Path
from typing import Sequence

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
        "variance": "#0b0b0b",
    },
    "dark": {
        "surface": "#1a1a19",
        "text_primary": "#f0efec",
        "text_secondary": "#c3c2b7",
        "grid": "#3a3a38",
        "series": "#3987e5",
        "summary": "#e66767",
        "variance": "#ffffff",
    },
}


def _wrap_to_width(text: str, fig_width_inches: float,
                   fontsize: float, left_margin: float = 0.03) -> str:
    """Wrap a caption to the figure's own width.

    Figure text is placed in figure coordinates and matplotlib does not wrap
    it, so a long footnote silently runs past the right edge -- and the part
    that disappears is the end of the sentence, which in this repository is
    reliably the caveat. The character budget comes from an average glyph
    width of ~0.55 em, which is approximate but only ever errs toward
    wrapping early.
    """
    import textwrap

    usable_inches = max(fig_width_inches * (1 - left_margin) - 0.2, 1.0)
    chars = max(int(usable_inches * 72 / (fontsize * 0.55)), 20)
    return "\n".join(
        line for para in text.split("\n")
        for line in textwrap.wrap(para, chars) or [""]
    )


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
    footnote: str | None = None,
    top_fraction: float = 0.10,
    theme: str = "light",
    max_points: int = 400,
    seed: int = 0,
    bin_width: float = 0.62,
    height: float = 5.0,
    min_for_summary: int = 5,
    trend: Sequence[float] | None = None,
    trend_label: str | None = None,
    summary_ink: str = "summary",
    highlight: Sequence[Sequence[Sequence[bool]]] | None = None,
    highlight_label: str = "",
) -> Path:
    """Draw one panel per group, one bin per structure within it.

    ``groups`` is ``[(group_name, [(bin_name, values), ...]), ...]``. Panels
    share a y-axis so bins are comparable across groups, which is the whole
    point of the figure.

    ``max_points`` subsamples the cloud per bin for legibility; the median and
    IQR are always computed on the **full** set, so the summary never depends
    on the subsample.

    ``summary_ink`` picks the palette key for the summary line. The default
    ``"summary"`` is the red used by the exploratory figures. The standard
    run figure passes ``"variance"``, which draws the median in the same ink
    as its spread bar -- black on the light surface, white on the dark one.
    "Black" cannot be literal in both: on the dark surface a black line is
    invisible, so the ink is the theme's, and the caption says which.

    ``highlight`` optionally marks a SUBSET of each bin's points to draw in
    the summary ink (red) on top of the rest, shaped like ``groups`` -- one
    boolean per value. It exists so a figure can show where a labelled class
    sits inside its own structure's distribution. Highlighted points bypass
    the ``max_points`` subsample: they are the reason the figure is being
    drawn, so losing them to a legibility cap would defeat it. The summary
    bars are unaffected and still describe the FULL bin, highlighted and not
    together, so the red marks can be read against the whole.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    palette = THEMES[theme]
    rng = np.random.default_rng(seed)

    # `bin_width` and `height` exist because the default inches-per-bin is
    # sized for a handful of structures; at 78 it produces a figure ten times
    # wider than it is tall, which is unreadable at any zoom.
    widths = [max(len(bins), 1) for _, bins in groups]
    fig, axes = plt.subplots(
        1, len(groups), figsize=(2.4 + bin_width * sum(widths), height),
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

        group_hl = None
        if highlight is not None:
            gi = [g for g, _ in groups].index(group_name)
            group_hl = highlight[gi]

        for i, (bin_name, values) in enumerate(bins):
            vals = np.asarray([v for v in values], dtype=float)
            if vals.size == 0:
                continue
            flags = None
            if group_hl is not None and i < len(group_hl):
                flags = np.asarray(group_hl[i], dtype=bool)
                if flags.size != vals.size:
                    raise ValueError(
                        f"highlight for bin {bin_name!r} has {flags.size} "
                        f"flags for {vals.size} values")
            # Subsample only the UNhighlighted points: the highlighted ones
            # are the subject of the figure and are always all drawn.
            base = vals if flags is None else vals[~flags]
            shown = base if base.size <= max_points else rng.choice(
                base, max_points, replace=False
            )
            # A bin holding a handful of values is a set of observations, not
            # a distribution: draw them as individual marks and skip the
            # summary bars below, which would otherwise draw a median and an
            # interquartile range over a single point.
            sparse = vals.size < min_for_summary
            # Jitter is cosmetic spread only; it carries no information.
            x = i + rng.uniform(-0.17, 0.17, size=shown.size)
            ax.scatter(x, shown,
                       s=34 if sparse else 5,
                       color=palette["summary"] if sparse else palette["series"],
                       alpha=0.95 if sparse else 0.30,
                       linewidths=0, zorder=4 if sparse else 2,
                       rasterized=not sparse)
            if flags is not None and flags.any():
                hv = vals[flags]
                hx = i + rng.uniform(-0.17, 0.17, size=hv.size)
                ax.scatter(hx, hv, s=18, color=palette["summary"],
                           alpha=0.95, linewidths=0, zorder=5,
                           label=(highlight_label
                                  if (highlight_label and i == 0
                                      and ax is axes[0]) else None))
            if sparse:
                continue

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
            # Thin white spread bar, drawn under the summary line. It reads
            # against the point cloud it overlays; where the cloud is sparse
            # on a light surface it is deliberately faint.
            ax.vlines(i, lo, hi, color=palette["variance"], linewidth=1.0,
                      alpha=0.95, zorder=3)
            # Caps mark where the span ends. Without them the line fades into
            # the point cloud and its extent has to be guessed.
            ax.hlines([lo, hi], i - 0.10, i + 0.10, color=palette["variance"],
                      linewidth=1.0, alpha=0.95, zorder=3)
            ax.hlines(bar, i - 0.30, i + 0.30, color=palette[summary_ink],
                      linewidth=2.4, zorder=4)

        # A fitted curve is drawn over the bins it was fitted on, so the
        # reader can see the model and its residuals in one place rather
        # than being told the fit is good.
        if trend is not None and len(trend) == len(bins):
            ax.plot(range(len(bins)), list(trend), color=palette["summary"],
                    linewidth=1.6, linestyle=(0, (5, 2)), zorder=5,
                    label=trend_label)
            if trend_label:
                ax.legend(loc="upper left", frameon=False, fontsize=8.5,
                          labelcolor=palette["text_secondary"])

        ax.set_xticks(range(len(bins)))
        ax.set_xticklabels(
            [b for b, _ in bins], rotation=45, ha="right", fontsize=9,
            color=palette["text_primary"],
        )
        ax.set_xlim(-0.6, len(bins) - 0.4)
        ax.set_title(group_name, fontsize=10.5, color=palette["text_primary"],
                     pad=9)
        counts = sum(len(v) for _, v in bins)
        ax.set_xlabel(f"n = {counts:,}", fontsize=8.5,
                      color=palette["text_secondary"], labelpad=7)

    axes[0].set_ylabel(y_label, fontsize=10, color=palette["text_primary"])

    if title:
        fig.text(0.012, 0.962, title, fontsize=13,
                 color=palette["text_primary"], ha="left", va="top")
    # The legend note goes on its own line under the subtitle, never
    # right-aligned on the same row: at that position a long subtitle and a
    # long note overprint each other, and which one wins depends on the
    # caller's text length rather than on anything the figure controls.
    lines = []
    if subtitle:
        lines.append(subtitle)
    # Colour-neutral wording: the spread bar is white on the dark surface and
    # black on the light one.
    ink = "red" if summary_ink == "summary" else (
        "black" if theme == "light" else "white")
    note = (f"{ink} line: median of the top {top_fraction:.0%} of motifs; "
            f"bar spans that decile"
            if summary == "top_decile"
            else f"{ink} line: median; vertical bar spans the IQR, with caps")
    if reference is not None:
        note = f"dashed line: {reference_label}   {note}"
    lines.append(note)
    if lines:
        fig.text(0.012, 0.902, "\n".join(lines), fontsize=9.0,
                 color=palette["text_secondary"], ha="left", va="top")

    if footnote:
        # Stated on the figure rather than only in the caption: a z-score
        # computed on a transformed quantity is not the same statistic as one
        # computed on the raw values, and a reader cannot tell which they are
        # looking at from the axis alone. Wrapped to the canvas, because an
        # unwrapped footnote runs off the right edge and the part that gets
        # cut is the caveat.
        footnote = _wrap_to_width(footnote, fig.get_figwidth(), fontsize=8.2)
        fig.text(0.012, 0.016, footnote, fontsize=8.2,
                 color=palette["text_secondary"], ha="left", va="bottom")

    fig.tight_layout(rect=(
        0, 0.034 * (footnote.count("\n") + 1) + 0.022 if footnote else 0, 1,
        0.845 if (title or subtitle or lines) else 1))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, facecolor=palette["surface"])
    plt.close(fig)
    return out_path


def paired_point_plot(
    bins: Sequence[tuple[str, Sequence[float], Sequence[float]]],
    out_path: Path,
    *,
    labels: tuple[str, str] = ("epitope", "other"),
    y_label: str = "score",
    title: str | None = None,
    subtitle: str | None = None,
    footnote: str | None = None,
    reference: float | None = None,
    theme: str = "light",
    min_for_median: int = 3,
) -> Path:
    """Two classes side by side within each bin, every point drawn.

    For a comparison between two labelled sets the faceted form hides the
    thing being asked: whether the distributions differ. Offsetting them
    inside one bin puts them against the same baseline, which is what makes
    an overlap -- or its absence -- immediately legible.

    Nothing is subsampled. The point of the figure is the shape of the two
    clouds, and thinning one of them would misrepresent it.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    palette = THEMES[theme]
    colours = ("#e34948", "#2a78d6") if theme == "light" else ("#e66767", "#3987e5")
    rng = np.random.default_rng(0)

    fig, ax = plt.subplots(figsize=(2.6 + 0.44 * len(bins), 5.6))
    fig.patch.set_facecolor(palette["surface"])
    ax.set_facecolor(palette["surface"])
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(palette["grid"])
    ax.tick_params(colors=palette["text_secondary"], length=0)
    ax.yaxis.grid(True, color=palette["grid"], linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    if reference is not None:
        ax.axhline(reference, color=palette["text_secondary"], linewidth=1.0,
                   linestyle=(0, (4, 3)), alpha=0.8, zorder=1)

    for i, (name, left, right) in enumerate(bins):
        for offset, vals, colour in ((-0.19, left, colours[0]),
                                     (0.19, right, colours[1])):
            v = np.asarray(list(vals), dtype=float)
            if v.size == 0:
                continue
            # A handful of values is a set of observations, not a
            # distribution: draw them as individual marks, and draw NO median
            # bar. A bold bar over two points lands at a value neither of
            # them has and reads as a summary of something that was measured.
            sparse = v.size < min_for_median
            x = i + offset + rng.uniform(-0.085, 0.085, size=v.size)
            ax.scatter(x, v,
                       s=30 if sparse else 5,
                       color=colour,
                       alpha=0.95 if sparse else 0.38,
                       linewidths=0, zorder=3 if sparse else 2,
                       rasterized=not sparse)
            if not sparse:
                ax.hlines(float(np.median(v)),
                          i + offset - 0.14, i + offset + 0.14,
                          color=palette["summary"], linewidth=1.8, zorder=4)

    ax.set_xticks(range(len(bins)))
    ax.set_xticklabels([b[0] for b in bins], rotation=90, fontsize=7,
                       color=palette["text_primary"])
    ax.set_xlim(-0.7, len(bins) - 0.3)
    ax.set_ylabel(y_label, fontsize=10, color=palette["text_primary"])

    handles = [plt.Line2D([], [], marker="o", linestyle="", markersize=6,
                          color=c, label=l) for c, l in zip(colours, labels)]
    handles.append(plt.Line2D([], [], color=palette["summary"], linewidth=1.8,
                              label="median"))
    ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=8.5,
              labelcolor=palette["text_secondary"], ncol=3)

    if title:
        fig.text(0.008, 0.975, title, fontsize=13,
                 color=palette["text_primary"], ha="left", va="top")
    if subtitle:
        fig.text(0.008, 0.936, subtitle, fontsize=9.5,
                 color=palette["text_secondary"], ha="left", va="top")
    if footnote:
        fig.text(0.008, 0.012, footnote, fontsize=8.2,
                 color=palette["text_secondary"], ha="left", va="bottom")

    # Reserve height per footnote line, so a wrapped method note is not
    # clipped off the bottom of the figure.
    fig.tight_layout(rect=(0, 0.035 * (footnote.count("\n") + 1) if footnote
                           else 0, 1, 0.90 if (title or subtitle) else 1))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, facecolor=palette["surface"])
    plt.close(fig)
    return out_path


def value_strip_plot(
    values: Sequence[tuple[str, float]],
    out_path: Path,
    *,
    reference: float = 0.5,
    reference_label: str = "chance",
    y_label: str = "value",
    title: str | None = None,
    subtitle: str | None = None,
    theme: str = "light",
) -> Path:
    """One value per item -- e.g. a per-structure AUC -- against a reference."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    palette = THEMES[theme]
    vals = sorted(values, key=lambda kv: kv[1])
    fig, ax = plt.subplots(figsize=(2.4 + 0.3 * len(vals), 4.6))
    fig.patch.set_facecolor(palette["surface"])
    ax.set_facecolor(palette["surface"])
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(palette["grid"])
    ax.tick_params(colors=palette["text_secondary"], length=0)
    ax.yaxis.grid(True, color=palette["grid"], linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)

    ax.axhline(reference, color="#e34948", linewidth=1.4,
               linestyle=(0, (4, 3)), zorder=1)
    xs = np.arange(len(vals))
    ys = np.array([v for _, v in vals])
    ax.scatter(xs, ys, s=34, color=palette["series"], alpha=0.9,
               linewidths=0, zorder=3)
    ax.set_xticks(xs)
    ax.set_xticklabels([k for k, _ in vals], rotation=90, fontsize=7,
                       color=palette["text_primary"])
    ax.set_xlim(-0.8, len(vals) - 0.2)
    ax.set_ylabel(y_label, fontsize=10, color=palette["text_primary"])
    ax.text(len(vals) - 0.4, reference, f"  {reference_label}", fontsize=8.5,
            color="#e34948", va="center", ha="left")

    above = sum(1 for _, v in vals if v > reference)
    if title:
        fig.text(0.008, 0.975, title, fontsize=13,
                 color=palette["text_primary"], ha="left", va="top")
    fig.text(0.008, 0.936,
             (subtitle + "   " if subtitle else "")
             + f"{above} of {len(vals)} above {reference_label}",
             fontsize=9.5, color=palette["text_secondary"], ha="left", va="top")

    fig.tight_layout(rect=(0, 0, 1, 0.90))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, facecolor=palette["surface"])
    plt.close(fig)
    return out_path


COUNT_COLORS = ("#d9f0d3", "#8fd4ae", "#36a9a0", "#1f6fb5", "#143d8f")
"""Low-to-high ramp for hit COUNTS: green -> teal -> blue.

A sequential ramp, deliberately not the diverging red/grey used for
z-scores, because the two quantities are not the same kind of thing. A
z-score has a meaningful centre to diverge from; a count has a floor at zero
and no centre, so its ramp runs one way.

Lightness decreases monotonically along it -- measured, not asserted:
0.884, 0.696, 0.437, 0.416, 0.320 in HLS. Two consequences are deliberate.
The order survives greyscale printing, and "many hits" is the visually
heavy end rather than the faint one, so the eye lands on the high counts.
A green-to-blue ramp also separates under the common forms of colour
blindness, where red-to-green does not. ``tests/test_visualization.py``
checks the monotonicity, so a later edit to these five values cannot
quietly break either property.
"""


def count_colormap(theme: str = "light"):
    """The green->teal->blue ramp as a matplotlib colormap.

    ``under`` is set to the surface colour so a count of zero reads as
    absent rather than as the low end of the scale: a residue assessed and
    matched by nothing is a different statement from a residue that merely
    scored low, and :mod:`count_scoring` keeps them apart.
    """
    from matplotlib.colors import LinearSegmentedColormap

    cmap = LinearSegmentedColormap.from_list("ibc_counts", COUNT_COLORS)
    return cmap.with_extremes(bad=THEMES[theme]["grid"])


def count_heatmap(
    rows: Sequence[tuple[str, Sequence[float]]],
    out_path: Path,
    *,
    value_label: str = "distinct hits",
    title: str | None = None,
    subtitle: str | None = None,
    footnote: str | None = None,
    theme: str = "light",
    max_width: int | None = None,
) -> Path:
    """Structures down, residue position across, colour = count.

    ``rows`` is ``[(structure_name, [count per residue in sequence order]),
    ...]``. Rows are ragged by nature -- chains differ in length -- and the
    short ones are padded with NaN, drawn in the grid colour. Padding is
    **not** drawn as zero: zero means "assessed, nothing matched", and a
    position that does not exist in a chain has no count at all.

    The x axis is residue **position within the row**, not a shared
    alignment. Nothing here aligns sequences, so column *i* of two rows is
    not a homologous site, and the axis label says so.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    palette = THEMES[theme]
    if not rows:
        raise ValueError("no rows to draw")

    width = max(len(v) for _, v in rows)
    if max_width is not None:
        width = min(width, max_width)
    grid = np.full((len(rows), width), np.nan)
    for i, (_, values) in enumerate(rows):
        vals = np.asarray(list(values)[:width], dtype=float)
        grid[i, :vals.size] = vals

    height = max(2.4, 0.22 * len(rows) + 1.9)
    fig, ax = plt.subplots(figsize=(12.0, height))
    fig.patch.set_facecolor(palette["surface"])
    ax.set_facecolor(palette["surface"])

    finite = grid[np.isfinite(grid)]
    vmax = float(finite.max()) if finite.size else 1.0
    im = ax.imshow(grid, aspect="auto", interpolation="nearest",
                   cmap=count_colormap(theme), vmin=0.0, vmax=max(vmax, 1.0))

    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([name for name, _ in rows], fontsize=7.5,
                       color=palette["text_primary"])
    ax.set_xlabel("residue position along the chain (NOT a sequence alignment: "
                  "column i of two rows is not a homologous site)",
                  fontsize=8.5, color=palette["text_secondary"], labelpad=7)
    ax.tick_params(colors=palette["text_secondary"], length=0)
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)

    bar = fig.colorbar(im, ax=ax, pad=0.012, fraction=0.022)
    bar.set_label(f"{value_label}  (green = few, blue = many)", fontsize=8.5,
                  color=palette["text_secondary"])
    bar.ax.tick_params(colors=palette["text_secondary"], length=0)
    bar.outline.set_visible(False)

    head = 0.0
    if title:
        fig.text(0.012, 0.985, title, fontsize=13,
                 color=palette["text_primary"], ha="left", va="top")
        head += 0.055
    if subtitle:
        subtitle = _wrap_to_width(subtitle, fig.get_figwidth(), fontsize=9.0)
        fig.text(0.012, 0.985 - head, subtitle, fontsize=9.0,
                 color=palette["text_secondary"], ha="left", va="top")
        head += 0.030 * (subtitle.count("\n") + 1)
    if footnote:
        footnote = _wrap_to_width(footnote, fig.get_figwidth(), fontsize=8.2)
        fig.text(0.012, 0.012, footnote, fontsize=8.2,
                 color=palette["text_secondary"], ha="left", va="bottom")

    fig.tight_layout(rect=(
        0, 0.030 * (footnote.count("\n") + 1) + 0.015 if footnote else 0,
        1, 1 - head - 0.02 if head else 1))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, facecolor=palette["surface"])
    plt.close(fig)
    return out_path


def scatter_plot(
    points: Sequence[tuple[str, float, float]],
    out_path: Path,
    *,
    x_label: str = "x",
    y_label: str = "y",
    title: str | None = None,
    subtitle: str | None = None,
    footnote: str | None = None,
    annotate: int = 0,
    theme: str = "light",
) -> Path:
    """One point per item, with Pearson and Spearman printed on the figure.

    Both coefficients, always: Pearson answers "is there a linear relation",
    Spearman "is there a monotone one at all". Reporting only the first turns
    a curved but real association into "no correlation", and reporting only
    the second hides that a strong rank correlation rests on three points.

    No trend line is drawn unless the relationship is worth one -- a fitted
    line through an uncorrelated cloud is the most persuasive way to show
    something that is not there.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    palette = THEMES[theme]
    xs = np.array([p[1] for p in points], dtype=float)
    ys = np.array([p[2] for p in points], dtype=float)

    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    fig.patch.set_facecolor(palette["surface"])
    ax.set_facecolor(palette["surface"])
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(palette["grid"])
    ax.tick_params(colors=palette["text_secondary"], length=0)
    ax.grid(True, color=palette["grid"], linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)

    ax.scatter(xs, ys, s=38, color=palette["series"], alpha=0.85,
               linewidths=0, zorder=3)
    ax.set_xlabel(x_label, fontsize=10, color=palette["text_primary"])
    ax.set_ylabel(y_label, fontsize=10, color=palette["text_primary"])

    # Label only the extremes, never every point: a name on all 78 is a
    # grey smear that hides the cloud it is meant to explain.
    if annotate and len(points) > annotate:
        for i in np.argsort(ys)[-annotate:]:
            ax.annotate(points[int(i)][0], (xs[int(i)], ys[int(i)]),
                        fontsize=7, color=palette["text_secondary"],
                        xytext=(4, 3), textcoords="offset points")

    stats_line = ""
    if len(xs) > 2 and xs.std() > 0 and ys.std() > 0:
        pearson = float(np.corrcoef(xs, ys)[0, 1])
        rx = np.argsort(np.argsort(xs)).astype(float)
        ry = np.argsort(np.argsort(ys)).astype(float)
        spearman = float(np.corrcoef(rx, ry)[0, 1])
        stats_line = (f"n = {len(xs)}   Pearson r = {pearson:+.3f}   "
                      f"Spearman rho = {spearman:+.3f}")

    if title:
        fig.text(0.008, 0.975, title, fontsize=13,
                 color=palette["text_primary"], ha="left", va="top")
    fig.text(0.008, 0.936, (subtitle + "   " if subtitle else "") + stats_line,
             fontsize=9.5, color=palette["text_secondary"], ha="left", va="top")
    if footnote:
        fig.text(0.008, 0.012, footnote, fontsize=8.2,
                 color=palette["text_secondary"], ha="left", va="bottom")

    # Reserve height per footnote line, so a multi-line method note is not
    # silently clipped off the bottom of the figure.
    lines = footnote.count("\n") + 1 if footnote else 0
    fig.tight_layout(rect=(0, 0.035 * lines, 1, 0.90))
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
    ap.add_argument("--y-label",
                    default="IDF z-score  (log10 scale, size-normalised)")
    ap.add_argument("--footnote", default=None,
                    help="method note rendered at the foot of the figure")
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
        footnote=args.footnote,
    )
    total = sum(len(v) for _, bins in groups for _, v in bins)
    print(f"wrote {out}  ({len(groups)} panels, "
          f"{sum(len(b) for _, b in groups)} bins, {total:,} motifs)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
