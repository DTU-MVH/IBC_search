"""Is an annotated residue better than the background WITHIN ITS OWN STRUCTURE?

Reads a ``few_balls_search`` run -- the method that places one ball on every
annotated residue and an equal number of dispersed background balls on the
same protein -- and answers that question per annotated residue, with a 1 or
a 0.

Why the comparison is per structure
-----------------------------------
Asking whether annotated motifs match more of a database than background
motifs across a whole panel was MEASURED to be dominated by things that have
nothing to do with the annotation: fold, resolution, deposition year, and
above all how many motifs each residue happened to sit in. Within one
structure none of those differ. Both classes come from the same protein, the
same radius, the same index and the same search, so the only thing separating
them is where their ball centres sit.

That is why the verdict is a within-structure comparison and why it is
reported per residue rather than pooled: pooling reintroduces between-
structure variation through the back door, by letting structures with many
annotated residues dominate the total.

The score being compared
------------------------
One number per ball: the count of **distinct matched sites** its query
returned, a site being one ``(reference structure, matched residue set)``
pair, deduplicated -- the same definition ``count_scoring.py`` uses. A ball is
a motif, so this is a per-motif score and needs no normalisation by search
effort: every ball is one query against one index.

The verdict
-----------
An annotated residue scores **1** when its own ball's count exceeds the
reference point of its structure's background balls, and **0** when it does
not. The reference point defaults to the background **median**, which is
robust to one or two unusual background balls in a set that may be as small
as five. A residue whose structure has no background balls has no verdict at
all and is reported as blank, not as 0: an unanswerable comparison is not a
failed one.

Ties count as 0, i.e. the test is strictly greater. With small integer counts
ties are common, so the rule is stated rather than left to rounding, and the
percentile is reported alongside so a 0 that was a near-miss can be told from
a 0 that was not.

What this does NOT establish
----------------------------
A higher count means the motif's geometry is **commoner in the searched
database**. It is not evidence of antigenicity, function, or that the
annotation is correct. A verdict of 1 says the annotated site's neighbourhood
is more widely shared than its own protein's background, which is a statement
about structural commonness and nothing else.
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "src"))

from count_scoring import load_motif_sites  # noqa: E402

SEED_CLASS_ANNOTATED = "annotated"
SEED_CLASS_BACKGROUND = "background"

VERDICT_HEADER = [
    "instance", "radius_angstrom", "seed_residue", "seed_class",
    "n_distinct_sites", "n_background_balls", "background_median",
    "background_mean", "percentile_within_background", "beats_background",
]


# ---------------------------------------------------------------------------
# reading a run
# ---------------------------------------------------------------------------

def motif_scores(run_dir: Path, inst: str) -> list[dict[str, Any]]:
    """One record per emitted ball of one structure, with its class and score.

    Balls whose query was never emitted are absent rather than zero: they were
    never assessed. Balls that were searched and matched nothing are 0,
    because they were.
    """
    d = run_dir / inst
    patches = d / "patches.tsv"
    if not patches.is_file():
        return []
    sites = load_motif_sites(d)
    out: list[dict[str, Any]] = []
    for r in csv.DictReader(patches.open(), delimiter="\t"):
        if r.get("emit_query", "true").lower() == "false":
            continue
        seed_class = (r.get("seed_class") or "").strip()
        if seed_class not in (SEED_CLASS_ANNOTATED, SEED_CLASS_BACKGROUND):
            # A run from a method that does not label its seeds cannot be
            # compared this way, and silently treating unlabelled balls as
            # background would invent the control.
            raise ValueError(
                f"{patches}: patch {r.get('patch_id')} has seed_class "
                f"{seed_class!r}; this analysis needs a few_balls_search run, "
                "whose balls are labelled 'annotated' or 'background'")
        try:
            radius = float(r["radius_angstrom"])
        except (KeyError, ValueError):
            continue
        out.append({
            "instance": inst,
            "radius": radius,
            "patch_id": r["patch_id"],
            "seed_residue": (r.get("seed_residue_labels") or "").strip(),
            "seed_class": seed_class,
            "n_distinct_sites": len(sites.get(r["patch_id"], ())),
        })
    return out


def read_run(run_dir: Path, panel: Path | None = None
             ) -> tuple[list[dict[str, Any]], list[str]]:
    """Every emitted ball of every panel member that produced output."""
    panel = panel or (run_dir / "panel.tsv")
    if not panel.is_file():
        raise SystemExit(f"no panel at {panel}")
    records: list[dict[str, Any]] = []
    missing: list[str] = []
    for row in csv.DictReader(panel.open(), delimiter="\t"):
        inst = row.get("instance") or row.get("pdb_id") or ""
        if not inst:
            continue
        got = motif_scores(run_dir, inst)
        if not got:
            missing.append(inst)
            continue
        records.extend(got)
    return records, missing


# ---------------------------------------------------------------------------
# the verdict
# ---------------------------------------------------------------------------

def annotated_beats_background(
    records: list[dict[str, Any]], *, reference: str = "median",
) -> list[dict[str, Any]]:
    """Per annotated ball: 1 if it beats its OWN structure's background.

    *records* is what :func:`read_run` returns. The comparison is made within
    each ``(instance, radius)`` group, so a run with several radii yields one
    verdict per radius rather than mixing scales that are not comparable.

    *reference* selects the background reference point: ``"median"`` (default,
    robust on the five-ball minimum) or ``"mean"``.

    Returns one row per annotated ball, with ``beats_background`` set to 1, 0,
    or ``None`` where the structure had no background ball to compare against.
    ``percentile_within_background`` is the fraction of background balls the
    annotated ball scored at least as high as, so a 0 that nearly passed is
    distinguishable from one that did not.
    """
    if reference not in ("median", "mean"):
        raise ValueError(f"reference must be 'median' or 'mean', got "
                         f"{reference!r}")

    groups: dict[tuple[str, float], list[dict[str, Any]]] = defaultdict(list)
    for rec in records:
        groups[(rec["instance"], rec["radius"])].append(rec)

    rows: list[dict[str, Any]] = []
    for (inst, radius), group in sorted(groups.items()):
        bg = [g["n_distinct_sites"] for g in group
              if g["seed_class"] == SEED_CLASS_BACKGROUND]
        ref = None
        if bg:
            ref = st.median(bg) if reference == "median" else st.mean(bg)
        for g in sorted(group, key=lambda x: x["patch_id"]):
            if g["seed_class"] != SEED_CLASS_ANNOTATED:
                continue
            value = g["n_distinct_sites"]
            if ref is None:
                verdict: int | None = None
                pct: float | None = None
            else:
                # Strictly greater: with small integer counts a tie is not
                # evidence of being better.
                verdict = int(value > ref)
                pct = round(sum(1 for b in bg if value >= b) / len(bg), 4)
            rows.append({
                "instance": inst,
                "radius_angstrom": radius,
                "seed_residue": g["seed_residue"],
                "seed_class": g["seed_class"],
                "n_distinct_sites": value,
                "n_background_balls": len(bg),
                "background_median": (round(st.median(bg), 4) if bg else None),
                "background_mean": (round(st.mean(bg), 4) if bg else None),
                "percentile_within_background": pct,
                "beats_background": verdict,
            })
    return rows


# ---------------------------------------------------------------------------
# statistics
# ---------------------------------------------------------------------------

def _describe(values: list[float]) -> dict[str, Any] | None:
    """Mean, variation and distinct counts of one category."""
    if not values:
        return None
    return {
        "n": len(values),
        "mean": round(st.mean(values), 4),
        "sd": round(st.stdev(values), 4) if len(values) > 1 else 0.0,
        "median": round(st.median(values), 4),
        "iqr": [round(q, 4) for q in (
            _percentile(values, 25), _percentile(values, 75))],
        "min": round(min(values), 4),
        "max": round(max(values), 4),
        # "Distinct counts" in two senses, both useful and easily confused:
        # how many different count VALUES occur, and the total of the counts.
        "n_distinct_values": len(set(values)),
        "sum_of_counts": round(sum(values), 4),
    }


def _percentile(values: list[float], pct: float) -> float:
    s = sorted(values)
    if len(s) == 1:
        return float(s[0])
    k = (len(s) - 1) * pct / 100.0
    lo = int(k)
    hi = min(lo + 1, len(s) - 1)
    return float(s[lo] + (s[hi] - s[lo]) * (k - lo))


def statistics(records: list[dict[str, Any]],
               verdicts: list[dict[str, Any]]) -> dict[str, Any]:
    """Per-class and per-verdict statistics, pooled and per structure."""
    by_class: dict[str, list[float]] = defaultdict(list)
    for rec in records:
        by_class[rec["seed_class"]].append(float(rec["n_distinct_sites"]))

    scored = [v for v in verdicts if v["beats_background"] is not None]
    by_verdict: dict[str, list[float]] = defaultdict(list)
    for v in scored:
        key = "beats_background" if v["beats_background"] else "does_not_beat"
        by_verdict[key].append(float(v["n_distinct_sites"]))

    # Per structure, so a handful of heavily annotated structures cannot set
    # the headline on their own.
    per_structure: list[dict[str, Any]] = []
    by_inst: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for v in verdicts:
        by_inst[v["instance"]].append(v)
    for inst in sorted(by_inst):
        rows = [v for v in by_inst[inst] if v["beats_background"] is not None]
        if not rows:
            per_structure.append({"instance": inst, "n_annotated": 0,
                                  "fraction_beating": None})
            continue
        per_structure.append({
            "instance": inst,
            "n_annotated": len(rows),
            "n_beating": sum(v["beats_background"] for v in rows),
            "fraction_beating": round(
                sum(v["beats_background"] for v in rows) / len(rows), 4),
            "background_median": rows[0]["background_median"],
        })

    fracs = [p["fraction_beating"] for p in per_structure
             if p["fraction_beating"] is not None]
    return {
        "score": "distinct matched sites per ball; a site is one (reference "
                 "structure, matched residue set) pair, deduplicated",
        "verdict_rule": "1 when an annotated ball's count is STRICTLY greater "
                        "than the median of its own structure's background "
                        "balls; 0 otherwise; blank where the structure had no "
                        "background ball",
        "by_seed_class": {k: _describe(v) for k, v in sorted(by_class.items())},
        "annotated_by_verdict": {k: _describe(v)
                                 for k, v in sorted(by_verdict.items())},
        "pooled_annotated_residues": len(verdicts),
        "pooled_annotated_with_a_verdict": len(scored),
        "pooled_n_beating": sum(v["beats_background"] for v in scored),
        "pooled_fraction_beating": (
            round(sum(v["beats_background"] for v in scored) / len(scored), 4)
            if scored else None),
        # The per-structure mean treats each structure as one observation,
        # which is the right unit of replication here; the pooled fraction
        # treats each residue as one and is reported beside it rather than
        # instead of it, because the two answer different questions.
        "mean_of_per_structure_fractions": (round(st.mean(fracs), 4)
                                            if fracs else None),
        "sd_of_per_structure_fractions": (round(st.stdev(fracs), 4)
                                          if len(fracs) > 1 else None),
        "n_structures_with_a_verdict": len(fracs),
        "per_structure": per_structure,
    }


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------

def draw_annotated_bins(records: list[dict[str, Any]], fig_dir: Path,
                        *, stem: str = "annotated_bins") -> list[str]:
    """The binned dot map, with the annotated residues drawn in red.

    One bin per structure, one point per ball, every ball of both classes in
    the cloud so the red marks can be read against the whole distribution
    they belong to. The median bar describes the full bin, not the red subset.
    """
    from ibc_search.visualization import binned_point_plot

    by_inst: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for rec in records:
        by_inst[rec["instance"]].append(rec)

    bins: list[tuple[str, list[float]]] = []
    flags: list[list[bool]] = []
    for inst in sorted(by_inst):
        rows = sorted(by_inst[inst], key=lambda r: r["patch_id"])
        bins.append((inst, [float(r["n_distinct_sites"]) for r in rows]))
        flags.append([r["seed_class"] == SEED_CLASS_ANNOTATED for r in rows])
    if not bins:
        return []

    written = []
    for theme in ("light", "dark"):
        out = fig_dir / (f"{stem}.png" if theme == "light"
                         else f"{stem}_dark.png")
        binned_point_plot(
            [("all structures", bins)], out,
            y_label="distinct matched sites per ball",
            title="Annotated residues against their own structure's background",
            subtitle=("one bin per structure, one point per ball; RED = a ball "
                      "centred on an annotated residue, blue = a background "
                      "ball placed on the same protein"),
            reference=None,
            summary="median", summary_ink="variance",
            highlight=[flags], highlight_label="annotated residue",
            theme=theme,
            footnote=("A COUNT, not a z-score: a high value means the ball's "
                      "geometry is COMMON in the searched database. The bar is "
                      "the median of the WHOLE bin, both classes together, so "
                      "the red points can be read against it. Nothing here is "
                      "evidence of antigenicity."),
        )
        written.append(out.name)
    return written


def draw_annotated_distribution(stats: dict[str, Any], fig_dir: Path,
                                *, stem: str = "annotated_distribution",
                                ) -> list[str]:
    """Mean, variation and distinct counts per category, as a figure.

    Two panels, because "category" means two different splits that are both
    worth seeing: the two seed CLASSES (annotated against background), and the
    annotated residues split by whether they beat their own background.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    from ibc_search.visualization import THEMES

    panels = [
        ("By seed class", stats["by_seed_class"]),
        ("Annotated residues, by verdict", stats["annotated_by_verdict"]),
    ]
    if not any(any(v for v in d.values()) for _, d in panels):
        return []

    written = []
    for theme in ("light", "dark"):
        p = THEMES[theme]
        fig, axes = plt.subplots(1, 2, figsize=(12.0, 6.0))
        fig.patch.set_facecolor(p["surface"])
        for ax, (name, cats) in zip(axes, panels):
            ax.set_facecolor(p["surface"])
            for side in ("top", "right", "left"):
                ax.spines[side].set_visible(False)
            ax.spines["bottom"].set_color(p["grid"])
            ax.tick_params(colors=p["text_secondary"], length=0)
            ax.yaxis.grid(True, color=p["grid"], linewidth=0.6, alpha=0.7)
            ax.set_axisbelow(True)

            items = [(k, v) for k, v in cats.items() if v]
            if not items:
                ax.set_title(f"{name}\n(no data)", color=p["text_primary"])
                continue
            xs = np.arange(len(items))
            means = [v["mean"] for _, v in items]
            sds = [v["sd"] for _, v in items]
            # Exact keys, not a substring test: "annotated" contains the
            # letters "not", which quietly made an earlier version colour the
            # annotated bar as if it were the "does_not_beat" category.
            RED_KEYS = {SEED_CLASS_ANNOTATED, "beats_background"}
            colours = [p["summary"] if k in RED_KEYS else p["series"]
                       for k, _ in items]
            ax.bar(xs, means, width=0.56, color=colours, alpha=0.85,
                   linewidth=0)
            # One standard deviation, in the ink the run figures use for
            # spread, labelled SD rather than left ambiguous between SD, SEM
            # and IQR. The lower arm is clamped at zero: these are counts, so
            # a whisker crossing zero would draw an impossible value and is
            # a sign the distribution is skewed, which the caption says.
            lower = [min(sd, m) for m, sd in zip(means, sds)]
            ax.errorbar(xs, means, yerr=[lower, sds], fmt="none",
                        ecolor=p["variance"], elinewidth=1.2, capsize=5,
                        capthick=1.2, zorder=3)
            for x, (k, v) in zip(xs, items):
                ax.annotate(
                    f"mean {v['mean']:.1f}\nSD {v['sd']:.1f}\n"
                    f"n {v['n']}\n{v['n_distinct_values']} distinct values",
                    (x, v["mean"] + v["sd"]),
                    textcoords="offset points", xytext=(0, 8),
                    ha="center", va="bottom", fontsize=8,
                    color=p["text_secondary"], annotation_clip=False)
            ax.set_xlim(-0.72, len(items) - 0.28)
            ax.set_xticks(xs)
            ax.set_xticklabels([k.replace("_", " ") for k, _ in items],
                               color=p["text_secondary"], fontsize=9)
            ax.set_ylabel("distinct matched sites per ball",
                          color=p["text_secondary"], fontsize=9)
            ax.set_title(name, color=p["text_primary"], fontsize=11,
                         loc="left")
            top = max(v["mean"] + v["sd"] for _, v in items)
            ax.set_ylim(0, top * 1.38 if top else 1)

        frac = stats.get("pooled_fraction_beating")
        rate = "n/a" if frac is None else f"{frac:.1%}"
        fig.suptitle(
            "Distinct matched sites per ball, by category",
            color=THEMES[theme]["text_primary"], fontsize=13, x=0.012,
            ha="left", y=0.985)
        fig.text(
            0.012, 0.925,
            f"bars are means, whiskers +/- 1 SD; {rate} of annotated residues "
            f"beat the median of their own structure's background",
            color=THEMES[theme]["text_secondary"], fontsize=9, ha="left")
        fig.text(
            0.012, 0.012,
            "A COUNT, not a z-score: a high value means the ball's geometry is "
            "COMMON in the searched database, which is not evidence of\n"
            "antigenicity. 'Distinct values' is how many different counts "
            "occur in the category, not the number of balls. A whisker is "
            "clamped at zero where\none SD reaches below it, which means the "
            "category is skewed rather than that it holds negative counts.",
            color=THEMES[theme]["text_secondary"], fontsize=8, ha="left",
            va="bottom", linespacing=1.5)
        fig.subplots_adjust(top=0.84, bottom=0.20, left=0.07, right=0.985,
                            wspace=0.22)
        out = fig_dir / (f"{stem}.png" if theme == "light"
                         else f"{stem}_dark.png")
        fig.savefig(out, dpi=200, facecolor=p["surface"])
        plt.close(fig)
        written.append(out.name)
    return written


# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", type=Path, help="a few_balls_search run")
    ap.add_argument("--panel", type=Path, default=None,
                    help="member list [<run_dir>/panel.tsv]")
    ap.add_argument("--reference", choices=("median", "mean"),
                    default="median",
                    help="the background reference point an annotated ball "
                         "must beat [median]")
    args = ap.parse_args(argv)

    fig_dir = args.run_dir / "00_figures"
    stat_dir = args.run_dir / "01_statistics"
    fig_dir.mkdir(parents=True, exist_ok=True)
    stat_dir.mkdir(parents=True, exist_ok=True)

    records, missing = read_run(args.run_dir, args.panel)
    if not records:
        print("error: no balls found; is this a few_balls_search run?",
              file=sys.stderr)
        return 2
    if missing:
        print(f"  {len(missing)} panel members produced no patches and were "
              f"skipped (e.g. {', '.join(missing[:5])})")

    verdicts = annotated_beats_background(records, reference=args.reference)
    stats = statistics(records, verdicts)
    stats["n_panel_members_without_output"] = len(missing)
    stats["background_reference"] = args.reference

    out_tsv = stat_dir / "annotated_vs_background.tsv"
    with out_tsv.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(VERDICT_HEADER)
        for v in verdicts:
            w.writerow(["" if v[c] is None else v[c] for c in VERDICT_HEADER])
    (stat_dir / "annotated_vs_background.json").write_text(
        json.dumps(stats, indent=2))

    figs = draw_annotated_bins(records, fig_dir)
    figs += draw_annotated_distribution(stats, fig_dir)

    print("statistics")
    ann = stats["by_seed_class"].get(SEED_CLASS_ANNOTATED)
    bg = stats["by_seed_class"].get(SEED_CLASS_BACKGROUND)
    if ann and bg:
        print(f"  annotated balls : n={ann['n']:5d}  mean={ann['mean']:8.2f} "
              f"SD={ann['sd']:7.2f}  distinct values={ann['n_distinct_values']}")
        print(f"  background balls: n={bg['n']:5d}  mean={bg['mean']:8.2f} "
              f"SD={bg['sd']:7.2f}  distinct values={bg['n_distinct_values']}")
    print(f"  annotated residues with a verdict: "
          f"{stats['pooled_annotated_with_a_verdict']} of "
          f"{stats['pooled_annotated_residues']}")
    print(f"  beating their own background: {stats['pooled_n_beating']} "
          f"({stats['pooled_fraction_beating']})")
    print(f"  mean of per-structure fractions: "
          f"{stats['mean_of_per_structure_fractions']} "
          f"(SD {stats['sd_of_per_structure_fractions']}, "
          f"n={stats['n_structures_with_a_verdict']} structures)")
    print(f"\n  {out_tsv}")
    print(f"  {stat_dir / 'annotated_vs_background.json'}")
    print("figures")
    for f in figs:
        print(f"  {fig_dir / f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
