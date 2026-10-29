"""Summary figures for a residue-scored run.

Which figures are drawn depends on whether the run has observed epitopes,
and the script decides that from the data rather than from a flag.

**With** epitope labels (an antibody-antigen complex), two figures:

``epitope_vs_other_scores.png``
    Do epitope residues score higher than the rest? Both classes are drawn
    side by side within each structure's bin, every point plotted -- the
    shape of the two clouds is the finding, so nothing is subsampled and
    nothing is bundled into a box.

``per_structure_auc.png``
    Is the ranking better than chance, structure by structure? One point per
    structure against the 0.5 line, which makes a split sample look like a
    split sample rather than a mean near chance.

**Without** labels (a designed binder predicted as a monomer), the same two
figures would be lies -- there is no positive class, so no AUC and no
"epitope vs other" exists. Three different figures are drawn instead, none
of which claims a correctness it cannot measure:

``residue_scores_by_structure.png``
    the per-residue score distribution, one bin per structure.
``calls_per_structure.png``
    how many residues each structure gets called at the threshold.
``elisa_vs_calls.png``
    those call counts against measured ADA reactivity, where the panel
    carries it. This is the only figure in the unlabelled set that can
    speak to whether the score tracks anything real, and it is drawn only
    when the measurement is present.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ibc_search.visualization import (  # noqa: E402
    binned_point_plot, paired_point_plot, scatter_plot, value_strip_plot,
)


def contact_distance(run_dir: Path) -> str:
    """How the labels were defined, read off the epitope file's name.

    A distance-defined label is named for its cutoff (``epitopes_4.5A.tsv``)
    and the figure can say "within 4.5A of". A label that is not a distance
    -- the PDB's own listed zinc sites, say -- has no cutoff to quote, and
    returning one would be a caption that lies about the label.
    """
    for path in sorted(run_dir.glob("epitopes_*.tsv")):
        tail = path.stem.split("_")[-1]
        if re.fullmatch(r"\d+(\.\d+)?A", tail):
            return tail
        return ""
    return ""


def unlabelled_figures(args, summary, by_inst, short, out_dir: Path) -> int:
    """Figures for a run with no observed epitope anywhere in it."""
    s = summary["summary"]
    per = {p["instance"]: p for p in summary["per_structure"]}

    # Ordered by how many residues each structure gets called, so the figure
    # answers "which binders carry the signal" by position, not by hunting.
    order = sorted(by_inst, key=lambda i: -per[i]["n_called"])
    binned_point_plot(
        [("", [(short(i), by_inst[i][1]) for i in order])],
        out_dir / "residue_scores_by_structure.png",
        y_label="max motif IDF z-score",
        title="Per-residue motif rarity across the designed binders",
        subtitle=(f"{s['structures']} binders, {s['residues_scored']:,} scored "
                  f"residues, {s['residues_called']} at or above "
                  f"+{s['threshold']:.1f} ({100 * s['call_rate']:.1f}%)"),
        footnote=("A residue's score is the highest z of any motif containing "
                  "it, against a natural-protein reference. Residues in no "
                  "emitted motif are absent, not zero. There is no antibody "
                  "in these structures: nothing here is a measured epitope. "
                  "Ordered by number of residues called."),
        reference=args.threshold,
        reference_label=f"+{s['threshold']:.1f} call threshold",
        max_points=10_000, theme=args.theme,
        bin_width=0.30, height=6.2,
    )

    # The title states the count rather than asserting a direction, so it
    # cannot drift out of agreement with the data the way a hardcoded claim
    # ("most binders carry none") silently does when the run changes.
    value_strip_plot(
        [(short(i), per[i]["n_called"]) for i in order],
        out_dir / "calls_per_structure.png",
        reference=0, reference_label="zero calls",
        y_label=f"residues at or above +{s['threshold']:.1f}",
        title=f"Residues called at +{s['threshold']:.1f}, per binder",
        subtitle=(f"median {s['median_calls_per_structure']:.0f}, "
                  f"max {s['max_calls_in_one_structure']}."),
        theme=args.theme,
    )

    # The ELISA column is the only external measurement available. Plotting
    # the score against it is the one honest test this dataset supports.
    elisa = {}
    panel = args.run_dir / "panel.tsv"
    if panel.is_file():
        for r in csv.DictReader(panel.open(), delimiter="\t"):
            try:
                elisa[r["instance"]] = float(r["elisa_mean"])
            except (KeyError, ValueError):
                continue
    points = [(short(i), elisa[i], float(per[i]["n_called"]))
              for i in order if i in elisa]
    if len(points) > 2:
        scatter_plot(
            points, out_dir / "elisa_vs_calls.png",
            x_label="ADA ELISA absorbance (blank-subtracted mean)",
            y_label=f"residues called at +{s['threshold']:.1f}",
            title="Called residues against measured ADA reactivity",
            subtitle="",
            footnote=("ELISA absorbance is pre-existing anti-drug-antibody "
                      "binding measured on the purified construct;\nthe call "
                      "count is this pipeline's structural prediction. "
                      "Neither was used to fit the other.\nBlank wells have "
                      "sd 0.131, so only 6 of 74 constructs clear 3 sd over "
                      "blank: this assay\ncannot validate or refute a "
                      "predictor, and the correlation is uninformative "
                      "either way."),
            annotate=5, theme=args.theme,
        )
        print(f"wrote {out_dir}/elisa_vs_calls.png")

    print(f"wrote {out_dir}/residue_scores_by_structure.png")
    print(f"wrote {out_dir}/calls_per_structure.png")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="[<run_dir>/figures]")
    ap.add_argument("--threshold", type=float, default=1.7)
    ap.add_argument("--theme", default="light", choices=("light", "dark"))
    ap.add_argument("--label", default="epitope",
                    help="what is_epitope==1 marks in this run; appears on "
                         "the figures [epitope]")
    ap.add_argument("--reference-feature", default="the antibody",
                    help="what the label distance is measured to "
                         "[the antibody]")
    args = ap.parse_args(argv)

    out_dir = args.out_dir or (args.run_dir / "00_figures")
    dist = contact_distance(args.run_dir)

    scores = args.run_dir / "residue_scores.tsv"
    summary_path = args.run_dir / "residue_score_summary.json"
    if not (scores.is_file() and summary_path.is_file()):
        print(f"error: run {args.run_dir} has no residue scores", file=sys.stderr)
        return 2
    summary = json.loads(summary_path.read_text())

    by_inst: dict[str, tuple[list[float], list[float]]] = defaultdict(
        lambda: ([], [])
    )
    for r in csv.DictReader(scores.open(), delimiter="\t"):
        z = float(r["max_motif_z"])
        by_inst[r["instance"]][0 if r["is_epitope"] == "1" else 1].append(z)

    # Short labels: the shared pdb_000 prefix, and the designed binders'
    # well/index/sample decorations, are noise on an axis.
    def short(inst: str) -> str:
        if "__" in inst:
            parts = inst.split("__")
            name = re.sub(r"^\d+_", "", parts[1]) if len(parts) > 1 else inst
            return f"{parts[0]} {name}"[:30]
        return inst.replace("pdb_0000", "").replace("_+", "")

    if not summary["summary"].get("ground_truth_epitopes", True):
        return unlabelled_figures(args, summary, by_inst, short, out_dir)

    # Ordered by how much the epitope cloud sits above the rest, so the
    # structures where the score did something are adjacent and comparable.
    def lift(kv):
        epi, oth = kv[1]
        if not epi or not oth:
            return 0.0
        return sum(epi) / len(epi) - sum(oth) / len(oth)

    bins = [(short(i), e, o)
            for i, (e, o) in sorted(by_inst.items(), key=lift, reverse=True)]

    s = summary["summary"]
    # The title states what was measured rather than asserting a direction.
    # A hardcoded "X do not score higher" was correct on the antibody panel
    # and flatly false on the zinc panel, where they do.
    lift_overall = s["epitope_mean_z"] - s["other_mean_z"]
    verdict = ("score higher than the rest" if lift_overall > 0.05
               else "score lower than the rest" if lift_overall < -0.05
               else "do not separate from the rest")
    paired_point_plot(
        bins, out_dir / "epitope_vs_other_scores.png",
        labels=(f"{args.label}{f' (<{dist})' if dist else ''}", "other"),
        y_label="max motif IDF z-score",
        title=f"{args.label.capitalize()} residues {verdict}",
        subtitle=(f"{s['structures']} structures, {s['residues_scored']:,} "
                  f"scored residues, {s['epitope_residues']} {args.label}.  "
                  f"{args.label} mean z {s['epitope_mean_z']:+.3f} vs "
                  f"other {s['other_mean_z']:+.3f}  "
                  f"(AUC {s['auc_mean']:.3f})"),
        footnote=("Every scored residue is drawn; a residue's score is the "
                  "highest z of any motif containing it.\nResidues in no "
                  "emitted motif have no score and are absent, not zero. "
                  f"Dashed line: the +{args.threshold:.1f} call threshold."
                  f"\nStructures ordered by {args.label}-minus-other mean."),
        reference=args.threshold, theme=args.theme,
    )

    value_strip_plot(
        [(short(p["instance"]), p["auc"]) for p in summary["per_structure"]],
        out_dir / "per_structure_auc.png",
        reference=0.5, reference_label="chance",
        y_label=f"AUC, {args.label} vs other",
        title=f"Per-structure ranking of {args.label} residues",
        subtitle=((f"within {dist} of {args.reference_feature}; "
                   if dist else f"{args.reference_feature}; ")
                  + f"mean AUC {s['auc_mean']:.3f}."),
        theme=args.theme,
    )

    print(f"wrote {out_dir}/epitope_vs_other_scores.png")
    print(f"wrote {out_dir}/per_structure_auc.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
