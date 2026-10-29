"""Place named query motifs against the reference population the z-score is fitted on.

The motifs handed to FoldDisco deliberately -- a zinc site, say -- are drawn
against the motif population whose mean is z = 0 and whose spread is the
z-score's unit. That population is the only honest comparison available.

Why NOT the structure's own motifs
----------------------------------
An earlier version of this figure compared each query motif against the other
motifs the surface scan extracted from the same structure, and reported a
percentile within it. That percentile measures nothing. The scan's other
motifs are not a negative class: any of them may be a genuinely rare,
genuinely interesting motif, and in a structure that happens to contain
several the query motif is pushed down the ranking for no reason connected to
the query. It also makes the number depend on how many motifs the scan
happened to emit, which is a property of the protein's size and surface, not
of the motif being asked about.

The reference population has neither problem. It is fixed, external, shared
by every structure, and it is what the z-score already means -- so a
percentile against it is just the z-score restated on a scale people read
more easily, and the two cannot disagree.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import json
import statistics as st
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ibc_search.visualization import binned_point_plot  # noqa: E402


def scores(path: Path) -> list[float]:
    if not path.is_file():
        return []
    out = []
    for r in csv.DictReader(path.open(), delimiter="\t"):
        try:
            out.append(float(r["idf_z"]))
        except (ValueError, KeyError):
            continue
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("query_run", type=Path, help="run of the named motifs")
    ap.add_argument("reference_run", type=Path,
                    help="the run the z-score parameters were fitted on")
    ap.add_argument("-o", "--out", type=Path, required=True)
    ap.add_argument("--label", default="zinc motif")
    ap.add_argument("--reference-label", default="reference population")
    ap.add_argument("--threshold", type=float, default=1.6)
    ap.add_argument("--theme", default="light", choices=("light", "dark"))
    args = ap.parse_args(argv)

    ref = sorted(z for d in args.reference_run.iterdir() if d.is_dir()
                 for z in scores(d / "idf_scored.tsv"))
    if not ref:
        print(f"error: no scored motifs under {args.reference_run}",
              file=sys.stderr)
        return 2

    def pct(z: float) -> float:
        return 100 * bisect.bisect_left(ref, z) / len(ref)

    rows = []
    for member in sorted(p for p in args.query_run.iterdir() if p.is_dir()):
        for z in scores(member / "idf_scored.tsv"):
            rows.append({"instance": member.name, "z": round(z, 4),
                         "percentile_vs_reference": round(pct(z), 2)})
    if not rows:
        print("error: no query motifs scored", file=sys.stderr)
        return 2
    rows.sort(key=lambda r: -r["z"])

    zs = [r["z"] for r in rows]
    above = sum(1 for z in zs if z >= args.threshold)
    binned_point_plot(
        [(args.reference_label, [(f"{len(ref):,} motifs", ref)]),
         (f"annotated {args.label}s", [(r["instance"], [r["z"]]) for r in rows])],
        args.out,
        y_label="motif IDF z-score",
        title=f"Annotated {args.label}s against the reference population",
        subtitle=(f"median z {st.median(zs):+.3f}, at the "
                  f"{st.median([r['percentile_vs_reference'] for r in rows]):.0f}th "
                  f"percentile of the reference; {above}/{len(rows)} at or "
                  f"above +{args.threshold:.1f}"),
        footnote=("Each query motif is exactly the residues the PDB entry "
                  "lists as its zinc site. Nothing was extracted by a scan;"
                  "\nself-hits are excluded. The reference is the motif set "
                  "the z-score was fitted on, so its mean is z = 0"
                  "\nand its spread is the y-axis unit. It replaces each "
                  "structure's own motifs, which are not a negative class"
                  "\nand give a percentile that measures nothing."),
        reference=args.threshold,
        reference_label=f"+{args.threshold:.1f} call threshold",
        theme=args.theme,
        bin_width=0.46, height=5.6, max_points=12000,
    )

    out_json = args.out.with_suffix(".json")
    out_json.write_text(json.dumps({
        "summary": {
            "query_motifs": len(rows),
            "reference_motifs": len(ref),
            "reference_median_z": round(st.median(ref), 4),
            "reference_p90": round(ref[int(0.90 * len(ref))], 4),
            "reference_p99": round(ref[int(0.99 * len(ref))], 4),
            "threshold": args.threshold,
            "at_or_above_threshold": above,
            "median_query_z": round(st.median(zs), 4),
            "median_percentile_vs_reference": round(
                st.median([r["percentile_vs_reference"] for r in rows]), 1),
        },
        "per_motif": rows,
    }, indent=2))

    print(f"wrote {args.out}\nwrote {out_json}")
    print(f"  {len(rows)} {args.label}s, median z {st.median(zs):+.3f}, "
          f"median percentile vs reference "
          f"{st.median([r['percentile_vs_reference'] for r in rows]):.1f}%")
    print(f"  {above} of {len(rows)} at or above +{args.threshold:.1f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
