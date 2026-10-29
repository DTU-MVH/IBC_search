"""Check that a fitted z-score reference actually standardises what it claims.

Fitting parameters and using them are two different acts, and nothing in the
fit verifies itself. This re-reads a reference run through its own fitted
parameters and asks the three questions that decide whether the z-score means
what its name says:

1. **Is the size trend real and linear in log space?** IDF is a sum over
   shared encodings, so it grows with motif size. If the trend is not linear
   the straight-line detrend leaves curvature behind, and the z-score then
   ranks motifs partly by size.
2. **Do the residuals centre on zero at every size?** A stratum whose
   residual mean is not zero is systematically scored high or low for no
   reason but its size.
3. **Is the result actually standard?** If the reference's own z-scores do
   not have mean 0 and sd 1, then "z = +1.6" does not mean 1.6 standard
   deviations and no threshold expressed in z is interpretable.

A reference that fails 3 on its own data cannot be trusted on anything else,
so this is run on the fitting set deliberately: it is a self-consistency
check, not a generalisation estimate.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ibc_search.visualization import binned_point_plot  # noqa: E402


def collect(run_dir: Path) -> list[tuple[int, float, float]]:
    """(motif size, log10 IDF aggregate, z) over every scored motif."""
    out = []
    for member in sorted(p for p in run_dir.iterdir() if p.is_dir()):
        path = member / "idf_scored.tsv"
        if not path.is_file():
            continue
        for r in csv.DictReader(path.open(), delimiter="\t"):
            try:
                out.append((int(r["motif_size"]), float(r["log10_idf"]),
                            float(r["idf_z"])))
            except (ValueError, KeyError):
                continue
    return out


def normal_quantile(p: float) -> float:
    """Inverse standard normal, Acklam's rational approximation."""
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    pl, ph = 0.02425, 1 - 0.02425
    if p < pl:
        q = math.sqrt(-2 * math.log(p))
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    if p > ph:
        q = math.sqrt(-2 * math.log(1 - p))
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
                ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    q = p - 0.5
    r = q * q
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / \
           (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_dir", type=Path, help="the fitted reference run")
    ap.add_argument("params", type=Path, help="its fitted parameter file")
    ap.add_argument("--out-dir", type=Path, default=None)
    ap.add_argument("--theme", default="light", choices=("light", "dark"))
    args = ap.parse_args(argv)

    out_dir = args.out_dir or (args.run_dir / "00_figures")
    out_dir.mkdir(parents=True, exist_ok=True)
    params = json.loads(args.params.read_text())
    a, b = params["intercept_a"], params["slope_b"]
    centre = {int(k): v for k, v in params.get("centre_by_size", {}).items()}
    sd_by = {int(k): v for k, v in params.get("sd_by_size", {}).items()}
    pooled = params["sd_pooled"]

    rows = collect(args.run_dir)
    if not rows:
        print(f"error: no scored motifs under {args.run_dir}", file=sys.stderr)
        return 2

    by_size: dict[int, list[tuple[float, float]]] = defaultdict(list)
    for size, logidf, z in rows:
        by_size[size].append((logidf, z))
    sizes = sorted(by_size)

    table = []
    for s in sizes:
        logs = [x for x, _ in by_size[s]]
        zs = [z for _, z in by_size[s]]
        expected = a + b * s + centre.get(s, 0.0)
        table.append({
            "motif_size": s,
            "n_motifs": len(logs),
            "mean_log10_idf": round(st.mean(logs), 4),
            "fitted_expected": round(expected, 4),
            "residual_mean": round(st.mean(logs) - expected, 4),
            "sd_used": round(sd_by.get(s, pooled), 4),
            "observed_sd": round(st.pstdev(logs), 4) if len(logs) > 1 else None,
            "z_mean": round(st.mean(zs), 4),
            "z_sd": round(st.pstdev(zs), 4) if len(zs) > 1 else None,
            "own_sd_fitted": s in sd_by,
        })

    allz = sorted(z for _, _, z in rows)
    n = len(allz)
    # Correlation of the sample quantiles with the normal ones: 1.0 would be
    # exactly normal. Reported instead of a test, because at n in the tens of
    # thousands every test rejects and the size of the departure is what
    # matters.
    nq = [normal_quantile((i + 0.5) / n) for i in range(n)]
    mz, mq = st.mean(allz), st.mean(nq)
    num = sum((x - mz) * (y - mq) for x, y in zip(allz, nq))
    den = math.sqrt(sum((x - mz) ** 2 for x in allz)
                    * sum((y - mq) ** 2 for y in nq))
    qq_r = num / den if den else float("nan")

    summary = {
        "run": str(args.run_dir), "params": str(args.params),
        "structures": sum(1 for p in args.run_dir.iterdir()
                          if (p / "idf_scored.tsv").is_file()),
        "motifs": n,
        "fit": {"intercept_a": a, "slope_b": b, "sd_pooled": pooled,
                "size_range": params.get("size_range"),
                "variance_explained_by_size":
                    params.get("variance_explained_by_size")},
        "z_mean": round(st.mean(allz), 4),
        "z_sd": round(st.pstdev(allz), 4),
        "z_median": round(st.median(allz), 4),
        "z_quantiles": {q: round(allz[min(n - 1, int(q / 100 * n))], 4)
                        for q in (1, 5, 25, 50, 75, 90, 95, 99, 99.9)},
        "fraction_above_1.6": round(sum(1 for z in allz if z >= 1.6) / n, 4),
        "normal_fraction_above_1.6": 0.0548,
        "qq_correlation_vs_normal": round(qq_r, 5),
        "worst_residual_mean_by_size": max(
            (abs(r["residual_mean"]), r["motif_size"]) for r in table),
        "worst_z_mean_by_size": max(
            (abs(r["z_mean"]), r["motif_size"]) for r in table),
    }

    tsv = out_dir / "calibration_by_motif_size.tsv"
    with tsv.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(list(table[0]))
        for r in table:
            w.writerow([r[k] for k in table[0]])
    (out_dir / "calibration_summary.json").write_text(
        json.dumps({"summary": summary, "by_motif_size": table}, indent=2))

    keep = [s for s in sizes if len(by_size[s]) >= 20]
    lab = [f"{s}  (n={len(by_size[s]):,})" for s in keep]

    binned_point_plot(
        [("", [(l, [x for x, _ in by_size[s]]) for l, s in zip(lab, keep)])],
        out_dir / "idf_by_motif_size.png",
        y_label="log10 of the per-motif IDF aggregate",
        title="IDF grows with motif size, and the fit removes it",
        subtitle=(f"{n:,} motifs from {summary['structures']} proteins.  "
                  f"fit: log10(IDF) = {a:.4f} + {b:.4f} x size, "
                  f"size explains "
                  f"{100 * (params.get('variance_explained_by_size') or 0):.0f}% "
                  f"of the variance"),
        footnote=("Dashed red: the fitted expectation per size, including the "
                  "per-size centre correction where a stratum was populated "
                  "enough to fit one.\nSizes with fewer than 20 motifs are "
                  "omitted from the figure; they appear in the table and fall "
                  "back to the pooled sd."),
        reference=None, theme=args.theme, summary="median",
        max_points=4000, bin_width=0.78, height=5.2,
        trend=[a + b * s + centre.get(s, 0.0) for s in keep],
        trend_label="fitted expectation",
    )

    binned_point_plot(
        [("", [(l, [x - (a + b * s + centre.get(s, 0.0))
                    for x, _ in by_size[s]]) for l, s in zip(lab, keep)])],
        out_dir / "size_trend_residuals.png",
        y_label="residual of log10 IDF after detrending",
        title="Residuals after removing the size trend",
        subtitle=("each stratum should centre on zero; the largest departure "
                  f"is {summary['worst_residual_mean_by_size'][0]:.3f} at size "
                  f"{summary['worst_residual_mean_by_size'][1]}"),
        footnote=("A stratum whose residuals do not centre on zero is scored "
                  "systematically high or low for no reason but its motif "
                  "size.\nDashed line: zero."),
        reference=0.0, reference_label="zero", theme=args.theme,
        summary="median", max_points=4000, bin_width=0.78, height=5.2,
    )

    binned_point_plot(
        [("", [(l, [z for _, z in by_size[s]]) for l, s in zip(lab, keep)])],
        out_dir / "z_by_motif_size.png",
        y_label="IDF z-score",
        title="The z-score after calibration, by motif size",
        subtitle=(f"overall mean {summary['z_mean']:+.4f}, sd "
                  f"{summary['z_sd']:.4f}; quantile correlation against a "
                  f"normal {qq_r:.4f}"),
        footnote=("If the reference's own z-scores are not mean 0 and sd 1, "
                  "then a threshold stated in z does not mean what it says."
                  "\nDashed line: zero."),
        reference=0.0, reference_label="zero", theme=args.theme,
        summary="median", max_points=4000, bin_width=0.78, height=5.2,
    )

    print(f"{summary['structures']} proteins, {n:,} motifs")
    print(f"  fit        log10(IDF) = {a:.4f} + {b:.4f} * size   "
          f"pooled sd {pooled:.4f}")
    print(f"  z          mean {summary['z_mean']:+.4f}  sd {summary['z_sd']:.4f}"
          f"  median {summary['z_median']:+.4f}")
    print(f"  normality  quantile correlation {qq_r:.4f}")
    print(f"  tail       {100 * summary['fraction_above_1.6']:.2f}% at or above "
          f"+1.6 (a normal would give 5.48%)")
    print(f"  worst per-size z mean {summary['worst_z_mean_by_size'][0]:.3f} "
          f"at size {summary['worst_z_mean_by_size'][1]}")
    print(f"\nwrote {tsv}")
    for name in ("idf_by_motif_size", "size_trend_residuals", "z_by_motif_size"):
        print(f"wrote {out_dir / (name + '.png')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
