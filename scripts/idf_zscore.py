"""Z-score normalise FoldDisco IDF, removing the motif-size dependence.

Why not a plain z-score on IDF
------------------------------
Raw IDF is heavily right-skewed -- measured on 883 motifs of 8YHU against the
whole-PDB index at ``--top 100``, mean 12613 against median 5288 (ratio 2.38)
with a standard deviation (22353) *larger than the mean*. Standardising that
directly lets a handful of extreme motifs set the scale and squashes
everything else into a narrow band. Taking log10 first makes the distribution
essentially symmetric (mean/median 1.03).

IDF also scales with motif size, because the score is a **sum** over the
encodings a candidate shares with the query and a larger motif contributes
more residue pairs. Across sizes 3-8 the mean rose 5.7x in raw terms, +0.69
decades in log terms. Left in, a z-score would mostly rank motifs by size.

The correction
--------------
In log space the size trend is linear and the detrending is clean::

    log10(IDF) = a + b * size            a = 3.169, b = 0.1333 on that run
    residual means by size: +0.03 -0.04 -0.00 +0.00 -0.02 +0.06

so::

    z = (log10(IDF_motif) - (a + b * size)) / s(size)

The denominator is size-dependent because residual spread grows with size
(sd 0.29 at size 3 to 0.47 at size 7); a single pooled sd would over-penalise
small motifs and under-penalise large ones.

The parameters are NOT universal
--------------------------------
IDF depends on the reference index (rarity is measured over *its* contents)
and on ``--top`` (which truncates the hit list the aggregate is taken from).
Parameters fitted against the whole-PDB index do not apply to the SAbDab2
antigen index, or to the same index queried with a different cap. ``fit``
records the index and top it was fitted on, and ``apply`` warns on a
mismatch rather than silently producing incomparable numbers.

Fit once on a reference run and reuse, rather than fitting per protein: a
self-fitted z-score centres every protein on zero, which is fine for ranking
motifs *within* one protein but makes two proteins incomparable -- and
comparing constructs is the point.
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
from ibc_search.config import load_config, resolve_database  # noqa: E402

MIN_STRATUM = 20
"""Motifs needed in a size stratum before its own sd is trusted."""


def per_motif(path: Path, aggregate: str, topk: int) -> dict[str, list[float]]:
    """Collect each motif's IDF values, keyed by its query string."""
    vals: dict[str, list[float]] = defaultdict(list)
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if "query_residues" not in (reader.fieldnames or []):
            raise SystemExit(f"{path}: no query_residues column")
        for r in reader:
            try:
                vals[r["query_residues"]].append(float(r["idf"]))
            except (ValueError, KeyError, TypeError):
                continue
    return vals


def aggregate_idf(values: list[float], how: str, topk: int) -> float | None:
    v = [x for x in values if x > 0]
    if not v:
        return None
    if how == "max":
        return max(v)
    if how == "median":
        return st.median(v)
    if how == "mean_topk":
        v.sort(reverse=True)
        return st.mean(v[:topk])
    raise SystemExit(f"unknown aggregate {how!r}")


def collect(paths: list[Path], how: str, topk: int) -> list[tuple[str, int, float]]:
    out = []
    for p in paths:
        for q, vals in per_motif(p, how, topk).items():
            a = aggregate_idf(vals, how, topk)
            if a is not None:
                out.append((q, len([t for t in q.split(",") if t]), a))
    return out


def cmd_fit(args) -> int:
    rows = collect(args.hits, args.aggregate, args.topk)
    if len(rows) < 50:
        raise SystemExit(f"only {len(rows)} motifs; too few to fit")
    pts = [(n, math.log10(a)) for _, n, a in rows]

    n = len(pts)
    sx = sum(x for x, _ in pts)
    sy = sum(y for _, y in pts)
    sxy = sum(x * y for x, y in pts)
    sxx = sum(x * x for x, _ in pts)
    denom = n * sxx - sx * sx
    if denom == 0:
        raise SystemExit("all motifs are the same size; cannot fit a size trend")
    b = (n * sxy - sx * sy) / denom
    a = (sy - b * sx) / n

    resid_by_size: dict[int, list[float]] = defaultdict(list)
    for x, y in pts:
        resid_by_size[x].append(y - (a + b * x))
    allr = [r for v in resid_by_size.values() for r in v]
    pooled = st.pstdev(allr)

    sd_by_size = {
        str(s): round(st.pstdev(v), 4)
        for s, v in sorted(resid_by_size.items())
        if len(v) >= MIN_STRATUM and st.pstdev(v) > 0
    }
    # The line is a good model of the size trend on a large index but not on
    # every one: against the 5951-structure antigen set the residual means
    # traced an inverted U (-0.33 at size 3, +0.21 at size 7, -0.70 at size
    # 10), because IDF saturates for large motifs when there are fewer
    # structures to be rare against. A per-size centre removes that curvature
    # where a stratum is well populated; sparse sizes still fall back to the
    # fitted line, which is why the line is kept.
    centre_by_size = {
        str(s): round(st.mean(v), 4)
        for s, v in sorted(resid_by_size.items())
        if len(v) >= MIN_STRATUM
    }

    # Min-max needs the observed range per size, in whichever space the
    # scaling is applied. Recorded for both so one fit serves either method.
    raw_by_size: dict[int, list[float]] = defaultdict(list)
    log_by_size: dict[int, list[float]] = defaultdict(list)
    for _, size_, agg_ in rows:
        raw_by_size[size_].append(agg_)
        log_by_size[size_].append(math.log10(agg_))
    range_by_size = {
        str(sz): {
            "raw_min": round(min(raw_by_size[sz]), 4),
            "raw_max": round(max(raw_by_size[sz]), 4),
            "log_min": round(min(log_by_size[sz]), 4),
            "log_max": round(max(log_by_size[sz]), 4),
            "n": len(raw_by_size[sz]),
        }
        for sz in sorted(raw_by_size)
        if len(raw_by_size[sz]) >= MIN_STRATUM
    }
    all_raw = [a for _, _, a in rows]
    range_global = {
        "raw_min": round(min(all_raw), 4), "raw_max": round(max(all_raw), 4),
        "log_min": round(math.log10(min(all_raw)), 4),
        "log_max": round(math.log10(max(all_raw)), 4),
    }

    params = {
        "aggregate": args.aggregate,
        "topk": args.topk if args.aggregate == "mean_topk" else None,
        "intercept_a": round(a, 4),
        "slope_b": round(b, 4),
        "sd_by_size": sd_by_size,
        "centre_by_size": centre_by_size,
        "range_by_size": range_by_size,
        "range_global": range_global,
        "sd_pooled": round(pooled, 4),
        "min_stratum": MIN_STRATUM,
        "n_motifs": n,
        "size_range": [min(x for x, _ in pts), max(x for x, _ in pts)],
        "variance_explained_by_size": round(
            1 - st.pvariance(allr) / st.pvariance([y for _, y in pts]), 4
        ),
        "fitted_on_index": args.index_label,
        "fitted_on_top": args.top_label,
        "fitted_on_radius_angstrom": args.radius_label,
        "source_files": [str(p) for p in args.hits],
    }
    args.out.write_text(json.dumps(params, indent=2))
    print(f"fitted on {n} motifs: log10(IDF) = {a:.4f} + {b:.4f} * size")
    print(f"  pooled residual sd {pooled:.4f}; per-size sd for "
          f"{len(sd_by_size)} strata")
    print(f"  size explains {100 * params['variance_explained_by_size']:.0f}% of variance")
    print(f"wrote {args.out}")
    return 0


def _params_path(args) -> Path:
    """Explicit --params wins; otherwise resolve from the database registry."""
    if args.params:
        return args.params
    entry = resolve_database(load_config(), args.database)
    path = entry.get("idf_params_path")
    if path is None:
        raise SystemExit(f"database {entry['name']!r} has no idf_params configured")
    if not path.is_file():
        raise SystemExit(
            f"database {entry['name']!r} points at {path}, which does not exist. "
            "Fit it first with `idf_zscore.py fit`."
        )
    print(f"database: {entry['name']}  params: {path.name}", file=sys.stderr)
    return path


def cmd_apply(args) -> int:
    params = json.loads(_params_path(args).read_text())
    a, b = params["intercept_a"], params["slope_b"]
    sd_by_size = {int(k): v for k, v in params["sd_by_size"].items()}
    centre_by_size = {int(k): v for k, v in params.get("centre_by_size", {}).items()}
    pooled = params["sd_pooled"]

    if args.index_label and params.get("fitted_on_index") not in (None, args.index_label):
        print(f"WARNING: parameters were fitted on index "
              f"{params.get('fitted_on_index')!r} but applying to "
              f"{args.index_label!r}. IDF rarity is defined over the index "
              f"contents, so these z-scores are not comparable.", file=sys.stderr)
    # `--top` ranks structures by IDF, so the MAXIMUM is retained at any cap:
    # verified identical for all 883 reference motifs at top 40 vs top 100.
    # Only cap-dependent aggregates need the warning.
    if (params["aggregate"] != "max"
            and args.top_label
            and params.get("fitted_on_top") not in (None, args.top_label)):
        print(f"WARNING: parameters were fitted at --top "
              f"{params.get('fitted_on_top')} but applying to {args.top_label}, "
              f"and aggregate={params['aggregate']!r} depends on where the hit "
              f"list is truncated.", file=sys.stderr)

    # Radius is checked as loudly as the index, because it was measured to
    # change the fit more than the protein set does: the size slope is 0.11
    # at 6 A and 0.057 at 5-5.5 A, so applying across radii reintroduces
    # exactly the size dependence the z-score exists to remove.
    fitted_radius = params.get("fitted_on_radius_angstrom")
    if (args.radius_label is not None and fitted_radius is not None
            and abs(float(fitted_radius) - args.radius_label) > 1e-9):
        print(f"WARNING: parameters were fitted at ball radius "
              f"{fitted_radius} A but applying to {args.radius_label} A. The "
              f"IDF size trend depends on radius, so these z-scores carry a "
              f"residual size dependence.", file=sys.stderr)

    rows = collect(args.hits, params["aggregate"], params["topk"] or 10)
    extrapolated = 0
    lo, hi = params["size_range"]
    with args.out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["query_residues", "motif_size", "idf_aggregate",
                    "log10_idf", "expected_log10", "sd_used", "idf_z"])
        for q, size, agg in sorted(rows, key=lambda r: -r[2]):
            y = math.log10(agg)
            # fitted line, corrected by the stratum's own residual mean where
            # one was measured
            exp = a + b * size + centre_by_size.get(size, 0.0)
            sd = sd_by_size.get(size, pooled)
            if not (lo <= size <= hi):
                extrapolated += 1
            w.writerow([q, size, f"{agg:.4f}", f"{y:.4f}", f"{exp:.4f}",
                        f"{sd:.4f}", f"{(y - exp) / sd:+.4f}"])
    print(f"wrote {len(rows)} motifs to {args.out}")
    if extrapolated:
        print(f"  note: {extrapolated} motifs fall outside the fitted size "
              f"range {lo}-{hi}; their trend value is extrapolated")
    return 0


def cmd_minmax(args) -> int:
    """Scale IDF strictly into [0, 1] per motif size.

    An alternative to the z-score that assumes nothing about the shape of the
    distribution -- no mean, no standard deviation, no normality. It is a
    linear rescaling, so within a size stratum it preserves the distribution
    exactly: only the axis changes.

    Scaling is per size stratum, because IDF scales multiplicatively with
    motif size (5.7x across sizes 3-8 on the PDB) and a single global range
    would simply rank motifs by size.

    ``--space raw`` preserves the raw distribution's shape, which is heavily
    right-skewed: on the PDB reference, mean 12613 against median 5288. Most
    motifs therefore land near 0 and only the extreme tail approaches 1,
    which is the intended behaviour if the tail is what you are hunting.
    ``--space log`` spreads the mid-range instead, at the cost of compressing
    the tail.

    Values outside the fitted range are clamped, and the count is reported:
    a new protein can legitimately exceed the reference maximum, and silently
    emitting 1.3 would break the [0, 1] guarantee the method exists to give.
    """
    params = json.loads(_params_path(args).read_text())
    ranges = {int(k): v for k, v in params.get("range_by_size", {}).items()}
    glob = params.get("range_global")
    if not ranges and not glob:
        raise SystemExit(
            "these parameters predate min-max support and carry no ranges; "
            "re-run `fit` to regenerate them"
        )

    lo_key, hi_key = f"{args.space}_min", f"{args.space}_max"
    rows = collect(args.hits, params["aggregate"], params["topk"] or 10)
    clamped = fellback = 0

    with args.out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["query_residues", "motif_size", "idf_aggregate",
                    "scaled_value", "range_min", "range_max", "idf_minmax"])
        for q, size, agg in sorted(rows, key=lambda r: -r[2]):
            entry = ranges.get(size)
            if entry is None:
                entry = glob
                fellback += 1
            lo, hi = entry[lo_key], entry[hi_key]
            value = agg if args.space == "raw" else math.log10(agg)
            if hi <= lo:
                scaled = 0.5          # degenerate stratum: no range to scale on
            else:
                scaled = (value - lo) / (hi - lo)
            if scaled < 0 or scaled > 1:
                clamped += 1
                scaled = min(1.0, max(0.0, scaled))
            w.writerow([q, size, f"{agg:.4f}", f"{value:.4f}",
                        f"{lo:.4f}", f"{hi:.4f}", f"{scaled:.4f}"])

    print(f"wrote {len(rows)} motifs to {args.out}  (space={args.space})")
    if fellback:
        print(f"  {fellback} motifs had no fitted stratum and used the global range")
    if clamped:
        print(f"  {clamped} motifs fell outside the reference range and were "
              f"clamped into [0, 1]")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    f = sub.add_parser("fit", help="fit size-trend parameters on reference hits")
    f.add_argument("hits", nargs="+", type=Path, help="folddisco_hits.tsv file(s)")
    f.add_argument("-o", "--out", type=Path, required=True)
    f.add_argument("--aggregate", choices=("max", "median", "mean_topk"), default="max")
    f.add_argument("--topk", type=int, default=10)
    f.add_argument("--index-label", default=None, help="which index these came from")
    f.add_argument("--top-label", type=int, default=None, help="--top used")
    f.add_argument("--radius-label", type=float, default=None,
                   help="ball radius the fitting runs were scanned at. "
                        "Recorded so `apply` can refuse to mix radii: the "
                        "size slope tracks the radius (0.11 at 6 A against "
                        "0.057 at 5-5.5 A), so a reference fitted at one "
                        "radius is a different population.")
    f.set_defaults(func=cmd_fit)

    p = sub.add_parser("apply", help="z-score motifs using fitted parameters")
    p.add_argument("hits", nargs="+", type=Path)
    p.add_argument("-p", "--params", type=Path, default=None,
                   help="explicit parameter file, overriding --database")
    p.add_argument("-d", "--database", default=None,
                   help="named database from [folddisco.databases]; its fitted "
                        "parameters are used (default: default_database)")
    p.add_argument("-o", "--out", type=Path, required=True)
    p.add_argument("--index-label", default=None)
    p.add_argument("--top-label", type=int, default=None)
    p.add_argument("--radius-label", type=float, default=None,
                   help="ball radius these hits were scanned at; warns on a "
                        "mismatch with the fitted radius")
    p.set_defaults(func=cmd_apply)

    m = sub.add_parser("minmax", help="scale IDF into [0,1] per size (no mean assumed)")
    m.add_argument("hits", nargs="+", type=Path)
    m.add_argument("-p", "--params", type=Path, default=None)
    m.add_argument("-d", "--database", default=None)
    m.add_argument("-o", "--out", type=Path, required=True)
    m.add_argument("--space", choices=("raw", "log"), default="raw",
                   help="scale the raw IDF (preserves its shape) or log10 "
                        "(spreads the mid-range) [raw]")
    m.set_defaults(func=cmd_minmax)

    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
