"""Do motifs seeded on ANNOTATED residues match more of the PDB than background motifs?

The question, stated carefully
------------------------------
"Do annotated residues get more hits than background?" cannot be answered by
comparing the two runs' per-residue counts, and this script does not try.
Two things differ between the runs that have nothing to do with the biology,
and both were MEASURED before this file was written:

**Motifs per residue.** In the background run every residue seeds a ball, so
a residue belongs to 7.8 motifs on average. In the seeded run only annotated
residues seed, so a residue belongs to 2.1. A per-residue count aggregates
over those motifs, so comparing per-residue counts would mostly measure how
many balls happened to contain the residue.

**Motif size.** Zinc motifs average 9.62 residues against the background's
7.88. A larger motif is more specific and so matches fewer places in the
database, which pushes the comparison the other way.

So the unit here is the **motif**, not the residue: one ball of one radius
against one index is directly comparable between runs. And the headline
comparison is **stratified by motif size**, because the size difference is
large enough to reverse a naive result.

What is compared
----------------
For every motif, the number of **distinct matched sites** its hits found --
the same ``(reference structure, matched residue set)`` deduplication
``count_scoring`` uses. Group A is the annotated-seeded motifs, group B the
background motifs.

Three views, because they can disagree and each is honest about something
different:

``pooled``
    every motif, both groups. Largest n, but motifs within a structure are
    not independent, so its interval is optimistic.
``per_structure``
    one median per structure, then the two sets of medians compared. Fewer
    points, but the unit of replication is right.
``by_motif_size``
    within each motif size separately. This is the one to read, because it
    is the only view in which the size difference cannot drive the result.

The statistic is a rank-based **AUC**: the probability that a randomly
chosen annotated motif matched more sites than a randomly chosen background
motif, with 0.5 meaning no difference. Rank-based because the count
distribution is right-tailed and a t-test's assumptions do not hold.

What a difference does NOT mean
-------------------------------
A higher count means the motif's geometry is **commoner in the PDB**, not
that it is more antigenic, more functional or more interesting. Zinc sites
are a structurally stereotyped, heavily-deposited motif, so finding them
common is close to a positive control for the method working, and is not a
discovery about zinc fingers. Nothing here establishes antigenicity.
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "src"))

from count_scoring import load_motif_sites  # noqa: E402


def search_top(run_dir: Path) -> int | None:
    """The ``--top`` retrieval cap the run used, from any run manifest."""
    for m in sorted(run_dir.glob("*/run_manifest.json"))[:1]:
        try:
            d = json.loads(m.read_text())
        except (json.JSONDecodeError, OSError):
            return None
        cfg = (d.get("resolved_config") or {}).get("folddisco") or {}
        top = cfg.get("search_top")
        return int(top) if top else None
    return None


def motif_site_counts(run_dir: Path) -> dict[str, list[tuple[int, int]]]:
    """``instance -> [(motif_size, distinct_sites)]`` for every emitted motif.

    A motif that was searched and matched nothing contributes ``0``, not an
    absent row: it was assessed. A motif that emitted no query is absent,
    because it never was.
    """
    out: dict[str, list[tuple[int, int]]] = defaultdict(list)
    panel = run_dir / "panel.tsv"
    if not panel.is_file():
        raise SystemExit(f"no panel.tsv in {run_dir}")

    for row in csv.DictReader(panel.open(), delimiter="\t"):
        inst = row["instance"]
        d = run_dir / inst
        patches = d / "patches.tsv"
        if not patches.is_file():
            continue
        sites = load_motif_sites(d)
        for r in csv.DictReader(patches.open(), delimiter="\t"):
            if r.get("emit_query", "true").lower() == "false":
                continue
            try:
                size = int(r["n_members"])
            except (KeyError, ValueError):
                continue
            out[inst].append((size, len(sites.get(r["patch_id"], ()))))
    return out


def auc(pos: list[float], neg: list[float]) -> float | None:
    """P(a random positive outranks a random negative); ties count half."""
    if not pos or not neg:
        return None
    merged = sorted([(v, 1) for v in pos] + [(v, 0) for v in neg])
    ranks: dict[int, float] = {}
    i = 0
    while i < len(merged):
        j = i
        while j < len(merged) and merged[j][0] == merged[i][0]:
            j += 1
        r = (i + j + 1) / 2
        for k in range(i, j):
            ranks[k] = r
        i = j
    rp = sum(ranks[k] for k, (_, g) in enumerate(merged) if g == 1)
    return (rp - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def describe(v: list[float]) -> dict:
    if not v:
        return {"n": 0}
    s = sorted(v)
    n = len(s)
    return {
        "n": n,
        "mean": round(st.fmean(s), 3),
        "sd": round(st.pstdev(s), 3) if n > 1 else 0.0,
        "min": s[0],
        "q25": s[n // 4],
        "median": round(st.median(s), 3),
        "q75": s[(3 * n) // 4],
        "max": s[-1],
    }


def bootstrap_auc_ci(pos: list[float], neg: list[float], *, n_boot: int,
                     seed: int) -> list[float] | None:
    """Percentile interval, resampling both groups.

    Reported because an AUC from a few hundred points carries a wide
    interval, and a bare point estimate invites reading a difference that
    the data does not support.
    """
    import random

    if not pos or not neg or n_boot <= 0:
        return None
    rng = random.Random(seed)
    vals = []
    for _ in range(n_boot):
        p = [pos[rng.randrange(len(pos))] for _ in range(len(pos))]
        q = [neg[rng.randrange(len(neg))] for _ in range(len(neg))]
        a = auc(p, q)
        if a is not None:
            vals.append(a)
    if not vals:
        return None
    vals.sort()
    return [round(vals[int(0.025 * len(vals))], 4),
            round(vals[min(int(0.975 * len(vals)), len(vals) - 1)], 4)]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("annotated_run", type=Path,
                    help="the seeded run (motifs centred on annotated residues)")
    ap.add_argument("background_run", type=Path)
    ap.add_argument("-o", "--out-dir", type=Path, default=None,
                    help="default: <annotated_run>/01_statistics")
    ap.add_argument("--bootstrap", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--min-per-size", type=int, default=30,
                    help="skip a size stratum with fewer than this in either "
                         "group; an AUC over a handful of motifs is noise [30]")
    args = ap.parse_args(argv)

    out_dir = args.out_dir or (args.annotated_run / "01_statistics")
    out_dir.mkdir(parents=True, exist_ok=True)

    ann = motif_site_counts(args.annotated_run)
    bg = motif_site_counts(args.background_run)
    if not ann or not bg:
        print("error: one of the runs has no emitted motifs", file=sys.stderr)
        return 2

    # --- censoring check, BEFORE any comparison --------------------------
    #
    # `--top N` truncates each query's hit list before the match stage, so a
    # motif commoner than N places in the PDB cannot report more than N
    # distinct sites. Its count is CENSORED, not measured. Where a large
    # share of one group is censored and the other is not, a difference
    # between the groups is an artefact of the cap and not a fact about the
    # motifs -- so the check runs first and the verdict is carried in the
    # output rather than left for a reader to infer.
    top = search_top(args.annotated_run) or search_top(args.background_run)
    def censored(pairs):
        if not top:
            return None
        n = sum(1 for _, c in pairs if c >= top)
        return {"cap": top, "n_at_or_above_cap": n, "n": len(pairs),
                "fraction": round(n / len(pairs), 4) if pairs else None}

    ann_all = [(s, c) for v in ann.values() for s, c in v]
    bg_all = [(s, c) for v in bg.values() for s, c in v]
    a_counts = [float(c) for _, c in ann_all]
    b_counts = [float(c) for _, c in bg_all]

    # --- pooled ----------------------------------------------------------
    pooled = {
        "annotated": describe(a_counts),
        "background": describe(b_counts),
        "auc": round(auc(a_counts, b_counts) or 0.0, 4),
        "auc_ci95": bootstrap_auc_ci(a_counts, b_counts,
                                     n_boot=min(args.bootstrap, 400),
                                     seed=args.seed),
        "caveat": ("motifs within a structure are not independent, so this "
                   "interval is optimistic; see per_structure"),
    }

    # --- per structure ---------------------------------------------------
    a_med = [st.median([c for _, c in v]) for v in ann.values() if v]
    b_med = [st.median([c for _, c in v]) for v in bg.values() if v]
    per_structure = {
        "annotated": describe(a_med),
        "background": describe(b_med),
        "auc": round(auc(a_med, b_med) or 0.0, 4),
        "auc_ci95": bootstrap_auc_ci(a_med, b_med, n_boot=args.bootstrap,
                                     seed=args.seed),
        "unit": "one median per structure, which is the honest unit of "
                "replication",
    }

    # --- by motif size ---------------------------------------------------
    a_by: dict[int, list[float]] = defaultdict(list)
    b_by: dict[int, list[float]] = defaultdict(list)
    for s, c in ann_all:
        a_by[s].append(float(c))
    for s, c in bg_all:
        b_by[s].append(float(c))

    strata = []
    for size in sorted(set(a_by) & set(b_by)):
        if len(a_by[size]) < args.min_per_size or \
                len(b_by[size]) < args.min_per_size:
            continue
        strata.append({
            "motif_size": size,
            "n_annotated": len(a_by[size]),
            "n_background": len(b_by[size]),
            "median_annotated": round(st.median(a_by[size]), 1),
            "median_background": round(st.median(b_by[size]), 1),
            "auc": round(auc(a_by[size], b_by[size]) or 0.0, 4),
        })

    cen_a, cen_b = censored(ann_all), censored(bg_all)
    heavily_censored = bool(
        cen_a and cen_b and max(cen_a["fraction"], cen_b["fraction"]) > 0.2
        and abs(cen_a["fraction"] - cen_b["fraction"]) > 0.1)

    payload = {
        "annotated_run": str(args.annotated_run),
        "background_run": str(args.background_run),
        "unit": "one MOTIF (one ball), not one residue",
        "value": "distinct matched sites per motif",
        "why_not_per_residue": (
            "the two runs differ in motifs per residue (background 7.8, "
            "seeded 2.1) because only annotated residues seed in the seeded "
            "run. A per-residue count aggregates over those motifs, so "
            "comparing per-residue counts would largely measure the seeding "
            "difference."),
        "confounder_motif_size": {
            "annotated_mean_size": round(
                st.fmean([s for s, _ in ann_all]), 2),
            "background_mean_size": round(
                st.fmean([s for s, _ in bg_all]), 2),
            "note": ("a larger motif is more specific and matches fewer "
                     "sites, so this difference alone can move the pooled "
                     "result. by_motif_size is the view in which it cannot."),
        },
        "censoring": {
            "annotated": cen_a,
            "background": cen_b,
            "comparison_valid": not heavily_censored,
            "note": (
                "--top N truncates each query's hit list before matching, so "
                "a motif commoner than N places cannot report more than N "
                "distinct sites: its count is CENSORED, not measured."),
            "verdict": (
                "INVALID as a measure of how common the motifs are. One "
                "group is largely at the cap and the other is not, so the "
                "AUC below mostly measures which group saturates. Re-run "
                "both sides with a higher --top, or with --no-top, before "
                "reading any magnitude from this."
                if heavily_censored else
                "censoring is low or similar in both groups; the comparison "
                "is readable."),
        },
        "pooled": pooled,
        "per_structure": per_structure,
        "by_motif_size": strata,
        "interpretation": (
            "AUC is P(a random annotated motif matched more distinct sites "
            "than a random background motif); 0.5 is no difference. A higher "
            "count means the geometry is COMMONER in the PDB -- not more "
            "antigenic, more functional or more interesting. Zinc sites are "
            "stereotyped and heavily deposited, so finding them common is "
            "closer to a positive control than a discovery."),
    }
    (out_dir / "annotated_vs_background.json").write_text(
        json.dumps(payload, indent=2))

    with (out_dir / "annotated_vs_background_by_size.tsv").open(
            "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["motif_size", "n_annotated", "n_background",
                    "median_annotated", "median_background", "auc"])
        for s in strata:
            w.writerow([s[k] for k in ("motif_size", "n_annotated",
                                       "n_background", "median_annotated",
                                       "median_background", "auc")])

    print(f"wrote {out_dir/'annotated_vs_background.json'}")
    if cen_a and cen_b:
        print(f"  CENSORING at --top {top}: annotated "
              f"{cen_a['fraction']:.1%} at the cap, background "
              f"{cen_b['fraction']:.1%}")
    if heavily_censored:
        print("  *** COMPARISON INVALID as a measure of motif commonness:")
        print("      one group saturates the retrieval cap and the other "
              "does not, so the AUC below mostly measures that.")
        print("      Re-run both sides with a higher --top before reading "
              "any magnitude from it.")
    print(f"  pooled        annotated median {pooled['annotated']['median']} "
          f"vs background {pooled['background']['median']}   "
          f"AUC {pooled['auc']} {pooled['auc_ci95']}")
    print(f"  per structure annotated median {per_structure['annotated']['median']} "
          f"vs background {per_structure['background']['median']}   "
          f"AUC {per_structure['auc']} {per_structure['auc_ci95']}")
    print(f"  size strata   {len(strata)} sizes with >= {args.min_per_size} "
          "motifs in both groups")
    for s in strata:
        print(f"     size {s['motif_size']:>3}: median "
              f"{s['median_annotated']:>7} vs {s['median_background']:>7}   "
              f"AUC {s['auc']}  (n {s['n_annotated']}/{s['n_background']})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
