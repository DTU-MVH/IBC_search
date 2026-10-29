"""Run-level statistics for an epitope-labelled panel, with honest intervals.

``residue_scores.py`` already writes a per-structure AUC and a mean over
structures. That mean is the wrong headline and this script exists to say so
with numbers: on a panel this size each structure carries few positives, its
AUC interval spans much of the range, and averaging forty noisy estimates
does not make a precise one. What is usable is the **pooled** figure with an
interval around it.

The four things computed here
-----------------------------
1. **Pooled AUC with a bootstrap interval resampled over STRUCTURES.**
   Residues within a structure are not independent -- they share a fold, a
   surface and a scan -- so resampling residues treats one protein's 200
   correlated residues as 200 independent facts and returns an interval that
   is far too narrow. The structure is the unit that was sampled, so it is
   the unit that is resampled. This is the number to quote.

2. **Per-structure AUC with its own interval,** reported precisely so the
   spread can be seen to be noise. Any ranking read off this column, or any
   story about which fold the method prefers, is reading sampling error.

3. **A threshold sweep** -- precision, recall, F1, MCC and enrichment over
   base rate at each cut -- because the call rate is extremely sensitive to
   where the threshold sits and a single reported number hides that.

4. **The coverage artifact (Defect 1).** A residue in *k* motifs takes the
   maximum of *k* draws, and the max of more draws is larger even when
   nothing differs. If epitope residues sit in more motifs than the rest,
   part of any separation is that counting and not rarity. The expected
   maximum E[max | k] is computed exactly from the run's own empirical motif
   distribution and subtracted, and the AUC is reported both ways.

What is deliberately not here
-----------------------------
No per-family or per-organism split. Partitioning a panel this size after
seeing the scores produces a significant-looking p-value from noise -- the
zinc panel's family split ranked first of all 210 possible partitions by gap
size, which is exactly what post-hoc selection manufactures. A split must be
declared before the scores are looked at, and none was.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path


def auc(pos: list[float], neg: list[float]) -> float | None:
    """Probability a positive outranks a negative; ties count half."""
    if not pos or not neg:
        return None
    merged = sorted([(v, 1) for v in pos] + [(v, 0) for v in neg])
    rank_sum = 0.0
    i = 0
    while i < len(merged):
        j = i
        while j < len(merged) and merged[j][0] == merged[i][0]:
            j += 1
        r = (i + j + 1) / 2
        rank_sum += r * sum(g for _, g in merged[i:j])
        i = j
    n_pos = len(pos)
    return (rank_sum - n_pos * (n_pos + 1) / 2) / (n_pos * len(neg))


def expected_max(sorted_pool: list[float], k: int) -> float:
    """E[max of k draws with replacement] from an empirical distribution.

    Exact rather than simulated: for sorted values z_(1)..z_(n), the maximum
    of k draws equals z_(i) with probability (i/n)^k - ((i-1)/n)^k.
    """
    n = len(sorted_pool)
    total = 0.0
    prev = 0.0
    for i, z in enumerate(sorted_pool, 1):
        p = (i / n) ** k
        total += z * (p - prev)
        prev = p
    return total


def mcc(tp: int, fp: int, fn: int, tn: int) -> float:
    d = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return ((tp * tn - fp * fn) / d) if d else 0.0


def load_rows(run_dir: Path) -> list[dict]:
    path = run_dir / "residue_scores.tsv"
    if not path.is_file():
        print(f"error: {path} missing; run residue_scores.py first",
              file=sys.stderr)
        raise SystemExit(2)
    rows = []
    for r in csv.DictReader(path.open(), delimiter="\t"):
        if r["is_epitope"] == "":
            continue  # unlabelled residue: no ground truth, excluded
        rows.append({
            "instance": r["instance"],
            "chain": r["chain"],
            "res_id": int(r["res_id"]),
            "is_epitope": int(r["is_epitope"]),
            "z": float(r["max_motif_z"]),
        })
    return rows


def load_coverage(run_dir: Path, instances: list[str]) -> dict:
    """(instance, chain, res_id) -> number of SCORED motifs containing it.

    Counts only motifs that carry a z, because an unscored motif contributes
    nothing to the maximum and must not inflate k.
    """
    cover: dict[tuple[str, str, int], int] = defaultdict(int)
    pool: list[float] = []
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
        patch_z = {}
        for r in csv.DictReader(index.open(), delimiter="\t"):
            z = by_query.get(r["query_string"])
            if z is not None:
                patch_z[r["patch_id"]] = z
        pool.extend(patch_z.values())
        for r in csv.DictReader(members.open(), delimiter="\t"):
            if r["patch_id"] in patch_z:
                cover[(inst, r["auth_chain_id"], int(r["auth_seq_id"]))] += 1
    return {"coverage": dict(cover), "pool": pool}


def bootstrap_pooled_auc(
    by_inst: dict[str, tuple[list[float], list[float]]],
    n_boot: int,
    seed: int,
) -> tuple[float | None, float | None]:
    """Percentile interval, resampling STRUCTURES rather than residues."""
    rng = random.Random(seed)
    names = list(by_inst)
    vals = []
    for _ in range(n_boot):
        pos: list[float] = []
        neg: list[float] = []
        for name in (rng.choice(names) for _ in names):
            p, n = by_inst[name]
            pos.extend(p)
            neg.extend(n)
        a = auc(pos, neg)
        if a is not None:
            vals.append(a)
    if len(vals) < 2:
        return None, None
    vals.sort()
    lo = vals[int(0.025 * len(vals))]
    hi = vals[min(int(0.975 * len(vals)), len(vals) - 1)]
    return lo, hi


def permutation_null(
    by_inst: dict[str, tuple[list[float], list[float]]],
    n_perm: int,
    seed: int,
) -> dict:
    """Label-permutation null for pooled AUC, permuting WITHIN each structure.

    Permuting within structures preserves each protein's own score
    distribution, so the test asks whether the epitope label picks out high
    scorers inside a protein -- not whether proteins differ from each other.

    **The null is not centred on 0.5, and that is the point.** Structures
    differ in both their score distributions and their epitope fractions, and
    pooling residues across them conflates that composition with real
    within-structure discrimination. The permuted mean measures exactly how
    much pooled AUC the composition alone buys, so it -- not 0.5 -- is what
    the observed value has to beat. Reporting a pooled AUC against an implied
    0.5 credits the score with the composition's work.
    """
    rng = random.Random(seed)
    pos = [z for p, _ in by_inst.values() for z in p]
    neg = [z for _, n in by_inst.values() for z in n]
    observed = auc(pos, neg)
    if observed is None:
        return {"p": float("nan"), "null_mean": None, "null_sd": None}
    vals = []
    for _ in range(n_perm):
        ppos: list[float] = []
        pneg: list[float] = []
        for p, n in by_inst.values():
            allz = p + n
            rng.shuffle(allz)
            ppos.extend(allz[: len(p)])
            pneg.extend(allz[len(p):])
        a = auc(ppos, pneg)
        if a is not None:
            vals.append(a)
    if not vals:
        return {"p": float("nan"), "null_mean": None, "null_sd": None}
    worse = sum(1 for a in vals if a >= observed)
    return {
        "p": (worse + 1) / (len(vals) + 1),
        "null_mean": round(st.mean(vals), 4),
        "null_sd": round(st.pstdev(vals), 4),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="[<run_dir>/01_statistics]")
    ap.add_argument("--threshold", type=float, default=1.6,
                    help="headline call threshold [1.6]")
    ap.add_argument("--sweep", type=float, nargs="+",
                    default=[1.0, 1.3, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0])
    ap.add_argument("--bootstrap", type=int, default=2000)
    ap.add_argument("--permutations", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(argv)

    out_dir = args.out_dir or (args.run_dir / "01_statistics")
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = load_rows(args.run_dir)
    if not rows:
        print("error: no labelled residues scored", file=sys.stderr)
        return 2

    by_inst: dict[str, tuple[list[float], list[float]]] = defaultdict(
        lambda: ([], []))
    for r in rows:
        by_inst[r["instance"]][0 if r["is_epitope"] else 1].append(r["z"])
    by_inst = {k: v for k, v in by_inst.items() if v[0] and v[1]}

    pos = [z for p, _ in by_inst.values() for z in p]
    neg = [z for _, n in by_inst.values() for z in n]
    pooled = auc(pos, neg)
    lo, hi = bootstrap_pooled_auc(by_inst, args.bootstrap, args.seed)
    perm = permutation_null(by_inst, args.permutations, args.seed)

    # ---- per-structure, with intervals that show the spread is noise ----
    rng = random.Random(args.seed)
    per_structure = []
    for inst, (p, n) in sorted(by_inst.items()):
        a = auc(p, n)
        boots = []
        for _ in range(500):
            bp = [rng.choice(p) for _ in p]
            bn = [rng.choice(n) for _ in n]
            b = auc(bp, bn)
            if b is not None:
                boots.append(b)
        boots.sort()
        per_structure.append({
            "instance": inst,
            "n_scored": len(p) + len(n),
            "n_epitope": len(p),
            "auc": round(a, 4),
            "ci_lo": round(boots[int(0.025 * len(boots))], 3) if boots else None,
            "ci_hi": round(boots[min(int(0.975 * len(boots)), len(boots) - 1)], 3)
                     if boots else None,
            "epitope_mean_z": round(st.mean(p), 4),
            "other_mean_z": round(st.mean(n), 4),
        })
    per_structure.sort(key=lambda d: d["auc"], reverse=True)

    with (out_dir / "per_structure_auc.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["instance", "n_scored", "n_epitope", "auc",
                    "ci_lo", "ci_hi", "epitope_mean_z", "other_mean_z"])
        for d in per_structure:
            w.writerow([d[k] for k in ("instance", "n_scored", "n_epitope",
                                       "auc", "ci_lo", "ci_hi",
                                       "epitope_mean_z", "other_mean_z")])

    # ---- threshold sweep ----
    base = len(pos) / (len(pos) + len(neg))
    sweep = []
    for t in args.sweep:
        tp = sum(1 for z in pos if z >= t)
        fn = len(pos) - tp
        fp = sum(1 for z in neg if z >= t)
        tn = len(neg) - fp
        prec = tp / (tp + fp) if tp + fp else float("nan")
        rec = tp / (tp + fn) if tp + fn else float("nan")
        f1 = (2 * prec * rec / (prec + rec)) if prec + rec else 0.0
        sweep.append({
            "threshold": t, "called": tp + fp, "tp": tp, "fp": fp,
            "precision": round(prec, 4) if prec == prec else None,
            "recall": round(rec, 4),
            "f1": round(f1, 4),
            "mcc": round(mcc(tp, fp, fn, tn), 4),
            "enrichment": round(prec / base, 3) if base and prec == prec else None,
        })
    with (out_dir / "threshold_sweep.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        cols = ["threshold", "called", "tp", "fp", "precision", "recall",
                "f1", "mcc", "enrichment"]
        w.writerow(cols)
        for d in sweep:
            w.writerow([d[c] for c in cols])

    # ---- Defect 1: the coverage artifact ----
    cov = load_coverage(args.run_dir, list(by_inst))
    coverage, cpool = cov["coverage"], sorted(cov["pool"])
    cov_stats = None
    if coverage and cpool:
        kp = [coverage.get((r["instance"], r["chain"], r["res_id"]), 0)
              for r in rows if r["is_epitope"]]
        kn = [coverage.get((r["instance"], r["chain"], r["res_id"]), 0)
              for r in rows if not r["is_epitope"]]
        emax = {k: expected_max(cpool, k) for k in range(1, max(kp + kn or [1]) + 1)}
        cpos, cneg, kpos, kneg = [], [], [], []
        for r in rows:
            k = coverage.get((r["instance"], r["chain"], r["res_id"]), 0)
            if k < 1:
                continue
            corrected = r["z"] - emax[k]
            (cpos if r["is_epitope"] else cneg).append(corrected)
            (kpos if r["is_epitope"] else kneg).append(float(k))
        cov_stats = {
            "epitope_mean_coverage": round(st.mean(kp), 3) if kp else None,
            "other_mean_coverage": round(st.mean(kn), 3) if kn else None,
            "auc_raw": round(pooled, 4) if pooled else None,
            "auc_coverage_alone": round(auc(kpos, kneg) or float("nan"), 4),
            "auc_corrected": round(auc(cpos, cneg) or float("nan"), 4),
        }
        with (out_dir / "coverage_artifact.tsv").open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["k_motifs", "n_epitope", "n_other", "expected_max_null"])
            for k in sorted(emax):
                w.writerow([k, kp.count(k), kn.count(k), round(emax[k], 4)])

    summary = {
        "run": str(args.run_dir),
        "structures": len(by_inst),
        "residues_scored": len(pos) + len(neg),
        "epitope_residues": len(pos),
        "epitope_base_rate": round(base, 4),
        "pooled_auc": round(pooled, 4) if pooled else None,
        "pooled_auc_ci95": [round(lo, 4), round(hi, 4)] if lo is not None else None,
        "pooled_auc_null_mean": perm["null_mean"],
        "pooled_auc_null_sd": perm["null_sd"],
        "pooled_auc_permutation_p": round(perm["p"], 4),
        "pooled_auc_excess_over_null": (
            round(pooled - perm["null_mean"], 4)
            if pooled is not None and perm["null_mean"] is not None else None),
        "bootstrap_unit": "structure",
        "auc_mean_over_structures": round(
            st.mean([d["auc"] for d in per_structure]), 4),
        "structures_above_chance": sum(1 for d in per_structure if d["auc"] > 0.5),
        "coverage_artifact": cov_stats,
        "headline_threshold": args.threshold,
    }
    (out_dir / "summary.json").write_text(json.dumps(
        {"summary": summary, "threshold_sweep": sweep,
         "per_structure": per_structure}, indent=2))

    print(f"wrote {out_dir}/")
    for k, v in summary.items():
        print(f"  {k:<32} {v}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
