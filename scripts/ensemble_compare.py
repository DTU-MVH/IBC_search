"""EXPERIMENTAL: quantify how much a perturbation ensemble changes results.

Companion to :mod:`ibc_search.ensemble`. It answers one question -- how much
of a single-structure result survives perturbation within the deposited
positional uncertainty -- at two stages:

1. **Motif extraction.** Which patches are even found? Comparable across
   members because perturbation moves coordinates but never renumbers
   residues, so a query string like ``A27,A28,A49`` denotes the same residues
   in every member.
2. **FoldDisco retrieval.** Which structures come back as hits, and is the
   top hit for a given motif stable?

Stability is reported as the number of members a thing appears in, which is
more informative than a mean: a motif found in 10/10 members and one found
in 3/10 can share an average yet mean very different things.

Usage::

    python scripts/ensemble_compare.py ensembles/8yhu
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path


def read_queries(run_dir: Path) -> set[str]:
    """The motif token strings a run extracted."""
    path = run_dir / "folddisco_queries.tsv"
    if not path.is_file():
        return set()
    out = set()
    for line in path.read_text().splitlines():
        if "\t" in line:
            out.add(line.split("\t", 1)[1])
    return out


def read_hits(run_dir: Path) -> tuple[dict[str, list[tuple[str, float, float]]], set[str]]:
    """Return ``{query_string: [(tid, idf, rmsd), ...]}`` and the set of tids."""
    path = run_dir / "folddisco_hits.tsv"
    by_query: dict[str, list[tuple[str, float, float]]] = {}
    tids: set[str] = set()
    if not path.is_file():
        return by_query, tids
    with path.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
        try:
            i_tid = header.index("tid")
            i_idf = header.index("idf")
            i_rmsd = header.index("rmsd")
            i_q = header.index("query_residues")
        except ValueError:
            return by_query, tids
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if len(f) != len(header):
                continue
            try:
                row = (f[i_tid], float(f[i_idf]), float(f[i_rmsd]))
            except ValueError:
                continue
            by_query.setdefault(f[i_q], []).append(row)
            tids.add(f[i_tid])
    return by_query, tids


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b)


def pct(n: int, d: int) -> str:
    return f"{100.0 * n / d:5.1f}%" if d else "    - "


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("ensemble_dir", type=Path)
    args = ap.parse_args(argv)

    root: Path = args.ensemble_dir
    baseline_dir = root / "run_baseline"
    control_dir = root / "run_00"
    member_dirs = sorted(d for d in root.glob("run_[0-9][0-9]") if d != control_dir)

    if not baseline_dir.is_dir():
        print(f"error: no {baseline_dir}", file=sys.stderr)
        return 2

    base_q = read_queries(baseline_dir)
    base_hits, base_tids = read_hits(baseline_dir)

    print("=" * 72)
    print("ENSEMBLE VARIANCE  (B-factor uncertainty ensemble, NOT conformational)")
    print("=" * 72)

    # ---- control -------------------------------------------------------
    if control_dir.is_dir():
        ctrl_q = read_queries(control_dir)
        _, ctrl_tids = read_hits(control_dir)
        same_q = ctrl_q == base_q
        same_t = ctrl_tids == base_tids
        print("\nCONTROL (unperturbed, same CIF round-trip as the members)")
        print(f"  identical motifs to baseline : {same_q}")
        print(f"  identical hit structures     : {same_t}")
        if not (same_q and same_t):
            print("  !! the round-trip itself changes results; variance below is NOT")
            print("     attributable to perturbation alone")

    # ---- stage 1: motif extraction -------------------------------------
    print("\n" + "-" * 72)
    print("STAGE 1 - MOTIF EXTRACTION")
    print("-" * 72)
    print(f"  baseline motifs: {len(base_q)}")
    print(f"\n  {'member':<10}{'motifs':>8}{'shared w/ base':>16}{'recovered':>11}{'new':>8}")
    member_qs = []
    for d in member_dirs:
        q = read_queries(d)
        member_qs.append(q)
        shared = len(q & base_q)
        print(
            f"  {d.name:<10}{len(q):>8}{shared:>16}"
            f"{pct(shared, len(base_q)):>11}{len(q - base_q):>8}"
        )

    if member_qs:
        counts = Counter()
        for q in member_qs:
            counts.update(q)
        n = len(member_qs)
        all_seen = set(counts)
        core = {k for k, v in counts.items() if v == n}
        once = {k for k, v in counts.items() if v == 1}
        print(f"\n  union over {n} members      : {len(all_seen)}")
        print(f"  in ALL {n} members (core)   : {len(core)}  ({pct(len(core), len(all_seen))} of union)")
        print(f"  in exactly 1 member        : {len(once)}  ({pct(len(once), len(all_seen))} of union)")
        base_core = len(base_q & core)
        print(f"  baseline motifs in core    : {base_core}  ({pct(base_core, len(base_q))} of baseline)")
        print("\n  motif occupancy spectrum (how many members each motif appears in):")
        spec = Counter(counts.values())
        for k in range(1, n + 1):
            bar = "#" * int(40 * spec.get(k, 0) / max(spec.values()))
            print(f"    {k:>2}/{n} members: {spec.get(k,0):>5}  {bar}")

    # ---- stage 1b: soft overlap ----------------------------------------
    # Exact motif identity is brittle: a 7-residue patch that gains or loses
    # one residue scores as completely different. These two measures
    # distinguish "the surface is being sampled differently" from "the same
    # surface, shifted by a residue".
    if member_qs:
        print("\n" + "-" * 72)
        print("STAGE 1b - SOFT OVERLAP (exact identity is brittle)")
        print("-" * 72)

        base_sets = [frozenset(q.split(",")) for q in base_q]
        print(f"\n  {'member':<10}{'mean best Jaccard':>19}{'>=0.8':>8}{'>=0.5':>8}")
        for d, q in zip(member_dirs, member_qs):
            msets = [frozenset(x.split(",")) for x in q]
            by_res: dict[str, list[frozenset]] = {}
            for s in msets:
                for r in s:
                    by_res.setdefault(r, []).append(s)
            best_scores = []
            for b in base_sets:
                cands = {id(s): s for r in b for s in by_res.get(r, [])}
                best_scores.append(
                    max((len(b & s) / len(b | s) for s in cands.values()), default=0.0)
                )
            hi = sum(1 for s in best_scores if s >= 0.8)
            mid = sum(1 for s in best_scores if s >= 0.5)
            mean = sum(best_scores) / len(best_scores) if best_scores else 0.0
            print(f"  {d.name:<10}{mean:>19.3f}{pct(hi,len(base_sets)):>8}{pct(mid,len(base_sets)):>8}")
        print("\n  (for each baseline motif, the most similar motif in that member)")

        # Residue-level: is the same surface being sampled at all?
        base_res = {r for q in base_q for r in q.split(",")}
        rc = Counter()
        for q in member_qs:
            rc.update({r for s in q for r in s.split(",")})
        n = len(member_qs)
        res_all = {r for r, v in rc.items() if v == n}
        print(f"\n  residues covered by baseline motifs : {len(base_res)}")
        print(f"  residues covered in ALL {n} members   : {len(res_all)}"
              f"  ({pct(len(res_all), len(base_res))} of baseline's)")
        print(f"  baseline residues never covered      : "
              f"{len(base_res - set(rc))}")

    # ---- stage 2: retrieval --------------------------------------------
    if base_tids:
        print("\n" + "-" * 72)
        print("STAGE 2 - FOLDDISCO RETRIEVAL")
        print("-" * 72)
        print(f"  baseline distinct hit structures: {len(base_tids)}")
        print(f"\n  {'member':<10}{'hit structs':>12}{'shared':>9}{'Jaccard':>9}{'top1 agree':>12}")
        member_tids = []
        for d in member_dirs:
            hits, tids = read_hits(d)
            member_tids.append(tids)
            shared = len(tids & base_tids)

            # Top-1 agreement on motifs both runs extracted.
            agree = total = 0
            for q, rows in hits.items():
                if q not in base_hits:
                    continue
                total += 1
                best = min(rows, key=lambda r: r[2])[0]
                base_best = min(base_hits[q], key=lambda r: r[2])[0]
                agree += best == base_best
            print(
                f"  {d.name:<10}{len(tids):>12}{shared:>9}"
                f"{jaccard(tids, base_tids):>9.3f}{pct(agree, total):>12}"
            )

        if member_tids:
            pair = [
                jaccard(member_tids[i], member_tids[j])
                for i in range(len(member_tids))
                for j in range(i + 1, len(member_tids))
            ]
            print(f"\n  pairwise Jaccard between members: "
                  f"min {min(pair):.3f}  mean {sum(pair)/len(pair):.3f}  max {max(pair):.3f}")

            tcounts = Counter()
            for t in member_tids:
                tcounts.update(t)
            n = len(member_tids)
            union = len(tcounts)
            robust = sum(1 for v in tcounts.values() if v == n)
            single = sum(1 for v in tcounts.values() if v == 1)
            print(f"  union of hit structures     : {union}")
            print(f"  found by ALL {n} members     : {robust}  ({pct(robust, union)} of union)")
            print(f"  found by exactly 1 member   : {single}  ({pct(single, union)} of union)")
    else:
        print("\n(no baseline hits file yet - run the searches first)")

    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
