"""Does sub-epitope geometry transfer? A discrimination test with ground truth.

The question
------------
Earlier measurements asked "of the reference residues a motif matched, how
many are epitope?" on an *unlabelled* query protein. That has no ground truth
on the query side, so nothing can be right or wrong about it.

This asks the discriminative question instead. For an antigen **A** whose
epitope is known from its own complex:

    do A's EPITOPE-derived motifs retrieve epitope regions in OTHER antigens
    more often than A's NON-EPITOPE-derived motifs do?

Both classes come from the same protein, so antigen size, fold, composition
and every per-protein confound cancel. The statistic is a difference between
classes, not an absolute rate against an assumed null.

Leakage control
---------------
A hit to A itself, to another deposition of A (same PDB ID), or to anything
sharing A's ``ag_cluster`` is discarded. Without that, A's epitope motifs
retrieve near-identical copies of A and find their own epitope -- a
guaranteed, meaningless positive. Exclusion is done here rather than with
``--exclude-id`` so a whole cluster can be dropped at once; queries therefore
run with a generous ``--top`` and are filtered afterwards.

Class definition
----------------
A motif is *epitope-derived* when at least ``--pos-coverage`` of its residues
lie in A's epitope, and *non-epitope-derived* when none do. Motifs that
straddle the boundary are excluded from the comparison, to contrast clean
positives against clean negatives rather than blur them.

Only matched residues that are **exposed** in the hit antigen are counted: a
buried residue cannot be an antibody epitope, so such a match is
uninformative and would only dilute the rate.
"""

from __future__ import annotations

import argparse
import csv
import random
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def _clusters(value: str) -> set[str]:
    return {c for c in (value or "").split("/") if c}


def scan_and_search(
    antigen: Path, index: Path, out_dir: Path, radius: float, top: int
) -> tuple[Path, Path] | None:
    """Run ibc-search scan then a folddisco query. Returns (members, hits)."""
    env = {**dict(__import__("os").environ), "PYTHONPATH": str(REPO / "src")}
    scan = subprocess.run(
        [str(REPO / ".venv/bin/python"), "-m", "ibc_search", "scan", str(antigen),
         "--radius", str(radius), "-o", str(out_dir), "--log-level", "ERROR"],
        capture_output=True, text=True, env=env,
    )
    members = out_dir / "patch_members.tsv"
    queries = out_dir / "folddisco_queries.tsv"
    if scan.returncode != 0 or not queries.is_file() or queries.stat().st_size == 0:
        print(f"      scan rc={scan.returncode} queries_exists={queries.is_file()} "
              f"stderr={scan.stderr.strip()[-200:]!r}", flush=True)
        return None

    hits = out_dir / "hits.tsv"
    with hits.open("w") as fh:
        q = subprocess.run(
            ["folddisco", "query", "-i", str(index), "-q", str(queries),
             "-t", "1", "--top", str(top), "--header", "--per-match",
             "--format-output", "tid,idf,rmsd,query_residues,matching_residues"],
            stdout=fh, stderr=subprocess.PIPE, text=True,
        )
    if q.returncode != 0:
        print(f"      folddisco query rc={q.returncode}: "
              f"{(q.stderr or '').strip()[-400:]!r}", flush=True)
        return None
    return members, hits


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("index_dir", type=Path)
    ap.add_argument("-o", "--out", type=Path, required=True, help="results TSV")
    ap.add_argument("-n", "--n-antigens", type=int, default=60)
    ap.add_argument("--radius", type=float, default=6.0)
    ap.add_argument("--top", type=int, default=50)
    ap.add_argument("--pos-coverage", type=float, default=0.5)
    ap.add_argument("--max-residues", type=int, default=400,
                    help="skip very large antigens to keep runtime sane [400]")
    ap.add_argument("--split", default="test", help="which ab_ag_split to sample")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(argv)

    root: Path = args.index_dir
    index = root / "index" / "sabdab2_antigen"
    structures = {r["instance"]: r for r in
                  csv.DictReader((root / "structures.tsv").open(), delimiter="\t")}

    pool = [
        i for i, r in structures.items()
        if r["ab_ag_split"] == args.split
        and 0 < int(r["n_antigen_residues"]) <= args.max_residues
        and int(r["n_epitope_residues"]) > 0
    ]
    random.Random(args.seed).shuffle(pool)
    chosen = pool[: args.n_antigens]
    print(f"{len(pool)} candidate antigens in split={args.split}; using {len(chosen)}",
          flush=True)

    # ---- phase 1: scan + search each antigen ---------------------------
    per_antigen = []
    with tempfile.TemporaryDirectory() as tmp:
        for n, inst in enumerate(chosen, 1):
            out_dir = Path(tmp) / inst
            res = scan_and_search(
                root / structures[inst]["antigen_file"], index, out_dir,
                args.radius, args.top,
            )
            if res is None:
                print(f"  [{n}/{len(chosen)}] {inst}: scan/search failed", flush=True)
                continue
            members, hits = res

            # motif -> its residues, from the scan
            motif_res: dict[str, set[tuple[str, int]]] = defaultdict(set)
            for r in csv.DictReader(members.open(), delimiter="\t"):
                motif_res[r["patch_id"]].add((r["auth_chain_id"], int(r["auth_seq_id"])))

            # query string -> motif residues, to join hits back
            qmap: dict[str, set[tuple[str, int]]] = {}
            for r in csv.DictReader(
                (out_dir / "folddisco_queries.index.tsv").open(), delimiter="\t"
            ):
                qmap[r["query_string"]] = motif_res.get(r["patch_id"], set())

            rows = []
            hdr = None
            for line in hits.open():
                f = line.rstrip("\n").split("\t")
                if f[0] == "tid":
                    hdr = f
                    continue
                if hdr is None or len(f) != len(hdr):
                    continue
                rows.append((f[0], f[3], f[4]))
            per_antigen.append((inst, qmap, rows))
            if n % 5 == 0:
                print(f"  [{n}/{len(chosen)}] done", flush=True)

    # ---- phase 2: labels, once ------------------------------------------
    print("\nloading reference labels...", flush=True)
    needed = {inst for inst, _, _ in per_antigen}
    for _, _, rows in per_antigen:
        for tid, _, _ in rows:
            needed.add(tid.split("/")[-1].removesuffix(".cif"))

    epi: dict[tuple[str, str, int], int] = {}
    for r in csv.DictReader((root / "epitope_residues.tsv").open(), delimiter="\t"):
        if r["instance"] in needed:
            epi[(r["instance"], r["chain"], int(r["res_id"]))] = int(r["is_epitope"])
    exposed: dict[tuple[str, str, int], int] = {}
    for r in csv.DictReader((root / "residue_sasa.tsv").open(), delimiter="\t"):
        if r["instance"] in needed:
            exposed[(r["instance"], r["chain"], int(r["res_id"]))] = int(r["is_exposed"])
    print(f"  labels for {len(needed)} instances: {len(epi)} residues", flush=True)

    # ---- phase 3: per-motif transfer rate, by class ---------------------
    with args.out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["query_antigen", "query_string", "motif_size", "motif_class",
                    "epitope_coverage_of_motif", "n_hits_kept", "n_matched_exposed",
                    "n_matched_epitope", "epitope_transfer"])

        for inst, qmap, rows in per_antigen:
            meta = structures[inst]
            own_clusters = _clusters(meta["ag_cluster"])
            own_pdb = meta["PDB_ID"] if "PDB_ID" in meta else meta.get("pdb_id", "")

            agg: dict[str, list] = defaultdict(lambda: [0, 0, 0])
            for tid, qres, mres in rows:
                hinst = tid.split("/")[-1].removesuffix(".cif")
                hmeta = structures.get(hinst)
                if hmeta is None:
                    continue
                # leakage control
                if hinst == inst:
                    continue
                if own_pdb and hmeta.get("pdb_id", hmeta.get("PDB_ID", "")) == own_pdb:
                    continue
                if own_clusters & _clusters(hmeta["ag_cluster"]):
                    continue

                a = agg[qres]
                a[0] += 1
                for t in mres.split(","):
                    if not t or t == "_":
                        continue
                    try:
                        key = (hinst, t[0], int(t[1:]))
                    except ValueError:
                        continue
                    if not exposed.get(key, 0):
                        continue  # buried: cannot be an epitope
                    a[1] += 1
                    a[2] += epi.get(key, 0)

            for qres, residues in qmap.items():
                own = [epi.get((inst, c, r)) for c, r in residues]
                known = [v for v in own if v is not None]
                if not known:
                    continue
                cov = sum(known) / len(known)
                if cov >= args.pos_coverage:
                    cls = "epitope"
                elif cov == 0:
                    cls = "non_epitope"
                else:
                    cls = "straddle"
                n_hits, n_exp, n_epi = agg.get(qres, [0, 0, 0])
                w.writerow([
                    inst, qres, len(residues), cls, f"{cov:.3f}",
                    n_hits, n_exp, n_epi,
                    f"{n_epi / n_exp:.4f}" if n_exp else "",
                ])

    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
