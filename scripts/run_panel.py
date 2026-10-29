"""Scan and FoldDisco-search every structure in a panel, in parallel.

Reads a ``panel.tsv`` naming one structure per row and runs, for each,
``ibc-search scan`` followed by ``folddisco_search.sh search``, writing into
``<run_dir>/<instance>/``. Each member's console output is kept beside its
directory as ``<instance>.log`` so a failure is diagnosable afterwards
without re-running it.

Parallelism is across structures, one FoldDisco process each, with each
process single-threaded. That is deliberate: the installed FoldDisco build
does not write stdout line-atomically, so at ``-t > 1`` a fraction of hit
rows come back torn (see ``docs/FOLDDISCO_CONTRACT.md``). Running N
independent single-threaded searches gets the cores busy without putting two
writers on one stream.

Already-complete members are skipped, so an interrupted run resumes by being
re-invoked. ``--force`` re-does them.
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SEARCH = REPO / "scripts" / "folddisco_search.sh"

STRUCTURE_COLUMNS = ("structure_file", "antigen_file", "path")
"""Columns naming a structure FILE, tried in order."""

ID_COLUMN = "pdb_id"
"""Fallback: a PDB entry ID, which ``scan`` fetches (and caches) itself."""


def _run_one(task: tuple) -> tuple[str, bool, str]:
    inst, structure, out_dir, radius, log_path, top, extra = task
    env = dict(os.environ, PYTHONPATH=str(REPO / "src"))
    with open(log_path, "w") as log:
        scan = subprocess.run(
            [str(REPO / ".venv/bin/python"), "-m", "ibc_search", "scan",
             structure, "--radius", str(radius), "-o", out_dir],
            stdout=log, stderr=subprocess.STDOUT, cwd=str(REPO), env=env,
        )
        if scan.returncode != 0:
            return inst, False, "scan failed"
        queries = Path(out_dir, "folddisco_queries.tsv")
        if not queries.is_file() or not queries.stat().st_size:
            return inst, False, "scan emitted no expressible queries"
        cmd = [str(SEARCH), "search", out_dir, "--top", str(top), *extra]
        search = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT,
                                cwd=str(REPO))
        if search.returncode != 0:
            return inst, False, "search failed"
    return inst, True, ""


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_dir", type=Path, help="directory holding panel.tsv")
    ap.add_argument("--radius", type=float, required=True,
                    help="ball radius in Angstroms for the scan")
    ap.add_argument("--top", type=int, default=40,
                    help="FoldDisco --top, structures kept per query [40]")
    ap.add_argument("-j", "--jobs", type=int, default=6)
    ap.add_argument("--force", action="store_true",
                    help="re-run members that already have hits")
    ap.add_argument("--search-arg", action="append", default=[],
                    help="extra argument for folddisco_search.sh; repeatable")
    args = ap.parse_args(argv)

    panel = args.run_dir / "panel.tsv"
    if not panel.is_file():
        print(f"error: no panel at {panel}", file=sys.stderr)
        return 2

    rows = list(csv.DictReader(panel.open(), delimiter="\t"))
    col = next((c for c in STRUCTURE_COLUMNS if rows and c in rows[0]), None)
    has_id = bool(rows) and ID_COLUMN in rows[0]
    if col is None and not has_id:
        print(f"error: {panel} has no structure column (looked for "
              f"{', '.join(STRUCTURE_COLUMNS)}, {ID_COLUMN})", file=sys.stderr)
        return 2

    tasks, skipped = [], 0
    for row in rows:
        # A path wins where the panel gives one -- a corrected or
        # chain-selected copy is there precisely to be used in place of the
        # deposited entry. Otherwise the ID goes to `scan`, which resolves it
        # from its cache or RCSB.
        structure = (row.get(col) or "").strip() if col else ""
        by_id = not structure
        if by_id:
            structure = (row.get(ID_COLUMN) or "").strip()
        inst = row.get("instance") or structure
        if not structure:
            print(f"  {inst}: no structure named, skipped")
            continue
        if not by_id and not Path(structure).is_file():
            print(f"  {inst}: structure missing, skipped")
            continue
        out_dir = args.run_dir / inst
        # "Already done" means a hits file with a header and at least one
        # row, not merely a file that exists. A run killed mid-search leaves
        # an empty or header-only file behind, and testing only for presence
        # silently adopts it as a finished member -- which then fails much
        # later, in whatever consumes the hits, or worse does not fail.
        hits = out_dir / "folddisco_hits.tsv"
        if not args.force and hits.is_file():
            with hits.open() as fh:
                complete = sum(1 for _, _ in zip(fh, range(2))) >= 2
            if complete:
                skipped += 1
                continue
            print(f"  {inst}: hits file is empty or truncated, re-running")
        out_dir.mkdir(parents=True, exist_ok=True)
        tasks.append((inst, structure, str(out_dir), args.radius,
                      str(args.run_dir / f"{inst}.log"), args.top,
                      tuple(args.search_arg)))

    if skipped:
        print(f"{skipped} members already searched, skipped (--force re-runs)")
    if not tasks:
        print("nothing to do")
        return 0
    print(f"scanning and searching {len(tasks)} structures at radius "
          f"{args.radius} A, {args.jobs} at a time", flush=True)

    started = time.time()
    ok = bad = 0
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(_run_one, t) for t in tasks]
        for n, fut in enumerate(as_completed(futures), 1):
            inst, good, detail = fut.result()
            if good:
                ok += 1
            else:
                bad += 1
                print(f"  {inst}: FAILED -- {detail} "
                      f"(see {args.run_dir / f'{inst}.log'})")
            if n % 10 == 0 or n == len(tasks):
                print(f"  {n}/{len(tasks)} done, {time.time() - started:.0f}s",
                      flush=True)

    print(f"\n{ok} searched, {bad} failed, {time.time() - started:.0f}s elapsed")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
