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

sys.path.insert(0, str(REPO / "src"))

from ibc_search.config import load_config, resolve_database  # noqa: E402

STRUCTURE_COLUMNS = ("structure_file", "antigen_file", "path")
"""Columns naming a structure FILE, tried in order."""

ID_COLUMN = "pdb_id"
"""Fallback: a PDB entry ID, which ``scan`` fetches (and caches) itself."""


def _run_one(task: tuple) -> tuple[str, bool, str]:
    (inst, structure, out_dir, radius, log_path, top, extra, method,
     seed_residues, seed_rounds, min_idf, background_mode,
     min_background, background_rsa) = task
    env = dict(os.environ, PYTHONPATH=str(REPO / "src"))
    with open(log_path, "w") as log:
        scan = subprocess.run(
            [str(REPO / ".venv/bin/python"), "-m", "ibc_search", method,
             structure, "--radius", str(radius), "-o", out_dir]
            + (["--seed-residues", seed_residues,
                # The panel knows the entry identity; the filename may not
                # carry it (biological assemblies are `<id>-assembly1`).
                "--seed-structure-id", inst] if seed_residues else [])
            # --seed-rounds belongs to dumb-search only: few-balls-search
            # has no expansion rounds and would reject the flag.
            + (["--seed-rounds", str(seed_rounds)]
               if seed_rounds and method == "dumb-search" else [])
            + (["--background-mode", background_mode]
               if background_mode and method == "few-balls-search" else [])
            + (["--min-background", str(min_background)]
               if min_background is not None
               and method == "few-balls-search" else [])
            + (["--background-rsa", str(background_rsa)]
               if background_rsa is not None
               and method == "few-balls-search" else []),
            stdout=log, stderr=subprocess.STDOUT, cwd=str(REPO), env=env,
        )
        if scan.returncode != 0:
            return inst, False, "scan failed"
        queries = Path(out_dir, "folddisco_queries.tsv")
        if not queries.is_file() or not queries.stat().st_size:
            return inst, False, "scan emitted no expressible queries"
        cmd = [str(SEARCH), "search", out_dir, "--top", str(top), *extra]
        if min_idf is not None:
            cmd += ["--min-idf", str(min_idf)]
        search = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT,
                                cwd=str(REPO))
        if search.returncode != 0:
            return inst, False, "search failed"
    return inst, True, ""


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_dir", type=Path, help="directory holding panel.tsv")
    ap.add_argument("--panel", type=Path, default=None, metavar="TSV",
                    help="read the member list from here instead of "
                         "<run_dir>/panel.tsv. Results still land in "
                         "run_dir, so several processes can share one run "
                         "directory as long as their panels are DISJOINT -- "
                         "which is how worker capacity is added to a run "
                         "already in progress without two workers racing on "
                         "one structure's output files")
    ap.add_argument("--radius", type=float, required=True,
                    help="ball radius in Angstroms for the scan")
    ap.add_argument("--method",
                    choices=("scan", "dumb-search", "few-balls-search"),
                    default="scan",
                    help="extraction method. `scan` is the filtered ball "
                         "search; `dumb-search` places one ball per residue "
                         "with no thresholds and so emits several times more "
                         "queries per structure; `few-balls-search` places "
                         "one ball per annotated residue plus an equal, "
                         "spatially dispersed background set on the same "
                         "structure, which needs --seed-residues and emits "
                         "far FEWER queries than either [scan]")
    ap.add_argument("--min-background", type=int, default=None, metavar="N",
                    help="few-balls-search only: floor on the number of "
                         "background balls per structure. The default of 5 "
                         "is a statistical minimum, not a good control: a "
                         "structure with 5 background balls cannot supply "
                         "an exposure-matched comparison, which is what "
                         "left the reviewed-enzyme matched null resting on "
                         "4 observations. Raise it for sparsely annotated "
                         "sets")
    ap.add_argument("--background-rsa", type=float, default=None,
                    metavar="FRAC",
                    help="few-balls-search only: the buried/exposed boundary "
                         "used by --background-mode surface and by the "
                         "burial stratification")
    ap.add_argument("--background-mode", choices=("protein", "surface"),
                    default=None,
                    help="few-balls-search only: `surface` restricts the "
                         "background balls to solvent-exposed residues, "
                         "which is what a surface annotation (interface, "
                         "epitope) needs so the control does not differ in "
                         "burial as well as in being annotated")
    ap.add_argument("--seed-rounds", type=int, default=None, metavar="N",
                    help="rounds of seeding with --seed-residues; 2 also "
                         "seeds every residue inside a round-1 ball")
    ap.add_argument("--seed-residues", default=None, metavar="TSV",
                    help="dumb-search only: place balls ONLY on the residues "
                         "this file names. Filtered per structure by its "
                         "pdb_id/instance column")
    ap.add_argument("--min-idf", type=float, default=None, metavar="X",
                    help="a hit counts only if its idf is >= X. Applied on "
                         "FoldDisco's idf column at the join, because this "
                         "build's own --score works on a different, "
                         "normalised scale and silently empties the output. "
                         "IDF is index-dependent, so choose per index: "
                         "median per-hit idf is ~478-787 on whole-PDB and "
                         "~43.7 on afdb_swissprot")
    ap.add_argument("--top", type=int, default=700,
                    help="FoldDisco --top, structures kept per query. 40 was "
                         "measured to censor 72-98%% of motifs, which makes "
                         "the hit count a measure of search effort rather "
                         "than of how common a motif is [700]")
    ap.add_argument("-d", "--database", default=None, metavar="NAME",
                    help="named reference database from "
                         "[folddisco.databases]; default: "
                         "[folddisco].default_database. List them with "
                         "`ibc-search list-databases`")
    ap.add_argument("-j", "--jobs", type=int, default=6)
    ap.add_argument("--force", action="store_true",
                    help="re-run members that already have hits")
    ap.add_argument("--search-arg", action="append", default=[],
                    help="extra argument for folddisco_search.sh; repeatable")
    ap.add_argument("--no-finalize", action="store_true",
                    help="do not run finalize_run.py afterwards. The standard "
                         "00_figures/01_statistics layout is produced by "
                         "default, so that every run has the same shape")
    ap.add_argument("--no-pymol", action="store_true",
                    help="finalize without the per-structure PyMOL renders")
    args = ap.parse_args(argv)

    config = load_config()
    # Resolve the database HERE rather than letting the shell script fall
    # back to its own default, so the run is searched against the index the
    # caller named and the choice is printed once, up front, instead of being
    # inferred afterwards from a path in a log.
    entry = resolve_database(config, args.database)
    if not entry["index_path"].exists():
        print(f"error: database {entry['name']!r} points at "
              f"{entry['index_path']}, which does not exist", file=sys.stderr)
        return 2
    print(f"database: {entry['name']}  ({entry['index_path']})")
    if entry["idf_params_path"] is None:
        print("  no fitted IDF parameters for this database: hits can be "
              "counted but not z-scored (see idf_zscore.py fit)")
    search_args = ["-i", str(entry["index_path"]), *args.search_arg]

    panel = args.panel or (args.run_dir / "panel.tsv")
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
                      tuple(search_args), args.method, args.seed_residues,
                      args.seed_rounds, args.min_idf,
                      args.background_mode, args.min_background,
                      args.background_rsa))

    if skipped:
        print(f"{skipped} members already searched, skipped (--force re-runs)")
    if not tasks:
        print("nothing to do")
        return 0
    print(f"{args.method}: {len(tasks)} structures at radius "
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

    # Finalise by default. The layout is only standard if producing it is
    # the default path rather than a step someone has to remember, which is
    # how runs drifted into having different figures in the first place.
    if args.no_finalize:
        print("skipping finalize (--no-finalize); run "
              "scripts/finalize_run.py yourself for the standard layout")
        return 0 if bad == 0 else 1
    if ok == 0:
        print("nothing searched, so nothing to finalise")
        return 1

    print("\nfinalising: standard figures, statistics and per-structure renders")
    cmd = [str(REPO / ".venv/bin/python"), str(REPO / "scripts/finalize_run.py"),
           str(args.run_dir), "-j", str(args.jobs)]
    if args.no_pymol:
        cmd.append("--no-pymol")
    final = subprocess.run(cmd, cwd=str(REPO))
    return 0 if (bad == 0 and final.returncode == 0) else 1


if __name__ == "__main__":
    sys.exit(main())
