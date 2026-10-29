"""Render a scored PyMOL figure for every structure in a run.

Produces, inside each structure's own directory, a ``00_figures/`` folder
holding the ray-traced PNG, the PDB carrying the score in its B-factor
column, the PyMOL script and a session file.

Epitope spheres appear only when the run actually has epitope labels. With
no ground truth -- a designed binder, say -- the figure is the red/grey
prediction alone and claims nothing about correctness, which is the honest
default for a structure whose binding site is unknown.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
RENDER = REPO / "scripts" / "render_pymol.py"


def _render(task: tuple) -> tuple[str, bool, str]:
    inst, structure, scores_tsv, out_png, threshold, mode, pymol, has_epitope = task
    cmd = [
        str(REPO / ".venv/bin/python"), str(RENDER), structure, scores_tsv,
        "--instance", inst, "-o", out_png,
        "--threshold", str(threshold), "--mode", mode,
    ]
    if pymol:
        cmd += ["--pymol", pymol]
    if not has_epitope:
        cmd.append("--no-epitope")
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=str(REPO))
    ok = proc.returncode == 0 and Path(out_png).is_file()
    lines = (proc.stderr or proc.stdout).strip().splitlines()
    return inst, ok, "" if ok else (lines[-1] if lines else "no output")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_dir", type=Path, help="run directory holding panel.tsv")
    ap.add_argument("--scores", type=Path, default=None,
                    help="residue_scores.tsv [<run_dir>/residue_scores.tsv]")
    ap.add_argument("--threshold", type=float, default=1.7)
    ap.add_argument("--mode", choices=("prediction", "confusion"),
                    default="prediction")
    ap.add_argument("--figures-dir", default="00_figures",
                    help="subdirectory created inside each structure folder")
    ap.add_argument("--pymol", default=None)
    ap.add_argument("-j", "--jobs", type=int, default=6)
    args = ap.parse_args(argv)

    scores_tsv = args.scores or (args.run_dir / "residue_scores.tsv")
    if not scores_tsv.is_file():
        print(f"error: no residue scores at {scores_tsv}", file=sys.stderr)
        return 2

    # Which structures actually carry epitope labels?
    labelled: set[str] = set()
    with scores_tsv.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("is_epitope") == "1":
                labelled.add(r.get("instance", ""))

    config = None
    tasks = []
    for row in csv.DictReader((args.run_dir / "panel.tsv").open(), delimiter="\t"):
        inst = row["instance"]
        struct = row.get("antigen_file") or row.get("structure_file")
        # A panel may name a PDB entry instead of a file, as the zinc-finger
        # and reference panels do. Resolve it through the same cache the scan
        # used, so the rendered coordinates are the scanned coordinates.
        if not struct and row.get("pdb_id"):
            if config is None:
                sys.path.insert(0, str(REPO / "src"))
                from ibc_search.config import load_config
                config = load_config()
            from ibc_search.structure_io import fetch_pdb_entry
            try:
                struct = str(fetch_pdb_entry(row["pdb_id"], config).path)
            except Exception as exc:
                print(f"  {inst}: cannot resolve {row['pdb_id']} ({exc})")
                continue
        if not struct or not Path(struct).is_file():
            print(f"  {inst}: structure file missing, skipped")
            continue
        fig_dir = args.run_dir / inst / args.figures_dir
        fig_dir.mkdir(parents=True, exist_ok=True)
        tasks.append((
            inst, struct, str(scores_tsv),
            str(fig_dir / f"{inst}_motif_score.png"),
            args.threshold, args.mode, args.pymol, inst in labelled,
        ))

    if not tasks:
        print("error: nothing to render", file=sys.stderr)
        return 2
    print(f"rendering {len(tasks)} structures, {args.jobs} at a time "
          f"({len(labelled)} with epitope labels)", flush=True)

    ok = bad = 0
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(_render, t) for t in tasks]
        for fut in as_completed(futures):
            inst, good, detail = fut.result()
            if good:
                ok += 1
            else:
                bad += 1
                print(f"  {inst}: FAILED {detail}")
    print(f"\n{ok} rendered, {bad} failed")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
