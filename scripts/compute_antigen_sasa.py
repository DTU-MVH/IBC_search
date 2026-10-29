"""Per-residue SASA/RSA for the SAbDab2 antigen reference, cached once.

Why this exists
---------------
``ibc-search scan`` computes exposure on the structure being scanned -- the
minibinder -- so motifs are built from exposed residues by construction. It
never parses the reference antigens a motif is searched against, so their
exposure is unknown, and two things depend on it:

1. **A buried antigen region cannot be an antibody epitope.** No antibody
   reaches it. Such a match carries no information either way, yet it dilutes
   the denominator of an ``epitope_transfer`` style feature -- and because the
   dilution depends on where a motif's matches happen to land, it adds
   per-motif noise rather than a constant offset. Measured on a 40-antigen
   sample: 12% of matched residues were buried, and excluding them moved the
   epitope rate from 9.2% to 10.5%.

2. **The base rate differs by exposure.** Epitope residues are necessarily
   exposed, so they concentrate in the exposed subset: on that same sample,
   11.5% of all residues were epitope against 14.0% of exposed ones. Comparing
   matched residues to the all-residue rate therefore understates the bar and
   flatters any enrichment.

The cost is one-off. The reference does not change per scan, so this is
computed once and reused for every protein that comes in.

RSA uses the same Tien et al. (2013) theoretical maxima as
``ibc_search.exposure``, so RSA means the same thing on both sides of a
search. Note ``--points`` defaults to 100 rather than the pipeline's 1000:
that is accurate enough to classify exposed/buried and to estimate a base
rate, and 10x cheaper over ~6000 structures. Absolute RSA values are
consequently a little noisier than the query side's.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from collections import defaultdict
from datetime import datetime, timezone
from multiprocessing import Pool
from pathlib import Path

import numpy as np

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx_io

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ibc_search.exposure import TIEN2013_THEORETICAL_MAX_ACC  # noqa: E402

EXPOSED_RSA = 0.05
"""RSA at or above which a residue counts as exposed.

Matches ``[ball_search].member_rsa_threshold``: the floor the pipeline already
uses to admit a residue into a motif at all. A reference residue below it
could not have been a motif member on the query side either.
"""


def _one(args: tuple[str, str, int]) -> tuple[str, list[tuple] | None, str]:
    """SASA for one antigen file. Returns (instance, rows, error)."""
    path, instance, points = args
    try:
        array = pdbx_io.get_structure(pdbx_io.CIFFile.read(path), model=1)
        array = array[array.element != "H"]
        if array.array_length() == 0:
            return instance, None, "no_heavy_atoms"
        sasa = np.nan_to_num(struc.sasa(array, point_number=points))

        acc: dict[tuple, float] = defaultdict(float)
        names: dict[tuple, str] = {}
        for chain, res_id, ins, name, value in zip(
            array.chain_id, array.res_id, array.ins_code, array.res_name, sasa
        ):
            key = (str(chain), int(res_id), str(ins or ""))
            acc[key] += float(value)
            names[key] = str(name)

        rows = []
        for key, area in acc.items():
            max_acc = TIEN2013_THEORETICAL_MAX_ACC.get(names[key])
            # A non-standard residue has no Tien maximum. Report the absolute
            # SASA and leave RSA empty rather than inventing a denominator;
            # exposure is then decided on area alone.
            if max_acc:
                rsa = area / max_acc
                exposed = int(rsa >= EXPOSED_RSA)
                rsa_s = f"{rsa:.4f}"
            else:
                rsa_s = ""
                exposed = int(area >= 10.0)
            rows.append(
                (instance, key[0], key[1], key[2], names[key],
                 f"{area:.2f}", rsa_s, exposed)
            )
        return instance, rows, ""
    except Exception as exc:  # noqa: BLE001 - one bad file must not kill the pass
        return instance, None, type(exc).__name__


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("index_dir", type=Path, help="sabdab2_antigen_index directory")
    ap.add_argument("-j", "--jobs", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    ap.add_argument("--points", type=int, default=100, help="SASA sphere points [100]")
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args(argv)

    root: Path = args.index_dir
    structures = list(csv.DictReader((root / "structures.tsv").open(), delimiter="\t"))
    if args.limit:
        structures = structures[: args.limit]

    tasks = [
        (str(root / r["antigen_file"]), r["instance"], args.points)
        for r in structures
        if (root / r["antigen_file"]).is_file()
    ]
    print(f"{len(tasks)} antigens, {args.jobs} processes, point_number={args.points}",
          flush=True)

    out_path = root / "residue_sasa.tsv"
    failures: dict[str, int] = {}
    n_res = n_exp = n_done = 0

    with out_path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["instance", "chain", "res_id", "ins_code", "res_name",
                    "sasa_a2", "rsa", "is_exposed"])
        with Pool(args.jobs) as pool:
            for instance, rows, err in pool.imap_unordered(_one, tasks, chunksize=8):
                n_done += 1
                if rows is None:
                    failures[err] = failures.get(err, 0) + 1
                    continue
                for row in rows:
                    w.writerow(row)
                    n_res += 1
                    n_exp += row[7]
                if n_done % 500 == 0:
                    print(f"  {n_done}/{len(tasks)}", flush=True)

    manifest = {
        "built_utc": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        "antigens_processed": n_done,
        "failures": failures,
        "residues": n_res,
        "exposed_residues": n_exp,
        "exposed_fraction": round(n_exp / n_res, 4) if n_res else None,
        "sasa_point_number": args.points,
        "exposed_rsa_threshold": EXPOSED_RSA,
        "rsa_reference": "tien2013_theoretical",
        "note": (
            "point_number=100 (pipeline scan uses 1000): adequate for "
            "exposed/buried classification and base rates, 10x cheaper here"
        ),
    }
    (root / "residue_sasa_manifest.json").write_text(json.dumps(manifest, indent=2))

    print(f"\nwrote {n_res} residues to {out_path}")
    print(f"  exposed (RSA >= {EXPOSED_RSA}): {n_exp} = {100 * n_exp / n_res:.1f}%"
          if n_res else "")
    if failures:
        print(f"  failures: {failures}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
