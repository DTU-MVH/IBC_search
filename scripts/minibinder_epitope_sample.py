"""Select minibinder-sized antigens and label their antibody binding sites.

Picks a sample of small antigens from the SAbDab2 reference -- one per
antigen cluster, so the sample is not dominated by the handful of antigens
the field has solved hundreds of times -- and records, for every residue,
whether it lies within ``--contact`` Angstroms of any antibody heavy atom.

Those labelled residues are the ground truth the motif score is tested
against: an epitope here is observed contact in a solved complex, not a
prediction.

Why these structures rather than the designed minibinders
---------------------------------------------------------
The ADA minibinder set is monomer predictions -- one chain, no antibody --
so an antibody contact distance cannot be computed from it at all. SAbDab2
carries the antibody and the antigen in one file, so it is the only source
where the binding site is observed rather than assumed. Restricting to
50-150 residues keeps the antigens in the size range of a designed binder.
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

import biotite.structure.io.pdbx as pdbx_io

WATER = frozenset({"HOH", "DOD", "WAT", "H2O"})


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("index_dir", type=Path, help="sabdab2_antigen_index directory")
    ap.add_argument("splits_dir", type=Path, help="splits_final directory")
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("-n", "--sample", type=int, default=40)
    ap.add_argument("--contact", type=float, default=5.5)
    ap.add_argument("--min-residues", type=int, default=50)
    ap.add_argument("--max-residues", type=int, default=150)
    ap.add_argument("--min-epitope", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(argv)

    structures = list(
        csv.DictReader((args.index_dir / "structures.tsv").open(), delimiter="\t")
    )
    eligible = [
        r for r in structures
        if args.min_residues <= int(r["n_antigen_residues"]) <= args.max_residues
        and int(r["n_epitope_residues"]) >= args.min_epitope
    ]
    # One per cluster: the field has solved some antigens hundreds of times,
    # and a sample weighted by that would measure those few antigens.
    by_cluster: dict[str, list] = defaultdict(list)
    for r in eligible:
        by_cluster[r["ag_cluster"]].append(r)
    rng = random.Random(args.seed)
    clusters = sorted(by_cluster)
    rng.shuffle(clusters)
    chosen = [rng.choice(by_cluster[c]) for c in clusters[: args.sample]]
    print(f"{len(eligible)} eligible antigens in {len(by_cluster)} clusters; "
          f"sampled {len(chosen)}", flush=True)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    panel = args.out_dir / "panel.tsv"
    epis = args.out_dir / f"epitopes_{args.contact}A.tsv"

    n_res = n_epi = 0
    kept = []
    with panel.open("w", newline="") as pf, epis.open("w", newline="") as ef:
        pw = csv.writer(pf, delimiter="\t", lineterminator="\n")
        ew = csv.writer(ef, delimiter="\t", lineterminator="\n")
        pw.writerow(["instance", "pdb_id", "antigen_chains", "antibody_chains",
                     "n_residues", "n_epitope", "epitope_fraction",
                     "ag_cluster", "ab_ag_split", "antigen_file"])
        ew.writerow(["instance", "chain", "res_id", "res_name", "is_epitope"])

        for row in chosen:
            inst = row["instance"]
            src = args.splits_dir / f"{inst}.cif"
            if not src.is_file():
                continue
            try:
                arr = pdbx_io.get_structure(pdbx_io.CIFFile.read(str(src)), model=1)
            except Exception:
                continue
            arr = arr[(arr.element != "H") & ~np.isin(arr.res_name, list(WATER))]

            ab = [c for c in (row["antibody_chains"] or "").split("/") if c]
            ag = [c for c in (row["antigen_chains"] or "").split("/") if c]
            abs_ = arr[np.isin(arr.chain_id, ab)]
            ags = arr[np.isin(arr.chain_id, ag)]
            if abs_.array_length() == 0 or ags.array_length() == 0:
                continue

            close = cKDTree(abs_.coord).query_ball_point(
                ags.coord, r=args.contact, return_length=True
            ) > 0
            keys = list(zip(ags.chain_id, ags.res_id, ags.res_name))
            epitope = {k[:2] for k, c in zip(keys, close) if c}

            seen, ordered = set(), []
            for k in keys:
                if k[:2] not in seen:
                    seen.add(k[:2]); ordered.append(k)
            for chain, res_id, res_name in ordered:
                ew.writerow([inst, str(chain), int(res_id), str(res_name),
                             int((chain, res_id) in epitope)])

            n_res += len(ordered); n_epi += len(epitope)
            kept.append(inst)
            pw.writerow([
                inst, row["pdb_id"], row["antigen_chains"], row["antibody_chains"],
                len(ordered), len(epitope),
                f"{len(epitope) / len(ordered):.4f}" if ordered else "",
                row["ag_cluster"], row["ab_ag_split"],
                str((args.index_dir / row["antigen_file"]).resolve()),
            ])

    (args.out_dir / "sample_manifest.json").write_text(json.dumps({
        "built_utc": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        "contact_distance_angstrom": args.contact,
        "size_range": [args.min_residues, args.max_residues],
        "min_epitope_residues": args.min_epitope,
        "one_per_cluster": True,
        "seed": args.seed,
        "structures": len(kept),
        "residues": n_res,
        "epitope_residues": n_epi,
        "epitope_fraction_overall": round(n_epi / n_res, 4) if n_res else None,
    }, indent=2))

    print(f"wrote {len(kept)} structures, {n_res} residues, "
          f"{n_epi} epitope ({100 * n_epi / n_res:.1f}%)")
    print(f"  {panel}\n  {epis}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
