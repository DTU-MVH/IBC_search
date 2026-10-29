"""Build an antigen-only structure set with epitope labels, for FoldDisco.

Turns the SAbDab2 ML release into a searchable reference in which **every
residue carries a known epitope label**, so a FoldDisco hit can be asked the
question that matters: did the motif match the part of that antigen an
antibody actually binds, or somewhere else on it?

What it does, per structure
---------------------------
1. Reads the complex (antibody cropped to the variable region + antigen
   chains in one file, per the SAbDab2 README).
2. Labels an antigen residue as **epitope** when any of its heavy atoms lies
   within ``--contact`` Angstroms of any antibody heavy atom.
3. Writes the **antigen chains only** to a new file. Stripping the antibody is
   not cosmetic: left in, motifs would match framework and CDR geometry, which
   is a different object and would pollute every downstream number.

Then one ``folddisco index`` is built over that directory.

Selection
---------
Only antigen chains typed ``PROTEIN`` are kept.

``PEPTIDE`` antigens are **deliberately excluded.** A 13-residue peptide is
almost entirely epitope, so its background epitope fraction approaches 1 and
including such entries would inflate any enrichment measured against
background -- manufacturing the signal the reference exists to test. SUGAR,
HAPTEN and ION antigens are excluded as non-polymer.

Not computed here
-----------------
The background for an enrichment should be the fraction of *solvent-exposed*
antigen residues that are epitope, since motifs only ever match exposed
surface. That needs a SASA pass per antigen and is left to the analysis step;
this script records total and epitope residue counts, which give only a
crude all-residue background.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

import biotite.structure.io.pdbx as pdbx_io

WATER = frozenset({"HOH", "DOD", "WAT", "H2O"})


def _heavy_polymer(array):
    """Heavy atoms of polymer residues: no hydrogens, no water."""
    keep = (array.element != "H") & ~np.isin(array.res_name, list(WATER))
    return array[keep]


def _residue_keys(array):
    return list(zip(array.chain_id, array.res_id, array.ins_code))


def build(
    splits_dir: Path,
    out_dir: Path,
    contact: float,
    limit: int | None,
    split_csv: str,
) -> int:
    rows = list(csv.DictReader((splits_dir / split_csv).open()))

    selected = []
    for row in rows:
        chains = (row.get("agchains") or "").strip()
        types = (row.get("agtypes") or "").strip()
        if not chains or not types:
            continue
        chain_list, type_list = chains.split("/"), types.split("/")
        if len(chain_list) != len(type_list):
            continue  # malformed pairing; refuse to guess
        protein = [c for c, t in zip(chain_list, type_list) if t == "PROTEIN"]
        if protein:
            selected.append((row, protein))

    if limit:
        selected = selected[:limit]
    print(f"{len(rows)} rows; {len(selected)} with >=1 PROTEIN antigen chain", flush=True)

    ag_dir = out_dir / "antigens"
    ag_dir.mkdir(parents=True, exist_ok=True)

    ep_path = out_dir / "epitope_residues.tsv"
    st_path = out_dir / "structures.tsv"
    skipped: dict[str, int] = {}
    n_ok = n_res = n_epi = 0

    with ep_path.open("w", newline="") as ep_f, st_path.open("w", newline="") as st_f:
        ep = csv.writer(ep_f, delimiter="\t", lineterminator="\n")
        st = csv.writer(st_f, delimiter="\t", lineterminator="\n")
        ep.writerow(["instance", "chain", "res_id", "ins_code", "res_name", "is_epitope"])
        st.writerow([
            "instance", "pdb_id", "antigen_chains", "antibody_chains",
            "n_antigen_residues", "n_epitope_residues", "epitope_fraction_all_residues",
            "ag_cluster", "ab_ag_split", "antigen_file",
        ])

        for i, (row, ag_chains) in enumerate(selected, 1):
            instance = row["INSTANCE"]
            src = splits_dir / f"{instance}.cif"
            if not src.is_file():
                skipped["missing_cif"] = skipped.get("missing_cif", 0) + 1
                continue

            ab_chains = [c for c in (row.get("Hchain"), row.get("Lchain")) if c]
            try:
                array = pdbx_io.get_structure(pdbx_io.CIFFile.read(str(src)), model=1)
            except Exception:
                skipped["parse_failed"] = skipped.get("parse_failed", 0) + 1
                continue

            array = _heavy_polymer(array)
            ab = array[np.isin(array.chain_id, ab_chains)] if ab_chains else array[:0]
            ag = array[np.isin(array.chain_id, ag_chains)]

            if ag.array_length() == 0:
                skipped["no_antigen_atoms"] = skipped.get("no_antigen_atoms", 0) + 1
                continue

            # FoldDisco cannot read a multi-character chain ID: indexing such a
            # file panics with "Chain name should be provided" (verified -- the
            # 68 affected structures killed the whole index build). That is
            # consistent with the query token grammar, where a chain is a
            # single character, so a hit on such a chain could not be expressed
            # or attributed anyway. Refuse them rather than ship an index that
            # cannot be built.
            if any(len(str(c)) != 1 for c in set(ag.chain_id)):
                skipped["multi_char_chain"] = skipped.get("multi_char_chain", 0) + 1
                continue
            if ab.array_length() == 0:
                # Without antibody atoms there is no epitope to label, so the
                # entry would enter the index unlabelled. Refuse it.
                skipped["no_antibody_atoms"] = skipped.get("no_antibody_atoms", 0) + 1
                continue

            # Epitope = antigen residues with a heavy atom within `contact` A
            # of any antibody heavy atom.
            tree = cKDTree(ab.coord)
            close = tree.query_ball_point(ag.coord, r=contact, return_length=True) > 0
            keys = _residue_keys(ag)
            epitope_keys = {k for k, c in zip(keys, close) if c}

            residue_keys: list = []
            seen = set()
            for k in keys:
                if k not in seen:
                    seen.add(k)
                    residue_keys.append(k)
            name_of = dict(zip(keys, ag.res_name))

            for k in residue_keys:
                ep.writerow([
                    instance, k[0], k[1], k[2] or "", name_of[k],
                    int(k in epitope_keys),
                ])

            dest = ag_dir / f"{instance}.cif"
            out = pdbx_io.CIFFile()
            pdbx_io.set_structure(out, ag)
            out.write(str(dest))

            n_total = len(residue_keys)
            n_ep = len(epitope_keys)
            st.writerow([
                instance, row.get("PDB_ID", ""), "/".join(ag_chains),
                "/".join(ab_chains), n_total, n_ep,
                f"{n_ep / n_total:.4f}" if n_total else "",
                row.get("agclusters", ""), row.get("ab_ag_split", ""),
                f"antigens/{instance}.cif",
            ])
            n_ok += 1
            n_res += n_total
            n_epi += n_ep

            if i % 250 == 0:
                print(f"  {i}/{len(selected)} processed, {n_ok} written", flush=True)

    manifest = {
        "built_utc": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        "source_splits_dir": str(splits_dir),
        "split_csv": split_csv,
        "contact_distance_angstrom": contact,
        "antigen_types_kept": ["PROTEIN"],
        "antigen_types_excluded": {
            "PEPTIDE": "near-100% epitope; would inflate background enrichment",
            "SUGAR/HAPTEN/ION": "non-polymer",
        },
        "also_excluded": {
            "multi_char_chain": (
                "FoldDisco panics on a multi-character chain ID, and its query "
                "token grammar cannot express one, so such a hit would be "
                "unattributable"
            ),
        },
        "rows_in_split_csv": len(rows),
        "rows_with_protein_antigen": len(selected),
        "structures_written": n_ok,
        "skipped": skipped,
        "antigen_residues_total": n_res,
        "epitope_residues_total": n_epi,
        "epitope_fraction_all_residues": round(n_epi / n_res, 4) if n_res else None,
    }
    (out_dir / "build_manifest.json").write_text(json.dumps(manifest, indent=2))

    print(f"\nwrote {n_ok} antigen-only structures to {ag_dir}")
    print(f"  antigen residues: {n_res}, epitope: {n_epi} "
          f"({100 * n_epi / n_res:.1f}% of all residues)" if n_res else "")
    if skipped:
        print(f"  skipped: {skipped}")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("splits_dir", type=Path, help="SAbDab2 splits_final directory")
    ap.add_argument("out_dir", type=Path, help="where to write the antigen set")
    ap.add_argument("--contact", type=float, default=5.0,
                    help="heavy-atom contact distance in Angstroms [5.0]")
    ap.add_argument("--split-csv", default="abag_split.csv",
                    help="use the antigen-aware split [abag_split.csv]")
    ap.add_argument("--limit", type=int, default=None, help="process only N (smoke test)")
    args = ap.parse_args(argv)
    return build(args.splits_dir, args.out_dir, args.contact, args.limit, args.split_csv)


if __name__ == "__main__":
    sys.exit(main())
