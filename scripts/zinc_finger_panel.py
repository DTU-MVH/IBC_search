"""Build a zinc-finger calibration panel and label its zinc sites.

Why this panel exists
---------------------
The +1.7 call threshold was never calibrated against anything observed. On
designed binders there is no ground truth at all, and on antibody-antigen
complexes the score turned out to be at chance, which leaves the threshold
unanchored in both directions.

A zinc finger is a useful anchor precisely because it is *not* an epitope.
It is a compact, structurally unambiguous, recurrent motif whose residues can
be identified from coordinates alone, with no annotation and no judgement
call: the residues around a bound Zn ion. If the motif score carries
structural information at all, a zinc site should separate from the rest of
the surface; if it cannot separate something this crisp, no threshold on it
will separate an epitope.

The expected direction is DOWNWARD. IDF rewards rarity, and zinc fingers are
among the most frequently deposited small motifs in the PDB, so a working
score should give them *low* z. Finding them unremarkable would mean the
score does not see them; finding them high would mean it is not measuring
rarity. Both outcomes are informative, which is what makes this a
calibration and not a fishing trip.

The label
---------
``is_epitope = 1`` for exactly the residues the PDB entry itself lists in a
zinc binding site -- the ``struct_site`` / ``struct_site_gen`` records for a
site whose ligand is ZN. Nothing is added by proximity. The column is named
``is_epitope`` only because that is the label column the scoring and
rendering code already reads -- **these residues are zinc sites, not
epitopes**, and nothing here says otherwise.

An earlier version labelled every residue within a radius of the metal. That
set was strictly larger (101 residues against 60) and included a second shell
the entry does not call part of the site, which makes the positive class a
judgement of this script rather than of the deposition. The radius is still
computed, but only as the reported ``distance_shell`` column, never as the
label.

**What this label is NOT.** Most of these ``struct_site`` records carry
``pdbx_evidence_code = Software``: the PDB generated them from proximity to
the ligand, so they are a curated-looking distance rule rather than an
author's assertion. Only 1VFY is marked ``Author``. The label is therefore
"what the entry lists as the zinc site", which is the most defensible
available definition, not an independent expert annotation.

It is also not the zinc-finger *domain*: the canonical motif is a fold of
~25-30 residues, while a listed site is the coordinating sphere. Residues are
matched on ``auth_asym_id``/``auth_seq_id``, the same numbering the scan
emits.

Structures are read from the same RCSB cache entry the scan parses, so the
labels cannot drift away from the coordinates that were scanned.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ibc_search.config import load_config  # noqa: E402
from ibc_search.structure_io import fetch_pdb_entry  # noqa: E402

import biotite.structure.io.pdbx as pdbx_io  # noqa: E402

WATER = frozenset({"HOH", "DOD", "WAT", "H2O"})

# 5 human, 5 non-human, chosen for distinct zinc-finger families and (on the
# non-human side) distinct organisms, all single-protein entries with Zn
# coordinates deposited and no nucleic acid.
PANEL: list[tuple[str, str, str, str]] = [
    # pdb_id, organism group, organism, zinc-finger family
    ("3IUF", "human", "Homo sapiens", "C2H2 (ubi-d4/DPF2)"),
    ("3O70", "human", "Homo sapiens", "PHD (PHF13)"),
    ("3Q1D", "human", "Homo sapiens", "B-box (TRIM54)"),
    ("3T7L", "human", "Homo sapiens", "FYVE (endofin/ZFYVE16)"),
    ("4Z0O", "human", "Homo sapiens", "CW-type (ZCWPW2)"),
    ("1PZW", "non_human", "Drosophila melanogaster", "ZAD"),
    ("1VFY", "non_human", "Saccharomyces cerevisiae", "FYVE (Vps27p)"),
    ("1LV3", "non_human", "Escherichia coli", "yacG zinc finger"),
    ("1WFH", "non_human", "Arabidopsis thaliana", "zf-AN1"),
    ("1U85", "non_human", "Mus musculus", "C2H2 (BKLF finger 3)"),
]


def pdb_zinc_sites(block) -> tuple[set[tuple[str, int]], str]:
    """Residues the entry itself lists in a ZN binding site, and the evidence.

    ``struct_site`` names each site and the ligand it is about;
    ``struct_site_gen`` lists that site's residues. Only sites whose ligand is
    ZN are taken, and the ligand and any waters are dropped from the residue
    list, leaving the polymer residues the deposition calls the zinc site.
    """
    if "struct_site" not in block or "struct_site_gen" not in block:
        return set(), "none"
    site = block["struct_site"]
    if "pdbx_auth_comp_id" not in site:
        return set(), "none"
    ids = site["id"].as_array(str)
    comps = site["pdbx_auth_comp_id"].as_array(str)
    evidence = (site["pdbx_evidence_code"].as_array(str)
                if "pdbx_evidence_code" in site else [])
    zinc_ids = {i for i, c in zip(ids, comps) if c.strip().upper() == "ZN"}
    codes = {e for i, e in zip(ids, evidence) if i in zinc_ids}
    if not zinc_ids:
        return set(), "none"

    gen = block["struct_site_gen"]
    out: set[tuple[str, int]] = set()
    for sid, chain, seq, comp in zip(
        gen["site_id"].as_array(str), gen["auth_asym_id"].as_array(str),
        gen["auth_seq_id"].as_array(str), gen["label_comp_id"].as_array(str),
    ):
        if sid not in zinc_ids or comp.strip().upper() in ({"ZN"} | WATER):
            continue
        try:
            out.add((str(chain), int(seq)))
        except ValueError:
            continue
    return out, "/".join(sorted(codes)) or "unspecified"


def label_structure(path: Path, radius: float, coordination: float):
    """-> (ordered residues, PDB site keys, coordinating keys, shell, n_zinc, evidence).

    The LABEL is the PDB's own site list. The distance sets are computed only
    so the report can say how the two relate; they never widen the label.
    """
    # The RCSB cache holds the download as served, which is gzipped.
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as fh:
            cif = pdbx_io.CIFFile.read(fh)
    else:
        cif = pdbx_io.CIFFile.read(str(path))
    array = pdbx_io.get_structure(cif, model=1)
    array = array[(array.element != "H") & ~np.isin(array.res_name, list(WATER))]

    zinc = array[array.element == "ZN"]
    protein = array[array.element != "ZN"]

    keys = [(str(c), int(r), str(n)) for c, r, n
            in zip(protein.chain_id, protein.res_id, protein.res_name)]
    seen, ordered = set(), []
    for k in keys:
        if k[:2] not in seen:
            seen.add(k[:2])
            ordered.append(k)

    listed, evidence = pdb_zinc_sites(cif.block)

    if zinc.array_length() == 0:
        return ordered, listed, set(), set(), 0, evidence

    tree = cKDTree(zinc.coord)
    near = tree.query_ball_point(protein.coord, r=radius, return_length=True) > 0
    coord_near = tree.query_ball_point(
        protein.coord, r=coordination, return_length=True) > 0
    shell = {k[:2] for k, c in zip(keys, near) if c}
    coordinating = {k[:2] for k, c in zip(keys, coord_near) if c}
    return ordered, listed, coordinating, shell, zinc.array_length(), evidence


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("--shell-radius", type=float, default=5.0,
                    help="reported only, never the label: distance from a Zn "
                         "at which a residue is called second shell [5.0]")
    ap.add_argument("--coordination", type=float, default=2.8,
                    help="distance defining direct Zn coordination, reported "
                         "separately [2.8]")
    args = ap.parse_args(argv)

    config = load_config()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    panel_path = args.out_dir / "panel.tsv"
    epi_path = args.out_dir / "epitopes_pdb_zinc_site.tsv"

    total_res = total_site = total_coord = total_shell = 0
    rows = []
    with epi_path.open("w", newline="") as ef:
        ew = csv.writer(ef, delimiter="\t", lineterminator="\n")
        ew.writerow(["instance", "chain", "res_id", "res_name", "is_epitope",
                     "is_zinc_coordinating", "in_distance_shell"])
        for pdb_id, group, organism, family in PANEL:
            source = fetch_pdb_entry(pdb_id, config)
            ordered, site, coordinating, shell, n_zn, evidence = label_structure(
                source.path, args.shell_radius, args.coordination)
            if not n_zn:
                print(f"  {pdb_id}: no Zn in the file, skipped", file=sys.stderr)
                continue
            if not site:
                print(f"  {pdb_id}: the entry lists no ZN site, skipped",
                      file=sys.stderr)
                continue
            for chain, res_id, res_name in ordered:
                ew.writerow([pdb_id, str(chain), int(res_id), str(res_name),
                             int((chain, res_id) in site),
                             int((chain, res_id) in coordinating),
                             int((chain, res_id) in shell)])
            total_res += len(ordered)
            total_site += len(site)
            total_coord += len(coordinating)
            total_shell += len(shell)
            rows.append([pdb_id, pdb_id, group, organism, family, n_zn,
                         evidence, len(ordered), len(site), len(coordinating),
                         len(shell),
                         f"{len(site) / len(ordered):.4f}" if ordered else ""])
            print(f"  {pdb_id}  {group:<9} {family:<22} {n_zn} Zn  "
                  f"{len(site):>3}/{len(ordered):<3} listed  "
                  f"({len(coordinating)} coordinating, {len(shell)} within "
                  f"{args.shell_radius} A)  evidence={evidence}")

    with panel_path.open("w", newline="") as pf:
        w = csv.writer(pf, delimiter="\t", lineterminator="\n")
        w.writerow(["instance", "pdb_id", "organism_group", "organism",
                    "zinc_finger_family", "n_zinc", "site_evidence",
                    "n_residues", "n_zinc_site", "n_coordinating",
                    "n_distance_shell", "zinc_site_fraction"])
        w.writerows(rows)

    (args.out_dir / "panel_manifest.json").write_text(json.dumps({
        "built_utc": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        "purpose": "threshold calibration against a known structural motif",
        "label_meaning": (
            "is_epitope marks exactly the residues the PDB entry lists in a "
            "ZN binding site (struct_site / struct_site_gen). Nothing is "
            "added by proximity. These are ZINC SITES, not epitopes; the "
            "column is named is_epitope only so the existing scoring and "
            "rendering code reads it."),
        "label_source": "pdb struct_site records, ligand ZN",
        "label_caveat": (
            "Most of these records carry pdbx_evidence_code=Software: the PDB "
            "generated them from proximity to the ligand, so they are not an "
            "author's assertion. The per-entry code is in panel.tsv."),
        "shell_radius_angstrom_reported_only": args.shell_radius,
        "coordination_angstrom": args.coordination,
        "structures": len(rows),
        "human": sum(1 for r in rows if r[2] == "human"),
        "non_human": sum(1 for r in rows if r[2] == "non_human"),
        "residues": total_res,
        "distance_shell_residues_not_labelled": total_shell,
        "zinc_site_residues": total_site,
        "coordinating_residues": total_coord,
    }, indent=2))

    print(f"\n{len(rows)} structures, {total_res} residues, "
          f"{total_site} listed by the PDB as zinc site "
          f"({100 * total_site / total_res:.1f}%), "
          f"of which {total_coord} directly coordinating")
    print(f"  for reference only, {total_shell} residues lie within "
          f"{args.shell_radius} A of a Zn; they are NOT labelled")
    print(f"  {panel_path}\n  {epi_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
