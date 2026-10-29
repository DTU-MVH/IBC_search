"""Query FoldDisco with the annotated zinc motif itself, not the whole surface.

The surface scan asks "do zinc-site residues end up in high-scoring patches",
which routes the answer through two extra steps -- which balls happened to
cover them, and the MAX projection back onto residues. Both add noise, and
the MAX step adds a measured order-statistic bias.

This asks the question directly: take the residues the PDB entry lists as the
zinc site, hand exactly those to FoldDisco as one motif, and compare that
motif's IDF z-score against the z-scores of the ordinary surface motifs of
the same structure. One number against a distribution, with nothing projected
and nothing aggregated.

Two query sets can be written. ``whole`` is the default and the one that
answers the question as posed: every residue the entry annotates, as one
motif (4 or 8 residues).

``site`` splits a two-zinc domain into one motif per ion, which keeps every
query at 4 residues. That was built on the worry that merging two sites makes
a motif larger than the background it is compared against. Measured, the
worry is backwards: the merged 8-residue motifs span 15.6-17.7 A against a
background of median 13.6 A and max 21.6 A, so they sit inside its range
(20% of background motifs are >= 15.6 A), while the split 4-residue motifs
span 8.1-10.7 A and are *smaller* than the background median. ``site`` is
kept because a per-ion view is still worth having, not because it is the
better comparison.

For the five single-zinc entries in the panel the two modes emit identical
queries and identical scores; only the five two-zinc entries differ.

Self-hits are excluded by the search script, which reads the entry ID from
the run manifest; this script copies that manifest across so the exclusion
works exactly as it does for a scan.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ibc_search.config import load_config  # noqa: E402
from ibc_search.structure_io import fetch_pdb_entry  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from zinc_finger_panel import pdb_zinc_sites  # noqa: E402

import biotite.structure.io.pdbx as pdbx_io  # noqa: E402


def nearest_zinc(path: Path, residues: list[tuple[str, int]]) -> dict:
    """Which Zn ion each annotated residue belongs to."""
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as fh:
            cif = pdbx_io.CIFFile.read(fh)
    else:
        cif = pdbx_io.CIFFile.read(str(path))
    arr = pdbx_io.get_structure(cif, model=1)
    arr = arr[arr.element != "H"]
    zinc = arr[arr.element == "ZN"]
    prot = arr[arr.element != "ZN"]

    want = set(residues)
    best: dict[tuple[str, int], tuple[float, int]] = {}
    for c, r, xyz in zip(prot.chain_id, prot.res_id, prot.coord):
        key = (str(c), int(r))
        if key not in want:
            continue
        d = np.linalg.norm(zinc.coord - xyz, axis=1)
        i = int(np.argmin(d))
        if key not in best or float(d[i]) < best[key][0]:
            best[key] = (float(d[i]), i)
    return {k: v[1] for k, v in best.items()}


def write_queries(out_dir: Path, structure: str, queries: list[tuple[str, str]]):
    """queries: (patch_id, query_string)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    with (out_dir / "folddisco_queries.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        for _, q in queries:
            w.writerow([structure, q])
    with (out_dir / "folddisco_queries.index.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["line_number", "patch_id", "structure_path",
                    "radius_angstrom", "scale", "n_members", "query_string"])
        for i, (pid, q) in enumerate(queries, 1):
            w.writerow([i, pid, structure, "", "zinc_site",
                        len([t for t in q.split(",") if t]), q])


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("--scan-dir", type=Path, default=None,
                    help="a completed surface-scan run; its epitopes_* label "
                         "file and manifests are reused")
    ap.add_argument("--entries", nargs="+", default=None,
                    help="PDB entry IDs to read zinc sites from directly, "
                         "with no prior scan. The site residues come from the "
                         "entry's own struct_site records.")
    ap.add_argument("--dedupe-sites", action="store_true",
                    help="keep one copy of each site that recurs across "
                         "chains, matched on its residue numbers. A complex "
                         "crystallised with several copies otherwise "
                         "contributes the same site many times and weights "
                         "the statistics by crystal packing.")
    ap.add_argument("--mode", choices=("whole", "site"), default="whole",
                    help="whole: every annotated residue of the structure as "
                         "one motif, which is the question 'is the annotated "
                         "zinc motif unusual'. site: one motif per Zn ion "
                         "[whole]")
    args = ap.parse_args(argv)

    if not args.scan_dir and not args.entries:
        print("error: pass --scan-dir or --entries", file=sys.stderr)
        return 2

    config = load_config()
    listed: dict[str, list[tuple[str, int]]] = defaultdict(list)
    epi = None
    if args.scan_dir:
        epi = next(iter(sorted(args.scan_dir.glob("epitopes_*.tsv"))), None)
        if epi is None:
            print(f"error: no epitope label file in {args.scan_dir}",
                  file=sys.stderr)
            return 2
        for r in csv.DictReader(epi.open(), delimiter="\t"):
            if r["is_epitope"] == "1":
                listed[r["instance"]].append((r["chain"], int(r["res_id"])))
    else:
        for entry in args.entries:
            src = fetch_pdb_entry(entry, config)
            if src.path.suffix == ".gz":
                with gzip.open(src.path, "rt") as fh:
                    block = pdbx_io.CIFFile.read(fh).block
            else:
                block = pdbx_io.CIFFile.read(str(src.path)).block
            site, _ = pdb_zinc_sites(block)
            if not site:
                print(f"  {entry}: entry lists no ZN site, skipped",
                      file=sys.stderr)
                continue
            listed[entry] = sorted(site)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    if args.scan_dir:
        shutil.copy(args.scan_dir / "panel.tsv", args.out_dir / "panel.tsv")

    total = 0
    for inst, residues in sorted(listed.items()):
        src = fetch_pdb_entry(inst, config)
        member = args.out_dir / inst
        member.mkdir(parents=True, exist_ok=True)
        # The manifest carries input.structure_id, which is what the search
        # script reads to exclude the entry's hits on itself.
        if args.scan_dir:
            shutil.copy(args.scan_dir / inst / "run_manifest.json",
                        member / "run_manifest.json")
        else:
            (member / "run_manifest.json").write_text(json.dumps(
                {"input": {"structure_id": inst, "pdb_id": inst,
                           "path": str(src.path)}}, indent=2))

        if args.mode == "whole":
            q = ",".join(f"{c}{r}" for c, r in sorted(residues))
            queries = [(f"{inst}_zinc_whole", q)]
        else:
            by_zn = nearest_zinc(src.path, residues)
            groups: dict[int, list[tuple[str, int]]] = defaultdict(list)
            for key, idx in by_zn.items():
                groups[idx].append(key)
            queries = [
                (f"{inst}_zinc_site{i}",
                 ",".join(f"{c}{r}" for c, r in sorted(g)))
                for i, (_, g) in enumerate(sorted(groups.items()))
            ]
            if args.dedupe_sites:
                # The same site in another chain has the same residue
                # numbers; keep the first copy only.
                seen_sig, unique = set(), []
                for pid_, q in queries:
                    sig = tuple(sorted(t[1:] for t in q.split(",") if t))
                    if sig in seen_sig:
                        continue
                    seen_sig.add(sig)
                    unique.append((pid_, q))
                queries = unique

        write_queries(member, str(src.path), queries)
        total += len(queries)
        sizes = [len([t for t in q.split(",") if t]) for _, q in queries]
        print(f"  {inst}  {len(queries)} query(ies), sizes {sizes}")

    (args.out_dir / "query_manifest.json").write_text(json.dumps({
        "mode": args.mode,
        "source_scan": str(args.scan_dir) if args.scan_dir else None,
        "label_file": str(epi) if epi else "pdb struct_site records (ZN)",
        "deduped_sites": bool(args.dedupe_sites),
        "structures": len(listed),
        "queries": total,
        "note": ("Each query is the residues the PDB entry lists as its zinc "
                 "site. Nothing was extracted by the surface scan."),
    }, indent=2))
    print(f"\n{len(listed)} structures, {total} zinc-motif queries -> {args.out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
