"""Are a run's FoldDisco hits enriched for human structures, or is that the PDB?

The PDB is not a uniform sample of life. Human chains are a large, and human
*protein* chains a very large, share of what has been deposited, so "most of
our hits are human" is the expected result of hitting anything at all. A
human fraction only means something against the fraction the search could
have returned by chance.

The denominator used here is therefore **the indexed universe**: every
structure in the FoldDisco index that was actually searched, read from its
``.lookup`` sidecar, not the PDB as a whole and not some external census.
Anything outside the index cannot be retrieved, so including it would deflate
the background and manufacture enrichment.

Unit of counting
----------------
A *chain*, not an entry. A motif matches particular chains, and an entry may
hold chains from several organisms, so resolving an entry to one organism
throws away the distinction that matters. Both numerator and denominator are
counted the same way, over ``(pdb id, chain)`` pairs in SIFTS.

Two numerators are reported, because they answer different questions:

``rows``
    every hit row. This is what the score "sees", but one popular entry
    matched by forty motifs counts forty times.
``distinct``
    each matched chain counted once per run. Immune to that inflation, and
    the honest figure for "how much of the PDB did we touch, and what was
    it".

``unknown`` is its own category throughout
------------------------------------------
SIFTS maps only chains it could resolve to UniProt; synthetic constructs, de
novo designs and unmapped chains are simply absent. A missing chain means
"not known", never "not human". Folding it into the non-human count would
dilute every fraction, so unknown chains are excluded from the fractions and
reported separately -- in the background exactly as in the hits, since the
same blind spot applies to both.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from top_motifs_by_organism import (  # noqa: E402
    DEFAULT_TAXONOMY, HUMAN_TAX_ID, matched_chains, pdb_id_from_tid,
)

DEFAULT_LOOKUP = Path(
    "/Users/mvh_mac/Documents/master_thesis_2026/folddisco_databases/"
    "index/pdb_folddisco.lookup"
)


def indexed_entries(lookup: Path) -> set[str]:
    """Every PDB entry the search could have returned."""
    out: set[str] = set()
    with lookup.open() as fh:
        for line in fh:
            parts = line.split("\t")
            if len(parts) > 1:
                pdb = pdb_id_from_tid(parts[1])
                if pdb:
                    out.add(pdb)
    return out


def human_chain_map(taxonomy: Path, keep: set[str]) -> dict[tuple[str, str], bool]:
    """(pdb, chain) -> is human, for chains of entries in ``keep``."""
    out: dict[tuple[str, str], bool] = {}
    with taxonomy.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or parts[0] == "PDB":
                continue
            pdb, chain, tax = parts[0].lower(), parts[1], parts[2]
            if pdb not in keep:
                continue
            key = (pdb, chain)
            # An entry can list a chain more than once (one row per mapped
            # UniProt segment). Human anywhere on the chain makes it human.
            out[key] = out.get(key, False) or (tax == HUMAN_TAX_ID)
    return out


def called_motifs(results_dir: Path, threshold: float) -> set[str]:
    """Query strings whose motif scored at or above the threshold."""
    scored = results_dir / "idf_scored.tsv"
    if not scored.is_file():
        return set()
    out = set()
    for r in csv.DictReader(scored.open(), delimiter="\t"):
        try:
            if float(r["idf_z"]) >= threshold:
                out.add(r["query_residues"])
        except (ValueError, KeyError):
            continue
    return out


def binomial_z(k: int, n: int, p: float) -> float | None:
    """Normal-approximation z for k successes of n against rate p."""
    if n == 0 or not 0 < p < 1:
        return None
    return (k - n * p) / math.sqrt(n * p * (1 - p))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--threshold", type=float, default=None,
                    help="restrict to hits of motifs scoring at or above this "
                         "z; omit to use every hit")
    ap.add_argument("--taxonomy", type=Path, default=DEFAULT_TAXONOMY)
    ap.add_argument("--lookup", type=Path, default=DEFAULT_LOOKUP)
    ap.add_argument("-o", "--out", type=Path, default=None)
    args = ap.parse_args(argv)

    for path in (args.taxonomy, args.lookup):
        if not path.is_file():
            print(f"error: not found: {path}", file=sys.stderr)
            return 2

    dirs = sorted(d for d in args.run_dir.iterdir()
                  if (d / "folddisco_hits.tsv").is_file())
    if not dirs:
        print(f"error: no searched members under {args.run_dir}", file=sys.stderr)
        return 2

    print(f"reading the indexed universe from {args.lookup.name}", flush=True)
    universe = indexed_entries(args.lookup)
    print(f"  {len(universe):,} entries in the index", flush=True)
    chains = human_chain_map(args.taxonomy, universe)
    bg_total = len(chains)
    bg_human = sum(chains.values())
    bg_entries = len({p for p, _ in chains})
    print(f"  {bg_total:,} SIFTS-mapped chains in {bg_entries:,} of those "
          f"entries; {bg_human:,} human ({100 * bg_human / bg_total:.1f}%)",
          flush=True)

    row_counts = Counter()
    seen: set[tuple[str, str]] = set()
    per_structure = []
    unmapped_entries: set[str] = set()

    for d in dirs:
        keep = called_motifs(d, args.threshold) if args.threshold is not None else None
        local = Counter()
        local_chains: set[tuple[str, str]] = set()
        for r in csv.DictReader((d / "folddisco_hits.tsv").open(), delimiter="\t"):
            if keep is not None and r.get("query_residues", "") not in keep:
                continue
            pdb = pdb_id_from_tid(r.get("tid", ""))
            if pdb is None:
                row_counts["unknown"] += 1
                local["unknown"] += 1
                continue
            hit_chains = matched_chains(r.get("matching_residues", ""))
            resolved = [(pdb, c) for c in hit_chains if (pdb, c) in chains]
            if not resolved:
                row_counts["unknown"] += 1
                local["unknown"] += 1
                unmapped_entries.add(pdb)
                continue
            # A hit row counts once, as human if any matched chain is human.
            label = "human" if any(chains[k] for k in resolved) else "other"
            row_counts[label] += 1
            local[label] += 1
            seen.update(resolved)
            local_chains.update(resolved)
        known = local["human"] + local["other"]
        per_structure.append({
            "instance": d.name,
            "hit_rows": sum(local.values()),
            "human_rows": local["human"], "other_rows": local["other"],
            "unknown_rows": local["unknown"],
            "human_row_fraction": round(local["human"] / known, 4) if known else None,
            "distinct_chains": len(local_chains),
        })

    distinct = Counter()
    for key in seen:
        distinct["human" if chains[key] else "other"] += 1

    def block(counts: Counter, label: str) -> dict:
        known = counts["human"] + counts["other"]
        frac = counts["human"] / known if known else None
        bg = bg_human / bg_total
        return {
            "unit": label,
            "human": counts["human"], "other": counts["other"],
            "unknown": counts["unknown"],
            "human_fraction_of_known": round(frac, 4) if frac else None,
            "background_human_fraction": round(bg, 4),
            "enrichment_over_background": round(frac / bg, 3) if frac else None,
            "binomial_z_vs_background": (
                round(binomial_z(counts["human"], known, bg), 1)
                if known else None
            ),
        }

    summary = {
        "run": str(args.run_dir),
        "threshold": args.threshold,
        "structures": len(dirs),
        "index_entries": len(universe),
        "background_chains_mapped": bg_total,
        "background_human_chains": bg_human,
        "background_entries_mapped": bg_entries,
        "background_entries_unmapped": len(universe) - bg_entries,
        "rows": block(row_counts, "hit rows"),
        "distinct": block(
            Counter({**distinct, "unknown": 0}), "distinct matched chains"),
        "distinct_chains_touched": len(seen),
        "entries_hit_but_unmapped": len(unmapped_entries),
    }

    out = args.out or (args.run_dir / "hit_organism_stats.json")
    out.write_text(json.dumps(
        {"summary": summary, "per_structure": per_structure}, indent=2))

    print(f"\nwrote {out}")
    for key in ("rows", "distinct"):
        b = summary[key]
        print(f"  {b['unit']:<24} human {b['human']:>7,}  other {b['other']:>7,}"
              f"  unknown {b['unknown']:>7,}")
        print(f"  {'':<24} human share {b['human_fraction_of_known']}"
              f"  vs background {b['background_human_fraction']}"
              f"  = {b['enrichment_over_background']}x"
              f"  (z {b['binomial_z_vs_background']})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
