"""Write each structure's top-scoring motifs with the organism each hit came from.

For every results directory it selects the highest-scoring motifs -- the top
``--fraction`` by IDF z-score -- and writes their hits to a separate file,
one row per hit, annotated with the source organism of the structure that
was hit.

Why hit-level rather than motif-level
-------------------------------------
An organism is a property of the structure a motif matched, not of the motif
itself, and a single motif typically matches structures from many organisms.
Collapsing that to one organism per motif would discard the distribution
that makes the column useful. Each row is therefore one (motif, hit) pair,
carrying the motif's rank and z-score so motif-level filtering is still a
one-liner.

How the organism is resolved
----------------------------
``matching_residues`` names the chains that actually matched, so the lookup
is per chain -- ``(pdb id, chain)`` against the SIFTS mapping -- rather than
per entry. In a complex of several proteins the matched chain is the one
whose organism is meaningful; the entry as a whole may contain several.

Three honest outcomes are distinguished in ``organism_source``:

``chain``
    the matched chains resolved directly.
``entry``
    the matched chains were not in SIFTS but every mapped chain of the entry
    agrees on one organism, so that is reported.
``unknown``
    no mapping. SIFTS covers only chains it could map to UniProt, so
    synthetic constructs and unmapped chains are simply absent. **A missing
    entry means "not known", never "no organism"** -- treating it as the
    latter would quietly dilute any fraction computed over these rows.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

DEFAULT_TAXONOMY = Path(
    "/Users/mvh_mac/documents_local_only/master_thesis_data_files/sifts/"
    "pdb_chain_taxonomy.tsv"
)

HUMAN_TAX_ID = "9606"

# A FoldDisco tid for the whole-PDB index looks like `pdb/mw/pdb2mwl.ent.gz`.
_PDB_FROM_TID = re.compile(r"pdb([0-9a-z]{4})\.ent", re.IGNORECASE)


def pdb_id_from_tid(tid: str) -> str | None:
    """Extract the 4-character entry ID a hit path refers to."""
    match = _PDB_FROM_TID.search(tid)
    if match:
        return match.group(1).lower()
    # Fall back to a bare 4-character basename, as other id schemes produce.
    stem = tid.split("/")[-1].split(".")[0].lower()
    return stem if re.fullmatch(r"[0-9][0-9a-z]{3}", stem) else None


def matched_chains(matching_residues: str) -> list[str]:
    """Chains named in a hit's matched tokens, in order, deduplicated.

    Read as a SET of tokens rather than positionally: this build omits some
    unmatched positions instead of padding them, so the field is not reliably
    aligned to the query (see docs/FOLDDISCO_CONTRACT.md).
    """
    seen: list[str] = []
    for token in (matching_residues or "").split(","):
        token = token.strip()
        if not token or token == "_":
            continue
        chain = token[0]
        if chain not in seen:
            seen.append(chain)
    return seen


def load_taxonomy(path: Path, wanted: set[str]) -> tuple[dict, dict]:
    """Return ``(by_chain, by_entry)`` restricted to the entries we need."""
    by_chain: dict[tuple[str, str], list[tuple[str, str]]] = defaultdict(list)
    by_entry: dict[str, set[tuple[str, str]]] = defaultdict(set)
    with path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4 or parts[0] == "PDB":
                continue
            pdb, chain, tax_id, name = parts[0].lower(), parts[1], parts[2], parts[3]
            if pdb not in wanted:
                continue
            by_chain[(pdb, chain)].append((tax_id, name))
            by_entry[pdb].add((tax_id, name))
    return by_chain, by_entry


def resolve_organism(pdb: str | None, chains: list[str], by_chain, by_entry):
    """-> (tax_ids, names, source). Never guesses when the mapping is absent."""
    if pdb is None:
        return "", "", "unknown"
    found: list[tuple[str, str]] = []
    for chain in chains:
        for pair in by_chain.get((pdb, chain), []):
            if pair not in found:
                found.append(pair)
    if found:
        return ("|".join(t for t, _ in found),
                "|".join(n for _, n in found), "chain")
    entry = by_entry.get(pdb) or set()
    # Only usable as a fallback when the entry is unambiguous.
    if len(entry) == 1:
        tax_id, name = next(iter(entry))
        return tax_id, name, "entry"
    return "", "", "unknown"


def top_motifs(scored: Path, fraction: float) -> dict[str, tuple[int, float]]:
    """Query string -> (rank, z) for the highest-scoring fraction of motifs."""
    rows = []
    for r in csv.DictReader(scored.open(), delimiter="\t"):
        try:
            rows.append((float(r["idf_z"]), r["query_residues"]))
        except (ValueError, KeyError):
            continue
    rows.sort(reverse=True)
    keep = max(1, round(fraction * len(rows))) if rows else 0
    return {q: (i + 1, z) for i, (z, q) in enumerate(rows[:keep])}


def process(
    results_dir: Path, fraction: float, by_chain, by_entry, out_name: str
) -> tuple[int, int, int] | None:
    scored = results_dir / "idf_scored.tsv"
    hits = results_dir / "folddisco_hits.tsv"
    if not scored.is_file() or not hits.is_file():
        return None

    selected = top_motifs(scored, fraction)
    if not selected:
        return None

    written = unknown = 0
    out = results_dir / out_name
    with out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "motif_rank", "idf_z", "query_residues", "motif_size",
            "patch_id", "hit_tid", "hit_pdb_id", "hit_chains",
            "idf", "rmsd", "node_count", "matching_residues",
            "organism_tax_ids", "organism_names", "organism_source", "is_human",
        ])
        for r in csv.DictReader(hits.open(), delimiter="\t"):
            q = r.get("query_residues", "")
            if q not in selected:
                continue
            rank, z = selected[q]
            pdb = pdb_id_from_tid(r.get("tid", ""))
            chains = matched_chains(r.get("matching_residues", ""))
            tax_ids, names, source = resolve_organism(pdb, chains, by_chain, by_entry)
            if source == "unknown":
                unknown += 1
            w.writerow([
                rank, z, q, len([t for t in q.split(",") if t]),
                r.get("patch_id", ""), r.get("tid", ""), pdb or "",
                "|".join(chains), r.get("idf", ""), r.get("rmsd", ""),
                r.get("node_count", ""), r.get("matching_residues", ""),
                tax_ids, names, source,
                "1" if HUMAN_TAX_ID in tax_ids.split("|") else "0",
            ])
            written += 1
    return len(selected), written, unknown


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("results_dirs", nargs="+", type=Path,
                    help="per-structure results directories")
    ap.add_argument("--fraction", type=float, default=0.05,
                    help="top fraction of motifs by z-score [0.05]")
    ap.add_argument("--taxonomy", type=Path, default=DEFAULT_TAXONOMY,
                    help="SIFTS pdb_chain_taxonomy.tsv")
    ap.add_argument("--out-name", default="top_motifs_organism.tsv")
    args = ap.parse_args(argv)

    dirs = [d for d in args.results_dirs if (d / "folddisco_hits.tsv").is_file()]
    if not dirs:
        print("error: no results directories with folddisco_hits.tsv", file=sys.stderr)
        return 2
    if not args.taxonomy.is_file():
        print(f"error: taxonomy file not found: {args.taxonomy}\n"
              "  Download pdb_chain_taxonomy.tsv from the SIFTS flatfiles, or "
              "pass --taxonomy.", file=sys.stderr)
        return 2

    # Only load the entries actually hit; the full file is 1.1M rows.
    wanted: set[str] = set()
    for d in dirs:
        for r in csv.DictReader((d / "folddisco_hits.tsv").open(), delimiter="\t"):
            pdb = pdb_id_from_tid(r.get("tid", ""))
            if pdb:
                wanted.add(pdb)
    print(f"{len(dirs)} structures; {len(wanted):,} distinct entries hit", flush=True)
    by_chain, by_entry = load_taxonomy(args.taxonomy, wanted)
    print(f"taxonomy: {len(by_chain):,} chain mappings for "
          f"{len(by_entry):,} entries", flush=True)

    tot_m = tot_r = tot_u = 0
    for d in sorted(dirs):
        res = process(d, args.fraction, by_chain, by_entry, args.out_name)
        if res is None:
            print(f"  {d.name}: skipped (no scored motifs)")
            continue
        motifs, rows, unknown = res
        tot_m += motifs; tot_r += rows; tot_u += unknown
        print(f"  {d.name:<10} top {motifs:>4} motifs -> {rows:>6} hit rows"
              f"  ({100 * unknown / rows if rows else 0:4.1f}% organism unknown)")
    print(f"\ntotal: {tot_m:,} motifs, {tot_r:,} rows, "
          f"{100 * tot_u / tot_r if tot_r else 0:.1f}% unknown organism")
    return 0


if __name__ == "__main__":
    sys.exit(main())
