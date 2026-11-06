"""Build ``background_dataset``: a deliberately mixed, non-redundant background.

What this set is for
--------------------
The three sets built by ``build_investigation_db.py`` are each selected
*for* something -- an interface, a zinc site, an EC number. This one is
selected for nothing. It exists to be the comparison: a spread of proteins
wide enough that a statistic computed on it describes proteins in general
rather than one family, one organism, one size or one experimental method.

Composition is therefore the product, not a by-product, and every axis of it
is recorded in the manifest so the spread can be checked rather than trusted.

How the spread is produced
--------------------------
**Non-redundancy.** RCSB is asked for one representative per **30%
sequence-identity cluster**. Raw PDB hits would be dominated by whatever the
field deposits hundreds of times -- lysozyme, trypsin, a few kinases -- and a
background made of those is not a background.

**Random, not best-first.** ``build_investigation_db.py`` takes cluster
representatives in ascending resolution, which is right when the goal is the
best members of a category. It is wrong here: sorting by resolution selects
small, well-ordered, heavily-studied X-ray structures, which is a bias on
three of the axes this set is meant to span. Candidates are shuffled under a
recorded seed instead, so the run is reproducible without being ordered.

**Size strata.** Random sampling alone would reproduce the PDB's own size
distribution, which is concentrated around 100-350 residues, and the large
end would be nearly absent. Candidates are therefore accepted into size bins
with per-bin quotas. A bin that cannot be filled redistributes its shortfall
and the manifest records which bin fell short, so a gap is visible instead of
being quietly smoothed over.

**Organism and function caps.** A per-organism cap stops *E. coli*, human and
*S. cerevisiae* from crowding out everything else. Function is not filtered
on -- there is no clean, complete functional label across the whole PDB -- but
the EC-class distribution of what was selected is reported, so the spread can
be inspected.

**Two sources.** Experimental structures come from the PDB; the rest are
AlphaFold DB predictions. They are not interchangeable, which is the point of
including both, and ``source`` distinguishes them in every row.

What "high quality" means, which differs by source
--------------------------------------------------
This is the part that cannot be collapsed into one column, and the code does
not try to:

* **X-ray and cryo-EM** -- resolution strictly below ``--resolution``.
* **NMR** -- *no resolution exists*. An NMR entry is admitted on the basis
  that it is NMR, not on a resolution it does not have. ``resolution_angstrom``
  is empty for those rows; it is not zero and not the threshold.
* **AlphaFold** -- no experiment happened at all. The quality filter is mean
  **pLDDT** at or above ``--min-plddt``, which is the model's own confidence
  in itself. It is **not** a resolution, it is not comparable to one, and it
  is reported in its own column.

Anything that mixes these three columns into a single "quality" number is
mixing a measurement, an absence of a measurement, and a prediction's
self-assessment. They are kept apart here and must stay apart downstream.

What a row does NOT mean
------------------------
Nothing in this set is labelled. There are no epitopes, no interfaces, no
sites -- that is what makes it a background. An AlphaFold row is a
*prediction*: its coordinates were never observed, and a high pLDDT is the
model's confidence, not evidence.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import random
import sys
import time
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx_io

REPO = Path(__file__).resolve().parents[1]

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL_URL = "https://data.rcsb.org/graphql"
RCSB_FILE_URL = "https://files.rcsb.org/download/{pdb}.cif.gz"
AFDB_API_URL = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"

USER_AGENT = "ibc-search background-dataset builder (academic structure analysis)"
"""Identifies this client to RCSB and AFDB.

Both reject requests that do not identify themselves. It names the tool
honestly rather than imitating a browser.
"""

WATER = frozenset({"HOH", "DOD", "WAT", "H2O"})

ALLOWED_METHODS = ("X-RAY DIFFRACTION", "ELECTRON MICROSCOPY", "SOLUTION NMR")
NO_RESOLUTION_METHODS = frozenset({"SOLUTION NMR", "SOLID-STATE NMR"})

SIZE_BINS: tuple[tuple[int, int], ...] = (
    (50, 100), (100, 200), (200, 350), (350, 550),
    (550, 800), (800, 1200), (1200, 2300),
)
"""Residue-count strata, half-open except the last, which includes 2300.

Seven bins spanning the requested 50-2300 range. The widths grow because the
PDB thins out sharply with size: equal-width bins would leave the top three
permanently unfillable and the bottom one saturated.
"""


# ---------------------------------------------------------------------------
# HTTP
# ---------------------------------------------------------------------------

def _get(url: str, timeout: int = 120) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as fh:
        return fh.read()


def _post_json(url: str, body: dict, timeout: int = 180) -> dict:
    req = urllib.request.Request(
        url, data=json.dumps(body).encode(),
        headers={"Content-Type": "application/json", "User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as fh:
        return json.load(fh)


def _retry(fn, *, tries: int = 3, note: str = ""):
    """Call ``fn``, backing off on failure. Returns None when all tries fail.

    Both services throttle bursts, so a failure is usually temporary and a
    retry is the right response; a permanent failure returns None and is
    recorded rather than aborting a run of hundreds of downloads.
    """
    for attempt in range(tries):
        try:
            return fn()
        except Exception as exc:  # noqa: BLE001 - recorded, not swallowed
            if attempt == tries - 1:
                if note:
                    print(f"    {note}: {type(exc).__name__}: {exc}"[:150])
                return None
            time.sleep(1.5 * (attempt + 1))
    return None


# ---------------------------------------------------------------------------
# Candidate search (PDB)
# ---------------------------------------------------------------------------

def _term(attr: str, op: str, val=None) -> dict:
    p = {"attribute": attr, "operator": op}
    if val is not None:
        p["value"] = val
    return {"type": "terminal", "service": "text", "parameters": p}


def search_cluster_representatives(resolution: float, identity: int,
                                   lo: int, hi: int) -> list[str]:
    """One polymer-entity id per sequence-identity cluster.

    The quality node is an OR, not an AND: an NMR entry carries no
    ``resolution_combined`` at all, so requiring a resolution would silently
    exclude every NMR structure while appearing to be a quality filter. The
    OR admits NMR on the basis of its method and holds everything else to the
    resolution limit.
    """
    quality = {"type": "group", "logical_operator": "or", "nodes": [
        _term("rcsb_entry_info.resolution_combined", "less", resolution),
        _term("exptl.method", "exact_match", "SOLUTION NMR"),
    ]}
    nodes = [
        _term("exptl.method", "in", list(ALLOWED_METHODS)),
        quality,
        _term("entity_poly.rcsb_entity_polymer_type", "exact_match", "Protein"),
        _term("rcsb_entry_info.deposited_polymer_monomer_count", "range",
              {"from": lo, "to": hi}),
    ]
    body = {
        "query": {"type": "group", "logical_operator": "and", "nodes": nodes},
        "return_type": "polymer_entity",
        "request_options": {
            "group_by": {"aggregation_method": "sequence_identity",
                         "similarity_cutoff": identity},
            "group_by_return_type": "representatives",
            "return_all_hits": True,
            "results_verbosity": "compact",
        },
    }
    payload = _post_json(SEARCH_URL, body)
    return list(payload.get("result_set", []))


ENTITY_QUERY = """query($ids:[String!]!){
  polymer_entities(entity_ids:$ids){
    rcsb_id
    rcsb_polymer_entity{pdbx_description rcsb_ec_lineage{id}}
    rcsb_polymer_entity_container_identifiers{uniprot_ids}
    rcsb_entity_source_organism{ncbi_scientific_name ncbi_taxonomy_id}
    entry{
      rcsb_id
      exptl{method}
      rcsb_entry_info{resolution_combined deposited_polymer_monomer_count
                      polymer_entity_count_protein}
      polymer_entities{
        rcsb_polymer_entity_container_identifiers{uniprot_ids}
      }
    }
  }
}"""


def entity_metadata(entity_ids: list[str]) -> list[dict]:
    """Method, resolution, size, organism, EC and UniProt, in one request.

    Batched because the alternative is one REST call per candidate, and the
    selection examines far more candidates than it keeps.
    """
    payload = _retry(
        lambda: _post_json(GRAPHQL_URL,
                           {"query": ENTITY_QUERY, "variables": {"ids": entity_ids}}),
        note="graphql")
    if not payload or "data" not in payload:
        return []
    return [e for e in (payload["data"].get("polymer_entities") or []) if e]


def describe_entity(entity: dict) -> dict | None:
    """Flatten one GraphQL entity into the fields selection needs."""
    entry = entity.get("entry") or {}
    info = entry.get("rcsb_entry_info") or {}
    pdb_id = entry.get("rcsb_id")
    if not pdb_id:
        return None
    methods = [m["method"] for m in (entry.get("exptl") or []) if m.get("method")]
    res_list = info.get("resolution_combined") or []
    resolution = float(res_list[0]) if res_list else None
    org = (entity.get("rcsb_entity_source_organism") or [{}])[0] or {}
    poly = entity.get("rcsb_polymer_entity") or {}
    ec = poly.get("rcsb_ec_lineage") or []
    ident = entity.get("rcsb_polymer_entity_container_identifiers") or {}
    # Every accession in the ENTRY, not only the clustered entity's. RCSB
    # clusters an entity but this script writes the whole entry to disk, so
    # redundancy has to be judged on the file: four thrombin-hirudin
    # complexes whose hirudin peptides fall in different 30% clusters are
    # four near-identical structures, and the entity-level accession alone
    # does not catch them.
    entry_accessions = {
        u
        for ent in (entry.get("polymer_entities") or [])
        for u in (((ent or {}).get("rcsb_polymer_entity_container_identifiers")
                   or {}).get("uniprot_ids") or [])
        if u
    }
    return {
        "entry_uniprot_ids": sorted(entry_accessions),
        "pdb_id": pdb_id,
        "entity_id": entity.get("rcsb_id"),
        "methods": methods,
        "resolution": resolution,
        "n_residues_reported": info.get("deposited_polymer_monomer_count"),
        "n_protein_entities": info.get("polymer_entity_count_protein"),
        "organism": org.get("ncbi_scientific_name") or "",
        "tax_id": org.get("ncbi_taxonomy_id"),
        "description": (poly.get("pdbx_description") or "").strip(),
        "ec_class": (ec[0]["id"] if ec else ""),
        "uniprot_ids": [u for u in (ident.get("uniprot_ids") or []) if u],
    }


def passes_quality(meta: dict, resolution_max: float) -> tuple[bool, str]:
    """The quality rule, stated once, applied to one candidate.

    A multi-method entry (X-ray plus neutron, say) is admitted if any of its
    methods is one this set accepts; the full list is kept in the row rather
    than reduced to the one that let it in.
    """
    if not any(m in ALLOWED_METHODS for m in meta["methods"]):
        return False, f"method {','.join(meta['methods']) or 'unknown'}"
    if meta["resolution"] is None:
        if not any(m in NO_RESOLUTION_METHODS for m in meta["methods"]):
            # No resolution and not an NMR entry: unexplained, so not
            # admitted. Treating a missing value as passing would let
            # unresolved entries in through the gap NMR opens.
            return False, "no resolution and not NMR"
        return True, ""
    if meta["resolution"] >= resolution_max:
        return False, f"resolution {meta['resolution']}"
    return True, ""


# ---------------------------------------------------------------------------
# Candidate pool (AlphaFold)
# ---------------------------------------------------------------------------

def afdb_candidates(lookup: Path, lo: int, hi: int,
                    min_plddt: float) -> list[tuple[str, int, float]]:
    """``(uniprot accession, residues, mean pLDDT)`` from the local AFDB lookup.

    Read from the FoldDisco index's own lookup table rather than queried over
    the network: it already carries the residue count and the mean pLDDT for
    every Swiss-Prot model, so the whole pool can be size- and
    confidence-filtered locally and only the chosen few are ever downloaded.
    """
    out: list[tuple[str, int, float]] = []
    with lookup.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            name = parts[1]
            try:
                n_res, plddt = int(parts[2]), float(parts[3])
            except ValueError:
                continue
            if not (lo <= n_res <= hi) or plddt < min_plddt:
                continue
            # AF-<accession>-F1-model_v4
            bits = name.split("-")
            if len(bits) < 2:
                continue
            out.append((bits[1], n_res, plddt))
    return out


def afdb_metadata(accession: str) -> dict | None:
    raw = _retry(lambda: _get(AFDB_API_URL.format(acc=accession), timeout=60),
                 tries=2, note=f"afdb api {accession}")
    if raw is None:
        return None
    try:
        records = json.loads(raw)
    except json.JSONDecodeError:
        return None
    if not records:
        return None
    r = records[0]
    return {
        "accession": accession,
        "cif_url": r.get("cifUrl"),
        "organism": r.get("organismScientificName") or "",
        "tax_id": r.get("taxId"),
        "description": (r.get("uniprotDescription") or "").strip(),
        "plddt": r.get("globalMetricValue"),
        "version": r.get("latestVersion"),
    }


# ---------------------------------------------------------------------------
# Selection
# ---------------------------------------------------------------------------

def bin_of(n_residues: int) -> int | None:
    for i, (lo, hi) in enumerate(SIZE_BINS):
        if lo <= n_residues < hi or (i == len(SIZE_BINS) - 1 and n_residues == hi):
            return i
    return None


def bin_quotas(total: int) -> list[int]:
    """Split a target evenly across the size bins, remainder to the small end.

    Even quotas, not proportional ones: a proportional split would reproduce
    the PDB's own size distribution, which is the bias the strata exist to
    counter. The remainder goes to the smaller bins because those are the
    ones that can actually supply it.
    """
    base, extra = divmod(total, len(SIZE_BINS))
    return [base + (1 if i < extra else 0) for i in range(len(SIZE_BINS))]


class Selector:
    """Accepts candidates into size bins under an organism cap.

    Holds the whole acceptance policy, so that the PDB and AlphaFold passes
    cannot drift into applying different rules to their candidates.
    """

    def __init__(self, quotas: list[int], organism_cap: int):
        self.quotas = list(quotas)
        self.organism_cap = organism_cap
        self.taken: list[dict] = []
        self.per_bin: Counter[int] = Counter()
        self.per_organism: Counter[str] = Counter()
        self.rejected: Counter[str] = Counter()

    @property
    def n(self) -> int:
        return len(self.taken)

    def offer(self, row: dict) -> bool:
        b = row["size_bin"]
        if self.per_bin[b] >= self.quotas[b]:
            self.rejected[f"bin {SIZE_BINS[b][0]}-{SIZE_BINS[b][1]} full"] += 1
            return False
        org = row.get("organism") or "unknown"
        if org != "unknown" and self.per_organism[org] >= self.organism_cap:
            self.rejected["organism cap"] += 1
            return False
        self.taken.append(row)
        self.per_bin[b] += 1
        self.per_organism[org] += 1
        return True

    def shortfall(self) -> dict[str, int]:
        """Which bins are below quota, and by how much.

        Reported rather than hidden: a bin that stays short means neither
        source had enough non-redundant structures of that size under the
        quality rule, which is a fact about the data and not a bug to paper
        over.
        """
        return {f"{SIZE_BINS[i][0]}-{SIZE_BINS[i][1]}":
                self.quotas[i] - self.per_bin[i]
                for i in range(len(SIZE_BINS))
                if self.per_bin[i] < self.quotas[i]}

    def relax(self, spare: int) -> None:
        """Let the bins that CAN still supply absorb an unfillable one's quota.

        Only ever called at the end, once both sources have been walked and
        the set is still short of its total. Calling it between passes would
        inflate quotas while the second source still had candidates for the
        bin that looked short, which is not the same thing as the bin being
        unfillable.
        """
        for i in range(len(SIZE_BINS)):
            self.quotas[i] += spare


# ---------------------------------------------------------------------------
# Structure files
# ---------------------------------------------------------------------------

def download_to(url: str, dest: Path, tries: int = 3) -> bool:
    if dest.is_file() and dest.stat().st_size > 0:
        return True
    data = _retry(lambda: _get(url, timeout=180), tries=tries,
                  note=f"download {dest.name}")
    if not data:
        return False
    tmp = dest.with_suffix(dest.suffix + ".part")
    if dest.name.endswith(".gz") and not data[:2] == b"\x1f\x8b":
        tmp.write_bytes(gzip.compress(data))
    else:
        tmp.write_bytes(data)
    tmp.replace(dest)
    return True


def count_residues(path: Path) -> int:
    """Amino-acid residues in model 1, counted from the coordinates.

    The authoritative number for this set. RCSB's
    ``deposited_polymer_monomer_count`` is used to shortlist candidates
    cheaply, but it counts the deposited polymer rather than what a parser
    sees, so the size bin and the 50-2300 bound are decided here, on the file
    that was actually written.
    """
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt") as fh:
        cif = pdbx_io.CIFFile.read(fh)
    arr = pdbx_io.get_structure(cif, model=1)
    arr = arr[(arr.element != "H") & ~np.isin(arr.res_name, list(WATER))]
    arr = arr[struc.filter_amino_acids(arr)]
    if arr.array_length() == 0:
        return 0
    return len({(str(c), int(s)) for c, s in zip(arr.chain_id, arr.res_id)})


def count_chains(path: Path) -> int:
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt") as fh:
        cif = pdbx_io.CIFFile.read(fh)
    arr = pdbx_io.get_structure(cif, model=1)
    arr = arr[struc.filter_amino_acids(arr)]
    return len(set(arr.chain_id.tolist())) if arr.array_length() else 0


ROW_HEADER = [
    "structure_id", "source", "method", "resolution_angstrom", "mean_plddt",
    "afdb_model_version", "n_residues", "size_bin", "n_chains", "organism",
    "tax_id", "ec_class", "uniprot_ids", "description", "file",
]


def _row(entry: dict) -> list:
    return [entry.get(k, "") for k in ROW_HEADER]


# ---------------------------------------------------------------------------
# Passes
# ---------------------------------------------------------------------------

def collect_pdb(sel: Selector, target: int, args, struct_dir: Path,
                stats: dict) -> None:
    print(f"searching RCSB for non-redundant cluster representatives "
          f"(<{args.resolution} A, {args.min_residues}-{args.max_residues} "
          f"residues, {args.identity}% identity)...", flush=True)
    reps = search_cluster_representatives(
        args.resolution, args.identity, args.min_residues, args.max_residues)
    stats["pdb_cluster_representatives"] = len(reps)
    print(f"  {len(reps):,} cluster representatives", flush=True)

    rng = random.Random(args.seed)
    rng.shuffle(reps)
    stats["pdb_candidates_examined"] = 0
    seen_pdb: set[str] = set()
    seen_uniprot: set[str] = set()

    for start in range(0, len(reps), args.batch):
        if sel.n >= target:
            break
        batch = reps[start:start + args.batch]
        stats["pdb_candidates_examined"] += len(batch)
        for entity in entity_metadata(batch):
            if sel.n >= target:
                break
            meta = describe_entity(entity)
            if meta is None:
                continue
            # One entry can contribute several entities; the file is the
            # entry, so a second entity of an entry already taken is a
            # duplicate structure, not a second data point.
            if meta["pdb_id"] in seen_pdb:
                sel.rejected["duplicate entry"] += 1
                continue
            # Any shared accession means the same protein is already in the
            # set, in some complex or other. 30% clustering is applied by
            # RCSB to one entity; this is the check on what gets written.
            if seen_uniprot.intersection(meta["entry_uniprot_ids"]):
                sel.rejected["uniprot already selected"] += 1
                continue
            ok, why = passes_quality(meta, args.resolution)
            if not ok:
                sel.rejected[why.split()[0] + " filter"] += 1
                continue
            reported = meta["n_residues_reported"] or 0
            b = bin_of(reported)
            if b is None:
                sel.rejected["size out of range"] += 1
                continue
            if sel.per_bin[b] >= sel.quotas[b]:
                sel.rejected[f"bin {SIZE_BINS[b][0]}-{SIZE_BINS[b][1]} full"] += 1
                continue

            pdb_id = meta["pdb_id"]
            dest = struct_dir / f"{pdb_id}.cif.gz"
            if not download_to(RCSB_FILE_URL.format(pdb=pdb_id), dest):
                sel.rejected["download failed"] += 1
                continue
            try:
                n_res = count_residues(dest)
                n_chains = count_chains(dest)
            except Exception as exc:  # noqa: BLE001
                sel.rejected["parse failed"] += 1
                print(f"    parse {pdb_id}: {type(exc).__name__}"[:120])
                dest.unlink(missing_ok=True)
                continue
            # Re-bin on the parsed count: the reported count got it onto the
            # shortlist, the file decides where it lands.
            b = bin_of(n_res)
            if b is None:
                sel.rejected["size out of range (parsed)"] += 1
                dest.unlink(missing_ok=True)
                continue

            row = {
                "structure_id": pdb_id,
                "source": "PDB",
                "method": "; ".join(meta["methods"]),
                "resolution_angstrom": ("" if meta["resolution"] is None
                                        else f"{meta['resolution']:.2f}"),
                "mean_plddt": "",
                "n_residues": n_res,
                "size_bin": b,
                "n_chains": n_chains,
                "organism": meta["organism"],
                "tax_id": meta["tax_id"] or "",
                "ec_class": meta["ec_class"],
                "uniprot_ids": ",".join(meta["entry_uniprot_ids"]),
                "description": meta["description"][:160],
                "file": dest.name,
            }
            if sel.offer(row):
                seen_pdb.add(pdb_id)
                seen_uniprot.update(meta["entry_uniprot_ids"])
                if sel.n % 25 == 0:
                    print(f"  {sel.n}/{target} experimental", flush=True)
            else:
                dest.unlink(missing_ok=True)


def collect_alphafold(sel: Selector, target: int, args, struct_dir: Path,
                      stats: dict, exclude_uniprot: set[str],
                      pool: list[tuple[str, int, float]], start: int = 0) -> int:
    """Fill up to ``target`` from the AlphaFold pool. Returns where it stopped.

    The caller owns the pool and the position in it, so a later top-up pass
    continues through candidates this one has not seen rather than
    re-examining -- and re-rejecting -- the ones it already declined.
    """
    for i in range(start, len(pool)):
        accession, n_reported, plddt = pool[i]
        if sel.n >= target:
            return i
        stats["afdb_candidates_examined"] = \
            stats.get("afdb_candidates_examined", 0) + 1
        # Cross-source redundancy: a prediction of a protein whose
        # experimental structure is already in the set is the same protein
        # twice, which is exactly what "no duplicates" rules out.
        if accession in exclude_uniprot:
            sel.rejected["uniprot already in PDB half"] += 1
            continue
        b = bin_of(n_reported)
        if b is None or sel.per_bin[b] >= sel.quotas[b]:
            sel.rejected["bin full (afdb)"] += 1
            continue
        meta = afdb_metadata(accession)
        if meta is None or not meta.get("cif_url"):
            sel.rejected["afdb metadata failed"] += 1
            continue
        # The local lookup is built from the v4 models; AFDB now serves v6,
        # and a re-prediction can move the confidence. The lookup is a cheap
        # pre-filter, but the floor has to be enforced on the pLDDT of the
        # model actually downloaded, or a row can carry a value below the
        # threshold the manifest claims. Same shortlist-then-verify shape as
        # the residue count.
        served_plddt = meta.get("plddt")
        if served_plddt is None or float(served_plddt) < args.min_plddt:
            sel.rejected["plddt below floor in served model"] += 1
            continue
        org = meta["organism"] or "unknown"
        if org != "unknown" and sel.per_organism[org] >= sel.organism_cap:
            sel.rejected["organism cap"] += 1
            continue

        dest = struct_dir / f"AF-{accession}.cif.gz"
        if not download_to(meta["cif_url"], dest):
            sel.rejected["download failed"] += 1
            continue
        try:
            n_res = count_residues(dest)
            n_chains = count_chains(dest)
        except Exception as exc:  # noqa: BLE001
            sel.rejected["parse failed"] += 1
            print(f"    parse AF-{accession}: {type(exc).__name__}"[:120])
            dest.unlink(missing_ok=True)
            continue
        b = bin_of(n_res)
        if b is None:
            sel.rejected["size out of range (parsed)"] += 1
            dest.unlink(missing_ok=True)
            continue

        row = {
            "structure_id": f"AF-{accession}",
            "source": "AlphaFold",
            "method": "ALPHAFOLD PREDICTION",
            # Deliberately empty: a prediction has no resolution, and writing
            # anything here would invite it being read as one.
            "resolution_angstrom": "",
            "mean_plddt": f"{float(served_plddt):.1f}",
            "afdb_model_version": meta.get("version") or "",
            "n_residues": n_res,
            "size_bin": b,
            "n_chains": n_chains,
            "organism": meta["organism"],
            "tax_id": meta["tax_id"] or "",
            "ec_class": "",
            "uniprot_ids": accession,
            "description": meta["description"][:160],
            "file": dest.name,
        }
        if sel.offer(row):
            if sel.n % 25 == 0:
                print(f"  {sel.n}/{target} with AlphaFold", flush=True)
        else:
            dest.unlink(missing_ok=True)
        time.sleep(args.afdb_delay)
    return len(pool)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--out-dir", type=Path,
                    default=REPO / "investigation_databases" / "background_dataset")
    ap.add_argument("--total", type=int, default=750,
                    help="structures in the finished set [750]")
    ap.add_argument("--alphafold", type=int, default=250,
                    help="minimum AlphaFold structures of that total [250]")
    ap.add_argument("--resolution", type=float, default=3.0,
                    help="strict upper bound, Angstroms; NMR is exempt [3.0]")
    ap.add_argument("--min-plddt", type=float, default=70.0,
                    help="mean pLDDT floor for AlphaFold models. NOT a "
                         "resolution: it is the model's confidence in itself "
                         "[70, AlphaFold's 'confident' cutoff]")
    ap.add_argument("--min-residues", type=int, default=50)
    ap.add_argument("--max-residues", type=int, default=2300)
    ap.add_argument("--identity", type=int, default=30,
                    help="sequence-identity clustering percentage [30]")
    ap.add_argument("--organism-cap", type=int, default=22,
                    help="most structures from any one organism [22, ~3%%]")
    ap.add_argument("--seed", type=int, default=20261030,
                    help="shuffling seed; recorded in the manifest")
    ap.add_argument("--batch", type=int, default=100,
                    help="candidates per GraphQL metadata request [100]")
    ap.add_argument("--afdb-delay", type=float, default=0.1,
                    help="seconds between AlphaFold requests [0.1]")
    ap.add_argument("--afdb-lookup", type=Path,
                    default=Path("~/documents_local_only/master_thesis_data_files/"
                                 "folddisco_databases/afdb_swissprot/"
                                 "afdb_swissprot_v4_folddisco.lookup").expanduser())
    args = ap.parse_args(argv)

    if args.alphafold > args.total:
        print("error: --alphafold exceeds --total", file=sys.stderr)
        return 2
    if not args.afdb_lookup.is_file():
        print(f"error: no AlphaFold lookup at {args.afdb_lookup}", file=sys.stderr)
        return 2

    out_dir = args.out_dir
    struct_dir = out_dir / "structures"
    struct_dir.mkdir(parents=True, exist_ok=True)

    stats: dict = {}
    started = time.time()
    n_experimental = args.total - args.alphafold
    sel = Selector(bin_quotas(args.total), args.organism_cap)

    # The experimental pass runs first and is capped at the experimental
    # share, so the AlphaFold minimum is guaranteed by construction rather
    # than by whatever is left over.
    collect_pdb(sel, n_experimental, args, struct_dir, stats)
    stats["pdb_selected"] = sel.n
    stats["pdb_bins_short_after_its_pass"] = sel.shortfall()
    print(f"experimental: {sel.n} selected", flush=True)

    print(f"reading the AlphaFold pool from {args.afdb_lookup}", flush=True)
    pool = afdb_candidates(args.afdb_lookup, args.min_residues,
                           args.max_residues, args.min_plddt)
    stats["afdb_pool_after_size_and_plddt"] = len(pool)
    print(f"  {len(pool):,} Swiss-Prot models at >= {args.min_plddt} pLDDT "
          f"and {args.min_residues}-{args.max_residues} residues", flush=True)
    random.Random(args.seed + 1).shuffle(pool)

    exclude = {u for row in sel.taken
               for u in (row["uniprot_ids"] or "").split(",") if u}
    stats["uniprot_accessions_in_pdb_half"] = len(exclude)
    at = collect_alphafold(sel, args.total, args, struct_dir, stats, exclude,
                           pool)
    stats["alphafold_selected"] = sel.n - stats["pdb_selected"]
    print(f"alphafold: {stats['alphafold_selected']} selected", flush=True)

    # Top-up. Reaching the requested total matters more than holding a bin
    # quota that the data cannot satisfy, so the strata are relaxed only
    # here, only by what is still missing, and the manifest records that it
    # happened and which bins were short when it did.
    stats["bins_short_before_topup"] = sel.shortfall()
    if sel.n < args.total:
        missing = args.total - sel.n
        print(f"topping up {missing} with relaxed size strata", flush=True)
        sel.relax(missing)
        at = collect_alphafold(sel, args.total, args, struct_dir, stats,
                               exclude, pool, start=at)
        stats["topped_up"] = sel.n - stats["pdb_selected"] - stats["alphafold_selected"]
        stats["alphafold_selected"] = sel.n - stats["pdb_selected"]
    stats["afdb_pool_position"] = at
    stats["bins_short_final"] = sel.shortfall()

    rows = sorted(sel.taken, key=lambda r: (r["source"], r["structure_id"]))
    tsv = out_dir / "structures.tsv"
    with tsv.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(ROW_HEADER)
        for row in rows:
            w.writerow(_row(row))

    # Remove any file left in structures/ that no row refers to: a run that
    # downloaded a candidate and then rejected it must not leave it behind
    # looking like a member of the set.
    keep = {r["file"] for r in rows}
    orphans = [p for p in struct_dir.iterdir()
               if p.is_file() and p.name not in keep]
    for p in orphans:
        p.unlink()

    by_size = Counter(r["size_bin"] for r in rows)
    manifest = {
        "built_utc": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        "category": "background_dataset",
        "purpose": ("an unlabelled, deliberately mixed background: spread "
                    "across size, organism, function, method and source, so "
                    "that a statistic computed on it describes proteins in "
                    "general rather than one family"),
        "sources": ["RCSB PDB (search API v2 + files.rcsb.org)",
                    "AlphaFold DB (api + files.alphafold.ebi.ac.uk)"],
        "random_seed": args.seed,
        "targets": {"total": args.total, "alphafold_minimum": args.alphafold,
                    "experimental": n_experimental},
        "filters": {
            "methods_accepted": list(ALLOWED_METHODS) + ["ALPHAFOLD PREDICTION"],
            "resolution_max_angstrom_exclusive": args.resolution,
            "resolution_caveat": (
                "applies to X-ray and cryo-EM only. NMR entries carry NO "
                "resolution and are admitted on method; their "
                "resolution_angstrom cell is EMPTY, not zero. AlphaFold "
                "models have no resolution at all and are filtered on mean "
                "pLDDT, which is a prediction's confidence in itself and is "
                "NOT comparable to a resolution."),
            "min_mean_plddt_alphafold": args.min_plddt,
            "plddt_basis": (
                "the local FoldDisco lookup is built from AFDB v4 and is used "
                "only to shortlist; the floor is enforced on the pLDDT of the "
                "model AFDB actually served, whose version is recorded per "
                "row in afdb_model_version"),
            "uniprot_non_redundancy": (
                "no UniProt accession appears in more than one row, across "
                "BOTH sources and counting every entity of a PDB entry rather "
                "than only the entity RCSB clustered: the file written is the "
                "whole entry, so redundancy is judged on that"),
            "residue_range": [args.min_residues, args.max_residues],
            "residue_count_basis": ("amino-acid residues parsed from model 1 "
                                    "of the written file, not the deposited "
                                    "monomer count used to shortlist"),
            "sequence_identity_cluster_percent": args.identity,
            "one_representative_per_cluster": True,
            "organism_cap": args.organism_cap,
            "size_bins": [list(b) for b in SIZE_BINS],
            "selection_order": ("shuffled under random_seed, NOT sorted by "
                                "resolution: sorting by quality biases the "
                                "set toward small, well-ordered, "
                                "heavily-studied X-ray structures"),
        },
        "counts": {
            "structures_written": len(rows),
            "by_source": dict(Counter(r["source"] for r in rows)),
            "by_method": dict(Counter(r["method"] for r in rows)),
            "by_size_bin": {f"{SIZE_BINS[i][0]}-{SIZE_BINS[i][1]}": by_size.get(i, 0)
                            for i in range(len(SIZE_BINS))},
            "distinct_organisms": len({r["organism"] for r in rows if r["organism"]}),
            "distinct_ec_classes": len({r["ec_class"] for r in rows if r["ec_class"]}),
            "residues_min": min((r["n_residues"] for r in rows), default=None),
            "residues_max": max((r["n_residues"] for r in rows), default=None),
            "residues_median": (int(np.median([r["n_residues"] for r in rows]))
                                if rows else None),
        },
        "top_organisms": dict(Counter(
            r["organism"] for r in rows if r["organism"]).most_common(20)),
        "ec_classes": dict(Counter(r["ec_class"] for r in rows if r["ec_class"])),
        "selection_stats": stats,
        "rejected_reasons": dict(sel.rejected),
        "orphan_files_removed": len(orphans),
        "seconds": round(time.time() - started, 1),
        "interpretation_notice": (
            "UNLABELLED by design: no epitopes, no interfaces, no sites. An "
            "AlphaFold row is a PREDICTION whose coordinates were never "
            "observed, and its pLDDT is the model's confidence, not "
            "evidence. Experimental and predicted rows must not be pooled "
            "into one 'quality' number."),
    }
    (out_dir / "build_manifest.json").write_text(json.dumps(manifest, indent=2))

    print(f"\nwrote {len(rows)} structures to {out_dir}")
    print(f"  by source   {manifest['counts']['by_source']}")
    print(f"  by size bin {manifest['counts']['by_size_bin']}")
    print(f"  organisms   {manifest['counts']['distinct_organisms']} distinct")
    print(f"  residues    {manifest['counts']['residues_min']}-"
          f"{manifest['counts']['residues_max']} "
          f"(median {manifest['counts']['residues_median']})")
    return 0 if len(rows) == args.total else 1


if __name__ == "__main__":
    sys.exit(main())
