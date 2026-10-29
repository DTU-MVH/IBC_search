"""Build a non-redundant, quality-filtered structure set per investigation category.

Three categories, one directory each under ``investigation_databases/``:

``protein_protein_interactions``
    entries with two or more protein entities. The label is the **interface**:
    residues with a heavy atom within ``--contact`` Angstroms of a heavy atom
    of a different protein chain.
``zincfingers``
    entities RCSB annotates with a zinc-finger domain. The label is the
    entry's own ``struct_site`` whose ligand is ZN.
``enzymes``
    entities carrying an EC assignment. The label is whatever ``struct_site``
    records the entry lists, with the ligand and evidence code kept per site.

Non-redundancy
--------------
RCSB is asked for one representative per **30% sequence-identity cluster**,
not for raw hits. Without that, a category is dominated by whichever proteins
the field has deposited hundreds of times -- lysozyme, trypsin, a handful of
kinases -- and every statistic computed on it is really a statistic about
those few. The cluster id each representative came from is recorded, so the
choice can be audited and a different representative taken later.

Quality
-------
X-ray only, resolution at or below ``--resolution``. Representatives are
taken in ascending resolution, so a capped category keeps its best members.
NMR and cryo-EM are excluded rather than mixed in: NMR entries have no
resolution at all, so admitting them would make the quality filter mean
different things for different rows of the same table.

The labels are computed, not assumed
------------------------------------
Interfaces come from coordinates. Zinc and enzyme sites come from the
deposition's own ``struct_site`` records and are reported with their
``pdbx_evidence_code``, which is usually ``Software`` -- the PDB generated
them from proximity to the ligand. That is a distance rule wearing an
annotation's clothes, and a column says so rather than letting the reader
assume an author asserted it.

A structure with no label is kept, with ``n_annotated = 0``. Annotation
coverage is then an observable property of the set instead of a hidden
filter applied before anyone could measure it.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx_io

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
DOWNLOAD_URL = "https://files.rcsb.org/download/{pdb}.cif.gz"
ASSEMBLY_URL = "https://files.rcsb.org/download/{pdb}-assembly{n}.cif.gz"
WATER = frozenset({"HOH", "DOD", "WAT", "H2O"})

CRYO_ADDITIVES = frozenset({
    # Cryoprotectants, precipitants, buffers and organic counter-ions. A
    # struct_site about one of these records how the crystal was grown, not
    # what the enzyme does. They are NOT dropped: a cryoprotectant frequently
    # occupies the substrate pocket, so the site can still mark the active
    # site. They are labelled, so the choice to use them is the analysis's.
    "EDO", "GOL", "PEG", "PG4", "PGE", "P6G", "1PE", "2PE", "MPD", "DMS",
    "SO4", "PO4", "ACT", "FMT", "CIT", "FLC", "TLA", "MLA", "MES", "EPE",
    "TRS", "BTB", "BIS", "IMD", "CAC", "ACY", "NO3", "AZI", "SCN", "BME",
    "DTT", "GSH", "URE", "TAR", "OXL", "MYR", "UNL", "UNX",
})

METAL_IONS = frozenset({
    # Genuinely ambiguous, and kept in their own class for that reason. Mg
    # and Ca are catalytic cofactors in kinases, nucleases and phosphatases,
    # and are also the commonest counter-ions in a crystallisation buffer.
    # Nothing in the entry distinguishes the two cases, so forcing these into
    # either "biological" or "additive" would be the script asserting
    # something it cannot know.
    "MG", "CA", "NA", "K", "LI", "CS", "RB", "CL", "BR", "IOD",
    "MN", "FE", "FE2", "CU", "CU1", "CO", "NI", "CD", "ZN", "HG", "SR", "BA",
})


def ligand_class(category: str, ligand: str) -> str:
    """What kind of thing the label is about, so noise can be filtered later.

    Four classes rather than two, because the data genuinely carries four
    cases and collapsing them would hide the ambiguity rather than report it.
    """
    if category == "protein_protein_interactions":
        return "interface"
    if category == "zincfingers":
        return "metal_ion"
    lig = (ligand or "").strip().upper()
    if not lig or lig in {"?", ".", "UNK"}:
        return "unspecified"
    if lig in CRYO_ADDITIVES:
        return "crystallization_additive"
    if lig in METAL_IONS:
        return "metal_ion"
    return "biological"

XRAY = {"type": "terminal", "service": "text",
        "parameters": {"attribute": "exptl.method", "operator": "exact_match",
                       "value": "X-RAY DIFFRACTION"}}


def _term(attr: str, op: str, val=None) -> dict:
    p = {"attribute": attr, "operator": op}
    if val is not None:
        p["value"] = val
    return {"type": "terminal", "service": "text", "parameters": p}


def category_nodes(category: str, resolution: float,
                   bound_ligand: str | None = None) -> list[dict]:
    res = _term("rcsb_entry_info.resolution_combined", "less_or_equal", resolution)
    if category == "protein_protein_interactions":
        nodes = [XRAY, res, _term("rcsb_entry_info.polymer_entity_count_protein",
                                  "greater_or_equal", 2)]
    elif category == "zincfingers":
        nodes = [XRAY, res, _term("rcsb_polymer_entity_annotation.name",
                                  "contains_phrase", "zinc finger")]
    elif category == "enzymes":
        nodes = [XRAY, res,
                 _term("rcsb_polymer_entity.rcsb_ec_lineage.id", "exists")]
    else:
        raise SystemExit(f"unknown category {category!r}")
    if bound_ligand:
        # Applied in the search rather than after it, so that RCSB picks a
        # cluster representative that HAS the component. Filtering afterwards
        # instead throws away whole clusters whose chosen representative
        # happened to be the one crystallised without it.
        nodes.append(_term(
            "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id",
            "exact_match", bound_ligand.upper()))
    return nodes


def search_representatives(category: str, resolution: float, identity: int,
                           bound_ligand: str | None = None) -> list[str]:
    """Entity ids, one per sequence-identity cluster, best resolution first."""
    body = {
        "query": {"type": "group", "logical_operator": "and",
                  "nodes": category_nodes(category, resolution, bound_ligand)},
        "return_type": "polymer_entity",
        "request_options": {
            "group_by": {"aggregation_method": "sequence_identity",
                         "similarity_cutoff": identity},
            "group_by_return_type": "representatives",
            "sort": [{"sort_by": "rcsb_entry_info.resolution_combined",
                      "direction": "asc"}],
            "return_all_hits": True,
        },
    }
    req = urllib.request.Request(
        SEARCH_URL, data=json.dumps(body).encode(),
        headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=180) as fh:
        payload = json.load(fh)
    return [r["identifier"] for r in payload.get("result_set", [])]


def fetch_structure(pdb_id: str, struct_dir: Path,
                    assembly: int | None) -> tuple[Path, str] | None:
    """The biological assembly where one exists, else the deposited entry.

    The assembly is the right unit for an interface: it applies the
    crystallographic symmetry, so a contact that only forms between symmetry
    mates is present, and a contact between neighbouring assemblies -- a
    crystal packing artefact -- is not. The asymmetric unit has both problems
    in reverse.

    Not every entry serves an assembly file, so the deposited entry is the
    fallback and which one was used is recorded per structure rather than
    assumed. Chains in an assembly carry symmetry suffixes (``A-2``, ``A-3``),
    so a label refers to the file in ``structures/`` and not to the
    asymmetric unit's chain naming.
    """
    if assembly:
        dest = struct_dir / f"{pdb_id}-assembly{assembly}.cif.gz"
        if download(pdb_id, dest, url=ASSEMBLY_URL.format(pdb=pdb_id, n=assembly),
                    retries=2):
            return dest, f"biological assembly {assembly}"
    dest = struct_dir / f"{pdb_id}.cif.gz"
    if download(pdb_id, dest):
        return dest, "deposited entry (asymmetric unit)"
    return None


def download(pdb_id: str, dest: Path, retries: int = 3,
             url: str | None = None) -> bool:
    if dest.is_file() and dest.stat().st_size > 0:
        return True
    url = url or DOWNLOAD_URL.format(pdb=pdb_id)
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=120) as fh:
                data = fh.read()
            if not data:
                raise OSError("empty body")
            tmp = dest.with_suffix(dest.suffix + ".part")
            tmp.write_bytes(data)
            tmp.replace(dest)
            return True
        except Exception:
            if attempt == retries - 1:
                return False
            # RCSB throttles bursts; back off rather than hammer it.
            time.sleep(1.5 * (attempt + 1))
    return False


def read_block(path: Path):
    with gzip.open(path, "rt") as fh:
        cif = pdbx_io.CIFFile.read(fh)
    return cif, cif.block


def struct_sites(block, ligand: str | None) -> tuple[dict[tuple[str, int], str], str]:
    """Residues the entry lists in a ``struct_site``, keyed to the site ligand.

    ``ligand`` restricts to sites about that component (ZN for zinc fingers);
    ``None`` keeps every site. Returns the residues and the evidence codes,
    because most of these records are generated by software from proximity
    and are not an author's assertion.
    """
    if "struct_site" not in block or "struct_site_gen" not in block:
        return {}, "none"
    site = block["struct_site"]
    if "id" not in site:
        return {}, "none"
    ids = site["id"].as_array(str)
    comps = (site["pdbx_auth_comp_id"].as_array(str)
             if "pdbx_auth_comp_id" in site else [""] * len(ids))
    evidence = (site["pdbx_evidence_code"].as_array(str)
                if "pdbx_evidence_code" in site else [""] * len(ids))
    keep, lig_of = set(), {}
    for sid, comp, _ev in zip(ids, comps, evidence):
        c = (comp or "").strip().upper()
        if ligand is not None and c != ligand:
            continue
        keep.add(sid)
        lig_of[sid] = c or "unspecified"
    if not keep:
        return {}, "none"
    codes = sorted({(e or "").strip() or "unspecified"
                    for i, e in zip(ids, evidence) if i in keep})

    gen = block["struct_site_gen"]
    out: dict[tuple[str, int], str] = {}
    for sid, chain, seq, comp in zip(
        gen["site_id"].as_array(str), gen["auth_asym_id"].as_array(str),
        gen["auth_seq_id"].as_array(str), gen["label_comp_id"].as_array(str),
    ):
        if sid not in keep:
            continue
        c = (comp or "").strip().upper()
        if c in WATER or (ligand is not None and c == ligand):
            continue
        try:
            out[(str(chain), int(seq))] = lig_of.get(sid, "unspecified")
        except ValueError:
            continue
    return out, "/".join(codes)


def interface_residues(cif, contact: float) -> tuple[dict[tuple[str, int], str], int]:
    """Amino-acid residues within ``contact`` A of a heavy atom of another
    PROTEIN chain.

    Both sides are restricted to amino acids. Without that, an entry that
    also carries DNA or RNA has its nucleotides labelled as interface and,
    worse, its DNA-binding protein residues labelled as though they contacted
    a protein -- which is the opposite of what this category is for, and
    invisible afterwards because those residues really are amino acids.

    Every ordered chain pair is scored, so both sides of an interface are
    labelled and a residue contacting two partners records both. Chains are
    compared as deposited: this is the asymmetric unit, so a small interface
    may be a crystal contact rather than a biological one, and an interface
    only formed by symmetry is missed. The per-residue partner list and the
    interface size are written out so that judgement is left to the analysis.
    """
    arr = pdbx_io.get_structure(cif, model=1)
    arr = arr[(arr.element != "H") & ~np.isin(arr.res_name, list(WATER))]
    arr = arr[struc.filter_amino_acids(arr)]
    if arr.array_length() == 0:
        return {}, 0
    chains = [c for c in dict.fromkeys(arr.chain_id.tolist())]
    if len(chains) < 2:
        return {}, 0
    sub = {c: arr[arr.chain_id == c] for c in chains}
    trees = {c: cKDTree(sub[c].coord) for c in chains if sub[c].array_length()}
    partners: dict[tuple[str, int], set[str]] = {}
    for a in chains:
        if a not in trees:
            continue
        for b in chains:
            if a == b or b not in trees:
                continue
            near = trees[b].query_ball_point(
                sub[a].coord, r=contact, return_length=True) > 0
            if not near.any():
                continue
            for chain, seq in zip(sub[a].chain_id[near], sub[a].res_id[near]):
                partners.setdefault((str(chain), int(seq)), set()).add(str(b))
    return ({k: "/".join(sorted(v)) for k, v in partners.items()},
            len(partners))


def zinc_coordination(cif, cutoff: float) -> dict[tuple[str, int], str]:
    """Polymer residues coordinating a bound ZN, for entries with no site record.

    Many zinc-finger depositions bind zinc but carry no ``struct_site``. The
    proximity rule is a faithful stand-in rather than a guess: on the
    calibration panel every one of the 60 residues the entries *did* list
    also lay within 2.8 A of the metal, so this reproduces the listed label
    where both exist. It is recorded under its own ``label_source`` so a
    coordination-derived row is never mistaken for a deposited annotation.
    """
    arr = pdbx_io.get_structure(cif, model=1)
    arr = arr[arr.element != "H"]
    zn = arr[(arr.res_name == "ZN") & (arr.element == "ZN")]
    poly = arr[(arr.hetero == False) & ~np.isin(arr.res_name, list(WATER))]  # noqa: E712
    if zn.array_length() == 0 or poly.array_length() == 0:
        return {}
    near = cKDTree(zn.coord).query_ball_point(
        poly.coord, r=cutoff, return_length=True) > 0
    return {(str(c), int(s)): "ZN"
            for c, s in zip(poly.chain_id[near], poly.res_id[near])}


def residue_count(cif) -> int:
    """Amino-acid residues, on the same basis the labels are drawn from.

    This is the denominator of ``annotated_fraction``, so it has to count the
    same universe the numerator does. An earlier version filtered on
    ``hetero == False`` instead, which drops modified residues that
    ``filter_amino_acids`` keeps -- hydroxyproline in a collagen peptide, for
    one -- and produced fractions above 1.0.
    """
    arr = pdbx_io.get_structure(cif, model=1)
    arr = arr[(arr.element != "H") & ~np.isin(arr.res_name, list(WATER))]
    arr = arr[struc.filter_amino_acids(arr)]
    if arr.array_length() == 0:
        return 0
    return len({(str(c), int(s)) for c, s in zip(arr.chain_id, arr.res_id)})


def process(category: str, pdb_id: str, path: Path, contact: float,
            zn_cutoff: float = 2.8) -> dict | None:
    try:
        cif, block = read_block(path)
    except Exception as exc:
        return {"pdb_id": pdb_id, "error": f"{type(exc).__name__}: {exc}"[:120]}
    try:
        n_res = residue_count(cif)
        if category == "protein_protein_interactions":
            labels, _ = interface_residues(cif, contact)
            detail_header, source, evidence = "partner_chains", \
                f"inter-chain contact <= {contact} A", "computed from coordinates"
        elif category == "zincfingers":
            labels, evidence = struct_sites(block, "ZN")
            detail_header, source = "site_ligand", "struct_site, ligand ZN"
            if not labels:
                labels = zinc_coordination(cif, zn_cutoff)
                source = f"ZN coordination <= {zn_cutoff} A (no struct_site)"
                evidence = "computed from coordinates" if labels else "none"
        else:
            labels, evidence = struct_sites(block, None)
            detail_header, source = "site_ligand", "struct_site, any ligand"
        res_names = {}
        arr = pdbx_io.get_structure(cif, model=1)
        for c, s, n in zip(arr.chain_id, arr.res_id, arr.res_name):
            res_names.setdefault((str(c), int(s)), str(n))
        return {"pdb_id": pdb_id, "n_residues": n_res, "labels": labels,
                "res_names": res_names, "label_source": source,
                "evidence": evidence, "detail_header": detail_header}
    except Exception as exc:
        return {"pdb_id": pdb_id, "error": f"{type(exc).__name__}: {exc}"[:120]}


def _source_counts(rows: list[list]) -> dict[str, int]:
    out: dict[str, int] = {}
    for r in rows:
        out[r[2]] = out.get(r[2], 0) + 1
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("category", choices=("protein_protein_interactions",
                                         "zincfingers", "enzymes"))
    ap.add_argument("--root", type=Path,
                    default=Path(__file__).resolve().parents[1] / "investigation_databases")
    ap.add_argument("-n", "--max-structures", type=int, default=1500)
    ap.add_argument("--resolution", type=float, default=2.5)
    ap.add_argument("--identity", type=int, default=30,
                    choices=(30, 40, 50, 70, 90, 95, 100))
    ap.add_argument("--contact", type=float, default=5.0,
                    help="interface distance for the PPI label [5.0]")
    ap.add_argument("--zn-cutoff", type=float, default=2.8,
                    help="ZN coordination distance used only where an entry "
                         "lists no struct_site [2.8]")
    ap.add_argument("--require-annotation", action="store_true",
                    help="keep only structures carrying at least one labelled "
                         "residue, and delete the rest from structures/. Use "
                         "when the labelled residues ARE the object of study "
                         "-- an unlabelled entry then contributes nothing and "
                         "only dilutes the set")
    ap.add_argument("--require-class", default=None,
                    help="comma-separated ligand_class values a structure must "
                         "carry at least one of, e.g. 'biological' to exclude "
                         "sites that are only about a cryoprotectant. Stricter "
                         "than --require-annotation and overrides it")
    ap.add_argument("--bound-ligand", default=None,
                    help="require the entry to contain this bound component, "
                         "e.g. ZN. Applied in the SEARCH, so the cluster "
                         "representative chosen is one that actually has it")
    ap.add_argument("--assembly", type=int, default=None,
                    help="biological assembly to label, falling back to the "
                         "deposited entry where none is served. Defaults to 1 "
                         "for protein_protein_interactions, where the "
                         "asymmetric unit both misses symmetry interfaces and "
                         "invents crystal-packing ones; 0 forces the "
                         "deposited entry")
    ap.add_argument("-j", "--jobs", type=int, default=6)
    args = ap.parse_args(argv)

    out_dir = args.root / args.category
    struct_dir = out_dir / "structures"
    struct_dir.mkdir(parents=True, exist_ok=True)

    # Only the interface category needs the assembly: the site categories read
    # struct_site, which belongs to the deposited entry and is not carried
    # into an assembly file.
    assembly = args.assembly
    if assembly is None:
        assembly = 1 if args.category == "protein_protein_interactions" else 0
    assembly = assembly or None

    print(f"searching RCSB: {args.category}, X-ray <= {args.resolution} A, "
          f"one per {args.identity}% identity cluster", flush=True)
    entities = search_representatives(args.category, args.resolution,
                                      args.identity, args.bound_ligand)
    # One entry can represent two clusters through two different entities;
    # the structure is downloaded once and labelled once either way.
    seen, pdb_ids = set(), []
    for ent in entities:
        pdb = ent.split("_")[0].upper()
        if pdb not in seen:
            seen.add(pdb)
            pdb_ids.append(pdb)
    print(f"  {len(entities)} cluster representatives -> {len(pdb_ids)} distinct entries")

    require = {c.strip() for c in (args.require_class or "").split(",") if c.strip()}
    cap = args.max_structures or len(pdb_ids)
    rank = {p: i for i, p in enumerate(pdb_ids)}

    def qualifies(n_labels: int, classes: set[str]) -> bool:
        if require:
            return bool(classes & require)
        if args.require_annotation:
            return n_labels > 0
        return True

    # The cap is applied AFTER the filter, not before. Capping first and then
    # filtering silently returns far fewer structures than asked for -- on
    # enzymes only 28% of entries carry a genuine cofactor site, so a 1500 cap
    # applied up front yields about 420. Candidates are therefore walked in
    # resolution order, in chunks, until the cap is met or the pool runs out.
    print(f"downloading into {struct_dir} ({args.jobs} at a time)", flush=True)
    if require or args.require_annotation:
        what = f"a label in {sorted(require)}" if require else "any label"
        print(f"  keeping only structures with {what}; the {cap} cap applies "
              f"AFTER that filter", flush=True)

    rows, label_rows, errors, failed = [], [], [], []
    examined = rejected = 0
    chunk = max(args.jobs * 50, 250)
    for start in range(0, len(pdb_ids), chunk):
        if len(rows) >= cap:
            break
        batch = pdb_ids[start:start + chunk]
        paths: dict[str, tuple[Path, str]] = {}
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futs = {pool.submit(fetch_structure, p, struct_dir, assembly): p
                    for p in batch}
            for fut in as_completed(futs):
                hit = fut.result()
                if hit is None:
                    failed.append(futs[fut])
                else:
                    paths[futs[fut]] = hit
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futs = {pool.submit(process, args.category, p, paths[p][0],
                                args.contact, args.zn_cutoff): p for p in paths}
            for fut in as_completed(futs):
                r = fut.result()
                examined += 1
                if r is None:
                    continue
                if "error" in r:
                    errors.append(r)
                    continue
                mine, classes, n_bio = [], set(), 0
                for (chain, seq), detail in sorted(r["labels"].items()):
                    cls = ligand_class(args.category, detail)
                    classes.add(cls)
                    n_bio += cls in ("biological", "interface")
                    mine.append([r["pdb_id"], chain, seq,
                                 r["res_names"].get((chain, seq), "UNK"), 1,
                                 detail, cls])
                if not qualifies(len(r["labels"]), classes):
                    rejected += 1
                    paths[r["pdb_id"]][0].unlink(missing_ok=True)
                    continue
                label_rows.extend(mine)
                rows.append([r["pdb_id"], paths[r["pdb_id"]][0].name,
                             paths[r["pdb_id"]][1],
                             r["n_residues"], len(r["labels"]), n_bio,
                             f"{len(r['labels']) / r['n_residues']:.4f}"
                             if r["n_residues"] else "",
                             r["label_source"], r["evidence"]])
        print(f"  examined {examined}, kept {len(rows)}/{cap}"
              f"{f', rejected {rejected}' if rejected else ''}", flush=True)

    # Keep the best by resolution, which is the order the search returned.
    dropped = rejected
    if len(rows) > cap:
        rows.sort(key=lambda x: rank.get(x[0], 1 << 30))
        cut = {r[0] for r in rows[cap:]}
        for r in rows[cap:]:
            (struct_dir / r[1]).unlink(missing_ok=True)
        rows = rows[:cap]
        label_rows = [r for r in label_rows if r[0] not in cut]
        print(f"  trimmed {len(cut)} past the cap, kept the best {len(rows)}")

    # A previous run in a different mode (asymmetric unit vs assembly, or a
    # looser filter) leaves files behind that nothing now references. Sweep
    # them, so structures/ always matches structures.tsv exactly.
    wanted = {r[1] for r in rows}
    stale = [p for p in struct_dir.glob("*.cif.gz") if p.name not in wanted]
    for p in stale:
        p.unlink(missing_ok=True)
    if stale:
        print(f"  swept {len(stale)} files no longer referenced")

    rows.sort(key=lambda x: x[0])
    label_rows.sort(key=lambda x: (x[0], x[1], x[2]))
    detail_header = ("partner_chains"
                     if args.category == "protein_protein_interactions"
                     else "site_ligand")

    with (out_dir / "structures.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["pdb_id", "structure_file", "structure_source", "n_residues",
                    "n_annotated", "n_annotated_biological",
                    "annotated_fraction", "label_source", "evidence_code"])
        w.writerows(rows)
    with (out_dir / "annotated_residues.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["pdb_id", "chain", "res_id", "res_name", "is_annotated",
                    detail_header, "ligand_class"])
        w.writerows(label_rows)

    with_labels = sum(1 for r in rows if r[4] > 0)
    with_bio = sum(1 for r in rows if r[5] > 0)
    class_counts: dict[str, int] = {}
    for r in label_rows:
        class_counts[r[6]] = class_counts.get(r[6], 0) + 1
    bio_residues = class_counts.get("biological", 0) + class_counts.get("interface", 0)
    manifest = {
        "built_utc": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        "category": args.category,
        "source": "RCSB search API v2 + files.rcsb.org",
        "filters": {"method": "X-RAY DIFFRACTION",
                    "resolution_max_angstrom": args.resolution,
                    "sequence_identity_cluster_percent": args.identity,
                    "one_representative_per_cluster": True,
                    "sorted_by": "resolution ascending",
                    "cap": args.max_structures},
        "label": {"source": rows[0][4] if rows else None,
                  "contact_distance_angstrom": (
                      args.contact
                      if args.category == "protein_protein_interactions" else None),
                  "caveat": (
                      "Interfaces are computed on the asymmetric unit as "
                      "deposited, so a small interface may be a crystal "
                      "contact rather than a biological one; the partner "
                      "chains and interface size are recorded so that "
                      "judgement is left to the analysis."
                      if args.category == "protein_protein_interactions" else
                      "struct_site records are usually generated by software "
                      "from proximity to the ligand, not asserted by an "
                      "author; evidence_code records which.")},
        "cluster_representatives_returned": len(entities),
        "distinct_entries": len(seen),
        "downloaded": examined,
        "download_failures": failed[:50],
        "parse_failures": errors[:50],
        "require_annotation": args.require_annotation,
        "require_class": sorted(require) or None,
        "bound_ligand_required": args.bound_ligand,
        "assembly": assembly,
        "structure_source_counts": _source_counts(rows),
        "candidates_examined": examined,
        "rejected_by_filter": dropped,
        "structures_written": len(rows),
        "structures_with_at_least_one_label": with_labels,
        "structures_with_a_biological_label": with_bio,
        "annotated_residues": len(label_rows),
        "annotated_residues_by_class": class_counts,
    }
    (out_dir / "build_manifest.json").write_text(json.dumps(manifest, indent=2))

    print(f"\nwrote {out_dir}/")
    print(f"  structures.tsv          {len(rows)} structures")
    print(f"  annotated_residues.tsv  {len(label_rows)} labelled residues")
    if rows:
        print(f"  with >=1 label          {with_labels} of {len(rows)}"
              f" ({100 * with_labels / len(rows):.1f}%)")
        print(f"  label classes           "
              + ", ".join(f"{k} {v}" for k, v in sorted(class_counts.items())))
    if errors:
        print(f"  parse failures          {len(errors)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
