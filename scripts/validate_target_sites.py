"""Audit biological-site overlap of retrieved FoldDisco target residues.

Three resumable stages: collect qualified hits, annotate local indexed PDB
structures, and report overlap against same-target-chain/residue-type controls.
No homology-transferred catalytic labels are used. Missing evidence is unknown,
not negative. Coordinate labels describe observed zinc/ligand contacts and the
first declared biological assembly, not proof of function or binding affinity.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import sqlite3
import sys
import urllib.request
import warnings

import numpy as np
from scipy.spatial import cKDTree
import biotite.structure as struc
from biotite.structure.io.pdb import PDBFile, list_assemblies

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))
from ibc_search.folddisco_contract import parse_query_token  # noqa: E402
from ibc_search.errors import ContractError  # noqa: E402
from build_investigation_db import ligand_class  # noqa: E402
from count_scoring import split_patch_ids  # noqa: E402
from top_motifs_by_organism import pdb_id_from_tid  # noqa: E402

MCSA_URL = "https://www.ebi.ac.uk/thornton-srv/m-csa/api/residues/?format=json"
RUNS = {
    "background": "background_dataset_dumb_r5.5",
    "zincfingers": "zincfingers_seeded_r5.5",
    "enzymes": "enzymes_seeded2_r5.5",
    "protein_protein_interactions": "protein_protein_interactions_seeded2_r5.5",
}
TRACKS = {
    "zincfingers": ("zinc",),
    "enzymes": ("catalytic", "ligand"),
    "protein_protein_interactions": ("interface",),
}
SCHEMA_VERSION = 1


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".part")
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    tmp.replace(path)


def rows(path: Path):
    with path.open() as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def qualified_hit(row: dict, min_size=6, max_size=14,
                  max_rmsd=2.0, min_coverage=0.8):
    """Return validated positional tokens, or a reason for exclusion."""
    q = row["query_residues"].split(",")
    t = row["matching_residues"].split(",")
    matched = [x for x in t if x != "_"]
    if not min_size <= len(q) <= max_size:
        return None, "size"
    if len(q) != len(t):
        return None, "unalignable"
    if len(set(q)) != len(q) or len(set(matched)) != len(matched):
        return None, "noninjective"
    try:
        for token in q + matched:
            parse_query_token(token)
        rmsd = float(row["rmsd"])
    except (ValueError, ContractError):
        return None, "invalid"
    if not math.isfinite(rmsd) or rmsd < 0 or rmsd > max_rmsd:
        return None, "rmsd"
    if len(matched) / len(q) < min_coverage:
        return None, "coverage"
    return (q, t, rmsd), "qualified"


def collect(args) -> None:
    """Store unique query-motif/target-site observations, preserving zero motifs."""
    dbpath = args.out / "hits.sqlite"
    if dbpath.exists():
        raise ValueError(f"{dbpath} exists; reuse annotate/report or choose a new output directory")
    args.out.mkdir(parents=True, exist_ok=True)
    db = sqlite3.connect(dbpath)
    db.executescript("""
        CREATE TABLE motifs(category TEXT, instance TEXT, patch TEXT, size INT,
          PRIMARY KEY(category,instance,patch));
        CREATE TABLE hits(category TEXT, instance TEXT, patch TEXT, size INT,
          target TEXT, tid TEXT, tokens TEXT, anchor TEXT, rmsd REAL,
          PRIMARY KEY(category,instance,patch,target,tokens));
        CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT);
    """)
    background_meta = {r["structure_id"]: r for r in rows(
        REPO / "investigation_databases/background_dataset/structures.tsv")}
    audit = {}; sources = []
    for category, run in RUNS.items():
        base = REPO / "results" / run
        annotations = defaultdict(dict)
        if category != "background":
            path = REPO / "investigation_databases" / category / "annotated_residues.tsv"
            sources.append({"path": str(path), "sha256": digest(path)})
            for r in rows(path):
                annotations[r["pdb_id"]][r["chain"] + r["res_id"]] = r["ligand_class"]
        counts = Counter()
        for item in rows(base / "panel.tsv"):
            inst = item["instance"]
            if category == "background" and "X-RAY" not in background_meta[inst]["method"]:
                counts["non_xray_background_structures_excluded"] += 1
                continue
            directory = base / inst
            selected = {}
            for p in rows(directory / "patches.tsv"):
                if p["emit_query"].lower() != "true":
                    continue
                seed = p["seed_residue_labels"]
                if category != "background":
                    cls = annotations[inst].get(seed)
                    if cls is None or (category == "enzymes" and cls != "biological"):
                        continue
                size = int(p["n_members"])
                if not args.min_size <= size <= args.max_size:
                    continue
                # Identical residue sets are one query, not independent seeds.
                query = p["query_string"]
                selected[p["patch_id"]] = (size, seed, query)
            canonical = {}
            aliases = {}
            for patch, (size, seed, query) in selected.items():
                canonical.setdefault(query, patch)
                aliases[patch] = canonical[query]
                if canonical[query] == patch:
                    db.execute("INSERT INTO motifs VALUES(?,?,?,?)", (category, inst, patch, size))
            counts["motifs"] += len(canonical)
            sources.append({"path": str(directory / "folddisco_hits.tsv"),
                            "sha256": digest(directory / "folddisco_hits.tsv")})
            for r in rows(directory / "folddisco_hits.tsv"):
                counts["rows_read"] += 1
                pids = [p for p in split_patch_ids(r["patch_id"]) if p in selected]
                if not pids:
                    continue
                parsed, reason = qualified_hit(r, args.min_size, args.max_size,
                                               args.max_rmsd, args.min_coverage)
                counts[reason] += 1
                if parsed is None:
                    continue
                target = pdb_id_from_tid(r["tid"])
                if target is None:
                    counts["unknown_target_id"] += 1
                    continue
                if target == inst.lower():
                    counts["self_hit_excluded"] += 1
                    continue
                q, t, rmsd = parsed
                tokens = ",".join(sorted(x for x in t if x != "_"))
                for pid in pids:
                    size, seed, query = selected[pid]
                    if q != query.split(",") or size != len(q):
                        raise ValueError(f"query contract mismatch: {directory}, {pid}")
                    anchor = t[q.index(seed)] if seed in q else "_"
                    db.execute("""INSERT INTO hits VALUES(?,?,?,?,?,?,?,?,?)
                        ON CONFLICT(category,instance,patch,target,tokens)
                        DO UPDATE SET rmsd=min(rmsd,excluded.rmsd)""",
                        (category, inst, aliases[pid], size, target, r["tid"], tokens,
                         anchor, rmsd))
            counts["structures"] += 1
            db.commit()
        audit[category] = dict(counts)
        print("collected", category, dict(counts), flush=True)
    db.executescript("CREATE INDEX target_idx ON hits(target); CREATE INDEX category_idx ON hits(category);")
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "filters": {k: getattr(args, k) for k in
                    ("min_size", "max_size", "max_rmsd", "min_coverage")},
        "query_selection": "original annotated seeds; biological ligand class for enzymes; X-ray background; identical queries collapsed",
        "audit": audit, "input_files": sources,
        "caveats": ["Existing top-40 retrieval, not exhaustive search.",
                    "Equal-length fields follow the existing positional contract; no independent coordinate remapping.",
                    "Self entries excluded; homologues are NOT excluded.",
                    "Query ligand-contact labels are not curated catalytic labels."],
    }
    db.execute("INSERT INTO metadata VALUES(?,?)", ("manifest", json.dumps(manifest)))
    db.commit(); db.close()
    write_json(args.out / "collection_manifest.json", manifest)


def load_mcsa(path: Path):
    """Use only explicit reference PDB/author residue mappings; never homologues."""
    data = json.loads(path.read_text())
    out = defaultdict(dict)
    for residue in data:
        for chain in residue.get("residue_chains", []):
            if chain.get("is_reference") is not True:
                continue
            number = chain.get("auth_resid")
            name = chain.get("chain_name")
            try:
                number = int(str(number))
            except (TypeError, ValueError):
                continue
            if not name or len(name) != 1:
                continue
            token = f"{name}{number}"
            out[chain["pdb_id"].lower()][token] = {
                "name": chain["code"].upper(), "mcsa_id": residue["mcsa_id"],
                "roles": residue.get("roles_summary", ""),
            }
    return dict(out)


def residue_records(array):
    """Only uniquely identifiable amino-acid residues; insertion collisions unknown."""
    aa = array[struc.filter_amino_acids(array) & ~np.isin(array.element, ["H", "D"])]
    identities = defaultdict(set)
    for c, n, ins, name in zip(aa.chain_id, aa.res_id, aa.ins_code, aa.res_name):
        identities[f"{c}{n}"].add((str(ins).strip(), str(name)))
    valid = {t: {"name": next(iter(ids))[1]} for t, ids in identities.items()
             if len(ids) == 1 and next(iter(ids))[0] == ""}
    return aa, valid, sorted(set(identities) - set(valid))


def contact_tokens(protein, other, cutoff):
    if not len(protein) or not len(other):
        return set()
    near = cKDTree(other.coord).query_ball_point(protein.coord, cutoff, return_length=True) > 0
    return {f"{c}{n}" for c, n in zip(protein.chain_id[near], protein.res_id[near])}


def assembly_interfaces(array, cutoff=5.0):
    """Separate symmetry copies, including homomeric copies with the same chain ID."""
    aa = array[struc.filter_amino_acids(array) & ~np.isin(array.element, ["H", "D"])]
    if "sym_id" not in aa.get_annotation_categories():
        raise ValueError("assembly has no symmetry-copy identifiers")
    groups = defaultdict(list)
    for i, key in enumerate(zip(aa.chain_id, aa.sym_id)):
        groups[(str(key[0]), int(key[1]))].append(i)
    labels = set()
    # A single neighbour graph is substantially cheaper than all chain pairs.
    if len(groups) > 1:
        tree = cKDTree(aa.coord)
        for indices in groups.values():
            index_set = set(indices)
            for i, neighbours in zip(indices, tree.query_ball_point(aa.coord[indices], cutoff)):
                if any(j not in index_set for j in neighbours):
                    labels.add(f"{aa.chain_id[i]}{aa.res_id[i]}")
    present = {f"{c}{n}" for c, n in zip(aa.chain_id, aa.res_id)}
    return labels, present


def annotate_structure(path: Path, catalytic: dict, max_atoms=150000):
    with gzip.open(path, "rt") as fh:
        pdb = PDBFile.read(fh)
    arr = pdb.get_structure(model=1)
    if len(arr) > max_atoms:
        raise ValueError(f"asymmetric unit exceeds atom safety limit {max_atoms}")
    protein, residues, ambiguous = residue_records(arr)
    heavy = arr[~np.isin(arr.element, ["H", "D"])]
    zn = heavy[(heavy.res_name == "ZN") & (heavy.element == "ZN")]
    biological_names = {str(n) for n in set(heavy.res_name[heavy.hetero])
                        if ligand_class("enzymes", str(n)) == "biological"
                        and str(n) not in {"HOH", "DOD", "WAT", "H2O"}}
    ligands = heavy[heavy.hetero & np.isin(heavy.res_name, list(biological_names))
                    & ~struc.filter_amino_acids(heavy)]
    zinc_labels = contact_tokens(protein, zn, 2.8)
    ligand_labels = contact_tokens(protein, ligands, 4.0)
    interface_labels = set(); interface_present = set(); assembly_note = "unknown"
    try:
        ids = list_assemblies(pdb)
        if ids:
            assembly = pdb.get_assembly(assembly_id=ids[0], model=1)
            if len(assembly) > max_atoms:
                raise ValueError(f"assembly exceeds atom safety limit {max_atoms}")
            interface_labels, interface_present = assembly_interfaces(assembly)
            assembly_note = f"first declared biological assembly {ids[0]}"
    except (ValueError, KeyError, IndexError) as exc:
        assembly_note = f"unknown: {type(exc).__name__}: {exc}"
    validated_catalytic = {t: v for t, v in catalytic.items()
                           if t in residues and residues[t]["name"] == v["name"]}
    known_catalytic_chains = {t[0] for t in validated_catalytic}
    for token, r in residues.items():
        r["zinc"] = int(token in zinc_labels) if len(zn) else None
        r["ligand"] = int(token in ligand_labels) if len(ligands) else None
        r["interface"] = int(token in interface_labels) if token in interface_present else None
        r["catalytic"] = int(token in validated_catalytic) if token[0] in known_catalytic_chains else None
    return {
        "residues": residues, "ambiguous_tokens": ambiguous,
        "assembly": assembly_note, "zinc_atoms": len(zn),
        "ligand_names": sorted(biological_names),
        "curated_residues": validated_catalytic,
        "curated_unmapped": sorted(set(catalytic) - set(validated_catalytic)),
    }


def annotate_one(task):
    target, path, catalytic, cache, max_atoms, mcsa_hash = task
    dest = Path(cache) / f"{target}.json"
    if dest.exists():
        existing = json.loads(dest.read_text())
        if existing.get("schema_version") == SCHEMA_VERSION and existing.get("mcsa_sha256") == mcsa_hash:
            return target, existing["status"]
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            value = annotate_structure(Path(path), catalytic, max_atoms)
        value.update(status="ok", source_path=str(path), source_sha256=digest(Path(path)))
    except Exception as exc:
        value = {"status": "unknown", "error": f"{type(exc).__name__}: {exc}", "residues": {}}
    value.update(schema_version=SCHEMA_VERSION, mcsa_sha256=mcsa_hash)
    write_json(dest, value)
    return target, value["status"]


def annotate(args) -> None:
    cache = args.out / "target_annotations"; cache.mkdir(exist_ok=True, parents=True)
    mcsa_path = args.out / "mcsa_reference_residues.json"
    if not mcsa_path.exists():
        with urllib.request.urlopen(MCSA_URL, timeout=120) as fh:
            data = json.load(fh)
        write_json(mcsa_path, data)
    mcsa = load_mcsa(mcsa_path); mcsa_hash = digest(mcsa_path)
    with sqlite3.connect(args.out / "hits.sqlite") as db:
        targets = [r[0] for r in db.execute("SELECT DISTINCT target FROM hits")]
    # Sampling depends only on target ID and seed, never on biological labels.
    targets.sort(key=lambda t: hashlib.sha256(f"{args.seed}:{t}".encode()).hexdigest())
    if args.curated_only:
        selected = [t for t in targets if t in mcsa]
    else:
        selected = targets[:args.max_targets] if args.max_targets else targets
    write_json(args.out / "annotation_manifest.json", {
        "target_universe": len(targets), "selected_targets": len(selected),
        "selection": ("Census of retrieved targets with explicit M-CSA reference PDB mappings"
                      if args.curated_only else
                      "SHA256(seed:target_id), independent of labels; shared across categories"),
        "seed": args.seed, "max_targets": args.max_targets,
        "max_atoms": args.max_atoms, "mcsa_url": MCSA_URL, "mcsa_sha256": mcsa_hash,
        "target_ids": selected,
    })
    tasks = [(t, str(args.pdb_root / "pdb" / t[1:3] / f"pdb{t}.ent.gz"),
              mcsa.get(t, {}), str(cache), args.max_atoms, mcsa_hash) for t in selected]
    counts = Counter()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for i, (_, status) in enumerate(pool.map(annotate_one, tasks, chunksize=5), 1):
            counts[status] += 1
            if i % 100 == 0 or i == len(tasks):
                print(f"annotated {i}/{len(tasks)}: {dict(counts)}", flush=True)


def overlap(tokens, annotation, track):
    """Observed fraction and exact expected fraction under type/chain controls.

    Every same-name residue in the same target chain is a possible control.
    This controls target, chain and amino-acid composition, NOT local geometry
    or solvent exposure. Entire sites with missing/ambiguous mappings are unknown.
    """
    residues = annotation.get("residues", {})
    if not tokens or any(t not in residues or residues[t].get(track) is None for t in tokens):
        return None
    means = annotation.setdefault("_control_means", {})
    if track not in means:
        pools = defaultdict(list)
        for t, r in residues.items():
            if r.get(track) is not None:
                pools[(t[0], r["name"])].append(r[track])
        means[track] = {k: float(np.mean(v)) for k, v in pools.items()}
    observed = float(np.mean([residues[t][track] for t in tokens]))
    expected = float(np.mean([means[track][(t[0], residues[t]["name"])] for t in tokens]))
    return observed, expected


def standardized_difference(annotated, background, bootstrap=2000, seed=5):
    """Query-structure-cluster bootstrap, standardized to background motif sizes."""
    sizes = sorted({r[1] for r in annotated} & {r[1] for r in background})
    sizes = [s for s in sizes if sum(r[1] == s for r in annotated) >= 10
             and sum(r[1] == s for r in background) >= 10]
    if not sizes:
        return {"status": "insufficient_common_size_strata"}
    aa = [r for r in annotated if r[1] in sizes]; bb = [r for r in background if r[1] in sizes]
    def matrix(values):
        ids = sorted({r[0] for r in values}); ii = {v: i for i, v in enumerate(ids)}
        ss = {v: i for i, v in enumerate(sizes)}
        n = np.zeros((len(ids), len(sizes))); y = n.copy()
        for inst, size, value in values:
            n[ii[inst], ss[size]] += 1; y[ii[inst], ss[size]] += value
        return n, y
    an, ay = matrix(aa); bn, by = matrix(bb); w = bn.sum(axis=0) / bn.sum()
    a = float(w @ (ay.sum(axis=0) / an.sum(axis=0)))
    b = float(w @ (by.sum(axis=0) / bn.sum(axis=0)))
    rng = np.random.default_rng(seed); differences = []
    for _ in range(bootstrap):
        ia = rng.integers(len(an), size=len(an)); ib = rng.integers(len(bn), size=len(bn))
        na = an[ia].sum(axis=0); nb = bn[ib].sum(axis=0)
        if (na == 0).any() or (nb == 0).any():
            continue
        differences.append(float(w @ (ay[ia].sum(axis=0) / na - by[ib].sum(axis=0) / nb)))
    return {"status": "exploratory", "sizes": sizes, "annotated": a, "background": b,
            "difference": a - b, "motifs": [len(aa), len(bb)],
            "query_structures": [len(an), len(bn)], "bootstrap_valid": len(differences),
            "difference_ci95": np.quantile(differences, [0.025, 0.975]).tolist() if differences else None,
            "bootstrap_warning": ("Many resamples lack a size stratum; interval is conditional on valid resamples. Treat as unstable."
                                  if len(differences) < 0.9 * bootstrap else None),
            "uncertainty_scope": "Query-structure bootstrap conditional on sampled target pool; NOT target/family-held-out uncertainty."}


def report(args) -> None:
    selection = json.loads((args.out / "annotation_manifest.json").read_text())
    annotations = {t: json.loads((args.out / "target_annotations" / f"{t}.json").read_text())
                   for t in selection["target_ids"]
                   if (args.out / "target_annotations" / f"{t}.json").exists()}
    accum = defaultdict(lambda: [0.0, 0.0, 0]); anchor_accum = defaultdict(lambda: [0.0, 0.0, 0])
    coverage = defaultdict(Counter); distinct_sites = defaultdict(dict)
    db = sqlite3.connect(args.out / "hits.sqlite")
    header = ["category", "query_instance", "patch_id", "motif_size", "target_pdb", "target_tokens",
              "track", "observed_fraction", "control_expected_fraction", "excess_fraction",
              "anchor_target", "anchor_observed", "anchor_expected", "rmsd", "evidence"]
    with (args.out / "target_site_overlap.tsv").open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n"); writer.writerow(header)
        for category, inst, patch, size, target, tid, token_field, anchor, rmsd in db.execute("SELECT * FROM hits"):
            tracks = ("zinc", "catalytic", "ligand", "interface") if category == "background" else TRACKS[category]
            for track in tracks:
                counter = coverage[(category, track)]; counter["qualified_motif_site_pairs"] += 1
                if target not in annotations:
                    counter["outside_annotated_target_sample"] += 1; continue
                counter["sampled_pairs"] += 1
                annotation = annotations[target]
                result = overlap(token_field.split(","), annotation, track)
                if result is None:
                    counter["unknown_pairs"] += 1; continue
                counter["assessable_pairs"] += 1
                obs, exp = result
                counter["pairs_with_any_site_overlap"] += obs > 0
                key = (category, track, inst, patch, size)
                a = accum[key]; a[0] += obs; a[1] += exp; a[2] += 1
                distinct_sites[(category, track)][(target, token_field)] = (obs, exp)
                ar = overlap([anchor], annotation, track) if anchor != "_" else None
                if ar is not None:
                    a = anchor_accum[key]; a[0] += ar[0]; a[1] += ar[1]; a[2] += 1
                evidence = {"zinc": "observed ZN heavy-atom contact <=2.8 A",
                            "ligand": "observed non-additive/non-metal ligand contact <=4 A",
                            "catalytic": "M-CSA manually curated reference residue; name-validated",
                            "interface": annotation.get("assembly", "unknown")}[track]
                writer.writerow([category, inst, patch, size, target, token_field, track,
                                 obs, exp, obs-exp, anchor,
                                 ar[0] if ar else "", ar[1] if ar else "", rmsd, evidence])
    db.close()
    comparisons = {}; summaries = {}
    for category, tracks in TRACKS.items():
        for track in tracks:
            label = f"{category}:{track}"
            metrics = {}
            for metric, column in [("observed_overlap", 0), ("control_expected_overlap", 1), ("excess_overlap", None)]:
                def values(cat):
                    return [(key[2], key[4], (v[0]-v[1])/v[2] if column is None else v[column]/v[2])
                            for key, v in accum.items() if key[:2] == (cat, track)]
                metrics[metric] = standardized_difference(values(category), values("background"), args.bootstrap, args.seed)
            comparisons[label] = metrics
    for (category, track), counter in sorted(coverage.items()):
        rr = [v for k, v in accum.items() if k[:2] == (category, track)]
        anchors = [v for k, v in anchor_accum.items() if k[:2] == (category, track)]
        ds = list(distinct_sites[(category, track)].values())
        summaries[f"{category}:{track}"] = {
            "coverage": dict(counter), "assessable_query_motifs": len(rr),
            "motif_mean_observed": float(np.mean([v[0]/v[2] for v in rr])) if rr else None,
            "motif_mean_expected": float(np.mean([v[1]/v[2] for v in rr])) if rr else None,
            "anchor_assessable_motifs": len(anchors),
            "anchor_mean_observed": float(np.mean([v[0]/v[2] for v in anchors])) if anchors else None,
            "anchor_mean_expected": float(np.mean([v[1]/v[2] for v in anchors])) if anchors else None,
            "distinct_target_sites": len(ds),
            "distinct_site_mean_observed": float(np.mean([v[0] for v in ds])) if ds else None,
            "distinct_site_mean_expected": float(np.mean([v[1] for v in ds])) if ds else None,
        }
    result = {
        "selection": selection, "summaries": summaries, "comparisons": comparisons,
        "annotation_status": dict(Counter(v["status"] for v in annotations.values())),
        "interpretation": [
            "Site overlap, not proof of shared function, mechanism, specificity or affinity.",
            "Unknowns excluded, explicitly counted; no Zn/ligand or no curated chain is NOT a biological negative.",
            "Interface negatives mean no observed contact in first declared assembly, not inability to interact.",
            "Control matches target chain and residue name, but NOT exposure, geometry or sequence family.",
            "Per-motif average avoids counting every target hit as an independent query replicate.",
            "CIs condition on the target sample; shared targets and homologues remain dependencies.",
            "Background and annotated query structures may differ in fold, quality and deposition biases.",
            "Curated catalytic overlap validates catalytic localization, NOT identical enzyme function.",
        ],
    }
    write_json(args.out / "target_site_validation.json", result)
    lines = ["# Retrieved target-site validation", "", f"Annotated targets: {len(annotations)} / {selection['target_universe']}", "",
             "| Track | Annotated overlap | Background overlap | Difference (95% query-cluster CI) |", "|---|---:|---:|---|"]
    for label, metrics in comparisons.items():
        m = metrics["observed_overlap"]
        if "difference" not in m:
            lines.append(f"| {label} | insufficient data | | |")
            continue
        ci = m["difference_ci95"]
        text = f"{m['difference']:+.3f}" + (f" ({ci[0]:+.3f}, {ci[1]:+.3f})" if ci else "")
        lines.append(f"| {label} | {m['annotated']:.3f} | {m['background']:.3f} | {text} |")
    lines += ["", "## Excess over same-chain / amino-acid controls", "",
              "| Track | Annotated excess | Background excess | Difference |", "|---|---:|---:|---:|"]
    for label, metrics in comparisons.items():
        m = metrics["excess_overlap"]
        if "difference" in m:
            lines.append(f"| {label} | {m['annotated']:+.3f} | {m['background']:+.3f} | {m['difference']:+.3f} |")
    lines += ["", "Values are fractions of matched residues, averaged per query motif and size-standardized; they are not percentages of whole motifs correctly classified."]
    for label, metrics in comparisons.items():
        m = metrics["observed_overlap"]
        if m.get("bootstrap_warning"):
            lines.append(f"\nWARNING ({label}): {m['bootstrap_warning']}")
    lines += ["", "## Limitations", ""] + [f"- {s}" for s in result["interpretation"]]
    (args.out / "REPORT.md").write_text("\n".join(lines) + "\n")
    print("report:", args.out / "REPORT.md", flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("stage", choices=["collect", "annotate", "report", "all"])
    ap.add_argument("--out", type=Path, default=REPO / "results/target_site_validation")
    ap.add_argument("--pdb-root", type=Path, default=Path(
        "~/documents_local_only/master_thesis_data_files/folddisco_databases/index").expanduser())
    ap.add_argument("--min-size", type=int, default=6)
    ap.add_argument("--max-size", type=int, default=14)
    ap.add_argument("--max-rmsd", type=float, default=2.0)
    ap.add_argument("--min-coverage", type=float, default=0.8)
    ap.add_argument("--max-targets", type=int, default=2000, help="0 means all targets")
    ap.add_argument("--curated-only", action="store_true",
                    help="annotate every retrieved M-CSA reference target, ignoring --max-targets; use a separate output directory")
    ap.add_argument("--max-atoms", type=int, default=150000)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--seed", type=int, default=5)
    ap.add_argument("--bootstrap", type=int, default=2000)
    args = ap.parse_args(argv)
    args.out = args.out.resolve(); args.pdb_root = args.pdb_root.resolve()
    if args.min_size < 1 or args.max_size < args.min_size or not 0 < args.min_coverage <= 1 or args.max_rmsd < 0:
        ap.error("invalid match filters")
    if args.max_targets < 0 or args.max_atoms < 1 or args.workers < 1 or args.bootstrap < 0:
        ap.error("invalid annotation/statistics limits")
    for stage, func in [("collect", collect), ("annotate", annotate), ("report", report)]:
        if args.stage in (stage, "all"):
            func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())