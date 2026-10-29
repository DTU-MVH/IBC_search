"""Diagnostic benchmark: does quality-controlled epitope transfer carry signal?

Run on a small, reproducible, family-diverse panel of SAbDab2 antigens. For
each antigen it scans surface motifs, searches them against the SAbDab2
antigen index, KEEPS the raw hits, and then reports:

1. Self-retrieval control. Does a motif retrieve its own coordinates with the
   right residues? A failure here is a search/mapping problem, not biology.
2. Stage accounting. Where evidence is lost: extracted -> searched ->
   retained after leakage exclusion -> identity-valid -> >= N families.
3. Candidate-budget sensitivity from independent real searches. The reference
   is training-only and excludes all query connected families and PDB entries.
4. Family-balanced epitope transfer, compared with two baselines, on the
   within-antigen contrast epitope-derived vs other-exposed-surface motifs.

Definitions (deliberately conservative)
---------------------------------------
* A "valid" match: node_count == query size, no repeated target token, every
  target token resolves to exactly ONE labelled residue (no insertion-code
  ambiguity), RMSD <= --max-rmsd.
* Per hit structure the best valid match (lowest RMSD) is used, then
  structures are averaged within an antigen-cluster family, then families are
  averaged. A hit row is not an independent observation; a family is.
* Only target residues that are exposed count. The score is the contact
  fraction minus the hit structure's own exposed-epitope base rate.
* Motif classes: "epitope" (>= --pos-coverage of residues are contacts),
  "other_exposed" (none are; NOT confirmed non-epitopes), "straddle"
  (excluded from AUC).
* No usable match means NO EVIDENCE. It is never scored as zero.

Caveats
-------
* The IDF baseline here is log10(max IDF) without the size correction in
  idf_zscore.py. It is only compared WITHIN an antigen, but motif size still
  confounds it. Treat it as a rough baseline.
* Residue-level label transfer is not done: matching_residues is not reliably
  positionally aligned to query_residues (docs/FOLDDISCO_CONTRACT.md), so only
  the patch-level contact fraction of the matched set is used.
* Development panel defaults to the TRAIN partition. Fix the procedure, then
  evaluate once on --split test.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import random
import statistics as st
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
UNMATCHED = "_"


# ---------------------------------------------------------------- utilities

def families(value: str) -> set[str]:
    return {c for c in (value or "").split("/") if c}


def family_of(value: str) -> str:
    """Only a single known cluster is an independent family in this pilot."""
    fs = families(value)
    if len(fs) != 1:
        raise ValueError(f"expected one known antigen cluster, got {value!r}")
    return next(iter(fs))


def family_components(meta: dict[str, dict]) -> dict[str, str]:
    """Conservative connected groups, including multi-antigen bridge entries."""
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]

    for row in meta.values():
        fs = sorted(families(row["ag_cluster"]))
        for f in fs:
            parent[find(f)] = find(fs[0])
    return {i: find(min(families(r["ag_cluster"])))
            for i, r in meta.items() if families(r["ag_cluster"])}


def inst_of_tid(tid: str) -> str:
    return tid.split("/")[-1].removesuffix(".cif")


def auc(pos: list[float], neg: list[float]) -> float | None:
    """P(positive outranks negative); ties count half."""
    if not pos or not neg:
        return None
    wins = 0.0
    for p in pos:
        for n in neg:
            wins += 1.0 if p > n else 0.5 if p == n else 0.0
    return wins / (len(pos) * len(neg))


def read_tsv(path: Path):
    with path.open() as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def parse_target(mres: str) -> list[tuple[str, int]] | None:
    """Target tokens -> [(chain, resnum)], dropping '_'. None if unparseable."""
    out = []
    for t in mres.split(","):
        t = t.strip()
        if not t or t == UNMATCHED:
            continue
        if len(t) < 2 or not t[1:].isascii() or not t[1:].isdigit():
            return None
        try:
            out.append((t[0], int(t[1:])))
        except ValueError:
            return None
    return out


# ------------------------------------------------------------------ panel

def select_panel(root: Path, n: int, split: str, lo: int, hi: int,
                 min_epi: int, seed: int) -> list[dict]:
    """One single-chain antigen per family, chosen without looking at scores."""
    rows = [r for r in read_tsv(root / "structures.tsv")
            if r["ab_ag_split"] == split
            and "/" not in r["antigen_chains"]
            and len(families(r["ag_cluster"])) == 1
            and lo <= int(r["n_antigen_residues"]) <= hi
            and int(r["n_epitope_residues"]) >= min_epi]
    all_meta = {r["instance"]: r for r in read_tsv(root / "structures.tsv")}
    components = family_components(all_meta)
    by_family: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_family[components[r["instance"]]].append(r)
    rng = random.Random(seed)
    fams = sorted(by_family)
    rng.shuffle(fams)
    return [rng.choice(by_family[f]) for f in fams[:n]]


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def file_identity(path: Path) -> dict:
    stat = path.stat()
    return {"path": str(path.resolve()), "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns, "sha256": digest(path)}


def cached_run(out_dir: Path, signature: dict) -> bool:
    """A success marker AND unchanged artifacts are required; zero hits is valid."""
    try:
        record = json.loads((out_dir / "completed.json").read_text())
        return record["signature"] == signature and all(
            digest(out_dir / name) == checksum
            for name, checksum in record["artifacts"].items()
        ) and bool(record["artifacts"])
    except (OSError, ValueError, KeyError, TypeError):
        return False


def prepare_reference(root: Path, out: Path, rows: list[dict]) -> Path:
    """Build an isolated, content-identified index with absolute source paths."""
    out.mkdir(parents=True, exist_ok=True)
    signature = {"structures": {r["instance"]: file_identity(root / r["antigen_file"])
                                for r in rows},
                 "version": subprocess.check_output(["folddisco", "version"], text=True).strip(),
                 "hash_type": "default", "distance_bins": 0, "angle_bins": 0}
    index = out / "index"
    if cached_run(out, signature):
        return index
    (out / "completed.json").unlink(missing_ok=True)
    structures = out / "structures"
    structures.mkdir(exist_ok=True)
    for old in structures.iterdir():
        if not old.is_symlink():
            raise ValueError(f"refusing to remove non-symlink {old}")
        old.unlink()
    for r in rows:
        (structures / f"{r['instance']}.cif").symlink_to((root / r["antigen_file"]).resolve())
    cmd = ["folddisco", "index", "-p", str(structures.resolve()), "-i", str(index.resolve()),
           "-t", "1", "--id", "abspath", "-y", "default", "-d", "0", "-a", "0"]
    run = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    (out / "build.log").write_text(run.stdout + run.stderr)
    (out / "command.json").write_text(json.dumps(cmd, indent=2))
    run.check_returncode()
    artifacts = {p.name: digest(p) for p in
                 [index, Path(str(index) + ".lookup"), Path(str(index) + ".type"),
                  Path(str(index) + ".offset")]}
    (out / "completed.json").write_text(json.dumps(
        {"signature": signature, "artifacts": artifacts}, indent=2))
    return index

# --------------------------------------------------------- scan and search

def scan_and_search(antigen: Path, index: Path, out_dir: Path,
                    radius: float, top: int, query_limit: int = 0) -> bool:
    """Scan + search. Raw outputs are KEPT; a finished antigen is reused."""
    out_dir.mkdir(parents=True, exist_ok=True)
    hits = out_dir / "hits.tsv"
    signature = {"input": file_identity(antigen), "radius": radius, "top": top,
                 "query_limit": query_limit,
                 "index": [file_identity(p) for p in
                           [index, Path(str(index) + ".lookup"), Path(str(index) + ".type"),
                            Path(str(index) + ".offset")]],
                 "version": subprocess.check_output(["folddisco", "version"], text=True).strip(),
                 "code": {str(p.relative_to(REPO)): digest(p) for p in
                          sorted((REPO / "src" / "ibc_search").rglob("*.py"))},
                 "defaults": digest(REPO / "src/ibc_search/data/default.toml"),
                 "diagnostic": digest(Path(__file__)), "distance": 0.5, "angle": 5.0}
    if cached_run(out_dir, signature):
        return True
    (out_dir / "completed.json").unlink(missing_ok=True)
    env = {**os.environ, "PYTHONPATH": str(REPO / "src")}
    scan = subprocess.run(
        [str(REPO / ".venv/bin/python"), "-m", "ibc_search", "scan", str(antigen),
         "--radius", str(radius), "-o", str(out_dir), "--log-level", "ERROR"],
        capture_output=True, text=True, env=env, timeout=120)
    (out_dir / "scan.log").write_text(scan.stdout + scan.stderr)
    queries = out_dir / "folddisco_queries.tsv"
    if scan.returncode != 0 or not queries.is_file() or queries.stat().st_size == 0:
        print(f"    scan failed rc={scan.returncode}: {scan.stderr.strip()[-200:]!r}")
        return False
    search_queries = out_dir / "search_queries.tsv"
    lines = queries.read_text().splitlines()
    search_queries.write_text("\n".join(lines[:query_limit] if query_limit else lines) + "\n")
    cmd = ["folddisco", "query", "-i", str(index.resolve()), "-q", str(search_queries.resolve()),
           "-t", "1", "--top", str(top), "-d", "0.5", "-a", "5.0", "--header",
           "--per-match", "--format-output",
           "tid,node_count,idf,rmsd,query_residues,matching_residues"]
    (out_dir / "search_command.json").write_text(json.dumps(cmd, indent=2))
    pending = out_dir / "hits.pending.tsv"
    with pending.open("w") as fh:
        q = subprocess.run(
            cmd, stdout=fh, stderr=subprocess.PIPE, text=True, timeout=600)
    (out_dir / "search.log").write_text(q.stderr or "")
    if q.returncode != 0:
        print(f"    folddisco rc={q.returncode}: {(q.stderr or '').strip()[-300:]!r}")
        return False
    pending.replace(hits)
    names = ["hits.tsv", "patch_members.tsv", "folddisco_queries.index.tsv", "search_queries.tsv"]
    (out_dir / "completed.json").write_text(json.dumps(
        {"signature": signature, "artifacts": {n: digest(out_dir / n) for n in names}}, indent=2))
    return True


def load_hits(path: Path) -> tuple[list[dict], int]:
    """Parsed hit rows; rows whose field count disagrees with a header are torn."""
    rows, torn, hdr = [], 0, None
    with path.open() as fh:
        for line in fh:
            f = line.rstrip("\r\n").split("\t")
            if f and f[0] == "tid":
                hdr = f
                continue
            if hdr is None or len(f) != len(hdr):
                torn += 1
                continue
            r = dict(zip(hdr, f))
            try:
                rows.append({
                    "inst": inst_of_tid(r["tid"]),
                    "node_count": int(r["node_count"]),
                    "idf": float(r["idf"]),
                    "rmsd": float(r["rmsd"]),
                    "q": r["query_residues"],
                    "m": r["matching_residues"],
                })
            except (ValueError, KeyError):
                torn += 1
    return rows, torn


# ----------------------------------------------------------------- labels

def load_labels(root: Path, needed: set[str]):
    """lab[(inst, chain, res)] -> list of [ins, is_epitope, is_exposed]."""
    lab: dict[tuple, list[list]] = defaultdict(list)
    idx: dict[tuple, list] = {}
    for r in read_tsv(root / "epitope_residues.tsv"):
        if r["instance"] in needed:
            rec = [r["ins_code"], int(r["is_epitope"]), None]
            lab[(r["instance"], r["chain"], int(r["res_id"]))].append(rec)
            idx[(r["instance"], r["chain"], int(r["res_id"]), r["ins_code"])] = rec
    for r in read_tsv(root / "residue_sasa.tsv"):
        if r["instance"] in needed:
            rec = idx.get((r["instance"], r["chain"], int(r["res_id"]), r["ins_code"]))
            if rec is not None:
                rec[2] = int(r["is_exposed"])
    return lab


def base_rates(lab, insts: set[str]) -> dict[str, float | None]:
    exp: dict[str, int] = defaultdict(int)
    epi: dict[str, int] = defaultdict(int)
    for (inst, _, _), recs in lab.items():
        if inst not in insts:
            continue
        for _, e, x in recs:
            if x:
                exp[inst] += 1
                epi[inst] += e
    return {i: (epi[i] / exp[i] if exp[i] else None) for i in insts}


# --------------------------------------------------------------- analysis

def classify_match(row: dict, nq: int, lab, max_rmsd: float):
    """-> (reason or None, exposed_total, exposed_epitope). None = valid."""
    if row["node_count"] != nq:
        return "partial" if row["node_count"] < nq else "excess_nodes", 0, 0
    tgt = parse_target(row["m"])
    if tgt is None:
        return "unparseable", 0, 0
    if len(tgt) != len(set(tgt)):
        return "repeated_target", 0, 0
    if len(tgt) != nq:
        return "target_count_mismatch", 0, 0
    if not math.isfinite(row["rmsd"]) or row["rmsd"] < 0:
        return "nonfinite_or_negative_rmsd", 0, 0
    if not math.isfinite(row["idf"]) or row["idf"] < 0:
        return "nonfinite_or_negative_idf", 0, 0
    if row["rmsd"] > max_rmsd:
        return "rmsd", 0, 0
    n_exp = n_epi = 0
    for chain, res in tgt:
        recs = lab.get((row["inst"], chain, res))
        if not recs:
            return "unresolved_label", 0, 0
        if len(recs) > 1:
            return "insertion_code_ambiguous", 0, 0
        ins, e, x = recs[0]
        if ins:
            return "insertion_code_unexpressible", 0, 0
        if e not in (0, 1) or x not in (0, 1):
            return "unresolved_label", 0, 0
        if x:
            n_exp += 1
            n_epi += e
    return None, n_exp, n_epi


def analyse_query(rows, nq, own, meta, lab, base, budget, args, counters):
    """Score one motif at one candidate budget. Returns a dict or None."""
    # These rows come from an actual --top=budget search, NOT a simulated cap.
    kept = {r["inst"] for r in rows}

    own_pdb = meta[own]["pdb_id"]
    own_fams = families(meta[own]["ag_cluster"])
    retained = set()
    for i in kept:
        m = meta.get(i)
        if m is None or i == own or m["pdb_id"] == own_pdb \
                or own_fams & families(m["ag_cluster"]) \
                or (m.get("family_group") is not None and
                    m.get("family_group") == meta[own].get("family_group")) \
                or m["ab_ag_split"] != "train" \
                or "/" in m["antigen_chains"] or len(families(m["ag_cluster"])) != 1:
            continue
        retained.add(i)

    best: dict[str, tuple] = {}
    n_rows = n_valid = 0
    for r in rows:
        if r["inst"] not in retained:
            continue
        n_rows += 1
        reason, n_exp, n_epi = classify_match(r, nq, lab, args.max_rmsd)
        if reason is None and n_exp >= args.min_exposed:
            n_valid += 1
            cur = best.get(r["inst"])
            if cur is None or r["rmsd"] < cur[0]:
                best[r["inst"]] = (r["rmsd"], n_epi / n_exp)
        else:
            counters[reason or "too_few_exposed"] += 1

    per_family: dict[str, list[tuple[float, float]]] = defaultdict(list)
    for i, (_, frac) in best.items():
        b = base.get(i)
        if b is not None:
            per_family[meta[i].get("family_group", family_of(meta[i]["ag_cluster"]))].append((frac, frac - b))
    fam_frac = [st.mean(f for f, _ in v) for v in per_family.values()]
    fam_enr = [st.mean(e for _, e in v) for v in per_family.values()]
    positive_idfs = [r["idf"] for r in rows if r["inst"] in retained
                     and math.isfinite(r["idf"]) and r["idf"] > 0]
    return {
        "kept": len(kept), "retained": len(retained), "rows": n_rows,
        "valid_rows": n_valid, "valid_structures": len(best),
        "families": len(fam_enr),
        "contact_fraction": st.mean(fam_frac) if fam_frac else None,
        "enrichment": st.mean(fam_enr) if fam_enr else None,
        "idf_log10_max": math.log10(max(positive_idfs)) if positive_idfs else None,
        "valid_row_fraction": n_valid / n_rows if n_rows else None,
    }


SCORE_FIELDS = ("enrichment", "contact_fraction", "idf_log10_max", "valid_row_fraction")


def evaluation_subset(qrows, inst, score, min_families, supported_only, paired):
    """Common support AND complete scores define the paired comparison cohort."""
    return [q for q in qrows if q["antigen"] == inst and q[score] is not None
            and (not supported_only or q["families"] >= min_families)
            and (not paired or (q["families"] >= min_families
                                and all(q[s] is not None for s in SCORE_FIELDS)))]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("index_dir", type=Path, help="sabdab2_antigen_index directory")
    ap.add_argument("out_dir", type=Path, help="kept raw outputs and reports")
    ap.add_argument("-n", "--families", type=int, default=12)
    ap.add_argument("--split", default="train", help="ab_ag_split to sample [train]")
    ap.add_argument("--min-residues", type=int, default=50)
    ap.add_argument("--max-residues", type=int, default=250)
    ap.add_argument("--min-epitope", type=int, default=8)
    ap.add_argument("--radius", type=float, default=6.0)
    ap.add_argument("--top", type=int, default=1000,
                    help="maximum allowed candidate budget [1000]")
    ap.add_argument("--budgets", default="50,200,1000",
                    help="actual candidate budgets searched independently, comma-separated")
    ap.add_argument("--query-limit", type=int, default=0,
                    help="pilot only: search first N emitted queries, 0 means all")
    ap.add_argument("--primary-budget", type=int, default=None,
                    help="budget used for AUC [largest of --budgets]")
    ap.add_argument("--max-rmsd", type=float, default=1.5)
    ap.add_argument("--min-exposed", type=int, default=2)
    ap.add_argument("--min-families", type=int, default=3,
                    help="families required to call a motif 'supported' [3]")
    ap.add_argument("--pos-coverage", type=float, default=0.5)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(argv)

    try:
        budgets = sorted({int(b) for b in args.budgets.split(",")})
    except ValueError:
        ap.error("--budgets must contain integers")
    primary = args.primary_budget or budgets[-1]
    if (not budgets or min(budgets) < 1 or primary not in budgets
            or args.families < 1 or args.query_limit < 0 or args.min_families < 1
            or args.min_exposed < 1 or args.min_residues < 1
            or args.max_residues < args.min_residues
            or not math.isfinite(args.radius) or args.radius <= 0
            or not math.isfinite(args.max_rmsd) or args.max_rmsd < 0
            or not 0 < args.pos_coverage <= 1):
        ap.error("invalid bounds, thresholds, or primary budget")
    if args.top < max(budgets):
        print(f"error: --top {args.top} < largest budget {max(budgets)}", file=sys.stderr)
        return 2

    root, out = args.index_dir.resolve(), args.out_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    meta = {r["instance"]: r for r in read_tsv(root / "structures.tsv")}
    groups = family_components(meta)
    for inst, row in meta.items():
        row["family_group"] = groups.get(inst, f"unknown:{inst}")

    panel = select_panel(root, args.families, args.split, args.min_residues,
                         args.max_residues, args.min_epitope, args.seed)
    if not panel:
        ap.error("no eligible panel antigens")
    held_groups = {meta[r["instance"]]["family_group"] for r in panel}
    held_pdbs = {r["pdb_id"] for r in panel}
    reference = [r for r in meta.values() if r["ab_ag_split"] == "train"
                 and "/" not in r["antigen_chains"]
                 and len(families(r["ag_cluster"])) == 1
                 and r["family_group"] not in held_groups and r["pdb_id"] not in held_pdbs]
    print(f"building/checking isolated reference: {len(reference)} training instances", flush=True)
    if not reference:
        ap.error("no eligible training reference after group exclusion")
    index = prepare_reference(root, out / "reference", reference)
    with (out / "panel.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["instance", "pdb_id", "ag_cluster", "n_residues", "n_epitope"])
        for r in panel:
            w.writerow([r["instance"], r["pdb_id"], r["ag_cluster"],
                        r["n_antigen_residues"], r["n_epitope_residues"]])
    print(f"panel: {len(panel)} antigens (split={args.split}, seed={args.seed})")

    # ---- phase 1: scan + search, raw outputs kept --------------------------
    done = []
    related_controls = {}
    for n, r in enumerate(panel, 1):
        inst = r["instance"]
        print(f"  [{n}/{len(panel)}] {inst}", flush=True)
        success = True
        for budget in budgets:
            success = scan_and_search(root / r["antigen_file"], index,
                                      out / inst / f"budget_{budget}",
                                      args.radius, budget, args.query_limit) and success
        # Exact-coordinate control uses a dedicated one-structure index.
        self_index = prepare_reference(root, out / inst / "self_index", [r])
        success = scan_and_search(root / r["antigen_file"], self_index,
                                  out / inst / "self_control", args.radius, 1,
                                  args.query_limit) and success
        related = sorted((m for m in meta.values()
                          if m["family_group"] == meta[inst]["family_group"]
                          and m["pdb_id"] != r["pdb_id"]
                          and m["ab_ag_split"] == "train"
                          and "/" not in m["antigen_chains"]
                          and len(families(m["ag_cluster"])) == 1),
                         key=lambda m: m["instance"])[:3]
        if related:
            related_index = prepare_reference(root, out / inst / "related_index", related)
            related_ok = scan_and_search(root / r["antigen_file"], related_index,
                                         out / inst / "related_control", args.radius,
                                         len(related), args.query_limit)
            related_controls[inst] = {"status": "completed" if related_ok else "failed",
                                      "reference_instances": [m["instance"] for m in related]}
            success = related_ok and success
        else:
            related_controls[inst] = {"status": "unavailable", "reference_instances": []}
        if success:
            done.append(inst)
    if not done:
        print("error: no antigen completed", file=sys.stderr)
        return 1

    # ---- phase 2: labels ----------------------------------------------------
    hits_by = {}
    torn_total = 0
    needed = set(done)
    for inst in done:
        for budget in budgets:
            rows, torn = load_hits(out / inst / f"budget_{budget}" / "hits.tsv")
            hits_by[(inst, budget)] = rows
            torn_total += torn
            needed |= {r["inst"] for r in rows}
        needed.update(related_controls[inst]["reference_instances"])
    print(f"loading labels for {len(needed)} structures...", flush=True)
    lab = load_labels(root, needed)
    base = base_rates(lab, needed)

    # ---- phase 3: analysis ---------------------------------------------------
    stage = defaultdict(int)
    rej: dict[int, dict[str, int]] = {b: defaultdict(int) for b in budgets}
    budget_stats = {b: {"queries": 0, "any_retained": 0, "supported": 0,
                        "retained_sum": 0} for b in budgets}
    self_stats = defaultdict(int)
    related_stats = {}
    qrows = []
    audit = (out / "match_audit.tsv").open("w", newline="")
    audit_writer = csv.writer(audit, delimiter="\t", lineterminator="\n")
    audit_writer.writerow(["query_antigen", "budget", "query", "target", "node_count",
                           "rmsd", "idf", "target_residues", "quality_rejection",
                           "exposed_count", "contact_count"])

    for inst in done:
        d = out / inst / f"budget_{primary}"
        pid_to_q = {r["patch_id"]: r["query_string"]
                    for r in read_tsv(d / "folddisco_queries.index.tsv")}
        members: dict[str, set] = defaultdict(set)
        n_patches = set()
        for r in read_tsv(d / "patch_members.tsv"):
            members[r["patch_id"]].add((r["auth_chain_id"], int(r["auth_seq_id"]),
                                        r["ins_code"]))
            n_patches.add(r["patch_id"])
        stage["motifs_extracted"] += len(n_patches)
        qmap = {q: members[p] for p, q in pid_to_q.items()}
        stage["motifs_emitted"] += len(qmap)
        selected = {line.split("\t", 1)[1] for line in
                    (d / "search_queries.tsv").read_text().splitlines()}
        qmap = {q: v for q, v in qmap.items() if q in selected}
        stage["motifs_searched"] += len(qmap)

        by_q = {b: defaultdict(list) for b in budgets}
        for b in budgets:
            for r in hits_by[(inst, b)]:
                by_q[b][r["q"]].append(r)
                nq = len(parse_target(r["q"]) or [])
                reason, nexp, nepi = classify_match(r, nq, lab, args.max_rmsd)
                if reason is None and nexp < args.min_exposed:
                    reason = "too_few_exposed"
                audit_writer.writerow([inst, b, r["q"], r["inst"], r["node_count"],
                                       r["rmsd"], r["idf"], r["m"], reason or "accepted",
                                       nexp, nepi])
        control_rows, control_torn = load_hits(out / inst / "self_control" / "hits.tsv")
        torn_total += control_torn
        control_by_q = defaultdict(list)
        for row in control_rows:
            control_by_q[row["q"]].append(row)
        related_by_q = defaultdict(list)
        info = related_controls[inst]
        related_stats[inst] = {**info, "queries": 0, "queries_with_valid_match": 0,
                               "queries_with_contact_overlap": 0}
        if info["status"] == "completed":
            related_rows, related_torn = load_hits(out / inst / "related_control" / "hits.tsv")
            torn_total += related_torn
            for row in related_rows:
                related_by_q[row["q"]].append(row)

        for qstr, residues in qmap.items():
            rows = by_q[primary].get(qstr, [])
            nq = len(residues)
            if rows:
                stage["motifs_with_any_hit"] += 1

            # own-label class (insertion-code aware; unresolved -> skipped)
            known = []
            for chain, res, ins in residues:
                m = [rec for rec in lab.get((inst, chain, res), []) if rec[0] == ins]
                if len(m) == 1:
                    known.append(m[0][1])
            if len(known) < nq:
                stage["motifs_label_unresolved"] += 1
                continue
            cov = sum(known) / nq
            cls = ("epitope" if cov >= args.pos_coverage
                   else "other_exposed" if cov == 0 else "straddle")
            if info["status"] == "completed":
                valid_related = [classify_match(row, nq, lab, args.max_rmsd)
                                 for row in related_by_q[qstr]]
                valid_related = [(nexp, nepi) for reason, nexp, nepi in valid_related
                                 if reason is None and nexp >= args.min_exposed]
                related_stats[inst]["queries"] += 1
                related_stats[inst]["queries_with_valid_match"] += bool(valid_related)
                related_stats[inst]["queries_with_contact_overlap"] += any(
                    nepi > 0 for _, nepi in valid_related)

            # self-retrieval control
            self_rows = [r for r in control_by_q[qstr] if r["inst"] == inst]
            self_stats["queries"] += 1
            if self_rows:
                self_stats["self_hit_present"] += 1
                want = {(c, n) for c, n, _ in residues}
                ok = False
                for r in self_rows:
                    t = parse_target(r["m"])
                    if (r["node_count"] == nq and r["rmsd"] < 0.1
                            and t is not None and set(t) == want and len(t) == nq):
                        ok = True
                        break
                self_stats["self_exact_recovery"] += ok

            for b in budgets:
                res = analyse_query(by_q[b].get(qstr, []), nq, inst, meta, lab, base,
                                    b, args, rej[b])
                bs = budget_stats[b]
                bs["queries"] += 1
                bs["retained_sum"] += res["retained"]
                bs["any_retained"] += res["retained"] > 0
                bs["supported"] += res["families"] >= args.min_families
                if b == primary:
                    stage["motifs_retained_after_exclusion"] += res["retained"] > 0
                    stage["motifs_with_valid_match"] += res["valid_structures"] > 0
                    stage["motifs_supported"] += res["families"] >= args.min_families
                    qrows.append({"antigen": inst, "query": qstr, "size": nq,
                                  "class": cls, "epitope_coverage": round(cov, 3), **res})

    audit.close()
    # ---- per-query table ------------------------------------------------------
    cols = ["antigen", "query", "size", "class", "epitope_coverage", "kept", "retained",
            "rows", "valid_rows", "valid_structures", "families", "contact_fraction",
            "enrichment", "idf_log10_max", "valid_row_fraction"]
    with (out / "per_query.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for q in qrows:
            w.writerow(["" if q[c] is None else
                        (f"{q[c]:.4f}" if isinstance(q[c], float) else q[c]) for c in cols])

    # ---- AUCs (within antigen, then averaged over antigens) ------------------
    def antigen_auc(score: str, supported_only: bool, paired: bool = False):
        per = []
        details = []
        for inst in done:
            sel = evaluation_subset(qrows, inst, score, args.min_families,
                                    supported_only, paired)
            a = auc([q[score] for q in sel if q["class"] == "epitope"],
                    [q[score] for q in sel if q["class"] == "other_exposed"])
            if a is not None:
                per.append(a)
            details.append({"antigen": inst, "auc": a,
                            "positive": sum(q["class"] == "epitope" for q in sel),
                            "other_exposed": sum(q["class"] == "other_exposed" for q in sel)})
        return {"antigens": len(per), "mean_auc": round(st.mean(per), 4) if per else None,
                "above_0.5": sum(a > 0.5 for a in per), "per_antigen": details}

    scores = {
        "epitope_enrichment (family-balanced)": ("enrichment", True),
        "contact_fraction (family-balanced)": ("contact_fraction", True),
        "baseline_idf_log10_max": ("idf_log10_max", False),
        "baseline_valid_row_fraction": ("valid_row_fraction", False),
    }
    auc_report = {k: antigen_auc(*v) for k, v in scores.items()}
    paired_report = {k: antigen_auc(*v, paired=True) for k, v in scores.items()}

    cls_counts = defaultdict(int)
    for q in qrows:
        cls_counts[q["class"]] += 1
    summary = {
        "settings": {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()},
        "antigens_searched": len(done), "torn_rows_dropped": torn_total,
        "failed_antigens": [r["instance"] for r in panel if r["instance"] not in done],
        "reference_instances": len(reference),
        "reference_partition": "train, single-chain, one known cluster; panel connected groups/PDBs excluded",
        "budget_method": "independent actual searches",
        "self_control_index": "dedicated exact-coordinate singleton per antigen",
        "related_controls": related_stats,
        "stage_accounting": dict(stage),
        "motif_classes_at_primary_budget": dict(cls_counts),
        "self_retrieval": dict(self_stats),
        "budget_sensitivity": {
            str(b): {**s, "mean_retained_structures":
                     round(s["retained_sum"] / s["queries"], 2) if s["queries"] else None}
            for b, s in budget_stats.items()},
        "match_rejections_by_budget": {str(b): dict(rej[b]) for b in budgets},
        "auc_within_antigen": auc_report,
        "auc_paired_same_motifs": paired_report,
        "not_implemented": ["family bootstrap confidence intervals",
                            "matched-composition surface background", "residue-level transfer"],
        "note": ("'other_exposed' motifs are NOT confirmed non-epitopes. A motif with "
                 "no usable match has no evidence and is excluded from AUC, not scored 0."),
    }
    (out / "diagnostic_summary.json").write_text(json.dumps(summary, indent=2))

    print("\n=== self-retrieval control ===")
    sq = self_stats["queries"] or 1
    print(f"  queries {self_stats['queries']}, self hit present "
          f"{self_stats['self_hit_present']} ({100 * self_stats['self_hit_present'] / sq:.1f}%), "
          f"exact residue recovery {self_stats['self_exact_recovery']} "
          f"({100 * self_stats['self_exact_recovery'] / sq:.1f}%)")
    print("=== stage accounting (primary budget) ===")
    for k, v in stage.items():
        print(f"  {k:<36} {v}")
    print("=== candidate-budget sensitivity ===")
    for b, s in budget_stats.items():
        print(f"  budget {b:>5}: mean retained structures "
              f"{s['retained_sum'] / max(s['queries'], 1):.1f}; motifs with any "
              f"{s['any_retained']}/{s['queries']}; supported {s['supported']}")
    print("=== within-antigen AUC, epitope vs other exposed ===")
    for k, v in auc_report.items():
        print(f"  {k:<40} {v}")
    print(f"\nwrote {out / 'diagnostic_summary.json'} and per_query.tsv")
    return 0 if len(done) == len(panel) else 1


if __name__ == "__main__":
    sys.exit(main())
