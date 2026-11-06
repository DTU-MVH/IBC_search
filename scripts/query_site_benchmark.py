"""Annotation-blind whole-query pilot: evaluate seed-centred FoldDisco scores.

prepare chooses proteins by ID hash and size, never site labels; run_panel.py
performs unseeded dumb-search. evaluate uses labels only after feature extraction.
No target biological annotations or fitted parameters enter the scores.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
import time

import numpy as np
from scipy.stats import mannwhitneyu

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "src"))
from count_scoring import split_patch_ids  # noqa: E402
from top_motifs_by_organism import pdb_id_from_tid  # noqa: E402
from ibc_search.folddisco_contract import parse_hit_residues, parse_query_token  # noqa: E402
from ibc_search.errors import ContractError  # noqa: E402

RUNS = {
    "zincfingers": "zincfingers_seeded_r5.5",
    "enzymes": "enzymes_seeded2_r5.5",
    "protein_protein_interactions": "protein_protein_interactions_seeded2_r5.5",
}
FEATURES = ("n_sites", "n_targets", "n_targets_cov80", "quality_support",
            "best_quality", "mean_coverage", "negative_median_rmsd",
            "size", "seed_rsa")


def rows(path):
    with Path(path).open() as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def write_table(path, header, values):
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(header); writer.writerows(values)


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def prepare(args):
    if (args.out / "panel.tsv").exists():
        raise ValueError("pilot panel already exists; use evaluate or a new output directory")
    args.out.mkdir(parents=True, exist_ok=True)
    chosen = []; candidates = {}
    for category, run in RUNS.items():
        eligible = []
        for item in rows(REPO / "results" / run / "panel.tsv"):
            source = REPO / "results" / run / item["instance"]
            count = sum(r["is_scanned_chain"] == "true" and r["residue_class"].startswith("polymer_")
                        for r in rows(source / "residues.tsv"))
            if args.min_residues <= count <= args.max_residues:
                rank = hashlib.sha256(f"query-pilot-{args.seed}:{item['instance']}".encode()).hexdigest()
                eligible.append((rank, item, count))
        candidates[category] = len(eligible)
        for _, item, count in sorted(eligible, key=lambda x: x[0])[:args.per_category]:
            chosen.append({"instance": item["instance"],
                           "structure_file": str((REPO / item["structure_file"]).resolve()),
                           "category": category, "previous_scanned_residues": count})
    write_table(args.out / "panel.tsv", ["instance", "structure_file", "category"],
                [[r[k] for k in ["instance", "structure_file", "category"]] for r in chosen])
    write_json(args.out / "pilot_manifest.json", {
        "purpose": "development feasibility pilot, NOT a held-out performance estimate",
        "selection": "SHA256(query-pilot-seed:ID), size-restricted; no site labels or hit scores used",
        "seed": args.seed, "size_range": [args.min_residues, args.max_residues],
        "eligible": candidates, "selected": chosen,
        "planned_search": {"method": "dumb-search", "seed_residues": None,
                           "radius_angstrom": 5.5, "top": 40, "num_residue": 50000,
                           "threads_per_process": 1},
    })
    print("prepared", args.out / "panel.tsv", flush=True)


def feature_rows(directory, instance):
    """No labels read here. One score vector for the ball's seed residue.

    Sites deduplicated per query. For quality_support each target protein
    contributes its maximum coverage*exp(-RMSD/2). No RMSD/coverage cutoff.
    Non-injective, malformed rows and same-entry hits are excluded. Compact
    unaligned target fields are usable: only sets and counts, never positional
    target-to-query correspondence, enter these features.
    """
    manifest = json.loads((directory / "run_manifest.json").read_text())
    if manifest["ball_search"].get("seeded", False):
        raise ValueError(f"annotation-guided scan cannot be evaluated: {directory}")
    residues = {}
    for r in rows(directory / "residues.tsv"):
        if r["is_scanned_chain"] == "true" and r["residue_class"].startswith("polymer_"):
            token = r["auth_chain_id"] + r["auth_seq_id"]
            if r["ins_code"]:
                continue
            if token in residues:
                raise ValueError(f"ambiguous seed token {token}")
            residues[token] = r
    patches = {r["patch_id"]: r for r in rows(directory / "patches.tsv")}
    sites = defaultdict(dict); audit = Counter()
    hits = directory / "folddisco_hits.tsv"
    if not hits.exists():
        raise ValueError(f"missing search output: {hits}")
    for h in rows(hits):
        audit["rows_total"] += 1
        target = pdb_id_from_tid(h["tid"])
        if target is None:
            audit["unknown_target"] += 1; continue
        if target == instance.lower():
            audit["self_excluded"] += 1; continue
        parsed = parse_hit_residues(h["query_residues"], h["matching_residues"])
        if not parsed.positionally_alignable:
            audit["compact_target_field"] += 1
        if not parsed.query_tokens or len(parsed.matched_tokens) > len(parsed.query_tokens):
            audit["invalid_lengths"] += 1; continue
        if not parsed.matched_tokens or len(set(parsed.matched_tokens)) != len(parsed.matched_tokens):
            audit["noninjective_or_empty"] += 1; continue
        try:
            rmsd = float(h["rmsd"])
            for t in parsed.matched_tokens:
                parse_query_token(t)
        except (ValueError, ContractError):
            audit["invalid"] += 1; continue
        if not math.isfinite(rmsd) or rmsd < 0:
            audit["invalid"] += 1; continue
        coverage = len(parsed.matched_tokens) / len(parsed.query_tokens)
        site = (target, tuple(sorted(parsed.matched_tokens)))
        for pid in split_patch_ids(h["patch_id"]):
            if pid not in patches or patches[pid]["emit_query"] != "true":
                raise ValueError(f"unknown/unsearched patch: {pid}")
            if patches[pid]["query_string"].split(",") != list(parsed.query_tokens):
                raise ValueError(f"query mismatch: {pid}")
            previous = sites[pid].get(site)
            if previous is None or rmsd < previous[1]:
                sites[pid][site] = (coverage, rmsd)
        audit["rows_usable"] += 1
    result = []
    for pid, p in patches.items():
        seed = p["seed_residue_labels"]
        if seed not in residues:
            continue
        emitted = p["emit_query"] == "true"
        features = {k: None for k in FEATURES}
        features["size"] = int(p["n_members"])
        features["seed_rsa"] = float(residues[seed]["rsa"]) if residues[seed]["rsa"] else None
        if emitted:
            found = sites.get(pid, {})
            target_quality = defaultdict(float); good_targets = set()
            for (target, _), (cov, rmsd) in found.items():
                target_quality[target] = max(target_quality[target], cov * math.exp(-rmsd / 2))
                if cov >= .8:
                    good_targets.add(target)
            features.update(n_sites=len(found), n_targets=len(target_quality),
                            n_targets_cov80=len(good_targets), quality_support=sum(target_quality.values()),
                            best_quality=max(target_quality.values(), default=0.),
                            mean_coverage=float(np.mean([v[0] for v in found.values()])) if found else 0.,
                            # -1e6 is an explicit below-all-valid-hits sentinel, not a measured RMSD.
                            negative_median_rmsd=-float(np.median([v[1] for v in found.values()])) if found else -1e6)
        result.append({"instance": instance, "seed": seed, "patch_id": pid,
                       "assessed": emitted, **features})
    return result, dict(audit), len(residues)


def ranking_metrics(labels, scores, fraction=.1):
    """Tie-aware AP and expected top-k recovery; no arbitrary ID tie breaking."""
    y = np.asarray(labels, dtype=int); s = np.asarray(scores, dtype=float)
    if not len(y) or not y.sum() or y.sum() == len(y):
        return None
    n = len(y); positives = int(y.sum()); auc = float(mannwhitneyu(s[y == 1], s[y == 0]).statistic / (positives * (n-positives)))
    tp = seen = 0; ap = 0.
    for value in sorted(set(s), reverse=True):
        group = y[s == value]; new = int(group.sum()); seen += len(group); tp += new
        ap += new / positives * tp / seen
    k = max(1, math.ceil(n * fraction)); threshold = np.sort(s)[-k]
    above = s > threshold; tied = s == threshold
    expected_tp = float(y[above].sum() + (k-int(above.sum())) * y[tied].mean())
    return {"n": n, "positives": positives, "prevalence": positives/n,
            "auc": auc, "average_precision": ap, "ap_over_prevalence": ap/(positives/n),
            "top_fraction": fraction, "top_k": k, "expected_precision_at_k": expected_tp/k,
            "expected_recall_at_k": expected_tp/positives,
            "enrichment_at_k": (expected_tp/k)/(positives/n)}


def evaluate(args):
    panel = list(rows(args.out / "panel.tsv")); annotations = {}
    for cat in RUNS:
        annotations[cat] = defaultdict(set)
        for r in rows(REPO / "investigation_databases" / cat / "annotated_residues.tsv"):
            if cat != "enzymes" or r["ligand_class"] == "biological":
                annotations[cat][r["pdb_id"]].add(r["chain"] + r["res_id"])
    table = []; evaluations = []; audits = {}; coverage = []
    # Extract features from all structures first, without passing site labels.
    for item in panel:
        inst = item["instance"]; directory = args.out / inst
        done = directory / "search_complete.json"
        if not done.exists():
            raise ValueError(f"missing explicit completion marker: {done}")
        expected = json.loads(done.read_text())["hits_sha256"]
        if hashlib.sha256((directory / "folddisco_hits.tsv").read_bytes()).hexdigest() != expected:
            raise ValueError(f"hits changed since completion: {directory}")
        f, audit, nres = feature_rows(directory, inst)
        audits[inst] = audit
        positives = annotations[item["category"]][inst]
        for r in f:
            r.update(category=item["category"], annotated=int(r["seed"] in positives))
        table.extend(f)
        coverage.append({"instance": inst, "category": item["category"], "parsed_seed_residues": nres,
                         "scored_seed_balls": sum(r["assessed"] for r in f),
                         "positive_seeds": sum(r["annotated"] for r in f),
                         "annotated_tokens_in_table": len(positives & {r['seed'] for r in f}),
                         "annotation_tokens_unmapped": sorted(positives - {r['seed'] for r in f})})
        for size_range in [(1, 1000000), (4, 7)]:
            eligible = [r for r in f if r["assessed"] and size_range[0] <= r["size"] <= size_range[1]]
            for feature in FEATURES:
                valid = [r for r in eligible if r[feature] is not None]
                metrics = ranking_metrics([r["annotated"] for r in valid], [r[feature] for r in valid])
                evaluations.append({"instance": inst, "category": item["category"],
                                    "size_range": list(size_range), "feature": feature, "metrics": metrics})
    header = ["instance", "category", "seed", "patch_id", "annotated", "assessed", *FEATURES]
    write_table(args.out / "query_scores.tsv", header, [[r.get(k) for k in header] for r in table])
    groups = defaultdict(list)
    for r in evaluations:
        if r["metrics"]:
            groups[(r["category"], tuple(r["size_range"]), r["feature"])].append(r["metrics"])
    summary = [{"category": k[0], "size_range": list(k[1]), "feature": k[2], "structures": len(v),
                **{name: float(np.mean([m[name] for m in v])) for name in
                   ["auc", "average_precision", "ap_over_prevalence", "enrichment_at_k", "expected_recall_at_k"]}}
               for k,v in groups.items()]
    write_json(args.out / "query_evaluation.json", {
        "status": "development pilot; no family holdout; no learned model or calibrated threshold",
        "label_definition": "seed belongs to known annotation; enzyme biological ligand contacts, NOT verified catalysis; unlabelled residues are annotation negatives, not proven nonfunctional",
        "feature_definition": "quality_support = sum across target proteins of max(coverage*exp(-RMSD/2)); no hard RMSD cutoff; all valid coverages retained",
        "coverage": coverage, "audits": audits, "per_structure": evaluations, "macro_summary": summary,
        "limitations": ["Top-40 truncation; counts are candidate-budget dependent.",
                        "Only same-entry self hits excluded; homologues remain.",
                        "Seed annotations are incomplete; absence is not proof of no function.",
                        "Small size-restricted development sample; no confidence intervals or generalization claim.",
                        "Full-background comparison and family-held-out validation remain to be completed.",
                        "Within-protein baseline comparisons are descriptive; no incremental predictive model fitted."]})
    lines = ["# Whole-query localization pilot", "", "Development feasibility only; no held-out validation.", "",
             "Scores are computed without site labels; evaluation compares annotated seeds against other seeds in the SAME protein.", "",
             "| Category | Sizes | Score | Proteins | Mean AUC | Mean AP/prevalence | Top-10% enrichment |", "|---|---|---|---:|---:|---:|---:|"]
    for r in summary:
        lines.append(f"| {r['category']} | {r['size_range']} | {r['feature']} | {r['structures']} | {r['auc']:.3f} | {r['ap_over_prevalence']:.2f} | {r['enrichment_at_k']:.2f} |")
    lines += ["", "See query_evaluation.json for exclusions, coverage, per-protein results, and limitations."]
    (args.out / "REPORT.md").write_text("\n".join(lines)+"\n")
    print("evaluated", args.out / "REPORT.md", flush=True)


def finish(args):
    """Wait for the pilot runner's successful join logs, then evaluate.

    A bounded wait avoids orphaned, indefinite pollers after a failed search.
    Successful shell join log lines are stronger evidence than file existence.
    """
    panel = list(rows(args.out / "panel.tsv"))
    deadline = time.monotonic() + args.wait_seconds
    while True:
        complete = []
        for item in panel:
            log = args.out / (item["instance"] + ".log")
            text = log.read_text() if log.exists() else ""
            if "establishes nothing about antibody binding or immunogenicity." in text:
                complete.append(item)
        if len(complete) == len(panel):
            break
        if time.monotonic() >= deadline:
            raise TimeoutError(f"only {len(complete)}/{len(panel)} searches completed; no evaluation produced")
        time.sleep(5)
    for item in complete:
        directory = args.out / item["instance"]
        write_json(directory / "search_complete.json", {
            "evidence": str(args.out / (item["instance"] + ".log")),
            "verification": "successful shell search join and final reminder logged",
            "hits_sha256": hashlib.sha256((directory / "folddisco_hits.tsv").read_bytes()).hexdigest(),
        })
    evaluate(args)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("stage", choices=["prepare", "evaluate", "finish"])
    ap.add_argument("--out", type=Path, default=REPO / "results/query_site_pilot")
    ap.add_argument("--per-category", type=int, default=2)
    ap.add_argument("--min-residues", type=int, default=60)
    ap.add_argument("--max-residues", type=int, default=150)
    ap.add_argument("--seed", type=int, default=5)
    ap.add_argument("--wait-seconds", type=int, default=3600)
    args = ap.parse_args(); args.out = args.out.resolve()
    if args.per_category < 1 or args.min_residues < 1 or args.max_residues < args.min_residues:
        ap.error("invalid pilot selection limits")
    {"prepare": prepare, "evaluate": evaluate, "finish": finish}[args.stage](args)


if __name__ == "__main__":
    main()