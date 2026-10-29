"""Project motif scores onto residues, and test them against known epitopes.

A motif score says something about a *set* of residues. To ask whether the
score finds antibody binding sites, it has to be expressed per residue, so
each residue takes the **maximum** z-score of any motif containing it.

Maximum, not mean: a residue belonging to one strongly-scoring motif and
several weak ones is implicated by the strong one. Averaging would let the
weak majority bury it, which is the wrong direction for a method whose
whole premise is that a few motifs carry the signal.

Residues in no emitted motif get no score and are excluded from the
comparison rather than scored zero -- they were never assessed, and scoring
them zero would silently stack the deck in whichever direction the epitope
happens to lie.

Runs with no epitopes at all
----------------------------
A designed binder predicted as a monomer has no antibody in the file, so no
residue can be labelled by contact and there is no positive class. That is a
normal state, not a missing input: the script detects it, writes the same
per-residue scores, and reports call counts only. **It does not report an
AUC, a precision or a recall**, because each of those needs a positive class
to be defined, and a number computed against an empty one would look like
performance while measuring nothing.
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path


def load_motif_residues(results_dir: Path) -> dict[str, list[tuple[str, int]]]:
    """patch_id -> its residues, from the scan."""
    out: dict[str, list[tuple[str, int]]] = defaultdict(list)
    path = results_dir / "patch_members.tsv"
    if not path.is_file():
        return out
    for r in csv.DictReader(path.open(), delimiter="\t"):
        out[r["patch_id"]].append((r["auth_chain_id"], int(r["auth_seq_id"])))
    return out


def load_motif_scores(results_dir: Path) -> dict[str, float]:
    """query string -> z, joined to patch ids through the query index."""
    scored = results_dir / "idf_scored.tsv"
    index = results_dir / "folddisco_queries.index.tsv"
    if not (scored.is_file() and index.is_file()):
        return {}
    by_query = {}
    for r in csv.DictReader(scored.open(), delimiter="\t"):
        try:
            by_query[r["query_residues"]] = float(r["idf_z"])
        except (ValueError, KeyError):
            continue
    out: dict[str, float] = {}
    for r in csv.DictReader(index.open(), delimiter="\t"):
        z = by_query.get(r["query_string"])
        if z is not None:
            # A query string can map to several patches; each gets the score.
            out[r["patch_id"]] = z
    return out


def residue_scores(results_dir: Path) -> dict[tuple[str, int], float]:
    members = load_motif_residues(results_dir)
    scores = load_motif_scores(results_dir)
    best: dict[tuple[str, int], float] = {}
    for patch_id, z in scores.items():
        for key in members.get(patch_id, ()):
            if key not in best or z > best[key]:
                best[key] = z
    return best


def auc(pos: list[float], neg: list[float]) -> float | None:
    """Probability a positive outranks a negative (ties count half)."""
    if not pos or not neg:
        return None
    merged = sorted([(v, 1) for v in pos] + [(v, 0) for v in neg])
    ranks: dict[int, float] = {}
    i = 0
    while i < len(merged):
        j = i
        while j < len(merged) and merged[j][0] == merged[i][0]:
            j += 1
        r = (i + j + 1) / 2
        for k in range(i, j):
            ranks[k] = r
        i = j
    rp = sum(ranks[k] for k, (_, g) in enumerate(merged) if g == 1)
    return (rp - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_dir", type=Path, help="directory holding panel.tsv")
    ap.add_argument("--epitopes", type=Path, default=None)
    ap.add_argument("--threshold", type=float, default=1.7,
                    help="z at or above which a residue is called [1.7]")
    ap.add_argument("-o", "--out", type=Path, default=None)
    args = ap.parse_args(argv)

    # No epitope file is a legitimate state, not an error: a designed binder
    # predicted as a monomer has no antibody in the file, so no residue can
    # be labelled by contact. The scores are then PREDICTIONS WITH NO GROUND
    # TRUTH, and the summary must say so rather than report an AUC computed
    # against an empty positive class.
    epi_path = args.epitopes or next(iter(sorted(
        args.run_dir.glob("epitopes_*.tsv"))), None)
    epitope: dict[tuple[str, str, int], int] = {}
    if epi_path is not None:
        for r in csv.DictReader(epi_path.open(), delimiter="\t"):
            epitope[(r["instance"], r["chain"], int(r["res_id"]))] = \
                int(r["is_epitope"])
    labelled = bool(epitope)

    out_path = args.out or (args.run_dir / "residue_scores.tsv")
    per_structure = []
    with out_path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["instance", "chain", "res_id", "is_epitope", "max_motif_z",
                    "called_at_threshold"])
        for row in csv.DictReader((args.run_dir / "panel.tsv").open(), delimiter="\t"):
            inst = row["instance"]
            scores = residue_scores(args.run_dir / inst)
            if not scores:
                continue
            pos, neg = [], []
            tp = fp = fn = tn = 0
            for (chain, res_id), z in sorted(scores.items()):
                lab = epitope.get((inst, chain, res_id))
                if lab is None and labelled:
                    continue  # residue not in the labelled antigen chains
                called = int(z >= args.threshold)
                w.writerow([inst, chain, res_id,
                            "" if lab is None else lab, f"{z:+.4f}", called])
                (pos if lab else neg).append(z)
                if lab and called: tp += 1
                elif lab: fn += 1
                elif called: fp += 1
                else: tn += 1
            if not labelled:
                per_structure.append({
                    "instance": inst, "n_scored": len(neg),
                    "n_called": fp, "call_rate": round(fp / len(neg), 4),
                    "mean_z": round(st.mean(neg), 4),
                    "max_z": round(max(neg), 4),
                })
                continue
            a = auc(pos, neg)
            if a is not None:
                per_structure.append({
                    "instance": inst, "n_scored": len(pos) + len(neg),
                    "n_epitope": len(pos), "auc": a,
                    "epitope_mean_z": st.mean(pos), "other_mean_z": st.mean(neg),
                    "tp": tp, "fp": fp, "fn": fn, "tn": tn,
                })

    if not per_structure:
        print("error: nothing scored", file=sys.stderr)
        return 2

    if not labelled:
        scored = sum(p["n_scored"] for p in per_structure)
        called = sum(p["n_called"] for p in per_structure)
        hit = [p for p in per_structure if p["n_called"] > 0]
        summary = {
            "ground_truth_epitopes": False,
            "note": ("No antibody in these structures: every call is a "
                     "prediction with nothing to score it against. No AUC, "
                     "precision or recall is defined here."),
            "structures": len(per_structure),
            "residues_scored": scored,
            "threshold": args.threshold,
            "residues_called": called,
            "call_rate": round(called / scored, 4) if scored else None,
            "structures_with_a_call": len(hit),
            "structures_with_no_call": len(per_structure) - len(hit),
            "mean_z": round(st.mean([p["mean_z"] for p in per_structure]), 4),
            "median_calls_per_structure": st.median(
                [p["n_called"] for p in per_structure]),
            "max_calls_in_one_structure": max(
                p["n_called"] for p in per_structure),
        }
        (args.run_dir / "residue_score_summary.json").write_text(json.dumps(
            {"summary": summary, "per_structure": per_structure}, indent=2))
        print(f"wrote {out_path}")
        for k, v in summary.items():
            print(f"  {k:<28} {v}")
        return 0

    aucs = [p["auc"] for p in per_structure]
    tp = sum(p["tp"] for p in per_structure); fp = sum(p["fp"] for p in per_structure)
    fn = sum(p["fn"] for p in per_structure); tn = sum(p["tn"] for p in per_structure)
    base = (tp + fn) / (tp + fp + fn + tn)
    prec = tp / (tp + fp) if tp + fp else float("nan")
    summary = {
        "structures": len(per_structure),
        "residues_scored": tp + fp + fn + tn,
        "epitope_residues": tp + fn,
        "epitope_base_rate": round(base, 4),
        "threshold": args.threshold,
        "auc_mean": round(st.mean(aucs), 4),
        "auc_median": round(st.median(aucs), 4),
        "auc_above_0.5": sum(1 for a in aucs if a > 0.5),
        "auc_below_0.5": sum(1 for a in aucs if a < 0.5),
        "precision_at_threshold": round(prec, 4),
        "recall_at_threshold": round(tp / (tp + fn), 4) if tp + fn else None,
        "enrichment_over_base_rate": round(prec / base, 3) if base else None,
        "epitope_mean_z": round(st.mean([p["epitope_mean_z"] for p in per_structure]), 4),
        "other_mean_z": round(st.mean([p["other_mean_z"] for p in per_structure]), 4),
    }
    (args.run_dir / "residue_score_summary.json").write_text(
        json.dumps({"summary": summary, "per_structure": per_structure}, indent=2)
    )

    print(f"wrote {out_path}")
    for k, v in summary.items():
        print(f"  {k:<28} {v}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
