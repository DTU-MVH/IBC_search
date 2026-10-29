"""Sweep the call threshold against labelled residues, and report where to put it.

Picking a cutoff by eye on an unlabelled run is circular: whatever it is, some
fraction of residues clears it. This reads a scored run that HAS labels and
reports, at every candidate threshold, what the call would actually buy --
precision, recall, F1, Matthews correlation and enrichment over the base rate.

Enrichment is the number to read first. A threshold whose precision merely
equals the base rate has bought nothing: it selects residues at the same rate
chance would. Below 1.0x it is worse than guessing.

The sweep also reports the AUC, which is threshold-free. **If the AUC is at
chance, no row in the table is a good threshold** -- the ranking carries no
signal, and precision at any cut is then base rate plus noise. The table
still prints, because seeing precision sit flat at the base rate across the
whole range is the clearest possible demonstration of that.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from residue_scores import auc  # noqa: E402


def sweep_row(pos: list[float], neg: list[float], t: float, base: float) -> dict:
    tp = sum(1 for v in pos if v >= t)
    fn = len(pos) - tp
    fp = sum(1 for v in neg if v >= t)
    tn = len(neg) - fp
    prec = tp / (tp + fp) if tp + fp else None
    rec = tp / (tp + fn) if tp + fn else None
    f1 = (2 * prec * rec / (prec + rec)
          if prec and rec and (prec + rec) else None)
    denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = ((tp * tn - fp * fn) / denom) if denom else None
    return {
        "threshold": round(t, 2),
        "called": tp + fp,
        "call_rate": round((tp + fp) / (len(pos) + len(neg)), 4),
        "tp": tp, "fp": fp, "fn": fn, "tn": tn,
        "precision": round(prec, 4) if prec is not None else None,
        "recall": round(rec, 4) if rec is not None else None,
        "f1": round(f1, 4) if f1 is not None else None,
        "mcc": round(mcc, 4) if mcc is not None else None,
        "enrichment": round(prec / base, 3) if prec is not None else None,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--scores", type=Path, default=None)
    ap.add_argument("--label", default="positive",
                    help="what is_epitope==1 means in this run, for the report")
    ap.add_argument("--start", type=float, default=-1.0)
    ap.add_argument("--stop", type=float, default=3.0)
    ap.add_argument("--step", type=float, default=0.1)
    ap.add_argument("-o", "--out", type=Path, default=None)
    args = ap.parse_args(argv)

    scores = args.scores or (args.run_dir / "residue_scores.tsv")
    if not scores.is_file():
        print(f"error: no residue scores at {scores}", file=sys.stderr)
        return 2

    pos: list[float] = []
    neg: list[float] = []
    per_inst: dict[str, tuple[list[float], list[float]]] = defaultdict(
        lambda: ([], []))
    for r in csv.DictReader(scores.open(), delimiter="\t"):
        try:
            z = float(r["max_motif_z"])
        except (KeyError, ValueError):
            continue
        if r.get("is_epitope") == "1":
            pos.append(z)
            per_inst[r["instance"]][0].append(z)
        elif r.get("is_epitope") == "0":
            neg.append(z)
            per_inst[r["instance"]][1].append(z)
    if not pos or not neg:
        print("error: this run has no labelled residues to sweep against",
              file=sys.stderr)
        return 2

    base = len(pos) / (len(pos) + len(neg))
    overall_auc = auc(pos, neg)
    per_auc = {i: auc(p, n) for i, (p, n) in per_inst.items()}
    per_auc = {i: a for i, a in per_auc.items() if a is not None}

    rows = []
    t = args.start
    while t <= args.stop + 1e-9:
        rows.append(sweep_row(pos, neg, t, base))
        t += args.step

    best_f1 = max((r for r in rows if r["f1"] is not None),
                  key=lambda r: r["f1"], default=None)
    best_mcc = max((r for r in rows if r["mcc"] is not None),
                   key=lambda r: r["mcc"], default=None)
    best_enr = max((r for r in rows if r["enrichment"] is not None
                    and r["called"] >= 10),
                   key=lambda r: r["enrichment"], default=None)

    summary = {
        "label_meaning": args.label,
        "residues": len(pos) + len(neg),
        "positive": len(pos),
        "base_rate": round(base, 4),
        "auc": round(overall_auc, 4) if overall_auc is not None else None,
        "auc_per_structure_mean": round(st.mean(per_auc.values()), 4) if per_auc else None,
        "auc_per_structure_above_half": sum(1 for a in per_auc.values() if a > 0.5),
        "structures": len(per_auc),
        "positive_mean_z": round(st.mean(pos), 4),
        "negative_mean_z": round(st.mean(neg), 4),
        "best_f1": best_f1,
        "best_mcc": best_mcc,
        "best_enrichment_min10_calls": best_enr,
    }

    out = args.out or (args.run_dir / "threshold_sweep.json")
    out.write_text(json.dumps(
        {"summary": summary, "sweep": rows,
         "per_structure_auc": {k: round(v, 4) for k, v in sorted(per_auc.items())}},
        indent=2))

    tsv = out.with_suffix(".tsv")
    with tsv.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(list(rows[0]))
        for r in rows:
            w.writerow([r[k] for k in rows[0]])

    print(f"{summary['residues']} residues, {summary['positive']} positive "
          f"({args.label}), base rate {base:.3f}")
    print(f"AUC {summary['auc']}   per-structure mean "
          f"{summary['auc_per_structure_mean']}  "
          f"({summary['auc_per_structure_above_half']}/{summary['structures']} "
          f"above 0.5)")
    print(f"positive mean z {summary['positive_mean_z']:+.3f} vs "
          f"negative {summary['negative_mean_z']:+.3f}")
    print()
    print("     z   called   call%   prec   recall      F1     MCC   enrich")
    for r in rows:
        if abs(r["threshold"] * 10 - round(r["threshold"] * 10)) < 1e-6 and \
                round(r["threshold"] * 10) % 2 == 0:
            print(f"  {r['threshold']:+5.1f}  {r['called']:>6}  "
                  f"{100 * r['call_rate']:5.1f}%  "
                  f"{(r['precision'] if r['precision'] is not None else float('nan')):.3f}  "
                  f"{(r['recall'] if r['recall'] is not None else float('nan')):.3f}  "
                  f"{(r['f1'] if r['f1'] is not None else float('nan')):6.3f}  "
                  f"{(r['mcc'] if r['mcc'] is not None else float('nan')):+6.3f}  "
                  f"{(r['enrichment'] if r['enrichment'] is not None else float('nan')):6.3f}")
    print(f"\nbest F1  at {best_f1['threshold']:+.1f} (F1 {best_f1['f1']}, "
          f"enrichment {best_f1['enrichment']}x)")
    print(f"best MCC at {best_mcc['threshold']:+.1f} (MCC {best_mcc['mcc']})")
    if best_enr:
        print(f"best enrichment at {best_enr['threshold']:+.1f} "
              f"({best_enr['enrichment']}x on {best_enr['called']} calls)")
    print(f"\nwrote {out}\nwrote {tsv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
