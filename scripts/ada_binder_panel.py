"""Assemble the panel of designed ADA binders to scan.

One row per construct: the predicted structure to scan, the design's
provenance, its prediction confidence, and the ELISA absorbance measured
against pre-existing anti-drug antibodies.

The ELISA column is an ANNOTATION, not a label
----------------------------------------------
Nothing downstream fits to it, thresholds on it or selects by it. It is
carried so the scan's output can afterwards be compared against the measured
ADA reactivity -- which is the only experimental evidence in this dataset
about whether a binder is actually recognised. Keeping it in the panel rather
than joining it later means the comparison cannot silently drift out of
register with the structures.

Why ``best/`` and not ``notebook_selection/``
---------------------------------------------
``best/`` is all 78 constructs at their highest-pLDDT sample.
``notebook_selection/`` is the 73 that a previous analysis kept, and it is
built by a sort that selects the LOWEST-pLDDT sample per construct (its own
README records the bug). Since the request is every binder, the 78-structure
folder is the one that answers it, and taking the best sample removes a known
defect rather than reproducing it.

These are single-chain predictions of the binder alone, with no target and no
antibody. There is therefore NO observed epitope anywhere in this dataset:
every residue-level call the pipeline makes here is a prediction with no
ground truth to score it against.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import statistics as st
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

DATA = Path("/Users/mvh_mac/documents_local_only/master_thesis_data_files/"
            "experimental_data_minibinders_Christian")


def elisa_by_well(path: Path) -> dict[str, dict]:
    """Mean blank-subtracted absorbance per ADA well, over its replicates."""
    by_well: dict[str, list[float]] = defaultdict(list)
    for r in csv.DictReader(path.open()):
        well = (r.get("ada_well") or "").strip()
        raw = r.get("absorbance_450_blank_subtracted")
        if not well or raw in ("", None):
            continue
        try:
            by_well[well].append(float(raw))
        except ValueError:
            continue
    return {
        w: {
            "elisa_mean": round(st.mean(v), 4),
            "elisa_sd": round(st.pstdev(v), 4) if len(v) > 1 else "",
            "elisa_n": len(v),
        }
        for w, v in by_well.items()
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("--data-dir", type=Path, default=DATA)
    ap.add_argument("--structures", default="best",
                    choices=("best", "notebook_selection"))
    args = ap.parse_args(argv)

    struct_dir = args.data_dir / "predicted_structures" / args.structures
    files = sorted(struct_dir.glob("*.cif"))
    if not files:
        print(f"error: no .cif under {struct_dir}", file=sys.stderr)
        return 2

    elisa = elisa_by_well(args.data_dir / "ADA_wells.csv")
    constructs = {
        r["ada_well"]: r
        for r in csv.DictReader((args.data_dir / "ADA_constructs.csv").open())
    }
    confidence = {
        (r["ada_well"], r["sample"]): r
        for r in csv.DictReader(
            (args.data_dir / "predicted_structures" /
             "structure_confidence.csv").open()
        )
    }

    args.out_dir.mkdir(parents=True, exist_ok=True)
    panel = args.out_dir / "panel.tsv"
    n_elisa = 0
    with panel.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["instance", "ada_well", "construct_name", "mb_target",
                    "mb_source", "design_strategy", "seq_len", "cys_count",
                    "plddt", "elisa_mean", "elisa_sd", "elisa_n",
                    "structure_file"])
        for path in files:
            # <well>__<construct>__best_sample<N>.cif
            stem = path.stem
            well = stem.split("__", 1)[0]
            sample = (re.search(r"sample(\d+)$", stem) or [None, ""])[1]
            con = constructs.get(well, {})
            conf = confidence.get((well, sample), {})
            e = elisa.get(well, {})
            if e:
                n_elisa += 1
            w.writerow([
                stem, well, con.get("construct_name", ""),
                con.get("mb_target", ""), con.get("mb_source", ""),
                con.get("design_strategy", ""), con.get("seq_len", ""),
                con.get("cys_count", ""), conf.get("plddt", ""),
                e.get("elisa_mean", ""), e.get("elisa_sd", ""),
                e.get("elisa_n", ""), str(path.resolve()),
            ])

    (args.out_dir / "panel_manifest.json").write_text(json.dumps({
        "built_utc": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        "structure_source": str(struct_dir),
        "structures": len(files),
        "with_elisa": n_elisa,
        "ground_truth_epitopes": False,
        "note": ("Single-chain binder predictions with no target and no "
                 "antibody: no observed epitope exists in this dataset."),
    }, indent=2))

    print(f"{len(files)} structures, {n_elisa} with ELISA absorbance")
    print(f"  {panel}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
