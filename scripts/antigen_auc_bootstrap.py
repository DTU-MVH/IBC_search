"""Offline, equal-antigen bootstrap of existing FoldDisco hit counts."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import mannwhitneyu

from count_scoring import site_key, split_patch_ids

ROOT = Path(__file__).resolve().parents[1]
THRESHOLDS = (50, 70, 80)


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_table(path, records):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]),
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def auc(annotated, background):
    a, b = np.asarray(annotated), np.asarray(background)
    if not len(a) or not len(b):
        raise ValueError("Both classes must be nonempty")
    return float(((a[:, None] > b).sum() + 0.5 * (a[:, None] == b).sum()) / (len(a) * len(b)))


def bootstrap(values, replicates=10000, seed=20261120):
    """Rows are antigens; columns are thresholds. Draw rows jointly."""
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or not len(values) or not np.isfinite(values).all():
        raise ValueError("Expected finite antigen-by-threshold matrix")
    draws = np.random.default_rng(seed).integers(0, len(values), (replicates, len(values)))
    return values[draws].mean(axis=1), draws


def counts(directory):
    patches = rows(directory / "patches.tsv")
    query_rows = rows(directory / "folddisco_queries.index.tsv")
    queried = {r["patch_id"] for r in query_rows}
    assert len(queried) == len(query_rows)
    assert queried == {r["patch_id"] for r in patches if r["emit_query"] == "true"}
    sites = {p: [set() for _ in THRESHOLDS] for p in queried}
    with (directory / "folddisco_hits.tsv").open() as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            idf = float(r["idf"])
            assert np.isfinite(idf) and idf >= 50
            ids = split_patch_ids(r["patch_id"])
            assert ids and set(ids) <= queried
            key = site_key(r["tid"], r["matching_residues"])
            if not key[1]:
                continue
            for p in ids:
                for j, threshold in enumerate(THRESHOLDS):
                    if idf >= threshold:
                        sites[p][j].add(key)
    result = []
    for p in patches:
        if p["patch_id"] not in queried:
            continue
        assert p["seed_class"] in ("annotated", "background")
        assert float(p["radius_angstrom"]) == 5.5
        values = [len(s) for s in sites[p["patch_id"]]]
        assert values == sorted(values, reverse=True)
        result.append(dict(instance=directory.name, patch_id=p["patch_id"],
                           seed_residue=p["seed_residue_labels"], seed_class=p["seed_class"],
                           **{f"idf_{t}": v for t, v in zip(THRESHOLDS, values)}))
    return result, len(patches) - len(queried)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, default=ROOT / "results/antibody_targets_fewballs_swissprot_r5.5")
    parser.add_argument("--out", type=Path, default=ROOT / "results/antibody_swissprot_antigen_auc_bootstrap")
    args = parser.parse_args()
    assert not args.out.exists(), "Use a new output directory; never overwrite an analysis"
    panel = rows(args.run / "panel.tsv")
    instances = sorted(r["instance"] for r in panel)
    assert len(instances) == len(set(instances)) == 70
    balls, per_antigen, matrix, unassessed = [], [], [], 0
    source_hashes = {}
    for inst in instances:
        directory = args.run / inst
        log = (args.run / (inst + ".log")).read_text()
        assert "queries searched:" in log and "wrote " in log
        for name in ("patches.tsv", "folddisco_queries.index.tsv", "folddisco_hits.tsv"):
            path = directory / name
            source_hashes[str(path.resolve())] = hashlib.sha256(path.read_bytes()).hexdigest()
        recs, skipped = counts(directory)
        balls.extend(recs)
        unassessed += skipped
        values = []
        for t in THRESHOLDS:
            a = [r[f"idf_{t}"] for r in recs if r["seed_class"] == "annotated"]
            b = [r[f"idf_{t}"] for r in recs if r["seed_class"] == "background"]
            value = auc(a, b)
            independent = mannwhitneyu(a, b, method="asymptotic").statistic / (len(a) * len(b))
            assert abs(value - independent) < 1e-12
            values.append(value)
            per_antigen.append(dict(instance=inst, min_idf=t, n_annotated=len(a), n_background=len(b),
                                   auc=value, annotated_median=float(np.median(a)), background_median=float(np.median(b))))
        matrix.append(values)
    original = rows(args.run / "01_statistics/annotated_vs_background.tsv")
    actual = {(r["instance"], r["seed_residue"]): r["idf_50"] for r in balls if r["seed_class"] == "annotated"}
    expected = {(r["instance"], r["seed_residue"]): int(r["n_distinct_sites"]) for r in original}
    assert len(actual) == len(expected) == len(original) == 1614
    assert actual == expected, "Baseline must exactly reproduce saved counts"
    matrix = np.array(matrix)
    boot, draws = bootstrap(matrix)
    summary, contrasts = [], []
    for j, t in enumerate(THRESHOLDS):
        a = np.array([r[f"idf_{t}"] for r in balls if r["seed_class"] == "annotated"])
        b = np.array([r[f"idf_{t}"] for r in balls if r["seed_class"] == "background"])
        assert np.median(a) == (46, 9, 5)[j]
        low, high = np.quantile(boot[:, j], [0.025, 0.975])
        summary.append(dict(min_idf=t, antigens=70, annotated_balls=len(a), background_balls=len(b),
                            annotated_median=float(np.median(a)), background_median=float(np.median(b)),
                            annotated_zero_fraction=float(np.mean(a == 0)), mean_auc=float(matrix[:, j].mean()),
                            ci95_low=float(low), ci95_high=float(high), illustrative_margin=0.60,
                            upper_bound_below_margin=bool(high < 0.60)))
        if j:
            delta = boot[:, j] - boot[:, 0]
            lo, hi = np.quantile(delta, [0.025, 0.975])
            contrasts.append(dict(comparison=f"IDF {t} minus IDF 50", mean_auc_change=float((matrix[:, j]-matrix[:, 0]).mean()),
                                  ci95_low=float(lo), ci95_high=float(hi)))
    args.out.mkdir(parents=True)
    write_table(args.out / "ball_counts.tsv", balls)
    write_table(args.out / "per_antigen_auc.tsv", per_antigen)
    write_table(args.out / "summary.tsv", summary)
    write_table(args.out / "paired_changes.tsv", contrasts)
    np.savez_compressed(args.out / "bootstrap.npz", draws=draws, mean_auc=boot, antigen_auc=matrix,
                        instances=np.array(instances), thresholds=np.array(THRESHOLDS))
    audit = dict(complete=True, source_run=str(args.run.resolve()), seed=20261120, replicates=10000,
                 weighting="equal antigen", interval="95% percentile bootstrap; pointwise, not simultaneous",
                 unassessed_patches=unassessed, baseline_1614_counts_exact=True, expected_medians_verified=True,
                 scipy_auc_crosscheck=True, source_sha256=source_hashes,
                 analysis_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (args.out / "validation.json").write_text(json.dumps(audit, indent=2) + "\n")
    text = ["# Antigen-level AUC bootstrap: Swiss-Prot antibody benchmark", "",
            "Retrospective analysis of saved searches; no new searches or parameter optimization.", "",
            "Score: distinct (target ID, matched target residue set) sites per emitted ball. Higher is predicted annotated.",
            "AUC compares annotated/background balls within an antigen with half credit for ties. All-zero antigens stay at 0.5.",
            "Mean AUC weights each of the 70 antigens equally. Bootstrap resamples whole antigens 10,000 times, seed 20261120.",
            "Identical antigen draws at all thresholds support paired differences. CIs are pointwise 95% percentile intervals.", "",
            "|IDF floor|Annotated median|BG median|Mean AUC|95% CI|Upper bound <0.60|",
            "|---:|---:|---:|---:|---|---|"]
    for r in summary:
        text.append(f"|{r['min_idf']}|{r['annotated_median']:g}|{r['background_median']:g}|{r['mean_auc']:.4f}|[{r['ci95_low']:.4f}, {r['ci95_high']:.4f}]|{r['upper_bound_below_margin']}|")
    text += ["", "## Paired changes", ""]
    for r in contrasts:
        text.append(f"- {r['comparison']}: {r['mean_auc_change']:+.4f}, 95% CI [{r['ci95_low']:+.4f}, {r['ci95_high']:+.4f}].")
    text += ["", "## Interpretation boundaries", "",
             "AUC 0.60 is an illustrative retrospective usefulness threshold, not a prespecified or biologically validated margin.",
             "An upper bound below 0.60 argues against that magnitude of average HIGHER-count advantage under this sampling model; it is not proof of impossibility or exact equivalence to chance.",
             "Intervals assume antigens are independent sampling units. Unmodeled family dependence and nonrandom benchmark selection limit generalization.",
             "This is the original unadjusted within-antigen comparison, not RSA- or motif-size-matched analysis. Background is not confirmed non-epitope.",
             "The scores are conditional on the original candidate retrieval limit; homology, label incompleteness, and other FoldDisco scores are not resolved here.",
             "Same observations at different IDF thresholds are sensitivity checks, not independent replications. Pointwise CIs do not account for threshold selection.",
             "These conclusions apply to average performance; individual proteins may differ. The result does not test a reversed score or all uses of FoldDisco.", "",
             "## Validation", "", f"All 1,614 baseline annotated counts exactly reproduced; medians 46/9/5 verified; {unassessed} unassessed patches excluded.",
             "Every per-antigen AUC cross-checked against SciPy Mann–Whitney U. Input hashes and bootstrap draws retained.", "",
             "## Reproduce", "", "```sh", f"{ROOT / '.venv/bin/python'} {Path(__file__).resolve()} --out /tmp/antigen_auc_bootstrap_reproduction", "```"]
    (args.out / "REPORT.md").write_text("\n".join(text) + "\n")
    print(json.dumps(dict(summary=summary, paired_changes=contrasts, unassessed=unassessed), indent=2))


if __name__ == "__main__":
    main()