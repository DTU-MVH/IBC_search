# Query-side whole-structure benchmark

`scripts/query_site_benchmark.py` prepares a small development panel and
evaluates annotation-blind, every-residue FoldDisco searches. This is distinct
from target-site localization: the prediction unit is a residue of the query,
scored only by its own centred motif. Annotations are not passed to search or
feature extraction. They enter evaluation only.

## Pilot

The default panel has two proteins per category, selected by a seeded SHA256
ordering of entry IDs among structures with 60–150 previously scanned amino-acid
residues. The restriction is computational, not a representative population
sample. No hit scores or annotated-residue positions enter selection.

Use the project's Python environment and absolute paths rooted at
`/Users/mvh_mac/documents_local_only/IBC_search`:

```sh
/Users/mvh_mac/documents_local_only/IBC_search/.venv/bin/python /Users/mvh_mac/documents_local_only/IBC_search/scripts/query_site_benchmark.py prepare
IBC_FOLDDISCO_NUM_RESIDUE=50000 /Users/mvh_mac/documents_local_only/IBC_search/.venv/bin/python /Users/mvh_mac/documents_local_only/IBC_search/scripts/run_panel.py /Users/mvh_mac/documents_local_only/IBC_search/results/query_site_pilot --method dumb-search --radius 5.5 --top 40 -j 3 --no-finalize
```

Before evaluating, each structure requires `search_complete.json`, an explicit
record of successful search completion. Do not create these for interrupted or
failed searches merely because a hits file exists. The `finish` command waits
up to one hour for successful shell search join logs, writes completion markers
with hit-file checksums, and evaluates. It can be started while searches run:

```sh
/Users/mvh_mac/documents_local_only/IBC_search/.venv/bin/python /Users/mvh_mac/documents_local_only/IBC_search/scripts/query_site_benchmark.py finish
```

## Scores and evaluation

- Distinct target sites and PDB entries; entries with a coverage ≥0.8 hit.
- Best and mean coverage-dependent quality; negative median RMSD.
- `quality_support = sum_target max_site(coverage * exp(-RMSD / 2))`.
  This is a fixed exploratory formula, not fitted or calibrated.
- Size and seed RSA are query-only baselines. No hard RMSD cutoff is imposed.
- Scores include all valid coverage levels. Duplicate sites are collapsed;
  same-entry self matches are removed again at evaluation.
- Evaluate all emitted motifs and, separately, sizes 4–7. Unemitted motifs are
  unassessed, not zero-hit negatives.
- Within-protein AUC, tie-group average precision, and tie-aware expected
  top-10% recovery; category summaries average proteins equally.

Enzyme positives are contacts to ligands classified as biological, **not curated
catalytic residues**. Unlabelled residues are not proven biological negatives.

## What remains before a biological prediction claim

This pilot does not resolve top-40 censoring, homologous-target inflation,
family leakage, or structural/compositional confounding. It has no held-out
test set. Background proteins and larger annotated panels need identical new
scans. Depth sensitivity, target-family counts, within-size/composition controls,
protein/family-held-out comparisons against stronger query-only baselines, and
structure-level uncertainty are required before claiming generalization.

Existing unrelated working-tree changes are deliberately left untouched.