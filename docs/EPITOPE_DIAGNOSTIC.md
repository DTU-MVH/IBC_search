# Quality-controlled SAbDab diagnostic

This is an exploratory **patch-level** benchmark, not a predictor of serum
binding and not a calibrated probability of an epitope. The existing IDF
scorer and its saved results are unchanged.

## Design

- Sample single-chain antigens with exactly one known antigen-cluster ID.
- Construct connected groups using **all** reference cluster memberships,
  including multichain bridge entries. Sample one query per connected group.
- Build a separate reference from training instances only, restricted to
  single-chain, single-cluster antigens. Exclude every selected query's PDB
  entry and connected group **before indexing**.
- Search each budget independently; do not approximate prefilter budgets
  from per-match output. `kept` counts structures with returned matches, not
  the unobserved prefilter candidate population.
- Use a singleton index per query for exact-coordinate self-retrieval.
- When available, index up to three training structures from different PDBs
  in the same connected group for a separate related-antigen control.
  Related matches are technical sensitivity controls, not held-out evidence.
- Retain only full-node matches with exactly as many unique, resolvable target
  residues as query residues, finite nonnegative scores/RMSD, and sufficient
  exposed target residues. Ambiguous and unexpressible insertion-code targets
  are rejected rather than guessed.
- Choose the lowest-RMSD valid match per target instance, average within
  connected reference groups, then average groups equally.
- Enrichment is target contact fraction minus that target instance's exposed
  contact base rate. This background does **not** match patch composition.
- Report method-specific coverage and paired AUCs using the same supported,
  complete-score query cohort. Other exposed surface is unlabelled, not a
  definitive negative epitope class. No evidence remains missing, not zero.

## Running

Use absolute paths for both positional arguments. A small pilot command:

```bash
~/documents_local_only/IBC_search/.venv/bin/python -B \
  ~/documents_local_only/IBC_search/scripts/epitope_diagnostic.py \
  ~/documents_local_only/master_thesis_data_files/sabdab2_antigen_index \
  ~/documents_local_only/IBC_search/results/epitope_diagnostic_pilot \
  -n 3 --budgets 10,50 --top 50 --query-limit 5
```

The first N emitted queries are a deterministic **technical** subset, not a
representative biological sample. Do not interpret its AUC as performance.
For complete surface scans omit `--query-limit`. Keep evaluation settings fixed
before selecting `--split test`; development searches never include test targets.

## Artifacts and resuming

The output directory retains index inputs as symlinks, index build commands and
logs, scan logs, exact search commands, raw hits for every actual budget, and
completion records. Completion records check input/index contents, extraction
radius, query limit, budget, Python source/defaults, FoldDisco version and output
checksums. Header-only successful searches can be reused; incomplete outputs
without matching completion records cannot. Hashing index contents on resume
is intentionally conservative and has a runtime cost.

`diagnostic_summary.json` contains stage counts, failed antigens, self-control
recovery, related-control availability/recovery, rejection counts and AUCs.
`per_query.tsv` includes evidence counts and scores. `match_audit.tsv` records
quality decisions for each parsed primary/reference hit at every budget; these
are quality checks, not residue-level correspondences. Raw malformed rows remain
in the original hit files and are counted, not used.

## Remaining limits

No positional residue-label transfer, antibody placement, family-bootstrap
intervals, matched-composition background, size-corrected IDF baseline, or
experimental minibinder localization is implemented. A related match having
contact overlap is not proof of equivalent antibody specificity. The reference
is deliberately narrower than the full SAbDab collection. Overlapping query
patches are dependent; motif counts are not independent sample sizes.

The first live smoke test recovered exact residues for 15/15 queried motifs.
Five had a full quality-qualified family-held-out match and none reached three
reference families. These small, first-five-query results validate plumbing,
not prediction quality. Read the latest saved summary for subsequent results.