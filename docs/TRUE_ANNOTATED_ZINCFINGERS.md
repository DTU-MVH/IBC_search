# true_annotated_zincfingers

Dataset root:
`/Users/mvh_mac/documents_local_only/IBC_search/investigation_databases/true_annotated_zincfingers`

Target: 100 unique reviewed UniProt/Swiss-Prot proteins with explicit Zn binding
residue annotations inside explicit zinc-finger domains. This database is
separate from the old `zincfingers` dataset and never imports its distance labels.

## Strict source rule

1. The source entry must be `UniProtKB reviewed (Swiss-Prot)`.
2. A domain must have feature type `Zinc finger` and exact, valid boundaries.
3. A positive residue must be an exact single-residue `Binding site` feature
   whose ligand identifier is `ChEBI:CHEBI:29105` (Zn2+).
4. That residue must fall inside exactly one accepted zinc-finger domain.
   Ambiguous overlapping domain assignments and duplicate residue labels reject
   the protein. Binding features outside zinc-finger domains are not included.
5. Every source domain and binding feature, including evidence, is retained.

No whole-domain positives, distance-derived labels, unrelated metal sites or
generic ligand pockets are introduced. This rule does not require complete
coordination shells; missing explicit residue annotations remain unknown.
Domain containment is an operational domain assignment, not an independent
experimental proof of metal coordination.

## Evidence interpretation

The requested name means **explicit reviewed-source annotations**, NOT that
every residue has been experimentally verified. UniProt evidence can include
homology and sequence-model inference or no feature-level evidence. Both domain
and binding evidence are exported. See `validation.json` for evidence counts.
Unlabelled residues are not established negatives.

## Structures and selection

AlphaFold DB predicted canonical full-length monomers are used only when the
sequence exactly matches UniProt. Every residue identity and number is checked;
no unverified PDB numbering transfer is performed. Models lack zinc atoms and
assembly context. No confidence threshold is applied: retain model confidence
and control for it downstream. Sequence accession order and availability drive
selection, not random sampling or sequence clustering. Family leakage remains
a downstream concern.

## Outputs and reproducibility

`panel.tsv`, `structures.tsv`, `annotated_residues.tsv`, `sequences.fasta`,
`provenance.json`, `structures/`, `build_manifest.json`, `validation.json`, and
`REPORT.md`. Join by `instance`/UniProt accession, not an invented PDB ID.
Raw responses, release headers and checksums are cached. Annotation tables retain
domain boundaries, domain type and per-residue evidence.

Build (resumes from cache if interrupted; finished datasets cannot be overwritten):

```sh
/Users/mvh_mac/documents_local_only/IBC_search/.venv/bin/python /Users/mvh_mac/documents_local_only/IBC_search/scripts/build_true_annotated_zincfingers.py
```

Re-run the final audit:

```sh
/Users/mvh_mac/documents_local_only/IBC_search/.venv/bin/python /Users/mvh_mac/documents_local_only/IBC_search/scripts/build_true_annotated_zincfingers.py --validate-only
```

The audit reproduces all labels from saved reviewed records, verifies model
checksums/sequences/numbering and checks annotation and panel membership. A quota
shortfall fails the build rather than relaxing labels or admitting unreviewed
entries. No FoldDisco searches are run.