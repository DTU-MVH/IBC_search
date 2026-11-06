# true_annotated_ppi

Dataset root:
`/Users/mvh_mac/documents_local_only/IBC_search/investigation_databases/true_annotated_ppi`

150 individual target proteins, each from a distinct experimental biological
assembly-1 source in the old PPI database. That database is not overwritten.
The inherited source filter is X-ray resolution <=2.5 A. No UniProt reviewed
requirement is imposed: labels are observed protein contacts, not UniProt features.

## Interface and selection

A positive target residue has an amino-acid heavy atom within 5.0 A inclusive of
an eligible partner protein amino-acid heavy atom. Water, non-protein ligands,
hydrogens and deuterium do not contribute. Partner identity and minimum distance
are preserved per residue pair. All touching eligible partner chains are retained
for labelling, but partner atoms are never exported to the target query file.

Both target and each partner have >=80 resolved amino-acid residues. Targets are
limited to <=1500 residues. This explicit size restriction excludes peptide-only
partners but does not independently prove a folded domain or interaction. In each
entry, the longest eligible contacting chain is selected, ties by source chain ID.
Homomeric as well as heteromeric contacts are allowed.

Source entries are ordered by a stable SHA256(`true-ppi-150:` + PDB ID) key. One
target is taken per source entry, and identical resolved target residue sequences
are excluded. No new family clustering is performed, so family leakage remains
a downstream concern. Antibody complexes are excluded by membership in the local
SAbDab2 metadata and antibody entity-name screening; this is not guaranteed to
recognize every unclassified antibody in all deposits.

## Residue identity correction

Every source residue key includes chain, residue number, insertion code and
residue identity. Unlike the old PPI builder, insertion codes are never merged.
Query structures are deliberately remapped to **chain A, residue numbers 1..N,
blank insertion code** for FoldDisco. Coordinates and atom/residue identities are
unchanged. `residue_mapping.tsv` is the reversible source-to-query mapping; use it
when interpreting source PDB positions. Every residue label also carries source
chain, number and insertion code. Modified amino acids remain identified.

## Files

- `structures/`: target-only heavy-atom amino-acid mmCIF query files.
- `source_complexes/`: retained full source assemblies, NOT query structures.
- `annotated_residues.tsv`: positive query residue annotations.
- `residue_labels.tsv`: all resolved target residues with contact state and
  minimum distance to selected eligible partners.
- `residue_mapping.tsv`: complete reversible query/source residue identities.
- `contacts.tsv`: per-partner residue contact pairs; target identifiers are query
  mapped, partner identifiers are source chain/number/insertion code.
- `panel.tsv`, `structures.tsv`, `provenance.json`, `build_manifest.json`,
  `validation.json`, `REPORT.md`.

Zero contact state means no observed contact to these partners, not proof of
nonfunction. Short-chain partners, unresolved regions, water/glycan-mediated
contacts and other conformations remain unassessed. These are bound conformations,
not experimentally isolated apo structures.

"True annotated" here means reproducible observed contacts in supplied biological
assemblies. It is NOT independent curated confirmation of biological relevance
or energetic hotspots. The old manifest's asymmetric-unit caveat is stale; actual
assembly-1 files and provenance are retained rather than silently rewriting it.

Build:

```sh
/Users/mvh_mac/documents_local_only/IBC_search/.venv/bin/python /Users/mvh_mac/documents_local_only/IBC_search/scripts/build_true_annotated_ppi.py
```

Final audit:

```sh
/Users/mvh_mac/documents_local_only/IBC_search/.venv/bin/python /Users/mvh_mac/documents_local_only/IBC_search/scripts/build_true_annotated_ppi.py --validate-only
```

The audit reconstructs target exports and tables from retained source assemblies,
checks full coordinate/identity equality and hashes, mapping and membership,
duplicate targets, and independently recomputes positive labels using dense
distances instead of KD-tree neighbor search. No FoldDisco searches are run.