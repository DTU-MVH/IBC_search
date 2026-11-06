# Reviewed UniProt enzyme database

Location: `/Users/mvh_mac/documents_local_only/IBC_search/investigation_databases/enzymes_reviewed_catalytic`.

This is separate from the original `enzymes` ligand-pocket dataset; nothing in
that dataset is overwritten. Target: 100 unique reviewed UniProt accessions in
each of nitrogenases, serine proteases, and metalloproteases. Nitrogenases are
individual subunit proteins, not 100 assembled nitrogenase complexes.

## Labels

- **Serine proteases:** exactly one explicit three-residue UniProt Active site
  charge-relay system, with Ser/His/Asp or Ser/His/Glu composition. Dyads,
  incomplete triads and multiple catalytic domains are excluded. No other
  residue features are annotated.
- **Nitrogenases:** exact single-residue Binding site features for nitrogenase
  catalytic metal clusters, P-clusters, or Fe-protein electron-transfer clusters.
  These three roles have distinct label types; electron-transfer ligands are
  NOT called catalytic residues or triads. ATP-binding and general pockets are
  excluded. Membership requires a nitrogenase EC assignment, not just a name.
- **Metalloproteases:** exact single-residue Active site features; additionally
  Binding site features explicitly marked with ligand note `catalytic`. Structural
  metals, substrate pockets and other ligand contacts are excluded.

Reviewed does not mean experimentally verified. Every annotation retains its
UniProt evidence, including homology/PROSITE inference or missing evidence.
The literature listed in the manifest supports the nitrogenase cofactor role;
it does not experimentally establish every residue in every selected homologue.
Full UniProt records, source features and AlphaFold metadata are preserved in
each class's `provenance.json`. No new sequence-based annotation transfer is done.

Nitrogenase context references: PMID 2500438 (Fe-protein cluster ligands), PMID
2287275 (MoFe alpha-chain conserved residues), DOI 10.1371/journal.pone.0072751
(MoFe P-cluster and FeMo-cofactor ligation). UniProt feature definitions:
`https://www.uniprot.org/help/sequence_annotation`.

## Structures and selection

Predicted AlphaFold DB full-length canonical monomer models are downloaded
through its API. The sequence must exactly match the reviewed UniProt entry;
every model residue identity and number is checked, not just annotated residues.
No PDB-to-UniProt alignment is assumed. These models lack cofactors and assembled
inter-subunit geometry. They are not substitutes for experimentally resolved
active complexes. No pLDDT threshold is applied; confidence is retained in models
and metadata and should be controlled in downstream analysis.

Selection is accession ascending among eligible entries with available matching
models; this is not random, taxonomically balanced or sequence-clustered. Family
redundancy must be addressed in any later benchmark. Download failures and label
exclusions are counted. A quota shortfall is reported as a failed/incomplete build,
not filled with unreviewed proteins or expanded labels.

Each class contains `panel.tsv`, `structures.tsv`, `annotated_residues.tsv`,
`sequences.fasta`, `provenance.json`, and `structures/`. The annotation join key
is `instance` (UniProt accession), not `pdb_id`: do not invent PDB IDs for models.
The existing panel runner accepts `instance` plus `structure_file` columns.
No FoldDisco searches are launched by the builder.

Reproduce with:

```sh
/Users/mvh_mac/documents_local_only/IBC_search/.venv/bin/python /Users/mvh_mac/documents_local_only/IBC_search/scripts/build_reviewed_enzymes.py --count 100 --jobs 6
```

Raw response caches enable interrupted runs to resume. A completed output cannot
be overwritten; choose a new `--out` for a different build. Release identifiers,
query URLs and hashes appear in `build_manifest.json`.

Run the independent final source/model audit with:

```sh
/Users/mvh_mac/documents_local_only/IBC_search/.venv/bin/python /Users/mvh_mac/documents_local_only/IBC_search/scripts/validate_reviewed_enzymes.py /Users/mvh_mac/documents_local_only/IBC_search/investigation_databases/enzymes_reviewed_catalytic
```

It writes `validation.json` only after reproducing all labels from the saved
reviewed source records and rechecking every model checksum, sequence, numbering,
class count and accession uniqueness.