# Experimental antibody target interfaces

Dataset:
`/Users/mvh_mac/documents_local_only/IBC_search/investigation_databases/true_annotated_antibody_targets`

70 experimentally resolved target-only structures, derived from the local SAbDab2
`abag_split.csv` and cropped complex files. These are not AlphaFold predictions
and are not required to be reviewed UniProt entries. "True annotated" denotes
observed structural interfaces, not curated energetic hotspots.

Selection requires holo X-ray complexes at resolution <=3.0 A, protein antigens,
valid separate antibody/antigen chain assignments, and a nonempty resolved
interface. Candidates are ordered by resolution then instance. No two selections
share a PDB entry, inherited antigen cluster or exact resolved target sequence
fingerprint. This is not random sampling; cluster inheritance is not an
independent family-level validation.

## Contact definition

Any resolved target amino-acid heavy atom within **5.0 A inclusive** of an
amino-acid heavy atom in the selected antibody variable region defines a positive
target residue. Both heavy and light variable chains contribute when present;
framework contacts are included, not just CDR contacts. The cropped SAbDab2
complex represents a selected antibody instance, not necessarily all antibodies
in the original deposition. No assumption of absence of other epitopes is made.

Target protein heavy atoms alone are exported into `structures/`. All resolved
protein target chains assigned to the instance are retained, including residues
outside the interface. No antibody atoms, glycans, water, ions or other ligands
are exported. Bound coordinates and source numbering/insertion codes are
preserved. Exported targets are not experimentally measured apo structures.
Sources may be constructs or unresolved fragments, not full-length proteins.
Resolved/expected sequences from source metadata are retained in structures.tsv.
Unresolved residues, water bridges and glycan-mediated interfaces are unassessed.

## Files

- `structures/`: target-only searchable mmCIF files.
- `source_complexes/`: input complexes retained separately for reproducibility;
  **do not scan or index this directory** as target-only structures.
- `annotated_residues.tsv`: positive target residue labels only.
- `residue_labels.tsv`: all resolved target residues with contact status and
  minimum antibody distance. Zero means no contact to this selected antibody
  under this definition, not proven nonfunctional or unable to bind antibodies.
- `contacts.tsv`: all contacting target/antibody residue pairs, with minimum
  heavy-atom distance. Insertion codes and partner identities are retained.
- `panel.tsv`, `structures.tsv`, `provenance.json`, `build_manifest.json`,
  `validation.json`, `REPORT.md`.

Build:

```sh
/Users/mvh_mac/documents_local_only/IBC_search/.venv/bin/python /Users/mvh_mac/documents_local_only/IBC_search/scripts/build_antibody_targets.py
```

Validate again:

```sh
/Users/mvh_mac/documents_local_only/IBC_search/.venv/bin/python /Users/mvh_mac/documents_local_only/IBC_search/scripts/build_antibody_targets.py --validate-only
```

Validation re-extracts targets and recomputes contacts from retained complexes,
checks coordinate and identity equality, annotation tables, checksums, panel
membership and redundancy controls. No FoldDisco index or search is launched.