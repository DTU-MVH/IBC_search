# Retrieved target-site validation

This is an **exploratory localization test**, not a functional annotation
predictor. The implementation is
`/Users/mvh_mac/documents_local_only/IBC_search/scripts/validate_target_sites.py`.

## Reproduce

```sh
/Users/mvh_mac/documents_local_only/IBC_search/.venv/bin/python \
  /Users/mvh_mac/documents_local_only/IBC_search/scripts/validate_target_sites.py \
  all --max-targets 2000 --workers 4
```

The default output directory is
`/Users/mvh_mac/documents_local_only/IBC_search/results/target_site_validation`.
Collection will not overwrite an existing hit database. After collection,
`annotate` and `report` can be run separately. Set `--max-targets 0` on the
annotation stage to extend the run to all qualified targets, then rerun
`report`. Annotation files are cached; the manifest records the selected
target pool. Sampling ranks target IDs by SHA256(seed:target), independently
of biological evidence and query category.

The completed initial run also has a separate `curated_target_census`
subdirectory. It uses the same hit database, restricted to **every** retrieved
target with an explicit M-CSA reference mapping (52 targets in this run), not
the random sample. For reproduction, create a separate output directory with
a `hits.sqlite` symlink to the collected database and a copy of the frozen
`mcsa_reference_residues.json`, then run `annotate --curated-only` and `report`
with that directory as `--out`. Do not merge that label-selected census into
the random target pool when interpreting the coordinate-derived tracks.

## Query/hit selection

- Existing completed panel entries only; no new searches.
- X-ray background; original annotated seeds, not expanded seeds.
- Enzymes: only query seeds with ligand class `biological`. These are
  ligand-contact labels, **not curated catalytic query annotations**.
- Identical query residue sets within a query structure are collapsed.
- Default motif size 6–14, RMSD <=2 Å, >=80% coverage, equal-length residue
  fields, unique matched target tokens; same-PDB self hits excluded.
- Existing top-40 candidate retrieval remains a limitation. Homologues and
  related families are **not** excluded.

## Target annotations and unknowns

The local indexed PDB coordinates are used directly, avoiding numbering
differences between the index and newly downloaded coordinate versions.
Insertion-code collisions and missing target residues cannot be resolved and
are unknown. Only model 1 and the parser's first alternate conformer are used.

1. **Zinc:** amino-acid heavy atom within 2.8 Å of observed Zn in the indexed
   asymmetric unit. Targets without observed Zn are unknown, not negatives.
   This is coordination evidence, **not zinc-finger domain classification**.
2. **Ligand:** amino-acid heavy atom within 4 Å of a non-water,
   non-additive/non-metal hetero ligand, using the existing project's ligand
   classification. This heuristic can include inhibitors or unrecognized
   additives; it is not a curated substrate/active-site assignment. Targets
   without eligible ligand atoms are unknown.
3. **Catalytic:** explicit manually curated M-CSA reference PDB residues,
   using API author residue numbers and chain names, validated against target
   residue names. No homologue-transferred labels. Only explicitly mapped
   reference chains are assessable. Other residues in those chains mean
   "not in this curated catalytic set", not proof of no biological function.
4. **Interface:** <=5 Å protein heavy-atom contact in the **first declared
   biological assembly**. Symmetry copies are separate partners even when
   author chain IDs are equal. Labels are projected to original author tokens;
   a token is positive if any assembly copy contacts another protein copy.
   Absence from that assembly, missing assembly metadata or annotation errors
   are unknown. An observed monomeric assembly has zero protein-interface
   labels; this does not establish inability to interact in another state.

Large structures/assemblies are skipped under a recorded atom safety limit.
No interface is inferred from asymmetric-unit crystal contacts as a fallback.
Annotation failures are cached and counted, never silently converted to zero.

Primary annotation references:

- https://www.ebi.ac.uk/thornton-srv/m-csa/download/
- https://www.ebi.ac.uk/thornton-srv/m-csa/api/residues/?format=json
- https://www.biotite-python.org/latest/apidoc/biotite.structure.io.pdb.get_assembly.html

The downloaded M-CSA data and SHA256 hash are retained with the results.

## Controls and statistics

For each assessable target site, observed overlap is the fraction of matched
target residues in the relevant biological-site set. The expected overlap is
the exact expectation obtained by replacing every matched residue with a
residue of the **same amino-acid name in the same target chain**. It controls
target identity, chain and residue composition, but does **not** control
exposure, secondary structure, local geometry or family similarity.

Site-level overlaps and expectations are averaged within each query motif.
Annotated versus background motifs are compared within common motif-size
strata (>=10 assessable query motifs per group), with weights from background
motif counts. The report includes observed, expected, and excess-over-control
comparisons. A query-structure bootstrap produces exploratory intervals.
These intervals condition on the selected targets: shared target structures,
families and cross-group dependencies are **not** handled by that bootstrap.
Do not interpret them as a family-held-out significance test.

The mapped annotated seed is additionally reported as an anchor-specific
check. It is unknown where that query position was unmatched. Distinct target
sites are also summarized to expose repeated retrieval of the same sites.

## Outputs

- `collection_manifest.json`: filters, input hashes and collection accounting.
- `hits.sqlite`: qualifying deduplicated matches and emitted query motifs.
- `annotation_manifest.json`: target-universe/sample sizes, IDs and provenance.
- `mcsa_reference_residues.json`: frozen curated API response.
- `target_annotations/`: per-target labels, provenance and failures.
- `target_site_overlap.tsv`: auditable per-query-motif/target-site observations.
- `target_site_validation.json`: coverage, unknowns, controls and comparisons.
- `REPORT.md`: concise results and limitations.

No numeric result by itself proves shared reaction mechanism, interface partner,
binding affinity or biological specificity. The correct initial claim is
**preferential localization to the relevant class of biological sites**, if
the observed and excess-over-control comparisons actually support it.