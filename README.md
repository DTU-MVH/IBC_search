# ibc-search

Surface ball-scanning extraction of **discontinuous structural motifs** from a
protein structure, emitted as FoldDisco queries.

This is the **first component** of a pipeline for identifying candidate
pre-existing anti-drug-antibody (ADA) epitopes on de novo protein binders. It
does three things and stops:

1. parses a structure (local file or PDB entry ID) into stable residue
   identities;
2. computes per-residue solvent accessibility;
3. extracts, deduplicates and tabulates candidate discontinuous surface
   motifs, and writes them in the FoldDisco query format.

## What the output does and does not mean

An extracted patch is a **geometric neighbourhood of solvent-exposed
residues**. It is a *candidate structural motif*.

It is **not** an epitope. It is **not** a validated antibody footprint. Nothing
in this component establishes antigenicity, antibody binding or
cross-reactivity, pre-existing ADA reactivity, induction of a new immune
response, or T-cell peptide presentation. These are six distinct,
non-interchangeable things and the code, column names and documentation keep
them apart.

**The `scan` command runs no FoldDisco search.** It stops at emitting
queries. Running them is a separate, explicit step —
`ibc-search folddisco-search`, which delegates to
[`scripts/folddisco_search.sh`](scripts/folddisco_search.sh) and leaves the
searching to the FoldDisco binary. Epitope labels, features, scoring and
benchmarking remain later stages; see
[`src/ibc_search/ext/`](src/ibc_search/ext/) for documented extension points
and [`docs/ASSUMPTIONS.md`](docs/ASSUMPTIONS.md) §5 for what blocks them.

## Install

Requires Python ≥ 3.11, NumPy, SciPy and biotite. No GPU, no large downloads.

```bash
cd ~/documents_local_only/IBC_search
python3 -m venv .venv
.venv/bin/python -m pip install numpy scipy biotite pytest
```

Then either install the package:

```bash
.venv/bin/python -m pip install -e .
```

…or skip the install and run from the source tree with `PYTHONPATH=src`, which
is how every command below is written.

> **Environment note (macOS + Python 3.14).** An editable install can appear
> to succeed yet leave the package unimportable. setuptools marks its
> `__editable__*.pth` file with the macOS `UF_HIDDEN` flag, and CPython 3.14
> skips hidden `.pth` files (`site.addpackage`), so the path is never added to
> `sys.path`. Symptom:
> `ModuleNotFoundError: No module named 'ibc_search'` straight after a
> successful `pip install -e .`. Confirm with
> `.venv/bin/python -v -c pass 2>&1 | grep -i pth` ("Skipping hidden .pth
> file"). Fix with
> `chflags nohidden .venv/lib/python3.14/site-packages/__editable__*.pth`, or
> just use `PYTHONPATH=src`. The test suite is unaffected: `pyproject.toml`
> sets pytest's `pythonpath = ["src"]`.

## Usage

All commands below have been run on this machine; none is hypothetical.

### Scan a local file

```bash
PYTHONPATH=src .venv/bin/python -m ibc_search scan \
    /Users/mvh_mac/Documents/master_thesis_2026/folddisco_test/pdb_files/4CHA.pdb \
    --chain A --chain B --chain C \
    -o results/4cha
```

### Scan a PDB entry ID (fetched from RCSB and cached)

Either spelling works — the legacy 4-character ID, or the extended
`pdb_XXXXXXXX` form the PDB is migrating to:

```bash
PYTHONPATH=src .venv/bin/python -m ibc_search scan 1CRN -o results/1crn
PYTHONPATH=src .venv/bin/python -m ibc_search scan pdb_00001crn -o results/1crn
```

Both resolve to the same download, cache entry and patch IDs (verified: the
two URLs serve byte-identical content). RCSB serves **no legacy PDB-format
file** for an entry outside the 4-character space, so with
`[input].fetch_format = "pdb"` such an entry is refused with a clear message;
the mmCIF default is unaffected.

mmCIF is preferred for fetched entries. The URL, UTC date and SHA-256 are
recorded in the run manifest. If the network is unavailable the run fails with
a clear error; it never silently falls back to another source.

### Useful flags

```bash
# one radius only, side-chain-centroid ball centres, a stricter seed threshold
PYTHONPATH=src .venv/bin/python -m ibc_search scan 1CRN -o results/1crn_sc \
    --radius 8 --center-mode sidechain_centroid --seed-rsa 0.30

# scan model 2 of a multi-model file
PYTHONPATH=src .venv/bin/python -m ibc_search scan model.pdb --model 2 -o results/m2

# override any configuration key (TOML scalar syntax)
PYTHONPATH=src .venv/bin/python -m ibc_search scan 1CRN -o results/x \
    --set connectivity.disconnected_policy="'flag'" \
    --set ball_search.radii_angstrom='[6.0, 9.0]'

# a user config file, merged over the shipped defaults
PYTHONPATH=src .venv/bin/python -m ibc_search scan 1CRN -c my_config.toml -o results/y
```

### `dumb_search`: one ball per residue, no thresholds

```bash
PYTHONPATH=src .venv/bin/python -m ibc_search dumb-search 1CRN --radius 6 -o results/1crn_dumb
```

The `scan` command shapes its motif set with five thresholds — a seed RSA
cutoff, a member RSA cutoff, a surface-connectivity check, deduplication of
identical member sets, and a minimum patch size. `dumb-search` removes all
five. It walks each scanned chain from its N-terminal residue to its
C-terminal residue, places one ball at every residue on the way, and keeps
whatever is inside:

> start at the terminal residue, draw the ball, move to the next residue, draw
> the ball, and so on until the final terminal.

Ball **geometry** is the same code: same `--radius`, same `--center-mode`,
same "a residue joins on its nearest heavy atom" membership rule, and the same
output files, so a `dumb-search` directory feeds `folddisco-search` and the
downstream scripts exactly as a `scan` directory does. One motif is emitted
per residue per radius — no more and no fewer — which makes the motif count a
property of the structure rather than of the parameters.

**The search input is deduplicated even though the motif set is not.** Two
balls centred on neighbouring residues of a helix can enclose the identical
residue set, and searching that query twice returns byte-identical hits for
twice the time. The query file therefore holds each distinct query once, the
query index keeps one row per patch, and `line_number` ties the sharing
patches to the line they share. Measured on a 5.5 Å background run, 0.16% of
queries were such repeats, in 32 of 60 structures sampled; a run that has any
says so in a warning.

It exists as a baseline. If a downstream result survives here, it does not
depend on the ball search's filtering; if it only appears under `scan`, the
filtering is doing the work. A `dumb_search` motif is a ball of residues, not
a surface patch: a ball centred on a buried residue is full of buried
residues, and nothing removes them. It is no more an epitope than anything
else this repository produces.

Motif IDs carry a `d` (`1CRN_m1_r6_d00042`) where `scan` writes a `p`, and the
run manifest records `ball_search.method = "dumb_search"`, so the two methods'
outputs can never be confused once separated from their directories.

### `few_balls_search`: annotated balls plus a matched background

```bash
PYTHONPATH=src .venv/bin/python -m ibc_search few-balls-search \
  investigation_databases/zincfingers/structures/12UM.cif.gz \
  --radius 5.5 -o results/12um_few \
  --seed-residues investigation_databases/zincfingers/annotated_residues.tsv \
  --seed-structure-id 12UM
```

Places one ball on every annotated residue, then the **same number** of
further balls elsewhere on the same protein. So eight annotated zinc-finger
residues give eight annotated balls and eight background balls, sixteen
motifs in total. A floor of five background balls applies however few
annotations a structure has.

The point is an **internal control**. Comparing a seeded run against a
separate background run confounds the annotation with fold, resolution,
deposition year and -- measured, and largest of all -- how many motifs each
residue happens to sit in. Here both classes come from one structure, one
radius, one index and one search, in equal number, so none of those can
explain a difference between them.

The background balls are chosen to be spread out *and* to represent the
protein rather than only its outside:

- **Spread out** by greedy farthest-point sampling: repeatedly take the
  residue whose nearest already-chosen centre is furthest away. The annotated
  centres are in the chosen set before the first pick, so the background is
  pushed away from the annotated sites as well as from itself.
- **Representative** by stratifying candidates into buried, exposed and
  RSA-unknown and giving each stratum a share proportional to its own share
  of the candidate pool. This is load-bearing, not a refinement: in a
  globular object the points furthest apart are the ones on its hull, so
  unstratified farthest-point sampling returns a *surface* sample. Annotated
  sites are often clefts and pockets, so that would confound annotation with
  burial.

Every motif is labelled in a `seed_class` column of `patches.tsv`, as
`annotated` or `background`. **Downstream analysis must read that column**
rather than infer the class from a patch ID or by re-deriving it from the
annotation table. Motif IDs carry an `f` (`12UM_m1_r5.5_f00007`) where `scan`
writes `p` and `dumb-search` writes `d`.

Nothing is filtered or deduplicated: a ball keeps every residue inside it,
annotated or not, buried or not. An annotated ball and a background ball may
overlap, and the run manifest reports how many background seeds landed inside
an annotated ball, together with the dispersion actually achieved and the
per-stratum request and delivery. Where the protein cannot supply the
requested background count, or where the floor makes the two classes unequal,
the manifest says so and the run warns rather than quietly returning a
smaller set.

A `few_balls_search` motif is a ball of residues. It is not a surface patch,
not an epitope, and establishes nothing about antigenicity.

### Is an annotated residue better than its own structure's background?

```bash
.venv/bin/python scripts/annotated_vs_background.py \
  results/zincfingers_fewballs_swissprot_r5.5
```

Reads a `few_balls_search` run and scores every annotated residue **1** or
**0**: 1 when its own ball matched more distinct sites than the median of the
background balls on the *same structure at the same radius*, 0 when it did
not, and blank where that structure had no background ball to compare
against. Blank is not 0 — an unanswerable comparison is not a failed one.

The test is strictly greater, since with small integer counts a tie is not
evidence of being better, and `percentile_within_background` is reported
beside the bit so a near miss can be told from a clear loss.

Writes to the run's standard layout:

| output | what it is |
|---|---|
| `01_statistics/annotated_vs_background.tsv` | one row per annotated residue, with its count, its structure's background median, the percentile and the 1/0 |
| `01_statistics/annotated_vs_background.json` | mean, SD, median, IQR and distinct counts per category, per structure and pooled |
| `00_figures/annotated_bins.png` | the binned dot map, one bin per structure, **annotated residues in red** against their own background in blue |
| `00_figures/annotated_distribution.png` | mean, ±1 SD and distinct counts per category — by seed class, and annotated residues split by verdict |

Two rates are reported and neither is given alone: the **pooled** fraction
treats each annotated residue as one observation and is dominated by heavily
annotated structures, while the **mean of per-structure fractions** treats
each structure as one.

A caution that applies to the means but not the verdict: the count is
censored by `--top`. Where a large share of annotated balls return exactly
the cap, the annotated mean is a lower bound. The verdict survives that — a
capped ball still beats a background median well below the cap — and only
becomes unreliable if both classes saturate.

### Filtering hits by IDF

`folddisco_search.sh --min-idf X` (and `run_panel.py --min-idf X`) drops hit
rows whose `idf` is below X, keeping them in a `.lowidf.tsv` side file so the
floor is auditable and a different one can be applied without re-searching.

This is **not** FoldDisco's own `--score`. Measured on the installed build
against `afdb_swissprot`, where per-hit `idf` ran 3.6 to 888 with a median of
43.7, `--score 0.5` kept every row and `--score 1` and above kept none — it
works on a different, normalised scale, so `--score 50` silently empties the
output instead of filtering at 50.

IDF is index-dependent, so a floor must be chosen per index rather than
inherited. Median per-hit IDF is roughly 478–787 against whole-PDB and 43.7
against `afdb_swissprot`, so a floor of 50 removes about 1% of PDB hits and
about 56% of Swiss-Prot hits.

### Inspect configuration and the FoldDisco contract

```bash
PYTHONPATH=src .venv/bin/python -m ibc_search show-config
PYTHONPATH=src .venv/bin/python -m ibc_search list-databases
PYTHONPATH=src .venv/bin/python -m ibc_search folddisco-info
```

`folddisco-info` reports the installed FoldDisco build and whether it still
matches the build the query contract was verified against. It runs no search.

### Run FoldDisco on the results

Searching is a separate step, performed by the FoldDisco binary itself. The
pipeline is `scan` → `folddisco-search`:

```bash
PYTHONPATH=src .venv/bin/python -m ibc_search scan 8YHU -o results/8yhu
PYTHONPATH=src .venv/bin/python -m ibc_search folddisco-search results/8yhu
```

That writes `results/8yhu/folddisco_hits.tsv`, which is FoldDisco's per-match
output with each row joined back to the `patch_id` that produced it.

The same thing from a plain shell, which is where the logic actually lives:

```bash
./scripts/folddisco_search.sh status                       # what is installed
./scripts/folddisco_search.sh extract                      # unpack the index (slow)
./scripts/folddisco_search.sh search results/8yhu
./scripts/folddisco_search.sh search results/8yhu -- --top 10 --rmsd 1.5
```

Anything after `--` goes straight to `folddisco query`.

> **The PDB database needs unpacking first.** The database directory holds
> the whole-PDB FoldDisco index as a single 65.7 GiB `pdb_folddisco.tar.lz4`.
> Until `extract` is run there is nothing to search, and the script says so
> rather than failing obscurely. The AFDB indexes below ship unpacked.
>
> Unpacked it is **66.31 GiB on disk** (65.87 GiB logical) in 230 659 files:
> an 18.7 GiB index table plus ~47 GiB of `index/pdb/**/pdb*.ent.gz`
> structure files. lz4 buys almost nothing here, so the unpack is barely
> larger than the archive, and keeping both needs ~132 GiB. Decompressing the
> full archive measured 19 s, so extraction costs minutes, not hours — it is
> bound by writing 66 GiB and creating 230k files.

### Choosing a reference database

Four are registered. List them, with whether each is ready to use:

```bash
PYTHONPATH=src .venv/bin/python -m ibc_search list-databases
```

Select one by name wherever a search runs:

```bash
PYTHONPATH=src .venv/bin/python -m ibc_search folddisco-search results/8yhu -d afdb_swissprot
.venv/bin/python scripts/run_panel.py results/myrun --radius 6 -d afdb_proteome
.venv/bin/python scripts/idf_zscore.py apply hits.tsv -d pdb
```

| name | contents | size | IDF parameters |
|---|---|---|---|
| `pdb` *(default)* | the whole PDB, 230 655 deposited structures | 132 GB | fitted and calibration-checked |
| `afdb_proteome` | AlphaFold DB proteome subset, 564 446 **predicted** models | 26 GB | **none fitted** |
| `afdb_swissprot` | AlphaFold DB Swiss-Prot subset, predicted models for reviewed UniProtKB entries | 24 GB | **none fitted** |
| `sabdab2_antigen` | SAbDab2 antigen chains, antibody stripped, epitope-labelled | small | fitted |

A name resolves to an index **and** the IDF parameters fitted against it, so
the two cannot be mixed by hand. The registry lives in
`[folddisco.databases]` in `src/ibc_search/data/default.toml`; add an entry
there to register another index.

**The two AFDB databases have no fitted IDF parameters**, so
`idf_zscore.py apply -d afdb_proteome` refuses rather than borrowing the
PDB's. That refusal is correct: IDF measures rarity over an index's own
contents, and measured on one 4-residue motif the same query scores IDF
**80.5** against the PDB index, **41.5** against `afdb_swissprot` and
**16.2** against `afdb_proteome`. Fit first with `idf_zscore.py fit`, then
add `idf_params` to the entry. Until then these indexes are usable for
counting hits — see `count_scoring.py`, which needs no fit — but not for
z-scores.

Two further cautions specific to AFDB. Its models are **predicted, not
observed**, each carrying a pLDDT confidence this pipeline does not read, so
a hit says a geometry occurs in a prediction rather than in a measured
structure. And the AFDB indexes carry no SIFTS chain mapping, so
`hit_organism_stats.py` does not apply to them; pass `--lookup` to point it
at whichever index was actually searched.

### Where the databases live

```
~/documents_local_only/master_thesis_data_files/folddisco_databases/
    index/pdb_folddisco           the PDB index (+ .lookup/.offset/.type)
    index/pdb/**/pdb*.ent.gz      its structure files
    pdb_folddisco.tar.lz4         the 65.7 GiB source archive
    afdb_proteome/afdb-proteome_folddisco
    afdb_swissprot/afdb_swissprot_v4_folddisco
```

233 GB in total. Every path inside an index — the `foldcomp_db` key in its
`.type` file, and the entries in its `.lookup` — is **relative to the index
directory**, so the whole tree can be moved without rebuilding anything;
only the paths registered in `default.toml` need updating.

Override the location with `IBC_FOLDDISCO_DB_DIR`, or an individual index
with `IBC_FOLDDISCO_INDEX` / `folddisco-search -i`.

> **Results are capped at the top 1000 structures per query, deliberately.**
> Measured against the extracted whole-PDB index (230 655 structures), one
> 3-residue query took 0.18 s at `--top 10`, 0.41 s at `--top 100`, 4.2 s at
> `--top 1000`, and produced **no output in 15 minutes uncapped**. The hash
> lookup is fast; the match stage loads and superposes every candidate, and a
> short motif matches much of the PDB. Without the cap the search does not
> finish.
>
> Treat `[folddisco].search_top` as a **retrieval cutoff, not a speed knob**:
> a true match ranked below it is never reported, and on a real scan the cap
> saturates — every query returned exactly 1000 rows, so all of them were
> truncated. How retrieval quality varies with this value is uncharacterised
> and belongs in the sensitivity study. `--no-top` removes it if you can
> afford the runtime.

> **Threads default to 1, deliberately.** The installed build does not write
> stdout line-atomically: at `-t 10`, 0–4 of ~1225 hit rows came back torn by
> interleaved writes, differing between otherwise identical runs, while `-t 1`
> was byte-identical across repeats. Torn rows are quarantined in
> `folddisco_hits.malformed.tsv` and dropped rather than parsed, and
> `--strict` makes them fatal. See
> [`docs/FOLDDISCO_CONTRACT.md`](docs/FOLDDISCO_CONTRACT.md).

A hit means a geometrically similar motif exists in a deposited structure. It
is not evidence of antigenicity: the PDB is a structural background, not a
collection of known epitopes, and `[folddisco].epitope_reference_index_path`
is still empty and still marked `blocked`.

### Tests

```bash
# 129 tests, fully offline and self-contained (synthetic fixtures)
.venv/bin/python -m pytest tests/ -m "not network and not folddisco"

# + 7 contract-regression tests that invoke the real FoldDisco executable
.venv/bin/python -m pytest tests/ -m folddisco

# everything, including 1 live RCSB fetch: 137 tests
.venv/bin/python -m pytest tests/
```

The `folddisco`-marked tests re-check the query contract against the installed
binary, including feeding it a freshly emitted query file. They skip
automatically if `folddisco` is not on `PATH`.

## Outputs

| file | contents |
|---|---|
| `residues.tsv` | one row per residue: identity, SASA, RSA (+ missing reasons), exposure flag, residue class, quality/missing-atom flags, ball-centre kind, FoldDisco expressibility |
| `patches.tsv` | one row per deduplicated patch: ID, seeds, radius, scale, centre mode, segment/gap structure, diameter, radius of gyration, compactness, cluster, connectivity counters, query string |
| `patch_members.tsv` | one row per (patch, member): the **selecting atom** and its distance to the ball centre, per-member RSA and flags |
| `folddisco_queries.tsv` | the handoff: `<structure_path>\t<tokens>` lines, the format verified against the installed FoldDisco build |
| `folddisco_queries.index.tsv` | maps each query line back to its `patch_id` (the query format has no field for one) |
| `run_manifest.json` | versions, exact command, input provenance + hashes, fully resolved config, seed, warnings, per-stage timings |

`folddisco-search` adds two more to the same directory:

| file | contents |
|---|---|
| `folddisco_hits.tsv` | FoldDisco's per-match rows, each prefixed with the `patch_id`, radius and scale it came from |
| `folddisco_hits.malformed.tsv` | written **only** if the build tore output rows; quarantined, never parsed |

`[run].table_format = "parquet"` writes Parquet instead of TSV (needs
`pyarrow`).

The query file was validated by feeding it to the real FoldDisco binary: all
545 queries emitted for 4CHA were consumed with exit 0, empty stderr and zero
panics, returning 4 871 hit rows in 2.98 s.

> **Reading hits back.** `matching_residues` is **not** reliably positionally
> aligned to `query_residues` in the installed build — only ~62% of real hit
> rows have equal-length fields, because unmatched positions are sometimes
> padded with `_` and sometimes omitted outright. Use
> `folddisco_contract.parse_hit_residues()` (set semantics, never raises) and
> see [`docs/FOLDDISCO_CONTRACT.md`](docs/FOLDDISCO_CONTRACT.md) for the
> measured counts. The `query_residues` side *is* reliable.

**Missing values** are written as an empty cell with the reason in a
companion `*_missing_reason` column — never as `0`. A SASA of 0.0 is a real
measurement for a buried residue, so overloading zero would corrupt the data.

## The standard run layout

Every FoldDisco run ends up with the same shape. `run_panel.py` finalises by
default; `--no-finalize` opts out, and `finalize_run.py` can be re-run on any
existing run directory at any time.

```
<run>/00_figures/              run-level figures
<run>/01_statistics/           run-level statistics
<run>/<instance>/00_figures/   per-structure figures, PyMOL session, PDB
```

```bash
.venv/bin/python scripts/run_panel.py results/myrun --radius 6 -d pdb   # finalises
.venv/bin/python scripts/finalize_run.py results/myrun                  # or re-finalise
```

This exists because the pieces — `residue_scores.py`, `count_scoring.py`,
`epitope_statistics.py`, the plotters, `render_all_structures.py` — used to
be run by hand in whatever order the moment called for, so two runs came out
with different figures and different thresholds and comparing them meant
first working out what each had been given. `finalize_run.py` is the single
answer to "what does a run always contain".

### The minimum, every run

**`00_figures/`**

| figure | what it shows |
|---|---|
| `score_bins.png` / `count_bins.png` | the **binned dot map**: one bin per structure, every assessed residue a point, the **median a black line** with a vertical black IQR bar and horizontal caps |
| `count_heatmap.png` | structures down, residue position across, **green → teal → blue** by distinct-hit count |
| `rate_bins.png` / `rate_heatmap.png` | the same two views of `hits_per_motif`, the count normalised by how many motifs the residue was searched in |

Both are written in light and dark (`_dark.png`). On the dark surface the
median ink is white — black would be invisible — and the caption names
whichever it drew rather than asserting a colour that isn't there.

The count ramp's lightness **decreases** monotonically (0.884 → 0.320 in
HLS, checked by a test), so the order survives greyscale printing, "many
hits" is the visually heavy end, and the ramp stays separable under the
common forms of colour blindness. Grey cells are positions a shorter chain
does not have — **not** a count of zero, which means assessed and unmatched.

**`01_statistics/`**

`distributions.json` and `<metric>_per_structure.tsv` always: per-structure
and pooled min/median/mean/sd/max, counts above threshold, and how many
residues scored zero. These are **descriptive**.

`set_statistics.json` / `.tsv` always: the whole-set view — the pooled
distribution over residues, the distribution of per-structure medians (two
different questions), and the **mean RSA of the residues a figure colours
red** against the rest. A count run also gets
`set_statistics_hits_per_motif.json` / `.tsv`, the same battery on the
effort-normalised score.

The two are kept **side by side, never merged**. The 3-D figures rank on
`n_distinct_hits`, so the normalised file's red set is a different set of
residues from the one on screen — it carries an explicit warning saying so.

When the run has epitope labels — detected from the data, never from a flag —
it also gets the full battery that `results/sabdab40_full_r6` established as
the baseline: `summary.json` with the pooled AUC, its bootstrap interval
resampled **over structures**, the permutation null and excess over it, the
coverage artifact, `threshold_sweep.tsv` and `per_structure_auc.tsv`, plus
their figures in `00_figures/`.

With no labels that battery is **skipped, not computed on an empty positive
class**, and `run_summary.json` says so in `no_ground_truth_note`. An AUC
against no positives would look like performance while measuring nothing.

**`<instance>/00_figures/`**

Per structure: the ray-traced still with its **see-through surface cloud**
over the cartoon, CA spheres on called or epitope residues where there are
any, the `.pse` PyMOL session, the `.pml` script, and a `.scored.pdb`
carrying the score in its **B-factor column** so any viewer — PyMOL,
ChimeraX, Mol\*, VMD — reproduces the colouring without this code.

3-D figures keep the **red/grey** scheme for both scorings. One convention in
3-D means a reader never has to work out which scale a structure was coloured
on; the green→blue ramp is the 2-D heatmap's alone.

## Downstream analysis scripts

`scan` and `folddisco-search` are the component this repository is *about*;
the scripts below consume their output. They are deliberately separate files
rather than subcommands, because each one makes an interpretive choice that
the core component refuses to make.

| script | what it does |
|---|---|
| `run_panel.py` | scans and searches every structure in a `panel.tsv`, N at a time, one single-threaded FoldDisco per structure; resumes by re-invocation |
| `idf_zscore.py` | `fit` / `apply` / `minmax` — size-detrended normalisation of FoldDisco IDF (see its docstring for why raw IDF cannot be used directly) |
| `residue_scores.py` | projects motif z-scores onto residues, each residue taking the **maximum** z of the motifs containing it; scores against epitope labels where they exist |
| `finalize_run.py` | produces the standard run layout above: scoring, `00_figures/`, `01_statistics/` and the per-structure renders. Run by `run_panel.py` automatically |
| `count_scoring.py` | the **count_scoring** alternative: each residue takes the number of **distinct matched sites** its motifs hit, deduplicated per residue; a count, never z-scored |
| `render_pymol.py` | one structure as a score-coloured cartoon in a translucent surface, via PyMOL; writes the score into the B-factor column so any viewer can reproduce it |
| `render_all_structures.py` | the above for a whole run, in parallel, into `<run>/<instance>/00_figures/` |
| `plot_epitope_scores.py` | the run-level summary figures, chosen by whether the run has epitope labels |
| `hit_organism_stats.py` | human vs non-human hits **against the composition of the searched index**, which is the only way that fraction means anything |
| `top_motifs_by_organism.py` | per-hit organism annotation for the top-scoring motifs |
| `minibinder_epitope_sample.py` | selects small SAbDab2 antigens and labels contact epitopes |
| `ada_binder_panel.py` | assembles the designed-binder panel with its ELISA annotation |

Two conventions hold across all of them, and both exist because the opposite
fails silently:

- **a residue in no emitted motif has no score** — it is absent from the
  output, never scored zero, because zero sits in the middle of the real
  range and would read as "assessed, unremarkable";
- **an unknown organism is `unknown`, never "not human"**, so any fraction
  computed over those rows has an honest denominator.

`count_scoring.py` sharpens the first convention rather than breaking it. A
residue whose motifs **were searched and matched nothing** scores `0`: it was
assessed and the database had nothing like it. A residue in no emitted query
is still absent. For a count, unlike a z-score, zero is at the end of the
range rather than the middle, so it can carry that meaning honestly.

### `count_scoring` vs `residue_scores`

The two answer different questions and will often disagree, by design:

| | `residue_scores.py` | `count_scoring.py` |
|---|---|---|
| per-residue rule | **maximum** IDF z of its motifs | **count** of distinct sites where the query residue explicitly matched |
| normalised? | yes — size-detrended, z-scored against a background | **no, and must not be** — a count is not a position in a null |
| high value means | the best motif here is **rare** | motifs here are **common** in the database |
| never assessed | absent | absent |
| assessed, found nothing | n/a (no hits means no score) | `0` |

A hit is one matched site: a (reference structure, matched residue set) pair.
Matched tokens are compared as a set and `_` placeholders are dropped, so
token order cannot split one site in two. Deduplication is **per residue**:
the same site reached through two of a residue's overlapping motifs counts
once for it, because it is one piece of evidence, while the same site counts
separately for each residue it informs. By default, a hit is attributed **only
to explicitly matched query residues**. A `_` target position gets no credit.
The scorer rejects unequal-length query/target fields rather than guessing
which query residues matched or silently reporting incomplete counts as zero.
The entire panel is validated before the output table is opened.

**Existing searches may require new, complete mappings:** the installed
FoldDisco sometimes omits unmatched positions, making those rows ambiguous.
Superposition columns alone are not decoded by this scorer. For reproducibility
only, `--attribution motif` explicitly restores the old whole-motif attribution;
the summary JSON records the chosen mode. Do not mix modes in a comparison.
`finalize_run.py` defaults to matched attribution too, accepts
`--count-attribution motif` for legacy runs, and stops if count scoring fails
instead of drawing figures from stale count tables.

```bash
.venv/bin/python scripts/count_scoring.py results/bg100_r6 --threshold 1
```

Writes `residue_hit_counts.tsv` and `residue_hit_count_summary.json`.

Three columns, two of them scores:

| column | meaning |
|---|---|
| `n_distinct_hits` | distinct matched sites, deduplicated per residue — the raw count |
| `n_motifs` | how many emitted queries contained the residue — the coverage |
| `hits_per_motif` | `n_distinct_hits / n_motifs` — the count **normalised for search effort** |

`dumb_search` places one ball per residue, so a mid-chain residue sits in
more balls than one near a terminus and is searched more often. The raw
count mixes that in; `hits_per_motif` divides it out. On the 187-structure
background run the difference is large — the mean RSA of the top 20% moves
from **0.071 raw to 0.230 normalised**, against 0.267 for everything else,
so most of the "high-count residues are buried" effect was coverage rather
than matching.

Two cautions. The ratio is **not** the mean hit count of the residue's
motifs: the numerator is deduplicated across them, so overlapping motifs
matching the same site count once and the ratio falls below that mean. And
it is bounded above by the search `--top` cap exactly as the count is — on a
censored run both pile up at the cap (pooled max was exactly 40.0 at
`--top 40`).

## The algorithm, in brief

For each solvent-exposed residue (a *seed*), place a ball at a representative
position for that residue, collect the exposed residues whose nearest heavy
atom lies inside the ball, discard members not surface-connected to the seed,
and emit the residue set as a candidate discontinuous motif.

| step | what happens | key parameter |
|---|---|---|
| A exposure | Shrake–Rupley SASA; RSA vs Tien et al. 2013 theoretical maxima | `probe_radius_angstrom` = 1.4 Å |
| B seeds | residues with RSA ≥ threshold | `seed_rsa_threshold` = 0.20 |
| C centres | Cβ atom, or side-chain heavy-atom centroid; glycine → **flagged virtual Cβ** | `center_mode` |
| D neighbours | min heavy-atom distance to centre ≤ radius; selecting atom + distance recorded | `radii_angstrom` |
| E discontinuity | segment count, lengths and author-numbering gaps annotated | — |
| F connectivity | keep the seed's contact-graph component | `contact_cutoff_angstrom` = 5.0 Å |
| G dedup | identical member sets merged, **all** seeds retained; Jaccard clustering | `jaccard_threshold` = 0.8 |
| H scales | several radii per run, relationship recorded via `member_set_hash` | — |

Sequence discontinuity is central: this is a 3D neighbourhood search, **not** a
sequence sliding window. On a hairpin fixture and on real structures the large
majority of patches are sequence-discontinuous (510 of 545 on 4CHA chains
A/B/C).

### Four parameters that are never conflated

Several default to similar magnitudes and are deliberately kept separate:

1. the SASA solvent probe radius — `[exposure].probe_radius_angstrom`
2. the patch-internal residue adjacency cutoff — `[connectivity].contact_cutoff_angstrom`
3. the patch extraction ball radius — `[ball_search].radii_angstrom`
4. the FoldDisco search tolerances — `[folddisco].distance_threshold_angstrom`, `.angle_threshold_degrees`

A fifth quantity, the ~5 Å antibody–antigen contact cutoff from the project
notes, is **not computed anywhere** in this component and must not be aliased
onto (2).

### Known limitations

A Euclidean ball is not a geodesic surface patch and not a physically
accessible antibody footprint. The connectivity check reduces, but does not
eliminate, merging across thin regions. Patches heavily overlap and are not
independent. The structure is treated as a single rigid conformer. See
[`docs/ASSUMPTIONS.md`](docs/ASSUMPTIONS.md) §4.

## Configuration

[`src/ibc_search/data/default.toml`](src/ibc_search/data/default.toml) is the
single source of defaults. It lives inside the package so that it ships with a
normal install; dump the resolved values any time with `ibc-search show-config`.
Every parameter carries a **unit**, a **provenance** note and a **status**
marker — `documented`, `verified`, `proposed` or `blocked`. A `proposed` value
is an engineering choice, not an established criterion, and must not be cited
as one. The fully resolved configuration is copied into every run manifest, so
a run is reproducible from its own output.

## Performance

Measured on this machine (Apple Silicon, arm64, Python 3.14.2), median of 3
runs, all four default radii (`[4, 5, 6, 7] Å`), single process:

| structure | residues | heavy atoms | seeds | patches | parse | SASA | ball search | **total** |
|---|---|---|---|---|---|---|---|---|
| 1CRN (crambin, minibinder scale) | 46 | 327 | 35 | 139 | 5 ms | 6 ms | 35 ms | **0.05 s** |
| 4CHA chains A/B/C | 562 | 3 591 | 110 | 435 | 16 ms | 35 ms | 117 ms | **0.17 s** |
| 6M0J chain E (SARS-CoV-2 RBD) | 878 | 6 558 | 109 | 433 | 68 ms | 31 ms | 118 ms | **0.22 s** |
| 3HFM (lysozyme–Fab, all chains) | 559 | 4 297 | 283 | 1 119 | 46 ms | 87 ms | 376 ms | **0.51 s** |

Well inside the "tens of seconds per binder" budget, and a few hundred
structures is a matter of minutes. Spatial queries use SciPy KD-trees, hot
paths are vectorised, parallelism is optional and bounded (default 1 worker),
and output is byte-identical regardless of worker count.

## Investigation databases

`investigation_databases/` holds curated structure sets, one directory each,
built by a script and described by a `build_manifest.json` that records every
filter that shaped it. Each directory has the same shape:

```
<set>/structures/*.cif.gz     the structure files
<set>/structures.tsv          one row per structure
<set>/build_manifest.json     what was asked for, what was found, what was rejected
```

| set | what it is | built by |
|---|---|---|
| `background_dataset` | 750 structures selected **for nothing** — a deliberately mixed background | `build_background_dataset.py` |
| `enzymes` | entities carrying an EC assignment, with their `struct_site` labels | `build_investigation_db.py` |
| `zincfingers` | zinc-finger entities, labelled by their ZN site | `build_investigation_db.py` |
| `protein_protein_interactions` | multi-entity entries, labelled by computed interface | `build_investigation_db.py` |

### `background_dataset`

The other three sets are selected *for* something. This one is the
comparison: a spread wide enough that a statistic computed on it describes
proteins in general rather than one family, one organism, one size or one
method. Composition is the product, so every axis of it is in the manifest.

```bash
.venv/bin/python scripts/build_background_dataset.py          # defaults below
.venv/bin/python scripts/build_background_dataset.py --total 300 --alphafold 100
```

**750 structures, at least 250 of them AlphaFold predictions**, 50–2300
residues, non-redundant at 30% sequence identity.

How the spread is produced, since random sampling alone would not produce it:

- **Non-redundancy.** One representative per 30% sequence-identity cluster.
  Raw PDB hits are dominated by whatever the field deposits hundreds of times,
  and a background made of lysozyme and trypsin is not a background.
- **Random, not best-first.** `build_investigation_db.py` takes representatives
  in ascending resolution, which is right for a category and wrong here:
  sorting by resolution selects small, well-ordered, heavily-studied X-ray
  structures, biasing three of the axes this set is meant to span. Candidates
  are shuffled under a recorded seed, so the run is reproducible without being
  ordered.
- **Size strata.** Seven bins from 50 to 2300 residues with equal quotas.
  Unstratified sampling would reproduce the PDB's own distribution, which is
  concentrated at 100–350 residues, and the large end would be nearly absent.
  Bin widths grow because the PDB thins out with size. A bin that cannot be
  filled is topped up from elsewhere and **named in the manifest**, so a gap
  is visible rather than smoothed over.
- **Organism cap.** At most 22 structures (~3%) from any one organism, so
  *E. coli*, human and yeast cannot crowd out the rest.

#### "High quality" means three different things here

This is the part that cannot be one column, and the table does not make it
one:

| source | quality filter | `resolution_angstrom` | `mean_plddt` |
|---|---|---|---|
| X-ray, cryo-EM | resolution **< 3.0 Å** | the value | empty |
| NMR | *no resolution exists* — admitted on method | **empty** | empty |
| AlphaFold | mean **pLDDT ≥ 70** | **empty** | the value |

An NMR entry's blank resolution is **not zero and not the threshold**: the
quantity does not exist for that method, and a search that required one would
silently exclude every NMR structure while looking like a quality filter.
An AlphaFold model has no resolution at all; pLDDT is the model's confidence
in *itself*, it is not comparable to a resolution, and it lives in its own
column. Pooling the three into one "quality" number mixes a measurement, the
absence of a measurement, and a prediction's self-assessment.

#### What a row does not mean

The set is **unlabelled by design** — no epitopes, no interfaces, no sites.
That is what makes it a background. An AlphaFold row is a **prediction**: its
coordinates were never observed, and a high pLDDT is confidence, not
evidence. Experimental and predicted rows should not be pooled without saying
so; `source` distinguishes them in every row.

Residue counts in `structures.tsv` are **parsed from model 1 of the written
file**, not taken from RCSB's deposited monomer count — the latter shortlists
candidates cheaply but counts a different thing, so the size bin and the
50–2300 bound are decided on the file that actually landed on disk.

## Repository layout

```
docs/ASSUMPTIONS.md          assumptions register + scientific limitations
docs/FOLDDISCO_CONTRACT.md   the verified FoldDisco interface and the evidence
src/ibc_search/
    data/default.toml        defaults with units, provenance and status markers
    cli.py                   scan / dumb-search / few-balls-search /
                             list-databases / show-config /
                             folddisco-info / folddisco-search
    config.py                loading, overrides, validation
    structure_io.py          input resolution, fetching, parsing, quality flags
    residue_id.py            stable identity + original/internal/query mapping
    exposure.py              SASA and RSA
    ball_search_algo.py      the extraction algorithm
    ball_placement.py        shared ball placement: member universe, centres,
                             patches -- used by every method below
    dumb_search.py           the dumb_search baseline: one ball per residue
    few_balls_search.py      annotated balls plus an equal, dispersed
                             background set on the same structure
    folddisco_contract.py    query emission and hit round-trip (no search)
    outputs.py               tables and the query file
    manifest.py              run manifest
    ext/                     documented extension points for later stages
investigation_databases/     curated structure sets, one directory each
    background_dataset/      the unlabelled mixed background (750)
    enzymes/ zincfingers/ protein_protein_interactions/
tests/                       122 offline tests + 1 live-network test
```

## Related

`../folddisco_test/` holds the earlier FoldDisco experiment (a 4CHA index and
`folddisco_test_motif_search.sh`) that was used as a live index to verify the
query contract. It is left untouched.
