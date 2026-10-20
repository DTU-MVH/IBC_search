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
cd /Users/mvh_mac/Documents/master_thesis_2026/IBC_search
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

```bash
PYTHONPATH=src .venv/bin/python -m ibc_search scan 1CRN -o results/1crn
```

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

### Inspect configuration and the FoldDisco contract

```bash
PYTHONPATH=src .venv/bin/python -m ibc_search show-config
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

> **The database needs unpacking first.** `../folddisco_databases/` holds the
> whole-PDB FoldDisco index as a single 65.7 GiB `pdb_folddisco.tar.lz4`.
> Until `extract` is run there is nothing to search, and the script says so
> rather than failing obscurely.
>
> Unpacked it is **66.31 GiB on disk** (65.87 GiB logical) in 230 659 files:
> an 18.7 GiB index table plus ~47 GiB of `index/pdb/**/pdb*.ent.gz`
> structure files. lz4 buys almost nothing here, so the unpack is barely
> larger than the archive, and keeping both needs ~132 GiB. Decompressing the
> full archive measured 19 s, so extraction costs minutes, not hours — it is
> bound by writing 66 GiB and creating 230k files.

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

## Repository layout

```
docs/ASSUMPTIONS.md          assumptions register + scientific limitations
docs/FOLDDISCO_CONTRACT.md   the verified FoldDisco interface and the evidence
src/ibc_search/
    data/default.toml        defaults with units, provenance and status markers
    cli.py                   scan / show-config / folddisco-info
    config.py                loading, overrides, validation
    structure_io.py          input resolution, fetching, parsing, quality flags
    residue_id.py            stable identity + original/internal/query mapping
    exposure.py              SASA and RSA
    ball_search_algo.py      the extraction algorithm
    folddisco_contract.py    query emission and hit round-trip (no search)
    outputs.py               tables and the query file
    manifest.py              run manifest
    ext/                     documented extension points for later stages
tests/                       122 offline tests + 1 live-network test
```

## Related

`../folddisco_test/` holds the earlier FoldDisco experiment (a 4CHA index and
`folddisco_test_motif_search.sh`) that was used as a live index to verify the
query contract. It is left untouched.
