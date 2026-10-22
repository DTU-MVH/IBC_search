# FoldDisco handoff contract

**This component runs no FoldDisco search.** It defines and validates the
interface, because the interface constrains the data model.

Everything below was obtained by invoking the FoldDisco build installed on
this machine. Nothing is inferred from documentation alone, and nothing is
assumed.

## Installed build

```
$ which folddisco
/usr/local/bin/folddisco

$ folddisco version
2a756d9522d5bab8f41bc73ec72ec2327e31cd97
```

The version is reported as a bare commit hash. That hash is recorded in
`ibc_search.folddisco_contract.FOLDDISCO_VERIFIED_COMMIT` and in every run
manifest, and `ibc-search folddisco-info` reports whether the currently
installed build still matches. **If it does not match, every rule on this page
must be re-verified.**

Pre-existing experiments in `../folddisco_test/` (a 4CHA index plus
`folddisco_test_motif_search.sh`) were the starting point and were reused as a
live index for these checks.

## Query token grammar

A query is a comma-separated list of tokens. A token is:

```
<single-character chain ID><non-negative integer residue number>
```

The parser takes the **first character** as the chain and parses the
**remainder as an integer**. Verified:

| Query | Result | Conclusion |
|---|---|---|
| `B57,B102,C195` | 2 hit rows | the canonical form works |
| `B57,B102` | 2 hit rows | 2 residues is enough |
| `B57` | header only, **zero rows** | 1 residue yields nothing |
| `B57-60` | 4 nodes: `B57,B58,B59,B60` | `-` is an inclusive range separator |
| `57,B102` | header only, zero rows, **no panic** | numeric chain IDs parse; chain `5` simply does not exist in 4CHA |
| `B57A` | panic, `Invalid residue: ParseIntError{InvalidDigit}` | **insertion codes cannot be expressed** |
| `AAA57` | panic, `Invalid residue: ParseIntError{InvalidDigit}` | **multi-character chain IDs cannot be expressed** |
| `B-5` | panic, `Invalid start residue: ParseIntError{Empty}` | **negative residue numbers cannot be expressed** |
| `B29,…,B58` (30 tokens) | hits in ~0.02 s | long queries are accepted and fast |

Reproduce, for example:

```bash
cd ../folddisco_test
folddisco query -i test_index -p pdb_files/4CHA.pdb -q B57,B102,C195 --header --per-match
folddisco query -i test_index -p pdb_files/4CHA.pdb -q B57 --header --per-match      # zero rows
folddisco query -i test_index -p pdb_files/4CHA.pdb -q B57A --header --per-match     # panics
```

### Consequences, enforced in code

1. **A minimum of 2 residues per query** is a property of the method, not a
   preference: FoldDisco hashes residue *pairs*, so a single node produces no
   hash. Encoded as the constant `MIN_QUERY_RESIDUES = 2`; the config
   validator refuses `min_patch_size < 2`.
2. **Residue identity is lossy in the token form.** A residue with an
   insertion code, a multi-character chain ID or a negative residue number
   **cannot be expressed**. Such residues are marked
   `folddisco_expressible = false` with a reason in `residues.tsv`, and any
   patch containing one is marked likewise in `patches.tsv`, gets
   `emit_query = false`, and is excluded from the query file. It is never
   mangled into a token the tool would reject or misread.
3. **Token collisions are reported, not resolved.** Residues `H10`, `H10A` and
   `H10B` all project to the token `H10`. `ResidueIndex.author_key_collisions()`
   finds these, the run manifest lists them, and `resolve_query_token` refuses
   to pick one arbitrarily.
4. **Ranges are never emitted.** Although `B57-60` works, this package always
   writes explicit tokens, so the `query_residues` column FoldDisco echoes
   back is directly comparable to what was sent.

## Query file format

`-q` accepts a path instead of an inline string. The file format is one query
per line:

```
<structure_path>\t<comma-separated tokens>
```

Verified: with this form `-p` is **not** required.

```bash
printf 'pdb_files/4CHA.pdb\tB57,B102,C195\n' > q.tsv
folddisco query -i test_index -q q.tsv --header --per-match
```

A first attempt with a file containing only `B57,B102,C195` (no path) failed
with `Failed to read structure from file: B57,B102,C195`, which is what
identified the required two-field layout.

`ibc-search` writes exactly this file as `folddisco_queries.tsv`. **It was
validated by feeding it to the real binary**, which consumed it and returned
hits, including an exact self-match at RMSD 0.0.

### The query file cannot express a model number

Each line is `<structure_path>\t<tokens>` and nothing else. There is no field
for a model, so **a multi-model file is read at its first model regardless of
which model was scanned.**

Verified the hard way. Sixteen runs were made over the sixteen models of the
NMR ensemble 2N74, each scanned with `--model N`, each query file pointing at
the one deposited file. Every run returned **byte-identical RMSDs** (e.g.
0.4072 / 0.2570 / 0.2882 for the same query in models 1 and 7): all sixteen
searches had silently used model 1. Splitting the ensemble into sixteen
single-model files and repeating the query returned entirely different hits
per model (model 1 → `7v5e`, model 7 → `2mk9`), confirming the cause.

Nothing in the output reveals this. The run manifest faithfully records
`model = 7`, the patches are genuinely extracted from model 7, and only the
*search* silently falls back to model 1 — so motif-level results are correct
while hit-level results are wrong, which is the worst possible failure shape.

`ibc-search scan` therefore warns whenever the input holds more than one
model, naming the model scanned and stating that a search will read the
first. To search a specific model, extract it to its own single-model file
and scan that.

### Patch identity

The query file format has no field for a patch ID, so
`folddisco_queries.index.tsv` carries the mapping:

| column | meaning |
|---|---|
| `line_number` | 1-based line in `folddisco_queries.tsv` |
| `patch_id` | the patch that produced the query |
| `structure_path` | path written on that line |
| `radius_angstrom`, `scale`, `n_members` | provenance of the patch |
| `query_string` | the tokens, duplicated for a direct join |

A later search stage joins hits back to patches through this file.

## Reading hits back

Relevant `--format-output` per-match keys: `qid`, `tid`, `nid`, `db_key`,
`node_count`, `idf`, `rmsd`, `matching_residues`, `query_residues`,
`tm_score`, `gdt_ts`, `gdt_ha`, `chamfer_distance`, `hausdorff_distance`,
`u_matrix`, `t_vector`, `matching_coordinates`.

### `matching_residues` is NOT reliably positionally aligned

An initial 30-token query suggested clean positional alignment with `_` for
unmatched positions:

```
matching_residues: _,_,_,_,G155,G156,_,_,_,_,_,_,_,G191,_,_,_,_,_,_,_,_,_,F31,_,_,_,_,F40,_
query_residues:    B29,B30,B31,...,B58
```

**That generalisation is false.** All 545 emitted queries were run against the
4CHA test index and all 4 871 resulting hit rows were classified:

| pattern | rows |
|---|---|
| same length as query, contains `_` | 2 221 |
| **shorter** than query, contains `_` | 1 577 |
| same length as query, fully matched | 786 |
| **shorter** than query, no `_` at all | 287 |

Only **62%** of rows have equal-length fields. Unmatched positions are
sometimes padded with `_` and sometimes **omitted outright, including leading
ones**. The same query produced both of these rows:

```
query:    B146,C217,C218,C219,C221
match:    C217,C218,C219,C221        <- first position dropped, not padded
match:    B75,B76,_,_                <- two padded, one omitted
```

**Consequence:** per-position residue correspondence **cannot be recovered
from these two columns alone** in this build.

Therefore:

- `parse_hit_residues()` is the safe reader. It returns the matched target
  tokens as a set, the placeholder count, and a `positionally_alignable` flag.
  It never raises. Verified on all 4 871 real rows: 0 failures.
- `align_hit_to_query()` gives per-position mapping but **refuses** rows where
  the lengths disagree, rather than inventing a correspondence.
- For a trustworthy per-residue mapping, use `--superpose` /
  `matching_coordinates` / `u_matrix` / `t_vector` instead of these two
  columns.

The query side is unaffected and fully reliable: all **111 930** query tokens
across those rows resolved back to a residue identity of the scanned structure
with zero failures, because `query_residues` is an exact echo of what was
sent.

Also note: with several queries in one file, FoldDisco prints a **header block
per query**. A parser must not assume a single header line.

## The full round trip

```
original structure numbering          (auth_chain_id, auth_seq_id, ins_code)
        |  ResidueId  -- never renumbered
        v
internal index                        (seq_index, 0-based, file order)
        |  to_query_token  -- validates representability, refuses if lossy
        v
FoldDisco query token                 "B57"
        |  folddisco query -q file
        v
hit row: query_residues (exact echo)  +  matching_residues (target residues)
        |
        |  query side: resolve_query_token + ResidueIndex.find_author
        |              -> reliable; 111930/111930 resolved on real output
        |
        |  target side: parse_hit_residues -> matched tokens as a SET
        |               (align_hit_to_query only for the 62% of rows whose
        |                field lengths agree; --superpose otherwise)
        v
back to the original residue identity
```

## Deliberately NOT conflated

Four parameters of similar magnitude are kept separate throughout, and a
fifth, from the project notes, is not computed at all:

| # | Parameter | Where | Default |
|---|---|---|---|
| 1 | SASA solvent probe radius | `[exposure].probe_radius_angstrom` | 1.4 Å |
| 2 | patch-internal residue adjacency cutoff | `[connectivity].contact_cutoff_angstrom` | 5.0 Å |
| 3 | patch extraction ball radius | `[ball_search].radii_angstrom` | 4/5/6/7 Å |
| 4 | FoldDisco search tolerances | `[folddisco].distance_threshold_angstrom`, `.angle_threshold_degrees` | 0.5 Å, 5.0° |
| — | antibody–antigen contact cutoff (~5 Å, notes) | **not computed anywhere** | — |

Parameter 2 is an intra-patch graph-adjacency parameter. It happens to default
to the same magnitude as the notes' antibody–antigen contact cutoff, but it is
a different quantity and must not be read as one: no antibody–antigen contact
is computed in this component.

Parameter 4 is recorded for the downstream stage only; `ball_search_algo` does
not use it.

## Running the search: `scripts/folddisco_search.sh`

The handoff is executed by `scripts/folddisco_search.sh`, either directly or
through `ibc-search folddisco-search <results_dir>`, which is a thin
delegation to it so the invocation has a single definition. The search is
still performed by the FoldDisco binary; nothing reimplements it in Python.

```bash
./scripts/folddisco_search.sh status                    # what is installed/extracted
./scripts/folddisco_search.sh extract                   # unpack the downloaded index
./scripts/folddisco_search.sh search results/8yhu_binder
```

### An uncapped per-match search over the whole PDB does not terminate

With the extracted whole-PDB index (230 655 structures), timings for a single
3-residue query:

| cap | time |
|---|---|
| `--top 10` | 0.18 s |
| `--top 100` | 0.41 s |
| `--top 1000` | 4.19 s |
| uncapped | **no output after 15 minutes** |

The hash lookup is not the bottleneck: the same query with
`--per-structure --skip-match` returns in **0.18 s**. The cost is the match
stage, which loads and superposes every candidate structure, and a short
motif matches a large fraction of the PDB.

The script therefore passes `--top 1000` by default
(`[folddisco].search_top`). This is a **retrieval cutoff**, not a performance
setting: a true match ranked below the cap is never reported. On a real scan
the cap saturates — all three smoke-test queries returned exactly 1000 rows —
so every query is truncated, and how retrieval quality varies with the cap is
uncharacterised. `--no-top` removes it.

Index parameters were checked for compatibility rather than assumed:
`pdb_folddisco.type` reports `hash_type = "PDBTrRosetta"`, `grid_width = 20.0`,
`num_bin_dist = 0`, `num_bin_angle = 0` — identical to the 4CHA test index the
contract was verified against.

The lookup table stores structure paths relative to the index directory
(`pdb/ij/pdb1ij0.ent.gz`), and they resolve correctly regardless of the
working directory the query is run from. Verified.

### `-t` above 1 corrupts output. Measured, not suspected.

**This build does not write stdout line-atomically.** With `-t > 1`, lines
from concurrent threads interleave mid-line and produce rows with the wrong
number of fields.

Five identical runs of the same 219 queries at `-t 10`, each over ~1225 hit
rows, produced **2, 0, 0, 2 and 4** torn rows — varying run to run. Three runs
at `-t 1` were **byte-identical** to one another with zero torn rows.

Worse, tearing is not always detectable: a few torn rows kept a plausible
field count but carried a mangled `query_residues` value, so they no longer
joined to any emitted query.

Consequently:

- the default is `-t 1`, in the script and in
  `[folddisco].search_threads`. This is a correctness default, not an
  oversight;
- rows whose field count disagrees with the header are written to
  `folddisco_hits.malformed.tsv` and **dropped**, never parsed;
- `--strict` turns any torn row into a non-zero exit.

Raising `-t` trades reproducibility for speed, and silently loses roughly
0.5% of hits. Do not use it for numbers you intend to report.

### `-o` is ignored

`folddisco query -o <path>` created no file on this build; output went to
stdout regardless. The script uses shell redirection instead.

### Joining hits back to patches

The script joins on the **query string**, not on output order. Order is not
usable: FoldDisco emits one header block per query, but at `-t > 1` those
blocks arrive in *completion* order, so counting header blocks to recover a
query's line number is correct single-threaded and silently wrong otherwise.

`query_residues` is a safe join key precisely because it is an exact echo of
what was sent (unlike `matching_residues`, above).

One wrinkle: a query string does **not** uniquely identify a patch. The same
member set extracted at two radii yields two patch records with one query
string — 219 queries from `8yhu_binder` carried 217 distinct strings. Those
hits are genuinely indistinguishable in FoldDisco's output, so the script
reports `patch_id` as `P_r5|P_r8` and says how many rows were affected, rather
than picking one.

## Known hazards for the search stage

- **The CLI panics** on a malformed token: a Rust panic, non-zero exit, no
  structured error. Queries from this package are pre-validated, so a panic
  means a contract regression and should be surfaced loudly.
- `--sampling-ratio` / `--sampling-count` trade sensitivity for speed on long
  queries. How retrieval quality varies with query length is
  **uncharacterised**.
- `idf` is a retrieval statistic over reference-index hash frequencies. It is
  not calibrated, not a probability, and not a measure of antibody
  recognition.
- An index must be built with the same hash type and bin counts
  (`-y`, `-d`, `-a`) that queries assume.
- `folddisco index --id` (`pdb`/`uniprot`/`afdb`/`relpath`/`abspath`)
  determines the `tid` in hit output, and therefore how a hit traces back to a
  reference entry. The default `relpath` is path-dependent and a poor choice
  for a long-lived reference index.
