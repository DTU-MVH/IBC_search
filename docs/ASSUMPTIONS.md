# Assumptions register

Status markers used throughout this repository:

| marker | meaning |
|---|---|
| **documented** | traceable to a named external reference or to the installed FoldDisco build |
| **verified** | empirically checked on this machine, with the check reproducible from this page |
| **proposed** | an engineering choice made by this implementation; not derived from literature; to be revisited by a sensitivity study |
| **blocked** | cannot be settled until an external input or decision arrives |

A **proposed** value must never be cited as an established criterion.

---

## 0. Status of the algorithm specification

**A-0.** No formal mathematical definition of the ball-scanning algorithm
("Algo 1") existed in this repository: the repository was **empty** before this
component (only `../folddisco_test/`, a 4CHA FoldDisco index and a shell
script, existed alongside it). `ball_search_algo` is therefore an explicit
**operational interpretation** of the prose brief, not an implementation of a
pre-existing spec. Status: **documented** (as an interpretation).

Divergences from the brief, and why:

| Brief | Implementation | Reason |
|---|---|---|
| "min/max patch size configurable, defaults justified" | `max_patch_size` suppresses **query emission** only; oversized patches still appear in `patches.tsv` | Dropping them would discard real extracted surface for a configuration reason. The size limit is a handoff policy, not an extraction rule. |
| "flag or drop disconnected components" | both, plus a third `keep`, selected by `[connectivity].disconnected_policy` | `keep` is needed to measure how often the connectivity check matters at all. |
| "deduplicate identical residue sets" | always on, not configurable | Emitting the identical motif twice has no interpretation; only the seed list differs, and that is retained. |

---

## 1. Input handling

**A-1.1** A token is treated as a PDB entry ID only if it matches either the
legacy form `^[0-9][A-Za-z0-9]{3}$` or the extended form
`^pdb_[0-9A-Za-z]{8}$` (case-insensitive) **and** no file of that name
exists. A local file always wins. Status: **proposed**.
*Consequence:* a 4-letter name not starting with a digit is never read as an
entry ID, so a four-letter filename cannot trigger a download.

**Both spellings are supported**, because the four-character space is close
to exhausted and the PDB is migrating to extended identifiers. Verified
against RCSB: `4HHB.cif.gz` and `pdb_00004hhb.cif.gz` decompress to
byte-identical content. An extended ID of the form `pdb_0000XXXX` therefore
canonicalises to its legacy `XXXX` equivalent for fetching, so either
spelling gives one URL, one cache entry, one SHA-256 and identical patch
IDs. An entry genuinely outside the four-character space keeps its extended
ID throughout.

*Format limitation, verified:* `pdb_00004hhb.pdb.gz` returns HTTP 404 while
the `.cif.gz` is served — the legacy PDB format cannot express an identifier
outside the four-character space, which is why extended IDs exist at all.
With `[input].fetch_format = "pdb"`, an entry having no legacy equivalent is
refused up front with a message naming mmCIF, rather than failing later on an
opaque 404. The default is mmCIF, so this does not arise by default.

*Not live-verified:* the genuinely-extended path (an ID **not** of the form
`pdb_0000XXXX`) is unit-tested but has never been exercised against a real
download, because the PDB has not yet issued such an entry. The parsing,
canonicalisation and format refusal are tested; the HTTP round trip is not.

**A-1.2** Fetched entries are retrieved as mmCIF (`.cif.gz`) from
`https://files.rcsb.org/download/<ID>.cif.gz`, with URL, UTC date, SHA-256 and
size recorded. Status: **documented**.
*Why mmCIF:* it carries both `label_asym_id` and `auth_asym_id`, and has
neither the 80-column limit nor the single-character chain limit of the PDB
format. Verified live for 1CRN.

**A-1.3** A network failure raises `FetchError`. There is **no fallback** to
another format, mirror or source, because a silent fallback would make the
manifest's provenance record untrue. Status: **documented**.

**A-1.4 Alternate conformers.** Per residue: atoms with a blank altloc are
always kept; among non-blank altlocs, the one with the greatest **summed**
occupancy wins; ties break to the **lexicographically smallest altloc ID**.
Status: **documented** (deterministic, and the only rule implemented).
*Why summed, not mean:* a conformer modelled with more atoms should not be
penalised. *Why the ID tie-break:* it removes any dependence on file order or
floating-point comparison order.

**A-1.5 Hydrogens** are discarded and counted per residue. All geometry and
SASA is heavy-atom based. Status: **proposed**.

**A-1.6 Residue classification** uses the residue *name* before the
HETATM/ATOM record type, so that modified polymer residues deposited as
HETATM (MSE and similar) are still treated as polymer. Classes:
`polymer_standard_aa`, `polymer_nonstandard_aa`, `nucleic`, `carbohydrate`,
`water`, `monoatomic_ion`, `other_hetero`. Status: **proposed**.

**A-1.7 Numbering is never altered.** Author chain ID, residue number and
insertion code are preserved exactly; the internal index is a separate,
0-based, file-order quantity. Non-contiguous, repeated and negative numbering
all survive. Status: **documented**; tested.

**A-1.8 Chain breaks** are reported as two independent facts: `seq_gap_before`
(author numbering is non-contiguous, i.e. residues missing from the model) and
`chain_break_before` (modelled residues are not peptide-bonded, C–N > 2.0 Å or
backbone atoms absent). Status: **proposed** (the 2.0 Å cutoff; a peptide C–N
bond is ~1.33 Å).

**A-1.9 Biological assembly.** Only the **deposited asymmetric unit** is
scanned. Assembly expansion is not applied. The declared assembly IDs are read
from mmCIF and recorded in the manifest. Status: **documented** limitation.
*Consequence:* a surface buried only in the biological assembly reads as
exposed here.

**A-1.10 Default model** is the first (model 1), configurable, and always
recorded along with the list of available models. Status: **proposed**.

---

## 2. Exposure (SASA / RSA)

**A-2.1** SASA is computed by Shrake–Rupley via
`biotite.structure.sasa` — an established implementation, not re-derived here.
Status: **documented**.

**A-2.2 Probe radius** defaults to **1.4 Å**, the conventional water probe
(Shrake & Rupley 1973; Lee & Richards 1971), configurable. Status:
**documented**. This is parameter (1) of four kept-separate radii.

**A-2.3 Van der Waals radii**: the ProtOr set (Tsai et al. 1999), biotite's
default, which assumes no hydrogens — consistent with A-1.5. If it fails, the
code falls back to the `Single` set **with a warning stating the values are no
longer comparable**. Status: **documented**.

**A-2.4 RSA normalisation** uses **Tien et al. 2013**, *PLoS ONE*
8(11):e80635, Table 1, **theoretical** maxima. Status: **documented**.
*Consequence:* RSA can slightly exceed 1.0 (observed max 1.06 on 4CHA),
because theoretical maxima come from an extended tripeptide context.

**A-2.5 Missing RSA is MISSING, never 0.** A residue type absent from the Tien
table (every non-standard residue) gets `rsa = MISSING` with the reason
recorded, and is **not** eligible as a seed. It is **not** counted as buried
either. Status: **documented**; tested.
*Rationale:* 0.0 is a legitimate RSA for a buried residue, so using it for
"unknown" would corrupt exposure statistics.

**A-2.6 SASA context** defaults to `scanned_chains`: SASA is computed using
only the chains being scanned. Status: **proposed**.
*This materially changes results.* With `scanned_chains`, a residue at an
interface with a non-scanned chain reads as **exposed**; with
`all_polymer_chains` it reads as buried. The choice is recorded in the
manifest and warned about at runtime.

**A-2.7 Waters are never SASA occluders.** They are bulk solvent, and a
crystallographic water shell would make the protein read as buried. Ligands,
glycans and ions are occluders only when
`[exposure].include_hetero_as_occluder = true` (default false, so a designed
binder's intrinsic surface is assessed without crystal additives). Status:
**proposed**.

---

## 3. `ball_search_algo`

**A-3.1 Seed threshold** `seed_rsa_threshold = 0.20`.
**This is an implementation choice, not an established epitope criterion.**
0.20 is a common exposed/buried split but carries no epitope meaning. Status:
**proposed**.

**A-3.2 Member threshold** `member_rsa_threshold = 0.05`, below the seed
threshold so partially occluded rim residues are retained. The validator
refuses `member > seed`, which would make every seed ineligible as a member of
its own patch. Status: **proposed**.

**A-3.3 Ball centre, `cb` mode.** Fallback ladder: observed CB → virtual CB
from backbone N/CA/C → CA. Glycine's virtual CB is labelled
`virtual_cb_glycine`; a *non-glycine* residue missing its CB is labelled
`virtual_cb_missing_cb`. These are deliberately different labels, because a
GLY without a CB is normal and a LEU without one is a defect. Status:
**proposed** (the ladder); **verified** (the construction, see A-3.4).

**A-3.4 Virtual CB construction** uses the standard ideal-geometry
coefficients (as used by trRosetta and related tools):

```
b = CA - N ;  c = C - CA ;  a = b x c
CB = -0.58273431*a + 0.56802827*b - 0.54067466*c + CA
```

Status: **verified**. The test suite rebuilds CB for every non-glycine residue
of 4CHA that has a real CB and compares: **mean deviation < 0.25 Å,
95th percentile < 0.5 Å** (asserted). A virtual centre is always flagged
`ball_center_is_virtual = true`; no coordinate is ever fabricated unflagged.

**A-3.5 Ball centre, `sidechain_centroid` mode.** Centroid of side-chain heavy
atoms, where "side chain" excludes `N, CA, C, O, OXT, OT1, OT2`. Terminal
oxygens are excluded so a chain terminus does not drag the centroid along the
backbone. A centroid built from fewer atoms than the residue type expects is
labelled `sidechain_heavy_atom_centroid_incomplete` and records the atom
count. Glycine and side-chain-free cases fall back to the A-3.3 ladder.
Status: **proposed**.

**A-3.6 Neighbour rule.** A residue joins the patch when its **minimum
eligible heavy-atom distance** to the ball centre is ≤ `patch_radius`. The
specific atom and its distance are recorded per member. All atoms of selected
residues remain available downstream (tested): a motif is a set of residues,
not a set of CB atoms. Status: **documented**.

**A-3.7 Ball radii** `[4, 5, 6, 7] Å`. **Engineering candidates for a
sensitivity study, not literature-derived optima.** Status: **proposed**.
This is parameter (3) of four kept-separate radii.

**Revision history**, the first two steps driven by A-3.7b (below):

1. The set was originally `[5, 8, 10, 12, 15] Å`. **15 Å was dropped and
   6.5 Å was added**, after A-3.7b showed 15 Å captured a majority of a real
   56-residue minibinder's surface in a single ball — not a local
   neighbourhood by any reading. 6.5 Å was added to keep four resolution
   points in the range that stays local down to minibinder scale.
2. **12 Å was then dropped too** — at 12 Å the same minibinder's median patch
   already covered 36% of the chain (max 55%), the same saturation pattern as
   15 Å, just less extreme. With it gone, every radius in the default set
   stays under 30% of the 56-residue chain even at its largest observed
   patch (see A-3.7b).
3. **8 Å and 10 Å were dropped and 4 Å and 6 Å added**, giving the current
   `[4, 5, 6, 7] Å`. Unlike steps 1 and 2 this was **not** forced by a
   saturation measurement — both 8 Å and 10 Å stayed under the 30% threshold
   that disqualified 12 Å and 15 Å. It is a deliberate narrowing of the grid
   onto the submotif scale, trading coverage of the intermediate scale for
   four closely spaced points between 4 and 7 Å. The intermediate scale is
   consequently **no longer sampled at all** by the defaults, which is a
   limitation of any sensitivity study run on this set: it can resolve
   behaviour within the submotif range but says nothing above 7 Å.

None of these changes is reflected in old output directories from earlier in
this project (`results/8yhu*`, produced under previous radius sets); re-run
`scan` to get output under the current set.

Measured patch sizes (4CHA chains A/B/C, `cb` mode, defaults), re-measured
when the set changed. The method was validated by reproducing the previous
set's figures exactly before the new ones were taken:

| radius | min | median | max |
|---|---|---|---|
| 4 Å | 3 | 4 | 5 |
| 5 Å | 3 | 4 | 8 |
| 6 Å | 3 | 6 | 10 |
| 7 Å | 5 | 7 | 13 |

For reference, the radii dropped in step 3 gave 8 Å: 5/9/16 and 10 Å: 8/14/23
on the same chains.

Every radius in the current set now sits inside or just above the 2–4 residue
"small local interaction motif" prior (medians 4–7 residues). The set no
longer contains anything approaching the ~15-residue whole-epitope prior —
the 10 Å ball that coincidentally matched it (median 14) is gone. That
coincidence was never evidence that a 10 Å ball recovers epitopes, and its
removal costs nothing in that respect.

Note that these sizes were measured on a 238-residue chain. They do **not**
transfer to a minibinder-sized input: see A-3.7b.

**A-3.7b The upper radii saturate on minibinder-sized inputs — this is why
15 Å and then 12 Å were dropped from the defaults.** Status: **verified**
(measured); this finding is the direct cause of A-3.7 revisions 1 and 2. It
is *not* the reason for revision 3, which narrowed the set to 4–7 Å by
choice rather than under measured pressure.

The original radius set was chosen as a generic exploration grid and measured
on 4CHA, a 238-residue chain. On a real minibinder it behaved very
differently. Measured on PDB **8YHU** ("hTLR3/minibinder 8.6"), author chain B
= minibinder 8.6, 56 residues, 457 heavy atoms:

* whole chain B: radius of gyration **11.0 Å**, maximum heavy-atom dimension
  **35.7 Å**.

Under the **current** radius set:

| radius | median patch | as fraction of the 56-residue chain | max patch |
|---|---|---|---|
| 4 Å | 4 res | 7% | 6 res (11%) |
| 5 Å | 5 res | 9% | 7 res (12%) |
| 6 Å | 7 res | 12% | 10 res (18%) |
| 7 Å | 9 res | 16% | 12 res (21%) |

The current set stays local on a minibinder by a wide margin: the largest
patch any default radius produces covers 21% of the chain, against the 30%
figure that the now-dropped 8 Å reached and the 55%/84% that disqualified
12 Å and 15 Å.

Under **earlier** radius sets, the same structure's larger balls gave:

| radius | median patch | as fraction of the 56-residue chain | max patch | status |
|---|---|---|---|---|
| 8 Å | 10 res | 18% | 17 res (30%) | dropped (step 3) |
| 10 Å | 15 res | 27% | 25 res (45%) | dropped (step 3) |
| 12 Å | 20 res | 36% | 31 res (55%) | dropped (step 2) |
| 15 Å | **32 res** | **57%** | **47 res (84%)** | dropped (step 1) |

The largest 15 Å patch (`8YHU_m1_r15_p00167`, seed B28, 41 residues in 9
segments) had a radius of gyration of **10.8 Å** and a spatial diameter of
**34.3 Å** — essentially the dimensions of the entire minibinder (11.0 Å /
35.7 Å). That measurement triggered dropping 15 Å; the 12 Å column above,
measured immediately after, showed the same saturation pattern one step down
and triggered dropping 12 Å as well. Both radii remain available via
`--radius 12` / `--radius 15` for anyone who wants them, with this caveat
applying to each.

Consequences, which must be respected when choosing radii:

1. A ball radius is only a *local* neighbourhood while it is small relative to
   the protein. Once the radius approaches the protein's own radius of
   gyration, the "patch" is most of the molecule. For a ~56-residue binder
   this starts to bite from about 10–12 Å upward (10 Å: 27% median / 45% max;
   12 Å: 36% median / 55% max), which is why the default set now tops out at
   10 Å rather than 12 Å or 15 Å.
2. The `scale` label (A-3.7 / `scale_labels`) is pure bookkeeping and becomes
   actively misleading at the top of a radius set: the 15 Å patch above,
   when it still carried the label `epitope_scale`, held 41 residues, roughly
   3× the ~15-residue whole-epitope prior. The current default set no longer
   has an `epitope_scale`-labelled radius at all, precisely because neither
   candidate for that label (12 Å, 15 Å) held up on a minibinder-sized input.
3. A query in the 30-55 residue range (12-15 Å territory on this binder)
   cannot correspond to an epitope-scale reference motif in any interpretable
   way, whatever retrieval score comes back. The useful range for a
   minibinder of this size is roughly **5–10 Å**, which is now exactly the
   default set's span.
4. Patch size is therefore not comparable across structures of different size
   at a fixed radius. Any cross-structure statistic must either normalise by
   protein size or fix the patch size rather than the radius.

Selecting a *production* radius within the current 5-10 Å range remains
blocked item **B-8**. Use `--radius` to probe beyond the defaults, e.g.
`--radius 12 --radius 15` to reproduce the saturated patches above.

**A-3.8 Minimum patch size** `min_patch_size = 3`. Status: **verified** floor
of 2, **proposed** value of 3. A 1-residue query returns zero rows from the
installed FoldDisco build; 3 sits above that hard floor and inside the 2–4
residue submotif prior. The validator refuses anything below 2.

**A-3.9 Maximum patch size** `max_patch_size = 60`. Status: **proposed**,
measured. A 30-residue query runs in ~0.02 s on the installed build, so query
length is not a tight tool limit; 60 keeps every radius in the current A-3.7
set comfortably emittable, with substantial headroom now that the largest
observed patch at the default set's largest radius is 23 residues (10 Å on
4CHA; see A-3.7). It was originally sized to also cover 12 Å and 15 Å radii
(observed max 32 and 55 residues respectively) that have since both been
dropped from the default set (A-3.7b), so the cap is now considerably looser
than it strictly needs to be; it was left at 60 rather than retightened,
since it is an engineering ceiling, not a sensitivity limit, and the headroom
costs nothing when `--radius` is used to probe beyond the defaults. **Open:**
how FoldDisco retrieval quality varies with query length is uncharacterised.

**A-3.10 Connectivity.** Patch members form a contact graph (edge when
minimum heavy-atom distance ≤ `contact_cutoff_angstrom`, default 5.0 Å); only
the connected component containing the seed is kept (`drop`), or the rest is
marked (`flag`), or connectivity is ignored (`keep`). Status: **proposed**.
This is parameter (2) of four kept-separate radii, and is **not** the notes'
antibody–antigen contact cutoff, which is not computed anywhere here.

**A-3.11 Seed-normal filter**, optional and **off by default**. The outward
normal is estimated as the direction from the centroid of surrounding heavy
atoms to the ball centre, and members below `seed_normal_min_cosine` (−0.2)
are removed. Status: **proposed**. Crude: reasonable on a convex surface, poor
in a concave pocket. The contact-graph check is primary.

**A-3.12 Deduplication** of identical member sets is unconditional, retaining
every contributing seed. Status: **documented**.

**A-3.13 Jaccard clustering** of strongly overlapping patches at
`jaccard_threshold = 0.8`, single-linkage, within one radius. Status:
**proposed**. Single linkage can chain — a known limitation, accepted because
the intent is only "these describe the same surface region".

**A-3.14 Segment structure.** Segments are maximal runs adjacent in the
**internal index** within one chain (i.e. "next modelled residue"), so they do
not depend on author numbering. Gaps between segments are reported in
**author numbering**, which is what a reader cares about. Status: **proposed**.

**A-3.15 Geometry metrics.** `spatial_diameter` is the largest pairwise
heavy-atom distance over the patch's residues and **can exceed 2 × radius**,
because membership is decided on a residue's *nearest* atom while its other
atoms may point outward. `radius_of_gyration` is the RMS distance of those
atoms from their centroid. `compactness_ratio = radius_of_gyration / radius`
is dimensionless and reads as how much of the ball the patch fills. Status:
**proposed** (definitions).

**A-3.16 Determinism.** Patches are canonically sorted before IDs are
assigned, so output is byte-identical regardless of `[performance].workers`.
Parallelism is confined to a single vectorised
`cKDTree.query_ball_point(..., workers=N)` call, whose results come back in
centre order. Status: **verified**; tested at 1, 2 and 4 workers.

---

## 3b. FoldDisco handoff

**A-3b.1 Query token grammar** is `<single-char chain><non-negative int>`.
Insertion codes, multi-character chain IDs and negative residue numbers
**cannot be expressed** (each verified to make the installed build panic).
Affected residues and any patch containing one are marked
`folddisco_expressible = false` with a reason and excluded from the query
file. Status: **verified**. See `docs/FOLDDISCO_CONTRACT.md`.

**A-3b.2 Query file format** is `<structure_path>\t<tokens>` per line.
Status: **verified** — and validated end to end by feeding all 545 queries
emitted for 4CHA to the real binary: exit 0, empty stderr, **zero panics**,
4 871 hit rows in 2.98 s.

**A-3b.3 `matching_residues` is NOT reliably positionally aligned to
`query_residues`.** Status: **verified**, and this **corrected an earlier
assumption made during this work**. An initial 30-residue probe showed clean
`_`-padded alignment, which suggested positional zipping was safe. Classifying
all 4 871 real hit rows disproved it: only 62% have equal-length fields, and
unmatched positions are sometimes padded with `_` and sometimes omitted
outright, including *leading* ones.
*Consequence:* `parse_hit_residues()` (set semantics, never raises, 0 failures
on all 4 871 rows) is the safe reader; `align_hit_to_query()` refuses
non-alignable rows rather than inventing a correspondence. A trustworthy
per-residue mapping needs `--superpose` / `matching_coordinates`.
The **query** side is reliable: all 111 930 query tokens resolved back to a
residue identity with zero failures.

**A-3b.4 Verified against one build only** — commit
`2a756d9522d5bab8f41bc73ec72ec2327e31cd97`. Recorded in every manifest;
`ibc-search folddisco-info` reports a mismatch. Status: **documented**.
*If the build changes, A-3b.1 to A-3b.3 must be re-verified.*

---

## 4. Scientific limitations of this extraction step

**L-1** A Euclidean ball is **not** a geodesic surface patch and **not** a
physically accessible antibody footprint. The connectivity check (A-3.10)
reduces, but does not eliminate, merging across thin regions of the protein.

**L-2** An extracted patch is a **candidate structural motif**. It is not an
epitope, and it carries no evidence of antigenicity, antibody binding or
cross-reactivity, pre-existing ADA reactivity, induction of a new immune
response, or T-cell peptide presentation.

**L-3** Patches are **massively overlapping and non-independent**: one seed per
exposed residue per radius (545 patches from 239 scanned residues across five
radii on 4CHA). Any downstream statistic must account for this; a naive hit
count will be inflated.

**L-4** The structure is treated as a **single rigid conformer**. No side-chain
flexibility, backbone dynamics or conformational ensemble is considered, and
surface residues are the most mobile part of a protein.

**L-5** Exposure is computed on the **deposited asymmetric unit** (A-1.9) in
the chosen chain context (A-2.6). Both choices change which residues are
seeds.

**L-6** Several patches matching separate references do **not** jointly
reconstitute one compatible antibody interface. They are independent
geometric matches and may be mutually inaccessible.

**L-7** Epitope chemistry spans polar, charged, hydrophobic and aromatic
residues. The notes' observation of polar/charged enrichment is **not** applied
as a filter anywhere, and must not become one.

**L-8** The seed RSA threshold (A-3.1) is the single most influential
unvalidated parameter: it alone decides which surface is scanned at all.

**L-9** **A patch is only "local" relative to the size of the protein.** A
fixed ball radius does not define a fixed biological scale: on a 56-residue
minibinder a 15 Å ball captured a median of 57% of the whole chain and a
12 Å ball 36% (A-3.7b), with the 15 Å ball's largest patch matching the
dimensions of the entire molecule. Both radii have since been removed from
the default set for exactly this reason. Patch sizes are therefore not
comparable between structures of different size at a fixed radius, and any
radius approaching a small input's own radius of gyration must not be
described as epitope-scale.

---

## 5. Blocked items

| # | Item | Blocked on |
|---|---|---|
| **B-1** | reference collection to search | Which dataset (IEDB-derived structures, SAbDab, other). Nothing selected, downloaded or inspected. |
| **B-2** | what a reference "entry" is | Whole antigen chain vs. interface-cropped vs. epitope residues only. This changes what a hit *means* and is the most consequential open decision. |
| **B-3** | `[folddisco].reference_index_path` | B-1/B-2. Marked `status: blocked` in the config. |
| **B-4** | FoldDisco search tolerances for this task | No reference index to tune against. Defaults (0.5 Å, 5.0°) recorded, unused. |
| **B-5** | negative/background set | Without a baseline, no hit count is interpretable. |
| **B-6** | redundancy control for references | Antibody–antigen data is heavily redundant; un-deduplicated references inflate any frequency statistic. |
| **B-7** | labelled evaluation set + positive-class definition | Needed before any benchmark number means anything. |
| **B-8** | which radius/radii to use in production | Requires the sensitivity study the radii set was chosen to enable. Partially informed by A-3.7b: 12 Å and 15 Å were both found non-local on a 56-residue binder and were dropped. The set was then narrowed again, by choice rather than measurement, to 4–7 Å (A-3.7 step 3), so the defaults now sample only the submotif scale and nothing above 7 Å. Still open: which of 4/5/6/7 Å (or a value in between) to use in production; whether dropping the intermediate scale entirely was right, since no default radius now tests it; and whether a size-relative rule (radius vs. the chain's radius of gyration, or fixing patch size instead of radius) should replace a fixed radius altogether. |
| **B-9** | a set of the project's own minibinder designs | Partially resolved: PDB **8YHU** (hTLR3/minibinder 8.6, author chain B, 56 residues) has been run end to end and is a genuine minibinder, alongside the reference structures 4CHA, 1CRN, 6M0J and 3HFM. Still outstanding: the project's actual de novo designs, and any binder supplied as a prediction rather than a crystal structure. |
| **B-10** | which conformational state to scan | 8YHU chain B is the TLR3-**bound** conformation, and its free-state surface may differ. Whether to scan complexed structures, free structures, or a predicted/ensemble model is undecided (see L-4, L-5 and A-2.6). |
