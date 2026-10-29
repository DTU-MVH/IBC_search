# Review of the per-residue motif score

What the score is, what it measures, the three defects measured in it,
and where its threshold belongs.
Everything here is computed on runs in this repository; no number is quoted
from elsewhere.

## The score

1. `scan` places a ball at each solvent-exposed residue and emits the exposed
   residues inside it as a candidate discontinuous motif.
2. FoldDisco searches each motif against the whole-PDB index and returns, per
   matched structure, an **IDF** — high when the motif's residue-pair
   encodings are rare across the index.
3. `idf_zscore.py` takes `log10` of the per-motif maximum IDF and removes the
   motif-size trend, giving a z against a fitted reference population.
4. `residue_scores.py` gives each residue the **maximum** z of the motifs
   containing it.

So a residue's score means: *the rarest surface patch this residue
participates in is this many standard deviations rarer than the reference*.
It does not mean antigenic. The PDB is a structural background, not a
collection of known epitopes, and `[folddisco].epitope_reference_index_path`
is still empty and still marked `blocked` in the shipped config.

## Defect 1 — the per-residue score is substantially an order statistic

A residue in *k* motifs takes the maximum of *k* draws, and `max` of more
draws is larger even when nothing differs. Residues sit in 1–8 motifs
(mean 3.6 on the 78 designed binders), and the score climbs with that count:

| motifs containing the residue | observed mean score | E[max] of *k* random draws from the same motif pool | difference |
|---|---|---|---|
| 1 | +0.234 | +0.479 | −0.245 |
| 2 | +0.765 | +0.837 | −0.072 |
| 3 | +0.960 | +0.996 | −0.036 |
| 4 | +1.055 | +1.091 | −0.036 |
| 5 | +1.150 | +1.163 | −0.013 |
| 6 | +1.199 | +1.224 | −0.025 |
| 7 | +1.300 | +1.271 | +0.029 |

From *k* = 3 the observed curve tracks the null to within 0.03 z. The
**+1.07 z rise from one motif to seven carries essentially no information
beyond the counting**, and it is most of the distance from the population
mean to the +1.7 call threshold. Pearson r(k, score) = +0.291.

Ball coverage is a geometric property — a convex, well-exposed residue falls
inside more balls — so the per-residue score is partly a convexity measure
wearing the clothes of motif rarity.

**What it is not.** On the 40-antigen labelled benchmark, coverage does *not*
differ between true epitope and non-epitope residues (4.80 vs 4.82 motifs),
so this defect does not fabricate an epitope signal there. It injects
variance unrelated to rarity, which dilutes rather than invents.

**Correcting it is not enough.** Subtracting the null E[max | k] per residue:

| residue score | AUC mean | AUC median | above chance |
|---|---|---|---|
| max (current) | 0.5018 | 0.4971 | 20/40 |
| mean | 0.4850 | 0.4850 | 18/40 |
| max − E[max \| k] | 0.5099 | 0.5185 | 21/40 |
| coverage *k* alone | 0.4860 | 0.5039 | 20/40 |

The correction is worth making — it is free and removes a known bias — but it
moves AUC from 0.502 to 0.510. It does not rescue the method.

## Defect 2 — the threshold, not the data, sets the headline number

The call rate is extremely sensitive to where the threshold sits. On the 78
designed binders, against the radius-matched natural-protein reference:

| threshold | residues called | % of 5 793 | binders with at least one |
|---|---|---|---|
| +1.0 | 2 756 | 47.6% | 78/78 |
| +1.3 | 1 140 | 19.7% | 72/78 |
| **+1.5** | **655** | **11.3%** | **65/78** |
| **+1.7** | **381** | **6.6%** | **52/78** |
| +2.0 | 181 | 3.1% | 31/78 |
| +2.5 | 58 | 1.0% | 11/78 |
| +3.0 | 8 | 0.1% | 2/78 |

Both thresholds in bold have been run end to end. Moving the cut from
+1.7 to +1.5 -- 0.2 z -- adds 274 residues (+72%) and 13 binders.
Nothing calibrates +1.7 against observed antibody binding, so the threshold
is a reporting convention, and any statement of the form "N residues are
likely epitopes" is really a statement about where the cut was placed.

**Residue counts also overstate the evidence.** A patch is called or not as
a unit and all its residues move together, so the residue count is roughly
the patch count times the patch size:

| threshold | called patches (of 4 442) | called residues | residues per patch | binders with a call |
|---|---|---|---|---|
| +1.5 | 154 (3.5%) | 655 | 4.3 | 65/78 |
| +1.7 | 87 (2.0%) | 381 | 4.4 | 52/78 |

That is why no binder has one or two calls — the distribution jumps straight
from 0 to 3. Report called *patches* when the question is how much
independent evidence there is.

## Defect 3 — the reference population defines the answer

The z-score is relative to whatever population `fit` was run on, and the
+1.7 threshold inherits that choice entirely:

- fitted **on the run being scored**, roughly 4.5% of motifs clear +1.7 by
  construction, whatever the structures are. The headline statistic is then a
  property of the normal distribution, not of the proteins.
- fitted on an **external natural-protein reference**, the number is a real
  measurement — it can legitimately be 0% or 20%.

Only the second answers "is this patch unusual". Any reported call rate must
name its reference population, and a reference fitted at a different ball
radius is in principle a different population, because radius sets the
motif-size mix.

Measured, that last effect turns out to be small. Fitting the same 20 natural
proteins at 5.0 Å and at 5.5 Å gives nearly the same parameters:

| reference | motifs | intercept a | slope b | pooled sd |
|---|---|---|---|---|
| natural 20, r = 5.0 Å | 11 016 | 3.2961 | 0.0568 | 0.4649 |
| natural 20, r = 5.5 Å | 11 188 | 3.2821 | 0.0574 | 0.5180 |

Scoring the 78 binders against each changed the call count from 367 to 381
of 5 793 residues (6.3% vs 6.6%) and the binders with a call from 49 to 52.
So radius matching is worth doing and cheap, but it was not what was holding
the method back.

## What the score has been shown not to do

On 40 SAbDab2 antigens with observed antibody contacts, over three epitope
definitions:

| contact distance | AUC mean | precision | enrichment over base rate |
|---|---|---|---|
| 4.5 Å | 0.502 | 0.187 | 0.99x |
| 5.5 Å | 0.465 | — | 0.649x |
| 6.0 Å | 0.502 | 0.216 | 0.98x |

Exactly 20 of 40 structures above chance at 4.5 Å. Leakage control was ruled
out as the cause: reconstructing the run without excluding the self entry or
its antigen cluster gave AUC 0.503 against 0.502.

On the 78 designed ADA binders there is no observed epitope at all, so the
only external check available is the ELISA: the correlation between called
residues and measured anti-drug-antibody absorbance is r = +0.023 at the
+1.7 threshold and r = −0.035 at +1.5, across 74 constructs.

**That null is uninformative, not evidence against the score.** Blank wells
(n = 66) have sd 0.131 against a maximum construct signal of 0.485, so only
**6 of 74** constructs clear three standard deviations over blank. With six
positives the assay has no power to confirm or refute any predictor, and the
apparent direction — the six ELISA-positive constructs average fewer calls
than the rest, 2.5 against 5.2 at +1.7 and 5.2 against 8.9 at +1.5 — rests
on those six points and should not be read as a result.

Hits are also not enriched for human structures once the composition of the
searched index is accounted for. Of the 230 655 indexed entries, SIFTS maps
701 876 chains, 28.05% of them human. Among the binders' hits:

| | human share | vs background |
|---|---|---|
| all hit rows | 27.68% | 0.99x |
| all distinct matched chains | 29.81% | 1.06x |
| rows from motifs called at +1.7 | 25.73% | 0.92x |
| distinct chains, called at +1.7 | 26.59% | 0.95x |
| rows from motifs called at +1.5 | 26.37% | 0.94x |
| distinct chains, called at +1.5 | 28.40% | 1.01x |

"Most of our hits are human" would have been true and meaningless; the PDB is
28% human to begin with.

## Calibration against a known structural motif

The threshold had never been anchored to anything observed. To anchor it,
10 small single-protein PDB entries with bound Zn were scanned at the same
5.0 Å ball radius and scored against the same reference: 5 human and 5
non-human, spanning the C2H2, PHD, B-box, FYVE, CW, ZAD, zf-AN1 and yacG
families and five non-human organisms.

**The label is the entry's own.** A residue is positive iff the PDB entry
lists it in a ``struct_site`` whose ligand is ZN. Nothing is added by
proximity. That gives 60 residues over the 10 entries, all 60 of which are
also within 2.8 Å of the metal, and all 60 Cys or His (50 Cys, 10 His) in
exactly the canonical per-family signatures — C2H2 as CCHH, PHD as C4HC3,
FYVE as C8. Two caveats travel with it:

- nine of the ten ``struct_site`` records carry
  ``pdbx_evidence_code = Software``, so the PDB generated them from proximity
  to the ligand. Only 1VFY is ``Author``. This is the deposition's own list,
  not an independent expert annotation;
- **23 of the 60 are never scored**, because they are buried and the scan
  emits surface motifs only. 22 of those 23 are Cys. The calibration is
  therefore over the 37 *solvent-exposed* members of the zinc sites, and the
  base rate is 7.2% of the 517 scored residues.

**The score finds zinc sites.**

| | antibody epitopes (40 antigens) | **zinc sites (10 entries)** |
|---|---|---|
| AUC, pooled | ~0.50 | **0.729**, 95% CI [0.644, 0.807] |
| AUC, per-structure mean | 0.502 | **0.669** |
| structures above chance | 20/40 | **9/10** |
| positive mean z | +0.764 | **+1.511** |
| other mean z | +0.783 | **+1.120** |
| best enrichment | 0.99x | **2.8x** |

This is the control the method was missing. The score is **not** broken and
the pipeline is **not** miscalibrated: it detects a real, compact structural
motif at AUC 0.73. Its failure on epitopes is specific — antibody footprints
are not a distinctive local geometry in the way a zinc site is, so there is
nothing for a rarity score to lock onto.

Note the direction: zinc sites score *higher*, not lower. The zinc-finger
*fold* is common in the PDB, but the tight tetrahedral Cys/His arrangement is
geometrically distinctive, and IDF scores the arrangement rather than the
fold.

### Where to put the threshold

Sweeping the cut against the PDB-listed zinc sites
(`scripts/threshold_sweep.py`):

| threshold | called | precision | recall | F1 | MCC | enrichment |
|---|---|---|---|---|---|---|
| +1.3 | 194 | 0.129 | 0.676 | 0.216 | +0.172 | 1.80x |
| +1.5 | 147 | 0.157 | 0.622 | 0.250 | +0.208 | 2.19x |
| **+1.6** | **126** | **0.175** | **0.595** | **0.270** | **+0.227** | **2.44x** |
| +1.7 | 110 | 0.173 | 0.513 | 0.259 | +0.204 | 2.41x |
| +1.8 | 86 | 0.186 | 0.432 | 0.260 | +0.198 | 2.60x |
| +1.9 | 65 | 0.200 | 0.351 | 0.255 | +0.189 | 2.80x |
| +2.0 | 41 | 0.195 | 0.216 | 0.205 | +0.141 | 2.73x |

**+1.6 is optimal on both F1 and MCC**, and the +1.5–+1.9 band is a broad
plateau rather than a peak — every cut in it buys 2.2-2.8x over base rate.
+1.7 was therefore very nearly right; +1.5 is slightly permissive. Above
+2.0 precision collapses on a handful of calls and the cut stops being
usable. Use **+1.6** for recall-leaning work and **+1.9** when precision
matters more than coverage.

The optimum did not move when the label was tightened from "within 5.0 Å of a
Zn" (101 residues) to "listed by the entry" (60 residues): +1.6 on both F1
and MCC either way. A threshold that survives that change of definition is
not an artefact of the definition.

This calibration transfers only as far as the assumption that an epitope
behaves like a zinc site, which the antibody benchmark says it does not. It
sets the threshold for *motif rarity*, not for antigenicity.

### Where it fails — and how much of that is measurable

The per-structure spread looks like a clean family story. It mostly is not,
and the reason is sample size: **each structure contributes 2-6 positives**,
because a zinc site is 4-8 listed residues and a third of those are buried
and never scored.

| pdb | n_pos | AUC | 95% CI (bootstrap) | family |
|---|---|---|---|---|
| 1VFY | 4 | 0.899 | [0.79, 0.98] | FYVE |
| 1LV3 | 3 | 0.888 | [0.78, 0.98] | yacG |
| 1WFH | 6 | 0.860 | [0.76, 0.94] | zf-AN1 |
| 3T7L | 3 | 0.697 | **[0.16, 0.99]** | FYVE |
| 3O70 | 5 | 0.693 | [0.43, 0.92] | PHD |
| 1PZW | 3 | 0.691 | **[0.14, 0.99]** | ZAD |
| 3Q1D | 5 | 0.566 | [0.27, 0.85] | B-box |
| 3IUF | 3 | 0.524 | [0.27, 0.79] | C2H2 |
| 1U85 | 3 | 0.506 | [0.22, 0.80] | C2H2 |
| 4Z0O | 2 | 0.365 | [0.20, 0.54] | CW |

Several intervals span almost the whole range. **No per-structure ranking
from this panel is supported**, and anything built on one — including a
family story — is reading noise. Only the pooled figure has usable
precision: AUC 0.729, 95% CI [0.644, 0.807], n_pos = 37.

**The organism split is not established.** Non-human minus human is +0.198.
Bootstrapping *residues* gives [+0.019, +0.363], which excludes zero — but
residues within a structure are not independent, so that interval is
anticonservative by pseudoreplication. At the correct unit, permuting the
organism label across the 10 structures, p = 0.104. With n = 5 a side,
nothing is demonstrated either way.

**The family split is post hoc.** Partitioning into cross-brace/loop sites
(FYVE, zf-AN1, yacG, ZAD, PHD) against classic C2H2 / B-box / CW gives mean
AUC 0.788 against 0.490, permutation p = 0.004. That p is not usable: of all
210 possible 6/4 splits of these structures, **this one ranks first by gap
size**, which is what choosing a partition after seeing the outcome
produces. It is a hypothesis for a larger panel, not a result.

### What the panel does show about the failure mode

One structure-level pattern is measurable rather than selected. Where the
score works, the scored site residues are markedly less exposed than the rest
of the surface (median RSA 0.19-0.22 against 0.43-0.60); where it fails, that
contrast disappears — in 3IUF the site residues are as exposed as everything
else (0.41 against 0.45). And in the small, zinc-dominated domains there is
barely any non-site left to contrast against: 70% of 3Q1D and 62% of 3O70 lie
within 8 Å of a metal, against 22% of 1PZW and 28% of 1LV3.

### How much of the zinc signal is Defect 1

Zinc-site residues sit in more emitted motifs than the rest — 4.30 against
3.78 — so the order-statistic bias of the MAX rule pushes them up for free:

| predictor | pooled AUC |
|---|---|
| max motif z (the score) | 0.729 |
| **ball coverage alone** | **0.620** |
| max z − E[max \| k] (artifact removed) | 0.705 |

Coverage alone separates zinc sites at 0.62, so a large part of the raw
number is reachable without looking at rarity at all. But removing the bias
costs only 0.024 AUC (0.729 → 0.705), so the rarity signal is real and
largely independent of the artifact. Both statements are needed: the score
genuinely detects zinc sites, *and* a naive reading of 0.729 credits rarity
with work that the counting did.

### Querying the zinc motif directly

The sweep above routes the question through the surface scan: which balls
happened to cover the site, then MAX back onto residues. Handing FoldDisco
the residues the entry annotates, *as the motif*, removes both steps. Ten
queries, one per structure, self-hits excluded, scored against the same
reference (`scripts/zinc_motif_query.py`).

**What they are compared against matters.** An earlier version of this figure
used each structure's own surface-scan motifs as the background and reported
a percentile within it. That number measures nothing: the scan's other motifs
are not a negative class — any of them may be a genuinely rare motif — so a
structure that happens to contain several pushes its zinc motif down the
ranking for reasons unconnected to the zinc motif, and the figure ends up
reporting a property of the protein's size and surface. The comparison used
now is the **reference population the z-score is fitted on**: 11 016 motifs
from 20 natural proteins, mean z = 0.001, sd 1.001. It is fixed, external and
shared by every structure, and a percentile against it is the z-score
restated on a scale that reads more easily.

| pdb | size | z | percentile of the reference population |
|---|---|---|---|
| 1WFH | 8 | +2.117 | 99.5 |
| 3O70 | 8 | +1.973 | 99.2 |
| 3T7L | 8 | +1.962 | 99.1 |
| 1VFY | 8 | +1.749 | 98.5 |
| 1U85 | 4 | +1.479 | 96.9 |
| 3Q1D | 8 | +1.438 | 96.5 |
| 1LV3 | 4 | +1.329 | 95.4 |
| 1PZW | 4 | +1.197 | 93.5 |
| 4Z0O | 4 | +0.398 | 63.8 |
| **3IUF** | 4 | **−3.900** | **0.6** |

Median z +1.458, at the 97th percentile of the reference; 4 of 10 clear
+1.6 and 3 clear the reference's own 99th percentile (+1.910). Eight of ten
sit above the 93rd percentile — and one sits at the 0.6th.

### Why the canonical sites score lowest

That bimodality is not noise, and it inverts the obvious reading. This part
of the analysis uses the per-ion queries (15 motifs of 4 residues each),
because match quality is a property of a single query and a two-zinc
structure would otherwise mix two sites into one number. Three of those
sites sit at z ≈ −3.9, and they are matched by **all 39 retrieved structures
at full node count with 0.14-0.52 Å RMSD** — near-perfect copies of the same geometry all
over the PDB. The high scorers are matched partially and loosely: 31-41% full
node count, RMSD up to 4.8 Å.

| | lowest 3 (z ≈ −3.9) | highest 3 (z ≥ +2.1) |
|---|---|---|
| hits at full node count | 100% | 31-77% |
| mean RMSD of hits | 0.14-0.52 Å | 0.42-1.23 Å |

Across all 15 sites, r(z, fraction of full-node hits) = **−0.496** and
r(z, mean RMSD) = **+0.411**.

**IDF rewards a motif for being hard to match, so a textbook-perfect zinc
finger is penalised precisely for being textbook-perfect.** This is the score
working as designed, not failing — but it means a high z is evidence of
geometric *unusualness*, and a canonical, highly conserved functional site
will score low. Any downstream use that reads "high z" as "important site"
has the relationship backwards for exactly the sites most likely to matter.

It also qualifies the earlier residue-level statement that zinc sites "score
higher". They do on average, through the surface scan; at motif level the
distribution is bimodal, and which mode a site falls into is set by how
common its geometry is.

## Recommendations

1. Subtract E[max | k] in `residue_scores.py`, or score patches rather than
   residues and stop projecting onto residues at all.
2. Never fit the z-score on the run being scored; keep fitted parameters in
   `scoring_params/` paired with their index, `--top` **and ball radius**.
3. Use **+1.6** as the default call threshold. It is the F1 and MCC optimum
   against observed zinc sites and sits mid-plateau, so it is not a knife
   edge. Say which reference and which radius any reported rate used.
4. Treat rarity as a feature, not a predictor of antigenicity. The zinc
   panel shows the score works on motif geometry; the antibody panel shows
   that is not what an epitope is. An epitope reference index remains the
   blocking dependency.
5. Do not read per-structure or per-family numbers off a panel this size.
   Each structure here carries 2-6 positives and its AUC interval can span
   half the range; only the pooled figure is usable. A panel built to compare
   families or organisms needs far more sites per group, chosen before the
   scores are looked at.
