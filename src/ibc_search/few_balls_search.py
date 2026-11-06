"""``few_balls_search`` -- annotated balls, plus a matched background set.

What it does, in one sentence
-----------------------------
Place one ball on every annotated residue, then place the **same number** of
further balls elsewhere on the same protein, chosen to be spread out in three
dimensions and to span buried and exposed positions alike, with a floor of
:data:`MIN_BACKGROUND_DEFAULT` background balls however few annotations there
are.

Why it exists
-------------
The methods next door answer the two ends of a question and neither answers
the middle. ``dumb_search`` unseeded places a ball on every residue, which is
a fair sample of the protein but thousands of queries per structure.
``dumb_search`` seeded places balls only on annotated residues, which is cheap
but has **no internal control**: to ask whether annotated motifs behave
differently you must compare against a different run, over different
structures, with a different number of balls per residue -- and that
comparison was measured to be dominated by those differences rather than by
the annotation.

This method puts the control *inside the same structure and the same run*. The
annotated balls and the background balls come from one protein, one radius,
one index and one search, in equal number. A difference between them cannot be
explained by fold, resolution, deposition year, motif count or search effort,
because those are shared. That is the entire point, and it is why the two
classes are labelled in the output rather than inferred afterwards.

How the background balls are chosen
-----------------------------------
Two requirements pull against each other: the background should be *spread
out*, and it should *represent the protein* rather than only its outside.

**Spread out** is greedy farthest-point sampling (the Gonzalez k-centre
heuristic): repeatedly take the candidate residue whose nearest already-chosen
centre is furthest away. The annotated centres are seeded into the chosen set
*before* the first pick, so the background is pushed away from the annotated
sites as well as away from itself. Ties break on residue index, so the result
is deterministic.

**Represents the protein** needs the stratification, because farthest-point
sampling alone is biased toward the surface: the points that are furthest
apart in a globular object are the ones on its hull. Left unstratified, a
"background" set would be a surface sample and the comparison against
annotated sites -- which are often clefts and pockets -- would confound
annotation with burial, which is exactly the artifact that normalising hit
counts by motif coverage already had to undo once. So candidates are split by
solvent accessibility into buried, exposed and RSA-unknown, each stratum is
given a share of the background balls *proportional to its own share of the
candidate pool*, and the farthest-point selection runs within each stratum
against the globally chosen set. A protein that is 70% buried gets a
background set that is about 70% buried.

What it does NOT do
-------------------
It does not filter members: a ball keeps every residue inside it, annotated or
not, buried or not, exactly as the other methods do. It does not deduplicate,
and an annotated ball and a background ball may overlap -- how much they do is
measured and reported rather than prevented, because suppressing overlap would
silently bias where the background can sit. It does not exclude residues that
are merely *inside* an annotated ball from being background centres; only the
annotated residues themselves are excluded. The number of background seeds
that landed inside an annotated ball is reported, so that choice is auditable
instead of invisible.

A ball is a ball of residues. It is not a surface patch, not an epitope, and
carries no evidence of antigenicity.
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np

from .ball_placement import (
    build_member_universe,
    build_patches,
    place_centers,
    seedable_residue_indices,
)
from . import background_selection as _bg
from .ball_search_algo import BallSearchResult, compute_ball_center
from .dumb_search import load_seed_residues
from .exposure import ExposureResult
from .structure_io import ParsedStructure

METHOD_NAME = "few_balls_search"

SEED_CLASS_ANNOTATED = "annotated"
SEED_CLASS_BACKGROUND = "background"

MIN_BACKGROUND_DEFAULT = 5
"""Floor on the background count, however few annotated residues there are.

A structure with two annotated residues would otherwise contribute two
background balls, and a two-against-two comparison within one structure is
not a comparison. The floor is a statistical minimum, not a geometric one.
"""

# The background selection lives in its own module: where the control is
# allowed to sit is a question with more than one defensible answer, and it
# is chosen per run by name. See ibc_search.background_selection.
#
# Re-exported here because callers and tests have referred to these through
# few_balls_search since before the split.
STRATUM_BURIED = _bg.STRATUM_BURIED
STRATUM_EXPOSED = _bg.STRATUM_EXPOSED
STRATUM_UNKNOWN = _bg.STRATUM_UNKNOWN
MODE_PROTEIN = _bg.MODE_PROTEIN
MODE_SURFACE = _bg.MODE_SURFACE
MODES = _bg.MODES

_strata_of = _bg.strata_of
_allocate = _bg.allocate
_farthest_point_select = _bg.farthest_point_select
select_background_seeds = _bg.select_background_seeds


def run_few_balls_search(
    parsed: ParsedStructure,
    exposure: ExposureResult,
    config: dict[str, Any],
) -> BallSearchResult:
    """Balls on the annotated residues, plus an equal dispersed background.

    Deterministic by construction: the candidate order, the stratification and
    the farthest-point tie-breaks all key on residue index, and no step
    depends on worker count.
    """
    bs = config["ball_search"]
    fbs = config.get("few_balls_search", {})
    min_background = int(fbs.get("min_background_balls",
                                 MIN_BACKGROUND_DEFAULT))
    if min_background < 0:
        raise ValueError(
            f"min_background_balls must be >= 0, got {min_background}")

    warnings: list[str] = []
    universe, timings = build_member_universe(parsed, config)

    # --- the annotated seeds ---------------------------------------------
    #
    # This method is DEFINED by an annotated set, so an absent list is an
    # error rather than a fallback to every residue: silently scanning the
    # whole protein would produce a table that looks like a few-balls run and
    # answers a different question.
    seed_structure_id = (bs.get("seed_structure_id")
                         or parsed.source.structure_id)
    seed_keys = load_seed_residues(bs.get("seed_residues_path") or "",
                                   seed_structure_id)
    if seed_keys is None:
        raise ValueError(
            "few_balls_search needs --seed-residues: the method places balls "
            "on annotated residues and an equal background set, so with no "
            "annotation there is no annotated class and nothing to match")

    present = {
        (rec.rid.auth_chain_id, rec.rid.auth_seq_id)
        for rec in parsed.residues
    }
    missing = set(seed_keys) - present
    if missing:
        sample = ", ".join(f"{c}{s}" for c, s in sorted(missing)[:5])
        warnings.append(
            f"{len(missing)} of {len(seed_keys)} annotated residues are not "
            f"present in this structure and were skipped (e.g. {sample}). "
            "They are NOT counted as seeds that found nothing."
        )

    # Centres for every seedable residue, computed once: the background count
    # must match the number of annotated balls actually PLACED, and a residue
    # that cannot take a centre places no ball.
    t0 = time.perf_counter()
    centers_by_index: dict[int, np.ndarray] = {}
    for i in seedable_residue_indices(parsed, config):
        c = compute_ball_center(parsed, i, bs["center_mode"])
        if c is not None:
            centers_by_index[i] = c.position
    timings["precompute_centers_s"] = time.perf_counter() - t0

    annotated_indices = [
        i for i, rec in enumerate(parsed.residues)
        if (rec.rid.auth_chain_id, rec.rid.auth_seq_id) in seed_keys
        and i in centers_by_index
    ]
    n_annotated = len(annotated_indices)
    if n_annotated == 0:
        warnings.append(
            "no annotated residue of this structure could take a ball centre, "
            "so this run has NO annotated class; the background balls below "
            "are not a control for anything and the run should not be pooled "
            "with structures that have both classes"
        )

    # --- the background seeds --------------------------------------------
    n_wanted = max(n_annotated, min_background)
    background_indices, audit = select_background_seeds(
        parsed, exposure, config,
        annotated_indices=annotated_indices,
        centers_by_index=centers_by_index,
        n_wanted=n_wanted,
    )
    if audit["n_background_shortfall"] > 0:
        warnings.append(
            f"asked for {n_wanted} background balls and placed "
            f"{audit['n_background_selected']}: the protein offered only "
            f"{audit['n_background_candidates']} eligible non-annotated "
            "residues with a placeable centre. The classes are therefore NOT "
            "equal in size for this structure"
        )
    if n_annotated and n_wanted > n_annotated:
        warnings.append(
            f"{n_annotated} annotated balls but {n_wanted} background balls: "
            f"the min_background_balls={min_background} floor applies, so the "
            "two classes are deliberately unequal here. Weight by class size "
            "rather than pooling raw counts"
        )

    # --- place the balls --------------------------------------------------
    # By INDEX, not by (chain, res_id). Those two are not the same selection
    # where insertion codes exist: A100 and A100A share an author number, so
    # a key-based walk picks up a residue this method never chose and which
    # therefore belongs to neither class. Measured on
    # true_annotated_antibody_targets, where it aborted a scan outright.
    centers, walk, n_no_center, place_warnings, place_timings = place_centers(
        parsed, config,
        seed_indices=set(annotated_indices) | set(background_indices),
    )
    warnings.extend(place_warnings)
    timings.update(place_timings)

    annotated_set = set(annotated_indices)
    background_set = set(background_indices)

    def seed_class(i: int) -> str:
        if i in annotated_set:
            return SEED_CLASS_ANNOTATED
        if i in background_set:
            return SEED_CLASS_BACKGROUND
        # Unreachable while the walk is built from the two sets above; a
        # blank here would quietly produce unlabelled motifs, so it raises.
        raise AssertionError(
            f"residue index {i} is in the walk but in neither seed class")

    # `f` marks the method in every patch_id, so a few-balls table cannot
    # collide with a scan (`p`) or a dumb-search (`d`) table.
    all_patches, per_radius, patch_warnings, patch_timings = build_patches(
        parsed, exposure, config, universe, centers, walk,
        id_marker="f", seed_class=seed_class,
    )
    warnings.extend(patch_warnings)
    timings.update(patch_timings)

    n_overlap = _n_background_inside_annotated_balls(
        all_patches, annotated_set, background_set)

    stats: dict[str, Any] = {
        "method": METHOD_NAME,
        "seeded": True,
        "seed_residues_path": bs.get("seed_residues_path") or None,
        "seed_structure_id": seed_structure_id,
        "n_annotated_requested": len(seed_keys),
        "n_annotated_not_in_structure": len(missing),
        "n_annotated_balls": n_annotated,
        "min_background_balls": min_background,
        "n_background_balls": len(background_indices),
        "classes_equal_size": n_annotated == len(background_indices),
        # How much the control overlaps what it controls for. Not prevented,
        # because preventing it would bias where the background may sit.
        "n_background_seeds_inside_an_annotated_ball": n_overlap,
        "n_seeds": len(walk),
        "n_seeds_skipped_no_center": n_no_center,
        "n_virtual_centers": sum(1 for c in centers.values() if c.is_virtual),
        "per_radius": per_radius,
        **audit,
    }
    stats["n_patches_total"] = len(all_patches)
    stats["n_patches_discontinuous"] = sum(
        1 for p in all_patches if p.is_discontinuous)
    stats["n_patches_query_emitted"] = sum(
        1 for p in all_patches if p.emit_query)
    for cls in (SEED_CLASS_ANNOTATED, SEED_CLASS_BACKGROUND):
        stats[f"n_patches_{cls}"] = sum(
            1 for p in all_patches if p.seed_class == cls)
        stats[f"n_patches_query_emitted_{cls}"] = sum(
            1 for p in all_patches if p.seed_class == cls and p.emit_query)

    if not walk:
        warnings.append(
            "no residue of a scanned chain could take a ball centre; no "
            "motifs were extracted")

    return BallSearchResult(
        patches=all_patches,
        centers=centers,
        seed_residue_indices=walk,
        warnings=warnings,
        timings=timings,
        stats=stats,
    )


def _n_background_inside_annotated_balls(
    patches: list, annotated: set[int], background: set[int]
) -> int:
    """How many background SEEDS fall inside some annotated ball.

    Counted over the smallest radius only, so the number does not multiply
    with the number of configured radii and means "at the tightest ball this
    run drew" rather than an unlabelled sum across scales.
    """
    if not patches or not background:
        return 0
    smallest = min(p.radius for p in patches)
    covered: set[int] = set()
    for p in patches:
        if p.radius != smallest or p.seed_class != SEED_CLASS_ANNOTATED:
            continue
        covered.update(m.residue_index for m in p.members)
    return len(covered & background)
