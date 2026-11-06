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

STRATUM_BURIED = "buried"
STRATUM_EXPOSED = "exposed"
STRATUM_UNKNOWN = "rsa_unknown"


def _strata_of(
    candidates: list[int],
    exposure: ExposureResult,
    threshold: float,
) -> dict[str, list[int]]:
    """Split candidate residues by solvent accessibility.

    RSA-unknown is a stratum of its own rather than being folded into either
    side or dropped: a residue whose exposure could not be computed is a
    different fact from a buried one, and dropping it would quietly remove
    whole regions (unresolved side chains, parse failures) from the background
    universe.
    """
    out: dict[str, list[int]] = {
        STRATUM_BURIED: [], STRATUM_EXPOSED: [], STRATUM_UNKNOWN: []}
    for i in candidates:
        rsa = exposure.rsa_value_or_none(i)
        if rsa is None:
            out[STRATUM_UNKNOWN].append(i)
        elif rsa >= threshold:
            out[STRATUM_EXPOSED].append(i)
        else:
            out[STRATUM_BURIED].append(i)
    return out


def _allocate(strata: dict[str, list[int]], total: int) -> dict[str, int]:
    """Share *total* balls between strata in proportion to their sizes.

    Largest-remainder apportionment, so the parts sum to *total* exactly and
    the rounding cannot drift. A non-empty stratum is guaranteed at least one
    ball when there is one to give, which is what makes "both internal and
    surface" a property of the method rather than of the particular protein:
    proportional shares alone would give a 96%-buried protein zero exposed
    background balls at small totals.
    """
    sizes = {k: len(v) for k, v in strata.items() if v}
    if not sizes or total <= 0:
        return {k: 0 for k in strata}
    pool = sum(sizes.values())
    alloc = {k: 0 for k in strata}

    # Floor of the fair share, then the largest remainders, then the
    # guaranteed minimum of one per non-empty stratum where budget allows.
    exact = {k: total * n / pool for k, n in sizes.items()}
    for k in sizes:
        alloc[k] = min(int(exact[k]), sizes[k])
    remaining = total - sum(alloc.values())
    order = sorted(sizes, key=lambda k: (-(exact[k] - int(exact[k])), k))
    while remaining > 0:
        progressed = False
        for k in order:
            if remaining == 0:
                break
            if alloc[k] < sizes[k]:
                alloc[k] += 1
                remaining -= 1
                progressed = True
        if not progressed:
            break          # every stratum exhausted: the shortfall is real

    if total >= len([k for k in sizes]):
        for k in sorted(sizes):
            if alloc[k] == 0:
                donor = max(
                    (d for d in sizes if alloc[d] > 1),
                    key=lambda d: (alloc[d], d), default=None)
                if donor is None:
                    break
                alloc[donor] -= 1
                alloc[k] += 1
    return alloc


def _farthest_point_select(
    positions: np.ndarray,
    candidate_rows: list[int],
    quota: int,
    chosen: np.ndarray | None,
) -> list[int]:
    """Greedy k-centre: repeatedly take the point furthest from those chosen.

    *positions* holds every candidate's centre, *candidate_rows* indexes the
    subset this call may pick from, and *chosen* is the already-selected
    coordinates to stay away from (annotated centres plus earlier picks).

    Returns rows of *positions*, in pick order. Deterministic: ties on the
    min-distance objective break toward the lower row, which is the lower
    residue index.
    """
    if quota <= 0 or not candidate_rows:
        return []
    rows = np.asarray(candidate_rows, dtype=np.int64)
    pts = positions[rows]

    if chosen is not None and len(chosen):
        mind = np.min(
            np.linalg.norm(pts[:, None, :] - chosen[None, :, :], axis=2),
            axis=1)
    else:
        # Nothing to stay away from yet. Start from the point furthest from
        # the candidate centroid -- a reproducible extremity, rather than an
        # arbitrary first row that would make the whole set depend on file
        # order.
        mind = np.linalg.norm(pts - pts.mean(axis=0), axis=1)

    picked: list[int] = []
    for _ in range(min(quota, len(rows))):
        best = int(np.argmax(mind))
        picked.append(int(rows[best]))
        d = np.linalg.norm(pts - pts[best], axis=1)
        mind = np.minimum(mind, d)
    return picked


def select_background_seeds(
    parsed: ParsedStructure,
    exposure: ExposureResult,
    config: dict[str, Any],
    *,
    annotated_indices: list[int],
    centers_by_index: dict[int, np.ndarray],
    n_wanted: int,
) -> tuple[list[int], dict[str, Any]]:
    """Choose *n_wanted* dispersed background centres away from the annotated.

    Returns the chosen residue indices and an audit record: the per-stratum
    request and delivery, the separation actually achieved, and the shortfall
    if the protein could not supply the request.
    """
    fbs = config.get("few_balls_search", {})
    threshold = float(fbs.get("background_rsa_threshold", 0.20))
    stratify = bool(fbs.get("background_stratify_by_burial", True))

    annotated = set(annotated_indices)
    candidates = [i for i in seedable_residue_indices(parsed, config)
                  if i in centers_by_index and i not in annotated]

    audit: dict[str, Any] = {
        "n_background_requested": n_wanted,
        "n_background_candidates": len(candidates),
        "background_rsa_threshold": threshold,
        "background_stratified_by_burial": stratify,
    }
    if not candidates or n_wanted <= 0:
        audit.update({"n_background_selected": 0,
                      "n_background_shortfall": max(n_wanted, 0),
                      "per_stratum": {}})
        return [], audit

    order = {i: n for n, i in enumerate(candidates)}
    positions = np.array([centers_by_index[i] for i in candidates])
    chosen_coords = (np.array([centers_by_index[i] for i in annotated_indices])
                     if annotated_indices else np.empty((0, 3)))

    strata = (_strata_of(candidates, exposure, threshold) if stratify
              else {STRATUM_EXPOSED: list(candidates)})
    quotas = _allocate(strata, n_wanted)

    # Round-robin across strata rather than draining one then the next, so
    # that global dispersion is shared out fairly instead of being spent
    # entirely on whichever stratum happened to go first.
    picked_rows: list[int] = []
    remaining = dict(quotas)
    rows_left = {k: [order[i] for i in v] for k, v in strata.items()}
    while any(remaining[k] > 0 and rows_left[k] for k in remaining):
        for k in sorted(remaining, key=lambda k: (-remaining[k], k)):
            if remaining[k] <= 0 or not rows_left[k]:
                continue
            against = (np.vstack([chosen_coords, positions[picked_rows]])
                       if picked_rows else chosen_coords)
            got = _farthest_point_select(positions, rows_left[k], 1, against)
            if not got:
                remaining[k] = 0
                continue
            picked_rows.append(got[0])
            rows_left[k] = [r for r in rows_left[k] if r != got[0]]
            remaining[k] -= 1
            break

    selected = [candidates[r] for r in picked_rows]
    delivered = {k: 0 for k in strata}
    for i in selected:
        for k, v in strata.items():
            if i in v:
                delivered[k] += 1
                break

    sep = None
    if len(picked_rows) > 1:
        pts = positions[picked_rows]
        d = np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=2)
        np.fill_diagonal(d, np.inf)
        sep = float(np.min(d))
    nearest_annotated = None
    if selected and len(chosen_coords):
        d = np.linalg.norm(
            positions[picked_rows][:, None, :] - chosen_coords[None, :, :],
            axis=2)
        nearest_annotated = float(np.min(d))

    audit.update({
        "n_background_selected": len(selected),
        "n_background_shortfall": n_wanted - len(selected),
        "per_stratum": {
            k: {"n_candidates": len(strata[k]), "n_requested": quotas[k],
                "n_selected": delivered[k]}
            for k in sorted(strata)
        },
        # The dispersion actually achieved, so "spread out" is a measured
        # property of the run and not a claim made by this docstring.
        "min_separation_between_background_centres_angstrom": sep,
        "min_distance_background_to_annotated_centre_angstrom":
            nearest_annotated,
    })
    return selected, audit


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
