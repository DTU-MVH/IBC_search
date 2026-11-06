"""Choosing the background balls: where the control is allowed to sit.

``few_balls_search`` places one ball on every annotated residue and an equal
number elsewhere on the same protein. *Elsewhere* is the whole question, and
there is no single right answer to it -- so the choice is a named mode,
recorded per run, rather than a constant buried in the extraction code.

The two modes answer different questions
----------------------------------------
``protein``
    the background is a fair sample of **the protein**: candidates are
    stratified into buried, exposed and RSA-unknown and each stratum gets a
    share proportional to its own share of the candidate pool. Asks "do the
    annotated sites look unlike this protein in general".

``surface``
    the background is a fair sample of **the protein's surface**: only
    residues at or above the exposure threshold may be chosen. Asks "do the
    annotated sites look unlike the rest of the surface".

Why ``surface`` exists
----------------------
For a label that is intrinsically a surface feature -- a protein-protein
interface, an antibody epitope -- ``protein`` is not a fair control. MEASURED
across the completed runs, the two classes land at opposite ends of the
exposure scale wherever the annotation is extreme:

    dataset              annotated mean RSA   background mean RSA
    true zinc fingers                 0.094                 0.581
    reviewed enzymes                  0.132                 0.420
    antibody targets                  0.358                 0.340

Zinc ligands and catalytic residues are buried while their background is
mostly exposed, so a difference in those runs confounds "is annotated" with
"is buried". The antibody targets happened to come out matched, but that was
luck of composition, not design. ``surface`` makes the match deliberate for
surface labels: both classes are then drawn from the solvent-accessible
residues and burial cannot explain a difference between them.

What ``surface`` does NOT do
----------------------------
It equalises the POOL, not the distribution. Annotated interface residues may
still be more exposed on average than a random surface residue, because
``surface`` is a threshold rather than a matched sample. It is a control for
"buried versus exposed", not for exposure as a continuous variable, and the
realised RSA of both classes is reported per run so the residual difference
is visible rather than assumed away.

It also cannot invent surface where there is none. A small or deeply buried
candidate pool yields fewer background balls than asked for, and that
shortfall is recorded and warned about rather than quietly filled from the
interior -- silently falling back would defeat the point of asking.
"""

from __future__ import annotations

import statistics as st
from typing import Any

import numpy as np

from .ball_placement import seedable_residue_indices
from .exposure import ExposureResult
from .structure_io import ParsedStructure

MODE_PROTEIN = "protein"
MODE_SURFACE = "surface"
MODES = (MODE_PROTEIN, MODE_SURFACE)

STRATUM_BURIED = "buried"
STRATUM_EXPOSED = "exposed"
STRATUM_UNKNOWN = "rsa_unknown"


def strata_of(candidates: list[int], exposure: ExposureResult,
              threshold: float) -> dict[str, list[int]]:
    """Split candidate residues by solvent accessibility.

    RSA-unknown is a stratum of its own rather than being folded into either
    side or dropped: a residue whose exposure could not be computed is a
    different fact from a buried one, and dropping it would quietly remove
    whole regions (unresolved side chains, parse failures) from the
    background universe.
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


def allocate(strata: dict[str, list[int]], total: int) -> dict[str, int]:
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

    if total >= len(sizes):
        for k in sorted(sizes):
            if alloc[k] == 0:
                donor = max((d for d in sizes if alloc[d] > 1),
                            key=lambda d: (alloc[d], d), default=None)
                if donor is None:
                    break
                alloc[donor] -= 1
                alloc[k] += 1
    return alloc


def farthest_point_select(positions: np.ndarray, candidate_rows: list[int],
                          quota: int, chosen: np.ndarray | None) -> list[int]:
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
        mind = np.minimum(mind, np.linalg.norm(pts - pts[best], axis=1))
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

    Returns the chosen residue indices and an audit record: the mode, the
    per-stratum request and delivery, the dispersion actually achieved, the
    realised exposure of both classes, and the shortfall if the protein could
    not supply the request.
    """
    fbs = config.get("few_balls_search", {})
    threshold = float(fbs.get("background_rsa_threshold", 0.20))
    stratify = bool(fbs.get("background_stratify_by_burial", True))
    mode = str(fbs.get("background_mode", MODE_PROTEIN))
    if mode not in MODES:
        raise ValueError(
            f"background_mode must be one of {MODES}, got {mode!r}")

    annotated = set(annotated_indices)
    candidates = [i for i in seedable_residue_indices(parsed, config)
                  if i in centers_by_index and i not in annotated]

    audit: dict[str, Any] = {
        "background_mode": mode,
        "n_background_requested": n_wanted,
        "background_rsa_threshold": threshold,
    }

    n_before = len(candidates)
    if mode == MODE_SURFACE:
        # Only residues CONFIRMED exposed. RSA-unknown is excluded rather
        # than assumed exposed: the mode's whole claim is that every
        # background ball sits on the surface, and an uncomputable RSA
        # cannot support that claim.
        candidates = [i for i in candidates
                      if (exposure.rsa_value_or_none(i) or -1.0) >= threshold]
        audit["n_candidates_before_surface_filter"] = n_before
        # Stratifying a pool that is exposed by construction would put every
        # ball in one stratum; saying so beats leaving a misleading `true`
        # in the manifest.
        stratify = False

    audit["n_background_candidates"] = len(candidates)
    audit["background_stratified_by_burial"] = stratify

    if not candidates or n_wanted <= 0:
        audit.update({"n_background_selected": 0,
                      "n_background_shortfall": max(n_wanted, 0),
                      "per_stratum": {}})
        return [], audit

    order = {i: n for n, i in enumerate(candidates)}
    positions = np.array([centers_by_index[i] for i in candidates])
    chosen_coords = (np.array([centers_by_index[i] for i in annotated_indices])
                     if annotated_indices else np.empty((0, 3)))

    strata = (strata_of(candidates, exposure, threshold) if stratify
              else {STRATUM_EXPOSED: list(candidates)})
    quotas = allocate(strata, n_wanted)

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
            got = farthest_point_select(positions, rows_left[k], 1, against)
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
        # property of the run and not a claim made by a docstring.
        "min_separation_between_background_centres_angstrom": sep,
        "min_distance_background_to_annotated_centre_angstrom":
            nearest_annotated,
        # The realised exposure of BOTH classes. `surface` equalises the pool
        # rather than the distribution, so the residual difference belongs in
        # the output where it can be read, not in an assumption.
        "mean_rsa_annotated": _mean_rsa(exposure, annotated_indices),
        "mean_rsa_background": _mean_rsa(exposure, selected),
    })
    return selected, audit


def _mean_rsa(exposure: ExposureResult, indices: list[int]) -> float | None:
    values = [exposure.rsa_value_or_none(i) for i in indices]
    known = [v for v in values if v is not None]
    return round(st.mean(known), 4) if known else None
