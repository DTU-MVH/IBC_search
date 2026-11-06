"""Shared ball placement: the member universe, the centres, and the patches.

Extracted from :mod:`ibc_search.dumb_search` so that more than one method can
place balls without reimplementing the geometry. Three methods now share it:
``dumb_search`` (one ball per residue), its seeded variants, and
``few_balls_search`` (annotated seeds plus a spatially dispersed background
set). What differs between them is **which residues become seeds**; everything
downstream of that choice -- the neighbour search, the membership rule, the
per-radius patch assembly, the FoldDisco contract -- is identical, and
identical here means literally the same code rather than three copies that
drift.

The membership rule, stated once: a residue joins a ball when any of its
eligible heavy atoms lies within the radius of the ball centre, and the atom
recorded against it is its closest such atom. Nothing about a residue's own
exposure, connectivity or identity enters that decision; methods that filter
on those do it by choosing seeds, not by changing membership.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Any, Callable

import numpy as np
from scipy.spatial import cKDTree

from .ball_search_algo import (
    BallCenter,
    Patch,
    PatchMember,
    _finalise_contract,
    _member_set_hash,
    _patch_metrics,
    _segment_structure,
    compute_ball_center,
)
from .config import scale_label
from .errors import Missing
from .exposure import ExposureResult
from .structure_io import (
    RESIDUE_CLASS_NONSTANDARD,
    RESIDUE_CLASS_STANDARD,
    ParsedStructure,
)

POLYMER_CLASSES = frozenset({RESIDUE_CLASS_STANDARD, RESIDUE_CLASS_NONSTANDARD})


@dataclass
class MemberUniverse:
    """The residues a ball may collect, and the tree used to find them.

    Built once per structure and reused across radii and methods. ``tree`` is
    indexed in the same order as ``eligible_coords``, so a neighbour index
    from the tree maps to a residue through ``eligible_res_of_atom`` and to a
    global atom through ``eligible_atom_idx``; keeping all three together is
    what makes that mapping impossible to get wrong at the call site.
    """

    tree: cKDTree
    eligible_coords: np.ndarray
    eligible_atom_idx: np.ndarray
    eligible_res_of_atom: np.ndarray
    eligible_res: np.ndarray
    """Boolean mask over ``parsed.residues``: may this residue be a member."""


def build_member_universe(
    parsed: ParsedStructure, config: dict[str, Any]
) -> tuple[MemberUniverse, dict[str, float]]:
    """Which residues can be ball members, plus a KD-tree over their atoms.

    The residue-class switches applied here say what counts as part of the
    protein, which is a parsing question rather than a filtering one: every
    method that places balls answers it the same way.
    """
    bs = config["ball_search"]
    t0 = time.perf_counter()

    eligible_res = np.zeros(len(parsed.residues), dtype=bool)
    for i, rec in enumerate(parsed.residues):
        if rec.residue_class not in POLYMER_CLASSES:
            continue
        if (
            rec.residue_class == RESIDUE_CLASS_NONSTANDARD
            and not bs["include_nonstandard_as_member"]
        ):
            continue
        if bs["members_restricted_to_scanned_chains"] and not rec.is_scanned_chain:
            continue
        eligible_res[i] = True

    eligible_atom_idx = np.flatnonzero(eligible_res[parsed.atom_res_index])
    if eligible_atom_idx.size == 0:
        raise ValueError(
            "no eligible heavy atoms for neighbour search; check chain "
            "selection and [ball_search] residue-class switches"
        )
    universe = MemberUniverse(
        tree=cKDTree(parsed.coords[eligible_atom_idx]),
        eligible_coords=parsed.coords[eligible_atom_idx],
        eligible_atom_idx=eligible_atom_idx,
        eligible_res_of_atom=parsed.atom_res_index[eligible_atom_idx],
        eligible_res=eligible_res,
    )
    return universe, {"build_kdtree_s": time.perf_counter() - t0}


def seedable_residue_indices(
    parsed: ParsedStructure, config: dict[str, Any]
) -> list[int]:
    """Residue indices eligible to CARRY a ball centre, before any seed list.

    Eligibility to seed is not the same question as eligibility to be a
    member: a seed must sit on a scanned chain and pass the seed-side
    residue-class switch. Whether a centre can actually be placed is a
    separate, geometric matter -- see :func:`place_centers`.
    """
    bs = config["ball_search"]
    out: list[int] = []
    for i, rec in enumerate(parsed.residues):
        if not rec.is_scanned_chain:
            continue
        if rec.residue_class not in POLYMER_CLASSES:
            continue
        if (
            rec.residue_class == RESIDUE_CLASS_NONSTANDARD
            and not bs["include_nonstandard_as_seed"]
        ):
            continue
        out.append(i)
    return out


def place_centers(
    parsed: ParsedStructure,
    config: dict[str, Any],
    *,
    seed_keys: set[tuple[str, int]] | None,
) -> tuple[dict[int, BallCenter], list[int], int, list[str], dict[str, float]]:
    """Put a ball centre on every seedable residue the seed list allows.

    Returns the centres, the walk order over them, how many seedable residues
    could take no centre, warnings, and timings.

    ``seed_keys`` of ``None`` means no restriction -- every seedable residue
    is a seed. An empty set means the opposite: nothing is a seed. Those are
    different requests and are kept distinct, because conflating them turns
    "this annotation matched nothing" into "scan the whole protein".

    The walk is in internal residue index order, which is file order within a
    chain, and does not depend on worker count. That is what makes a run
    reproducible.
    """
    t0 = time.perf_counter()
    centers: dict[int, BallCenter] = {}
    walk: list[int] = []
    warnings: list[str] = []
    n_no_center = 0

    for i in seedable_residue_indices(parsed, config):
        rec = parsed.residues[i]
        if seed_keys is not None and (
            rec.rid.auth_chain_id,
            rec.rid.auth_seq_id,
        ) not in seed_keys:
            continue
        center = compute_ball_center(parsed, i, config["ball_search"]["center_mode"])
        if center is None:
            # No CB, no complete backbone, no CA. The walk steps over this
            # residue rather than inventing a position for it; it can still
            # be a member of its neighbours' balls.
            n_no_center += 1
            warnings.append(
                f"residue {rec.rid.label()} ({rec.rid.res_name}) has no CB, no "
                "complete backbone and no CA; no ball centre could be placed "
                "and the walk skipped it (no coordinates were fabricated)"
            )
            continue
        centers[i] = center
        walk.append(i)

    n_virtual = sum(1 for c in centers.values() if c.is_virtual)
    if n_virtual:
        warnings.append(
            f"{n_virtual} of {len(centers)} ball centres are VIRTUAL "
            "(constructed from backbone geometry, mostly glycine); they are "
            "flagged per patch in the output"
        )
    return (centers, walk, n_no_center, warnings,
            {"place_centers_s": time.perf_counter() - t0})


def build_patches(
    parsed: ParsedStructure,
    exposure: ExposureResult,
    config: dict[str, Any],
    universe: MemberUniverse,
    centers: dict[int, BallCenter],
    walk: list[int],
    *,
    id_marker: str,
    seed_class: Callable[[int], str] | None = None,
) -> tuple[list[Patch], dict[str, Any], list[str], dict[str, float]]:
    """One ball per centre, per configured radius, as :class:`Patch` records.

    *id_marker* is the single letter that goes into every ``patch_id`` of this
    run, so tables produced by different methods can never collide on an id.

    *seed_class* optionally labels each patch by what kind of seed it was
    placed on, given that seed's residue index. Methods that place only one
    kind of seed leave it unset and the column stays blank; the point of the
    label is methods that place several, where the whole analysis depends on
    telling them apart afterwards.
    """
    bs = config["ball_search"]
    workers = int(config["performance"]["workers"])
    max_size = int(bs["max_patch_size"])
    warnings: list[str] = []
    timings: dict[str, float] = {}
    per_radius: dict[str, Any] = {}
    all_patches: list[Patch] = []

    if not walk:
        return all_patches, per_radius, warnings, timings

    center_array = np.array([centers[i].position for i in walk])

    for radius in [float(r) for r in bs["radii_angstrom"]]:
        t_radius = time.perf_counter()
        neighbour_lists = universe.tree.query_ball_point(
            center_array, r=radius, workers=workers
        )
        radius_patches: list[Patch] = []
        n_oversized = 0

        for step, seed_i in enumerate(walk):
            hits = np.asarray(neighbour_lists[step], dtype=np.int64)
            if hits.size == 0:
                continue
            center = center_array[step]
            dists = np.linalg.norm(universe.eligible_coords[hits] - center, axis=1)
            res_of_hit = universe.eligible_res_of_atom[hits]

            # Nearest hit atom per residue: sort by (residue, distance) and
            # take the first entry of each residue block. A residue joins on
            # its closest heavy atom.
            order = np.lexsort((dists, res_of_hit))
            sorted_res = res_of_hit[order]
            first = np.concatenate(([True], sorted_res[1:] != sorted_res[:-1]))
            sel_positions = order[first]
            sel_res = sorted_res[first]
            sel_atom_global = universe.eligible_atom_idx[hits[sel_positions]]
            sel_dist = dists[sel_positions]

            members = tuple(
                PatchMember(
                    residue_index=int(res_i),
                    selecting_atom_name=parsed.atom_name_at(int(atom_g)),
                    selecting_atom_distance=float(dist),
                    rsa=exposure.rsa_value_or_none(int(res_i)),
                    rsa_missing_reason=missing_rsa_reason(exposure, int(res_i)),
                    is_seed=(int(res_i) == seed_i),
                )
                for res_i, atom_g, dist in sorted(
                    zip(sel_res, sel_atom_global, sel_dist),
                    key=lambda t: int(t[0]),
                )
            )

            residue_indices = [m.residue_index for m in members]
            n_seg, seg_lengths, gaps, chains = _segment_structure(
                parsed, residue_indices
            )
            diameter, rg, compactness = _patch_metrics(
                parsed, residue_indices, radius
            )
            if len(members) > max_size:
                n_oversized += 1
            radius_patches.append(
                Patch(
                    patch_id="",
                    radius=radius,
                    scale=scale_label(config, radius),
                    center_mode=bs["center_mode"],
                    members=members,
                    seed_residue_indices=(seed_i,),
                    seed_center_kinds=(centers[seed_i].kind,),
                    seed_class="" if seed_class is None else seed_class(seed_i),
                    member_set_hash=_member_set_hash(parsed, residue_indices),
                    n_segments=n_seg,
                    segment_lengths=seg_lengths,
                    author_numbering_gaps=gaps,
                    chains=chains,
                    spatial_diameter=diameter,
                    radius_of_gyration=rg,
                    compactness_ratio=compactness,
                    # Nothing is dropped by these methods, so these are
                    # structurally zero rather than coincidentally zero.
                    n_dropped_buried=0,
                    n_dropped_disconnected=0,
                    n_components_before_connectivity=0,
                )
            )

        per_radius[f"{radius:g}"] = {
            "n_balls_placed": len(walk),
            "n_patches": len(radius_patches),
            "n_oversized_flagged": n_oversized,
            "seconds": time.perf_counter() - t_radius,
        }
        if n_oversized:
            warnings.append(
                f"radius {radius:g} A: {n_oversized} motif(s) exceed "
                f"max_patch_size={max_size}; they are kept in the table but no "
                "FoldDisco query is emitted for them"
            )
        all_patches.extend(radius_patches)
        timings[f"balls_r{radius:g}_s"] = per_radius[f"{radius:g}"]["seconds"]

    # IDs and the FoldDisco contract. No canonical re-sort: the patches are
    # already in (radius, walk order), which is the order the method is
    # defined by, and that order does not depend on worker count.
    t0 = time.perf_counter()
    structure_id = parsed.source.structure_id
    for n, patch in enumerate(all_patches):
        patch.patch_id = (
            f"{structure_id}_m{parsed.model}_r{patch.radius:g}"
            f"_{id_marker}{n:05d}"
        )
        _finalise_contract(parsed, patch, max_size)
    timings["finalise_s"] = time.perf_counter() - t0

    return all_patches, per_radius, warnings, timings


def missing_rsa_reason(exposure: ExposureResult, res_i: int) -> str:
    """Why this residue has no RSA, or ``""`` when it has one.

    A missing RSA changes nothing about membership -- no threshold is applied
    there -- but the reason is still carried into the table, because a blank
    RSA column and an untestable one are different facts.
    """
    if exposure.rsa_value_or_none(res_i) is not None:
        return ""
    value = exposure.rsa[res_i]
    return value.reason if isinstance(value, Missing) else "unknown"
