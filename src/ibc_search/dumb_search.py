"""``dumb_search`` -- the simplest possible motif enumeration.

What it does, in one sentence
-----------------------------
Walk the scanned chains from the N-terminal residue to the C-terminal residue,
place one ball at every residue on the way, and emit the residues inside each
ball as a motif.

Why it exists
-------------
:mod:`ibc_search.ball_search_algo` makes a series of defensible-but-arguable
choices -- an RSA threshold for seeds, an RSA threshold for members, a
contact-graph connectivity check, deduplication of identical member sets, a
minimum patch size. Each is a knob, and each knob is a place where the motif
set could have been shaped by the parameter rather than by the structure.

This module removes every one of them. There is one rule and no thresholds:
**one ball per residue, in sequence order, keeping whatever falls inside.**
That makes it a baseline. If a downstream result holds here, it does not
depend on the ball search's filtering; if it only holds there, the filtering is
doing the work, and that is worth knowing.

It is deliberately dumb and should stay that way. Resist adding filters here:
the filtered version already exists next door.

What is reused and what is not
------------------------------
Ball **geometry** is shared with the real ball search, not reimplemented:
:func:`~ibc_search.ball_search_algo.compute_ball_center` places the centre,
the same ``[ball_search].radii_angstrom`` and ``center_mode`` apply, and
membership is still "a residue's nearest heavy atom lies within the radius".
The output objects are the same :class:`~ibc_search.ball_search_algo.Patch`
records, so every table, query file and downstream script reads a dumb-search
run exactly as it reads a scan.

What differs is only *which* balls are drawn and *what survives* them:

===========================  ==================  ====================
step                         ball_search         dumb_search
===========================  ==================  ====================
seeds                        RSA >= threshold    every residue
member exposure rule         RSA >= threshold    none
surface connectivity         component of seed   none
deduplication                identical sets      none
min/max patch size           dropped/flagged     neither
===========================  ==================  ====================

RSA is still computed and still recorded per member, because the output tables
carry it and it remains useful to read. It simply decides nothing here.

What a dumb-search motif is NOT
-------------------------------
It is a ball of residues. It is **not** a surface patch -- a ball centred on a
buried residue is full of buried residues, and nothing here filters those out.
It is **not** an epitope, and it carries no more evidence of antigenicity than
anything else in this repository: none.
"""

from __future__ import annotations

import csv
import time
from pathlib import Path
from typing import Any

import numpy as np
from scipy.spatial import cKDTree

from .ball_search_algo import (
    BallCenter,
    BallSearchResult,
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

METHOD_NAME = "dumb_search"


def load_seed_residues(path: str, structure_id: str
                       ) -> set[tuple[str, int]] | None:
    """``{(chain, res_id)}`` to seed balls on, or ``None`` for every residue.

    The file is a TSV. Chain and residue number are read from
    ``chain``/``res_id``, or ``auth_chain_id``/``auth_seq_id``, whichever it
    carries. Two optional columns are honoured when present:

    ``pdb_id`` or ``instance``
        rows are filtered to the structure being scanned, so one annotation
        table can drive a whole panel without being split per structure.
    ``is_annotated``
        only rows equal to ``1`` are used, so a table that lists every
        residue with a flag works as-is.

    Returning ``None`` for an empty path is what keeps the unseeded method
    the default: an absent list means "no restriction", never "no seeds".
    """
    if not path:
        return None
    src = Path(path).expanduser()
    if not src.is_file():
        raise ValueError(f"seed residue file not found: {src}")

    keys: set[tuple[str, int]] = set()
    with src.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = set(reader.fieldnames or ())
        chain_col = "chain" if "chain" in fields else "auth_chain_id"
        res_col = "res_id" if "res_id" in fields else "auth_seq_id"
        id_col = next((c for c in ("pdb_id", "instance", "structure_id")
                       if c in fields), None)
        if chain_col not in fields or res_col not in fields:
            raise ValueError(
                f"{src} needs chain/res_id (or auth_chain_id/auth_seq_id) "
                f"columns; found {sorted(fields)}")
        for row in reader:
            if id_col and row.get(id_col, "").strip().upper() != \
                    structure_id.strip().upper():
                continue
            if "is_annotated" in fields and row.get("is_annotated") != "1":
                continue
            try:
                keys.add((row[chain_col], int(row[res_col])))
            except (KeyError, ValueError):
                continue
    return keys


def expand_seeds_one_round(
    parsed: ParsedStructure,
    seed_keys: set[tuple[str, int]],
    *,
    radius: float,
    center_mode: str,
    tree: cKDTree,
    eligible_res_of_atom: np.ndarray,
    workers: int,
) -> tuple[set[tuple[str, int]], int]:
    """One more round of seeding: the members of the current seeds' balls.

    Returns the seed set GROWN by every eligible residue that falls inside a
    ball placed on a current seed, together with how many current seeds could
    not be expanded because no centre could be placed on them.

    Why a second round is a different question
    ------------------------------------------
    Seeding on annotated residues alone asks "what does the neighbourhood of
    each annotated residue look like". The neighbourhoods it returns are
    centred on the annotation, so the motifs it scores are the annotation's
    own, and a residue one step outside the annotated set is scored only
    through balls centred elsewhere -- or, at a chain terminus, not at all.

    Adding a round asks instead "what do the neighbourhoods of the annotated
    SITE look like", the site being the annotated residues plus everything
    packed against them. Every residue of that site gets a ball of its own
    and is therefore scored from its own centre.

    The expansion is MEMBERSHIP-based, using the same rule as a ball's
    membership: a residue joins when any of its eligible heavy atoms lies
    inside the radius. (The patches additionally record WHICH atom was
    nearest, but that cannot change who is a member.) So the residues added
    in round 2 are exactly the residues round 1's motifs already contained
    -- nothing is invented here, the existing members are promoted from
    members to seeds.

    The radius is the LARGEST configured radius. Balls at one centre are
    nested as the radius grows, so the largest is the union over all of them
    and no smaller radius can add a residue this misses.

    It is ONE round, applied once. Iterating to a fixed point would walk
    outward across the whole protein a shell at a time and end at "every
    residue", which is the unseeded method with extra steps.
    """
    current = sorted(
        (i, rec) for i, rec in enumerate(parsed.residues)
        if (rec.rid.auth_chain_id, rec.rid.auth_seq_id) in seed_keys
    )
    centers: list[np.ndarray] = []
    n_no_center = 0
    for i, _rec in current:
        center = compute_ball_center(parsed, i, center_mode)
        if center is None:
            # Cannot be expanded from, exactly as it cannot be walked from.
            # Counted rather than silently ignored: it is a seed the caller
            # asked for whose neighbourhood was never collected.
            n_no_center += 1
            continue
        centers.append(center.position)
    if not centers:
        return set(seed_keys), n_no_center

    grown = set(seed_keys)
    neighbour_lists = tree.query_ball_point(
        np.array(centers), r=radius, workers=workers)
    for hits in neighbour_lists:
        idx = np.asarray(hits, dtype=np.int64)
        if idx.size == 0:
            continue
        # The same membership rule the patches use, so a residue is added
        # here if and only if it would appear as a member there.
        for res_i in np.unique(eligible_res_of_atom[idx]):
            rec = parsed.residues[int(res_i)]
            grown.add((rec.rid.auth_chain_id, rec.rid.auth_seq_id))
    return grown, n_no_center


def run_dumb_search(
    parsed: ParsedStructure,
    exposure: ExposureResult,
    config: dict[str, Any],
) -> BallSearchResult:
    """Place one ball at every residue of every scanned chain.

    Residues are visited in internal index order, which is file order within a
    chain: the walk therefore starts at the first modelled residue of a chain
    and ends at its last. Author numbering is not consulted, so numbering gaps
    and insertion codes cannot reorder or interrupt the walk.

    Deterministic by construction: the output order is the walk order, and no
    step depends on worker count.
    """
    bs = config["ball_search"]
    center_mode = bs["center_mode"]
    workers = int(config["performance"]["workers"])
    max_size = int(bs["max_patch_size"])

    warnings: list[str] = []
    timings: dict[str, float] = {}
    polymer_classes = {RESIDUE_CLASS_STANDARD, RESIDUE_CLASS_NONSTANDARD}

    # --- member universe: every polymer residue -------------------------
    #
    # The same residue-class switches as the ball search, because they say
    # what counts as part of the protein, which is a parsing question rather
    # than one of the filtering knobs this module exists to remove.
    t0 = time.perf_counter()
    eligible_res = np.zeros(len(parsed.residues), dtype=bool)
    for i, rec in enumerate(parsed.residues):
        if rec.residue_class not in polymer_classes:
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
    eligible_coords = parsed.coords[eligible_atom_idx]
    eligible_res_of_atom = parsed.atom_res_index[eligible_atom_idx]
    tree = cKDTree(eligible_coords)
    timings["build_kdtree_s"] = time.perf_counter() - t0

    # --- optional seed restriction ---------------------------------------
    #
    # With no seed list this is the method as described: one ball per
    # residue. With one, balls are placed ONLY on the named residues, which
    # is a different question -- "what does the neighbourhood of these
    # particular residues look like" rather than "what does every
    # neighbourhood look like" -- and the manifest says which was asked.
    #
    # Only the SEEDS are restricted. Membership is unchanged: a ball still
    # collects every eligible residue within the radius, named or not, so a
    # motif is the seed's real neighbourhood and not a subgraph of the list.
    # The id the seed table is filtered on. It defaults to the parsed
    # structure id, which is the file stem -- but a biological-assembly file
    # is named `11KP-assembly1` while its annotation table says `11KP`, and
    # silently matching nothing there would seed no balls at all. So the
    # caller may name the id explicitly, and the panel does.
    seed_structure_id = (bs.get("seed_structure_id")
                         or parsed.source.structure_id)
    seed_keys = load_seed_residues(bs.get("seed_residues_path") or "",
                                   seed_structure_id)
    # Not `or 1`: that would read an explicit 0 as "unset" and silently run
    # one round, when 0 is a request that cannot be honoured.
    raw_rounds = bs.get("seed_rounds")
    seed_rounds = 1 if raw_rounds is None else int(raw_rounds)
    if seed_rounds < 1:
        raise ValueError(f"seed_rounds must be >= 1, got {seed_rounds}")
    if seed_rounds > 1 and seed_keys is None:
        raise ValueError(
            "seed_rounds > 1 needs a seed residue list: with no list every "
            "residue is already a seed, so there is nothing for a further "
            "round to add")
    # Captured BEFORE any expansion round, so "requested" keeps meaning the
    # list the caller handed in rather than what the rounds grew it into.
    n_seed_requested = len(seed_keys) if seed_keys is not None else None
    n_seed_unmatched = 0
    n_seed_round1_present = None
    n_seed_added_by_expansion = None
    n_seed_not_expandable = 0
    if seed_keys is not None:
        wanted = set(seed_keys)
        found = {
            (rec.rid.auth_chain_id, rec.rid.auth_seq_id)
            for rec in parsed.residues
        }
        missing = wanted - found
        if missing:
            n_seed_unmatched = len(missing)
            sample = ", ".join(f"{c}{s}" for c, s in sorted(missing)[:5])
            warnings.append(
                f"{len(missing)} of {len(wanted)} requested seed residues are "
                f"not present in this structure and were skipped (e.g. "
                f"{sample}). They are NOT counted as seeds that found nothing."
            )
        if not (wanted & found):
            warnings.append(
                "none of the requested seed residues is present in this "
                "structure; no motifs were extracted"
            )

        # --- further rounds of seeding --------------------------------
        #
        # Round 1's seeds are the named residues. Each further round
        # promotes the MEMBERS of the current seeds' balls to seeds too, so
        # the walk covers the annotated site rather than only the annotated
        # residues. See expand_seeds_one_round for why this is a different
        # question and why it is not iterated to a fixed point.
        present = wanted & found
        n_seed_round1_present = len(present)
        if seed_rounds > 1:
            radii = [float(r) for r in bs["radii_angstrom"]]
            expansion_radius = max(radii)
            grown = set(seed_keys)
            for _round in range(seed_rounds - 1):
                grown, n_no_expand = expand_seeds_one_round(
                    parsed, grown,
                    radius=expansion_radius,
                    center_mode=center_mode,
                    tree=tree,
                    eligible_res_of_atom=eligible_res_of_atom,
                    workers=workers,
                )
                n_seed_not_expandable += n_no_expand
            n_seed_added_by_expansion = len(grown) - n_seed_requested
            seed_keys = grown
            if len(radii) > 1:
                warnings.append(
                    f"seed expansion used the largest configured radius "
                    f"({expansion_radius:g} A) for all {len(radii)} radii: "
                    "balls at one centre are nested, so the largest is the "
                    "union over them"
                )

    # --- the walk: one centre per residue, N-terminus to C-terminus ------
    t0 = time.perf_counter()
    centers: dict[int, BallCenter] = {}
    walk: list[int] = []
    n_no_center = 0
    for i, rec in enumerate(parsed.residues):
        if not rec.is_scanned_chain:
            continue
        if seed_keys is not None and (
            rec.rid.auth_chain_id, rec.rid.auth_seq_id
        ) not in seed_keys:
            continue
        if (
            rec.residue_class == RESIDUE_CLASS_NONSTANDARD
            and not bs["include_nonstandard_as_seed"]
        ):
            continue
        if rec.residue_class not in polymer_classes:
            continue
        center = compute_ball_center(parsed, i, center_mode)
        if center is None:
            # No CB, no complete backbone, no CA. The walk steps over this
            # residue rather than inventing a position for it; it can still be
            # a member of its neighbours' balls.
            n_no_center += 1
            warnings.append(
                f"residue {rec.rid.label()} ({rec.rid.res_name}) has no CB, no "
                "complete backbone and no CA; no ball centre could be placed "
                "and the walk skipped it (no coordinates were fabricated)"
            )
            continue
        centers[i] = center
        walk.append(i)
    timings["place_centers_s"] = time.perf_counter() - t0

    n_virtual = sum(1 for c in centers.values() if c.is_virtual)
    if n_virtual:
        warnings.append(
            f"{n_virtual} of {len(centers)} ball centres are VIRTUAL "
            "(constructed from backbone geometry, mostly glycine); they are "
            "flagged per patch in the output"
        )

    stats: dict[str, Any] = {
        "method": METHOD_NAME,
        # Recorded because it changes what the run MEANS: unseeded is one
        # ball per residue, seeded is one ball per named residue only.
        "seeded": seed_keys is not None,
        "n_seed_residues_requested": n_seed_requested,
        "n_seed_residues_not_in_structure": n_seed_unmatched,
        "seed_residues_path": bs.get("seed_residues_path") or None,
        "seed_structure_id": seed_structure_id if seed_keys is not None else None,
        # Rounds of seeding. 1 = the named residues only. 2 = those plus
        # every residue inside a ball placed on one of them, so the walk
        # covers the annotated SITE and not just the annotation.
        "seed_rounds": seed_rounds if seed_keys is not None else None,
        "n_seed_residues_round1_present": n_seed_round1_present,
        "n_seed_residues_added_by_expansion": n_seed_added_by_expansion,
        "n_seed_residues_after_expansion": (len(seed_keys)
                                            if seed_keys is not None else None),
        "n_round1_seeds_without_a_centre": n_seed_not_expandable,
        "n_seeds": len(walk),
        "n_seeds_skipped_no_center": n_no_center,
        "n_virtual_centers": n_virtual,
        "per_radius": {},
    }

    if not walk:
        warnings.append(
            "no residue of a scanned chain could take a ball centre; no "
            "motifs were extracted"
        )
        return BallSearchResult(
            patches=[],
            centers=centers,
            seed_residue_indices=[],
            warnings=warnings,
            timings=timings,
            stats=stats,
        )

    center_array = np.array([centers[i].position for i in walk])

    # --- one ball per centre, per radius ---------------------------------
    all_patches: list[Patch] = []
    for radius in [float(r) for r in bs["radii_angstrom"]]:
        t_radius = time.perf_counter()
        neighbour_lists = tree.query_ball_point(
            center_array, r=radius, workers=workers
        )

        radius_patches: list[Patch] = []
        n_oversized = 0
        for step, seed_i in enumerate(walk):
            hits = np.asarray(neighbour_lists[step], dtype=np.int64)
            if hits.size == 0:
                continue
            center = center_array[step]
            dists = np.linalg.norm(eligible_coords[hits] - center, axis=1)
            res_of_hit = eligible_res_of_atom[hits]

            # Nearest hit atom per residue: sort by (residue, distance) and
            # take the first entry of each residue block. Same membership rule
            # as the ball search -- a residue joins on its closest heavy atom.
            order = np.lexsort((dists, res_of_hit))
            sorted_res = res_of_hit[order]
            first = np.concatenate(([True], sorted_res[1:] != sorted_res[:-1]))
            sel_positions = order[first]
            sel_res = sorted_res[first]
            sel_atom_global = eligible_atom_idx[hits[sel_positions]]
            sel_dist = dists[sel_positions]

            members = tuple(
                PatchMember(
                    residue_index=int(res_i),
                    selecting_atom_name=parsed.atom_name_at(int(atom_g)),
                    selecting_atom_distance=float(dist),
                    rsa=exposure.rsa_value_or_none(int(res_i)),
                    rsa_missing_reason=_missing_reason(exposure, int(res_i)),
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
                    center_mode=center_mode,
                    members=members,
                    seed_residue_indices=(seed_i,),
                    seed_center_kinds=(centers[seed_i].kind,),
                    member_set_hash=_member_set_hash(parsed, residue_indices),
                    n_segments=n_seg,
                    segment_lengths=seg_lengths,
                    author_numbering_gaps=gaps,
                    chains=chains,
                    spatial_diameter=diameter,
                    radius_of_gyration=rg,
                    compactness_ratio=compactness,
                    # Nothing is dropped here, so these are structurally zero
                    # rather than coincidentally zero.
                    n_dropped_buried=0,
                    n_dropped_disconnected=0,
                    n_components_before_connectivity=0,
                )
            )

        stats["per_radius"][f"{radius:g}"] = {
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
        timings[f"dumb_search_r{radius:g}_s"] = (
            stats["per_radius"][f"{radius:g}"]["seconds"]
        )

    # --- IDs and the FoldDisco contract ----------------------------------
    #
    # No canonical re-sort: the patches are already in (radius, walk order),
    # which is the order the method is defined by, and the walk order does not
    # depend on worker count. IDs carry a `d` to mark the method, so a
    # dumb-search table and a scan table never collide.
    t0 = time.perf_counter()
    structure_id = parsed.source.structure_id
    for n, patch in enumerate(all_patches):
        patch.patch_id = (
            f"{structure_id}_m{parsed.model}_r{patch.radius:g}_d{n:05d}"
        )
        _finalise_contract(parsed, patch, max_size)
    timings["finalise_s"] = time.perf_counter() - t0

    stats["n_patches_total"] = len(all_patches)
    stats["n_patches_discontinuous"] = sum(
        1 for p in all_patches if p.is_discontinuous
    )
    stats["n_patches_query_emitted"] = sum(1 for p in all_patches if p.emit_query)

    return BallSearchResult(
        patches=all_patches,
        centers=centers,
        seed_residue_indices=walk,
        warnings=warnings,
        timings=timings,
        stats=stats,
    )


def _missing_reason(exposure: ExposureResult, res_i: int) -> str:
    """Why this residue has no RSA, or ``""`` when it has one.

    A missing RSA changes nothing about membership here -- no threshold is
    applied -- but the reason is still carried into the table, because a blank
    RSA column and an untestable one are different facts.
    """
    if exposure.rsa_value_or_none(res_i) is not None:
        return ""
    value = exposure.rsa[res_i]
    return value.reason if isinstance(value, Missing) else "unknown"
