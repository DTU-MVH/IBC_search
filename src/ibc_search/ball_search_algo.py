"""``ball_search_algo`` -- surface ball-scanning extraction of discontinuous motifs.

Status of this specification
----------------------------
No formal mathematical definition of this algorithm exists in this repository
(the repository was empty before this component). What follows is therefore an
**operational interpretation**, written to be simple, fast and auditable. Every
choice below is recorded in ``docs/ASSUMPTIONS.md`` with its status, and every
threshold is configurable.

What it does, in one sentence
-----------------------------
For each sufficiently solvent-exposed residue (a *seed*), place a ball at a
representative position for that residue, collect the exposed residues whose
nearest heavy atom falls inside the ball, discard members that are not
surface-connected to the seed, and emit the resulting residue set as a
candidate discontinuous structural motif.

Pipeline
--------
A. exposure        -- per-residue SASA/RSA (:mod:`ibc_search.exposure`)
B. seeds           -- residues with RSA >= ``seed_rsa_threshold``
C. ball centres    -- one per seed, Cbeta or side-chain centroid
D. neighbours      -- min heavy-atom distance to centre <= ``patch_radius``
E. discontinuity   -- segments/gaps annotated, never collapsed
F. connectivity    -- contact-graph component containing the seed
G. deduplication   -- identical member sets merged, seeds retained
H. scales          -- several radii per run, their relationship recorded

What a patch is NOT
-------------------
A patch is a Euclidean neighbourhood of exposed residues. It is **not** a
geodesic surface patch, **not** a physically accessible antibody footprint, and
**not** an epitope. A ball of radius r centred near the surface can still reach
across a thin region of the protein, which is why step F exists -- and step F is
a cheap graph heuristic, not a proof of accessibility. Nothing computed here
speaks to antigenicity, antibody cross-reactivity, pre-existing anti-drug
antibody reactivity, induction of a new immune response, or T-cell peptide
presentation.
"""

from __future__ import annotations

import hashlib
import time
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from scipy.spatial.distance import pdist

from .config import scale_label
from .errors import Missing
from .exposure import ExposureResult
from .folddisco_contract import (
    MIN_QUERY_RESIDUES,
    check_representable,
    motif_query_string,
)
from .residue_id import ResidueId
from .structure_io import (
    RESIDUE_CLASS_NONSTANDARD,
    RESIDUE_CLASS_STANDARD,
    ParsedStructure,
)

SIDECHAIN_EXCLUDED_ATOM_NAMES = frozenset({"N", "CA", "C", "O", "OXT", "OT1", "OT2"})
"""Atoms treated as backbone when forming a side-chain centroid.

Terminal oxygens (OXT/OT1/OT2) are excluded so that a chain terminus does not
drag its centroid along the backbone direction.
"""

# Virtual C-beta construction constants.
#
# Given backbone N, CA, C, a C-beta position is reconstructed as
#     b = CA - N ;  c = C - CA ;  a = b x c
#     CB = -0.58273431*a + 0.56802827*b - 0.54067466*c + CA
# These are the standard ideal-geometry coefficients used by trRosetta and
# related tools. They are validated in the test suite against *real* C-beta
# atoms of non-glycine residues, which is the only honest way to justify them.
_VCB_A = -0.58273431
_VCB_B = 0.56802827
_VCB_C = -0.54067466

CENTER_KIND_CB = "cb_atom"
CENTER_KIND_VIRTUAL_CB_GLY = "virtual_cb_glycine"
CENTER_KIND_VIRTUAL_CB_MISSING = "virtual_cb_missing_cb"
CENTER_KIND_CA_FALLBACK = "ca_fallback"
CENTER_KIND_SIDECHAIN_CENTROID = "sidechain_heavy_atom_centroid"
CENTER_KIND_SIDECHAIN_PARTIAL = "sidechain_heavy_atom_centroid_incomplete"


@dataclass(frozen=True)
class BallCenter:
    """A scanning-ball centre for one seed residue."""

    residue_index: int
    position: np.ndarray
    kind: str
    is_virtual: bool
    """True when the position was constructed rather than observed.

    A virtual centre is always flagged. Coordinates are never fabricated
    without this flag being set and carried into the output tables.
    """

    n_atoms_used: int
    note: str = ""


@dataclass(frozen=True)
class PatchMember:
    """One residue inside a patch, with the evidence that put it there."""

    residue_index: int
    selecting_atom_name: str
    """The specific atom whose distance to the ball centre was smallest."""

    selecting_atom_distance: float
    """That atom's distance to the ball centre, in Angstroms."""

    rsa: float | None
    rsa_missing_reason: str
    is_seed: bool
    buried_flagged: bool = False
    """Retained despite RSA below ``member_rsa_threshold`` (policy='flag')."""

    disconnected_flagged: bool = False
    """Retained despite not being in the seed's contact component."""


@dataclass
class Patch:
    """A deduplicated candidate discontinuous structural motif."""

    patch_id: str
    radius: float
    scale: str
    center_mode: str
    members: tuple[PatchMember, ...]
    seed_residue_indices: tuple[int, ...]
    seed_center_kinds: tuple[str, ...]
    member_set_hash: str
    n_segments: int
    segment_lengths: tuple[int, ...]
    author_numbering_gaps: tuple[int, ...]
    chains: tuple[str, ...]
    spatial_diameter: float
    radius_of_gyration: float
    compactness_ratio: float
    n_dropped_buried: int
    n_dropped_disconnected: int
    n_components_before_connectivity: int
    cluster_id: str = ""
    folddisco_expressible: bool = True
    folddisco_block_reason: str = ""
    query_string: str = ""
    emit_query: bool = True
    emit_block_reason: str = ""

    @property
    def size(self) -> int:
        return len(self.members)

    @property
    def is_discontinuous(self) -> bool:
        return self.n_segments > 1

    def residue_indices(self) -> tuple[int, ...]:
        return tuple(m.residue_index for m in self.members)


@dataclass
class BallSearchResult:
    """Everything one ball-search run produced."""

    patches: list[Patch]
    centers: dict[int, BallCenter]
    seed_residue_indices: list[int]
    warnings: list[str] = field(default_factory=list)
    timings: dict[str, float] = field(default_factory=dict)
    stats: dict[str, Any] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# C. Ball centres
# ---------------------------------------------------------------------------

def _residue_atom_lookup(
    parsed: ParsedStructure, res_i: int
) -> dict[str, np.ndarray]:
    rec = parsed.residues[res_i]
    sl = slice(rec.atom_start, rec.atom_stop)
    return {
        str(name): parsed.coords[rec.atom_start + offset]
        for offset, name in enumerate(parsed.atom_names[sl])
    }


def virtual_cb(n: np.ndarray, ca: np.ndarray, c: np.ndarray) -> np.ndarray:
    """Construct an ideal-geometry C-beta position from backbone N, CA, C."""
    b = ca - n
    c_vec = c - ca
    a = np.cross(b, c_vec)
    return _VCB_A * a + _VCB_B * b + _VCB_C * c_vec + ca


def compute_ball_center(
    parsed: ParsedStructure, res_i: int, center_mode: str
) -> BallCenter | None:
    """Compute the scanning-ball centre for one residue.

    Returns ``None`` when no defensible position exists (no CB, no complete
    backbone and no CA), rather than inventing one.

    Fallback order, documented and deterministic:

    ``center_mode="cb"``
        real CB atom -> virtual CB from N/CA/C -> CA atom.
        For glycine the virtual CB is the *expected* path, not an error, and is
        labelled ``virtual_cb_glycine``.

    ``center_mode="sidechain_centroid"``
        centroid of side-chain heavy atoms -> (glycine, or no side-chain atom
        present) the ``cb`` ladder above. A centroid built from fewer atoms
        than the residue type should have is labelled
        ``sidechain_heavy_atom_centroid_incomplete`` and records how many atoms
        it used, so a truncated side chain is never mistaken for a complete one.
    """
    rec = parsed.residues[res_i]
    atoms = _residue_atom_lookup(parsed, res_i)

    if center_mode == "sidechain_centroid":
        sc_names = [
            name for name in atoms if name not in SIDECHAIN_EXCLUDED_ATOM_NAMES
        ]
        if sc_names:
            coords = np.array([atoms[name] for name in sorted(sc_names)])
            expected = rec.expected_sidechain_heavy_atoms
            incomplete = expected is not None and len(sc_names) < expected
            return BallCenter(
                residue_index=res_i,
                position=coords.mean(axis=0),
                kind=(
                    CENTER_KIND_SIDECHAIN_PARTIAL
                    if incomplete
                    else CENTER_KIND_SIDECHAIN_CENTROID
                ),
                is_virtual=False,
                n_atoms_used=len(sc_names),
                note=(
                    f"centroid of {len(sc_names)} of {expected} expected "
                    "side-chain heavy atoms"
                    if incomplete
                    else f"centroid of {len(sc_names)} side-chain heavy atoms"
                ),
            )
        # Glycine and side-chain-free cases fall through to the CB ladder.

    if "CB" in atoms:
        return BallCenter(
            residue_index=res_i,
            position=atoms["CB"],
            kind=CENTER_KIND_CB,
            is_virtual=False,
            n_atoms_used=1,
            note="observed CB atom",
        )

    if all(name in atoms for name in ("N", "CA", "C")):
        kind = (
            CENTER_KIND_VIRTUAL_CB_GLY
            if rec.rid.res_name == "GLY"
            else CENTER_KIND_VIRTUAL_CB_MISSING
        )
        return BallCenter(
            residue_index=res_i,
            position=virtual_cb(atoms["N"], atoms["CA"], atoms["C"]),
            kind=kind,
            is_virtual=True,
            n_atoms_used=3,
            note=(
                "CB constructed from backbone N/CA/C with ideal geometry; "
                "position is VIRTUAL, not observed"
            ),
        )

    if "CA" in atoms:
        return BallCenter(
            residue_index=res_i,
            position=atoms["CA"],
            kind=CENTER_KIND_CA_FALLBACK,
            is_virtual=False,
            n_atoms_used=1,
            note=(
                "no CB and incomplete backbone; fell back to the CA atom, which "
                "sits further from the side-chain surface than a CB would"
            ),
        )

    return None


# ---------------------------------------------------------------------------
# E. Discontinuity annotation
# ---------------------------------------------------------------------------

def _segment_structure(
    parsed: ParsedStructure, residue_indices: list[int]
) -> tuple[int, tuple[int, ...], tuple[int, ...], tuple[str, ...]]:
    """Describe the sequence-segment structure of a residue set.

    Segments are maximal runs that are adjacent in the **internal index**
    within one chain. The internal index follows file order, so it expresses
    "next modelled residue of this chain" without depending on author
    numbering, which may contain gaps, insertion codes or non-contiguous runs.

    Gaps between consecutive segments of the same chain are reported in
    **author numbering**, because that is the quantity a reader cares about
    (how many residues of the deposited sequence are skipped).
    """
    by_chain: dict[str, list[int]] = {}
    for res_i in residue_indices:
        rid = parsed.residues[res_i].rid
        by_chain.setdefault(rid.auth_chain_id, []).append(res_i)

    n_segments = 0
    lengths: list[int] = []
    gaps: list[int] = []
    for chain in sorted(by_chain):
        members = sorted(by_chain[chain])
        runs: list[list[int]] = [[members[0]]]
        for prev, cur in zip(members, members[1:]):
            if cur == prev + 1:
                runs[-1].append(cur)
            else:
                runs.append([cur])
        n_segments += len(runs)
        lengths.extend(len(r) for r in runs)
        for prev_run, next_run in zip(runs, runs[1:]):
            prev_last = parsed.residues[prev_run[-1]].rid.auth_seq_id
            next_first = parsed.residues[next_run[0]].rid.auth_seq_id
            gaps.append(int(next_first - prev_last - 1))

    return n_segments, tuple(lengths), tuple(gaps), tuple(sorted(by_chain))


# ---------------------------------------------------------------------------
# F. Connectivity
# ---------------------------------------------------------------------------

def _patch_atom_indices(parsed: ParsedStructure, residue_indices: list[int]) -> np.ndarray:
    """All heavy-atom indices of the given residues.

    Extraction keeps every atom of a selected residue available: the motif is
    a set of residues, not a set of C-beta atoms.
    """
    slices = [
        np.arange(parsed.residues[i].atom_start, parsed.residues[i].atom_stop)
        for i in residue_indices
    ]
    return np.concatenate(slices) if slices else np.empty(0, dtype=np.int64)


def _connected_component_of_seed(
    parsed: ParsedStructure,
    residue_indices: list[int],
    seed_index: int,
    contact_cutoff: float,
) -> tuple[list[int], list[int], int]:
    """Split patch members into the seed's contact component and the rest.

    Two residues are adjacent when their closest heavy atoms are within
    ``contact_cutoff``. This is an *intra-patch graph adjacency* parameter. It
    is deliberately NOT the ~5 A antibody--antigen contact cutoff from the
    project notes: no antibody--antigen contact is computed anywhere in this
    component.

    Returns ``(kept, rejected, n_components)``.
    """
    if len(residue_indices) <= 1:
        return list(residue_indices), [], 1

    local_of_res = {res_i: k for k, res_i in enumerate(residue_indices)}
    atom_idx = _patch_atom_indices(parsed, residue_indices)
    atom_local = np.array(
        [local_of_res[int(parsed.atom_res_index[a])] for a in atom_idx]
    )
    coords = parsed.coords[atom_idx]

    tree = cKDTree(coords)
    pairs = tree.query_pairs(r=contact_cutoff, output_type="ndarray")
    if pairs.size == 0:
        rows = np.empty(0, dtype=np.int64)
        cols = np.empty(0, dtype=np.int64)
    else:
        rows = atom_local[pairs[:, 0]]
        cols = atom_local[pairs[:, 1]]
        keep = rows != cols
        rows, cols = rows[keep], cols[keep]

    n = len(residue_indices)
    data = np.ones(len(rows), dtype=np.int8)
    graph = coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()
    n_comp, labels = connected_components(graph, directed=False)

    seed_label = labels[local_of_res[seed_index]]
    kept = [res_i for res_i in residue_indices if labels[local_of_res[res_i]] == seed_label]
    rejected = [
        res_i for res_i in residue_indices if labels[local_of_res[res_i]] != seed_label
    ]
    return kept, rejected, int(n_comp)


def _seed_normal(
    center: np.ndarray, coords: np.ndarray, neighbour_idx: np.ndarray
) -> np.ndarray | None:
    """Estimate an outward surface normal at a ball centre.

    Taken as the direction from the centroid of the surrounding heavy atoms to
    the ball centre. This is a crude estimate: it is reasonable on a convex
    surface and poor in a concave pocket, which is why it is off by default and
    the contact-graph check is primary.
    """
    if neighbour_idx.size < 3:
        return None
    centroid = coords[neighbour_idx].mean(axis=0)
    normal = center - centroid
    norm = float(np.linalg.norm(normal))
    if norm < 1e-6:
        return None
    return normal / norm


# ---------------------------------------------------------------------------
# Geometry metrics
# ---------------------------------------------------------------------------

def _patch_metrics(
    parsed: ParsedStructure, residue_indices: list[int], radius: float
) -> tuple[float, float, float]:
    """Return ``(spatial_diameter, radius_of_gyration, compactness_ratio)``.

    * ``spatial_diameter`` -- the largest distance between any two heavy atoms
      of the patch's residues, i.e. the patch's true extent. It can exceed
      ``2 * radius``, because membership is decided on a residue's *nearest*
      atom while its other atoms may point outward.
    * ``radius_of_gyration`` -- RMS distance of those heavy atoms from their
      centroid.
    * ``compactness_ratio`` -- ``radius_of_gyration / radius``: dimensionless,
      and reads as how much of the scanning ball the patch actually fills.
    """
    atom_idx = _patch_atom_indices(parsed, residue_indices)
    coords = parsed.coords[atom_idx]
    if len(coords) < 2:
        return 0.0, 0.0, 0.0
    centroid = coords.mean(axis=0)
    rg = float(np.sqrt(((coords - centroid) ** 2).sum(axis=1).mean()))
    # Condensed pairwise distances rather than a full (n, n, 3) broadcast:
    # patch atom counts reach a few hundred, and the broadcast form allocated
    # a large temporary per patch, which profiling showed to dominate the
    # whole ball search.
    diameter = float(pdist(coords).max())
    return diameter, rg, (rg / radius if radius > 0 else 0.0)


def _member_set_hash(parsed: ParsedStructure, residue_indices: list[int]) -> str:
    """Stable hash of a member set, keyed on author identity.

    Author identity rather than internal index, so the same residue set in two
    runs (or in two files of the same entry) hashes identically.
    """
    tokens = sorted(
        f"{parsed.residues[i].rid.auth_chain_id}:"
        f"{parsed.residues[i].rid.auth_seq_id}:"
        f"{parsed.residues[i].rid.ins_code}"
        for i in residue_indices
    )
    return hashlib.sha256("|".join(tokens).encode()).hexdigest()[:16]


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_ball_search(
    parsed: ParsedStructure,
    exposure: ExposureResult,
    config: dict[str, Any],
) -> BallSearchResult:
    """Run the full ball search over every configured radius.

    Deterministic: patches are canonically sorted before IDs are assigned, so
    the output is byte-identical regardless of ``[performance].workers``.
    """
    bs = config["ball_search"]
    conn = config["connectivity"]
    center_mode = bs["center_mode"]
    seed_threshold = float(bs["seed_rsa_threshold"])
    member_threshold = float(bs["member_rsa_threshold"])
    buried_policy = bs["buried_member_policy"]
    disconnected_policy = conn["disconnected_policy"]
    contact_cutoff = float(conn["contact_cutoff_angstrom"])
    use_normal = bool(conn["use_seed_normal_filter"])
    normal_min_cos = float(conn["seed_normal_min_cosine"])
    min_size = int(bs["min_patch_size"])
    max_size = int(bs["max_patch_size"])
    workers = int(config["performance"]["workers"])

    warnings: list[str] = []
    timings: dict[str, float] = {}
    polymer_classes = {RESIDUE_CLASS_STANDARD, RESIDUE_CLASS_NONSTANDARD}

    # --- eligible member universe ---------------------------------------
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

    eligible_atom_mask = eligible_res[parsed.atom_res_index]
    eligible_atom_idx = np.flatnonzero(eligible_atom_mask)
    if eligible_atom_idx.size == 0:
        raise ValueError(
            "no eligible heavy atoms for neighbour search; check chain "
            "selection and [ball_search] residue-class switches"
        )
    eligible_coords = parsed.coords[eligible_atom_idx]
    eligible_res_of_atom = parsed.atom_res_index[eligible_atom_idx]
    tree = cKDTree(eligible_coords)
    timings["build_kdtree_s"] = time.perf_counter() - t0

    # --- B. seeds --------------------------------------------------------
    t0 = time.perf_counter()
    seed_indices: list[int] = []
    n_seed_rsa_missing = 0
    for i, rec in enumerate(parsed.residues):
        if not rec.is_scanned_chain:
            continue
        if (
            rec.residue_class == RESIDUE_CLASS_NONSTANDARD
            and not bs["include_nonstandard_as_seed"]
        ):
            continue
        rsa = exposure.rsa_value_or_none(i)
        if rsa is None:
            n_seed_rsa_missing += 1
            continue
        if rsa >= seed_threshold:
            seed_indices.append(i)

    if n_seed_rsa_missing:
        warnings.append(
            f"{n_seed_rsa_missing} scanned residue(s) could not be tested as "
            "seeds because their RSA is MISSING (no reference maximum for the "
            "residue type); they were skipped as seeds, not treated as buried"
        )

    # --- C. ball centres -------------------------------------------------
    centers: dict[int, BallCenter] = {}
    n_no_center = 0
    for res_i in seed_indices:
        center = compute_ball_center(parsed, res_i, center_mode)
        if center is None:
            n_no_center += 1
            warnings.append(
                f"seed {parsed.residues[res_i].rid.label()} "
                f"({parsed.residues[res_i].rid.res_name}) has no CB, no complete "
                "backbone and no CA; no ball centre could be placed and it was "
                "skipped as a seed (no coordinates were fabricated)"
            )
            continue
        centers[res_i] = center
    seed_indices = [i for i in seed_indices if i in centers]
    timings["select_seeds_and_centers_s"] = time.perf_counter() - t0

    n_virtual = sum(1 for c in centers.values() if c.is_virtual)
    if n_virtual:
        warnings.append(
            f"{n_virtual} of {len(centers)} ball centres are VIRTUAL "
            "(constructed from backbone geometry, mostly glycine); they are "
            "flagged per patch in the output"
        )

    if not seed_indices:
        warnings.append(
            f"no residue reached the seed RSA threshold of {seed_threshold}; "
            "no patches were extracted"
        )
        return BallSearchResult(
            patches=[],
            centers=centers,
            seed_residue_indices=[],
            warnings=warnings,
            timings=timings,
            stats={"n_seeds": 0},
        )

    center_array = np.array([centers[i].position for i in seed_indices])

    # --- D-H per radius --------------------------------------------------
    all_patches: list[Patch] = []
    stats: dict[str, Any] = {
        "n_seeds": len(seed_indices),
        "n_seeds_skipped_no_center": n_no_center,
        "n_virtual_centers": n_virtual,
        "per_radius": {},
    }

    for radius in [float(r) for r in bs["radii_angstrom"]]:
        t_radius = time.perf_counter()
        # One vectorised ball query for every seed. scipy parallelises this
        # internally; results are returned per centre in input order, so the
        # worker count cannot change the outcome.
        neighbour_lists = tree.query_ball_point(
            center_array, r=radius, workers=workers
        )

        raw: dict[frozenset[int], dict[str, Any]] = {}
        n_oversized = 0
        n_undersized = 0
        n_dropped_disc_total = 0
        n_dropped_buried_total = 0

        for seed_pos, seed_i in enumerate(seed_indices):
            hits = np.asarray(neighbour_lists[seed_pos], dtype=np.int64)
            if hits.size == 0:
                continue
            center = center_array[seed_pos]
            dists = np.linalg.norm(eligible_coords[hits] - center, axis=1)
            res_of_hit = eligible_res_of_atom[hits]

            # Nearest hit atom per residue: sort by (residue, distance) and
            # take the first entry of each residue block.
            order = np.lexsort((dists, res_of_hit))
            sorted_res = res_of_hit[order]
            first = np.concatenate(([True], sorted_res[1:] != sorted_res[:-1]))
            sel_positions = order[first]
            sel_res = sorted_res[first]
            sel_atom_global = eligible_atom_idx[hits[sel_positions]]
            sel_dist = dists[sel_positions]

            # Optional seed-normal direction filter.
            if use_normal:
                normal = _seed_normal(center, eligible_coords, hits)
                if normal is not None:
                    vectors = parsed.coords[sel_atom_global] - center
                    norms = np.linalg.norm(vectors, axis=1)
                    safe = norms > 1e-6
                    cosines = np.ones(len(vectors))
                    cosines[safe] = (
                        vectors[safe] @ normal
                    ) / norms[safe]
                    keep_dir = (cosines >= normal_min_cos) | (
                        sel_res == seed_i
                    )
                    sel_res = sel_res[keep_dir]
                    sel_atom_global = sel_atom_global[keep_dir]
                    sel_dist = sel_dist[keep_dir]

            # Exposure rule for members.
            candidates: list[tuple[int, int, float, float | None, str, bool]] = []
            n_dropped_buried = 0
            for res_i, atom_g, dist in zip(sel_res, sel_atom_global, sel_dist):
                res_i = int(res_i)
                rsa = exposure.rsa_value_or_none(res_i)
                missing_reason = ""
                if rsa is None:
                    value = exposure.rsa[res_i]
                    missing_reason = (
                        value.reason if isinstance(value, Missing) else "unknown"
                    )
                    # MISSING RSA is not evidence of burial. Such a residue is
                    # kept on geometry alone and flagged, never silently
                    # dropped by a threshold it cannot be tested against.
                    candidates.append(
                        (res_i, int(atom_g), float(dist), None, missing_reason, False)
                    )
                    continue
                if rsa < member_threshold and res_i != seed_i:
                    if buried_policy == "drop":
                        n_dropped_buried += 1
                        continue
                    candidates.append(
                        (res_i, int(atom_g), float(dist), rsa, "", True)
                    )
                    continue
                candidates.append((res_i, int(atom_g), float(dist), rsa, "", False))

            if not any(c[0] == seed_i for c in candidates):
                # The seed must be in its own patch; if the exposure rule
                # removed it the patch is not interpretable.
                continue

            member_ids = [c[0] for c in candidates]

            # F. connectivity.
            if disconnected_policy == "keep":
                kept_ids, rejected_ids = member_ids, []
                n_comp = 1
            else:
                kept_ids, rejected_ids, n_comp = _connected_component_of_seed(
                    parsed, member_ids, seed_i, contact_cutoff
                )
            if disconnected_policy == "drop":
                final_ids = set(kept_ids)
                flagged_disc: set[int] = set()
                n_dropped_disc = len(rejected_ids)
            else:
                final_ids = set(member_ids)
                flagged_disc = set(rejected_ids)
                n_dropped_disc = 0

            final = [c for c in candidates if c[0] in final_ids]
            if len(final) < min_size:
                n_undersized += 1
                continue
            if len(final) > max_size:
                n_oversized += 1
            n_dropped_disc_total += n_dropped_disc
            n_dropped_buried_total += n_dropped_buried

            key = frozenset(c[0] for c in final)
            entry = raw.get(key)
            if entry is None:
                members = tuple(
                    PatchMember(
                        residue_index=res_i,
                        selecting_atom_name=parsed.atom_name_at(atom_g),
                        selecting_atom_distance=dist,
                        rsa=rsa,
                        rsa_missing_reason=reason,
                        is_seed=(res_i == seed_i),
                        buried_flagged=buried,
                        disconnected_flagged=(res_i in flagged_disc),
                    )
                    for res_i, atom_g, dist, rsa, reason, buried in sorted(
                        final, key=lambda c: c[0]
                    )
                )
                raw[key] = {
                    "members": members,
                    "seeds": [seed_i],
                    "n_dropped_buried": n_dropped_buried,
                    "n_dropped_disconnected": n_dropped_disc,
                    "n_components": n_comp,
                }
            else:
                # G. deduplication: identical residue set, so retain the extra
                # seed rather than discarding this ball. Per-member selecting
                # atoms are kept from the first seed that produced the set;
                # which seed that is becomes deterministic after the canonical
                # sort below, because `seeds` is sorted.
                entry["seeds"].append(seed_i)
                entry["n_components"] = max(entry["n_components"], n_comp)

        # Build Patch objects for this radius.
        radius_patches: list[Patch] = []
        for key, entry in raw.items():
            residue_indices = sorted(key)
            n_seg, seg_lengths, gaps, chains = _segment_structure(
                parsed, residue_indices
            )
            diameter, rg, compactness = _patch_metrics(
                parsed, residue_indices, radius
            )
            seeds = sorted(entry["seeds"])
            # Re-mark is_seed across every seed that produced this set.
            seed_set = set(seeds)
            members = tuple(
                PatchMember(
                    residue_index=m.residue_index,
                    selecting_atom_name=m.selecting_atom_name,
                    selecting_atom_distance=m.selecting_atom_distance,
                    rsa=m.rsa,
                    rsa_missing_reason=m.rsa_missing_reason,
                    is_seed=(m.residue_index in seed_set),
                    buried_flagged=m.buried_flagged,
                    disconnected_flagged=m.disconnected_flagged,
                )
                for m in entry["members"]
            )
            patch = Patch(
                patch_id="",
                radius=radius,
                scale=scale_label(config, radius),
                center_mode=center_mode,
                members=members,
                seed_residue_indices=tuple(seeds),
                seed_center_kinds=tuple(centers[s].kind for s in seeds),
                member_set_hash=_member_set_hash(parsed, residue_indices),
                n_segments=n_seg,
                segment_lengths=seg_lengths,
                author_numbering_gaps=gaps,
                chains=chains,
                spatial_diameter=diameter,
                radius_of_gyration=rg,
                compactness_ratio=compactness,
                n_dropped_buried=int(entry["n_dropped_buried"]),
                n_dropped_disconnected=int(entry["n_dropped_disconnected"]),
                n_components_before_connectivity=int(entry["n_components"]),
            )
            radius_patches.append(patch)

        stats["per_radius"][f"{radius:g}"] = {
            "n_seeds_scanned": len(seed_indices),
            "n_patches_after_dedup": len(radius_patches),
            "n_undersized_dropped": n_undersized,
            "n_oversized_flagged": n_oversized,
            "n_members_dropped_buried": n_dropped_buried_total,
            "n_members_dropped_disconnected": n_dropped_disc_total,
            "seconds": time.perf_counter() - t_radius,
        }
        if n_oversized:
            warnings.append(
                f"radius {radius:g} A: {n_oversized} patch(es) exceed "
                f"max_patch_size={max_size}; they are kept in the patches table "
                "but no FoldDisco query is emitted for them"
            )
        all_patches.extend(radius_patches)
        timings[f"ball_search_r{radius:g}_s"] = (
            stats["per_radius"][f"{radius:g}"]["seconds"]
        )

    # --- canonical ordering, IDs, clustering, contract check -------------
    t0 = time.perf_counter()
    all_patches.sort(
        key=lambda p: (
            p.radius,
            tuple(
                parsed.residues[i].rid.auth_chain_id for i in p.residue_indices()
            ),
            tuple(parsed.residues[i].rid.auth_seq_id for i in p.residue_indices()),
            p.residue_indices(),
            p.seed_residue_indices,
        )
    )
    structure_id = parsed.source.structure_id
    for n, patch in enumerate(all_patches):
        patch.patch_id = (
            f"{structure_id}_m{parsed.model}_r{patch.radius:g}_p{n:05d}"
        )
        _finalise_contract(parsed, patch, max_size)

    if config["dedup"]["cluster_overlapping"]:
        _assign_clusters(all_patches, float(config["dedup"]["jaccard_threshold"]))
    timings["finalise_s"] = time.perf_counter() - t0

    stats["n_patches_total"] = len(all_patches)
    stats["n_patches_discontinuous"] = sum(1 for p in all_patches if p.is_discontinuous)
    stats["n_patches_query_emitted"] = sum(1 for p in all_patches if p.emit_query)

    return BallSearchResult(
        patches=all_patches,
        centers=centers,
        seed_residue_indices=seed_indices,
        warnings=warnings,
        timings=timings,
        stats=stats,
    )


def _finalise_contract(parsed: ParsedStructure, patch: Patch, max_size: int) -> None:
    """Decide whether a patch can be emitted as a FoldDisco query.

    Two independent reasons can block emission, and they are reported
    separately because they mean different things:

    * ``folddisco_expressible = False`` -- a residue identity cannot be written
      in the verified token grammar (insertion code, multi-character chain,
      negative residue number). This is a hard limit of the installed build.
    * ``emit_query = False`` with a size reason -- a configuration limit
      (``min_patch_size`` / ``max_patch_size``), not a limit of FoldDisco.
    """
    residues: list[ResidueId] = [
        parsed.residues[m.residue_index].rid for m in patch.members
    ]
    blockers = []
    for rid in residues:
        rep = check_representable(rid)
        if not rep.ok:
            blockers.append(f"{rid.label()}: {rep.reason}")
    if blockers:
        patch.folddisco_expressible = False
        patch.folddisco_block_reason = "; ".join(blockers[:3])
        patch.emit_query = False
        patch.emit_block_reason = "not expressible in the FoldDisco token grammar"
        patch.query_string = ""
        return

    patch.folddisco_expressible = True
    if patch.size < MIN_QUERY_RESIDUES:
        patch.emit_query = False
        patch.emit_block_reason = (
            f"patch has {patch.size} residues, below the verified FoldDisco "
            f"floor of {MIN_QUERY_RESIDUES}"
        )
        patch.query_string = ""
        return
    if patch.size > max_size:
        patch.emit_query = False
        patch.emit_block_reason = (
            f"patch has {patch.size} residues, above the configured "
            f"max_patch_size={max_size}"
        )
        # The query string is still computed, so a later stage can choose to
        # use it after raising the limit.
        patch.query_string = motif_query_string(residues)
        return

    patch.emit_query = True
    patch.query_string = motif_query_string(residues)


def _assign_clusters(patches: list[Patch], threshold: float) -> None:
    """Single-linkage cluster patches of equal radius by Jaccard overlap.

    Deduplication of *identical* member sets has already happened. This step
    only groups strongly overlapping but distinct patches, so that a later
    stage can treat a family of near-identical balls as one piece of evidence
    instead of many.

    Jaccard index ``|A n B| / |A u B|`` on member residue sets. Single linkage
    is used because the intent is "these describe the same surface region";
    it can chain, which is recorded as a known limitation.
    """
    by_radius: dict[float, list[int]] = {}
    for i, patch in enumerate(patches):
        by_radius.setdefault(patch.radius, []).append(i)

    parent = list(range(len(patches)))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            # Union by smaller index so the representative is deterministic.
            if ra < rb:
                parent[rb] = ra
            else:
                parent[ra] = rb

    for radius, idxs in by_radius.items():
        sets = {i: set(patches[i].residue_indices()) for i in idxs}
        for pos, i in enumerate(idxs):
            for j in idxs[pos + 1 :]:
                a, b = sets[i], sets[j]
                union_size = len(a | b)
                if union_size == 0:
                    continue
                if len(a & b) / union_size >= threshold:
                    union(i, j)

    # Name clusters by their representative patch's canonical position.
    for i, patch in enumerate(patches):
        root = find(i)
        patch.cluster_id = f"cl_{patches[root].radius:g}_{root:05d}"
