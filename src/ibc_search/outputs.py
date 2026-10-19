"""Machine-readable output tables and the FoldDisco query file.

Three tables plus one query file are written per run:

``residues.tsv``
    One row per residue in the parsed model: identity, SASA, RSA, exposure
    flag, residue class and quality/missing-atom flags.

``patches.tsv``
    One row per deduplicated patch: identity, seeds, radius, representation
    mode, segment/gap structure, geometry metrics, dedup/cluster membership,
    connectivity counters and the FoldDisco query string.

``patch_members.tsv``
    One row per (patch, member residue): the *selecting atom* and its distance
    to the ball centre, plus per-member flags. Kept separate from
    ``patches.tsv`` so that per-member evidence is queryable rather than
    buried in a packed string.

``folddisco_queries.tsv``
    The handoff file: ``<structure_path>\\t<comma-separated tokens>`` lines,
    the form verified against the installed FoldDisco build. A companion
    ``folddisco_queries.index.tsv`` maps each line back to its patch ID,
    because the query file format itself has no room for one.

Missing values
--------------
A missing quantity is written as an empty cell, with the reason in the
corresponding ``*_missing_reason`` column. It is never written as ``0``.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from .ball_search_algo import BallSearchResult
from .errors import Missing
from .exposure import ExposureResult
from .folddisco_contract import query_file_line
from .structure_io import ParsedStructure

RESIDUE_COLUMNS = [
    "structure_id", "model", "auth_chain_id", "label_chain_id", "auth_seq_id",
    "ins_code", "res_name", "seq_index", "residue_class", "is_scanned_chain",
    "sasa_angstrom2", "sasa_missing_reason", "rsa", "rsa_missing_reason",
    "exposed_at_seed_threshold", "n_atoms_heavy", "missing_backbone_atoms",
    "has_cb", "n_sidechain_heavy_atoms", "expected_sidechain_heavy_atoms",
    "incomplete_sidechain", "had_altloc", "altloc_chosen", "altloc_candidates",
    "n_hydrogens_discarded", "seq_gap_before", "chain_break_before",
    "is_seed", "ball_center_kind", "ball_center_is_virtual",
    "folddisco_expressible", "folddisco_block_reason", "quality_flags",
]

PATCH_COLUMNS = [
    "patch_id", "structure_id", "model", "radius_angstrom", "scale",
    "center_mode", "n_members", "n_segments", "segment_lengths",
    "author_numbering_gaps", "is_discontinuous", "chains",
    "seed_residue_labels", "seed_center_kinds", "n_seeds",
    "spatial_diameter_angstrom", "radius_of_gyration_angstrom",
    "compactness_ratio", "member_set_hash", "cluster_id",
    "n_members_dropped_buried", "n_members_dropped_disconnected",
    "n_components_before_connectivity", "has_buried_flagged_member",
    "has_disconnected_flagged_member", "has_virtual_center",
    "folddisco_expressible", "folddisco_block_reason", "emit_query",
    "emit_block_reason", "query_string", "member_labels",
]

MEMBER_COLUMNS = [
    "patch_id", "structure_id", "model", "radius_angstrom", "residue_label",
    "auth_chain_id", "auth_seq_id", "ins_code", "res_name", "seq_index",
    "selecting_atom_name", "selecting_atom_distance_angstrom", "rsa",
    "rsa_missing_reason", "is_seed", "buried_flagged", "disconnected_flagged",
    "folddisco_query_token",
]


def _fmt(value: Any) -> str:
    """Render one cell. ``None``/:class:`Missing` become an empty cell."""
    if value is None or isinstance(value, Missing):
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return f"{value:.4f}"
    if isinstance(value, (list, tuple)):
        return ";".join(str(v) for v in value)
    return str(value)


def _write_tsv(path: Path, columns: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for row in rows:
            writer.writerow([_fmt(row.get(col)) for col in columns])


def _write_parquet(path: Path, columns: list[str], rows: list[dict[str, Any]]) -> None:
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise RuntimeError(
            "table_format='parquet' requires pyarrow; install the 'parquet' "
            "extra or set [run].table_format = 'tsv'"
        ) from exc
    table = pa.table(
        {col: [_fmt(row.get(col)) or None for row in rows] for col in columns}
    )
    pq.write_table(table, path)


def residue_rows(
    parsed: ParsedStructure,
    exposure: ExposureResult,
    result: BallSearchResult,
    config: dict[str, Any],
) -> list[dict[str, Any]]:
    """Build the residues table."""
    from .folddisco_contract import check_representable

    threshold = float(config["ball_search"]["seed_rsa_threshold"])
    seeds = set(result.seed_residue_indices)
    rows: list[dict[str, Any]] = []
    for i, rec in enumerate(parsed.residues):
        sasa_value = exposure.sasa[i]
        rsa_value = exposure.rsa[i]
        center = result.centers.get(i)
        rep = check_representable(rec.rid)
        row: dict[str, Any] = dict(rec.rid.as_dict())
        row.update(
            {
                "residue_class": rec.residue_class,
                "is_scanned_chain": rec.is_scanned_chain,
                "sasa_angstrom2": None if isinstance(sasa_value, Missing) else sasa_value,
                "sasa_missing_reason": (
                    sasa_value.reason if isinstance(sasa_value, Missing) else ""
                ),
                "rsa": None if isinstance(rsa_value, Missing) else rsa_value,
                "rsa_missing_reason": (
                    rsa_value.reason if isinstance(rsa_value, Missing) else ""
                ),
                # A residue with MISSING RSA is reported as not-exposed here but
                # the reason column says why; it is never counted as buried in
                # the statistics.
                "exposed_at_seed_threshold": exposure.is_exposed(i, threshold),
                "n_atoms_heavy": rec.atom_stop - rec.atom_start,
                "missing_backbone_atoms": rec.missing_backbone_atoms,
                "has_cb": rec.has_cb,
                "n_sidechain_heavy_atoms": rec.n_sidechain_heavy_atoms,
                "expected_sidechain_heavy_atoms": rec.expected_sidechain_heavy_atoms,
                "incomplete_sidechain": rec.incomplete_sidechain,
                "had_altloc": rec.had_altloc,
                "altloc_chosen": rec.altloc_chosen,
                "altloc_candidates": rec.altloc_candidates,
                "n_hydrogens_discarded": rec.n_hydrogens_discarded,
                "seq_gap_before": rec.seq_gap_before,
                "chain_break_before": rec.geometric_break_before,
                "is_seed": i in seeds,
                "ball_center_kind": center.kind if center else "",
                "ball_center_is_virtual": center.is_virtual if center else None,
                "folddisco_expressible": rep.ok,
                "folddisco_block_reason": rep.reason,
                "quality_flags": rec.quality_flags(),
            }
        )
        rows.append(row)
    return rows


def patch_rows(
    parsed: ParsedStructure, result: BallSearchResult
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Build the patches table and the patch-members table."""
    patch_out: list[dict[str, Any]] = []
    member_out: list[dict[str, Any]] = []
    for patch in result.patches:
        member_labels = [
            parsed.residues[m.residue_index].rid.label() for m in patch.members
        ]
        seed_labels = [
            parsed.residues[i].rid.label() for i in patch.seed_residue_indices
        ]
        patch_out.append(
            {
                "patch_id": patch.patch_id,
                "structure_id": parsed.source.structure_id,
                "model": parsed.model,
                "radius_angstrom": patch.radius,
                "scale": patch.scale,
                "center_mode": patch.center_mode,
                "n_members": patch.size,
                "n_segments": patch.n_segments,
                "segment_lengths": patch.segment_lengths,
                "author_numbering_gaps": patch.author_numbering_gaps,
                "is_discontinuous": patch.is_discontinuous,
                "chains": patch.chains,
                "seed_residue_labels": seed_labels,
                "seed_center_kinds": patch.seed_center_kinds,
                "n_seeds": len(patch.seed_residue_indices),
                "spatial_diameter_angstrom": patch.spatial_diameter,
                "radius_of_gyration_angstrom": patch.radius_of_gyration,
                "compactness_ratio": patch.compactness_ratio,
                "member_set_hash": patch.member_set_hash,
                "cluster_id": patch.cluster_id,
                "n_members_dropped_buried": patch.n_dropped_buried,
                "n_members_dropped_disconnected": patch.n_dropped_disconnected,
                "n_components_before_connectivity": (
                    patch.n_components_before_connectivity
                ),
                "has_buried_flagged_member": any(
                    m.buried_flagged for m in patch.members
                ),
                "has_disconnected_flagged_member": any(
                    m.disconnected_flagged for m in patch.members
                ),
                "has_virtual_center": any(
                    result.centers[s].is_virtual for s in patch.seed_residue_indices
                ),
                "folddisco_expressible": patch.folddisco_expressible,
                "folddisco_block_reason": patch.folddisco_block_reason,
                "emit_query": patch.emit_query,
                "emit_block_reason": patch.emit_block_reason,
                "query_string": patch.query_string,
                "member_labels": member_labels,
            }
        )
        for member in patch.members:
            rid = parsed.residues[member.residue_index].rid
            from .folddisco_contract import check_representable

            rep = check_representable(rid)
            member_out.append(
                {
                    "patch_id": patch.patch_id,
                    "structure_id": parsed.source.structure_id,
                    "model": parsed.model,
                    "radius_angstrom": patch.radius,
                    "residue_label": rid.label(),
                    "auth_chain_id": rid.auth_chain_id,
                    "auth_seq_id": rid.auth_seq_id,
                    "ins_code": rid.ins_code,
                    "res_name": rid.res_name,
                    "seq_index": rid.seq_index,
                    "selecting_atom_name": member.selecting_atom_name,
                    "selecting_atom_distance_angstrom": (
                        member.selecting_atom_distance
                    ),
                    "rsa": member.rsa,
                    "rsa_missing_reason": member.rsa_missing_reason,
                    "is_seed": member.is_seed,
                    "buried_flagged": member.buried_flagged,
                    "disconnected_flagged": member.disconnected_flagged,
                    "folddisco_query_token": (
                        f"{rid.auth_chain_id}{rid.auth_seq_id}" if rep.ok else ""
                    ),
                }
            )
    return patch_out, member_out


def write_tables(
    out_dir: Path,
    parsed: ParsedStructure,
    exposure: ExposureResult,
    result: BallSearchResult,
    config: dict[str, Any],
) -> dict[str, Path]:
    """Write all tables and the FoldDisco query file."""
    out_dir.mkdir(parents=True, exist_ok=True)
    fmt = config["run"]["table_format"]
    writer = _write_parquet if fmt == "parquet" else _write_tsv
    ext = "parquet" if fmt == "parquet" else "tsv"

    res_rows = residue_rows(parsed, exposure, result, config)
    pat_rows, mem_rows = patch_rows(parsed, result)

    paths: dict[str, Path] = {}
    paths["residues"] = out_dir / f"residues.{ext}"
    writer(paths["residues"], RESIDUE_COLUMNS, res_rows)
    paths["patches"] = out_dir / f"patches.{ext}"
    writer(paths["patches"], PATCH_COLUMNS, pat_rows)
    paths["patch_members"] = out_dir / f"patch_members.{ext}"
    writer(paths["patch_members"], MEMBER_COLUMNS, mem_rows)

    paths.update(write_query_file(out_dir, parsed, result))
    return paths


def write_query_file(
    out_dir: Path, parsed: ParsedStructure, result: BallSearchResult
) -> dict[str, Path]:
    """Write the FoldDisco query file plus its patch-ID index.

    The query file holds only lines FoldDisco can actually consume. The index
    file is what preserves the patch identity that the query format cannot
    carry, and it is what a later search stage joins hits back on.
    """
    query_path = out_dir / "folddisco_queries.tsv"
    index_path = out_dir / "folddisco_queries.index.tsv"
    structure_path = str(parsed.source.path)

    lines: list[str] = []
    index_rows: list[dict[str, Any]] = []
    for patch in result.patches:
        if not patch.emit_query:
            continue
        residues = [parsed.residues[m.residue_index].rid for m in patch.members]
        lines.append(query_file_line(structure_path, residues))
        index_rows.append(
            {
                "line_number": len(lines),
                "patch_id": patch.patch_id,
                "structure_path": structure_path,
                "radius_angstrom": patch.radius,
                "scale": patch.scale,
                "n_members": patch.size,
                "query_string": patch.query_string,
            }
        )

    query_path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
    _write_tsv(
        index_path,
        [
            "line_number", "patch_id", "structure_path", "radius_angstrom",
            "scale", "n_members", "query_string",
        ],
        index_rows,
    )
    return {"folddisco_queries": query_path, "folddisco_query_index": index_path}


def write_manifest(out_dir: Path, manifest: dict[str, Any]) -> Path:
    """Write the run manifest as JSON."""
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "run_manifest.json"
    path.write_text(json.dumps(manifest, indent=2, sort_keys=False), encoding="utf-8")
    return path
