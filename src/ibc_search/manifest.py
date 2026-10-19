"""Run manifest assembly.

The manifest is the record that makes a run reproducible and auditable:
versions, the exact command, input provenance and hashes, the fully resolved
configuration, seeds, warnings and per-stage timings.

It also records the state of the *downstream* FoldDisco dependency -- which
executable was found and what version it reports -- even though this component
runs no search, because the query contract was derived from that build and
would need re-verifying against a different one.
"""

from __future__ import annotations

import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from typing import Any

from . import __version__
from .ball_search_algo import BallSearchResult
from .exposure import ExposureResult
from .folddisco_contract import FOLDDISCO_VERIFIED_COMMIT, MIN_QUERY_RESIDUES
from .structure_io import ParsedStructure


def library_versions() -> dict[str, str]:
    versions: dict[str, str] = {
        "ibc_search": __version__,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "machine": platform.machine(),
    }
    for name in ("numpy", "scipy", "biotite", "pyarrow"):
        try:
            module = __import__(name)
            versions[name] = getattr(module, "__version__", "unknown")
        except ImportError:
            versions[name] = "not installed"
    return versions


def folddisco_status(executable: str) -> dict[str, Any]:
    """Probe the FoldDisco executable without running any search.

    Recorded because the query contract in
    :mod:`ibc_search.folddisco_contract` was verified against one specific
    build; a mismatch here means those rules need re-checking.
    """
    path = shutil.which(executable)
    status: dict[str, Any] = {
        "configured_executable": executable,
        "resolved_path": path,
        "found": path is not None,
        "contract_verified_against_commit": FOLDDISCO_VERIFIED_COMMIT,
        "min_query_residues_verified": MIN_QUERY_RESIDUES,
        "search_executed_by_this_component": False,
    }
    if path is None:
        status["note"] = (
            "FoldDisco executable not found on PATH. This component runs no "
            "search, so the run is still valid, but the emitted query files "
            "could not be checked against a live build."
        )
        return status
    try:
        completed = subprocess.run(  # noqa: S603 - fixed argv, no shell
            [path, "version"],
            capture_output=True,
            text=True,
            timeout=15,
            check=False,
        )
        reported = completed.stdout.strip().splitlines()
        status["reported_version"] = reported[-1].strip() if reported else ""
        status["version_matches_verified_contract"] = (
            status["reported_version"] == FOLDDISCO_VERIFIED_COMMIT
        )
    except (OSError, subprocess.SubprocessError) as exc:
        status["reported_version"] = None
        status["probe_error"] = str(exc)
    return status


def build_manifest(
    *,
    command: list[str],
    config: dict[str, Any],
    config_source: str,
    parsed: ParsedStructure,
    exposure: ExposureResult,
    result: BallSearchResult,
    timings: dict[str, float],
    output_paths: dict[str, Any],
) -> dict[str, Any]:
    """Assemble the full run manifest."""
    all_warnings = (
        [{"stage": "parse", "message": w} for w in parsed.warnings]
        + [{"stage": "exposure", "message": w} for w in exposure.warnings]
        + [{"stage": "ball_search", "message": w} for w in result.warnings]
    )

    collisions = parsed.residue_index.author_key_collisions()

    return {
        "schema_version": 1,
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "command": command,
        "command_line": " ".join(command),
        "versions": library_versions(),
        "input": parsed.source.as_dict(),
        "structure": {
            "structure_id": parsed.source.structure_id,
            "model_scanned": parsed.model,
            "models_available": parsed.available_models,
            "chains_available": parsed.available_chains,
            "chains_scanned": parsed.scanned_chains,
            "n_residues_parsed": len(parsed.residues),
            "n_heavy_atoms_parsed": int(len(parsed.coords)),
            "biological_assembly": parsed.assembly_note,
            "author_key_collisions": {
                f"{chain}{seq}": [r.label() for r in residues]
                for (chain, seq), residues in collisions.items()
            },
            "author_key_collision_note": (
                "residues sharing a (chain, residue number) pair differ only by "
                "insertion code and therefore collide in the FoldDisco token "
                "grammar; such patches are marked not folddisco_expressible"
            ),
        },
        "exposure": {
            "method": "Shrake-Rupley via biotite.structure.sasa",
            "probe_radius_angstrom": exposure.probe_radius,
            "sphere_points": exposure.sphere_points,
            "vdw_radii_set": exposure.vdw_radii,
            "rsa_reference": exposure.rsa_reference,
            "rsa_reference_citation": (
                "Tien MZ, Meyer AG, Sydykova DK, Spielman SJ, Wilke CO (2013) "
                "PLoS ONE 8(11):e80635, theoretical maxima"
            ),
            "context": exposure.context_description,
            "n_context_atoms": exposure.n_context_atoms,
        },
        "ball_search": result.stats,
        "folddisco": folddisco_status(config["folddisco"]["executable"]),
        "resolved_config": config,
        "config_source": config_source,
        "random_seed": config["run"]["random_seed"],
        "timings_seconds": timings | result.timings,
        "outputs": {k: str(v) for k, v in output_paths.items()},
        "warnings": all_warnings,
        "n_warnings": len(all_warnings),
        "interpretation_notice": (
            "Patches in this run are candidate STRUCTURAL MOTIFS: geometric "
            "neighbourhoods of solvent-exposed residues. They are not epitopes, "
            "not validated antibody footprints, and they carry no evidence of "
            "antigenicity, antibody binding or cross-reactivity, pre-existing "
            "anti-drug-antibody reactivity, induction of a new immune response, "
            "or T-cell peptide presentation. No FoldDisco search was run by "
            "this component."
        ),
    }
