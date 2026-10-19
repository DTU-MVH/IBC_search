"""Per-residue solvent accessibility (SASA) and relative accessibility (RSA).

SASA is computed with biotite's Shrake--Rupley implementation
(:func:`biotite.structure.sasa`), which is an established implementation
rather than a re-derivation here.

Four separate parameters of similar magnitude exist in this pipeline. The one
used in this module is **the solvent probe radius** (default 1.4 A). It is not
the antibody--antigen contact cutoff, not the patch ball radius, and not a
FoldDisco tolerance. See ``ibc_search/data/default.toml``.

RSA normalisation
-----------------
``RSA = SASA / max_acc(residue type)``, with ``max_acc`` taken from
Tien et al. (2013) "Maximum allowed solvent accessibilites of residues in
proteins", PLoS ONE 8(11):e80635 -- the **theoretical** maxima column.
A residue type absent from that table (any non-standard residue) yields
``MISSING``, never ``0.0``: zero is a legitimate RSA for a buried residue, so
using it to mean "unknown" would corrupt downstream exposure statistics.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

import biotite.structure as struc

from .errors import Missing
from .structure_io import (
    RESIDUE_CLASS_NONSTANDARD,
    RESIDUE_CLASS_STANDARD,
    RESIDUE_CLASS_WATER,
    ParsedStructure,
)

TIEN2013_THEORETICAL_MAX_ACC: dict[str, float] = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLU": 223.0, "GLN": 225.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}
"""Theoretical maximum SASA per residue type, in A^2.

Source: Tien MZ, Meyer AG, Sydykova DK, Spielman SJ, Wilke CO (2013),
PLoS ONE 8(11):e80635, Table 1, "Theoretical" column.
"""

RSA_REFERENCE_NAME = "tien2013_theoretical"


@dataclass
class ExposureResult:
    """Per-residue SASA/RSA over the parsed structure, in residue order."""

    sasa: list[float | Missing]
    rsa: list[float | Missing]
    context_description: str
    probe_radius: float
    sphere_points: int
    vdw_radii: str
    rsa_reference: str
    n_context_atoms: int
    warnings: list[str]

    def rsa_value_or_none(self, res_i: int) -> float | None:
        value = self.rsa[res_i]
        return None if isinstance(value, Missing) else float(value)

    def is_exposed(self, res_i: int, threshold: float) -> bool:
        """True when RSA is known and at/above *threshold*.

        A residue with MISSING RSA is **not** treated as exposed. It is also
        not treated as buried: callers that care about the distinction read
        :meth:`rsa_value_or_none` and the residue's flags.
        """
        value = self.rsa_value_or_none(res_i)
        return value is not None and value >= threshold


def _build_context_mask(
    parsed: ParsedStructure, config: dict[str, Any]
) -> tuple[np.ndarray, str, list[str]]:
    """Select the atoms that form the SASA context (the occluding set).

    Waters are always excluded: they are bulk solvent, and including a
    crystallographic water shell would make the protein surface read as
    buried. Ligands/glycans/ions are included only when configured, because a
    designed binder's intrinsic surface is normally assessed without crystal
    additives.
    """
    warnings: list[str] = []
    context_setting = config["exposure"]["context"]
    include_hetero = bool(config["exposure"]["include_hetero_as_occluder"])

    polymer_classes = {RESIDUE_CLASS_STANDARD, RESIDUE_CLASS_NONSTANDARD}
    keep_residue = np.zeros(len(parsed.residues), dtype=bool)
    for i, rec in enumerate(parsed.residues):
        cls = rec.residue_class
        if cls in polymer_classes:
            if context_setting == "scanned_chains":
                keep_residue[i] = rec.rid.auth_chain_id in parsed.scanned_chains
            else:
                keep_residue[i] = True
        elif cls == RESIDUE_CLASS_WATER:
            keep_residue[i] = False
        else:
            keep_residue[i] = include_hetero

    n_hetero = sum(
        1
        for rec in parsed.residues
        if rec.residue_class not in polymer_classes
        and rec.residue_class != RESIDUE_CLASS_WATER
    )
    if n_hetero and not include_hetero:
        warnings.append(
            f"{n_hetero} ligand/glycan/ion residue(s) were excluded from the "
            "SASA context ([exposure].include_hetero_as_occluder = false), so "
            "any surface they cover is reported as exposed"
        )
    if context_setting == "scanned_chains" and len(parsed.available_chains) > len(
        parsed.scanned_chains
    ):
        warnings.append(
            "SASA context is restricted to the scanned chains, so residues at "
            "an interface with a non-scanned chain are reported as exposed"
        )

    atom_mask = keep_residue[parsed.atom_res_index]
    description = (
        f"context={context_setting}; chains={parsed.scanned_chains}; "
        f"waters excluded; hetero_as_occluder={include_hetero}; "
        "hydrogens removed before SASA"
    )
    return atom_mask, description, warnings


def _atom_array_for(parsed: ParsedStructure, mask: np.ndarray) -> struc.AtomArray:
    """Rebuild a biotite AtomArray for the masked atom subset.

    ``res_name`` and ``atom_name`` are populated because the ProtOr radius set
    is keyed on them; ``res_id``/``chain_id``/``ins_code`` are populated so
    biotite's own residue-level helpers stay usable on the result.
    """
    idx = np.flatnonzero(mask)
    array = struc.AtomArray(len(idx))
    array.coord = parsed.coords[idx].astype(np.float32)
    array.element = np.array(
        [str(e).upper() for e in parsed.elements[idx]], dtype="U2"
    )
    array.atom_name = np.array([str(n) for n in parsed.atom_names[idx]], dtype="U6")
    res_i = parsed.atom_res_index[idx]
    array.res_name = np.array(
        [parsed.residues[i].rid.res_name for i in res_i], dtype="U5"
    )
    array.chain_id = np.array(
        [parsed.residues[i].rid.auth_chain_id for i in res_i], dtype="U4"
    )
    array.res_id = np.array(
        [parsed.residues[i].rid.auth_seq_id for i in res_i], dtype=int
    )
    array.ins_code = np.array(
        [parsed.residues[i].rid.ins_code for i in res_i], dtype="U1"
    )
    array.hetero = np.zeros(len(idx), dtype=bool)
    return array


def compute_exposure(
    parsed: ParsedStructure, config: dict[str, Any]
) -> ExposureResult:
    """Compute per-residue SASA and RSA for *parsed*."""
    warnings: list[str] = []
    probe = float(config["exposure"]["probe_radius_angstrom"])
    points = int(config["exposure"]["sphere_points"])
    radii_set = str(config["exposure"]["vdw_radii"])

    mask, description, ctx_warnings = _build_context_mask(parsed, config)
    warnings.extend(ctx_warnings)
    n_context = int(mask.sum())
    if n_context == 0:
        raise ValueError(
            "the SASA context contains no atoms; check [input].chains and "
            "[exposure].context"
        )

    array = _atom_array_for(parsed, mask)
    try:
        atom_sasa = struc.sasa(
            array,
            probe_radius=probe,
            point_number=points,
            vdw_radii=radii_set,
            ignore_ions=False,
        )
    except Exception as exc:  # noqa: BLE001 - fall back with an explicit warning
        if radii_set != "Single":
            warnings.append(
                f"vdw_radii={radii_set!r} failed ({exc}); recomputed with the "
                "'Single' radius set. SASA values are therefore NOT comparable "
                "to a ProtOr-based run"
            )
            radii_set = "Single"
            atom_sasa = struc.sasa(
                array,
                probe_radius=probe,
                point_number=points,
                vdw_radii="Single",
                ignore_ions=False,
            )
        else:
            raise

    # Aggregate atom SASA to residues. biotite returns NaN for atoms it did not
    # evaluate; NaN is propagated as MISSING rather than silently read as 0.
    atom_idx = np.flatnonzero(mask)
    res_of_atom = parsed.atom_res_index[atom_idx]
    n_res = len(parsed.residues)
    totals = np.zeros(n_res, dtype=np.float64)
    nan_counts = np.zeros(n_res, dtype=np.int64)
    counts = np.zeros(n_res, dtype=np.int64)
    finite = np.isfinite(atom_sasa)
    np.add.at(totals, res_of_atom[finite], atom_sasa[finite])
    np.add.at(nan_counts, res_of_atom[~finite], 1)
    np.add.at(counts, res_of_atom, 1)

    sasa: list[float | Missing] = []
    rsa: list[float | Missing] = []
    n_missing_rsa = 0
    for i, rec in enumerate(parsed.residues):
        if counts[i] == 0:
            reason = (
                f"residue is not in the SASA context ({description})"
            )
            sasa.append(Missing(reason))
            rsa.append(Missing(reason))
            continue
        if nan_counts[i] == counts[i]:
            reason = "SASA implementation returned no finite value for any atom"
            sasa.append(Missing(reason))
            rsa.append(Missing(reason))
            continue
        if nan_counts[i]:
            warnings.append(
                f"residue {rec.rid.label()} had {nan_counts[i]}/{counts[i]} "
                "atoms without a finite SASA; the residue total is a lower bound"
            )
        total = float(totals[i])
        sasa.append(total)
        max_acc = TIEN2013_THEORETICAL_MAX_ACC.get(rec.rid.res_name)
        if max_acc is None:
            n_missing_rsa += 1
            rsa.append(
                Missing(
                    f"no Tien et al. 2013 theoretical maximum for residue type "
                    f"{rec.rid.res_name!r}"
                )
            )
        else:
            rsa.append(total / max_acc)

    if n_missing_rsa:
        warnings.append(
            f"{n_missing_rsa} residue(s) have no RSA because their residue type "
            "is absent from the Tien et al. 2013 table; their RSA is MISSING, "
            "not 0, and they cannot be selected as seeds by an RSA threshold"
        )

    return ExposureResult(
        sasa=sasa,
        rsa=rsa,
        context_description=description,
        probe_radius=probe,
        sphere_points=points,
        vdw_radii=radii_set,
        rsa_reference=RSA_REFERENCE_NAME,
        n_context_atoms=n_context,
        warnings=warnings,
    )
