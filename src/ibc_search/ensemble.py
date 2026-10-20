"""EXPERIMENTAL: generate a perturbation ensemble from one structure.

**This is a probe, not a feature.** It exists to answer one question before
any ensemble support is designed: *does ensembling change the extracted
motifs and their FoldDisco hits at all, and by how much?* If the answer is
"barely", ensemble support is not worth building; if it is "a lot", the
pipeline's single-structure results carry an uncharacterised variance that
has to be reported.

What this does and does not generate
------------------------------------
This produces an **uncertainty ensemble**, not a conformational one. Each
member displaces every residue as a rigid body by a random vector drawn from
that residue's own deposited B-factors, which under the usual harmonic
reading are an isotropic displacement variance::

    B = 8 * pi^2 * <u^2>   =>   sigma = sqrt(B / (8 * pi^2))

So a member answers "where might these residues actually be, given the
experiment's stated positional uncertainty?" It does **not** answer "what
other conformations does this protein adopt?" There is no backbone dihedral
sampling, no side-chain repacking, no energy function and no dynamics. A
real conformational ensemble needs MD, normal modes or multiple deposited
models, and none of those is approximated here.

Why whole residues rather than individual atoms
-----------------------------------------------
Displacing atoms independently would scramble bond geometry: with a median
B of 68 (sigma ~0.93 A) two bonded atoms would move ~1.3 A relative to each
other, distorting a 1.5 A covalent bond beyond recognition. The resulting
"structure" would not be a plausible position of the molecule, and the
motif differences it produced would be artefacts of broken chemistry.

Rigid per-residue translation keeps every residue's internal geometry exact
and perturbs only where residues sit relative to one another, which is the
quantity the ball search actually consumes (residue centres, residue-residue
distances, per-residue SASA). It is still an approximation -- real
displacements are correlated between neighbouring residues, and these are
drawn independently, so this likely *overstates* the achievable variance
between adjacent residues. Treat the measured spread as a loose upper
bound, not an estimate.

Usage::

    python -m ibc_search.ensemble 8YHU --members 10 --out-dir ensembles/8yhu
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

import biotite.structure.io.pdbx as pdbx_io

from .config import load_config
from .errors import IbcSearchError
from .structure_io import resolve_input

EIGHT_PI_SQUARED = 8.0 * np.pi**2

DEFAULT_MEMBERS = 10
"""Ensemble size. Ten is enough to see whether a spread exists at all."""


@dataclass(frozen=True)
class EnsembleMember:
    index: int
    path: Path
    seed: int
    rmsd_to_reference: float
    """Heavy-atom RMSD to the unperturbed structure, in Angstroms."""


def _bfactor_sigma(b_factors: np.ndarray, floor: float, cap: float) -> np.ndarray:
    """Convert B-factors to isotropic displacement sigmas, in Angstroms.

    Values are clamped: a B-factor of 0 usually means "not refined" rather
    than "perfectly determined", and a very large one would throw a residue
    across the structure.
    """
    with np.errstate(invalid="ignore"):
        sigma = np.sqrt(np.maximum(b_factors, 0.0) / EIGHT_PI_SQUARED)
    sigma = np.nan_to_num(sigma, nan=floor)
    return np.clip(sigma, floor, cap)


def make_ensemble(
    input_spec: str,
    out_dir: Path,
    members: int = DEFAULT_MEMBERS,
    seed: int = 0,
    sigma_floor: float = 0.05,
    sigma_cap: float = 2.0,
    config: dict[str, Any] | None = None,
) -> list[EnsembleMember]:
    """Write ``members`` perturbed copies of ``input_spec`` into ``out_dir``.

    Member 1 is **not** the original: every member is perturbed, so no member
    is privileged. The unperturbed structure stays available as the existing
    single-structure run, which is the baseline the comparison uses.
    """
    if members < 2:
        raise IbcSearchError("an ensemble needs at least 2 members")

    config = config or load_config()
    source = resolve_input(input_spec, config)
    out_dir.mkdir(parents=True, exist_ok=True)

    if source.file_format != "mmcif":
        raise IbcSearchError(
            f"this probe reads mmCIF only; {source.path.name} is "
            f"{source.file_format}. Fetch the entry (mmCIF is the default) "
            "or pass a .cif file."
        )
    # Fetched entries are cached gzipped, so open accordingly.
    if source.gzipped:
        with gzip.open(source.path, "rt", encoding="utf-8") as handle:
            cif = pdbx_io.CIFFile.read(handle)
    else:
        cif = pdbx_io.CIFFile.read(str(source.path))
    reference = pdbx_io.get_structure(cif, model=1, extra_fields=["b_factor"])

    coords = reference.coord.astype(np.float64)
    sigma_atom = _bfactor_sigma(
        np.asarray(reference.b_factor, dtype=np.float64), sigma_floor, sigma_cap
    )

    # Group atoms by residue so each residue moves as one rigid body.
    keys = [
        f"{c}|{r}|{i}"
        for c, r, i in zip(reference.chain_id, reference.res_id, reference.ins_code)
    ]
    _, residue_index = np.unique(np.asarray(keys), return_inverse=True)
    n_residues = int(residue_index.max()) + 1

    # One sigma per residue: the mean of its atoms'.
    sigma_residue = np.zeros(n_residues)
    counts = np.bincount(residue_index, minlength=n_residues)
    np.add.at(sigma_residue, residue_index, sigma_atom)
    sigma_residue /= np.maximum(counts, 1)

    written: list[EnsembleMember] = []
    for member in range(1, members + 1):
        member_seed = seed + member
        rng = np.random.default_rng(member_seed)
        shifts = rng.normal(0.0, 1.0, size=(n_residues, 3)) * sigma_residue[:, None]

        perturbed = reference.copy()
        perturbed.coord = (coords + shifts[residue_index]).astype(np.float32)

        rmsd = float(
            np.sqrt(np.mean(np.sum((perturbed.coord - reference.coord) ** 2, axis=1)))
        )

        path = out_dir / f"member_{member:02d}.cif"
        out = pdbx_io.CIFFile()
        pdbx_io.set_structure(out, perturbed)
        out.write(str(path))

        written.append(
            EnsembleMember(
                index=member, path=path, seed=member_seed, rmsd_to_reference=rmsd
            )
        )

    manifest = {
        "EXPERIMENTAL": (
            "uncertainty ensemble by rigid per-residue displacement from "
            "B-factors; NOT a conformational ensemble"
        ),
        "source": str(source.path),
        "source_sha256": source.sha256,
        "structure_id": source.structure_id,
        "members": members,
        "base_seed": seed,
        "sigma_floor_angstrom": sigma_floor,
        "sigma_cap_angstrom": sigma_cap,
        "residues_perturbed": n_residues,
        "sigma_residue_angstrom": {
            "min": float(sigma_residue.min()),
            "median": float(np.median(sigma_residue)),
            "max": float(sigma_residue.max()),
        },
        "member_rmsd_to_reference_angstrom": {
            m.path.name: round(m.rmsd_to_reference, 4) for m in written
        },
    }
    (out_dir / "ensemble_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    return written


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m ibc_search.ensemble",
        description=(
            "EXPERIMENTAL probe: write a B-factor uncertainty ensemble. This "
            "is not a conformational ensemble and is not a supported feature."
        ),
    )
    parser.add_argument("input", help="structure file or PDB entry ID")
    parser.add_argument("-o", "--out-dir", type=Path, required=True)
    parser.add_argument("-n", "--members", type=int, default=DEFAULT_MEMBERS)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--sigma-cap", type=float, default=2.0)
    args = parser.parse_args(argv)

    try:
        written = make_ensemble(
            args.input,
            args.out_dir,
            members=args.members,
            seed=args.seed,
            sigma_cap=args.sigma_cap,
        )
    except IbcSearchError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    for m in written:
        print(f"{m.path}  seed={m.seed}  rmsd_to_reference={m.rmsd_to_reference:.3f} A")
    return 0


if __name__ == "__main__":
    sys.exit(main())
