"""Shared fixtures: hand-built synthetic structures with known geometry.

Synthetic fixtures are used rather than downloaded entries so that the tests
are offline, fast and -- most importantly -- have *predictable* geometry, so a
test can assert which residues a ball must contain instead of asserting
whatever the code happened to produce.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pytest

ELEMENT_OF = {
    "N": "N", "CA": "C", "C": "C", "O": "O", "CB": "C", "CG": "C",
    "CD": "C", "OG": "O", "SG": "S", "OXT": "O", "NZ": "N", "CE": "C",
}


@dataclass
class Atom:
    name: str
    xyz: tuple[float, float, float]
    altloc: str = " "
    occupancy: float = 1.00


@dataclass
class Residue:
    chain: str
    seq: int
    name: str
    atoms: list[Atom]
    icode: str = " "
    hetero: bool = False


@dataclass
class Model:
    residues: list[Residue] = field(default_factory=list)


def _atom_line(serial: int, atom: Atom, residue: Residue) -> str:
    record = "HETATM" if residue.hetero else "ATOM  "
    name = atom.name
    # PDB atom-name padding: one-or-two-character element names start in
    # column 14, so that "CA" (calcium) and " CA " (alpha carbon) differ.
    name_field = f" {name:<3}" if len(name) < 4 else name
    element = ELEMENT_OF.get(atom.name, atom.name[0])
    x, y, z = atom.xyz
    return (
        f"{record}{serial:>5} {name_field}{atom.altloc}{residue.name:>3}"
        f" {residue.chain}{residue.seq:>4}{residue.icode}   "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{atom.occupancy:>6.2f}{0.0:>6.2f}"
        f"          {element:>2}"
    )


def write_pdb(path: Path, models: list[Model]) -> Path:
    """Write a multi-model PDB file from the given model list."""
    lines: list[str] = []
    serial = 1
    multi = len(models) > 1
    for model_i, model in enumerate(models, start=1):
        if multi:
            lines.append(f"MODEL     {model_i:>4}")
        for residue in model.residues:
            for atom in residue.atoms:
                lines.append(_atom_line(serial, atom, residue))
                serial += 1
        if multi:
            lines.append("ENDMDL")
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")
    return path


def backbone_atoms(origin: np.ndarray, with_cb: bool = True) -> list[Atom]:
    """A minimal residue: N, CA, C, O and optionally CB, in plausible geometry.

    Offsets are small and roughly peptide-like; exact ideality is not needed
    because the tests assert on relationships the algorithm computes, and the
    virtual-CB formula is validated separately against real CB atoms.
    """
    ox, oy, oz = origin
    atoms = [
        Atom("N", (ox - 1.2, oy, oz)),
        Atom("CA", (ox, oy, oz)),
        Atom("C", (ox + 1.2, oy + 0.4, oz)),
        Atom("O", (ox + 1.3, oy + 1.6, oz)),
    ]
    if with_cb:
        atoms.append(Atom("CB", (ox - 0.1, oy - 1.5, oz + 0.2)))
    return atoms


@pytest.fixture
def linear_chain_pdb(tmp_path: Path) -> Path:
    """Ten ALA residues spaced 4 A apart along x, in one chain.

    Spacing is chosen so that a 5 A ball centred on a residue reaches its
    immediate sequence neighbours but not residues two positions away, which
    makes patch membership predictable.
    """
    residues = [
        Residue("A", i + 1, "ALA", backbone_atoms(np.array([4.0 * i, 0.0, 0.0])))
        for i in range(10)
    ]
    return write_pdb(tmp_path / "linear.pdb", [Model(residues)])


@pytest.fixture
def hairpin_pdb(tmp_path: Path) -> Path:
    """A hairpin: residues 1-6 run outward, residues 7-12 run back alongside.

    Residue i and residue 13-i are spatial neighbours but far apart in
    sequence, so a ball anywhere along the hairpin must produce a
    sequence-discontinuous patch. This is the fixture that proves the
    extractor is a 3D neighbourhood search and not a sequence window.
    """
    residues: list[Residue] = []
    for i in range(6):
        residues.append(
            Residue("A", i + 1, "ALA", backbone_atoms(np.array([4.0 * i, 0.0, 0.0])))
        )
    for i in range(6):
        residues.append(
            Residue(
                "A",
                i + 7,
                "ALA",
                backbone_atoms(np.array([4.0 * (5 - i), 5.0, 0.0])),
            )
        )
    return write_pdb(tmp_path / "hairpin.pdb", [Model(residues)])


@pytest.fixture
def icode_pdb(tmp_path: Path) -> Path:
    """A chain containing insertion codes and a numbering gap.

    Residue 10 is followed by 10A and 10B, then numbering jumps to 20. This is
    the fixture for "never silently renumber" and for the FoldDisco token
    collision that insertion codes create.
    """
    residues = [
        Residue("H", 10, "ALA", backbone_atoms(np.array([0.0, 0.0, 0.0]))),
        Residue(
            "H", 10, "SER", backbone_atoms(np.array([3.5, 0.0, 0.0])), icode="A"
        ),
        Residue(
            "H", 10, "CYS", backbone_atoms(np.array([7.0, 0.0, 0.0])), icode="B"
        ),
        Residue("H", 20, "ALA", backbone_atoms(np.array([10.5, 0.0, 0.0]))),
    ]
    return write_pdb(tmp_path / "icode.pdb", [Model(residues)])


@pytest.fixture
def altloc_pdb(tmp_path: Path) -> Path:
    """Two residues with alternate conformers.

    Residue 1 has altlocs A (occupancy 0.3) and B (occupancy 0.7): B must win
    on occupancy. Residue 2 has altlocs A and B both at 0.5: A must win on the
    documented lexicographic tie-break.
    """
    res1 = Residue(
        "A", 1, "SER",
        [
            Atom("N", (-1.2, 0.0, 0.0)),
            Atom("CA", (0.0, 0.0, 0.0)),
            Atom("C", (1.2, 0.4, 0.0)),
            Atom("O", (1.3, 1.6, 0.0)),
            Atom("CB", (0.0, -1.5, 0.0), altloc="A", occupancy=0.30),
            Atom("OG", (0.0, -2.9, 0.0), altloc="A", occupancy=0.30),
            Atom("CB", (0.0, -1.5, 1.0), altloc="B", occupancy=0.70),
            Atom("OG", (0.0, -2.9, 1.0), altloc="B", occupancy=0.70),
        ],
    )
    res2 = Residue(
        "A", 2, "SER",
        [
            Atom("N", (2.8, 0.0, 0.0)),
            Atom("CA", (4.0, 0.0, 0.0)),
            Atom("C", (5.2, 0.4, 0.0)),
            Atom("O", (5.3, 1.6, 0.0)),
            Atom("CB", (4.0, -1.5, 0.0), altloc="A", occupancy=0.50),
            Atom("OG", (4.0, -2.9, 0.0), altloc="A", occupancy=0.50),
            Atom("CB", (4.0, -1.5, 1.0), altloc="B", occupancy=0.50),
            Atom("OG", (4.0, -2.9, 1.0), altloc="B", occupancy=0.50),
        ],
    )
    return write_pdb(tmp_path / "altloc.pdb", [Model([res1, res2])])


@pytest.fixture
def glycine_pdb(tmp_path: Path) -> Path:
    """A chain of ALA with a glycine in the middle and a CB-less LEU.

    Residue 3 is GLY (no CB by definition -> virtual CB is the expected path).
    Residue 5 is LEU with its side chain absent (-> virtual CB, flagged as a
    missing-CB reconstruction, which is a different condition from glycine).
    Residue 7 has only a CA (-> CA fallback).
    """
    residues: list[Residue] = []
    for i in range(8):
        origin = np.array([4.0 * i, 0.0, 0.0])
        if i == 2:
            residues.append(
                Residue("A", i + 1, "GLY", backbone_atoms(origin, with_cb=False))
            )
        elif i == 4:
            residues.append(
                Residue("A", i + 1, "LEU", backbone_atoms(origin, with_cb=False))
            )
        elif i == 6:
            residues.append(
                Residue("A", i + 1, "ALA", [Atom("CA", tuple(origin))])
            )
        else:
            residues.append(Residue("A", i + 1, "ALA", backbone_atoms(origin)))
    return write_pdb(tmp_path / "glycine.pdb", [Model(residues)])


@pytest.fixture
def multi_model_pdb(tmp_path: Path) -> Path:
    """Two models whose coordinates differ by a 10 A shift in z."""
    models = []
    for shift in (0.0, 10.0):
        residues = [
            Residue(
                "A", i + 1, "ALA",
                backbone_atoms(np.array([4.0 * i, 0.0, shift])),
            )
            for i in range(6)
        ]
        models.append(Model(residues))
    return write_pdb(tmp_path / "models.pdb", models)


@pytest.fixture
def two_chain_pdb(tmp_path: Path) -> Path:
    """Two well-separated chains, A and B, 40 A apart."""
    residues: list[Residue] = []
    for i in range(6):
        residues.append(
            Residue("A", i + 1, "ALA", backbone_atoms(np.array([4.0 * i, 0.0, 0.0])))
        )
    for i in range(6):
        residues.append(
            Residue("B", i + 1, "ALA", backbone_atoms(np.array([4.0 * i, 40.0, 0.0])))
        )
    return write_pdb(tmp_path / "two_chain.pdb", [Model(residues)])


@pytest.fixture
def bridged_pdb(tmp_path: Path) -> Path:
    """Two residue clusters close enough to share a ball but not in contact.

    Cluster 1 sits at x=0, cluster 2 at x=9. With a 10 A ball centred in
    cluster 1 both clusters fall inside the ball, but with a 5 A contact
    cutoff they are in different contact components, so the connectivity check
    must separate them. This is the "ball merged two unrelated surfaces" case.
    """
    residues: list[Residue] = []
    for i in range(2):
        residues.append(
            Residue("A", i + 1, "ALA", backbone_atoms(np.array([0.0, 3.0 * i, 0.0])))
        )
    for i in range(2):
        residues.append(
            Residue("A", i + 10, "ALA", backbone_atoms(np.array([9.0, 3.0 * i, 0.0])))
        )
    return write_pdb(tmp_path / "bridged.pdb", [Model(residues)])


@pytest.fixture
def waters_and_ligand_pdb(tmp_path: Path) -> Path:
    """A short chain plus waters, an ion and a sugar, all as HETATM."""
    residues = [
        Residue("A", i + 1, "ALA", backbone_atoms(np.array([4.0 * i, 0.0, 0.0])))
        for i in range(5)
    ]
    residues.append(
        Residue("A", 100, "HOH", [Atom("O", (0.0, 6.0, 0.0))], hetero=True)
    )
    residues.append(
        Residue("A", 101, "HOH", [Atom("O", (4.0, 6.0, 0.0))], hetero=True)
    )
    residues.append(
        Residue("B", 1, "ZN", [Atom("ZN", (8.0, 6.0, 0.0))], hetero=True)
    )
    residues.append(
        Residue(
            "B", 2, "NAG",
            [Atom("C1", (12.0, 6.0, 0.0)), Atom("O5", (13.0, 6.0, 0.0))],
            hetero=True,
        )
    )
    return write_pdb(tmp_path / "hetero.pdb", [Model(residues)])


@pytest.fixture
def base_config(tmp_path: Path) -> dict:
    """Default config, pointed at a temporary output/cache directory."""
    from ibc_search.config import load_config

    config = load_config()
    config["run"]["output_dir"] = str(tmp_path / "out")
    config["input"]["cache_dir"] = str(tmp_path / "cache")
    return config


def rotate_translate(
    coords: np.ndarray, angle: float = 0.7, axis: int = 2,
    shift: tuple[float, float, float] = (13.0, -7.5, 4.25),
) -> np.ndarray:
    """Apply a rigid rotation about *axis* plus a translation."""
    c, s = math.cos(angle), math.sin(angle)
    if axis == 2:
        rot = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    elif axis == 1:
        rot = np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])
    else:
        rot = np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])
    return coords @ rot.T + np.array(shift)


def write_rotated_pdb(source: Path, dest: Path) -> Path:
    """Rewrite a PDB file with every coordinate rigidly transformed."""
    out: list[str] = []
    for line in source.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            xyz = np.array(
                [[float(line[30:38]), float(line[38:46]), float(line[46:54])]]
            )
            new = rotate_translate(xyz)[0]
            out.append(
                line[:30]
                + f"{new[0]:>8.3f}{new[1]:>8.3f}{new[2]:>8.3f}"
                + line[54:]
            )
        else:
            out.append(line)
    dest.write_text("\n".join(out) + "\n")
    return dest
