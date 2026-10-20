"""Structure input: local file or PDB entry ID, parsed into stable identities.

Responsibilities
----------------
* Resolve an input specification to either a local file or a cached RCSB
  download, recording URL, UTC fetch date and SHA-256 for the manifest.
* Parse PDB / mmCIF (gzip-aware) into a flat atom array.
* Resolve alternate conformers by a single documented deterministic rule.
* Group atoms into residues **without ever renumbering**, preserving author
  numbering, insertion codes and non-contiguous numbering exactly.
* Classify every residue (standard polymer / non-standard polymer / nucleic /
  carbohydrate / water / ion / other hetero) and attach quality flags for
  missing backbone atoms, missing side-chain atoms, hydrogens, altlocs and
  chain breaks.

Nothing here is silent: every condition the specification lists is either
handled explicitly or recorded as a flag and a warning.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import re
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

import biotite.structure.info as struc_info
import biotite.structure.io.pdb as pdb_io
import biotite.structure.io.pdbx as pdbx_io

from .errors import FetchError, InputError, StructureParseError
from .residue_id import ResidueId, ResidueIndex

PDB_ID_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
"""A legacy PDB entry ID: four characters, first a digit.

Deliberately strict so that an ambiguous four-character *filename* is not
mistaken for an entry ID. Every legacy ID begins with a digit, so a
four-letter word never resolves to a download.
"""

PDB_EXTENDED_ID_RE = re.compile(r"^pdb_[0-9A-Za-z]{8}$", re.IGNORECASE)
"""An extended PDBx identifier: ``pdb_`` followed by eight characters.

The PDB is migrating to this form because the four-character space is close
to exhausted, so both spellings must be accepted. Matching is
case-insensitive; ``pdb_00004hhb`` is the canonical lowercase form.
"""

_LEGACY_IN_EXTENDED_RE = re.compile(r"^pdb_0000([0-9][0-9A-Za-z]{3})$", re.IGNORECASE)
"""Extended IDs that encode a legacy entry: ``pdb_0000XXXX`` -> ``XXXX``."""

WATER_NAMES = frozenset({"HOH", "DOD", "WAT", "H2O"})

BACKBONE_ATOMS = ("N", "CA", "C", "O")

SIDECHAIN_EXCLUDED_BACKBONE_NAMES = frozenset(
    {"N", "CA", "C", "O", "OXT", "OT1", "OT2"}
)
"""Atom names not counted towards a residue's side-chain heavy-atom total.

Terminal oxygens are included here so a chain terminus is not mistaken for an
extra side-chain atom.
"""

SIDECHAIN_HEAVY_ATOM_COUNTS: dict[str, int] = {
    "GLY": 0, "ALA": 1, "SER": 2, "CYS": 2, "THR": 3, "VAL": 3, "PRO": 3,
    "LEU": 4, "ILE": 4, "ASN": 4, "ASP": 4, "MET": 4, "GLN": 5, "GLU": 5,
    "LYS": 5, "HIS": 6, "ARG": 7, "PHE": 7, "TYR": 8, "TRP": 10,
}
"""Expected side-chain heavy-atom count per canonical residue (OXT excluded).

Used only to *flag* an incomplete side chain. A deficit never triggers
fabrication of coordinates.
"""

RESIDUE_CLASS_STANDARD = "polymer_standard_aa"
RESIDUE_CLASS_NONSTANDARD = "polymer_nonstandard_aa"
RESIDUE_CLASS_NUCLEIC = "nucleic"
RESIDUE_CLASS_CARBOHYDRATE = "carbohydrate"
RESIDUE_CLASS_WATER = "water"
RESIDUE_CLASS_ION = "monoatomic_ion"
RESIDUE_CLASS_OTHER = "other_hetero"


@dataclass
class InputSource:
    """Provenance of the structure that was parsed."""

    kind: str
    """``"file"`` or ``"pdb_id"``."""

    path: Path
    """Path actually parsed (the cached download for a PDB ID)."""

    file_format: str
    """``"pdb"`` or ``"mmcif"``."""

    gzipped: bool
    sha256: str
    size_bytes: int
    structure_id: str
    pdb_id: str | None = None
    fetch_url: str | None = None
    fetch_date_utc: str | None = None
    from_cache: bool = False

    def as_dict(self) -> dict[str, Any]:
        return {
            "kind": self.kind,
            "path": str(self.path),
            "file_format": self.file_format,
            "gzipped": self.gzipped,
            "sha256": self.sha256,
            "size_bytes": self.size_bytes,
            "structure_id": self.structure_id,
            "pdb_id": self.pdb_id,
            "fetch_url": self.fetch_url,
            "fetch_date_utc": self.fetch_date_utc,
            "served_from_cache": self.from_cache,
        }


@dataclass
class ResidueRecord:
    """One residue, its identity, classification and quality flags."""

    rid: ResidueId
    residue_class: str
    atom_start: int
    atom_stop: int
    """Half-open slice into the parsed (altloc-resolved) atom arrays."""

    is_scanned_chain: bool
    missing_backbone_atoms: tuple[str, ...] = ()
    has_cb: bool = False
    n_sidechain_heavy_atoms: int = 0
    expected_sidechain_heavy_atoms: int | None = None
    incomplete_sidechain: bool = False
    had_altloc: bool = False
    altloc_chosen: str = ""
    altloc_candidates: tuple[str, ...] = ()
    n_hydrogens_discarded: int = 0
    seq_gap_before: bool = False
    geometric_break_before: bool = False

    @property
    def is_glycine(self) -> bool:
        return self.rid.res_name == "GLY"

    def quality_flags(self) -> list[str]:
        """Compact list of flags for the residues table."""
        flags: list[str] = []
        if self.missing_backbone_atoms:
            flags.append("missing_backbone:" + "/".join(self.missing_backbone_atoms))
        if self.incomplete_sidechain:
            flags.append(
                f"incomplete_sidechain:{self.n_sidechain_heavy_atoms}"
                f"of{self.expected_sidechain_heavy_atoms}"
            )
        if not self.has_cb and not self.is_glycine:
            flags.append("no_cb")
        if self.is_glycine:
            flags.append("glycine")
        if self.had_altloc:
            flags.append(
                "altloc_resolved:" + self.altloc_chosen
                + "(from " + "/".join(self.altloc_candidates) + ")"
            )
        if self.n_hydrogens_discarded:
            flags.append(f"hydrogens_discarded:{self.n_hydrogens_discarded}")
        if self.seq_gap_before:
            flags.append("seq_gap_before")
        if self.geometric_break_before:
            flags.append("chain_break_before")
        if self.residue_class != RESIDUE_CLASS_STANDARD:
            flags.append(f"class:{self.residue_class}")
        return flags


@dataclass
class ParsedStructure:
    """An altloc-resolved, hydrogen-free model with residue bookkeeping."""

    source: InputSource
    model: int
    coords: np.ndarray
    """``(n_atoms, 3)`` float64 coordinates, hydrogens already removed."""

    atom_names: np.ndarray
    elements: np.ndarray
    atom_res_index: np.ndarray
    """``(n_atoms,)`` int: index into :attr:`residues` owning each atom."""

    occupancies: np.ndarray
    residues: list[ResidueRecord]
    residue_index: ResidueIndex
    scanned_chains: list[str]
    available_chains: list[str]
    available_models: list[int]
    warnings: list[str] = field(default_factory=list)
    assembly_note: str = ""

    def residue_atom_slice(self, res_i: int) -> slice:
        rec = self.residues[res_i]
        return slice(rec.atom_start, rec.atom_stop)

    def atom_name_at(self, atom_i: int) -> str:
        return str(self.atom_names[atom_i])


# --------------------------------------------------------------------------
# Input resolution and fetching
# --------------------------------------------------------------------------

def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _detect_format(path: Path) -> tuple[str, bool]:
    """Infer ``(file_format, gzipped)`` from the filename suffixes."""
    name = path.name.lower()
    gzipped = name.endswith(".gz")
    stem = name[:-3] if gzipped else name
    if stem.endswith((".cif", ".mmcif")):
        return "mmcif", gzipped
    if stem.endswith((".pdb", ".ent")):
        return "pdb", gzipped
    raise InputError(
        f"cannot determine structure format from filename {path.name!r}; "
        "expected .pdb, .ent, .cif or .mmcif (optionally .gz)"
    )


@dataclass(frozen=True)
class EntryId:
    """A PDB entry identifier in both spellings it may legitimately have.

    ``legacy`` is the four-character form, or ``None`` for an entry outside
    that space. ``extended`` is the ``pdb_XXXXXXXX`` form, which every entry
    has.

    Verified against RCSB: ``4HHB.cif.gz`` and ``pdb_00004hhb.cif.gz``
    decompress to byte-identical content, so the legacy form is preferred for
    fetching whenever one exists. That keeps one URL, one cache entry and one
    SHA-256 per entry however the user spelled it -- which matters because the
    manifest records all three as provenance.
    """

    legacy: str | None
    extended: str

    @property
    def fetch_id(self) -> str:
        """The identifier used in the download URL."""
        return self.legacy if self.legacy is not None else self.extended

    @property
    def structure_id(self) -> str:
        """Identifier used in patch IDs and output tables.

        The legacy form when the entry has one, so that ``4CHA`` and
        ``pdb_00004cha`` -- the same entry -- yield identical patch IDs.
        """
        return self.fetch_id


def parse_entry_id(spec: str) -> EntryId | None:
    """Parse either accepted entry-ID spelling, or return ``None``.

    Accepts the legacy four-character form and the extended
    ``pdb_XXXXXXXX`` form, the latter case-insensitively.
    """
    token = spec.strip()
    if PDB_ID_RE.match(token):
        legacy = token.upper()
        return EntryId(legacy=legacy, extended=f"pdb_0000{legacy.lower()}")
    if PDB_EXTENDED_ID_RE.match(token):
        match = _LEGACY_IN_EXTENDED_RE.match(token)
        return EntryId(
            legacy=match.group(1).upper() if match else None,
            extended=token.lower(),
        )
    return None


def looks_like_pdb_id(spec: str) -> bool:
    """True for either accepted entry-ID spelling."""
    return parse_entry_id(spec) is not None


def resolve_input(spec: str, config: dict[str, Any]) -> InputSource:
    """Resolve an input spec to a concrete local file plus provenance.

    A spec is treated as a PDB entry ID only when it matches
    :data:`PDB_ID_RE` *and* no file of that name exists, so a local file
    always wins and the choice is never surprising.
    """
    candidate = Path(spec).expanduser()
    if candidate.exists():
        return _source_from_file(candidate)

    entry = parse_entry_id(spec)
    if entry is not None:
        return fetch_pdb_entry(entry, config)

    raise InputError(
        f"input {spec!r} is neither an existing file nor a PDB entry ID "
        "(expected a 4-character ID such as '4CHA', an extended ID such as "
        "'pdb_00004cha', or a local file path)"
    )


def _source_from_file(path: Path) -> InputSource:
    if not path.is_file():
        raise InputError(f"input path is not a regular file: {path}")
    file_format, gzipped = _detect_format(path)
    stem = path.name
    for suffix in (".gz", ".pdb", ".ent", ".cif", ".mmcif"):
        if stem.lower().endswith(suffix):
            stem = stem[: -len(suffix)]
    return InputSource(
        kind="file",
        path=path.resolve(),
        file_format=file_format,
        gzipped=gzipped,
        sha256=sha256_of(path),
        size_bytes=path.stat().st_size,
        structure_id=stem or path.name,
    )


def fetch_pdb_entry(entry: EntryId | str, config: dict[str, Any]) -> InputSource:
    """Fetch a PDB entry from RCSB, caching it on disk.

    Accepts either spelling of an entry ID; a bare string is parsed. mmCIF is
    preferred (see config provenance): it carries both ``label_asym_id`` and
    ``auth_asym_id`` and has no 80-column or single-character-chain limits.

    On any network failure this raises :class:`FetchError`. It never falls
    back to another format or source, because a silent fallback would make the
    manifest's provenance record untrue.
    """
    if isinstance(entry, str):
        parsed = parse_entry_id(entry)
        if parsed is None:
            raise InputError(f"{entry!r} is not a PDB entry ID in either form")
        entry = parsed

    fmt = config["input"]["fetch_format"]
    ext = "cif.gz" if fmt == "mmcif" else "pdb.gz"
    file_format = "mmcif" if fmt == "mmcif" else "pdb"

    # Verified against RCSB: pdb_00004hhb.pdb.gz returns HTTP 404 while the
    # .cif.gz is served. The legacy format cannot express an identifier
    # outside the four-character space -- that exhaustion is precisely why
    # extended IDs exist -- so there is nothing to fall back to, and refusing
    # here beats an opaque 404 from the fetch below.
    if file_format == "pdb" and entry.legacy is None:
        raise InputError(
            f"entry {entry.extended!r} lies outside the four-character ID "
            "space, and RCSB serves no legacy PDB-format file for it, only "
            'mmCIF. Set [input].fetch_format = "mmcif" (the default).'
        )

    pdb_id = entry.fetch_id
    url = f"https://files.rcsb.org/download/{pdb_id}.{ext}"

    cache_dir = Path(config["input"]["cache_dir"]).expanduser()
    cache_dir.mkdir(parents=True, exist_ok=True)
    cached = cache_dir / f"{pdb_id}.{ext}"
    meta_path = cache_dir / f"{pdb_id}.{ext}.fetch.json"

    if cached.is_file() and cached.stat().st_size > 0:
        meta: dict[str, Any] = {}
        if meta_path.is_file():
            try:
                meta = json.loads(meta_path.read_text())
            except json.JSONDecodeError:
                meta = {}
        return InputSource(
            kind="pdb_id",
            path=cached.resolve(),
            file_format=file_format,
            gzipped=True,
            sha256=sha256_of(cached),
            size_bytes=cached.stat().st_size,
            structure_id=entry.structure_id,
            pdb_id=pdb_id,
            fetch_url=meta.get("fetch_url", url),
            fetch_date_utc=meta.get("fetch_date_utc"),
            from_cache=True,
        )

    timeout = float(config["input"]["fetch_timeout_seconds"])
    request = urllib.request.Request(
        url, headers={"User-Agent": "ibc-search/0.1 (academic research)"}
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            payload = response.read()
    except urllib.error.HTTPError as exc:
        raise FetchError(
            f"RCSB returned HTTP {exc.code} for {url}. If {pdb_id!r} is not a "
            "real entry ID, pass a local file path instead."
        ) from exc
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        raise FetchError(
            f"could not fetch {url}: {exc}. The network appears unavailable; "
            "this tool does not fall back to any other source. Pass a local "
            f"file, or pre-populate the cache at {cached}."
        ) from exc

    if not payload:
        raise FetchError(f"RCSB returned an empty response for {url}")

    # Write via a temporary file so an interrupted fetch cannot leave a
    # truncated file that a later run would trust as a cache hit.
    tmp = cached.with_suffix(cached.suffix + ".part")
    tmp.write_bytes(payload)
    tmp.replace(cached)
    fetch_date = datetime.now(timezone.utc).isoformat(timespec="seconds")
    digest = sha256_of(cached)
    meta_path.write_text(
        json.dumps(
            {
                "pdb_id": pdb_id,
                "fetch_url": url,
                "fetch_date_utc": fetch_date,
                "sha256": digest,
                "size_bytes": cached.stat().st_size,
            },
            indent=2,
        )
    )
    return InputSource(
        kind="pdb_id",
        path=cached.resolve(),
        file_format=file_format,
        gzipped=True,
        sha256=digest,
        size_bytes=cached.stat().st_size,
        structure_id=entry.structure_id,
        pdb_id=pdb_id,
        fetch_url=url,
        fetch_date_utc=fetch_date,
        from_cache=False,
    )


# --------------------------------------------------------------------------
# Parsing
# --------------------------------------------------------------------------

def _open_text(source: InputSource) -> io.TextIOBase:
    if source.gzipped:
        return gzip.open(source.path, "rt")
    return source.path.open("rt")


def _read_atom_array(source: InputSource, model: int) -> tuple[Any, list[int], str]:
    """Read one model as an AtomArray with every altloc retained.

    Returns ``(atom_array, available_models, assembly_note)``.
    """
    extra = ["occupancy", "b_factor"]
    with _open_text(source) as fh:
        if source.file_format == "mmcif":
            cif = pdbx_io.CIFFile.read(fh)
            block = cif.block
            model_numbers = _mmcif_model_numbers(block)
            if model not in model_numbers:
                raise StructureParseError(
                    f"model {model} not present in {source.path.name}; "
                    f"available models: {model_numbers}"
                )
            try:
                array = pdbx_io.get_structure(
                    cif,
                    model=model,
                    altloc="all",
                    extra_fields=[*extra, "label_asym_id"],
                    use_author_fields=True,
                )
            except Exception as exc:  # noqa: BLE001 - re-raised with context
                raise StructureParseError(
                    f"biotite failed to parse mmCIF {source.path.name}: {exc}"
                ) from exc
            note = _assembly_note(block)
            return array, model_numbers, note

        pdb_file = pdb_io.PDBFile.read(fh)
        model_count = pdb_file.get_model_count()
        model_numbers = list(range(1, model_count + 1))
        if model not in model_numbers:
            raise StructureParseError(
                f"model {model} not present in {source.path.name}; "
                f"available models: {model_numbers}"
            )
        try:
            array = pdb_file.get_structure(
                model=model, altloc="all", extra_fields=extra
            )
        except Exception as exc:  # noqa: BLE001 - re-raised with context
            raise StructureParseError(
                f"biotite failed to parse PDB {source.path.name}: {exc}"
            ) from exc
        return array, model_numbers, (
            "scanned the deposited asymmetric unit; PDB-format input carries no "
            "machine-readable assembly metadata, so no assembly was read. "
            "Assembly expansion is NOT applied, so a surface buried only in the "
            "biological assembly will read as exposed here."
        )


def _mmcif_model_numbers(block: Any) -> list[int]:
    atom_site = block["atom_site"]
    if "pdbx_PDB_model_num" in atom_site:
        values = np.asarray(atom_site["pdbx_PDB_model_num"].as_array(int))
        return sorted({int(v) for v in values})
    return [1]


def _assembly_note(block: Any) -> str:
    """Describe biological-assembly ambiguity rather than resolving it.

    This component scans the **deposited asymmetric unit**. If the entry
    declares assemblies that differ from it, the surface being scanned may not
    be the biologically relevant one. That ambiguity is recorded, never
    silently resolved.
    """
    try:
        if "pdbx_struct_assembly" in block:
            assembly = block["pdbx_struct_assembly"]
            ids = [str(v) for v in np.atleast_1d(assembly["id"].as_array(str))]
            details = []
            if "oligomeric_details" in assembly:
                details = [
                    str(v)
                    for v in np.atleast_1d(
                        assembly["oligomeric_details"].as_array(str)
                    )
                ]
            return (
                "scanned the deposited asymmetric unit; entry declares "
                f"assembly id(s) {ids} ({'; '.join(details) if details else 'no details'}). "
                "Assembly expansion is NOT applied, so a surface buried only in "
                "the biological assembly will read as exposed here."
            )
    except Exception:  # noqa: BLE001 - metadata is best-effort provenance
        pass
    return (
        "scanned the deposited asymmetric unit; no assembly metadata was read. "
        "Assembly expansion is not applied."
    )


def _residue_group_starts(
    chain_id: np.ndarray,
    res_id: np.ndarray,
    ins_code: np.ndarray,
    res_name: np.ndarray,
) -> np.ndarray:
    """Indices at which a new residue begins, in file order.

    A new residue starts wherever any of chain, residue number, insertion code
    or residue name changes between consecutive atoms. Grouping by *change*
    rather than by sorted numbering is what keeps insertion codes and
    non-contiguous or repeated numbering intact, and is why this package never
    renumbers anything.
    """
    n = len(res_id)
    if n == 0:
        return np.empty(0, dtype=np.int64)
    changed = (
        (chain_id[1:] != chain_id[:-1])
        | (res_id[1:] != res_id[:-1])
        | (ins_code[1:] != ins_code[:-1])
        | (res_name[1:] != res_name[:-1])
    )
    return np.concatenate(([0], np.flatnonzero(changed) + 1)).astype(np.int64)


def _classify_residue(res_name: str, is_hetero: bool, n_heavy_atoms: int) -> str:
    """Assign a residue class from its name, HETATM status and atom count.

    Name is checked before HETATM status on purpose: selenomethionine (MSE)
    and similar modified residues are deposited as HETATM but are genuine
    polymer residues, so trusting the record type alone would wrongly exclude
    them from the scanned polymer.
    """
    name = res_name.strip().upper()
    if name in WATER_NAMES:
        return RESIDUE_CLASS_WATER
    if name in SIDECHAIN_HEAVY_ATOM_COUNTS:
        return RESIDUE_CLASS_STANDARD
    if name in struc_info.amino_acid_names():
        return RESIDUE_CLASS_NONSTANDARD
    if name in struc_info.nucleotide_names():
        return RESIDUE_CLASS_NUCLEIC
    if name in struc_info.carbohydrate_names():
        return RESIDUE_CLASS_CARBOHYDRATE
    if n_heavy_atoms == 1 and is_hetero:
        return RESIDUE_CLASS_ION
    return RESIDUE_CLASS_OTHER


def _resolve_altlocs(
    altloc_id: np.ndarray, occupancy: np.ndarray, group_starts: np.ndarray, n_atoms: int
) -> tuple[np.ndarray, dict[int, tuple[str, tuple[str, ...]]]]:
    """Pick one alternate conformer per residue, deterministically.

    Rule (documented, and the only rule implemented):
      1. Atoms with a blank altloc ID are always kept -- they are shared by
         every conformer.
      2. Among the non-blank altloc IDs present in a residue, choose the one
         with the greatest **summed occupancy** over its atoms.
      3. Ties are broken by the lexicographically smallest altloc ID.

    Summed rather than mean occupancy is used so that a conformer modelled
    with more atoms is not penalised; ties are then broken by ID so the result
    cannot depend on file ordering or on floating-point comparison order.
    """
    keep = np.ones(n_atoms, dtype=bool)
    notes: dict[int, tuple[str, tuple[str, ...]]] = {}
    stops = np.append(group_starts[1:], n_atoms)
    for res_i, (start, stop) in enumerate(zip(group_starts, stops)):
        ids = altloc_id[start:stop]
        blank = (ids == "") | (ids == " ") | (ids == ".") | (ids == "?")
        present = sorted({str(v) for v in ids[~blank]})
        if not present:
            continue
        best_id = ""
        best_sum = -np.inf
        for candidate in present:
            mask = ids == candidate
            total = float(np.nansum(occupancy[start:stop][mask]))
            # `>` keeps the first (lexicographically smallest) candidate on a tie.
            if total > best_sum:
                best_sum = total
                best_id = candidate
        local_keep = blank | (ids == best_id)
        keep[start:stop] = local_keep
        notes[res_i] = (best_id, tuple(present))
    return keep, notes


def parse_structure(
    source: InputSource, config: dict[str, Any]
) -> ParsedStructure:
    """Parse *source* into a :class:`ParsedStructure`.

    Hydrogens are discarded (counted per residue), alternate conformers are
    resolved by the documented rule, and residues are classified and flagged.
    Waters, ions, carbohydrates, nucleic acids and other heteroatoms are
    retained in the arrays -- so they remain available as SASA occluders and
    for inspection -- but are excluded from the scanned polymer by
    classification rather than by deletion.
    """
    model = int(config["input"]["model"])
    array, available_models, assembly_note = _read_atom_array(source, model)
    warnings: list[str] = []

    if array.array_length() == 0:
        raise StructureParseError(
            f"{source.path.name} model {model} contains no atoms"
        )

    chain_id = np.asarray(array.chain_id, dtype=object).astype(str)
    res_id = np.asarray(array.res_id, dtype=np.int64)
    ins_code = np.asarray(array.ins_code, dtype=object).astype(str)
    ins_code = np.array([c.strip() for c in ins_code], dtype=object)
    res_name = np.asarray(array.res_name, dtype=object).astype(str)
    atom_name = np.asarray(array.atom_name, dtype=object).astype(str)
    element = np.asarray(array.element, dtype=object).astype(str)
    hetero = np.asarray(array.hetero, dtype=bool)
    occupancy = np.asarray(
        getattr(array, "occupancy", np.ones(array.array_length())), dtype=np.float64
    )
    altloc_id = np.asarray(
        getattr(array, "altloc_id", np.full(array.array_length(), "")), dtype=object
    ).astype(str)
    label_asym = None
    if "label_asym_id" in array.get_annotation_categories():
        label_asym = np.asarray(array.label_asym_id, dtype=object).astype(str)
    coords = np.asarray(array.coord, dtype=np.float64)

    group_starts = _residue_group_starts(chain_id, res_id, ins_code, res_name)

    # --- alternate conformers -------------------------------------------
    altloc_keep, altloc_notes = _resolve_altlocs(
        altloc_id, occupancy, group_starts, len(res_id)
    )
    n_altloc_dropped = int((~altloc_keep).sum())
    if n_altloc_dropped:
        warnings.append(
            f"resolved alternate conformers: discarded {n_altloc_dropped} atoms "
            f"from non-selected altlocs in {len(altloc_notes)} residues "
            "(rule: highest summed occupancy, ties by smallest altloc ID)"
        )

    # --- hydrogens -------------------------------------------------------
    is_hydrogen = np.isin(np.char.upper(element.astype(str)), ["H", "D"])
    n_hydrogen = int(is_hydrogen.sum())
    if n_hydrogen:
        warnings.append(
            f"discarded {n_hydrogen} hydrogen/deuterium atoms; all geometry and "
            "SASA here is heavy-atom based"
        )

    keep = altloc_keep & ~is_hydrogen

    # Count discarded hydrogens per residue before compaction.
    stops = np.append(group_starts[1:], len(res_id))
    hydrogens_per_group = np.array(
        [int(is_hydrogen[a:b].sum()) for a, b in zip(group_starts, stops)]
    )

    # --- build residue records on the compacted arrays -------------------
    atom_group = np.zeros(len(res_id), dtype=np.int64)
    for res_i, (start, stop) in enumerate(zip(group_starts, stops)):
        atom_group[start:stop] = res_i

    kept_idx = np.flatnonzero(keep)
    c_coords = coords[kept_idx]
    c_atom_name = atom_name[kept_idx]
    c_element = element[kept_idx]
    c_occupancy = occupancy[kept_idx]
    c_group = atom_group[kept_idx]

    groups = _RawResidueGroups(
        starts=group_starts,
        stops=stops,
        chain=[str(chain_id[s]) for s in group_starts],
        residue_class=[
            _classify_residue(
                res_name[s], bool(hetero[s]), int((~is_hydrogen[s:e]).sum())
            )
            for s, e in zip(group_starts, stops)
        ],
        res_id=res_id,
        ins_code=ins_code,
        res_name=res_name,
        label_asym=label_asym,
        altloc_notes=altloc_notes,
        hydrogens_per_group=hydrogens_per_group,
    )

    scanned_chains = _select_scanned_chains(groups, source, model, config, warnings)

    records, residue_ids, atom_res_index = _build_residue_records(
        groups=groups,
        source=source,
        model=model,
        n_kept_atoms=len(kept_idx),
        compacted_group=c_group,
        compacted_atom_name=c_atom_name,
        scanned_chains=scanned_chains,
        warnings=warnings,
    )

    _flag_chain_breaks(records, c_coords, c_atom_name, warnings)
    _summarise_quality(records, warnings)

    return ParsedStructure(
        source=source,
        model=model,
        coords=c_coords,
        atom_names=c_atom_name,
        elements=c_element,
        atom_res_index=atom_res_index,
        occupancies=c_occupancy,
        residues=records,
        residue_index=ResidueIndex(residue_ids),
        scanned_chains=scanned_chains,
        available_chains=sorted(set(groups.chain)),
        available_models=available_models,
        warnings=warnings,
        assembly_note=assembly_note,
    )


@dataclass
class _RawResidueGroups:
    """Per-residue bookkeeping over the *uncompacted* atom arrays.

    Groups are indexed in original file order. Compaction (altloc and hydrogen
    removal) happens afterwards, so these arrays and the compacted ones are
    related only through :attr:`starts`/:attr:`stops` and the group index.
    """

    starts: np.ndarray
    stops: np.ndarray
    chain: list[str]
    residue_class: list[str]
    res_id: np.ndarray
    ins_code: np.ndarray
    res_name: np.ndarray
    label_asym: np.ndarray | None
    altloc_notes: dict[int, tuple[str, tuple[str, ...]]]
    hydrogens_per_group: np.ndarray

    def __len__(self) -> int:
        return len(self.starts)


_POLYMER_CLASSES = frozenset({RESIDUE_CLASS_STANDARD, RESIDUE_CLASS_NONSTANDARD})


def _select_scanned_chains(
    groups: _RawResidueGroups,
    source: InputSource,
    model: int,
    config: dict[str, Any],
    warnings: list[str],
) -> list[str]:
    """Decide which polymer chains are scanned, and warn about the choice."""
    requested = [str(c) for c in config["input"]["chains"]]
    polymer_chains = sorted(
        {
            chain
            for chain, cls in zip(groups.chain, groups.residue_class)
            if cls in _POLYMER_CLASSES
        }
    )

    if requested:
        unknown = [c for c in requested if c not in polymer_chains]
        if unknown:
            raise InputError(
                f"requested chain(s) {unknown} are not polymer chains in "
                f"{source.structure_id} model {model}; available polymer "
                f"chains: {polymer_chains}"
            )
        scanned = [c for c in polymer_chains if c in requested]
    else:
        scanned = list(polymer_chains)

    if not scanned:
        raise StructureParseError(
            f"{source.structure_id} model {model} contains no amino-acid "
            "polymer chain to scan"
        )
    if len(polymer_chains) > 1 and not requested:
        warnings.append(
            f"input has {len(polymer_chains)} polymer chains {polymer_chains} and "
            "no --chain selection was given, so all of them are scanned "
            "together; for a single designed binder this is usually not what "
            "you want, because inter-chain interfaces will read as buried"
        )
    return scanned


def _build_residue_records(
    *,
    groups: _RawResidueGroups,
    source: InputSource,
    model: int,
    n_kept_atoms: int,
    compacted_group: np.ndarray,
    compacted_atom_name: np.ndarray,
    scanned_chains: list[str],
    warnings: list[str],
) -> tuple[list[ResidueRecord], list[ResidueId], np.ndarray]:
    """Build residue records over the compacted atom arrays.

    Assigns the internal ``seq_index`` in file order, independently of author
    numbering, and maps every kept atom to its owning residue. A residue whose
    every atom was removed (all hydrogens, or an unselected altloc) is omitted
    with a warning rather than left as an empty record.
    """
    records: list[ResidueRecord] = []
    residue_ids: list[ResidueId] = []
    atom_res_index = np.full(n_kept_atoms, -1, dtype=np.int64)
    seq_index = 0

    for group_i in range(len(groups)):
        start = int(groups.starts[group_i])
        selected = np.flatnonzero(compacted_group == group_i)
        if selected.size == 0:
            warnings.append(
                f"residue {groups.chain[group_i]}{groups.res_id[start]}"
                f"{groups.ins_code[start]} ({groups.res_name[start]}) retained "
                "no heavy atoms after altloc/hydrogen handling and is omitted"
            )
            continue

        first, last = int(selected[0]), int(selected[-1]) + 1
        if not np.array_equal(selected, np.arange(first, last)):
            raise StructureParseError(
                "internal error: atoms of one residue are not contiguous after "
                "compaction; refusing to continue with ambiguous indexing"
            )

        rid = ResidueId(
            structure_id=source.structure_id,
            model=model,
            auth_chain_id=groups.chain[group_i],
            auth_seq_id=int(groups.res_id[start]),
            ins_code=str(groups.ins_code[start]),
            res_name=str(groups.res_name[start]),
            seq_index=seq_index,
            label_chain_id=(
                str(groups.label_asym[start])
                if groups.label_asym is not None
                else None
            ),
        )

        names = {str(n) for n in compacted_atom_name[first:last]}
        residue_class = groups.residue_class[group_i]
        is_polymer = residue_class in _POLYMER_CLASSES
        expected_sidechain = SIDECHAIN_HEAVY_ATOM_COUNTS.get(rid.res_name)
        n_sidechain = len(
            {n for n in names if n not in SIDECHAIN_EXCLUDED_BACKBONE_NAMES}
        )
        altloc_chosen, altloc_candidates = groups.altloc_notes.get(
            group_i, ("", ())
        )

        records.append(
            ResidueRecord(
                rid=rid,
                residue_class=residue_class,
                atom_start=first,
                atom_stop=last,
                is_scanned_chain=(
                    rid.auth_chain_id in scanned_chains and is_polymer
                ),
                missing_backbone_atoms=(
                    tuple(a for a in BACKBONE_ATOMS if a not in names)
                    if is_polymer
                    else ()
                ),
                has_cb=("CB" in names),
                n_sidechain_heavy_atoms=n_sidechain,
                expected_sidechain_heavy_atoms=expected_sidechain,
                incomplete_sidechain=(
                    expected_sidechain is not None
                    and n_sidechain < expected_sidechain
                ),
                had_altloc=group_i in groups.altloc_notes,
                altloc_chosen=altloc_chosen,
                altloc_candidates=altloc_candidates,
                n_hydrogens_discarded=int(groups.hydrogens_per_group[group_i]),
            )
        )
        residue_ids.append(rid)
        atom_res_index[first:last] = seq_index
        seq_index += 1

    if np.any(atom_res_index < 0):
        raise StructureParseError(
            "internal error: some atoms were not assigned to a residue"
        )
    return records, residue_ids, atom_res_index


def _flag_chain_breaks(
    records: list[ResidueRecord],
    coords: np.ndarray,
    atom_names: np.ndarray,
    warnings: list[str],
) -> None:
    """Flag numbering gaps and geometric chain breaks.

    Both are recorded, because they are different facts: a numbering gap means
    residues are *missing from the model*, while a geometric break means the
    modelled residues are not peptide-bonded. Either can occur without the
    other.
    """
    n_seq_gap = 0
    n_geom_break = 0
    for i in range(1, len(records)):
        prev, cur = records[i - 1], records[i]
        if cur.rid.auth_chain_id != prev.rid.auth_chain_id:
            continue
        if cur.residue_class not in (
            RESIDUE_CLASS_STANDARD, RESIDUE_CLASS_NONSTANDARD
        ) or prev.residue_class not in (
            RESIDUE_CLASS_STANDARD, RESIDUE_CLASS_NONSTANDARD
        ):
            continue
        # Numbering gap: next author number is not previous + 1 and there is no
        # insertion code explaining it.
        if cur.rid.auth_seq_id != prev.rid.auth_seq_id + 1 and not cur.rid.ins_code:
            if cur.rid.auth_seq_id != prev.rid.auth_seq_id:
                cur.seq_gap_before = True
                n_seq_gap += 1
        c_prev = _atom_coord(prev, "C", coords, atom_names)
        n_cur = _atom_coord(cur, "N", coords, atom_names)
        if c_prev is None or n_cur is None:
            cur.geometric_break_before = True
            n_geom_break += 1
            continue
        # A peptide C-N bond is ~1.33 A. 2.0 A is a generous cutoff that still
        # separates bonded neighbours from a genuine break.
        if float(np.linalg.norm(c_prev - n_cur)) > 2.0:
            cur.geometric_break_before = True
            n_geom_break += 1
    if n_seq_gap:
        warnings.append(
            f"{n_seq_gap} residue(s) follow a gap in author numbering "
            "(missing residues in the model); numbering was NOT altered"
        )
    if n_geom_break:
        warnings.append(
            f"{n_geom_break} residue(s) are not peptide-bonded to their "
            "sequence predecessor (chain break or missing backbone atoms); "
            "sequence-segment statistics treat these as breaks"
        )


def _atom_coord(
    rec: ResidueRecord, name: str, coords: np.ndarray, atom_names: np.ndarray
) -> np.ndarray | None:
    sl = slice(rec.atom_start, rec.atom_stop)
    hits = np.flatnonzero(atom_names[sl] == name)
    if hits.size == 0:
        return None
    return coords[rec.atom_start + int(hits[0])]


def _summarise_quality(records: list[ResidueRecord], warnings: list[str]) -> None:
    n_missing_bb = sum(1 for r in records if r.missing_backbone_atoms)
    n_incomplete_sc = sum(1 for r in records if r.incomplete_sidechain)
    n_nonstd = sum(1 for r in records if r.residue_class == RESIDUE_CLASS_NONSTANDARD)
    n_water = sum(1 for r in records if r.residue_class == RESIDUE_CLASS_WATER)
    n_carb = sum(1 for r in records if r.residue_class == RESIDUE_CLASS_CARBOHYDRATE)
    n_other = sum(1 for r in records if r.residue_class == RESIDUE_CLASS_OTHER)
    n_nuc = sum(1 for r in records if r.residue_class == RESIDUE_CLASS_NUCLEIC)
    n_ion = sum(1 for r in records if r.residue_class == RESIDUE_CLASS_ION)
    if n_missing_bb:
        warnings.append(
            f"{n_missing_bb} polymer residue(s) are missing >=1 backbone atom "
            "of N/CA/C/O"
        )
    if n_incomplete_sc:
        warnings.append(
            f"{n_incomplete_sc} residue(s) have fewer side-chain heavy atoms "
            "than their canonical complement; no coordinates were fabricated"
        )
    if n_nonstd:
        warnings.append(
            f"{n_nonstd} non-standard amino-acid residue(s) present; these have "
            "no Tien et al. 2013 reference maximum, so their RSA is reported as "
            "MISSING rather than 0"
        )
    for count, label in (
        (n_water, "water"), (n_carb, "carbohydrate/glycan"),
        (n_other, "other heteroatom (ligand etc.)"), (n_nuc, "nucleic-acid"),
        (n_ion, "monoatomic ion"),
    ):
        if count:
            warnings.append(
                f"{count} {label} residue(s) present; excluded from the scanned "
                "polymer and from seed selection"
            )
