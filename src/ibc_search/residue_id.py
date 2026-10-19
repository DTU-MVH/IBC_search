"""Stable residue identity and the three-way numbering map.

Three distinct residue namespaces exist in this pipeline and must never be
confused:

1. **Original structure numbering** -- what the depositor wrote:
   author chain ID, author residue number, insertion code. Not necessarily
   contiguous, not necessarily starting at 1, may be negative, may repeat with
   different insertion codes. This package NEVER renumbers it.

2. **Internal index** (``seq_index``) -- a 0-based position into this run's
   residue arrays, assigned in file order over the selected model. Used for
   NumPy addressing and for sequence-adjacency reasoning only.

3. **FoldDisco query token** -- ``"<chain><resnum>"``, e.g. ``"B57"``. A
   lossy projection of (1): the installed FoldDisco build's token grammar has
   no slot for an insertion code, cannot hold a multi-character chain ID and
   cannot hold a negative residue number. See
   :mod:`ibc_search.folddisco_contract`.

Because (3) is lossy, round-tripping a FoldDisco hit back to a residue needs an
explicit map rather than string parsing alone, and a motif containing a residue
that is not representable in (3) must be reported as such instead of being
quietly mangled.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True, order=True)
class ResidueId:
    """Immutable identity of one residue in one model of one structure.

    Hashable and totally ordered, so it can key dictionaries and sort
    deterministically. Ordering is by ``(structure_id, model, auth_chain_id,
    auth_seq_id, ins_code, seq_index)`` through the field declaration order
    below, which puts ``seq_index`` last so it never overrides the
    depositor-facing keys.
    """

    structure_id: str
    """Identifier of the source structure (PDB ID, or the input file stem)."""

    model: int
    """Deposited model number that was scanned (1-based)."""

    auth_chain_id: str
    """Author chain ID (mmCIF ``auth_asym_id`` / PDB column 22)."""

    auth_seq_id: int
    """Author residue number. May be negative or non-contiguous."""

    ins_code: str
    """Insertion code, or ``""`` when absent. Never ``None``."""

    res_name: str
    """Three-letter (or longer, for some mmCIF components) residue name."""

    seq_index: int
    """0-based internal array index for this run."""

    label_chain_id: str | None = field(default=None, compare=False)
    """mmCIF ``label_asym_id`` when the source provided one, else ``None``.

    Excluded from comparison/hashing: it is provenance, and a PDB-format input
    has none, so including it would make the same residue compare unequal
    across formats.
    """

    def label(self) -> str:
        """Human-readable label, e.g. ``"B57"`` or ``"B57A"`` with an icode.

        This is for logs and tables. It is NOT the FoldDisco query token; use
        :func:`ibc_search.folddisco_contract.to_query_token` for that, which
        validates representability instead of silently producing a string
        FoldDisco would reject.
        """
        return f"{self.auth_chain_id}{self.auth_seq_id}{self.ins_code}"

    def as_dict(self) -> dict[str, object]:
        """Flat mapping for table emission."""
        return {
            "structure_id": self.structure_id,
            "model": self.model,
            "auth_chain_id": self.auth_chain_id,
            "label_chain_id": self.label_chain_id if self.label_chain_id is not None else "",
            "auth_seq_id": self.auth_seq_id,
            "ins_code": self.ins_code,
            "res_name": self.res_name,
            "seq_index": self.seq_index,
        }


class ResidueIndex:
    """Bidirectional map between residue identities and internal indices.

    Holds the authoritative original -> internal -> query mapping for a run.
    Built once by the parser and then treated as read-only.
    """

    def __init__(self, residues: list[ResidueId]) -> None:
        self._residues = list(residues)
        self._by_seq_index: dict[int, ResidueId] = {r.seq_index: r for r in self._residues}
        # Keyed on depositor-facing identity, which is what a FoldDisco hit or
        # an external annotation will refer to.
        self._by_author_key: dict[tuple[str, int, str], ResidueId] = {
            (r.auth_chain_id, r.auth_seq_id, r.ins_code): r for r in self._residues
        }
        if len(self._by_author_key) != len(self._residues):
            # Duplicate (chain, resnum, icode) within one model means the
            # structure is ambiguous; surfacing it is better than guessing.
            seen: dict[tuple[str, int, str], int] = {}
            dupes = []
            for r in self._residues:
                key = (r.auth_chain_id, r.auth_seq_id, r.ins_code)
                if key in seen:
                    dupes.append(key)
                seen[key] = r.seq_index
            raise ValueError(
                f"duplicate author residue keys within one model: {dupes[:5]}"
            )

    def __len__(self) -> int:
        return len(self._residues)

    def __iter__(self):
        return iter(self._residues)

    @property
    def residues(self) -> list[ResidueId]:
        """Residues in internal index order."""
        return list(self._residues)

    def by_seq_index(self, seq_index: int) -> ResidueId:
        """Internal index -> identity."""
        return self._by_seq_index[seq_index]

    def by_author(self, chain: str, seq_id: int, ins_code: str = "") -> ResidueId:
        """Original structure numbering -> identity.

        This is the entry point for reading a FoldDisco hit back: a hit token
        yields ``(chain, seq_id)``, and ``ins_code`` defaults to ``""``
        because the token grammar cannot carry one.
        """
        return self._by_author_key[(chain, seq_id, ins_code)]

    def find_author(
        self, chain: str, seq_id: int, ins_code: str = ""
    ) -> ResidueId | None:
        """Like :meth:`by_author` but returns ``None`` when absent."""
        return self._by_author_key.get((chain, seq_id, ins_code))

    def author_key_collisions(self) -> dict[tuple[str, int], list[ResidueId]]:
        """Residues sharing a (chain, resnum) pair but differing by icode.

        These are exactly the residues whose FoldDisco query tokens would
        collide, since the token grammar drops the insertion code. Reported in
        the run manifest so the ambiguity is never silent.
        """
        groups: dict[tuple[str, int], list[ResidueId]] = {}
        for r in self._residues:
            groups.setdefault((r.auth_chain_id, r.auth_seq_id), []).append(r)
        return {k: v for k, v in groups.items() if len(v) > 1}
