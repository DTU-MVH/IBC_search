"""The FoldDisco handoff contract. NO SEARCH IS RUN HERE.

This module encodes *only* what was empirically verified against the FoldDisco
build installed on this machine, plus the round-trip machinery that contract
implies. It deliberately contains no search execution: running queries and
building reference indexes are later stages (see :mod:`ibc_search.ext.search`).

Installed build
---------------
``/usr/local/bin/folddisco``, ``folddisco version`` reports commit
``2a756d9522d5bab8f41bc73ec72ec2327e31cd97``.

Verified facts (each reproduced by a real invocation; see
docs/FOLDDISCO_CONTRACT.md for the exact commands and captured output)
----------------------------------------------------------------------
* A query is passed either inline as ``-q B57,B102,C195`` or as a file via
  ``-q <path>``, where the file holds one
  ``"<structure_path>\\t<comma-separated tokens>"`` line per query. The
  tab-separated file form was confirmed to run without ``-p``.
* A token is ``<chain><resnum>``. The chain field is a **single character**;
  the remainder is parsed as an integer.
* Rejected, with a Rust panic rather than a graceful error:
  - ``B57A``  -> ``Invalid residue: ParseIntError{InvalidDigit}``
    => **insertion codes cannot be expressed**.
  - ``AAA57`` -> ``Invalid residue: ParseIntError{InvalidDigit}``
    => **multi-character chain IDs cannot be expressed**.
  - ``B-5``   -> ``Invalid start residue: ParseIntError{Empty}``
    => **negative residue numbers cannot be expressed**; ``-`` is the
    inclusive range separator (``B57-60`` expands to B57,B58,B59,B60).
* Output columns are selected with ``--format-output``; ``query_residues``
  echoes the query tokens back and ``matching_residues`` carries the hit
  tokens, both in the same ``<chain><resnum>`` form. That echo is the
  round-trip key.
* ``matching_residues`` reports the matched target residues in query order and
  uses ``_`` for *some* unmatched query positions -- but **it is not reliably
  positionally aligned** to ``query_residues``. Measured over 4871 real hit
  rows (545 emitted queries against the 4CHA test index):

  =====================================================  =====
  pattern                                                 rows
  =====================================================  =====
  same length as query, contains ``_``                    2221
  **shorter** than query, contains ``_``                  1577
  same length as query, fully matched                      786
  **shorter** than query, no ``_`` at all                  287
  =====================================================  =====

  Unmatched positions are sometimes padded with ``_`` and sometimes omitted
  outright, including *leading* ones: a query ``B146,C217,C218,C219,C221``
  produced ``C217,C218,C219,C221`` (first position dropped, not padded) and
  also ``B75,B76,_,_`` (two padded, one omitted). Therefore **per-position
  residue identity cannot be recovered from these two columns alone** in this
  build. Use :func:`parse_hit_residues` and treat the matched tokens as a set;
  use :func:`align_hit_to_query` only on rows where the lengths agree, and
  obtain a trustworthy per-residue correspondence from ``--superpose`` /
  ``matching_coordinates`` instead.
* A 30-residue query is accepted and returns in ~0.02 s against a
  single-structure index, so query length is not a tight constraint of the
  tool itself.

Consequences for the data model
-------------------------------
The token form is a **lossy** projection of a residue identity. Therefore:

* Representability is validated up front. A motif containing a residue with an
  insertion code, a multi-character chain ID or a negative residue number is
  marked ``folddisco_expressible = False`` with a reason, and is NOT silently
  mangled into a token FoldDisco would mis-parse or reject.
* Round-tripping a hit uses :class:`ibc_search.residue_id.ResidueIndex`, not
  string parsing alone, and a token that maps to more than one residue (same
  chain and number, different insertion codes) is reported as ambiguous.
* Because tokens are emitted explicitly and never as ranges, the echoed
  ``query_residues`` column is directly comparable to what was sent.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from .errors import ContractError
from .residue_id import ResidueId, ResidueIndex

FOLDDISCO_VERIFIED_COMMIT = "2a756d9522d5bab8f41bc73ec72ec2327e31cd97"
"""The build these rules were verified against. Re-verify if this changes."""

MIN_QUERY_RESIDUES = 2
"""Hard floor on query length.

Verified: a single-residue query (``-q B57``) returns zero result rows from the
installed build, because FoldDisco hashes residue *pairs* -- one node produces
no hash. Two residues return rows. This is a property of the search method, not
a configuration choice, so it is a constant here rather than in the config.
"""

UNMATCHED_TOKEN = "_"
"""Placeholder FoldDisco writes in ``matching_residues`` for an unmatched node.

Verified, but **not written for every unmatched position**: see the module
docstring. Some unmatched positions are omitted from the field entirely.
"""

_TOKEN_RE = re.compile(r"^(?P<chain>.)(?P<resnum>\d+)$")


@dataclass(frozen=True)
class Representability:
    """Whether a residue can be written as a FoldDisco query token."""

    ok: bool
    reason: str = ""


def check_representable(residue: ResidueId) -> Representability:
    """Validate one residue against the verified token grammar."""
    if residue.ins_code not in ("", " "):
        return Representability(
            False,
            f"insertion code {residue.ins_code!r} has no slot in the FoldDisco "
            "token grammar (verified: 'B57A' is rejected with ParseIntError)",
        )
    if len(residue.auth_chain_id) != 1:
        return Representability(
            False,
            f"chain ID {residue.auth_chain_id!r} is not a single character "
            "(verified: 'AAA57' is rejected with ParseIntError)",
        )
    if residue.auth_seq_id < 0:
        return Representability(
            False,
            f"negative residue number {residue.auth_seq_id} cannot be expressed "
            "('-' is the inclusive range separator; verified: 'B-5' is rejected)",
        )
    return Representability(True)


def to_query_token(residue: ResidueId) -> str:
    """Render one residue as a FoldDisco query token.

    Raises :class:`ContractError` rather than emitting a token the installed
    build would mis-parse. Callers that want a soft check should call
    :func:`check_representable` first.
    """
    rep = check_representable(residue)
    if not rep.ok:
        raise ContractError(
            f"residue {residue.label()} is not expressible as a FoldDisco "
            f"query token: {rep.reason}"
        )
    return f"{residue.auth_chain_id}{residue.auth_seq_id}"


def motif_query_string(residues: list[ResidueId]) -> str:
    """Render an ordered motif as a comma-separated query string.

    Tokens are emitted explicitly, never collapsed into ``start-end`` ranges,
    so that the ``query_residues`` column echoed back by FoldDisco is directly
    comparable to what was sent.
    """
    if len(residues) < MIN_QUERY_RESIDUES:
        raise ContractError(
            f"a FoldDisco query needs at least {MIN_QUERY_RESIDUES} residues "
            f"(verified: a 1-residue query returns zero rows); got "
            f"{len(residues)}"
        )
    return ",".join(to_query_token(r) for r in residues)


def parse_query_token(token: str) -> tuple[str, int]:
    """Parse one ``<chain><resnum>`` token into its two fields.

    Mirrors the installed build's grammar: first character is the chain, the
    rest is an integer. Used when reading ``matching_residues`` /
    ``query_residues`` back from a future search stage.
    """
    match = _TOKEN_RE.match(token.strip())
    if match is None:
        raise ContractError(
            f"token {token!r} does not match the verified FoldDisco grammar "
            "'<single-char chain><non-negative integer>'"
        )
    return match.group("chain"), int(match.group("resnum"))


def resolve_query_token(
    token: str, index: ResidueIndex
) -> tuple[ResidueId | None, str]:
    """Map a FoldDisco token back to a residue identity.

    Returns ``(residue, note)``. ``residue`` is ``None`` when the token cannot
    be resolved. ``note`` is ``""`` on a clean resolution, and otherwise states
    why resolution was unclean -- in particular, token collisions created by
    the grammar's missing insertion-code slot are reported rather than being
    resolved by an arbitrary pick.
    """
    chain, seq_id = parse_query_token(token)
    collisions = index.author_key_collisions().get((chain, seq_id))
    if collisions:
        return None, (
            f"token {token!r} is ambiguous: it maps to "
            f"{len(collisions)} residues differing only by insertion code "
            f"({', '.join(r.label() for r in collisions)})"
        )
    residue = index.find_author(chain, seq_id, "")
    if residue is None:
        return None, f"token {token!r} does not correspond to any scanned residue"
    return residue, ""


@dataclass(frozen=True)
class AlignedHitPosition:
    """One query position and what a FoldDisco hit put there.

    ``query_residue`` is the residue of *our* structure that occupied this
    query position; ``hit_token`` is the token FoldDisco reported in the target
    structure, or ``None`` when the position was unmatched (``_``).
    """

    position: int
    query_token: str
    query_residue: ResidueId | None
    hit_token: str | None
    note: str = ""

    @property
    def matched(self) -> bool:
        return self.hit_token is not None


@dataclass(frozen=True)
class HitResidues:
    """A hit's residue fields, parsed without assuming positional alignment.

    This is the safe reading of a hit. ``matched_tokens`` is the set of target
    residues FoldDisco reported, which is well defined; a per-query-position
    correspondence is not (see the module docstring).
    """

    query_tokens: tuple[str, ...]
    matched_tokens: tuple[str, ...]
    """Target residue tokens, with ``_`` placeholders removed."""

    n_unmatched_placeholders: int
    positionally_alignable: bool
    """True when the two fields have equal length, so positions line up.

    Measured on real output: true for about 62% of hit rows. When false, the
    unmatched positions were omitted rather than padded and no per-position
    mapping can be derived from these columns.
    """

    raw_matching_field: str

    @property
    def n_matched(self) -> int:
        return len(self.matched_tokens)


def parse_hit_residues(
    query_residues_field: str, matching_residues_field: str
) -> HitResidues:
    """Parse a hit's two residue fields safely.

    Prefer this over :func:`align_hit_to_query`: it makes no assumption that
    the fields are positionally aligned, because measurement showed they often
    are not.
    """
    query_tokens = tuple(
        t for t in query_residues_field.strip().split(",") if t
    )
    raw_hit = tuple(t for t in matching_residues_field.strip().split(",") if t)
    matched = tuple(t for t in raw_hit if t != UNMATCHED_TOKEN)
    return HitResidues(
        query_tokens=query_tokens,
        matched_tokens=matched,
        n_unmatched_placeholders=sum(1 for t in raw_hit if t == UNMATCHED_TOKEN),
        positionally_alignable=len(raw_hit) == len(query_tokens),
        raw_matching_field=matching_residues_field,
    )


def align_hit_to_query(
    query_residues_field: str,
    matching_residues_field: str,
    index: ResidueIndex,
) -> list[AlignedHitPosition]:
    """Zip a hit's ``matching_residues`` against its ``query_residues``.

    The round-trip step: original structure -> query token -> hit, read back to
    a residue identity of the scanned structure.

    **Only valid when the two fields have equal length.** Raises
    :class:`ContractError` otherwise, because in that case the installed build
    omitted unmatched positions instead of padding them and the positional
    correspondence is genuinely unrecoverable from these columns. Callers that
    must handle every row should use :func:`parse_hit_residues` and set
    semantics, or request ``--superpose``/``matching_coordinates`` for a
    trustworthy per-residue mapping.
    """
    parsed = parse_hit_residues(query_residues_field, matching_residues_field)
    hit_tokens = [t for t in matching_residues_field.strip().split(",") if t]
    if not parsed.positionally_alignable:
        raise ContractError(
            f"query_residues has {len(parsed.query_tokens)} tokens but "
            f"matching_residues has {len(hit_tokens)}; this build omits some "
            "unmatched positions instead of padding them with "
            f"{UNMATCHED_TOKEN!r}, so no positional alignment exists for this "
            "row. Use parse_hit_residues() instead. (Contract verified against "
            f"FoldDisco commit {FOLDDISCO_VERIFIED_COMMIT}.)"
        )
    aligned: list[AlignedHitPosition] = []
    for position, (q_token, h_token) in enumerate(
        zip(parsed.query_tokens, hit_tokens)
    ):
        residue, note = resolve_query_token(q_token, index)
        aligned.append(
            AlignedHitPosition(
                position=position,
                query_token=q_token,
                query_residue=residue,
                hit_token=None if h_token == UNMATCHED_TOKEN else h_token,
                note=note,
            )
        )
    return aligned


def query_file_line(structure_path: str, residues: list[ResidueId]) -> str:
    """Build one line of a FoldDisco ``-q <file>`` query file.

    Format verified against the installed build:
    ``"<structure_path>\\t<comma-separated tokens>"``.
    """
    if "\t" in structure_path or "\n" in structure_path:
        raise ContractError(
            f"structure path {structure_path!r} contains a tab or newline and "
            "cannot be written into a tab-separated FoldDisco query file"
        )
    return f"{structure_path}\t{motif_query_string(residues)}"
