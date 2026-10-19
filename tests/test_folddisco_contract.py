"""The FoldDisco handoff contract.

These tests encode behaviour that was verified by invoking the installed
FoldDisco build (see docs/FOLDDISCO_CONTRACT.md). They assert that this
package refuses to emit anything that build would reject, and that a hit can
be read back to a residue identity.
"""

from __future__ import annotations

import pytest

from ibc_search.errors import ContractError
from ibc_search.folddisco_contract import (
    MIN_QUERY_RESIDUES,
    UNMATCHED_TOKEN,
    align_hit_to_query,
    parse_hit_residues,
    check_representable,
    motif_query_string,
    parse_query_token,
    query_file_line,
    resolve_query_token,
    to_query_token,
)
from ibc_search.residue_id import ResidueId, ResidueIndex


def rid(chain="B", seq=57, icode="", name="ALA", index=0, structure="X"):
    return ResidueId(
        structure_id=structure, model=1, auth_chain_id=chain, auth_seq_id=seq,
        ins_code=icode, res_name=name, seq_index=index,
    )


# --- token grammar ----------------------------------------------------------

def test_plain_residue_becomes_a_token():
    assert to_query_token(rid("B", 57)) == "B57"
    assert to_query_token(rid("C", 195)) == "C195"


def test_motif_query_string_is_explicit_never_a_range():
    residues = [rid("B", 57, index=0), rid("B", 58, index=1), rid("B", 59, index=2)]
    # Consecutive residues could be written "B57-59", but the explicit form is
    # emitted so that the echoed query_residues column is directly comparable.
    assert motif_query_string(residues) == "B57,B58,B59"


def test_insertion_code_is_rejected_not_mangled():
    residue = rid("B", 57, icode="A")
    rep = check_representable(residue)
    assert rep.ok is False
    assert "insertion code" in rep.reason
    # Verified against the installed build: "B57A" panics with a ParseIntError.
    with pytest.raises(ContractError, match="insertion code"):
        to_query_token(residue)


def test_multi_character_chain_is_rejected():
    rep = check_representable(rid("AAA", 57))
    assert rep.ok is False
    assert "single character" in rep.reason
    with pytest.raises(ContractError):
        to_query_token(rid("AAA", 57))


def test_negative_residue_number_is_rejected():
    rep = check_representable(rid("B", -5))
    assert rep.ok is False
    assert "range separator" in rep.reason
    with pytest.raises(ContractError):
        to_query_token(rid("B", -5))


def test_query_below_verified_floor_is_rejected():
    assert MIN_QUERY_RESIDUES == 2
    # Verified: a single-residue query returns zero rows, because FoldDisco
    # hashes residue pairs.
    with pytest.raises(ContractError, match="at least 2 residues"):
        motif_query_string([rid("B", 57)])


def test_parse_query_token_matches_the_verified_grammar():
    assert parse_query_token("B57") == ("B", 57)
    assert parse_query_token(" C195 ") == ("C", 195)

    # A numeric chain ID is legal: the installed build takes the first
    # character as the chain and parses the remainder as an integer, so "57"
    # means chain "5", residue 7. Mirroring that exactly matters more than
    # looking tidy -- guessing "residue 57, chain unknown" would silently
    # disagree with the tool.
    assert parse_query_token("57") == ("5", 7)

    for bad in ("B57A", "AAA57", "B-5", "", "BB"):
        with pytest.raises(ContractError):
            parse_query_token(bad)


def test_query_file_line_uses_the_verified_tab_format():
    residues = [rid("B", 57, index=0), rid("C", 195, index=1)]
    line = query_file_line("/data/model.pdb", residues)
    assert line == "/data/model.pdb\tB57,C195"
    assert line.count("\t") == 1


def test_query_file_line_refuses_a_path_containing_a_tab():
    residues = [rid("B", 57, index=0), rid("C", 195, index=1)]
    with pytest.raises(ContractError, match="tab or newline"):
        query_file_line("/data/we\tird.pdb", residues)


# --- round-trip -------------------------------------------------------------

def test_token_resolves_back_to_a_residue():
    index = ResidueIndex([rid("B", 57, index=0), rid("C", 195, index=1)])
    residue, note = resolve_query_token("B57", index)
    assert residue is not None and residue.res_name == "ALA"
    assert residue.seq_index == 0
    assert note == ""


def test_unknown_token_resolves_to_none_with_a_reason():
    index = ResidueIndex([rid("B", 57, index=0)])
    residue, note = resolve_query_token("Z999", index)
    assert residue is None
    assert "does not correspond" in note


def test_insertion_code_collision_is_reported_not_guessed():
    """Two residues differing only by insertion code share one token.

    The grammar cannot distinguish them, so resolution must report the
    ambiguity rather than silently returning one of them.
    """
    index = ResidueIndex(
        [
            rid("H", 10, icode="", index=0),
            rid("H", 10, icode="A", name="SER", index=1),
        ]
    )
    residue, note = resolve_query_token("H10", index)
    assert residue is None
    assert "ambiguous" in note
    assert "insertion code" in note


def test_align_hit_to_query_handles_unmatched_placeholder():
    """``matching_residues`` pads unmatched positions with ``_``.

    Verified with a 30-residue query whose hit row read ``_,_,_,_,G155,...``.
    """
    index = ResidueIndex(
        [
            rid("B", 57, index=0),
            rid("B", 102, name="SER", index=1),
            rid("C", 195, name="CYS", index=2),
        ]
    )
    aligned = align_hit_to_query("B57,B102,C195", f"F57,{UNMATCHED_TOKEN},G195", index)

    assert [a.position for a in aligned] == [0, 1, 2]
    assert [a.matched for a in aligned] == [True, False, True]
    assert aligned[0].hit_token == "F57"
    assert aligned[1].hit_token is None
    assert aligned[2].hit_token == "G195"
    # Each position still maps back to the residue of our own structure.
    assert [a.query_residue.seq_index for a in aligned] == [0, 1, 2]


def test_align_hit_to_query_refuses_rows_it_cannot_align():
    """Length mismatch means the correspondence is genuinely unrecoverable.

    Measured on real output: this build omits some unmatched positions rather
    than padding them, so refusing is the only honest behaviour.
    """
    index = ResidueIndex([rid("B", 57, index=0), rid("B", 102, index=1)])
    with pytest.raises(ContractError, match="no positional alignment exists"):
        align_hit_to_query("B57,B102", "B57", index)


# --- safe hit parsing, against real observed rows ---------------------------

@pytest.mark.parametrize(
    ("query", "matching", "n_matched", "n_placeholders", "alignable"),
    [
        # Real rows captured from `folddisco query` output over 545 emitted
        # queries against the 4CHA test index.
        ("A1,A3,A4", "A1,A3,A4", 3, 0, True),
        ("A9,A10,A11", "A9,A10", 2, 0, False),
        ("A7,A8,A9,A10,B23", "F30,_,_,F33,_", 2, 3, True),
        ("B146,C217,C218,C219,C221", "C217,C218,C219,C221", 4, 0, False),
        ("B146,C217,C218,C219,C221", "B75,B76,_,_", 2, 2, False),
        ("A10,A11,B135,C159", "A10,B135,C159", 3, 0, False),
    ],
)
def test_parse_hit_residues_handles_every_observed_pattern(
    query, matching, n_matched, n_placeholders, alignable
):
    parsed = parse_hit_residues(query, matching)
    assert parsed.n_matched == n_matched
    assert parsed.n_unmatched_placeholders == n_placeholders
    assert parsed.positionally_alignable is alignable
    # Placeholders never leak into the matched set.
    assert UNMATCHED_TOKEN not in parsed.matched_tokens
    # Every matched token is a parseable residue token.
    for token in parsed.matched_tokens:
        parse_query_token(token)


def test_parse_hit_residues_never_raises_on_real_patterns():
    """The safe parser must cope with every row, unlike the aligner."""
    rows = [
        ("A1,A3,A4", "A1,A3,A4"),
        ("A9,A10,A11", "A9,A10"),
        ("B146,C217,C218,C219,C221", "B75,B76,_,_"),
        ("A7,A8", "_,_"),
    ]
    for query, matching in rows:
        parsed = parse_hit_residues(query, matching)
        assert parsed.query_tokens
        assert parsed.raw_matching_field == matching

    # A fully unmatched row yields zero matched residues, not an error.
    empty = parse_hit_residues("A7,A8", "_,_")
    assert empty.n_matched == 0
    assert empty.n_unmatched_placeholders == 2
