"""Contract regression tests against the real FoldDisco executable.

Opt-in: ``pytest -m folddisco``. Skipped automatically when the executable is
absent, so the default suite stays offline and self-contained.

These exist because the query contract is derived from one specific build.
If that build changes, the rules in :mod:`ibc_search.folddisco_contract` may
silently stop holding, and the only way to notice is to ask the binary.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

from ibc_search.cli import main
from ibc_search.folddisco_contract import (
    FOLDDISCO_VERIFIED_COMMIT,
    MIN_QUERY_RESIDUES,
    parse_hit_residues,
)

pytestmark = pytest.mark.folddisco

FOLDDISCO = shutil.which("folddisco")
TEST_INDEX = Path(
    "/Users/mvh_mac/Documents/master_thesis_2026/folddisco_test/test_index"
)
REF_PDB = Path(
    "/Users/mvh_mac/Documents/master_thesis_2026/folddisco_test/pdb_files/4CHA.pdb"
)

requires_folddisco = pytest.mark.skipif(
    FOLDDISCO is None, reason="folddisco executable not on PATH"
)
requires_index = pytest.mark.skipif(
    not TEST_INDEX.with_suffix(".lookup").is_file(),
    reason="folddisco_test index not available",
)


def run_folddisco(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(  # noqa: S603 - fixed argv, no shell
        [FOLDDISCO, *args], capture_output=True, text=True, timeout=120, check=False
    )


@requires_folddisco
def test_installed_build_still_matches_the_verified_contract():
    """A build change invalidates every rule derived from it."""
    result = run_folddisco("version")
    reported = result.stdout.strip().splitlines()[-1].strip()
    assert reported == FOLDDISCO_VERIFIED_COMMIT, (
        f"installed FoldDisco is {reported!r} but the query contract was "
        f"verified against {FOLDDISCO_VERIFIED_COMMIT!r}. Re-verify "
        "docs/FOLDDISCO_CONTRACT.md before trusting emitted queries."
    )


@requires_folddisco
@requires_index
def test_single_residue_query_still_returns_nothing():
    """The basis of MIN_QUERY_RESIDUES = 2."""
    assert MIN_QUERY_RESIDUES == 2
    one = run_folddisco(
        "query", "-i", str(TEST_INDEX), "-p", str(REF_PDB), "-q", "B57", "--per-match"
    )
    assert one.returncode == 0
    assert one.stdout.strip() == "", "a 1-residue query unexpectedly returned rows"

    two = run_folddisco(
        "query", "-i", str(TEST_INDEX), "-p", str(REF_PDB), "-q", "B57,B102",
        "--per-match",
    )
    assert two.returncode == 0
    assert two.stdout.strip(), "a 2-residue query must return rows"


@requires_folddisco
@requires_index
@pytest.mark.parametrize("bad_query", ["B57A", "AAA57", "B-5"])
def test_inexpressible_tokens_are_still_rejected(bad_query):
    """These rejections are why representability is validated before emitting."""
    result = run_folddisco(
        "query", "-i", str(TEST_INDEX), "-p", str(REF_PDB), "-q", bad_query,
        "--per-match",
    )
    assert result.returncode != 0
    assert "panicked" in result.stderr or "Invalid" in result.stderr


@requires_folddisco
@requires_index
def test_emitted_query_file_is_consumed_without_error(tmp_path):
    """End-to-end: every query this package emits must be accepted.

    This is the test that would have caught a bad token format, and it is the
    reason the query file format is trusted.
    """
    out = tmp_path / "out"
    assert main(
        ["scan", str(REF_PDB), "-o", str(out),
         "--chain", "A", "--chain", "B", "--chain", "C", "--radius", "8"]
    ) == 0

    query_file = out / "folddisco_queries.tsv"
    lines = query_file.read_text().splitlines()
    assert lines, "no queries were emitted"

    result = run_folddisco(
        "query", "-i", str(TEST_INDEX), "-q", str(query_file), "--per-match"
    )
    assert result.returncode == 0, result.stderr[:500]
    assert "panicked" not in result.stderr
    assert result.stderr.strip() == ""
    assert result.stdout.strip(), "the emitted queries returned no hits at all"


@requires_folddisco
@requires_index
def test_hit_rows_parse_and_query_tokens_round_trip(tmp_path):
    """Every hit row must be readable, and every query token must resolve.

    Uses `parse_hit_residues`, which makes no positional-alignment assumption
    -- the measured behaviour is that only ~62% of rows are alignable.
    """
    from ibc_search.config import load_config
    from ibc_search.folddisco_contract import resolve_query_token
    from ibc_search.structure_io import parse_structure, resolve_input

    config = load_config()
    config["input"]["chains"] = ["A", "B", "C"]
    config["ball_search"]["radii_angstrom"] = [8.0]
    config["run"]["output_dir"] = str(tmp_path / "out")

    parsed = parse_structure(resolve_input(str(REF_PDB), config), config)
    index = parsed.residue_index

    out = tmp_path / "out2"
    assert main(
        ["scan", str(REF_PDB), "-o", str(out),
         "--chain", "A", "--chain", "B", "--chain", "C", "--radius", "8"]
    ) == 0
    result = run_folddisco(
        "query", "-i", str(TEST_INDEX), "-q",
        str(out / "folddisco_queries.tsv"), "--per-match",
    )
    assert result.returncode == 0

    n_rows = 0
    n_alignable = 0
    for line in result.stdout.splitlines():
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 6 or fields[0] == "tid":
            continue
        n_rows += 1
        hit = parse_hit_residues(fields[5], fields[4])
        n_alignable += hit.positionally_alignable
        assert hit.query_tokens
        # The query side is an exact echo, so it must always resolve.
        for token in hit.query_tokens:
            residue, note = resolve_query_token(token, index)
            assert residue is not None, f"{token!r} did not resolve: {note}"

    assert n_rows > 0
    # Documented reality, asserted so a build change that "fixes" this is noticed.
    assert 0 < n_alignable < n_rows, (
        f"{n_alignable}/{n_rows} rows alignable; the contract documents a mix "
        "of padded and truncated rows"
    )
