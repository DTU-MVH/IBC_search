"""Tests for scripts/folddisco_search.sh.

These run offline. Rather than invoke the real binary, they put a stub named
``folddisco`` first on ``PATH`` that replays canned output, which makes the
interesting cases reproducible: a query returning nothing, two patches sharing
one query string, and a row torn by the interleaved-write bug that this build
exhibits above one thread.

The contract-regression tests that do exercise the real binary live in
tests/test_folddisco_live.py.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "folddisco_search.sh"

HEADER = "tid\tidf\trmsd\tnode_count\tquery_residues\tmatching_residues"


def test_script_exists_and_is_executable():
    assert SCRIPT.is_file(), f"missing {SCRIPT}"
    assert os.access(SCRIPT, os.X_OK), f"{SCRIPT} is not executable"


def test_script_is_syntactically_valid():
    completed = subprocess.run(
        ["bash", "-n", str(SCRIPT)], capture_output=True, text=True
    )
    assert completed.returncode == 0, completed.stderr


def _make_results_dir(
    tmp_path: Path, rows: list[tuple[str, str]], structure_id: str | None = None
) -> Path:
    """Write the two files a scan leaves behind that the script consumes.

    ``rows`` is a list of ``(patch_id, query_string)``.
    """
    results = tmp_path / "results"
    results.mkdir(parents=True, exist_ok=True)
    (results / "folddisco_queries.tsv").write_text(
        "".join(f"/fake/structure.pdb\t{q}\n" for _, q in rows), encoding="utf-8"
    )
    lines = ["line_number\tpatch_id\tstructure_path\tradius_angstrom\tscale\tn_members\tquery_string"]
    for i, (patch_id, query) in enumerate(rows, start=1):
        lines.append(
            f"{i}\t{patch_id}\t/fake/structure.pdb\t5.0000\tsubmotif\t"
            f"{len(query.split(','))}\t{query}"
        )
    (results / "folddisco_queries.index.tsv").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    if structure_id is not None:
        (results / "run_manifest.json").write_text(
            json.dumps({"input": {"structure_id": structure_id}}), encoding="utf-8"
        )
    return results


def _stub_folddisco(tmp_path: Path, output: str) -> dict[str, str]:
    """Put a fake `folddisco` first on PATH and return the patched env."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    stub = bin_dir / "folddisco"
    payload = output.replace("'", "'\\''")
    stub.write_text(
        "#!/usr/bin/env bash\n"
        'if [[ "$1" == "version" ]]; then\n'
        # Report the verified commit so the script does not warn.
        "  echo 2a756d9522d5bab8f41bc73ec72ec2327e31cd97\n"
        "  exit 0\n"
        "fi\n"
        f"printf '%s' '{payload}'\n",
        encoding="utf-8",
    )
    stub.chmod(0o755)
    env = dict(os.environ)
    env["PATH"] = f"{bin_dir}{os.pathsep}{env['PATH']}"
    return env


def _run(results: Path, index: Path, env: dict[str, str], *extra: str):
    return subprocess.run(
        [str(SCRIPT), "search", str(results), "-i", str(index), *extra],
        capture_output=True,
        text=True,
        env=env,
    )


@pytest.fixture()
def index_file(tmp_path: Path) -> Path:
    path = tmp_path / "index"
    path.write_text("not a real index\n", encoding="utf-8")
    return path


def test_hits_are_joined_back_to_their_patch(tmp_path: Path, index_file: Path):
    results = _make_results_dir(
        tmp_path, [("P_a", "A1,A2,A3"), ("P_b", "B7,B8,B9")]
    )
    # Second query returns nothing: its header block appears with no rows.
    output = (
        f"{HEADER}\n"
        "4cha.pdb\t0.5\t0.17\t2\tA1,A2,A3\tB96,B97,_\n"
        f"{HEADER}\n"
    )
    completed = _run(results, index_file, _stub_folddisco(tmp_path, output))
    assert completed.returncode == 0, completed.stderr

    lines = (results / "folddisco_hits.tsv").read_text().strip().split("\n")
    assert lines[0].split("\t")[:3] == ["patch_id", "radius_angstrom", "scale"]
    assert len(lines) == 2, "one hit row expected"
    assert lines[1].startswith("P_a\t5.0000\tsubmotif\t")
    assert "queries searched: 2 of 2, with hits: 1, with none: 1" in completed.stderr
    assert "PARTIAL" not in completed.stderr


def test_join_does_not_depend_on_output_order(tmp_path: Path, index_file: Path):
    """Blocks arrive in completion order at -t > 1, so order must not matter."""
    results = _make_results_dir(
        tmp_path, [("P_first", "A1,A2,A3"), ("P_second", "B7,B8,B9")]
    )
    # The SECOND query's block is emitted first.
    output = (
        f"{HEADER}\n"
        "4cha.pdb\t0.5\t0.17\t2\tB7,B8,B9\tC1,C2,_\n"
        f"{HEADER}\n"
        "4cha.pdb\t0.5\t0.31\t2\tA1,A2,A3\tD1,D2,_\n"
    )
    completed = _run(results, index_file, _stub_folddisco(tmp_path, output))
    assert completed.returncode == 0, completed.stderr

    rows = [
        line.split("\t")
        for line in (results / "folddisco_hits.tsv").read_text().strip().split("\n")[1:]
    ]
    by_patch = {r[0]: r for r in rows}
    assert by_patch["P_second"][-2] == "B7,B8,B9"
    assert by_patch["P_first"][-2] == "A1,A2,A3"


def test_a_query_string_shared_by_two_patches_lists_both(
    tmp_path: Path, index_file: Path
):
    """Identical member sets at two radii are indistinguishable in the output.

    The script must say so rather than pick one.
    """
    results = _make_results_dir(
        tmp_path, [("P_r5", "A1,A2,A3"), ("P_r8", "A1,A2,A3")]
    )
    output = f"{HEADER}\n4cha.pdb\t0.5\t0.17\t2\tA1,A2,A3\tB96,B97,_\n"
    completed = _run(results, index_file, _stub_folddisco(tmp_path, output))
    assert completed.returncode == 0, completed.stderr

    row = (results / "folddisco_hits.tsv").read_text().strip().split("\n")[1]
    assert row.split("\t")[0] == "P_r5|P_r8"
    assert "shared by >1 patch" in completed.stderr


def test_torn_rows_are_quarantined_not_parsed(tmp_path: Path, index_file: Path):
    """A row with the wrong field count is unparseable and must not be emitted."""
    results = _make_results_dir(tmp_path, [("P_a", "A1,A2,A3")])
    output = (
        f"{HEADER}\n"
        "4cha.pdb\t0.5\t0.17\t2\tA1,A2,A3\tB96,B97,_\n"
        "4cha.pdb\t0.5\t0.2\tA1,A2\n"              # torn: too few fields
    )
    completed = _run(results, index_file, _stub_folddisco(tmp_path, output))
    assert completed.returncode == 0, completed.stderr

    body = (results / "folddisco_hits.tsv").read_text().strip().split("\n")[1:]
    assert len(body) == 1, "the torn row must not reach the hits table"
    assert "1 rows were torn" in completed.stderr
    quarantine = results / "folddisco_hits.malformed.tsv"
    assert quarantine.is_file() and "A1,A2" in quarantine.read_text()


def test_strict_fails_on_a_torn_row(tmp_path: Path, index_file: Path):
    results = _make_results_dir(tmp_path, [("P_a", "A1,A2,A3")])
    output = f"{HEADER}\n4cha.pdb\t0.5\t0.2\tA1,A2\n"
    completed = _run(
        results, index_file, _stub_folddisco(tmp_path, output), "--strict"
    )
    assert completed.returncode != 0
    assert "--strict" in completed.stderr


def test_missing_index_is_a_clear_error(tmp_path: Path):
    results = _make_results_dir(tmp_path, [("P_a", "A1,A2,A3")])
    env = _stub_folddisco(tmp_path, "")
    completed = _run(results, tmp_path / "nope", env)
    assert completed.returncode == 2
    assert "index not found" in completed.stderr
    assert "extract" in completed.stderr, "should point at the unpack step"


def test_a_scan_with_no_queries_is_refused(tmp_path: Path, index_file: Path):
    results = _make_results_dir(tmp_path, [])
    completed = _run(results, index_file, _stub_folddisco(tmp_path, ""))
    assert completed.returncode == 2
    assert "emitted no expressible queries" in completed.stderr


def test_format_output_without_query_residues_is_refused(
    tmp_path: Path, index_file: Path
):
    """Dropping the join key must fail loudly, not produce unjoinable hits."""
    results = _make_results_dir(tmp_path, [("P_a", "A1,A2,A3")])
    output = "tid\trmsd\n4cha.pdb\t0.17\n"
    completed = _run(
        results,
        index_file,
        _stub_folddisco(tmp_path, output),
        "--format-output",
        "tid,rmsd",
    )
    assert "must include query_residues" in completed.stderr


@pytest.mark.skipif(shutil.which("folddisco") is None, reason="folddisco not on PATH")
def test_status_reports_the_installed_build():
    completed = subprocess.run(
        [str(SCRIPT), "status"], capture_output=True, text=True
    )
    assert completed.returncode == 0, completed.stderr
    assert "folddisco binary" in completed.stdout
    assert "index" in completed.stdout


def test_partial_output_is_not_reported_as_empty_queries(
    tmp_path: Path, index_file: Path
):
    """An interrupted search must not look like queries that found nothing.

    A long search can be stopped part-way and the completed blocks are still
    worth joining, but "never ran" and "ran and matched nothing" are different
    facts and conflating them would silently understate coverage.
    """
    results = _make_results_dir(
        tmp_path,
        [("P_a", "A1,A2,A3"), ("P_b", "B7,B8,B9"), ("P_c", "C1,C2,C3")],
    )
    # Only the first query's block is present: the run was cut short.
    output = f"{HEADER}\n4cha.pdb\t0.5\t0.17\t2\tA1,A2,A3\tB96,B97,_\n"
    completed = _run(results, index_file, _stub_folddisco(tmp_path, output))
    assert completed.returncode == 0, completed.stderr

    assert "queries searched: 1 of 3" in completed.stderr
    assert "PARTIAL -- 2 of 3 queries were never searched" in completed.stderr


def test_join_runs_standalone_on_an_existing_raw_file(
    tmp_path: Path, index_file: Path
):
    """`join` must work without re-searching, for interrupted runs."""
    results = _make_results_dir(tmp_path, [("P_a", "A1,A2,A3")])
    raw = results / "folddisco_hits.raw.tsv"
    raw.write_text(
        f"{HEADER}\n4cha.pdb\t0.5\t0.17\t2\tA1,A2,A3\tB96,B97,_\n", encoding="utf-8"
    )
    completed = subprocess.run(
        [str(SCRIPT), "join", str(results)], capture_output=True, text=True
    )
    assert completed.returncode == 0, completed.stderr
    row = (results / "folddisco_hits.tsv").read_text().strip().split("\n")[1]
    assert row.startswith("P_a\t")


# --- self-hit exclusion ---------------------------------------------------

SELF_HEADER = "tid\tidf\trmsd\tnode_count\tquery_residues\tmatching_residues"


def _self_output() -> str:
    """Canned output mixing the scanned entry with two other structures."""
    return (
        f"{SELF_HEADER}\n"
        "pdb/n7/pdb2n74.ent.gz\t9.9\t0.00\t3\tA1,A2,A3\tA1,A2,A3\n"
        "pdb/ab/pdb1abc.ent.gz\t5.0\t0.21\t3\tA1,A2,A3\tB7,B8,B9\n"
        "pdb/cd/pdb9xyz.ent.gz\t4.0\t0.33\t2\tA1,A2,A3\tC1,_,C3\n"
    )


def test_self_hits_are_excluded_by_default(tmp_path: Path, index_file: Path):
    """Searching a deposited entry against the whole PDB finds it matching itself.

    Those rows are real but circular: on 2N74 the entry was its own top hit
    for 56 of 87 motifs.
    """
    results = _make_results_dir(tmp_path, [("P_a", "A1,A2,A3")], structure_id="2N74")
    completed = _run(results, index_file, _stub_folddisco(tmp_path, _self_output()))
    assert completed.returncode == 0, completed.stderr

    body = (results / "folddisco_hits.tsv").read_text()
    assert "pdb2n74" not in body, "the entry itself must not remain in the hits"
    assert "pdb1abc" in body and "pdb9xyz" in body, "other structures must survive"
    assert "excluded 1 self-hit" in completed.stderr
    # Excluded rows are auditable, not silently dropped.
    side = results / "folddisco_hits.selfhits.tsv"
    assert side.is_file() and "pdb2n74" in side.read_text()


def test_no_exclude_self_keeps_them(tmp_path: Path, index_file: Path):
    results = _make_results_dir(tmp_path, [("P_a", "A1,A2,A3")], structure_id="2N74")
    completed = _run(
        results, index_file, _stub_folddisco(tmp_path, _self_output()),
        "--no-exclude-self",
    )
    assert completed.returncode == 0, completed.stderr
    assert "pdb2n74" in (results / "folddisco_hits.tsv").read_text()


def test_exclude_id_drops_additional_entries(tmp_path: Path, index_file: Path):
    """Homologues cannot be inferred from an entry ID, so they are opt-in."""
    results = _make_results_dir(tmp_path, [("P_a", "A1,A2,A3")], structure_id="2N74")
    completed = _run(
        results, index_file, _stub_folddisco(tmp_path, _self_output()),
        "--exclude-id", "1ABC",
    )
    assert completed.returncode == 0, completed.stderr
    body = (results / "folddisco_hits.tsv").read_text()
    assert "pdb2n74" not in body and "pdb1abc" not in body
    assert "pdb9xyz" in body


def test_a_non_entry_structure_id_excludes_nothing(tmp_path: Path, index_file: Path):
    """A design or local model is not a deposited entry; nothing is guessed."""
    results = _make_results_dir(
        tmp_path, [("P_a", "A1,A2,A3")], structure_id="member_01"
    )
    completed = _run(results, index_file, _stub_folddisco(tmp_path, _self_output()))
    assert completed.returncode == 0, completed.stderr
    assert "pdb2n74" in (results / "folddisco_hits.tsv").read_text()
    assert "excluded" not in completed.stderr


# --- parallel ensemble ----------------------------------------------------


def _ensemble(tmp_path: Path, n: int) -> Path:
    root = tmp_path / "ens"
    root.mkdir()
    for i in range(1, n + 1):
        d = _make_results_dir(tmp_path / f"m{i}", [(f"P{i}", "A1,A2,A3")])
        d.rename(root / f"run_{i:02d}")
    return root


def test_search_ensemble_runs_every_member(tmp_path: Path, index_file: Path):
    root = _ensemble(tmp_path, 3)
    env = _stub_folddisco(tmp_path, _self_output())
    completed = subprocess.run(
        [str(SCRIPT), "search-ensemble", str(root), "-i", str(index_file), "-j", "3"],
        capture_output=True, text=True, env=env,
    )
    assert completed.returncode == 0, completed.stderr
    assert "3 members" in completed.stderr
    for i in (1, 2, 3):
        assert (root / f"run_{i:02d}" / "folddisco_hits.tsv").is_file()


def test_search_ensemble_declines_to_parallelise_a_single_run(
    tmp_path: Path, index_file: Path
):
    """Parallelism is gated on actually finding an ensemble."""
    root = _ensemble(tmp_path, 1)
    env = _stub_folddisco(tmp_path, _self_output())
    completed = subprocess.run(
        [str(SCRIPT), "search-ensemble", str(root), "-i", str(index_file)],
        capture_output=True, text=True, env=env,
    )
    assert completed.returncode == 0, completed.stderr
    assert "not an ensemble" in completed.stderr
    assert (root / "run_01" / "folddisco_hits.tsv").is_file()


def test_search_ensemble_needs_scanned_members(tmp_path: Path, index_file: Path):
    empty = tmp_path / "nothing"
    empty.mkdir()
    completed = subprocess.run(
        [str(SCRIPT), "search-ensemble", str(empty), "-i", str(index_file)],
        capture_output=True, text=True,
    )
    assert completed.returncode == 2
    assert "no run_* directories" in completed.stderr
