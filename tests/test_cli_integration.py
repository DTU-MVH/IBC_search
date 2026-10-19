"""End-to-end CLI runs: output tables, manifest, query file, determinism."""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

import pytest

from ibc_search.cli import main

from .conftest import write_rotated_pdb


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def run_scan(args: list[str]) -> int:
    return main(["scan", *args])


@pytest.fixture
def scanned(hairpin_pdb, tmp_path) -> Path:
    out = tmp_path / "out"
    assert run_scan([str(hairpin_pdb), "-o", str(out), "--radius", "8"]) == 0
    return out


# --- outputs exist and are well formed -------------------------------------

def test_scan_writes_all_expected_outputs(scanned):
    for name in (
        "residues.tsv", "patches.tsv", "patch_members.tsv",
        "folddisco_queries.tsv", "folddisco_queries.index.tsv",
        "run_manifest.json",
    ):
        assert (scanned / name).is_file(), name


def test_residue_table_has_identity_and_exposure_columns(scanned):
    rows = read_tsv(scanned / "residues.tsv")
    assert rows
    required = {
        "structure_id", "model", "auth_chain_id", "auth_seq_id", "ins_code",
        "res_name", "seq_index", "sasa_angstrom2", "rsa",
        "exposed_at_seed_threshold", "quality_flags", "residue_class",
    }
    assert required <= set(rows[0])
    for row in rows:
        assert row["seq_index"].isdigit()
        assert row["res_name"]
        assert row["exposed_at_seed_threshold"] in ("true", "false")


def test_patch_table_carries_the_required_evidence(scanned):
    patches = read_tsv(scanned / "patches.tsv")
    members = read_tsv(scanned / "patch_members.tsv")
    assert patches and members

    required = {
        "patch_id", "radius_angstrom", "scale", "center_mode", "n_members",
        "n_segments", "segment_lengths", "author_numbering_gaps",
        "is_discontinuous", "seed_residue_labels",
        "spatial_diameter_angstrom", "radius_of_gyration_angstrom",
        "compactness_ratio", "member_set_hash", "cluster_id",
        "n_members_dropped_disconnected", "query_string",
    }
    assert required <= set(patches[0])

    # Every member row joins to a patch, and each records its selecting atom
    # and that atom's distance -- the per-member evidence of inclusion.
    patch_ids = {row["patch_id"] for row in patches}
    for row in members:
        assert row["patch_id"] in patch_ids
        assert row["selecting_atom_name"]
        assert float(row["selecting_atom_distance_angstrom"]) >= 0.0

    # Member counts agree between the two tables.
    counts: dict[str, int] = {}
    for row in members:
        counts[row["patch_id"]] = counts.get(row["patch_id"], 0) + 1
    for row in patches:
        assert counts[row["patch_id"]] == int(row["n_members"])


def test_patch_ids_are_unique(scanned):
    patches = read_tsv(scanned / "patches.tsv")
    ids = [row["patch_id"] for row in patches]
    assert len(ids) == len(set(ids))


def test_query_file_matches_the_verified_folddisco_format(scanned):
    lines = (scanned / "folddisco_queries.tsv").read_text().splitlines()
    index = read_tsv(scanned / "folddisco_queries.index.tsv")
    assert lines
    assert len(lines) == len(index)

    for line_number, (line, row) in enumerate(zip(lines, index), start=1):
        path, query = line.split("\t")
        assert Path(path).is_file()
        assert query == row["query_string"]
        assert int(row["line_number"]) == line_number
        tokens = query.split(",")
        # The verified floor: at least two residues per query.
        assert len(tokens) >= 2
        for token in tokens:
            assert token[0].isalnum()
            assert token[1:].isdigit()


def test_query_index_links_every_emitted_patch_back_to_its_id(scanned):
    patches = {row["patch_id"]: row for row in read_tsv(scanned / "patches.tsv")}
    index = read_tsv(scanned / "folddisco_queries.index.tsv")

    emitted = {row["patch_id"] for row in index}
    expected = {
        pid for pid, row in patches.items() if row["emit_query"] == "true"
    }
    assert emitted == expected
    for row in index:
        assert row["query_string"] == patches[row["patch_id"]]["query_string"]


# --- manifest ---------------------------------------------------------------

def test_manifest_records_provenance_config_timings_and_warnings(scanned, hairpin_pdb):
    manifest = json.loads((scanned / "run_manifest.json").read_text())

    assert manifest["input"]["kind"] == "file"
    assert manifest["input"]["path"] == str(hairpin_pdb.resolve())
    assert len(manifest["input"]["sha256"]) == 64
    assert manifest["command"][0] == "ibc-search"
    assert "scan" in manifest["command"]

    # Versions of the tools that produced the numbers.
    for key in ("ibc_search", "python", "numpy", "scipy", "biotite", "platform"):
        assert manifest["versions"][key]

    # The fully resolved configuration, so the run reproduces from its output.
    assert manifest["resolved_config"]["ball_search"]["radii_angstrom"] == [8.0]
    assert manifest["random_seed"] == 0

    # Per-stage timings, so a regression is visible.
    for stage in (
        "parse_structure_s", "compute_exposure_s", "ball_search_total_s",
        "write_tables_s",
    ):
        assert manifest["timings_seconds"][stage] >= 0.0

    assert "warnings" in manifest and isinstance(manifest["n_warnings"], int)
    assert manifest["folddisco"]["search_executed_by_this_component"] is False
    assert "not epitopes" in manifest["interpretation_notice"]


def test_manifest_records_exposure_method_and_reference(scanned):
    manifest = json.loads((scanned / "run_manifest.json").read_text())
    exposure = manifest["exposure"]
    assert "Shrake-Rupley" in exposure["method"]
    assert exposure["probe_radius_angstrom"] == 1.4
    assert "Tien" in exposure["rsa_reference_citation"]
    assert "e80635" in exposure["rsa_reference_citation"]


def test_manifest_records_models_chains_and_assembly_caveat(scanned):
    manifest = json.loads((scanned / "run_manifest.json").read_text())
    structure = manifest["structure"]
    assert structure["model_scanned"] == 1
    assert structure["models_available"] == [1]
    assert structure["chains_scanned"] == ["A"]
    assert "asymmetric unit" in structure["biological_assembly"]


# --- determinism ------------------------------------------------------------

def test_repeated_runs_are_byte_identical(hairpin_pdb, tmp_path):
    first, second = tmp_path / "a", tmp_path / "b"
    assert run_scan([str(hairpin_pdb), "-o", str(first), "--radius", "8"]) == 0
    assert run_scan([str(hairpin_pdb), "-o", str(second), "--radius", "8"]) == 0

    for name in ("residues.tsv", "patches.tsv", "patch_members.tsv",
                 "folddisco_queries.tsv"):
        assert (first / name).read_bytes() == (second / name).read_bytes(), name


@pytest.mark.parametrize("workers", [2, 4])
def test_worker_count_does_not_change_output_bytes(hairpin_pdb, tmp_path, workers):
    single, parallel = tmp_path / "w1", tmp_path / f"w{workers}"
    assert run_scan(
        [str(hairpin_pdb), "-o", str(single), "--radius", "8", "--workers", "1"]
    ) == 0
    assert run_scan(
        [str(hairpin_pdb), "-o", str(parallel), "--radius", "8",
         "--workers", str(workers)]
    ) == 0

    for name in ("patches.tsv", "patch_members.tsv", "folddisco_queries.tsv"):
        assert (single / name).read_bytes() == (parallel / name).read_bytes(), name


# --- input routes -----------------------------------------------------------

def test_scan_accepts_a_gzipped_file(hairpin_pdb, tmp_path):
    gz = tmp_path / "hairpin.pdb.gz"
    gz.write_bytes(gzip.compress(hairpin_pdb.read_bytes()))
    out = tmp_path / "out_gz"
    assert run_scan([str(gz), "-o", str(out), "--radius", "8"]) == 0
    assert read_tsv(out / "patches.tsv")


def test_scan_from_cached_pdb_id_without_network(
    hairpin_pdb, tmp_path, monkeypatch
):
    """The PDB-ID route, exercised offline from a pre-populated cache."""
    import urllib.request

    cache = tmp_path / "cache"
    cache.mkdir()
    (cache / "1ABC.pdb.gz").write_bytes(gzip.compress(hairpin_pdb.read_bytes()))

    def explode(*args, **kwargs):
        raise AssertionError("cached entry must not trigger a download")

    monkeypatch.setattr(urllib.request, "urlopen", explode)

    out = tmp_path / "out_id"
    assert run_scan(
        [
            "1ABC", "-o", str(out), "--radius", "8",
            "--set", f"input.cache_dir='{cache}'",
            "--set", "input.fetch_format='pdb'",
        ]
    ) == 0

    manifest = json.loads((out / "run_manifest.json").read_text())
    assert manifest["input"]["kind"] == "pdb_id"
    assert manifest["input"]["pdb_id"] == "1ABC"
    assert manifest["input"]["served_from_cache"] is True
    # Residue identities are stamped with the entry ID, not the filename.
    rows = read_tsv(out / "residues.tsv")
    assert {row["structure_id"] for row in rows} == {"1ABC"}


def test_unreachable_pdb_id_exits_with_an_error(tmp_path, monkeypatch, capsys):
    import urllib.error
    import urllib.request

    def unreachable(*args, **kwargs):
        raise urllib.error.URLError("simulated offline")

    monkeypatch.setattr(urllib.request, "urlopen", unreachable)
    code = run_scan(
        ["9XYZ", "-o", str(tmp_path / "out"),
         "--set", f"input.cache_dir='{tmp_path / 'c'}'"]
    )
    assert code == 2
    assert "does not fall back" in capsys.readouterr().err


def test_malformed_input_exits_with_an_error(tmp_path, capsys):
    bad = tmp_path / "bad.pdb"
    bad.write_text("not a structure\n")
    code = run_scan([str(bad), "-o", str(tmp_path / "out")])
    assert code == 2
    assert "error:" in capsys.readouterr().err


def test_nonexistent_input_exits_with_an_error(tmp_path, capsys):
    code = run_scan(["/no/such/file.pdb", "-o", str(tmp_path / "out")])
    assert code == 2
    assert "error:" in capsys.readouterr().err


# --- multi-chain / multi-model through the CLI ------------------------------

def test_chain_selection_through_the_cli(two_chain_pdb, tmp_path):
    out = tmp_path / "out_chain"
    assert run_scan(
        [str(two_chain_pdb), "-o", str(out), "--chain", "B", "--radius", "8"]
    ) == 0
    members = read_tsv(out / "patch_members.tsv")
    assert members
    assert {row["auth_chain_id"] for row in members} == {"B"}

    manifest = json.loads((out / "run_manifest.json").read_text())
    assert manifest["structure"]["chains_scanned"] == ["B"]
    assert manifest["structure"]["chains_available"] == ["A", "B"]


def test_model_selection_through_the_cli(multi_model_pdb, tmp_path):
    out = tmp_path / "out_model"
    assert run_scan(
        [str(multi_model_pdb), "-o", str(out), "--model", "2", "--radius", "8"]
    ) == 0
    manifest = json.loads((out / "run_manifest.json").read_text())
    assert manifest["structure"]["model_scanned"] == 2
    assert manifest["structure"]["models_available"] == [1, 2]
    rows = read_tsv(out / "residues.tsv")
    assert {row["model"] for row in rows} == {"2"}


# --- insertion codes end to end --------------------------------------------

def test_insertion_code_residues_are_not_emitted_as_queries(icode_pdb, tmp_path):
    """A residue with an insertion code cannot be written as a query token.

    It must still appear in the tables, with the reason recorded -- the motif
    is real, it simply cannot be handed to this FoldDisco build.
    """
    out = tmp_path / "out_icode"
    assert run_scan(
        [str(icode_pdb), "-o", str(out), "--radius", "8",
         "--set", "ball_search.min_patch_size=2"]
    ) == 0

    residues = read_tsv(out / "residues.tsv")
    icode_rows = [row for row in residues if row["ins_code"]]
    assert icode_rows
    for row in icode_rows:
        assert row["folddisco_expressible"] == "false"
        assert "insertion code" in row["folddisco_block_reason"]

    patches = read_tsv(out / "patches.tsv")
    assert patches
    blocked = [row for row in patches if row["folddisco_expressible"] == "false"]
    assert blocked, "patches containing icode residues must be marked"
    for row in blocked:
        assert row["emit_query"] == "false"
        assert row["query_string"] == ""

    # Those patches must not leak into the query file.
    queries = (out / "folddisco_queries.tsv").read_text()
    assert "10A" not in queries and "10B" not in queries

    manifest = json.loads((out / "run_manifest.json").read_text())
    assert manifest["structure"]["author_key_collisions"]


# --- other commands ---------------------------------------------------------

def test_show_config_prints_valid_json(capsys):
    assert main(["show-config"]) == 0
    config = json.loads(capsys.readouterr().out)
    assert config["exposure"]["probe_radius_angstrom"] == 1.4


def test_folddisco_info_reports_the_verified_contract(capsys):
    assert main(["folddisco-info"]) == 0
    info = json.loads(capsys.readouterr().out)
    contract = info["verified_query_contract"]
    assert contract["min_query_residues"] == 2
    assert contract["unmatched_position_token"] == "_"
    assert any("insertion codes" in item for item in contract["cannot_express"])
    assert info["note"] == "This component runs no FoldDisco search."


# --- rigid-motion invariance through the full CLI ---------------------------

def test_cli_output_is_invariant_under_rigid_motion(hairpin_pdb, tmp_path):
    """Patch membership and queries must not depend on the frame of reference."""
    rotated = write_rotated_pdb(hairpin_pdb, tmp_path / "rot.pdb")
    out_a, out_b = tmp_path / "ra", tmp_path / "rb"
    assert run_scan([str(hairpin_pdb), "-o", str(out_a), "--radius", "8"]) == 0
    assert run_scan([str(rotated), "-o", str(out_b), "--radius", "8"]) == 0

    def queries(path: Path) -> list[str]:
        return [
            line.split("\t")[1]
            for line in (path / "folddisco_queries.tsv").read_text().splitlines()
        ]

    assert queries(out_a) == queries(out_b)
