"""Query-side ranking tests (no FoldDisco index or network required)."""
import importlib.util
from pathlib import Path
import json
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("query_site_benchmark", ROOT / "scripts/query_site_benchmark.py")
benchmark = importlib.util.module_from_spec(spec)
spec.loader.exec_module(benchmark)


def test_ties_are_chance():
    m = benchmark.ranking_metrics([1, 0, 0, 0], [0, 0, 0, 0])
    assert m["auc"] == .5
    assert m["average_precision"] == .25
    assert m["enrichment_at_k"] == 1


def test_perfect_ranking():
    m = benchmark.ranking_metrics([1, 0, 0, 0], [1, 0, 0, 0])
    assert m["auc"] == 1
    assert m["average_precision"] == 1
    assert m["enrichment_at_k"] == 4
    assert benchmark.ranking_metrics([0, 0], [1, 2]) is None


def test_seeded_run_rejected(tmp_path):
    (tmp_path / "run_manifest.json").write_text(json.dumps({"ball_search": {"seeded": True}}))
    with pytest.raises(ValueError, match="annotation-guided"):
        benchmark.feature_rows(tmp_path, "1ABC")


def test_features_are_seed_centred_self_filtered_and_target_deduplicated(tmp_path):
    (tmp_path / "run_manifest.json").write_text(json.dumps({"ball_search": {"seeded": False}}))
    (tmp_path / "residues.tsv").write_text(
        "auth_chain_id\tauth_seq_id\tins_code\tresidue_class\tis_scanned_chain\trsa\n"
        "A\t1\t\tpolymer_standard_aa\ttrue\t0.5\n"
        "A\t2\t\tpolymer_standard_aa\ttrue\t0.2\n")
    (tmp_path / "patches.tsv").write_text(
        "patch_id\tseed_residue_labels\tn_members\temit_query\tquery_string\n"
        "p1\tA1\t4\ttrue\tA1,A2,A3,A4\n"
        "p2\tA2\t1\tfalse\t\n")
    (tmp_path / "folddisco_hits.tsv").write_text(
        "patch_id\ttid\tquery_residues\tmatching_residues\trmsd\n"
        "p1\t1abc.pdb\tA1,A2,A3,A4\tA1,A2,A3,A4\t0\n"
        "p1\t2abc.pdb\tA1,A2,A3,A4\tB1,B2,B3,B4\t0\n"
        "p1\t2abc.pdb\tA1,A2,A3,A4\tB1,B2,B3,B4\t1\n"
        "p1\t2abc.pdb\tA1,A2,A3,A4\tB5,B6\t2\n")
    f, audit, n = benchmark.feature_rows(tmp_path, "1ABC")
    assert n == 2
    assert f[0]["n_sites"] == 2
    assert f[0]["n_targets"] == 1
    assert f[0]["quality_support"] == 1
    assert f[0]["mean_coverage"] == .75
    assert f[1]["assessed"] is False
    assert f[1]["n_sites"] is None
    assert audit["self_excluded"] == 1
    assert audit["compact_target_field"] == 1


def test_finish_requires_completion_log(tmp_path):
    (tmp_path / "panel.tsv").write_text("instance\tcategory\n1ABC\tzincfingers\n")
    with pytest.raises(TimeoutError):
        benchmark.finish(SimpleNamespace(out=tmp_path, wait_seconds=0))
    assert not (tmp_path / "query_evaluation.json").exists()


def test_finish_records_hash_before_evaluation(tmp_path, monkeypatch):
    (tmp_path / "panel.tsv").write_text("instance\tcategory\n1ABC\tzincfingers\n")
    (tmp_path / "1ABC").mkdir()
    (tmp_path / "1ABC/folddisco_hits.tsv").write_text("header\n")
    (tmp_path / "1ABC.log").write_text("establishes nothing about antibody binding or immunogenicity.")
    called = []
    monkeypatch.setattr(benchmark, "evaluate", lambda args: called.append(args))
    benchmark.finish(SimpleNamespace(out=tmp_path, wait_seconds=0))
    marker = json.loads((tmp_path / "1ABC/search_complete.json").read_text())
    assert len(marker["hits_sha256"]) == 64
    assert len(called) == 1