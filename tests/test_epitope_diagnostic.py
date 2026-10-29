"""Regression checks for the exploratory SAbDab diagnostic, not biological claims."""

import importlib.util
import json
import math
from pathlib import Path
from types import SimpleNamespace

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/epitope_diagnostic.py"
spec = importlib.util.spec_from_file_location("epitope_diagnostic", SCRIPT)
d = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d)


def match(**changes):
    return dict(inst="target", node_count=2, m="A1,A2", rmsd=0.2, idf=10, **changes)


def labels():
    return {("target", "A", 1): [["", 1, 1]],
            ("target", "A", 2): [["", 0, 1]]}


@pytest.mark.parametrize("field,value,reason", [
    ("node_count", 1, "partial"), ("node_count", 3, "excess_nodes"),
    ("m", "A1,A1", "repeated_target"),
    ("m", "A1", "target_count_mismatch"),
    ("m", "A1,A2,A3", "target_count_mismatch"),
    ("m", "A1,A-2", "unparseable"),
    ("rmsd", math.nan, "nonfinite_or_negative_rmsd"),
    ("rmsd", math.inf, "nonfinite_or_negative_rmsd"),
    ("rmsd", -1, "nonfinite_or_negative_rmsd"),
    ("rmsd", 2, "rmsd"),
    ("idf", math.nan, "nonfinite_or_negative_idf"),
])
def test_invalid_match(field, value, reason):
    row = match()
    row[field] = value
    assert d.classify_match(row, 2, labels(), 1.5)[0] == reason


def test_valid_and_ambiguous_labels():
    lab = labels()
    assert d.classify_match(match(), 2, lab, 1.5) == (None, 2, 1)
    lab[("target", "A", 1)].append(["A", 0, 1])
    assert d.classify_match(match(), 2, lab, 1.5)[0] == "insertion_code_ambiguous"
    lab[("target", "A", 1)] = [["", 1, None]]
    assert d.classify_match(match(), 2, lab, 1.5)[0] == "unresolved_label"


def test_family_components_include_transitive_bridges():
    meta = {"a": {"ag_cluster": "P1/P2"}, "b": {"ag_cluster": "P2/P3"},
            "c": {"ag_cluster": "P3"}, "d": {"ag_cluster": "P9"}}
    groups = d.family_components(meta)
    assert groups["a"] == groups["b"] == groups["c"]
    assert groups["c"] != groups["d"]
    with pytest.raises(ValueError):
        d.family_of("P1/P2")


def test_cache_requires_matching_signature_and_intact_artifacts(tmp_path):
    artifact = tmp_path / "hits.tsv"
    artifact.write_text("tid\n")  # Successful zero-hit search is cacheable.
    signature = {"radius": 6, "top": 50}
    record = {"signature": signature, "artifacts": {"hits.tsv": d.digest(artifact)}}
    (tmp_path / "completed.json").write_text(json.dumps(record))
    assert d.cached_run(tmp_path, signature)
    assert not d.cached_run(tmp_path, {**signature, "radius": 5})
    artifact.write_text("truncated")
    assert not d.cached_run(tmp_path, signature)


def test_missing_marker_is_not_complete(tmp_path):
    (tmp_path / "hits.tsv").write_text("some output\n")
    assert not d.cached_run(tmp_path, {})


def test_no_evidence_and_test_partition_exclusion():
    meta = {
        "own": dict(pdb_id="p1", ag_cluster="P1", family_group="P1",
                    ab_ag_split="train", antigen_chains="A"),
        "target": dict(pdb_id="p2", ag_cluster="P2", family_group="P2",
                       ab_ag_split="test", antigen_chains="A"),
    }
    result = d.analyse_query([match()], 2, "own", meta, labels(), {"target": .5},
                             50, SimpleNamespace(max_rmsd=1.5, min_exposed=2), {})
    assert result["retained"] == 0
    assert result["enrichment"] is None
    assert result["families"] == 0


def test_crlf_repeated_headers_and_malformed_rows(tmp_path):
    path = tmp_path / "hits.tsv"
    header = "tid\tnode_count\tidf\trmsd\tquery_residues\tmatching_residues\r\n"
    path.write_bytes((header + "target.cif\t2\t10\t.2\tA1,A2\tA1,A2\r\n" +
                      header + "broken\r\n").encode())
    rows, rejected = d.load_hits(path)
    assert len(rows) == 1 and rejected == 1
    assert rows[0]["m"] == "A1,A2"


def test_paired_evaluation_uses_identical_queries():
    def query(name, support, enrichment):
        return dict(antigen="a", query=name, families=support, enrichment=enrichment,
                    contact_fraction=.5, idf_log10_max=3, valid_row_fraction=.1)
    rows = [query("supported", 3, .2), query("unsupported", 1, .1),
            query("missing", 3, None)]
    for field in d.SCORE_FIELDS:
        selected = d.evaluation_subset(rows, "a", field, 3, False, True)
        assert [q["query"] for q in selected] == ["supported"]
    assert len(d.evaluation_subset(rows, "a", "idf_log10_max", 3, False, False)) == 3


def test_family_balancing_and_rejections():
    meta = {"own": dict(pdb_id="own", ag_cluster="P0", family_group="P0",
                         ab_ag_split="train", antigen_chains="A")}
    lab = {}
    rows = []
    for inst, family, contact in [("t1", "P1", 1), ("t2", "P1", 1), ("t3", "P2", 0)]:
        meta[inst] = dict(pdb_id=inst, ag_cluster=family, family_group=family,
                          ab_ag_split="train", antigen_chains="A")
        for n in (1, 2):
            lab[(inst, "A", n)] = [["", contact, 1]]
        row = match()
        row["inst"] = inst
        rows.extend([row, row.copy()])
    result = d.analyse_query(rows, 2, "own", meta, lab,
                             {i: .25 for i in ("t1", "t2", "t3")}, 50,
                             SimpleNamespace(max_rmsd=1.5, min_exposed=2), {})
    assert result["families"] == 2
    assert result["contact_fraction"] == .5
    assert result["enrichment"] == .25