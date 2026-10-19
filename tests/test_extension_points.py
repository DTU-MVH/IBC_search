"""Extension points must be explicit placeholders, not silent no-ops.

A stub that returned ``None`` would let a later caller believe a stage ran.
Each must raise :class:`NotImplementedError` with a reason.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from ibc_search.ext import benchmark, features, reference_index, scoring, search


def test_search_stage_is_an_explicit_placeholder():
    """Only the *in-process* search is a placeholder.

    Searching itself is done by scripts/folddisco_search.sh, so the error has
    a second job beyond refusing: it must send the caller to the path that
    works, rather than implying no search exists.
    """
    with pytest.raises(NotImplementedError, match="in-process") as excinfo:
        search.run_folddisco_search(Path("q.tsv"), Path("idx"), {})
    assert "folddisco_search.sh" in str(excinfo.value)
    with pytest.raises(NotImplementedError):
        search.parse_folddisco_output(Path("out.tsv"))


def test_reference_index_stage_is_blocked_on_dataset_selection():
    with pytest.raises(NotImplementedError, match="dataset selection"):
        reference_index.build_reference_index(Path("a"), Path("b"), {})
    with pytest.raises(NotImplementedError):
        reference_index.ingest_reference_metadata(Path("meta"))


def test_feature_and_scoring_stages_are_placeholders():
    with pytest.raises(NotImplementedError):
        features.compute_patch_features(None, None, {})
    with pytest.raises(NotImplementedError):
        scoring.aggregate_patch_evidence(None, {})
    with pytest.raises(NotImplementedError):
        benchmark.evaluate(None, None, {})


def test_scoring_module_documents_forbidden_label_names():
    """The naming prohibition must be written down where it is needed."""
    doc = scoring.__doc__ or ""
    for forbidden in (
        "probability_of_immunogenicity", "ada_risk", "antigenicity_score",
        "cross_reactivity",
    ):
        assert forbidden in doc, f"{forbidden} should be listed as forbidden"
    assert "does not establish" in doc


def test_search_module_records_the_verified_invocation_shape():
    doc = search.__doc__ or ""
    assert "folddisco query" in doc
    assert "matching_residues" in doc
    assert "parse_hit_residues" in doc
    # The panic hazard is real and must be documented for the next stage.
    assert "panic" in doc.lower()


def test_no_speculative_implementation_exists():
    """No extension module may ship a working implementation by accident."""
    import inspect

    for module in (search, reference_index, features, scoring, benchmark):
        source = inspect.getsource(module)
        assert "NotImplementedError" in source
        # No stage may quietly import a search/dataset dependency yet.
        assert "subprocess" not in source, module.__name__
