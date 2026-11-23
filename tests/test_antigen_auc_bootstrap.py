import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from antigen_auc_bootstrap import auc, bootstrap, counts, write_table


@pytest.mark.parametrize("a,b,expected", [([0], [0], .5), ([2], [1], 1), ([0], [1], 0),
                                         ([0, 2], [0, 1], .625)])
def test_auc(a, b, expected):
    assert auc(a, b) == expected


def test_empty_class_rejected():
    with pytest.raises(ValueError):
        auc([], [1])


def test_paired_bootstrap_equal_weight_reproducible():
    values = np.array([[0., .1], [1., 1.1]])
    result, draws = bootstrap(values, 100)
    again, again_draws = bootstrap(values, 100)
    assert np.array_equal(result, again)
    assert np.array_equal(draws, again_draws)
    assert np.allclose(result, values[draws].mean(axis=1))
    assert np.allclose(result[:, 1] - result[:, 0], .1)


def test_counts_zero_unassessed_dedup_and_shared_query(tmp_path):
    write_table(tmp_path / "patches.tsv", [dict(patch_id=p, emit_query=e, seed_class=c,
                 radius_angstrom=5.5, seed_residue_labels=p) for p,e,c in
                 [("a","true","annotated"),("b","true","background"),
                  ("zero","true","background"),("unassessed","false","annotated")]])
    write_table(tmp_path / "folddisco_queries.index.tsv", [dict(patch_id=p) for p in ("a","b","zero")])
    write_table(tmp_path / "folddisco_hits.tsv", [dict(patch_id=p,tid="target",idf=i,matching_residues=m)
                for p,i,m in [("a|b",70,"A1,_,A2"),("a",80,"A2,A1"),("a",50,"A3"),("zero",90,"_,_")]])
    result, skipped = counts(tmp_path)
    by_id = {r['patch_id']: r for r in result}
    assert skipped == 1 and set(by_id) == {"a","b","zero"}
    assert [by_id['a'][f'idf_{t}'] for t in (50,70,80)] == [2,1,1]
    assert [by_id['b'][f'idf_{t}'] for t in (50,70,80)] == [1,1,0]
    assert by_id['zero']['idf_50'] == 0


def test_missing_hits_not_inferred_as_zero(tmp_path):
    write_table(tmp_path / "patches.tsv", [dict(patch_id='a',emit_query='true')])
    write_table(tmp_path / "folddisco_queries.index.tsv", [dict(patch_id='a')])
    with pytest.raises(FileNotFoundError):
        counts(tmp_path)