"""Tests for the IDF z-score normalisation.

The properties worth pinning are the ones the measurements forced:

* the size trend is **multiplicative**, so it is removed in log space -- a
  linear detrend on raw IDF leaves residual size bias (measured: correlation
  with size -0.063 raw vs -0.003 log);
* the transform is monotonic within a size class, so ranking motifs of equal
  size is unaffected;
* parameters are index- and ``--top``-specific, so applying them elsewhere
  warns rather than silently emitting incomparable numbers.
"""

from __future__ import annotations

import csv
import json
import math
import random
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "idf_zscore.py"

HEADER = "tid\tidf\trmsd\tnode_count\tquery_residues\tmatching_residues"


def _hits(path: Path, motifs: list[tuple[str, float]]) -> Path:
    """Write a minimal hits file: one row per motif."""
    lines = [HEADER]
    for query, idf in motifs:
        lines.append(f"x.cif\t{idf}\t0.1\t2\t{query}\tA1,A2")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _synthetic(n_per_size: int = 60, seed: int = 0) -> list[tuple[str, float]]:
    """Motifs whose IDF scales multiplicatively with size, as real ones do."""
    rng = random.Random(seed)
    out = []
    for size in range(3, 9):
        for i in range(n_per_size):
            # geometric growth with size, lognormal scatter
            idf = (10 ** (3.0 + 0.13 * size)) * math.exp(rng.gauss(0, 0.8))
            out.append((",".join(f"A{i}{j}" for j in range(size)), idf))
    return out


def _run(args: list[str]):
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args], capture_output=True, text=True
    )


def _fit(tmp_path: Path, motifs, **kw) -> Path:
    hits = _hits(tmp_path / "hits.tsv", motifs)
    params = tmp_path / "p.json"
    extra = []
    for k, v in kw.items():
        extra += [f"--{k.replace('_', '-')}", str(v)]
    r = _run(["fit", str(hits), "-o", str(params), *extra])
    assert r.returncode == 0, r.stderr
    return params


def test_fit_recovers_a_positive_size_slope(tmp_path: Path):
    params = json.loads(_fit(tmp_path, _synthetic()).read_text())
    assert params["slope_b"] > 0, "IDF grows with motif size; slope must be positive"
    assert params["sd_pooled"] > 0
    assert params["n_motifs"] > 0


def test_zscore_removes_the_size_dependence(tmp_path: Path):
    """The point of the whole exercise."""
    motifs = _synthetic()
    params = _fit(tmp_path, motifs)
    hits = _hits(tmp_path / "h2.tsv", motifs)
    out = tmp_path / "z.tsv"
    r = _run(["apply", str(hits), "-p", str(params), "-o", str(out)])
    assert r.returncode == 0, r.stderr

    rows = list(csv.DictReader(out.open(), delimiter="\t"))
    assert rows
    by_size: dict[int, list[float]] = {}
    for row in rows:
        by_size.setdefault(int(row["motif_size"]), []).append(float(row["idf_z"]))

    # every stratum centred near zero, none left systematically high or low
    for size, zs in by_size.items():
        mean = sum(zs) / len(zs)
        assert abs(mean) < 0.25, f"size {size} mean z = {mean:.3f}, not detrended"


def test_ranking_is_preserved_within_a_size_class(tmp_path: Path):
    """log is monotonic: it reorders nothing among motifs of equal size."""
    motifs = _synthetic()
    params = _fit(tmp_path, motifs)
    hits = _hits(tmp_path / "h3.tsv", motifs)
    out = tmp_path / "z3.tsv"
    assert _run(["apply", str(hits), "-p", str(params), "-o", str(out)]).returncode == 0

    rows = [r for r in csv.DictReader(out.open(), delimiter="\t")
            if int(r["motif_size"]) == 5]
    assert len(rows) > 5
    by_idf = sorted(rows, key=lambda r: -float(r["idf_aggregate"]))
    by_z = sorted(rows, key=lambda r: -float(r["idf_z"]))
    assert [r["query_residues"] for r in by_idf] == [r["query_residues"] for r in by_z]


def test_applying_parameters_from_another_index_warns(tmp_path: Path):
    """IDF rarity is defined over the index contents; mixing them is invalid."""
    motifs = _synthetic()
    params = _fit(tmp_path, motifs, index_label="pdb_folddisco", top_label=100)
    hits = _hits(tmp_path / "h4.tsv", motifs)
    r = _run(["apply", str(hits), "-p", str(params), "-o", str(tmp_path / "z4.tsv"),
              "--index-label", "sabdab2_antigen", "--top-label", "100"])
    assert r.returncode == 0
    assert "not comparable" in r.stderr, r.stderr


def test_too_few_motifs_is_refused(tmp_path: Path):
    """Fitting a trend on a handful of motifs would be noise."""
    hits = _hits(tmp_path / "h6.tsv", [("A1,A2,A3", 1000.0), ("A1,A2,A3,A4", 2000.0)])
    r = _run(["fit", str(hits), "-o", str(tmp_path / "p6.json")])
    assert r.returncode != 0
    assert "too few" in r.stderr


def test_zero_and_negative_idf_are_dropped_not_logged(tmp_path: Path):
    """log10 of a non-positive IDF is undefined; such rows must not crash."""
    motifs = _synthetic(n_per_size=25)
    motifs += [("B1,B2,B3", 0.0), ("B1,B2,B3,B4", -5.0)]
    params = _fit(tmp_path, motifs)
    hits = _hits(tmp_path / "h7.tsv", motifs)
    out = tmp_path / "z7.tsv"
    r = _run(["apply", str(hits), "-p", str(params), "-o", str(out)])
    assert r.returncode == 0, r.stderr
    queries = {row["query_residues"] for row in csv.DictReader(out.open(), delimiter="\t")}
    assert "B1,B2,B3" not in queries and "B1,B2,B3,B4" not in queries


# --- database registry ----------------------------------------------------


def test_registry_pairs_an_index_with_its_own_parameters():
    """An index and its IDF parameters must travel together.

    IDF measures rarity over the index contents, so a fit against one
    database is meaningless against another -- measured: the intercept is
    3.17 for the 230655-structure PDB and 0.99 for the 5951-structure antigen
    set. Selecting by name is what stops the two being mixed by hand.
    """
    from ibc_search.config import database_names, load_config, resolve_database

    config = load_config()
    names = database_names(config)
    assert "pdb" in names

    for name in names:
        entry = resolve_database(config, name)
        assert entry["name"] == name
        assert entry["index_path"].is_absolute()
        assert entry["idf_params_path"] is not None

    # the default resolves to a real registry entry
    assert resolve_database(config)["name"] in names


def test_unknown_database_names_the_alternatives():
    from ibc_search.config import load_config, resolve_database
    from ibc_search.errors import ConfigError

    with pytest.raises(ConfigError, match="available"):
        resolve_database(load_config(), "no_such_database")


def test_max_aggregate_does_not_warn_on_a_different_top(tmp_path: Path):
    """`--top` ranks by IDF, so the maximum survives any cap.

    Verified on all 883 reference motifs: max IDF was identical at top 40 and
    top 100. Warning about a cap mismatch would be noise for aggregate=max.
    """
    motifs = _synthetic()
    params = _fit(tmp_path, motifs, index_label="pdb", top_label=40)
    hits = _hits(tmp_path / "h8.tsv", motifs)
    r = _run(["apply", str(hits), "-p", str(params), "-o", str(tmp_path / "z8.tsv"),
              "--index-label", "pdb", "--top-label", "100"])
    assert r.returncode == 0
    assert "truncated" not in r.stderr and "depends on where" not in r.stderr


def test_cap_dependent_aggregate_still_warns(tmp_path: Path):
    motifs = _synthetic()
    params = _fit(tmp_path, motifs, index_label="pdb", top_label=40,
                  aggregate="mean_topk")
    hits = _hits(tmp_path / "h9.tsv", motifs)
    r = _run(["apply", str(hits), "-p", str(params), "-o", str(tmp_path / "z9.tsv"),
              "--index-label", "pdb", "--top-label", "100"])
    assert r.returncode == 0
    assert "depends on where" in r.stderr, r.stderr


# --- min-max scaling ------------------------------------------------------


def _minmax(tmp_path: Path, motifs, params: Path, space: str = "raw") -> list[dict]:
    hits = _hits(tmp_path / f"mm_{space}.tsv", motifs)
    out = tmp_path / f"mmout_{space}.tsv"
    r = _run(["minmax", str(hits), "-p", str(params), "-o", str(out),
              "--space", space])
    assert r.returncode == 0, r.stderr
    return list(csv.DictReader(out.open(), delimiter="\t"))


def test_minmax_output_is_strictly_within_unit_interval(tmp_path: Path):
    """The whole point of the method is the [0,1] guarantee."""
    motifs = _synthetic()
    params = _fit(tmp_path, motifs)
    for space in ("raw", "log"):
        rows = _minmax(tmp_path, motifs, params, space)
        vals = [float(r["idf_minmax"]) for r in rows]
        assert vals, space
        assert min(vals) >= 0.0 and max(vals) <= 1.0, (space, min(vals), max(vals))


def test_minmax_preserves_order_within_a_size_class(tmp_path: Path):
    """It is a linear rescale, so it reorders nothing of equal size."""
    motifs = _synthetic()
    params = _fit(tmp_path, motifs)
    rows = [r for r in _minmax(tmp_path, motifs, params) if int(r["motif_size"]) == 6]
    assert len(rows) > 5
    by_idf = [r["query_residues"] for r in
              sorted(rows, key=lambda r: -float(r["idf_aggregate"]))]
    by_mm = [r["query_residues"] for r in
             sorted(rows, key=lambda r: -float(r["idf_minmax"]))]
    assert by_idf == by_mm


def test_minmax_scales_per_size_not_globally(tmp_path: Path):
    """IDF grows with size; one global range would just rank motifs by size."""
    motifs = _synthetic()
    params = _fit(tmp_path, motifs)
    rows = _minmax(tmp_path, motifs, params)
    by_size: dict[int, list[float]] = {}
    for r in rows:
        by_size.setdefault(int(r["motif_size"]), []).append(float(r["idf_minmax"]))
    means = {s: sum(v) / len(v) for s, v in by_size.items() if len(v) >= 20}
    assert len(means) >= 3
    # no stratum systematically at the top or bottom of the scale
    assert max(means.values()) - min(means.values()) < 0.35, means


def test_minmax_clamps_values_beyond_the_reference_range(tmp_path: Path):
    """A new protein may exceed the reference maximum; [0,1] must still hold."""
    motifs = _synthetic()
    params = _fit(tmp_path, motifs)
    beyond = [("Z1,Z2,Z3", 1e12), ("Z1,Z2,Z3,Z4", 1e-6)]
    hits = _hits(tmp_path / "mm_out.tsv", motifs + beyond)
    out = tmp_path / "mm_clamped.tsv"
    r = _run(["minmax", str(hits), "-p", str(params), "-o", str(out)])
    assert r.returncode == 0, r.stderr
    assert "clamped" in r.stdout, r.stdout
    vals = [float(x["idf_minmax"]) for x in csv.DictReader(out.open(), delimiter="\t")]
    assert min(vals) >= 0.0 and max(vals) <= 1.0


def test_raw_and_log_spaces_differ_in_spread(tmp_path: Path):
    """Raw keeps the right skew; log spreads the mid-range."""
    motifs = _synthetic()
    params = _fit(tmp_path, motifs)
    raw = sorted(float(r["idf_minmax"]) for r in _minmax(tmp_path, motifs, params, "raw"))
    log = sorted(float(r["idf_minmax"]) for r in _minmax(tmp_path, motifs, params, "log"))
    med_raw, med_log = raw[len(raw) // 2], log[len(log) // 2]
    assert med_log > med_raw, (med_raw, med_log)


def test_minmax_on_parameters_without_ranges_is_refused(tmp_path: Path):
    """Parameters fitted before min-max support carry no ranges."""
    motifs = _synthetic()
    params = _fit(tmp_path, motifs)
    data = json.loads(params.read_text())
    data.pop("range_by_size", None)
    data.pop("range_global", None)
    stale = tmp_path / "stale.json"
    stale.write_text(json.dumps(data))
    hits = _hits(tmp_path / "mm2.tsv", motifs)
    r = _run(["minmax", str(hits), "-p", str(stale), "-o", str(tmp_path / "o.tsv")])
    assert r.returncode != 0
    assert "re-run `fit`" in r.stderr, r.stderr
