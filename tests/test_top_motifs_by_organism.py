"""Tests for the top-motif organism annotation.

The behaviours worth pinning are the ones that could silently corrupt a
fraction computed over these rows: that a missing SIFTS mapping is reported
as unknown rather than assumed, that chimeric chains are not collapsed to a
single organism, and that the chain actually matched drives the lookup.
"""

from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "top_motifs_by_organism.py"
sys.path.insert(0, str(SCRIPT.parent))

from top_motifs_by_organism import (  # noqa: E402
    matched_chains,
    pdb_id_from_tid,
    resolve_organism,
)


def test_pdb_id_is_extracted_from_a_folddisco_path():
    assert pdb_id_from_tid("pdb/mw/pdb2mwl.ent.gz") == "2mwl"
    assert pdb_id_from_tid("/abs/path/pdb1abc.ent.gz") == "1abc"
    assert pdb_id_from_tid("antigens/pdb_00009fgy_C_+.cif") is None


def test_matched_chains_ignores_unmatched_placeholders():
    """`_` marks a query position the hit did not match."""
    assert matched_chains("A14,_,A16,B20") == ["A", "B"]
    assert matched_chains("") == []
    assert matched_chains("_,_") == []


def test_a_missing_mapping_is_unknown_not_absent():
    """SIFTS covers only chains it could map to UniProt.

    Reporting an unmapped chain as having no organism would quietly dilute
    any foreign/human fraction computed over these rows.
    """
    ids, names, source = resolve_organism("9zzz", ["A"], {}, {})
    assert (ids, names, source) == ("", "", "unknown")
    assert resolve_organism(None, ["A"], {}, {})[2] == "unknown"


def test_the_matched_chain_decides_the_organism():
    """In a multi-species complex the matched chain is the meaningful one."""
    by_chain = {("1abc", "A"): [("9606", "Homo sapiens")],
                ("1abc", "B"): [("10665", "Tequatrovirus T4")]}
    by_entry = {"1abc": {("9606", "Homo sapiens"), ("10665", "Tequatrovirus T4")}}
    assert resolve_organism("1abc", ["B"], by_chain, by_entry)[1] == "Tequatrovirus T4"
    assert resolve_organism("1abc", ["A"], by_chain, by_entry)[1] == "Homo sapiens"


def test_a_chimeric_chain_keeps_every_organism():
    by_chain = {("1abc", "A"): [("9606", "Homo sapiens"), ("10090", "Mus musculus")]}
    ids, names, source = resolve_organism("1abc", ["A"], by_chain, {})
    assert names == "Homo sapiens|Mus musculus" and ids == "9606|10090"
    assert source == "chain"


def test_entry_fallback_only_when_the_entry_is_unambiguous():
    """Falling back on a mixed-species entry would invent an answer."""
    single = {"1abc": {("9606", "Homo sapiens")}}
    assert resolve_organism("1abc", ["Z"], {}, single)[2] == "entry"
    mixed = {"1abc": {("9606", "Homo sapiens"), ("10090", "Mus musculus")}}
    assert resolve_organism("1abc", ["Z"], {}, mixed)[2] == "unknown"


def _fixture(tmp_path: Path) -> tuple[Path, Path]:
    d = tmp_path / "1XYZ"
    d.mkdir()
    (d / "idf_scored.tsv").write_text(
        "query_residues\tmotif_size\tidf_z\n"
        "A1,A2,A3\t3\t+3.0\nB1,B2,B3\t3\t+2.0\n"
        "C1,C2,C3\t3\t+1.0\nD1,D2,D3\t3\t-1.0\n"
    )
    (d / "folddisco_hits.tsv").write_text(
        "patch_id\tradius_angstrom\tscale\ttid\tidf\trmsd\tnode_count\t"
        "query_residues\tmatching_residues\n"
        "P0\t5.5\ts\tpdb/ab/pdb1abc.ent.gz\t100\t0.2\t3\tA1,A2,A3\tA7,A8,_\n"
        "P1\t5.5\ts\tpdb/de/pdb2def.ent.gz\t90\t0.3\t3\tD1,D2,D3\tA1,A2,A3\n"
    )
    tax = tmp_path / "tax.tsv"
    tax.write_text("# header line\nPDB\tCHAIN\tTAX_ID\tSCIENTIFIC_NAME\n"
                   "1abc\tA\t9606\tHomo sapiens\n")
    return d, tax


def test_only_the_top_fraction_is_written(tmp_path: Path):
    d, tax = _fixture(tmp_path)
    r = subprocess.run(
        [sys.executable, str(SCRIPT), str(d), "--taxonomy", str(tax),
         "--fraction", "0.25"],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr
    rows = list(csv.DictReader((d / "top_motifs_organism.tsv").open(), delimiter="\t"))
    # 4 motifs, top 25% = 1: the z=+3.0 one. The z=-1.0 hit must not appear.
    assert [x["query_residues"] for x in rows] == ["A1,A2,A3"]
    assert rows[0]["organism_names"] == "Homo sapiens"
    assert rows[0]["is_human"] == "1"
    assert rows[0]["motif_rank"] == "1"


def test_a_missing_taxonomy_file_is_refused(tmp_path: Path):
    d, _ = _fixture(tmp_path)
    r = subprocess.run(
        [sys.executable, str(SCRIPT), str(d), "--taxonomy", str(tmp_path / "nope.tsv")],
        capture_output=True, text=True,
    )
    assert r.returncode == 2
    assert "taxonomy file not found" in r.stderr


def test_rows_are_sorted_by_motif_rank_then_strongest_hit(tmp_path: Path):
    """The file should open on the motifs worth reading.

    FoldDisco's own output order is completion order and carries no meaning,
    so the rows are reordered: best motif first, and within a motif its
    highest-IDF hit first.
    """
    d = tmp_path / "2ABC"
    d.mkdir()
    (d / "idf_scored.tsv").write_text(
        "query_residues\tmotif_size\tidf_z\n"
        "LOW,L2,L3\t3\t+1.0\nTOP,T2,T3\t3\t+9.0\nMID,M2,M3\t3\t+5.0\n"
    )
    # deliberately written worst-first, and with hits out of IDF order
    (d / "folddisco_hits.tsv").write_text(
        "patch_id\tradius_angstrom\tscale\ttid\tidf\trmsd\tnode_count\t"
        "query_residues\tmatching_residues\n"
        "P\t5.5\ts\tpdb/ab/pdb1abc.ent.gz\t10\t0.5\t3\tLOW,L2,L3\tA1,A2,A3\n"
        "P\t5.5\ts\tpdb/ab/pdb1abc.ent.gz\t50\t0.4\t3\tMID,M2,M3\tA1,A2,A3\n"
        "P\t5.5\ts\tpdb/ab/pdb1abc.ent.gz\t20\t0.3\t3\tTOP,T2,T3\tA1,A2,A3\n"
        "P\t5.5\ts\tpdb/ab/pdb1abc.ent.gz\t99\t0.2\t3\tTOP,T2,T3\tA1,A2,A3\n"
    )
    tax = tmp_path / "tax.tsv"
    tax.write_text("PDB\tCHAIN\tTAX_ID\tSCIENTIFIC_NAME\n1abc\tA\t9606\tHomo sapiens\n")

    r = subprocess.run(
        [sys.executable, str(SCRIPT), str(d), "--taxonomy", str(tax),
         "--fraction", "1.0"],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr
    rows = list(csv.DictReader((d / "top_motifs_organism.tsv").open(), delimiter="\t"))

    ranks = [int(x["motif_rank"]) for x in rows]
    assert ranks == sorted(ranks), ranks
    assert rows[0]["query_residues"] == "TOP,T2,T3"
    # within the top motif, the stronger hit leads
    assert float(rows[0]["idf"]) == 99.0
    assert float(rows[1]["idf"]) == 20.0
    assert rows[-1]["query_residues"] == "LOW,L2,L3"
