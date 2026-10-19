"""Input parsing, residue identity and the "never silently" requirements."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from ibc_search.errors import FetchError, InputError, StructureParseError
from ibc_search.structure_io import (
    RESIDUE_CLASS_CARBOHYDRATE,
    RESIDUE_CLASS_ION,
    RESIDUE_CLASS_WATER,
    looks_like_pdb_id,
    parse_structure,
    resolve_input,
)


def _parse(path: Path, config: dict):
    return parse_structure(resolve_input(str(path), config), config)


# --- residue identity -------------------------------------------------------

def test_insertion_codes_preserved_and_not_renumbered(icode_pdb, base_config):
    parsed = _parse(icode_pdb, base_config)
    labels = [r.rid.label() for r in parsed.residues]
    assert labels == ["H10", "H10A", "H10B", "H20"]

    # Author numbering is untouched: three residues share number 10, and the
    # jump from 10 to 20 is preserved rather than compacted.
    assert [r.rid.auth_seq_id for r in parsed.residues] == [10, 10, 10, 20]
    assert [r.rid.ins_code for r in parsed.residues] == ["", "A", "B", ""]
    # Internal indices are contiguous and independent of author numbering.
    assert [r.rid.seq_index for r in parsed.residues] == [0, 1, 2, 3]


def test_original_to_internal_to_query_roundtrip(icode_pdb, base_config):
    parsed = _parse(icode_pdb, base_config)
    index = parsed.residue_index

    # original -> identity
    residue = index.by_author("H", 10, "A")
    assert residue.res_name == "SER"
    # identity -> internal -> identity
    assert index.by_seq_index(residue.seq_index) is residue

    # The insertion-code residues collide in the FoldDisco token grammar, and
    # that collision must be reported rather than silently resolved.
    collisions = index.author_key_collisions()
    assert ("H", 10) in collisions
    assert len(collisions[("H", 10)]) == 3


def test_numbering_gap_is_flagged(icode_pdb, base_config):
    parsed = _parse(icode_pdb, base_config)
    by_label = {r.rid.label(): r for r in parsed.residues}
    assert by_label["H20"].seq_gap_before is True
    assert any("gap in author numbering" in w for w in parsed.warnings)


# --- alternate conformers ---------------------------------------------------

def test_altloc_resolved_by_occupancy_then_by_id(altloc_pdb, base_config):
    parsed = _parse(altloc_pdb, base_config)
    assert len(parsed.residues) == 2
    res1, res2 = parsed.residues

    assert res1.had_altloc and res1.altloc_candidates == ("A", "B")
    # Residue 1: B has occupancy 0.7 against A's 0.3.
    assert res1.altloc_chosen == "B"
    # Residue 2: both 0.5, so the lexicographic tie-break picks A.
    assert res2.altloc_chosen == "A"

    # Exactly one conformer of each atom survives.
    names1 = [
        str(parsed.atom_names[i])
        for i in range(res1.atom_start, res1.atom_stop)
    ]
    assert names1.count("CB") == 1 and names1.count("OG") == 1

    # The chosen conformer's geometry is the one kept: residue 1 CB at z=1.0
    # (altloc B), residue 2 CB at z=0.0 (altloc A).
    def coord_of(rec, name):
        for i in range(rec.atom_start, rec.atom_stop):
            if str(parsed.atom_names[i]) == name:
                return parsed.coords[i]
        raise AssertionError(f"{name} not found")

    assert coord_of(res1, "CB")[2] == pytest.approx(1.0, abs=1e-3)
    assert coord_of(res2, "CB")[2] == pytest.approx(0.0, abs=1e-3)
    assert any("alternate conformers" in w for w in parsed.warnings)


# --- missing atoms ----------------------------------------------------------

def test_missing_atoms_flagged_not_fabricated(glycine_pdb, base_config):
    parsed = _parse(glycine_pdb, base_config)
    by_label = {r.rid.label(): r for r in parsed.residues}

    gly = by_label["A3"]
    assert gly.is_glycine and not gly.has_cb
    assert gly.incomplete_sidechain is False  # GLY expects zero side-chain atoms

    leu = by_label["A5"]
    assert leu.rid.res_name == "LEU" and not leu.has_cb
    assert leu.incomplete_sidechain is True
    assert "incomplete_sidechain:0of4" in leu.quality_flags()

    ca_only = by_label["A7"]
    assert set(ca_only.missing_backbone_atoms) == {"N", "C", "O"}
    assert "missing_backbone:N/C/O" in ca_only.quality_flags()
    assert any("missing >=1 backbone atom" in w for w in parsed.warnings)


# --- models and chains ------------------------------------------------------

def test_default_model_is_first_and_selection_is_honoured(
    multi_model_pdb, base_config
):
    parsed1 = _parse(multi_model_pdb, base_config)
    assert parsed1.model == 1
    assert parsed1.available_models == [1, 2]
    z1 = parsed1.coords[:, 2].mean()

    base_config["input"]["model"] = 2
    parsed2 = _parse(multi_model_pdb, base_config)
    assert parsed2.model == 2
    # Model 2 is the 10 A z-shifted copy.
    assert parsed2.coords[:, 2].mean() == pytest.approx(z1 + 10.0, abs=1e-3)


def test_missing_model_raises_clearly(multi_model_pdb, base_config):
    base_config["input"]["model"] = 7
    with pytest.raises(StructureParseError, match="model 7 not present"):
        _parse(multi_model_pdb, base_config)


def test_chain_selection_and_scanned_chain_record(two_chain_pdb, base_config):
    parsed_all = _parse(two_chain_pdb, base_config)
    assert parsed_all.scanned_chains == ["A", "B"]
    assert any("2 polymer chains" in w for w in parsed_all.warnings)

    base_config["input"]["chains"] = ["B"]
    parsed_b = _parse(two_chain_pdb, base_config)
    assert parsed_b.scanned_chains == ["B"]
    assert parsed_b.available_chains == ["A", "B"]
    # Chain A residues are still parsed, just not scanned.
    assert {r.rid.auth_chain_id for r in parsed_b.residues} == {"A", "B"}
    assert all(
        r.is_scanned_chain == (r.rid.auth_chain_id == "B") for r in parsed_b.residues
    )


def test_unknown_chain_selection_raises(two_chain_pdb, base_config):
    base_config["input"]["chains"] = ["Z"]
    with pytest.raises(InputError, match="not polymer chains"):
        _parse(two_chain_pdb, base_config)


# --- heteroatoms ------------------------------------------------------------

def test_waters_ligands_glycans_classified_and_flagged(
    waters_and_ligand_pdb, base_config
):
    parsed = _parse(waters_and_ligand_pdb, base_config)
    classes = {r.rid.res_name: r.residue_class for r in parsed.residues}
    assert classes["HOH"] == RESIDUE_CLASS_WATER
    assert classes["ZN"] == RESIDUE_CLASS_ION
    assert classes["NAG"] == RESIDUE_CLASS_CARBOHYDRATE

    # None of them is scanned, and each is reported.
    assert all(
        not r.is_scanned_chain
        for r in parsed.residues
        if r.residue_class != "polymer_standard_aa"
    )
    joined = " ".join(parsed.warnings)
    assert "water residue(s) present" in joined
    assert "carbohydrate/glycan" in joined
    assert "monoatomic ion" in joined


def test_biological_assembly_ambiguity_is_recorded(linear_chain_pdb, base_config):
    parsed = _parse(linear_chain_pdb, base_config)
    assert "asymmetric unit" in parsed.assembly_note
    assert "not applied" in parsed.assembly_note.lower()


# --- input resolution -------------------------------------------------------

def test_gzip_input_is_transparent(linear_chain_pdb, base_config, tmp_path):
    gz_path = tmp_path / "linear.pdb.gz"
    gz_path.write_bytes(gzip.compress(linear_chain_pdb.read_bytes()))

    plain = _parse(linear_chain_pdb, base_config)
    gzipped = _parse(gz_path, base_config)
    assert [r.rid.label() for r in plain.residues] == [
        r.rid.label() for r in gzipped.residues
    ]
    assert resolve_input(str(gz_path), base_config).gzipped is True


def test_pdb_id_detection():
    assert looks_like_pdb_id("4CHA")
    assert looks_like_pdb_id("1crn")
    # First character must be a digit, so a 4-letter word is not an entry ID.
    assert not looks_like_pdb_id("ABCD")
    assert not looks_like_pdb_id("4CH")
    assert not looks_like_pdb_id("4CHAX")


def test_local_file_wins_over_pdb_id_interpretation(tmp_path, base_config, linear_chain_pdb):
    # A file literally named like an entry ID must be read as a file.
    named = tmp_path / "1ABC"
    named.write_bytes(linear_chain_pdb.read_bytes())
    with pytest.raises(InputError, match="cannot determine structure format"):
        resolve_input(str(named), base_config)


def test_unknown_extension_raises(tmp_path, base_config):
    bad = tmp_path / "model.xyz"
    bad.write_text("nonsense\n")
    with pytest.raises(InputError, match="cannot determine structure format"):
        resolve_input(str(bad), base_config)


def test_nonexistent_input_raises(base_config):
    with pytest.raises(InputError, match="neither an existing file"):
        resolve_input("/no/such/path/model.pdb", base_config)


def test_extended_pdb_id_rejected_explicitly(base_config):
    with pytest.raises(InputError, match="not supported yet"):
        resolve_input("pdb_00001abc", base_config)


def test_malformed_structure_file_raises(tmp_path, base_config):
    bad = tmp_path / "broken.pdb"
    bad.write_text("this is not a PDB file at all\nATOM nonsense\n")
    with pytest.raises((StructureParseError, InputError)):
        _parse(bad, base_config)


def test_empty_structure_file_raises(tmp_path, base_config):
    empty = tmp_path / "empty.pdb"
    empty.write_text("END\n")
    with pytest.raises(StructureParseError):
        _parse(empty, base_config)


# --- PDB ID fetch path ------------------------------------------------------

def test_pdb_id_served_from_cache_without_network(
    tmp_path, base_config, monkeypatch, linear_chain_pdb
):
    """A pre-populated cache must be used without touching the network.

    The network call is replaced by one that fails loudly, so a cache hit is
    proven rather than assumed.
    """
    import urllib.request

    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    base_config["input"]["cache_dir"] = str(cache_dir)
    base_config["input"]["fetch_format"] = "pdb"
    cached = cache_dir / "1ABC.pdb.gz"
    cached.write_bytes(gzip.compress(linear_chain_pdb.read_bytes()))

    def explode(*args, **kwargs):
        raise AssertionError("network must not be used when the cache has the entry")

    monkeypatch.setattr(urllib.request, "urlopen", explode)

    source = resolve_input("1ABC", base_config)
    assert source.kind == "pdb_id"
    assert source.pdb_id == "1ABC"
    assert source.from_cache is True
    assert source.sha256 and source.size_bytes > 0

    parsed = parse_structure(source, base_config)
    assert len(parsed.residues) == 10
    assert parsed.source.structure_id == "1ABC"


def test_unreachable_network_raises_fetch_error_with_no_fallback(
    tmp_path, base_config, monkeypatch
):
    import urllib.error
    import urllib.request

    base_config["input"]["cache_dir"] = str(tmp_path / "empty_cache")

    def unreachable(*args, **kwargs):
        raise urllib.error.URLError("simulated: no route to host")

    monkeypatch.setattr(urllib.request, "urlopen", unreachable)

    with pytest.raises(FetchError) as excinfo:
        resolve_input("9XYZ", base_config)
    message = str(excinfo.value)
    assert "does not fall back" in message
    assert "9XYZ" in message or "9xyz" in message.lower()


def test_http_error_is_reported_clearly(tmp_path, base_config, monkeypatch):
    import urllib.error
    import urllib.request

    base_config["input"]["cache_dir"] = str(tmp_path / "empty_cache2")

    def not_found(*args, **kwargs):
        raise urllib.error.HTTPError(
            url="https://files.rcsb.org/download/9ZZZ.cif.gz",
            code=404, msg="Not Found", hdrs=None, fp=None,
        )

    monkeypatch.setattr(urllib.request, "urlopen", not_found)
    with pytest.raises(FetchError, match="HTTP 404"):
        resolve_input("9ZZZ", base_config)


@pytest.mark.network
def test_real_rcsb_fetch_records_provenance(tmp_path, base_config):
    """Live fetch. Deselect with -m 'not network'."""
    base_config["input"]["cache_dir"] = str(tmp_path / "live_cache")
    source = resolve_input("1CRN", base_config)
    assert source.file_format == "mmcif"
    assert source.fetch_url == "https://files.rcsb.org/download/1CRN.cif.gz"
    assert source.fetch_date_utc is not None
    assert len(source.sha256) == 64
    assert source.from_cache is False

    # Second resolution must hit the cache and record it.
    again = resolve_input("1CRN", base_config)
    assert again.from_cache is True
    assert again.sha256 == source.sha256
