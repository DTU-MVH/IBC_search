"""Tests for how the background set is selected.

Composition is this set's whole product, so the selection rules are the thing
worth pinning. Every test here is offline: the network parts are the easy
parts, and the parts that would silently produce a biased set are these.

The quality rule gets the most attention because it is the one that fails
invisibly: "resolution below 3 Angstroms" silently excludes every NMR
structure, which has no resolution at all, while still looking like a quality
filter.
"""

from __future__ import annotations

import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "build_background_dataset.py"
sys.path.insert(0, str(SCRIPT.parent))

from build_background_dataset import (  # noqa: E402
    SIZE_BINS,
    Selector,
    afdb_candidates,
    bin_of,
    bin_quotas,
    describe_entity,
    passes_quality,
)


def meta(methods, resolution, **kw):
    out = {"methods": list(methods), "resolution": resolution,
           "n_residues_reported": 200, "organism": "Homo sapiens"}
    out.update(kw)
    return out


# --- the quality rule, which means three different things -------------------

def test_xray_is_held_to_the_resolution_limit():
    assert passes_quality(meta(["X-RAY DIFFRACTION"], 2.4), 3.0)[0]
    assert not passes_quality(meta(["X-RAY DIFFRACTION"], 3.4), 3.0)[0]


def test_the_resolution_bound_is_strict():
    """'less than 3 A' excludes exactly 3.00, rather than rounding it in."""
    assert not passes_quality(meta(["X-RAY DIFFRACTION"], 3.0), 3.0)[0]
    assert passes_quality(meta(["X-RAY DIFFRACTION"], 2.99), 3.0)[0]


def test_nmr_is_admitted_although_it_has_no_resolution():
    """The failure this guards: NMR silently excluded by a resolution filter."""
    ok, why = passes_quality(meta(["SOLUTION NMR"], None), 3.0)
    assert ok, why


def test_a_missing_resolution_does_not_admit_a_non_nmr_entry():
    """Absent must not read as passing, or the NMR exemption becomes a hole."""
    assert not passes_quality(meta(["X-RAY DIFFRACTION"], None), 3.0)[0]


def test_cryo_em_is_held_to_the_limit_like_xray():
    assert passes_quality(meta(["ELECTRON MICROSCOPY"], 2.8), 3.0)[0]
    assert not passes_quality(meta(["ELECTRON MICROSCOPY"], 3.1), 3.0)[0]


def test_an_unlisted_method_is_rejected():
    assert not passes_quality(meta(["NEUTRON DIFFRACTION"], 2.0), 3.0)[0]


def test_a_multi_method_entry_passes_on_any_accepted_method():
    """X-ray plus neutron is an X-ray structure; the full list is kept."""
    ok, _ = passes_quality(meta(["X-RAY DIFFRACTION", "NEUTRON DIFFRACTION"], 2.0), 3.0)
    assert ok


# --- size strata ------------------------------------------------------------

def test_the_bins_cover_the_whole_requested_range_without_gaps():
    assert SIZE_BINS[0][0] == 50 and SIZE_BINS[-1][1] == 2300
    for (_, hi), (lo, _) in zip(SIZE_BINS, SIZE_BINS[1:]):
        assert hi == lo, "a gap between bins would drop sizes silently"


def test_every_size_in_range_lands_in_exactly_one_bin():
    for n in (50, 99, 100, 349, 350, 1199, 1200, 2299, 2300):
        assert bin_of(n) is not None, n
    assert bin_of(49) is None
    assert bin_of(2301) is None


def test_the_top_bin_includes_its_upper_bound():
    """2300 is in range, so the last bin has to be closed at both ends."""
    assert bin_of(2300) == len(SIZE_BINS) - 1


def test_quotas_sum_to_the_target():
    for total in (7, 21, 750, 751):
        q = bin_quotas(total)
        assert sum(q) == total
        assert max(q) - min(q) <= 1, "quotas must be even, not proportional"


# --- the selector -----------------------------------------------------------

def row(n_residues: int, organism: str = "Homo sapiens") -> dict:
    return {"size_bin": bin_of(n_residues), "organism": organism,
            "n_residues": n_residues}


def test_a_full_bin_stops_accepting():
    sel = Selector([1] * len(SIZE_BINS), organism_cap=99)
    assert sel.offer(row(120))
    assert not sel.offer(row(150)), "same bin, quota of 1 already used"
    assert sel.offer(row(600)), "a different bin is unaffected"


def test_the_organism_cap_stops_one_species_dominating():
    sel = Selector([99] * len(SIZE_BINS), organism_cap=2)
    assert sel.offer(row(120, "Escherichia coli"))
    assert sel.offer(row(130, "Escherichia coli"))
    assert not sel.offer(row(140, "Escherichia coli"))
    assert sel.offer(row(150, "Homo sapiens")), "the cap is per organism"


def test_unknown_organisms_are_not_capped_as_if_they_were_one_species():
    """'unknown' is an absence of a label, not a species they share."""
    sel = Selector([99] * len(SIZE_BINS), organism_cap=1)
    assert sel.offer(row(120, ""))
    assert sel.offer(row(130, ""))


def test_shortfall_names_the_bins_that_could_not_be_filled():
    sel = Selector([2] * len(SIZE_BINS), organism_cap=99)
    sel.offer(row(120))
    short = sel.shortfall()
    assert short["100-200"] == 1
    assert short["1200-2300"] == 2
    assert "50-100" in short


def test_relax_only_raises_quotas_and_never_drops_a_selection():
    sel = Selector([1] * len(SIZE_BINS), organism_cap=99)
    sel.offer(row(120))
    before = list(sel.taken)
    sel.relax(3)
    assert sel.taken == before
    assert all(q == 4 for q in sel.quotas)
    assert sel.offer(row(150)), "the relaxed bin accepts again"


def test_rejections_are_counted_rather_than_discarded():
    sel = Selector([1] * len(SIZE_BINS), organism_cap=99)
    sel.offer(row(120))
    sel.offer(row(150))
    assert sum(sel.rejected.values()) == 1


# --- the AlphaFold pool -----------------------------------------------------

def test_the_afdb_pool_is_filtered_on_size_and_plddt(tmp_path):
    lookup = tmp_path / "x.lookup"
    lookup.write_text(
        "0\tAF-P00001-F1-model_v4\t120\t91.0\t0\n"   # keep
        "1\tAF-P00002-F1-model_v4\t40\t95.0\t1\n"    # too small
        "2\tAF-P00003-F1-model_v4\t300\t55.0\t2\n"   # pLDDT too low
        "3\tAF-P00004-F1-model_v4\t3000\t99.0\t3\n"  # too large
        "4\tAF-P00005-F1-model_v4\t2300\t70.0\t4\n"  # both bounds inclusive
    )
    got = afdb_candidates(lookup, 50, 2300, 70.0)

    assert [a for a, _, _ in got] == ["P00001", "P00005"]
    assert got[0] == ("P00001", 120, 91.0)


def test_a_malformed_lookup_line_is_skipped_not_fatal(tmp_path):
    lookup = tmp_path / "x.lookup"
    lookup.write_text(
        "0\tAF-P00001-F1-model_v4\t120\t91.0\t0\n"
        "junk\n"
        "2\tAF-P00002-F1-model_v4\tnot-a-number\t91.0\t2\n"
        "3\tAF-P00003-F1-model_v4\t200\t88.0\t3\n"
    )
    assert [a for a, _, _ in afdb_candidates(lookup, 50, 2300, 70.0)] == \
        ["P00001", "P00003"]


# --- metadata flattening ----------------------------------------------------

def test_describe_entity_survives_the_fields_rcsb_leaves_null():
    """Every nested block here is one RCSB really does return as null."""
    entity = {
        "rcsb_id": "1ABC_1",
        "rcsb_polymer_entity": None,
        "rcsb_polymer_entity_container_identifiers": None,
        "rcsb_entity_source_organism": None,
        "entry": {"rcsb_id": "1ABC", "exptl": [{"method": "SOLUTION NMR"}],
                  "rcsb_entry_info": {"resolution_combined": None,
                                      "deposited_polymer_monomer_count": 88,
                                      "polymer_entity_count_protein": 1}},
    }
    got = describe_entity(entity)

    assert got["pdb_id"] == "1ABC"
    assert got["resolution"] is None
    assert got["organism"] == "" and got["uniprot_ids"] == []
    assert got["n_residues_reported"] == 88
    assert passes_quality(got, 3.0)[0], "an NMR entry must still get through"


def test_describe_entity_takes_the_first_of_several_resolutions():
    entity = {
        "rcsb_id": "1ABC_1",
        "rcsb_polymer_entity": {"pdbx_description": "a protein",
                                "rcsb_ec_lineage": [{"id": "3"}]},
        "rcsb_polymer_entity_container_identifiers": {"uniprot_ids": ["P12345"]},
        "rcsb_entity_source_organism": [{"ncbi_scientific_name": "Homo sapiens",
                                         "ncbi_taxonomy_id": 9606}],
        "entry": {"rcsb_id": "1ABC", "exptl": [{"method": "X-RAY DIFFRACTION"}],
                  "rcsb_entry_info": {"resolution_combined": [1.9, 2.4],
                                      "deposited_polymer_monomer_count": 250,
                                      "polymer_entity_count_protein": 2}},
    }
    got = describe_entity(entity)

    assert got["resolution"] == 1.9
    assert got["ec_class"] == "3"
    assert got["uniprot_ids"] == ["P12345"]
    assert got["tax_id"] == 9606


# --- redundancy at the level of what actually gets written ------------------

def test_entry_accessions_cover_every_entity_not_just_the_clustered_one():
    """The regression this guards, found in a real 750-structure build.

    RCSB clusters one ENTITY but the builder writes the whole ENTRY. Four
    thrombin-hirudin complexes were selected because their small hirudin
    peptides fell into different 30% clusters, while the files were
    near-identical. Judging redundancy on the entry's full accession set is
    what catches them.
    """
    entity = {
        "rcsb_id": "1D3Q_2",
        "rcsb_polymer_entity": {"pdbx_description": "HIRUGEN"},
        # The clustered entity is the small peptide...
        "rcsb_polymer_entity_container_identifiers": {"uniprot_ids": ["P01050"]},
        "rcsb_entity_source_organism": [{"ncbi_scientific_name": "Hirudo medicinalis"}],
        "entry": {
            "rcsb_id": "1D3Q",
            "exptl": [{"method": "X-RAY DIFFRACTION"}],
            "rcsb_entry_info": {"resolution_combined": [2.0],
                                "deposited_polymer_monomer_count": 290,
                                "polymer_entity_count_protein": 2},
            # ...but the file also holds thrombin, which is most of it.
            "polymer_entities": [
                {"rcsb_polymer_entity_container_identifiers":
                    {"uniprot_ids": ["P00734"]}},
                {"rcsb_polymer_entity_container_identifiers":
                    {"uniprot_ids": ["P01050"]}},
            ],
        },
    }
    got = describe_entity(entity)

    assert got["uniprot_ids"] == ["P01050"], "the clustered entity's own id"
    assert got["entry_uniprot_ids"] == ["P00734", "P01050"], \
        "the whole file's content is what redundancy must be judged on"


def test_entry_accessions_are_empty_rather_than_absent_when_unmapped():
    entity = {
        "rcsb_id": "1ABC_1",
        "rcsb_polymer_entity": None,
        "rcsb_polymer_entity_container_identifiers": None,
        "rcsb_entity_source_organism": None,
        "entry": {"rcsb_id": "1ABC", "exptl": [{"method": "SOLUTION NMR"}],
                  "rcsb_entry_info": {"resolution_combined": None,
                                      "deposited_polymer_monomer_count": 88,
                                      "polymer_entity_count_protein": 1},
                  "polymer_entities": [
                      {"rcsb_polymer_entity_container_identifiers": None}]},
    }
    assert describe_entity(entity)["entry_uniprot_ids"] == []
