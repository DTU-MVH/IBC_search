"""Configuration loading, overriding and validation."""

from __future__ import annotations

from pathlib import Path

import pytest

from ibc_search.config import (
    default_config_path,
    load_config,
    scale_label,
    validate_config,
)
from ibc_search.errors import ConfigError


def test_shipped_defaults_load_and_validate():
    config = load_config()
    assert config["exposure"]["probe_radius_angstrom"] == 1.4
    assert config["ball_search"]["center_mode"] == "cb"
    assert config["ball_search"]["radii_angstrom"] == [5.0, 6.5, 8.0, 10.0]
    assert config["performance"]["workers"] == 1
    validate_config(config)


def test_every_parameter_line_carries_a_status_marker():
    """The config is also documentation: each parameter states its status."""
    text = default_config_path().read_text()
    markers = ("documented", "proposed", "verified", "blocked")

    missing: list[str] = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if "=" not in stripped or stripped.startswith("["):
            continue
        if "status:" not in stripped:
            missing.append(stripped[:60])
            continue
        assert any(f"status: {m}" in stripped for m in markers), stripped[:80]

    # scale_labels is a lookup table documented by the radii entry above it.
    assert missing == [] or all("scale_labels" in m for m in missing), missing


def test_user_config_is_merged_over_defaults(tmp_path: Path):
    user = tmp_path / "user.toml"
    user.write_text(
        "[ball_search]\ncenter_mode = 'sidechain_centroid'\nseed_rsa_threshold = 0.35\n"
    )
    config = load_config(user)
    assert config["ball_search"]["center_mode"] == "sidechain_centroid"
    assert config["ball_search"]["seed_rsa_threshold"] == 0.35
    # Untouched keys keep their default.
    assert config["exposure"]["probe_radius_angstrom"] == 1.4


def test_dotted_overrides_use_toml_scalar_rules():
    config = load_config(
        None,
        {
            "ball_search.radii_angstrom": "[6.0, 9.0]",
            "performance.workers": "4",
            "connectivity.use_seed_normal_filter": "true",
            "run.table_format": "'tsv'",
        },
    )
    assert config["ball_search"]["radii_angstrom"] == [6.0, 9.0]
    assert config["performance"]["workers"] == 4
    assert config["connectivity"]["use_seed_normal_filter"] is True
    assert config["run"]["table_format"] == "tsv"


def test_override_of_unknown_key_is_refused():
    """Refusing unknown keys prevents a silently ineffective override."""
    with pytest.raises(ConfigError, match="unknown key"):
        load_config(None, {"ball_search.not_a_real_knob": "1"})
    with pytest.raises(ConfigError, match="no such config section"):
        load_config(None, {"nosuchsection.key": "1"})


def test_missing_config_file_raises():
    with pytest.raises(ConfigError, match="config file not found"):
        load_config(Path("/no/such/config.toml"))


@pytest.mark.parametrize(
    ("section", "key", "value", "pattern"),
    [
        ("ball_search", "center_mode", "magic", "not one of"),
        ("connectivity", "disconnected_policy", "maybe", "not one of"),
        ("run", "table_format", "xlsx", "not one of"),
        ("exposure", "context", "everything", "not one of"),
        ("exposure", "probe_radius_angstrom", 0.0, "plausible range"),
        ("exposure", "probe_radius_angstrom", 50.0, "plausible range"),
        ("ball_search", "seed_rsa_threshold", 1.4, "fraction"),
        ("dedup", "jaccard_threshold", 0.0, "must be in"),
        ("performance", "workers", 0, "must be >= 1"),
    ],
)
def test_invalid_values_are_rejected(section, key, value, pattern):
    config = load_config()
    config[section][key] = value
    with pytest.raises(ConfigError, match=pattern):
        validate_config(config)


def test_member_threshold_above_seed_threshold_is_refused():
    """A seed must be eligible as a member of its own patch."""
    config = load_config()
    config["ball_search"]["seed_rsa_threshold"] = 0.10
    config["ball_search"]["member_rsa_threshold"] = 0.50
    with pytest.raises(ConfigError, match="every seed would then be ineligible"):
        validate_config(config)


def test_min_patch_size_below_verified_folddisco_floor_is_refused():
    config = load_config()
    config["ball_search"]["min_patch_size"] = 1
    with pytest.raises(ConfigError, match="verified FoldDisco floor"):
        validate_config(config)


def test_max_patch_size_below_min_is_refused():
    config = load_config()
    config["ball_search"]["min_patch_size"] = 10
    config["ball_search"]["max_patch_size"] = 5
    with pytest.raises(ConfigError, match="below"):
        validate_config(config)


def test_empty_radius_list_is_refused():
    config = load_config()
    config["ball_search"]["radii_angstrom"] = []
    with pytest.raises(ConfigError, match="non-empty list"):
        validate_config(config)


def test_scale_labels_cover_the_default_radii():
    config = load_config()
    for radius in config["ball_search"]["radii_angstrom"]:
        assert scale_label(config, radius) != "unlabelled"
    # An unlisted radius is labelled explicitly rather than mislabelled.
    assert scale_label(config, 7.3) == "unlabelled"


def test_epitope_reference_index_is_marked_blocked():
    """The *epitope* reference collection is still an open decision.

    `reference_index_path` was unblocked once the whole-PDB FoldDisco index
    was adopted as the search background, but that index answers a different
    question: structural similarity to any deposited structure, not
    resemblance to a known epitope. The curated collection remains unchosen,
    and keeping it as its own blocked key stops the PDB index from quietly
    standing in for it.
    """
    config = load_config()
    assert config["folddisco"]["epitope_reference_index_path"] == ""
    text = default_config_path().read_text()
    assert "status: blocked" in text


def test_search_defaults_are_single_threaded():
    """-t 1 is a correctness default, not a performance oversight.

    The installed build does not write stdout line-atomically; above one
    thread a fraction of hit rows come back torn by interleaved writes. See
    the comment in default.toml and docs/FOLDDISCO_CONTRACT.md.
    """
    config = load_config()
    assert config["folddisco"]["search_threads"] == 1
    # query_residues is the join key from a hit back to its patch.
    assert "query_residues" in config["folddisco"]["search_format_output"]
