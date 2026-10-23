"""Configuration loading, overriding and validation.

The shipped ``ibc_search/data/default.toml`` is the single source of truth for
defaults
and carries a unit, a provenance note and a status marker per parameter. A user
config is merged over it key-by-key, and CLI overrides are merged last. The
fully resolved mapping is written verbatim into the run manifest, so a run can
always be reproduced from its own output.
"""

from __future__ import annotations

import copy
import tomllib
from importlib import resources
from pathlib import Path
from typing import Any

from .errors import ConfigError

_DEFAULT_IN_TREE = Path(__file__).resolve().parent / "data" / "default.toml"

VALID_CENTER_MODES = ("cb", "sidechain_centroid")
VALID_BURIED_POLICIES = ("drop", "flag")
VALID_DISCONNECTED_POLICIES = ("drop", "flag", "keep")
VALID_TABLE_FORMATS = ("tsv", "parquet")
VALID_EXPOSURE_CONTEXTS = ("scanned_chains", "all_polymer_chains")
VALID_ALTLOC_POLICIES = ("highest_occupancy",)
VALID_FETCH_FORMATS = ("mmcif", "pdb")


def default_config_path() -> Path:
    """Locate the shipped default config.

    The file lives *inside* the package (``ibc_search/data/default.toml``) so
    that it ships with a normal install; a copy at the repository root would
    not be installed and would leave the console script unable to start.
    """
    if _DEFAULT_IN_TREE.is_file():
        return _DEFAULT_IN_TREE
    packaged = resources.files("ibc_search") / "data" / "default.toml"
    if packaged.is_file():
        return Path(str(packaged))
    raise ConfigError(
        "cannot locate the shipped default config "
        "(ibc_search/data/default.toml)"
    )


def _deep_merge(base: dict[str, Any], over: dict[str, Any]) -> dict[str, Any]:
    """Recursively merge *over* into a copy of *base*."""
    out = copy.deepcopy(base)
    for key, value in over.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = _deep_merge(out[key], value)
        else:
            out[key] = copy.deepcopy(value)
    return out


def _coerce_scalar(text: str) -> Any:
    """Interpret a CLI override value with TOML scalar rules."""
    try:
        parsed = tomllib.loads(f"v = {text}")
        return parsed["v"]
    except tomllib.TOMLDecodeError:
        return text


def load_config(
    user_path: Path | None = None,
    overrides: dict[str, str] | None = None,
) -> dict[str, Any]:
    """Load defaults, merge a user config and dotted CLI overrides, validate.

    *overrides* maps a dotted path (``"ball_search.center_mode"``) to a string
    value, which is interpreted with TOML scalar rules so that
    ``--set ball_search.radii_angstrom='[8.0]'`` works.
    """
    with default_config_path().open("rb") as fh:
        config = tomllib.load(fh)

    if user_path is not None:
        if not user_path.is_file():
            raise ConfigError(f"config file not found: {user_path}")
        with user_path.open("rb") as fh:
            config = _deep_merge(config, tomllib.load(fh))

    for dotted, raw in (overrides or {}).items():
        parts = dotted.split(".")
        if len(parts) < 2:
            raise ConfigError(
                f"override {dotted!r} must be a dotted path like "
                "'section.key'"
            )
        cursor: dict[str, Any] = config
        for part in parts[:-1]:
            nxt = cursor.get(part)
            if not isinstance(nxt, dict):
                raise ConfigError(f"override {dotted!r}: no such config section")
            cursor = nxt
        if parts[-1] not in cursor:
            raise ConfigError(
                f"override {dotted!r}: unknown key (refusing to invent "
                "configuration that has no documented default)"
            )
        cursor[parts[-1]] = _coerce_scalar(raw)

    validate_config(config)
    return config


def validate_config(config: dict[str, Any]) -> None:
    """Fail loudly on an inconsistent configuration."""

    def req(section: str, key: str) -> Any:
        if section not in config or key not in config[section]:
            raise ConfigError(f"missing required config key [{section}].{key}")
        return config[section][key]

    def enum(section: str, key: str, allowed: tuple[str, ...]) -> None:
        value = req(section, key)
        if value not in allowed:
            raise ConfigError(
                f"[{section}].{key} = {value!r} is not one of {allowed}"
            )

    enum("run", "table_format", VALID_TABLE_FORMATS)
    enum("input", "altloc_policy", VALID_ALTLOC_POLICIES)
    enum("input", "fetch_format", VALID_FETCH_FORMATS)
    enum("exposure", "context", VALID_EXPOSURE_CONTEXTS)
    enum("ball_search", "center_mode", VALID_CENTER_MODES)
    enum("ball_search", "buried_member_policy", VALID_BURIED_POLICIES)
    enum("connectivity", "disconnected_policy", VALID_DISCONNECTED_POLICIES)

    probe = float(req("exposure", "probe_radius_angstrom"))
    if not 0.0 < probe <= 5.0:
        raise ConfigError(
            f"[exposure].probe_radius_angstrom = {probe} is outside the "
            "plausible range (0, 5] A"
        )

    radii = req("ball_search", "radii_angstrom")
    if not isinstance(radii, list) or not radii:
        raise ConfigError("[ball_search].radii_angstrom must be a non-empty list")
    for radius in radii:
        if not isinstance(radius, (int, float)) or radius <= 0:
            raise ConfigError(
                f"[ball_search].radii_angstrom contains invalid radius {radius!r}"
            )

    seed_rsa = float(req("ball_search", "seed_rsa_threshold"))
    member_rsa = float(req("ball_search", "member_rsa_threshold"))
    for name, value in (("seed_rsa_threshold", seed_rsa),
                        ("member_rsa_threshold", member_rsa)):
        if not 0.0 <= value <= 1.0:
            raise ConfigError(
                f"[ball_search].{name} = {value} must be a fraction in [0, 1]"
            )
    if member_rsa > seed_rsa:
        raise ConfigError(
            f"[ball_search].member_rsa_threshold ({member_rsa}) exceeds "
            f"seed_rsa_threshold ({seed_rsa}): every seed would then be "
            "ineligible as a member of its own patch"
        )

    min_size = int(req("ball_search", "min_patch_size"))
    max_size = int(req("ball_search", "max_patch_size"))
    # Imported here to keep the verified floor in one place.
    from .folddisco_contract import MIN_QUERY_RESIDUES

    if min_size < MIN_QUERY_RESIDUES:
        raise ConfigError(
            f"[ball_search].min_patch_size = {min_size} is below the verified "
            f"FoldDisco floor of {MIN_QUERY_RESIDUES} residues (a 1-residue "
            "query returns zero rows from the installed build)"
        )
    if max_size < min_size:
        raise ConfigError(
            f"[ball_search].max_patch_size ({max_size}) is below "
            f"min_patch_size ({min_size})"
        )

    jaccard = float(req("dedup", "jaccard_threshold"))
    if not 0.0 < jaccard <= 1.0:
        raise ConfigError(
            f"[dedup].jaccard_threshold = {jaccard} must be in (0, 1]"
        )

    workers = int(req("performance", "workers"))
    cap = int(req("performance", "max_workers_cap"))
    if workers < 1:
        raise ConfigError(f"[performance].workers = {workers} must be >= 1")
    if workers > cap:
        raise ConfigError(
            f"[performance].workers ({workers}) exceeds max_workers_cap ({cap})"
        )


def scale_label(config: dict[str, Any], radius: float) -> str:
    """Descriptive scale label for a radius.

    Bookkeeping only: the label records which exploration scale a patch came
    from and asserts nothing about epitope dimensions.
    """
    labels = config["ball_search"].get("scale_labels", {})
    for key, value in labels.items():
        if abs(float(key) - float(radius)) < 1e-9:
            return str(value)
    return "unlabelled"


REPO_ROOT = Path(__file__).resolve().parents[2]
"""Repository root, which registry paths may be written relative to."""


def _resolve_path(value: str) -> Path:
    """Expand ``~`` and resolve a registry path against the repository root."""
    path = Path(value).expanduser()
    return path if path.is_absolute() else (REPO_ROOT / path)


def database_names(config: dict[str, Any]) -> list[str]:
    return sorted((config.get("folddisco", {}).get("databases") or {}))


def resolve_database(config: dict[str, Any], name: str | None = None) -> dict[str, Any]:
    """Look up a named reference database and resolve its paths.

    Returns the registry entry plus ``name``, ``index_path`` and
    ``idf_params_path`` as :class:`Path`. Raises :class:`ConfigError` naming
    the available databases rather than failing later on a missing file.

    The pairing matters: IDF measures rarity over an index's contents and the
    per-motif aggregate is taken over a ``--top``-truncated hit list, so
    parameters fitted against one database are meaningless against another.
    Resolving both from one entry is what keeps them together.
    """
    registry = (config.get("folddisco", {}) or {}).get("databases") or {}
    if not registry:
        raise ConfigError("no [folddisco.databases] are configured")

    chosen = name or config.get("folddisco", {}).get("default_database")
    if not chosen:
        raise ConfigError(
            "no database given and [folddisco].default_database is unset; "
            f"available: {', '.join(sorted(registry))}"
        )
    if chosen not in registry:
        raise ConfigError(
            f"unknown reference database {chosen!r}; available: "
            f"{', '.join(sorted(registry))}"
        )

    entry = dict(registry[chosen])
    entry["name"] = chosen
    entry["index_path"] = _resolve_path(entry.get("index_path", ""))
    params = entry.get("idf_params")
    entry["idf_params_path"] = _resolve_path(params) if params else None
    return entry
