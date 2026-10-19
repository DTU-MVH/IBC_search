"""Exception types and the MISSING sentinel.

A central rule of this package: a quantity that could not be computed is
represented as :data:`MISSING` carrying a human-readable reason, and never as
``0``. Zero is a legitimate SASA value for a fully buried residue, so using it
as "unknown" would silently corrupt downstream exposure statistics.
"""

from __future__ import annotations

from dataclasses import dataclass


class IbcSearchError(Exception):
    """Base class for all errors raised by this package."""


class InputError(IbcSearchError):
    """The requested input could not be interpreted or does not exist."""


class StructureParseError(IbcSearchError):
    """A structure file exists but could not be parsed into an atom array."""


class FetchError(IbcSearchError):
    """A PDB entry could not be fetched from RCSB.

    Raised instead of falling back to any other source. A silent fallback
    would make the run manifest's provenance record false.
    """


class ConfigError(IbcSearchError):
    """The configuration is malformed or internally inconsistent."""


class ContractError(IbcSearchError):
    """A motif cannot be expressed under the FoldDisco query contract."""


@dataclass(frozen=True)
class Missing:
    """A value that is unavailable, with a stated reason.

    Serialises to an empty table cell plus a populated ``*_missing_reason``
    column, so the reason survives into the output tables.
    """

    reason: str

    def __bool__(self) -> bool:  # pragma: no cover - guard against truthiness
        raise TypeError(
            "Missing has no truth value; test with `isinstance(x, Missing)` "
            "or `is_missing(x)` instead of a boolean check."
        )

    def __str__(self) -> str:
        return ""


def is_missing(value: object) -> bool:
    """Return True if *value* is a :class:`Missing` sentinel."""
    return isinstance(value, Missing)


MISSING = Missing
"""Alias so callers read ``MISSING("no reference value")`` at use sites."""
