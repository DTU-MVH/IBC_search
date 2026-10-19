"""Extension point: per-patch and per-residue features. NOT IMPLEMENTED.

The patch table already carries the geometric descriptors the extractor can
justify: size, segment/gap structure, spatial diameter, radius of gyration,
compactness ratio, per-member selecting atom and distance, connectivity
counters and quality flags.

A later stage may add physicochemical and structural descriptors, for example
amino-acid composition summaries, charge and polarity summaries, secondary
structure context, or per-residue curvature.

Constraints this stage must respect
-----------------------------------
* **Chemistry is a descriptor, never a filter.** Project notes mention
  polar/charged enrichment in epitopes. That must not become a hard filter at
  any point: epitope chemistry spans polar, charged, hydrophobic and aromatic
  residues, and filtering on composition would discard real surface and bias
  every downstream statistic.
* Buried surface area is **not** computed here. It requires a bound partner,
  which a single minibinder structure does not have.
* Any feature derived from an unvalidated model must be named for what it
  measures, not for what it is hoped to predict.
"""

from __future__ import annotations

from typing import Any

_NOT_IMPLEMENTED = (
    "the extended feature set is out of scope for this component; the "
    "geometric descriptors that are justified are already in patches.tsv"
)


def compute_patch_features(patch: Any, parsed: Any, config: dict[str, Any]) -> None:
    """Compute extended per-patch features. NOT IMPLEMENTED."""
    raise NotImplementedError(_NOT_IMPLEMENTED)
