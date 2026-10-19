"""Extension point: aggregate hits into evidence. NOT IMPLEMENTED.

Naming rules, which are not negotiable
--------------------------------------
Any quantity produced here must be named for what it actually measures. A
score computed from FoldDisco hits against a collection of observed epitopes
is a **structural-motif-similarity** statistic. Permitted names are of the
form ``motif_similarity_score`` or
``n_hits_vs_observed_bcell_epitope_reference``.

Forbidden, because each asserts something unmeasured:

* ``probability_of_immunogenicity``
* ``ada_risk`` / ``ada_probability``
* ``antigenicity_score``
* ``epitope_score`` (a hit is similarity to an epitope *reference*, not an
  epitope)
* ``cross_reactivity``

A FoldDisco hit establishes that two residue constellations are geometrically
similar under a hash-based tolerance. It does not establish that an antibody
binds either one, that any circulating antibody recognises the query, or that
the query would induce a response.

Also required of this stage
---------------------------
* A per-patch score needs a stated null model. ``idf`` is a retrieval
  statistic computed from the reference index's hash frequencies; it is not
  calibrated and is not a probability.
* Aggregating several patches must not imply that separate matches jointly
  reconstitute one compatible antibody interface. They are independent
  geometric matches, possibly mutually inaccessible, and the patch radii that
  produced them are an exploration grid, not an interface model.
* Multiple-testing control is needed: this component emits hundreds of
  overlapping patches per structure (545 for a 239-residue, 3-chain input
  across five radii), and they are not independent.
"""

from __future__ import annotations

from typing import Any

_NOT_IMPLEMENTED = (
    "scoring and aggregation are out of scope for this component; see the "
    "module docstring for the naming rules any implementation must follow"
)


def aggregate_patch_evidence(hits: Any, config: dict[str, Any]) -> None:
    """Aggregate hits to patch-level evidence. NOT IMPLEMENTED."""
    raise NotImplementedError(_NOT_IMPLEMENTED)
