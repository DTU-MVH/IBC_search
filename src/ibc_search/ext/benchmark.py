"""Extension point: benchmarking and comparison. NOT IMPLEMENTED.

Intended scope for a later stage: evaluate motif-similarity evidence against
labelled data, and compare against existing predictors (e.g. DiscoTope for
B-cell epitope prediction).

Open requirements before any benchmark number is meaningful
-----------------------------------------------------------
* A labelled evaluation set with a stated definition of a positive. "Epitope
  residue" is defined differently across sources (contact-distance based,
  buried-surface based, curated), and the choice changes every metric.
* A redundancy-aware split. Antibody-antigen data is highly redundant; random
  splits leak.
* A baseline. Surface-exposed residues are not uniformly distributed, so any
  predictor that merely prefers exposed, protruding residues will score above
  chance. Exposure-matched baselines are required.
* A comparison protocol that is fair to the other method: DiscoTope predicts
  per-residue B-cell epitope propensity from structure, which is a different
  output type from patch-level motif similarity, so the two are not directly
  comparable without an explicit, stated mapping.

Out of scope here and not implemented anywhere in this component: ProteinMPNN,
ESM-IF, MaSIF, docking, clinical risk prediction and T-cell epitope
prediction.
"""

from __future__ import annotations

from typing import Any

_NOT_IMPLEMENTED = (
    "benchmarking is out of scope for this component and is blocked on a "
    "labelled evaluation set and a stated positive-class definition"
)


def evaluate(predictions: Any, labels: Any, config: dict[str, Any]) -> None:
    """Evaluate predictions against labels. NOT IMPLEMENTED."""
    raise NotImplementedError(_NOT_IMPLEMENTED)
