"""Extension points for later pipeline stages.

Every module here is a **documented placeholder**. Each states the inputs it
will consume, the contract it must honour, and the open decisions that block
it. None contains a speculative implementation: writing one now would mean
inventing FoldDisco behaviour, dataset semantics or biological labels that have
not been established.

Stage order, and where each stage attaches:

1. :mod:`ibc_search.ext.reference_index` -- build/curate the reference
   collections to search against. Blocked on dataset selection.
2. :mod:`ibc_search.ext.search` -- run FoldDisco over the emitted query file
   and parse hits back to residue identities. Unblocked technically; depends
   on stage 1 for anything meaningful to search.
3. :mod:`ibc_search.ext.features` -- per-patch and per-residue features.
4. :mod:`ibc_search.ext.scoring` -- aggregate hits into interpretable
   evidence. Must not produce a label implying antigenicity.
5. :mod:`ibc_search.ext.benchmark` -- evaluate against labelled data.

Terminology discipline applies across all of them. These are five distinct,
non-interchangeable things and no column, variable or report may blur them:

* structural motif similarity;
* similarity to experimentally observed B-cell epitopes;
* antibody binding / cross-reactivity;
* pre-existing anti-drug-antibody (ADA) reactivity;
* induction of a new immune response;
* T-cell peptide presentation.
"""

from __future__ import annotations

__all__ = [
    "benchmark",
    "features",
    "reference_index",
    "scoring",
    "search",
]
