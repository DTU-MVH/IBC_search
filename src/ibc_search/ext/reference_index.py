"""Extension point: curate and build reference collections. NOT IMPLEMENTED.

This stage must assemble the reference structures that extracted motifs are
searched against, and build a FoldDisco index over them.

Verified index-building shape
-----------------------------
``folddisco index -p <pdb_dir|foldcomp_db> -i <index_out> [-t N] [-r]``
``[--id pdb|uniprot|afdb|relpath|abspath] [-y default|pdb|trrosetta|ppf|3di]``
``[-d <dist_bins>] [-a <angle_bins>] [--multiple-bins 16-4,8-3]``

The ``--id`` choice determines the ``tid`` value that appears in hit output,
so it determines how a hit is traced back to a reference entry. ``relpath`` is
the default and is path-dependent, which makes it a poor choice for a
long-lived reference index.

Open decisions that BLOCK this stage
------------------------------------
* **Which reference collection.** Candidates discussed for the project are
  experimentally observed B-cell epitopes / antibody-antigen interfaces
  (e.g. IEDB-derived structural entries, SAbDab complexes). Nothing has been
  selected, downloaded or inspected. No dataset is assumed here.
* **What a reference "entry" is.** A whole antigen chain, an antigen chain
  cropped to its interface, or the epitope residues only. This changes what a
  hit means and is the single most consequential open decision: searching
  whole antigens answers "does this motif recur in proteins that happen to be
  antigens", which is a different question from "does this motif resemble an
  observed epitope".
* **Negative/background set.** Without one, a hit count has no baseline and no
  score computed from it is interpretable.
* **Redundancy control.** Antibody-antigen databases are heavily redundant
  (the same antigen bound by many antibodies). Un-deduplicated references will
  inflate any hit-frequency statistic.
* **Indexing parameters.** The hash type and bin counts must match between
  index and query.

Terminology constraint
----------------------
A reference collection of observed epitopes supports statements about
*similarity to observed B-cell epitopes*. It does not support statements about
antibody binding, cross-reactivity, pre-existing ADA reactivity or
immunogenicity.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

_NOT_IMPLEMENTED = (
    "reference index construction is out of scope for this component and is "
    "blocked on dataset selection; see the module docstring"
)


def build_reference_index(
    structure_dir: Path, index_out: Path, config: dict[str, Any]
) -> None:
    """Build a FoldDisco index over a reference collection. NOT IMPLEMENTED."""
    raise NotImplementedError(_NOT_IMPLEMENTED)


def ingest_reference_metadata(source: Path) -> None:
    """Ingest epitope/interface annotations. NOT IMPLEMENTED.

    Deliberately absent: no epitope label, contact label or buried-surface-area
    computation is implemented anywhere in this component, so no such column
    exists in any output table and none may be invented downstream.
    """
    raise NotImplementedError(_NOT_IMPLEMENTED)
