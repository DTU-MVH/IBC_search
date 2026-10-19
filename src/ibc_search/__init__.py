"""ibc-search: discontinuous surface motif extraction for FoldDisco queries.

This package extracts candidate **structural motifs** from the solvent-exposed
surface of a protein structure and prepares them as FoldDisco queries.

Scope and interpretation
------------------------
An extracted patch is a *geometric neighbourhood of exposed residues*. It is
NOT an epitope, not a validated antibody footprint, and carries no claim about
antigenicity, antibody cross-reactivity, pre-existing anti-drug-antibody
reactivity, induction of a new immune response, or T-cell peptide
presentation. Those are separate, non-interchangeable concepts; see
docs/ASSUMPTIONS.md.
"""

from __future__ import annotations

__version__ = "0.1.0"

__all__ = ["__version__"]
