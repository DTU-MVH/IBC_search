"""Extension point: run FoldDisco searches *from Python*. NOT IMPLEMENTED.

A search can already be run -- by ``scripts/folddisco_search.sh``, reachable
as ``ibc-search folddisco-search <results_dir>``. That script invokes the
FoldDisco binary, drops rows torn by the build's non-atomic multithreaded
writes, and joins each hit back to its ``patch_id``. Prefer it. What remains
unimplemented here is an *in-process* equivalent returning typed records, and
nothing in the pipeline needs one yet.

If that changes, note that the shell script is the reference implementation of
the invocation, and anything written here must reproduce its three
hard-won behaviours (all documented in docs/FOLDDISCO_CONTRACT.md):

* default to ``-t 1`` -- above one thread this build interleaves concurrent
  writes and tears roughly 0.5% of hit rows, nondeterministically;
* never join hits to queries by output order, which is completion order once
  threads are involved; join on the echoed ``query_residues`` instead;
* treat a query string as possibly mapping to several patches, since the same
  member set can be extracted at two radii.

What is already fixed, and what this stage must honour, is the contract in
:mod:`ibc_search.folddisco_contract`, all of it verified against the installed
build (commit ``2a756d95...``):

Invocation shape (verified)
---------------------------
``folddisco query -i <index> -q <query_file> [-t N] [--header] [--per-match]``

* ``-q`` accepts either an inline token string or a path to a file of
  ``"<structure_path>\\t<tokens>"`` lines -- which is exactly what
  :func:`ibc_search.outputs.write_query_file` writes. When ``-q`` is a file,
  ``-p`` is not required.
* Output columns are chosen with ``--format-output``; the per-match keys
  relevant here are ``tid``, ``node_count``, ``idf``, ``rmsd``,
  ``matching_residues``, ``query_residues`` (and optional ``tm_score``,
  ``gdt_ts``, ``gdt_ha``, ``chamfer_distance``, ``hausdorff_distance``).
* With several queries in one file, FoldDisco prints a header block **per
  query**, so a parser must not assume a single header line.

Reading hits back
-----------------
Use :func:`ibc_search.folddisco_contract.parse_hit_residues`, **not** naive
positional zipping. Measured over 4871 real hit rows, ``matching_residues`` is
equal in length to ``query_residues`` for only ~62% of rows: unmatched
positions are sometimes padded with ``_`` and sometimes omitted entirely,
including leading ones. ``align_hit_to_query`` is available for the alignable
rows and refuses the rest rather than inventing a correspondence; for a
trustworthy per-residue mapping request ``--superpose`` /
``matching_coordinates`` instead.

The ``query_residues`` column *is* reliable -- it echoes the emitted tokens
exactly -- so the query side of the round trip is sound.

Joining a hit to the patch that produced it goes through
``folddisco_queries.index.tsv``, because the query file format itself has no
field for a patch ID.

Known hazards to handle when this is written
--------------------------------------------
* The CLI **panics** (Rust panic, non-zero exit, no structured error) on a
  malformed token. Queries emitted by this package are pre-validated, so a
  panic indicates a contract regression and should be surfaced, not retried.
* ``--sampling-ratio`` / ``--sampling-count`` trade sensitivity for speed on
  long queries. How retrieval quality varies with query length is
  **uncharacterised** and must be measured before any threshold is trusted.
* ``idf`` is an inverse-document-frequency style score over hash matches. It
  is a *retrieval* statistic. It is not a probability, not a similarity to an
  epitope, and not a measure of antibody recognition.

Open decisions blocking this stage
----------------------------------
* which reference collection to search (see :mod:`.reference_index`);
* which hashing parameters (``-y/-d/-a``) the reference index is built with,
  since a query must be searched against a compatibly built index;
* what tolerance values (``-d``/``-a`` at query time) are appropriate -- these
  are FoldDisco's own distance/angle tolerances and are a *fourth*, separate
  parameter from the SASA probe radius, the patch ball radius and the
  patch-internal contact cutoff.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

_NOT_IMPLEMENTED = (
    "No in-process FoldDisco search is implemented. Use "
    "`ibc-search folddisco-search <results_dir>`, or scripts/folddisco_search.sh "
    "directly, which runs the search and joins hits back to their patches. See "
    "the module docstring before writing a Python equivalent."
)


def run_folddisco_search(
    query_file: Path,
    index_path: Path,
    config: dict[str, Any],
) -> None:
    """Run a FoldDisco search over an emitted query file. NOT IMPLEMENTED."""
    raise NotImplementedError(_NOT_IMPLEMENTED)


def parse_folddisco_output(output_path: Path) -> None:
    """Parse FoldDisco per-match output into records. NOT IMPLEMENTED.

    Must handle the per-query header blocks and the ``_`` placeholder, and
    must route residue resolution through
    :func:`ibc_search.folddisco_contract.align_hit_to_query` rather than
    parsing tokens ad hoc.
    """
    raise NotImplementedError(_NOT_IMPLEMENTED)
