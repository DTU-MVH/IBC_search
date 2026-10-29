"""Every TSV this project writes must use LF line endings.

Python's `csv.writer` defaults to CRLF. `csv.DictReader` consumes that
transparently, so the problem is invisible from Python and surfaces only
when a file is parsed with shell tools -- where the trailing `\\r` silently
joins the last field. It cost a whole 40-structure run: every path read out
of a generated panel with `awk` came back as `...cif\\r`, and all 40 scans
failed with "neither an existing file nor a PDB entry ID".

src/ibc_search/outputs.py always set `lineterminator="\\n"`. This test
exists so the scripts cannot drift from that again.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SOURCES = sorted(ROOT.glob("scripts/*.py")) + sorted(ROOT.glob("src/ibc_search/*.py"))


def _csv_writer_calls(path: Path):
    """Yield every `csv.writer(...)` call node in a source file."""
    tree = ast.parse(path.read_text())
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr == "writer"
                and isinstance(node.func.value, ast.Name)
                and node.func.value.id == "csv"):
            yield node


@pytest.mark.parametrize("path", SOURCES, ids=lambda p: p.name)
def test_csv_writers_set_lf_line_endings(path: Path):
    for call in _csv_writer_calls(path):
        kwargs = {k.arg: k.value for k in call.keywords}
        assert "lineterminator" in kwargs, (
            f"{path.name}:{call.lineno} csv.writer without lineterminator; "
            "it defaults to CRLF, which breaks awk/cut parsing"
        )
        value = kwargs["lineterminator"]
        assert isinstance(value, ast.Constant) and value.value == "\n", (
            f"{path.name}:{call.lineno} lineterminator must be '\\n'"
        )


def test_generated_tsvs_on_disk_have_no_carriage_returns():
    """Spot-check real output, not just the source that wrote it."""
    checked = 0
    for tsv in list(ROOT.glob("results/*/*.tsv")) + list(ROOT.glob("results/*/*/*.tsv")):
        head = tsv.open("rb").read(8192)
        assert b"\r\n" not in head, f"{tsv} has CRLF line endings"
        checked += 1
        if checked >= 40:
            break
    if checked == 0:
        pytest.skip("no generated TSVs present to check")
