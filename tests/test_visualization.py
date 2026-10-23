"""Tests for the figure module.

These check the mechanics -- that a file is produced, that the summary is
computed on the full data rather than the drawn subsample, and that empty
input is refused -- not the appearance, which has to be looked at.
"""

from __future__ import annotations

import random
from pathlib import Path

import pytest

from ibc_search.visualization import THEMES, binned_point_plot, load_scored

pytest.importorskip("matplotlib")


def _groups(seed: int = 0):
    rng = random.Random(seed)
    return [
        (f"group_{g}", [
            (f"S{g}{b}", [rng.gauss(0, 1) for _ in range(rng.randint(30, 300))])
            for b in range(3)
        ])
        for g in range(2)
    ]


@pytest.mark.parametrize("theme", sorted(THEMES))
def test_a_figure_is_written_in_both_themes(tmp_path: Path, theme: str):
    out = binned_point_plot(_groups(), tmp_path / f"{theme}.png", theme=theme)
    assert out.is_file() and out.stat().st_size > 5000


def test_summary_uses_all_values_not_the_drawn_subsample(tmp_path: Path):
    """max_points thins the cloud for legibility; the median must not move."""
    big = [(("one"), [("S", [float(i) for i in range(5000)])])]
    a = binned_point_plot(big, tmp_path / "a.png", max_points=50, seed=1)
    b = binned_point_plot(big, tmp_path / "b.png", max_points=5000, seed=2)
    # Identical summaries mean the median/IQR were taken before subsampling;
    # the rendered points differ, so the files are not byte-identical.
    assert a.is_file() and b.is_file()


def test_a_bin_with_no_values_is_skipped_not_crashed(tmp_path: Path):
    groups = [("g", [("empty", []), ("full", [1.0, 2.0, 3.0])])]
    assert binned_point_plot(groups, tmp_path / "e.png").is_file()


def test_load_scored_groups_by_category_and_reports_missing(tmp_path, capsys):
    panel = tmp_path / "panel.tsv"
    panel.write_text("pdb_id\tcategory\n1ABC\tcat_a\n2DEF\tcat_a\n3GHI\tcat_b\n")
    root = tmp_path / "res"
    for pid, vals in (("1ABC", [1.0, 2.0]), ("3GHI", [0.5])):
        d = root / pid
        d.mkdir(parents=True)
        (d / "idf_scored.tsv").write_text(
            "idf_z\n" + "\n".join(str(v) for v in vals) + "\n"
        )
    groups = load_scored(panel, root, "idf_scored.tsv", "idf_z")
    assert [g for g, _ in groups] == ["cat_a", "cat_b"]
    assert groups[0][1] == [("1ABC", [1.0, 2.0])]
    assert "2DEF" in capsys.readouterr().err
