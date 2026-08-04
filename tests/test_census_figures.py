"""Guards for the census figure (fig28).

An earlier version of this file recomputed the figure's numbers in a PARALLEL
re-implementation and compared that to the prose. It therefore passed while the
figure plotted the opposite finding, while the log axis was removed, and while
the PNG was replaced with an unrelated chart -- because nothing here ever ran
the figure code.

These call `fig28_census_capture` itself, into a temporary directory, and assert
on what it reports having plotted. The parallel-reimplementation checks are kept
for the prose-vs-data comparison, which is a different question.
"""

import json
import re
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from atlas_landscape import PHARMACOLOGICAL, PHYSICAL  # noqa: E402

LANDSCAPE_JSON = REPO_ROOT / "analysis" / "atlas-landscape.json"
LANDSCAPE_MD = REPO_ROOT / "analysis" / "atlas-landscape.md"
FIG_DIR = REPO_ROOT / "article" / "figures"


def _rows():
    return [r for r in json.loads(LANDSCAPE_JSON.read_text())["rows"]
            if r["mesh_census"] > 0 and r["mesh_frozen"] > 0]


@pytest.fixture(scope="module")
def plotted(tmp_path_factory):
    """Run the real figure code into a temp dir and return what it plotted."""
    matplotlib = pytest.importorskip("matplotlib")
    import generate_census_figures as gcf
    out = tmp_path_factory.mktemp("fig28")
    original = gcf.FIG_DIR
    gcf.FIG_DIR = out
    try:
        result = gcf.fig28_census_capture()
    finally:
        gcf.FIG_DIR = original
    result["_files"] = sorted(p.name for p in out.iterdir())
    return result


def test_figure_files_exist():
    for ext in ("pdf", "png"):
        f = FIG_DIR / f"fig28_census_capture.{ext}"
        assert f.exists() and f.stat().st_size > 5000, f


def test_the_figure_code_emits_both_formats(plotted):
    assert plotted["_files"] == ["fig28_census_capture.pdf", "fig28_census_capture.png"]


def test_panel_a_uses_a_log_axis(plotted):
    """The encoding is load-bearing, not cosmetic.

    A 213-fold spread on a linear axis collapses every mechanism under 10% onto
    the axis, which is the whole finding. Bars were also replaced with dots for
    the same reason: a log axis has no zero, so bar LENGTH would encode the
    auto-chosen left limit rather than the data.
    """
    assert plotted["xscale"] == "log"


def test_the_figure_plots_the_spread_the_analysis_states(plotted):
    """Panel A's headline, read off the figure code rather than recomputed."""
    stated = re.search(r"a (\d+)-fold spread", LANDSCAPE_MD.read_text())
    assert stated, "atlas-landscape.md no longer states a fold spread"
    assert round(plotted["spread"]) == int(stated.group(1)), (
        f"the figure plots {plotted['spread']:.0f}x, the analysis says {stated.group(1)}x")


def test_the_figure_plots_all_three_arms_in_the_right_order(plotted):
    """Panel B must show the method arm, not just the first and last.

    Showing only 9.1 and 17.6 attributes the whole move to corpus selection and
    does not reconcile: 9.1 x 3.3 is 30, not 17.6. The missing factor is the
    method effect, and it runs the other way.
    """
    ratios = plotted["ratios"]
    assert len(ratios) == 3, f"expected three arms, got {len(ratios)}"
    kw, mesh_frozen, census = ratios
    md = LANDSCAPE_MD.read_text()
    for v in ratios:
        assert f"{v:.1f} : 1" in md, f"the figure plots {v:.1f}:1, absent from the analysis"
    # The load-bearing shape: method CUTS the ratio, corpus RAISES it past the start.
    assert mesh_frozen < kw, "the method arm no longer cuts the ratio"
    assert census > kw, "the census arm no longer exceeds the manuscript's figure"


def test_the_figure_drops_only_mechanisms_it_cannot_place(plotted):
    """0/0 plotted as 0% would read as 'never captured', a different claim, and
    a zero capture has no position on a log axis."""
    all_rows = json.loads(LANDSCAPE_JSON.read_text())["rows"]
    expected = sorted(r["mechanism"] for r in all_rows
                      if not (r["mesh_census"] > 0 and r["mesh_frozen"] > 0))
    assert sorted(plotted["dropped"]) == expected
    assert expected, "fixture assumption gone: nothing left to exclude"
    assert all(v > 0 for v in plotted["capture"].values())


def test_panel_b_ratios_match_the_committed_analysis():
    """Independent recomputation, as a cross-check on the figure's own arithmetic."""
    R = {r["mechanism"]: r for r in _rows()}
    tot = lambda ks, c: sum(R[k][c] for k in ks if k in R and R[k].get(c))  # noqa: E731
    kw = tot(PHARMACOLOGICAL, "keyword_frozen") / tot(PHYSICAL, "keyword_frozen")
    census = tot(PHARMACOLOGICAL, "mesh_census") / tot(PHYSICAL, "mesh_census")
    md = LANDSCAPE_MD.read_text()
    assert f"{kw:.1f} : 1" in md and f"{census:.1f} : 1" in md
    assert census > kw


def test_pharmacological_set_is_curated_not_the_complement():
    """The distinction that moved the ratio from 12.5:1 to 9.1:1."""
    others = {r["mechanism"] for r in _rows()} - PHYSICAL
    assert PHARMACOLOGICAL < others, "pharmacological set must be a strict subset"
    for platform in ("nanoparticle", "crispr", "oncolytic-virus", "mrna-vaccine"):
        assert platform not in PHARMACOLOGICAL, (
            f"{platform} is a delivery platform or genetic tool, not a drug modality")
