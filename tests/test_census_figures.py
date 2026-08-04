"""Guards for the census figure (fig28).

A figure that quotes numbers is a second copy of them, so it can drift away from
the prose it illustrates without anything failing. These recompute the figure's
numbers from the committed atlas JSON and assert they still match what
`analysis/atlas-landscape.md` says, so the figure and the analysis cannot
disagree silently.

The pharmacological/physical class definitions are imported from
`atlas_landscape`, not restated, for the same reason: a hand-written list beside
the real one is how the 12.5:1-vs-9.1:1 discrepancy arose in the first place.
"""

import json
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from atlas_landscape import PHARMACOLOGICAL, PHYSICAL  # noqa: E402

LANDSCAPE_JSON = REPO_ROOT / "analysis" / "atlas-landscape.json"
LANDSCAPE_MD = REPO_ROOT / "analysis" / "atlas-landscape.md"
FIG_DIR = REPO_ROOT / "article" / "figures"


def _rows():
    return [r for r in json.loads(LANDSCAPE_JSON.read_text())["rows"]
            if r["mesh_census"] > 0]


def test_figure_files_exist():
    for ext in ("pdf", "png"):
        f = FIG_DIR / f"fig28_census_capture.{ext}"
        assert f.exists() and f.stat().st_size > 5000, f


def test_capture_spread_matches_the_committed_analysis():
    """Panel A's headline: the 213-fold spread."""
    caps = sorted(r["mesh_frozen"] / r["mesh_census"] for r in _rows()
                  if r["mesh_frozen"])
    spread = caps[-1] / caps[0]
    stated = re.search(r"a (\d+)-fold spread", LANDSCAPE_MD.read_text())
    assert stated, "atlas-landscape.md no longer states a fold spread"
    assert round(spread) == int(stated.group(1)), (
        f"figure would plot {spread:.0f}x, the analysis says {stated.group(1)}x")


def test_panel_b_ratios_match_the_committed_analysis():
    """Panel B's two bars: 9.1:1 by the manuscript's method, 17.6:1 on the census."""
    R = {r["mechanism"]: r for r in _rows()}
    tot = lambda ks, c: sum(R[k][c] for k in ks if k in R and R[k].get(c))  # noqa: E731
    kw = tot(PHARMACOLOGICAL, "keyword_frozen") / tot(PHYSICAL, "keyword_frozen")
    census = tot(PHARMACOLOGICAL, "mesh_census") / tot(PHYSICAL, "mesh_census")

    md = LANDSCAPE_MD.read_text()
    assert f"{kw:.1f} : 1" in md, f"figure would plot {kw:.1f}:1, absent from the analysis"
    assert f"{census:.1f} : 1" in md, f"figure would plot {census:.1f}:1, absent"
    # The direction is the load-bearing part: the census ratio must exceed the
    # manuscript's, which is what "the manuscript understates its own case" means.
    assert census > kw


def test_pharmacological_set_is_curated_not_the_complement():
    """The distinction that moved the ratio from 12.5:1 to 9.1:1.

    If PHARMACOLOGICAL ever becomes "everything not physical", the figure and the
    analysis would both silently start counting delivery platforms and genetic
    tools as drug modalities.
    """
    others = {r["mechanism"] for r in _rows()} - PHYSICAL
    assert PHARMACOLOGICAL < others, "pharmacological set must be a strict subset"
    for platform in ("nanoparticle", "crispr", "oncolytic-virus", "mrna-vaccine"):
        assert platform not in PHARMACOLOGICAL, (
            f"{platform} is a delivery platform or genetic tool, not a drug modality")


def test_mechanisms_without_census_articles_are_excluded_not_zeroed():
    """0/0 plotted as 0% would read as 'never captured', which is a different claim."""
    all_rows = json.loads(LANDSCAPE_JSON.read_text())["rows"]
    assert any(r["mesh_census"] == 0 for r in all_rows), (
        "fixture assumption gone: no zero-census mechanism left to exclude")
    assert all(r["mesh_census"] > 0 for r in _rows())
