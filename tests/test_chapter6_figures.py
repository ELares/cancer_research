"""Guards for Chapter 6's figures.

The chapter had ZERO numbered manuscript figures against 24 for the rest of
the book, while citing two supplementary ones and quoting their numbers -- so
the deficit was real and the two it did cite were already doing manuscript
work. These guards exist because a figure is a CLAIM: it can go stale against
the artifact it draws, it can contradict its own caption, and it can present
an uncalibrated placeholder as a result.
"""
import json
import re
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parent.parent
MANUSCRIPT = REPO / "article/drafts/v1.md"
GEN = REPO / "scripts/generate_conceptual_diagrams.py"
FIGURES = REPO / "FIGURES.yaml"
CHAPTER6_FIGURES = {27: "fig32_modality_tme", 28: "fig33_adoptive_barriers"}


@pytest.fixture(scope="module")
def entries():
    return yaml.safe_load(FIGURES.read_text())["figures"]


def _chapter6() -> str:
    t = MANUSCRIPT.read_text()
    return t[t.index("## Chapter 6:"):t.index("## Chapter 7:")]


def test_chapter_six_now_has_numbered_manuscript_figures(entries):
    """The measured gap this closes."""
    by_num = {e["manuscript_figure"]: e for e in entries
              if e.get("manuscript_figure")}
    ch6 = _chapter6()
    cited = {int(n) for n in re.findall(r"Figure (\d+)", ch6)}
    numbered = {n for n in cited if n in by_num}
    assert numbered, (
        "Chapter 6 cites no numbered manuscript figure, which is the deficit "
        f"this file exists to prevent returning; it cites {sorted(cited)}")
    for num, stem in CHAPTER6_FIGURES.items():
        assert num in by_num, f"Figure {num} is no longer registered"
        assert by_num[num]["filename"] == stem
        assert by_num[num]["status"] == "manuscript", (
            f"Figure {num} was demoted to supplementary while the chapter "
            "still cites it as a numbered figure")
        assert num in cited, f"Chapter 6 no longer cites Figure {num}"


def test_the_figure_files_exist_for_every_registered_chapter6_figure(entries):
    for num, stem in CHAPTER6_FIGURES.items():
        for ext in ("pdf", "png"):
            f = REPO / "article/figures" / f"{stem}.{ext}"
            assert f.exists() and f.stat().st_size > 0, f"{f} missing or empty"


def test_the_waterfall_quotes_the_live_collapse():
    """The chapter states the collapse the figure draws.

    Both read `analysis/modality-panel.json`, so a run that moves the panel
    must move both. An earlier version of the decomposition beside this was a
    residual that could not disagree with itself; this is the check that the
    prose cannot disagree with the artifact either.
    """
    ab = json.loads((REPO / "analysis/modality-panel.json").read_text())[
        "adoptive_barriers"]
    collapse = ab["leukaemia_kill_fraction"] / ab["solid_tumour_kill_fraction"]
    assert f"{collapse:,.0f}-fold" in _chapter6(), (
        f"Chapter 6 does not quote the live {collapse:,.0f}-fold collapse")


def test_the_diverging_map_is_not_reversed():
    """The first draft drew every loss BLUE while its caption said red.

    A figure contradicting its own legend is the same defect class as prose
    contradicting the artifact beside it, and it is invisible in any check
    that only asks whether the figure was produced.
    """
    src = GEN.read_text()
    fn = src[src.index("def fig32_modality_tme"):src.index("def fig33_adoptive_barriers")]
    assert 'cmap="RdBu"' in fn, "fig32 no longer uses the RdBu map"
    assert 'cmap="RdBu_r"' not in fn, (
        "fig32 uses the REVERSED map, so losses draw blue while its caption "
        "says red is a loss")
    assert "red is a loss, blue a GAIN" in fn, (
        "the caption no longer states the colour convention it is checked "
        "against")


def test_both_figures_carry_their_refusal():
    """Each caption must say what the figure is NOT.

    Every arm drawn is an uncalibrated placeholder, and a chart of kill
    fractions invites exactly the clinical reading the analysis pages refuse.
    """
    src = GEN.read_text()
    for start, end, must in (
        ("def fig32_modality_tme", "def fig33_adoptive_barriers",
         ("NOT a ranking", "placeholder")),
        ("def fig33_adoptive_barriers", "\nif __name__",
         ("NOT a clinical comparison", "uncalibrated placeholder")),
    ):
        fn = src[src.index(start):src.index(end)]
        for phrase in must:
            assert phrase in fn, f"{start} lost its refusal: {phrase}"
