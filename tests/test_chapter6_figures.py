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


def test_every_figure_reference_in_the_manuscript_resolves():
    """A dangling `Figure N` is a reader following a pointer to nothing.

    This PR renumbered two Chapter 6 citations and missed a THIRD citation of
    the same image in Chapter 11, so one figure was cited under two numbers
    and one of them was above the highest that exists. Nothing caught it,
    because the chapter guard slices Chapter 6 and no test looked at the rest
    of the manuscript.
    """
    entries = yaml.safe_load(FIGURES.read_text())["figures"]
    registered = {e["manuscript_figure"] for e in entries
                  if e.get("manuscript_figure")}
    text = MANUSCRIPT.read_text()
    cited = {int(n) for n in re.findall(r"Figure (\d+)", text)}
    dangling = sorted(n for n in cited if n not in registered)
    assert not dangling, (
        f"the manuscript cites figures that are not registered: {dangling}; "
        f"registered are {min(registered)}-{max(registered)}")


def test_the_manuscript_quotes_the_panel_numbers_it_draws():
    """Chapter 6 states arm kill fractions beside the figure that draws them.

    The ADC row gained a bystander term and moved 1.71% -> 1.84%; the figure
    was regenerated and the sentence beside it was not, so the promoted figure
    and the adjacent prose disagreed in the same commit. No test read panel
    numbers out of the manuscript.
    """
    panel = json.loads((REPO / "analysis/modality-panel.json").read_text())
    arms = {a["arm"]: a["kill_fraction"] for a in panel["arms"]}
    text = MANUSCRIPT.read_text()
    for arm in ("AntibodyDrugConjugate", "SDT"):
        live = f"{arms[arm]:.2%}"
        assert live in text, (
            f"the manuscript does not quote {arm}'s live kill fraction {live}")
    # and the stale pair must be gone as a live claim
    assert "1.71% against sonodynamic" not in text


def test_the_dominant_axis_sentence_matches_the_sweep():
    """Chapter 6 names which axis dominates in each cell state.

    It named the wrong axis for one state and quoted 121% for the other --
    which is the RETRACTED absolute-value artefact the same section calls
    impossible two paragraphs below. The figure drew +124 in that cell.
    """
    tme = json.loads((REPO / "analysis/modality-tme.json").read_text())
    text = MANUSCRIPT.read_text()
    for state, (axis, value) in tme["dominant_axis"].items():
        claim = f"{axis} in the {state} state ({value:.0%})"
        assert claim in text, (
            f"Chapter 6 does not state the live dominant axis: expected "
            f"{claim!r}")


def test_the_published_captions_state_the_colour_convention_correctly():
    """The GENERATOR's caption was guarded and the PUBLISHED ones were not.

    A reviewer inverted the manuscript and LaTeX captions to "blue is a loss
    and red a gain", regenerated, and the whole suite passed -- because the
    guard scanned `generate_conceptual_diagrams.py`. What a reader sees is the
    published caption, so that is what must be checked.
    """
    md = MANUSCRIPT.read_text()
    tex = (REPO / "scripts/generate_latex.py").read_text()
    for name, src in (("manuscript", md), ("LaTeX", tex)):
        i = src.find("fig32") if name == "LaTeX" else src.find("[FIGURE 27:")
        assert i >= 0, f"the {name} caption for Figure 27 is missing"
        cap = src[i:i + 900]
        assert "Red is a loss and blue a gain" in cap or \
               "red is a loss and blue a gain" in cap, (
            f"the {name} caption does not state the colour convention the "
            f"figure draws: {cap[:200]}")
        assert "lue is a loss" not in cap, (
            f"the {name} caption has the colour convention REVERSED")


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


def test_both_figures_actually_DRAW_their_refusal():
    """Read the PDF, not the generator.

    The first version substring-scanned `generate_conceptual_diagrams.py`, so
    a reviewer deleted the `fig.text(...)` refusal from BOTH figures, left a
    Python COMMENT carrying the same words, regenerated, and the guard passed
    while the drawn figures contained no refusal at all. What a reader sees is
    the PDF, so that is what is checked.
    """
    fitz = pytest.importorskip("fitz")
    for stem, must in (
        ("fig32_modality_tme", ("NOT a ranking", "placeholder")),
        ("fig33_adoptive_barriers",
         ("NOT a clinical comparison", "uncalibrated placeholder")),
    ):
        doc = fitz.open(REPO / "article/figures" / f"{stem}.pdf")
        text = " ".join(pg.get_text("text") for pg in doc)
        doc.close()
        flat = " ".join(text.split())
        for phrase in must:
            assert phrase in flat, (
                f"{stem}.pdf does not DRAW its refusal: {phrase!r}")


def test_the_drawn_numbers_are_the_artifact_numbers():
    """Close the hole the freshness docstring claimed was already closed.

    That docstring asserted every number on these figures was pinned by a
    test. It was false: no test named either figure, and multiplying every
    drawn cell label in fig32 by a thousand left the whole suite green. The
    same hole let fig31 sit stale against its own committed input for an
    entire PR cycle with a 4.3x error on one arm.
    """
    fitz = pytest.importorskip("fitz")
    tme = json.loads((REPO / "analysis/modality-tme.json").read_text())
    doc = fitz.open(REPO / "article/figures/fig32_modality_tme.pdf")
    drawn = " ".join(pg.get_text("text") for pg in doc)
    doc.close()
    tokens = set(re.findall(r"[+-]\d+", drawn))
    expected, missing = [], []
    for ph, axmap in tme["effects_by_phenotype"].items():
        for axis, arms in axmap.items():
            for arm, v in arms.items():
                if abs(v) < 0.005:
                    continue
                expected.append(f"{v * 100:+.0f}")
    for e in set(expected):
        if e not in tokens:
            missing.append(e)
    assert not missing, (
        f"fig32 does not draw these measured effects: {sorted(missing)[:6]}; "
        f"it drew {sorted(tokens)[:10]}")

    panel = json.loads((REPO / "analysis/modality-panel.json").read_text())
    ab = panel["adoptive_barriers"]
    doc = fitz.open(REPO / "article/figures/fig33_adoptive_barriers.pdf")
    w = " ".join(" ".join(pg.get_text("text").split()) for pg in doc)
    doc.close()
    collapse = ab["leukaemia_kill_fraction"] / ab["solid_tumour_kill_fraction"]
    assert f"{collapse:,.0f}x collapse" in w, (
        f"fig33 does not draw the live {collapse:,.0f}x collapse")
    for v in (ab["leukaemia_kill_fraction"], ab["solid_tumour_kill_fraction"]):
        assert f"{v * 100:.4g}%" in w, f"fig33 does not draw {v * 100:.4g}%"


def test_the_inert_step_is_labelled_and_counted_from_the_data():
    """`no effect` and the title's inert count were both unguarded.

    Changing the threshold that decides them removed the label and made the
    title say zero while the last two bars stayed identical, with a green
    suite.
    """
    fitz = pytest.importorskip("fitz")
    ab = json.loads((REPO / "analysis/modality-panel.json").read_text())[
        "adoptive_barriers"]
    doc = fitz.open(REPO / "article/figures/fig33_adoptive_barriers.pdf")
    w = " ".join(" ".join(pg.get_text("text").split()) for pg in doc)
    doc.close()
    binds = ab["antigen_ceiling_binds"]
    if binds:
        assert "no effect" not in w, (
            "the ceiling binds but the figure labels its step 'no effect'")
        assert "1 of the three doing nothing" not in w
    else:
        assert "no effect" in w, (
            "the ceiling does not bind and the figure does not say so")
        assert "1 of the three doing nothing here" in w, (
            "the title's inert count disagrees with the data")


def test_the_compiled_pdf_prints_the_number_the_prose_cites():
    """LaTeX numbers floats in DOCUMENT ORDER, not by our numbering.

    The four Chapter 6 figures printed as 13-16 in the compiled PDF while the
    prose beside them said 25-28, so this PR's whole deliverable -- "Chapter 6
    now has numbered manuscript figures" -- was true of the markdown reading
    and false of the artifact a reader downloads. Worse, inserting any float
    shifted every later figure's printed number, silently.
    """
    tex = (REPO / "article/drafts/v1.tex").read_text()
    entries = yaml.safe_load(FIGURES.read_text())["figures"]
    want = {e["filename"]: e["manuscript_figure"] for e in entries
            if e.get("manuscript_figure")}
    floats = re.findall(r"\\label\{fig:([^}]+)\}", tex)
    assert floats, "the manuscript has no figure floats"
    pinned = {fn: int(n) + 1 for n, fn in re.findall(
        r"\\setcounter\{figure\}\{(\d+)\}.*?\\label\{fig:([^}]+)\}", tex, re.S)}
    unpinned = [f for f in floats if f not in pinned]
    assert not unpinned, (
        f"these floats do not pin their printed number, so LaTeX will number "
        f"them by position: {unpinned}")
    wrong = {f: (pinned[f], want.get(f)) for f in floats
             if want.get(f) != pinned[f]}
    assert not wrong, (
        f"the compiled number disagrees with FIGURES.yaml (printed, expected): "
        f"{wrong}")
