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
    ch6 = _chapter6()
    assert f"{collapse:,.0f}-fold" in ch6, (
        f"Chapter 6 does not quote the live {collapse:,.0f}-fold collapse")
    # ATTACHED TO THE RIGHT COMPARISON. "collapses 633-fold between two doses
    # of the same drug" passed the presence check, and the whole point of the
    # number is that it is ONE construct against TWO DISEASES.
    i = ch6.index(f"{collapse:,.0f}-fold")
    around = ch6[max(0, i - 220):i + 220]
    assert "same construct" in around and "between the two diseases" in around, (
        "the collapse figure is not attached to the one-construct/two-diseases "
        f"comparison it measures: ...{around[:200]}...")


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
    # AND EACH CHAPTER-6 FIGURE MUST BE CITED IN CHAPTER 6. Membership in the
    # valid range is satisfied by citing any figure anywhere, so swapping a
    # citation for an unrelated but valid number passed.
    ch6 = _chapter6()
    for num in CHAPTER6_FIGURES:
        assert f"Figure {num}" in ch6, (
            f"Figure {num} is a Chapter 6 figure and Chapter 6 does not cite it")
    # AND THE CITATION THIS PR RENUMBERED MUST STILL POINT AT ITS OWN IMAGE.
    # A number inside the valid range says nothing about WHICH figure is meant:
    # replacing Chapter 11's landscape citation with an unrelated valid number
    # left every other assertion green. These are the sentences whose subject
    # is fixed by what the figure draws, so each is pinned to its figure.
    by_file = {e["filename"]: e["manuscript_figure"] for e in entries
               if e.get("manuscript_figure")}
    for stem, phrase in (
        ("fig30_modality_landscape", "shows the shape of that"),
        ("fig30_modality_landscape", "Nine treatment arms exist where four did"),
        ("fig31_modality_panel", "runs every applicable arm against the identical tumour"),
    ):
        n = by_file[stem]
        i = text.find(phrase)
        assert i >= 0, f"the sentence citing {stem} is gone: {phrase!r}"
        near = text[max(0, i - 260):i + 260]
        assert f"Figure {n}" in near, (
            f"the sentence {phrase!r} no longer cites Figure {n} ({stem}); "
            f"it reads: ...{near[:180]}...")


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
    # ATTACHED TO THE RIGHT ARM. Presence alone let a reviewer swap the two
    # numbers -- "AntibodyDrugConjugate kills 87.24% against sonodynamic
    # therapy's 1.84%" -- inverting the chapter's whole delivery argument
    # while every assertion passed.
    adc, sdt = arms["AntibodyDrugConjugate"], arms["SDT"]
    claim = (f"`AntibodyDrugConjugate` kills {adc:.2%} against sonodynamic "
             f"therapy's {sdt:.2%}")
    assert claim in text, (
        f"Chapter 6 does not attach the panel numbers to their arms: "
        f"expected {claim!r}")
    assert adc < sdt, "the delivery argument depends on the ADC being smaller"
    assert "1.71% against sonodynamic" not in text


def test_the_dominant_axis_sentence_matches_the_sweep():
    """Chapter 6 names which axis dominates in each cell state.

    It named the wrong axis for one state and quoted 121% for the other --
    which is the RETRACTED absolute-value artefact the same section calls
    impossible two paragraphs below. The figure drew +124 in that cell.
    """
    tme = json.loads((REPO / "analysis/modality-tme.json").read_text())
    text = MANUSCRIPT.read_text()
    for state, r in tme["dominant_axis"].items():
        # WITH ITS SIGN AND ITS ARM. The first version compared a magnitude
        # against an absolute-valued field, so flipping a total loss to a
        # doubling left it green -- and the sentence it guarded called a GAIN
        # a "pressure" for exactly that reason.
        claim = (f"the {state} state it is {r['axis']}, a {r['direction']} of "
                 f"{abs(r['value']):.0%} for `{r['arm']}`")
        assert claim in text, (
            f"Chapter 6 does not state the live dominant axis with its "
            f"direction: expected {claim!r}")
        if r["tied_with"]:
            assert "is arbitrary" in text, (
                f"{state}'s axis ties with {r['tied_with']} and the chapter "
                "presents it as a unique winner")


def test_the_published_captions_state_the_colour_convention_correctly():
    """The GENERATOR's caption was guarded and the PUBLISHED ones were not.

    A reviewer inverted the manuscript and LaTeX captions to "blue is a loss
    and red a gain", regenerated, and the whole suite passed -- because the
    guard scanned `generate_conceptual_diagrams.py`. What a reader sees is the
    published caption, so that is what must be checked.
    """
    # BOTH published artifacts. The first version's "LaTeX" arm read
    # `generate_latex.py` -- the generator source -- which is exactly the
    # defect its own docstring says it retired. `v1.tex` is committed and is
    # what compiles.
    md = MANUSCRIPT.read_text()
    tex = (REPO / "article/drafts/v1.tex").read_text()
    for name, src in (("manuscript", md), ("v1.tex", tex)):
        i = src.find("fig32") if name == "v1.tex" else src.find("[FIGURE 27:")
        assert i >= 0, f"the {name} caption for Figure 27 is missing"
        cap = src[i:i + 900]
        assert "Red is a loss and blue a gain" in cap or \
               "red is a loss and blue a gain" in cap, (
            f"the {name} caption does not state the colour convention the "
            f"figure draws: {cap[:200]}")
        assert "lue is a loss" not in cap, (
            f"the {name} caption has the colour convention REVERSED")


def test_the_drawn_colours_match_the_stated_convention():
    """Read the colour out of the PDF, not the colormap name out of the source.

    A reviewer regenerated fig32 with the reversed map and then reverted only
    the SOURCE. The committed figure drew every loss blue, every caption said
    red, and all twelve guards passed -- because none of them looked at the
    artifact's pixels.
    """
    fitz = pytest.importorskip("fitz")
    tme = json.loads((REPO / "analysis/modality-tme.json").read_text())
    flat = [(v, ph, ax, arm)
            for ph, axmap in tme["effects_by_phenotype"].items()
            for ax, arms in axmap.items() for arm, v in arms.items()]
    worst = min(flat)[0]
    best = max(flat)[0]
    assert worst < 0 < best, "the sweep has no signed spread to check"
    doc = fitz.open(REPO / "article/figures/fig32_modality_tme.pdf")
    pix = doc[0].get_pixmap(dpi=100)
    doc.close()
    px = [(pix.pixel(x, y)) for y in range(0, pix.height, 3)
          for x in range(0, pix.width, 3)]
    reds = sum(1 for r, g, b in px if r > g + 40 and r > b + 40)
    blues = sum(1 for r, g, b in px if b > r + 40 and b > g + 40)
    assert reds > 0 and blues > 0, (
        f"fig32 draws no diverging colour at all (red px {reds}, blue {blues})")
    # Losses outnumber gains in this sweep, so red must outnumber blue. If the
    # map is reversed the inequality flips.
    n_loss = sum(1 for v, *_ in flat if v < -0.005)
    n_gain = sum(1 for v, *_ in flat if v > 0.005)
    assert n_loss > n_gain, "the sweep is no longer loss-dominated"
    assert reds > blues, (
        f"fig32 draws more BLUE than red ({blues} vs {reds}) while the sweep "
        f"holds {n_loss} losses against {n_gain} gains and every caption says "
        "red is a loss -- the colormap is reversed in the committed artifact")


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
        page = doc[0]
        spans = [s0 for b in page.get_text("dict")["blocks"]
                 for l in b.get("lines", []) for s0 in l["spans"]]
        flat = " ".join(" ".join(pg.get_text("text").split()) for pg in doc)
        doc.close()
        for phrase in must:
            assert phrase in flat, (
                f"{stem}.pdf does not DRAW its refusal: {phrase!r}")
        # AND IT MUST BE VISIBLE. Setting the refusal to `fontsize=0.0,
        # color="white"` left it extractable by PyMuPDF and unreadable by a
        # human, which passed the text-only check.
        head = must[0].split()[0]
        cands = [s0 for s0 in spans if head in s0["text"]]
        assert cands, f"{stem}: no span carries the refusal"
        for s0 in cands:
            assert s0["size"] >= 5.0, (
                f"{stem}: the refusal is drawn at {s0['size']}pt, which no "
                "reader can see")
            assert s0["color"] != 0xFFFFFF, (
                f"{stem}: the refusal is drawn white on white")


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
    # CELL IDENTITY, not a bag of tokens. The first version collected the
    # drawn numbers into a SET, so reversing every row -- putting SDT's -77 on
    # Ablation -- passed, and a suppressed label passed too because the
    # colorbar's extreme tick is computed from the same maximum and reprints
    # the value. Positions are read from the PDF and matched to the grid.
    doc = fitz.open(REPO / "article/figures/fig32_modality_tme.pdf")
    page = doc[0]
    spans = [(round(s0["bbox"][0], 1), round(s0["bbox"][1], 1), s0["text"].strip())
             for b in page.get_text("dict")["blocks"]
             for l in b.get("lines", []) for s0 in l["spans"]]
    doc.close()
    labels = {t for _, _, t in spans}
    arms = sorted({a for ax in tme["effects_by_phenotype"]["glycolytic"].values()
                   for a in ax})
    # Row order is alphabetical by arm, so a row-reversal moves every arm
    # label's y relative to its values. Anchor on the ARM NAME's y and require
    # each of its drawn values to share that row.
    ys = {t: y for x, y, t in spans if t in arms}
    assert len(ys) == len(arms), (
        f"fig32 does not label every arm row: missing {sorted(set(arms) - set(ys))}")
    missing = []
    for ph, axmap in tme["effects_by_phenotype"].items():
        for axis, per_arm in axmap.items():
            for arm, v in per_arm.items():
                if abs(v) < 0.005 or v == 0.0:
                    continue
                want = f"{v * 100:+.0f}"
                row_y = ys[arm]
                on_row = [t for x, y, t in spans if abs(y - row_y) < 4.0]
                if want not in on_row:
                    missing.append(f"{ph}/{axis}/{arm}={want}")
    assert not missing, (
        f"fig32 draws these measured effects on the wrong row, or not at all: "
        f"{sorted(missing)[:6]}")
    assert f"{max(abs(v) for ax in tme['effects_by_phenotype'].values() for a in ax.values() for v in a.values()) * 100:+.0f}" in labels

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
    # BEFORE THE CAPTION. `\\caption` is what consumes the counter, so a
    # setcounter placed after it -- but still before the label -- printed 6.26
    # and passed a regex that only required setcounter-then-label.
    pinned = {fn: int(n) + 1 for n, fn in re.findall(
        r"\\setcounter\{figure\}\{(\d+)\}[^\\]*\\includegraphics"
        r".*?\\caption\{.*?\\label\{fig:([^}]+)\}", tex, re.S)}
    unpinned = [f for f in floats if f not in pinned]
    assert not unpinned, (
        f"these floats do not pin their printed number, so LaTeX will number "
        f"them by position: {unpinned}")
    wrong = {f: (pinned[f], want.get(f)) for f in floats
             if want.get(f) != pinned[f]}
    assert not wrong, (
        f"the compiled number disagrees with FIGURES.yaml (printed, expected): "
        f"{wrong}")


def test_the_latex_cannot_silently_stop_compiling():
    """Two one-character faults took the whole book out, and no test saw it.

    A caption carried `0\\%%`: `\\%` sets a percent and the SECOND `%` opens a
    TeX comment that swallows the rest of the line INCLUDING the closing brace,
    so the build aborted in Chapter 5 and Chapters 6-12 -- every figure this
    file guards -- never rendered. Separately `report.cls` numbers figures
    `\\thechapter.\\arabic{figure}`, so the per-float `\\setcounter` printed 6.25
    where the prose said 25 and NOT ONE of the 29 citations resolved.

    Neither is detectable without compiling, and there is no LaTeX toolchain in
    this environment, so these check the two specific faults rather than
    pretending to check the build.
    """
    tex = (REPO / "article/drafts/v1.tex").read_text()
    bad = re.findall(r"\\%%", tex)
    assert not bad, (
        f"{len(bad)} occurrence(s) of `\\\\%%` in v1.tex: the second % comments "
        "out the rest of the line, including the caption's closing brace, and "
        "the manuscript stops compiling there")
    assert r"\renewcommand{\thefigure}{\arabic{figure}}" in tex, (
        "v1.tex does not flatten figure numbering, so report.cls will print "
        "chapter-prefixed numbers (6.25) while the prose says 25")
