"""What the simulation figures DRAW, checked against the committed fixtures.

`tests/test_quantitative_figure_data.py` guards these figures' CAPTIONS against
`tests/fixtures/*.json`, and `test_flagship_figure_data.py` does the same for
fig27's. Neither opens a PDF. So the caption could be right and the figure
wrong -- which is exactly what happened to fig17 in #790, where the panels came
from one simulation scenario and the numbers printed on them from another, and
every caption guard passed.

These figures read `simulations/output/`, which is gitignored, so CI cannot
regenerate them (#788). But it does not have to: the numbers they draw are
already committed as fixtures for the caption guards, so comparing the two
gates the artifact a reader actually looks at, in CI, with no new data and
without tracking simulation output.

WHAT THIS DOES NOT DO. It cannot notice a figure whose PLOT is wrong while its
annotations are right -- a mis-drawn curve, a swapped panel, a wrong colour
map. It compares the numbers rendered as text. That is a smaller claim than
"the figure is correct", and it is the claim the committed data can support.

Four narrower limits, each measured rather than assumed:

- **A change too small to survive rounding is invisible.** fig24's four bars
  are drawn `:.1f`, and RSL3's hypoxic value is 0.1% whether it is the mean
  over the four lambdas or the single 120um condition. Making that swap for
  RSL3 is caught, but only by the collapse ANNOTATION, which is `:.0f` and
  moves. SDT has no such annotation, so the same swap made for SDT alone would
  leave the PDF text byte-identical and nothing here would fire.
- **Bar HEIGHTS are not read.** Swapping the two height lists on fig24's
  `axA.bar` calls, so the blue "Normoxic" bars are drawn at hypoxic heights
  while every label and annotation stays put, changes no text and passes. That
  is the plot-level limit above, stated concretely because it is reachable by
  a one-line edit.
- **fig25 panel (b) is backed for one pair.** `tests/fixtures/bliss_synergy.json`
  holds `rsl3_fsp1i` only, so the other two scores are compared against
  nothing; dropping or altering `FSP1i+HDACi` cannot be detected here.
- **fig24 panel (b) is ungated**, as is fig26 panel (b)'s GPX4 right-hand axis,
  whose tick labels are data-derived numbers no assertion reads.
- **`closes ~day 3` is a presence check.** It is a hardcoded string in the
  generator, so it would still read day 3 if the window moved. It is the ONLY
  such string left in these three figures: the two cohort footnotes read as
  hardcoded too, and both turned out to be real fields of the simulation
  output, so they are compared against the fixtures rather than listed here.

BINDING. matplotlib writes a bar's annotation and its axis label into the text
stream as separate runs with nothing tying the two together, so for fig24 and
fig25's panel (a) a value is bound to its bar by ORDER. Order alone is not
enough: with only the values pinned, swapping fig24's two group labels or
reordering fig25's four bar labels leaves every number exactly where it was and
passes. Both tests therefore also pin the LABELS in the order they are drawn,
which is what makes the positional binding say whose bar is whose.

ORDER IN THE TEXT STREAM IS DRAW ORDER, NOT LAYOUT, and for fig24 the two
differ. Its annotations are placed at hand-computed `xi -/+ w/2` offsets with
no connection to the bars they describe, so swapping those two offsets moves
each number over the other bar -- RSL3 re-captioned as 0.1% normoxic and 3.7%
hypoxic, the same inversion the legend check was added to stop -- while the
text stream stays byte-identical. An earlier version of this paragraph claimed
pinning the labels "is what makes the positional binding say whose bar is
whose"; it does not, and the mutation that proved it passed every assertion in
this file. fig24's numbers are therefore ALSO compared left-to-right by
bounding box, which is the check that actually binds a value to a bar.
fig25's panel (a) needs no such check: it anchors each annotation to
`bar.get_x() + bar.get_width() / 2`, so its text cannot drift off its bar.

fig25's panel (b) is the one place a label sits beside its own value, so there
each score is bound to its pair directly. A set comparison is not enough --
two scores swapped between pairs would satisfy one -- and that is the defect
the fig17 guard was fixed for in #790.
"""
import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
FIG_DIR = REPO / "article/figures"
FIXTURES = REPO / "tests/fixtures"


def _drawn(stem):
    try:
        import pymupdf
    except ImportError:
        try:
            import fitz as pymupdf
        except ImportError:
            pytest.skip("no PDF reader available")
    path = FIG_DIR / f"{stem}.pdf"
    assert path.exists(), f"{stem}.pdf is not committed"
    doc = pymupdf.open(path)
    try:
        return " ".join(" ".join(p.get_text().split()) for p in doc)
    finally:
        doc.close()


def _words(stem):
    """Every word in the PDF with its bounding box, page 0.

    `_drawn` returns the text STREAM, which is draw order, not layout. For a
    matplotlib figure those two coincide only when the artist positions are
    themselves in order -- and fig24's are not: it places each annotation at a
    hand-computed `xi -/+ w/2` offset, entirely independent of the bar it is
    meant to sit over. Swapping the two offsets moves both numbers onto the
    wrong bars and leaves the text stream byte-identical, so no assertion over
    `_drawn` can see it. Geometry can.
    """
    try:
        import pymupdf
    except ImportError:
        try:
            import fitz as pymupdf
        except ImportError:
            pytest.skip("no PDF reader available")
    path = FIG_DIR / f"{stem}.pdf"
    assert path.exists(), f"{stem}.pdf is not committed"
    doc = pymupdf.open(path)
    try:
        return [(w[0], w[1], w[4]) for w in doc[0].get_text("words")]
    finally:
        doc.close()


def _pct(value):
    """A rate rendered the way the generators render it."""
    return f"{value * 100:.1f}%"


def test_fig25_draws_its_bliss_numbers():
    """Panel (a): RSL3 alone, FSP1i alone, Bliss expected, observed combination.

    Bar order is the order the generator appends them, which is also the order
    the caption reads.
    """
    fx = json.loads((FIXTURES / "bliss_synergy.json").read_text())["rsl3_fsp1i"]
    text = _drawn("fig25_bliss_synergy")
    expected = [_pct(fx[k]) for k in
                ("rate_a", "rate_b", "bliss_prediction", "rate_combo")]
    found = re.findall(r"\d+\.\d%", text)
    assert found[:4] == expected, (
        f"fig25's panel (a) draws {found[:4]} and bliss_synergy.json gives "
        f"{expected}. Regenerate the figure, or the fixture if the engine moved")
    # THE BAR LABELS, IN ORDER, for the same reason as fig24's: reordering
    # them without the values captioned the Bliss expectation as the observed
    # combination and passed.
    bars = ("RSL3 alone", "FSP1i alone", "Bliss expected",
            "Observed combination")
    for label in bars:
        assert label in text, f"fig25 no longer labels a bar {label!r}"
    order = [text.index(l) for l in bars]
    assert order == sorted(order), (
        "fig25's panel (a) labels are not in the order its values are "
        f"compared in; the bars may be mislabelled. positions={order}")
    # SCOPED TO PANEL (a). A bare `in text` was satisfied by panel (b), which
    # draws the same score as a bar label -- so replacing panel (a)'s
    # annotation with a wrong value, or deleting it, both passed.
    assert f"{fx['synergy_score']:.2f}× synergy" in text, (
        f"fig25's panel (a) does not annotate "
        f"{fx['synergy_score']:.2f}x synergy; panel (b) drawing that number "
        "as a bar label is not the same claim")
    # THE COHORT SIZE. It is drawn as a footnote and was compared against
    # nothing, so the figure could claim any n. It is a real field of the
    # simulation output, so it is now a fixture value rather than a scope
    # limit -- the same move that turned fig26's timepoints from "nothing to
    # compare" into a comparison.
    whole = json.loads((FIXTURES / "bliss_synergy.json").read_text())
    n = whole["n_cells_per_condition"]
    assert f"{n:,} persister cells/condition" in text, (
        f"fig25's footnote does not say {n:,} persister cells/condition, "
        "which is what combo_summary.json records as n_cells_per_condition")


def test_fig25_binds_each_pair_to_its_own_score():
    """Panel (b) labels each bar, so a set comparison is not good enough: two
    scores swapped between pairs would pass one."""
    fx = json.loads((FIXTURES / "bliss_synergy.json").read_text())["rsl3_fsp1i"]
    text = _drawn("fig25_bliss_synergy")
    # PAIR NAMES DERIVED FROM THE TEXT, not listed. A hardcoded alternation of
    # the three current pairs fails the moment the engine emits a fourth --
    # measured: adding RSL3+ACSL4i made a correct figure fail with "draws 1
    # pairs and 4 scores" -- and fails again if the sim ever emits the pair
    # with drug_a and drug_b the other way round, which the generator's own
    # lookup deliberately tolerates.
    #
    # THE CHARACTER CLASS IS `[^\s+]`, NOT `[A-Za-z0-9]`. Deriving the names
    # is only half the fix: an alphanumeric class cannot span a hyphen, and
    # ferroptosis drug names are full of them. Adding a legitimate `Fer-1 +
    # HDACi` pair to the fixture made a CORRECT figure fail -- and worse, the
    # class silently matched the tail `1+HDACi`, so the guard reported "2 pairs
    # and 4 scores" rather than saying it could not read the name.
    #
    # The pair added was `Fer-1`+`HDACi`, and it was added to
    # `simulations/output/combo-mech/combo_summary.json`, which is what panel
    # (b) actually reads. An earlier version of this comment said "to the
    # fixture"; adding a pair to `tests/fixtures/bliss_synergy.json` changes
    # nothing at all, because that file is read only as `["rsl3_fsp1i"]` and
    # panel (b) never touches it. The mutation was right and the sentence
    # pointed at the wrong file.
    labels = re.search(r"((?:[^\s+]+\+[^\s+]+\s+)+)"
                       r"((?:\d+\.\d+×\s*)+)", text)
    assert labels, "fig25's panel (b) no longer lists its pairs and scores"
    pairs = labels.group(1).split()
    scores = re.findall(r"\d+\.\d+×", labels.group(2))
    assert len(pairs) == len(scores), (
        f"fig25 draws {len(pairs)} pairs and {len(scores)} scores")
    mapping = dict(zip(pairs, scores))
    key = next((p for p in pairs if set(p.split("+")) == {"RSL3", "FSP1i"}), None)
    assert key, f"fig25's panel (b) no longer draws the RSL3/FSP1i pair: {pairs}"
    assert mapping.get(key) == f"{fx['synergy_score']:.2f}×", (
        f"fig25 gives {key} {mapping.get(key)} and the fixture "
        f"says {fx['synergy_score']:.2f}x -- the label and the score it is "
        "drawn beside disagree with the data")


def test_fig24_draws_its_hypoxia_numbers():
    """Panel (a): RSL3 normoxic, RSL3 hypoxic, SDT normoxic, SDT hypoxic.

    THE HYPOXIC BAR IS THE MEAN OVER ALL FOUR LAMBDAS, not `gradient_120um`.
    A first version bound it to the single 120um condition and PASSED, because
    at one decimal the two collide on today's data (RSL3 0.0009841 against
    0.0010249, both "0.1%"; SDT 0.8780347 against 0.8781699, both "87.8%") --
    so the guard pinned a rounding coincidence, missed a generator changed to
    draw a single lambda, and would have failed a legitimate regeneration whose
    four lambdas spread. The derivation is imported from the caption guard next
    door, which had it right, so the two cannot diverge.
    """
    import sys as _sys

    _sys.path.insert(0, str(Path(__file__).resolve().parent))
    from test_quantitative_figure_data import _hypoxia_kills

    rows = json.loads(
        (FIXTURES / "hypoxia_killcurve_rows.json").read_text())["conditions"]
    kills = _hypoxia_kills(rows)          # {treatment: (normoxic, hypoxic)} in %
    expected = [f"{kills['RSL3'][0]:.1f}%", f"{kills['RSL3'][1]:.1f}%",
                f"{kills['SDT'][0]:.1f}%", f"{kills['SDT'][1]:.1f}%"]
    text = _drawn("fig24_hypoxia_killcurve")
    found = re.findall(r"\d+\.\d%", text)
    assert found[:4] == expected, (
        f"fig24 draws {found[:4]} and hypoxia_killcurve_rows.json gives "
        f"{expected}. NOTE this compares DRAW ORDER, so a deliberate reordering "
        "of the two treatment groups fails here too and reports as data drift; "
        "the bar order is fixed by the generator's append order on purpose")

    # BY POSITION, not by draw order -- this is the assertion that makes the
    # binding real. The generator places each annotation at a hand-computed
    # `xi - w/2` / `xi + w/2` offset that has no connection to the bar it sits
    # over, so swapping the two offsets re-captions RSL3 as 0.1% normoxic and
    # 3.7% hypoxic, inverting the figure's claim, and leaves the text stream
    # BYTE-IDENTICAL. Every assertion above passed on exactly that mutation.
    #
    # Panel (a) is the left subplot; the collapse annotation draws two more
    # percentages, and they are the two sharing a row with the arrow, so the
    # arrow's y locates the row to drop rather than a hardcoded coordinate.
    words = _words("fig24_hypoxia_killcurve")
    panel_a = [w for w in words if w[0] < 330]
    arrow_y = [y for x, y, t in panel_a if t.strip() in {"→", "\u2192"}]
    assert len(arrow_y) == 1, (
        f"expected exactly one arrow in fig24 panel (a), found {len(arrow_y)}; "
        "the collapse annotation is what identifies the row to exclude")
    bars = sorted(((x, t) for x, y, t in panel_a
                   if re.fullmatch(r"\d+\.\d%", t)
                   and abs(y - arrow_y[0]) > 3), key=lambda w: w[0])
    assert [t for _, t in bars] == expected, (
        f"fig24's bar annotations read {[t for _, t in bars]} LEFT TO RIGHT "
        f"and the data gives {expected}. The numbers are drawn at fixed "
        "offsets, so this is the only check that says each one sits over the "
        "bar it describes")
    # THE GROUP LABELS, IN ORDER. Without this every bar can be relabelled --
    # swapping the two group labels puts RSL3's 3.7%/0.1% under the SDT
    # heading and passed, which is the #790 defect this file exists to catch.
    labels = re.search(r"(RSL3 \(GPX4 inhibitor\)) (SDT \(exogenous ROS\))", text)
    assert labels, (
        "fig24's group labels are not RSL3 then SDT in that order, so the "
        "positional binding above no longer says which bar is whose")
    # AND THE SERIES LEGEND, which is the other half of a 2x2 and was the
    # weaker half: pinning only the group axis left "which of these two bars
    # is the hypoxic one" unstated. Swapping the two `label=` strings on
    # `axA.bar` re-captions 3.7%/91.9% as the HYPOXIC bars and 0.1%/87.8% as
    # the normoxic ones -- the thesis of the figure inverted -- and passed.
    # The legend is emitted in bar-call order, the same order as the values.
    series = re.search(r"Normoxic \(uniform O2\)\s+Hypoxic \(O2 gradient\)", text)
    assert series, (
        "fig24's series legend is not `Normoxic (uniform O2)` then `Hypoxic "
        "(O2 gradient)`. The bar VALUES are compared in that order, so if the "
        "legend disagrees the figure says hypoxia raises the kill rate")
    # And the collapse annotation must be the ratio of the two it names, not a
    # number carried along beside them.
    collapse = re.search(r"(\d+\.\d%) → (\d+\.\d%) \(~(\d+)× collapse\)", text)
    assert collapse, "fig24 no longer annotates the collapse ratio"
    assert [collapse.group(1), collapse.group(2)] == expected[:2], (
        f"fig24's collapse annotation names {collapse.group(1)} → "
        f"{collapse.group(2)} while its bars draw {expected[:2]}")
    # THE TOLERANCE IS DERIVED FROM THE FORMAT STRING, not from the magnitude.
    # The generator renders this ratio with `:.0f`, so a CORRECT figure can
    # differ from the fixture by at most 0.5 -- at any magnitude, since
    # rounding error does not grow with the value. 0.6 leaves that half-unit
    # of headroom and nothing more.
    #
    # Two earlier attempts were looser for reasons that do not survive being
    # checked. An absolute 1.0 admits a genuinely wrong annotation: adding 0.9
    # before formatting draws `~38x` against a data ratio of 37.18 and passed.
    # A `max(1.0, 0.02 * ratio)` was no better -- at 37.18 the relative term is
    # 0.744, so the floor always won and the relative branch never ran.
    #
    # The clamp matches the generator's own `max(hyp, 0.01)`; both operate on
    # percentages, so a hypoxic rate at or below 0.01% is pinned identically on
    # both sides instead of one of them dividing by zero.
    ratio = kills["RSL3"][0] / max(kills["RSL3"][1], 0.01)
    drawn_ratio = int(collapse.group(3))
    assert abs(drawn_ratio - ratio) <= 0.6, (
        f"fig24 says ~{drawn_ratio}x collapse and the fixture gives "
        f"{ratio:.1f}x. The generator formats with `:.0f`, so a correct "
        "figure is within 0.5 of the data")


def test_fig26_draws_the_timepoints_its_fixture_holds():
    """fig26 DOES draw data values: its x tick labels are the fixture's
    `timepoint_days`, rendered by the generator as `f"{d:g}"`.

    An earlier version of this test asserted the opposite -- "it prints no data
    value, so there is nothing to compare" -- and left behind an unused
    `axis_ticks` set listing those very strings, which is the tell: the check
    that would have consumed it was dropped and the sentence claiming it was
    unnecessary stayed. Nothing else gates them, so dropping a timepoint from
    the sweep changes the figure with no test noticing.

    The death rates themselves are plotted as curves, not annotated, so they
    stay ungated here; `test_quantitative_figure_data.py` pins the four the
    caption quotes.
    """
    rows = json.loads((FIXTURES / "vulnerability_window.json").read_text())["rows"]
    days = sorted({r["timepoint_days"] for r in rows})
    expected = [f"{d:g}" for d in days]
    text = _drawn("fig26_vulnerability_window")
    # Both panels share the axis, so the run appears twice; check the first.
    run = re.search(r"\b" + r"\s+".join(re.escape(d) for d in expected) + r"\b",
                    text)
    assert run, (
        f"fig26 does not draw its timepoints {expected} in order; "
        "vulnerability_window.json is what it plots against")
    assert text.count(run.group(0)) >= 2, (
        "fig26 draws the timepoint axis once; both panels share it, so a "
        "single occurrence means one panel lost its labels")
    # No data value beyond the axes is annotated, and if one appears it needs
    # a comparison rather than passing unnoticed.
    values = re.findall(r"\d+\.\d%|\d+\.\d+×", text)
    assert not values, (
        f"fig26 now annotates {values}; add them to this comparison against "
        "vulnerability_window.json rather than leaving them ungated")
    # THE CURVE LEGEND. The timepoints above say WHEN, and nothing said WHICH
    # CURVE IS WHICH. Swapping the two `label=` strings on `axA.plot` labels
    # the curve that holds near 100% out to day 28 as RSL3 and the collapsing
    # one as SDT -- "days for RSL3, weeks for SDT" reversed, under a suptitle
    # still asserting the original claim -- and passed. Legend entries are
    # emitted in plot-call order, so pinning the order pins the pairing.
    legend = re.search(r"SDT \(exogenous ROS\)\s+RSL3 \(GPX4 inhibitor\)", text)
    assert legend, (
        "fig26's panel (a) legend is not `SDT (exogenous ROS)` then `RSL3 "
        "(GPX4 inhibitor)`; the curves may be labelled the wrong way round, "
        "which reverses the window claim the figure is drawn to make")
    whole = json.loads((FIXTURES / "vulnerability_window.json").read_text())
    n = whole["n_cells"]
    assert f"{n:,} cells/condition" in text, (
        f"fig26's footnote does not say {n:,} cells/condition, which is the "
        "n_cells its own rows record")
    # THE UNIT. The timepoints above are the only numbers this test gates, and
    # this is the only text saying what they are measured in and from. Pinning
    # the values while leaving the axis caption free is the half-bound shape
    # the legend checks were added to close.
    assert "Time post-chemotherapy (days)" in text, (
        "fig26's x axis no longer says `Time post-chemotherapy (days)`, so "
        "the timepoints this test compares are drawn without a unit or an "
        "origin")
    assert "closes ~day 3" in text, (
        "fig26's window annotation is gone. NOTE this is a hardcoded string in "
        "the generator, not a derived one -- if the window moved, the "
        "annotation would still read day 3 and this assertion would still "
        "pass, so it is a presence check and not a data guard")
