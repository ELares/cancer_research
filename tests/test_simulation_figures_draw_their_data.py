"""What the simulation figures DRAW, checked against the committed fixtures.

`tests/test_quantitative_figure_data.py` guards these figures' CAPTIONS against
`tests/fixtures/*.json`, and `test_flagship_figure_data.py` does the same for
fig27's. Neither opens a PDF. So the caption could be right and the figure
wrong -- which is exactly what happened to fig17 in #790, where the panels came
from one simulation scenario and the numbers printed on them from another, and
every caption guard passed.

These figures read `simulations/output/`, which is gitignored, so CI cannot
regenerate them (#788). But it does not have to: the numbers they draw are
committed as fixtures, so comparing the two gates the artifact a reader
actually looks at, in CI, without tracking simulation output.

Most of those fixture values were already there for the caption guards. Two
were not -- the cohort sizes fig25 and fig26 state in their footnotes -- and
they are the reason `test_quantitative_figure_data.py`'s live-agreement guards
were extended to cover them. A fixture value that no live guard reads lets the
figure and the fixture agree with each other while both drift from the run.

WHAT THIS DOES NOT DO. It cannot notice a figure whose PLOT is wrong while its
annotations are right -- a mis-drawn curve, a swapped panel, a wrong colour
map. It compares the numbers rendered as text. That is a smaller claim than
"the figure is correct", and it is the claim the committed data can support.

Five narrower limits, each measured rather than assumed. READ THEM AS LIMITS
OF THE SEMANTIC CHECKS ONLY. The fingerprint at the end of this file hashes
as much of the drawing as it has been taught to read, so every mutation named below now fails there -- verified,
each one -- but it fails as "something moved and nothing explained it", which
is a weaker statement than the assertion that would name the defect. These are
the places where no assertion would name it:

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
  holds one combination, `rsl3_fsp1i` (alongside the cohort size), so the other
  scores are compared against nothing; dropping or altering `FSP1i+HDACi`
  cannot be detected here.
- **fig26 panel (b)'s GPX4 right-hand axis is ungated** -- its tick labels are
  data-derived numbers no assertion reads. fig24 panel (b) was in this list
  and is not any more: its two dashed reference lines are bound to its own
  legend below. The rest of that panel -- the four kill curves themselves --
  is still unread except by the fingerprint.
- **`closes ~day 3` is a presence check.** It is a hardcoded string in the
  generator, so it would still read day 3 if the window moved.

  The two cohort footnotes were in the same class and are no longer: an
  earlier version of this bullet claimed they were "compared against the
  fixtures", which was false in the way that matters. They were string
  literals in the generator, so halving n in the simulation output left both
  captions unchanged and every guard green -- each figure then overstated its
  cohort by 2x, with nothing able to say so. (An earlier version of this
  sentence said "2x and 200x"; 200x is what a different mutation gives, not
  the halving described beside it.) Comparing a drawn
  constant to a fixture cannot fail. The GENERATOR now derives both numbers
  from the run it plots, which is what makes the comparison mean anything.

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
An earlier version of this paragraph exempted fig25's panel (a) on the grounds
that it anchors each annotation to `bar.get_x() + bar.get_width() / 2`. That
anchors a label to A bar, not to ITS bar -- the pairing is the `zip`, and
reversing it puts every number neatly on top of the wrong bar with the stream
unchanged, inverting the synergy claim. Panel (a) is checked by position too.

fig25's panel (b) draws each score on the same ROW as its pair label, and that
row is the binding -- they are 226 to 305 points apart in x, so "beside" was never
true and reading them as two parallel lists was reading draw order again.
Mirroring the y coordinate in the generator's single drawing loop re-pairs
every score with a different pair, leaving the text stream identical; it
credited the flagship RSL3+FSP1i pair with a barely-additive 1.21x and gave
its 1.99x to FSP1i+HDACi, with this file green. Scores are now attributed by
shared row.

EVERY POSITIONAL CHECK HERE EXISTS BECAUSE THE STREAM ORDER OF SOME ELEMENT
TURNED OUT NOT TO BE ITS LAYOUT, in every figure and for every element class
here: bar annotations (anchored at a fixed offset on fig24, and to the wrong
bar of a reversed `zip` on fig25), tick labels (fig24's groups, fig25's bars,
fig26's timepoints), and panel (b)'s scores (one loop, mirrored y).

Tick labels are emitted in TICK-ARRAY order -- the order of the sequence handed
to `set_xticks`, not the left-to-right order of the resulting positions. That
is exactly why reversing the array leaves the text stream untouched while the
drawing inverts. An earlier version of this paragraph said "tick-location
order", which is the opposite and refutes its own conclusion: were it true, the
stream WOULD flip and the existing stream-order checks would already catch it.

Each of these was found after a previous round wrote a sentence claiming the
binding was already established, so the rule this file now follows is that
anything a reader takes as identifying a value -- a label, a unit, a position
in a series -- is compared by geometry, and the text stream is treated as
evidence of content only, never of arrangement.
"""
import hashlib
import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
FIG_DIR = REPO / "article/figures"
FIXTURES = REPO / "tests/fixtures"


def _drawn(stem):
    pymupdf = _reader()
    path = FIG_DIR / f"{stem}.pdf"
    assert path.exists(), f"{stem}.pdf is not committed"
    doc = pymupdf.open(path)
    try:
        return " ".join(" ".join(p.get_text().split()) for p in doc)
    finally:
        doc.close()


def _reader():
    """The PDF reader, or a FAILURE -- never a skip.

    Both readers used to `pytest.skip` when neither import worked, which turns
    the entire gate off silently: every figure in this file would report as
    fine on a machine that cannot open a PDF. A gate that disables itself when
    its instrument is missing is the failure it exists to prevent, and it is
    the same reasoning as the freshness suite's positive control failing rather
    than skipping when its pinned blob is unreachable. PyMuPDF is pinned in
    `requirements-lock.txt`, so this is a hard requirement, not an optional
    extra.
    """
    try:
        import pymupdf
        return pymupdf
    except ImportError:
        pass
    try:
        import fitz
        return fitz
    except ImportError as exc:
        raise AssertionError(
            "no PDF reader available, so none of these figures can be "
            "checked. PyMuPDF is pinned in requirements-lock.txt; install it "
            "rather than letting this gate pass by not running"
        ) from exc


def _filled_rects(stem):
    """Filled rectangles as `(x0, x1, y0, y1, fill)`, page 0.

    Words are not enough. fig24 draws its bar annotations at fixed offsets AND
    its bars at fixed offsets, so swapping the two `axA.bar` x offsets moves
    the RECTANGLES while every word stays exactly where it was -- the whole
    `(x, y, text)` dump is identical, in text and in geometry. The red
    "Hypoxic" bar then carries 91.9% and the blue "Normoxic" bar 87.8%,
    hypoxia raising the kill rate, and no check that reads words can ever see
    it. The drawing layer can: the rectangles carry their fill colour.
    """
    pymupdf = _reader()
    path = FIG_DIR / f"{stem}.pdf"
    assert path.exists(), f"{stem}.pdf is not committed"
    doc = pymupdf.open(path)
    try:
        out = []
        for d in doc[0].get_drawings():
            fill = d.get("fill")
            if fill is None:
                continue
            for item in d["items"]:
                if item[0] == "re":
                    r = item[1]
                    out.append((r.x0, r.x1, r.y0, r.y1,
                                tuple(round(c, 4) for c in fill)))
        return out
    finally:
        doc.close()


def _vertical_labels(words):
    """Rotated axis labels, as `x -> "the words bottom to top"`.

    A y-axis label is drawn rotated, so its words share an x and differ in y.
    Grouping on x recovers each panel's label separately, which counting
    occurrences cannot: fig24 gives both panels the same y-axis text, so
    `text.count(...) == 2` was satisfied while panel (a) had NO label at all
    (its words moved to the title) and again while panel (a) was relabelled and
    a second copy added to panel (b). A count cannot say where anything is.
    """
    cols = {}
    for x, y, t in words:
        cols.setdefault(round(x, 1), []).append((y, t))
    # BOTTOM TO TOP. A 90-degree rotated label puts its first word at the
    # LARGEST y, so reading in ascending y gives "(%) kill tumor Overall".
    return {x: " ".join(t for _, t in sorted(ws, reverse=True))
            for x, ws in cols.items() if len(ws) > 1}


def _panels(stem):
    """Each panel's x range, from the axis spines: `[(x0, x1), ...]` left to right.

    The spines are the long black horizontal strokes, and they are the only
    thing on the page that states where a panel actually ENDS. Earlier versions
    of this file guessed that boundary from panel titles (which broke when a
    title wrapped) and from the rightmost bar annotation (which is 53pt inside
    panel (a)'s real edge, so a stray label between the two panels counted as
    belonging to the second one).
    """
    pymupdf = _reader()
    doc = pymupdf.open(FIG_DIR / f"{stem}.pdf")
    try:
        spans = set()
        for d in doc[0].get_drawings():
            col = d.get("color")
            if col is None or tuple(round(c, 2) for c in col) != (0.0, 0.0, 0.0):
                continue
            for item in d["items"]:
                if item[0] != "l":
                    continue
                a, b = item[1], item[2]
                if abs(a.y - b.y) < 0.5 and abs(b.x - a.x) > 100:
                    spans.add((round(min(a.x, b.x), 1), round(max(a.x, b.x), 1)))
    finally:
        doc.close()
    assert spans, f"{stem}: no axis spines found, so panels cannot be located"
    return sorted(spans)


def _dashed_lines(stem):
    """Long DASHED horizontal strokes as `(y, x0, x1, colour)`.

    DASHED IS THE DISCRIMINATOR. An earlier version took every stroke with
    |dy| < 0.5 and |dx| > 100, which on fig24 also matches both axis spines,
    panel (a)'s legend frame, and -- the reason it was wrong rather than merely
    loose -- the near-flat segments of panel (b)'s own DATA CURVES, which run
    108.7pt with a dy of 0.05 to 0.09. Taking `min()` per colour over that set
    picked a curve segment whenever a curve happened to sit above the reference
    line, so deleting an `axhline` outright still passed: a curve answered for
    it, and the failure message named a y that was never a reference line.

    The two reference lines carry `dashes='[ 3.7 1.6 ] 0'`; every curve and
    spine is solid.
    """
    pymupdf = _reader()
    doc = pymupdf.open(FIG_DIR / f"{stem}.pdf")
    try:
        out = []
        for d in doc[0].get_drawings():
            dashes = (d.get("dashes") or "").strip()
            col = d.get("color")
            if col is None or not dashes or dashes == "[] 0":
                continue
            for item in d["items"]:
                if item[0] != "l":
                    continue
                a, b = item[1], item[2]
                if abs(a.y - b.y) < 0.5:
                    out.append((round(a.y, 2), min(a.x, b.x), max(a.x, b.x),
                                tuple(round(c, 4) for c in col)))
        return out
    finally:
        doc.close()


def _dashed_verticals(stem):
    """Long DASHED vertical strokes as `(x, y0, y1, colour)`."""
    pymupdf = _reader()
    doc = pymupdf.open(FIG_DIR / f"{stem}.pdf")
    try:
        out = []
        for d in doc[0].get_drawings():
            dashes = (d.get("dashes") or "").strip()
            col = d.get("color")
            if col is None or not dashes or dashes == "[] 0":
                continue
            for item in d["items"]:
                if item[0] != "l":
                    continue
                a, b = item[1], item[2]
                if abs(a.x - b.x) < 0.5 and abs(b.y - a.y) > 50:
                    out.append((round(a.x, 2), min(a.y, b.y), max(a.y, b.y),
                                tuple(round(c, 4) for c in col)))
        return out
    finally:
        doc.close()


def _hlines_any(stem):
    """Short horizontal stroked segments as `(x0, y, colour)` -- legend samples."""
    pymupdf = _reader()
    doc = pymupdf.open(FIG_DIR / f"{stem}.pdf")
    try:
        out = []
        for d in doc[0].get_drawings():
            col = d.get("color")
            if col is None:
                continue
            for item in d["items"]:
                if item[0] != "l":
                    continue
                a, b = item[1], item[2]
                if abs(a.y - b.y) < 0.5 and 3 < abs(b.x - a.x) < 45:
                    out.append((min(a.x, b.x), round(a.y, 2),
                                tuple(round(c, 4) for c in col)))
        return sorted(out)
    finally:
        doc.close()


def _stroke_points(stem):
    """Every stroked segment endpoint as `(x, y, colour)`."""
    pymupdf = _reader()
    doc = pymupdf.open(FIG_DIR / f"{stem}.pdf")
    try:
        out = []
        for d in doc[0].get_drawings():
            col = d.get("color")
            if col is None:
                continue
            c = tuple(round(v, 4) for v in col)
            for item in d["items"]:
                if item[0] != "l":
                    continue
                for pt in (item[1], item[2]):
                    out.append((pt.x, pt.y, c))
        return out
    finally:
        doc.close()


def _rows(words):
    """Group words by y, rounded, preserving `(x, text)` per row.

    Three checks need it: two panels drawn side by side share y coordinates for
    their tick labels and axis captions, so "same row" is what identifies a
    horizontal series -- and, when both panels land in one row, what makes a
    per-panel check have to work at word level rather than row level.
    """
    out = {}
    for x, y, t in words:
        out.setdefault(round(y, 1), []).append((x, t))
    return out


def _words(stem):
    """Every word on page 0 as `(x0, y0, text)`.

    Only page 0, and only the top-left corner of each box -- these figures are
    single-page and every check here needs left-to-right or same-row ordering,
    not extents. `_drawn` reads all pages; this deliberately does not.

    `_drawn` returns the text STREAM, which is draw order, not layout. For a
    matplotlib figure those two coincide only when the artist positions are
    themselves in order -- and fig24's are not: it places each annotation at a
    hand-computed `xi -/+ w/2` offset, entirely independent of the bar it is
    meant to sit over. Swapping the two offsets moves both numbers onto the
    wrong bars and leaves the text stream byte-identical, so no assertion over
    `_drawn` can see it. Geometry can.
    """
    pymupdf = _reader()
    path = FIG_DIR / f"{stem}.pdf"
    assert path.exists(), f"{stem}.pdf is not committed"
    doc = pymupdf.open(path)
    try:
        return [(w[0], w[1], w[4]) for w in doc[0].get_text("words")]
    finally:
        doc.close()


def _word_centres(stem):
    """`(centre_x, y0, text)` per word.

    The bar matcher needs the annotation's CENTRE, because the generator
    anchors these with `ha="center"`. Keying on the left edge only worked
    while the bars were much wider than their labels: at a bar width of 0.13 --
    a pure style constant, figure otherwise correct -- a label's left edge
    falls outside its own bar and the guard reported a data inversion.
    """
    pymupdf = _reader()
    path = FIG_DIR / f"{stem}.pdf"
    assert path.exists(), f"{stem}.pdf is not committed"
    doc = pymupdf.open(path)
    try:
        return [((w[0] + w[2]) / 2, w[1], w[4])
                for w in doc[0].get_text("words")]
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

    # BOTH BY POSITION, for the reasons fig24 needed it. The two checks above
    # read the text stream, and for panel (a) the stream is neither the bar
    # order nor the label order:
    #
    #  - Reversing the `zip` that pairs bars with values puts every number
    #    neatly on top of the WRONG bar -- RSL3 alone captioned 84.1% and the
    #    combination 40.0%, the synergy claim inverted -- with the stream
    #    unchanged. `bar.get_x() + bar.get_width()/2` anchors each label to A
    #    bar, not to ITS bar; the pairing is the zip, and nothing read it.
    #  - `set_xticks(...[::-1])` ALONE does not reverse this axis: it is
    #    categorical, the labels stay at exactly the same x, and the run above
    #    catches it at the `text.index` check. Adding `set_xticklabels(labels)`
    #    after it is the mutation that matters -- stream-identical, labels
    #    genuinely reversed on the page -- and only the row check below sees
    #    it. An earlier version of this bullet named the wrong mutation and
    #    claimed `text.index` could not see the one it named.
    wds = _words("fig25_bliss_synergy")
    drawn = sorted((x, t) for x, y, t in wds if re.fullmatch(r"\d+\.\d%", t))
    assert [t for _, t in drawn] == expected, (
        f"fig25's panel (a) bar annotations read {[t for _, t in drawn]} LEFT "
        f"TO RIGHT and bliss_synergy.json gives {expected}. Each number is "
        "anchored to the centre of some bar, so only position says which")
    # The labels are two-line, so the FIRST word of each carries its position.
    heads = [b.split()[0] for b in bars]
    rows = [r for r in _rows(wds).values()
            if sorted(t for _, t in r) == sorted(heads)]
    assert len(rows) == 1, (
        f"expected one row holding the four bar-label heads {heads}, found "
        f"{len(rows)}")
    assert [t for _, t in sorted(rows[0])] == heads, (
        f"fig25's panel (a) bar labels read "
        f"{[t for _, t in sorted(rows[0])]} LEFT TO RIGHT and its values are "
        f"compared as {heads}; the bars are labelled in the wrong order")
    # SCOPED TO PANEL (a). A bare `in text` was satisfied by panel (b), which
    # draws the same score as a bar label -- so replacing panel (a)'s
    # annotation with a wrong value, or deleting it, both passed.
    # IN PANEL (a), BY POSITION. `f"{score:.2f}× synergy" in text` is stream
    # adjacency, not scope: deleting panel (a)'s annotation outright and
    # rewording panel (b)'s title to begin `synergy, pairwise ...` recreates
    # the substring out of panel (b)'s bar label and title, and the check that
    # exists to require the annotation passes with the annotation gone. Panel
    # (b)'s pair labels are the rightmost thing panel (a) must be left of.
    words_a = _words("fig25_bliss_synergy")
    pair_x = [x for x, y, t in words_a
              if re.fullmatch(r"[^\s+]+\+[^\s+]+", t)]
    assert pair_x, "fig25 draws no pair labels, so panel (b) cannot be located"
    want = f"{fx['synergy_score']:.2f}×"
    hits = [(x, y) for x, y, t in words_a if t == want and x < min(pair_x)]
    assert hits, (
        f"fig25's panel (a) does not annotate {want} anywhere left of panel "
        f"(b)'s pair labels. Panel (b) drawing that number as a bar label is "
        "not the same claim")
    # The two words sit 0.87pt apart vertically (57.115 and 57.985), so an
    # exact row key splits them; the annotation is one visual line.
    ax, ay = hits[0]
    # THE WHOLE VISUAL LINE, both directions. Filtering to `x >= ax - 1`
    # discarded anything drawn to the LEFT of the number, so prefixing the
    # annotation with `no ` negated the panel's claim and passed the check
    # written to pin that exact annotation.
    line = sorted((x, t) for x, y, t in words_a
                  if abs(y - ay) < 2 and x < min(pair_x))
    assert " ".join(t for _, t in line) == f"{want} synergy", (
        f"fig25's panel (a) synergy annotation reads "
        f"{' '.join(t for _, t in line)!r}, expected {want + ' synergy'!r}")
    # BY POSITION, for the reason fig24's needed it: a containment check does
    # not say WHICH panel carries the label, or that any panel does.
    cols = _vertical_labels(_words("fig25_bliss_synergy"))
    kill_axes = sorted(x for x, lbl in cols.items()
                       if lbl == "Persister kill (%)")
    assert len(kill_axes) == 1, (
        "fig25 does not draw `Persister kill (%)` as a y-axis label exactly "
        f"once; vertical labels found: {sorted(cols.values())}. The four "
        "values compared above are drawn as percentages")
    # AND IT LABELS PANEL (a). Counting was the same defect fig24's version
    # had: removing panel (a)'s y-label and dropping one rotated copy past the
    # right edge of panel (b) left the count at one, 475pt from the bars it is
    # supposed to describe.
    fig25_panels = _panels("fig25_bliss_synergy")
    assert kill_axes[0] < fig25_panels[0][0], (
        f"fig25's `Persister kill (%)` label is at x={kill_axes[0]:.0f}, not "
        f"left of panel (a)'s axis at x={fig25_panels[0][0]:.0f}; it does not "
        "label the bars whose values are compared above")

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
    key = next((q for q in pairs if set(q.split("+")) == {"RSL3", "FSP1i"}), None)
    assert key, f"fig25's panel (b) no longer draws the RSL3/FSP1i pair: {pairs}"

    # PAIRED BY ROW, not by stream order. `dict(zip(pairs, scores))` binds the
    # nth name to the nth score, which is draw order -- and the generator draws
    # both in one loop, so mirroring the y coordinate (`len(scores) - 1 - i`)
    # re-pairs every score against a different drug pair while leaving the
    # stream identical. That published a figure crediting the 1.99x synergy to
    # FSP1i+HDACi and giving the flagship RSL3+FSP1i a barely-additive 1.21x,
    # with this test green. The docstring called this binding "direct"; it was
    # the same defect the fig17 guard was written for.
    #
    # A label and its score share a row and are 226 to 305 points apart in x,
    # so the row is the only thing tying them together.
    words = _words("fig25_bliss_synergy")
    # NO PANEL FILTER HERE, deliberately. Panel (b)'s pair labels are its
    # y-tick labels, drawn to the LEFT of its own axes -- far enough left to
    # fall on panel (a)'s side of any title-derived boundary, which is why the
    # split used for fig24 is wrong for this figure. The row IS the filter:
    # panel (a)'s synergy annotation is the only other `N.NN×` word and no
    # pair label shares its row, so it is excluded by having no partner
    # rather than by a coordinate.
    names = [(x, y, t) for x, y, t in words
             if re.fullmatch(r"[^\s+]+\+[^\s+]+", t)]
    vals = [(x, y, t) for x, y, t in words if re.fullmatch(r"\d+\.\d+×", t)]
    assert len(names) == len(pairs), (
        f"fig25 draws {len(names)} pair labels by geometry and {len(pairs)} "
        "in the text stream; the two readings disagree")
    by_row = {}
    for x, y, t in names:
        near = [v for v in vals if abs(v[1] - y) < 3]
        assert len(near) == 1, (
            f"{t} shares a row with {len(near)} scores, so no score can be "
            "attributed to it")
        by_row[t] = near[0][2]
    # THE SCORE'S DEFINITION. Panel (b)'s values are pinned and the caption
    # saying what they MEAN was not: flipping it to "(expected / observed)"
    # inverts every score on the panel and passed.
    # Panel (b)'s caption is horizontal, so position is its row. Bind it to
    # the panel by requiring it to sit to the RIGHT of every pair label, which
    # is what makes it panel (b)'s caption rather than a string anywhere.
    caption = [(x, y) for x, y, t in words if t == "score"]
    assert len(caption) == 1, (
        f"fig25 draws the word `score` {len(caption)} times; the panel (b) "
        "x-axis caption cannot be identified")
    # SCOPED TO PANEL (b) BY X. Panel (b)'s caption sits 0.35pt from panel
    # (a)'s second-line bar labels (`alone alone expected combination`), and
    # only `round(y, 1)` landing on 319.8 against 319.5 kept them apart. One
    # point of `axes.labelsize` closes that gap to 0.006pt, merges the two
    # groups into one row, and a CORRECT figure fails with a message accusing
    # the panel of drawing the reciprocal. The two groups are 76 points apart
    # in x, which is the difference that is actually real -- the same
    # nearest-neighbour lesson fig26's unit check already learned.
    _b0 = _panels("fig25_bliss_synergy")[-1][0]
    cap_row = [t for x, t in sorted(_rows(words)[round(caption[0][1], 1)])
               if x >= _b0 - 20]
    assert " ".join(cap_row) == "Bliss synergy score (observed / expected)", (
        f"fig25's panel (b) x axis reads {' '.join(cap_row)!r}; the values "
        "compared here are the observed-over-expected ratio, and the "
        "reciprocal is a different claim about every bar on the panel")
    fig25_panels = _panels("fig25_bliss_synergy")
    # THE ADDITIVE LINE, on panel (b)'s own scale. `axvline(1.0)` is what
    # turns a score into a synergy claim, and nothing read it: moving it to
    # 1.9 with its label unchanged makes two of the three pairs read
    # sub-additive, and only the fingerprint noticed. Panel (b)'s x tick
    # labels give the scale, so the line's position is checked in DATA units
    # rather than points.
    # CENTRES. A tick label is centred on its tick, so its left edge is offset
    # by half its own width -- using x0 put the additive line at 1.08.
    # `\d+(?:\.\d+)?`, NOT `\d+\.\d`. Panel (b)'s x limit is
    # `max(scores) * 1.25`, so once the top score passes ~3.2 matplotlib
    # switches to INTEGER tick labels and a one-decimal pattern matches none of
    # them. Measured: raising RSL3+FSP1i's synergy to 3.5 -- a stronger result,
    # the direction this figure argues for -- made a correct figure fail with
    # "0 numeric x tick labels". This file already paid for that exact mistake
    # once, in the pair-name regex a few blocks up.
    ticks = sorted((cx, float(t))
                   for cx, y, t in _word_centres("fig25_bliss_synergy")
                   if re.fullmatch(r"\d+(?:\.\d+)?", t)
                   and cx > fig25_panels[-1][0] - 20)
    assert len(ticks) >= 2, (
        f"fig25's panel (b) draws {len(ticks)} numeric x tick labels; its "
        "scale cannot be recovered")
    (tx0, tv0), (tx1, tv1) = ticks[0], ticks[-1]
    scale = (tx1 - tx0) / (tv1 - tv0)
    vlines = [x for x, y0, y1, col in _dashed_verticals("fig25_bliss_synergy")
              if fig25_panels[-1][0] <= x <= fig25_panels[-1][1]]
    assert len(vlines) == 1, (
        f"fig25's panel (b) draws {len(vlines)} dashed vertical lines, "
        "expected one additive threshold")
    drawn_at = (vlines[0] - tx0) / scale + tv0
    assert abs(drawn_at - 1.0) < 0.02, (
        f"fig25's additive threshold line is drawn at {drawn_at:.2f} on its "
        "own x scale, not 1.0. Every score on the panel is read against that "
        "line, so moving it changes which pairs read as synergistic")
    assert "additive" in text, (
        "fig25's panel (b) no longer labels its threshold line `additive`")

    assert caption[0][0] > max(x for x, _, t in names), (
        "fig25's score caption is not to the right of its pair labels, so it "
        "is not panel (b)'s x axis")
    assert by_row[key] == f"{fx['synergy_score']:.2f}×", (
        f"fig25 draws {key} on the same row as {by_row[key]} and the fixture "
        f"says {fx['synergy_score']:.2f}x -- the pair and the score beside it "
        "disagree with the data")


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
    # NO PANEL FILTER. The only `N.N%` words fig24 draws are its four bar
    # annotations and the collapse annotation's two; panel (b) plots curves
    # against a plain axis and draws none. A panel split was used here and is
    # gone: derived from the titles it false-failed the moment a title wrapped
    # onto two lines, and every version of it was a coordinate standing in for
    # a fact the words already carry. If panel (b) ever does draw a percentage
    # the comparison below fails loudly with the extra value listed, which is
    # the right way to find out.
    panel_a = _words("fig24_hypoxia_killcurve")
    arrow_y = [y for x, y, t in panel_a if t.strip() in {"→", "\u2192"}]
    assert len(arrow_y) == 1, (
        f"expected exactly one arrow in fig24 panel (a), found {len(arrow_y)}; "
        "the collapse annotation is what identifies the row to exclude")
    # THE COLLAPSE PAIR IS THE TWO FLANKING THE ARROW, not everything on its
    # row. Dropping the whole row false-failed on correct data: a bar whose top
    # lands near the collapse annotation puts its label on that row, and with
    # RSL3's uniform kill anywhere in roughly 56.6-59.2% the guard reported
    # three bars and blamed the figure for numbers it had drawn correctly.
    arrow_x = [x for x, y, t in panel_a
               if t.strip() in {"→", "\u2192"}][0]
    on_row = sorted((x, t) for x, y, t in panel_a
                    if re.fullmatch(r"\d+\.\d%", t) and abs(y - arrow_y[0]) <= 3)
    left = [w for w in on_row if w[0] < arrow_x]
    right = [w for w in on_row if w[0] > arrow_x]
    assert left and right, (
        "fig24's collapse annotation does not read `N.N% -> N.N%`; the two "
        f"percentages flanking its arrow are what identify it. row={on_row}")
    bars = sorted(((x, t) for x, y, t in panel_a
                   if re.fullmatch(r"\d+\.\d%", t)
                   and not (abs(y - arrow_y[0]) <= 3
                            and (x, t) in (left[-1], right[0]))),
                  key=lambda w: w[0])
    assert [t for _, t in bars] == expected, (
        f"fig24's bar annotations read {[t for _, t in bars]} LEFT TO RIGHT "
        f"and the data gives {expected}. The numbers are drawn at fixed "
        "offsets, so this is the only check that says each one sits over the "
        "bar it describes")

    # AND THE GROUP LABELS BY POSITION TOO. Pinning the numbers to their bars
    # while leaving the labels bound to the text stream binds nothing: tick
    # labels are emitted in tick-ARRAY order, so `set_xticks(x[::-1])`
    # moves `RSL3 (GPX4 inhibitor)` over the SDT bars and leaves the stream
    # untouched. The four numbers stay over their own bars and the figure now
    # captions 3.7%/0.1% as SDT -- the #790 inversion again, reached by a
    # different route, and the suite was green on it.
    #
    # The tick row is the one y where both names appear: `SDT` also occurs in
    # the footnote, so bottom-most is not the right rule.
    rows = {}
    for x, y, t in panel_a:
        if t in ("RSL3", "SDT"):
            rows.setdefault(round(y, 1), []).append((x, t))
    tick_rows = [r for r in rows.values()
                 if {t for _, t in r} == {"RSL3", "SDT"}]
    assert len(tick_rows) == 1, (
        f"expected exactly one row carrying both group labels, found "
        f"{len(tick_rows)}; that row is what identifies the tick labels")
    order = [t for _, t in sorted(tick_rows[0])]
    assert order == ["RSL3", "SDT"], (
        f"fig24's group labels read {order} LEFT TO RIGHT while its bars are "
        f"compared as {expected} -- RSL3's pair first. The labels sit over the "
        "wrong groups, so the figure captions each treatment with the other's "
        "numbers")
    # THE GROUP LABELS, IN ORDER. Without this every bar can be relabelled --
    # swapping the two group labels puts RSL3's 3.7%/0.1% under the SDT
    # heading and passed, which is the #790 defect this file exists to catch.
    labels = re.search(r"(RSL3 \(GPX4 inhibitor\)) (SDT \(exogenous ROS\))", text)
    assert labels, (
        "fig24's group labels are not RSL3 then SDT in the text stream")
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
    # THE BARS THEMSELVES, BY COLOUR AND POSITION. Everything above reads
    # words, and fig24 draws its bars at fixed offsets exactly as it draws its
    # annotations -- so swapping the two `axA.bar` x offsets slides the
    # RECTANGLES under the numbers and leaves the entire word dump identical,
    # text and geometry alike. The red "Hypoxic" bar then carries 91.9% and the
    # blue "Normoxic" bar 87.8%: hypoxia raising the kill rate, which is the
    # claim this figure exists to deny.
    #
    # The chain closed here is legend text -> swatch colour -> bar position.
    # The four bars are the series-coloured rectangles standing on the axis
    # baseline, which is the y1 they share and the legend swatches do not.
    rects = _filled_rects("fig24_hypoxia_killcurve")
    # EXACTLY ONE ENTRY EACH. The scan is page-wide, and panel (b) draws
    # `normoxic` and `(hypoxic)` in its own legend -- only their case and
    # parentheses keep them out. Capitalising them there rebinds panel (a)'s
    # colours through last-write-wins, so the count is checked rather than
    # assumed.
    for want in ("Normoxic", "Hypoxic"):
        n = sum(1 for _, _, t in panel_a if t == want)
        assert n == 1, (
            f"{want!r} is drawn {n} times on fig24; the legend entry it names "
            "can no longer be identified unambiguously")
    legend = {}
    for x, y, t in panel_a:
        if t in ("Normoxic", "Hypoxic"):
            # MATCH ON THE TOP EDGE, and let `sorted` pick. Measured on the
            # committed figure: each label's own swatch is 4.20pt below its
            # top, and the OTHER swatch is 8.80pt away for `Hypoxic` (the
            # entry above it) and 17.20pt for `Normoxic`. So the margin this
            # has to work in is 4.20 against 8.80 -- it is `sorted(...)`
            # taking the nearest that makes it right, not the `< 8` bound,
            # which only rejects far-away candidates. Comparing instead to the
            # rect's BOTTOM edge (y1) puts the correct swatch 9.80pt away,
            # outside the bound, and the entry resolves to nothing at all.
            # (Two earlier versions of this comment were wrong here: one said
            # both entries resolved to the same colour, the other quoted 13pt,
            # which is the label-to-label pitch and not any swatch distance.)
            near = sorted((abs(r[2] - y), r) for r in rects
                          if abs(r[2] - y) < 8 and r[1] - r[0] < 40)
            assert near, f"fig24's {t!r} legend entry has no swatch beside it"
            legend[t] = near[0][1][4]
    assert set(legend) == {"Normoxic", "Hypoxic"}, (
        f"fig24's panel (a) legend entries are {sorted(legend)}, expected "
        "Normoxic and Hypoxic")
    assert legend["Normoxic"] != legend["Hypoxic"], (
        "fig24 draws both legend swatches in the same colour, so colour "
        "cannot say which bar is which")
    series = [legend["Normoxic"], legend["Hypoxic"]]
    baseline = max(r[3] for r in rects if r[4] in series)
    bars_drawn = sorted((r for r in rects
                         if r[4] in series and abs(r[3] - baseline) < 0.5),
                        key=lambda r: r[0])
    # A ZERO-HEIGHT BAR DRAWS NO RECTANGLE AT ALL. matplotlib emits no `re`
    # item for it, so requiring four false-failed on correct data: setting
    # RSL3's hypoxic rate to 0.0 in the run AND the fixture together left every
    # value, label and legend assertion passing and died here on "found 3". A
    # zero rate is anticipated elsewhere in this figure -- the generator clamps
    # its collapse denominator at 0.01 for exactly that -- and today's value is
    # 0.1%. So each bar is matched to the annotation above it, and an
    # annotation with no bar under it is a bar of height zero, not a defect.
    assert bars_drawn, "fig24 draws no bars on its baseline at all"
    assert len(bars_drawn) <= 4, (
        f"fig24 draws {len(bars_drawn)} bars on its baseline, more than the "
        "four its data has; something else is being read as a bar")
    centres = {(round(y, 2), t): cx
               for cx, y, t in _word_centres("fig24_hypoxia_killcurve")}
    got = []
    for x, t in bars:
        cx = next((c for (yy, tt), c in centres.items()
                   if tt == t and abs(c - x) < 40), x)
        under = [r for r in bars_drawn if r[0] <= cx <= r[1]]
        if not under:
            got.append(None)                        # a zero-height bar
            continue
        assert len(under) == 1, (
            f"the annotation {t} at x={x:.0f} sits over {len(under)} bars")
        got.append("Normoxic" if under[0][4] == legend["Normoxic"] else "Hypoxic")
    want = ["Normoxic", "Hypoxic", "Normoxic", "Hypoxic"]
    paired = [(g, w) for g, w in zip(got, want) if g is not None]
    assert len(paired) >= 2, (
        f"only {len(paired)} of fig24's four annotations have a bar beneath "
        "them, so colour cannot say which series they belong to")
    assert [g for g, _ in paired] == [w for _, w in paired], (
        f"fig24's bars read {got} under its annotations, which are compared "
        f"as {expected} -- normoxic then hypoxic within each treatment. The "
        "numbers are sitting over bars of the other series")

    # THE UNIT THE BARS ARE IN. Four percentages are pinned above and nothing
    # said they were percentages: relabelling the axis "Overall tumor kill
    # (fraction surviving)" left every value in place and passed, inverting
    # what the bars mean. The same omission covered fig25's two axes.
    # BY POSITION. fig24 gives both panels the same y-axis label, so
    # `"..." in text` is satisfied by panel (b)'s copy -- the third time in
    # this file a bare containment check was answered by the other panel --
    # after fig25's synergy annotation and fig26's `(days)` caption. The
    # count that replaced it was no better: `== 2` held while panel (a) had no
    # y-axis label at all (the words moved into its title), and held again
    # while panel (a) was relabelled `(fraction surviving)` and a second copy
    # added to panel (b). A count says how many, never where.
    px = _panels("fig24_hypoxia_killcurve")
    assert len(px) >= 2, f"fig24 has {len(px)} panels by its spines, expected 2"
    columns = _vertical_labels(panel_a)
    kill_axes = sorted(x for x, lbl in columns.items()
                       if lbl == "Overall tumor kill (%)")
    assert len(kill_axes) == 2, (
        f"fig24 draws `Overall tumor kill (%)` as a y-axis label at "
        f"{len(kill_axes)} positions, expected one per panel. Columns found: "
        f"{sorted(columns.values())}. The values compared above are drawn as "
        "percentages")
    # AND ONE OF THEM LABELS PANEL (a). Counting columns is still counting:
    # stripping panel (a)'s y-label and dropping one stray rotated copy at the
    # page edge leaves two columns and panel (a) with no label at all, while
    # the four percentages it is supposed to describe stay where they are. A
    # y-axis label sits to the LEFT of the data it labels, so panel (a)'s must
    # be left of panel (a)'s own bars and panel (b)'s to their right.
    # BOUNDED ON THE PANEL SPINES. Using panel (a)'s rightmost bar annotation
    # as the divider put the boundary 53pt inside panel (a)'s real edge, so a
    # rotated label dropped ON TOP of panel (a)'s own bars counted as
    # labelling panel (b), and the message asserted a conclusion the bound
    # could not support.
    (a_x0, a_x1), (b_x0, b_x1) = px[0], px[-1]
    assert kill_axes[0] < a_x0, (
        f"fig24's leftmost `Overall tumor kill (%)` label is at "
        f"x={kill_axes[0]:.0f}, inside panel (a) rather than left of its axis "
        f"at x={a_x0:.0f}. Panel (a)'s values are drawn without the axis "
        "label that says what they are")
    # SLACK KEYED TO PANEL WIDTH, and a claim the bound can support. A y-axis
    # label is drawn just outside its own axes, and panel (b)'s sits only
    # 7.6pt clear of panel (a)'s right spine -- so a larger font or a narrower
    # figure slides it a couple of points INSIDE that spine while it still
    # plainly labels panel (b). The bound now allows an intrusion of 5% of a
    # panel's width and says what it actually detects; a label dropped in the
    # middle of panel (a) (round 9's mutation put one at x=308 against a
    # spine at 330) is still far outside it.
    slack = (a_x1 - a_x0) * 0.05
    assert a_x1 - slack < kill_axes[1] < b_x0, (
        f"fig24's second `Overall tumor kill (%)` label is at "
        f"x={kill_axes[1]:.0f}, which is not just outside panel (b)'s left "
        f"edge ({b_x0:.0f}) -- it sits more than {slack:.0f}pt inside panel "
        f"(a), whose axis ends at {a_x1:.0f}, so it is labelling panel (a)'s "
        "data rather than panel (b)'s")

    # PANEL (b)'S REFERENCE LINES, bound to PANEL (b)'S OWN LEGEND. Each
    # treatment's normoxic kill is drawn there as a dashed horizontal line,
    # and SDT's sits far above RSL3's -- that gap is the panel's claim.
    #
    # Two things this must not do, both of which an earlier version did. It
    # must not take any long flat stroke as a reference line, because panel
    # (b)'s own data curves flatten out and answered for a deleted `axhline`.
    # And it must not name the treatments by panel (a)'s series colours: those
    # happen to be the same two hex values, so recolouring panel (a) alone --
    # a style change leaving panel (b) correct -- failed with a message about
    # panel (b)'s legend, which it was not reading. Panel (b) captions its own
    # lines, and each caption has a dashed sample beside it.
    b_x0, b_x1 = px[-1]
    dashed = _dashed_lines("fig24_hypoxia_killcurve")
    ref_colour = {}
    for x, y, t in panel_a:
        if t != "normoxic":
            continue
        owner = [w for w in panel_a
                 if abs(w[1] - y) < 1 and w[0] < x and w[2] in ("SDT", "RSL3")]
        assert len(owner) == 1, (
            f"fig24's panel (b) legend row at y={y:.0f} names "
            f"{len(owner)} treatments")
        sample = [d for d in dashed
                  if 0 < d[0] - y < 10 and d[2] - d[1] < 40]
        assert sample, (
            f"fig24's `{owner[0][2]} normoxic` legend entry has no dashed "
            "sample beside it, so its colour cannot be read from the legend")
        ref_colour[owner[0][2]] = sample[0][3]
    assert set(ref_colour) == {"SDT", "RSL3"}, (
        f"fig24's panel (b) legend captions {sorted(ref_colour)} as normoxic "
        "references, expected SDT and RSL3")
    assert ref_colour["SDT"] != ref_colour["RSL3"], (
        "fig24's panel (b) draws both normoxic references in one colour")
    refs = {}
    for y, x0, x1, col in dashed:
        if x0 >= b_x0 - 1 and x1 <= b_x1 + 1 and (x1 - x0) > (b_x1 - b_x0) * 0.8:
            refs.setdefault(col, []).append(y)
    for name in ("SDT", "RSL3"):
        assert ref_colour[name] in refs, (
            f"fig24's panel (b) draws no full-width dashed {name} reference "
            "line in the colour its own legend gives that treatment")
        assert len(refs[ref_colour[name]]) == 1, (
            f"fig24's panel (b) draws {len(refs[ref_colour[name]])} {name} "
            "reference lines")
    sdt_y = refs[ref_colour["SDT"]][0]
    rsl3_y = refs[ref_colour["RSL3"]][0]
    assert sdt_y < rsl3_y, (
        f"fig24's panel (b) draws SDT's normoxic reference at y={sdt_y:.0f} "
        f"and RSL3's at y={rsl3_y:.0f}, so RSL3's reads as the higher kill "
        f"rate. The data gives SDT {kills['SDT'][0]:.1f}% against RSL3 "
        f"{kills['RSL3'][0]:.1f}% -- the two lines are exchanged")

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
        f"fig26's footnote does not say {n:,} cells/condition. The fixture "
        "carries n_cells as a top-level scalar, hoisted from the live run's "
        "rows, which the generator now reads rather than hardcoding")
    # THE UNIT. The timepoints above are the only numbers this test gates, and
    # this is the only text saying what they are measured in and from. Pinning
    # the values while leaving the axis caption free is the half-bound shape
    # the legend checks were added to close.
    wds = _words("fig26_vulnerability_window")

    # THE TIMEPOINTS BY POSITION. The run check above reads the text stream,
    # which is tick-ARRAY order: `set_xticks(x[::-1])` reverses the time axis
    # on both panels -- day 28 at the left, the window shading and the
    # `closes ~day 3` arrow pointing at the tick labelled 1, the collapse running
    # backwards -- and leaves the stream byte-identical. These labels are the
    # only data values this test gates, so stream order was gating nothing.
    #
    # Both panels share the tick row, so the row holds 2N labels: the left N
    # are panel (a)'s and the right N panel (b)'s, and each must read in
    # ascending order.
    # ROWS ACCUMULATED UPWARD FROM THE CAPTION until 2N labels are in hand.
    # Neither simpler rule survives rotation, and the two failures are
    # different, which is why this took three attempts. Measured, with
    # `set_xticklabels(labels, rotation=45)`:
    #
    #   unrotated          one row of 18, 14.7 above the caption
    #   rotated on axA     TWO rows of 9, at 13.5 and 14.7
    #   rotated on both    one row of 18, at 27.6 (the caption moves too)
    #
    # Requiring a single row of 2N fails the middle case; a fixed band of 25
    # fails the last. Taking rows nearest the caption first and stopping at 2N
    # handles all three, and stops well before the y-axis ticks, which are
    # single-element rows 32 or more above it.
    cap_y = min(y for x, y, t in wds if t == "post-chemotherapy")
    numeric = [(x, y, t) for x, y, t in wds
               if re.fullmatch(r"\d+(?:\.\d+)?", t) and 0 < cap_y - y < 60]
    by_row = {}
    for x, y, t in numeric:
        by_row.setdefault(round(y, 1), []).append((x, t))
    band = []
    for y in sorted(by_row, reverse=True):          # nearest the caption first
        if len(band) >= 2 * len(expected):
            break
        band.extend(by_row[y])
    assert len(band) == 2 * len(expected), (
        f"expected {2 * len(expected)} numeric tick labels above fig26's "
        f"x-axis caption (both panels share the axis), found {len(band)}: "
        f"{sorted(band)}")
    ticks = [t for _, t in sorted(band)]
    assert ticks[:len(expected)] == expected and ticks[len(expected):] == expected, (
        f"fig26's tick labels read {ticks} LEFT TO RIGHT; each panel should "
        f"read {expected} in ascending order. The axis may be reversed, which "
        "the text stream cannot show")

    # THE UNIT, PER PANEL AND BY WORDS. Both x-labels are drawn at the SAME y,
    # so grouping the caption by row merges them into one string -- an earlier
    # version of this check did exactly that and `"(days)" in cap` was
    # satisfied by panel (b)'s copy while panel (a) said hours. Split the row
    # at each `Time` and require every caption to end in days.
    # EACH CAPTION'S UNIT IS ITS NEAREST ONE, by distance rather than by row.
    # Both panels' captions are drawn at the same y, so grouping by row merges
    # them and a containment check is answered by the other panel's copy --
    # that shipped once already, as the fix for the round before it. Splitting
    # the row at each `Time` closed that and false-failed the moment the
    # caption wrapped onto two lines. Nearest-unit handles both: the unit is
    # adjacent on one row, and directly below when wrapped.
    units = [(x, y, t) for x, y, t in wds if re.fullmatch(r"\(\w+\)", t)]
    subjects = [(x, y, t) for x, y, t in wds if t == "post-chemotherapy"]
    assert len(subjects) == 2, (
        f"fig26 draws {len(subjects)} x-axis captions, expected one per panel")
    for x, y, _ in subjects:
        assert units, "fig26 draws no parenthesised unit beside its x axis"
        nearest = min(units, key=lambda u: (u[0] - x) ** 2 + (u[1] - y) ** 2)
        assert nearest[2] == "(days)", (
            f"the x-axis caption at x={x:.0f} is nearest the unit "
            f"{nearest[2]!r}. The ticks are the fixture's timepoint_days, so a "
            "panel saying anything else declares a unit its own numbers are "
            "not in")
    # THE CURVE LEGEND, BY COLOUR AND POSITION. The stream-order check above
    # binds nothing on its own: `axA.legend(handles[::-1], labels, ...)` gives
    # the SDT label the RSL3 line's colour and marker -- "days for RSL3, weeks
    # for SDT" inverted -- and passed every semantic check here, caught only by
    # the fingerprint as an unexplained change. fig24's two legends were both
    # given geometric bindings for exactly this reason; this one was left out.
    #
    # Each legend entry has its line sample ~6pt below its text, and the claim
    # is which curve is still high at the last timepoint.
    p26 = _panels("fig26_vulnerability_window")
    a_lo, a_hi = p26[0]
    segs = []
    doc_lines = _hlines_any("fig26_vulnerability_window")
    legend_col = {}
    for x, y, t in wds:
        if t not in ("SDT", "RSL3") or not (a_lo <= x <= a_hi):
            continue
        sample = [seg for seg in doc_lines
                  if 0 < seg[1] - y < 10 and abs(seg[0] - x) < 60]
        if sample:
            legend_col.setdefault(t, sample[0][2])
    assert set(legend_col) == {"SDT", "RSL3"}, (
        f"fig26's panel (a) legend captions {sorted(legend_col)}, expected "
        "SDT and RSL3, each with a line sample beside it")
    assert legend_col["SDT"] != legend_col["RSL3"], (
        "fig26 draws both panel (a) legend samples in one colour")
    ends = {}
    for x, y, col in _stroke_points("fig26_vulnerability_window"):
        if not (a_lo <= x <= a_hi) or col not in legend_col.values():
            continue
        if x > ends.get(col, (-1, 0))[0]:
            ends[col] = (x, y)
    for name in ("SDT", "RSL3"):
        assert legend_col[name] in ends, (
            f"fig26's panel (a) draws no curve in the colour its legend gives "
            f"{name}")
    sdt_end = ends[legend_col["SDT"]][1]
    rsl3_end = ends[legend_col["RSL3"]][1]
    assert sdt_end < rsl3_end, (
        f"fig26's panel (a) ends with the SDT curve at y={sdt_end:.0f} and "
        f"RSL3 at y={rsl3_end:.0f}, so RSL3 is the one still killing at the "
        "last timepoint. The figure's title is `days for RSL3, weeks for "
        "SDT`; the curves or the legend are exchanged")

    assert "closes ~day 3" in text, (
        "fig26's window annotation is gone. NOTE this is a hardcoded string in "
        "the generator, not a derived one -- if the window moved, the "
        "annotation would still read day 3 and this assertion would still "
        "pass, so it is a presence check and not a data guard")


# ---------------------------------------------------------------------------
# The backstop
# ---------------------------------------------------------------------------

# Recorded fingerprints of what each figure DRAWS: every text span with its
# box, colour, size, font and opacity, and every path item with its geometry
# and every paint property, in paint order.
#
# WHY THIS EXISTS. Successive adversarial rounds each found a different
# element whose ARRANGEMENT was unchecked -- annotations at fixed offsets, annotations
# anchored to the wrong bar of a reversed zip, three sets of tick labels, a
# legend, a label/value list, and finally the bar rectangles themselves, which
# move under their own labels leaving every word untouched. Each round closed
# the element it found and the next round found another. Closing them one at a
# time is not converging, because the list is not the mutation SPACE.
#
# It is much wider than the semantic checks and it is still not "the space".
# FOUR rounds have now found elements it did not cover, each after it had been
# widened for the last one: it read only rectangle items, missing every line
# and curve; then no paint property beyond colour, so a confidence band turned
# opaque enough to hide the curves left it unmoved; then no text opacity and
# no paint order, so a value could be transparent or buried behind a bar; and
# then no clip, so a path could be in the stream and scissored off the page.
#
# It now reads every path item with its geometry, every paint property, the
# clip in force and the nesting level, transparency groups, and text spans
# with box, colour, size, font and opacity -- in paint order. That is a claim
# about what it READS, not a promise that nothing can slip past it, and the
# history above is the reason to state it that way. Treat a green fingerprint
# as "nothing I know how to look at moved", and keep adding semantic checks,
# which say what a change MEANS.
#
# It is deliberately a BACKSTOP, not a replacement. It says "something in the
# drawing moved and nothing above noticed"; it cannot say what the change
# means. That is the semantic checks' job, and they stay. When this fails,
# either add the check that explains the change or update the hash on purpose.
#
# PORTABLE BY CONSTRUCTION: computed from the COMMITTED PDF, not from a
# regeneration, so every platform hashes the same bytes with the same pinned
# reader. That is why it can run in CI where the figures cannot be rebuilt.
DRAWING_FINGERPRINTS = {
    "fig24_hypoxia_killcurve": "7822f5da7e31d9d8",
    "fig25_bliss_synergy": "0e6e769553b76f80",
    "fig26_vulnerability_window": "456a227169e9af16",
}


def _fingerprint(stem):
    """A hash of everything drawn: word positions and path geometry.

    EVERY ITEM TYPE AND EVERY PAINT PROPERTY, over spans rather than words. The first version hashed only `re`
    items -- 14 of fig24's 112 drawing items, and 24 of fig26's 277, where the
    plotted CURVES are the data. Swapping the two `axhline` VALUES in fig24's
    panel (b), keeping every label and colour, moved the two reference lines
    (y=72.06 and y=285.91 exchanged), inverted what panel (b) claims, and
    passed the whole file. A backstop that reads an eighth of the drawing is
    not a backstop.

    COLOUR IS PART OF IT, and has to be. Dropping it looked attractive because
    PyMuPDF 1.24.14 reports fig25's `#999999` bar as 0.664 where 1.28.2 reports
    0.5999 -- but 1.24.14 is below the `>=1.28.0` floor requirements.txt
    declares, so that divergence is outside the supported range. And without
    colour the hash is blind to the mutation that motivated widening it in the
    first place: swapping which VALUE each of fig24's two `axhline` calls draws
    exchanges the two reference lines, and since the pair of line geometries is
    unchanged as a SET, only the colour attached to each one moves. Geometry
    alone cannot see a swap; the swap is the defect.

    The hash is therefore bound to the reader as well as the file. That is
    stated in the failure message rather than engineered away, because the
    alternative -- a backstop that cannot see a swap -- is worse.
    """
    pymupdf = _reader()
    path = FIG_DIR / f"{stem}.pdf"
    assert path.exists(), f"{stem}.pdf is not committed"
    doc = pymupdf.open(path)
    try:
        page = doc[0]
        # SPANS, NOT WORDS, because a word carries no colour or size. Drawing
        # fig24's four bar annotations in white leaves them invisible on a
        # white ground -- the four values this whole file exists to gate gone
        # from the artifact -- while every word is still at its old position
        # and every check, this one included, passed.
        spans = []
        for block in page.get_text("dict")["blocks"]:
            for line in block.get("lines", ()):
                for span in line["spans"]:
                    # ALPHA TOO. Colour was closed and the property next to
                    # it reaching the identical outcome was not: `alpha=0.0`
                    # on fig24's four bar annotations renders 0 non-white
                    # pixels where the committed figure renders 1176, with
                    # every span's bbox, colour, size and font unchanged and
                    # all sixteen tests green. Paths already had their two
                    # opacities hashed; text did not.
                    spans.append((tuple(round(v, 2) for v in span["bbox"]),
                                  span["text"], span["color"],
                                  round(span["size"], 2), span.get("font"),
                                  span.get("alpha")))
        # DOCUMENT ORDER, NOT SORTED. Sorting made the hash blind to what is
        # drawn on top of what: `zorder=0` on fig25's additive threshold line
        # hides 77% of it behind the bars it is the reference for (709 visible
        # scanlines down to 164) without moving a single coordinate. Order in
        # the content stream IS paint order, and it is deterministic for a
        # given file, so keeping it costs nothing and closes occlusion.
        words = spans
        items = []
        # extended=True ADDS THE CLIP STATE. Without it a path can be in the
        # content stream and invisible on the page: setting a clip box on
        # panel (b)'s four kill curves erases 60% of them -- 8507 pixels at
        # 200 dpi -- while every coordinate and paint property stays exactly
        # as it was, so the hash did not move. The extended list carries
        # `clip` entries with the scissor rect in force and a `level` for each
        # entry, which is what says WHICH clip applies to which path.
        for d in page.get_drawings(extended=True):
            if d.get("type") == "group":
                # A transparency group: no path items, but its bbox and
                # blend/opacity are part of what the page shows.
                r = d.get("rect")
                items.append(("group", d.get("level"),
                              None if r is None else
                              tuple(round(v, 2) for v in (r.x0, r.y0, r.x1, r.y1)),
                              d.get("blendmode"), d.get("opacity")))
                continue
            if d.get("type") == "clip":
                sc = d.get("scissor")
                items.append(("clip", d.get("level"),
                              None if sc is None else
                              tuple(round(v, 2) for v in
                                    (sc.x0, sc.y0, sc.x1, sc.y1))))
                continue
            # EVERY PAINT PROPERTY, not just the two colours. Width, dash
            # pattern and opacity are all invisible to a colour-only hash, and
            # each is enough on its own to change what the figure SAYS:
            # turning fig26's confidence bands from alpha 0.15 to 0.95 covers
            # the plotted curves entirely -- the data of the figure -- and the
            # hash did not move.
            paint = (
                tuple(None if d.get(k) is None
                      else tuple(round(c, 2) for c in d[k])
                      for k in ("fill", "color"))
                + tuple(round(d[k], 3) if isinstance(d.get(k), float)
                        else d.get(k)
                        for k in ("width", "dashes", "stroke_opacity",
                                  "fill_opacity", "even_odd", "closePath")))
            for item in d["items"]:
                kind, rest = item[0], item[1:]
                coords = []
                for part in rest:
                    if hasattr(part, "ul"):                  # Quad
                        for corner in (part.ul, part.ur, part.ll, part.lr):
                            coords.extend([round(corner.x, 2),
                                           round(corner.y, 2)])
                    elif hasattr(part, "x0"):                # Rect
                        coords.extend(round(getattr(part, a), 2)
                                      for a in ("x0", "y0", "x1", "y1"))
                    elif hasattr(part, "x"):                 # Point
                        coords.extend([round(part.x, 2), round(part.y, 2)])
                    else:
                        coords.append(part)                  # scalar (radius)
                # AN UNRECOGNISED SHAPE MUST NOT HASH AS NOTHING. The first
                # version read x0/y0/x1/y1 or x/y and fell through silently
                # otherwise -- a `qu` item has neither, so a quad would have
                # been recorded as an empty tuple and any change to it would
                # have been invisible while the docstring claimed every item
                # type was covered. None of these three figures draws one
                # today, which is exactly why it would have gone unnoticed.
                assert coords, (
                    f"{stem}: a `{kind}` drawing item produced no coordinates, "
                    "so it would be hashed as nothing. Teach this extractor "
                    "the shape rather than letting it hash an empty tuple")
                items.append((kind, tuple(coords), paint, d.get("level")))
    finally:
        doc.close()
    blob = json.dumps({"words": words, "items": items}, sort_keys=True)
    return hashlib.sha256(blob.encode()).hexdigest()[:16]


@pytest.mark.parametrize("stem", sorted(DRAWING_FINGERPRINTS))
def test_nothing_in_the_drawing_moved_unnoticed(stem):
    """The check that does not need to know what to look for.

    Every assertion above had to be written for a specific element after a
    specific defect. This one covers the elements nobody has thought about
    yet, which is where every review round has found its defects.
    """
    got = _fingerprint(stem)
    reader = _reader()
    version = getattr(reader, "version", ("?",))[0]
    assert got == DRAWING_FINGERPRINTS[stem], (
        f"{stem}'s drawing changed: fingerprint {got}, recorded "
        f"{DRAWING_FINGERPRINTS[stem]}. Something moved, recoloured, appeared "
        "or disappeared, and no check above noticed -- which is the case this "
        "exists for. Work out what changed and either add the assertion that "
        "explains it, or update the hash here deliberately if the new drawing "
        "is correct. Do not update it to make the suite green.\n"
        f"\nPDF reader: PyMuPDF {version}. The hash is read out of the "
        "committed file, so the same bytes give the same value everywhere "
        "with the same reader -- but the reader is part of it. CI installs "
        "requirements-lock.txt, which pins 1.28.2; requirements.txt allows "
        ">=1.28.0, so a newer one here can move this hash without the figure "
        "having changed at all. Check the version before hunting the drawing.")
