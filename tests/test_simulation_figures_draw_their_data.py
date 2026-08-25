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

    Takes `(x0, x1, y, text)` tuples -- see `_boxed_row_words`.

    A y-axis label is drawn rotated, so its words share an x and differ in y.
    Grouping on x recovers each panel's label separately, which counting
    occurrences cannot: fig24 gives both panels the same y-axis text, so
    `text.count(...) == 2` was satisfied while panel (a) had NO label at all
    (its words moved to the title) and again while panel (a) was relabelled and
    a second copy added to panel (b). A count cannot say where anything is.
    """
    # KEYED ON THE FULL X EXTENT. Grouping on x0 alone let a HORIZONTAL word
    # join a rotated column whenever it happened to start at the same x: at
    # one style setting the footnote's word `upper` shares x0=489.629 with
    # panel (b)'s label and the column read "upper Overall tumor kill (%)",
    # failing a correct figure. A rotated label's words all have the same
    # x-extent -- their width is the text HEIGHT -- so x1 separates them from
    # a horizontal word of a DIFFERENT width.
    #
    # It is not a general separation: any words of equal width stacked at the
    # same x group together, whether or not they are a rotated label. Measured
    # on the committed fig24, this returns SIX columns -- the two real labels
    # at (3.0, 20.5) and (337.8, 355.3), the two y-tick stacks `20 40 60 80`
    # at (30.2, 44.2) and (365.0, 379.0), and panel (b)'s legend read
    # vertically as `SDT SDT` and `RSL3 RSL3`, both starting at x0 595.6.
    # That is harmless because callers select by the label TEXT.
    #
    # It is still why the returned key is the full extent rather than x0: the
    # two legend columns share an x0 with each other, so collapsing onto x0
    # would merge them last-write-wins. (Two earlier versions of this note
    # said "one such column" and then "TWO", and neither counted the tick
    # stacks. The conclusion held; the count did not.)
    cols = {}
    for x0, x1, y, t in words:
        cols.setdefault((round(x0, 1), round(x1, 1)), []).append((y, t))
    # BOTTOM TO TOP. A 90-degree rotated label puts its first word at the
    # LARGEST y, so reading in ascending y gives "(%) kill tumor Overall".
    # Returns text AND the column's own vertical extent, so a caller checking
    # where the label sits does not have to re-find its words by x0 -- which
    # re-opens the very hole this function was fixed for.
    # KEYED ON THE FULL EXTENT ON THE WAY OUT TOO. Grouping on (x0, x1) and
    # then returning a dict keyed on x0 alone put two columns that share a
    # rounded x0 back into one slot, last-write-wins -- undoing the grouping
    # it had just done. fig24 already carries a spurious horizontal column
    # (`RSL3 RSL3`, its panel (b) legend, at x=595.6), so this is not
    # hypothetical; x1 keeps it separate only if x1 survives.
    return {key: (" ".join(t for _, t in sorted(ws, reverse=True)),
                  min(y for y, _ in ws), max(y for y, _ in ws))
            for key, ws in cols.items() if len(ws) > 1}


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
            # ANY GREY, NOT ONLY BLACK. Requiring exactly (0,0,0) meant
            # `axes.edgecolor="0.3"` -- seaborn's default and most journal
            # templates -- aborted every geometric reader in this file at
            # once, on figures that render perfectly. An axis frame is
            # greyscale; a data series is not.
            if col is None or max(col) - min(col) > 0.1 or min(col) > 0.75:
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


def _panel_y(stem):
    """The vertical extent of the plotting area, from the axis spines."""
    pymupdf = _reader()
    doc = pymupdf.open(FIG_DIR / f"{stem}.pdf")
    try:
        ys = set()
        for d in doc[0].get_drawings():
            col = d.get("color")
            # ANY GREY, NOT ONLY BLACK. Requiring exactly (0,0,0) meant
            # `axes.edgecolor="0.3"` -- seaborn's default and most journal
            # templates -- aborted every geometric reader in this file at
            # once, on figures that render perfectly. An axis frame is
            # greyscale; a data series is not.
            if col is None or max(col) - min(col) > 0.1 or min(col) > 0.75:
                continue
            for item in d["items"]:
                if item[0] != "l":
                    continue
                a, b = item[1], item[2]
                if abs(a.y - b.y) < 0.5 and abs(b.x - a.x) > 100:
                    ys.add(round(a.y, 2))
    finally:
        doc.close()
    # ONE SPINE IS ENOUGH WHEN THE TICKS GIVE THE OTHER END. `axes.spines.top:
    # False` is the commonest publication restyle and it removed a whole row,
    # aborting fig24 on a correct figure. The tick marks span the axis, so
    # they close the range when a spine is gone.
    if len(ys) < 2:
        tick_ys = {round(a.y, 2) for _, _, _, _, _, a, _ in _axis_ticks(stem)}
        ys |= tick_ys
    assert len(ys) >= 2, (
        f"{stem}: fewer than two spine rows and no tick marks to close the "
        "vertical extent")
    return min(ys), max(ys)


def _axis_ticks_vertical(stem):
    """Short black tick marks on a horizontal axis (vertical strokes)."""
    pymupdf = _reader()
    doc = pymupdf.open(FIG_DIR / f"{stem}.pdf")
    try:
        out = []
        for d in doc[0].get_drawings():
            col = d.get("color")
            # ANY GREY, NOT ONLY BLACK. Requiring exactly (0,0,0) meant
            # `axes.edgecolor="0.3"` -- seaborn's default and most journal
            # templates -- aborted every geometric reader in this file at
            # once, on figures that render perfectly. An axis frame is
            # greyscale; a data series is not.
            if col is None or max(col) - min(col) > 0.1 or min(col) > 0.75:
                continue
            for item in d["items"]:
                if item[0] != "l":
                    continue
                a, b = item[1], item[2]
                if abs(a.x - b.x) < 0.2 and 0.5 <= abs(b.y - a.y) <= 20:
                    out.append((a.x, min(a.y, b.y), b.x, max(a.y, b.y),
                                tuple(round(c, 4) for c in col), a, b))
        return out
    finally:
        doc.close()


def _axis_ticks(stem):
    """Short black axis tick marks as `(x0, y0, x1, y1, colour, p1, p2)`."""
    pymupdf = _reader()
    doc = pymupdf.open(FIG_DIR / f"{stem}.pdf")
    try:
        out = []
        for d in doc[0].get_drawings():
            col = d.get("color")
            # ANY GREY, NOT ONLY BLACK. Requiring exactly (0,0,0) meant
            # `axes.edgecolor="0.3"` -- seaborn's default and most journal
            # templates -- aborted every geometric reader in this file at
            # once, on figures that render perfectly. An axis frame is
            # greyscale; a data series is not.
            if col is None or max(col) - min(col) > 0.1 or min(col) > 0.75:
                continue
            for item in d["items"]:
                if item[0] != "l":
                    continue
                a, b = item[1], item[2]
                # 0.5 TO 20pt. The window was `1 < len < 8`, which failed at
                # both ends on correct figures: `xtick.major.size: 8` and
                # `1.0` each left zero ticks found, silently disabling the
                # scale readings that depend on them. matplotlib's default is
                # 3.5 and journal styles routinely set 1 or 8.
                if abs(a.y - b.y) < 0.2 and 0.5 <= abs(b.x - a.x) <= 20:
                    out.append((min(a.x, b.x), a.y, max(a.x, b.x), b.y,
                                tuple(round(c, 4) for c in col), a, b))
        return out
    finally:
        doc.close()


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
    """Every STROKE-ONLY segment endpoint as `(x, y, colour)`.

    `fill is None` IS THE CURVE FILTER. A `fill_between` band is reported with
    a stroke colour too (type `fs`, 81 of fig26's 103 coloured drawings), and
    it is drawn in the series colour and spans the same x -- so taking every
    coloured segment made "where does this curve end" read the BAND. On the
    committed figure both endpoints already came from bands, and swapping only
    the two `axA.plot` data arguments -- so the curve labelled SDT traces
    RSL3's data and collapses to zero by day 28, inverting the panel's whole
    claim -- passed, because the untouched bands still ended where the
    assertion expected. The mirror image was worse: swapping only the two
    `fill_between` calls rejected a correct pair of curves and blamed the
    curves for it.
    """
    pymupdf = _reader()
    doc = pymupdf.open(FIG_DIR / f"{stem}.pdf")
    try:
        out = []
        for d in doc[0].get_drawings():
            col = d.get("color")
            if col is None or d.get("fill") is not None:
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


def _word_bboxes(stem):
    """Every word as `(x0, y0, x1, y1, text)` -- the full box.

    A tick label is centred on its tick, so mapping a drawn y to a data value
    needs the label's vertical CENTRE. Using its top edge put fig24's SDT
    reference at 88.6% against a true 91.9%.
    """
    pymupdf = _reader()
    doc = pymupdf.open(FIG_DIR / f"{stem}.pdf")
    try:
        return [(w[0], w[1], w[2], w[3], w[4])
                for w in doc[0].get_text("words")]
    finally:
        doc.close()


def _boxed_row_words(stem):
    """Every word as `(x0, x1, y0, text)`."""
    pymupdf = _reader()
    doc = pymupdf.open(FIG_DIR / f"{stem}.pdf")
    try:
        return [(w[0], w[2], w[1], w[4]) for w in doc[0].get_text("words")]
    finally:
        doc.close()


def _boxed_rows(stem):
    """Rows of `(x0, x1, text)` keyed by rounded y, with word EXTENTS.

    `_rows` keeps only each word's left edge, so a "gap" computed from it is
    start-to-start and carries the width of the preceding word. That inflated
    the typical gap enough to swallow the real one: at label size 13 the
    between-group gap is 63.3pt and a start-to-start median put the threshold
    just above it.
    """
    pymupdf = _reader()
    doc = pymupdf.open(FIG_DIR / f"{stem}.pdf")
    try:
        ws = [(w[1], w[0], w[2], w[4]) for w in doc[0].get_text("words")]
        return {round(min(y for y, _, _, _ in row), 1):
                sorted((x0, x1, t) for _, x0, x1, t in row)
                for row in _cluster_rows(ws, key=lambda w: w[0])}
    finally:
        doc.close()


def _cluster_rows(items, key, tol=1.0):
    """Group `items` into rows by `key(item)`, within `tol` of each other.

    ROUNDING IS NOT GROUPING. `round(y, 1)` puts two words on different rows
    whenever they straddle a boundary, however close they are -- and glyph
    extents make that a coin flip, not a property of the figure. Measured:
    fig25's four bar-label heads sit at y0 339.0485, 339.0485, 339.0504,
    339.0504 AT FIGURE HEIGHT 5 -- three of nine sampled heights hit this --
    a spread of 0.0019pt, and `round` split them into 339.0 and 339.1 because
    `RSL3`/`FSP1i` are cap-height while `Bliss`/`Observed` carry ascenders.
    On the committed figure they are at 306.2885/306.2905, which `round`
    happens NOT to split; an earlier version of this note quoted the height-5
    coordinates without saying so, which reads as though the committed figure
    were affected. A correct figure then failed being told its bars were
    mislabelled. This file had already diagnosed the same knife edge for the
    caption row and fixed it only there.
    """
    out = []
    for item in sorted(items, key=key):
        if out and abs(key(item) - key(out[-1][-1])) <= tol:
            out[-1].append(item)
        else:
            out.append([item])
    return out


def _rows(words):
    """Group words into rows by y, tolerantly, as `representative_y -> [(x, text)]`.

    Two panels drawn side by side share y coordinates for their tick labels and
    axis captions, so "same row" is what identifies a horizontal series -- and,
    when both panels land in one row, what makes a per-panel check have to work
    at word level rather than row level.
    """
    return {round(min(y for _, y, _ in row), 1): [(x, t) for x, _, t in row]
            for row in _cluster_rows(words, key=lambda w: w[1])}


def _row_containing(rows, y, tol=1.5):
    """The row whose key is nearest `y`, or raise if none is close."""
    near = sorted(rows, key=lambda k: abs(k - y))
    assert near and abs(near[0] - y) <= tol, (
        f"no row within {tol}pt of y={y:.2f}; rows are at "
        f"{sorted(rows)[:8]}")
    return rows[near[0]]


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
    # SCOPED TO PANEL (a) BY X, not by requiring the row to hold only these
    # four words. Panel (b)'s x tick labels sit just above the heads, and the
    # separation is 0.1 x `xtick.labelsize` -- 1.101pt at the committed size,
    # so tolerance clustering at 1.0 cleared it by 0.101pt and merged the two
    # groups at tick size 9, at `font.size` 8, and under Helvetica (0.164pt).
    # The previous `round(y, 1)` keying had the same knife edge in a different
    # place; clustering moved it rather than removing it. The two groups are
    # ~70pt apart in x and on opposite sides of a panel spine, which is the
    # separation that does not depend on a font.
    heads = [b.split()[0] for b in bars]
    a_hi = _panels("fig25_bliss_synergy")[0][1]
    rows = [[w for w in r if w[0] <= a_hi] for r in _rows(wds).values()]
    rows = [r for r in rows if sorted(t for _, t in r) == sorted(heads)]
    assert len(rows) == 1, (
        f"expected one row of panel (a) holding the four bar-label heads "
        f"{heads}, found {len(rows)}")
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
    # BOUNDED ON PANEL (b)'S SPINE, not on its pair labels. Those labels are
    # y-ticks drawn OUTSIDE the spine, so a longer drug name pushes them left
    # towards panel (a)'s annotation: with real compound names the margin goes
    # from 81pt today to 1.0pt for `Liproxstatin-1+Ferrostatin-1` and to below
    # zero for `Liproxstatin-1+Deferoxamine`, failing a correct figure. The
    # spine moves the other way -- longer tick labels shrink panel (b), so the
    # bound gains room exactly when the labels need it. This file has already
    # been fixed twice for assuming today's pair names.
    b_spine = _panels("fig25_bliss_synergy")[-1][0]
    hits = [(x, y) for x, y, t in words_a if t == want and x < b_spine]
    assert hits, (
        f"fig25's panel (a) does not annotate {want} anywhere left of panel "
        f"(b)'s axis at x={b_spine:.0f}. Panel (b) drawing that number as a "
        "bar label is not the same claim")
    # The two words sit 0.87pt apart vertically (57.115 and 57.985), so an
    # exact row key splits them; the annotation is one visual line.
    _ax, ay = hits[0]
    ay_mid = next((y0 + y1) / 2
                  for x0, y0, x1, y1, t in _word_bboxes("fig25_bliss_synergy")
                  if t == want and abs(x0 - _ax) < 0.01 and abs(y0 - ay) < 0.01)
    # THE WHOLE VISUAL LINE, both directions. Filtering to `x >= ax - 1`
    # discarded anything drawn to the LEFT of the number, so prefixing the
    # annotation with `no ` negated the panel's claim and passed the check
    # written to pin that exact annotation.
    # THE CONTIGUOUS RUN CONTAINING THE NUMBER. Bounding this by panel (b)'s
    # pair labels held only while those labels stayed put: at a figure width
    # of 8 or less they move left of the annotation's own trailing word, and a
    # correct figure failed reading `1.99×` with its `synergy` dropped. The
    # annotation is one visual run, so the run is what defines it.
    # BY MIDLINE. `y` here is each word's TOP edge, and a number and the word
    # beside it share a baseline but not a top: 0.87pt apart on the committed
    # font, 2.70pt under Helvetica, and 2.26pt at annotation fontsize 26 --
    # each of which dropped `synergy` from the run and failed a correct
    # figure. Bottom edges are no better (0.73pt apart here). The vertical
    # centre differs by 0.1pt on the committed figure and stays close as the
    # font changes, because both words are set on one line. This file states
    # the same lesson forty lines below for the score matcher and had applied
    # it only there.
    # SCOPED TO PANEL (a). `hits` was bounded by panel (b)'s spine but the run
    # was built over the whole page, so the 20pt gap-walk absorbed panel (b)'s
    # top y-tick label once the panel grew: at 17 pairs the annotation read
    # `1.99× synergy RSL3+FSP1i` and a correct figure failed. Panel (b)'s pair
    # labels are drawn LEFT of its spine, so the spine cannot bound this --
    # panel (a)'s right edge can, and `_panels` is already called here.
    a_hi = _panels("fig25_bliss_synergy")[0][1]
    same_row = sorted(
        (x0, x1, t)
        for x0, y0, x1, y1, t in _word_bboxes("fig25_bliss_synergy")
        if abs((y0 + y1) / 2 - ay_mid) < 4 and x1 <= a_hi)
    at = next(i for i, (x0, _, t) in enumerate(same_row) if t == want)
    lo = hi = at
    while lo > 0 and same_row[lo][0] - same_row[lo - 1][1] < 20:
        lo -= 1
    while hi + 1 < len(same_row) and same_row[hi + 1][0] - same_row[hi][1] < 20:
        hi += 1
    line = [(x0, t) for x0, _, t in same_row[lo:hi + 1]]
    assert " ".join(t for _, t in line) == f"{want} synergy", (
        f"fig25's panel (a) synergy annotation reads "
        f"{' '.join(t for _, t in line)!r}, expected {want + ' synergy'!r}")
    # BY POSITION, for the reason fig24's needed it: a containment check does
    # not say WHICH panel carries the label, or that any panel does.
    # THE SUPTITLE, which states this figure's whole claim. Rewriting it to
    # say the dual-pathway depletion is ANTAGONISTIC -- the reverse of what
    # the synergy score and the additive reference line show -- passed every
    # check here. The same omission covered fig24's and fig26's headlines.
    assert "Dual-pathway (GPX4 + FSP1) depletion is synergistic" in text, (
        "fig25's suptitle no longer states that dual-pathway depletion is "
        "synergistic, which is the claim its synergy score and its additive "
        "reference line exist to support")

    cols = _vertical_labels(_boxed_row_words("fig25_bliss_synergy"))
    kill_axes = sorted(k[0] for k, (lbl, _, _) in cols.items()
                       if lbl == "Persister kill (%)")
    assert len(kill_axes) == 1, (
        "fig25 does not draw `Persister kill (%)` as a y-axis label exactly "
        f"once; vertical labels found: {sorted(v[0] for v in cols.values())}. "
        f"The four "
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
    # NEAREST, NOT WITHIN 3pt. `y` here is each word's TOP edge, and the pair
    # labels are y-ticks at `font.size` while the scores are drawn at
    # fontsize=8 -- vertically centred on the same row, their tops diverge by
    # about 0.73pt per point of size difference. A fixed 3pt bound therefore
    # broke on `font.size: 15`, on `tick_params(labelsize=16)`, and on a
    # larger score annotation, each a correct figure. In every one of those
    # the right partner was still nearest by 7x or more, which is the property
    # that actually holds -- the same nearest-neighbour lesson fig26's unit
    # check learned, left un-applied here.
    # RANKED WITHIN PANEL (b). The comment above says panel (a)'s synergy
    # annotation is excluded "by having no partner rather than by a
    # coordinate" -- true at three pairs, false as soon as panel (b) gains
    # rows and its top bar climbs toward that annotation. Measured, correct
    # figures failed at 12 pairs with font.size 12, at 9 pairs with 15, at 5
    # with 16, and at 25 pairs even at the committed size: the runner-up was
    # panel (a)'s annotation, 174pt left of panel (b)'s spine. x separates
    # them unambiguously and the spine is already computed above.
    b_spine = _panels("fig25_bliss_synergy")[-1][0]
    scoped = [v for v in vals if v[0] >= b_spine]
    assert scoped, (
        f"fig25 draws no synergy scores right of panel (b)'s axis at "
        f"x={b_spine:.0f}")
    by_row = {}
    for x, y, t in names:
        ranked = sorted(scoped, key=lambda v: abs(v[1] - y))
        assert ranked, f"fig25 draws no scores at all beside {t}"
        best = ranked[0]
        runner = ranked[1] if len(ranked) > 1 else None
        assert runner is None or abs(runner[1] - y) > 3 * abs(best[1] - y) + 3, (
            f"{t}'s nearest score is {abs(best[1] - y):.1f}pt away and the "
            f"next is {abs(runner[1] - y):.1f}pt; the rows are too close "
            "together to attribute a score to a pair")
        by_row[t] = best[2]
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
    # SPLIT ON THE GAP, not at a fixed x. Panel (b)'s caption shares a row
    # with panel (a)'s second-line bar labels -- 0.35pt apart at the committed
    # label size, 0.006pt at one point larger, so `round(y, 1)` alone stops
    # separating them. An absolute cut near panel (b)'s spine fixed that and
    # broke something worse: the caption is CENTRED, so once it is wider than
    # its panel its leading word falls left of the cut and is dropped --
    # measured, every `axes.labelsize` from 16 up and every figure width of 8
    # or less failed on a correct figure, where the fixed cut had failed at
    # exactly one label size. A knife edge became an open-ended one.
    #
    # What is real is the space between the two groups: measured edge to edge,
    # `combination` ends at 325.4 and `Bliss` starts at 399.1, a gap of 73.8pt
    # at the committed size and 63.3pt at label size 13, against word spacing
    # of a few points inside either group. So the row is split at its widest
    # gap when that gap dwarfs the typical one, and the right-hand group is
    # panel (b)'s caption. No width ceiling.
    #
    # (An earlier version of this comment said 76pt and called the approach
    # nearest-neighbour. Neither was true: the gap is 73.8pt, and
    # nearest-neighbour is what fig26's unit check does -- this is a split,
    # and the difference is exactly why the fixed cut it replaced broke.)
    row_words = _row_containing(_boxed_rows("fig25_bliss_synergy"),
                                caption[0][1])
    gaps = [(row_words[i + 1][0] - row_words[i][1], i)      # EDGE to edge
            for i in range(len(row_words) - 1)]
    cap_row = [t for _, _, t in row_words]
    if gaps:
        widest, at = max(gaps)
        typical = sorted(g for g, _ in gaps)[len(gaps) // 2]
        if widest > max(20.0, 3 * typical):
            cap_row = [t for _, _, t in row_words[at + 1:]]
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
    # THE LABEL MUST NAME THE POSITION THE LINE IS AT. `"additive" in text`
    # is a presence check: relabelling the entry `additive (2.0×)` while the
    # line stays at 1.0 leaves the legend contradicting the line every bar on
    # the panel is read against, and passed. The label is derived from the
    # same number the line is checked at.
    assert f"additive ({1.0:.1f}×)" in text.replace("$\\times$", "×"), (
        "fig25's panel (b) threshold legend does not read `additive (1.0×)`. "
        f"The line itself is drawn at {drawn_at:.2f} on the panel's own "
        "scale, so a legend naming any other value contradicts it")

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
    # BY MIDLINE, not top edge -- the third site with this assumption, and
    # the file said so while fixing only the other two. Under cmr10, the LaTeX
    # default, the arrow's top sits 3.86pt from the percentage beside it
    # (bound was 3) while their midlines are 1.08pt apart, and a correct
    # annotation reading `3.7% -> 0.1% (~37x collapse)` was rejected.
    boxes = {(round(x0, 2), t): (y0 + y1) / 2
             for x0, y0, x1, y1, t in _word_bboxes("fig24_hypoxia_killcurve")}
    arrow = [(x, boxes[(round(x, 2), t)]) for x, y, t in panel_a
             if t.strip() in {"→", "\u2192"}]
    assert arrow, "fig24's collapse annotation draws no arrow"
    arrow_x, arrow_mid = arrow[0]
    on_row = sorted((x, t) for x, y, t in panel_a
                    if re.fullmatch(r"\d+\.\d%", t)
                    and abs(boxes[(round(x, 2), t)] - arrow_mid) <= 3)
    left = [w for w in on_row if w[0] < arrow_x]
    right = [w for w in on_row if w[0] > arrow_x]
    assert left and right, (
        "fig24's collapse annotation does not read `N.N% -> N.N%`; the two "
        f"percentages flanking its arrow are what identify it. row={on_row}")
    bars = sorted(((x, t) for x, y, t in panel_a
                   if re.fullmatch(r"\d+\.\d%", t)
                   and not (abs(boxes[(round(x, 2), t)] - arrow_mid) <= 3
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
    # CLUSTERED, NOT ROUNDED -- this file's own lesson, still unapplied here.
    # Under Courier New the two group labels sit 0.17pt apart and straddle a
    # rounding boundary, so `round(y, 1)` put them on separate rows and a
    # clean, fully legible fig24 was rejected.
    named = [(x, y, t) for x, y, t in panel_a if t in ("RSL3", "SDT")]
    tick_rows = [[(x, t) for x, _, t in row]
                 for row in _cluster_rows(named, key=lambda w: w[1])]
    tick_rows = [r for r in tick_rows if {t for _, t in r} == {"RSL3", "SDT"}]
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
            # LEFT OF THE TEXT, AND CLOSE TO IT. There was no x term at all,
            # so `sorted` picked the nearest candidate ANYWHERE on the row --
            # and once the bars are narrower than 40pt the SDT normoxic BAR,
            # 87pt to the right, sat 3.28pt from `Hypoxic` and beat that
            # entry's own swatch at 4.20pt. Both entries then resolved to the
            # same colour and a correct figure was told it drew one colour
            # twice. A legend swatch is immediately left of its label; the
            # real one is 22.4pt left, the impostor 87pt right.
            near = sorted((abs(r[2] - y), r) for r in rects
                          if abs(r[2] - y) < 8 and r[1] - r[0] < 52
                          and 0 < x - r[1] < 60)
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

    # THE SUPTITLE. It states the figure's whole claim, and nothing read it:
    # rewriting it to say hypoxia SPARES pharmacologic ferroptosis and
    # collapses exogenous ROS -- the exact reverse of what the bars show --
    # passed every check here. This file already gates comparable non-numeric
    # claims (the score's definition, the additive line's label, the axis
    # unit), so leaving the headline unread was an omission rather than a
    # scope boundary.
    assert ("Hypoxia collapses pharmacologic ferroptosis but not exogenous ROS"
            in text), (
        "fig24's suptitle no longer states that hypoxia collapses "
        "pharmacologic ferroptosis but not exogenous ROS, which is what its "
        "bars and its panel (b) references are drawn to show")

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
    columns = _vertical_labels(_boxed_row_words("fig24_hypoxia_killcurve"))
    kill_cols = sorted(k for k, (lbl, _, _) in columns.items()
                       if lbl == "Overall tumor kill (%)")
    kill_axes = [k[0] for k in kill_cols]
    assert len(kill_axes) == 2, (
        f"fig24 draws `Overall tumor kill (%)` as a y-axis label at "
        f"{len(kill_axes)} positions, expected one per panel. Columns found: "
        f"{sorted(v[0] for v in columns.values())}. The values compared above "
        f"are drawn as "
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
    (a_x0, _a_x1), (b_x0, b_x1) = px[0], px[-1]
    assert kill_axes[0] < a_x0, (
        f"fig24's leftmost `Overall tumor kill (%)` label is at "
        f"x={kill_axes[0]:.0f}, inside panel (a) rather than left of its axis "
        f"at x={a_x0:.0f}. Panel (a)'s values are drawn without the axis "
        "label that says what they are")
    # EACH LABEL SITS THE SAME DISTANCE LEFT OF ITS OWN AXIS. That offset is
    # 48.2pt for both panels on the committed figure -- identical, because it
    # is set by the tick labels and the font, not by where the panel is. So
    # the invariant is that the two offsets AGREE, which needs no slack and
    # does not care about the figure's width.
    #
    # The bound this replaces was `a_x1 - 5% of panel width < x < b_x0`, and
    # its slack scaled the wrong way: the label's intrusion into panel (a) is
    # font-derived and absolute, so as the figure narrows the intrusion grows
    # while the slack shrinks. Measured, a correct fig24 at width 7 intruded
    # 12.6pt against 8.9pt of slack and failed.
    off_a = a_x0 - kill_axes[0]
    off_b = b_x0 - kill_axes[1]
    assert off_a > 0 and off_b > 0, (
        f"fig24 draws a `Overall tumor kill (%)` label at x={kill_axes} that "
        f"is not left of its own axis (panels start at {a_x0:.0f} and "
        f"{b_x0:.0f}); a y-axis label sits outside the data it labels")
    # AND EACH LABEL MUST SIT BESIDE ITS PANEL VERTICALLY. The offsets above
    # compare x only, so a rotated label floating ABOVE panel (a) -- with no
    # label beside its bars at all -- satisfied them. The spines give the
    # panel's y range as well as its x range.
    y_lo, y_hi = _panel_y("fig24_hypoxia_killcurve")
    for key in kill_cols:
        _, col_lo, col_hi = columns[key]
        x = key[0]
        # OVERHANG IS NORMAL. A centred y-label longer than its axis extends
        # past both ends -- at `axes.labelsize` 24 fig24's spans 45 to 221
        # against an axis of 53 to 302 -- and that is ordinary matplotlib
        # output, not a mislabelled panel. What matters is that the label is
        # CENTRED on its axis, so the test is on its midpoint.
        col_mid = (col_lo + col_hi) / 2
        axis_mid = (y_lo + y_hi) / 2
        assert abs(col_mid - axis_mid) <= (y_hi - y_lo) * 0.25, (
            f"fig24's `Overall tumor kill (%)` label at x={x:.0f} is "
            f"centred at y={col_mid:.0f} while its panel's axis is centred at "
            f"{axis_mid:.0f} (spanning {y_lo:.0f} to {y_hi:.0f}); a y-axis "
            "label sits beside the data it labels")

    tol = max(4.0, 0.15 * max(off_a, off_b))
    assert abs(off_a - off_b) <= tol, (
        f"fig24's two `Overall tumor kill (%)` labels sit {off_a:.0f}pt and "
        f"{off_b:.0f}pt left of their panels' axes. A y-axis label is offset "
        "by its own tick labels and font, so the two should agree -- one of "
        "these is not labelling the panel beside it")

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
        # Same correction as fig26's, for the same reason: the dashed sample
        # is centred on its label, so a band measured from the label's top
        # loses it once the legend font grows past about 13pt.
        y_mid = next((yy0 + yy1) / 2
                     for xx0, yy0, xx1, yy1, tt in _word_bboxes(
                         "fig24_hypoxia_killcurve")
                     if tt == t and abs(xx0 - x) < 0.01 and abs(yy0 - y) < 0.01)
        sample = sorted((d for d in dashed
                         if abs(d[0] - y_mid) < 12 and d[2] - d[1] < 40),
                        key=lambda d: abs(d[0] - y_mid))
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
    # READ IN DATA UNITS, off panel (b)'s own y-axis ticks. Comparing only
    # which line is higher left the VALUES unread: drawing SDT's reference at
    # 0.6x its real rate -- 55% while its own bar says 91.9% -- kept the
    # ordering and passed every semantic check. And the ordering itself was
    # hardcoded, so a run where RSL3's uniform kill genuinely exceeded SDT's
    # would be told its lines were exchanged when they were correct. The tick
    # labels give the scale, the same way fig25's additive line is checked.
    # SCALE FROM THE TICK MARKS, NOT THE TICK LABELS. A label's bbox centre
    # carries that font's digit-glyph asymmetry as a constant offset, so the
    # mapping was biased by the typeface: 0.006pp on the committed serif but
    # +0.243 on monospace, +0.386 on STIXGeneral and +1.143 on Hoefler Text.
    # Widening the bound to absorb that was the wrong repair -- it let a real
    # 0.45pp error through and still rejected correct figures in five journal
    # serifs. The tick MARKS are byte-identical across every font tried
    # (52.560, 101.003, 149.445, 197.888, 246.331, 294.774), because they are
    # geometry rather than glyphs, so the bias disappears and the bound can be
    # tighter than the one it replaces.
    marks = sorted({round(a.y, 3)
                    for x0, y0, x1, y1, col, a, b in _axis_ticks(
                        "fig24_hypoxia_killcurve")
                    if abs(min(a.x, b.x) - b_x0) < 8})
    labels = sorted(
        ((y0 + y1) / 2, float(t))
        for x0, y0, x1, y1, t in _word_bboxes("fig24_hypoxia_killcurve")
        if re.fullmatch(r"\d+", t) and b_x0 - 40 < x1 < b_x0 + 2)
    assert len(marks) >= 2 and len(marks) == len(labels), (
        f"fig24's panel (b) shows {len(marks)} y tick marks against "
        f"{len(labels)} numeric labels; its scale cannot be recovered")
    # Marks and labels are both in y order, and the tick pitch (48.4pt) dwarfs
    # any glyph offset, so pairing by order is unambiguous.
    ty0, tv0 = marks[0], labels[0][1]
    ty1, tv1 = marks[-1], labels[-1][1]
    scale = (ty1 - ty0) / (tv1 - tv0)
    for name in ("SDT", "RSL3"):
        drawn = (refs[ref_colour[name]][0] - ty0) / scale + tv0
        want = kills[name][0]
        # 0.1pp, AGAINST A MAPPING THAT NO LONGER HAS A FONT BIAS. Reading
        # the tick MARKS rather than the tick labels drops the worst residual
        # to 0.0011pp across fourteen style variations, so this bound rejects
        # none of them while killing a +0.15pp error. (The block that stood
        # here argued for 0.5pp and described the label-centre mapping the
        # same commit had already replaced -- it also quoted a 0.013pp
        # residual, which was neither mapping's.)
        assert abs(drawn - want) <= 0.1, (
            f"fig24's panel (b) draws {name}'s normoxic reference at "
            f"{drawn:.1f}% on its own y scale and the data gives "
            f"{want:.1f}%. That line is what the depth curves are compared "
            "against, so it has to be the value it claims")

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
    # NO STREAM-ORDER RUN CHECK HERE. There was one, asserting the nine
    # labels appear in order in the extracted text, and it added no coverage
    # the positional check below does not already have -- while failing FIRST,
    # with a message about a panel losing its labels, in cases where PyMuPDF
    # merely merges two crowded tick words into one (a narrow figure, or more
    # timepoints). Stream order is not layout, which is why the positional
    # check exists; keeping a weaker duplicate in front of it only made the
    # diagnosis worse.

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
    # CLUSTERED, NOT ROUNDED. Found by searching for the shape rather than
    # waiting for it to be reported: this accumulator still keyed rows on
    # `round(y, 1)` after that was corrected at three other sites in this
    # file. A row that straddles a rounding boundary splits, which both
    # changes what `[:2]` below selects and lets a word fall outside the
    # collision scan's row filter.
    rows26 = _cluster_rows(numeric, key=lambda w: w[1])
    rows26.sort(key=lambda r: -min(y for _, y, _ in r))   # nearest the caption
    band = []
    for row in rows26:
        if len(band) >= 2 * len(expected):
            break
        band.extend((x, t) for x, _, t in row)
    # DIAGNOSE A COLLISION RATHER THAN JUST COUNTING. When two tick labels
    # overlap on the page PyMuPDF returns them as ONE word -- at font.size 15
    # panel (a) draws `0.250.5` -- and the figure really is defective, but a
    # bare count reports it as missing labels and sends the reader looking for
    # the wrong thing.
    # SCANNED ON THE RAW WORDS, not on the band. A collided pair reads
    # `0.250.5` -- two decimal points -- so the numeric pattern that builds
    # the band rejects it, panel (a) contributes one label too few, and the
    # accumulator reaches its count by pulling in a y-axis tick row. The
    # symptom then looks like the wrong labels rather than like two labels on
    # top of each other.
    # ADJACENT PAIRS ONLY, AND ONLY ON THE TICK ROWS. Testing `t in joined`
    # -- a substring of every day label concatenated -- let a y-axis tick
    # impersonate a collision: with a 0.75-day timepoint added, `joined`
    # contains "0.50", which is panel (b)'s GPX4 axis tick 36.5pt above the
    # caption at a completely different x. A correct figure failed being told
    # its x axis was unreadable. Only two labels drawn side by side can run
    # together, so only adjacent concatenations count, and only on the rows
    # the band above already identified as tick rows.
    tick_ys = {round(y, 2) for row in rows26[:2] for _, y, _ in row}
    # RUNS OF ANY LENGTH, not just pairs. Three or more labels can merge --
    # `00.250.5` at one figure size, `00.250.50.751` with sixteen timepoints --
    # and a pairs-only set stayed silent, letting the failure revert to the
    # count message that drags in panel (b)'s GPX4 tick. The verdict was still
    # a correct FAIL; the diagnosis regressed.
    adjacent = {"".join(expected[i:j])
                for i in range(len(expected))
                for j in range(i + 2, len(expected) + 1)}
    merged = sorted({t for x, y, t in wds
                     if any(abs(y - ty) < 1.0 for ty in tick_ys)
                     and t in adjacent})
    assert not merged, (
        f"fig26 draws {merged} as single words: adjacent tick labels have "
        "collided and been run together, so the axis is unreadable at this "
        "size. The expected labels are " + ", ".join(expected))
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
    stem26 = "fig26_vulnerability_window"
    p26 = _panels(stem26)
    a_lo, a_hi = p26[0]
    doc_lines = _hlines_any("fig26_vulnerability_window")
    legend_col = {}
    for x, y, t in wds:
        if t not in ("SDT", "RSL3") or not (a_lo <= x <= a_hi):
            continue
        # NEAREST, not lowest-x. `_hlines_any` returns sorted by x, so
        # `sample[0]` took whichever candidate started furthest left rather
        # than the one belonging to this entry. fig24's equivalent already
        # ranks by distance.
        # BY MIDLINE, AND NEAREST. The sample sits at the label's vertical
        # CENTRE, so a one-sided 10pt band from the label's TOP edge grows
        # out of range as the legend font grows: at fontsize 16 the offset is
        # 12.16pt and a correct figure was told its legend had no samples.
        # This is the third site with the top-edge assumption; the previous
        # round fixed one and said so.
        y_mid = next((yy0 + yy1) / 2
                     for xx0, yy0, xx1, yy1, tt in _word_bboxes(stem26)
                     if tt == t and abs(xx0 - x) < 0.01 and abs(yy0 - y) < 0.01)
        sample = sorted((seg for seg in doc_lines
                         if abs(seg[1] - y_mid) < 12 and abs(seg[0] - x) < 60),
                        key=lambda seg: (abs(seg[1] - y_mid), abs(seg[0] - x)))
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

    # THE WINDOW SHADING, in tick units. `axvspan(-0.3, win_end)` is what shows
    # WHERE the window is, and nothing read it: shading the whole 28 days, or
    # only days 14-28, passed. `win_end` is the last timepoint at or under day
    # 3, derived from the fixture, and the tick marks give the x scale.
    xmarks = sorted({round(a.x, 2)
                     for x0, y0, x1, y1, col, a, b in _axis_ticks_vertical(stem26)
                     if a.x <= a_hi})
    assert len(xmarks) == len(expected), (
        f"fig26 panel (a) shows {len(xmarks)} x tick marks against "
        f"{len(expected)} timepoints")
    win_end = max(i for i, d in enumerate(days) if d <= 3.0)
    shade = [r for r in _filled_rects(stem26)
             if r[0] >= a_lo - 1 and r[1] <= a_hi + 1
             and (r[1] - r[0]) > 20 and (r[3] - r[2]) > 100
             and r[4] not in ((1.0, 1.0, 1.0),)]
    assert len(shade) == 1, (
        f"fig26's panel (a) draws {len(shade)} shaded spans, expected one "
        "marking the RSL3 window")
    # BOTH EDGES. Reading only the right one let the window be made to OPEN
    # late: `axvspan(win_end - 1, win_end)` shades days 2 to 3 instead of 0 to
    # 3, so the figure claims the ferroptosis-sensitive window does not open
    # until day 2 -- the inverse of the panel's argument -- and passed. The
    # generator opens it at `-0.3`, i.e. just left of the first tick, so the
    # left edge belongs to day 0.
    left_at = min(range(len(xmarks)), key=lambda i: abs(xmarks[i] - shade[0][0]))
    assert left_at == 0, (
        f"fig26 opens the RSL3 window at the tick for day {expected[left_at]}, "
        "and the data has it open from the first timepoint. The shaded span is "
        "what shows where the window is")
    right_at = min(range(len(xmarks)), key=lambda i: abs(xmarks[i] - shade[0][1]))
    assert right_at == win_end, (
        f"fig26 shades the RSL3 window out to the tick for day "
        f"{expected[right_at]}, and the data closes it at day "
        f"{expected[win_end]}. The shaded span is what shows where the window "
        "is; the annotation beside it is a fixed string")

    assert ("The ferroptosis-sensitive window: days for RSL3, weeks for SDT"
            in text), (
        "fig26's suptitle no longer says `days for RSL3, weeks for SDT`. The "
        "curve check below compares against exactly that claim, and quoted it "
        "in its own failure message while never reading it")

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
    "fig24_hypoxia_killcurve": "bc6fb1c80d9f3a30",
    "fig25_bliss_synergy": "5508eb163d0d74e7",
    "fig26_vulnerability_window": "41d03195ab16c649",
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
        page_rect = page.rect
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
        # (kept as `words` because that is the JSON key the hash uses)
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
                # `isolated` and `knockout` belong HERE. They were in the
                # path paint tuple, where PyMuPDF never sets them -- the group
                # branch returns before that tuple is built -- so replacing
                # both with names that do not exist left all three hashes
                # unchanged. Inert keys read as coverage.
                items.append(("group", d.get("level"),
                              None if r is None else
                              tuple(round(v, 2) for v in (r.x0, r.y0, r.x1, r.y1)),
                              d.get("blendmode"), d.get("opacity"),
                              d.get("isolated"), d.get("knockout")))
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
                                  "fill_opacity", "even_odd", "closePath",
                                  "lineCap", "lineJoin")))
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
    # RASTER IMAGES TOO. `get_drawings` returns no image XObjects and
    # `get_text` cannot see them, so an entire layer was unread: painting a
    # solid black raster over 18.9% of fig25 -- covering the panel -- left the
    # hash byte-identical and all seven tests green. It is the one layer that
    # can hide a whole figure, and it is reachable in this very file, since
    # fig17 already draws its heatmaps with `imshow`. Both the xref metadata
    # and where each image LANDS are recorded, because the same image moved
    # is a different page.
    doc = _reader().open(FIG_DIR / f"{stem}.pdf")
    try:
        page = doc[0]
        images = []
        for info in page.get_images(full=True):
            rects = sorted(tuple(round(v, 2) for v in r)
                           for r in page.get_image_rects(info[0]))
            images.append((info[1:], rects))
    finally:
        doc.close()

    # THE PAGE BOX TOO. Every coordinate above is measured from the top-left,
    # so extending the MediaBox leaves all of them unchanged and the hash
    # unmoved while the rendered page is a different shape. Not reachable from
    # a generator edit today -- matplotlib derives the box from the content --
    # but it is one line, and "a layer it did not read" is the recurring
    # shape of this backstop's history.
    box = tuple(round(v, 2) for v in (page_rect.x0, page_rect.y0,
                                      page_rect.x1, page_rect.y1))
    blob = json.dumps({"words": words, "items": items, "images": images,
                       "box": box}, sort_keys=True)
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
