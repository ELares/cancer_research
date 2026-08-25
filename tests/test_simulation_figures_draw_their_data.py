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
- **fig26 panel (b) is ungated except for its axis captions and titles.** Its
  GPX4 right-hand axis draws data-derived tick labels no assertion reads, and
  -- said explicitly here because the parallel fig24 bullet says it and this
  one did not -- its CURVES are unread too, exactly as fig24 panel (b)'s are.
  Only the fingerprint would notice either moving.

  fig24 panel (b) was in this list and is partly out of it now: its two dashed
  reference lines are bound to its own legend below, its lambda tick labels and
  both its axis captions ARE read (an earlier version of this sentence listed
  those as unread, in the same commit that gated them). Its two kill curves
  remain unread except by the fingerprint.
- **The collapse arrow's glyph is pinned, so `mathtext.fontset: cm` fails
  here** -- and so do `plt.style.use("bmh")` and `plt.style.use("classic")`,
  both of which set that fontset. cm
  embeds cmsy10 with no Unicode map, so PyMuPDF reads `$\to$` as `!`. This
  was briefly not a limit: the pattern matched "one non-space character" in
  the arrow's place so that any font would do, and that silently removed the
  ONLY check on which way the collapse runs -- `3.7% <- 0.1%` reads as
  hypoxia RAISING the kill rate from 0.1% to 3.7%, the #790 inversion, and it
  passed. Both values stay in order when only the arrow flips, so nothing
  else can see it. A loud failure on an unusual font is the better half of
  that trade, and it is stated rather than silent. `_digits` and `_mult`
  handle the font-dependent glyphs that carry no direction (`$O_2$` under
  `stixsans`, `$\times$`), which is why they are normalised and this is not.

- **A title that COLLIDES with the next panel's is reported, not tolerated.**
  At `axes.titlesize: 16` fig26 draws `(b)` on top of `stays` and `open` on
  top of `Why:`; the figure is unreadable there and the check says so. The
  remedy is to wrap the title, and that now works -- wrapping both fig26
  titles at that size passes.

- **A panel needs at least one horizontal spine DRAWN.** `_panels` locates
  panels from long horizontal axis rules, so turning off BOTH
  `axes.spines.top` and `axes.spines.bottom` removes what it looks for -- and
  so does `axes.linewidth: 0`, which leaves the spines enabled and draws
  nothing. That is how `seaborn-v0_8`, `-dark` and `-darkgrid` reach this
  limit; naming only the `spines.*` route would have read as though a style
  that sets neither were safe. Any single spine removed is
  fine, including `bottom`, whose absence used to put the axis 6.8pt out and
  silently emptied fig26's x scale.

  The failure is NOT uniform, and an earlier version of this bullet said it
  was ("every reader here aborts"). Measured with both off: fig25 aborts with
  "no horizontal axis rules found", but fig24 returns `[(56.8, 176.4)]` and
  fig26 returns `[(330.0, 446.4), (624.3, 781.2)]` -- LEGEND FRAMES, because a
  legend frame is also a long horizontal rule and with the real ones gone
  nothing outranks it. The suite still fails overall, so nothing ships, but
  the mechanism is `_panels` mistaking a legend for a panel rather than
  declining to answer.

- **The figure must draw visible tick marks.** THREE scale readings -- fig24
  panel (b)'s reference lines, fig26's window shading, and fig25's additive
  threshold line -- locate themselves from the tick MARKS. fig25's was the
  last on labels and joined them because label bias is not only the font's
  glyph-centre offset (up to 1.14pp under Hoefler Text): `xtick.alignment:
  left` shifts every label relative to its tick, and the additive line then
  read 0.92 instead of 1.0, accusing a correct figure of moving the line that
  every score on the panel is compared against. A style that sets
  `xtick.major.size: 0` or `ytick.major.size: 0` therefore fails loudly with
  "shows 0 tick marks" or
  "its scale cannot be recovered" rather than silently. This project's own
  rcParams set neither, and the failure names its cause, so it is a stated
  requirement rather than a defect to chase: a label-centre fallback was tried
  and reintroduced the bias it was built to remove.

  Which styles those are, measured on the pinned matplotlib 3.11.1 rather
  than characterised: six of the 28 available -- `fivethirtyeight`,
  `seaborn-v0_8`, `-dark`, `-darkgrid`, `-white` and `-whitegrid`. An earlier
  version of this bullet said "every `seaborn-v0_8` grid preset", which is
  wrong in both directions: `-white` and `-dark` are not grid presets and do
  set it, while `-ticks` (6.0) and `-talk` (3.5) are seaborn and do not, and
  `fivethirtyeight` is not seaborn at all.

  This bullet was WRONG about fig24 until the fallback was deleted. fig24 kept
  one, widening its bound 0.1pp -> 1.5pp when no marks were drawn, so under
  `ytick.major.size: 0` that reading did not fail at all -- it silently ran at
  a fifteenth of its precision and passed a reference line 1.1pp off its
  value. Only fig26 behaved as described here.

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
import unicodedata
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
    # FULL EXTENTS. This took `(x0, x1, y0, text)` and reported min/max over
    # TOP edges, so any midpoint computed from it was biased upward by half a
    # line height -- at `axes.labelsize` 24 the y-label read 44.3pt
    # off-centre when its true offset was 0.0, consuming 71% of the
    # centredness budget with an artifact. It takes full boxes now.
    cols = {}
    for x0, y0, x1, y1, t in words:
        cols.setdefault((round(x0, 1), round(x1, 1)), []).append((y0, y1, t))
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
    return {key: (" ".join(t for _, _, t in sorted(ws, reverse=True)),
                  min(y0 for y0, _, _ in ws), max(y1 for _, y1, _ in ws))
            for key, ws in cols.items() if len(ws) > 1}


def _panels(stem):
    """Each panel's x range, from the axis rules: `[(x0, x1), ...]` left to right.

    NESTING IS THE DISCRIMINATOR, not colour, not pairing, not width. Three
    earlier versions each keyed on something that is true of the committed
    figure and not of a figure in general:

      * exactly black -- `axes.edgecolor="0.3"` broke it;
      * greyscale plus a brightness cut -- matplotlib's own legend frame is
        greyscale 0.8, so the cut was doing the real work and any edgecolor
        at or above 0.76 aborted every reader here;
      * paired rules at least 90% of the widest -- `axes.spines.top: False`
        leaves ONE rule per panel and aborted all three figures, and
        `width_ratios=[1.2, 1]` dropped the narrower panel outright.

    What is actually true is structural, in two steps. A legend frame drawn
    INSIDE its panel is contained by it and a panel is contained by nothing,
    so discard every x range another range contains. And a legend frame drawn
    OUTSIDE its panel -- `bbox_to_anchor=(1.02, 0.5)` -- is contained by
    nothing either, so nesting alone reported it as a third panel; what
    separates it there is that a row of subplots SHARES one vertical extent
    while a legend has its own, so the largest such group is the panels, and
    an equal-sized group is settled on total width.
    No constant, and it survives a despined, recoloured or unequally sized
    figure with its legend inside or out.
    """
    pymupdf = _reader()
    doc = pymupdf.open(FIG_DIR / f"{stem}.pdf")
    try:
        rules = set()
        for d in doc[0].get_drawings():
            if d.get("color") is None:
                continue
            # A DASHED FULL-WIDTH RULE IS A REFERENCE LINE, NOT A SPINE.
            # `axhline` spans the whole axes, so fig24 panel (b)'s two dashed
            # normoxic references carry that panel's exact x range. While
            # only the x range was read that was harmless -- they duplicated a
            # span the spine already gave. Reading their y is not: they put
            # panel (b) in a different vertical-extent group from panel (a),
            # and under `axes.spines.top: False` the two panels stopped being
            # a group at all.
            dashes = (d.get("dashes") or "").strip()
            if dashes and dashes != "[] 0":
                continue
            for item in d["items"]:
                if item[0] != "l":
                    continue
                a, b = item[1], item[2]
                if abs(a.y - b.y) < 0.5 and abs(b.x - a.x) > 100:
                    rules.add((round(min(a.x, b.x), 1),
                               round(max(a.x, b.x), 1), round(a.y, 1)))
    finally:
        doc.close()
    assert rules, f"{stem}: no horizontal axis rules found, so panels cannot be located"
    spans = {(x0, x1) for x0, x1, _ in rules}
    outer = [k for k in spans
             if not any(o != k and o[0] <= k[0] and k[1] <= o[1] for o in spans)]
    # NESTING IS NOT ENOUGH ONCE THE LEGEND LEAVES THE AXES. Placed with
    # `bbox_to_anchor=(1.02, 0.5)` a legend frame sits BESIDE its panel rather
    # than inside it, so nothing contains it and fig26 reported THREE panels
    # -- which then mis-assigned the axis captions and reported a correct
    # figure as drawing no unit.
    #
    # What still separates them is exact and needs no threshold: subplots in a
    # row share one vertical extent, and a legend has its own. So group the
    # candidates by the (top, bottom) their own rules span, and take the
    # LARGEST GROUP -- a row of panels, against a legend that is one box.
    #
    # Size before height, because height alone is not enough: under
    # `axes.spines.top: False` a panel's rules collapse to a single y and its
    # extent is ZERO while the legend keeps its whole box, so "the tallest
    # group" chose the legend and reported one panel where there are two.
    #
    # AND THE TIE IS BROKEN ON TOTAL WIDTH, NOT HEIGHT. An earlier version
    # said height, and justified it with "which is what a single-panel figure
    # with an outside legend gives" -- a claim about when ties arise that was
    # simply untrue. Two panels with two outside legends tie on size just as
    # readily, and combined with the despining above (where the panels' height
    # is zero and the legends' is not) `_panels` returned the two LEGEND
    # FRAMES. Width does not depend on the spines, which is why it is the
    # better tie-break.
    #
    # IT IS A TIE-BREAK, NOT AN INVARIANT, and the sentence that stood here
    # ("panels are the wide things on a figure and a legend is a box beside
    # them") overstates it. A legend can be made wider than its panel:
    # `legend.fontsize: 20` with `handlelength: 8` and two entries outside the
    # axes measures 418pt against a 279pt panel, and `_panels` then returns
    # the legend. On fig26 the same legend CONTAINS panel (b), so the nesting
    # filter deletes the real panel before this code runs at all. Both fail
    # loudly rather than silently, and both are reachable at the parent too,
    # so this is a known bound on the helper and not a claim that it cannot
    # happen.
    extent = {}
    for k in outer:
        ys = [y for x0, x1, y in rules if (x0, x1) == k]
        extent[k] = (round(min(ys), 1), round(max(ys), 1))
    groups = {}
    for k, ext in extent.items():
        groups.setdefault(ext, []).append(k)
    best = max(groups.items(),
               key=lambda kv: (len(kv[1]), sum(hi - lo for lo, hi in kv[1])))
    return sorted(best[1])


def _axes_box(stem):
    """The plotting area's (top, bottom), from the CLIP matplotlib applies.

    THE CLIP IS THE ANSWER, and three rounds of hunting for it in the spines
    were the wrong question. matplotlib clips a panel's artists to its axes
    rectangle, so `get_drawings(extended=True)` carries a `clip` entry whose
    `scissor` IS the plotting area -- measured identical to `_panel_y` on all
    three committed figures, (52.56, 302.04) in every case.

    What that removes, all at once: it does not care whether a spine is drawn,
    so no despining breaks it; it does not care how far a confidence band
    reaches, because the band is the thing being clipped; and it does not care
    about `axes.xmargin`, `ytick.alignment` or any other placement rcParam,
    because a scissor is geometry the renderer emitted rather than a property
    inferred from what happened to be drawn.

    Returns None when no clip spans a panel, and the spine and tick paths stay
    behind it for that case.
    """
    pymupdf = _reader()
    panel_x = _panels(stem)
    doc = pymupdf.open(FIG_DIR / f"{stem}.pdf")
    try:
        boxes = []
        for d in doc[0].get_drawings(extended=True):
            if d.get("type") != "clip":
                continue
            sc = d.get("scissor")
            if sc is None:
                continue
            if any(abs(sc.x0 - lo) < 1 and abs(sc.x1 - hi) < 1
                   for lo, hi in panel_x):
                boxes.append((sc.y0, sc.y1))
    finally:
        doc.close()
    if not boxes:
        return None
    return min(b[0] for b in boxes), max(b[1] for b in boxes)


def _long_vertical_ends(stem):
    """Both y ends of every long vertical rule STANDING ON A PANEL EDGE.

    A left or right spine runs the full height of the axes, so its ends ARE
    the plotting area's top and bottom. `_panel_y` and `_axes_top` both need
    that when a horizontal spine has been turned off, and each was reaching
    for it separately.

    THE LONGEST STROKE AT EACH PANEL EDGE, because "long and vertical" is not
    enough to be a spine and neither is "at a panel edge". A `fill_between`
    band is drawn with a stroke colour, so its edges survive every other
    filter here: widening fig26's confidence band to +/-40 puts EIGHT more
    values in this list, one of them at 385.96, below the axis at 302.04, and
    two ABOVE the axes entirely. Filtering to a panel edge removes the
    interior ones but not the outermost, because `axes.xmargin: 0` puts the
    first and last data points ON the spine -- a fact this file documents in
    two other places -- so the band's own outer edges land exactly there.

    "Longest at each panel edge" was the next attempt and it is FALSE in the
    configuration the paragraph above names. PyMuPDF returns PRE-CLIP
    geometry, so a band edge running from above the axes to below them is
    LONGER than the spine: at fig26's shared edge with `axes.xmargin: 0` and
    a band wider than the y range, the spine is 249.48 and two band edges are
    272.36 and 272.17, so the reduction picked a band and `_axes_top` returned
    -66.61 for a true 52.56.

    This helper is now a FALLBACK behind `_axes_box`, which reads the clip
    matplotlib actually applies and is not guessing at all. What is kept here
    is the part that does hold: a stroke enclosed by a clip spanning the whole
    panel is that panel's DATA, whatever its length, and the reduction is
    keyed on the panel EDGE rather than on the stroke's own x -- an earlier
    version keyed on `round(a.x, 1)`, so a stroke admitted by the 1pt
    tolerance but more than 0.05pt off got its own key and was never reduced
    against the spine at all.
    """
    pymupdf = _reader()
    panel_x = _panels(stem)
    edges = {x for panel in panel_x for x in panel}
    doc = pymupdf.open(FIG_DIR / f"{stem}.pdf")
    try:
        clipped_to_panel = False
        longest = {}
        for d in doc[0].get_drawings(extended=True):
            if d.get("type") == "clip":
                sc = d.get("scissor")
                clipped_to_panel = sc is not None and any(
                    abs(sc.x0 - lo) < 1 and abs(sc.x1 - hi) < 1
                    for lo, hi in panel_x)
                continue
            if d.get("type") == "group" or d.get("color") is None:
                continue
            dashes = (d.get("dashes") or "").strip()
            if dashes and dashes != "[] 0":
                continue
            if clipped_to_panel:      # this is a panel's data, not its frame
                continue
            for item in d["items"]:
                if item[0] != "l":
                    continue
                a, b = item[1], item[2]
                if abs(a.x - b.x) >= 0.5 or abs(b.y - a.y) <= 100:
                    continue
                near = [e for e in edges if abs(a.x - e) < 1]
                if not near:
                    continue
                key = round(near[0], 1)          # THE EDGE, not the stroke
                span = (min(a.y, b.y), max(a.y, b.y))
                if key not in longest or (span[1] - span[0]) > (
                        longest[key][1] - longest[key][0]):
                    longest[key] = span
        return [v for span in longest.values() for v in span]
    finally:
        doc.close()


def _panel_y(stem):
    """The vertical extent of the plotting area, from the axis spines."""
    pymupdf = _reader()
    doc = pymupdf.open(FIG_DIR / f"{stem}.pdf")
    try:
        rules = []
        for d in doc[0].get_drawings():
            col = d.get("color")
            # ANY COLOUR. Requiring black, then greyscale-and-dark, each
            # broke a shipped style: `axes.edgecolor` is freely settable and
            # matplotlib's own legend frame is greyscale 0.8. A tick mark is
            # identified by its SHAPE -- a short stroke perpendicular to an
            # axis -- not by how it is painted.
            if col is None:
                continue
            # A DASHED FULL-WIDTH RULE IS A REFERENCE LINE, NOT A SPINE --
            # the THIRD site to need this guard. `_panels` and `_axes_top`
            # got it and this one did not, so under `axes.spines.bottom:
            # False` fig24 panel (b)'s dashed normoxic references were the
            # only rules left at that x range and the axis "bottom" came back
            # as 285.91 instead of 302.04. That also kept the tick fallback
            # below from firing, since three rules looked like plenty.
            dashes = (d.get("dashes") or "").strip()
            if dashes and dashes != "[] 0":
                continue
            for item in d["items"]:
                if item[0] != "l":
                    continue
                a, b = item[1], item[2]
                if abs(a.y - b.y) < 0.5 and abs(b.x - a.x) > 100:
                    rules.append((round(min(a.x, b.x), 1),
                                  round(max(a.x, b.x), 1), round(a.y, 2)))
    finally:
        doc.close()
    # ONLY RULES THAT BELONG TO A PANEL. Taking every long horizontal rule
    # included a legend placed BELOW the axes -- `bbox_to_anchor=(0.5, -0.22)`
    # -- whose frame put the bottom at 375.3 instead of the axis at 302.0, so
    # no tick touched the "axis" and all nine were discarded. A panel's rules
    # span exactly that panel's x range; a legend's do not.
    panel_x = _panels(stem)
    ys = {y for x0, x1, y in rules
          if any(abs(x0 - lo) < 1 and abs(x1 - hi) < 1 for lo, hi in panel_x)}
    # ONE SPINE IS ENOUGH WHEN SOMETHING ELSE GIVES THE OTHER END.
    # `axes.spines.top: False` is the commonest publication restyle and it
    # removed a whole row, aborting fig24 on a correct figure.
    #
    # THE SIDE SPINES FIRST, THEN THE TICKS. A vertical spine runs the full
    # height of the axes, so it gives the missing end exactly; the y-tick
    # positions only bracket it, and they stop at the outermost TICK. Under
    # `axes.spines.bottom: False` on fig26 -- ylim(-3, 107), so the lowest
    # tick is 7pt above the axis -- the tick fallback returned 295.24 for an
    # axis at 302.04, and every x tick then missed it and the panel read as
    # having no scale at all. The ticks stay as the last resort, for a figure
    # despined on every side.
    if len(ys) < 2:
        box = _axes_box(stem)
        if box is not None:
            ys |= {round(v, 2) for v in box}
    if len(ys) < 2:
        ys |= {round(v, 2) for v in _long_vertical_ends(stem)}
    if len(ys) < 2:
        tick_ys = {round(a.y, 2) for _, _, _, _, _, a, _ in _axis_ticks(stem)}
        ys |= tick_ys
    assert len(ys) >= 2, (
        f"{stem}: fewer than two spine rows and no tick marks to close the "
        "vertical extent")
    return min(ys), max(ys)


def _axes_top(stem):
    """The y of the top of the plotting area.

    `_panel_y` answers this from horizontal rules and falls back to the TICK
    positions when a spine is missing -- which returns the topmost tick, well
    inside the axes, not the top of them. That is fine for its own callers and
    wrong as a title-block boundary: under `axes.spines.top: False` the band
    then swallowed the y-axis labels and the legend, and the panel title read
    as `(a) Kill collapse under hypoxia 100 Normoxic (uniform O2) ...`.

    Either spine gives the same answer, so both are read: the top horizontal
    rule spanning a panel, and the upper end of a long vertical rule.

    "Anything drawn INSIDE the axes has a larger y than their top, so taking
    the minimum is safe against data strokes" stood here and is FALSE, which
    the body of this function goes on to say: PyMuPDF reports pre-clip
    geometry, so a confidence band reaches above the axes with no spine
    missing at all. `_axes_box` is consulted first for exactly that reason.

    Returns None when the top genuinely cannot be established -- one lone
    horizontal rule and no vertical spine, which is what `axes.spines.left`,
    `.top` and `.right` all set False leaves. GUESSING THERE IS WORSE THAN
    DECLINING: the single remaining rule is the BOTTOM of the axes, and using
    it put every tick label and the whole legend inside the title block. The
    caller falls back to the one-row reading, which cannot see a wrapped
    title but also cannot invent one.
    """
    pymupdf = _reader()
    panel_x = _panels(stem)
    doc = pymupdf.open(FIG_DIR / f"{stem}.pdf")
    try:
        horiz = []
        for d in doc[0].get_drawings():
            if d.get("color") is None:
                continue
            # A DASHED FULL-WIDTH RULE IS A REFERENCE LINE, NOT A SPINE.
            # `axhline` spans the whole axes, so fig24 panel (b)'s two dashed
            # normoxic references match a panel's x range exactly. With the
            # top spine present that is harmless -- the real top is smaller --
            # but under `axes.spines.top: False` the upper reference BECAME
            # the top, and the title block then swallowed the tick labels and
            # the legend.
            dashes = (d.get("dashes") or "").strip()
            if dashes and dashes != "[] 0":
                continue
            for item in d["items"]:
                if item[0] != "l":
                    continue
                a, b = item[1], item[2]
                if abs(a.y - b.y) < 0.5 and abs(b.x - a.x) > 100:
                    lo, hi = min(a.x, b.x), max(a.x, b.x)
                    if any(abs(lo - x0) < 1 and abs(hi - x1) < 1
                           for x0, x1 in panel_x):
                        horiz.append(a.y)
    finally:
        doc.close()
    # THROUGH THE SHARED HELPER, which is what the helper's own docstring said
    # was happening and was not: `_long_vertical_ends` had exactly one caller
    # and this function kept an unfiltered copy of the same scan. That left
    # the panel-edge and longest-at-each-edge filters protecting `_panel_y`
    # only -- and this is the WORSE site to leave open, because it takes
    # `min`: a confidence band whose edges reach ABOVE the axes needs no
    # missing spine at all. Widening fig26's band to +/-23pp returned an axes
    # top of 16.28 against a true 52.56, and every figure reported having no
    # title, which is round 23's symptom reached by a second route.
    box = _axes_box(stem)
    if box is not None:
        return box[0]
    vert = _long_vertical_ends(stem)
    if vert:
        return min(vert + horiz)
    # TWO ROWS OF RULES MEAN A TOP AND A BOTTOM; one row is only a bottom.
    if len({round(y, 1) for y in horiz}) >= 2:
        return min(horiz)
    return None


def _digits(text):
    """Normalise non-ASCII digit glyphs to ASCII.

    Exactly the problem `_mult` solves for the multiplication sign, at the
    sites that were left literal. `mathtext.fontset: stixsans` sets `$O_2$`
    with U+1D7E4 MATHEMATICAL SANS-SERIF DIGIT TWO, so `Normoxic (uniform
    O2)` never matches although the figure is correct and legible. The glyph
    a font chooses moves; the character it denotes does not.
    """
    # `unicodedata.digit`, NOT `.decimal`: `.decimal` RAISES on a character
    # that is `isdigit()` without being decimal -- a superscript two is the
    # everyday example -- so the guard would have crashed the reader on a
    # figure it was meant to read.
    def one(ch):
        if ch.isascii():
            return ch
        value = unicodedata.digit(ch, None)
        return ch if value is None else str(value)

    return "".join(one(ch) for ch in text)


def _mult(text):
    """Normalise the multiplication sign, whatever the font extracted it as.

    `mathtext.fontset: cm` embeds cmsy10 with no Unicode map, so PyMuPDF reads
    `$\\times$` as a pound sign and `$\\to$` as an exclamation mark. The figure
    renders correctly; only the extraction moves. Every score in fig25 is a
    number followed by that one sign, so the sign is normalised rather than
    matched literally.
    """
    # TOKEN-WISE, and only where the whole token is a decimal followed by one
    # sign. A blanket `(?<=\d)[^\d\s\w]` also rewrote `RSL3+FSP1i` into
    # `RSL3×FSP1i`, because the `+` in a pair name follows a digit too.
    return " ".join(
        re.sub(r"^(\(?)(\d+\.\d+)[^\d\s\w)](\)?)$",
               "\\1\\2" + "\u00d7" + "\\3", tok)
        for tok in text.split(" "))


def _assert_titled(panels, boxes, stem, titles, top):
    """Each title in `titles` must be centred over the panel of the same index.

    `top` is the y of the top of the plotting area (`_axes_top(...)`, which
    may be None when no spine establishes it); the title block is everything
    between the topmost title line and it. `stem` is a display label for the
    messages, so the caller supplies `top` rather than reading it back off a
    file name.

    A containment check on the extracted text cannot see a SWAP: exchanging
    two panels' titles leaves both strings present while `(b)` now captions
    the left panel and `(a)` the right, so every reference to "panel (a)"
    points at the wrong one.

    CENTRED OVER, not contained by. A title is centred on its axes and a long
    one overhangs both ends -- fig26's panel (a) title starts 33pt left of its
    own panel -- so the test is which panel centre it is nearest, which needs
    no tolerance at all.
    """
    assert len(panels) >= len(titles), (
        f"{stem} has {len(panels)} panels by its axis rules and "
        f"{len(titles)} titles to place")
    heads = [want.split()[0] for want in titles]           # `(a)` / `(b)`
    marks = [b for b in boxes if b[4] in heads]
    assert len(marks) == len(titles), (
        f"{stem} draws {len(marks)} panel-letter words, expected "
        f"{len(titles)}: {sorted(b[4] for b in marks)}")
    # THE TITLE BLOCK, NOT ONE ROW. Taking the words within 2pt of the
    # topmost panel letter rejected a WRAPPED title. matplotlib anchors a
    # title at its BOTTOM, so a two-line title puts its first line -- the line
    # carrying the panel letter -- above the other panel's letter, and the
    # single row then held one title where two were expected. That matters
    # more than an ordinary false failure: this same check correctly reports a
    # genuine title collision under a large `axes.titlesize`, and wrapping is
    # the natural fix for it, so the guard was rejecting its own remedy. The
    # block is every word between the topmost title line and the top of the
    # plotting area.
    row_y = min(b[1] for b in marks)
    # HOW FAR A WORD REACHES INTO THE AXES, AS A SHARE OF ITSELF. Requiring
    # the whole glyph box to clear the axes top by 1pt left the reader 2.49pt
    # from the edge on the committed figures, and `axes.titlepad` is a plain
    # rcParam defaulting to 6.0: at 3 the panel letter's parenthesis descender
    # legitimately crosses the top (`(a)` spans 35.01-53.95 against a top of
    # 52.56), every title word fell outside the band, and all three figures
    # were reported as having no title at all.
    #
    # The midpoint alone is too fine a line, and it fails in the direction
    # that matters: fig26's topmost tick label `50` is CENTRED on a tick
    # sitting exactly at the axes top, so it reaches 49.89% of its own height
    # in and its midpoint misses the boundary by 0.018pt. That near-half is
    # structural rather than lucky -- a label centred ON the boundary
    # straddles it -- so anything at or above about half is in-axes text.
    #
    # The line is drawn between the two measured extremes, each the deepest
    # panel-letter reach over all three figures:
    #
    #   committed                                    -7.86%  (clear of the axes)
    #   axes.titlepad: 3                              7.33%
    #   axes.titlepad: 0                              7.33%  (SATURATES here --
    #       matplotlib clamps the title's ink bottom to the axes top, so
    #       reducing the pad further changes nothing)
    #   axes.titlepad: 0 + axes.titley: 1.0          23.82%
    #   Courier New + titlepad: 0                    28.88%
    #   Courier New + titlepad: 0 + titley: 1.0      39.98%  (deepest seen)
    #   ---- 45% ----
    #   fig26's `50` tick label                      49.89%  (nearest in-axes)
    #
    # The two Courier rows are an UPPER BOUND on how deep a title can reach,
    # not figures that must pass: monospace widens fig26's titles until `(b)`
    # (x 549.1-572.5) is drawn over `open` (549.5-580.7), so those figures
    # fail the genuine-collision check anyway. `stays` is clear by 7.4pt here
    # -- `(b)` over `stays` is the `axes.titlesize: 16` case, and an earlier
    # version of this comment reused that pairing for this configuration. Bounding the cut with them errs toward accepting, which is the
    # safe direction here -- the thing on the other side is at 49.89%.
    #
    # An earlier version cut at a third and justified it with "at
    # `titlepad: 0` still only 23%", which is two errors: 0 gives 7.33%, and
    # 23.82% needs `axes.titley` set as well. Its replacement then said a
    # third "left 4.4 points of margin", which is 33.33 minus the 28.88 row
    # and does not follow from the 39.98 named in the same sentence -- against
    # that, a third leaves MINUS 6.65, i.e. it already rejected the deepest
    # title measured.
    #
    # None of this table decides the question any more. Depth cannot separate
    # the two classes at all, because a tick label's reach is a property of
    # `ytick.alignment` -- 49.89% at the default, 41.64% at `center`, 23.82%
    # at `baseline`, which is shallower than titles in this very table. The
    # cut is kept as a cheap upper bound on the band; POSITION is what
    # separates a title from a tick label, and that is done per RUN above.
    band = ([b for b in boxes
             if b[1] >= row_y - 0.5 and b[3] - top < 0.45 * (b[3] - b[1])]
            if top is not None
            else [b for b in boxes if abs(b[1] - row_y) < 2])
    rows = _cluster_rows(band, key=lambda b: b[1])
    groups, loose = [], []
    for r in rows:
        ws = sorted(r, key=lambda w: w[0])
        # A ROW HOLDS RUNS, AND A GAP ENDS ONE. Accumulating every word after
        # a panel letter into that letter's title assumed the only stray text
        # on a row comes BEFORE the first letter. It does not: wrap only the
        # RIGHT panel's title and its second line shares the row with the left
        # panel's whole title, sitting 173pt to its RIGHT, and was appended to
        # the left title -- which then reported panel (a) as mistitled when
        # panel (a) was byte-perfect. Words inside one title are a space
        # apart; two titles are a panel apart, so the gap separates them and
        # is measured against the row's own glyph height rather than a
        # constant, which keeps it true at any font size.
        heights = sorted(w[3] - w[1] for w in ws)
        gap = 1.5 * heights[len(heights) // 2]
        runs, run = [], [ws[0]]
        for prev, b in zip(ws, ws[1:]):
            if b[0] - prev[2] > gap:
                runs.append(run)
                run = [b]
            else:
                run.append(b)
        runs.append(run)
        # A RUN MUST OVERLAP A PANEL. Depth alone cannot separate a title from
        # a y-tick label, and the table below only made that look possible
        # because it was measured at one `ytick.alignment`. The label's reach
        # is a property of THAT setting, not of sitting on the boundary:
        # `center` (which `plt.style.use("classic")` sets) puts fig26's
        # topmost label at 41.64% and `baseline` at 23.82% -- shallower than
        # several legitimate titles, so no fraction exists that splits the two
        # classes. Position does: a y-tick label is drawn OUTSIDE its panel,
        # left of the spine it labels, while a title is centred ON its panel.
        #
        # By RUN rather than by word, because a title may overhang -- fig26's
        # panel (a) title starts 33pt left of its own panel, so its first word
        # overlaps nothing while the run it belongs to plainly does.
        runs = [r for r in runs
                if any(r[0][0] <= hi and r[-1][2] >= lo for lo, hi in panels)]
        for run in runs:
            cur = None
            for b in run:
                if b[4] in heads:
                    if cur is not None:
                        groups.append(cur)
                    cur = [b]
                elif cur is None:
                    # A WRAPPED TITLE'S SECOND LINE SHARES A ROW WITH THE
                    # NEXT PANEL'S LETTER, because a title is anchored at its
                    # BOTTOM: wrap panel (a)'s title and `under hypoxia`
                    # lands on the row holding `(b)`, to the left of it.
                    # Skipping rows carrying no letter at all was not enough.
                    loose.append(b)
                else:
                    cur.append(b)
            if cur is not None:
                groups.append(cur)
    assert len(groups) == len(titles), (
        f"{stem}'s title block splits into {len(groups)} titles, expected "
        f"{len(titles)}")
    # ORDERED BY PANEL LETTER, not by where the fragment landed: with a
    # wrapped title the letters sit on different rows, so reading order is no
    # longer left-to-right. Which panel each lettered title then captions is
    # still decided below, by the centre comparison, so this does not decide
    # the swap question -- it only says which expected title to compare with.
    groups.sort(key=lambda g: heads.index(g[0][4]))
    # A CONTINUATION LINE CARRIES NO PANEL LETTER, so each of its words joins
    # the title whose lettered fragment it sits nearest to horizontally.
    spans = [(min(b[0] for b in g), max(b[2] for b in g)) for g in groups]
    for b in loose:
        c = (b[0] + b[2]) / 2
        j = min(range(len(spans)),
                key=lambda k: max(spans[k][0] - c, 0.0, c - spans[k][1]))
        groups[j].append(b)
    groups = [sorted(g, key=lambda b: (round(b[1], 1), b[0])) for g in groups]
    groups = [[(b[0], b[2], b[4]) for b in g] for g in groups]
    centres = [(min(g_[0] for g_ in g) + max(g_[1] for g_ in g)) / 2
               for g in groups]
    panel_centres = [(lo + hi) / 2 for lo, hi in panels[:len(titles)]]
    for i, (want, c) in enumerate(zip(titles, centres)):
        nearest = min(range(len(panel_centres)),
                      key=lambda j: abs(panel_centres[j] - c))
        got = " ".join(t for _, _, t in groups[i])
        assert nearest == i, (
            f"{stem}'s title {got!r} is centred at x={c:.0f}, nearest panel "
            f"{nearest} rather than panel {i}. The titles are swapped, so "
            "every reference to that panel letter points at the other panel")
        # EXACT, NOT A PREFIX. `want.startswith(got)` admitted any truncation:
        # `(a) Kill collapse under hypoxia` cut to `(a) Kill`, or to `(a)`
        # alone, passed -- the whole panel claim deleted with only the
        # fingerprint firing. The branch existed to tolerate a wrapped title's
        # first line, but the group-count assertion above rejects a wrapped
        # title before reaching here, so it bought nothing.
        assert got == want, (
            f"{stem}'s panel {i} is titled {got!r}, expected {want!r}")


def _major_only(marks, axis_at=None):
    """From `(position, span_lo, span_hi, length)`, the major tick positions.

    TOUCHING THE AXIS IS THE FILTER, and only then is length used to separate
    major ticks from minor ones. Selecting purely by longest length was wrong
    twice over: the candidate set is every short stroke on the panel, not just
    ticks, so a `fill_between` band edge, the `closes ~day 3` arrow shaft or a
    legend frame's side could outrank a real tick and delete all nine.

    That was reachable from DATA alone, with no style change: widening RSL3's
    day-0 confidence interval from a half-width of 0.31 to 0.8 percentage
    points -- a purely stochastic quantity no fixture pins, 2.6x today's
    0.3063pp half-width (0.6126pp full) -- made the band's leftmost edge the
    longest stroke and the x axis read as empty. HALF-width in both figures,
    since the "2.6x" only holds on that reading and the bare numbers did not
    say which they were. It
    also re-broke `xtick.major.size: 1`, which the tick window had been
    widened to 0.5 specifically to accept, one commit earlier.

    A tick starts on the axis line; a band edge and an arrow do not. Among the
    strokes that do touch it, the majors are the longest -- minor ticks are
    drawn shorter, which is the only thing length is asked to decide.
    """
    if not marks:
        return []
    if axis_at is not None:
        # TOUCHING MEANS THE AXIS LIES IN THE STROKE'S SPAN, not that an
        # endpoint sits on it. `xtick.direction: inout` -- a plain rcParam --
        # straddles the axis, so at the default `major.size` 3.5 each endpoint
        # is 1.75pt away and an endpoint test discarded all nine ticks of a
        # perfectly correct figure. Span containment covers `out`, `in` and
        # `inout` alike, at any tick size.
        marks = [m for m in marks if m[1] - 1 <= axis_at <= m[2] + 1]
        if not marks:
            return []
    # TICKS ARE EVENLY SPACED, and that is what keeps the wider touch test
    # honest. Admitting a stroke that merely touches or crosses the axis lets
    # in data: under `axes.xmargin: 0` the first point sits ON the spine, so
    # the first segment of each near-flat curve starts there -- two strokes of
    # 79.71pt against six ticks of 3.5, and "the longest of the strokes that
    # touch" then read the curves as the scale. Two coincidental strokes are
    # not a tick series; six equally pitched ones are. Minor ticks are evenly
    # spaced too, which is why length still decides between the series that
    # qualify.
    #
    # A non-uniform axis (a log scale) has no regular series, and then this
    # falls back to "appears more than once", which is where it started.
    by_len = {}
    for m in marks:
        by_len.setdefault(round(m[3], 1), []).append(m)

    def _regular(ms):
        pos = sorted({round(m[0], 2) for m in ms})
        if len(pos) < 3:
            return False
        gaps = [b - a for a, b in zip(pos, pos[1:])]
        return max(gaps) - min(gaps) <= 0.05 * max(gaps)

    regular = {length for length, ms in by_len.items() if _regular(ms)}
    series = regular or {length for length, ms in by_len.items() if len(ms) > 1}
    longest = max(series or by_len)
    return sorted({m[0] for m in marks if round(m[3], 1) == longest})


def _axis_ticks_vertical(stem):
    """Short black tick marks on a horizontal axis (vertical strokes)."""
    pymupdf = _reader()
    doc = pymupdf.open(FIG_DIR / f"{stem}.pdf")
    try:
        cand = []
        for d in doc[0].get_drawings():
            col = d.get("color")
            # ANY COLOUR. Requiring black, then greyscale-and-dark, each
            # broke a shipped style: `axes.edgecolor` is freely settable and
            # matplotlib's own legend frame is greyscale 0.8. A tick mark is
            # identified by its SHAPE -- a short stroke perpendicular to an
            # axis -- not by how it is painted.
            if col is None:
                continue
            for item in d["items"]:
                if item[0] != "l":
                    continue
                a, b = item[1], item[2]
                if abs(a.x - b.x) < 0.2 and abs(b.y - a.y) >= 0.5:
                    cand.append((abs(b.y - a.y), a, b,
                                 tuple(round(c, 4) for c in col)))
        # SHORT RELATIVE TO THE AXIS, not shorter than 20pt. The absolute cap
        # was a measured constant sitting one tick-size away from ordinary
        # settings. Measured on all three figures: plain `xtick.major.size:
        # 20` draws strokes of exactly 20.0 (and 19.999996185302734), which
        # the inclusive cap accepted -- so `major.size: 20` alone was NOT the
        # false failure, and an earlier version of this comment said it was,
        # quoting a length no plain-size-20 figure produces. The reachable
        # cases are `major.size: 21` and up, and `xtick.direction: inout` at
        # size 20, which draws 20.000030517578125 and lost every tick.
        #
        # What is structurally true is that a tick decorates an axis and is a
        # fraction of it, so the longest stroke on the figure is never a tick
        # -- and `0.5 * longest` can only ever be MORE permissive than the old
        # 20, since it is floored there, so nothing this used to collect can
        # stop being collected.
        cap = max(20.0, 0.5 * max((c[0] for c in cand), default=0.0))
        return [(a.x, min(a.y, b.y), b.x, max(a.y, b.y), col, a, b)
                for length, a, b, col in cand if length <= cap]
    finally:
        doc.close()


def _axis_ticks(stem):
    """Short black axis tick marks as `(x0, y0, x1, y1, colour, p1, p2)`."""
    pymupdf = _reader()
    doc = pymupdf.open(FIG_DIR / f"{stem}.pdf")
    try:
        cand = []
        for d in doc[0].get_drawings():
            col = d.get("color")
            # ANY COLOUR. Requiring black, then greyscale-and-dark, each
            # broke a shipped style: `axes.edgecolor` is freely settable and
            # matplotlib's own legend frame is greyscale 0.8. A tick mark is
            # identified by its SHAPE -- a short stroke perpendicular to an
            # axis -- not by how it is painted.
            if col is None:
                continue
            for item in d["items"]:
                if item[0] != "l":
                    continue
                a, b = item[1], item[2]
                if abs(a.y - b.y) < 0.2 and abs(b.x - a.x) >= 0.5:
                    cand.append((abs(b.x - a.x), a, b,
                                 tuple(round(c, 4) for c in col)))
        # SHORT RELATIVE TO THE AXIS, not shorter than 20pt. The window was
        # `1 < len < 8` and failed at both ends on correct figures
        # (`ytick.major.size: 8` and `1.0` each left zero ticks and silently
        # disabled the scale readings); widening it to an absolute 20 left the
        # cap one tick-size from ordinary settings. See the vertical
        # extractor for the measurement: plain size 20 draws exactly 20.0 and
        # was accepted, so the reachable cases are size 21 and up, and
        # `direction: inout` at size 20. What is structurally true is that a
        # tick decorates an axis and is a fraction of it, so the longest
        # stroke on the figure is never a tick. Floored at the old 20, so
        # this can only be MORE permissive -- nothing it used to collect can
        # stop being collected.
        cap = max(20.0, 0.5 * max((c[0] for c in cand), default=0.0))
        return [(min(a.x, b.x), a.y, max(a.x, b.x), b.y, col, a, b)
                for length, a, b, col in cand if length <= cap]
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
    # +1pt, BECAUSE THE PANEL EDGE IS ROUNDED. `_panels` rounds each x to
    # one decimal, so an artist drawn exactly on the edge lands a few
    # ten-thousandths outside it. Under `axes.xmargin: 0` that is where the
    # first and last points sit, and the bare `<=` dropped fig26's ninth
    # tick -- eight marks against nine timepoints on a correct figure.
    rows = [[w for w in r if w[0] <= a_hi + 1] for r in _rows(wds).values()]
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
    want = f"{fx['synergy_score']:.2f}\u00d7"
    # BOUNDED ON PANEL (b)'S SPINE, not on its pair labels. Those labels are
    # y-ticks drawn OUTSIDE the spine, so a longer drug name pushes them left
    # towards panel (a)'s annotation: with real compound names the margin goes
    # from 81pt today to 1.0pt for `Liproxstatin-1+Ferrostatin-1` and to below
    # zero for `Liproxstatin-1+Deferoxamine`, failing a correct figure. The
    # spine moves the other way -- longer tick labels shrink panel (b), so the
    # bound gains room exactly when the labels need it. This file has already
    # been fixed twice for assuming today's pair names.
    b_spine = _panels("fig25_bliss_synergy")[-1][0]
    hits = [(x, y) for x, y, t in words_a
            if _mult(t) == want and x < b_spine]
    assert hits, (
        f"fig25's panel (a) does not annotate {want} anywhere left of panel "
        f"(b)'s axis at x={b_spine:.0f}. Panel (b) drawing that number as a "
        "bar label is not the same claim")
    # The two words sit 0.87pt apart vertically (57.115 and 57.985), so an
    # exact row key splits them; the annotation is one visual line.
    _ax, ay = hits[0]
    ay_mid = next((y0 + y1) / 2
                  for x0, y0, x1, y1, t in _word_bboxes("fig25_bliss_synergy")
                  if _mult(t) == want and abs(x0 - _ax) < 0.01
                  and abs(y0 - ay) < 0.01)
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
        # WHERE THE WORD'S MIDDLE IS. Neither edge separates these two cases,
        # because in both of them a word starts inside panel (a) and ends
        # outside it:
        #
        #   * bounding on the RIGHT edge clips the annotation's own text the
        #     moment it overhangs -- `annotate(..., fontsize=18)` pushes
        #     `synergy` past the edge and the run read `1.99x` alone;
        #   * bounding on the LEFT edge re-opened the case the comment above
        #     records, where panel (b)'s top pair label `RSL3+FSP1i` starts
        #     23.6pt INSIDE panel (a) at 17 pairs and the 20pt gap-walk
        #     absorbed it, reading `1.99x synergy RSL3+FSP1i`.
        #
        # It cannot widen to the inter-panel gap the way the fig24 and fig26
        # legend scans did, because that gap is where panel (b)'s labels are
        # drawn. The midpoint decides both: at 17 pairs the annotation's
        # `synergy` centres at 275.19 and panel (b)'s `RSL3+FSP1i` at 342.80,
        # against a panel edge of 330.2 -- 55pt of separation, not the 9pt an
        # earlier version of this comment implied by quoting 321.4, which is
        # `synergy`'s centre in a DIFFERENT figure (three pairs, annotation
        # fontsize 18). The other two numbers in that sentence were right.
        #
        # IT IS A TRADE, NOT AN IMPROVEMENT ON EVERY AXIS. Measured on the
        # same figure by sweeping the annotation's fontsize, first failure:
        #
        #   right edge (x1)   16      17 pairs: passes
        #   left edge  (x0)   32      17 pairs: FAILS
        #   midpoint          20      17 pairs: passes
        #
        # So the midpoint buys the pair-count case and gives back annotation
        # headroom, 32 down to 20. An earlier version of this comment called
        # it "strictly better" than the left-edge bound and cited 16 as that
        # bound's ceiling -- 16 belongs to the RIGHT-edge bound two rounds
        # earlier, so the number was attached to the wrong thing and the
        # conclusion it supported was the reverse of the measurement.
        #
        # Past fontsize 20 the annotation is more than half outside the panel
        # it annotates, which is a figure worth fixing rather than a reading
        # worth widening for.
        if abs((y0 + y1) / 2 - ay_mid) < 4 and (x0 + x1) / 2 <= a_hi + 1)
    at = next(i for i, (x0, _, t) in enumerate(same_row) if _mult(t) == want)
    lo = hi = at
    while lo > 0 and same_row[lo][0] - same_row[lo - 1][1] < 20:
        lo -= 1
    while hi + 1 < len(same_row) and same_row[hi + 1][0] - same_row[hi][1] < 20:
        hi += 1
    line = [(x0, t) for x0, _, t in same_row[lo:hi + 1]]
    assert _mult(" ".join(t for _, t in line)) == f"{want} synergy", (
        f"fig25's panel (a) synergy annotation reads "
        f"{' '.join(t for _, t in line)!r}, expected {want + ' synergy'!r}")
    # BY POSITION, for the reason fig24's needed it: a containment check does
    # not say WHICH panel carries the label, or that any panel does.
    # EACH OVER ITS OWN PANEL. Containment leaves a swap invisible: exchange
    # the two titles and both strings are still present while `(b)` captions
    # the left panel.
    _assert_titled(_panels("fig25_bliss_synergy"),
                   _word_bboxes("fig25_bliss_synergy"), "fig25",
                   ["(a) RSL3 + FSP1i: dual-pathway synergy",
                    "(b) Pairwise synergy (SDT pairs excluded: ceiling)"],
                   _axes_top("fig25_bliss_synergy"))

    # THE SUPTITLE, which states this figure's whole claim. Rewriting it to
    # say the dual-pathway depletion is ANTAGONISTIC -- the reverse of what
    # the synergy score and the additive reference line show -- passed every
    # check here. The same omission covered fig24's and fig26's headlines.
    assert "Dual-pathway (GPX4 + FSP1) depletion is synergistic" in text, (
        "fig25's suptitle no longer states that dual-pathway depletion is "
        "synergistic, which is the claim its synergy score and its additive "
        "reference line exist to support")

    cols = _vertical_labels(_word_bboxes("fig25_bliss_synergy"))
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
                       r"((?:\d+\.\d+\u00d7\s*)+)", _mult(text))
    assert labels, "fig25's panel (b) no longer lists its pairs and scores"
    pairs = labels.group(1).split()
    scores = re.findall(r"\d+\.\d+\u00d7", labels.group(2))
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
    vals = [(x, y, _mult(t)) for x, y, t in words
            if re.fullmatch(r"\d+\.\d+\u00d7", _mult(t))]
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
    # SPLIT AT THE GAP THAT STRADDLES THE PANEL BOUNDARY, not at the widest
    # one. "Widest" is a property of the committed layout: with
    # `gridspec_kw={"width_ratios": [1.6, 1]}` -- a plain option -- panel (a)
    # grows, a gap inside ITS OWN labels becomes the widest in the row, and
    # the split landed there, reading the caption as `alone expected
    # combination Bliss synergy score (observed / expected)` on a correct
    # figure. Which gap separates the two panels' text is not a matter of
    # size: it is the one spanning the space between the panels, and
    # `_panels` already says where that is.
    row_words = _row_containing(_boxed_rows("fig25_bliss_synergy"),
                                caption[0][1])
    fig25_panels = _panels("fig25_bliss_synergy")
    # BOTH CONDITIONS, because either alone is wrong. A gap must DOMINATE the
    # row's word spacing to be a separation at all -- at `axes.labelsize: 16`
    # panel (a)'s labels move to another row entirely and this one holds only
    # the caption, so a rule that always splits drops the caption's first
    # word. And among the gaps that do dominate, the separating one is the
    # one at the PANEL BOUNDARY, not the widest: at `width_ratios: [1.6, 1]`
    # panel (a) grows until a gap inside its own labels is the widest in the
    # row.
    cap_row = [t for _, _, t in row_words]
    gaps = [(row_words[i + 1][0] - row_words[i][1], i)      # EDGE to edge
            for i in range(len(row_words) - 1)]
    if gaps and len(fig25_panels) > 1:
        typical = sorted(g for g, _ in gaps)[len(gaps) // 2]
        wide = [(g, i) for g, i in gaps if g > max(20.0, 3 * typical)]
        if wide:
            between = (fig25_panels[0][1] + fig25_panels[1][0]) / 2
            _, at = min(wide, key=lambda gi: max(
                row_words[gi[1]][1] - between, 0.0,
                between - row_words[gi[1] + 1][0]))
            cap_row = [t for _, _, t in row_words[at + 1:]]
    assert " ".join(cap_row) == "Bliss synergy score (observed / expected)", (
        f"fig25's panel (b) x axis reads {' '.join(cap_row)!r}; the values "
        "compared here are the observed-over-expected ratio, and the "
        "reciprocal is a different claim about every bar on the panel")
    # THE ADDITIVE LINE, on panel (b)'s own scale. `axvline(1.0)` is what
    # turns a score into a synergy claim, and nothing read it: moving it to
    # 1.9 with its label unchanged makes two of the three pairs read
    # sub-additive, and only the fingerprint noticed. Panel (b)'s x tick
    # labels give the scale, so the line's position is checked in DATA units
    # rather than points.
    # THE CENTRE IS NOW ONLY AN X FILTER. This block used to take the label's
    # centre as its POSITION, because a label is centred on its tick and using
    # `x0` put the additive line at 1.08. Position comes from the tick marks
    # now, so what survives here is `cx` deciding which labels belong to panel
    # (b) and, below, which mark each label sits nearest.
    # `\d+(?:\.\d+)?`, NOT `\d+\.\d`. Panel (b)'s x limit is
    # `max(scores) * 1.25`, so once the top score passes ~3.2 matplotlib
    # switches to INTEGER tick labels and a one-decimal pattern matches none of
    # them. Measured: raising RSL3+FSP1i's synergy to 3.5 -- a stronger result,
    # the direction this figure argues for -- made a correct figure fail with
    # "0 numeric x tick labels". This file already paid for that exact mistake
    # once, in the pair-name regex a few blocks up.
    # SCALE FROM THE TICK MARKS, NOT THE TICK LABELS -- the third scale
    # reading to need this and the one that was left behind. fig24's y scale
    # was moved to the marks because a label's bbox centre carries the font's
    # digit-glyph asymmetry; fig25's x scale kept reading labels, and the bias
    # is not only the font's. `xtick.alignment: left` shifts every label
    # relative to its tick, and the additive line then read 0.92 instead of
    # 1.0 -- a correct figure accused of moving the line every score on the
    # panel is compared against. Marks are geometry and do not move.
    b_lo, b_hi = fig25_panels[-1]
    bottom = _panel_y("fig25_bliss_synergy")[1]
    mark_xs = _major_only(
        [(round(a.x, 2), min(a.y, b.y), max(a.y, b.y), round(abs(b.y - a.y), 2))
         for x0, y0, x1, y1, col, a, b
         in _axis_ticks_vertical("fig25_bliss_synergy")
         if b_lo - 1 <= a.x <= b_hi + 1],
        axis_at=bottom)
    labels = sorted((cx, float(t))
                    for cx, y, t in _word_centres("fig25_bliss_synergy")
                    if re.fullmatch(r"\d+(?:\.\d+)?", t) and cx > b_lo - 20)
    assert len(mark_xs) >= 2 and len(mark_xs) == len(labels), (
        f"fig25's panel (b) shows {len(mark_xs)} x tick marks against "
        f"{len(labels)} numeric labels; its scale cannot be recovered")
    # PAIR EACH MARK WITH THE LABEL NEAREST IT, NOT BY RANK. Taking the mark
    # positions was right -- a label's own x carries the font's glyph bias and
    # moves with `xtick.alignment`. Sorting the VALUES and zipping them onto
    # the marks was not: it manufactures the ascending order that the comment
    # here then claimed to have observed, so a label sitting on the wrong tick
    # became invisible. Relabelling panel (b)'s ticks `0.0 0.5 2.0 1.5 1.0`
    # left the additive line under a label reading 2.0 -- every score on the
    # panel compared against a line the axis says is somewhere else -- and
    # passed every semantic check, where the label-based reading it replaced
    # had caught it.
    #
    # Nearest-in-x keeps the association the figure actually draws while still
    # taking position from the marks.
    ticks = [(mx, min(labels, key=lambda lv: abs(lv[0] - mx))[1])
             for mx in mark_xs]
    assert len({min(labels, key=lambda lv: abs(lv[0] - mx))[0]
                for mx in mark_xs}) == len(mark_xs), (
        "fig25's panel (b) has two tick marks claiming the same label; the "
        "axis scale cannot be attributed")
    drawn_values = [v for _, v in ticks]
    assert all(b > a for a, b in zip(drawn_values, drawn_values[1:])), (
        f"fig25's panel (b) x tick labels read {drawn_values} left to right "
        "and must ascend. A label on the wrong tick moves the additive line "
        "relative to the scores without moving either")
    (tx0, tv0), (tx1, tv1) = ticks[0], ticks[-1]
    scale = (tx1 - tx0) / (tv1 - tv0)
    dashed = sorted(x for x, y0, y1, col
                    in _dashed_verticals("fig25_bliss_synergy")
                    if fig25_panels[-1][0] <= x <= fig25_panels[-1][1])
    # A GRID DRAWS ONE LINE AT EVERY TICK; the threshold draws one wherever
    # its value falls. `axes.grid: True` with `grid.linestyle: '--'` -- and
    # every style that presets it, `bmh` among them -- therefore added five
    # more dashed verticals and a correct figure read as having five spurious
    # thresholds. The threshold's own value is 1.0, which IS a tick here, so
    # "drop the ones at tick positions" would drop the threshold too; what
    # separates them is that a grid covers ALL the ticks, so it is only
    # subtracted when all of them are covered.
    tick_xs = [cx for cx, _ in ticks]
    gridded = all(any(abs(x - cx) < 1 for x in dashed) for cx in tick_xs)
    vlines = list(dashed)
    if gridded:
        for cx in tick_xs:
            hit = next((x for x in vlines if abs(x - cx) < 1), None)
            if hit is not None:
                vlines.remove(hit)
    assert len(vlines) == 1, (
        f"fig25's panel (b) draws {len(vlines)} dashed vertical lines that "
        f"are not grid lines, expected one additive threshold")
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
    assert f"additive ({1.0:.1f}\u00d7)" in _mult(
            text.replace("$\\times$", "\u00d7")), (
        "fig25's panel (b) threshold legend does not read `additive (1.0×)`. "
        f"The line itself is drawn at {drawn_at:.2f} on the panel's own "
        "scale, so a legend naming any other value contradicts it")

    assert caption[0][0] > max(x for x, _, t in names), (
        "fig25's score caption is not to the right of its pair labels, so it "
        "is not panel (b)'s x axis")
    assert by_row[key] == f"{fx['synergy_score']:.2f}\u00d7", (
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
    # THE COLLAPSE PAIR IS THE PERCENTAGES FLANKING THE ARROW. Identifying
    # the pair instead as "the row that holds two of them" avoided naming a
    # glyph, and put a property of today's four DISTINCT bar values in place
    # of a property of the figure. Give SDT's four gradient rates its uniform
    # rate -- which is what fig24's own footnote says the model assumes,
    # "SDT is modeled as O$_2$-independent" -- and two bar annotations sit at
    # one height, so a correct, fully legible figure was rejected. Values
    # within ~0.04pp of each other collided the same way.
    #
    # The arrow is pinned rather than avoided, for the reason it is pinned in
    # the collapse-sentence check below: it is the only thing in the figure
    # carrying the DIRECTION of the collapse. `mathtext.fontset: cm` extracts
    # it as `!` and fails here loudly; that is the stated limit, and a loud
    # failure on an unusual font beats being blind to a reversed arrow.
    boxes_a = _word_bboxes("fig24_hypoxia_killcurve")
    arrows = [b for b in boxes_a if b[4] == "\u2192"]
    assert len(arrows) == 1, (
        f"expected fig24 panel (a) to draw exactly one collapse arrow, "
        f"found {len(arrows)}")
    ax0, ay0, ax1, ay1 = arrows[0][:4]
    # VERTICAL OVERLAP, NOT A SHARED ROW. Clustering by the words' top edges
    # separates the arrow from the very percentages it sits between, because
    # a glyph box's height is a property of the FONT: under `seaborn-v0_8-
    # ticks` the arrow spans y 152.44-166.01 while `3.7%` and `0.1%` span
    # 154.26-164.91, so their tops are 1.82pt apart and any row tolerance
    # tight enough to be useful splits them. Their extents nevertheless
    # overlap -- the arrow's box CONTAINS the digits' -- and overlap is what
    # "on the same line" actually means.
    same_line = [b for b in boxes_a if b[1] < ay1 and b[3] > ay0]
    pcts = [b for b in same_line if re.fullmatch(r"\d+\.\d%", b[4])]
    # FLANKING, not merely on the line: a bar annotation that happens to land
    # at the annotation's height must not displace either half of the pair.
    before = [b for b in pcts if b[2] <= ax0]
    after = [b for b in pcts if b[0] >= ax1]
    assert before and after, (
        f"fig24's collapse arrow is not flanked by two percentages; its line "
        f"reads {[b[4] for b in sorted(pcts)]}")
    lb, rb = max(before, key=lambda b: b[0]), min(after, key=lambda b: b[0])
    left, right = (lb[0], lb[4]), (rb[0], rb[4])

    bars = sorted(((x, t) for x, y, t in panel_a
                   if re.fullmatch(r"\d+\.\d%", t)
                   and (x, t) not in (left, right)),
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
    series = re.search(r"Normoxic \(uniform O2\)\s+Hypoxic \(O2 gradient\)",
                       _digits(text))
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
    # SCOPED TO PANEL (a). The scan is page-wide and panel (b) draws its
    # `o-`/`s-` markers as filled rectangles in the SAME two series colours,
    # at the bottom of a panel that shares fig24's y range. `lines.markersize:
    # 14` alone -- no generator edit -- made it report six bars on the
    # baseline, and `axB.set_ylim(0, 100)` put four panel-(b) markers 2.81pt
    # BELOW the real bars, so `max(y1)` picked them and every real bar was
    # excluded. The text scan beside this was scoped for exactly this reason
    # and the rectangle scan was not.
    # UP TO WHERE PANEL (b) STARTS, for the same reason as fig26's legend
    # scan: a legend moved out of the axes with `bbox_to_anchor=(1.02, 0.5)`
    # draws its swatches in the gap between the panels, and a scope ending at
    # panel (a)'s own right edge lost them -- a correct figure reported that
    # its `Normoxic` entry had no swatch. The gap holds nothing else, and the
    # bars are still selected by standing on the axis.
    fig24_panels = _panels("fig24_hypoxia_killcurve")
    a_right = (fig24_panels[1][0] - 1 if len(fig24_panels) > 1
               else fig24_panels[0][1] + 1)
    rects = [r for r in _filled_rects("fig24_hypoxia_killcurve")
             if r[1] <= a_right
             # THE AXES BACKGROUND IS NOT A BAR. matplotlib fills the plotting
             # area with a rectangle spanning the panel EXACTLY, and its
             # bottom edge is the axis -- so it silently became "the lowest
             # thing standing on the axis" and the bars, which stand on data
             # zero, were all 11.34pt above it and matched nothing. Nothing
             # that spans a whole panel is a bar or a swatch.
             and not any(r[0] <= lo + 1 and r[1] >= hi - 1
                         for lo, hi in fig24_panels)]
    # THE AXIS ITSELF, not the bottom-most filled edge. The bars stand on the
    # axis and a legend swatch does not, which is the distinction being drawn
    # -- but deriving the baseline from the rectangles assumes nothing is
    # drawn BELOW the bars, and a legend placed under the axes
    # (`bbox_to_anchor=(0.5, -0.22)`, the ordinary way to put one there) is
    # exactly that. With `ncol=2` both swatches then share the bottom-most
    # edge, so both were excluded as "the baseline" and neither entry
    # resolved, on a correct figure. The bottom spine is where the bars
    # stand, whatever else the figure draws.
    # AT OR ABOVE THE AXIS, not ON it. Deriving the baseline from `max(r[3])`
    # broke when a legend was placed BELOW the axes; replacing it with the
    # axis itself broke the opposite way, because bars stand on DATA ZERO and
    # zero is only on the axis when the y limit starts there. `ylim(-5, 105)`
    # -- the convention this same generator already uses on fig24 panel (b)
    # and both fig26 panels -- puts all four bars 11.34pt above the spine and
    # the guard reported no bars at all.
    #
    # What holds in both cases: the bars share one bottom edge, and nothing
    # BELOW the axis is a bar. So take the lowest edge that is not below the
    # axis -- which is the bars' own, whether or not it touches it.
    _axis_y = _panel_y("fig24_hypoxia_killcurve")[1]
    _standing = [r[3] for r in rects if r[3] <= _axis_y + 0.5]
    assert _standing, "fig24 panel (a) draws no filled rectangle above its axis"
    bar_baseline = max(_standing)
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
                          # NOT ON THE BASELINE, rather than under a width
                          # bound. `< 52` was a measured midpoint between the
                          # swatch (16pt) and the bars (53pt), and it rejected
                          # `legend.handlelength` 8, whose swatch is 80pt.
                          # What actually separates them is that a bar stands
                          # on the axis and a swatch floats.
                          if abs(r[2] - y) < 8
                          and abs(r[3] - bar_baseline) > 1
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
    # THE SAME BASELINE, for the same reason as `bar_baseline` above. This
    # site had the identical shape -- a baseline derived from the rectangles
    # -- and the round that fixed the swatch site did not look for the
    # sibling: with the legend below the axes its two swatches ARE in the
    # series colours and sit lower than the bars, so `max` made them the
    # baseline and every real bar was excluded.
    bars_drawn = sorted((r for r in rects
                         if r[4] in series
                         and abs(r[3] - bar_baseline) < 0.5),
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

    # THE PANEL TITLES, EACH OVER ITS OWN PANEL. A panel caption can
    # contradict the suptitle this file asserts and the check beside it, so
    # `(a) Kill collapse under hypoxia` rewritten to `Kill preserved` needed
    # gating. But `in text` alone was not enough either: SWAPPING the two
    # titles between the axes leaves both strings present, puts `(b)` over the
    # left panel and `(a)` over the right, and every manuscript reference to
    # "panel (a)" then points at the wrong figure. This file's own rule is
    # that anything identifying a value is compared by geometry.
    _assert_titled(_panels("fig24_hypoxia_killcurve"),
                   _word_bboxes("fig24_hypoxia_killcurve"), "fig24",
                   ["(a) Kill collapse under hypoxia",
                    "(b) Robust across penetration length"],
                   _axes_top("fig24_hypoxia_killcurve"))

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
    columns = _vertical_labels(_word_bboxes("fig24_hypoxia_killcurve"))
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
        assert abs(col_mid - axis_mid) <= (y_hi - y_lo) * 0.08, (
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
    marks_raw = [(round(a.y, 3), min(a.x, b.x), max(a.x, b.x),
                  round(abs(b.x - a.x), 2))
                 for x0, y0, x1, y1, col, a, b in _axis_ticks(
                     "fig24_hypoxia_killcurve")
                 # THE SPINE LIES IN THE STROKE, rather than the stroke's left
                 # end sitting within 8pt of it. `ytick.major.size: 20` puts
                 # that end 20pt out, so every tick was discarded -- and the
                 # test still passed, because the label fallback deleted above
                 # silently caught it at a fifteenth of the precision. Same
                 # rule as `_major_only`'s, so `in` / `out` / `inout` and any
                 # tick size all read alike.
                 if min(a.x, b.x) - 1 <= b_x0 <= max(a.x, b.x) + 1]
    labels = sorted(
        ((y0 + y1) / 2, float(t))
        for x0, y0, x1, y1, t in _word_bboxes("fig24_hypoxia_killcurve")
        if re.fullmatch(r"\d+", t) and b_x0 - 40 < x1 < b_x0 + 2)
    # MAJOR TICKS ONLY. `xtick.minor.visible: True` is a plain rcParam and
    # gave 21 marks -- 15 minor and 6 major -- against 6 numeric labels on
    # this axis. (The "46 against 9 timepoints" an earlier version quoted here
    # is fig26's x-axis measurement; fig24 panel (b) has no timepoints.) Minor
    # ticks are drawn shorter, so the majors are the longest of the strokes
    # that touch the axis.
    marks = _major_only(marks_raw, axis_at=b_x0)
    # NO LABEL FALLBACK. This used to swap in the label centres when no marks
    # were drawn, widening the bound 0.1pp -> 1.5pp to absorb the font bias
    # that swap reintroduces. The widening is a 15x loss of precision and it
    # was applied SILENTLY: under `ytick.major.size: 0` a reference line drawn
    # 1.1pp off its value -- 91.9% rendered as 90.8% -- passed with the whole
    # suite green. The bullet at the top of this file said BOTH scale readings
    # fail loudly on a zero tick size, and only fig26's did. A guard that
    # keeps running at a fifteenth of its precision without saying so is worse
    # than one that stops, so this one stops.
    tolerance = 0.1
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
        assert abs(drawn - want) <= tolerance, (
            f"fig24's panel (b) draws {name}'s normoxic reference at "
            f"{drawn:.1f}% on its own y scale and the data gives "
            f"{want:.1f}%. That line is what the depth curves are compared "
            "against, so it has to be the value it claims")

    # PANEL (b)'S LAMBDA AXIS, in ascending order. Reversing those tick labels
    # -- so the axis reads 150/120/100/80 while the curves are plotted 80 to
    # 150 -- survived every check here, which is the same tick-array inversion
    # already closed on panel (a) and on fig26. The values are in the
    # fixture's own condition names.
    lam = sorted(int(re.search(r"gradient_(\d+)um", r["o2_condition"]).group(1))
                 for r in {r["o2_condition"]: r for r in rows
                           if r["o2_condition"].startswith("gradient_")}.values())
    assert len(lam) >= 2, (
        f"hypoxia_killcurve_rows.json names {len(lam)} gradient conditions; "
        "fig24's panel (b) plots one point per lambda")
    lam_labels = [t for _, t in sorted(
        ((x0 + x1) / 2, t)
        for x0, y0, x1, y1, t in _word_bboxes("fig24_hypoxia_killcurve")
        # BELOW THE PANEL MIDLINE, not below the spine. `xtick.direction: in`
        # lifts the labels 3.5pt and put them above a spine they clear by only
        # 3.15pt outward, so a correct axis read as empty and was reported
        # reversed.
        # BY ITS CENTRE, not its left edge. A tick label is centred on its
        # tick and every tick lies on the axis, so the centre is in the
        # panel by construction while the left edge is not: under
        # `axes.xmargin: 0` the first point sits ON the spine and `80`
        # overhangs it, so the label was dropped and the axis read as
        # reversed on a correct figure.
        if re.fullmatch(r"\d+", t) and b_x0 - 1 <= (x0 + x1) / 2 <= b_x1 + 1
        and y0 > (y_lo + y_hi) / 2)]
    assert lam_labels == [str(v) for v in lam], (
        f"fig24's panel (b) x axis reads {lam_labels} LEFT TO RIGHT and the "
        f"fixture's gradient conditions are {lam}. The curves are plotted in "
        "ascending lambda, so a reversed axis reads the depth trend backwards")

    # And the collapse annotation must be the ratio of the two it names, not a
    # number carried along beside them.
    # THE ARROW IS PINNED, BECAUSE IT CARRIES THE DIRECTION. A previous
    # version matched `\S` here so that any font's rendering would do -- and
    # that silently removed the only check on which way the collapse runs:
    # `3.7% <- 0.1%` reads as hypoxia RAISING the kill rate from 0.1% to 3.7%,
    # the #790 inversion, and it passed. The two values are already compared
    # in order, but both stay in order when only the arrow flips.
    #
    # So the rightward arrow is required literally. `mathtext.fontset: cm`
    # extracts it as `!` and will fail here; that is a stated limit below,
    # and a loud failure on an unusual font is the right trade against being
    # blind to a reversed arrow.
    collapse = re.search(
        r"(\d+\.\d%)\s*\u2192\s*(\d+\.\d%)\s*\(~(\d+)\s*"
        r"\u00d7\s*collapse\)", text)
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
    # SCOPED BY PANEL, NOT BY LINE. Two earlier rules each failed one way:
    # nearest-anywhere let a panel title beat the caption's own unit at large
    # label sizes, and restricting to the caption's own LINE re-opened the
    # very hole this block exists to close -- wrap panel (a)'s caption onto
    # two lines and change its unit to `(hours)` and the search skipped it,
    # answering from panel (b)'s `(days)`, which is the cross-panel answer.
    # A caption and its unit belong to the same PANEL, which is an x range
    # `_panels` already gives us, and that holds whether the caption wraps or
    # not.
    stem26 = "fig26_vulnerability_window"
    p26 = _panels(stem26)
    panel_of = lambda xx: next((i for i, (lo, hi) in enumerate(p26)
                                if lo - 30 <= xx <= hi + 30), None)
    units = [(x, y, t) for x, y, t in wds if re.fullmatch(r"\(\w+\)", t)]
    subjects = [(x, y, t) for x, y, t in wds if t == "post-chemotherapy"]
    assert len(subjects) == 2, (
        f"fig26 draws {len(subjects)} x-axis captions, expected one per panel")
    for x, y, _ in subjects:
        panel = panel_of(x)
        assert panel is not None, (
            f"fig26's x-axis caption at x={x:.0f} is not inside either panel")
        # BELOW the caption's own row and within its own panel: a wrapped
        # caption puts the unit on the next line down, never on the other
        # panel.
        here = [u for u in units
                if panel_of(u[0]) == panel and u[1] >= y - 2]
        assert here, (
            f"fig26's panel {panel} x-axis caption draws no parenthesised "
            "unit; every panel states the unit its ticks are in")
        nearest = min(here, key=lambda u: (u[0] - x) ** 2 + (u[1] - y) ** 2)
        assert nearest[2] == "(days)", (
            f"fig26's panel {panel} x-axis caption is nearest the unit "
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
    a_lo, a_hi = p26[0]
    # A PANEL'S LEGEND BELONGS TO IT EVEN WHEN DRAWN BESIDE IT. Scoping this
    # to panel (a)'s own x range assumed the legend is inside the axes, and
    # `bbox_to_anchor=(1.02, 0.5)` -- the standard way to move one out -- puts
    # both entries in the gap between the panels, so a correct figure reported
    # no legend captions at all. The region that belongs to panel (a) runs up
    # to where panel (b) starts; only the LEGEND scan is widened, because the
    # curve, tick and shading scans read things that are inside the axes by
    # construction and have nothing to gain from the gap.
    legend_hi = p26[1][0] - 1 if len(p26) > 1 else a_hi + 1
    doc_lines = _hlines_any("fig26_vulnerability_window")
    legend_col = {}
    for x, y, t in wds:
        if t not in ("SDT", "RSL3") or not (a_lo - 1 <= x <= legend_hi):
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
        # NEAREST WINS, and a candidate must be genuinely beside its label.
        # The loop takes any `SDT`/`RSL3` word inside panel (a) -- which
        # includes the in-panel `RSL3 window open` annotation -- and
        # `setdefault` locked whichever came first. At legend fontsize 10 and
        # up the annotation's `RSL3` sat 6.02pt from the SDT sample and
        # claimed red, while the real entry 0.51pt from its own blue sample
        # was ignored, so a clean figure was told it drew one colour twice.
        if sample:
            dist = abs(sample[0][1] - y_mid)
            if t not in legend_col or dist < legend_col[t][0]:
                legend_col[t] = (dist, sample[0][2])
    legend_col = {k: v[1] for k, v in legend_col.items()}
    assert set(legend_col) == {"SDT", "RSL3"}, (
        f"fig26's panel (a) legend captions {sorted(legend_col)}, expected "
        "SDT and RSL3, each with a line sample beside it")
    assert legend_col["SDT"] != legend_col["RSL3"], (
        "fig26 draws both panel (a) legend samples in one colour")
    ends = {}
    for x, y, col in _stroke_points("fig26_vulnerability_window"):
        if not (a_lo - 1 <= x <= a_hi + 1) or col not in legend_col.values():
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
    _, a_bottom = _panel_y(stem26)
    xmarks = _major_only([(round(a.x, 2), min(a.y, b.y), max(a.y, b.y),
                           round(abs(b.y - a.y), 2))
                          for x0, y0, x1, y1, col, a, b
                          in _axis_ticks_vertical(stem26)
                          if a.x <= a_hi + 1],
                         axis_at=a_bottom)
    assert len(xmarks) == len(expected), (
        f"fig26 panel (a) shows {len(xmarks)} x tick marks against "
        f"{len(expected)} timepoints")
    win_end = max(i for i, d in enumerate(days) if d <= 3.0)
    shade = [r for r in _filled_rects(stem26)
             if r[0] >= a_lo - 1 and r[1] <= a_hi + 1
             and (r[1] - r[0]) > 20 and (r[3] - r[2]) > 100
             # NOT THE AXES PATCH. Excluding white by literal meant any
             # `axes.facecolor` -- ggplot, Solarize_Light2, bmh,
             # dark_background -- counted the panel's own background as a
             # second span. The patch spans the panel EXACTLY; the window
             # shade is strictly inside it.
             and not (abs(r[0] - a_lo) < 1 and abs(r[1] - a_hi) < 1)]
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

    # EACH TITLE OVER ITS OWN PANEL. Panel (a)'s title names which treatment
    # closes, which is the claim the curve check above compares against -- and
    # a containment check cannot see the titles SWAPPED, which puts that claim
    # over the wrong panel while both strings remain present.
    _assert_titled(p26, _word_bboxes(stem26), "fig26",
                   ["(a) Treatment window: RSL3 closes, SDT stays open",
                    "(b) Why: GPX4 re-expression closes the window"],
                   _axes_top(stem26))

    # THE Y-AXIS LABELS, BOUND TO THEIR AXES. fig24's and fig25's equivalents
    # are gated with an explicit "the unit the bars are in" rationale; fig26's
    # were not, so `Persister kill (%)` could become `survival (%)` and pass.
    # Containment is not enough here either: panel (b) carries TWO y axes, and
    # exchanging `RSL3 kill (%)` with `Mean GPX4 (recovered fraction)` leaves
    # the blue 0-50 kill axis labelled as a recovered fraction and the green
    # 0.5-0.8 GPX4 axis labelled as a kill percentage, with both strings still
    # on the page. A rotated label is a column at one x, so each is required
    # at its own.
    fig26_cols = _vertical_labels(_word_bboxes(stem26))
    by_text = {}
    for (cx0, _cx1), (lbl, _lo, _hi) in fig26_cols.items():
        by_text.setdefault(lbl, []).append(cx0)
    for want in ("Persister kill (%)", "RSL3 kill (%)",
                 "Mean GPX4 (recovered fraction)"):
        assert want in by_text, (
            f"fig26 no longer draws {want!r} as a y-axis label; those labels "
            f"say what the plotted values are. Found: {sorted(by_text)}")
    # Left to right: panel (a)'s kill axis, panel (b)'s kill axis, then the
    # twin GPX4 axis on panel (b)'s right.
    order = sorted((min(by_text[w]), w) for w in
                   ("Persister kill (%)", "RSL3 kill (%)",
                    "Mean GPX4 (recovered fraction)"))
    assert [w for _, w in order] == ["Persister kill (%)", "RSL3 kill (%)",
                                     "Mean GPX4 (recovered fraction)"], (
        f"fig26's y-axis labels read {[w for _, w in order]} LEFT TO RIGHT. "
        "Panel (b)'s left axis carries the RSL3 kill curve and its right axis "
        "the GPX4 fraction; exchanged, each axis is labelled with the other's "
        "quantity")

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
            # THE PIXELS, NOT JUST THE METADATA. The raster's dimensions come
            # from the figure geometry rather than from the data, so a heatmap
            # redrawn from DIFFERENT numbers keeps the same size, colourspace
            # and placement and hashed identically -- verified with two random
            # seeds. Latent for these three figures, which draw no raster, but
            # fig17 draws its heatmaps with `imshow`, which is the case this
            # backstop names as reachable.
            try:
                raw = doc.extract_image(info[0])["image"]
                digest = hashlib.sha256(raw).hexdigest()[:16]
            except Exception:
                digest = None
            images.append((info[1:], rects, digest))
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
