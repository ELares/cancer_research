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

BINDING IS POSITIONAL where the artifact gives no alternative. A matplotlib bar
annotation carries no label in the text stream, so for fig24 and fig25's panel
(a) the order is the bar order, documented per figure below. Where labels ARE
adjacent -- fig25's panel (b) -- each value is bound to its own label, because
a set comparison passes when two values are swapped between subjects, which is
the defect the fig17 guard was fixed for.
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
    assert f"{fx['synergy_score']:.2f}×" in text, (
        f"fig25 does not draw the synergy score {fx['synergy_score']:.2f}x "
        "the fixture holds")


def test_fig25_binds_each_pair_to_its_own_score():
    """Panel (b) labels each bar, so a set comparison is not good enough: two
    scores swapped between pairs would pass one."""
    fx = json.loads((FIXTURES / "bliss_synergy.json").read_text())["rsl3_fsp1i"]
    text = _drawn("fig25_bliss_synergy")
    labels = re.search(
        r"((?:FSP1i\+HDACi|RSL3\+HDACi|RSL3\+FSP1i)(?:\s+(?:FSP1i\+HDACi|"
        r"RSL3\+HDACi|RSL3\+FSP1i))*)\s+((?:\d+\.\d+×\s*)+)", text)
    assert labels, "fig25's panel (b) no longer lists its pairs and scores"
    pairs = labels.group(1).split()
    scores = re.findall(r"\d+\.\d+×", labels.group(2))
    assert len(pairs) == len(scores), (
        f"fig25 draws {len(pairs)} pairs and {len(scores)} scores")
    mapping = dict(zip(pairs, scores))
    assert mapping.get("RSL3+FSP1i") == f"{fx['synergy_score']:.2f}×", (
        f"fig25 gives RSL3+FSP1i {mapping.get('RSL3+FSP1i')} and the fixture "
        f"says {fx['synergy_score']:.2f}x -- the label and the score it is "
        "drawn beside disagree with the data")


def test_fig24_draws_its_hypoxia_numbers():
    """Panel (a): RSL3 normoxic, RSL3 hypoxic, SDT normoxic, SDT hypoxic.

    `uniform` is the normoxic arm and `gradient_120um` the hypoxic one -- the
    same pair the caption quotes.
    """
    rows = json.loads(
        (FIXTURES / "hypoxia_killcurve_rows.json").read_text())["conditions"]
    by = {(r["treatment"], r["o2_condition"]): r["overall_kill_rate"]
          for r in rows}
    expected = [_pct(by[("RSL3", "uniform")]), _pct(by[("RSL3", "gradient_120um")]),
                _pct(by[("SDT", "uniform")]), _pct(by[("SDT", "gradient_120um")])]
    text = _drawn("fig24_hypoxia_killcurve")
    found = re.findall(r"\d+\.\d%", text)
    assert found[:4] == expected, (
        f"fig24 draws {found[:4]} and hypoxia_killcurve_rows.json gives "
        f"{expected}")
    # And the collapse annotation must be the ratio of the two it names, not a
    # number carried along beside them.
    collapse = re.search(r"(\d+\.\d%) → (\d+\.\d%) \(~(\d+)× collapse\)", text)
    assert collapse, "fig24 no longer annotates the collapse ratio"
    assert [collapse.group(1), collapse.group(2)] == expected[:2], (
        f"fig24's collapse annotation names {collapse.group(1)} → "
        f"{collapse.group(2)} while its bars draw {expected[:2]}")
    ratio = by[("RSL3", "uniform")] / by[("RSL3", "gradient_120um")]
    assert abs(int(collapse.group(3)) - round(ratio)) <= 1, (
        f"fig24 says ~{collapse.group(3)}x collapse and the fixture gives "
        f"{ratio:.1f}x")


def test_fig26_annotates_no_values_and_so_is_not_gated_here():
    """Stated rather than left as a silent gap.

    fig26 plots curves and labels axes; it prints no data value, so there is
    nothing for this file to compare against `vulnerability_window.json`. Its
    caption numbers are guarded by `test_quantitative_figure_data.py`. If it
    ever gains an annotation, this fails and the annotation gets a guard.
    """
    text = _drawn("fig26_vulnerability_window")
    axis_ticks = {"0", "0.25", "0.5", "1", "2", "3", "7", "14", "28",
                  "10", "20", "30", "40", "50", "60", "80", "100",
                  "0.50", "0.55", "0.60", "0.65", "0.70", "0.75"}
    values = [n for n in re.findall(r"\d+\.\d%|\d+\.\d+×", text)]
    assert not values, (
        f"fig26 now draws data values {values}; add them to this file's "
        "comparison against vulnerability_window.json rather than leaving "
        "them ungated")
    assert "closes ~day 3" in text, (
        "fig26's window annotation is gone; it is the one qualitative claim "
        "the figure draws")
