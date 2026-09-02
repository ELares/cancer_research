"""Guards for the oncolytic percolation threshold (#844).

This page makes the strongest external claim in the chapter -- an emergent
threshold reproducing a lattice constant from computational physics -- so the
guards are aimed at the ways such a claim goes hollow: a bracket so wide it
would contain any constant, a classification threshold picked to make the
answer come out, a comparator quietly edited, and the closed-form contrast
disappearing so the page's whole reason for existing goes with it.
"""
import json
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "scripts" / "validate_oncolytic_percolation.py"
MD = REPO / "analysis" / "calibration" / "oncolytic-percolation-validation.md"
JSON_ = REPO / "analysis" / "calibration" / "oncolytic-percolation-validation.json"
SWEEP = REPO / "analysis" / "calibration" / "oncolytic_percolation_sweep.txt"
BIN = REPO / "simulations" / "sim-tme-3d" / "src" / "main.rs"


def _d():
    return json.loads(JSON_.read_text())


def test_the_committed_page_is_what_the_script_produces_now():
    before_md, before_json = MD.read_text(), JSON_.read_text()
    try:
        r = subprocess.run([sys.executable, str(SCRIPT)], cwd=REPO,
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-2000:]
        assert MD.read_text() == before_md, "the page is stale"
        assert JSON_.read_text() == before_json, "the JSON is stale"
    finally:
        MD.write_text(before_md)
        JSON_.write_text(before_json)


def test_the_bracket_is_tight_enough_to_exclude_the_other_constant():
    """A bracket wide enough to contain both thresholds would confirm nothing.

    The whole claim is that the front behaves like 26-neighbour percolation
    SPECIFICALLY, and that is only a claim while the measurement can tell the
    two neighbourhoods apart.
    """
    d = _d()
    assert d["bracket_contains_moore"]
    assert not d["bracket_contains_von_neumann"], (
        "the measured bracket contains BOTH lattice constants, so it does not "
        "discriminate between neighbourhoods and confirms nothing")
    lo, hi = d["certain_transmission_bracket"]
    assert lo is not None and hi is not None and hi > lo
    assert hi / lo < 2.0, (
        f"the bracket [{lo}, {hi}] spans more than a factor of two, which is "
        "too loose to call agreement with a constant")


def test_the_crossing_line_sits_in_an_empty_gap():
    """The classification must not depend on where the line is drawn.

    If some run stopped near the threshold, moving it slightly would move the
    answer and the bracket would be a property of the choice.
    """
    d = _d()
    cut = d["cross_threshold_um"]
    seeded = [p for r in d["by_replication"] for p in r["points"] if p["seeded"]]
    non_crossing = [p["front_radius_um"] for p in seeded if not p["crossed"]]
    crossing = [p["front_radius_um"] for p in seeded if p["crossed"]]
    assert non_crossing and crossing
    top, bottom = max(non_crossing), min(crossing)
    # The measured GAP, not a band around the cut. The first version of this
    # guard asked whether any run sat within a factor of two of the line and
    # flagged runs at 538 µm against a 270 µm cut -- which are as far onto the
    # crossing side as the data goes. It was measuring the wrong thing.
    assert bottom / top > 2.0, (
        f"the largest non-crossing front is {top:.0f} µm and the smallest "
        f"crossing one {bottom:.0f} µm; without a clear gap the classification "
        "is a property of where the line was drawn")
    assert top < cut < bottom, (
        f"the {cut:.0f} µm crossing line does not sit inside the measured gap "
        f"[{top:.0f}, {bottom:.0f}]")


def test_the_closed_form_contrast_is_the_reason_the_page_exists():
    d = _d()
    assert d["closed_form_positive_everywhere"], (
        "the closed form no longer predicts spread everywhere, which is the "
        "contrast this page is built on")
    assert d["rows_where_closed_form_is_wrong"] >= 10, (
        "too few conditions where the closed form predicts spread and the "
        "front does not cross; the contrast has gone")
    md = MD.read_text()
    assert "Fisher-KPP" in md and "no term for how much of the tumour" in md


def test_the_threshold_rises_as_replication_falls():
    """Site-BOND percolation, and the model must show it or the framing that
    the top row approaches the pure SITE constant is unearned."""
    d = _d()
    assert d["threshold_rises_as_replication_falls"]
    rows = d["by_replication"]
    assert len(rows) >= 3
    assert rows[0]["replication"] > rows[-1]["replication"]
    assert rows[-1]["lowest_crossing"] > rows[0]["lowest_crossing"], (
        "the lowest replication crosses at the same permissive fraction as the "
        "highest, so transmission is not acting as a bond probability")


def test_the_lattice_constants_are_the_published_ones():
    d = _d()
    assert d["pc_moore_26"] == 0.0976
    assert d["pc_von_neumann_6"] == 0.3116
    src = (REPO / "scripts" / "validate_oncolytic_percolation.py").read_text()
    assert "PC_MOORE_26 = 0.0976" in src and "PC_VON_NEUMANN_6 = 0.3116" in src
    # The neighbourhood claim must be true of the engine, not just asserted.
    grid = (REPO / "simulations" / "ferroptosis-core" / "src" / "grid.rs").read_text()
    assert "26-Moore" in grid and "[(usize, usize, usize); 26]" in grid, (
        "the grid no longer returns a 26-Moore neighbourhood, so the Moore "
        "constant is no longer the right comparator")


def test_a_row_that_never_seeded_is_not_reported_as_a_failure_to_spread():
    """Two different claims that render as the same zero."""
    d = _d()
    unseeded = [(r["replication"], f) for r in d["by_replication"]
                for f in r["unseeded_fractions"]]
    assert unseeded, (
        "no row exercises the never-seeded marker, so the distinction is "
        "untested")
    for r in d["by_replication"]:
        for p in r["points"]:
            if not p["seeded"]:
                assert not p["crossed"]
    assert "*never seeded*" in MD.read_text()


def test_the_sweep_is_the_binarys_own_output():
    txt = SWEEP.read_text()
    for l in txt.splitlines():
        assert l.startswith("ONCOLYTIC_PERC"), f"stray line: {l[:60]!r}"
    assert '"--oncolytic-percolation-sweep"' in BIN.read_text()
    assert "fn run_oncolytic_percolation_sweep" in BIN.read_text()


def test_the_page_states_what_it_does_not_establish():
    md = MD.read_text()
    for phrase in ("Not a claim about real oncolytic viruses",
                   "not evidence that tumours percolate",
                   "contact process, not a PDE"):
        assert phrase in md, f"the page no longer states its limit: {phrase!r}"
