"""Guards for the emergent oxygen enhancement ratio (#844).

The trap this page invites is specific and this project has fallen into it
before: `radiation::dna_channel_dose_modifying_factor` returns a single number
that IS the published hyperbola restated, and a page reporting *that* against
the published band would be a guard computing its own expectation. So the
guards below check that the measured factor is genuinely a different quantity
-- that it varies with a spatial parameter the formula does not contain -- and
that the shortfall is reported with its structural cause rather than rounded
up to a pass.
"""
import json
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "scripts" / "validate_radiation_oer.py"
MD = REPO / "analysis" / "calibration" / "radiation-oer-validation.md"
JSON = REPO / "analysis" / "calibration" / "radiation-oer-validation.json"
SWEEP = REPO / "analysis" / "calibration" / "radiation_oer_sweep.txt"
BIN = REPO / "simulations" / "sim-tme-3d" / "src" / "main.rs"


def _d():
    return json.loads(JSON.read_text())


def test_the_committed_page_is_what_the_script_produces_now():
    before_md, before_json = MD.read_text(), JSON.read_text()
    try:
        r = subprocess.run([sys.executable, str(SCRIPT)], cwd=REPO,
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-2000:]
        assert MD.read_text() == before_md, "the page is stale; re-run the generator"
        assert JSON.read_text() == before_json, "the JSON is stale"
    finally:
        MD.write_text(before_md)
        JSON.write_text(before_json)


def test_the_measured_factor_is_not_the_formula_restated():
    """The load-bearing check.

    If the spatial factor merely reproduced the single-cell restatement, the
    sweep would be measuring the formula and this page's whole argument would
    be false. It must vary with the O2 gradient, which the formula has no term
    for, and by a margin far larger than numerical noise.
    """
    d = _d()
    vals = [r["dmf_mean"] for r in d["by_lambda"] if r["dmf_mean"] is not None]
    assert len(vals) >= 4, "too few gradients to say the factor varies"
    assert max(vals) / min(vals) > 1.5, (
        f"the factor spans only {max(vals)/min(vals):.2f}x across the swept "
        "gradients; if it were flat it would be the formula restated")
    restated = d["restated_single_cell_dmf"]
    assert not any(abs(v - restated) < 0.05 for v in vals), (
        f"a measured factor coincides with the restated single-cell value "
        f"({restated}), which is what the page claims it is not")


def test_the_optimum_is_interior_and_the_shortfall_keeps_its_reason():
    d = _d()
    assert d["optimum_is_interior"], (
        "the best gradient sits on the edge of the swept range, so 'the model "
        "has a best gradient' is a statement about the scan")
    assert d["arms_below_band"], (
        "no gradient falls below the band any more; the page's central "
        "explanation -- that one lambda sets rim and core together -- has to "
        "be re-derived rather than left standing")
    md = MD.read_text()
    for phrase in ("lower bound", "one λ sets the rim and the core together",
                   "category error"):
        assert phrase in md, f"the page no longer states: {phrase!r}"


def test_the_band_is_the_published_one_and_the_verdict_follows_from_it():
    d = _d()
    assert d["published_oer_band"] == [2.5, 3.0]
    lo, hi = d["published_oer_band"]
    best = d["best_dmf"]
    assert d["best_reaches_band"] == (lo <= best <= hi), (
        "the verdict flag disagrees with the arithmetic behind it")
    assert ("PARTIAL" in d["verdict"]) == (
        d["best_reaches_band"] and bool(d["arms_below_band"]))


def test_every_reported_factor_has_a_populated_core():
    """The empty-zone trap, which the first sweep ran straight into.

    `zone_kill_rates_3d` returns 0.0 for a zone with no cells, so below grid 40
    the hypoxic kill rate is exactly zero at every dose -- total
    radioresistance, from a division by nothing.
    """
    for r in _d()["by_lambda"]:
        assert r["core_n"] > 500, (
            f"lambda {r['lambda_um']} reports a factor over a core of only "
            f"{r['core_n']} cells")
        assert r["rim_n"] > 500
    md = MD.read_text()
    assert "empty" in md.lower() and "division by nothing" in md, (
        "the page no longer records the empty-zone trap it fell into")


def test_the_sweep_is_the_binarys_output_and_not_a_hand_written_table():
    """Every row must parse as the binary's own print, so the artifact cannot
    be edited into agreement."""
    txt = SWEEP.read_text()
    sweep_rows = [l for l in txt.splitlines() if l.startswith("RADIATION_OER_SWEEP")]
    zone_rows = [l for l in txt.splitlines() if l.startswith("RADIATION_OER_ZONES")]
    assert len(sweep_rows) >= 50 and len(zone_rows) >= 4
    for l in txt.splitlines():
        assert l.startswith("RADIATION_OER_"), f"stray line in the sweep: {l[:60]!r}"
    assert "fn run_radiation_oer_sweep" in BIN.read_text()
    assert '"--radiation-oer-sweep"' in BIN.read_text(), (
        "the flag that regenerates this artifact no longer exists")


def test_only_the_dna_channel_is_wired_and_the_page_says_so():
    """The second channel's `ros_per_gy` has no gray-to-ROS conversion in the
    literature, so enabling it would make this measurement a function of an
    unanchored knob. That has to stay visible, or a later reader will take the
    factor for a two-channel result."""
    assert "ros_per_gy: 0.0" in BIN.read_text()
    md = MD.read_text()
    assert "Only the DNA channel is wired" in md
    assert "unanchored knob" in md
