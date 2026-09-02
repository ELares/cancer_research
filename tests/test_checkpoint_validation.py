"""Guards for the checkpoint arm's ratio-based calibration.

The claim this page makes is unusual and easy to overstate: that a RATIO
between two strata of one trial constrains something an absolute response rate
could not. It is true, it is narrow, and the ways it could quietly become an
overclaim are what these guards watch:

  * the ratio argument depends on the two strata differing in ONE thing. If a
    stratification that also moves the brake were used, the mapping would stop
    cancelling and the comparison would be invalid.
  * "constrained" must not drift into "identified". One equation, two unknowns.
  * the representative-burden choice moves the model's answer by about as much
    as the target band is wide, and the page has to keep saying so.
"""
import json
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "scripts" / "validate_checkpoint.py"
MD = REPO / "analysis" / "calibration" / "checkpoint-validation.md"
JSON = REPO / "analysis" / "calibration" / "checkpoint-validation.json"
CSV = REPO / "analysis" / "calibration" / "checkpoint_strata.csv"
RUST = REPO / "simulations" / "ferroptosis-core" / "src" / "checkpoint.rs"


def _d():
    return json.loads(JSON.read_text())


def test_the_committed_result_is_what_the_script_produces_now():
    before_md, before_json = MD.read_text(), JSON.read_text()
    try:
        r = subprocess.run([sys.executable, str(SCRIPT)], cwd=REPO,
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-2000:]
        assert MD.read_text() == before_md, "stale; run scripts/validate_checkpoint.py"
        assert JSON.read_text() == before_json, "the committed JSON is stale"
    finally:
        MD.write_text(before_md)
        JSON.write_text(before_json)


def test_the_target_is_a_ratio_from_one_trial_and_one_endpoint():
    """The whole argument depends on this. Two response rates from different
    trials would not share a mapping constant, and nothing would cancel."""
    text = CSV.read_text()
    rows = [l for l in text.splitlines() if l and not l.startswith("#")][1:]
    assert len(rows) >= 1
    for row in rows:
        cells = row.split(",")
        assert cells[1].isdigit(), f"no PMID on {row[:50]}"
    d = _d()
    t = d["trial"]
    assert t["pmid"] == "32919526", "the trial behind the ratio has changed"
    # High and low must come from the SAME trial row, which is how the file is
    # shaped -- one row carries both strata precisely so they cannot drift
    # apart into two sources.
    assert t["stratum_high"] and t["stratum_low"]
    assert float(t["orr_high_pct"]) > float(t["orr_low_pct"])


def test_the_measured_ratio_is_recomputed_not_copied():
    d = _d()
    t = d["trial"]
    expected = float(t["orr_high_pct"]) / float(t["orr_low_pct"])
    assert abs(expected - d["measured_ratio"]) < 0.01, (
        f"the page's ratio {d['measured_ratio']} is not "
        f"{t['orr_high_pct']}/{t['orr_low_pct']}")
    lo, hi = d["measured_band"]
    assert lo < d["measured_ratio"] < hi, "the point estimate is outside its own band"


def test_the_model_lands_in_the_band_and_the_page_says_where():
    d = _d()
    lo, hi = d["measured_band"]
    model = d["model_ratio_at_defaults"]
    assert lo <= model <= hi, (
        f"the model's ratio {model} is outside the published band {lo}-{hi}; "
        "that is the check failing, and the prose has to change with it")
    assert model < d["measured_ratio"], (
        "the model no longer under-predicts the point ratio; the page says it "
        "does and must be re-derived")
    assert "below the point estimate" in MD.read_text()


def test_constrained_has_not_drifted_into_identified():
    """One ratio is one equation. The page may say the antigenicity SHAPE is
    constrained; it may not say the arm is fitted."""
    d = _d()
    assert d["verdict"]["brake_identified"].startswith("NO")
    assert d["verdict"]["shape_constrained"] in (
        "CONSTRAINED", "PARTIALLY CONSTRAINED", "UNCONSTRAINED")
    md = MD.read_text()
    assert "is not an identification" in md
    assert "still unidentified" in md or "still unidentified" in json.dumps(d)
    # And the admitted fraction has to be reported, because a verdict word
    # without the number behind it is the thing this repository keeps retracting.
    assert f"{d['admissible_fraction']:.0%}" in md


def test_the_representative_burden_choice_is_reported_as_a_limit():
    """It moves the answer by about as much as the band is wide, and a target
    that moves with an unstated choice is not a target."""
    d = _d()
    spread = d["verdict"]["representative_tmb_spread"]
    assert spread[1] > spread[0], "the sensitivity sweep produced no spread"
    band = d["measured_band"]
    assert (spread[1] - spread[0]) > 0.2 * (band[1] - band[0]), (
        "the sensitivity is now small relative to the band; the page's central "
        "caveat has changed and should be re-derived rather than kept")
    md = MD.read_text()
    assert "comparable in size to the target" in md
    assert len(d["sensitivity_to_representative_tmb"]) >= 4


def test_the_structural_check_is_structural():
    d = _d()
    assert d["b2m_null_response_index"] == 0.0
    src = RUST.read_text()
    assert "a_tumour_that_cannot_present_does_not_respond_at_any_burden" in src, (
        "the crate no longer asserts the B2M-null structural property")
    assert "the_ratio_cancels_everything_the_two_strata_share" in src, (
        "the crate no longer asserts that the ratio cancels the brake, which "
        "is the argument the trial comparison rests on")


def test_the_page_states_what_it_cannot_establish():
    md = MD.read_text().lower()
    for phrase in ("the brake is still unidentified", "not a response rate",
                   "first order only"):
        assert phrase in md, f"the page no longer states: {phrase!r}"
