"""Guards for the ADC drug-loading result.

This is the campaign's strongest calibration: an optimum that falls out of
somebody else's measured clearance ratios rather than being chosen. Two things
would quietly destroy it -- smoothing the curve until it fits one ratio and
misses the other, and letting "delivered payload" drift into being read as
efficacy.
"""
import json
import math
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "scripts" / "validate_adc.py"
MD = REPO / "analysis" / "calibration" / "adc-validation.md"
JSON = REPO / "analysis" / "calibration" / "adc-validation.json"
CSV = REPO / "analysis" / "calibration" / "adc_drug_loading.csv"
RUST = REPO / "simulations" / "ferroptosis-core" / "src" / "adc.rs"


def _d():
    return json.loads(JSON.read_text())


def test_the_committed_result_is_what_the_script_produces_now():
    before_md, before_json = MD.read_text(), JSON.read_text()
    try:
        r = subprocess.run([sys.executable, str(SCRIPT)], cwd=REPO,
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-2000:]
        assert MD.read_text() == before_md, "stale; run scripts/validate_adc.py"
        assert JSON.read_text() == before_json, "stale JSON"
    finally:
        MD.write_text(before_md)
        JSON.write_text(before_json)


def test_the_anchor_points_are_the_measured_ratios():
    """Parsed from the crate, and checked against the committed data file --
    which is the paper's numbers and nobody's arithmetic."""
    d = _d()
    rows = {r["conjugate"]: r for r in d["measured_rows"]}
    assert set(rows) == {"E2", "E4", "E8"}
    for row in rows.values():
        assert row["pmid"] == "15501986"
    e2, e4, e8 = (float(rows[k]["clearance_relative_to_dar2"])
                  for k in ("E2", "E4", "E8"))
    assert abs(e8 / e4 - 3.0) < 0.01, "the E8/E4 clearance ratio is not 3"
    assert abs(e8 / e2 - 5.0) < 0.01, "the E8/E2 clearance ratio is not 5"
    anchors = {p["dar"]: p["clearance_relative_to_dar2"] for p in d["anchor_points"]}
    assert abs(anchors[4.0] - e4) < 0.01 and abs(anchors[8.0] - e8) < 0.01, (
        "the crate's anchor points disagree with the committed measurement")


def test_a_single_power_law_still_misses_one_ratio():
    """The reason the curve is piecewise. If this ever stops being true the
    justification for the shape has changed and the page must be re-derived."""
    d = _d()
    assert d["verdict"]["single_power_law_fits_both_ratios"] == "NO"
    predicted = d["single_power_law_predicts_c8_over_c4"]
    assert abs(predicted - 3.0) > 0.5, (
        f"a single power law now predicts {predicted} against a measured 3, so "
        "the acceleration this result rests on has gone")
    # Recomputed here, independently of the script.
    exponent = math.log(5.0) / math.log(4.0)
    expected = (4.0 ** exponent) / (2.0 ** exponent)
    assert abs(expected - predicted) < 0.01


def test_the_optimum_is_emergent_and_where_the_field_settled():
    d = _d()
    assert abs(d["optimal_dar"] - 4.0) < 0.3, (
        f"the optimum moved to DAR {d['optimal_dar']}")
    assert d["delivered_8_over_4"] < 1.0, (
        "a DAR-8 conjugate now delivers more than a DAR-4 one, which is the "
        "opposite of what the study reports")
    assert d["delivered_8_over_4"] > 0.5
    # It has to be an INTERIOR maximum, or it is not an optimum at all.
    curve = {row["dar"]: row["delivered_per_dose"] for row in d["curve"]}
    peak = d["optimal_dar"]
    assert curve[min(curve)] < curve[peak] and curve[max(curve)] < curve[peak], (
        "the delivered-payload curve is monotonic, so there is no optimum")


def test_the_two_orderings_still_disagree():
    d = _d()
    assert d["verdict"]["in_vitro_and_in_vivo_orderings_agree"] == "NO"
    curve = {row["dar"]: row for row in d["curve"]}
    assert curve[8.0]["in_vitro_potency"] > curve[4.0]["in_vitro_potency"]
    assert curve[8.0]["delivered_per_dose"] < curve[4.0]["delivered_per_dose"]
    assert "opposite ordering" in MD.read_text()


def test_the_page_states_what_it_cannot_establish():
    md = MD.read_text().lower()
    for phrase in ("the interpolation is an assumption",
                   "one conjugate, one payload",
                   "delivered payload is not efficacy"):
        assert phrase in md, f"the page no longer states: {phrase!r}"
