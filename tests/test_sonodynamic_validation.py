"""Guards for the sonodynamic frequency-optimum validation.

The load-bearing risk here is not that the model is wrong -- one of its three
claims IS wrong and the page says so. It is that the refutation quietly turns
into a confirmation on some later regeneration, because a refuted claim is the
uncomfortable one and the document has to keep carrying it.

So these guards pin the REFUTATION as hard as the confirmations, and pin the
things that would make either hollow: that the page is regenerated rather than
committed stale, that the closed form is checked against a scan rather than
against itself, and that the two mismatches the page refuses to reason across
are still stated.
"""
import json
import math
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "scripts" / "validate_sonodynamic.py"
MD = REPO / "analysis" / "calibration" / "sonodynamic-validation.md"
JSON = REPO / "analysis" / "calibration" / "sonodynamic-validation.json"
RUST = REPO / "simulations" / "ferroptosis-core" / "src" / "sonodynamic.rs"
PARAMS = REPO / "simulations" / "ferroptosis-core" / "src" / "params.rs"


def _d():
    return json.loads(JSON.read_text())


def test_the_committed_page_is_what_the_script_produces_now():
    before_md, before_json = MD.read_text(), JSON.read_text()
    try:
        r = subprocess.run([sys.executable, str(SCRIPT)], cwd=REPO,
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-2000:]
        assert MD.read_text() == before_md, (
            "sonodynamic-validation.md is stale; run scripts/validate_sonodynamic.py")
        assert JSON.read_text() == before_json, "the JSON is stale"
    finally:
        MD.write_text(before_md)
        JSON.write_text(before_json)


def test_the_python_closed_form_is_the_rust_one():
    """A drift guard, parsed from source rather than restated.

    Two implementations of one formula is the arrangement this repository
    uses to make a validator independent; it only works while they agree.
    """
    src = RUST.read_text()
    assert "10.0 / denom" in src, "the Rust closed form has changed shape"
    assert "alpha_db_cm_mhz * std::f64::consts::LN_10 * depth_cm" in src, (
        "the Rust denominator has changed; the Python re-implementation in "
        "validate_sonodynamic.py must move with it")
    # And the attenuation constant the page quotes comes from params.rs, not
    # from a literal in the validator.
    alpha = float(re.search(r"sdt_alpha:\s*([0-9.]+)", PARAMS.read_text()).group(1))
    assert _d()["alpha_db_cm_mhz_from_params_rs"] == alpha


def test_the_closed_form_is_checked_against_a_scan_and_not_against_itself():
    """The check that makes claim 1 a derivation.

    A guard whose expected value comes from the thing it is guarding proves
    nothing; here the expectation is an independent brute scan.
    """
    for row in _d()["interior"]:
        assert row["agrees"], f"closed form and scan disagree at {row['depth_cm']} cm"
        assert row["beats_low_end"] and row["beats_high_end"], (
            f"the optimum at {row['depth_cm']} cm is not interior -- a "
            "monotonic model would satisfy the closed form and still be wrong")
        assert row["scanned_mhz"] > 0


def test_the_depth_scaling_is_still_reported_as_refuted():
    """The uncomfortable claim, pinned.

    If the comparator's own optimum ever stops being flat, or the model's
    stops sliding, this must be re-derived rather than left saying REFUTED --
    and equally, a regeneration must not quietly upgrade it.
    """
    d = _d()
    assert d["claim_depth_scaling"] == "REFUTED"
    assert d["model_depth_span_ratio"] > 2.0, (
        "the model no longer predicts a large depth swing, so the refutation "
        "no longer describes what the model does")
    assert all(v < 1.6 for v in d["reported_depth_span_ratios"].values()), (
        "the comparator's optimum is no longer flat with depth")
    md = MD.read_text()
    assert "no value of" in md and "flat" in md, (
        "the page no longer states that this is a disagreement rather than a "
        "tuning gap")
    assert "near-field" in md.lower(), (
        "the page no longer names the missing term, which is the only useful "
        "thing a refutation leaves behind")


def test_the_two_confirmed_claims_are_confirmed_for_a_reason():
    d = _d()
    assert d["claim_interior_optimum"] == "CONFIRMED"
    assert d["claim_falls_with_attenuation"] == "CONFIRMED"
    # Not just the model's own direction: the comparator must independently
    # show it too, or claim 2 is the model agreeing with itself.
    assert d["reported_optimum_falls_with_attenuation"], (
        "the comparator no longer shows the attenuation direction, so claim 2 "
        "rests on the model alone")
    for row in d["attenuation"]:
        assert row["falls"]
        assert abs(row["ratio"] - 2.0) < 1e-6, (
            "doubling alpha no longer halves the optimum exactly; that exact "
            "factor is the sharper half of the claim")


def test_the_page_refuses_the_numerical_comparison_and_says_why():
    md = MD.read_text()
    for phrase in ("No numerical agreement", "opposite sign", "Np/m/MHz",
                   "factor of two"):
        assert phrase in md, (
            f"the page no longer states the mismatch {phrase!r}; a numeric "
            "match across either of them would be uninterpretable")
    assert d_verdict() == "PARTIAL"


def d_verdict():
    return _d()["verdict"]


def test_the_comparator_rows_are_the_published_ones():
    """The anchor cannot drift into being convenient.

    These six values are what PMID 26233216 reports; if any of them is edited
    the claim changes and the edit must be deliberate.
    """
    rows = {(e["alpha_np_m_mhz"], e["depth_mm"]): e["reported_khz"]
            for e in _d()["ellens"]}
    assert rows == {
        (5.0, 50): 750, (5.0, 100): 750, (5.0, 150): 750,
        (10.0, 50): 750, (10.0, 100): 500, (10.0, 150): 500,
    }


def test_the_model_answers_inside_the_band_the_comparator_scanned():
    """Weaker than a fit and stated as such: a model whose optimum fell
    outside 250-1500 kHz everywhere would not be answering the same question,
    and 'the directions agree' would mean nothing."""
    band = _d()["band"]
    assert band and all(r["in_scanned_band"] for r in band)
    # Straddling: the model must be above the comparator at the shallow end
    # and below it at the deep end, which is what a sliding optimum against a
    # flat one MEANS. Both on one side would be an offset, not a scaling gap.
    khz = {r["depth_mm"]: r["model_khz"] for r in band}
    assert khz[50] > 750 > khz[150], (
        f"the model no longer straddles the comparator's 750 kHz: {khz}")


def test_the_engine_module_carries_the_threshold_distinction():
    """SDT's limit is a threshold, which is what makes it structurally unlike
    the dose-response arms. If `cavitates` ever becomes a gradient the whole
    section's argument changes."""
    src = RUST.read_text()
    assert "pub fn cavitates(" in src
    assert re.search(r"index >= threshold", src), (
        "cavitation is no longer a threshold comparison")
    assert "MI_DIAGNOSTIC_CAP" in src
    assert "REGULATORY LIMIT, NOT A MEASUREMENT" in src, (
        "the diagnostic cap no longer says it is a regulatory ceiling rather "
        "than the pressure at which tissue cavitates")


def test_the_closed_form_reproduces_by_hand():
    """One arithmetic check the guards do not take from the artifact."""
    # f* = 10 / (alpha * ln10 * z); at alpha=1, z=1 that is 10/ln(10).
    from importlib import util
    spec = util.spec_from_file_location("vs", SCRIPT)
    m = util.module_from_spec(spec)
    spec.loader.exec_module(m)
    assert abs(m.optimal_frequency_mhz(1.0, 1.0) - 10.0 / math.log(10.0)) < 1e-12
    assert math.isinf(m.optimal_frequency_mhz(0.0, 1.0))
