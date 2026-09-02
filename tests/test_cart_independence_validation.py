"""Guards for the CAR-T barrier-independence check (#844).

The whole design rests on one control: with antigen independent of position the
product must be right, and the spatial run must agree with it. Without that,
every disagreement elsewhere is indistinguishable from a bug -- and the first
version of this measurement HAD one, double-counting the infiltration barrier,
which the control caught. So the guards protect the control first.
"""
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "scripts" / "validate_cart_independence.py"
MD = REPO / "analysis" / "calibration" / "cart-independence-validation.md"
JSON_ = REPO / "analysis" / "calibration" / "cart-independence-validation.json"
SWEEP = REPO / "analysis" / "calibration" / "cart_independence_sweep.txt"
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


def test_the_control_holds_and_the_sweep_contains_it():
    """Zero correlation must be IN the sweep and must agree with the product.

    A sweep that never visits zero could report any divergence it liked with
    nothing to check it against.
    """
    d = _d()
    corrs = [p["correlation"] for p in d["points"]]
    assert 0.0 in corrs, "the sweep has no zero-correlation control"
    assert d["control_holds"], (
        f"the spatial run and the product disagree by {d['control_error']:.1%} "
        "where antigen is independent of position -- that is a bug in the "
        "wiring, not a finding about geometry")
    assert d["control_error"] < d["control_tolerance"]
    md = MD.read_text()
    assert "double-counted the same barrier" in md, (
        "the page no longer records the defect its own control caught")


def test_the_sweep_straddles_zero_in_both_directions():
    d = _d()
    corrs = [p["correlation"] for p in d["points"]]
    assert min(corrs) < -0.3 and max(corrs) > 0.3, (
        "the sweep does not reach far enough either side of zero to show both "
        "directions of the effect")
    assert d["product_optimistic_when_antigen_is_deep"]
    assert d["product_pessimistic_when_antigen_is_rim"]


def test_the_divergence_is_large_enough_to_matter():
    d = _d()
    assert d["max_overstatement"] > 1.10, (
        f"the product only over-states by {d['max_overstatement']:.3f}x at "
        "worst, which is within the noise of an 80k-cell Monte Carlo")
    assert d["max_understatement"] < 0.93


def test_saturation_at_the_extreme_is_marked_not_smoothed():
    """The positive arm stops growing once rim antigen is certain. That is
    saturation, and reporting it as a reversal -- or hiding it -- would both be
    wrong."""
    d = _d()
    md = MD.read_text()
    if not d["positive_arm_is_monotone"]:
        assert d["saturating_correlations"], (
            "the positive arm is not monotone and no correlation is named as "
            "saturating")
        assert "saturation, not a reversal" in md
    else:
        assert "stops growing at the positive extreme" not in md


def test_the_both_ways_bucket_exists_and_is_reported():
    """A product of fractions has no slot for a cell that failed both ways,
    and the two failure modes have different remedies."""
    d = _d()
    lo, hi = d["lost_to_both_range"]
    assert lo > 0 and hi >= lo
    for p in d["points"]:
        assert p["lost_to_both"] > 0
        assert p["unreached"] > 0
    assert "no slot for" in MD.read_text()


def test_the_sweep_is_the_binarys_own_output_and_normalises_the_field():
    txt = SWEEP.read_text()
    for l in txt.splitlines():
        assert l.startswith("CART_INDEP"), f"stray line: {l[:60]!r}"
    src = BIN.read_text()
    assert '"--cart-independence-sweep"' in src
    assert "fn run_cart_independence_sweep" in src
    # The normalisation is what makes the two models differ in DISTRIBUTION
    # only. Without it the comparison is not a comparison.
    assert "NORMALISED SO ITS TUMOUR MEAN IS THE POINT MODEL" in src
    assert "double-counts the same barrier" in src
    # And the achieved mean must be reported, since the clamp moves it.
    assert "adoptive_reach_mean" in src
    for p in _d()["points"]:
        assert 0.0 < p["reach_mean"] <= 1.0


def test_the_page_states_what_it_does_not_establish():
    md = MD.read_text()
    for phrase in ("Not a measurement of any real correlation",
                   "placeholder", "radial-depth proxy", "cap, not a barrier"):
        assert phrase in md, f"the page no longer states its limit: {phrase!r}"
