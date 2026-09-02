"""Guards for the CAR-T arm's discrimination, and for its refusal to fit.

The interesting thing about this page is what it declines to do. It has a
famous number available -- 81% remission in ELIANA -- and does not fit to it,
because a remission is not a kill fraction and, unlike the checkpoint arm, no
ratio is available that cancels the mapping. A page in that position drifts in
one direction: toward quoting the number as though it validated something.
"""
import json
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "scripts" / "validate_adoptive.py"
MD = REPO / "analysis" / "calibration" / "adoptive-validation.md"
JSON = REPO / "analysis" / "calibration" / "adoptive-validation.json"
RUST = REPO / "simulations" / "ferroptosis-core" / "src" / "adoptive.rs"


def _d():
    return json.loads(JSON.read_text())


def test_the_committed_result_is_what_the_script_produces_now():
    before_md, before_json = MD.read_text(), JSON.read_text()
    try:
        r = subprocess.run([sys.executable, str(SCRIPT)], cwd=REPO,
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-2000:]
        assert MD.read_text() == before_md, "stale; run scripts/validate_adoptive.py"
        assert JSON.read_text() == before_json, "stale JSON"
    finally:
        MD.write_text(before_md)
        JSON.write_text(before_json)


def test_the_famous_number_is_present_and_not_used_as_a_target():
    """It would be easy, and wrong, to quote 81% as though the model
    reproduced it."""
    d = _d()
    assert d["anchor"]["remission_rate_pct"] == 81
    assert d["verdict"]["anchor_used_as_target"].startswith("NO")
    md = MD.read_text()
    assert "cannot be used to fit this model" in md
    assert "This page therefore fits nothing." in md
    # And the reason the checkpoint trick does not transfer has to be stated,
    # because the obvious next move is to try it.
    assert "does not transfer here" in md.lower()


def test_the_discrimination_is_a_real_separation():
    d = _d()
    v = d["verdict"]
    assert v["discrimination_holds"] == "YES", (
        "escalating the dose no longer separates a delivery-limited failure "
        "from a density-limited one, which is this arm's only prediction")
    assert v["delivery_gain"] > 5.0
    assert v["density_gain"] < 1.1
    # The two cases must share their barriers, or the comparison is between
    # two different tumours rather than two failure modes.
    cases = {c["case"].split(" (")[0]: c for c in d["discrimination"]}
    assert (cases["delivery-limited"]["delivery_fraction"]
            == cases["density-limited"]["delivery_fraction"]), (
        "the two cases differ in delivery as well as density, so the "
        "comparison does not isolate the failure mode")
    assert (cases["delivery-limited"]["density"]
            > cases["density-limited"]["density"])


def test_the_engagement_curve_is_a_threshold():
    d = _d()
    t = d["constants"]["ANTIGEN_DENSITY_THRESHOLD"]
    curve = dict((row[0], row[1]) for row in d["engagement_curve"])
    at = curve.get(int(t))
    assert at is not None and abs(at - 0.5) < 1e-6, (
        f"engagement at the threshold is {at}, not one half")
    low = [v for k, v in curve.items() if k <= t * 0.2]
    high = [v for k, v in curve.items() if k >= t * 5]
    assert all(v < 0.02 for v in low), "the curve does not switch off below"
    assert all(v > 0.98 for v in high), "the curve does not saturate above"


def test_the_constants_come_from_the_crate():
    src = RUST.read_text()
    c = _d()["constants"]
    m = re.search(r"pub const ANTIGEN_DENSITY_THRESHOLD: f64 = ([0-9.]+);", src)
    assert m and abs(float(m.group(1)) - c["ANTIGEN_DENSITY_THRESHOLD"]) < 1e-9
    for field in ("growth_per_day", "contraction_per_day", "max_fold",
                  "memory_fraction"):
        fm = re.search(rf"{field}: ([0-9.]+),", src)
        assert fm and abs(float(fm.group(1)) - c[field]) < 1e-9, field


def test_the_page_states_what_it_cannot_establish():
    md = MD.read_text().lower()
    for phrase in ("nothing is fitted", "prediction about experiments",
                   "graded in reality"):
        assert phrase in md, f"the page no longer states: {phrase!r}"
