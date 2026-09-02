"""Guards for the ablation arm's spatial prediction.

This arm's calibration is UNCONSTRAINED and will stay that way: a threshold
observable cannot identify a threshold parameter. What replaces it is a
prediction about WHERE the failures are, and the two ways that could quietly
become an overclaim are a sleeve size read as though it were measured, and a
degenerate row read as though it were a sleeve.
"""
import json
import math
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "scripts" / "validate_ablation.py"
MD = REPO / "analysis" / "calibration" / "ablation-validation.md"
JSON = REPO / "analysis" / "calibration" / "ablation-validation.json"
RUST = REPO / "simulations" / "ferroptosis-core" / "src" / "ablation.rs"


def _d():
    return json.loads(JSON.read_text())


def test_the_committed_result_is_what_the_script_produces_now():
    before_md, before_json = MD.read_text(), JSON.read_text()
    try:
        r = subprocess.run([sys.executable, str(SCRIPT)], cwd=REPO,
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-2000:]
        assert MD.read_text() == before_md, "stale; run scripts/validate_ablation.py"
        assert JSON.read_text() == before_json, "stale JSON"
    finally:
        MD.write_text(before_md)
        JSON.write_text(before_json)


def test_a_total_failure_is_not_reported_as_a_sleeve():
    """The same class of defect the oncolytic page had to fix: a degenerate
    row whose value looks like the quantity being reported and means something
    else."""
    d = _d()
    assert d["n_total_failure"] > 0, (
        "no row fails everywhere any more, so the distinction this guard "
        "protects has stopped being exercised and the caveat should be "
        "re-derived rather than kept")
    for row in d["by_applicator_temperature"]:
        if row["total_failure"]:
            assert row["thermal_sleeve_mm"] >= 10.0, (
                "a total failure should return the scan limit, not a plausible "
                "sleeve")
    assert "fails everywhere" in MD.read_text()


def test_the_sleeve_shrinks_with_temperature_over_the_live_rows():
    d = _d()
    live = [r for r in d["by_applicator_temperature"] if not r["total_failure"]]
    sleeves = [r["thermal_sleeve_mm"] for r in live]
    assert len(sleeves) >= 4
    assert all(a >= b for a, b in zip(sleeves, sleeves[1:])), (
        f"the sleeve does not shrink monotonically with temperature: {sleeves}")
    assert d["verdict"]["sleeve_shrinks_with_temperature"] == "YES"
    lo, hi = d["verdict"]["sleeve_range_mm"]
    assert 0.0 < lo < hi < 20.0
    # Millimetre scale, which is what makes it a clinical problem rather than a
    # rounding error.
    assert hi < 10.0, f"the sleeve is {hi} mm, which is not a margin question"


def test_the_electroporation_contrast_is_exactly_zero_everywhere():
    d = _d()
    for row in d["by_applicator_temperature"]:
        assert row["electroporation_sleeve_mm"] == 0.0
    assert d["verdict"]["electroporation_contrast"].startswith("YES")
    src = RUST.read_text()
    assert "electroporation_has_no_perivascular_sleeve_and_that_is_the_contrast" in src, (
        "the crate no longer asserts the contrast that this arm's prediction is")


def test_the_size_is_not_claimed_to_be_anchored():
    d = _d()
    assert d["verdict"]["sleeve_size_anchored"] == "NO"
    assert d["anchor"]["kind"].startswith("DIRECTION")
    md = MD.read_text()
    assert "Direction only" in md
    # The sensitivity to the placeholder has to be shown, because it is what
    # makes the size a restatement of an assumption.
    assert len(d["sensitivity_to_cooling_length"]) >= 4
    sizes = [r["sleeve_at_60c_mm"] for r in d["sensitivity_to_cooling_length"]]
    assert max(sizes) > 3 * min(sizes), (
        "the sleeve is no longer strongly dependent on the cooling-length "
        "placeholder, so the page's central caveat has changed")


def test_the_dose_criterion_is_a_dose():
    """A brief exposure at a lethal temperature is not lethal; the crate has a
    guard for that because a mutation replacing the dose with a bare
    temperature survived everything else."""
    src = RUST.read_text()
    assert "the_kill_criterion_is_a_thermal_DOSE_and_not_a_temperature" in src
    md = MD.read_text()
    assert "thermal DOSE and not a temperature" in md
