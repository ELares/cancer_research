"""Guards for the chemotherapy arm's two structural predictions.

The page these check is unusual: its headline is that the measurement is out of
reach. That makes it exactly the page most likely to drift into claiming
something it cannot support, so the guards are aimed at the three ways it
could:

  * the "no reachable target" statement quietly disappearing once somebody
    forgets why it is there
  * the two-sided dose-density window collapsing to the one-sided story that
    was expected before the model was run
  * the two implementations agreeing because they share code
"""
import json
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "scripts" / "validate_chemo.py"
MD = REPO / "analysis" / "calibration" / "chemo-validation.md"
JSON = REPO / "analysis" / "calibration" / "chemo-validation.json"
RUST = REPO / "simulations" / "ferroptosis-core" / "src" / "chemo.rs"


def _d():
    return json.loads(JSON.read_text())


def test_the_committed_result_is_what_the_script_produces_now():
    before_md, before_json = MD.read_text(), JSON.read_text()
    try:
        r = subprocess.run([sys.executable, str(SCRIPT)], cwd=REPO,
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-2000:]
        assert MD.read_text() == before_md, "the page is stale; run scripts/validate_chemo.py"
        assert JSON.read_text() == before_json, "the JSON is stale"
    finally:
        MD.write_text(before_md)
        JSON.write_text(before_json)


def test_the_two_implementations_pin_each_other():
    """Neither side reads the other, so agreement means something.

    The Rust test pins two points of the dose-response to four decimals; the
    Python artifact carries the same two numbers, derived from a table parsed
    out of the Rust. Deleting the Rust pin fails here, and changing either
    implementation fails there.
    """
    src = RUST.read_text()
    assert "the_dose_response_matches_the_independent_python_mirror" in src, (
        "the crate no longer pins the cross-implementation points")
    for literal in ("0.4531", "0.0275"):
        assert literal in src, (
            f"the crate's cross-implementation test no longer pins {literal}")
    curves = _d()["dose_response"]
    at8 = {cls: dict(pts)[8] for cls, pts in curves.items()}
    assert abs(at8["SPhaseSpecific"] - 0.4531) < 5e-5, at8
    assert abs(at8["PhaseNonspecific"] - 0.0275) < 5e-5, at8


def test_the_sensitivity_table_is_read_from_the_crate():
    """A Python copy of the table would let the two implementations agree
    while both drifted away from the code they describe."""
    table = _d()["sensitivity_table"]
    assert table, "the artifact carries no sensitivity table"
    src = RUST.read_text()
    for key, value in table.items():
        cls, phase = key.split("/")
        assert re.search(
            rf"\({re.escape(cls)}, {re.escape(phase)}\) => {re.escape(str(value))}", src), (
            f"{key} = {value} is not what chemo.rs says")


def test_the_dose_density_window_has_two_ends():
    """The result the model produced and the author did not expect.

    A one-sided story -- "faster is better when the tumour regrows" -- is the
    one everybody starts with, and it is what this guard exists to stop the
    page relaxing back into.
    """
    dd = _d()["dose_density"]
    assert dd["peak_advantage"] > 2.0, (
        f"no dose-density advantage anywhere: peak {dd['peak_advantage']}")
    assert dd["advantage_at_zero_regrowth"] < 1.05, (
        "an advantage appears with nothing to outrun, which would be an "
        f"artifact: {dd['advantage_at_zero_regrowth']}")
    assert dd["advantage_at_fast_regrowth"] < 1.2, (
        "the advantage does not vanish at fast regrowth, so the two-sided "
        f"finding has changed: {dd['advantage_at_fast_regrowth']}")
    assert dd["window_lo"] is not None and dd["window_hi"] > dd["window_lo"]
    assert 0.0 < dd["peak_at_regrowth_per_day"] < 0.5


def test_the_residue_gap_narrows_out_of_cycle():
    r = _d()["residue_ratio_at_dose_8"]
    for cls in ("SPhaseSpecific", "MPhaseSpecific"):
        assert r["proliferating"][cls] > 1.0, cls
        assert r["quiescent_rich"][cls] < r["proliferating"][cls], (
            f"{cls}: the gap did not narrow out of cycle "
            f"({r['quiescent_rich'][cls]} vs {r['proliferating'][cls]})")


def test_the_page_still_says_the_target_is_unreachable():
    """The most important sentence on the page, and the easiest to lose."""
    d = _d()
    assert d["dose_response_target"]["status"] == "UNREACHABLE"
    assert d["verdicts"]["dose_response_magnitude"] == "NO TARGET REACHABLE"
    md = MD.read_text()
    assert "UNREACHABLE" in md
    assert "verification page" in md, (
        "the page no longer says WHY the target is unreachable, which is what "
        "makes it a recoverable block rather than an excuse")
    for phrase in ("placeholder", "Not a reproduction of CALGB 9741",
                   "population-level"):
        assert phrase in md, f"the page no longer states: {phrase!r}"
