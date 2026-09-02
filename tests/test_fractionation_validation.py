"""Guards for the fractionation calibration leg.

The page these check makes a strong claim -- that two numbers out of a trial
protocol predict an α/β somebody else measured -- and a weaker one that matters
more: that the leg which DISAGREES is reported at the same weight as the leg
that agrees.

Both are easy to lose. A generator can drift until every verdict is
REPRODUCES; a renderer can print the passing row and drop the failing one; a
Python copy of a Rust constant can go stale and take the whole check with it.
Each of those is a separate guard below.
"""
import json
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "scripts" / "validate_fractionation.py"
MD = REPO / "analysis" / "calibration" / "fractionation-validation.md"
JSON = REPO / "analysis" / "calibration" / "fractionation-validation.json"
CSV = REPO / "analysis" / "calibration" / "fractionation_trials.csv"
RUST = REPO / "simulations" / "ferroptosis-core" / "src" / "radiation.rs"


def _d():
    return json.loads(JSON.read_text())


def test_the_committed_result_is_what_the_script_produces_now():
    before_md, before_json = MD.read_text(), JSON.read_text()
    try:
        r = subprocess.run([sys.executable, str(SCRIPT)], cwd=REPO,
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-2000:]
        assert MD.read_text() == before_md, (
            "the committed validation page is stale; run "
            "scripts/validate_fractionation.py")
        assert JSON.read_text() == before_json, "the committed JSON is stale"
    finally:
        MD.write_text(before_md)
        JSON.write_text(before_json)


def test_the_constants_come_from_the_rust_and_not_from_a_python_copy():
    """A second source of truth for a number is how a validation goes stale
    without failing. The script parses radiation.rs; this checks it parsed the
    values that are actually there."""
    src = RUST.read_text()
    for name, value in _d()["constants"].items():
        if isinstance(value, list):          # the published band, a tuple
            m = re.search(rf"pub const {name}: \(f64, f64\) = "
                          rf"\(([0-9.]+), ([0-9.]+)\);", src)
            assert m, f"{name} is no longer in radiation.rs"
            assert [float(m.group(1)), float(m.group(2))] == value, (
                f"{name} disagrees with the crate")
        else:
            m = re.search(rf"pub const {name}: f64 = ([0-9.]+);", src)
            assert m, f"{name} is no longer in radiation.rs"
            assert abs(float(m.group(1)) - value) < 1e-12, (
                f"{name} in the artifact is {value}, in the crate "
                f"{m.group(1)} -- the validation is describing code that has "
                "changed under it")


def test_the_isoeffect_inversion_is_arithmetic_anyone_can_repeat():
    """Recomputed here from the committed schedules, independently of both the
    script and the crate. Three implementations of a four-line formula is not
    redundancy: it is the only way a sign error in one of them is visible."""
    d = _d()
    for row in d["isoeffect"]:
        n1, d1 = (float(x) for x in row["arm_a"].split(" x ")[0:1] + [
            row["arm_a"].split(" x ")[1].split(" ")[0]])
        n2, d2 = (float(x) for x in row["arm_b"].split(" x ")[0:1] + [
            row["arm_b"].split(" x ")[1].split(" ")[0]])
        total1, total2 = n1 * d1, n2 * d2
        expected = (total2 * d2 - total1 * d1) / (total1 - total2)
        assert abs(expected - row["implied_alpha_beta_gy"]) < 0.01, (
            f"{row['site']}: recomputing the inversion gives {expected:.3f}, "
            f"the page says {row['implied_alpha_beta_gy']}")


def test_one_leg_disagrees_and_the_page_says_so():
    """THE GUARD THAT MATTERS.

    A calibration page reporting only successes is marketing. The breast leg
    fails -- START-B's shorter arm delivers less EQD2 and was still not
    inferior, so no positive α/β makes the two equivalent -- and that row is
    the more informative of the two. If it ever starts passing, the finding has
    changed and the prose has to be rewritten rather than the guard relaxed.
    """
    d = _d()
    verdicts = {r["site"]: r["verdict"] for r in d["isoeffect"]}
    assert verdicts.get("breast") == "DISAGREES", (
        "the breast leg now reproduces; the page's central observation about "
        "what EQD2 leaves out has changed and must be re-derived")
    assert verdicts.get("prostate") == "REPRODUCES", (
        "the prostate leg stopped reproducing, which is the check failing")
    md = MD.read_text()
    assert "DISAGREES" in md, "the failing verdict is not on the page"
    # And it must be argued, not merely tabulated: the reason is the whole
    # value of the row.
    assert "still not inferior" in md or "still was not inferior" in md, (
        "the page tabulates the disagreement without saying why it happens")


def test_the_repopulation_prediction_reports_its_own_sensitivity():
    d = _d()["repopulation"]
    assert d["verdict"] == "REPRODUCES", (
        "D_prolif has left the published band; that is the prediction failing")
    assert d["gbm_verdict"] == "DISAGREES", (
        "the glioblastoma α now lands inside the head-and-neck band, so the "
        "sensitivity caveat this page carries has stopped being true")
    assert d["d_prolif_head_and_neck_gy_per_day"] < d["d_prolif_with_gbm_alpha_gy_per_day"]


def test_the_page_states_what_it_cannot_establish():
    md = MD.read_text().lower()
    for phrase in ("not a trial analysis", "non-inferiority",
                   "uncalibrated", "not in use"):
        assert phrase in md, f"the page no longer states: {phrase!r}"


def test_the_trial_data_carries_its_own_provenance():
    """Every row is a claim about a published protocol, so every row names the
    paper it came from and the independent source of the α/β it is checked
    against."""
    text = CSV.read_text()
    assert text.startswith("#"), "the data file has no header comment"
    rows = [l for l in text.splitlines() if l and not l.startswith("#")][1:]
    assert rows, "no trials in the committed data"
    for row in rows:
        cells = row.split(",")
        assert cells[2].isdigit() and len(cells[2]) >= 7, (
            f"row has no PMID: {row[:60]}")
        assert cells[-1].isdigit(), (
            f"row does not name the source of its published α/β: {row[:60]}")
