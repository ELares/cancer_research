"""Guards for the ablation superposition measurement (#844).

This page qualifies a REGISTERED prediction (P20), so the guards protect the
qualification from softening in either direction: the control must stay at the
wide end where the analytic model is genuinely right, the divergence must stay
large enough to matter, and the saturated rows must stay excluded rather than
averaged into a headline that would then read as agreement.
"""
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "scripts" / "validate_ablation_superposition.py"
MD = REPO / "analysis" / "calibration" / "ablation-superposition-validation.md"
JSON_ = REPO / "analysis" / "calibration" / "ablation-superposition-validation.json"
SWEEP = REPO / "analysis" / "calibration" / "ablation_superposition_sweep.txt"
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


def test_the_control_is_at_the_wide_end_where_the_analytic_model_is_right():
    """Agreement has to happen where the physics says it should.

    At wide spacing the second-nearest vessel contributes a factor of ~1, so
    the product reduces to the nearest term and the analytic model IS correct.
    Agreement anywhere else would be a coincidence rather than a control.
    """
    d = _d()
    assert d["control_holds"], "no spacing reproduces the analytic model"
    assert d["control_is_the_wide_end"], (
        "the models agree at a TIGHTER spacing than where they diverge, which "
        "is backwards: the control must sit where the extra cooling factors "
        "are negligible")
    assert len(d["control_points"]) >= 2, (
        "only one spacing agrees, which is too thin to call a control")


def test_the_divergence_is_large_and_in_the_right_direction():
    """The one-vessel model must UNDER-state the survivors. Over-stating would
    mean cooling from extra vessels was somehow making tissue hotter."""
    d = _d()
    assert d["max_understatement"] > 2.0, (
        f"the analytic model under-states by only {d['max_understatement']:.2f}x "
        "at worst, which does not qualify P20 by anything worth saying")
    for p in d["points"]:
        if p["total_failure"]:
            continue
        assert p["understatement"] >= 1.0 - d["control_tolerance"], (
            f"at {p['inter_vessel_um']} um the all-vessel model leaves FEWER "
            "survivors than the nearest-vessel one, which would mean extra "
            "cooling made the tissue hotter")


def test_the_saturated_rows_are_excluded_not_averaged():
    """A ratio across an arm that left the tumour entirely intact is not a
    comparison, and at the tight end it drifts back toward 1 for exactly that
    reason. Reading it as renewed agreement would invert the finding."""
    d = _d()
    assert d["total_failure_points"], (
        "no row exercises the total-failure marker, so the exclusion rule is "
        "untested")
    for p in d["points"]:
        if p["total_failure"]:
            assert p["survivors_all_vessels"] >= p["total_tumor"] - 1
            assert p["kills"] == 0
    # The headline must come from a live row.
    live = [p for p in d["points"] if not p["total_failure"]]
    assert d["max_understatement"] in [p["understatement"] for p in live]
    md = MD.read_text()
    assert "not a comparison" in md and "total failure" in md


def test_it_qualifies_p20_rather_than_refuting_it():
    md = MD.read_text()
    assert "qualifies P20 rather than refuting it" in md
    assert "lower bound" in md
    prereg = (REPO / "PREREGISTRATION.md").read_text()
    assert "**P20." in prereg, "P20 is gone; this page qualifies nothing"
    # P20's direction must still be the registered one, or the qualification
    # is attached to a claim that no longer exists.
    assert "electroporation" in prereg


def test_the_scale_problem_is_recorded_where_it_bit():
    """The default 20 um grid makes the tumour SMALLER than the cooling length,
    and the first sweep reported zero kills everywhere because of it."""
    src = BIN.read_text()
    assert "smaller than the cooling length" in src
    assert "mm_per_cell" in src
    md = MD.read_text()
    assert "zero kills at every vessel density" in md
    assert _d()["meta"]["mm_per_cell"] > 0.05, (
        "the sweep is back on a spheroid-scale cell size, where this arm "
        "cannot ablate anything")


def test_the_sweep_is_the_binarys_own_output():
    txt = SWEEP.read_text()
    for l in txt.splitlines():
        assert l.startswith("ABLATION_SUPER"), f"stray line: {l[:60]!r}"
    src = BIN.read_text()
    assert '"--ablation-superposition-sweep"' in src
    assert "fn run_ablation_superposition_sweep" in src
    # Both models must be computed on the same cells, or they are not a pair.
    assert "survivors_nearest_only" in src and "survivors_all_vessels" in src


def test_the_page_states_what_it_does_not_establish():
    md = MD.read_text()
    for phrase in ("UNCONSTRAINED", "cooling length is still a placeholder",
                   "modelling choice", "Vessels are points"):
        assert phrase in md, f"the page no longer states its limit: {phrase!r}"
