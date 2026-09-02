"""Guards for the oncolytic arm's conditional finding.

The claim is a CONDITION rather than a result: an interior optimum in immune
competence exists above a crossover in priming efficiency and not below it. A
page in that position fails in one specific way -- it stops reporting the
regime where the optimum is at the boundary, and the conditional becomes an
assertion that immunity always helps.
"""
import json
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "scripts" / "validate_oncolytic.py"
MD = REPO / "analysis" / "calibration" / "oncolytic-validation.md"
JSON = REPO / "analysis" / "calibration" / "oncolytic-validation.json"
RUST = REPO / "simulations" / "ferroptosis-core" / "src" / "oncolytic.rs"


def _d():
    return json.loads(JSON.read_text())


def test_the_committed_result_is_what_the_script_produces_now():
    before_md, before_json = MD.read_text(), JSON.read_text()
    try:
        r = subprocess.run([sys.executable, str(SCRIPT)], cwd=REPO,
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-2000:]
        assert MD.read_text() == before_md, "stale; run scripts/validate_oncolytic.py"
        assert JSON.read_text() == before_json, "stale JSON"
    finally:
        MD.write_text(before_md)
        JSON.write_text(before_json)


def test_both_regimes_are_present_or_the_finding_is_not_conditional():
    """THE GUARD THAT MATTERS. If every efficiency gives an interior optimum,
    the model is asserting that immunity helps rather than deriving when."""
    d = _d()
    # SATURATED rows are excluded, and that exclusion is itself the finding
    # of a bug: once the outcome reaches 1 the scan returns the lowest
    # competence that gets there, which reads as a boundary optimum and means
    # the trade-off has stopped operating. Mixed in, the two regimes
    # interleave and the crossover looks unstable.
    rows = [r for r in d["optimum_by_efficiency"] if not r["saturated"]]
    assert d["n_saturated"] > 0, (
        "no row saturates any more, so the exclusion this guard relies on has "
        "stopped being necessary and the page's caveat should be re-derived")
    boundary = [r for r in rows if not r["interior"]]
    interior = [r for r in rows if r["interior"]]
    assert boundary, (
        "no priming efficiency puts the optimum at full suppression; the "
        "finding has stopped being conditional and is now an assertion")
    assert interior, "no priming efficiency gives an interior optimum"
    assert d["verdict"]["interior_optimum_is_conditional"] == "YES"
    # And the two regimes must be ORDERED -- low efficiency at the boundary,
    # high efficiency interior -- or the crossover is not a crossover.
    assert max(r["priming_efficiency"] for r in boundary) \
        < min(r["priming_efficiency"] for r in interior) + 1e-9, (
        "the regimes interleave, so there is no single crossover")


def test_the_crossover_is_recomputed_not_stored():
    d = _d()
    rows = d["optimum_by_efficiency"]
    first_interior = min(r["priming_efficiency"] for r in rows if r["interior"])
    assert abs(first_interior - d["crossover_priming_efficiency"]) < 1e-9, (
        "the reported crossover is not the first interior row")


def test_the_interior_definition_is_stated_because_it_moves_the_answer():
    d = _d()
    assert 0.0 < d["interior_tolerance_fraction"] < 0.5
    md = MD.read_text()
    assert "`interior` is a definition" in md
    assert f"{d['interior_tolerance_fraction']:.0%}" in md


def test_the_trial_is_direction_only_and_the_page_says_why():
    d = _d()
    assert d["trial"]["pmid"] == "26014293"
    assert d["verdict"]["trial_anchor"].startswith("DIRECTION")
    md = MD.read_text()
    assert "cannot anchor this model" in md
    assert "two different agents" in md, (
        "the page no longer says WHY the checkpoint arm's ratio trick does not "
        "transfer, which is the whole reason this anchor is weaker")


def test_the_two_implementations_agree_on_the_crossover():
    """The Python re-implements the crate's stepper rather than importing it,
    so agreement is evidence."""
    src = RUST.read_text()
    assert "priming_efficiency_for_interior_optimum" in src
    assert "the_interior_optimum_is_a_condition_and_not_a_built_in" in src, (
        "the crate no longer asserts the conditional finding")
    d = _d()
    assert 0.5 <= d["crossover_priming_efficiency"] <= 10.0


def test_the_page_states_what_it_cannot_establish():
    md = MD.read_text().lower()
    for phrase in ("not a measurement", "saturating regime",
                   "nothing here is fitted"):
        assert phrase in md, f"the page no longer states: {phrase!r}"
