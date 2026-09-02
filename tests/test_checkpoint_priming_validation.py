"""Guards for the checkpoint-blockade priming measurement (#844).

The finding here is a NEGATIVE about the point model's own output -- that its
fold-change is the quantity which cannot distinguish the two cases -- so the
guards protect the comparison that makes it a negative: the fold must stay
agent-independent while the share does not. If either half drifted, the claim
would invert without anything failing.
"""
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "scripts" / "validate_checkpoint_priming.py"
MD = REPO / "analysis" / "calibration" / "checkpoint-priming-validation.md"
JSON_ = REPO / "analysis" / "calibration" / "checkpoint-priming-validation.json"
SWEEP = REPO / "analysis" / "calibration" / "checkpoint_priming_sweep.txt"
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


def test_a_cold_tumour_gets_nothing_at_any_blockade_strength():
    """Checkpoint blockade cannot START a response.

    Not 'gets little' -- gets exactly nothing, at every strength tested. If
    the control ever produced a kill the claim would need rewriting rather
    than the threshold moving.
    """
    d = _d()
    assert d["cold_tumour_gets_nothing"]
    cold = next(a for a in d["arms"] if a["treatment"] == "Control")
    assert len(cold["points"]) >= 3, "too few blockade strengths to say 'any'"
    assert all(p["immune_kills"] == 0 for p in cold["points"])
    assert cold["ferroptosis_kills"] < 10, (
        "the control is no longer a cold tumour -- it is killing cells, so it "
        "does not test what this row claims")


def test_the_fold_is_agent_independent_and_the_share_is_not():
    """The whole negative, in one guard.

    A multiplier does not know what it multiplies, so the fold must be the
    same for both active treatments. The share must NOT be, or the point
    model's output would have been sufficient after all.
    """
    d = _d()
    assert d["fold_is_agent_independent"], (
        f"the fold benefit differs by {d['fold_spread']:.2f}x between "
        "treatments; if the multiplier varied by agent the page's central "
        "claim would be wrong")
    assert d["share_differs_by_an_order_of_magnitude"], (
        f"the immune share differs by only {d['share_ratio']:.1f}x, so the "
        "denominator no longer distinguishes the cases and the finding is gone")
    assert all(f > 1.5 for f in d["fold_benefits"]), (
        "blockade barely moves immune killing at all, so there is no effect "
        "whose denominator is worth arguing about")


def test_more_danger_signal_goes_with_a_smaller_share():
    """The counterintuitive part, and the one a reader would most like to
    forget. Ranking treatments by immunogenic signal and picking the top one
    as a blockade partner gets it backwards here."""
    d = _d()
    assert d["more_damp_smaller_share"]
    assert d["damp_ratio"] > 2.0, (
        "the two active treatments no longer differ much in DAMP, so 'more "
        "signal, smaller share' is not a contrast")
    md = MD.read_text()
    assert "gets this backwards" in md


def test_the_response_is_monotone_in_blockade_strength():
    d = _d()
    assert d["monotone_in_blockade"]
    for a in d["arms"]:
        if a["treatment"] == "Control":
            continue
        brakes = [p["combined_brake"] for p in a["points"]]
        assert all(x >= y for x, y in zip(brakes, brakes[1:])), (
            "the combined brake does not fall as blockade rises, so the panel "
            "is not being driven")


def test_the_sweep_is_the_binarys_own_output_and_needs_the_immune_layer():
    txt = SWEEP.read_text()
    for l in txt.splitlines():
        assert l.startswith("CHECKPOINT_PRIMING"), f"stray line: {l[:60]!r}"
    src = BIN.read_text()
    assert '"--checkpoint-priming-sweep"' in src
    assert "fn run_checkpoint_priming_sweep" in src
    # The arm's defining property, in the code rather than only the prose.
    assert "checkpoint blockade cannot START a response" in src
    assert "immune: Some(SpatialImmuneConfig::for_3d())" in src


def test_the_arm_is_a_modifier_and_the_parity_table_says_so():
    """Selecting `Treatment::Immunotherapy` does nothing; the brake acts
    through the immune layer whatever treatment is chosen. Rounding that up to
    'wired' would overstate the campaign."""
    parity = json.loads((REPO / "analysis" / "arm-parity.json").read_text())
    assert "Immunotherapy" in parity["spatial_as_modifier"]
    assert "Immunotherapy" not in parity["spatial_on_demand"]
    assert "Immunotherapy" not in parity["spatial_by_default"]
    row = next(r for r in parity["arms"] if r["arm"] == "Immunotherapy")
    assert row["spatial"] == "as modifier"
    src = (REPO / "simulations" / "sim-tme-3d" / "src" / "main.rs").read_text()
    assert "belongs in the on-demand tier" in src, (
        "the Rust probe no longer checks that selecting the Treatment does "
        "nothing, so the modifier label rests on nobody's measurement")


def test_the_page_states_what_it_does_not_establish():
    md = MD.read_text()
    for phrase in ("104:1", "counts are small", "placeholders",
                   "INADMISSIBLE", "three points"):
        assert phrase in md, f"the page no longer states its limit: {phrase!r}"
