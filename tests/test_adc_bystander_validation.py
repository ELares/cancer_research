"""Guards for the spatial ADC bystander effect (#844).

The prediction here came with a condition the sweep FOUND rather than the
author choosing, and the temptation this page has to resist is reporting the
reach that agreed. So the guards require the disagreeing arm to stay visible,
require the two opposing columns to stay side by side, and require the scalar
the point model would have returned to keep sitting there not moving.
"""
import json
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "scripts" / "validate_adc_bystander.py"
MD = REPO / "analysis" / "calibration" / "adc-bystander-validation.md"
JSON_ = REPO / "analysis" / "calibration" / "adc-bystander-validation.json"
SWEEP = REPO / "analysis" / "calibration" / "adc_bystander_sweep.txt"
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


def test_the_arm_that_disagrees_stays_on_the_page():
    """The condition is the finding, not a caveat on it.

    At a payload reach of one cell the advantage is flat in penetration,
    because the payload lands on cells the conjugate reached anyway. Dropping
    that row would turn a conditional result into an unconditional one.
    """
    d = _d()
    assert d["reaches_where_fold_falls"], "no reach shows the predicted fall"
    assert d["reaches_where_fold_is_flat"], (
        "every reach now shows the fall, so the page's central condition -- "
        "that the payload must travel beyond one cell -- is no longer true and "
        "has to be re-derived rather than quietly dropped")
    assert d["flat_arm_span_ratio"] < 1.25, (
        "the 'flat' arm spans more than 25%, which is not flat and makes the "
        "word wrong")
    assert "PARTIAL" in d["verdict"]
    md = MD.read_text()
    assert "The sweep found that condition; it was not chosen" in md
    assert "picking the arm that agreed" in md


def test_both_opposing_columns_are_reported():
    """The ratio falls and the absolute count rises. Publishing one is
    publishing the more flattering of two true answers."""
    d = _d()
    md = MD.read_text()
    assert "absolute bystander count moves the other way" in md
    assert "more flattering of two" in md
    for r in d["by_reach"]:
        for p in r["points"]:
            assert p["bystander_kills"] >= 0
            assert p["fold_advantage"] is not None
            assert p["with_payload"] >= p["direct_kills"], (
                "a cleavable linker killed FEWER cells than a non-cleavable "
                "one at the same settings, which cannot be right")


def test_the_fall_is_large_enough_to_be_a_claim():
    d = _d()
    falling = [c for c in d["checks"] if c["fold_falls_with_penetration"]]
    assert falling
    for c in falling:
        lo, hi = c["fold_range"]
        assert hi / lo > 1.4, (
            f"at reach {c['payload_reach_cells']} the advantage spans only "
            f"{hi/lo:.2f}x across a sixteenfold penetration range, which is "
            "too little to call a trend")


def test_the_scalar_the_point_model_returns_never_moves():
    """The whole reason the arm needed a 'where'. If this ever varied with
    penetration the point model would already have expressed the effect."""
    d = _d()
    assert d["scalar_never_moves"]
    for r in d["by_reach"]:
        vals = {p["scalar_bystander_fraction"] for p in r["points"]}
        assert len(vals) == 1
    src = (REPO / "simulations" / "ferroptosis-core" / "src" / "adc.rs").read_text()
    assert "neighbours_in_reach" in src
    assert "An INPUT, not derived" in src, (
        "adc.rs no longer records that neighbours_in_reach is an input rather "
        "than a geometric consequence, which is the gap this page fills")


def test_the_bystander_is_one_hop_and_the_page_says_why():
    """A bystander-killed cell must release nothing. The payload came from the
    conjugate, and a cell that never took up conjugate has none to give;
    letting it release would be a chain reaction, not a bystander effect."""
    src = BIN.read_text()
    assert "chain reaction rather than a" in src and "bystander effect" in src
    assert "Payload release is one hop, deliberately" in MD.read_text()


def test_the_sweep_is_the_binarys_own_output():
    txt = SWEEP.read_text()
    for l in txt.splitlines():
        assert l.startswith("ADC_BYSTANDER"), f"stray line: {l[:60]!r}"
    assert '"--adc-bystander-sweep"' in BIN.read_text()
    assert "fn run_adc_bystander_sweep" in BIN.read_text()
    # Both linker arms must be present, or the difference is not a difference.
    assert "cleavable=true" in txt and "cleavable=false" in txt


def test_the_page_states_what_it_does_not_establish():
    md = MD.read_text()
    for phrase in ("No calibration", "NO TARGET", "payload reach is a parameter",
                   "binding-site consumption"):
        assert phrase in md, f"the page no longer states its limit: {phrase!r}"
