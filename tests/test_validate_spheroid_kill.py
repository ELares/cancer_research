"""Tests for the #333 spheroid kill-vs-size analysis (directional validation +
falsifiable prediction). The sweep CSV is produced by the compiled binary
(sim-tme-3d --spheroid-size-sweep) and committed; these pure-stdlib tests run in
CI against that committed CSV.
"""

import json
import sys
from types import SimpleNamespace
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import validate_spheroid_kill as vsk  # noqa: E402

SWEEP_CSV = REPO_ROOT / "analysis" / "calibration" / "spheroid_kill_vs_size.csv"


def test_monotone_helper():
    assert vsk._is_monotone_decreasing([0.036, 0.016, 0.008, 0.003, 0.0019]) is True
    assert vsk._is_monotone_decreasing([0.0, 0.0, 0.0006, 0.0029, 0.0019]) is False


def test_committed_sweep_loads_both_modes():
    rows = vsk.load_sweep(SWEEP_CSV)
    assert len(rows) == 10  # 5 sizes x {fixed, size-aware}
    assert {r["size_aware"] for r in rows} == {True, False}
    # every size has both modes
    radii = sorted({r["radius_um"] for r in rows})
    assert radii == [144.0, 216.0, 288.0, 432.0, 540.0]


def test_fixed_thresholds_reproduce_bigger_resists_more():
    rows = vsk.load_sweep(SWEEP_CSV)
    r = vsk.analyze(rows)
    sg = r["supply_gradient_direction"]
    # RSL3 kill falls monotonically with size (supply-gradient direction)
    assert sg["monotone_decreasing_with_size"] is True
    # fold-drop is substantial and lands in the measured cytotoxic size-resistance range
    assert sg["fold_drop_small_to_large"] > 5.0
    assert sg["fold_in_measured_range"] is True


def test_size_aware_persister_targeting_prediction():
    rows = vsk.load_sweep(SWEEP_CSV)
    r = vsk.analyze(rows)
    pt = r["persister_targeting_prediction"]
    # small all-proliferating spheroids resist the persister-targeting inducer (~0 kill)
    assert pt["small_spheroids_resist"] is True
    # vulnerability emerges as the persister core appears at larger size
    assert pt["vulnerability_emerges_with_core"] is True


def test_run_writes_outputs_and_is_deterministic():
    a = vsk.run(SimpleNamespace(csv=SWEEP_CSV))
    assert vsk.OUT_MD.exists() and vsk.OUT_JSON.exists()
    b = vsk.run(SimpleNamespace(csv=SWEEP_CSV))
    assert json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True)
