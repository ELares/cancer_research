"""Tests for the #330 kill-switch calibration (fit leg).

The model-dependent fit (sim_batch grid search) needs the compiled
`ferroptosis_core` extension and is NOT re-run in CI, so these tests cover the
pure helpers (dose->intensity map, logistic, error metrics, empirical median) and
the internal consistency of the committed calibration result.
"""

import json
import math
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import calibrate_kill_switch as ck  # noqa: E402

CALIB_JSON = REPO_ROOT / "analysis" / "calibration" / "kill-switch-calibration.json"


# --------------------------------------------------------------------------
# Pure helpers
# --------------------------------------------------------------------------


def test_dose_to_inhib_saturating_map():
    assert ck.dose_to_inhib(0.0, 0.5) == 0.0
    assert ck.dose_to_inhib(0.5, 0.5) == pytest.approx(0.5)  # dose == K -> half
    assert ck.dose_to_inhib(1e6, 0.5) == pytest.approx(1.0, abs=1e-5)  # saturates to 1
    # strictly increasing
    xs = [ck.dose_to_inhib(d, 0.5) for d in (0.01, 0.1, 1.0, 10.0)]
    assert all(b > a for a, b in zip(xs, xs[1:]))
    assert all(0.0 <= x < 1.0 for x in xs)


def test_ctrp_viability_decreasing_and_overflow_safe():
    v_low = ck.ctrp_viability(1e-3, 0.05, 1.0, 0.5, -6.0)
    v_high = ck.ctrp_viability(33.0, 0.05, 1.0, 0.5, -6.0)
    assert v_low > v_high
    assert v_low == pytest.approx(1.0, abs=0.02)
    assert v_high == pytest.approx(0.05, abs=0.02)
    # non-responder (huge EC50) must not overflow
    assert 0.5 <= ck.ctrp_viability(10.0, 0.5, 1.01, 4.75e8, -0.45) <= 1.01


def test_error_metrics():
    a = [1.0, 0.5, 0.0]
    b = [1.0, 0.5, 0.0]
    assert ck.sse(a, b) == 0.0
    assert ck.rmse(a, b) == 0.0
    assert ck.sse([1.0, 0.0], [0.0, 0.0]) == 1.0
    assert ck.rmse([1.0, 0.0], [0.0, 0.0]) == pytest.approx(math.sqrt(0.5))


def test_empirical_median_viability_synthetic():
    # Two cell lines with EC50 0.3 and 3.0; at dose 1.0 median is between them.
    rows = [
        {"LowerAsymptote": "0.0", "UpperAsymptote": "1.0", "EC50": "0.3", "Slope": "-4"},
        {"LowerAsymptote": "0.0", "UpperAsymptote": "1.0", "EC50": "3.0", "Slope": "-4"},
    ]
    med = ck.empirical_median_viability(rows, doses=(0.01, 1.0, 100.0))
    assert med[0] == pytest.approx(1.0, abs=0.05)   # both ~alive at low dose
    assert med[2] == pytest.approx(0.0, abs=0.05)   # both ~dead at high dose
    assert 0.0 < med[1] < 1.0                        # mid dose between


# --------------------------------------------------------------------------
# Committed calibration result: structure + internal consistency
# --------------------------------------------------------------------------


def test_committed_calibration_structure():
    r = json.loads(CALIB_JSON.read_text())
    assert r["fit_compound"] == "ML162" and r["heldout_compound"] == "ML210"
    p = r["calibrated_params"]
    assert set(p) == {"lp_propagation", "lp_rate", "k_um"}
    assert 0.0 < p["lp_propagation"] <= 1.0
    assert 0.0 < p["lp_rate"] <= 1.0
    assert p["k_um"] > 0.0


def test_committed_calibration_rmse_self_consistent():
    r = json.loads(CALIB_JSON.read_text())
    c = r["curves"]
    # reported RMSEs recompute from the stored (rounded) curve arrays
    assert ck.rmse(c["model_fit"], c["empirical_fit"]) == pytest.approx(r["fit_rmse"], abs=3e-3)
    assert ck.rmse(c["model_heldout"], c["empirical_heldout"]) == pytest.approx(r["heldout_rmse"], abs=3e-3)
    assert ck.rmse(c["default_uncalibrated_model"], c["empirical_fit"]) == pytest.approx(
        r["default_uncalibrated_rmse"], abs=3e-3
    )


def test_calibration_beats_default_and_generalizes():
    r = json.loads(CALIB_JSON.read_text())
    # the calibrated fit is dramatically better than the uncalibrated default
    assert r["fit_rmse"] < 0.15
    assert r["default_uncalibrated_rmse"] > 0.4
    assert r["default_uncalibrated_rmse"] > 3 * r["fit_rmse"]
    # held-out (same GPX4i mechanism) generalizes: comparable to the fit error
    assert r["heldout_rmse"] < 3 * r["fit_rmse"]
    # cross-mechanism (erastin) is worse than the same-mechanism held-out (as expected)
    assert r["cross_mechanism_rmse"] > r["heldout_rmse"]
