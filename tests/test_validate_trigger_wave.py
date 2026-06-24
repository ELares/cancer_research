"""Tests for the ferroptotic trigger-wave validation (#482).

Pure Python stdlib (no compiled extension), so these run in CI. They cover the
Nagumo front-speed formula, the iron square-root scaling, the numeric/analytic
self-consistency, the drift-guard against the Rust trigger_wave.rs baseline(),
and the overall validation outcome against the measured Co 2024 speeds.
"""

import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import validate_trigger_wave as vt  # noqa: E402


def test_analytical_nagumo_speed_formula():
    # c = sqrt(D*k/2)*(1-2a). D=30, k=8.13, a=0.25 -> sqrt(30*8.13/2)*0.5.
    expected = math.sqrt(30.0 * 8.13 / 2.0) * (1.0 - 2.0 * 0.25)
    assert abs(vt.analytical_front_speed(30.0, 8.13, 0.25) - expected) < 1e-9
    # At a = 0.5 the front stalls (speed exactly 0).
    assert vt.analytical_front_speed(30.0, 8.13, 0.5) == 0.0
    # Above 0.5 the healthy state re-invades (negative speed).
    assert vt.analytical_front_speed(30.0, 8.13, 0.6) < 0.0


def test_speed_scales_as_sqrt_iron():
    """c ~ sqrt(iron_level): doubling iron multiplies the speed by sqrt(2)."""
    base = vt.model_speed(1.0)
    double = vt.model_speed(2.0)
    assert abs(double / base - math.sqrt(2.0)) < 1e-9
    # The measured loaded/baseline ratio (9.40/5.52) corresponds to iron ~2.9.
    loaded = vt.model_speed(2.9)
    assert abs(loaded / base - math.sqrt(2.9)) < 1e-9


def test_baseline_calibrated_to_measured_speed():
    """The baseline model speed lands on the measured 5.52 um/min."""
    assert abs(vt.model_speed(1.0) - 5.52) < 0.6


def test_iron_dose_matches_measured_speeds():
    """At the implied iron fold-changes, the model reproduces DFO 2.33 and
    iron-loaded 9.40 um/min, and the ordering DFO < control < loaded."""
    dfo = vt.model_speed(0.18)
    control = vt.model_speed(1.0)
    loaded = vt.model_speed(2.9)
    assert dfo < control < loaded
    assert abs(dfo - 2.33) < 0.6
    assert abs(loaded - 9.40) < 1.0


def test_gpx4_defense_slows_the_front():
    assert vt.model_speed(1.0, gpx4_defense=0.15) < vt.model_speed(1.0)


def test_numeric_solver_matches_analytical():
    """The pure-Python explicit-Euler solve agrees with the closed form (the
    cross-language self-consistency check)."""
    numeric = vt.numeric_front_speed(1.0)
    analytic = vt.model_speed(1.0)
    assert numeric > 0
    assert abs(numeric - analytic) / analytic < 0.06


def test_drift_guard_matches_rust_baseline():
    """The Python BASELINE constants must equal the Rust trigger_wave.rs
    baseline() (D, base_reaction_rate, ignition_threshold)."""
    checks = vt.drift_guard()
    assert checks["diffusion_um2_per_min"] == 30.0
    assert checks["base_reaction_rate"] == 8.13
    assert checks["ignition_threshold"] == 0.25


def test_full_validation_passes():
    result = vt.validate()
    assert result["all_passed"], result["checks"]
    assert all(result["checks"].values())


def test_measured_data_cites_the_paper():
    rows = vt.load_measured()
    conds = {r["condition"] for r in rows}
    assert conds == {"baseline", "iron_chelation_DFO", "iron_loaded"}
    assert all("38987590" in r["source"] for r in rows)
