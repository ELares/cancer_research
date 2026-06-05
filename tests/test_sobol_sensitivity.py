"""Validate the #331 Sobol estimator and guard the headline sensitivity result.

The estimator math is checked against the analytic Ishigami function (whose
first-order and total-effect Sobol indices are known in closed form); the
model-specific claim (the ferroptosis kill switch is driven by lp_propagation)
is guarded by a small binding-based run.
"""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
SPEC = importlib.util.spec_from_file_location(
    "sobol_sensitivity", REPO / "scripts" / "sobol_sensitivity.py"
)
sob = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(sob)


def test_estimator_matches_ishigami_analytic_indices():
    """Ishigami (a=7, b=0.1) has closed-form Sobol indices:
    S1 ≈ [0.314, 0.442, 0.0], ST ≈ [0.557, 0.442, 0.244]."""
    a, b = 7.0, 0.1

    def ishigami(rows):
        x1, x2, x3 = rows[:, 0], rows[:, 1], rows[:, 2]
        return np.sin(x1) + a * np.sin(x2) ** 2 + b * x3**4 * np.sin(x1)

    lows = [-np.pi, -np.pi, -np.pi]
    highs = [np.pi, np.pi, np.pi]
    s1, st, _, _ = sob.sobol_indices(ishigami, lows, highs, n_base=32768, rng_seed=7)

    s1_true = np.array([0.3139, 0.4424, 0.0])
    st_true = np.array([0.5574, 0.4424, 0.2436])
    assert np.allclose(s1, s1_true, atol=0.04), f"S1 {s1} vs {s1_true}"
    assert np.allclose(st, st_true, atol=0.04), f"ST {st} vs {st_true}"
    # x3 has zero first-order effect but a real total effect (pure interaction).
    assert s1[2] < 0.05 < st[2], f"x3 should be interaction-only: S1={s1[2]}, ST={st[2]}"


def test_estimator_invariants():
    """ST >= S1 >= ~0 for every parameter, on a simple additive model."""

    def additive(rows):
        return 2.0 * rows[:, 0] + rows[:, 1]  # no interactions

    s1, st, _, _ = sob.sobol_indices(
        additive, [0, 0, 0], [1, 1, 1], n_base=4096, rng_seed=3
    )
    assert np.all(st >= s1 - 0.02), f"ST must be >= S1: S1={s1}, ST={st}"
    assert np.all(s1 >= -0.02) and np.all(st >= -0.02)
    # Purely additive ⇒ ST ≈ S1 (no interactions), and x2 (unused) ≈ 0.
    assert abs(st[2]) < 0.03 and abs(s1[2]) < 0.03


def test_lp_propagation_drives_the_ferroptosis_switch():
    """Small binding-based run: the autocatalytic propagation rate dominates the
    kill switch, and the ROS-driven GPX4 degradation is kill-rate-insensitive.
    Skipped where the compiled `ferroptosis_core` extension is not built (the
    Python CI does not build it; the estimator math is covered by the Ishigami
    test, which needs no binding)."""
    pytest.importorskip("ferroptosis_core")
    s1, st, _, _ = sob.saltelli_indices(n_base=256)
    names = [p[0] for p in sob.PARAMS]
    rank = {n: st[i] for i, n in enumerate(names)}
    top = max(rank, key=rank.get)
    assert top == "lp_propagation", f"expected lp_propagation to dominate; got {top}"
    assert rank["lp_propagation"] > rank["gpx4_rate"] > 0.05
    # gpx4_degradation_by_ros barely moves the kill rate ⇒ non-identifiable here.
    assert rank["gpx4_degradation_by_ros"] < 0.05
