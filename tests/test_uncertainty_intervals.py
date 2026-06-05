"""Validate the #332 prior-predictive interval estimator and guard the harness.

The estimator math is checked against an analytic case (identity model over a
uniform prior, whose quantiles are known in closed form); the model-specific path
(real death-rate intervals are ordered and in [0, 1]) is guarded by a small
binding-based run that is skipped where the compiled `ferroptosis_core` extension
is not built (the Python CI does not build it).
"""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
SPEC = importlib.util.spec_from_file_location(
    "uncertainty_intervals", REPO / "scripts" / "uncertainty_intervals.py"
)
ui = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ui)


def test_estimator_matches_uniform_analytic_quantiles():
    """Identity model over a single Uniform(2, 10) prior: the output IS the
    parameter, so its 2.5/50/97.5 percentiles are the uniform's closed-form
    quantiles lo + q*(hi-lo), mean = (lo+hi)/2, std = (hi-lo)/sqrt(12)."""
    lo, hi = 2.0, 10.0

    def identity(draws):
        return draws[:, 0]  # (n,) single output

    quantiles = [0.025, 0.5, 0.975]
    stats, out = ui.prior_predictive_intervals(
        identity, [lo], [hi], n_samples=200000, rng_seed=7, quantiles=quantiles
    )
    # stats row: [q2.5, q50, q97.5, mean, std]
    q025, q50, q975, mean, std = stats[0]
    span = hi - lo
    assert abs(q025 - (lo + 0.025 * span)) < 0.05
    assert abs(q50 - (lo + 0.5 * span)) < 0.05
    assert abs(q975 - (lo + 0.975 * span)) < 0.05
    assert abs(mean - (lo + hi) / 2) < 0.05
    assert abs(std - span / np.sqrt(12)) < 0.05
    assert out.shape == (1, 200000)


def test_estimator_invariants_multi_output():
    """Quantiles are ordered (2.5 <= 50 <= 97.5) per output, and a monotone scale
    of the input scales the output interval proportionally."""

    def two_outputs(draws):
        # output 0 = x; output 1 = 3*x, over the same prior.
        x = draws[:, 0]
        return np.column_stack([x, 3.0 * x])  # (n, 2)

    quantiles = [0.025, 0.5, 0.975]
    stats, out = ui.prior_predictive_intervals(
        two_outputs, [0.0], [1.0], n_samples=50000, rng_seed=3, quantiles=quantiles
    )
    assert out.shape == (2, 50000)
    for row in stats:
        assert row[0] <= row[1] <= row[2], "quantiles must be ordered"
    # output 1 interval width is ~3x output 0's (linear scaling).
    w0 = stats[0][2] - stats[0][0]
    w1 = stats[1][2] - stats[1][0]
    assert abs(w1 - 3.0 * w0) < 0.02


def test_death_rate_intervals_are_ordered_and_in_unit_range():
    """Small binding-based run: every condition's prior-predictive interval is
    bounded in [0, 1] and ordered (2.5 <= median <= 97.5). Skipped where the
    compiled `ferroptosis_core` extension is not built; the estimator math is
    covered by the analytic tests above, which need no binding."""
    pytest.importorskip("ferroptosis_core")
    lows = np.array([p[1] for p in ui.PARAMS], float)
    highs = np.array([p[2] for p in ui.PARAMS], float)
    stats, out = ui.prior_predictive_intervals(
        ui.evaluate, lows, highs, n_samples=24, rng_seed=1, quantiles=ui.QUANTILES
    )
    assert out.shape == (len(ui.CONDITIONS), 24)
    assert np.all(out >= 0.0) and np.all(out <= 1.0), "death rates must be in [0, 1]"
    for row in stats:
        assert row[0] <= row[1] <= row[2], "interval quantiles must be ordered"
        assert 0.0 <= row[0] and row[2] <= 1.0
