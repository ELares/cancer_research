"""Validate the self-contained Morris estimator in scripts/headline_sensitivity.py.

The estimator math is tested on an analytic function with a KNOWN sensitivity
structure (no simulation binary needed, so this runs in the normal Python CI):
a linear + interaction function whose elementary effects are exact, so the
mu_star ranking and the sigma (interaction) signal are predictable.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

hs = pytest.importorskip("headline_sensitivity")


def test_param_set_matches_prcc():
    """The screened set is exactly the 11 PRCC rate constants (apples-to-apples)."""
    assert len(hs.PARAM_NAMES) == 11
    for name in ("lp_propagation", "gpx4_rate", "lp_rate", "sdt_ros", "rsl3_gpx4_inhib"):
        assert name in hs.PARAM_NAMES
    assert hs.LOWS.shape == hs.HIGHS.shape == (11,)
    assert np.all(hs.HIGHS > hs.LOWS)


def test_morris_trajectory_changes_one_coord_at_a_time():
    rng = np.random.default_rng(0)
    k, levels = 5, 4
    delta = levels / (2.0 * (levels - 1.0))
    traj = hs.morris_trajectory(k, levels, delta, rng)
    assert traj.shape == (k + 1, k)
    changed = set()
    for r in range(k):
        diff = np.abs(traj[r + 1] - traj[r])
        # Exactly one coordinate moves on each step...
        nonzero = np.where(diff > 1e-9)[0]
        assert len(nonzero) == 1, f"step {r} changed {len(nonzero)} coords"
        # ...by the Morris step size delta...
        assert abs(diff[nonzero[0]] - delta) < 1e-9
        changed.add(int(nonzero[0]))
    # ...and every coordinate is perturbed exactly once across the trajectory.
    assert changed == set(range(k))
    assert traj.min() >= -1e-9 and traj.max() <= 1.0 + 1e-9


def test_morris_recovers_known_sensitivity_ranking():
    # f(x) = 5*x0 + 3*x1 + 0*x2 + 2*x0*x3 on the unit cube.
    #   EE_0 = 5 + 2*x3  (large; varies with x3 -> sigma > 0)
    #   EE_1 = 3         (linear -> sigma ~ 0)
    #   EE_2 = 0         (inert -> mu_star ~ 0)
    #   EE_3 = 2*x0      (interaction-only -> sigma > 0)
    def f(rows):
        x = np.asarray(rows)
        return 5 * x[:, 0] + 3 * x[:, 1] + 2 * x[:, 0] * x[:, 3]

    lows = np.zeros(4)
    highs = np.ones(4)
    mu_star, sigma = hs.morris_indices(f, lows, highs, n_traj=60, levels=4, rng_seed=7)

    order = list(np.argsort(-mu_star))
    # x0 is the strongest, x2 the weakest (inert).
    assert order[0] == 0, f"expected x0 top, got {order} (mu*={mu_star})"
    assert order[-1] == 2, f"expected x2 (inert) last, got {order}"
    # x1 (linear, coeff 3) outranks x3 (interaction-only, mean effect ~1).
    assert mu_star[1] > mu_star[3]
    # The inert parameter has ~zero importance.
    assert mu_star[2] < 0.2
    # Interaction/nonlinear params show sigma; the purely-linear x1 does not.
    # (EE_0 = 5 + 2*x3 has std ~ 2*std(x3) ~ 0.75 over the Morris grid; the linear
    # EE_1 = 3 is constant so its sigma is ~0.)
    assert sigma[1] < 0.3, f"linear param should have ~0 sigma, got {sigma[1]}"
    assert sigma[0] > 0.4, f"interaction param x0 should have sigma>0, got {sigma[0]}"
    assert sigma[3] > sigma[1], "interaction-only x3 should out-sigma the linear x1"


def test_default_binary_lookup_is_none_when_absent(tmp_path, monkeypatch):
    # _default_binary returns None when neither release nor debug binary exists,
    # so main() can print a build hint rather than crash.
    monkeypatch.setattr(hs, "REPO", tmp_path)
    # The function reads REPO at call time via the module global.
    assert hs._default_binary() is None
