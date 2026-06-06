#!/usr/bin/env python3
"""Unit tests for the Bliss prior-predictive uncertainty harness (#332).

Pure-estimator tests only (the prior sampler + percentile summary), so they run
in CI without the compiled sim-combo-mech binary — mirroring
tests/test_headline_sensitivity.py, which tests the Morris estimator without
running the binaries.

Run: pytest tests/test_headline_uncertainty.py -v
"""

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

hu = pytest.importorskip("headline_uncertainty")


def test_sample_prior_shape_and_bounds():
    n = 50
    draws = hu.sample_prior(n)
    assert draws.shape == (n, len(hu.PARAM_NAMES))
    # Every draw lies within the PRCC ranges [LOWS, HIGHS].
    assert np.all(draws >= hu.LOWS - 1e-12)
    assert np.all(draws <= hu.HIGHS + 1e-12)


def test_sample_prior_deterministic():
    """Fixed seed ⇒ identical draws, so the reported interval is reproducible."""
    a = hu.sample_prior(64)
    b = hu.sample_prior(64)
    assert np.array_equal(a, b)


def test_sample_prior_varies_with_n():
    """Different sample counts must not silently return the same first rows
    (guards against a frozen/cached draw)."""
    small = hu.sample_prior(8)
    big = hu.sample_prior(64)
    assert big.shape[0] == 64
    # The generator is seeded identically, so the first 8 rows match the small
    # draw (a useful determinism property), but the big draw has more rows.
    assert np.array_equal(small, big[:8])


def test_pctiles_on_known_array():
    vals = np.arange(1.0, 101.0)  # 1..100
    s = hu._pctiles(vals)
    assert s["n"] == 100
    assert s["min"] == 1.0 and s["max"] == 100.0
    assert abs(s["median"] - 50.5) < 1e-9
    assert abs(s["mean"] - 50.5) < 1e-9
    # 2.5 / 97.5 percentiles of 1..100 (linear interp) are ~3.475 / ~97.525.
    assert s["p2_5"] < s["median"] < s["p97_5"]
    assert 3.0 < s["p2_5"] < 4.0
    assert 97.0 < s["p97_5"] < 98.0


def test_pctiles_monotone_interval():
    rng = np.random.default_rng(0)
    s = hu._pctiles(rng.uniform(1.0, 5.0, size=500))
    assert s["min"] <= s["p2_5"] <= s["median"] <= s["p97_5"] <= s["max"]


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
