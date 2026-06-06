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


def test_partition_drops_nan_and_failures():
    """Regression for the blocker: a NaN synergy (sim-combo-mech emits NaN when
    Bliss <= 0.001, i.e. both single agents ~0%) must be dropped as
    undefined-Bliss, not pass through into np.percentile (where one NaN poisons
    every quantile). Subprocess failures (None) are a distinct category."""
    results = [1.0, float("nan"), None, 2.0, float("inf"), 1.5, None]
    finite, n_failed, n_undefined = hu._partition(results)
    assert sorted(finite.tolist()) == [1.0, 1.5, 2.0]
    assert n_failed == 2  # the two None
    assert n_undefined == 2  # nan + inf
    # The surviving array is all-finite, so the summary is never poisoned.
    assert np.all(np.isfinite(finite))
    s = hu._pctiles(finite)
    assert np.isfinite(s["median"]) and np.isfinite(s["p2_5"]) and np.isfinite(s["p97_5"])


def test_partition_all_undefined_yields_empty():
    finite, n_failed, n_undefined = hu._partition([float("nan"), float("nan"), None])
    assert finite.size == 0 and n_undefined == 2 and n_failed == 1


def test_partition_tme_drops_nan_per_observable_and_counts_failures():
    """sim-tme partitioning: a run dict with a non-finite observable is dropped
    from THAT observable's array (independently), and None is a run failure.
    Guards the same NaN-poisons-percentile contract as the Bliss path, per
    observable, without the costly binary."""
    results = [
        {"hypoxia": 0.9, "immune": 0.04},
        {"hypoxia": float("nan"), "immune": 0.05},  # hypoxia dropped, immune kept
        None,  # run failure
        {"hypoxia": 0.8, "immune": float("inf")},  # immune dropped, hypoxia kept
        {"hypoxia": 0.7, "immune": 0.03},
    ]
    hyp, imm, n_failed = hu._partition_tme(results)
    assert n_failed == 1
    assert sorted(hyp.tolist()) == [0.7, 0.8, 0.9]  # the NaN-hypoxia row dropped
    assert sorted(imm.tolist()) == [0.03, 0.04, 0.05]  # the inf-immune row dropped
    assert np.all(np.isfinite(hyp)) and np.all(np.isfinite(imm))


def test_partition_tme_all_failed_yields_empty():
    hyp, imm, n_failed = hu._partition_tme([None, None])
    assert hyp.size == 0 and imm.size == 0 and n_failed == 2


def test_partition_penetration_drops_nan_per_tissue_and_counts_failures():
    """Penetration partitioning: per-tissue non-finite death rates are dropped
    from THAT tissue's array independently; None is a run failure. Same NaN
    contract as the other headlines, across the three tissue scenarios, without
    the binary."""
    k0, k1, k2 = (k for k, _ in hu.PENETRATION_TISSUES)
    results = [
        {k0: 0.12, k1: 0.03, k2: 0.02},
        {k0: float("nan"), k1: 0.04, k2: 0.03},  # k0 dropped, k1/k2 kept
        None,  # run failure
        {k0: 0.20, k1: 0.05, k2: float("inf")},  # k2 dropped, k0/k1 kept
    ]
    per_tissue, n_failed = hu._partition_penetration(results)
    assert n_failed == 1
    assert sorted(per_tissue[k0].tolist()) == [0.12, 0.20]
    assert sorted(per_tissue[k1].tolist()) == [0.03, 0.04, 0.05]
    assert sorted(per_tissue[k2].tolist()) == [0.02, 0.03]
    for k in (k0, k1, k2):
        assert np.all(np.isfinite(per_tissue[k]))


def test_partition_penetration_missing_key_treated_as_nan():
    """A result dict missing a tissue key must not KeyError; that tissue's value
    is treated as non-finite (dropped) for that draw."""
    k0, k1, k2 = (k for k, _ in hu.PENETRATION_TISSUES)
    per_tissue, n_failed = hu._partition_penetration([{k0: 0.1, k1: 0.2}])  # k2 absent
    assert n_failed == 0
    assert per_tissue[k0].tolist() == [0.1]
    assert per_tissue[k2].size == 0  # missing key -> dropped, no crash


def test_ordering_preserved_fraction_is_within_draw_paired():
    """The per-draw ordering test (well >= poorly >= CNS) must count ONLY draws
    with all three finite, and detect a per-draw inversion that overlapping
    marginals would hide."""
    k0, k1, k2 = (k for k, _ in hu.PENETRATION_TISSUES)
    results = [
        {k0: 0.9, k1: 0.4, k2: 0.1},  # ordered
        {k0: 0.2, k1: 0.5, k2: 0.1},  # INVERTED (poorly > well) -> not counted
        {k0: 0.3, k1: 0.3, k2: 0.3},  # ties ok (>=)
        None,  # failure -> excluded from valid
        {k0: float("nan"), k1: 0.2, k2: 0.1},  # incomplete -> excluded from valid
    ]
    frac, n_valid = hu._ordering_preserved_fraction(results)
    assert n_valid == 3  # the two ordered/tie draws + the inverted one
    assert abs(frac - 2 / 3) < 1e-9  # 2 of 3 valid draws monotone


def test_ordering_preserved_fraction_empty_is_nan():
    import math

    frac, n_valid = hu._ordering_preserved_fraction([None, None])
    assert n_valid == 0 and math.isnan(frac)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
