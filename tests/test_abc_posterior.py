"""Tests for the #332 ABC posterior (data-conditioned, in-vitro).

The ABC run itself needs the compiled `ferroptosis_core` extension and is not
re-run in CI, so these tests cover the prior-design guard and the committed
posterior result (structure, ordering, the load-bearing in-vivo/in-vitro
disjunction, and the held-out coverage).
"""

import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import abc_posterior as abc  # noqa: E402

ABC_JSON = REPO_ROOT / "analysis" / "calibration" / "abc-posterior.json"


def test_priors_span_the_330_point_calibration():
    """The in-vitro priors must bracket the #330 point fit (lp_propagation=0.7,
    lp_rate=0.4, K=0.25); otherwise the posterior would be biased by truncation."""
    bounds = {name: (lo, hi) for name, lo, hi in abc.PRIORS}
    for name, val in (("lp_propagation", 0.7), ("lp_rate", 0.4), ("k_um", 0.25)):
        lo, hi = bounds[name]
        assert lo < val < hi, f"{name}={val} not strictly inside prior [{lo},{hi}]"


def test_priors_exceed_invivo_prcc_ranges():
    """The in-vitro priors must extend ABOVE the in-vivo PRCC ranges (the whole
    point: the in-vitro regime is out of the in-vivo prior's reach)."""
    bounds = {name: (lo, hi) for name, lo, hi in abc.PRIORS}
    for name, (lo, hi) in abc.INVIVO_PRCC.items():
        assert bounds[name][1] > hi


def test_committed_posterior_structure_and_ordering():
    r = json.loads(ABC_JSON.read_text())
    assert r["n_accepted"] < r["n_draws"]
    assert r["fit_compound"] == "ML162" and r["heldout_compound"] == "ML210"
    for name in ("lp_propagation", "lp_rate", "k_um", "gpx4_rate"):
        p = r["posterior"][name]
        assert p["q2_5"] <= p["median"] <= p["q97_5"]  # quantiles ordered


def test_invivo_invitro_disjunction_is_the_result():
    """The load-bearing finding: the in-vitro posterior for the lp cascade lies
    ENTIRELY above the in-vivo PRCC ranges, so the in-vivo priors cannot be
    conditioned on the in-vitro data."""
    r = json.loads(ABC_JSON.read_text())
    dj = r["invivo_prior_disjunction"]
    for name in ("lp_propagation", "lp_rate"):
        assert dj[name]["posterior_above_invivo_range"] is True
        assert dj[name]["posterior_q2_5"] > dj[name]["invivo_prcc_range"][1]


def test_heldout_posterior_predictive_coverage_reports_strict_and_tolerant():
    """Coverage is reported BOTH strictly (inside the 95% band) and within a small
    viability tolerance, so the strict number is never hidden behind the tolerant
    one. The tolerant band should bracket most points; the strict count is lower
    (single-cell model vs cell-line-median curves), and that is documented."""
    r = json.loads(ABC_JSON.read_text())
    strict, total_s = map(int, r["heldout_coverage_strict"].split("/"))
    tol, total_t = map(int, r["heldout_coverage_tolerant"].split("/"))
    assert total_s == total_t == 7
    assert tol >= strict  # tolerance can only add coverage
    assert tol >= 5        # the band is in the right place
    assert strict >= 1     # at least some points are covered without any tolerance
    assert 0.0 < r["heldout_tolerance"] <= 0.1


def test_posterior_in_invitro_regime():
    """Posterior medians sit in the in-vitro regime (above the in-vivo max),
    consistent with the #330 point calibration."""
    r = json.loads(ABC_JSON.read_text())
    assert r["posterior"]["lp_propagation"]["median"] > 0.2
    assert r["posterior"]["lp_rate"]["median"] > 0.12
