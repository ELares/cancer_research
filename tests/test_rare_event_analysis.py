"""Guards for the rare-event sweep's statistics and its classification rule.

WHY THESE AND NOT OTHERS
------------------------
The sweep's whole output is two things a reader takes on trust: an exact Poisson
interval computed from a hand-rolled gamma quantile, and a label saying whether
a condition's number is a measurement or an artifact of the sample size. Both
were checked once in a shell during development, which is exactly the kind of
verification that does not survive the session it happened in.

The Poisson interval matters more than it looks. At 42% death a normal
approximation is fine and nobody would notice a slightly wrong interval. At one
event in a million the interval spans two orders of magnitude, the approximation
is meaningless, and a bug in the gamma quantile would produce a plausible-looking
number with no way to eyeball it. So it is checked against values that come from
outside this repo.
"""

import importlib.util
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))


def _load(name):
    """Import a script by path; `scripts/` is not a package."""
    spec = importlib.util.spec_from_file_location(
        name, REPO_ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


sweep = _load("rare_event_sweep")
analysis = _load("rare_event_analysis")


# --- the interval ---------------------------------------------------------
# Reference values are the standard exact (Garwood) Poisson interval, which is
# defined through the chi-square quantiles:
#     lower = chi2(0.025, 2k)/2,  upper = chi2(0.975, 2k+2)/2
# These come from the textbook definition, not from running this code.

def test_zero_events_gives_the_rule_of_three():
    """k=0 has no lower bound, and the upper is ~3/n -- the rule of three."""
    lo, hi = sweep.poisson_ci(0, 1_000_000)
    assert lo == 0.0
    # chi2(0.975, 2)/2 = 3.6889, so the exact bound is 3.6889e-6 and the
    # familiar 3/n is the round-number approximation of it.
    assert abs(hi - 3.6889e-6) < 1e-9, hi
    assert 3.0 / 1_000_000 < hi, "3/n must be inside the exact bound, not above it"


def test_a_single_event_spans_two_orders_of_magnitude():
    """The case the normal approximation cannot handle at all."""
    lo, hi = sweep.poisson_ci(1, 1_000_000)
    # chi2(0.025,2)/2 = 0.02532, chi2(0.975,4)/2 = 5.5716
    assert abs(lo - 2.532e-8) < 1e-10, lo
    assert abs(hi - 5.572e-6) < 1e-9, hi
    assert hi / lo > 200, "a single event really is this uninformative"


def test_a_large_count_gives_a_tight_symmetric_interval():
    """The sanity end: with hundreds of events it should look ordinary."""
    lo, hi = sweep.poisson_ci(451, 1_000_000)
    assert abs(lo - 4.104e-4) < 1e-6, lo
    assert abs(hi - 4.947e-4) < 1e-6, hi
    rate = 451 / 1_000_000
    assert lo < rate < hi
    # Near-symmetry on a log scale is what a large count should produce.
    assert abs((rate - lo) - (hi - rate)) / rate < 0.10


def test_the_interval_always_contains_the_point_estimate():
    for k in (0, 1, 2, 7, 50, 451, 11_375):
        lo, hi = sweep.poisson_ci(k, 1_000_000)
        assert lo <= k / 1_000_000 <= hi, k


def test_the_bound_falls_exactly_one_decade_per_decade_of_n():
    """The property the whole figure is built to expose.

    With zero events the bound is proportional to 1/n, so it is a straight line
    of slope -1 on log-log axes -- a statement about the sample, not the
    biology. If this ever stopped holding, the figure's reference line would be
    wrong and every 'resolution-limited' reading with it.
    """
    prev = None
    for n in (10**6, 10**7, 10**8, 10**9):
        _, hi = sweep.poisson_ci(0, n)
        if prev is not None:
            assert abs(prev / hi - 10.0) < 1e-6, (n, prev, hi)
        prev = hi


# --- the classification ---------------------------------------------------

def _row(n, dead):
    return {"n_cells": n, "n_dead": dead,
            "death_rate": dead / n,
            "zero_event_upper_bound_95": 3.0 / n,
            "poisson_ci_low": 0.0, "poisson_ci_high": 3.0 / n}


def test_all_zero_is_resolution_limited():
    rows = [_row(10**6, 0), _row(10**9, 0)]
    assert analysis.classify(rows) == "resolution-limited"
    assert analysis.tracks_rule_of_three(rows)


def test_zero_then_events_is_emergent():
    """The most informative outcome: it locates where the tail begins."""
    rows = [_row(10**6, 0), _row(10**9, 4)]
    assert analysis.classify(rows) == "emergent"
    assert not analysis.tracks_rule_of_three(rows)


def test_events_throughout_is_resolved():
    rows = [_row(10**6, 451), _row(10**9, 451_000)]
    assert analysis.classify(rows) == "resolved"
    assert not analysis.tracks_rule_of_three(rows)


def test_a_single_event_at_the_top_still_counts_as_emergent():
    """One death in a hundred billion is a real departure from the floor.

    This is the boundary the sweep exists to find, so it must not be rounded
    away into 'resolution-limited'.
    """
    rows = [_row(10**6, 0), _row(10**9, 0), _row(10**11, 1)]
    assert analysis.classify(rows) == "emergent"
    assert not analysis.tracks_rule_of_three(rows)


# --- the count parser -----------------------------------------------------

def test_the_conditions_swept_are_the_ones_at_the_resolution_limit():
    """Each swept condition must carry its reason, since spending a hundred
    billion cell-simulations on a condition needs one."""
    assert sweep.CONDITIONS, "no conditions defined"
    for pheno, tx, why in sweep.CONDITIONS:
        assert pheno and tx
        assert len(why) > 40, f"{pheno}/{tx} has no stated rationale"


# --- the nesting witness --------------------------------------------------
# The prose claims the samples are nested, and the evidence it offers is a
# real count sequence. That sequence used to be typed out by hand next to the
# generated table it described, which is the shape that goes stale silently.

def test_the_witness_is_computed_and_reports_the_real_counts():
    by_cond = {("PersisterNrf2", "Control"):
               [_row(10**6, 1), _row(10**7, 7), _row(10**8, 76)]}
    w = analysis._witness(by_cond)
    assert "1, 7, 76" in w
    assert "PersisterNrf2 + Control" in w


def test_the_witness_refuses_to_claim_monotonicity_it_does_not_have():
    """If nesting ever broke, the sentence asserting it must not still print.

    A larger sample containing a smaller one cannot lose a death, so a
    decreasing count means the seeding assumption is wrong -- and that is
    exactly when a hand-written 'the counts are monotone' would be a lie.
    """
    broken = {("X", "Y"): [_row(10**6, 9), _row(10**7, 2)]}
    w = analysis._witness(broken)
    assert "NOT monotone" in w
    assert "investigated" in w


def test_the_witness_says_so_when_nothing_has_events_yet():
    assert "no condition has events" in analysis._witness(
        {("X", "Y"): [_row(10**6, 0)]})


# --- the count parser's zero boundary ------------------------------------

def test_zero_cells_is_rejected_by_both_parse_branches():
    """`--cells 0` used to parse as an exact usize and skip the guard.

    The guard covered only the float-shorthand branch, so "0" returned Ok(0),
    the run reported death_rate and all three means as null (0/0 is NaN, which
    serde_json writes as null) with exit 0, and the sweep driver then failed
    reading the null -- a bad argument surfacing as a parse error one layer
    away from its cause.
    """
    import subprocess
    root = Path(__file__).resolve().parent.parent
    exe = root / "simulations" / "target" / "release" / "sim-scale"
    if not exe.exists():
        import pytest
        pytest.skip("sim-scale not built")
    for bad in ("0", "0.0", "-1", "0e9"):
        r = subprocess.run([str(exe), "--cells", bad, "--phenotype", "Glycolytic",
                            "--treatment", "RSL3"], capture_output=True, text=True)
        assert r.returncode != 0, f"--cells {bad} was accepted"
        assert "null" not in r.stdout, f"--cells {bad} emitted a record of nulls"
