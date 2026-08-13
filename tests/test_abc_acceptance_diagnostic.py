"""Guards for the ABC acceptance diagnostic.

THE FINDING
-----------
The joint ABC returned a posterior whose median fits worse than a vector already
committed in the repository. The cause is not a truncated prior and not a small
sample: acceptance is a fixed 2% FRACTION, so the run always keeps its best 2%
however bad they are, and the reported epsilon is an output rather than a
criterion. With uniform draws almost never reaching the good region in seven
dimensions, that cut lands well above what is achievable.

These guards pin the three claims separately, because they fail in different
ways: the numbers come from the committed artifact, the "not a truncated prior"
conclusion has to follow from its own test rather than be asserted, and the
report must keep saying what the finding does NOT overturn.
"""

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
CAL = REPO_ROOT / "analysis" / "calibration"
DIAG_JSON = CAL / "abc-acceptance-diagnostic.json"
DIAG_MD = CAL / "abc-acceptance-diagnostic.md"
JOINT = CAL / "joint-posterior.json"
GENERATOR = REPO_ROOT / "scripts" / "abc_joint_posterior.py"


def diag() -> dict:
    return json.loads(DIAG_JSON.read_text())


def test_the_fix_holds_and_the_document_reports_it_as_resolved():
    """The defect is fixed; this now guards the FIX rather than the bug.

    Written first as "the committed vector beats the posterior median", which was
    the observation at the time. Changing the acceptance rule inverted it, and
    the guard fired -- correctly, because a document describing a live defect
    must not survive the defect's removal unchanged.
    """
    d = diag()
    assert d["posterior_median_distance"] < d["committed_distance"], (
        "the posterior median is worse than the reference again; the acceptance "
        "fix has regressed")
    assert d["reported_epsilon"] <= d["committed_distance"] * 1.15, (
        f"epsilon {d['reported_epsilon']} is no longer anchored near the "
        f"achievable {d['committed_distance']}; it has drifted back toward being "
        "an output rather than a criterion")
    assert "Status: RESOLVED" in DIAG_MD.read_text()


def test_the_prior_is_not_what_truncates_the_fit():
    """The conclusion must FOLLOW from the test, not sit beside it.

    The committed vector sits on two prior bounds, which looks like clipping.
    Stepping outside has to show the fit not improving, or the document's
    "it is not a truncated prior" section is wrong.
    """
    t = diag()["prior_truncation_test"]
    at_bound = t["k_erastin"]["3.0"]
    outside = [v for k, v in t["k_erastin"].items() if float(k) < 3.0]
    assert outside, "the k_erastin truncation test no longer probes outside the bound"
    assert all(v >= at_bound for v in outside), (
        "k_erastin improves outside its prior bound, so the prior IS truncating "
        "and the diagnostic's central negative claim is wrong")


def test_hill_is_inert_rather_than_merely_unidentified():
    """A parameter that changes nothing is a stronger statement than a wide CrI."""
    hl = diag()["prior_truncation_test"]["hill"]
    assert len(set(hl.values())) == 1, (
        f"hill now changes the joint distance ({hl}); the report says it is inert, "
        "which is why the information-content analysis found it uninformed")


def test_the_good_region_is_rare_but_reachable():
    """A RATE, not an absolute zero.

    The first version asserted 0 draws beat the reference, taken from a 300-draw
    sample that was simply too small to see the region. At 20,000 draws nine do.
    "Never" and "about 5 in 10,000" support the same argument -- a 1,500-draw run
    lands in the region essentially never -- but only one of them is true.
    """
    s_ = diag()["sampling"]
    rate = s_["n_beating_committed"] / s_["draws"]
    assert rate < 1e-3, (
        f"the good region is now reached at {rate:.1e} per draw; a plain "
        "rejection sampler would find it readily and the diagnosis needs redoing")
    assert s_["draws"] >= 10000, (
        "the rate is measured on too few draws to distinguish 'rare' from "
        "'absent' -- the earlier 300-draw sample reported zero and was wrong")


def test_epsilon_is_now_anchored_rather_than_floating():
    """Epsilon must sit NEAR the achievable distance, not far above it.

    This asserted the gap was LARGE, which was the defect's signature. Under the
    fix epsilon is the reference times the tolerance factor, so the same quantity
    is now small and positive -- the assertion had to invert with the thing it
    measures.
    """
    d = diag()
    gap = d["epsilon_excess_over_committed"]
    assert gap is not None, "the epsilon-to-reference gap is no longer recorded"
    assert 0 < gap <= 0.15, (
        f"epsilon is {100*gap:.0f}% above the achievable distance; anchored "
        "acceptance should keep it within the tolerance factor")


def test_the_report_says_what_it_does_not_overturn():
    """It must not read as 'the posterior is worthless'."""
    txt = DIAG_MD.read_text().lower()
    assert "does not claim" in txt or "does not overturn" in txt, (
        "the report no longer bounds what the fix does and does not establish")
    assert "not a validated one" in txt, (
        "the report no longer says a better-fitting posterior is not a validated "
        "one -- the overclaim this whole arc is most exposed to")
    assert "inadmissible" in txt, (
        "the report no longer points at the in-vivo regime where substituting "
        "these values remains inadmissible")


def test_the_acceptance_rule_is_a_tolerance_not_a_quantile():
    """The fix itself. This guard used to assert the OPPOSITE.

    It was written as "the rule is still a fixed fraction -- if that ever changes
    the document must be reworded", the right shape for a guard on a known-live
    defect. Fixing the rule fired it, the document was reworded, and it now pins
    the fix.
    """
    src = GENERATOR.read_text()
    assert "n_accept = max(10, int(args.n_draws * ACCEPT_FRAC))" not in src, (
        "the generator is back to a fixed-fraction acceptance")
    assert "distances <= eps" in src, "acceptance is no longer a tolerance"
    assert "REFERENCE_VECTOR" in src and "TOLERANCE_FACTOR" in src, (
        "the tolerance is no longer anchored to a reachable reference distance")
    assert "underpowered" in src, (
        "the run no longer reports a shortfall; padding back to a quota is the "
        "behaviour that produced the original defect")
