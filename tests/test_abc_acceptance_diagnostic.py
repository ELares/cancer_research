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


def test_the_committed_vector_beats_the_posterior_median():
    """The observation the whole diagnostic exists to explain."""
    d = diag()
    assert d["committed_distance"] < d["posterior_median_distance"], (
        "the posterior median now fits at least as well as the committed vector; "
        "the diagnostic's premise no longer holds and its prose must be revisited")
    assert d["committed_distance"] < d["reported_epsilon"], (
        "the committed vector is outside the acceptance region, which would be a "
        "different problem than the one documented")


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


def test_uniform_sampling_does_not_reach_the_good_region():
    d = diag()["sampling"]
    assert d["n_beating_committed"] == 0, (
        f"{d['n_beating_committed']} of {d['draws']} uniform draws now beat the "
        "committed vector; the sampling argument needs rechecking")
    assert d["best"] > diag()["committed_distance"]


def test_the_epsilon_gap_is_quantified_not_asserted():
    d = diag()
    gap = d["epsilon_excess_over_committed"]
    assert gap is not None and gap > 0.1, (
        "the reported epsilon is no longer meaningfully above the achievable "
        "distance, so the 'shell' framing does not apply")
    assert f"{100*gap:.0f}%" in DIAG_MD.read_text(), (
        "the epsilon gap in the prose is not the one in the artifact")


def test_the_report_says_what_it_does_not_overturn():
    """It must not read as 'the posterior is worthless'."""
    txt = DIAG_MD.read_text()
    assert "does not overturn" in txt.lower()
    assert "informed" in txt, (
        "the report no longer records that five of seven parameters are still "
        "genuinely informed by the data")


def test_the_acceptance_rule_is_still_a_fixed_fraction():
    """If this is ever fixed, the diagnostic becomes historical and must say so."""
    src = GENERATOR.read_text()
    fixed = "n_accept = max(10, int(args.n_draws * ACCEPT_FRAC))" in src
    assert fixed, (
        "the generator no longer accepts a fixed fraction. That is the fix this "
        "diagnostic recommends — good — but the document now describes the past "
        "and must be reworded as historical")
