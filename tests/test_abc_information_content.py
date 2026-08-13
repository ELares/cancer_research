"""The "unconstrained" flag must be judged against a null, not a constant.

WHAT WAS WRONG
--------------
`abc_joint_posterior.py` flagged a parameter unconstrained when its 95% posterior
interval still spanned at least 0.6 of the prior width. That threshold ignores
how many draws were accepted, and that is precisely what sets how narrow an
UNINFORMATIVE posterior looks: with 30 accepted draws, samples taken from the
prior and nothing else still span a median of ~0.90 of the prior width, because
30 points rarely reach the corners.

So 0.6 sat far below anything noise produces, and the flag fired on WELL
determined parameters. It labelled `lp_propagation`, `lp_rate` and `gpx4_rate`
unconstrained while they sit at the 0th percentile of that null — and those are
the cascade parameters the manuscript quotes credible intervals for. The
mislabel had propagated into `joint-posterior.md` ("the parameters the in-vitro
dose-response panel does not identify") and into CLAUDE.md ("the LP cascade ...
stay loosely constrained"), both asserting the opposite of the truth.

WHY THE NULL IS THE RIGHT REFERENCE
-----------------------------------
"Posterior width as a fraction of prior width" has no fixed scale. Its
uninformative value depends only on the number of accepted draws, so the honest
comparison is against that distribution rather than against a number someone
picked. These guards pin the null itself, so a change to the sampler cannot
quietly move what counts as informed.
"""

import json
import math
import random
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
CAL = REPO_ROOT / "analysis" / "calibration"
INFO_JSON = CAL / "abc-information-content.json"
INFO_MD = CAL / "abc-information-content.md"
JOINT_MD = CAL / "joint-posterior.md"
GENERATOR = REPO_ROOT / "scripts" / "abc_joint_posterior.py"


def reports() -> list:
    return json.loads(INFO_JSON.read_text())


def test_the_null_is_what_an_uninformative_posterior_looks_like():
    """Independent re-derivation, so the reported null is not self-certifying.

    Computed here with a different seed and a different generator than the
    analysis uses. If the two disagreed, the null would be an artifact of one
    RNG rather than a property of the statistic.
    """
    r = reports()[0]
    n = r["n_accepted"]
    rng = random.Random(12345)
    widths = []
    for _ in range(4000):
        u = sorted(rng.random() for _ in range(n))
        k_hi = (len(u) - 1) * 0.975
        k_lo = (len(u) - 1) * 0.025
        hi = u[int(k_hi)] + (u[min(int(k_hi) + 1, len(u) - 1)] - u[int(k_hi)]) * (k_hi - int(k_hi))
        lo = u[int(k_lo)] + (u[min(int(k_lo) + 1, len(u) - 1)] - u[int(k_lo)]) * (k_lo - int(k_lo))
        widths.append(hi - lo)
    widths.sort()
    median = widths[len(widths) // 2]
    assert abs(median - r["null_median_width"]) < 0.02, (
        f"independent null median {median:.3f} disagrees with the reported "
        f"{r['null_median_width']}")
    # The property that makes 0.6 wrong: noise does NOT look like 0.6.
    assert median > 0.8, (
        f"an uninformative posterior of {n} draws spans {median:.3f} of the "
        "prior; if this ever fell near 0.6 the original threshold would have "
        "been defensible and this analysis would need revisiting")


def test_the_legacy_threshold_mislabels_informed_parameters():
    """The finding. If it ever stops holding, the corrections must be revisited."""
    mis = [(r["artifact"], n) for r in reports()
           for n, v in r.get("parameters", {}).items()
           if v["informed"] and v["legacy_flag_unconstrained"]]
    assert mis, (
        "no parameter is both informed and flagged unconstrained any more, so "
        "the corrections written into joint-posterior.md and CLAUDE.md no "
        "longer describe the data and should be re-checked")
    names = {n for _, n in mis}
    for expected in ("lp_propagation", "lp_rate"):
        assert expected in names, (
            f"{expected} is no longer among the mislabelled parameters; the "
            "manuscript quotes its credible interval, so this matters")


def test_the_cascade_parameters_are_informed_by_the_data():
    """The substantive claim the corrections assert."""
    joint = [r for r in reports() if "joint-posterior" in r["artifact"]][0]
    for p in ("lp_propagation", "lp_rate", "gpx4_rate"):
        v = joint["parameters"][p]
        assert v["informed"], f"{p} is no longer informed by the data"
        # The decision threshold is 5%, not 1%. This asserted <= 1.0 because on
        # the 1,500-draw run all three sat at 0.0-0.1. The tolerance-anchored
        # 40,000-draw run moved gpx4_rate to 1.5 -- still decisively informed,
        # but a guard pinned to the incidental value rather than the criterion
        # fired on an improvement.
        assert v["null_percentile"] <= 5.0, (
            f"{p} sits at the {v['null_percentile']}th percentile of the null, "
            "above the 5% threshold this analysis uses to call a parameter "
            "informed")


def test_the_genuinely_uninformed_parameters_are_still_reported_as_such():
    """The correction must not overshoot into claiming everything is informed."""
    joint = [r for r in reports() if "joint-posterior" in r["artifact"]][0]
    # `hill` only. k_erastin WAS uninformed on the 1,500-draw quantile run
    # (0.813 of prior width, 8.5th percentile) and became informed under the
    # tolerance-anchored 40,000-draw run (0.391, 0.0th) -- a real gain from the
    # fix, not a drift, so the guard follows it rather than pinning the old set.
    assert not joint["parameters"]["hill"]["informed"], (
        "hill is now reported as informed; it is the parameter measured to be "
        "INERT (the joint distance is identical at hill 6, 8 and 10), so this "
        "would contradict abc-acceptance-diagnostic.md")
    assert joint["parameters"]["hill"]["null_percentile"] > 25, (
        "hill's width is no longer unremarkable for noise")
    assert joint["parameters"]["k_erastin"]["informed"], (
        "k_erastin is uninformed again; the denser tolerance-anchored run had "
        "brought it inside, so this suggests the acceptance fix regressed")


def test_the_generator_no_longer_uses_a_bare_constant():
    src = GENERATOR.read_text()
    assert 'post[n]["posterior_width_frac_of_prior"] >= 0.6' not in src, (
        "the generator is back to the fixed 0.6 threshold")
    assert "null_p5" in src, (
        "the generator does not compute a null-calibrated threshold")


def test_the_history_survives_regeneration():
    """The corrections were written INTO a generated file and were erased.

    joint-posterior.md is produced by abc_joint_posterior.py. Two correction
    notes were added directly to the markdown; the next run of the generator
    overwrote both without warning. The history now lives in the generator's own
    template, so regenerating reproduces it -- which is what this checks.
    """
    txt = JOINT_MD.read_text()
    assert "History of this run's acceptance rule" in txt, (
        "the acceptance-rule history is missing from joint-posterior.md; if it "
        "was written into the markdown directly it will be erased again on the "
        "next regeneration -- it belongs in the generator")
    src = (REPO_ROOT / "scripts" / "abc_joint_posterior.py").read_text()
    assert "History of this run's acceptance rule" in src, (
        "the history is in the generated document but not in the generator, so "
        "it will not survive the next run")


def test_the_analysis_states_what_it_does_not_measure():
    txt = INFO_MD.read_text().lower()
    assert "not whether the model" in txt or "not whether the" in txt, (
        "the report no longer distinguishes 'the data moved this parameter' "
        "from 'the fit is right'")


# --- the identifiability report's closing section --------------------------

IDENT_MD = REPO_ROOT / "analysis" / "identifiability-report.md"
IDENT_GEN = REPO_ROOT / "scripts" / "identifiability_report.py"


def test_the_closing_section_does_not_still_promise_a_route_already_taken():
    """It said #500 + #502 "would condition" the constants. Both landed.

    Neither made any headline point-estimable, so the promise had become a
    deferred note describing where the author stopped looking rather than what
    turned out to be true.
    """
    txt = IDENT_MD.read_text()
    tail = txt[txt.index("## What would make a headline point-estimable"):]
    assert "has been taken, and it is closed" in tail, (
        "the closing section no longer records that the route it used to name "
        "was taken and did not work")
    assert "INADMISSIBLE" in tail, (
        "the closing section does not record the demonstrated reason the "
        "substitution route is closed")


def test_the_closing_figures_are_derived_from_the_artifacts():
    """Every number in that section must come from a committed artifact.

    They were hand-typed beside the two files that compute them, which is the
    shape this repository keeps rediscovering: the artifact moves, the sentence
    does not, and the stale figure reads as freshly checked.
    """
    import re
    txt = IDENT_MD.read_text()
    tail = txt[txt.index("## What would make a headline point-estimable"):]

    h = json.loads((CAL.parent / "headline-at-fitted-cascade.json").read_text())
    worst = max(("ctrpv2_point", "posterior_median"),
                key=lambda k: h[k]["admissibility"]["worst_rate"])
    a = h[worst]["admissibility"]
    assert f"{a['worst_rate']*100:.2f}%" in tail, (
        "the inadmissible untreated-death rate in the prose is not the one in "
        "headline-at-fitted-cascade.json")
    assert f"{h['default']['admissibility']['worst_rate']*100:.2f}%" in tail

    joint = [r for r in reports() if "joint-posterior" in r["artifact"]][0]
    informed = sum(1 for v in joint["parameters"].values() if v["informed"])
    assert f"{informed} of its {len(joint['parameters'])} parameters" in tail, (
        "the informed-parameter count in the prose is not the one in "
        "abc-information-content.json")

    # The generator must CALL the helper, not merely define it. Asserting the
    # name appears anywhere passed a mutation that deleted the call site, because
    # the function's own `def` line contains the name -- a guard satisfied by the
    # thing it is checking for.
    src = IDENT_GEN.read_text()
    body = src[src.index("def write_report("):]
    assert "_fitted_cascade_facts()" in body, (
        "write_report no longer calls _fitted_cascade_facts(); the closing "
        "figures are not being derived even if the helper still exists")
    for literal in ("99.96%", "5 of its 7 parameters are informed"):
        assert literal not in src, (
            f"{literal!r} is hardcoded in the generator again")


def test_the_disjunction_framing_carries_its_qualification():
    """The in-vivo PRCC ranges are +/-50% of the defaults under question.

    So the disjunction restates the falsification rather than supplying
    independent grounds to discount the fit. Without that, "DISJOINT" reads as
    stronger evidence than it is.
    """
    txt = IDENT_MD.read_text()
    assert "DISJOINT" in txt, "the disjunction is no longer reported at all"
    # NOT an `or`. The first version accepted either the explanation or the bare
    # "+/-50%", and the bare number appears elsewhere in the document -- so
    # deleting the sentence that does the explaining still passed.
    assert "restates the falsification" in txt, (
        "the report asserts DISJOINT without the circularity qualification that "
        "analysis/calibration/in-vivo-prior-provenance.md establishes: the "
        "in-vivo PRCC ranges are +/-50% bands around the defaults under "
        "question, so the disjunction restates the falsification rather than "
        "supplying independent grounds to discount the fit")
