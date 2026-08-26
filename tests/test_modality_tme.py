"""Guards for `analysis/modality-tme.{md,json}`.

This is the page that turns "the engine can express radiation" into "the
engine has something to say about radiation", so its guards are about the
three ways such a page misleads:

* an axis that moved nothing is listed as tested;
* an effect is attributed to an axis that co-varied with another;
* an ordering is reported as a result when it was put in by hand.
"""
import importlib.util
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
MD = REPO / "analysis/modality-tme.md"
JSON_ = REPO / "analysis/modality-tme.json"


def _load():
    spec = importlib.util.spec_from_file_location(
        "modality_tme_report", REPO / "scripts/modality_tme_report.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


TM = _load()


@pytest.fixture(scope="module")
def d():
    return json.loads(JSON_.read_text())


@pytest.fixture(scope="module")
def md():
    return MD.read_text()


def test_the_sweep_is_a_full_factorial(d):
    """Every combination of the three axes, or the paired comparisons below
    have no partner to compare against."""
    conds = d["conditions"]
    seen = {(c["hypoxic"], c["stroma"], c["acidic"]) for c in conds}
    assert len(seen) == 8, f"{len(seen)} of 8 combinations present: {sorted(seen)}"
    assert len(conds) == 8
    arms = {a for c in conds for a in c["arms"]}
    for c in conds:
        assert set(c["arms"]) == arms, (
            "an arm is missing from some conditions, so its effect sizes are "
            "computed over a different set than the others")


def test_effects_are_paired_so_an_axis_cannot_borrow_anothers(d):
    """The failure this page would otherwise be prone to.

    Comparing "all hypoxic rows" against "all normoxic rows" attributes to
    hypoxia anything that co-varied with it. `_effect` pairs each ON condition
    with the OFF condition identical in the other two axes, and this
    recomputes it to make sure it still does.
    """
    conds = d["conditions"]
    for label, key in (("hypoxia", "hypoxic"), ("stromal shielding", "stroma"),
                       ("acidic pH", "acidic")):
        for arm in d["arms"]:
            assert abs(TM._effect(conds, arm, key) - d["effects"][label][arm]) < 1e-12
    # And the pairing must be real: a mutated sweep where an axis co-varies
    # with another must NOT show up as that axis's effect.
    doctored = [dict(c) for c in conds]
    for c in doctored:
        c["arms"] = dict(c["arms"])
        if c["stroma"]:                      # make STROMA rows differ wildly
            c["arms"]["SDT"] = 0.01
    stroma_effect = TM._effect(doctored, "SDT", "stroma")
    hypoxia_effect = TM._effect(doctored, "SDT", "hypoxic")
    assert stroma_effect > 0.5, "the pairing cannot see a real stroma effect"
    # Hypoxia's paired effect must be computed within matched stroma states,
    # so it is not inflated by the stroma change.
    assert hypoxia_effect < stroma_effect, (
        "a change confined to the stroma axis leaked into hypoxia's effect; "
        "the comparisons are not paired")


def test_an_axis_that_moved_nothing_is_called_inert(d, md):
    """A table of 0.08% differences invites a reader to believe an axis was
    tested when the configuration could not see it."""
    thr = d["inert_threshold"]
    assert 0 < thr < 0.5
    # THE THRESHOLD ITSELF MUST BE MEANINGFUL, or this guard just agrees with
    # whatever it is set to. Lowering it to 1e-7 reclassified two inert axes
    # as live and every assertion below still passed, because they were all
    # derived from the same number. So a LIVE axis must actually move
    # something substantially -- an axis promoted by a tiny threshold moves
    # nothing and fails here.
    for label in d["live_axes"]:
        worst = max(d["effects"][label].values(), default=0.0)
        assert worst >= 0.10, (
            f"{label} is reported LIVE and moves every arm by at most "
            f"{worst:.4f}. Either the threshold has been lowered until an "
            "axis that does nothing counts as tested, or the axis genuinely "
            "stopped biting and the page's ordering paragraph is stale.")
    for label, effs in d["effects"].items():
        worst = max(effs.values(), default=0.0)
        inert = worst < thr
        assert (label in d["inert_axes"]) == inert, (
            f"{label} moves arms by at most {worst:.4f} and is called "
            f"{'inert' if label in d['inert_axes'] else 'live'}")
    if d["inert_axes"]:
        assert "INERT in this configuration" in md
        # The page must say WHY each is invisible, not merely that it is.
        assert "nothing for the trapping to scale" in md
        assert "overwhelm it by an order of magnitude" in md
        # And must not claim the arms are robust to it.
        assert "they were not visible" in md


def test_the_ordering_follows_the_mechanisms(d, md):
    """The result is the ORDERING, so it must be derived and it must match
    what the mechanisms imply -- a threshold arm loses nothing, an
    oxygen-dependent arm loses most."""
    live = d["live_axes"]
    assert live, "no axis discriminates at all; the sweep shows nothing"
    worst = {a: max(d["effects"][l][a] for l in live) for a in d["arms"]}
    assert d["robustness_order"] == sorted(d["arms"], key=lambda a: worst[a])
    # Ablation is a threshold: it must be unaffected, and if it ever is not
    # the model has stopped being self-consistent.
    assert worst["Ablation"] == 0.0, (
        f"ablation lost {worst['Ablation']:.3f} of its kill to an axis; a "
        "destroyed cell does not care about its oxygen tension")
    # The oxygen-dependent arm must lose more than the dose-modified one.
    assert worst["SDT"] > worst["Radiation"] > 0.0, (
        f"SDT {worst['SDT']:.3f} vs Radiation {worst['Radiation']:.3f}: an "
        "arm whose lethality DEPENDS on oxygen should lose more than one "
        "whose dose is merely modified by it")
    assert "That ordering was not tuned for" in md


def test_the_page_refuses_the_three_over_readings(md):
    """Each names a limit that would change a claim if lifted."""
    for frag in ("they were not visible",
                 "is a PREDICTION, not a measurement",
                 "the model being consistent, not a result",
                 "the ORDERING is the result, and the numbers are not"):
        assert frag in md, f"the page no longer says: {frag}"
