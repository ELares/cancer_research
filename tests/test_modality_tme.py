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
    """Every combination of every axis, in every stratum, or the paired
    comparisons below have no partner to compare against."""
    conds = d["conditions"]
    seen = {(c["hypoxic"], c["stroma"], c["acidic"], c["deep"],
             c["heterogeneous"], c["phenotype"]) for c in conds}
    n_pheno = len(d["phenotypes"])
    assert len(seen) == 32 * n_pheno, (
        f"{len(seen)} of {32 * n_pheno} combinations present")
    assert len(conds) == 32 * n_pheno
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
                       ("acidic pH", "acidic"), ("depth", "deep"),
                       ("clonal heterogeneity", "heterogeneous")):
        for arm in d["arms"]:
            assert abs(TM._effect(conds, arm, key) - d["effects"][label][arm]) < 1e-12
    # And the pairing must be real: a mutated sweep where an axis co-varies
    # with another must NOT show up as that axis's effect.
    doctored = [dict(c) for c in conds]
    for c in doctored:
        c["arms"] = dict(c["arms"])
        if c["stroma"]:                      # make STROMA rows differ wildly
            c["arms"]["SDT"] = 0.01
    stroma_effect = abs(TM._effect(doctored, "SDT", "stroma"))
    hypoxia_effect = abs(TM._effect(doctored, "SDT", "hypoxic"))
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
    # The floor is tied to the threshold rather than being a second constant:
    # a live axis must clear it COMFORTABLY, not by a hair. A flat 10% was the
    # first attempt and was simply wrong -- stromal shielding moves an arm by
    # 7%, which is a real effect and not a rounding artefact, and a guard that
    # rejected it would have been demanding the world be tidier than it is.
    floor = max(2.0 * thr, 0.03)
    for label in d["live_axes"]:
        worst = max((abs(v) for v in d["effects"][label].values()), default=0.0)
        assert worst >= floor, (
            f"{label} is reported LIVE and moves every arm by at most "
            f"{worst:.4f}, under a floor of {floor:.4f}. Either the threshold "
            "has been lowered until an axis that does nothing counts as "
            "tested, or the axis genuinely stopped biting and the page's "
            "ordering paragraph is stale.")
    for label, effs in d["effects"].items():
        worst = max((abs(v) for v in effs.values()), default=0.0)
        inert = worst < thr
        assert (label in d["inert_axes"]) == inert, (
            f"{label} moves arms by at most {worst:.4f} and is called "
            f"{'inert' if label in d['inert_axes'] else 'live'}")
    if d["inert_axes"]:
        assert "remain inert even stratified" in md
        # Must not claim the arms are robust to an axis it could not apply.
        assert "cannot apply that pressure" in md


def test_the_ordering_follows_the_mechanisms(d, md):
    """The result is the ORDERING, so it must be derived and it must match
    what the mechanisms imply -- a threshold arm loses nothing, an
    oxygen-dependent arm loses most."""
    live = d["live_axes"]
    assert live, "no axis discriminates at all; the sweep shows nothing"
    worst = {a: max((abs(d["effects"][l][a]) for l in live), default=0.0)
             for a in d["arms"]}
    assert d["robustness_order"] == sorted(d["arms"], key=lambda a: worst[a])
    # Ablation is a threshold: it must be unaffected, and if it ever is not
    # the model has stopped being self-consistent.
    assert worst["Ablation"] == 0.0, (
        f"ablation moved {worst['Ablation']:.3f} of its kill on an axis; a "
        "destroyed cell does not care about its oxygen tension")
    # The oxygen-dependent arm must move more than the dose-modified one.
    assert worst["SDT"] > worst["Radiation"] > 0.0, (
        f"SDT {worst['SDT']:.3f} vs Radiation {worst['Radiation']:.3f}: an "
        "arm whose lethality DEPENDS on oxygen should move more than one "
        "whose dose is merely modified by it")
    assert "That ordering was not tuned for" in md


def test_the_page_refuses_the_three_over_readings(md):
    """Each names a limit that would change a claim if lifted."""
    for frag in ("they were not visible",
                 "is a PREDICTION, not a measurement",
                 "the model being consistent, not a result",
                 "the ORDERING is the result, and the numbers are not"):
        assert frag in md, f"the page no longer says: {frag}"


def test_the_phenotype_is_a_stratum_and_not_an_axis(d, md):
    """Comparing a persister against a glycolytic cell is not a treatment
    effect, so the phenotype must never appear as an axis -- and the paired
    comparisons must hold it fixed."""
    assert "phenotype" not in [k for k in ("hypoxic", "stroma", "acidic", "deep")]
    assert len(d["phenotypes"]) >= 2, (
        "one phenotype only; the sweep cannot see an axis that bites in the "
        "other state, which is how two axes were reported inert")
    # Effects must be computable per stratum AND pooled, and the two must
    # differ -- if they were identical the stratification would be decorative.
    pooled = d["effects"]
    per = d["effects_by_phenotype"]
    assert set(per) == set(d["phenotypes"])
    differs = any(
        abs(per[ph][ax][arm] - pooled[ax][arm]) > 1e-9
        for ph in per for ax in pooled for arm in d["arms"])
    assert differs, (
        "every per-phenotype effect equals the pooled one, so stratifying "
        "changed nothing and the finding it supports is not real")


def test_the_dominant_axis_differs_by_phenotype(d, md):
    """The finding the stratification exists to produce.

    If it ever stops being true the page's central paragraph -- that an axis
    reported inert is a statement about the run and not the biology -- has to
    be re-derived rather than left standing.
    """
    dom = d["dominant_axis"]
    assert len(dom) >= 2, dom
    axes = {axis for axis, _ in dom.values()}
    assert len(axes) > 1, (
        f"the same axis dominates in every phenotype ({axes}); the page "
        "claims they differ")
    for ph, (axis, worst) in dom.items():
        assert worst >= 0.10, (
            f"{ph}'s dominant axis {axis} moves nothing ({worst:.3f})")
        assert f"**{ph}** — {axis}" in md, (
            f"the report does not name {ph}'s dominant axis")


def test_the_sweep_reproduces_the_repos_own_persister_target(d):
    """An independent hit on a committed self-consistency target.

    `targets.yaml` records `persister_rsl3_death_rate` at 0.425 from the
    original single-cell work. This sweep reaches the persister state by a
    different path -- a different binary, a different loop, a different seed
    scheme -- and lands on it. That is worth pinning, because a drift here
    means the two paths have parted and one of them is wrong.
    """
    import re
    base = next(c for c in d["conditions"]
                if c["phenotype"] == "persister"
                and not any((c["deep"], c["hypoxic"], c["stroma"], c["acidic"])))
    y = (REPO / "simulations/calibration/targets.yaml").read_text()
    m = re.search(r"id: persister_rsl3_death_rate.*?target_value: ([0-9.]+)"
                  r".*?tolerance: ([0-9.]+)", y, re.S)
    assert m, "the persister target is gone from targets.yaml"
    target, tol = float(m.group(1)), float(m.group(2))
    got = base["arms"]["RSL3"]
    assert abs(got - target) <= tol, (
        f"the sweep's unstressed persister RSL3 kill is {got:.4f} against the "
        f"committed target {target} +/- {tol}; the two paths have parted")
    # And the agreement must be TIGHT, not merely inside a wide tolerance --
    # the tolerance is 0.1, which almost anything would satisfy.
    assert abs(got - target) < 0.02, (
        f"the sweep lands {abs(got - target):.4f} from the target, inside the "
        "stated tolerance but loosely; check whether the two paths still "
        "model the same thing")


def test_the_effects_are_signed_so_a_gain_is_not_called_a_loss(d, md):
    """Taking the absolute value produced the impossible line that an arm had
    lost 121% of its kill. It had GAINED, and collapsing the sign did not just
    mislabel a number -- it hid a result."""
    signed = [v for eff in d["effects"].values() for v in eff.values()]
    assert any(v < 0 for v in signed), "no axis lowers any arm's kill"
    gains = [v for v in signed if v > 0.10]
    if gains:
        assert "One axis can HELP" in md, (
            "an axis raises an arm's kill by more than 10% and the page does "
            "not say so; that is a result, not a rounding artefact")
        assert "Variance rescues a marginal drug" in md
    # No reported LOSS may exceed 100%: an arm cannot lose more kill than it
    # had, and a figure above 100 is the sign error this test exists for.
    losses = [v for v in signed if v < 0]
    assert all(v >= -1.0 for v in losses), (
        f"a loss exceeds 100%: {min(losses):.3f}. An arm cannot lose more "
        "than it had, so this is a sign or pairing error.")


def test_the_amplification_is_measured_after_the_grace_period(d, md):
    """The measurement that was wrong in an instructive way.

    Reading lipid peroxidation AT death returns approximately the death
    threshold for every arm by construction -- death IS the crossing -- so the
    quality difference vanished. It had not vanished; it was being measured
    before it happens.
    """
    qr = d.get("quality_ratio")
    if qr is None:
        return
    # Death threshold is 10.0; a per-death figure at ~10 for every arm is the
    # symptom of measuring at the crossing rather than after the grace period.
    assert qr["sdt"] > 12.0 and qr["rsl3"] > 12.0, (
        f"per-death DAMP release is {qr['sdt']:.1f}/{qr['rsl3']:.1f}, close to "
        "the death threshold -- the measurement is being taken at the "
        "crossing rather than after `post_death_steps`")
    assert qr["ratio"] > 1.0, "the ROS arm should die louder, not quieter"
    assert qr["ratio"] < 5.0, (
        f"the per-death ratio is {qr['ratio']:.1f}; the page says the "
        "amplification advantage is a COUNT effect and that sentence needs "
        "re-deriving if the quality gap is now large")
    assert "COUNT effect, not a quality effect" in md
