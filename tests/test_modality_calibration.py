"""Guards for `analysis/modality-calibration.{md,json}`.

A calibration page is the easiest artifact in a repository to make dishonest,
and the dishonesty is almost never a wrong number. It is a WIDE fit reported as
a fit, a target the form could not have missed, or a clinical endpoint quietly
equated with a model output. So the guards here are about the three refusals
rather than about the fitted values:

* a fit that admits most of the search range is UNCONSTRAINED and must say so;
* an INADMISSIBLE row must carry a diagnosis, not just a verdict;
* every clinical target must carry its mapping, because that is the weak link
  and it is the thing a reader skims past.
"""
import importlib.util
import json
import math
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
MD = REPO / "analysis/modality-calibration.md"
JSON_ = REPO / "analysis/modality-calibration.json"
CORE = REPO / "simulations/ferroptosis-core/src"


def _load():
    spec = importlib.util.spec_from_file_location(
        "calibrate_modality_arms", REPO / "scripts/calibrate_modality_arms.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


CM = _load()


@pytest.fixture(scope="module")
def d():
    return json.loads(JSON_.read_text())


@pytest.fixture(scope="module")
def md():
    return MD.read_text()


def test_every_arm_names_a_source_that_is_in_the_corpus(d):
    """A target with no citation is a preference. Every PMID quoted must be a
    file this repository actually holds, so a number cannot be carried in from
    memory -- which this project has done and retracted before."""
    import re
    for a in d["arms"]:
        assert a["source"].strip(), a["arm"]
        for pmid in re.findall(r"PMID (\d+)", a["source"]):
            f = REPO / "corpus" / "by-pmid" / f"{pmid}.md"
            assert f.exists(), (
                f"{a['arm']} cites PMID {pmid}, which is not in the frozen "
                "corpus; a target quoted from outside it cannot be checked")
        # An arm may ALSO name an external comparator, and two arms do -- the
        # frozen corpus was retrieved from ferroptosis and photo/sonodynamic
        # queries and structurally holds no acoustics-physics literature, so
        # requiring every citation to be corpus-held would forbid comparing
        # against the only relevant study rather than making anything
        # checkable. The split is the point: `source` stays corpus-checked and
        # a comparator must SAY it is outside, so an unverifiable citation
        # cannot be moved into the checked field to escape this guard.
        ext = a.get("comparator_source")
        if ext:
            assert "OUTSIDE the frozen corpus" in ext or "outside the frozen corpus" in ext, (
                f"{a['arm']} names a comparator without saying it is outside "
                "the corpus, which is the whole reason it lives in a separate "
                f"field: {ext[:120]}")
            assert re.search(r"\b(19|20)\d\d\b", ext), (
                f"{a['arm']}'s comparator has no year; an external citation "
                "must carry enough to be located")


def test_a_wide_fit_is_called_unconstrained_rather_than_fitted(d):
    """The failure mode this page exists to avoid.

    A fit admitting most of the search range has been given a target that
    cannot discriminate. Reporting it as calibrated is the same error as
    reporting a p-value without an effect size.
    """
    cap = d["unconstrained_width"]
    assert 0.0 < cap < 1.0, cap
    for a in d["arms"]:
        f = a["fit"]
        # A DIRECTIONAL arm is exempt from the width test and NOT from
        # scrutiny: it has no fitted range because its target constrains a
        # sign rather than a value, so what is checked instead is that it
        # actually scored some directions and did not simply decline to be
        # measured. An exemption that asserted nothing would be the easiest
        # place in this file to hide an unfitted arm.
        if a.get("directional") is not None:
            assert a["verdict"] in ("DIRECTIONAL", "PARTLY REFUTED"), a["arm"]
            dd = a["directional"]
            assert dd["confirmed"] + dd["refuted"] >= 2, (
                f"{a['arm']} claims a directional verdict on fewer than two "
                "scored claims, which is not a verdict")
            assert (a["verdict"] == "PARTLY REFUTED") == (dd["refuted"] > 0), (
                f"{a['arm']}'s verdict and its refuted count disagree -- a "
                "refuted direction must show in the headline, since that is "
                "the one a reader would otherwise never see")
            assert f is None
            continue
        if a["target"] is None:
            # A row with NO TARGET is a different outcome from one whose
            # target nothing satisfies, and collapsing them would hide the
            # more interesting admission: that this project has no number.
            assert a["verdict"] == "NO TARGET", a["arm"]
            assert f is None
            continue
        if f is None:
            assert a["verdict"] == "INADMISSIBLE", a["arm"]
            continue
        wide = f["width_fraction"] > cap
        assert (a["verdict"] == "UNCONSTRAINED") == wide, (
            f"{a['arm']} has width {f['width_fraction']:.2f} against a cap of "
            f"{cap} and is called {a['verdict']}")
        assert 0.0 <= f["width_fraction"] <= 1.0
        assert f["lo"] <= f["hi"]


def test_the_fits_are_recomputed_not_stored(d):
    """The artifact must be what the script produces, not what it produced
    once. Recomputed live here, so a hand-edited number fails."""
    live = CM.assemble(CM.scan())
    assert [a["arm"] for a in live["arms"]] == [a["arm"] for a in d["arms"]]
    for got, want in zip(live["arms"], d["arms"]):
        assert got["verdict"] == want["verdict"], got["arm"]
        if got["fit"] is None:
            assert want["fit"] is None
            continue
        for k in ("lo", "hi", "width_fraction"):
            assert math.isclose(got["fit"][k], want["fit"][k], rel_tol=1e-9), (
                f"{got['arm']}.{k} drifted")


def test_an_inadmissible_row_carries_a_diagnosis_not_just_a_verdict(d, md):
    """INADMISSIBLE is the outcome the page calls worth wanting, so it has to
    earn that by explaining itself. A verdict alone would be an unexplained
    failure dressed as a finding."""
    bad = [a for a in d["arms"] if a["verdict"] == "INADMISSIBLE"]
    if not bad:
        # Every arm fitted. That is allowed, and the page must then not be
        # claiming an inadmissible row it does not have.
        assert "The inadmissible row is the informative one" not in md
        return
    assert "The inadmissible row is the informative one" in md
    for a in bad:
        diag = a.get("diagnosis", "")
        assert len(diag) > 120, (
            f"{a['arm']} is INADMISSIBLE with a {len(diag)}-character "
            "diagnosis; a verdict without an explanation is not a finding")
        assert diag in md, f"{a['arm']}'s diagnosis is not in the report"


def test_the_checkpoint_diagnosis_is_arithmetic_that_still_holds(d):
    """The identifiability finding is the campaign's most valuable result, so
    it is recomputed from the crate rather than trusted.

    If the parameters ever change so that the arm CAN reach its target, this
    fails and the paragraph must be rewritten -- which is the intent.
    """
    import re
    arm = next((a for a in d["arms"] if a["arm"] == "Checkpoint blockade"), None)
    assert arm, "the checkpoint row is gone"
    params = (CORE / "params.rs").read_text()
    g = lambda f: float(re.search(rf"\b{f}: ([0-9.]+),", params).group(1))
    coef = (g("dc_maturation_rate") * g("tcell_priming_rate") * 0.02
            * (1.0 - g("pd1_brake") * (1.0 - g("anti_pd1_efficacy"))))
    assert math.isclose(arm["ceiling_at_full_presentation"], min(coef, 1.0),
                        rel_tol=1e-9)
    assert math.isclose(arm["antigenicity_required"], 0.20 / coef, rel_tol=1e-9)
    if arm["verdict"] == "INADMISSIBLE":
        assert arm["antigenicity_required"] > 1.0, (
            "the required antigenicity is now within its own bound, so the "
            "row should fit and the diagnosis is stale")
        assert arm["ceiling_at_full_presentation"] < 0.20


def test_a_row_with_no_target_says_so_rather_than_inventing_one(d, md):
    """The outcome a calibration page is most tempted to hide.

    Inventing a target that a flexible form then satisfies looks exactly like
    calibration and constrains nothing, so a row with no number must SAY it
    has none and say why -- and must not carry a fit.
    """
    none = [a for a in d["arms"] if a["verdict"] == "NO TARGET"]
    if not none:
        assert "The rows with nothing to fit to" not in md
        return
    assert "The rows with nothing to fit to" in md
    assert "no target at all" in md
    for a in none:
        assert a["target"] is None and a["fit"] is None, a["arm"]
        reason = a.get("no_target_reason", "")
        assert len(reason) > 120, (
            f"{a['arm']} has no target and a {len(reason)}-character reason; "
            "an unexplained absence is not an admission")
        assert reason in md
        # The SOURCE must still exist -- the layer landed under the
        # layer-freeze policy on a mechanism anchor, and losing that would
        # mean the layer should not have landed.
        assert a["source"].strip(), a["arm"]


def test_an_unconstrained_row_is_named_and_explained(d, md):
    """A 99%-wide window is consistent with the model and says nothing about
    it. That has to read as a fact about the TARGET, not as a fit."""
    unc = [a for a in d["arms"] if a["verdict"] == "UNCONSTRAINED"]
    if not unc:
        assert "The rows whose target excludes almost nothing" not in md
        return
    assert "The rows whose target excludes almost nothing" in md
    for a in unc:
        assert a["fit"]["width_fraction"] > d["unconstrained_width"]
        assert f"**{a['arm']}**" in md


def test_every_clinical_target_carries_its_mapping(d, md):
    """The weak link, and the one a reader skims past.

    A clinical endpoint and a lattice kill fraction are different quantities.
    Any row whose target is clinical must say what it assumed to compare them.
    """
    clinical = [a for a in d["arms"] if "clinical" in a["target_kind"]]
    assert clinical, "no clinical targets at all; check the classification"
    for a in clinical:
        assert a["mapping"], f"{a['arm']} has a clinical target and no mapping"
        assert a["mapping"] in md, f"{a['arm']}'s mapping is not in the report"
    # And a NON-clinical row must not claim one, or the distinction is noise.
    for a in d["arms"]:
        if "clinical" not in a["target_kind"]:
            assert not a["mapping"], (
                f"{a['arm']} is not a clinical target and carries a mapping")
    # The CAR-T mapping must carry the caveat from its own source, because a
    # fit to the headline band describes the indication these therapies were
    # approved for and not the setting this engine simulates.
    cart = next((a for a in d["arms"] if a["arm"].startswith("CAR-T")), None)
    if cart:
        assert "not transferred to solid" in cart["mapping"] or \
               "NOT transferred to solid" in cart["mapping"], (
            "the CAR-T row fits the leukaemia band without carrying the "
            "caveat its own source attaches to it")


def test_the_page_refuses_to_call_a_fit_a_validation(md):
    """Three sentences stop this being a victory lap, and each names a limit
    that would change a claim if lifted."""
    for frag in (
        "A fitted parameter is not a validated model",
        "Hitting a target is weak evidence when nothing could miss it",
        "None of these fits is wired into a binary",
        "used in any reported number",
    ):
        assert frag in md, f"the page no longer says: {frag}"


def test_only_one_parameter_moves_per_arm(d):
    """Fitting two would leave most of these targets underdetermined, and an
    underdetermined fit reporting a point estimate is worse than no fit."""
    for a in d["arms"]:
        assert isinstance(a["parameter"], str) and "," not in a["parameter"], (
            f"{a['arm']} fits more than one parameter: {a['parameter']}")
        if a["range_scanned"] is None:
            assert a["verdict"] in ("NO TARGET", "DIRECTIONAL",
                                    "PARTLY REFUTED"), a["arm"]
            # NO TARGET means there is nothing to fit; DIRECTIONAL means there
            # is a target and it constrains a sign. Collapsing them would hide
            # the more interesting admission, which is the same reason the
            # page splits NO TARGET from UNCONSTRAINED.
            assert (a["verdict"] == "NO TARGET") == (a["target"] is None), (
                f"{a['arm']}: a directional arm must NAME its target and an "
                "arm with no target must not claim a direction")
            continue
        lo, hi = a["range_scanned"]
        assert lo < hi and lo >= 0.0


def test_the_cart_fit_reads_the_rust_default_and_depends_on_it():
    """The module claimed this default was load-bearing while it was not.

    `adoptive.rs` documented -- and `CALIBRATION_STATUS.md` repeated -- that a
    solid-tumour default would move this fit. It could not have: the fit
    reimplemented the kill in Python with no barrier term, and never opened
    `adoptive.rs`. A reviewer found it by tracing callers rather than grepping
    for the type name, and re-ran the generator with the default replaced to
    show the artifact regenerated byte-identically.

    Two things are checked, because reading the file is not the same as
    depending on it: the parsed default must be the leukaemia case the
    published band was measured in, and substituting solid-tumour barriers
    must actually change the outcome.
    """
    mod = _load()
    barriers = mod._rust_default_barriers()
    assert barriers["trafficking"] == 1.0
    assert barriers["infiltration"] == 1.0
    assert barriers["activation"] == 1.0
    assert barriers["exhaustion_rate"] == 0.0
    assert barriers["antigen_positive_fraction"] == 1.0, (
        "the B-ALL band is measured with every barrier open; a default that "
        "applies solid-tumour barriers makes this fit describe something else")

    baseline = mod._cart_arm()["fit"]
    assert baseline is not None, "the leukaemia fit no longer lands"

    real = mod._rust_default_barriers
    try:
        mod._rust_default_barriers = lambda: {
            "trafficking": 0.3, "infiltration": 0.4, "activation": 0.5,
            "exhaustion_rate": 0.02, "antigen_positive_fraction": 0.8,
        }
        moved = mod._cart_arm()["fit"]
    finally:
        mod._rust_default_barriers = real
    assert moved != baseline, (
        "the fitted band is unchanged under solid-tumour barriers, so the "
        "default is not load-bearing and the claim that it is must be "
        "retracted rather than repeated")
