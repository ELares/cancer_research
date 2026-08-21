"""Guards for the shared-axis modality comparison and the artifact it found.

The analysis answers half of what Section 9.4 asks and found, in passing, that
a number this project had put in its own manuscript was measuring a different
modality: `Ultrasonic Therapy` gives sonodynamic therapy a 4.54% clinical-trial
share, and only 2 of those 114 trials mention sonodynamic therapy at all. The
rest are 1980s ultrasound hyperthermia.

THREE THINGS ARE GUARDED, in descending order of how badly they fail silently.

1. The correction must reach every site. A retraction that lives in one
   analysis while the manuscript still quotes the old figure is worse than not
   finding it, because the figure now has a citation.
2. The comparison must be drawn on ONE arm. Radioligand therapy -- the
   alternative the section turns on -- has no descriptor, so a mixed-arm
   comparison would measure the alternative one way and the thesis modality
   another, exactly where the argument is.
3. The two axes it CANNOT measure must stay named. Delivery burden and
   applicability to residual disease reduce to nothing in a bibliographic
   record, and a reader who forgets that will take a half-comparison for a
   whole one.

OFFLINE: these read only committed artifacts.
"""
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
JSON = REPO / "analysis/census-modality-comparison.json"
MD = REPO / "analysis/census-modality-comparison.md"
MANUSCRIPT = REPO / "article/drafts/v1.md"
# Pinned literals, not read from the artifact.
DISAGREE_AT = 3.0
MIN_TRIALS = 10


@pytest.fixture(scope="module")
def d():
    return json.loads(JSON.read_text())


def test_shares_and_arm_ratios_recompute(d):
    for r in d["rows"]:
        for arm in ("mesh", "text"):
            art, tri = r[f"{arm}_articles"], r[f"{arm}_trials"]
            if art:
                assert r[f"{arm}_trial_share"] == pytest.approx(
                    100 * tri / art, abs=0.02), (r["modality"], arm)
            assert r[f"{arm}_interpretable"] == (tri >= MIN_TRIALS)
            ph = r[f"{arm}_phases"]
            assert r[f"{arm}_phased"] == sum(ph.values())
            assert r[f"{arm}_late_phase"] == ph["III"] + ph["IV"]
        a, b = r["text_trial_share"], r["mesh_trial_share"]
        if a and b:
            assert r["arm_ratio"] == pytest.approx(max(a, b) / min(a, b), abs=0.06)
            assert r["arms_disagree"] == (r["arm_ratio"] >= DISAGREE_AT)
        else:
            assert r["arm_ratio"] is None and not r["arms_disagree"]


def test_the_comparison_is_drawn_on_an_arm_that_covers_every_modality(d):
    """A mixed-arm comparison would fall exactly where the argument is."""
    assert d["common_arm"] == "text", (
        "the comparison is not drawn on the arm covering every modality; "
        "radioligand therapy has no descriptor, so a descriptor-arm comparison "
        "would measure the named alternative one way and the thesis modality "
        "another")
    for r in d["rows"]:
        assert r["text_measurable"], (
            f"{r['modality']} is not measurable on the common arm, so the "
            "comparison cannot include it")
    assert d["mesh_missing"], (
        "no modality lacks a descriptor, which would mean the design note "
        "about instrument mismatch no longer describes this data")


def test_the_descriptor_artifact_is_derived_and_named(d):
    """The finding that corrected the manuscript.

    It must be reached by the test rather than asserted: a modality whose two
    arms disagree by three-fold or more is flagged, and sonodynamic must be
    among them or the artifact has silently gone away.
    """
    flagged = set(d["arms_disagree"])
    assert flagged == {r["modality"] for r in d["rows"] if r["arms_disagree"]}
    assert "sonodynamic" in flagged, (
        "sonodynamic no longer fails the arm-agreement test. Either the data "
        "changed or the test weakened; do not remove the manuscript correction "
        "without establishing which")
    assert set(d["validity_failed"]) <= flagged
    # The other modalities must pass, or the test is flagging its own method
    # rather than a property of one descriptor.
    passing = [r for r in d["rows"] if r["mesh_trial_share"] and not r["arms_disagree"]]
    assert len(passing) >= 4, (
        "most modalities fail the arm-agreement test, which means the two arms "
        "disagree generally and the sonodynamic finding is not specific")
    md = MD.read_text()
    assert "exactly **2** mention sonodynamic" in md
    assert "hyperthermia" in md.lower()


def test_every_failure_is_traced_rather_than_labelled(d):
    """Three failures with three different causes and two directions.

    A test that flags mismatch does not say which arm is right. Reporting them
    under one label -- "descriptor problems" -- would hide that one is a
    descriptor too broad, one too narrow, and one too young, and that only the
    first warranted a manuscript correction.
    """
    md = MD.read_text()
    for m in d["validity_failed"]:
        assert f"**`{m}`**" in md, f"{m} fails the test and is not diagnosed"
    # The three causes must be distinguished, not merged.
    assert "too BROAD" in md and "too NARROW" in md and "too YOUNG" in md, (
        "the failures are reported under a single label; they have different "
        "causes and only one of them implies a correction")
    assert "does not tell you which arm is right" in md


def test_the_validity_floor_is_on_articles_not_trials(d):
    """The blind spot that would have hidden the motivating case.

    A descriptor that inflates a share does it by piling trials onto the
    descriptor arm while the term's own arm stays thin -- so requiring trials
    on BOTH arms excludes exactly the worst artifacts. Sonodynamic has 4 text
    trials against 114 descriptor trials and would have been skipped.
    """
    assert "sonodynamic" in d["validity_tested"], (
        "sonodynamic is not in the tested set, so the sweep cannot see the "
        "artifact it was built from")
    assert "sonodynamic" in d["hidden_by_a_trial_floor"], (
        "sonodynamic no longer has a thin arm, which would mean the blind-spot "
        "note no longer describes this data")
    for m in d["validity_tested"]:
        r = next(x for x in d["rows"] if x["modality"] == m)
        assert min(r["mesh_articles"], r["text_articles"]) >= d["min_articles_for_validity"]
        assert max(r["mesh_trials"], r["text_trials"]) >= d["min_trials"]
    md = MD.read_text()
    assert "floor is on ARTICLES, not trials" in md
    assert "blind to the thing that motivated it" in md


def test_most_mechanisms_pass_or_the_test_measures_its_own_method(d):
    tested, failed = len(d["validity_tested"]), len(d["validity_failed"])
    assert tested >= 12, f"only {tested} mechanisms testable; the sweep is thin"
    assert failed < tested / 2, (
        f"{failed} of {tested} fail, so the two arms disagree generally and a "
        "single failure is not evidence about a single descriptor")


def test_the_correction_reached_the_manuscript():
    """A retraction living in one analysis while the manuscript quotes the old
    figure is worse than not finding it -- the figure now has a citation."""
    txt = " ".join(MANUSCRIPT.read_text().split())
    assert "0.29%" in txt, (
        "the manuscript does not carry the corrected sonodynamic trial share")
    # The old figure may still appear, but only inside its own correction.
    for m in ("4.54%",):
        if m in txt:
            assert "98% borrowed" in txt or "only 2 of those 114" in txt, (
                f"the manuscript still quotes {m} without the correction "
                "beside it")
    # And the comparison that changed by a factor of 15 must not still be
    # stated at its old size.
    assert "a factor of 1.6 below HIFU" not in txt, (
        "the manuscript still states the HIFU-to-sonodynamic gap at the size "
        "the descriptor artifact produced")


def test_the_two_unmeasurable_axes_stay_named(d):
    """A reader who forgets which half is missing takes a half-comparison for
    a whole one."""
    md = MD.read_text()
    assert "delivery burden" in md.lower()
    assert "applicability to residual disease" in md.lower()
    assert "half of the comparison" in md.lower() or "half the job" in md.lower()


def test_it_does_not_read_a_trial_share_as_evidence_of_efficacy(d):
    md = MD.read_text()
    assert "## What a trial share is not" in md
    assert "not trials that read out well" in md or "not trials that succeeded" in md
    for overclaim in ("proves radioligand", "shows sonodynamic does not work",
                      "more effective than"):
        assert overclaim not in md.lower()


def test_the_section_argument_is_not_overstated(d):
    """Being further along a development path is not meeting the same clinical
    need, and the report must not let the first stand in for the second."""
    md = MD.read_text()
    assert "does not establish meeting the same need" in md
    assert "targetable receptor" in md
