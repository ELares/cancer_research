"""Guards for the metric-vocabulary measurement behind P1's falsification rule.

THE HAZARD HERE IS THE OPPOSITE OF THE USUAL ONE. Most analyses in this repo
risk claiming more than their counts support; this one started by doing
exactly that -- an 8-versus-0 split inside the ferroptosis subject arm was
rendered as a ranking, which is not a ranking. The fix was a control arm over
the wider census where the counts separate, and the guards exist to keep the
ranking attached to the arm that can carry it.

So the properties pinned are about POWER rather than about the answer: the
subject arm must declare itself unrankable when it is, the ranking must come
from the control, and the two arms' agreement must be stated as agreement
rather than promoted into evidence.

The second hazard is the absolute reading. Counting title and abstract
undercounts every metric heavily, and a table of small numbers invites being
read as a census of practice. The report must keep saying it is not one.

OFFLINE: reads only committed artifacts.
"""
import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
JSON = REPO / "analysis/census-synergy-metrics.json"
MD = REPO / "analysis/census-synergy-metrics.md"
PROTOCOL = REPO / "analysis/p1-wetlab-protocol.md"


@pytest.fixture(scope="module")
def d():
    return json.loads(JSON.read_text())


def test_the_subject_arm_declares_its_own_power(d):
    """Derived from the counts, not asserted.

    If the ferroptosis literature ever names enough metrics to rank them, this
    must flip on its own rather than keep deferring to a control it no longer
    needs.
    """
    assert d["subject_rankable"] == (d["any_metric"] >= d["rankable_min"])
    md = " ".join(MD.read_text().split())
    if not d["subject_rankable"]:
        assert "cannot rank these metrics" in md
        assert f"below the {d['rankable_min']}" in md, (
            "the report does not say what threshold the subject arm failed, "
            "so 'too thin' reads as a judgement rather than a rule")


def test_the_ranking_comes_from_the_control_arm(d):
    """The published order must be the control's order.

    Sorting by the subject arm would reinstate the 8-versus-0 ranking this
    analysis exists to refuse.
    """
    order = [r["control_articles"] for r in d["rows"]]
    assert order == sorted(order, reverse=True), (
        "rows are not ordered by the control arm, so the published ranking is "
        "not the one the counts support")
    assert d["leading_metric"] == d["rows"][0]["metric"]
    assert d["leading_control_articles"] == d["rows"][0]["control_articles"]


def test_the_lead_is_wider_than_the_counting_noise(d):
    """A leader inside Poisson noise of the runner-up is not a leader.

    These are counts of independent occurrences, so the standard error on each
    is about its square root; requiring the gap to exceed two combined errors
    is the weakest bar that still refuses a coin flip.
    """
    a, b = d["leading_control_articles"], d["runner_up_control_articles"]
    se = (a + b) ** 0.5
    assert a - b > 2 * se, (
        f"lead {a} over {b} is within counting noise (2 s.e. = {2 * se:.1f}), "
        "so the control cannot rank these two either")


def test_the_arms_agreement_is_derived_and_stated_as_agreement(d):
    subj = max(d["rows"], key=lambda r: r["combination_articles"])
    expect = subj["metric"] if subj["combination_articles"] else None
    assert d["subject_top_metric"] == expect
    assert d["arms_agree"] == (d["subject_top_metric"] == d["leading_metric"])
    md = " ".join(MD.read_text().split())
    assert ("AGREES" in md) == d["arms_agree"]
    assert "the most it can support" in md, (
        "the report does not mark the agreement as the limit of what the thin "
        "arm establishes, which is how agreement gets read as confirmation")


def test_the_protocol_metrics_are_read_from_the_protocol(d):
    txt = PROTOCOL.read_text().lower()
    named = {r["metric"] for r in d["rows"] if r["in_protocol"]}
    assert named == set(d["protocol_named"])
    assert named, "no protocol metric was detected, so the P1 column is inert"
    for m in named:
        # Every named metric must be findable in the protocol by some spelling
        # its own pattern accepts -- the mapping is many-to-one, so this checks
        # the claim rather than re-running the matcher's exact regex.
        head = re.split(r"[ (]", m)[0].lower()
        assert head in txt, f"{m} is marked in-protocol and is not in it"


def test_the_control_column_is_not_presented_as_a_share(d):
    """It counts records of all kinds, so it is not comparable to the column
    beside it and would be over-read as prevalence if left unexplained."""
    md = " ".join(MD.read_text().split())
    assert "to ORDER the metrics, not to size them" in md
    assert f"{d['control_records']:,} records" in md
    assert f"1-in-{d['control_stride']}" in md, (
        "the control's sampling rate is not stated, so its counts cannot be "
        "related to the census they came from")


def test_the_undercount_and_its_direction_stay_stated(d):
    md = " ".join(MD.read_text().split())
    assert "UNDERCOUNTS every metric" in md
    assert "upper bound on how bad the real gap is" in md, (
        "the claim-to-metric ratio appears without the note that abstract-only "
        "counting flatters it, so it reads as a measured practice rate")
    assert "says nothing about whether a reported synergy is real" in md


def test_the_claim_to_metric_gap_recomputes(d):
    assert d["claim_with_metric_share"] == pytest.approx(
        100 * d["any_metric_combination"] / d["claim_combination"], abs=0.02)
    assert d["claim_only_combination"] == (
        d["claim_combination"] - d["any_metric_combination"])
    assert d["any_metric_combination"] <= d["any_metric"]
    assert d["combination_articles"] <= d["ferroptosis_articles"]
