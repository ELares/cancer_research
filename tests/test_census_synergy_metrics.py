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
import gzip
import re
from pathlib import Path
import sys

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
import census_synergy_metrics as scanner
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


def test_the_observed_control_lead_is_descriptive(d):
    """Shared-record matches and chronological shards are not independent counts."""
    a, b = d["leading_control_articles"], d["runner_up_control_articles"]
    assert a > b
    markdown = MD.read_text()
    assert "descriptive ordering" in markdown
    assert "do not supply an independent-count uncertainty model" in markdown


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


def test_mentions_are_not_promoted_to_computation_rates_or_inferred_bias(d):
    md = " ".join(MD.read_text().split())
    assert "Metric-specific reporting rates are unknown" in md
    assert "neither the true computation rate nor a bound or bias direction" in md
    assert "undercount applies to every metric alike" not in md
    assert "output is a Bliss excess" not in md
    assert "says nothing about whether a reported synergy is real" in md


def test_the_claim_to_metric_gap_recomputes(d):
    assert d["claim_with_metric_share"] == pytest.approx(
        100 * (d["claim_combination"] - d["claim_only_combination"])
        / d["claim_combination"], abs=0.02)
    assert d["claim_combination"] - d["claim_only_combination"] <= d["any_metric_combination"]
    assert d["any_metric_combination"] <= d["any_metric"]
    assert d["combination_articles"] <= d["ferroptosis_articles"]


def test_the_manuscript_does_not_offer_the_metrics_as_an_either_or(d):
    """Preserve compatible metric reporting without inflating mention counts."""
    manuscript = (REPO / "article/drafts/v1.md").read_text()
    txt = " ".join(manuscript.split())
    either_or = "Chou-Talalay analysis or compare against Bliss independence"
    if either_or in txt:
        assert "should be an *and*" in txt or "Report both" in txt, (
            "the manuscript offers the two synergy metrics as alternatives "
            "with nothing saying why both are needed")
        assert f"{d['control_records']:,} records" in txt, (
            "the correction is present without the measurement behind it")
    paragraph = next(p for p in manuscript.split("\n\n") if "analysis/census-synergy-metrics.md" in p)
    assert f"{d['control_records']:,} records" in paragraph
    assert f"Chou-Talalay appears in {d['control_metric']['Chou-Talalay combination index']} titles or abstracts" in paragraph
    assert f"Bliss in {d['control_metric']['Bliss independence']}" in paragraph
    for phrase in (
        "1.99× output is a Bliss ratio",
        "mention counts, not measurements of which metrics were computed",
        "differential reporting could change their ordering",
        "do not establish that Bliss-based results lack published comparators",
        "Report both metrics when the dose-response data support them",
        "publish the raw dose-response matrix",
    ):
        assert phrase in paragraph, f"the manuscript dropped the qualification: {phrase}"
    assert "nothing published beside it" not in paragraph


def _records(tmp_path, monkeypatch, records):
    with gzip.open(tmp_path / "records.jsonl.gz", "wt") as stream:
        for record in records:
            stream.write(json.dumps(record) + "\n")
    monkeypatch.setattr(scanner, "RECORDS", tmp_path)
    return scanner.assemble(scanner.scan())


def test_control_keeps_its_own_sorted_stride_when_subject_stride_changes(tmp_path, monkeypatch):
    for index in reversed(range(81)):
        with gzip.open(tmp_path / f"{index:03}.jsonl.gz", "wt") as stream:
            stream.write(json.dumps({"mesh": ["Ferroptosis"],
                                     "title": "Bliss in combination"}) + "\n")
    monkeypatch.setattr(scanner, "RECORDS", tmp_path)
    result = scanner.scan(3)
    assert result["ferroptosis_articles"] == 27
    assert result["metric_combination"] == {"Bliss independence": 27}
    assert result["control_records"] == 3
    assert result["control_stride"] == 40
    assert result["control_metric"] == {"Bliss independence": 3}


def test_zero_matches_do_not_fabricate_a_leader_or_arm_agreement(tmp_path, monkeypatch):
    result = _records(tmp_path, monkeypatch, [{"title": "unrelated observation"}])
    assert result["control_records"] == 1
    assert result["leading_metric"] is None
    assert result["subject_top_metric"] is None
    assert result["protocol_names_the_leader"] is False
    assert result["arms_agree"] is False
    markdown = scanner.render(result)
    assert "no unique leading metric" in markdown
    assert "share is unavailable" in markdown
    assert "None" not in markdown
    assert "the control can" not in markdown
    assert "Computing both indices is the right call" not in markdown


def test_tied_positive_counts_have_no_unique_leader(tmp_path, monkeypatch):
    result = _records(tmp_path, monkeypatch, [{"mesh": ["Ferroptosis"],
        "title": "Bliss and Loewe in combination"}])
    assert result["leading_control_articles"] == 1
    assert result["leading_metric"] is None
    assert result["subject_top_metric"] is None
    assert result["arms_agree"] is False


def test_claim_metric_share_counts_the_intersection_not_all_metric_articles(tmp_path, monkeypatch):
    result = _records(tmp_path, monkeypatch, [
        {"mesh": ["Ferroptosis"], "title": "synergy in combination"},
        {"mesh": ["Ferroptosis"], "title": "Bliss in combination"},
        {"mesh": ["Ferroptosis"], "title": "Loewe in combination"},
    ])
    assert result["any_metric_combination"] == 2
    assert result["claim_combination"] == result["claim_only_combination"] == 1
    assert result["claim_with_metric_share"] == 0
    assert "0 of them (0.0%)" in scanner.render(result)


def test_report_derives_rankability_and_does_not_assume_protocol_names_leader(tmp_path, monkeypatch):
    result = _records(tmp_path, monkeypatch, [
        {"mesh": ["Ferroptosis"], "title": "Loewe synergy in combination"}
        for _ in range(40)
    ])
    assert result["subject_rankable"] is True
    assert result["protocol_names_the_leader"] is False
    markdown = scanner.render(result)
    assert "meets the 30-article threshold" in markdown
    assert "protocol does not name the metric" in markdown
    assert "subject arm cannot rank" not in markdown
    assert "The protocol already computes both" not in markdown


def test_existing_numerical_artifact_is_unchanged_by_reassembly(d):
    assert scanner.assemble(d) == d
