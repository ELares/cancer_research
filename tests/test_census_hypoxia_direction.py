"""Guards for historical hypoxia evidence and fresh, explicitly bound labels.

This analysis set out to check whether "contested" describes the literature on
the manuscript's most disputed leg. Its keyword classifier FAILED -- 34%
agreement with adjudication, and seven outright direction reversals -- so the
regex is demoted to candidate generation. The historical adjudication lacks
record identifiers and cannot be automatically applied to fresh census records.

THE FAILURE MODE IS THE VALUABLE PART and is guarded so it cannot be quietly
"fixed" by tuning: **"hypoxia-induced ferroptosis resistance" contains
"hypoxia-induced ferroptosis"**. A phrase asserting protection contains, as a
substring, the phrase asserting sensitisation. No proximity rule separates
them, and the construction concentrates in tumour biology.

The hazard specific to this file is that the adjudicated answer favours the
direction this project's simulation already assumes. So the guards check that
the adjudication is committed with reasons, that off-topic articles are
excluded rather than assigned, and that the interval is reported rather than
the point estimate alone.

OFFLINE: reads only committed artifacts.
"""
import csv
import importlib.util
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
JSON = REPO / "analysis/census-hypoxia-direction.json"
MD = REPO / "analysis/census-hypoxia-direction.md"
CSV = REPO / "analysis/hypoxia-direction-adjudication.csv"
SCRIPT = REPO / "scripts/census_hypoxia_direction.py"


@pytest.fixture(scope="module")
def d():
    return json.loads(JSON.read_text())


@pytest.fixture(scope="module")
def rows():
    return list(csv.DictReader(CSV.open(encoding="utf-8")))


@pytest.fixture
def scanner(monkeypatch):
    monkeypatch.syspath_prepend(str(REPO / "scripts"))
    spec = importlib.util.spec_from_file_location("hypoxia_unit_tests", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_every_adjudicated_row_carries_a_reason(rows):
    """A label without a reason cannot be disagreed with, which is the only
    thing that makes a single-adjudicator judgement usable by anyone else."""
    assert len(rows) >= 30
    for r in rows:
        assert r["title"].strip(), "a row has no title"
        assert r["adjudicated"] in ("protects", "sensitises", "off-topic",
                                    "ambiguous"), r
        assert len(r["reason"].strip()) > 5, (
            f"no reason given for {r['title'][:50]!r}")


def test_historical_summary_preserves_title_only_adjudication_numbers(d, rows):
    a = d["adjudication"]
    assert a["n"] == len(rows)
    from collections import Counter

    adj = Counter(r["adjudicated"] for r in rows)
    assert a["protects"] == adj["protects"]
    assert a["sensitises"] == adj["sensitises"]
    assert a["directional"] == adj["protects"] + adj["sensitises"]
    # Off-topic and ambiguous are EXCLUDED, not assigned to a side. Assigning
    # them would let a dozen irrelevant articles decide a direction.
    assert a["directional"] < a["n"], (
        "every article is directional, which means off-topic and ambiguous "
        "articles are being assigned to a side")
    assert d["adj_protects_share"] == pytest.approx(
        100 * a["protects"] / a["directional"], abs=0.1)


def test_the_classifier_failure_is_recorded_with_its_mechanism(d):
    """A retired instrument whose failure is not explained gets rebuilt."""
    a = d["adjudication"]
    assert a["regex_agreement"] < 60, (
        "the classifier now agrees with the adjudication most of the time; if "
        "it has genuinely improved, the demotion should be revisited "
        "deliberately rather than left in place")
    assert a["regex_reversed"] > 0
    md = MD.read_text()
    assert "CONTAINS" in md and "hypoxia-induced ferroptosis resistance" in md, (
        "the substring mechanism is not stated, so a future editor will try to "
        "fix the classifier by widening its patterns")
    assert "demoted to generating candidates" in md


def test_the_direction_reversal_is_reported(d):
    """The classifier and the adjudication disagree about WHICH WAY the
    literature leans, and burying that would leave the earlier number quotable."""
    assert d["regex_direction_reversed_by_adjudication"] is True
    md = MD.read_text()
    assert "adjudication reverses the recorded majority direction" in md


def test_the_interval_is_reported_and_wide_enough_to_matter(d):
    """18 directional articles cannot settle a ratio, and the favourable
    direction is the one this project's simulation assumes -- so the point
    estimate must not travel without its interval."""
    lo, hi = d["adj_ci"]
    assert lo < d["adj_protects_share"] < hi
    md = MD.read_text()
    assert f"{lo}-{hi}%" in md
    if lo < 50 < hi:
        assert "wide enough to contain an even split" in md, (
            "the interval spans an even split and the report does not say so")
        assert "do not establish their ratio in the broader" in md


def test_it_does_not_settle_the_biology(d):
    md = MD.read_text()
    assert "## What this cannot do" in md
    assert "Settle the biology" in md
    for overclaim in ("proves hypoxia", "settles the dispute",
                      "confirms the simulation", "the field agrees"):
        assert overclaim not in md.lower()


def test_the_limits_of_a_single_adjudicator_are_stated(d):
    md = MD.read_text()
    assert "## Limits of the adjudication" in md
    assert "one adjudicator" in md
    assert "hypoxia-direction-adjudication.csv" in md, (
        "the report does not point at the committed labels, so a reader cannot "
        "disagree with any of them individually")


def test_historical_report_does_not_claim_verified_record_linkage(d, scanner):
    markdown = scanner.render(d)
    assert markdown == MD.read_text()
    assert "Historical snapshot; record linkage is not verifiable" in markdown
    assert "no record identifiers or complete candidate inventory" in markdown
    assert "All 32 articles the classifier labelled" not in markdown
    assert "'contested' framing therefore reflects the literature" not in markdown
    assert "Examples of regex candidates, not adjudicated directions" in markdown


def test_assembly_drops_unbound_historical_adjudication(d, scanner):
    result = scanner.assemble(d)
    assert result["adjudication"] == {}
    assert result["adj_protects_share"] is None
    for field in ("adj_ci", "protects_ci", "dominance_threshold", "min_classified",
                  "interpretable", "verdict", "verdict_survives_interval"):
        assert field not in result
    assert result["regex_direction_reversed_by_adjudication"] is False
    assert "No adjudications are bound to this selected cohort" in scanner.render(result)


@pytest.mark.parametrize("decision", ["protects", "sensitises", "off-topic", "ambiguous"])
def test_fresh_report_uses_explicit_decisions_and_actual_denominator(scanner, decision):
    raw = {
        "total": 1,
        "counts": {"protects": 1, "sensitises": 0, "both": 0, "neither": 0},
        "sample": [], "hypoxia_descriptors": sorted(scanner.HYPOXIA),
        "cohort": {"cohort_sha256": "a" * 64},
    }
    decisions = [{"regex_label": "protects", "adjudicated": decision,
                  "reason": "Explicit reviewed direction"}]
    result = scanner.assemble(raw, adjudication_rows=decisions)
    summary = result["adjudication"]
    assert summary["mode"] == "complete-current-cohort"
    assert summary["cohort_sha256"] == "a" * 64
    assert summary["decisions"] == decisions
    assert summary["n"] == 1
    assert "regex_reversed_cancer" not in summary
    assert "adj_ci" not in result
    assert "protects_ci" not in result
    assert "verdict" not in result
    markdown = scanner.render(result)
    assert "All 1 singly classified candidates have validated labels" in markdown
    assert "18 directional articles" not in markdown
    assert "Historical snapshot" not in markdown
    assert "Wilson" not in markdown
    assert "95%" not in markdown
    assert "Both directions occur" not in markdown
    if decision in {"off-topic", "ambiguous"}:
        assert summary["directional"] == 0
        assert "There are no directional adjudications" in markdown
        assert "reverses the recorded majority" not in markdown
    else:
        assert summary["directional"] == 1
        assert result["adj_protects_share"] == (100.0 if decision == "protects" else 0.0)
        assert result["regex_direction_reversed_by_adjudication"] == (decision == "sensitises")
        assert f"**{result['adj_protects_share']}%** over 1 directional articles" in markdown
        assert "direct count of the completed candidate adjudication" in markdown


@pytest.mark.parametrize("decisions", [None, []])
def test_zero_candidates_render_without_undefined_percentages(scanner, decisions):
    raw = {
        "total": 0,
        "counts": {"protects": 0, "sensitises": 0, "both": 0, "neither": 0},
        "sample": [], "hypoxia_descriptors": sorted(scanner.HYPOXIA),
        "cohort": {"cohort_sha256": "a" * 64},
    }
    result = scanner.assemble(raw, adjudication_rows=decisions)
    markdown = scanner.render(result)
    assert "None%" not in markdown
    assert "no hypoxia/ferroptosis intersection articles" in markdown
    assert result["adj_protects_share"] is None
    assert "adj_ci" not in result
    assert "verdict" not in result


def test_complete_labels_clear_stale_historical_intervals_and_verdicts(d, scanner):
    raw = {**d, "cohort": {"cohort_sha256": "a" * 64}}
    decisions = [{"regex_label": "sensitises", "adjudicated": "protects",
                  "reason": "Cancer mentioned without an inferred context label"}]
    result = scanner.assemble(raw, adjudication_rows=decisions)
    for field in ("adj_ci", "protects_ci", "dominance_threshold", "min_classified",
                  "interpretable", "verdict", "verdict_survives_interval"):
        assert field not in result
    assert "regex_reversed_cancer" not in result["adjudication"]
    markdown = scanner.render(result)
    assert "Wilson" not in markdown
    assert "Both directions occur" not in markdown
    assert "**100.0%** over 1 directional articles" in markdown
