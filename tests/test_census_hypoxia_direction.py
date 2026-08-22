"""Guards for the hypoxia-direction measurement and the classifier it retired.

This analysis set out to check whether "contested" describes the literature on
the manuscript's most disputed leg. Its keyword classifier FAILED -- 34%
agreement with adjudication, and seven outright direction reversals -- so the
regex is demoted to candidate generation and the committed adjudication is the
measurement.

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


def test_the_adjudication_and_not_the_regex_is_the_measurement(d, rows):
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
    assert "It reverses the answer" in md


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
        assert "NOT their ratio" in md


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
