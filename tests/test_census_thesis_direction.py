"""Keep historical evidence and candidate-level interpretation distinct.

THE SPECIFIC HAZARD IS FITTING. The exploit vocabulary IS the thesis
vocabulary, so widening the pattern set after seeing the split would tune the
result toward the answer this project wants, and nothing in the output would
show it. The patterns are pinned here, in a second file, so changing them takes
a deliberate edit in two places.

Offline fixtures exercise the renderer and preserve the published calculation.
CLI identity/coverage behavior is tested separately with isolated gzip inputs.
"""
import importlib.util
import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
JSON = REPO / "analysis/census-thesis-direction.json"
MD = REPO / "analysis/census-thesis-direction.md"
MANUSCRIPT = REPO / "article/drafts/v1.md"
SCRIPT = REPO / "scripts/census_thesis_direction.py"
# Pinned INDEPENDENTLY of the generator. Counting the alternations rather than
# copying the patterns: the point is that the vocabulary cannot quietly grow.
EXPLOIT_TERMS = 10
OBSTACLE_TERMS = 7


@pytest.fixture(scope="module")
def d():
    return json.loads(JSON.read_text())


@pytest.fixture
def module(monkeypatch):
    monkeypatch.syspath_prepend(str(REPO / "scripts"))
    spec = importlib.util.spec_from_file_location("thesis_unit_test", SCRIPT)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def test_the_pattern_sets_have_not_grown(d):
    """Widening the exploit set after seeing the split would fit the answer.

    Counted rather than compared string-for-string, so a rewording that keeps
    the same coverage passes while an addition fails -- which is the change
    that would move the result.
    """
    src = SCRIPT.read_text()

    def top_level_terms(name: str) -> int:
        body = re.search(rf"{name} = re\.compile\(\s*(.*?)\)\n", src, re.S).group(1)
        # Inner non-capturing groups carry their own alternations -- the
        # `(?:drug |chemo)` inside one exploit term is not a separate term --
        # so they are stripped before counting. A raw pipe count would move
        # when a term was rephrased, which is not the change worth blocking.
        body = re.sub(r"\(\?:[^)]*\)", "X", body)
        return body.count("|") + 1

    assert top_level_terms("EXPLOIT") == EXPLOIT_TERMS, (
        f"the exploit vocabulary now has {top_level_terms('EXPLOIT')} terms "
        f"where this guard pins {EXPLOIT_TERMS}. That vocabulary IS the thesis "
        "vocabulary; growing it after seeing the split fits the answer, so the "
        "change must be deliberate in both files.")
    assert top_level_terms("OBSTACLE") == OBSTACLE_TERMS
    assert "FIXED BEFORE THE RESULT WAS READ" in src, (
        "the no-tuning commitment is no longer stated in the generator")


def test_the_ratio_is_computed_over_singly_classified_articles_only(d):
    """An article carrying BOTH framings is not evidence for either.

    Folding those into the larger side would let ambiguous cases inflate
    whichever direction already leads -- and here that is the direction this
    project wants.
    """
    c = d["counts"]
    assert d["classified"] == c["exploit"] + c["obstacle"]
    assert d["unclassified"] == c["neither"] + c["both"], (
        "articles carrying both framings are being counted as classified")
    assert d["exploit_share_of_classified"] == pytest.approx(
        100 * c["exploit"] / d["classified"], abs=0.1)
    assert sum(c.values()) == d["total"]

    # INTERNAL CONSISTENCY IS NOT ENOUGH, and a mutation proved it: moving the
    # `both` count into `exploit` and zeroing it keeps every equation above
    # true while inflating the direction this project wants. The artifact has
    # no independent record of what `both` should be, so the check has to come
    # from outside it -- the SOURCE must not fold ambiguous articles into
    # either side, and a leg this size cannot plausibly have zero overlap.
    src = SCRIPT.read_text()
    import re as _re

    # Anchored to the END of the expression, not to a prefix: adding `+
    # c["both"]` leaves the prefix intact, so a substring check passes on the
    # exact edit it exists to block.
    assert _re.search(r'single = c\["exploit"\] \+ c\["obstacle"\]\s*\n', src), (
        "the ratio is no longer computed over singly classified articles")
    assert 'c["neither"] + c["both"]' in src, (
        "articles carrying both framings are no longer counted as unclassified")
    assert c["both"] > 0, (
        f"{d['total']} articles and zero carry both framings, which two "
        "overlapping vocabularies over a literature this size do not produce. "
        "The ambiguous cases have most likely been folded into one side.")


def test_historical_verdict_does_not_become_a_current_claim(d, module):
    """Preserve the snapshot without presenting its boolean as new evidence."""
    assert d["points_the_projects_way"] is True
    md = MD.read_text()
    assert md == module.render(d)
    assert "historical snapshot" in md
    assert "does point the way the manuscript needs" not in md
    assert "should stop" not in md
    fresh = module.assemble(d)
    assert not fresh.get("adjudication")
    assert "points_the_projects_way" not in fresh


def test_the_recall_limit_is_shown_rather_than_asserted(d):
    """A large unclassified share with no examples reads as a formality.

    The sample is what makes it concrete: the first unmatched title is the
    thesis direction stated in words the exploit set does not contain.
    """
    assert d["unclassified_share"] > 0
    assert d["unclassified_sample"], "no unclassified examples are shown"
    md = MD.read_text()
    assert f"{d['unclassified_share']}%" in md
    for s in d["unclassified_sample"]:
        assert s["title"][:40] in md, "an unclassified example is not shown"
    assert "Recall remains unmeasured" in md
    assert "Their biological direction" in md
    assert "9-to-1" not in md
    assert "The first of those is the thesis direction" not in md


def test_it_does_not_read_attention_as_endorsement(d):
    """A field framing something as exploitable is evidence about what it is
    TRYING. This is the caveat most likely to be dropped, because the result
    is favourable."""
    md = MD.read_text()
    assert "## What this does not establish" in md
    assert "what the field is TRYING" in md
    for overclaim in ("confirms the thesis", "supports the hypothesis",
                      "validates the project", "evidence that ferroptosis works"):
        assert overclaim not in md.lower()


# --- the classifier errs in both directions, and the raw split overstates ----

def test_historical_adjudication_counts_are_preserved(d):
    """The frozen correction stays auditable despite its linkage limitations."""
    a = d["adjudication"]
    assert a, "the adjudication is missing, so the raw split stands uncorrected"
    ex, ob = a["by_label"]["exploit"], a["by_label"]["obstacle"]
    assert ob["adjudicated_exploit"] > 0, (
        "no obstacle-labelled article was found to be an exploit paper, which "
        "would mean the classifier no longer has the failure this correction "
        "exists for")
    assert ob["sampled"] == 29 and ex["sampled"] == 20
    assert a["rows"] == 49
    assert 0 < ex["precision"] < 1 and 0 < ob["precision"] < 1
    # The correction must MOVE the answer, or it is decoration.
    assert abs(a["corrected_exploit_share"] - d["exploit_share_of_classified"]) > 5, (
        "the corrected share matches the raw one; either the classifier "
        "improved or the correction is not being applied")


def test_historical_range_is_explicitly_conditional(d):
    """A sampled precision interval cannot cover identity and sampling bias."""
    a = d["adjudication"]
    lo, hi = a["corrected_range"]
    assert lo <= a["corrected_exploit_share"] <= hi
    assert a["direction_survives"] == (lo > 50)
    md = MD.read_text()
    assert f"{a['corrected_exploit_share']}%" in md
    assert f"{lo}-{hi}%" in md
    assert "Historical conditional correction" in md
    assert "representativeness is unverified" in md
    assert "verified article-identity link or complete coverage" in md
    assert "every plausible correction" not in md
    assert "roughly 3 to 1" not in md


def test_the_manuscript_keeps_the_historical_figure_with_its_limitations():
    """A retained historical result must not regain an unconditional verdict."""
    txt = " ".join(MANUSCRIPT.read_text().split())
    a = json.loads(JSON.read_text())["adjudication"]
    assert f"{a['corrected_exploit_share']}%" in txt, (
        "the manuscript does not carry the corrected exploit share")
    assert "a ratio of 8.2 to 1" not in txt, (
        "the manuscript still states the uncorrected ratio, which the "
        "adjudication showed overstates the lead by roughly threefold")
    paragraph = next(line for line in MANUSCRIPT.read_text().splitlines()
                     if "analysis/census-thesis-direction.md" in line)
    assert "historical" in paragraph.lower()
    assert "representativeness" in paragraph.lower()
    assert "direction survives every plausible correction" not in paragraph


@pytest.mark.parametrize("counts", [
    {"exploit": 0, "obstacle": 0, "both": 0, "neither": 0},
    {"exploit": 0, "obstacle": 0, "both": 1, "neither": 2},
    {"exploit": 1, "obstacle": 0, "both": 0, "neither": 0},
    {"exploit": 0, "obstacle": 1, "both": 0, "neither": 0},
    {"exploit": 1, "obstacle": 1, "both": 0, "neither": 0},
])
def test_unadjudicated_candidates_cannot_assert_direction(module, counts):
    data = module.assemble({"total": sum(counts.values()), "counts": counts,
                            "unclassified_sample": []})
    text = module.render(data)
    assert not data.get("adjudication")
    assert "points_the_projects_way" not in data
    assert "No validated adjudication" in text
    assert "do not establish a biological direction" in text
    assert "None" not in text and "9-to-1" not in text
    assert "should stop" not in text


@pytest.mark.parametrize("decisions,direction,share", [
    (["exploit", "exploit", "exploit", "obstacle", "ambiguous"], "exploit", 75.0),
    (["exploit", "obstacle", "obstacle", "obstacle", "ambiguous"], "obstacle", 25.0),
    (["exploit", "obstacle", "ambiguous"], "tie", 50.0),
    (["exploit"], "exploit", 100.0),
    (["obstacle"], "obstacle", 0.0),
    (["ambiguous"], None, None),
    ([], None, None),
])
def test_complete_decisions_bound_direction_to_decided_candidates(
        module, decisions, direction, share):
    rows = [{"adjudicated": decision} for decision in decisions]
    raw = {"total": len(rows) + 1,
           "counts": {"exploit": len(rows), "obstacle": 0,
                      "both": 0, "neither": 1},
           "unclassified_sample": [], "cohort": {"cohort_sha256": "fixture"}}
    data = module.assemble(raw, adjudication_rows=rows)
    adjudication = data["adjudication"]
    assert adjudication["direction"] == direction
    assert adjudication["exploit_share_of_decided"] == share
    assert adjudication["decisions"] == rows
    assert adjudication["cohort_sha256"] == "fixture"
    assert "corrected_range" not in adjudication
    text = module.render(data)
    assert "historical" not in text.lower()
    assert "None" not in text
    if direction is None:
        assert "no direction comparison is available" in text
    elif direction == "tie":
        assert "counts are equal" in text
    else:
        assert f"adjudicated {direction} framing leads within this candidate set" in text
    if share is not None:
        assert f"{share}%" in text
        assert "without extrapolation to unmatched articles" in text
