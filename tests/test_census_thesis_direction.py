"""Guards for the direction check on this project's strongest leg.

The analysis asks whether the 479-article leg the manuscript leans on points
the way the manuscript needs. It came back favourable, which is the situation
in which a guard earns its place: a check on one's own thesis that returns good
news needs to be one that could have returned bad news.

THE SPECIFIC HAZARD IS FITTING. The exploit vocabulary IS the thesis
vocabulary, so widening the pattern set after seeing the split would tune the
result toward the answer this project wants, and nothing in the output would
show it. The patterns are pinned here, in a second file, so changing them takes
a deliberate edit in two places.

OFFLINE: reads only the committed artifact.
"""
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


def test_the_verdict_is_derived_and_could_have_gone_the_other_way(d):
    """A check on one's own thesis that cannot fail is not a check."""
    assert d["points_the_projects_way"] == (
        d["exploit_share_of_classified"] > 50)
    md = MD.read_text()
    if d["points_the_projects_way"]:
        assert "does point the way the manuscript needs" in md
        assert "should stop" not in md
    else:
        assert "does NOT point the project's way" in md
        assert "should stop" in md
    # The generator must contain the branch that reports the bad news, or the
    # favourable verdict is the only thing it can say.
    assert "does NOT point the project's way" in SCRIPT.read_text()


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
    assert "MAGNITUDE is not claimed" in md


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

def test_the_adjudication_corrects_the_raw_split(d):
    """The raw 89% is not the measurement.

    Adjudicating found the classifier wrong in BOTH directions: exploit papers
    labelled obstacle because they contain the obstacle phrase ("counteracts
    ferroptosis resistance"), and obstacle papers labelled exploit because they
    contain the induction phrase ("suppress erastin-induced ferroptosis"). The
    same substring failure the sibling hypoxia analysis names, twice over.
    """
    a = d["adjudication"]
    assert a, "the adjudication is missing, so the raw split stands uncorrected"
    ex, ob = a["by_label"]["exploit"], a["by_label"]["obstacle"]
    assert ob["adjudicated_exploit"] > 0, (
        "no obstacle-labelled article was found to be an exploit paper, which "
        "would mean the classifier no longer has the failure this correction "
        "exists for")
    assert ob["sampled"] == 29 or ob["sampled"] >= ob["decided"], ob
    assert 0 < ex["precision"] < 1 and 0 < ob["precision"] < 1
    # The correction must MOVE the answer, or it is decoration.
    assert abs(a["corrected_exploit_share"] - d["exploit_share_of_classified"]) > 5, (
        "the corrected share matches the raw one; either the classifier "
        "improved or the correction is not being applied")


def test_the_corrected_range_separates_direction_from_magnitude(d):
    """Direction can survive a correction that destroys the magnitude, and
    reporting only one of those is how a weakened claim keeps its old force."""
    a = d["adjudication"]
    lo, hi = a["corrected_range"]
    assert lo <= a["corrected_exploit_share"] <= hi
    assert a["direction_survives"] == (lo > 50)
    md = MD.read_text()
    assert f"{a['corrected_exploit_share']}%" in md
    assert f"{lo}-{hi}%" in md
    if a["direction_survives"]:
        assert "DIRECTION survives" in md and "MAGNITUDE does not" in md
        assert "should not be quoted" in md


def test_the_manuscript_carries_the_corrected_figure_not_the_raw_one():
    """The raw figure was in the manuscript before this correction existed."""
    txt = " ".join(MANUSCRIPT.read_text().split())
    a = json.loads(JSON.read_text())["adjudication"]
    assert f"{a['corrected_exploit_share']}%" in txt, (
        "the manuscript does not carry the corrected exploit share")
    assert "a ratio of 8.2 to 1" not in txt, (
        "the manuscript still states the uncorrected ratio, which the "
        "adjudication showed overstates the lead by roughly threefold")
    assert "errs in BOTH directions" in txt or "errs in both directions" in txt
