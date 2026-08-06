"""Seed-specific signal adds nothing over a candidate-only prior, and that
result must stay honest about its own decision rule.

The discovery leaderboard compares rankers in ISOLATION, which cannot answer
whether a seed-aware method is COMPLEMENTARY to the candidate prior -- a signal
can be weaker alone and still add on top. This measures that directly.

The guards here exist because the first version of the analysis had two ways to
overclaim, and one of them fired: a weight sweep whose verdict read the SIGN of
each cell rather than testing it, which credited +0.03 hits out of 20 as a lead
worth pursuing. Every positive cell is now put through the same paired bootstrap
as the headline, and the prose is conditional on the outcome.
"""

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
RAW = REPO_ROOT / "analysis" / "atlas-discovery-headroom.json"
DOC = REPO_ROOT / "analysis" / "atlas-discovery-headroom.md"
DUMP_SIDE = REPO_ROOT / "analysis" / "atlas-discovery-candidates-seeds.json"
EVAL = REPO_ROOT / "analysis" / "atlas-discovery-eval.json"


def _raw():
    return json.loads(RAW.read_text())


def test_the_benchmark_reproduces_the_popularity_ranking():
    """The candidate-only prior IS popularity, so it must score like it.

    `deg_c` does not vary with the seed, so ranking a pool by it is exactly the
    published `popularity` ranking. If the benchmark's precision drifts from
    that number, the blend is being compared against something else.
    """
    d = _raw()
    ev = json.loads(EVAL.read_text())["headline"]
    assert abs(d["benchmark_precision"] - ev["precision"]["popularity"]) < 1e-9, (
        f"benchmark scores {d['benchmark_precision']:.4f}, the published "
        f"popularity ranking {ev['precision']['popularity']:.4f}")


def test_it_describes_the_same_build_as_the_evaluation():
    d = _raw()
    ev = json.loads(EVAL.read_text())["headline"]
    assert d.get("pairs_before") == ev["pairs_before"], (
        "the headroom test and the precisions it compares against come from "
        "different graph builds")


def test_every_positive_sweep_cell_was_tested():
    """A positive mean with no interval is a sign, not a finding."""
    d = _raw()
    positives = [(w, m) for w, row in d["sweep"].items()
                 for m, v in row.items() if v > 0]
    tested = {(str(t["weight"]), t["method"]) for t in d["sweep_positive_cells_tested"]}
    missing = [c for c in positives if c not in tested]
    assert not missing, (
        f"these sweep cells are positive but were never put through the paired "
        f"bootstrap: {missing}. Reporting them without one is how a +0.03 "
        "difference becomes a lead.")


def test_the_verdict_follows_the_decision_rule_not_the_signs():
    d, txt = _raw(), DOC.read_text()
    survivors = d["sweep_survivors"]
    if survivors:
        assert "survive the decision rule" in txt
    else:
        assert "No cell survives the decision rule" in txt, (
            "no cell passed the test, but the document does not say so")
        # and it must not imply the positives mean something
        assert "worth pursuing" not in txt, (
            "the retracted framing that treated positive-but-untested cells as "
            "leads is back")


def test_the_headline_matches_the_measured_headroom():
    """Both directions, and neither satisfiable by the other's wording.

    The positive branch used to assert "adds measurable precision", which is a
    SUBSTRING of the negative headline "No method adds measurable precision" --
    so the branch that exists to catch "the result flipped, the prose did not"
    could not catch it. Proved by setting any_headroom truthy against the
    negative document: it passed.
    """
    d, txt = _raw(), DOC.read_text()
    if d["any_headroom"]:
        assert "No method adds measurable precision" not in txt, (
            "a method now has headroom but the document still says none does")
        assert "adds measurable precision" in txt
    else:
        assert "No method adds measurable precision" in txt


def test_the_identifiability_limit_is_stated():
    """One step past this result is a claim the corpus cannot support."""
    txt = DOC.read_text()
    assert "## What this cannot settle" in txt, "the limits section is gone"
    flat = " ".join(txt[txt.index("## What this cannot settle"):].split())
    # Scoped to the limits SECTION. The previous version searched the whole
    # document for "statement about the" and "evaluation", both of which appear
    # in earlier sections another guard already pins -- so deleting this entire
    # section left it green.
    assert "not identifiable" in flat, (
        "without this, a reader takes 'the prior exhausts what is measurable' "
        "as evidence that popularity is the right prior")
    assert "ONE combination family" in flat, (
        "the result is stated over all combinations when one family was tested")


def test_the_findings_page_no_longer_asserts_the_second_half():
    """'A good candidate generator and a bad ranker' was half-wrong.

    The first half stands -- the candidate set beats random severalfold. The
    second does not follow, and this is the page where that claim travels
    furthest, so it must carry the qualification rather than the slogan.
    """
    import re

    page = REPO_ROOT / "analysis" / "census-findings.md"
    assert page.exists(), "the findings page is missing"
    txt = page.read_text()
    anchor = "Literature-based discovery does not work"
    assert anchor in txt, (
        "the discovery finding is gone from the findings page, or its heading "
        "was reworded -- this guard used to return silently in that case")
    start = txt.index(anchor)
    section = txt[start:txt.index("*Source:*", start)]
    d = _raw()
    if d["any_headroom"]:
        return                      # the claim would be defensible again
    # Punctuation-independent. The previous version checked for "a bad ranker."
    # with a trailing period, so "a bad ranker, plainly." re-asserted the exact
    # retracted claim and passed.
    assert not re.search(r"\ba bad ranker\b", section), (
        "the findings page still asserts the ranker is bad; the headroom test "
        "shows a degree-correcting ranker and a bad one are indistinguishable "
        "on this metric")
    assert "does not follow" in section


def test_the_mission_statement_carries_the_same_correction():
    """MISSION.md quotes this finding, and its wording outlives the analysis.

    It ended with the same "a good generator, a bad ranker" summary the findings
    page carried. The first half stands; the second does not follow. Guarded in
    the same conditional way, so if a method ever shows headroom the original
    summary becomes defensible again rather than being frozen out.
    """
    import re

    mission = REPO_ROOT / "MISSION.md"
    assert mission.exists()
    txt = mission.read_text()
    anchor = "shipped ABC ranking scores"
    assert anchor in txt, (
        "the discovery finding is gone from MISSION.md, or its wording changed "
        "-- this guard must not pass silently in that case")
    section = txt[txt.index(anchor):]
    section = section[:section.index("atlas-discovery-headroom.md") + 40] \
        if "atlas-discovery-headroom.md" in section else section[:2500]
    d = _raw()
    if d["any_headroom"]:
        return
    assert not re.search(r"^\s*>?\s*A good generator, a bad ranker", section, re.M), (
        "MISSION.md still ends on the retracted summary")
    assert "does not follow" in section
    assert "not identifiable" in section, (
        "the mission statement asserts the metric's limit without the stronger "
        "limit that hub-selection's correctness is unidentifiable here")
