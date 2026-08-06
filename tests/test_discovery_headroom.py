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
    d, txt = _raw(), DOC.read_text()
    if d["any_headroom"]:
        assert "adds measurable precision" in txt
    else:
        assert "No method adds measurable precision" in txt


def test_the_identifiability_limit_is_stated():
    """One step past this result is a claim the corpus cannot support."""
    flat = " ".join(DOC.read_text().split())
    assert "not identifiable" in flat, (
        "without this, a reader takes 'the prior exhausts what is measurable' "
        "as evidence that popularity is the right prior")
    assert "statement about the" in flat and "evaluation" in flat
