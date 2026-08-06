"""What the discovery rankers select, pinned against what they score.

`atlas-discovery-eval.md` observed that the harder a ranking corrects for
candidate degree, the worse it does -- ordering the methods by a hand-written
column describing what each formula does. `atlas-discovery-degree-bias.md`
measures that instead: `L` is the mean degree of a ranker's top-k over the mean
degree of the pool it drew them from, so `L = 1` is degree-neutral and
popularity is the ceiling by construction.

The measured ordering confirms the hand-written one and adds a magnitude. These
guards pin the relationship, not the numbers, and pin the limits with it --
because the interesting reading (the metric rewards not correcting) is one step
from an overclaim (hub-selection is wrong), and that step is not identifiable
from this corpus.
"""

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
RAW = REPO_ROOT / "analysis" / "atlas-discovery-degree-bias.json"
DOC = REPO_ROOT / "analysis" / "atlas-discovery-degree-bias.md"
EVAL = REPO_ROOT / "analysis" / "atlas-discovery-eval.json"


def _raw():
    return json.loads(RAW.read_text())


def test_the_measure_is_calibrated_on_its_own_control():
    """`random` must land at degree-neutral, or L means nothing.

    It is the only method whose expected L is known a priori: sampling the pool
    uniformly gives the pool mean. If it drifts off 1, the normalisation is
    wrong and every other row is uninterpretable.
    """
    L = _raw()["median_L"]
    assert "random" in L, "the control ranking is missing"
    assert 0.8 <= L["random"] <= 1.25, (
        f"random sits at {L['random']:.2f}x, not ~1x, so L is not measuring "
        "selectivity relative to the pool")


def test_popularity_is_the_ceiling_and_jaccard_is_below_neutral():
    """The two structural endpoints, which follow from the formulas."""
    L = _raw()["median_L"]
    assert L["popularity"] == max(L.values()), (
        "popularity is degree, so nothing can be more degree-selective")
    assert L["jaccard"] < 1.0, (
        "Jaccard normalises by both degrees and should over-correct past "
        "neutral; if it no longer does, the spread this analysis rests on is gone")


def test_selectivity_predicts_precision():
    """The relationship the document is about."""
    d = _raw()
    rho = d["spearman_L_vs_precision_excluding_random"]
    assert rho > 0.8, (
        f"L no longer predicts precision (rho={rho:.2f}); the document's central "
        "claim does not hold and its reading must be revisited")


def test_the_two_documents_describe_the_same_build():
    """L and precision are quoted side by side, so they must share a split.

    Quoting a selectivity computed on one graph beside a precision computed on
    another is the error that let a sibling analysis in this repo invert its own
    finding. Cheap to check, so checked.
    """
    d = _raw()
    ev = json.loads(EVAL.read_text())["headline"]
    assert d["split"] == ev["split_year"], (
        f"degree bias is computed at split {d['split']} but the precisions it "
        f"quotes come from {ev['split_year']}")
    assert d["top"] == ev["top_k"], "different k"
    for m, p in d["precision_at_k"].items():
        assert abs(p - ev["precision"][m]) < 1e-9, (
            f"the precision quoted for {m} is not the evaluation's")


def _flat() -> str:
    """The document with whitespace collapsed.

    The generator wraps prose, so a phrase can arrive split across a newline --
    `"not identifiable"` is stored as `"not\nidentifiable"`. Matching raw text
    makes a guard fail for a reason that has nothing to do with what it checks,
    and (worse) makes it PASS if a later rewrap happens to join the words.
    """
    return " ".join(DOC.read_text().split()).lower()


def test_the_document_does_not_claim_hub_selection_is_wrong():
    """One step past the finding is an unidentifiable claim."""
    txt = _flat()
    assert "does not show hub-selection is wrong" in txt, (
        "the document draws the identifiability limit nowhere")
    assert "not identifiable" in txt, (
        "the document does not say the separation is unidentifiable, so a "
        "reader may take the correlation as evidence that correcting for "
        "degree is the right thing to do")


def test_the_document_credits_the_prior_observation():
    """The direction was already in the evaluation; only the magnitude is new."""
    txt = DOC.read_text()
    assert "already suspected this" in txt, (
        "the document presents a measurement of an existing observation as a "
        "new observation")
