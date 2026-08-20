"""The manuscript's claim about WHICH WAY the keyword arm errs must be derived.

Section 3.6 and Section 3.13 both state a direction: the keyword tagger's
45.6% co-occurrence rate is more likely INFLATED than deflated, so the
instrument step of the decomposition is the generous reading rather than the
conservative one.

THAT SENTENCE WAS WRONG WHEN FIRST WRITTEN, in the opposite direction, and it
was wrong for the reason this repo keeps rediscovering: it reasoned from ONE of
two measurements that sit side by side. Mechanism RECALL (~90%) says the tagger
misses tags, which deflates a co-occurrence rate. Mechanism PRECISION (82.5%
strict, measured separately over 325 adjudicated pairs) says it adds tags,
which inflates one. Reading only the first gives "under-detects, so the gap is
a lower bound" -- the reverse of what the pair supports.

The asymmetry is what decides it, and it is a property of the STATISTIC rather
than of either measurement: a rate requiring an article to carry TWO tags is
created outright by one spurious tag, while a missed tag only removes one. So
precision error is amplified and recall error is not.

These guards derive the direction from the two committed measurements. If
either is remeasured such that the direction flips, the manuscript sentence
fails rather than quietly becoming false.
"""
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
MANUSCRIPT = REPO / "article/drafts/v1.md"
PRECISION_MD = REPO / "analysis/mechanism-precision-report.md"
RECALL_JSON = REPO / "analysis/mechanism-recall.json"


def _strict_precision() -> float:
    """Overall strict mechanism-tag precision, read from the committed report."""
    m = re.search(r"\*\*Overall strict precision: \d+/\d+ = ([\d.]+)%\*\*",
                  PRECISION_MD.read_text())
    assert m, "the precision report no longer states an overall strict figure"
    return float(m.group(1))


def _volume_weighted_recall() -> float:
    d = json.loads(RECALL_JSON.read_text())
    agg = d["aggregates"]
    # The real key. A first version guessed three plausible names, matched
    # none, and SKIPPED -- so the guard reported green while checking nothing,
    # which is the shape it exists to prevent one level up.
    key = "volume_weighted_recall_leakage_free"
    assert key in agg, (
        f"mechanism-recall.json no longer carries {key!r}; aggregates are "
        f"{sorted(agg)}. Fix the key rather than letting this skip, because a "
        "skipped direction check is indistinguishable from a passing one.")
    v = agg[key]
    return v * 100 if v <= 1 else v


def test_the_two_measurements_are_both_present_and_disagree_in_sign():
    """The whole point: one measure alone gives the wrong answer.

    If precision ever exceeds recall, the tagger's dominant error becomes
    omission rather than over-tagging, and the manuscript's direction would
    have to be rewritten rather than kept.
    """
    p, r = _strict_precision(), _volume_weighted_recall()
    assert 0 < p < 100 and 0 < r < 100
    assert p < r, (
        f"mechanism precision ({p}%) is no longer below recall ({r}%), so the "
        "tagger's dominant error is no longer over-tagging and the "
        "manuscript's 'more likely inflated' direction does not follow. "
        "Rewrite the sentence rather than leaving it standing.")


def test_the_manuscript_states_the_direction_the_measurements_support():
    txt = " ".join(MANUSCRIPT.read_text().split())
    p = _strict_precision()
    assert f"{p}%" in txt, (
        f"the manuscript does not quote the measured mechanism precision "
        f"({p}%), which is the measurement its direction rests on")
    assert "more likely inflated than deflated" in txt, (
        "the manuscript no longer states which way the keyword arm errs, so a "
        "reader cannot tell whether the instrument gap is the generous or the "
        "conservative reading")
    # The banned sentence is the one that was actually wrong, kept banned by
    # its reasoning rather than its wording: 'under-detects' is a claim from
    # recall alone.
    assert "the keyword arm under-detects" not in txt, (
        "the manuscript has reverted to reasoning from recall alone, which "
        "gives the opposite direction to what the pair of measurements "
        "supports")


def test_the_asymmetry_argument_is_stated_not_assumed():
    """Without it the direction looks like a preference between two numbers.

    The argument is that the statistic requires TWO tags, so one spurious tag
    creates a co-occurrence while one missed tag only removes one. A reader
    who does not see that has no way to check the conclusion.
    """
    txt = " ".join(MANUSCRIPT.read_text().split())
    # ANCHORED TO THE ARGUMENT, not to a word. A first version asserted
    # `"amplified" in txt`, which the manuscript satisfies three times over --
    # "amplified by checkpoint blockade" among them -- so deleting the actual
    # reasoning left the guard green. The substring trap, in the guard written
    # to stop a reasoning error.
    for phrase in ("precision error is amplified and the recall error is not",
                   "amplifies a precision error while a recall error only "
                   "removes one"):
        assert phrase in txt, (
            f"the asymmetry argument is missing: {phrase!r}. Without it the "
            "direction reads as a preference between two numbers rather than "
            "a property of a statistic that requires two tags.")
    assert "TWO tags" in txt or "two tags" in txt


def test_both_sites_that_carry_the_direction_agree():
    """It appears in Section 3.6 and again in Section 3.13.

    A corrected claim's summary phrase travels further than its analysis, and
    this one was wrong in both places before it was right in either.
    """
    txt = " ".join(MANUSCRIPT.read_text().split())
    assert txt.count("generous") >= 2, (
        "only one site states the direction; the decomposition in Section 3.13 "
        "and the limitation in Section 3.6 must agree, or a reader meeting one "
        "of them takes away the wrong reading")
