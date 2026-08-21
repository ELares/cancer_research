"""Guards for the first-literature-to-first-trial measurement.

Three design decisions carry this analysis, and each fails silently if it is
weakened. They are guarded rather than commented.

1. RIGHT-CENSORING. A mechanism with a literature and no indexed trial has no
   lag. Scoring it as a large one would make "never reached a trial"
   indistinguishable from "reached it slowly", and would quietly pull the
   median. The censored set must be derived, and censored rows must carry no
   lag.
2. TWO ARMS, KEPT SEPARATE. The MeSH arm cannot see a concept before its
   descriptor existed, so a lag measured on descriptors alone is partly a
   measurement of when NLM minted a term. The text arm must therefore read
   title and abstract ONLY -- folding MeSH into it would make the two arms
   partly one instrument and collapse the very gap being measured.
3. THE START THRESHOLD IS A COUNT, so it is NOT sample-invariant. A strided
   run reaches the threshold later on both arms and compresses every lag, which
   means the committed artifact has to come from a full pass.
4. WORD BOUNDARIES. Unbounded substrings dated `car-t` to 1947 by matching
   inside "s(car t)issue" and `electrolysis` to 1950 via "li(echt)enstein".
   A first-appearance statistic is decided by its single worst false positive
   across four million records, so it has none of the error averaging that
   makes the same vocabulary adequate for a prevalence estimate.

OFFLINE: these read only the committed artifact.
"""
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
JSON = REPO / "analysis/census-translation-lag.json"
MD = REPO / "analysis/census-translation-lag.md"
SCRIPT = REPO / "scripts/census_translation_lag.py"
# The indexed census. A committed artifact built from fewer records than this
# was strided, and a count threshold does not survive striding.
CENSUS_RECORDS = 4_403_994


@pytest.fixture(scope="module")
def d():
    if not JSON.exists():
        pytest.skip("translation-lag artifact not built")
    return json.loads(JSON.read_text())


def test_the_artifact_came_from_a_full_pass(d):
    """A count threshold is not sample-invariant.

    At stride 40 the 5-article start threshold is reached years later on both
    arms, every lag compresses, and the descriptor-delay median collapses
    toward zero -- a strided run produces a plausible table that means
    something else.
    """
    assert d["census"] == CENSUS_RECORDS, (
        f"the artifact was built from {d['census']:,} records, not the full "
        f"{CENSUS_RECORDS:,}. The start threshold counts articles, so a sample "
        "shifts every start year later and compresses every lag.")


def test_censored_mechanisms_carry_no_lag_and_are_derived(d):
    for r in d["rows"]:
        for arm in ("mesh", "text"):
            if r[f"{arm}_censored"]:
                assert r[f"{arm}_lag"] is None, (
                    f"{r['mechanism']} is censored on the {arm} arm but carries "
                    "a lag, which would score never-reached-a-trial like "
                    "reached-it-slowly")
                assert r[f"{arm}_start"] and not r[f"{arm}_first_trial"]
    assert sorted(d["censored_text"]) == sorted(
        r["mechanism"] for r in d["rows"] if r["text_censored"])


def test_the_median_excludes_censored_mechanisms(d):
    """The failure a censored row causes if it is included anyway."""
    import statistics

    lags = [r["text_lag"] for r in d["rows"] if r["text_lag"] is not None]
    assert d["median_text_lag"] == statistics.median(lags)
    assert all(not r["text_censored"] for r in d["rows"]
               if r["text_lag"] is not None)


def test_the_robust_subset_is_derived_and_leads_the_report(d):
    """Only a minority of mechanisms support the measurement, and the report
    must lead with that rather than with a median over all of them."""
    import statistics

    robust = [r for r in d["rows"] if r["text_lag"] is not None
              and (r["text_start_fragility"] or 0) <= 5]
    assert sorted(r["mechanism"] for r in robust) == sorted(d["robust"])
    assert not set(d["robust"]) & set(d["fragile"])
    if robust:
        assert d["median_robust_lag"] == statistics.median(
            [r["text_lag"] for r in robust])
    md = MD.read_text()
    assert f"Only {len(d['robust'])} of " in md
    assert "should not be quoted" in md, (
        "the all-mechanism median is presented without the warning that it "
        "averages durations with numbers set by one false positive")


def test_fragile_lags_are_called_upper_bounds(d):
    """The direction, which I first got backwards.

    A false positive sets the START earlier than the truth, so the computed
    lag (trial - start) comes out LONGER than the real duration. The reported
    number is therefore an upper bound. Reasoning from "the error makes the
    start too early" to "the lag is a lower bound" skips the subtraction.
    """
    md = MD.read_text()
    if not d["fragile"]:
        pytest.skip("no fragile rows")
    assert "UPPER BOUNDS" in md
    assert "lags are LOWER BOUNDS" not in md
    assert "sets the start EARLIER than the truth" in md


def test_the_text_arm_is_word_bounded(d):
    """The fix that moved car-t by 53 years."""
    src = SCRIPT.read_text()
    assert r"\b{re.escape(t.lower())}\b" in src or "rf\"\\b" in src, (
        "the text arm no longer word-bounds its terms, so `car t` matches "
        "inside 'scar tissue' and the start years become substring artifacts")
    md = MD.read_text()
    assert "s(car t)issue" in md and "li(echt)enstein" in md, (
        "the report no longer names the substring failures it was corrected "
        "for, so a future editor may 'simplify' the patterns back")


def test_the_remaining_polysemy_failure_is_named(d):
    """Word boundaries do not fix a real word used for another subject, and a
    report that only described the fixed defects would imply the rows are now
    sound."""
    md = MD.read_text()
    assert "POLYSEMY" in md
    assert "cuproptosis" in md and "copper ionophore" in md
    assert "electrolysis" in md
    assert "no error averaging" in md, (
        "the transferable point -- a minimum has no error averaging while a "
        "prevalence estimate does -- is missing")


def test_every_lag_recomputes_from_its_two_years(d):
    for r in d["rows"]:
        for arm in ("mesh", "text"):
            s, t, lag = r[f"{arm}_start"], r[f"{arm}_first_trial"], r[f"{arm}_lag"]
            if lag is not None:
                assert s is not None and t is not None
                assert lag == t - s
                assert lag >= 0, (
                    f"{r['mechanism']} has a trial before its literature starts "
                    f"on the {arm} arm, which means the start threshold is "
                    "hiding early articles")


def test_the_text_arm_reads_no_mesh(d):
    """The gap between the arms IS the measurement.

    If the text arm folded MeSH in, the two would share an instrument and the
    descriptor-delay column would understate itself.
    """
    src = SCRIPT.read_text()
    assert "Title and abstract ONLY" in src
    body = src[src.index("blob = f\""):src.index("return {", src.index("blob = f\""))]
    assert "mesh" not in body.lower(), (
        "the text arm's matching block references mesh, so the two arms are "
        "not independent and the descriptor delay is not what it claims")


def test_the_descriptor_delay_is_reported_per_mechanism_not_only_as_a_median(d):
    """A median delay invites applying one correction to every row. The delay
    is a property of each descriptor's own history and is not uniform."""
    md = MD.read_text()
    # Named `arm gap`, not `descriptor delay`. The first name asserted a CAUSE
    # -- when NLM minted a term -- and the sign refuted it: a negative value
    # means the descriptor is older or broader than the word, which is a
    # different effect entirely.
    assert "arm gap" in md
    assert "NOT a measure of when a descriptor was minted" in md
    assert "read the per-mechanism column, not the median" in md
    have = [r for r in d["rows"] if r["arm_start_gap"] is not None]
    assert have, "no mechanism has both arms, so no delay is measurable"
    spread = {r["arm_start_gap"] for r in have}
    if len(spread) > 1:
        assert d["median_arm_start_gap"] is not None


def test_text_only_mechanisms_are_marked_on_their_rows(d):
    """TTFields, bioelectric and CAP are reported as not measurable everywhere
    else here. They ARE measurable on one arm, and a row that did not say so
    would look like an ordinary two-arm result with empty cells."""
    for r in d["rows"]:
        if r["text_measurable"] and not r["mesh_measurable"]:
            assert r["mesh_start"] is None and r["mesh_lag"] is None
    md = MD.read_text()
    if d["text_only"]:
        assert "measurable on the TEXT arm only" in md
        for m in d["text_only"]:
            assert f"`{m}`" in md


def test_it_does_not_claim_a_lag_measures_success(d):
    """An indexed trial says a trial happened. A short lag can mean a low
    barrier to a first-in-human study rather than a strong result."""
    md = MD.read_text()
    assert "## What a lag is not" in md
    assert "not evidence that a mechanism translated WELL" in md
    assert "82.5%" in md, (
        "the report does not carry the text arm's measured mechanism "
        "precision, which is what bounds an early stray match starting the "
        "clock early")
    for overclaim in ("translated fastest", "most successful mechanism",
                      "proves that"):
        assert overclaim not in md.lower()


# --- the lesson has to reach the sites that inherit it --------------------

def test_the_manuscript_carries_the_statistic_limit():
    """The finding constrains other claims, so it belongs where they are made.

    Only 4 of 25 mechanisms survived this measurement, and the reason is a
    property of the INSTRUMENT rather than of translation: a minimum has no
    error averaging. Left only in this analysis, a reader has no way to know
    which of the manuscript's other numbers inherit the limit.
    """
    md = (REPO / "article/drafts/v1.md").read_text()
    d = json.loads(JSON.read_text())
    n_rob = len(d["robust"])
    n_all = n_rob + len(d["fragile"])
    assert f"{n_rob} of {n_all} mechanisms" in md, (
        "the manuscript does not state how many mechanisms survived, which is "
        "the measurement that makes the limit concrete rather than cautionary")
    assert "census-translation-lag.md" in md
    for phrase in ("averages its errors", "MINIMUM"):
        assert phrase in md, (
            f"the manuscript states the limit without {phrase!r}, so a reader "
            "cannot tell which statistics it applies to")


def test_the_replication_analysis_carries_the_same_limit():
    """It dates 2.37 million pairs by first assertion from an ~79.6-F1
    extractor, so it inherits this defect exactly.

    A newly found defect class has to reach the analyses that share it, or the
    lesson stays local to where it happened to be discovered.
    """
    repl = REPO / "analysis/atlas-replication.md"
    if not repl.exists():
        pytest.skip("replication analysis not built")
    md = repl.read_text()
    assert "no error averaging" in md, (
        "atlas-replication dates every pair by a MINIMUM over an imperfect "
        "extractor and does not say so")
    # ANCHORED TO THE SENTENCE, not to the words. `biased` and `DOWN` both
    # already appear in an unrelated paragraph about MeSH indexing lag, so a
    # two-word check passes on a document that never mentions this defect --
    # the substring trap, in the guard written to propagate a finding.
    assert "falls OUTSIDE the window and the pair scores unreplicated" in md, (
        "the limit is stated without the mechanism that gives it a direction: "
        "a false-early date starts the equal window early, so a genuine second "
        "paper can arrive after it closes")
    assert "replication rate is therefore biased DOWN" in md, (
        "the direction is missing, which is the part a reader needs to know "
        "whether the rate is optimistic or pessimistic")
