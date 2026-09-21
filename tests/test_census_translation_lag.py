"""Guards for the first-literature-to-first-trial measurement.

Three design decisions carry this analysis, and each fails silently if it is
weakened. They are guarded rather than commented.

1. RIGHT-CENSORING. A mechanism with a literature and no indexed trial has no
   lag. Scoring it as a large one would make "never reached a trial"
   indistinguishable from "reached it slowly", and would quietly pull the
   median. The censored set must be derived, and censored rows must carry no
   lag.
2. TWO ARMS, KEPT SEPARATE. Descriptor and keyword coverage can differ, and
   first-match ordering alone does not identify why. The text arm must read
   title and abstract ONLY -- folding MeSH into it would make the two arms
   partly one instrument and collapse the very gap being measured.
3. THE START THRESHOLD IS A COUNT, so it is NOT sample-invariant. A strided
   run reaches the threshold later on both arms and compresses every lag, which
   means the committed artifact has to come from a full pass.
4. WORD BOUNDARIES. Unbounded substrings dated `car-t` to 1947 by matching
   inside "s(car t)issue" and `electrolysis` to 1950 via "li(echt)enstein".
   A first-appearance statistic depends on its earliest match, so a false
   positive can move it earlier and a missed match can move it later.
   Aggregate prevalence accuracy also needs both precision and recall.

OFFLINE: these read only the committed artifact.
"""
import json
import gzip
from pathlib import Path
import sys

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
import census_translation_lag as scanner
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
              and r["text_start_fragility"] is not None
              and r["text_start_fragility"] <= 5]
    assert sorted(r["mechanism"] for r in robust) == sorted(d["robust"])
    assert not set(d["robust"]) & set(d["fragile"])
    if robust:
        assert d["median_robust_lag"] == statistics.median(
            [r["text_lag"] for r in robust])
    md = MD.read_text()
    assert f"Only {len(d['robust'])} of " in md
    assert "should not be quoted" in md, (
        "the all-mechanism median is presented without the warning that it "
        "includes durations sensitive to the article-count threshold")


def test_threshold_sensitive_lags_are_not_claimed_to_bound_true_duration(d):
    """Neither first-match accuracy follows from an article-count threshold."""
    if not d["fragile"]:
        pytest.skip("no fragile rows")
    for md in (MD.read_text(), scanner.render(scanner.assemble(d))):
        assert "threshold-sensitive, not established bounds" in md
        assert "without checking both the first literature and first trial matches" in md
        assert "82.5% precision alone does not establish the accuracy of prevalence estimates" in md
        assert "sensitivity to the article-count threshold, not whether individual matches are correct" in md
        assert "Both stable and threshold-sensitive lags need source-level verification" in md
        for overclaim in (
            "upper bound on the true duration",
            "off by a predictable fraction",
            "supports one statistic and not the other",
            "wrong in a way the fragility column flags",
            "which is what the fragility column measures",
        ):
            assert overclaim not in md


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
        "the transferable point -- a minimum has no error averaging -- is missing")


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


def test_the_arm_gap_reports_ordering_without_attributing_a_cause(d):
    """Date differences must not diagnose descriptor history or breadth."""
    reports = [MD.read_text(), scanner.render(scanner.assemble(d))]
    for mesh_year, text_year in ((2000, 2020), (2020, 2000)):
        sparse = scanner.assemble({
            "census": 2, "first_min": 1, "stable_min": 5, "last_full_year": 2025,
            "mesh_measurable": ["example"], "text_measurable": ["example"],
            "arms": {"mesh": {"example": {str(mesh_year): 1}},
                     "text": {"example": {str(text_year): 1}}},
            "trials": {"mesh": {}, "text": {}},
        })
        assert sparse["rows"][0]["arm_start_gap"] == mesh_year - text_year
        reports.append(scanner.render(sparse))
    for md in reports:
        assert "arm gap" in md
        assert "NOT a measure of when a descriptor was minted" in md
        assert "read the per-mechanism column, not the median" in md
        assert "The sign alone does not identify the cause" in md
        assert "source-level checks of the matches and descriptor history" in md
        for overclaim in (
            "measurement of MeSH's own history",
            "which is the descriptor-introduction effect",
            "NEGATIVE gap means the DESCRIPTOR is older or broader",
            "For these the descriptor is the broader instrument",
        ):
            assert overclaim not in md
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
        "precision and its limits for interpreting aggregate counts")
    for overclaim in ("translated fastest", "most successful mechanism",
                      "proves that"):
        assert overclaim not in md.lower()


# --- the lesson has to reach the sites that inherit it --------------------

def test_the_manuscript_carries_the_statistic_limit():
    """Keep the measured stability result separate from verified date accuracy."""
    md = (REPO / "article/drafts/v1.md").read_text()
    d = json.loads(JSON.read_text())
    n_rob = len(d["robust"])
    n_all = n_rob + len(d["fragile"])
    paragraph = next(p for p in md.split("\n\n") if "analysis/census-translation-lag.md" in p)
    assert f"{n_rob} of {n_all} mechanisms met the declared stability screen" in paragraph
    for phrase in (
        "counts depend on both precision and recall",
        "82.5% precision alone does not establish the accuracy of prevalence estimates",
        "not the accuracy of the earliest literature or trial matches",
        "not established bounds on true durations",
        "stable lags still need source-level verification",
    ):
        assert phrase in paragraph, f"the manuscript dropped the qualification: {phrase}"
    assert "survived it" not in paragraph


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


def test_shared_reader_preserves_sorted_shard_stride_and_arm_separation(tmp_path, monkeypatch):
    for index in (4, 1, 3, 0, 2):
        with gzip.open(tmp_path / f"{index:03}.jsonl.gz", "wt") as stream:
            json.dump({"year": 2010 + index, "title": "immunotherapy",
                       "mesh": ["Autophagy"], "pub_types": []}, stream)
            stream.write("\n")
    monkeypatch.setattr(scanner, "RECORDS", tmp_path)
    result = scanner.scan(2)
    assert result["census"] == 3
    assert result["arms"]["text"]["immunotherapy"] == {
        "2010": 1, "2012": 1, "2014": 1}
    assert "immunotherapy" not in result["arms"]["mesh"]


def test_valid_census_without_matches_reports_unavailable_lags(tmp_path, monkeypatch):
    with gzip.open(tmp_path / "records.jsonl.gz", "wt") as stream:
        stream.write(json.dumps({"year": 2020, "title": "unrelated observation"}) + "\n")
    monkeypatch.setattr(scanner, "RECORDS", tmp_path)
    result = scanner.assemble(scanner.scan())
    assert result["census"] == 1
    assert result["rows"] == []
    assert result["median_text_lag"] is None
    markdown = scanner.render(result)
    assert "No mechanisms matched" in markdown
    assert "None" not in markdown


def test_missing_stability_threshold_is_not_evidence_of_robustness():
    result = scanner.assemble({
        "census": 1, "first_min": 1, "stable_min": 5, "last_full_year": 2025,
        "mesh_measurable": [], "text_measurable": ["example"],
        "arms": {"mesh": {}, "text": {"example": {"2020": 1}}},
        "trials": {"mesh": {}, "text": {"example": {"2020": 1}}},
    })
    assert result["rows"][0]["text_lag"] == 0
    assert result["robust"] == result["fragile"] == []
    assert result["median_robust_lag"] is None
    markdown = scanner.render(result)
    assert "not classified as robust or fragile" in markdown
    assert "None" not in markdown


def test_existing_numerical_artifact_is_unchanged_by_reassembly(d):
    assert scanner.assemble(d) == d


def test_all_stable_observed_lags_do_not_inherit_the_historical_failure_headline():
    result = scanner.assemble({
        "census": 5, "first_min": 1, "stable_min": 5, "last_full_year": 2025,
        "mesh_measurable": [], "text_measurable": ["example"],
        "arms": {"mesh": {}, "text": {"example": {"2020": 5}}},
        "trials": {"mesh": {}, "text": {"example": {"2020": 1}}},
    })
    assert result["robust"] == ["example"]
    markdown = scanner.render(result)
    assert "All 1 of 1 mechanisms" in markdown
    assert "mostly does not work" not in markdown
    assert "not answerable with the vocabulary" not in markdown
    assert "should not be quoted" not in markdown
