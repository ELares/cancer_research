"""Guards for the open-access ranking-sensitivity test (Section 3.3.1).

The analysis this replaced compared 4,830 full-text records against 5,586
abstract-only ones with the full-text side 98.7% open access -- which is not a
contrast at all, since both arms were drawn from a retrieval that had already
selected for availability. The census version splits 936,347 against 3,467,647
by PMC identifier presence while using identical expert descriptors.

THE CONFOUND IS THE POINT OF THE GUARDS, not the ranking. Identifier coverage
may reflect publication era as well as subject and access. Without the year
column and qualification, the identifier column can read as a fact about
access, an attribution this comparison has no design to support.
"""
import copy
import gzip
import importlib.util
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
JSON = REPO / "analysis/census-oa-bias.json"
MD = REPO / "analysis/census-oa-bias.md"
MANUSCRIPT = REPO / "article/drafts/v1.md"
# PINNED HERE AS A LITERAL, not read from the artifact. A first version took
# the threshold from the JSON, which made every guard that used it compare the
# generator to itself: raising it from 3 to 6 emptied the "moved" set, changed
# the headline, and passed the whole file. A judgement that can be loosened
# until the result disappears is exactly the kind that needs two files to
# change.
BIG_SHIFT = 3


@pytest.fixture(scope="module")
def d():
    return json.loads(JSON.read_text())


def _rows(d):
    return d["mechanisms"]


def test_both_arms_partition_the_same_articles(d):
    """The comparison is only meaningful if the two arms are complements."""
    for r in _rows(d):
        assert r["with_fulltext"] + r["without"] == r["total"], (
            f"{r['mechanism']}: the two arms do not sum to its census count, so "
            "they are not a partition and the ranking contrast is meaningless")
        # Stored as a FRACTION; the render multiplies by 100. Comparing against
        # a percentage here would have passed only for a mechanism at 1%.
        assert r["oa_rate"] == pytest.approx(
            r["with_fulltext"] / r["total"], abs=0.0005)


def test_the_shift_column_is_derived_from_the_two_rankings(d):
    rows = _rows(d)
    oa = {r["mechanism"]: 1 + sum(other["with_fulltext"] > r["with_fulltext"]
                                 for other in rows) for r in rows}
    non = {r["mechanism"]: 1 + sum(other["without"] > r["without"]
                                  for other in rows) for r in rows}
    for r in rows:
        assert r["rank_oa"] == oa[r["mechanism"]], r["mechanism"]
        assert r["rank_non_oa"] == non[r["mechanism"]], r["mechanism"]
        assert r["shift"] == non[r["mechanism"]] - oa[r["mechanism"]], r["mechanism"]


def test_the_headline_count_matches_the_shifts_it_counts(d):
    """An extremum or a count stated over a set nobody enumerated is how this
    repo has been wrong before.

    The moved SET is now stored in the artifact, so the count in the prose, the
    set in the JSON and the bolded rows in the table are three renderings of one
    derivation rather than three chances to disagree.
    """
    assert d["big_shift_threshold"] == BIG_SHIFT, (
        f"the generator now calls a shift 'large' at {d['big_shift_threshold']} "
        f"places where this guard pins {BIG_SHIFT}. That threshold decides how "
        "many mechanisms the headline reports as moving, so changing it must be "
        "a deliberate edit in both places rather than a quiet one here")
    rows = _rows(d)
    big = [r for r in rows if abs(r["shift"]) >= BIG_SHIFT]
    assert sorted(r["mechanism"] for r in big) == sorted(d["moved"])
    md = MD.read_text()
    assert f"**{len(big)} of {len(rows)} mechanisms shift" in md, (
        f"{len(big)} mechanisms shift {BIG_SHIFT}+ places and the report does "
        "not say so")
    for r in big:
        assert f"**{r['shift']:+d}**" in md, (
            f"{r['mechanism']} shifts {r['shift']} and is not marked in the table")


def test_the_era_confound_is_printed_beside_the_oa_column(d):
    """Read alone, the OA column reads as a fact about access. It is not one.

    Every row must carry a median year, and the report must say the confound is
    uncontrolled -- otherwise a reader takes an association for a cause.
    """
    for r in _rows(d):
        assert isinstance(r.get("median_year"), (int, float)), (
            f"{r['mechanism']} has no median year, so its OA rate has nothing "
            "to be read against")
    md = MD.read_text()
    assert "not controlled here" in md
    assert "median year" in md.lower()
    # The worked example: the lowest-OA mechanism should also be among the
    # oldest, which is what makes the confound concrete rather than abstract.
    rows = _rows(d)
    lowest = min(rows, key=lambda r: r["oa_rate"])
    median_years = sorted(r["median_year"] for r in rows)
    assert lowest["median_year"] <= sorted(median_years)[len(median_years) // 2], (
        "the lowest-OA mechanism is no longer among the older half, so the "
        "report's worked example of the confound no longer holds -- check the "
        "prose before trusting it")


def test_it_does_not_claim_the_ranking_is_unaffected(d):
    md = MD.read_text()
    for overclaim in ("ranking is stable", "no effect on the ranking",
                      "availability does not matter"):
        assert overclaim not in md.lower()


def test_the_manuscript_quotes_this_analysis_not_the_superseded_one():
    txt = " ".join(MANUSCRIPT.read_text().split())
    d = json.loads(JSON.read_text())
    rows = _rows(d)
    big = [r for r in rows if abs(r["shift"]) >= BIG_SHIFT]
    assert "census-oa-bias.md" in txt, (
        "Section 3.3.1 does not cite the census version of this test")
    assert f"{len(big)} of {len(rows)} mechanisms shift" in txt, (
        "the manuscript's rank-shift count is not the one the artifact derives")
    # The superseded contrast may still be NAMED in a retraction; what it must
    # not do is stand as the section's measurement. So the guard requires the
    # census figures to be the ones quoted, rather than banning a string that a
    # legitimate retraction also contains.
    assert f"{d['with_pmcid']:,}" in txt, (
        "Section 3.3.1 does not quote the census availability split, so it may "
        "still be resting on the superseded 4,830-vs-5,586 contrast whose two "
        "arms were both drawn from a retrieval that selected for availability")


def _module():
    spec = importlib.util.spec_from_file_location(
        "oa_bias_sparse", REPO / "scripts/census_oa_bias.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _counts(pairs, with_pmcid=None, census=None):
    rows = [dict(mechanism=name, with_fulltext=have, without=lack,
                 total=have + lack, oa_rate=have / (have + lack), median_year=None)
            for name, have, lack in pairs]
    have = sum(row["with_fulltext"] for row in rows) if with_pmcid is None else with_pmcid
    n = sum(row["total"] for row in rows) if census is None else census
    return dict(census=n, with_pmcid=have, census_oa_rate=have / n, mechanisms=rows)


def test_equal_counts_share_ranks_independent_of_input_order():
    module = _module()
    data = _counts([("a", 5, 4), ("b", 5, 5), ("c", 1, 4), ("d", 0, 1)])
    expected = {"a": (1, 2, 1), "b": (1, 1, 0),
                "c": (3, 2, -1), "d": (4, 4, 0)}
    for ordered in (data["mechanisms"], list(reversed(data["mechanisms"]))):
        result = module.derive({**data, "mechanisms": copy.deepcopy(ordered)})
        actual = {r["mechanism"]: (r["rank_oa"], r["rank_non_oa"], r["shift"])
                  for r in result["mechanisms"]}
        assert actual == expected
        assert result["moved"] == []


@pytest.mark.parametrize("pairs,have,n,empty_rank", [
    ([("a", 0, 3), ("b", 0, 1)], 0, 4, "rank_oa"),
    ([("a", 3, 0), ("b", 1, 0)], 4, 4, "rank_non_oa"),
    # Both arms contain articles, but one has no mapped mechanisms.
    ([("a", 0, 3), ("b", 0, 1)], 5, 9, "rank_oa"),
    ([("a", 3, 0), ("b", 1, 0)], 4, 9, "rank_non_oa"),
])
def test_absent_mechanism_arm_has_no_ranks_or_shift_claim(pairs, have, n, empty_rank):
    module = _module()
    data = module.derive(_counts(pairs, have, n))
    assert data["moved"] is None
    assert all(r[empty_rank] is None and r["shift"] is None for r in data["mechanisms"])
    report = module.render(data)
    assert "Ranking comparison unavailable" in report
    assert "N/A" in report
    assert "**0 of" not in report


def test_no_mapped_mechanisms_is_a_valid_zero_match_report():
    module = _module()
    data = module.derive(_counts([], with_pmcid=0, census=1))
    assert data["mechanisms"] == [] and data["moved"] is None
    report = module.render(data)
    assert "no mapped mechanism matched" in report
    assert "0.0%" in report


def test_observed_zero_shifts_remain_distinct_from_an_unavailable_comparison():
    module = _module()
    data = module.derive(_counts([("a", 5, 5), ("b", 2, 2)]))
    assert data["moved"] == []
    assert all(row["shift"] == 0 for row in data["mechanisms"])
    assert "**0 of 2 mechanisms shift" in module.render(data)


def test_historical_counts_and_derived_results_remain_unchanged(d):
    assert _module().derive(copy.deepcopy(d)) == d


def test_report_identifies_the_identifier_proxy_without_claiming_verified_access(d):
    report = _module().render(copy.deepcopy(d))
    assert "PMC identifier presence" in report
    assert "does not verify current access" in report
    assert "local full-text copy" in report


def test_recent_records_do_not_imply_high_identifier_coverage_or_a_causal_explanation():
    module = _module()
    data = _counts([("recent", 0, 1)])
    data["mechanisms"][0]["median_year"] = 2025
    report = module.render(data)
    assert "0.0%" in report and "2025" in report
    assert "only availability differs" not in report
    assert "high PMC identifier rate" not in report
    assert "does not separate these effects or quantify an era contribution" in report


def test_broad_ultrasound_bucket_is_not_presented_as_sonodynamic_specific(
        tmp_path, monkeypatch):
    module = _module()
    records = tmp_path / "records"
    records.mkdir()
    with gzip.open(records / "part.jsonl.gz", "wt", encoding="utf-8") as target:
        target.write(json.dumps({"pmid": "1", "title": "General ultrasound therapy",
                                 "mesh": ["Ultrasonic Therapy"], "year": 2020}) + "\n")
    monkeypatch.setattr(module, "RECORDS", records)
    data = module.scan()
    bucket = next(row for row in data["mechanisms"] if row["mechanism"] == "sonodynamic")
    assert bucket["total"] == 1
    report = module.render(data)
    manuscript = " ".join(MANUSCRIPT.read_text().split())
    for text in (report, manuscript):
        assert "`Ultrasonic Therapy` descriptor" in text
        assert "does not isolate sonodynamic therapy" in text
