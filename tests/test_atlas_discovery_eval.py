"""Reconstruct discovery comparisons from retained seed counts without raw data."""

import copy
import hashlib
import json
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import atlas_discovery_eval as evaluation

HISTORICAL = ROOT / "analysis/atlas-discovery-eval.json"
HISTORICAL_SHA = "2f174509b5c7ee84efa3fa423e141ca671ebd2f64bd13007018b755e16ebf843"
METHODS = ("abc", "popularity", "adamic_adar", "resource_alloc",
           "jaccard", "bridges", "random")


def seed(identity, abc, popularity, candidates=10, **counts):
    row = {"seed": identity, "seed_name": f"Entity {identity}",
           "degree": 40, "candidates": candidates,
           **{method: max(popularity - 1, 0) for method in METHODS}}
    row.update(abc=abc, popularity=popularity, random=0)
    row.update(counts)
    return row


def split(rows=None, year=2018, top=5):
    return {"split_year": year, "top_k": top,
            "pairs_before": 100, "pairs_after": 50,
            "per_seed": rows if rows is not None else
            [seed("a", 1, 3), seed("b", 2, 4)]}


def raw(rows=None, robustness=None):
    return {"headline": split(rows),
            "robustness": robustness if robustness is not None else []}


def inventory():
    return {"fingerprint_kind": "filesystem-metadata-v1",
            "inventory_sha256": "a" * 64,
            "shards_per_stream": {name: int(name == "records")
                                  for name in evaluation.YEAR_STREAMS}}


def date_support():
    return {"schema_version": 1, "year_policy": "earliest-year-per-pmid",
            "dated_pmids": 2, "record_year_min": 2010, "record_year_max": 2024,
            "dated_pairs": 150, "pair_first_year_min": 2010,
            "pair_first_year_max": 2024, "source_inventory": inventory()}


def flat(text):
    return " ".join(text.lower().split())


@pytest.fixture
def outputs(monkeypatch, tmp_path):
    md, stored = tmp_path / "report.md", tmp_path / "report.json"
    md.write_bytes(b"Previously published Markdown.\n")
    stored.write_bytes(b"Previously published JSON.\n")
    monkeypatch.setattr(evaluation, "OUT_MD", md)
    monkeypatch.setattr(evaluation, "OUT_JSON", stored)
    return md, stored


def assert_preserved(outputs, before):
    assert tuple(path.read_bytes() for path in outputs) == before


def forbid_readers(monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail("raw atlas input was accessed")

    for name in ("atlas_root", "load_index", "pmid_years", "load_pmid_years", "_scan_pmid_years",
                 "load_corrections", "pair_first_year", "evaluate"):
        monkeypatch.setattr(evaluation, name, forbidden)


def synthetic_scan(monkeypatch, tmp_path, result=None):
    source = raw() if result is None else result
    monkeypatch.setattr(evaluation, "atlas_root", lambda: tmp_path / "atlas")
    monkeypatch.setattr(evaluation, "load_index", lambda *args: {"canon": {}})
    monkeypatch.setattr(evaluation, "load_pmid_years",
                        lambda *args: ({"1": 2010, "2": 2024}, inventory()))
    monkeypatch.setattr(evaluation, "load_corrections", lambda *args: {})
    monkeypatch.setattr(evaluation, "pair_first_year",
                        lambda *args: {**{("a", f"old{i}"): 2010 for i in range(100)},
                                      **{("a", f"new{i}"): 2024 for i in range(50)}})
    monkeypatch.setattr(evaluation, "evaluate",
                        lambda *args, **kwargs: copy.deepcopy(source["headline"]))


def test_historical_counts_reconstruct_all_paired_intervals_and_published_report():
    assert hashlib.sha256(HISTORICAL.read_bytes()).hexdigest() == HISTORICAL_SHA
    stored = json.loads(HISTORICAL.read_text())
    original = copy.deepcopy(stored)
    reconstructed = evaluation.assemble(stored)
    assert len(reconstructed["robustness"]) == len(stored["robustness"])
    for expected, actual in zip([stored["headline"], *stored["robustness"]],
                                [reconstructed["headline"], *reconstructed["robustness"]]):
        for key in ("seeds_evaluated", "hits", "predictions", "precision",
                    "abc_over_popularity", "paired", "paired_all", "per_seed"):
            assert actual[key] == expected[key], (expected["split_year"], key)
    report = evaluation.render(stored)
    assert stored == original
    assert report == (ROOT / "analysis/atlas-discovery-eval.md").read_text()
    for method in ("adamic_adar", "bridges"):
        assert reconstructed["headline"]["paired_all"][method]["ci95"][1] < 0
        table_row = next(line for line in report.splitlines()
                         if line.startswith(f"| {method} |"))
        assert "spans 0" not in table_row
        assert "lower than popularity" in table_row
    prose = flat(report)
    assert "tie with it (their intervals span zero)" not in prose
    assert "candidate set is doing real work" not in prose
    assert "good candidate generator" not in prose


def test_stale_derived_fields_cannot_change_comparisons_and_raw_is_not_mutated():
    source = raw()
    expected = evaluation.assemble(source)
    stale = copy.deepcopy(source)
    stale["headline"].update(
        seeds_evaluated=999, hits={"abc": 999}, predictions={"abc": 1},
        precision={"abc": 999}, abc_over_popularity=float("inf"),
        paired={"mean_diff": 999, "ci95": [999, 999], "decided": False},
        paired_all={"abc": {"mean_diff": 999}})
    before = copy.deepcopy(stale)
    assert evaluation.assemble(stale) == expected
    assert evaluation.render(stale) == evaluation.render(source)
    assert stale == before
    # Reusing the same counts must not expose mutable cached result objects.
    expected["headline"]["paired_all"]["abc"]["ci95"][0] = 999
    assert evaluation.assemble(source)["headline"]["paired_all"]["abc"]["ci95"] == [-2, -2]


def test_short_candidate_pools_use_actual_predictions_in_precision():
    short = seed("short", 2, 2, candidates=2)
    short.update(dict.fromkeys(METHODS, 2))
    source = raw([short, seed("long", 3, 4, candidates=10)])
    result = evaluation.assemble(source)["headline"]
    assert result["seeds_evaluated"] == 2
    assert result["predictions"] == dict.fromkeys(METHODS, 7)
    assert result["hits"]["abc"] == 5
    assert result["hits"]["popularity"] == 6
    assert result["precision"]["abc"] == pytest.approx(5 / 7)
    assert result["paired"]["mean_diff"] == -0.5
    assert result["paired"]["ci95"] == [-1, 0]
    assert result["paired"]["abc_ahead"] == 0
    assert result["paired"]["abc_behind"] == 1


@pytest.mark.parametrize("candidates,abc,popularity", [
    (2, 2, 1), (5, 3, 2), (6, 3, 1), (8, 5, 1),
], ids=["short-pool", "exact-full-pool", "one-omitted", "three-omitted"])
def test_shared_pool_rejects_impossible_differences_between_ranking_hits(
        candidates, abc, popularity):
    row = seed("a", abc, popularity, candidates=candidates)
    row.update(dict.fromkeys(METHODS, popularity))
    row["abc"] = abc
    with pytest.raises(ValueError, match="shared candidate pool"):
        evaluation.assemble(raw([row]))


@pytest.mark.parametrize("candidates,abc,popularity", [
    (2, 1, 1), (5, 3, 3), (6, 3, 2), (8, 4, 1),
], ids=["short-pool", "exact-full-pool", "one-omitted", "three-omitted"])
def test_shared_pool_accepts_realisable_ranking_hit_counts(
        candidates, abc, popularity):
    row = seed("a", abc, popularity, candidates=candidates)
    row.update(dict.fromkeys(METHODS, popularity))
    row["abc"] = abc
    result = evaluation.assemble(raw([row]))["headline"]
    assert result["hits"]["abc"] == abc
    assert result["hits"]["popularity"] == popularity


def test_evaluate_preserves_pre_split_ranking_sampling_and_paired_results():
    """Freeze a small graph's outputs from the evaluator before report refactoring.

    The graph includes withheld edges and gives every seed eligible bridge
    candidates. Its different rankings and shuffled control expose changes to
    graph timing, seed order, tie handling, candidate selection, and hit scoring.
    """
    first = {}
    for i in range(120):
        for j in range(i + 1, 120):
            pair = (f"{i:03}", f"{j:03}")
            if (i * 37 + j * 13) % 7 < 3:
                first[pair] = 2010
            elif (i * 17 + j * 23) % 5 == 0:
                first[pair] = 2020
    result = evaluation.evaluate(first, {"canon": {}}, 2018, 8, 5, log=False)
    assert result["pairs_before"] == 3060
    assert result["pairs_after"] == 789
    assert result["seeds_evaluated"] == 8
    assert result["predictions"] == dict.fromkeys(METHODS, 40)
    rows = result["per_seed"]
    assert [row["seed"] for row in rows] == ["110", "018", "030", "067", "080", "045", "008", "102"]
    assert [row["degree"] for row in rows] == [50, 52, 51, 52, 52, 52, 50, 52]
    assert [row["candidates"] for row in rows] == [68, 67, 68, 67, 67, 67, 64, 67]
    for method in ("abc", "adamic_adar", "resource_alloc", "jaccard", "bridges"):
        assert [row[method] for row in rows] == [0, 0, 1, 0, 0, 0, 1, 0]
        assert result["precision"][method] == 0.05
        assert result["paired_all"][method] == {
            "mean_diff": -0.625, "ci95": [-0.875, -0.25], "decided": True,
            "ahead": 0, "behind": 5}
    assert [row["popularity"] for row in rows] == [1, 0, 1, 1, 1, 1, 1, 1]
    assert [row["random"] for row in rows] == [1, 0, 2, 1, 0, 0, 0, 2]
    assert result["precision"]["popularity"] == 0.175
    assert result["precision"]["random"] == 0.15
    assert result["paired_all"]["random"] == {
        "mean_diff": -0.125, "ci95": [-0.625, 0.375], "decided": False,
        "ahead": 2, "behind": 3}


@pytest.mark.parametrize("counts,mean,interval,decided,ahead,behind", [
    ([(4, 2), (3, 1)], 2, [2, 2], True, 2, 0),
    ([(1, 3), (2, 4)], -2, [-2, -2], True, 0, 2),
    ([(3, 2), (3, 2), (1, 2)], 1 / 3, [-1, 1], False, 2, 1),
    ([(1, 2), (1, 2), (3, 2)], -1 / 3, [-1, 1], False, 1, 2),
    ([(2, 2), (2, 2)], 0, [0, 0], False, 0, 0),
    ([(3, 2), (1, 2)], 0, [-1, 1], False, 1, 1),
    ([(0, 0), (0, 0)], 0, [0, 0], False, 0, 0),
], ids=["abc-wins", "abc-loses", "uncertain-abc-leads", "uncertain-popularity-leads",
        "same-counts", "equal-means-different-seeds", "zero-hits"])
def test_paired_comparisons_distinguish_direction_uncertainty_and_observed_ties(
        counts, mean, interval, decided, ahead, behind):
    source = raw([seed(str(index), a, p) for index, (a, p) in enumerate(counts)])
    paired = evaluation.assemble(source)["headline"]["paired"]
    assert paired["mean_diff"] == pytest.approx(mean)
    assert paired["ci95"] == interval
    assert paired["decided"] is decided
    assert paired["abc_ahead"] == ahead
    assert paired["abc_behind"] == behind
    report = evaluation.render(source)
    verdict = next(line for line in report.splitlines() if line.startswith("### Verdict:"))
    row = next(line for line in report.splitlines() if line.startswith("| abc |"))
    assert f"{mean:+.2f} [{interval[0]:+.2f}, {interval[1]:+.2f}]" in row
    if decided:
        assert ("beats popularity" if mean > 0 else "performs worse than popularity") in verdict
        assert ("higher than popularity" if mean > 0 else "lower than popularity") in row
    else:
        assert "inconclusive" in verdict
        assert "interval includes zero" in row
        prose = flat(report)
        expected = ("point estimate favors abc" if mean > 0 else
                    "point estimate favors popularity" if mean < 0 else
                    "point estimates are equal")
        assert expected in prose
        assert "does not establish that two rankings are equal" in prose


@pytest.mark.parametrize("other_counts,other_conclusion,unanimous", [
    ([(1, 3), (2, 4)], "lower than popularity", True),
    ([(4, 2), (3, 1)], "higher than popularity", False),
    ([(3, 2), (1, 2)], "interval includes zero", False),
], ids=["all-losses", "mixed-directions", "one-unresolved"])
def test_robustness_claims_use_every_requested_split(other_counts, other_conclusion, unanimous):
    other = split([seed(str(i), a, p) for i, (a, p) in enumerate(other_counts)], year=2021)
    source = raw(robustness=[split(year=2015), other])
    result = evaluation.assemble(source)
    assert [item["split_year"] for item in result["robustness"]] == [2015, 2021]
    report = evaluation.render(source)
    row = next(line for line in report.splitlines() if line.startswith("| 2021 |"))
    assert other_conclusion in row
    prose = flat(report)
    if unanimous:
        assert "abc is lower than popularity in every evaluated split" in prose
    else:
        assert "do not all resolve the abc comparison in the same direction" in prose
        assert "abc is lower than popularity in every evaluated split" not in prose


def test_beating_within_pool_random_does_not_establish_candidate_generator_value():
    source = raw([seed("a", 5, 4, random=1), seed("b", 5, 4, random=1)])
    result = evaluation.assemble(source)["headline"]
    assert result["precision"]["abc"] == 1.0
    assert result["precision"]["popularity"] == 0.8
    assert result["precision"]["random"] == 0.2
    prose = flat(evaluation.render(source))
    assert ("outperforming the within-pool random baseline does **not** establish that "
            "the candidate generator improves on unrestricted candidate selection") in prose
    assert "no such control is evaluated here" in prose
    assert "good candidate generator" not in prose


@pytest.mark.parametrize("abc", [0, 2])
def test_zero_popularity_baseline_is_valid_with_no_ratio(abc):
    source = raw([seed("a", abc, 0), seed("b", abc, 0)])
    result = evaluation.assemble(source)
    assert result["headline"]["precision"]["popularity"] == 0
    assert result["headline"]["abc_over_popularity"] is None
    json.dumps(result, allow_nan=False)
    report = evaluation.render(source)
    assert "infinity" not in flat(report)
    assert "nan" not in flat(report)


@pytest.mark.parametrize("field,value", [
    ("abc", -1), ("abc", True), ("abc", 1.0), ("abc", "1"),
    ("abc", None), ("abc", 6), ("popularity", float("nan")),
    ("random", float("inf")), ("candidates", 0), ("candidates", True),
    ("degree", -1), ("degree", False), ("seed", ""), ("seed", 1),
    ("seed_name", None),
])
def test_invalid_seed_counts_and_types_are_rejected(field, value):
    source = raw()
    source["headline"]["per_seed"][0][field] = value
    with pytest.raises(ValueError):
        evaluation.assemble(source)


@pytest.mark.parametrize("field,value", [
    ("split_year", 0), ("split_year", True), ("split_year", "2018"),
    ("top_k", 0), ("top_k", -1), ("top_k", 2.0),
    ("pairs_before", -1), ("pairs_after", False),
    ("per_seed", []), ("per_seed", {}), ("per_seed", [None]),
])
def test_invalid_split_metadata_is_rejected(field, value):
    source = raw()
    source["headline"][field] = value
    with pytest.raises(ValueError):
        evaluation.assemble(source)


def test_a_hit_count_cannot_exceed_a_short_candidate_pool():
    source = raw([seed("a", 3, 1, candidates=2)])
    with pytest.raises(ValueError):
        evaluation.assemble(source)


@pytest.mark.parametrize("pairs_after,method,hits", [
    (0, "abc", 1), (1, "random", 2),
], ids=["no-observed-future-pairs", "more-hits-than-future-pairs"])
def test_seed_hits_cannot_exceed_observed_post_split_pairs(pairs_after, method, hits):
    source = raw([seed("a", 0, 0)])
    source["headline"]["pairs_after"] = pairs_after
    source["headline"]["per_seed"][0][method] = hits
    with pytest.raises(ValueError, match="observed post-split pairs"):
        evaluation.assemble(source)


@pytest.mark.parametrize("pairs_after", [0, 2], ids=["no-future-pairs", "all-future-pairs-hit"])
def test_seed_can_hit_every_observed_post_split_pair(pairs_after):
    source = raw([seed("a", pairs_after, 0)])
    source["headline"]["pairs_after"] = pairs_after
    result = evaluation.assemble(source)["headline"]
    assert result["hits"]["abc"] == pairs_after


@pytest.mark.parametrize("missing", METHODS)
def test_missing_ranking_counts_are_not_interpreted_as_zero(missing):
    source = raw()
    del source["headline"]["per_seed"][0][missing]
    with pytest.raises(ValueError):
        evaluation.assemble(source)


@pytest.mark.parametrize("duplicate", ["seed", "headline-year", "robustness-year"])
def test_duplicate_evaluation_identities_are_rejected(duplicate):
    source = raw()
    if duplicate == "seed":
        source["headline"]["per_seed"][1]["seed"] = "a"
    elif duplicate == "headline-year":
        source["robustness"] = [split(year=2018)]
    else:
        source["robustness"] = [split(year=2015), split(year=2015)]
    with pytest.raises(ValueError):
        evaluation.assemble(source)


def test_render_only_is_offline_and_preserves_historical_json_format(outputs, monkeypatch):
    md, stored = outputs
    source = raw()
    # Preserve spacing and stale derived data as well as the scientific counts.
    source["headline"]["precision"] = {"abc": 8765}
    stored.write_text(json.dumps(source, indent=4) + "\n\n")
    before = stored.read_bytes()
    forbid_readers(monkeypatch)
    assert evaluation.main(["--render-only"]) == 0
    assert stored.read_bytes() == before
    assert md.read_text() == evaluation.render(source)
    assert "8765" not in md.read_text()


@pytest.mark.parametrize("problem", ["bad-json", "bad-count", "empty-seeds", "impossible-pool",
                                     "insufficient-future-pairs"])
def test_invalid_offline_input_preserves_both_reports(outputs, monkeypatch, problem):
    _, stored = outputs
    source = raw()
    if problem == "bad-count":
        source["headline"]["per_seed"][0]["abc"] = 6
    elif problem == "empty-seeds":
        source["headline"]["per_seed"] = []
    elif problem == "impossible-pool":
        row = source["headline"]["per_seed"][0]
        row.update(candidates=2, **dict.fromkeys(METHODS, 0))
        row["abc"] = 2
    elif problem == "insufficient-future-pairs":
        source["headline"]["pairs_after"] = 0
    stored.write_text("{broken" if problem == "bad-json" else json.dumps(source))
    before = tuple(path.read_bytes() for path in outputs)
    forbid_readers(monkeypatch)
    assert evaluation.main(["--render-only"]) != 0
    assert_preserved(outputs, before)


@pytest.mark.parametrize("offline", [False, True])
@pytest.mark.parametrize("stage", ["render", "serialize"])
def test_preparation_failures_preserve_both_reports(
        outputs, monkeypatch, tmp_path, offline, stage):
    _, stored = outputs
    if offline:
        stored.write_text(json.dumps(raw()))
        forbid_readers(monkeypatch)
    else:
        synthetic_scan(monkeypatch, tmp_path)
    before = tuple(path.read_bytes() for path in outputs)

    def fail(*args, **kwargs):
        raise TypeError("synthetic report preparation failure")

    if stage == "render":
        monkeypatch.setattr(evaluation, "render", fail)
    else:
        monkeypatch.setattr(evaluation.json, "dumps", fail)
    with pytest.raises(TypeError, match="synthetic report preparation failure"):
        evaluation.main(["--render-only"] if offline else ["--also-years"])
    assert_preserved(outputs, before)


@pytest.mark.parametrize("problem", ["no-years", "no-seeds", "invalid-count",
                                     "insufficient-future-pairs"])
def test_invalid_fresh_inputs_preserve_both_reports(
        outputs, monkeypatch, tmp_path, problem):
    synthetic_scan(monkeypatch, tmp_path)
    if problem == "no-years":
        monkeypatch.setattr(evaluation, "load_pmid_years",
                            lambda *args: ({}, inventory()))
    elif problem == "no-seeds":
        monkeypatch.setattr(evaluation, "evaluate", lambda *args, **kwargs: None)
    else:
        source = raw()
        if problem == "invalid-count":
            source["headline"]["per_seed"][0]["abc"] = -1
        else:
            source["headline"]["pairs_after"] = 0
        synthetic_scan(monkeypatch, tmp_path, source)
    before = tuple(path.read_bytes() for path in outputs)
    assert evaluation.main(["--also-years"]) != 0
    assert_preserved(outputs, before)


def test_an_unevaluable_requested_robustness_split_is_not_silently_omitted(
        outputs, monkeypatch, tmp_path):
    synthetic_scan(monkeypatch, tmp_path)

    def evaluate(first, index, year, *args, **kwargs):
        return split(year=year) if year == 2018 else None

    monkeypatch.setattr(evaluation, "evaluate", evaluate)
    before = tuple(path.read_bytes() for path in outputs)
    assert evaluation.main(["--also-years", "2015"]) != 0
    assert_preserved(outputs, before)


def test_a_requested_year_outside_dated_support_fails_before_pair_loading(
        outputs, monkeypatch, tmp_path):
    synthetic_scan(monkeypatch, tmp_path)
    monkeypatch.setattr(evaluation, "pair_first_year",
                        lambda *args: pytest.fail("unsupported split loaded pairs"))
    before = tuple(path.read_bytes() for path in outputs)
    assert evaluation.main(["--also-years", "2030"]) != 0
    assert_preserved(outputs, before)


@pytest.mark.parametrize("args", [
    ["--seeds", "0"], ["--seeds", "-1"], ["--top", "0"],
    ["--top", "-2"], ["--split-year", "0"], ["--split-year", "-1"],
    ["--also-years", "2015", "0"], ["--also-years", "-1"],
])
def test_cli_rejects_invalid_parameters_before_reading_inputs(outputs, monkeypatch, args):
    forbid_readers(monkeypatch)
    before = tuple(path.read_bytes() for path in outputs)
    with pytest.raises(SystemExit) as error:
        evaluation.main(args)
    assert error.value.code == 2
    assert_preserved(outputs, before)


def test_fresh_zero_hit_results_serialize_strictly_and_render(outputs, monkeypatch, tmp_path):
    source = raw([seed("a", 0, 0), seed("b", 0, 0)])
    synthetic_scan(monkeypatch, tmp_path, source)
    assert evaluation.main(["--also-years"]) == 0
    md, stored = outputs
    result = json.loads(stored.read_text(), parse_constant=lambda value: pytest.fail(value))
    assert result["headline"]["hits"] == dict.fromkeys(METHODS, 0)
    assert result["headline"]["abc_over_popularity"] is None
    assert md.read_text() == evaluation.render(result)


def test_fresh_date_support_is_retained_and_rendered_without_claiming_a_cutoff(
        outputs, monkeypatch, tmp_path):
    synthetic_scan(monkeypatch, tmp_path)
    assert evaluation.main(["--also-years"]) == 0
    md, stored = outputs
    published = json.loads(stored.read_text())
    assert published["date_support"] == date_support()
    original = copy.deepcopy(published)
    rebuilt = evaluation.assemble(published)
    assert rebuilt == published
    rebuilt["date_support"]["source_inventory"]["shards_per_stream"]["records"] = 9
    assert published == original
    prose = flat(md.read_text())
    assert "2 dated pmids spanning **2010-2024**" in prose
    assert "150 dated pairs" in prose
    assert "not a complete observation cutoff" in prose
    assert "filesystem metadata, not article contents" in prose
    before_json = stored.read_bytes()
    forbid_readers(monkeypatch)
    assert evaluation.main(["--render-only"]) == 0
    assert stored.read_bytes() == before_json
    assert md.read_text() == evaluation.render(published)


@pytest.mark.parametrize("key,value", [
    ("schema_version", True), ("schema_version", 2),
    ("year_policy", "last-seen-year"), ("dated_pmids", 0),
    ("dated_pmids", True), ("record_year_min", None),
    ("record_year_max", "2024"), ("record_year_min", 2025),
    ("record_year_max", 2009), ("dated_pairs", 149),
    ("pair_first_year_min", 2009), ("pair_first_year_max", 2025),
    ("pair_first_year_min", 2025), ("pair_first_year_min", 2018),
    ("pair_first_year_max", 2017), ("source_inventory", None),
])
def test_invalid_retained_date_support_preserves_reports(
        outputs, monkeypatch, key, value):
    md, stored = outputs
    source = raw()
    source["date_support"] = date_support()
    source["date_support"][key] = value
    stored.write_text(json.dumps(source))
    before = tuple(path.read_bytes() for path in outputs)
    forbid_readers(monkeypatch)
    assert evaluation.main(["--render-only"]) != 0
    assert_preserved(outputs, before)


@pytest.mark.parametrize("key,value", [
    ("fingerprint_kind", "sha256-content"),
    ("inventory_sha256", None), ("inventory_sha256", "a" * 63),
    ("inventory_sha256", "Z" * 64), ("shards_per_stream", []),
    ("shards_per_stream", {"records": 1}),
    ("shards_per_stream", dict.fromkeys(evaluation.YEAR_STREAMS, 0)),
    ("shards_per_stream", dict.fromkeys(evaluation.YEAR_STREAMS, True)),
    ("shards_per_stream", dict.fromkeys(evaluation.YEAR_STREAMS, -1)),
])
def test_invalid_retained_inventory_provenance_is_rejected(key, value):
    source = raw()
    source["date_support"] = date_support()
    source["date_support"]["source_inventory"][key] = value
    with pytest.raises(ValueError):
        evaluation.assemble(source)


@pytest.mark.parametrize("support", [None, [], "unrecorded"])
def test_present_but_invalid_support_is_not_treated_as_legacy(support):
    source = raw()
    source["date_support"] = support
    with pytest.raises(ValueError, match="date_support"):
        evaluation.render(source)


def test_all_splits_must_reconcile_with_the_same_dated_pair_population():
    source = raw(robustness=[split(year=2021)])
    source["date_support"] = date_support()
    assert evaluation.assemble(source)["date_support"] == date_support()
    source["robustness"][0]["pairs_before"] += 1
    with pytest.raises(ValueError, match="dated_pairs"):
        evaluation.assemble(source)


@pytest.mark.parametrize("year,before", [(2021, 50), (2015, 125)])
def test_decreasing_pre_split_pair_counts_preserve_reports(
        outputs, monkeypatch, year, before):
    source = raw(robustness=[split(year=year)])
    source["robustness"][0].update(pairs_before=before, pairs_after=150 - before)
    source["date_support"] = date_support()
    source["date_support"]["dated_pmids"] = 4
    with pytest.raises(ValueError, match="nondecreasing"):
        evaluation.assemble(source)
    outputs[1].write_text(json.dumps(source))
    previous = tuple(path.read_bytes() for path in outputs)
    forbid_readers(monkeypatch)
    assert evaluation.main(["--render-only"]) != 0
    assert_preserved(outputs, previous)


@pytest.mark.parametrize("before_counts", [(100, 100, 100), (50, 100, 125)])
def test_monotonic_pair_counts_preserve_scrambled_split_order(before_counts):
    source = raw(robustness=[split(year=2021), split(year=2015)])
    by_year = dict(zip((2015, 2018, 2021), before_counts))
    for item in [source["headline"], *source["robustness"]]:
        before = by_year[item["split_year"]]
        item.update(pairs_before=before, pairs_after=150 - before)
    source["date_support"] = date_support()
    source["date_support"]["dated_pmids"] = 4
    original = copy.deepcopy(source)
    result = evaluation.assemble(source)
    assert source == original
    assert result["date_support"] == source["date_support"]
    assert [(item["split_year"], item["pairs_before"])
            for item in [result["headline"], *result["robustness"]]] == [
                (2018, by_year[2018]), (2021, by_year[2021]), (2015, by_year[2015])]


def test_legacy_summaries_without_date_support_keep_split_validation_contract():
    source = raw(robustness=[split(year=2021)])
    source["robustness"][0].update(pairs_before=50, pairs_after=100)
    result = evaluation.assemble(source)
    assert "date_support" not in result
    assert result["robustness"][0]["pairs_before"] == 50


def test_record_span_may_extend_beyond_pair_dates_and_zero_future_pairs_are_valid():
    rows = [seed("zero", 0, 0)]
    rows[0].update(dict.fromkeys(METHODS, 0))
    source = raw(rows)
    source["headline"].update(pairs_before=150, pairs_after=0)
    support = date_support()
    support["pair_first_year_max"] = 2017
    source["date_support"] = support
    assert evaluation.assemble(source)["date_support"] == support
    assert "**2010-2017**" in evaluation.render(source)


def test_added_revision_shard_changes_pair_from_future_hit_to_prior_knowledge(tmp_path):
    import gzip

    records, updates, relations = (tmp_path / name for name in
                                   ("records", "records_updates", "relations"))
    for path in (records, updates, relations):
        path.mkdir()
    with gzip.open(records / "base.jsonl.gz", "wt") as stream:
        stream.write(json.dumps({"pmid": "1", "year": 2025}) + "\n")
    with gzip.open(relations / "relations.tsv.gz", "wt") as stream:
        stream.write("1\trelation\tGene|a\tGene|b\n")
    years = evaluation.pmid_years(tmp_path)
    assert evaluation.pair_first_year(tmp_path, years, {}) == {("a", "b"): 2025}
    with gzip.open(updates / "revision.jsonl.gz", "wt") as stream:
        stream.write(json.dumps({"pmid": "1", "year": 2010}) + "\n")
    years = evaluation.pmid_years(tmp_path)
    assert evaluation.pair_first_year(tmp_path, years, {}) == {("a", "b"): 2010}


@pytest.mark.parametrize("count", ["dated_pmids", "dated_pairs"])
def test_single_dated_observation_cannot_claim_multiple_years(count):
    support = date_support()
    support[count] = 1
    with pytest.raises(ValueError, match="single"):
        evaluation._assemble_date_support(support, [])
