"""Units, sparse outcomes, reproducibility and output preservation for diagnostics."""

import copy
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import atlas_contradiction_quality as quality

HISTORICAL = ROOT / "analysis/atlas-contradiction-quality.json"
HISTORICAL_SHA = "f02453b3f16e918b7b15ec074cb984cc587d1e602706b951bbaf3298bba12e9b"


def ambiguity():
    return {"by_type": {name: {"sense_rows": []}
                        for name in ("gene", "chemical", "disease")}}


def relation(pmid="1", predicate="positive_correlate", a="Gene|101", b="Gene|102"):
    return f"{pmid}\t{predicate}\t{a}\t{b}\n"


@pytest.fixture
def outputs(monkeypatch, tmp_path):
    atlas = tmp_path / "external data"
    (atlas / "relations").mkdir(parents=True)
    monkeypatch.setenv("FERRO_ATLAS_ROOT", str(atlas))
    scan = tmp_path / "ambiguity.json"
    scan.write_text(json.dumps(ambiguity()))
    monkeypatch.setattr(quality, "SCAN", scan)
    md, raw = tmp_path / "report.md", tmp_path / "report.json"
    md.write_text("previous markdown")
    raw.write_text("previous JSON")
    monkeypatch.setattr(quality, "OUT_MD", md)
    monkeypatch.setattr(quality, "OUT_JSON", raw)
    return atlas, md, raw


def write_relations(root, text):
    with gzip.open(root / "relations/relations.tsv.gz", "wt", encoding="utf-8") as stream:
        stream.write(text)


def failed_main():
    try:
        result = quality.main([])
    except (ValueError, OSError) as error:
        assert str(error)
    except SystemExit as error:
        assert error.code not in (None, 0)
    else:
        assert result != 0, "invalid input must not succeed"


def test_historical_json_is_preserved_and_uses_pair_units():
    assert hashlib.sha256(HISTORICAL.read_bytes()).hexdigest() == HISTORICAL_SHA
    raw = json.loads(HISTORICAL.read_text())
    before = copy.deepcopy(raw)
    md = quality.render(raw)
    assert raw == before
    assert "196,363" in md and "196,364" in md and "6,373" in md
    assert "pair-PMID" in md
    assert "unique-paper counts were not retained" in md.lower()
    assert "1.45x" in md
    assert "115,024" not in md
    assert "conflicts really are" not in md
    assert "biases the ratio DOWN" not in md
    assert "RISES with assertion count" not in md


def test_render_recomputes_point_estimates_but_preserves_historical_interval():
    raw = json.loads(HISTORICAL.read_text())
    raw["crude_risk_ratio"] = 9876.0
    raw["mantel_haenszel"] = 6543.0
    md = quality.render(raw)
    assert "1.46x" in md and "1.45x" in md
    assert "9876" not in md and "6543" not in md
    assert "1.39" in md and "1.50" in md
    assert "historical" in md.lower()
    assert "independent" in md.lower()


def test_capped_highest_support_bin_is_open_ended():
    md = quality.render(json.loads(HISTORICAL.read_text()))
    assert "1024+" in md
    assert "1024-2047" not in md


def test_duplicate_and_reversed_rows_count_one_pmid_per_direction(tmp_path):
    (tmp_path / "relations").mkdir()
    write_relations(tmp_path, relation() + relation() + relation(a="Gene|102", b="Gene|101")
                    + relation(predicate="negative_correlate") + relation(predicate="treat"))
    pos, neg = quality.load_directional(tmp_path)
    assert dict(pos) == {("101", "102"): {"1"}}
    assert dict(neg) == {("101", "102"): {"1"}}


@pytest.mark.parametrize("content", ["", "\n", "bad\n", relation() + "broken later row\n",
                                    relation(pmid="not-a-pmid"), relation(a="Gene|"),
                                    relation(b="")])
def test_empty_or_malformed_relations_fail_without_replacing_reports(outputs, content):
    atlas, md, raw = outputs
    write_relations(atlas, content)
    failed_main()
    assert md.read_text() == "previous markdown"
    assert raw.read_text() == "previous JSON"


@pytest.mark.parametrize("problem", ["missing", "corrupt", "utf8"])
def test_unreadable_relation_dump_preserves_both_outputs(outputs, problem):
    atlas, md, raw = outputs
    path = atlas / "relations/relations.tsv.gz"
    if problem == "corrupt":
        path.write_bytes(b"not gzip")
    elif problem == "utf8":
        path.write_bytes(gzip.compress(relation().encode() + b"\xff\n"))
    failed_main()
    assert md.read_text() == "previous markdown"
    assert raw.read_text() == "previous JSON"


@pytest.mark.parametrize("content", [relation(predicate="treat"), relation()])
def test_readable_zero_eligible_population_is_valid_with_null_comparisons(outputs, content):
    atlas, md, raw = outputs
    write_relations(atlas, content)
    assert quality.main([]) == 0
    result = json.loads(raw.read_text(), parse_constant=lambda value: pytest.fail(value))
    assert result["eligible_pairs"] == 0
    assert result["crude_risk_ratio"] is None
    assert result["mantel_haenszel"] is None
    assert result["mh_ci95"] is None
    assert "not estimable" in md.read_text().lower()
    assert "nan" not in md.read_text().lower()


def test_same_papers_on_multiple_pairs_are_not_a_unique_paper_denominator():
    pos = {("a", "b"): set("12345"), ("c", "d"): set("12345")}
    neg = {("a", "b"): set("167"), ("c", "d"): set("167")}
    result = quality.analyze(pos, neg, {"a"}, bootstrap=64)
    assert result["conflicting_pairs"] == 2
    assert result["pairs_with_self_contradiction"] == 2
    assert result["self_contradicting_assertions"] == 2
    assert result["total_assertions_in_conflicts"] == 16
    assert result["total_assertions_in_conflicts"] - result["self_contradicting_assertions"] == 14
    assert result["unique_papers_in_conflicts"] == 7
    assert result["unique_papers_with_overlap"] == 1
    assert result["crude_risk_ratio"] == 1.0
    assert result["mantel_haenszel"] == 1.0


def test_no_clean_conflicts_has_no_finite_risk_ratio_or_interval():
    pos = {("a", "b"): set("12345"), ("c", "d"): set("12345678")}
    neg = {("a", "b"): set("678")}
    result = quality.analyze(pos, neg, {"a"}, bootstrap=64)
    assert result["eligible_pairs"] == 2
    assert result["ambiguous"] == {"n": 1, "conflicted": 1}
    assert result["clean"] == {"n": 1, "conflicted": 0}
    assert result["crude_risk_ratio"] is None
    assert result["mantel_haenszel"] is None
    assert result["mh_ci95"] is None
    assert result["bootstrap"]["finite_resamples"] == 0
    assert result["bootstrap"]["undefined_resamples"] == 64
    json.dumps(result, allow_nan=False)


def test_one_stratum_reduces_to_the_ordinary_risk_ratio():
    pos, neg, contested = {}, {}, set()
    for group, flagged_count in (("amb", 25), ("clean", 10)):
        for index in range(50):
            key = (f"{group}{index}", f"partner{index}")
            flagged = index < flagged_count
            pos[key] = set("12345" if flagged else "12345678")
            neg[key] = set("678") if flagged else set()
            if group == "amb":
                contested.add(key[0])
    result = quality.analyze(pos, neg, contested, bootstrap=128)
    assert result["crude_risk_ratio"] == 2.5
    assert result["mantel_haenszel"] == 2.5
    assert result["bootstrap"]["finite_resamples"] == 128
    assert result["bootstrap"]["undefined_resamples"] == 0
    assert result["mh_ci95"][0] < 2.5 < result["mh_ci95"][1]


@pytest.mark.parametrize("mutation", ["negative", "boolean", "margin", "overlap", "bucket",
                                     "nan_interval", "reversed_interval", "unique_papers",
                                     "support_lower_bound", "support_upper_bound"])
def test_invalid_offline_counts_preserve_the_existing_reports(outputs, mutation):
    _, md, raw_path = outputs
    raw = json.loads(HISTORICAL.read_text())
    if mutation == "negative":
        raw["eligible_pairs"] = -1
    elif mutation == "boolean":
        raw["self_contradicting_assertions"] = True
    elif mutation == "margin":
        raw["strata"]["3"]["amb"][0] += 1
    elif mutation == "overlap":
        raw["self_contradicting_assertions"] = raw["total_assertions_in_conflicts"]
    elif mutation == "bucket":
        raw["strata"]["11"] = raw["strata"].pop("10")
    elif mutation == "nan_interval":
        raw["mh_ci95"][1] = float("nan")
    elif mutation == "reversed_interval":
        raw["mh_ci95"].reverse()
    elif mutation == "unique_papers":
        raw["unique_papers_in_conflicts"] = 1
        raw["unique_papers_with_overlap"] = 2
    elif mutation == "support_lower_bound":
        raw["total_assertions_in_conflicts"] = raw["conflicting_pairs"] * quality.MIN_TOTAL
    elif mutation == "support_upper_bound":
        raw = quality.analyze({("a", "b"): set("12345")},
                              {("a", "b"): set("678")}, set(), bootstrap=0)
        raw["total_assertions_in_conflicts"] = 16  # the retained 8–15 stratum excludes this
    raw_path.write_text(json.dumps(raw))
    before = raw_path.read_bytes()
    assert quality.main(["--render-only"]) != 0
    assert raw_path.read_bytes() == before
    assert md.read_text() == "previous markdown"


@pytest.mark.parametrize("supports, field, lowered_count", [
    ([(10, 5, 0)], "unique_papers_in_conflicts", 8),
    ([(5, 3, 2), (5, 3, 0)], "unique_papers_with_overlap", 1),
    ([(13, 3, 0)] + [(5, 3, 0)] * 4, "unique_papers_in_conflicts", 10),
    ([(12, 4, 4)] + [(5, 3, 0)] * 4, "unique_papers_in_conflicts", 9),
    ([(1021, 3, 0)] + [(5, 3, 0)] * 128, "unique_papers_in_conflicts", 16),
], ids=["pair-pmid-total", "overlap-total", "largest-bin-no-overlap",
        "largest-bin-with-overlap", "capped-largest-bin"])
def test_impossible_unique_paper_lower_bounds_preserve_reports(
        outputs, supports, field, lowered_count):
    _, md, raw_path = outputs
    pos, neg = {}, {}
    for index, (positive, negative, overlap) in enumerate(supports):
        pair = (f"entity{index}", f"partner{index}")
        pos[pair] = {str(p) for p in range(positive)}
        neg[pair] = {str(p) for p in range(positive - overlap,
                                          positive - overlap + negative)}
    raw = quality.analyze(pos, neg, set(), bootstrap=0)
    # These feasible examples share PMIDs across pairs, including at bin edges.
    assert quality.render(raw)
    raw[field] = lowered_count
    raw_path.write_text(json.dumps(raw))
    before = raw_path.read_bytes()
    assert quality.main(["--render-only"]) != 0
    assert raw_path.read_bytes() == before
    assert md.read_text() == "previous markdown"


@pytest.mark.parametrize("scan", [None, {}, {"by_type": {}},
                                 {"by_type": {"gene": {"sense_rows": None}}}])
def test_invalid_ambiguity_input_preserves_both_outputs(outputs, scan):
    atlas, md, raw = outputs
    write_relations(atlas, relation())
    quality.SCAN.write_text(json.dumps(scan))
    failed_main()
    assert md.read_text() == "previous markdown"
    assert raw.read_text() == "previous JSON"


@pytest.mark.parametrize("stage", ["render", "serialize"])
def test_report_preparation_failure_preserves_both_outputs(outputs, monkeypatch, stage):
    atlas, md, raw = outputs
    write_relations(atlas, relation())

    def fail(*args, **kwargs):
        raise ValueError("synthetic preparation failure")

    if stage == "render":
        monkeypatch.setattr(quality, "render", fail)
    else:
        monkeypatch.setattr(quality.json, "dumps", fail)
    failed_main()
    assert md.read_text() == "previous markdown"
    assert raw.read_text() == "previous JSON"


def test_offline_cli_needs_no_raw_inputs_and_never_rewrites_json(outputs, monkeypatch):
    atlas, md, raw = outputs
    raw.write_bytes(HISTORICAL.read_bytes())
    quality.SCAN.unlink()
    monkeypatch.setattr(quality, "load_directional", lambda *args: pytest.fail("offline scan"))
    before = raw.read_bytes()
    assert quality.main(["--render-only"]) == 0
    assert raw.read_bytes() == before
    assert md.read_text() == quality.render(json.loads(before))


def test_seeded_pair_resampling_is_independent_of_python_hash_seed():
    code = '''
import json, sys
sys.path.insert(0, sys.argv[1])
import atlas_contradiction_quality as q
pos, neg, contested = {}, {}, set()
for i in range(32):
    pair = (f"a{i}", f"b{i}")
    pos[pair] = {str(p) for p in range(5 + i % 7)}
    neg[pair] = {str(100+p) for p in range(i % 6)}
    if i % 3:
        contested.add(pair[0])
print(json.dumps(q.analyze(pos, neg, contested, bootstrap=256), sort_keys=True, allow_nan=False))
'''
    results = []
    for seed in ("1", "8675309"):
        result = subprocess.run([sys.executable, "-c", code, str(ROOT / "scripts")],
                                env=dict(os.environ, PYTHONHASHSEED=seed),
                                capture_output=True, text=True, timeout=30, check=True)
        results.append(json.loads(result.stdout))
    assert results[0] == results[1]
