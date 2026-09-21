"""Offline paired comparisons must preserve articles, ties, and support limits.

The synthetic counts below have known answers. They test the analysis without
requiring any particular scientific result from the frozen corpus.
"""

import copy
import json
import math
import sys
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import paired_partner_agreement as agreement  # noqa: E402


def _snapshot(records, vocabulary=("f", "p", "q", "r", "s", "z")):
    return {
        "schema_version": 1,
        "vocabulary": sorted(vocabulary),
        "records": records,
        "descriptors": {},
        "descriptor_metadata": {},
        "provenance": {},
    }


def _record(pmid, keyword, mesh):
    return {"pmid": str(pmid), "keyword": keyword, "mesh": mesh}


def _selection_snapshot():
    records = []
    for partner, n in (("p", 15), ("q", 10), ("r", 5)):
        for _ in range(n):
            records.append(_record(len(records) + 1, ["f", partner], ["f", partner]))
    for _ in range(30):
        records.append(_record(len(records) + 1, ["f", "s"], []))
        records.append(_record(len(records) + 1, [], ["f", "z"]))
    records.append(_record("998", [], []))
    records.append(_record("999", ["f", "s"], None))
    # Labels outside the common vocabulary cannot add paired edges.
    records[0]["keyword"] = ["f", "outside", "p"]
    return _snapshot(records)


def _row(result, mechanism="f"):
    return next(row for row in result["primary"]["rows"] if row["mechanism"] == mechanism)


def test_average_ranks_and_spearman_use_average_ties():
    assert agreement.average_ranks([4, 1, 1, 3, 0]) == [5, 2.5, 2.5, 4, 1]
    assert agreement.spearman([1, 1, 2, 3], [1, 2, 3, 3]) == pytest.approx(8 / 9)
    assert agreement.spearman([1, 1, 2], [1, 2, 3]) == pytest.approx(math.sqrt(3) / 2)
    assert agreement.spearman([1, 2, 3], [3, 2, 1]) == pytest.approx(-1)


@pytest.mark.parametrize("left,right", [([], []), ([1], [2]), ([1, 1, 1], [1, 2, 3]),
                                         ([1, 2, 3], [0, 0, 0])])
def test_spearman_is_undefined_for_short_or_constant_vectors(left, right):
    assert agreement.spearman(left, right) is None


def test_spearman_retains_zero_partners():
    # Omitting partners absent in one arm would turn this into perfect agreement.
    assert agreement.spearman([9, 4, 0, 0, 0], [0, 4, 9, 0, 0]) == pytest.approx(-1 / 8)


def test_top_three_shares_boundary_ties_without_alphabetic_selection():
    counts = {"a": 10, "b": 7, "c": 5, "d": 5, "e": 0}
    expected = {"a": 1, "b": 1, "c": 0.5, "d": 0.5, "e": 0}
    assert agreement.top_k_membership(counts) == expected
    assert agreement.top_k_membership(dict(reversed(list(counts.items())))) == expected
    renaming = dict(zip(counts, reversed(list(counts))))
    renamed = agreement.top_k_membership({renaming[key]: value for key, value in counts.items()})
    assert {key: renamed[renaming[key]] for key in counts} == expected
    assert agreement.top_k_overlap(counts, {"a": 10, "b": 7, "c": 6, "d": 4, "e": 0}) == pytest.approx(5 / 6)
    assert agreement.top_k_overlap(counts, counts) == pytest.approx(1)


def test_top_three_never_fills_slots_with_zero_count_partners():
    sparse = {"a": 9, "b": 2, "c": 0, "d": 0}
    assert agreement.top_k_membership(sparse) is None
    assert agreement.top_k_overlap(sparse, {"a": 9, "b": 2, "c": 1, "d": 0}) is None
    assert agreement.top_k_membership(dict.fromkeys("abcd", 2)) == dict.fromkeys("abcd", 0.75)


def test_fractional_overlap_is_exactly_reproducible_when_partner_names_change():
    # Shared membership is 1/3 + 1 + 1, hence overlap 7/9. Ordinary float
    # accumulation on Python 3.10 changes its last bit when the fractional
    # term sorts before the two unit terms. Persisted JSON must not depend on
    # label ordering or Python's version-specific built-in sum algorithm.
    left = {"a": 1, "b": 5, "c": 4, "d": 1, "e": 1}
    right = {"a": 3, "b": 5, "c": 4, "d": 0, "e": 0}
    assert agreement.top_k_overlap(left, right) == 7 / 9
    renamed_left = {("z" if key == "a" else key): value for key, value in left.items()}
    renamed_right = {("z" if key == "a" else key): value for key, value in right.items()}
    assert agreement.top_k_overlap(renamed_left, renamed_right) == 7 / 9


def test_primary_pairs_use_same_articles_common_labels_and_all_denominators():
    result = agreement.assemble(_selection_snapshot())
    primary = result["primary"]
    assert primary["n"] == 91
    assert result["missing_mesh"]["n"] == 1
    for arm in ("keyword", "mesh"):
        summary = primary["arms"][arm]
        assert summary["n"] == 91
        assert summary["tagged"] == 60
        assert summary["multi_tagged"] == 60
        assert summary["tagged_denominator"] == 91
        assert summary["multi_tagged_denominator"] == 60
        assert summary["multi_tagged_fraction"] == 1
    pairs = {tuple(pair["mechanisms"]): pair for pair in primary["pairs"]}
    assert len(pairs) == len(primary["pairs"])
    assert all(a < b for a, b in pairs)
    assert all("outside" not in pair for pair in pairs)
    assert pairs[("f", "p")]["keyword"] == pairs[("f", "p")]["mesh"] == 15
    assert pairs[("f", "s")]["keyword"] == 30
    assert pairs[("f", "s")]["mesh"] == 0
    row = _row(result)
    assert row["target_counts"] == {"keyword": 60, "mesh": 60, "both": 30}
    assert row["partner_counts"] == {
        "keyword": {"p": 15, "q": 10, "r": 5, "s": 30, "z": 0},
        "mesh": {"p": 15, "q": 10, "r": 5, "s": 0, "z": 30},
    }
    assert row["positive_partners"] == {"keyword": 4, "mesh": 4}
    assert row["zero_partners"] == {"keyword": 1, "mesh": 1}


def test_shared_focal_sensitivity_changes_selection_before_comparing_partners():
    row = _row(agreement.assemble(_selection_snapshot()))
    assert row["scored"]
    assert row["spearman"] == pytest.approx(-0.6)
    assert row["top_k_overlap"] == pytest.approx(2 / 3)
    shared = row["shared_focal"]
    assert shared["n"] == 30
    assert shared["scored"]
    assert shared["partner_counts"]["keyword"] == shared["partner_counts"]["mesh"] == {
        "p": 15, "q": 10, "r": 5, "s": 0, "z": 0,
    }
    assert shared["spearman"] == pytest.approx(1)
    assert shared["top_k_overlap"] == pytest.approx(1)


def test_each_arm_must_reach_target_support_and_raw_counts_remain_visible():
    records = [_record(i, ["f", "pqr"[i % 3]], ["f", "pqr"[i % 3]]) for i in range(30)]
    records[-1]["keyword"] = ["r"]
    row = _row(agreement.assemble(_snapshot(records)))
    assert agreement.MIN_TARGET == 30
    assert agreement.TOP_K == 3
    assert row["target_counts"]["keyword"] == 29
    assert row["target_counts"]["mesh"] == 30
    assert row["partner_counts"]["keyword"]["r"] == 9
    assert not row["scored"]
    assert row["spearman"] is None and row["top_k_overlap"] is None
    assert row["unscored_reasons"]


def test_three_positive_partners_are_required_even_with_thirty_focal_articles():
    records = [_record(i, ["f", "pq"[i % 2]], ["f", "pq"[i % 2]]) for i in range(30)]
    row = _row(agreement.assemble(_snapshot(records)))
    assert row["target_counts"]["keyword"] == row["target_counts"]["mesh"] == 30
    assert row["positive_partners"] == {"keyword": 2, "mesh": 2}
    assert not row["scored"]
    assert row["spearman"] is None and row["top_k_overlap"] is None
    assert row["unscored_reasons"]
    assert row["shared_focal"]["unscored_reasons"]


def test_shared_focal_support_is_checked_on_its_own_subset():
    snapshot = _selection_snapshot()
    del snapshot["records"][0]
    row = _row(agreement.assemble(snapshot))
    assert row["scored"]
    shared = row["shared_focal"]
    assert shared["n"] == 29
    assert not shared["scored"]
    assert shared["spearman"] is None and shared["top_k_overlap"] is None
    assert shared["partner_counts"]["keyword"]["p"] == 14
    assert shared["unscored_reasons"]


def test_supported_constant_profiles_keep_overlap_and_explain_undefined_rho():
    records = [_record(i, ["f", "pqr"[i % 3]], ["f", "pqr"[i % 3]]) for i in range(30)]
    row = _row(agreement.assemble(_snapshot(records, vocabulary=("f", "p", "q", "r"))))
    assert row["scored"]
    assert row["spearman"] is None
    assert "constant" in row["spearman_reason"].lower()
    assert row["top_k_overlap"] == pytest.approx(1)


@pytest.fixture
def corpus_files(tmp_path):
    index = tmp_path / "INDEX.jsonl"
    articles = tmp_path / "articles"
    articles.mkdir()
    mapping = tmp_path / "map.yaml"
    descriptors = {
        "car-t": ["Receptors, Chimeric Antigen"],
        "immunotherapy": ["Immune Checkpoint Inhibitors"],
        "mrna-vaccine": ["mRNA Vaccines"],
        "nanoparticle": ["Nanoparticles"],
        "frequency-therapy": [],
    }
    mapping.write_text(yaml.safe_dump({"mechanisms": {
        name: {"descriptors": terms} for name, terms in descriptors.items()
    }}))
    rows = [
        {"pmid": "202", "mechanisms": ["immunotherapy", "mRNA-vaccine"]},
        {"pmid": "101", "mechanisms": ["car-t", "car-t", "bioelectric", "frequency-therapy"]},
        {"pmid": "404", "mechanisms": []},
        {"pmid": "303", "mechanisms": ["nanoparticle"]},
    ]
    mesh = {
        "101": ["Receptors, Chimeric Antigen", "Receptors, Chimeric Antigen",
                "Immune Checkpoint Inhibitors"],
        "202": ["Humans"],
        "303": [],
        "404": ["mRNA Vaccines"],
    }
    index.write_text("".join(json.dumps(row) + "\n" for row in rows))
    for row in rows:
        data = {**row, "mesh_terms": mesh[row["pmid"]]}
        (articles / f"{row['pmid']}.md").write_text(
            "---\n" + yaml.safe_dump(data) + "---\nIgnored article body.\n")
    return index, articles, mapping


def test_loader_joins_pmids_and_distinguishes_missing_mesh_from_observed_no_hits(corpus_files):
    snapshot = agreement.load_snapshot(*corpus_files)
    records = {row["pmid"]: row for row in snapshot["records"]}
    assert snapshot["vocabulary"] == ["car-t", "immunotherapy", "mrna-vaccine", "nanoparticle"]
    assert records["101"]["keyword"] == ["bioelectric", "car-t", "frequency-therapy"]
    assert records["101"]["mesh"] == ["car-t", "immunotherapy"]
    assert records["202"]["keyword"] == ["immunotherapy", "mrna-vaccine"]
    assert records["202"]["mesh"] == []
    assert records["303"]["mesh"] is None
    assert records["404"]["keyword"] == []
    assert records["404"]["mesh"] == ["mrna-vaccine"]
    result = agreement.assemble(snapshot)
    assert result["primary"]["n"] == 3
    assert result["missing_mesh"]["n"] == 1
    assert result["missing_mesh"]["mechanism_counts"]["nanoparticle"] == 1
    assert result["primary"]["arms"]["keyword"]["tagged"] == 2
    assert result["primary"]["arms"]["keyword"]["multi_tagged"] == 1
    assert result["primary"]["arms"]["mesh"]["tagged"] == 2
    assert result["primary"]["arms"]["mesh"]["multi_tagged"] == 1
    pairs = {tuple(pair["mechanisms"]): pair for pair in result["primary"]["pairs"]}
    assert pairs[("car-t", "immunotherapy")]["mesh"] == 1


@pytest.mark.parametrize("fault,reason", [
    ("duplicate_index", "duplicate.*PMID"),
    ("missing_article", "missing article"),
    ("extra_article", "extra article"),
    ("wrong_frontmatter_pmid", "PMID mismatch"),
    ("duplicate_frontmatter_pmid", "duplicate.*PMID"),
    ("wrong_frontmatter_labels", "mechanism mismatch"),
    ("bad_yaml", "invalid metadata"),
    ("unclosed_frontmatter", "closing delimiter"),
])
def test_invalid_inputs_fail_without_replacing_existing_reports(corpus_files, tmp_path, monkeypatch, fault, reason):
    index, articles, mapping = corpus_files
    article = articles / "202.md"
    if fault == "duplicate_index":
        index.write_text(index.read_text() + index.read_text().splitlines()[0] + "\n")
    elif fault == "missing_article":
        article.unlink()
    elif fault == "extra_article":
        (articles / "505.md").write_text(article.read_text())
    elif fault in {"wrong_frontmatter_pmid", "duplicate_frontmatter_pmid"}:
        replacement = "505" if fault == "wrong_frontmatter_pmid" else "101"
        article.write_text(article.read_text().replace("'202'", f"'{replacement}'"))
    elif fault == "wrong_frontmatter_labels":
        article.write_text(article.read_text().replace("- immunotherapy", "- sonodynamic"))
    elif fault == "bad_yaml":
        article.write_text("---\npmid: [\n---\n")
    else:
        article.write_text("---\npmid: '202'\n")
    with pytest.raises(ValueError, match=reason):
        agreement.load_snapshot(index, articles, mapping)
    json_output, md_output = tmp_path / "report.json", tmp_path / "report.md"
    json_output.write_bytes(b"existing JSON\n")
    md_output.write_bytes(b"existing Markdown\n")
    monkeypatch.setattr(agreement, "OUT_JSON", json_output)
    monkeypatch.setattr(agreement, "OUT_MD", md_output)
    monkeypatch.setattr(agreement, "INDEX", index)
    monkeypatch.setattr(agreement, "ARTICLES", articles)
    monkeypatch.setattr(agreement, "MAP", mapping)
    assert agreement.main([]) == 2
    assert json_output.read_bytes() == b"existing JSON\n"
    assert md_output.read_bytes() == b"existing Markdown\n"


def test_render_only_recomputes_stale_derived_values_without_loading_corpus(tmp_path, monkeypatch):
    expected = agreement.assemble(_selection_snapshot())
    stale = copy.deepcopy(expected)
    stale["primary"]["n"] = -1
    _row(stale)["top_k_overlap"] = -1
    _row(stale)["shared_focal"]["n"] = -1
    json_output, md_output = tmp_path / "report.json", tmp_path / "report.md"
    json_output.write_text(json.dumps(stale))
    md_output.write_text("stale report\n")
    monkeypatch.setattr(agreement, "OUT_JSON", json_output)
    monkeypatch.setattr(agreement, "OUT_MD", md_output)
    monkeypatch.setattr(agreement, "INDEX", tmp_path / "absent-index.jsonl")
    monkeypatch.setattr(agreement, "ARTICLES", tmp_path / "absent-articles")
    monkeypatch.setattr(agreement, "MAP", tmp_path / "absent-map.yaml")

    def no_corpus_access(*args, **kwargs):
        pytest.fail("--render-only must reconstruct from stored article labels")

    monkeypatch.setattr(agreement, "load_snapshot", no_corpus_access)
    agreement.main(["--render-only"])
    assert json.loads(json_output.read_text()) == expected
    assert md_output.read_text() == agreement.render(expected)


def test_committed_artifacts_reconstruct_and_match_current_frozen_inputs():
    """Guard freshness, without pinning agreement to a preferred outcome."""
    stored = json.loads(agreement.OUT_JSON.read_text())
    assert stored["snapshot"] == agreement.load_snapshot()
    rebuilt = agreement.assemble(stored["snapshot"])
    assert stored == rebuilt
    assert agreement.OUT_MD.read_text() == agreement.render(rebuilt)
