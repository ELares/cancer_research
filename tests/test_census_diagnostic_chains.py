"""Guards for the census diagnostic-therapy chain analysis.

Aggregate agreement is a reproducibility check, not proof that matchers or
cohorts are identical. Annotation effects and rankings are measured separately
on the selected records, including ties and overlapping chain memberships.
"""
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
JSON = REPO / "analysis/census-diagnostic-chains.json"
MD = REPO / "analysis/census-diagnostic-chains.md"
MANUSCRIPT = REPO / "article/drafts/v1.md"
PUBLISHED_CORPUS_MATCHES = 240


@pytest.fixture(scope="module")
def d():
    return json.loads(JSON.read_text())


def test_the_corpus_arm_reproduces_the_published_figure(d):
    """Without this the census column is being compared against nothing."""
    assert d["corpus_matched_production_text"] == PUBLISHED_CORPUS_MATCHES, (
        f"the corpus arm gives {d['corpus_matched_production_text']} where the "
        f"manuscript published {PUBLISHED_CORPUS_MATCHES}; the matcher "
        "or comparator has changed and the reproduction must be checked")


def test_the_annotation_channel_cost_is_measured_not_assumed(d):
    """The one channel the census lacks. Its cost must be a number, and the
    report must state which of the two cases it is in -- an earlier draft said
    the census "reads slightly low" in a sentence that shipped unchanged when
    the measured cost was zero."""
    cost = d["annotation_channel_cost_records"]
    assert cost == (d["corpus_matched_production_text"]
                    - d["corpus_matched_without_annotations"])
    md = MD.read_text()
    if cost == 0:
        assert "does not change the count of records matching at least one chain" in md
    else:
        assert f"a difference of {cost:,}" in md
    assert "do not estimate an annotation effect in the census" in md
    assert "reads low by about that much" not in md


def test_a_census_zero_and_a_corpus_zero_are_reported_apart(d):
    """They are different claims. A corpus zero says the retrieval never
    reached those papers; a census zero would be a statement about the indexed
    literature. Collapsing them is the error the section exists to correct."""
    for chain in d["chains_zero_on_corpus"]:
        row = next(r for r in d["rows"] if r["chain"] == chain)
        assert row["corpus_production_text"] == 0
    for chain in d["chains_zero_on_census"]:
        row = next(r for r in d["rows"] if r["chain"] == chain)
        assert row["census"] == 0
    md = MD.read_text()
    assert ("## The zeros" in md) and ("property of a retrieval" in md
                                       or "claim about the indexed literature" in md)


def test_the_ordering_disagreement_is_derived(d):
    """The section's actual finding is that the ordering moves, and an earlier
    draft would have left a reader to spot that from two columns."""
    # Average the occupied positions for each tied count.
    def positions(key):
        ordered = sorted(d["rows"], key=lambda r: -r[key])
        return {r["chain"]: sum(i for i, peer in enumerate(ordered, 1)
                                if peer[key] == r[key])
                / sum(peer[key] == r[key] for peer in ordered)
                for r in ordered}

    census, corpus = positions("census"), positions("corpus_production_text")
    moved = [c for c in census if abs(census[c] - corpus[c]) >= 3]
    md = MD.read_text()
    assert f"disagrees on {len(moved)} of {len(d['rows'])}" in md
    for c in moved:
        assert (f"`{c}`: corpus rank {corpus[c]:g}, "
                f"census rank {census[c]:g}") in md


def test_the_manuscript_quotes_the_census_column_not_the_corpus_one():
    txt = " ".join(MANUSCRIPT.read_text().split())
    d = json.loads(JSON.read_text())
    assert f"{d['census_matched']:,} records match at least one chain" in txt
    top = max(d["rows"], key=lambda r: r["census"])
    assert f"{top['census']:,}" in txt, (
        "the manuscript does not quote the leading chain's census count")


def test_the_layers_limit_survives_the_scale_up(d):
    """Scale does not turn co-occurrence into causation, and the report must
    keep saying so -- a bigger number reads as a stronger claim unless the
    limit travels with it."""
    md = MD.read_text()
    assert "does not say the paper USES one to select the other" in md


@pytest.fixture
def scanner(tmp_path, monkeypatch):
    import importlib.util
    import sys

    monkeypatch.setenv("FERRO_ATLAS_ROOT", str(tmp_path / "atlas"))
    monkeypatch.syspath_prepend(str(REPO / "scripts"))
    spec = importlib.util.spec_from_file_location(
        "diagnostic_input_test", REPO / "scripts/census_diagnostic_chains.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.RECORDS == tmp_path / "atlas/records"
    monkeypatch.setattr(module, "OUT_JSON", tmp_path / "report.json")
    monkeypatch.setattr(module, "OUT_MD", tmp_path / "report.md")
    module.OUT_JSON.write_text('{"published": true}\n')
    module.OUT_MD.write_text("Published interpretation.\n")
    monkeypatch.setattr(sys, "argv", ["census_diagnostic_chains.py"])
    import tag_articles

    monkeypatch.setattr(tag_articles, "PMID_DIR", tmp_path / "frozen-corpus")
    return module


def _write_input(path, records):
    import gzip

    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        for record in records:
            fh.write(json.dumps(record) + "\n")


def test_missing_census_precedes_corpus_scan(scanner, monkeypatch):
    monkeypatch.setattr(scanner, "scan_corpus_without_annotations",
                        lambda: pytest.fail("corpus scan reached"))
    before = scanner.OUT_JSON.read_bytes(), scanner.OUT_MD.read_bytes()
    with pytest.raises(SystemExit, match="No census records"):
        scanner.main()
    assert (scanner.OUT_JSON.read_bytes(), scanner.OUT_MD.read_bytes()) == before


def test_selected_shards_preserve_chain_counts_and_abstract_channel(scanner):
    record = {"title": "HER2 testing and trastuzumab", "abstract": "Observed response"}
    _write_input(scanner.RECORDS / "a.jsonl.gz", [record])
    _write_input(scanner.RECORDS / "b.jsonl.gz", [record] * 4)
    _write_input(scanner.RECORDS / "c.jsonl.gz", [{"title": "Unmatched article"}])
    result = scanner.scan_census(2)
    assert result["shards"] == 2
    assert result["records"] == 2
    assert result["with_abstract"] == 1
    assert result["matched"] == 1
    assert result["per_chain"] == {"her2-testing-to-trastuzumab": 1}


@pytest.mark.parametrize("state", ["missing", "empty", "unreadable"])
def test_unavailable_corpus_preserves_reports(scanner, state):
    import tag_articles

    _write_input(scanner.RECORDS / "a.jsonl.gz", [{"title": "Unmatched article"}])
    corpus = Path(tag_articles.PMID_DIR)
    if state != "missing":
        corpus.mkdir()
    if state == "unreadable":
        (corpus / "article.md").write_text("No article frontmatter.\n")
    before = scanner.OUT_JSON.read_bytes(), scanner.OUT_MD.read_bytes()
    with pytest.raises(SystemExit, match="No readable frozen-corpus articles"):
        scanner.main()
    assert (scanner.OUT_JSON.read_bytes(), scanner.OUT_MD.read_bytes()) == before


def test_valid_zero_matches_do_not_claim_historical_reproduction(scanner):
    import tag_articles

    _write_input(scanner.RECORDS / "a.jsonl.gz", [{"title": "Unmatched article"}])
    corpus = Path(tag_articles.PMID_DIR)
    corpus.mkdir()
    (corpus / "article.md").write_text("---\ntitle: Unmatched article\n---\n\n")
    assert scanner.main() == 0
    result = json.loads(scanner.OUT_JSON.read_text())
    assert result["census_records"] == result["corpus_records"] == 1
    assert result["census_matched"] == result["corpus_matched_production_text"] == 0
    markdown = scanner.OUT_MD.read_text()
    assert "does not reproduce the historical corpus result" in markdown
    assert "ordering comparison is unavailable" in markdown
    assert "rank the chains the same way" not in markdown
    assert "does not establish absence elsewhere in the indexed literature" in markdown
    assert "this one is a claim about the indexed literature" not in markdown


def test_frozen_corpus_comparison_keeps_the_annotation_channel(scanner):
    import tag_articles

    corpus = Path(tag_articles.PMID_DIR)
    corpus.mkdir()
    (corpus / "article.md").write_text(
        "---\ntitle: HER2 testing\ndrugs:\n- trastuzumab\n---\n\n")
    result = scanner.scan_corpus_without_annotations()
    assert result["records"] == 1
    assert result["matched_production_text"] == 1
    assert result["matched_without_annotations"] == 0
    _write_input(scanner.RECORDS / "a.jsonl.gz", [{"title": "Unmatched article"}])
    markdown = scanner.render(scanner.assemble(scanner.scan_census(1), result))
    assert "from 1 to 0 (a difference of 1, 100.0%" in markdown
    assert "do not estimate an annotation effect in the census" in markdown
    assert "reads low by about that much" not in markdown


def test_render_failure_preserves_both_diagnostic_reports(scanner, monkeypatch):
    import tag_articles

    _write_input(scanner.RECORDS / "a.jsonl.gz", [{"title": "Unmatched article"}])
    corpus = Path(tag_articles.PMID_DIR)
    corpus.mkdir()
    (corpus / "article.md").write_text("---\ntitle: Unmatched article\n---\n\n")

    def broken_render(_):
        raise ValueError("render failure")

    monkeypatch.setattr(scanner, "render", broken_render)
    before = scanner.OUT_JSON.read_bytes(), scanner.OUT_MD.read_bytes()
    with pytest.raises(ValueError, match="render failure"):
        scanner.main()
    assert (scanner.OUT_JSON.read_bytes(), scanner.OUT_MD.read_bytes()) == before


def test_aggregate_agreement_does_not_certify_matcher_or_cohort(d):
    for markdown in (MD.read_text(), MANUSCRIPT.read_text()):
        assert "does not establish an identical matcher or article cohort" in markdown
        assert "so the matcher is the one the manuscript reported" not in markdown


def test_unchanged_any_chain_total_does_not_hide_per_chain_loss(scanner):
    import tag_articles

    corpus = Path(tag_articles.PMID_DIR)
    corpus.mkdir()
    (corpus / "article.md").write_text(
        "---\ntitle: HER2 testing trastuzumab and EGFR mutation\n"
        "drugs:\n- erlotinib\n---\n\n")
    _write_input(scanner.RECORDS / "a.jsonl.gz", [{"title": "Unmatched article"}])
    assert scanner.main() == 0
    result = json.loads(scanner.OUT_JSON.read_text())
    assert result["corpus_matched_production_text"] == 1
    assert result["corpus_matched_without_annotations"] == 1
    assert result["annotation_channel_cost_records"] == 0
    markdown = scanner.OUT_MD.read_text()
    assert "does not change the count of records matching at least one chain (1)" in markdown
    assert "changes 1 per-chain count(s)" in markdown
    assert "`egfr-mutation-to-egfr-inhibitor` 1 to 0" in markdown
    assert "None of the per-chain counts changes" not in markdown


@pytest.mark.parametrize("census,corpus,expected", [
    ([10, 9, 8, 7, 6, 5, 4, 3, 2, 1], [9, 10, 8, 7, 6, 5, 4, 3, 2, 1],
     "The rankings differ, but no chain differs by three or more places."),
    ([1] * 10, [1] * 10, "The two columns rank the chains the same way."),
    ([1] * 10, [2] + [1] * 9, "corpus rank 1, census rank 5.5"),
    ([2] + [1] * 9, [1] * 10, "corpus rank 5.5, census rank 1"),
])
def test_rank_comparisons_handle_small_changes_and_ties(scanner, d, census, corpus, expected):
    from copy import deepcopy

    result = deepcopy(d)
    for row, cen, cor in zip(result["rows"], census, corpus):
        row.update(census=cen, corpus_production_text=cor, corpus_without_annotations=cor)
    result.update(census_matched=sum(census), corpus_matched_production_text=sum(corpus),
                  corpus_matched_without_annotations=sum(corpus), annotation_channel_cost_records=0,
                  chains_zero_on_census=[], chains_zero_on_corpus=[])
    for rows in (result["rows"], list(reversed(result["rows"]))):
        result["rows"] = rows
        markdown = scanner.render(result)
        assert expected in markdown
        assert "equal counts receive their average tied rank" in markdown
        assert "do not by themselves identify their cause" in markdown
        if census != corpus:
            assert "rank the chains the same way" not in markdown
