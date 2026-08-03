#!/usr/bin/env python3
"""Unit guards for the atlas pipeline (#ATLAS).

The atlas ingests millions of records from network sources into a data root
that is gitignored, so none of it can be exercised in CI directly. What CAN be
guarded, and is guarded here, is the pure logic that decides what enters the
census and how it is interpreted:

  * the cancer definition (C04 vs the adjacent-descriptor extension),
  * the XML record parser,
  * the text-match fallback for articles MeSH has not indexed,
  * entity resolution and the identifier-vs-symbol hazard,
  * the contradiction detector's thresholds.

Every test uses synthetic fixtures. No network, no atlas data, no PubMed.

Run: pytest tests/test_atlas.py -v
"""

import gzip
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import atlas_baseline as ab  # noqa: E402
import atlas_contradictions as ac  # noqa: E402
import atlas_graph as ag  # noqa: E402
from atlas_unindexed import is_cancer_text  # noqa: E402


# --------------------------------------------------------------------------
# The cancer definition
# --------------------------------------------------------------------------

def test_adjacent_descriptors_are_experimental_context_not_topical():
    """The extension exists to catch mechanism papers, not to widen the topic.

    Broad process terms would pull in most of cell biology; Apoptosis is the
    canonical example and must stay out.
    """
    assert "D045744" in ab.ADJACENT_DESCRIPTORS      # Cell Line, Tumor
    assert "D000079403" in ab.ADJACENT_DESCRIPTORS   # Ferroptosis
    assert "D023041" in ab.ADJACENT_DESCRIPTORS      # Xenograft Model Antitumor Assays
    assert "D017209" not in ab.ADJACENT_DESCRIPTORS, "Apoptosis must not be adjacent"


def _article_xml(pmid: str, mesh_uis, title="T", abstract="A", year="2020",
                 pmcid=None, doi=None):
    heads = "".join(
        f'<MeshHeading><DescriptorName UI="{u}" MajorTopicYN="N">D{u}</DescriptorName></MeshHeading>'
        for u in mesh_uis)
    mesh_block = f"<MeshHeadingList>{heads}</MeshHeadingList>" if mesh_uis else ""
    ids = ""
    if doi:
        ids += f'<ArticleId IdType="doi">{doi}</ArticleId>'
    if pmcid:
        ids += f'<ArticleId IdType="pmc">{pmcid}</ArticleId>'
    return (
        "<PubmedArticle><MedlineCitation>"
        f"<PMID Version='1'>{pmid}</PMID>"
        "<Article>"
        "<Journal><Title>J Test</Title><JournalIssue><PubDate>"
        f"<Year>{year}</Year></PubDate></JournalIssue></Journal>"
        f"<ArticleTitle>{title}</ArticleTitle>"
        f"<Abstract><AbstractText>{abstract}</AbstractText></Abstract>"
        "<PublicationTypeList><PublicationType UI='D016428'>Journal Article</PublicationType>"
        "</PublicationTypeList>"
        "</Article>"
        f"{mesh_block}"
        "</MedlineCitation>"
        f"<PubmedData><ArticleIdList>{ids}</ArticleIdList></PubmedData>"
        "</PubmedArticle>")


def _write_baseline(tmp_path: Path, articles) -> Path:
    xml = ("<?xml version='1.0'?><PubmedArticleSet>" + "".join(articles)
           + "</PubmedArticleSet>")
    p = tmp_path / "pubmed26n9999.xml.gz"
    with gzip.open(p, "wt", encoding="utf-8") as fh:
        fh.write(xml)
    return p


def test_parser_keeps_c04_and_adjacent_and_drops_neither_nor(tmp_path):
    c04 = {"D009369": "Neoplasms"}
    path = _write_baseline(tmp_path, [
        _article_xml("1", ["D009369"], doi="10.1/a", pmcid="PMC1"),   # true C04
        _article_xml("2", ["D045744"]),                               # adjacent only
        _article_xml("3", ["D006801"]),                               # neither
        _article_xml("4", []),                                        # unindexed
    ])
    got = {r["pmid"]: r for r in ab.parse_articles(path, c04)}
    assert set(got) == {"1", "2"}, "kept exactly the C04 and adjacent records"
    assert got["1"]["cancer_basis"] == "C04"
    assert got["2"]["cancer_basis"] == "adjacent"
    assert got["1"]["cancer_ui"] == ["D009369"]
    assert got["2"]["cancer_ui"] == [] and got["2"]["adjacent_ui"] == ["D045744"]
    # metadata round-trips
    assert got["1"]["doi"] == "10.1/a" and got["1"]["pmcid"] == "PMC1"
    assert got["1"]["year"] == 2020


def test_parser_reports_unindexed_separately(tmp_path):
    """Articles with no MeSH cannot be classified and must not silently vanish."""
    path = _write_baseline(tmp_path, [
        _article_xml("1", ["D009369"]),
        _article_xml("2", []),
        _article_xml("3", []),
    ])
    total, nomesh = ab.count_articles(path)
    assert (total, nomesh) == (3, 2)


# --------------------------------------------------------------------------
# The text fallback for un-indexed articles
# --------------------------------------------------------------------------

@pytest.mark.parametrize("title,abstract,expected", [
    ("A phase 2 trial in breast cancer", "", True),
    ("Carcinoma of the lung", "", True),
    ("Tumour microenvironment remodelling", "", True),
    ("Ferroptosis is driven by lipid peroxidation", "in tumour cells", True),
    ("Chemotherapy outcomes", "", True),
    ("Cardiac remodelling after myocardial infarction", "left ventricular function", False),
    ("A study of pulsed electromagnetic fields on bone healing", "osteoblasts", False),
    ("mRNA vaccine against SARS-CoV-2", "immunogenicity in healthy adults", False),
])
def test_cancer_text_matcher(title, abstract, expected):
    assert is_cancer_text(title, abstract) is expected


def test_cancer_matcher_excludes_deliberately_broad_words():
    """Words that would raise recall at a real precision cost stay out."""
    for t in ("Cell growth regulation", "A benign lesion", "Tissue mass measurement"):
        assert is_cancer_text(t, "") is False


# --------------------------------------------------------------------------
# Entity resolution
# --------------------------------------------------------------------------

def _idx():
    return {
        "alias": {"gpx4": "2879", "fsp1": "51062", "aifm2": "84883",
                  "ferroptosis suppressor protein 1": "51062"},
        "canon": {"2879": "GPX4", "51062": "ATL1", "84883": "AIFM2"},
        "edges": {}, "pmids": {},
    }


def test_resolve_accepts_symbol_or_identifier():
    idx = _idx()
    assert ag.resolve(idx, "GPX4") == "2879"
    assert ag.resolve(idx, "gpx4") == "2879"
    assert ag.resolve(idx, "2879") == "2879", "an identifier resolves to itself"
    assert ag.resolve(idx, "not-a-gene") is None


def test_fsp1_symbol_collision_is_reproduced_not_silently_fixed():
    """The hazard that motivated scripts/atlas_entity_audit.py.

    PubTator3 maps `FSP1` and `ferroptosis suppressor protein 1` to gene 51062,
    which NCBI calls ATL1 (atlastin GTPase 1). The real FSP1 is 84883 (AIFM2).
    The resolver deliberately reports what the data says rather than patching
    it, so downstream code that queries by SYMBOL gets the wrong gene and the
    audit is what catches it. If this test ever fails because resolution
    changed, re-run scripts/atlas_entity_audit.py before trusting any result.
    """
    idx = _idx()
    assert ag.resolve(idx, "FSP1") == "51062"
    assert idx["canon"]["51062"] == "ATL1"
    assert ag.resolve(idx, "AIFM2") == "84883"
    assert ag.resolve(idx, "FSP1") != ag.resolve(idx, "AIFM2")


def test_support_is_order_independent():
    idx = _idx()
    idx["edges"][("2879", "84883")] = {"negative_correlate": 4}
    idx["pmids"][("2879", "84883")] = ["111", "222"]
    a = ag.support(idx, "GPX4", "AIFM2")
    b = ag.support(idx, "AIFM2", "GPX4")
    assert a["total"] == b["total"] == 4
    assert a["predicates"] == {"negative_correlate": 4}


# --------------------------------------------------------------------------
# Contradiction detection
# --------------------------------------------------------------------------

def test_contradiction_ranks_balanced_disputes_above_lopsided_ones():
    """A 12-vs-9 split is a live disagreement; 50-vs-1 is a settled claim with
    one outlier. Ranking by the WEAKER side is what encodes that."""
    idx = {
        "canon": {"a": "A", "b": "B", "c": "C", "d": "D"},
        "pmids": {},
        "edges": {
            ("a", "b"): {"positive_correlate": 12, "negative_correlate": 9},
            ("c", "d"): {"positive_correlate": 50, "negative_correlate": 1},
        },
    }
    direction, _valence = ac.scan(idx)
    assert [r["a_name"] for r in direction] == ["A"], \
        "the lopsided pair falls below the MIN_WEAK threshold entirely"
    assert direction[0]["weaker"] == 9
    assert direction[0]["balance"] == pytest.approx(9 / 12)


def test_contradiction_thresholds_suppress_thin_evidence():
    idx = {
        "canon": {"a": "A", "b": "B"},
        "pmids": {},
        # both sides present but the pair barely discussed
        "edges": {("a", "b"): {"positive_correlate": 2, "negative_correlate": 2}},
    }
    direction, _ = ac.scan(idx)
    assert direction == [], "below MIN_WEAK and MIN_TOTAL, so not a reportable dispute"


def test_valence_conflict_detects_treats_and_causes():
    idx = {
        "canon": {"x": "drugX", "y": "diseaseY"},
        "pmids": {},
        "edges": {("x", "y"): {"treat": 20, "cause": 11, "associate": 4}},
    }
    _direction, valence = ac.scan(idx)
    assert len(valence) == 1
    assert valence[0]["treat"] == 20 and valence[0]["cause"] == 11
    assert valence[0]["weaker"] == 11
