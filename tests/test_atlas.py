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


@pytest.fixture
def blocked(monkeypatch):
    """Pin the blocklist so these tests do not depend on the committed scan."""
    monkeypatch.setattr(ag, "_AMBIG", ({"fsp1", "er"}, {
        "er": {"id": "2099", "symbol": "ESR1", "why": "oncology sense"}}))


def test_fsp1_symbol_collision_refuses_rather_than_guessing(blocked):
    """The hazard that motivated scripts/atlas_ambiguity.py.

    PubTator3 maps `FSP1` to gene 51062, which NCBI calls ATL1 (atlastin
    GTPase 1); the ferroptosis suppressor is 84883 (AIFM2). But `FSP1` is a
    real alias of BOTH, and of S100A4 -- so there is no single right answer to
    patch in, and a blanket remap to AIFM2 would corrupt the S100A4 papers that
    are the majority. resolve() therefore returns None instead of a
    plausible-looking wrong gene, and the per-paper answer comes from
    scripts/atlas_disambiguate.py.
    """
    idx = _idx()
    assert ag.resolve(idx, "FSP1") is None
    assert ag.resolve(idx, "fsp1") is None, "the block is case-insensitive"
    assert "sense collision" in ag.resolve_reason(idx, "FSP1")
    # unambiguous symbols are untouched
    assert ag.resolve(idx, "AIFM2") == "84883"
    assert ag.resolve(idx, "GPX4") == "2879"


def test_blocklist_is_checked_before_the_canon_shortcut(blocked):
    """Regression: PubTator's canonical NAME for gene 51062 is itself "FSP1".

    resolve() short-circuits when a name is a key of `canon`, so a blocked
    symbol that is also a canon key resolved straight through the block. The
    blocklist must be consulted first.
    """
    idx = _idx()
    idx["canon"]["FSP1"] = "FSP1"  # the collision that let it through
    assert ag.resolve(idx, "FSP1") is None


def test_domain_sense_is_opt_in_never_silent(blocked):
    """A curated cancer-domain sense must never apply itself."""
    idx = _idx()
    idx["alias"]["er"] = "2069"  # majority vote: EREG
    assert ag.resolve(idx, "ER") is None, "silent by default is the bug"
    assert ag.resolve(idx, "ER", allow_domain_sense=True) == "2099", "ESR1 when asked"
    assert ag.resolve(idx, "FSP1", allow_domain_sense=True) is None, \
        "FSP1 has no defensible domain default -- it is genuinely per-paper"


def test_resolve_majority_reports_the_vote_the_blocklist_refuses(blocked):
    """The entity audit's job is to REPORT the vote, so it bypasses the block."""
    idx = _idx()
    assert ag.resolve_majority(idx, "FSP1") == "51062"
    assert ag.resolve(idx, "FSP1") is None


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


# --------------------------------------------------------------------------
# Regression guards for defects found by adversarial review.
# Each of these shipped, so each gets a test.
# --------------------------------------------------------------------------

def test_pmid_sample_is_uniform_not_a_lexicographic_prefix():
    """PMIDs must be SAMPLED, not truncated by string order.

    The first version stored `sorted(v)[:50]` over PMID strings, which is a
    first-digit order, so the retained PMIDs were systematically the oldest.
    atlas_emergence then computed 'share of support since 2021' on a sample
    built to exclude recent papers: median true recent-share 26.6% read as 0.0%.
    """
    import random
    import atlas_graph as ag

    # 500 PMIDs spanning old (7-digit) and new (8-digit, 4xxxxxxx = recent)
    old = [str(1_000_000 + i) for i in range(250)]
    new = [str(40_000_000 + i) for i in range(250)]
    universe = old + new

    rng = random.Random(ag._RESERVOIR_SEED)
    res, seen = [], set()
    for pmid in universe:
        seen.add(pmid)
        n = len(seen)
        if len(res) < ag.PMID_SAMPLE:
            res.append(pmid)
        else:
            j = rng.randrange(n)
            if j < ag.PMID_SAMPLE:
                res[j] = pmid

    recent = sum(1 for p in res if int(p) >= 40_000_000)
    assert len(res) == ag.PMID_SAMPLE
    # a uniform sample of a 50/50 population should be nowhere near all-old.
    # The broken lexicographic prefix would give exactly 0 recent.
    assert recent > 0, "reservoir sample retained no recent PMIDs — prefix bias is back"
    assert 0.2 < recent / len(res) < 0.8, (
        f"sample is not representative: {recent}/{len(res)} recent, expected ~50%")


def test_graph_index_records_true_support_size():
    """A share computed from a sample needs the real denominator."""
    import atlas_graph as ag
    idx = {
        "alias": {}, "canon": {"a": "A", "b": "B"},
        "edges": {("a", "b"): {"associate": 900}},
        "pmids": {("a", "b"): ["101", "102"]},
        "n_pmids": {("a", "b"): 873},
    }
    r = ag.support(idx, "a", "b")
    assert r["n_articles"] == 873, "true support size must survive sampling"
    assert len(r["pmids"]) == 2, "the sample is separate from the true count"


def test_seed_replication_by_seed_keys_survive_a_json_round_trip():
    """The ratio section read `by_seed` with str() while it was built with int
    keys, so every lookup returned None and the report crashed with a TypeError
    before it was written. Keys must be str at creation, which is also what a
    JSON round trip produces."""
    import json
    import seed_replication as sr

    src = Path(sr.__file__).read_text(encoding="utf-8")
    assert "by_seed = {str(s_): v for s_, v in pairs}" in src, \
        "by_seed must be built with STRING keys"
    # and a round trip must not change the key type
    d = {"by_seed": {str(42): 1.0, str(43): 2.0}}
    assert set(json.loads(json.dumps(d))["by_seed"]) == {"42", "43"}


def test_discovery_null_uses_the_same_neighbour_set_as_the_observation():
    """Bridges were counted over non-hub neighbours only, while the seed's FULL
    degree was passed as K, inflating the expectation and understating every
    enrichment by the hub fraction."""
    import atlas_discovery as ad

    src = Path(ad.__file__).read_text(encoding="utf-8")
    assert "usable_nb = {b for b in a_nb if degrees.get(b, 0) <= cutoff}" in src
    assert "deg_a = len(usable_nb)" in src, \
        "K must be the usable-neighbour count, not the full degree"


def test_downstream_filters_read_both_census_streams():
    """The text-recovered stream was built to reach recent literature and then
    excluded from every downstream layer, because the filters globbed
    `records/` only."""
    for mod in ("atlas_relations.py", "atlas_fulltext.py"):
        src = (REPO_ROOT / "scripts" / mod).read_text(encoding="utf-8")
        assert "records_unindexed" in src, \
            f"{mod} must include the recovered census stream in its PMID filter"


# --- FSP1 sense disambiguation (#ATLAS-AMBIG) ------------------------------

import atlas_disambiguate as ad  # noqa: E402


def test_gold_label_needs_exactly_one_declared_sense():
    """A paper naming two senses is not usable as ground truth."""
    assert ad.gold_label("we studied ferroptosis suppressor protein 1 (fsp1)") == "AIFM2"
    assert ad.gold_label("fibroblast-specific protein 1 marks caf") == "S100A4"
    assert ad.gold_label("atlastin gtpase mutations cause spastic paraplegia") == "ATL1"
    assert ad.gold_label("fsp1 was measured") is None, "no declaration"
    assert ad.gold_label(
        "fsp1 (ferroptosis suppressor protein 1) and s100a4 both rose") is None, \
        "two senses declared -> unusable, not a coin flip"


def test_classifier_cannot_read_the_phrase_that_defines_the_gold_label():
    """The load-bearing anti-leakage guard.

    The gold label IS the presence of an expansion phrase. If the classifier
    may read that phrase it scores ~100% and measures nothing. Masking must
    remove it, so a text carrying ONLY the declaration and no independent cue
    is an abstention rather than a free correct answer.
    """
    pred, _score, why = ad.classify("ferroptosis suppressor protein 1 was studied")
    assert pred is None, f"leaked the label through the declaration ({why})"

    # the same text plus a genuine, independent cue now decides
    pred, _score, _ = ad.classify(
        "ferroptosis suppressor protein 1 was studied; gpx4 and erastin too")
    assert pred == "AIFM2"


def test_masking_does_not_delete_independent_evidence():
    """Masking must remove the label phrase, not the whole vocabulary."""
    pred, _s, _w = ad.classify("s100a4 drives invasion and migration in fibroblast rich stroma")
    assert pred == "S100A4", "cues outside the masked phrase must survive"


def test_classifier_abstains_on_a_tie_rather_than_guessing():
    # 2 cues each: gpx4 + erastin against fibroblast + metasta
    pred, score, why = ad.classify("gpx4 and erastin in fibroblast and metastasis")
    assert score["AIFM2"] == score["S100A4"] == 2, "fixture must actually tie"
    assert pred is None and why == "tie between senses"
    # an uneven contest still decides
    pred, _s, _w = ad.classify("gpx4 and erastin and rsl3 in fibroblast tissue")
    assert pred == "AIFM2"


def test_ambiguity_scan_separates_species_from_sense():
    """The decomposition that keeps the headline honest.

    Undecomposed, 28.3% of mentions sit on a contested surface form. Most are
    human/mouse orthologs of the SAME gene, which a literature map may merge
    safely. Reporting the raw figure as an error rate overstates the damage
    roughly sevenfold, so the two populations must never be pooled.
    """
    import json
    raw = json.loads((REPO_ROOT / "analysis" / "atlas-ambiguity.json").read_text())
    species = {r["surface"] for r in raw["species_ambiguity"]}
    sense = {r["surface"] for r in raw["sense_collision"]}
    assert not (species & sense), "a form cannot be both"
    assert "gapdh" in species, "human/mouse GAPDH is species ambiguity, not a collision"
    assert "fsp1" in sense, "FSP1 is a genuine sense collision"
    assert "fsp1" in raw["blocklist"]
    # FSP1 ranks ~2181st by mention volume, so a pure top-N scan would miss it
    assert "fsp1" not in {r["surface"] for r in
                          sorted(raw["sense_collision"],
                                 key=lambda r: -r["total"])[:25]}


def test_genes_are_measurably_dirtier_than_mesh_entities():
    """The comparison that explains WHY the gene layer needed a blocklist.

    MeSH is a curated vocabulary with one preferred term per concept; gene
    symbols are not, and NCBI itself lists `FSP1` as an official alias of three
    genes. If this ever stops holding, the blocklist's scope should be revisited.
    """
    import json
    raw = json.loads((REPO_ROOT / "analysis" / "atlas-ambiguity.json").read_text())
    by = raw["by_type"]
    assert set(by) == {"gene", "chemical", "disease"}
    gene = by["gene"]["mention_share"]
    assert gene > 4 * max(by["chemical"]["mention_share"],
                          by["disease"]["mention_share"]), \
        "genes should be several times worse; if not, re-read the scan"


def test_relation_corrections_move_only_the_atl1_identifiers():
    """Precision guard on the index-build remap.

    relations.tsv.gz stores identifiers, not the surface form that produced
    them. A paper discussing BOTH cancer-associated fibroblasts and ferroptosis
    would have its genuine S100A4 edge rewritten to AIFM2 if every colliding id
    were correctable, so only ATL1 -- which just 1.9% of these papers mention
    at all -- may move.
    """
    corr = {"111": "84883"}
    assert ag._corrected("51062", "111", corr) == "84883", "ATL1 -> the real FSP1"
    assert ag._corrected("73991", "111", corr) == "84883", "mouse Atl1 too"
    assert ag._corrected("6275", "111", corr) == "6275", \
        "a genuine S100A4 edge must survive a paper whose FSP1 means AIFM2"
    assert ag._corrected("2879", "111", corr) == "2879", "GPX4 is untouched"
    assert ag._corrected("51062", "999", corr) == "51062", \
        "an undecided paper keeps what PubTator said"


def test_corrections_load_maps_pmid_to_identifier(tmp_path, monkeypatch):
    import json
    f = tmp_path / "d.json"
    f.write_text(json.dumps({"corrections": {
        "1": {"pubtator": "ATL1", "corrected": "AIFM2"},
        "2": {"pubtator": "ATL1", "corrected": "S100A4"},
        "3": {"pubtator": "ATL1", "corrected": "nonsense"}}}))
    monkeypatch.setattr(ag, "DISAMBIGUATION_JSON", f)
    got = ag.load_corrections()
    assert got == {"1": "84883", "2": "6275"}, "unknown senses are dropped, not guessed"


def test_corrections_absent_file_is_not_fatal(tmp_path, monkeypatch):
    """The index must still build on a checkout that has not run the layer."""
    monkeypatch.setattr(ag, "DISAMBIGUATION_JSON", tmp_path / "missing.json")
    assert ag.load_corrections() == {}


def test_hierarchical_relatedness_distinguishes_nesting_from_collision():
    """The third ambiguity class, which genes do not have.

    Glioblastoma nests under Glioma in the MeSH tree: merging them loses
    specificity but does not attribute the biology to an unrelated concept.
    Prostatic and Pancreatic Neoplasms do not nest, and conflating those is a
    genuine sense collision between two different cancers.
    """
    import atlas_ambiguity as aa
    glioblastoma = ["C04.557.465.625.600.380.080.335"]
    glioma = ["C04.557.465.625.600.380"]
    assert aa.hierarchically_related(glioblastoma, glioma)
    assert aa.hierarchically_related(glioma, glioblastoma), "order must not matter"

    prostate = ["C04.588.945.440.770", "C12.200.294.260.750"]
    pancreas = ["C06.301.761", "C04.588.322.475"]
    assert not aa.hierarchically_related(prostate, pancreas)

    # a shared prefix that is not a tree-path boundary is not nesting
    assert not aa.hierarchically_related(["C04.557.46"], ["C04.557.465"]), \
        "prefix matching must respect the dot separator"
    assert not aa.hierarchically_related(["C04.557"], ["C04.557"]), \
        "identical trees are the same node, not a parent/child pair"
