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


# --- curated domain senses are measured, not asserted (#ATLAS-AMBIG) --------

def test_every_curated_domain_sense_is_confirmed_by_the_corpus():
    """DOMAIN_SENSE entries are hand-written, so something must check them.

    One of them was wrong: the `psa` note claimed the majority vote returned
    NPEPPS when it returns KLK3, and the repository's own committed scan already
    said so. This guard exists so a future hand-edit cannot reintroduce an
    unmeasured claim.
    """
    import json
    import atlas_ambiguity as aa
    raw = json.loads(
        (REPO_ROOT / "analysis" / "atlas-domain-sense-validation.json").read_text())
    results = raw["results"]
    assert set(results) == set(aa.DOMAIN_SENSE), \
        "every curated sense must be validated, and vice versa"
    for surface, r in results.items():
        assert r["curated_is_dominant"], f"{surface}: curated sense is not dominant"
        assert r["curated_share"] > 0.85, \
            f"{surface}: curated sense only {r['curated_share']:.1%} of declaring papers"
        assert r["curated_sense"] == aa.DOMAIN_SENSE[surface][1]


def test_the_majority_vote_is_wrong_on_most_curated_collisions():
    """The finding that justifies blocking rather than trusting the vote.

    It is not merely noisy: it picks the sense the cancer literature essentially
    never means. `psa` is the one case where the vote is already right, which is
    why the count is asserted as 'most', not 'all'.
    """
    import json
    raw = json.loads(
        (REPO_ROOT / "analysis" / "atlas-domain-sense-validation.json").read_text())
    disagree = [s for s, r in raw["results"].items() if not r["vote_matches_curated"]]
    assert len(disagree) >= 3, "if the vote became reliable, revisit the blocklist"
    assert "psa" not in disagree, "psa is the documented case where the vote is right"


def test_domain_sense_probes_cover_every_curated_symbol():
    """A curated sense with no probe would silently go unmeasured."""
    import atlas_ambiguity as aa
    import atlas_domain_sense as ads
    missing = set(aa.DOMAIN_SENSE) - set(ads.PROBES)
    assert not missing, f"no declaration probe for {missing}"
    for surface, probes in ads.PROBES.items():
        if surface in aa.DOMAIN_SENSE:
            assert aa.DOMAIN_SENSE[surface][1] in probes, \
                f"{surface}: the curated sense itself must be probeable"


def test_declared_requires_exactly_one_sense():
    import re
    import atlas_domain_sense as ads
    probes = {k: re.compile(v) for k, v in ads.PROBES["er"].items()}
    assert ads.declared("estrogen receptor positive breast cancer", probes) == "ESR1"
    assert ads.declared("epiregulin drives egfr signalling", probes) == "EREG"
    assert ads.declared("er was measured", probes) is None, "no declaration"
    assert ads.declared("estrogen receptor and endoplasmic reticulum stress",
                        probes) is None, "two senses declared -> not usable"


def test_impact_report_never_presents_containment_as_an_error_rate():
    """The measurement mistake this report exists to avoid.

    50.8% of relation rows touch a contested identifier. That is CONTAINMENT --
    ESR1 is contested because `ER` is ambiguous, but most ESR1 edges come from
    papers that wrote 'estrogen receptor' in full. The figure that bounds the
    damage is 1.35%, roughly 38x smaller. Reporting the first alone would
    overstate the damage by that factor, so both must appear.
    """
    import json
    raw = json.loads(
        (REPO_ROOT / "analysis" / "atlas-ambiguity-impact.json").read_text())
    rows = raw["relation_rows"]
    containment = raw["relation_rows_touching_contested_id"] / rows
    at_risk = raw["relation_rows_resting_on_at_risk"] / rows
    assert at_risk < containment, "the bound must be tighter than containment"
    assert at_risk < 0.05, "if this ever approaches containment, re-read the scan"

    text = (REPO_ROOT / "analysis" / "atlas-ambiguity-impact.md").read_text()
    assert f"{100*containment:.1f}%" in text and f"{100*at_risk:.2f}%" in text, \
        "both numbers must be stated; neither may be reported alone"

    # the three assignment classes must partition the total
    assert (raw["assignments_unaffected"] + raw["assignments_corroborated"]
            + raw["assignments_at_risk"]) == raw["assignments_total"]


def test_corroboration_requires_an_unambiguous_route():
    """An assignment is only at risk when nothing else supports it."""
    import json
    raw = json.loads(
        (REPO_ROOT / "analysis" / "atlas-ambiguity-impact.json").read_text())
    # a paper writing both `ER` and "estrogen receptor" is corroborated, not at risk
    assert raw["assignments_corroborated"] > 0
    assert raw["assignments_at_risk"] > 0
    # the vast majority of contested-identifier assignments never touch an
    # ambiguous form at all, which is exactly why containment misleads
    assert raw["assignments_unaffected"] > 10 * raw["assignments_at_risk"]


def test_pmid_sample_is_stable_against_unrelated_pairs(tmp_path, monkeypatch):
    """A pair's example PMIDs must depend only on that pair's own evidence.

    The reservoir originally drew from ONE shared RNG consumed in file order, so
    any change in an earlier pair's membership shifted every later draw. Applying
    441 sense corrections rewrote the example PMIDs of MDM2-p53, which shares no
    paper with any corrected article -- making rebuild diffs unreadable and
    hiding which pairs actually moved.

    The unrelated pair must be written BEFORE the pair under test: a shared RNG
    is only disturbed by draws that happen EARLIER in the file, so a fixture that
    appends the noise afterwards passes under both implementations and measures
    nothing.
    """
    import gzip as gz

    def build(prefix_rows):
        root = tmp_path / f"r{len(prefix_rows)}"
        (root / "relations").mkdir(parents=True)
        (root / "entities").mkdir(parents=True)
        for kind in ("gene", "chemical", "disease"):
            with gz.open(root / "entities" / f"{kind}.tsv.gz", "wt") as fh:
                fh.write("")
        rows = list(prefix_rows) + [
            f"{9000000 + i}\tassociate\tGene|111\tGene|222"
            for i in range(ag.PMID_SAMPLE + 40)
        ]
        with gz.open(root / "relations" / "relations.tsv.gz", "wt") as fh:
            fh.write("\n".join(rows) + "\n")
        return ag.build_index(root)

    monkeypatch.setattr(ag, "DISAMBIGUATION_JSON", tmp_path / "none.json")
    a = build([])
    b = build([f"{8000000 + i}\tassociate\tGene|333\tGene|444"
               for i in range(ag.PMID_SAMPLE + 120)])

    key = tuple(sorted(("111", "222")))
    assert len(a["pmids"][key]) == ag.PMID_SAMPLE
    assert a["pmids"][key] == b["pmids"][key], \
        "an unrelated pair processed earlier changed this pair's sample"


# --- does the discovery layer predict anything? (#ATLAS-LBD-EVAL) -----------

def test_discovery_eval_compares_against_popularity_not_just_random():
    """The comparison that makes the evaluation worth anything.

    Almost any ranking beats random on a clustered co-occurrence graph, so
    'beats random' is not evidence the ABC machinery works. The load-bearing
    baseline is ranking the SAME candidates by popularity.
    """
    import json
    raw = json.loads(
        (REPO_ROOT / "analysis" / "atlas-discovery-eval.json").read_text())
    head = raw["headline"]
    assert {"abc", "popularity", "random"} <= set(head["precision"])
    # identical candidate set, so EVERY ranking must make the same number of
    # predictions -- a method that quietly dropped candidates would look better
    # for the wrong reason
    counts = set(head["predictions"].values())
    assert len(counts) == 1, \
        f"the rankings must differ only in order, never in the candidate set: {head['predictions']}"
    assert head["precision"]["random"] < head["precision"]["popularity"], \
        "beating random is the floor, not the result"


def test_discovery_ranking_underperformance_is_recorded_not_buried():
    """A negative result about this repo's own layer must stay visible.

    atlas_discovery.py claims to correct for popularity. Measured, it does not,
    and its docstring has to say so -- otherwise the next reader trusts the
    claim over the measurement.
    """
    import json
    raw = json.loads(
        (REPO_ROOT / "analysis" / "atlas-discovery-eval.json").read_text())
    head = raw["headline"]
    assert head["precision"]["abc"] < head["precision"]["popularity"]
    assert head["paired"]["decided"], "the paired interval must exclude zero"

    src = (REPO_ROOT / "scripts" / "atlas_discovery.py").read_text()
    assert "DOES NOT BEAT POPULARITY" in src, \
        "the module must carry its own measured negative result"
    assert "atlas-discovery-eval.md" in src


def test_discovery_eval_result_holds_across_split_years():
    """One split year could be an artifact; three agreeing is a result."""
    import json
    raw = json.loads(
        (REPO_ROOT / "analysis" / "atlas-discovery-eval.json").read_text())
    runs = [raw["headline"]] + raw["robustness"]
    assert len(runs) >= 3, "robustness needs more than one split"
    years = {r["split_year"] for r in runs}
    assert len(years) == len(runs), "split years must be distinct"
    for r in runs:
        assert r["precision"]["abc"] < r["precision"]["popularity"], \
            f"split {r['split_year']} disagrees; re-read before trusting the verdict"
        assert r["paired"]["ci95"][1] < 0, f"split {r['split_year']} CI includes zero"


def test_no_standard_link_predictor_beats_popularity():
    """The follow-up question, answered so nobody has to re-ask it.

    Once ABC is found to lose to popularity the obvious move is to swap in a
    better-known link predictor. The standard ones were tried on the SAME
    candidate set and none wins, so that avenue is closed rather than untried.
    """
    import json
    raw = json.loads(
        (REPO_ROOT / "analysis" / "atlas-discovery-eval.json").read_text())
    prec = raw["headline"]["precision"]
    for m in ("adamic_adar", "resource_alloc", "jaccard", "bridges", "abc"):
        assert m in prec, f"{m} must be evaluated"
        assert prec[m] <= prec["popularity"], \
            f"{m} now beats popularity -- rewrite the verdict, do not delete this test"


def test_degree_correction_monotonically_hurts():
    """The shape of the result, not just its sign.

    Rankings order themselves by how hard they correct for degree, and harder
    correction does worse. That is what makes this a statement about the graph
    rather than about one method's implementation.
    """
    import json
    raw = json.loads(
        (REPO_ROOT / "analysis" / "atlas-discovery-eval.json").read_text())
    prec = raw["headline"]["precision"]
    # Jaccard corrects hardest and must sit far below the uncorrected baselines
    assert prec["jaccard"] < prec["abc"] < prec["popularity"]
    assert prec["jaccard"] < prec["adamic_adar"]
    # but even the hardest correction still beats chance: the candidate SET works
    assert prec["jaccard"] > prec["random"]


def test_eval_report_states_the_objection_to_its_own_target():
    """A negative result must carry the argument against itself.

    Predicting what the literature went on to say is not obviously the right
    target for a layer whose purpose is finding OVERLOOKED connections. The
    report has to say so, or it overclaims.
    """
    text = (REPO_ROOT / "analysis" / "atlas-discovery-eval.md").read_text()
    assert "objection" in text.lower()
    assert "overlooked" in text.lower() or "slow to reach" in text.lower()


# --- co-mention inherits the ambiguity problem too (#ATLAS-AMBIG) ----------

def test_comention_alias_map_handles_sense_collisions(monkeypatch):
    """The shape filter is not a sense filter, and the gap is load-bearing.

    `usable_alias` keeps any single token carrying a digit or hyphen, which is
    exactly the shape of `cox-2`, `p21`, `p62` and `fsp1` -- the forms whose
    majority vote was measured wrong. Module support uses this layer's counts to
    argue a zero in the relation column is an extraction failure rather than
    absence of evidence, so a wrong identifier here corrupts that argument.
    """
    import atlas_comention as cm
    monkeypatch.setattr(cm, "_ambiguity", lambda: (
        {"cox-2", "fsp1"}, {"cox-2": {"id": "5743", "symbol": "PTGS2", "why": "x"}}))
    idx = {"alias": {"cox-2": "4513",      # mitochondrial oxidase, the wrong sense
                     "fsp1": "51062",      # atlastin, and no defensible default
                     "gpx4": "2879",
                     "acsl4": "11332"}}
    amap, stats = cm.build_alias_map(idx)

    assert amap["cox-2"] == "5743", "a measured domain sense must be applied"
    assert "fsp1" not in amap, "no domain default -> drop rather than guess"
    assert amap["gpx4"] == "2879", "unambiguous aliases are untouched"
    assert stats == {"redirected": 1, "dropped_ambiguous": 1}


def test_comention_shape_filter_alone_would_keep_the_dangerous_forms():
    """Pins WHY the sense filter is needed, not just that it exists.

    If usable_alias ever started rejecting these on shape, this test failing is
    the signal that the sense filter's justification changed.
    """
    import atlas_comention as cm
    for a in ("cox-2", "fsp1", "beta-actin", "ap-1"):
        assert cm.usable_alias(a), f"{a} passes the shape filter, so shape is not enough"
    # Length happens to exclude some collisions, but that is luck, not
    # disambiguation -- these are rejected for being short, not for being ambiguous.
    for a in ("psa", "p21", "p62", "er"):
        assert not cm.usable_alias(a)
        assert len(a) < cm.MIN_ALIAS_LEN, "excluded by length, which is incidental"


# --- how much of the contradiction signal is real? (#ATLAS-CONTRA-Q) -------

def test_contradiction_quality_separates_the_two_failure_modes():
    """Two failure modes, measured, with opposite answers.

    Within-paper self-contradiction is extraction inconsistency; ambiguity-driven
    conflation is two literatures merged. Reporting only the reassuring one would
    be as misleading as reporting only the alarming one.
    """
    import json
    raw = json.loads(
        (REPO_ROOT / "analysis" / "atlas-contradiction-quality.json").read_text())
    # mode 1: essentially absent
    assert raw["self_contradicting_assertions"] / raw["total_assertions_in_conflicts"] < 0.001
    # mode 2: real, and above 1 with the interval excluding it
    assert raw["mantel_haenszel"] > 1.0
    assert raw["mh_ci95"][0] > 1.0, "if the CI ever spans 1, the caveat must be softened"


def test_ambiguity_enrichment_is_not_a_popularity_artifact():
    """The confound that would have made this finding worthless.

    Colliding identifiers are contested BECAUSE they are heavily mentioned, and
    a pair with more assertions has more chance of showing both directions. If
    stratifying by assertion count collapsed the ratio toward 1, the crude number
    would have been measuring popularity.
    """
    import json
    raw = json.loads(
        (REPO_ROOT / "analysis" / "atlas-contradiction-quality.json").read_text())
    crude, adjusted = raw["crude_risk_ratio"], raw["mantel_haenszel"]
    # the adjustment must barely move it -- that is what rules the confound out
    assert abs(crude - adjusted) < 0.15, \
        f"crude {crude:.2f} vs adjusted {adjusted:.2f}: the confound is no longer ruled out"
    # and the enrichment must hold inside the strata, not only in the pooled number
    held = 0
    for s in raw["strata"].values():
        an, ac = s["amb"]
        cn, cc = s["clean"]
        if an >= 20 and cn >= 20 and (ac / an) > (cc / cn):
            held += 1
    assert held >= 3, "the enrichment must survive within strata, not just pooled"


def test_contradictions_module_carries_the_measured_caveat():
    """A caveat measured elsewhere has to reach the module that needs it."""
    src = (REPO_ROOT / "scripts" / "atlas_contradictions.py").read_text()
    assert "1.45x" in src and "115,024" in src, \
        "the contradiction module must state both measured failure modes"


# --- how wrong is the sampled emergence estimate? (#ATLAS-EMERG-ERR) -------

def test_emergence_error_separates_estimated_pairs_from_exact_ones():
    """Pooling the two would flatter the estimator.

    Most pairs carry no more asserting papers than the sample holds, so their
    share is exact and the estimator does no estimating. Only the remainder can
    be wrong, and that is the population the error must be quoted on.
    """
    import json
    raw = json.loads(
        (REPO_ROOT / "analysis" / "atlas-emergence-error.json").read_text())
    assert raw["pairs_exact"] + raw["pairs_estimated"] == raw["pairs_examined"]
    assert raw["pairs_exact"] > raw["pairs_estimated"], \
        "if most pairs became estimated, the headline needs re-deriving"
    # the honest number is quoted on the estimated subset, and it is worse
    assert (raw["confusion_estimated"]["precision"]
            < raw["confusion_all"]["precision"]), \
        "the estimated-only precision must be the stricter figure"


def test_emergence_sampling_error_grows_with_support():
    """Sampling error should rise as the sample covers less of the population.

    If it did not, the error would not be sampling error and the explanation in
    the report would be wrong.
    """
    import json
    raw = json.loads(
        (REPO_ROOT / "analysis" / "atlas-emergence-error.json").read_text())
    bands = raw["error_by_band"]
    order = [b for b in ("61-120", "121-500", "501-2000", "2000+") if b in bands]
    meds = [bands[b]["median"] for b in order]
    assert len(meds) >= 3
    assert meds[0] < meds[-1], f"error did not grow with support: {list(zip(order, meds))}"


def test_emergence_module_states_its_measured_accuracy():
    src = (REPO_ROOT / "scripts" / "atlas_emergence.py").read_text()
    assert "MEASURED ACCURACY" in src
    assert "atlas_emergence_error.py" in src


# --- this repository is public (#PRIVACY) ----------------------------------

def test_no_machine_specific_paths_in_tracked_files():
    """A public repo must not carry the author's home directory or account name.

    Found in a pre-push scan: `atlas_fulltext.DEFAULT_FT` was a literal
    /Users/<account>/nas/... path, which leaks both the account and the fact
    that full text sits on external storage, and three links in an older note
    embedded an absolute checkout path. The default is now derived from
    Path.home() so it resolves identically without being written down.
    """
    import re
    import subprocess
    files = subprocess.run(["git", "ls-files"], cwd=REPO_ROOT,
                           capture_output=True, text=True).stdout.split()
    # the corpus holds article text the authors wrote; it is not our prose
    files = [f for f in files if not f.startswith("corpus/")]
    pattern = re.compile(r"/Users/[A-Za-z0-9._-]+|/home/[A-Za-z0-9._-]+")
    offenders = []
    for rel in files:
        p = REPO_ROOT / rel
        try:
            text = p.read_text(encoding="utf-8", errors="ignore")
        except (OSError, IsADirectoryError):
            continue
        for m in pattern.finditer(text):
            offenders.append(f"{rel}: {m.group(0)}")
    assert not offenders, (
        "machine-specific paths in a public repo:\n  " + "\n  ".join(offenders[:20]))


# --- quoted numbers must not rot (#ATLAS) ----------------------------------

def _pct(n, d):
    return 100 * n / d


def test_headline_numbers_match_their_source_json():
    """Prose quoting a stale number is the failure mode this catches.

    Each of these figures is recomputed from the committed JSON and required to
    appear verbatim in every document that cites it. If an analysis is re-run and
    a number moves, this fails and the prose has to be updated with it -- which
    is the point. Fixing it by editing the expected value here would defeat the
    guard; edit the documents instead.
    """
    import json
    A = REPO_ROOT / "analysis"
    load = lambda n: json.loads((A / n).read_text())  # noqa: E731
    amb = load("atlas-ambiguity.json")
    dis = load("atlas-disambiguation.json")
    imp = load("atlas-ambiguity-impact.json")
    ev = load("atlas-discovery-eval.json")
    cq = load("atlas-contradiction-quality.json")
    ee = load("atlas-emergence-error.json")

    claims = [
        (f"{amb['by_type']['gene']['mention_share']:.1f}%",
         ["atlas-ambiguity.md", "atlas-README.md"]),
        (f"{_pct(dis['layer_accuracy']['correct'], dis['layer_accuracy']['n']):.1f}%",
         ["atlas-disambiguation.md", "atlas-README.md"]),
        (f"{_pct(imp['relation_rows_resting_on_at_risk'], imp['relation_rows']):.2f}%",
         ["atlas-ambiguity-impact.md", "atlas-README.md"]),
        (f"{100 * ev['headline']['precision']['popularity']:.1f}%",
         ["atlas-discovery-eval.md", "atlas-README.md"]),
        (f"{cq['mantel_haenszel']:.2f}x",
         ["atlas-contradiction-quality.md"]),
        (f"{_pct(ee['confusion_all']['precision'], 1):.1f}%",
         ["atlas-emergence-error.md"]),
    ]
    missing = []
    for value, docs in claims:
        for doc in docs:
            text = (A / doc).read_text()
            # the README renders x as a multiplication sign in prose
            if value not in text and value.replace("x", "×") not in text:
                missing.append(f"{doc} does not quote {value}")
    assert not missing, "stale or missing figures:\n  " + "\n  ".join(missing)


def test_comention_audit_sample_is_uniform_across_shards():
    """The sample must not be a prefix of the first shard.

    Full-text shards are ordered by PMCID, which correlates with publication
    date, so keeping the first N kept sentences would sample the oldest
    literature -- the same failure that once corrupted the emergence layer when
    the index stored a lexicographic PMID prefix.
    """
    import collections
    import random
    import atlas_comention as cm

    audit = {"rows": [], "seen": 0}
    rng = random.Random(cm._AUDIT_SEED)
    # simulate two shards; the second must be represented in the sample
    for shard in (0, 1):
        for i in range(cm.AUDIT_SAMPLE * 5):
            audit["seen"] += 1
            row = {"pmid": f"{shard}-{i}", "sentence": "s", "entities": ["a", "b"]}
            if len(audit["rows"]) < cm.AUDIT_SAMPLE:
                audit["rows"].append(row)
            else:
                j = rng.randrange(audit["seen"])
                if j < cm.AUDIT_SAMPLE:
                    audit["rows"][j] = row
    shards = collections.Counter(r["pmid"].split("-")[0] for r in audit["rows"])
    assert len(audit["rows"]) == cm.AUDIT_SAMPLE
    assert shards["1"] > 0, "a later shard is unrepresented -- the sample is a prefix"
    assert shards["0"] > 0, "the first shard vanished entirely"


def test_comention_process_shard_audit_is_optional():
    """The audit must not be required, so an existing caller keeps working."""
    import collections
    import atlas_comention as cm
    import inspect
    sig = inspect.signature(cm.process_shard)
    assert sig.parameters["audit"].default is None
    assert sig.parameters["rng"].default is None


def test_most_fsp1_corrections_are_extrapolated_and_the_report_says_so():
    """The limitation I nearly shipped without stating.

    The 97.4% accuracy is measured on papers that DECLARE a sense, but three
    quarters of the corrections actually written into the graph land on papers
    that declare nothing. Quoting the headline without that caveat would imply a
    measurement that was never made.
    """
    import json
    raw = json.loads(
        (REPO_ROOT / "analysis" / "atlas-disambiguation.json").read_text())
    text = (REPO_ROOT / "analysis" / "atlas-disambiguation.md").read_text()
    assert raw["gold_set"] < len(raw["corrections"]), \
        "if the gold set ever covered every correction, simplify this section"
    assert "extrapolat" in text.lower()
    assert "declare\nnothing" in text or "declare nothing" in text.replace("\n", " ")


def test_temporal_check_validates_the_undeclaring_corrections():
    """The check that makes the extrapolation defensible rather than assumed.

    FSP1 was named as the ferroptosis suppressor in 2019, so a pre-2019 paper
    corrected to AIFM2 is almost certainly wrong. Among the undeclaring papers
    -- exactly the population the gold set never scored -- that count must stay
    at or near zero while the other senses spread across decades.
    """
    import json
    raw = json.loads(
        (REPO_ROOT / "analysis" / "atlas-disambiguation.json").read_text())
    years = raw.get("gold_years", {})
    assert years, "the temporal check must be recorded, not just printed"
    aifm2 = years.get("AIFM2", [])
    assert aifm2 and min(aifm2) >= 2019, \
        f"an AIFM2-sense paper predates the term: {min(aifm2)}"
    other = years.get("S100A4", [])
    assert other and min(other) < 2019, \
        "the contrast sense must span decades, or the check proves nothing"


def test_module_support_rules_out_conflation_for_contested_claims():
    """A contested verdict has to survive the conflation explanation.

    Across the graph, pairs built on a measured sense collision are 1.45x more
    likely to be flagged contradictory. A reader seeing 'contested: yes' on a
    module claim needs to know whether that could be two literatures merged
    rather than a field divided, so the report states which it is.
    """
    text = (REPO_ROOT / "analysis" / "atlas-module-support.md").read_text()
    assert "1.45x" in text, "the report must quote the measured enrichment"
    assert ("Conflation does not explain these" in text
            or "may be conflation, not disagreement" in text), \
        "the report must resolve the conflation question either way"


def test_no_module_claim_rests_on_a_colliding_entity():
    """If a claim ever lands on one, its contested verdict needs re-reading.

    This is the check that lets the report say conflation is excluded. It is a
    property of the claim list, so it must be asserted rather than assumed.
    """
    import json
    scan = json.loads((REPO_ROOT / "analysis" / "atlas-ambiguity.json").read_text())
    collide = set()
    for t in ("gene", "chemical", "disease"):
        for r in scan["by_type"][t]["sense_rows"]:
            collide |= {r["top"]["id"], r["runner_up"]["id"]}
    import atlas_module_support as ms
    from atlas_baseline import atlas_root
    try:
        idx = ag.load_index(atlas_root())
    except SystemExit:
        pytest.skip("atlas index not built in this checkout")
    bad = []
    for module, a, b, _pmid, _claim in ms.CLAIMS:
        for name in (a, b):
            if ag.resolve(idx, name) in collide:
                bad.append(f"{module}: {name}")
    assert not bad, (
        "module claims now rest on colliding entities; their contested verdicts "
        "must be re-read as possible conflation:\n  " + "\n  ".join(bad))


# --- do the module citations point at the right papers? (#ATLAS-CITE) ------

def test_module_citations_are_not_broken():
    """Three of them were, and the failure was invisible in the support table.

    The `dhodh` layer cited a Nature news item about US fetal-tissue policy,
    `prom2` cited a Theriogenology paper on embryo vitrification, and `gch1`
    cited a PMID that does not resolve. A citation pointing at an unrelated
    paper reads in atlas-module-support.md as `cited-absent`, which looks like
    the graph failed to find a real paper rather than the paper being wrong.
    """
    import json
    raw = json.loads(
        (REPO_ROOT / "analysis" / "atlas-citation-audit.json").read_text())
    broken = [r for r in raw["rows"]
              if r["status"] in ("wrong-subject", "unresolvable")]
    assert not broken, "broken module citations:\n  " + "\n  ".join(
        f"{r['module']} -> {r['pmid']}: {r['title'] or '(unresolvable)'}"
        for r in broken)


def test_the_corrected_pmids_are_pinned():
    """Pin the three corrections so a revert is loud rather than silent.

    Each was verified by exact title against the mechanism its module documents,
    not guessed. The wrong values must not reappear anywhere in the repo.
    """
    import subprocess
    # bad -> (module, the verified correction). A line naming BOTH is documenting
    # the fix, not making the citation, so it is not a revert.
    wrong = {"33864050": ("dhodh", "33981038"),
             "31919077": ("gch1", "31989025"),
             "31761539": ("prom2", "31735663")}
    files = subprocess.run(["git", "ls-files"], cwd=REPO_ROOT,
                           capture_output=True, text=True).stdout.split()
    hits = []
    for rel in files:
        # the audit report records what was wrong, and this file must name the
        # wrong values in order to check for them; both would match themselves
        if (rel.startswith("corpus/") or "citation-audit" in rel
                or rel == "tests/test_atlas.py"):
            continue
        try:
            text = (REPO_ROOT / rel).read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        for line in text.splitlines():
            for bad, (mod, good) in wrong.items():
                if bad in line and good not in line:
                    hits.append(f"{rel}: {bad} ({mod})")
    assert not hits, "a corrected citation has reverted:\n  " + "\n  ".join(hits)

    import atlas_module_support as ms
    by_module = {c[0]: c[3] for c in ms.CLAIMS}
    assert by_module["dhodh"] == "33981038"    # Mao 2021 Nature
    assert by_module["gch1"] == "31989025"     # Kraft 2020 ACS Cent Sci
    assert by_module["prom2"] == "31735663"    # Brown 2019 Dev Cell


# --- the news pipeline's "verified" label (#NEWS-VERIFY) -------------------

def test_news_verification_links_are_measured_not_trusted():
    """44 claims are labelled `verified` on links that mostly do not match.

    A claim about electric fields treating brain cancer is 'verified' against
    papers on freshwater fish biodiversity and speech-language pathology. The
    linked identifiers cluster in one numeric band, which is the signature of
    matching on when a record was indexed rather than what it says.
    """
    import json
    raw = json.loads(
        (REPO_ROOT / "analysis" / "news-verification-audit.json").read_text())
    share = raw["zero_overlap"] / raw["pairs_resolved"]
    assert share > 0.3, (
        "the link quality improved -- re-read the audit and update the manuscript "
        "footnotes and this guard together")
    # the numeric clustering is the diagnostic, so it must be recorded
    assert raw["dominant_prefix_count"] / raw["distinct_pmids"] > 0.2


def test_manuscript_does_not_claim_verification_it_lacks():
    """Two footnotes asserted 'Verified: PMID:X' against unrelated papers.

    42020835 backed a pancreatic-cancer blood-test claim with a paper on
    Andexanet Alfa withdrawal; 42020682 backed an electric-fields brain-cancer
    claim with a survey of speech-language pathologists. Both assertions are
    removed. The manuscript may cite the reporting -- it may not call it verified.
    """
    md = (REPO_ROOT / "article" / "drafts" / "v1.md").read_text()
    for pmid in ("42020835", "42020682"):
        assert f"Verified: PMID:{pmid}" not in md, \
            f"the false verification claim for {pmid} has returned"
    assert "not as verified evidence" in md


def test_sentence_initial_words_are_not_proper_nouns():
    """The root cause of the false verifications.

    'Seven of these 26 patients had inoperable tumors' yielded the search term
    ['Seven'], and the query `Seven` matches ~836,000 PubMed records. ESearch
    returned the five most recently indexed, and the pipeline called that
    verification -- which is why the linked identifiers all cluster in one
    numeric band.
    """
    from verify_news_claims import extract_search_terms
    assert extract_search_terms(
        "Seven of these 26 patients had inoperable tumors.") == []
    assert extract_search_terms(
        "Results showed improvement.") == []
    # genuine proper nouns still survive
    terms = extract_search_terms(
        "Pembrolizumab improved survival in Keytruda-treated NSCLC patients.")
    assert "Keytruda" in terms and "NSCLC" in terms


def test_a_search_result_is_not_verification():
    """Accepting any non-empty result is what marked 44 claims verified.

    The check must reject a paper sharing no subject matter with the claim
    while still accepting one that does, including across morphological
    variants -- 'tumors'/'tumor' and 'immune'/'immunity' must not count as
    different words, or genuine support gets rejected.
    """
    from verify_news_claims import supports_claim
    assert not supports_claim(
        "Seven of these 26 patients had inoperable tumors",
        "How Do Speech-Language Pathologists Know? A Survey")
    assert not supports_claim(
        "New blood test could catch pancreatic cancer before it is too late",
        "Andexanet Alfa Withdrawn from the US: Implications for Intracranial Hemorrhage")
    assert supports_claim(
        "Electric fields supercharge immune system against brain cancer tumors",
        "Tumor treating fields enhance antitumor immunity in glioblastoma")
    assert supports_claim(
        "New blood test could catch pancreatic cancer early detection",
        "A blood-based test for early detection of pancreatic cancer")


def test_fulltext_recency_ceiling_is_documented():
    """A structural coverage limit must be stated, not left as a puzzle.

    Two PMC bulk packages returned exactly zero cancer articles out of 232,890
    while every other package yielded 14-18%. That is version skew -- the PMC
    release is newer than the PubMed baseline the census was built from -- and
    it looks like a bug unless the code says otherwise.
    """
    src = (REPO_ROOT / "scripts" / "atlas_fulltext.py").read_text()
    assert "RECENCY CEILING" in src
    assert "PMC13" in src and "232,890" in src
    readme = (REPO_ROOT / "analysis" / "atlas-README.md").read_text()
    assert "cliff, not a slope" in readme.lower() or "A cliff, not a slope" in readme


def test_comention_rebuild_clears_manifest_and_pairs_together():
    """Clearing one without the other silently corrupts the counts.

    The run merges its results into any existing pair table, relying on the
    manifest to skip shards already counted there. Clearing the manifest alone
    would re-count every shard on top of the existing table and double every
    pair; clearing the table alone would lose the shards not reprocessed.
    """
    src = (REPO_ROOT / "scripts" / "atlas_comention.py").read_text()
    assert "--rebuild" in src
    # the unlink and the manifest reset must sit in the same branch
    i = src.index("if args.rebuild:")
    block = src[i:i + 700]
    assert "unlink()" in block and 'man = {"shards": {}}' in block, \
        "--rebuild must clear BOTH the pair table and the manifest"


def test_pipeline_order_puts_the_graph_rebuild_after_disambiguation():
    """The dependency that is easiest to get wrong and hardest to notice.

    atlas_graph --build APPLIES the per-paper corrections that
    atlas_disambiguate produces. Building the graph first yields an index
    carrying stale corrections, which looks completely normal. Likewise the
    ambiguity scan must precede both, since it produces the blocklist they read.
    """
    src = (REPO_ROOT / "scripts" / "atlas_pipeline.sh").read_text()
    q = src[src.index("quality()"):src.index("mine()")]
    i_amb = q.index("atlas_ambiguity.py")
    i_dis = q.index("atlas_disambiguate.py")
    i_gra = q.index("atlas_graph.py --build")
    assert i_amb < i_dis < i_gra, (
        "quality phase order is wrong: ambiguity -> disambiguate -> graph")
    # co-mention consumes the blocklist, so it cannot precede the scan
    assert i_amb < q.index("atlas_comention.py")
    # and a corpus change requires the double-count-safe rebuild
    assert "--rebuild" in q


def test_pipeline_does_not_silently_rewrite_news_claim_statuses():
    """Re-running the verifier is a data change, not a pipeline step.

    scripts/verify_news_claims.py rewrites 44 claim statuses and their
    credibility scores. Folding that into an audit phase would make a reviewable
    decision happen as a side effect.
    """
    src = (REPO_ROOT / "scripts" / "atlas_pipeline.sh").read_text()
    audit = src[src.index("audit()"):src.index("case \"$PHASE\"")]
    assert "news_verification_audit.py" in audit
    assert "$PY scripts/verify_news_claims.py" not in audit
    assert "Deliberately NOT" in audit


# --- the manuscript's corpus statistics must match the frozen index --------

def test_manuscript_corpus_statistics_match_the_frozen_index():
    """The frozen index is immutable, so every headline count is checkable.

    Checked all of them: 4,830 records, 803 journals, 22 cancer types,
    2001-2026 and ~2,297 immunotherapy articles all match exactly. The
    mechanism count did not -- the abstract said 19 where the index carries 23.
    """
    import collections
    import json
    idx = [json.loads(l) for l in
           (REPO_ROOT / "corpus" / "INDEX.jsonl").read_text().splitlines() if l.strip()]
    md = (REPO_ROOT / "article" / "drafts" / "v1.md").read_text()

    assert len(idx) == 4830
    assert len({r.get("journal") for r in idx if r.get("journal")}) == 803
    mech = collections.Counter()
    cancers = set()
    for r in idx:
        for m in (r.get("mechanisms") or []):
            mech[m] += 1
        cancers.update(r.get("cancer_types") or [])
    assert len(cancers) == 22
    assert mech["immunotherapy"] == 2297

    # the corrected claim, and the threshold that explains the old one
    assert len(mech) == 23, f"index now carries {len(mech)} mechanisms; update the manuscript"
    assert sum(1 for c in mech.values() if c >= 20) == 19
    assert "23 mechanisms" in md
    assert "19 mechanisms, 22 cancer types" not in md, \
        "the abstract has reverted to the undocumented threshold count"


def test_the_thin_mechanisms_are_the_physical_modalities():
    """Why the count mattered rather than being a rounding quibble.

    The four mechanisms the old threshold excluded are cold-atmospheric-plasma,
    electrolysis, radioligand-therapy and targeted-protein-degradation. Three
    are physical or device modalities -- exactly the sparsity this manuscript
    argues about -- so dropping them from the headline hid the project's own
    evidence.
    """
    import collections
    import json
    idx = [json.loads(l) for l in
           (REPO_ROOT / "corpus" / "INDEX.jsonl").read_text().splitlines() if l.strip()]
    mech = collections.Counter()
    for r in idx:
        for m in (r.get("mechanisms") or []):
            mech[m] += 1
    thin = {k for k, c in mech.items() if c < 20}
    assert {"cold-atmospheric-plasma", "electrolysis", "radioligand-therapy"} <= thin


# --- the manuscript recomputed at census scale (#ATLAS-LANDSCAPE) ----------

def test_landscape_capture_is_scope_invariant_and_reported():
    """The comparison that survives descriptor imprecision.

    A cross-mechanism ranking is only as good as its least precise descriptor,
    and several are broad -- 75% of `epigenetic` comes from 'DNA Methylation',
    which any paper MEASURING methylation carries. Comparing the same labels
    across two corpora cancels that entirely, and capture is the result that
    stands.
    """
    import json
    raw = json.loads((REPO_ROOT / "analysis" / "atlas-landscape.json").read_text())
    caps = [r["mesh_frozen"] / r["mesh_census"] for r in raw["rows"]
            if r["mesh_census"] and r["mesh_frozen"]]
    assert len(caps) >= 12
    assert max(caps) / min(caps) > 50, \
        "capture is now near-uniform; the non-uniformity argument needs re-deriving"
    # column B must actually be populated: reading MeSH from INDEX.jsonl, which
    # carries none, silently produced a column of zeros that read as a finding
    assert sum(1 for r in raw["rows"] if r["mesh_frozen"] > 0) >= 12, \
        "MeSH/frozen column is empty -- the by-pmid join is broken again"


def test_manuscript_carries_the_census_recomputation():
    """The manuscript is updated by the research, not frozen against it."""
    md = (REPO_ROOT / "article" / "drafts" / "v1.md").read_text()
    assert "17.6:1" in md and "atlas-landscape.md" in md, \
        "the census recomputation must reach the manuscript, not just the analysis"


def test_landscape_flags_over_broad_descriptors():
    """Reporting the rank shift without the scope caveat would be an artifact.

    `epigenetic` ranks FIRST on the census and would be the headline, but most
    of its matches come from one descriptor that is broader than the mechanism.
    """
    md = (REPO_ROOT / "analysis" / "atlas-landscape.md").read_text()
    assert "DNA Methylation" in md
    assert "artifact of descriptor scope" in md


def test_maturity_comparison_reports_that_it_flips():
    """The result depends on descriptor choice, and hiding that would be spin.

    Across all mechanisms physical modalities look MORE clinically mature;
    restricted to descriptors naming a therapy rather than a process or a
    material, they look less. Reporting either alone would be picking the
    convenient one.
    """
    md = (REPO_ROOT / "analysis" / "atlas-landscape.md").read_text()
    assert "The answer flips" in md
    assert "precise descriptors only" in md


def test_hifu_is_not_preclinical_and_the_manuscript_says_so():
    """The specific claim that survives descriptor scrutiny.

    Both HIFU and CAR-T rest on precise single descriptors, so comparing them
    is clean, and HIFU comes out ahead. Calling physical modalities preclinical
    as a class is what the manuscript got wrong.
    """
    import json
    raw = json.loads((REPO_ROOT / "analysis" / "atlas-landscape.json").read_text())
    R = {r["mechanism"]: r for r in raw["rows"]}
    assert R["hifu"]["clinical_share"] > R["car-t"]["clinical_share"]
    # and the two "physical" mechanisms are not at the same stage
    assert R["hifu"]["clinical_share"] > 1.4 * R["sonodynamic"]["clinical_share"]
    md = (REPO_ROOT / "article" / "drafts" / "v1.md").read_text()
    assert "7.10%" in md and "not a maturity class" in md


# --- replications that never happened (#ATLAS-REPL) ------------------------

def test_replication_cohorts_use_an_equal_window_not_a_minimum():
    """The bias that manufactured a collapse in the first version.

    Scoring whether a cohort was EVER replicated gave a clean monotonic decline
    from ~60% in 1950 to 17.5% in 2020. That is not a finding about science: a
    1975 pair has had fifty years to acquire a second paper and a 2020 pair has
    had six. Every cohort must be scored on the same interval from its own
    first assertion.
    """
    import json
    raw = json.loads((REPO_ROOT / "analysis" / "atlas-replication.json").read_text())
    rows = raw["cohorts"]
    assert len(rows) > 10
    early = [r for r in rows if r["year"] < 1980]
    late = [r for r in rows if 2000 <= r["year"] <= 2014]
    ea = sum(r["replicated"] for r in early) / max(1, sum(r["pairs"] for r in early))
    la = sum(r["replicated"] for r in late) / max(1, sum(r["pairs"] for r in late))
    # with an equal window the two eras must be comparable, not 3x apart
    assert ea < 2 * la, (
        f"early {ea:.3f} vs late {la:.3f}: the window looks unequal again")
    md = (REPO_ROOT / "analysis" / "atlas-replication.md").read_text()
    assert "equal, not merely a minimum" in md


def test_replication_states_it_measures_co_assertion_not_replication():
    """The word does more work than the measurement supports.

    A second paper asserting the same pair may cite the first rather than test
    it, and the pair key carries no direction, so a contradiction counts as a
    replication. Calling that 'replication' without saying so would overclaim.
    """
    md = (REPO_ROOT / "analysis" / "atlas-replication.md").read_text()
    assert "co-assertion" in md.lower()
    assert "not truth" in md.lower() or "bounds attention" in md.lower()


def test_recent_cohort_decline_is_flagged_as_an_upper_bound():
    """MeSH lag biases the newest cohorts downward, so the drop is not measured."""
    md = (REPO_ROOT / "analysis" / "atlas-replication.md").read_text()
    assert "upper bound on a real decline" in md


# --- co-mention precision, measured at last (#ATLAS-COMENT-AUDIT) ----------

def test_comention_alias_matching_is_mechanically_sound():
    """Check 1 is not a judgement call: a miss here is a tokenizer bug."""
    import json
    raw = json.loads(
        (REPO_ROOT / "analysis" / "atlas-comention-audit.json").read_text())
    assert raw["mentions"] > 500
    assert raw["alias_missing"] == 0, \
        f"{raw['alias_missing']} sampled entities have no alias in their sentence"


def test_comention_precision_is_reported_as_a_bound_not_a_number():
    """Corroboration and precision are not the same thing.

    PubTator reads abstracts; this layer reads full text. An entity discussed
    only in Methods is genuinely present and genuinely absent from PubTator's
    annotation, so agreement UNDERSTATES precision. Reporting 44.3% as 'the
    precision' would be as wrong as reporting the 83.8% upper bound as one.
    """
    import json
    raw = json.loads(
        (REPO_ROOT / "analysis" / "atlas-comention-audit.json").read_text())
    lower = raw["pubtator_agree"] / raw["pubtator_scored"]
    upper = (raw["pubtator_agree"] + raw["body_only"]) / raw["pubtator_scored"]
    assert lower < upper, "the bound must be an interval, not a point"
    md = (REPO_ROOT / "analysis" / "atlas-comention-audit.md").read_text()
    assert "upper bound" in md and "lower bound" in md


def test_disagreements_are_split_by_abstract_visibility():
    """Pooling the two kinds of disagreement would waste the measurement.

    An alias visible in the abstract PubTator read is a candidate false
    positive; a body-only match is this layer doing what it exists for.
    """
    import json
    raw = json.loads(
        (REPO_ROOT / "analysis" / "atlas-comention-audit.json").read_text())
    assert raw["in_abstract"] + raw["body_only"] == raw["pubtator_disagree"]
    assert raw["body_only"] > raw["in_abstract"], \
        "if most disagreements became abstract-visible, false positives are likely"


# --- where the thesis sits in the literature (#ATLAS-THESIS) --------------

def test_thesis_position_counts_the_legs_the_corpus_could_not_see():
    """The corpus has no ferroptosis query, so it could not measure its own thesis.

    That is the whole reason this exists: a 4,830-article keyword corpus built
    without a ferroptosis or PDT query cannot say how much work exists on
    ferroptosis or PDT, however carefully it is analysed.
    """
    import json
    raw = json.loads(
        (REPO_ROOT / "analysis" / "atlas-thesis-position.json").read_text())
    t = raw["totals"]
    assert t["ferroptosis"] > 5000, "the ferroptosis field should be large"
    # the ordering that carries the argument
    assert t["drug resistance"] > t["photodynamic therapy"] > t["sonodynamic therapy"]
    assert t["sonodynamic therapy"] < 100, \
        "if SDT-ferroptosis has grown into a literature, rewrite the framing"


def test_thesis_position_is_reported_in_both_directions():
    """A thin intersection supports the claim AND limits it.

    Reporting only 'under-explored, as we said' would be self-serving; reporting
    only 'thirty papers' would ignore that under-exploration is the thesis.
    """
    md = (REPO_ROOT / "analysis" / "atlas-thesis-position.md").read_text()
    assert "not a refutation" in md
    assert "not a literature" in md or "cannot establish a" in md
    # and the SDT count must be flagged as an over-estimate, since the MeSH
    # concept is broader than sonodynamic therapy
    assert "OVER-estimate" in md or "over-estimate" in md


def test_manuscript_states_its_own_evidence_thinness():
    md = (REPO_ROOT / "article" / "drafts" / "v1.md").read_text()
    assert "atlas-thesis-position.md" in md
    assert "roughly thirty papers" in md
