"""Tests for the non-circular mechanism-recall measurement (#412).

The load-bearing guards here are:

* `test_leak_classifier_*` — the leak classifier (`is_keyword_substring`) must
  keep the exact independent/leaky split the headline depends on. If a descriptor
  silently flips from independent to leaky (or back), the headline recall would
  either lose a clean signal or admit a tautological one. This is the CI tripwire.

* `test_leakage_free_path_excludes_mesh` / `test_mesh_only_signal_not_tagged` —
  prove the non-circularity guard itself: the MeSH reference label is invisible to
  the matcher whose recall we measure. If someone reverts the `include_metadata`
  switch in `tag_articles.get_searchable_text`, these fail loudly.

The rest are pure-logic tests of the aggregation, plus light structural checks of
the shipped YAML map and generated artifacts.
"""

import json
import sys
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

from config import MECHANISM_KEYWORDS  # noqa: E402
from tag_articles import get_searchable_text, match_mechanisms  # noqa: E402
import mechanism_recall as mr  # noqa: E402

MAP_PATH = REPO_ROOT / "analysis" / "mesh-mechanism-map.yaml"
REPORT_PATH = REPO_ROOT / "analysis" / "mechanism-recall-report.md"
JSON_PATH = REPO_ROOT / "analysis" / "mechanism-recall.json"


# ---------------------------------------------------------------------------
# Leak classifier — the load-bearing CI guard for the independent/leaky split.
# ---------------------------------------------------------------------------

# (descriptor, mechanism, expected_is_leaky). Hand-verified against the live
# keyword set. A change here means the headline pool composition changed.
_LEAK_CASES = [
    # Independent: descriptor string does NOT contain a mechanism keyword.
    ("Receptors, Chimeric Antigen", "car-t", False),  # word order differs from "chimeric antigen receptor"
    ("Immunotherapy, Adoptive", "car-t", False),
    ("Immunoconjugates", "antibody-drug-conjugate", False),
    ("Antibodies, Bispecific", "bispecific-antibody", False),  # vs "bispecific antibody"
    ("B7-H1 Antigen", "immunotherapy", False),
    ("CTLA-4 Antigen", "immunotherapy", False),  # vs "anti-ctla-4"
    ("DNA Methylation", "epigenetic", False),  # vs "dna methylation cancer"
    ("Glycolysis", "metabolic-targeting", False),  # vs "glycolysis inhibit"
    ("Poly(ADP-ribose) Polymerase Inhibitors", "synthetic-lethality", False),  # vs "parp inhibitor"
    ("Electroporation", "electrochemical-therapy", False),  # vs "electroporation cancer/therapy"
    ("Ultrasonic Therapy", "sonodynamic", False),
    # Leaky: descriptor string DOES contain a mechanism keyword (near-tautological).
    ("Oncolytic Virotherapy", "oncolytic-virus", True),
    ("Immune Checkpoint Inhibitors", "immunotherapy", True),  # contains "immune checkpoint"
    ("Histone Deacetylase Inhibitors", "epigenetic", True),  # contains "histone deacetylase"
    ("Azacitidine", "epigenetic", True),
    ("CRISPR-Cas Systems", "crispr", True),  # contains "crispr"
    ("Gene Editing", "crispr", True),
    ("Nanoparticles", "nanoparticle", True),  # contains "nanoparticle"
    ("Fecal Microbiota Transplantation", "microbiome", True),  # contains "fecal microbiota"
    ("Radiofrequency Ablation", "frequency-therapy", True),
    ("CD47 Antigen", "phagocytosis-checkpoint", True),  # contains "cd47"
    ("mRNA Vaccines", "mrna-vaccine", True),  # contains "mrna vaccine"
]


@pytest.mark.parametrize("descriptor,mechanism,expected_leaky", _LEAK_CASES)
def test_leak_classifier(descriptor, mechanism, expected_leaky):
    assert mr.is_keyword_substring(descriptor, mechanism) is expected_leaky


def test_leak_classifier_uses_live_keywords():
    """The classifier must read the LIVE keyword set, not a frozen copy: adding a
    keyword that matches a descriptor flips it to leaky on the next run."""
    descriptor = "Immunoconjugates"
    assert mr.is_keyword_substring(descriptor, "antibody-drug-conjugate") is False
    saved = MECHANISM_KEYWORDS["antibody-drug-conjugate"]
    try:
        MECHANISM_KEYWORDS["antibody-drug-conjugate"] = saved + ["immunoconjugate"]
        assert mr.is_keyword_substring(descriptor, "antibody-drug-conjugate") is True
    finally:
        MECHANISM_KEYWORDS["antibody-drug-conjugate"] = saved


# ---------------------------------------------------------------------------
# Non-circularity: the MeSH reference label is invisible to the matcher.
# ---------------------------------------------------------------------------


def test_leakage_free_path_excludes_mesh():
    fm = {
        "title": "A cancer study",
        "mesh_terms": ["Oncolytic Virotherapy"],
        "diseases_annotated": ["zztoken_disease"],
        "genes": ["ZZGENE"],
        "drugs": ["zzdrug"],
    }
    body = "## Abstract\n\nThis tumor work is unrelated text.\n"
    with_meta = get_searchable_text(fm, body, include_metadata=True)
    without_meta = get_searchable_text(fm, body, include_metadata=False)
    # MeSH + annotations present with metadata, absent without it.
    for token in ("oncolytic virotherapy", "zztoken_disease", "zzgene", "zzdrug"):
        assert token in with_meta
        assert token not in without_meta
    # Title + abstract survive in both.
    assert "cancer study" in without_meta and "unrelated text" in without_meta


def test_mesh_only_signal_not_tagged():
    """A record whose ONLY mechanism signal is in MeSH must NOT be tagged by the
    leakage-free matcher, but IS tagged by the metadata path — the exact circularity
    the measurement removes."""
    fm = {"title": "A tumor treatment study", "mesh_terms": ["Oncolytic Virotherapy"]}
    body = "## Abstract\n\nWe treated the tumor; no mechanism keyword appears here.\n"
    lf = match_mechanisms(get_searchable_text(fm, body, include_metadata=False), fm["title"].lower())
    meta = match_mechanisms(get_searchable_text(fm, body, include_metadata=True), fm["title"].lower())
    assert "oncolytic-virus" not in lf
    assert "oncolytic-virus" in meta


# ---------------------------------------------------------------------------
# Pure-logic: canonicalization, classification, recall, aggregation.
# ---------------------------------------------------------------------------


def test_canonical_mechanism_resolves_case():
    assert mr.canonical_mechanism("mrna-vaccine") == "mRNA-vaccine"
    assert mr.canonical_mechanism("car-t") == "car-t"
    assert mr.canonical_mechanism("unknown-x") == "unknown-x"


def test_classify_descriptors_splits():
    m = {
        "mechanisms": {
            "car-t": {"descriptors": ["Receptors, Chimeric Antigen", "Immunotherapy, Adoptive"]},
            "oncolytic-virus": {"descriptors": ["Oncolytic Virotherapy"]},
        }
    }
    c = mr.classify_descriptors(m)
    assert c["car-t"]["independent"] == ["Receptors, Chimeric Antigen", "Immunotherapy, Adoptive"]
    assert c["car-t"]["leaky"] == []
    assert c["oncolytic-virus"]["independent"] == []
    assert c["oncolytic-virus"]["leaky"] == ["Oncolytic Virotherapy"]


def test_recall_helper():
    assert mr.recall(3, 4) == 0.75
    assert mr.recall(0, 0) is None
    assert mr.recall(0, 5) == 0.0


def _synthetic_results():
    classified = {
        "big": {"canon": "big", "independent": ["D1"], "leaky": [], "note": "n", "proxy_confounded": False},
        "small": {"canon": "small", "independent": ["D2"], "leaky": [], "note": "", "proxy_confounded": False},
        "proxy": {"canon": "proxy", "independent": ["D3"], "leaky": [], "note": "", "proxy_confounded": True},
        "leakyonly": {"canon": "leakyonly", "independent": [], "leaky": ["D4"], "note": "", "proxy_confounded": False},
    }
    acc = {
        "big": dict(indep_pool=100, indep_hit_lf=90, indep_hit_idx=95, leaky_pool=0, leaky_hit_lf=0),
        "small": dict(indep_pool=10, indep_hit_lf=8, indep_hit_idx=8, leaky_pool=0, leaky_hit_lf=0),
        "proxy": dict(indep_pool=200, indep_hit_lf=50, indep_hit_idx=50, leaky_pool=0, leaky_hit_lf=0),
        "leakyonly": dict(indep_pool=0, indep_hit_lf=0, indep_hit_idx=0, leaky_pool=50, leaky_hit_lf=48),
    }
    m = {"unmeasurable": {"foo": "reason"}}
    return mr.build_results(classified, acc, m, n_records=500, n_with_mesh=400)


def test_build_results_partitions_measurable():
    r = _synthetic_results()
    # Only "big" clears N_MIN and is non-proxy.
    assert set(r["measurable"]) == {"big"}
    assert r["measurable"]["big"]["recall_leakage_free"] == 0.90
    assert r["measurable"]["big"]["recall_production_index"] == 0.95


def test_build_results_routes_subfloor_and_proxy():
    r = _synthetic_results()
    assert "small" in r["leaky_only"] and "N_MIN" in r["leaky_only"]["small"]["reason"]
    assert "proxy" in r["leaky_only"] and "proxy" in r["leaky_only"]["proxy"]["reason"]
    assert "leakyonly" in r["leaky_only"]
    assert r["leaky_only"]["leakyonly"]["recall_leaky_near_tautological"] == 48 / 50


def test_build_results_aggregates_over_measurable_only():
    r = _synthetic_results()
    agg = r["aggregates"]
    # Only "big" is measurable, so both aggregates equal its recall.
    assert agg["n_measurable_mechanisms"] == 1
    assert agg["total_independent_pool"] == 100
    assert agg["volume_weighted_recall_leakage_free"] == 0.90
    assert agg["macro_recall_leakage_free"] == 0.90


def test_build_results_passes_through_unmeasurable():
    r = _synthetic_results()
    assert r["unmeasurable"] == {"foo": "reason"}


# ---------------------------------------------------------------------------
# Structural: shipped YAML map + generated artifacts.
# ---------------------------------------------------------------------------


def _load_map():
    return yaml.safe_load(MAP_PATH.read_text(encoding="utf-8"))


def test_map_mechanisms_resolve_to_real_tags():
    m = _load_map()
    canon_keys = {k.lower() for k in MECHANISM_KEYWORDS}
    for mech in m["mechanisms"]:
        assert mech.lower() in canon_keys, f"{mech} is not a known mechanism tag"


def test_map_no_mechanism_listed_twice():
    m = _load_map()
    overlap = set(m["mechanisms"]) & set(m.get("unmeasurable") or {})
    assert not overlap, f"mechanisms in both measurable and unmeasurable: {overlap}"


def test_committed_json_is_self_consistent():
    data = json.loads(JSON_PATH.read_text(encoding="utf-8"))
    agg = data["aggregates"]
    measurable = data["measurable"]
    assert agg["n_measurable_mechanisms"] == len(measurable)
    # Every headline mechanism clears the floor and is a valid fraction.
    for mech, e in measurable.items():
        assert e["n_independent"] >= data["n_min"]
        assert 0.0 <= e["recall_leakage_free"] <= 1.0
        assert not e["proxy_confounded"]
    # Volume-weighted aggregate recomputes from the per-mechanism counts.
    tot_pool = sum(e["n_independent"] for e in measurable.values())
    tot_hit = sum(round(e["recall_leakage_free"] * e["n_independent"]) for e in measurable.values())
    assert agg["total_independent_pool"] == tot_pool
    if tot_pool:
        assert abs(agg["volume_weighted_recall_leakage_free"] - tot_hit / tot_pool) < 1e-9


def test_committed_report_headline_matches_json():
    data = json.loads(JSON_PATH.read_text(encoding="utf-8"))
    report = REPORT_PATH.read_text(encoding="utf-8")
    n = data["aggregates"]["n_measurable_mechanisms"]
    assert f"Measurable mechanisms (independent pool >= {data['n_min']}): **{n}**" in report
