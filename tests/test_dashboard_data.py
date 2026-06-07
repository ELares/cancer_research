"""Tests for the dashboard's data layer (#354).

Covers the pure aggregation/filter helpers (CI-safe; the Streamlit app itself is a
UI-only, non-pinned dependency not exercised here) plus a smoke check against the
committed corpus index.
"""

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import dashboard_data as dd  # noqa: E402

SYNTH = [
    {"pmid": "1", "year": 2020, "mechanisms": ["immunotherapy", "car-t"], "cancer_types": ["lung"],
     "evidence_level": "phase3-clinical"},
    {"pmid": "2", "year": 2021, "mechanisms": ["immunotherapy"], "cancer_types": ["breast", "lung"],
     "evidence_level": "preclinical-invivo"},
    {"pmid": "3", "year": 2021, "mechanisms": ["nanoparticle"], "cancer_types": ["breast"],
     "evidence_level": ""},
    {"pmid": "4", "year": None, "mechanisms": [], "cancer_types": [], "evidence_level": None},
]


def test_value_counts_list_and_scalar():
    mech = dd.value_counts(SYNTH, "mechanisms")
    assert mech["immunotherapy"] == 2 and mech["car-t"] == 1 and mech["nanoparticle"] == 1
    # ordered by frequency (most common first)
    assert list(mech)[0] == "immunotherapy"
    ev = dd.value_counts(SYNTH, "evidence_level")
    assert ev == {"phase3-clinical": 1, "preclinical-invivo": 1}  # empty/None excluded


def test_year_histogram_sorted_ints_only():
    h = dd.year_histogram(SYNTH)
    assert h == {2020: 1, 2021: 2}  # None year excluded; sorted
    assert list(h) == sorted(h)


def test_mechanism_cancer_matrix():
    m = dd.mechanism_cancer_matrix(SYNTH)
    assert m[("immunotherapy", "lung")] == 2
    assert m[("car-t", "lung")] == 1
    assert m[("nanoparticle", "breast")] == 1
    # top-N restriction keeps only the most frequent mechanism
    m1 = dd.mechanism_cancer_matrix(SYNTH, top_mech=1, top_cancer=1)
    assert set(mech for mech, _ in m1) == {"immunotherapy"}


def test_filter_records():
    # mechanism filter (OR within filter)
    assert {r["pmid"] for r in dd.filter_records(SYNTH, mechanisms=["car-t"])} == {"1"}
    # cancer + year AND across filters
    f = dd.filter_records(SYNTH, cancer_types=["lung"], year_range=(2021, 2026))
    assert {r["pmid"] for r in f} == {"2"}
    # evidence-level filter
    assert {r["pmid"] for r in dd.filter_records(SYNTH, evidence_levels=["phase3-clinical"])} == {"1"}
    # year filter drops records with no parseable year (record 4)
    assert all(r["pmid"] != "4" for r in dd.filter_records(SYNTH, year_range=(2000, 2030)))


def test_summary_stats():
    s = dd.summary_stats(SYNTH)
    assert s["n_records"] == 4
    assert s["n_mechanisms"] == 3
    assert s["n_cancer_types"] == 2
    assert s["n_evidence_tagged"] == 2
    assert s["year_min"] == 2020 and s["year_max"] == 2021


# --- committed-corpus smoke ---


def test_load_real_index_and_summary():
    recs = dd.load_index()
    assert len(recs) > 4000  # ~4830
    s = dd.summary_stats(recs)
    assert s["n_mechanisms"] >= 19
    assert s["year_min"] >= 2000 and s["year_max"] <= 2030
    # immunotherapy is the most-studied mechanism (a stable corpus fact)
    assert list(dd.value_counts(recs, "mechanisms"))[0] == "immunotherapy"
