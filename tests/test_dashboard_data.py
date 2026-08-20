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


# --- census layer (#RETIRE-FROZEN) ----------------------------------------

def test_census_headline_returns_both_denominators_or_neither():
    """Quoting one denominator alone misleads in a KNOWN direction.

    Against the whole census the trial share understates; against classifiable
    records alone it overstates. So the function must never return a shape a
    caller can render half of -- either both, or None.
    """
    design = {"census": 1000, "classifiable": 400, "classes": {"trial": 40}}
    h = dd.census_headline(design)
    assert h["share_of_census"] == 4.0
    assert h["share_of_classifiable"] == 10.0
    assert h["undetermined"] == 600
    # every field a caller needs is present together
    for key in ("census", "classifiable", "undetermined", "trials",
                "share_of_census", "share_of_classifiable"):
        assert key in h


def test_census_headline_refuses_a_partial_artifact():
    """A half-written artifact must render a notice, not a plausible number.

    An artifact missing `classifiable` would let a caller show the census-wide
    share alone, which is exactly the misleading half.
    """
    assert dd.census_headline(None) is None
    assert dd.census_headline({}) is None
    assert dd.census_headline({"census": 10, "classes": {"trial": 1}}) is None
    assert dd.census_headline({"census": 10, "classifiable": 5}) is None


def test_census_rows_order_by_trial_share_not_volume():
    """The ordering is a finding, not a display preference.

    Descriptor breadth varies enormously between mechanisms, so a volume
    ranking is substantially a ranking of how broad each descriptor is. This
    fixture makes the two orderings DISAGREE: the largest mechanism has the
    lowest share, so a volume sort would put it first.
    """
    profile = {"rows": [
        {"mechanism": "big-and-preclinical", "census": 36788, "trials": 182,
         "trial_share": 0.49, "growth": 2.49, "top_sites": [], "top_partners": []},
        {"mechanism": "small-and-clinical", "census": 1352, "trials": 96,
         "trial_share": 7.10, "growth": 1.16,
         "top_sites": [{"site": "cervix/uterus", "enrichment": 6.48, "n": 338}],
         "top_partners": [{"mechanism": "nanoparticle", "n": 37}]},
    ]}
    rows = dd.census_mechanism_rows(profile)
    assert [r["mechanism"] for r in rows] == ["small-and-clinical", "big-and-preclinical"]
    assert rows[0]["top site"] == "cervix/uterus"
    assert rows[0]["top partner"] == "nanoparticle"
    # a mechanism with no site or partner data renders as None, not as a crash
    assert rows[1]["top site"] is None


def test_census_loader_fails_soft_per_artifact():
    """A missing artifact must not take the others down with it.

    Fail-soft per key rather than fail-open: the caller gets None for what is
    absent and real data for what is present, so a panel can say which half it
    is missing instead of silently rendering an empty census as a complete one.
    """
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        out = dd.load_census(base=d)
        assert set(out) == set(dd.CENSUS_ARTIFACTS)
        assert all(v is None for v in out.values())
        # a corrupt file is treated the same as an absent one
        p = Path(d) / dd.CENSUS_ARTIFACTS["design"]
        p.write_text("{not json", encoding="utf-8")
        assert dd.load_census(names=["design"], base=d)["design"] is None


def test_the_committed_census_artifacts_are_loadable():
    """The artifacts the dashboard ships against are actually there.

    Unlike the census records themselves, these are COMMITTED, so this is not
    an offline-contract skip -- if it fails, the front door is broken.
    """
    c = dd.load_census()
    missing = sorted(k for k, v in c.items() if v is None)
    assert not missing, (
        f"committed census aggregates missing or unreadable: {missing}. "
        "Regenerate with the scripts/census_*.py generators.")
    assert dd.census_headline(c["design"]) is not None
    assert dd.census_mechanism_rows(c["profile"])
