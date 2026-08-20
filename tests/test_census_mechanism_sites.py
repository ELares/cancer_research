"""Guards for the mechanism-class-by-site analysis.

THE ANALYSIS EXISTS TO OVERTURN HALF OF A MANUSCRIPT CLAIM and confirm the
other half, so the guards protect the thing that makes either verdict
possible: the two classes are measured over the SAME sites with the SAME
denominator. Whatever inflates or deflates a site in the tagged literature
generally is then common-mode, and only a sign DISAGREEMENT carries.

The caveat is load-bearing on exactly one row and is guarded as such: this
physical class holds three mechanisms and omits radiotherapy, which is central
to brain practice, so the brain/CNS row must never be read as a statement about
physically delivered treatment.
"""
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
JSON = REPO / "analysis/census-mechanism-sites.json"
MD = REPO / "analysis/census-mechanism-sites.md"


@pytest.fixture(scope="module")
def d():
    return json.loads(JSON.read_text())


def test_both_classes_are_measured_over_the_same_sites(d):
    """The property the whole comparison rests on."""
    sites = {r["site"] for r in d["rows"]}
    assert sites == set(d["site_totals"])
    for cls in ("physical", "pharmacological"):
        assert set(d["class_by_site"][cls]) <= sites


def test_enrichment_is_recomputable_from_the_stored_counts(d):
    """A derived column that cannot be re-derived is a stored opinion."""
    base_tot = sum(d["site_totals"].values())
    for cls, key in (("physical", "physical_enrichment"),
                     ("pharmacological", "pharmacological_enrichment")):
        ctot = sum(d["class_by_site"][cls].values())
        for r in d["rows"]:
            n = d["class_by_site"][cls].get(r["site"], 0)
            expect = ((n / ctot) / (d["site_totals"][r["site"]] / base_tot)) if n else 0.0
            assert r[key] == pytest.approx(expect, abs=0.005), (
                f"{cls} enrichment for {r['site']} does not recompute")


def test_the_opposed_set_is_derived(d):
    """Sign disagreement is the reading that does not depend on the base rate,
    so it must be computed rather than listed."""
    expect = sorted(r["site"] for r in d["rows"]
                    if (r["physical_enrichment"] - 1)
                    * (r["pharmacological_enrichment"] - 1) < 0)
    assert sorted(d["opposed_sites"]) == expect
    md = MD.read_text()
    if expect:
        assert f"{len(expect)} site(s) where the two classes move in OPPOSITE" in md


def test_the_haematologic_verdict_rests_on_the_contrast_not_the_base_rate(d):
    """If the pharmacological class were ALSO depleted there, the finding would
    be about how those sites are written about rather than about modality."""
    md = MD.read_text()
    if "haematologic half SURVIVES" not in md:
        pytest.skip("the report no longer makes the haematologic claim")
    for site in ("leukaemia", "lymphoma"):
        r = next(x for x in d["rows"] if x["site"] == site)
        assert r["physical_enrichment"] < 1.0, f"{site} is not depleted for physical"
        assert r["pharmacological_enrichment"] > 1.0, (
            f"{site} is not ENRICHED for pharmacological, so the contrast the "
            "report leans on does not hold and the claim reduces to a base-rate "
            "observation")


def test_the_brain_row_carries_the_radiotherapy_caveat(d):
    """The class omits radiotherapy, and brain is where that bites hardest."""
    md = MD.read_text()
    if "neuroectodermal half does NOT survive" in md:
        assert "radiotherapy is outside this physical class" in md
        r = next(x for x in d["rows"] if x["site"] == "brain/CNS")
        assert abs(r["physical_enrichment"] - r["pharmacological_enrichment"]) < 0.15, (
            "the report calls the two classes indistinguishable at brain/CNS "
            f"but they are {r['physical_enrichment']} and "
            f"{r['pharmacological_enrichment']}")


def test_the_class_lists_are_imported_not_restated():
    """A hand-written copy beside the real list is how the #ATLAS-LANDSCAPE
    discrepancy arose."""
    src = (REPO / "scripts/census_mechanism_sites.py").read_text()
    assert "al.PHYSICAL" in src and "al.PHARMACOLOGICAL" in src
    assert "PHYSICAL = {" not in src and "PHARMACOLOGICAL = {" not in src


def test_the_small_denominator_is_reported(d):
    """The physical class is ~50x smaller, so a single row is weak evidence
    even where the ordering is not."""
    assert d["physical_total"] < d["pharmacological_total"] / 5
    md = MD.read_text()
    assert "the ordering is better determined than any single row" in md
