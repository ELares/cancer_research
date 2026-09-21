"""Guards for the mechanism-class-by-site analysis.

The two classes use the same sites and census assignment denominator.
Opposite enrichment directions depend on that baseline; sharing it does not
adjust for indexing or selection bias. The ratio between class enrichments
cancels the baseline mathematically, while the opposite-direction flag does not.

The caveat is load-bearing on exactly one row and is guarded as such: this
physical class holds three mechanisms and omits radiotherapy, which is central
to brain practice, so the brain/CNS row must never be read as a statement about
physically delivered treatment.
"""
import gzip
import importlib.util
import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
JSON = REPO / "analysis/census-mechanism-sites.json"
MD = REPO / "analysis/census-mechanism-sites.md"


def _load_sites():
    spec = importlib.util.spec_from_file_location(
        "census_sites_units_test", REPO / "scripts/census_mechanism_sites.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_overlapping_sites_and_classes_keep_assignment_denominators(tmp_path, monkeypatch):
    """Multiple mechanisms in one class do not multiply that class's count."""
    sites = _load_sites()
    records = tmp_path / "records"
    records.mkdir()
    monkeypatch.setattr(sites, "RECORDS", records)
    monkeypatch.setattr(sites, "load_sites", lambda: {
        "site-a": {"site a"}, "site-b": {"site b"},
    })
    monkeypatch.setattr(sites, "load_mechanisms", lambda: {
        "physical-a": {"physical a"}, "physical-b": {"physical b"},
        "drug": {"drug"},
    })
    monkeypatch.setattr(sites, "load_classes", lambda: (
        {"physical-a", "physical-b"}, {"drug"},
    ))
    with gzip.open(records / "part.jsonl.gz", "wt", encoding="utf-8") as fh:
        for i in range(2):
            fh.write(json.dumps({
                "pmid": str(i),
                "mesh": ["Site A", "Site B", "Physical A", "Physical B"]
                + (["Drug"] if i == 0 else []),
            }) + "\n")

    result = sites.assemble(sites.scan())
    assert result["census"] == 2
    assert result["site_totals"] == {"site-a": 2, "site-b": 2}
    assert result["class_by_site"] == {
        "physical": {"site-a": 2, "site-b": 2},
        "pharmacological": {"site-a": 1, "site-b": 1},
    }
    assert result["site_assigned_records"] == result["physical_total"] == 4
    assert result["pharmacological_total"] == 2
    assert sum(row["base_share"] for row in result["rows"]) == 1.0
    for row in result["rows"]:
        assert row["site_records"] == row["physical"] == 2
        assert row["pharmacological"] == 1
        assert row["base_share"] == 0.5
        assert row["physical_enrichment"] == row["pharmacological_enrichment"] == 1.0
    report = sites.render(result)
    assert "site assignments (4)" in report
    assert ("physical class contributes 4 site assignments against the "
            "pharmacological class's 2") in report
    assert "a site's share of a class's site assignments" in report
    assert "share of census site assignments |" in report
    assert "summed site assignments are not unique article counts" in report
    assert "site-assigned records" not in report
    assert "sums past the census" not in report


@pytest.mark.parametrize("present_class", [None, "physical", "pharmacological"])
def test_missing_class_assignments_use_assignment_units(present_class):
    sites = _load_sites()
    counts = {"physical": {}, "pharmacological": {}}
    if present_class:
        counts[present_class] = {"site-a": 1, "site-b": 1}
    result = sites.assemble({
        "census": 1,
        "site_totals": {"site-a": 1, "site-b": 1} if present_class else {},
        "class_by_site": counts,
        "physical_members": ["physical-a"],
        "pharmacological_members": ["drug"],
    })
    assert result["site_assigned_records"] == (2 if present_class else 0)
    assert result["physical_total"] == (2 if present_class == "physical" else 0)
    assert result["pharmacological_total"] == (2 if present_class == "pharmacological" else 0)
    report = sites.render(result)
    assert f"{result['site_assigned_records']} site assignments" in report
    assert "comparison of class enrichment is unavailable" in report
    assert "site-assigned records" not in report


@pytest.fixture(scope="module")
def d():
    return json.loads(JSON.read_text())


def test_committed_report_uses_current_assignment_units(d):
    sites = _load_sites()
    report = MD.read_text()
    assert report == sites.render(sites.assemble(d))
    assert f"site assignments ({d['site_assigned_records']:,})" in report
    assert f"{d['physical_total']:,} site assignments" in report
    assert "site-assigned records" not in report
    assert "sums past the census" not in report


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
    """Opposite directions are derived relative to the current census baseline."""
    expect = sorted(r["site"] for r in d["rows"]
                    if (r["physical_enrichment"] - 1)
                    * (r["pharmacological_enrichment"] - 1) < 0)
    assert sorted(d["opposed_sites"]) == expect
    md = MD.read_text()
    if expect:
        assert f"At {len(expect)} site(s) the two classes move in OPPOSITE" in md
    assert "Opposite-direction flags depend on the chosen baseline" in md
    assert "not bias-adjusted effects" in md
    assert "does not depend on the base rate" not in md


def test_opposed_flags_depend_on_baseline_but_class_enrichment_ratios_do_not():
    """Moving the baseline alone can remove every opposite-direction flag."""
    sites = _load_sites()
    class_counts = {
        "physical": {"a": 60, "b": 40},
        "pharmacological": {"a": 30, "b": 70},
    }
    results = [sites.assemble({"site_totals": baseline, "class_by_site": class_counts})
               for baseline in ({"a": 450, "b": 550}, {"a": 650, "b": 350})]
    assert set(results[0]["opposed_sites"]) == {"a", "b"}
    assert results[1]["opposed_sites"] == []
    assert results[0]["class_by_site"] == results[1]["class_by_site"] == class_counts
    for result in results:
        ratios = {row["site"]: row["physical_enrichment"] / row["pharmacological_enrichment"]
                  for row in result["rows"]}
        assert ratios == pytest.approx({"a": 2.0, "b": 4 / 7})


def test_haematologic_description_reports_current_enrichment_values(d):
    """The descriptive contrast must quote both classes at the current baseline."""
    md = MD.read_text()
    for site in ("leukaemia", "lymphoma"):
        r = next(x for x in d["rows"] if x["site"] == site)
        assert r["physical_enrichment"] < 1.0, f"{site} is not depleted for physical"
        assert r["pharmacological_enrichment"] > 1.0, f"{site} is not enriched for pharmacological"
        assert (f"{site} {r['physical_enrichment']:.2f}x physical against "
                f"{r['pharmacological_enrichment']:.2f}x pharmacological") in md


def test_the_brain_row_carries_the_radiotherapy_caveat(d):
    """The class omits radiotherapy, and brain is where that bites hardest."""
    md = MD.read_text()
    assert "Radiotherapy is outside this physical class" in md
    r = next(x for x in d["rows"] if x["site"] == "brain/CNS")
    assert (f"brain/CNS sits at {r['physical_enrichment']:.2f}x for the physical class "
            f"against {r['pharmacological_enrichment']:.2f}x for the pharmacological one") in md
    assert "indistinguishable" not in md


def test_the_class_lists_are_imported_not_restated():
    """A hand-written copy beside the real list is how the #ATLAS-LANDSCAPE
    discrepancy arose."""
    src = (REPO / "scripts/census_mechanism_sites.py").read_text()
    assert "al.PHYSICAL" in src and "al.PHARMACOLOGICAL" in src
    assert "PHYSICAL = {" not in src and "PHARMACOLOGICAL = {" not in src


def test_assignment_totals_do_not_establish_coverage_or_precision(d):
    """Different assignment totals alone establish neither coverage nor precision."""
    assert d["physical_total"] < d["pharmacological_total"] / 5
    md = MD.read_text()
    assert (f"physical class contributes {d['physical_total']:,} site assignments "
            f"against the pharmacological class's {d['pharmacological_total']:,}") in md
    assert "totals do not measure unique article coverage" in md
    assert "or establish the precision of the ordering" in md
    assert "the ordering is better determined than any single row" not in md
    # A correct report cannot rescue an old article-count claim in the book.
    manuscript = (REPO / "article/drafts/v1.md").read_text()
    section = re.search(r"^### 4\.2 [^\n]+\n(.*?)(?=^### |\Z)",
                        manuscript, re.M | re.S).group(1)
    assert f"{d['physical_total']:,} site assignments" in section
    assert f"pharmacological class's {d['pharmacological_total']:,}" in section
    assert "not counts of unique articles" in section
    assert "site-assigned records" not in section
    caption = re.search(r"\[FIGURE 11: ([^\n]+)\]", section).group(1)
    assert "share of class site assignments divided by its share of all census site assignments" in caption
