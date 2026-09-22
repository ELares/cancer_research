"""Guards for the per-mechanism census profile.

Six manuscript sections read this one artifact. That is the point -- six
separate scans is how one quantity ends up quoted three ways -- but it also
means a defect here reaches six sections at once, so the columns the prose
leans on are pinned to their derivations rather than to their stored values.

The profile declines to rank mechanisms by volume (descriptor breadth varies
enormously) or interpret a co-tagging fraction as a tested-combination rate
(Section 3.13 distinguishes labels, denominators and population). Its report and generator state these
limits. A synthetic expansion of only the target's descriptors also guards
the caveat that trial share, site enrichment and partner ordering describe
the selected articles and can change when descriptor coverage changes.
"""
import gzip
import importlib.util
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
JSON = REPO / "analysis/census-mechanism-profile.json"
MD = REPO / "analysis/census-mechanism-profile.md"
MANUSCRIPT = REPO / "article/drafts/v1.md"


def _load_profile():
    spec = importlib.util.spec_from_file_location(
        "census_profile_input_test", REPO / "scripts/census_mechanism_profile.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("has_sites", [True, False], ids=["overlapping-sites", "no-sites"])
def test_site_assignment_units_from_scan_to_report(tmp_path, monkeypatch, has_sites):
    """Twenty articles can contribute forty assignments but only five trials."""
    profile = _load_profile()
    records = tmp_path / "records"
    records.mkdir()
    monkeypatch.setattr(profile, "RECORDS", records)
    monkeypatch.setattr(profile, "load_sites", lambda: {
        "site-a": {"site a"}, "site-b": {"site b"},
    })
    mech_map = tmp_path / "mechanisms.yaml"
    mech_map.write_text("mechanisms:\n  mechanism-a:\n    descriptors: [Mechanism A]\n")
    monkeypatch.setattr(profile, "MECH_MAP", mech_map)
    with gzip.open(records / "part.jsonl.gz", "wt", encoding="utf-8") as fh:
        for i in range(20):
            fh.write(json.dumps({
                "pmid": str(i),
                "mesh": ["Mechanism A"] + (["Site A", "Site B"] if has_sites else []),
                "pub_types": ["Clinical Trial"] if i < 5 else [],
            }) + "\n")

    result = profile.assemble(profile.scan())
    row, = result["rows"]
    assert result["census"] == row["census"] == 20
    assert row["trials"] == 5
    assert row["trial_share"] == 25.0  # Articles remain the trial-share denominator.
    expected_sites = {"site-a": 20, "site-b": 20} if has_sites else {}
    assert result["site_totals"] == expected_sites
    assert result["by_site"].get("mechanism-a", {}) == expected_sites
    assert row["site_assigned"] == (40 if has_sites else 0)
    assert row["top_sites"] == [
        {"site": site, "n": n, "enrichment": 1.0}
        for site, n in expected_sites.items()
    ]
    report = profile.render(result)
    assert f"{row['site_assigned']} site assignments" in report
    assert "5 carrying a clinical-trial publication type (25.0%)" in report
    assert "a site's share of a mechanism's site assignments" in report
    assert "by its share of all census site assignments" in report
    assert "assignment totals are not unique article counts" in report
    assert "are assignable to a site" not in report
    assert "site-assigned records" not in report


def test_sparse_site_profile_describes_reporting_limits_without_claiming_no_concentration():
    profile = _load_profile()
    result = profile.assemble({
        "census": 100, "site_totals": {"site-a": 19, "site-b": 81},
        "count": {"target": 19}, "trials": {}, "by_year": {},
        "by_site": {"target": {"site-a": 19}}, "partners": {},
    })

    row, = result["rows"]
    assert row["site_assigned"] == 19
    assert row["top_sites"] == []
    # Every target article belongs to a site holding only 19% of the census;
    # the reporting floor must not turn that concentration into an absence.
    report = profile.render(result)
    assert "20-article reporting threshold" in report
    assert "no measurable anatomical concentration" not in report
    assert "other mapped mechanisms" in report
    assert "No mechanism co-occurs with it in the census" not in report


def test_equal_mechanism_counts_have_stable_rows_and_rendering():
    profile = _load_profile()
    raw = {
        "census": 20, "site_totals": {}, "count": {"alpha": 10, "beta": 10},
        "trials": {}, "by_year": {}, "by_site": {}, "partners": {},
    }
    forward = profile.assemble(raw)
    backward = profile.assemble({**raw, "count": {"beta": 10, "alpha": 10}})

    assert [row["mechanism"] for row in forward["rows"]] == ["alpha", "beta"]
    assert backward["rows"] == forward["rows"]
    assert profile.render(backward) == profile.render(forward)


def test_descriptor_expansion_changes_profile_with_fixed_census(tmp_path, monkeypatch):
    """Broader target coverage changes all three summaries within a fixed census."""
    profile = _load_profile()
    records = tmp_path / "records"
    records.mkdir()
    monkeypatch.setattr(profile, "RECORDS", records)
    monkeypatch.setattr(profile, "load_sites", lambda: {
        "site-a": {"site a"}, "site-b": {"site b"},
    })
    mech_map = tmp_path / "mechanisms.json"
    monkeypatch.setattr(profile, "MECH_MAP", mech_map)
    mechanisms = {
        "target": {"descriptors": ["Narrow"]},
        "partner-x": {"descriptors": ["Partner X"]},
        "partner-y": {"descriptors": ["Partner Y"]},
    }
    with gzip.open(records / "part.jsonl.gz", "wt", encoding="utf-8") as fh:
        for i in range(60):
            mesh = ["Narrow" if i < 20 else "Additional",
                    "Site A" if i < 20 else "Site B"]
            if i < 12:
                mesh.append("Partner X")
            if 12 <= i < 16 or 20 <= i < 56:
                mesh.append("Partner Y")
            fh.write(json.dumps({
                "pmid": str(i), "mesh": mesh,
                "pub_types": ["Clinical Trial"] if i < 10 else [],
            }) + "\n")

    results = []
    for descriptors in (["Narrow"], ["Narrow", "Additional"]):
        mechanisms["target"]["descriptors"] = descriptors
        mech_map.write_text(json.dumps({"mechanisms": mechanisms}))
        result = profile.assemble(profile.scan())
        assert result["census"] == 60
        assert result["site_totals"] == {"site-a": 20, "site-b": 40}
        assert result["count"]["partner-x"] == 12
        assert result["count"]["partner-y"] == 40
        assert result["trials"] == {"target": 10, "partner-x": 10}
        results.append(result)

    narrow, broader = results
    assert narrow["count"]["target"] == 20
    assert broader["count"]["target"] == 60
    assert narrow["by_site"]["target"] == {"site-a": 20}
    assert broader["by_site"]["target"] == {"site-a": 20, "site-b": 40}
    assert narrow["partners"]["target"] == {"partner-x": 12, "partner-y": 4}
    assert broader["partners"]["target"] == {"partner-x": 12, "partner-y": 40}
    narrow_row, broader_row = [
        next(row for row in result["rows"] if row["mechanism"] == "target")
        for result in results
    ]
    assert (narrow_row["census"], broader_row["census"]) == (20, 60)
    assert (narrow_row["trials"], broader_row["trials"]) == (10, 10)
    assert (narrow_row["trial_share"], broader_row["trial_share"]) == (50.0, 16.67)
    assert (narrow_row["site_assigned"], broader_row["site_assigned"]) == (20, 60)
    assert narrow_row["top_sites"] == [
        {"site": "site-a", "n": 20, "enrichment": 3.0},
    ]
    assert broader_row["top_sites"] == [
        {"site": "site-a", "n": 20, "enrichment": 1.0},
        {"site": "site-b", "n": 40, "enrichment": 1.0},
    ]
    assert narrow_row["top_partners"] == [
        {"mechanism": "partner-x", "n": 12},
        {"mechanism": "partner-y", "n": 4},
    ]
    assert broader_row["top_partners"] == [
        {"mechanism": "partner-y", "n": 40},
        {"mechanism": "partner-x", "n": 12},
    ]


@pytest.fixture(scope="module")
def d():
    return json.loads(JSON.read_text())


def test_committed_report_uses_current_assignment_units(d):
    profile = _load_profile()
    report = MD.read_text()
    assert report == profile.render(profile.assemble(d))
    for row in d["rows"]:
        assert f"{row['site_assigned']:,} site assignments" in report
    assert "`site_assigned` stores this assignment total" in report
    assert "are assignable to a site" not in report
    assert "site-assigned records" not in report
    assert "curated descriptor map" in report
    assert "selected roots in NLM's C04" in report
    assert "none of the three assigned by this project" not in report


def test_trial_share_recomputes_from_its_own_counts(d):
    for r in d["rows"]:
        assert r["trial_share"] == pytest.approx(
            100 * r["trials"] / r["census"], abs=0.01), r["mechanism"]
        assert r["trials"] <= r["census"], (
            f"{r['mechanism']} has more trials than articles")


def test_site_enrichment_recomputes_and_uses_the_mechanisms_own_denominator(d):
    """Enrichment is against the SITE's weight, normalised by the MECHANISM's
    assigned total -- mixing those two denominators is the easiest way to get a
    plausible-looking wrong number here."""
    base_tot = sum(d["site_totals"].values())
    for r in d["rows"]:
        assigned = r["site_assigned"]
        if not assigned:
            continue
        counted = sum(d["by_site"].get(r["mechanism"], {}).values())
        assert counted == assigned, r["mechanism"]
        for sr in r["top_sites"]:
            n = d["by_site"][r["mechanism"]][sr["site"]]
            assert sr["n"] == n
            assert sr["enrichment"] == pytest.approx(
                (n / assigned) / (d["site_totals"][sr["site"]] / base_tot), abs=0.01)


def test_top_sites_are_ordered_by_enrichment_not_by_raw_count(d):
    """The defect this prevents was live in a draft: by raw volume the census
    puts breast top for sonodynamic therapy, and breast is the top target of
    cancer research generally. Reading the raw column names the wrong site."""
    for r in d["rows"]:
        e = [s["enrichment"] for s in r["top_sites"]]
        assert e == sorted(e, reverse=True), r["mechanism"]


def test_partners_are_symmetric(d):
    """Co-occurrence is a symmetric relation; an asymmetry means the counter
    double-counted or skipped a direction."""
    part = d["partners"]
    for a, others in part.items():
        for b, n in others.items():
            assert part.get(b, {}).get(a) == n, f"{a}<->{b} disagree: {n}"


def test_the_report_refuses_a_volume_ranking_and_a_cooccurrence_rate(d):
    """The report scopes its summaries; the generator also refuses a rate."""
    md = MD.read_text()
    assert "Volume is NOT comparable across mechanisms" in md
    assert "how broad each descriptor is" in md
    assert "articles selected by each mechanism's descriptors" in md
    assert "Partner ordering ranks raw co-occurrence counts" in md
    assert "can change with descriptor coverage and indexing" in md
    assert "Normalization does not establish comparability across mechanism definitions" in md
    assert "survive that objection" not in md
    assert "does not report a co-occurrence RATE" in (
        REPO / "scripts/census_mechanism_profile.py").read_text()


def test_the_unmeasurable_mechanisms_are_absent_not_zero(d):
    """TTFields and bioelectric modulation have no MeSH descriptor. If either
    ever appears here with a count of 0, a reader will take it for a field
    nobody works on -- and TTFields has FDA approval in two indications."""
    for absent in ("ttfields", "bioelectric"):
        assert absent not in d["count"], (
            f"{absent} has no MeSH descriptor; a row for it would report "
            "unmeasurable as zero")
    txt = " ".join(MANUSCRIPT.read_text().split())
    assert "not measurable at census scale" in txt or "unmeasurable" in txt


def test_every_mechanism_the_manuscript_quotes_is_in_the_artifact(d):
    """Catches a figure typed from a scratch run rather than from the artifact."""
    txt = " ".join(MANUSCRIPT.read_text().split())
    rows = {r["mechanism"]: r for r in d["rows"]}
    quoted = [m for m in ("immunotherapy", "nanoparticle", "sonodynamic", "hifu",
                          "microbiome", "metabolic-targeting")
              if f"{rows[m]['census']:,}" in txt]
    assert len(quoted) >= 5, (
        f"only {len(quoted)} of the mechanisms this artifact profiles have "
        "their census volume quoted in the manuscript; the prose may be "
        "carrying figures from somewhere else")
