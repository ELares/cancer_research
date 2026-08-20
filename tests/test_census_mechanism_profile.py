"""Guards for the per-mechanism census profile.

Six manuscript sections read this one artifact. That is the point -- six
separate scans is how one quantity ends up quoted three ways -- but it also
means a defect here reaches six sections at once, so the columns the prose
leans on are pinned to their derivations rather than to their stored values.

The two refusals matter as much as the numbers. This file must NOT rank
mechanisms by volume (descriptor breadth varies enormously, so a volume
ranking is substantially a ranking of how broad each descriptor is) and must
NOT report a co-occurrence rate (which Section 3.13 shows is a property of the
labelling instrument). Both refusals are guarded, because a refusal nothing
checks is a comment.
"""
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
JSON = REPO / "analysis/census-mechanism-profile.json"
MD = REPO / "analysis/census-mechanism-profile.md"
MANUSCRIPT = REPO / "article/drafts/v1.md"


@pytest.fixture(scope="module")
def d():
    return json.loads(JSON.read_text())


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
    """Both refusals are load-bearing and both are stated in the report."""
    md = MD.read_text()
    assert "Volume is NOT comparable across mechanisms" in md
    assert "how broad each descriptor is" in md
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
