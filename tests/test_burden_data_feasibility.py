"""Guards for the burden-denominator feasibility probe.

This artifact exists to STOP work, so the failure that matters is a verdict
reading FEASIBLE while a requirement is unmet -- which would license computing
a research-effort-per-death ratio on a 2004 country-level alcohol-attributable
subset, and that number would look like an answer.

The probe's own discriminator is the interesting part and is guarded: every GCO
route answers HTTP 200 with the single-page-app shell, so a status check or an
exit code reads them as a working API. Only the body tells the truth.

OFFLINE: these read the committed artifact and never make a request.
"""
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
JSON = REPO / "analysis/burden-data-feasibility.json"
MD = REPO / "analysis/burden-data-feasibility.md"
SCRIPT = REPO / "scripts/burden_data_feasibility.py"
SITE_MAP = REPO / "analysis/site-descriptor-map.tsv"


@pytest.fixture(scope="module")
def d():
    return json.loads(JSON.read_text())


def test_the_verdict_is_derived_from_every_requirement(d):
    """A verdict that can be true while a requirement is false is a preference."""
    assert d["feasible"] == all(d["requirements"].values())
    if not d["feasible"]:
        unmet = [k for k, v in d["requirements"].items() if not v]
        assert unmet, "not feasible but every requirement is met"
        assert "NOT feasible" in MD.read_text()
    else:
        assert "**Verdict: FEASIBLE" in MD.read_text()


def test_each_requirement_recomputes_from_the_probe_data(d):
    req = d["requirements"]
    covered = d["gho_sites_covered"]
    assert req["site_resolved_at_census_granularity"] == (
        len(covered) >= len(d["census_sites"]) * 0.8)
    assert req["recent_enough"] == bool(
        d["gho_latest_year"] and d["gho_latest_year"] >= 2015)
    assert req["global_not_national"] == bool(d["gho_has_global"])


def test_a_200_with_an_html_body_is_not_counted_as_usable(d):
    """The whole reason the probe records the BODY.

    Every GCO route returns HTTP 200 and serves the app shell. A checker that
    trusted the status would report the API as available and the analysis as
    unblocked.
    """
    for c in d["candidates"]:
        for p in c["probes"]:
            if p["kind"] != "json":
                assert not p["usable"], (
                    f"{p['url']} returned {p['kind']} and is counted usable")
        assert c["any_json"] == any(p["usable"] for p in c["probes"])
    gco = next(c for c in d["candidates"] if "GLOBOCAN" in c["source"])
    if not gco["any_json"]:
        assert any(p["status"] == 200 and p["kind"] != "json"
                   for p in gco["probes"]), (
            "the GCO probes no longer show the 200-with-HTML case, which is the "
            "reason this probe checks bodies rather than statuses")
        assert "HTTP 200 with the single-page-app" in MD.read_text()


def test_the_site_list_matches_the_one_the_census_actually_assigns(d):
    """A feasibility answer about 18 sites is only about them if it names the
    same 18. Read from the committed site map rather than trusted."""
    sites = set()
    for ln in SITE_MAP.read_text(encoding="utf-8").splitlines():
        if ln.startswith("#") or not ln.strip():
            continue
        parts = ln.split("\t")
        if len(parts) >= 3:
            sites.add(parts[0])
    assert set(d["census_sites"]) == sites, (
        "the feasibility probe's site list differs from the census site map: "
        f"only in probe {sorted(set(d['census_sites']) - sites)}, "
        f"only in map {sorted(sites - set(d['census_sites']))}")


def test_the_covered_and_missing_sites_partition_the_list(d):
    assert not set(d["gho_sites_covered"]) & set(d["gho_sites_missing"])
    assert (set(d["gho_sites_covered"]) | set(d["gho_sites_missing"])
            == set(d["census_sites"]))
    md = MD.read_text()
    for s in d["gho_sites_missing"]:
        assert f"`{s}`" in md, f"{s} has no burden data and the report omits it"


def test_the_report_names_what_would_unblock_it():
    """A negative result that does not say what would change it is a shrug."""
    md = MD.read_text()
    assert "## What would unblock it" in md
    assert "committed as a derived artifact" in md
    # The precedent matters: this repo already handles two registration-gated
    # sources this way, so the route is known to work here.
    assert "CTRPv2" in md


def test_the_report_does_not_endorse_the_ratio_it_gates():
    """Establishing that the data is missing says nothing about whether the
    measure would be sound, and conflating the two would leave a future reader
    thinking only the data stands in the way."""
    md = MD.read_text()
    assert "## What this does NOT say" in md
    # SCOPED TO THE ENDORSING HALF. A first version banned the phrase "a good
    # measure" anywhere, and the report says "It does not say the ratio would be
    # a good measure" -- a NEGATION containing the banned phrase. Substring bans
    # cannot read negation, which is the same trap that makes a GitHub closing
    # keyword fire inside "not fixed: #N". The disclaimer section is excluded
    # and the rest of the document is checked.
    body = md.split("## What this does NOT say")[0]
    for endorsement in ("would show which cancers are neglected",
                        "a good measure", "would answer definitively",
                        "reveals which cancers are neglected"):
        assert endorsement not in body, (
            f"the report endorses the ratio outside its disclaimer: "
            f"{endorsement!r}")
    assert "high-income research systems" in md, (
        "the report does not name the confound that makes a low ratio "
        "ambiguous between neglect and where research happens")


def test_the_probe_stays_offline_for_ci():
    """It reaches the network by design, so nothing may import it at test time
    beyond reading its artifact."""
    src = SCRIPT.read_text()
    assert "OFFLINE CONTRACT" in src
    assert "urllib.request" in src
    assert "--render-only" in src, (
        "the probe cannot re-render without re-probing, so a documentation fix "
        "would require network access")
