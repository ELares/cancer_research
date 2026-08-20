#!/usr/bin/env python3
"""Can the census's site counts be divided by cancer burden? Probed, not assumed.

THE QUESTION THIS GATES is the most directly useful one a literature census can
ask: which cancers are under-researched relative to how much harm they do. The
census now supports the numerator -- 2,546,944 records assignable to 18 sites
by NLM's own C04 tree, with the assignment measured rather than asserted
(`analysis/atlas-site-coverage.md`) -- and `atlas_site_coverage.py` deliberately
computes no ratio, with a guard forbidding mortality terms as identifiers,
precisely so the question could not be answered before the denominator was.

So the remaining gate is DATA, and this establishes which side of it we are on.
The repo has a precedent for exactly this: `calibration_feasibility.py` asked
whether a calibration target existed for four proposed layers, found none, and
recorded four "proposed and NOT built" rows. Asking first is cheaper than
building and discovering afterwards, and a negative answer is a result: it says
what would unblock the work, which is more actionable than a vague intention.

WHAT A BURDEN DENOMINATOR HAS TO BE, for this join:
  * site-resolved at roughly the granularity of the 18-site list;
  * global, since the census is the world's literature and a national series
    would compare one country's deaths against everyone's papers;
  * recent enough that a 2026 literature count divided by it means something;
  * reachable without registration, because CI reads only committed artifacts
    and a source needing a login cannot be re-derived by a reader.

OFFLINE CONTRACT: this probes the network and writes a committed artifact. CI
reads the artifact and never makes a request.
"""
import argparse
import json
import urllib.error
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
OUT_MD = REPO / "analysis/burden-data-feasibility.md"
OUT_JSON = REPO / "analysis/burden-data-feasibility.json"
UA = {"User-Agent": "cancer-research-feasibility/1.0", "Accept": "application/json"}

# The 18 sites the census assigns to. A source is only useful here if it
# resolves cancer at roughly this granularity.
CENSUS_SITES = [
    "breast", "lung", "colorectal", "skin/melanoma", "liver", "brain/CNS",
    "prostate", "cervix/uterus", "head and neck", "leukaemia", "stomach",
    "lymphoma", "pancreas", "ovary", "kidney", "bladder", "oesophagus",
    "thyroid",
]

CANDIDATES = [
    {
        "source": "IARC Global Cancer Observatory (GLOBOCAN)",
        "why": "the standard global site-resolved incidence and mortality series",
        "urls": [
            "https://gco.iarc.fr/gco-api/v1/cancers",
            "https://gco.iarc.who.int/gco-api/analysis/today/data-population/900/",
            "https://gco.iarc.fr/today/api/data/populations/900/",
            "https://gco.iarc.fr/overtime/api/data/",
        ],
    },
    {
        "source": "WHO Global Health Observatory (OData)",
        "why": "login-free and well documented; already used for other WHO series",
        "urls": [
            "https://ghoapi.azureedge.net/api/Indicator?$select=IndicatorCode,IndicatorName",
        ],
    },
]

# Site-resolved cancer mortality indicators GHO actually carries, found by
# listing all 3,090 indicators rather than by guessing codes.
GHO_SITE_INDICATORS = {
    "SA_0000001438": "breast",
    "SA_0000001439": "colorectal",
    "SA_0000001445": "liver",
    "SA_0000001448": "head and neck",
    "SA_0000001449": "oesophagus",
}


def probe(url: str) -> dict:
    """One request. Records what came back, not whether it 200'd.

    A 200 is not success here: every GCO endpoint tried returns HTTP 200 with
    the single-page-app shell, which an exit-code check or a status check would
    read as a working API. The discriminator is whether the body parses as
    JSON.
    """
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers=UA),
                                    timeout=25) as r:
            head = r.read(400)
            status = r.status
    except urllib.error.HTTPError as e:
        return {"url": url, "status": e.code, "kind": "http-error", "usable": False}
    except Exception as e:
        return {"url": url, "status": None, "kind": type(e).__name__,
                "usable": False}
    is_json = head.lstrip()[:1] in (b"{", b"[")
    return {
        "url": url, "status": status,
        "kind": "json" if is_json else "html-or-other",
        # An HTML body from an API path means the route no longer exists and
        # the framework is serving the app instead.
        "usable": is_json,
    }


def gho_site_coverage() -> dict:
    """What GHO's site-resolved mortality indicators actually contain."""
    out = {}
    for code, site in GHO_SITE_INDICATORS.items():
        try:
            with urllib.request.urlopen(
                    urllib.request.Request(
                        f"https://ghoapi.azureedge.net/api/{code}", headers=UA),
                    timeout=40) as r:
                rows = json.load(r)["value"]
        except Exception as e:
            out[site] = {"code": code, "error": f"{type(e).__name__}"}
            continue
        years = sorted({x.get("TimeDim") for x in rows if x.get("TimeDim")})
        spatial = sorted({x.get("SpatialDimType") for x in rows})
        out[site] = {
            "code": code, "rows": len(rows),
            "years": [years[0], years[-1]] if years else [],
            "spatial_types": spatial,
            "has_global_aggregate": any(x.get("SpatialDim") == "GLOBAL"
                                        for x in rows),
        }
    return out


def scan() -> dict:
    results = []
    for c in CANDIDATES:
        probes = [probe(u) for u in c["urls"]]
        results.append({**c, "probes": probes,
                        "any_json": any(p["usable"] for p in probes)})
    return {"census_sites": CENSUS_SITES, "candidates": results,
            "gho_site_indicators": gho_site_coverage()}


def assemble(d: dict) -> dict:
    gho = d["gho_site_indicators"]
    covered = [s for s in gho if s in d["census_sites"] and "error" not in gho[s]]
    latest = [y for s in covered for y in (gho[s].get("years") or []) if y]
    d = dict(d)
    d["gho_sites_covered"] = sorted(covered)
    d["gho_sites_missing"] = sorted(set(d["census_sites"]) - set(covered))
    d["gho_latest_year"] = max(latest) if latest else None
    d["gho_has_global"] = all(gho[s].get("has_global_aggregate")
                              for s in covered) if covered else False
    d["globocan_reachable"] = next(
        (c["any_json"] for c in d["candidates"] if "GLOBOCAN" in c["source"]),
        False)
    # The verdict is DERIVED from the four requirements in the docstring, so it
    # cannot say "usable" while a requirement is unmet.
    d["requirements"] = {
        "site_resolved_at_census_granularity":
            len(covered) >= len(d["census_sites"]) * 0.8,
        "global_not_national": bool(d["gho_has_global"]),
        "recent_enough": bool(d["gho_latest_year"] and d["gho_latest_year"] >= 2015),
        "reachable_without_registration": bool(d["globocan_reachable"] or covered),
    }
    d["feasible"] = all(d["requirements"].values())
    return d


def render(d: dict) -> str:
    L = ["# Is a burden denominator reachable?\n"]
    L.append(
        f"Generated by `scripts/burden_data_feasibility.py`. The census assigns "
        f"{len(d['census_sites'])} anatomical sites; dividing those counts by "
        f"cancer burden would answer which cancers are under-researched relative "
        f"to the harm they do. This asks whether the denominator exists in a "
        f"form that join can use.\n"
    )
    L.append(f"**Verdict: {'FEASIBLE' if d['feasible'] else 'NOT feasible'} "
             f"from the sources probed.**\n")
    L.append("| requirement | met |")
    L.append("|---|---|")
    for k, v in d["requirements"].items():
        L.append(f"| {k.replace('_', ' ')} | {'yes' if v else '**no**'} |")
    L.append("")
    L.append("## What each source returns\n")
    for c in d["candidates"]:
        L.append(f"**{c['source']}** — {c['why']}\n")
        for p in c["probes"]:
            L.append(f"- `{p['url']}` → HTTP {p['status']}, {p['kind']}")
        L.append("")
    L.append(
        "Every GCO route tried answers **HTTP 200 with the single-page-app "
        "shell**. That matters for how this was tested: a status check or an "
        "exit code would read those as a working API. The discriminator is "
        "whether the body parses as JSON, and none does. The documented API "
        "path has moved, and guessing further is not a method -- this repo has "
        "already hit the same failure with a DepMap catalogue that began "
        "serving HTML where it had served CSV.\n"
    )
    gho = d["gho_site_indicators"]
    L.append("## What WHO GHO covers\n")
    L.append("| census site | indicator | rows | years | global aggregate |")
    L.append("|---|---|--:|---|---|")
    for site in sorted(gho):
        g = gho[site]
        if "error" in g:
            L.append(f"| {site} | {g['code']} | - | *{g['error']}* | - |")
            continue
        yr = "-".join(str(y) for y in g["years"]) if g["years"] else "-"
        L.append(f"| {site} | `{g['code']}` | {g['rows']:,} | {yr} | "
                 f"{'yes' if g['has_global_aggregate'] else '**no**'} |")
    L.append("")
    L.append(
        f"{len(d['gho_sites_covered'])} of {len(d['census_sites'])} census sites "
        f"have a site-resolved mortality indicator, the most recent data is "
        f"{d['gho_latest_year']}, and there is no global aggregate -- the series "
        f"is country-level. Missing entirely: "
        + ", ".join(f"`{s}`" for s in d["gho_sites_missing"]) + ".\n"
    )
    L.append(
        "These are also not a general cancer-mortality panel. They belong to an "
        "alcohol-attributable-burden series, which is why the covered sites are "
        "the alcohol-associated ones and why the panel stops where it does. "
        "Dividing 2026 literature counts by a 2004 country-level "
        "alcohol-attributable subset, and presenting the result as research "
        "effort per death, would be worse than not asking -- the number would "
        "look like an answer.\n"
    )
    L.append("## What would unblock it\n")
    L.append(
        "One thing: a site-resolved global cancer mortality table that can be "
        "committed as a derived artifact. GLOBOCAN is the right source and its "
        "data are downloadable through the web interface; what is missing is a "
        "documented machine route to them. Any of these closes the gap:\n"
    )
    L.append(
        "- a GLOBOCAN export saved by hand and committed, with its retrieval "
        "date and the exact selection recorded, exactly as the CTRPv2 and "
        "cBioPortal legs are handled here;\n"
        "- the current GCO API path, if IARC documents one;\n"
        "- an IHME GBD extract, which is site-resolved and global but needs "
        "registration, so it would have to be fetched once and committed.\n"
    )
    L.append("## What this does NOT say\n")
    L.append(
        "It does not say the ratio would be a good measure once the data "
        "arrived. Literature is global but concentrated in high-income "
        "research systems while mortality is not, so a low ratio could mean "
        "neglect or could mean the disease's burden falls where research does "
        "not happen -- a distinction the ratio alone cannot make. Incidence and "
        "mortality also diverge sharply by site, and a cancer that is already "
        "curable rationally attracts less new work. Those confounds are "
        "arguments for how to report the ratio, not reasons to skip it, but "
        "they should be settled before the number is computed rather than "
        "after.\n"
    )
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--render-only", action="store_true")
    a = ap.parse_args()
    d = assemble(json.loads(OUT_JSON.read_text()) if a.render_only else scan())
    OUT_JSON.write_text(json.dumps(d, indent=1) + "\n")
    OUT_MD.write_text(render(d))
    print(f"wrote {OUT_MD}")
    print(f"  feasible={d['feasible']}  gho sites "
          f"{len(d['gho_sites_covered'])}/{len(d['census_sites'])}  "
          f"latest {d['gho_latest_year']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
