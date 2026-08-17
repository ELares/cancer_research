#!/usr/bin/env python3
"""Can the census assign articles to cancer sites well enough to weight by burden? (#729)

WHY THIS COMES FIRST
--------------------
#729 wants literature-per-death: attention per cancer site divided by how many
people that site kills. That is the question a reader who cares about helping
people would ask, and nothing in this project has ever asked it -- every
coverage and capture figure here is weighted by publication volume.

But a per-site rate is only as good as the site assignment underneath it, and
nobody has measured that. A reviewer's gate on the issue was explicit: measure
site-assignment completeness first and report it before any ratio. If assignment
is sparse or uneven across sites, a burden ratio is not computable at the
precision the claim would need, and computing one anyway would produce a
confident number whose error structure nobody knows.

So this measures the denominator's denominator, and stops there.

WHAT IT MEASURES
----------------
For each major cancer site, how many census articles carry a MeSH descriptor
that identifies it, and what share of the census is assignable to any site at
all. The interesting quantity is not the total -- it is the SPREAD, because an
attention-per-death ratio compares sites against each other, so uneven
assignability biases the comparison even when total coverage looks fine.

WHAT IT DELIBERATELY DOES NOT DO
---------------------------------
It does not compute a burden ratio, join to GLOBOCAN, or rank sites by neglect.
Doing that here would answer the question before establishing whether it can be
answered. The issue's own framing also warned that mortality and publication
counts have different denominators and different geography -- global mortality
is dominated by regions whose research output is not proportional to their
burden -- so a literature-per-death figure partly measures where science is
funded. That is a real finding too, and a different one.

Usage:
    python scripts/atlas_site_coverage.py
"""

import argparse
import gzip
import json
from collections import Counter
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
ATLAS = PROJECT_ROOT / "corpus" / "atlas"
OUT_MD = PROJECT_ROOT / "analysis" / "atlas-site-coverage.md"
OUT_JSON = PROJECT_ROOT / "analysis" / "atlas-site-coverage.json"

# Major sites by MeSH descriptor. Deliberately the top-level neoplasm
# descriptors rather than a deep subtree walk: the question is whether a site is
# identifiable at all, and a shallow, checkable list makes a low number a
# statement about the census rather than about a mapping nobody can audit.
SITES = {
    "lung": {"lung neoplasms", "carcinoma, non-small-cell lung",
             "small cell lung carcinoma"},
    "breast": {"breast neoplasms", "triple negative breast neoplasms",
               "carcinoma, ductal, breast"},
    "colorectal": {"colorectal neoplasms", "colonic neoplasms",
                   "rectal neoplasms"},
    "prostate": {"prostatic neoplasms", "prostatic neoplasms, castration-resistant"},
    "stomach": {"stomach neoplasms"},
    "liver": {"liver neoplasms", "carcinoma, hepatocellular"},
    "oesophagus": {"esophageal neoplasms", "esophageal squamous cell carcinoma"},
    "pancreas": {"pancreatic neoplasms", "carcinoma, pancreatic ductal"},
    "cervix/uterus": {"uterine cervical neoplasms", "uterine neoplasms",
                      "endometrial neoplasms"},
    "ovary": {"ovarian neoplasms"},
    "bladder": {"urinary bladder neoplasms"},
    "kidney": {"kidney neoplasms", "carcinoma, renal cell"},
    "brain/CNS": {"brain neoplasms", "glioblastoma", "glioma",
                  "central nervous system neoplasms"},
    "leukaemia": {"leukemia", "leukemia, myeloid, acute",
                  "leukemia, lymphocytic, chronic, b-cell"},
    "lymphoma": {"lymphoma", "lymphoma, non-hodgkin", "hodgkin disease"},
    "head and neck": {"head and neck neoplasms",
                      "squamous cell carcinoma of head and neck",
                      "mouth neoplasms", "laryngeal neoplasms"},
    "skin/melanoma": {"melanoma", "skin neoplasms"},
    "thyroid": {"thyroid neoplasms"},
}


def scan() -> dict:
    lut = {}
    for site, descs in SITES.items():
        for d in descs:
            lut.setdefault(d, set()).add(site)

    per_site = Counter()
    assigned = 0
    multi = 0
    total = 0
    no_mesh = 0
    for f in sorted((ATLAS / "records").glob("*.jsonl.gz")):
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                total += 1
                mesh = {m.lower() for m in (r.get("mesh") or [])}
                if not mesh:
                    no_mesh += 1
                    continue
                hits = set()
                for m in mesh & lut.keys():
                    hits |= lut[m]
                if hits:
                    assigned += 1
                    if len(hits) > 1:
                        multi += 1
                    for s in hits:
                        per_site[s] += 1
    return {"census": total, "no_mesh": no_mesh, "assigned": assigned,
            "multi_site": multi, "sites": dict(per_site.most_common()),
            "n_sites": len(SITES)}


def render(d: dict) -> str:
    n, a = d["census"], d["assigned"]
    sites = d["sites"]
    L = ["# Can the census assign articles to cancer sites?", ""]
    L += ["*Generated by `scripts/atlas_site_coverage.py`. This measures the "
          "denominator a burden-weighted analysis would need. It does NOT "
          "compute a burden ratio.*", ""]

    L += ["| | count | share of census |", "|---|--:|--:|"]
    L += [f"| census articles | {n:,} | |",
          f"| carry no MeSH at all | {d['no_mesh']:,} | {100*d['no_mesh']/n:.1f}% |",
          f"| **assignable to a site** | **{a:,}** | **{100*a/n:.1f}%** |",
          f"| assigned to more than one | {d['multi_site']:,} | "
          f"{100*d['multi_site']/max(a,1):.1f}% of assigned |", ""]

    L += [f"Across {d['n_sites']} major sites:", ""]
    L += ["| site | articles | share of census |", "|---|--:|--:|"]
    for s, c in sites.items():
        L.append(f"| {s} | {c:,} | {100*c/n:.2f}% |")
    L += [""]

    lo = min(sites.values()) if sites else 0
    hi = max(sites.values()) if sites else 0
    L += ["## What this says about the burden question", ""]
    L += [f"**{100*a/n:.1f}% of the census is assignable to one of these "
          f"sites.** The remainder is not a failure of the census: much cancer "
          f"literature is about biology, methods or cancer in general rather "
          f"than a site, and a site-weighted analysis simply cannot speak to "
          f"it.", ""]
    L += [f"The spread across sites is {lo:,} to {hi:,} articles, a factor of "
          f"{hi/max(lo,1):.0f}. That spread is the thing a burden ratio would "
          f"divide into mortality, so it carries directly into every "
          f"literature-per-death figure.", ""]
    L += [f"Multi-site assignment is {100*d['multi_site']/max(a,1):.1f}% of "
          f"assigned articles. Those are counted once per site here, so the "
          f"per-site column sums to more than the assigned total -- correct for "
          f"'how much literature touches this site', wrong for a partition. A "
          f"burden ratio has to state which it wants.", ""]

    L += ["## What is still missing before a ratio is defensible", ""]
    L += ["* GLOBOCAN site definitions do not map one-to-one onto MeSH "
          "descriptors, and the mismatch is not uniform: sites like "
          "`cervix/uterus` merge in MeSH where mortality data separates them.",
          "* Mortality and publication counts have different denominators and "
          "different geography. Global mortality is dominated by regions whose "
          "research output is not proportional to their burden, so a "
          "literature-per-death ratio partly measures where science is funded. "
          "That is a real finding, and a different one from neglect.",
          "* The site list here is shallow and checkable by design. A deep "
          "subtree walk would raise assignability, at the cost of a mapping "
          "nobody can audit -- and an unauditable denominator is worse than a "
          "conservative one.",
          ""]
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--render-only", action="store_true")
    args = ap.parse_args()
    if args.render_only:
        d = json.loads(OUT_JSON.read_text())
    else:
        d = scan()
        if d["assigned"] == 0:
            raise SystemExit(
                "no article was assigned to a site, which is not a finding -- "
                "it is what a descriptor-case mismatch looks like.")
        OUT_JSON.write_text(json.dumps(d, indent=1, sort_keys=True) + "\n",
                            encoding="utf-8")
    OUT_MD.write_text(render(d), encoding="utf-8")
    print(f"wrote {OUT_MD}")
    print(f"wrote {OUT_JSON}")
    print(f"  assignable: {d['assigned']:,} of {d['census']:,} "
          f"({100*d['assigned']/d['census']:.1f}%)")
    for s, c in list(d["sites"].items())[:6]:
        print(f"    {s:16s} {c:>8,}")


if __name__ == "__main__":
    main()
