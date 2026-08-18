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
import re
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


# ---------------------------------------------------------------------------
# THE DEEP MAP. The 18-site list above is deliberately shallow, and the page
# used to justify that with a sentence nobody had measured: "a deep subtree
# walk would raise assignability, at the cost of a mapping nobody can audit".
# Both halves were wrong. The benefit is large and is measured below, and the
# repo ALREADY commits its whole cancer definition as a 704-line descriptor
# file -- an enumerated per-site list is no less auditable than that.
#
# THE REAL PROBLEM WITH THE SHALLOW LIST IS NOT ITS SIZE, IT IS THAT ITS DEPTH
# IS NOT UNIFORM. `stomach`, `ovary`, `bladder` and `thyroid` get ONE
# descriptor while `brain/CNS` and `head and neck` get four, and the C04 tree
# holds ~40 leukaemia descriptors against 3 for prostate. So the per-site
# column is understated by very different factors for different sites, and it
# is the column a burden ratio would divide into mortality.
#
# The rules are patterns over the COMMITTED C04 descriptor file, and the
# resolved map is written out beside the report so a placement can be
# disputed the same way the partitions in #724 can.
# ---------------------------------------------------------------------------
C04_TSV = ATLAS / "mesh" / "c04-descriptors.tsv"
OUT_MAP = PROJECT_ROOT / "analysis" / "site-descriptor-map.tsv"

DEEP_RULES = {
    "lung": r"\blung\b|pulmonary blastoma|pancoast|bronchial neoplasms",
    "breast": r"\bbreast\b|phyllodes",
    "colorectal": r"colorectal|colonic neoplasms|rectal neoplasms|anus neoplasms|"
                  r"sigmoid neoplasms|cecal neoplasms|adenomatous polyposis coli",
    "prostate": r"prostat",
    "stomach": r"stomach neoplasms|gastrointestinal stromal|linitis plastica",
    "liver": r"liver neoplasms|hepatocellular|hepatoblastoma|liver cell adenoma",
    "oesophagus": r"esophag",
    "pancreas": r"pancrea|insulinoma|glucagonoma|somatostatinoma|vipoma|"
                r"gastrinoma|zollinger",
    "cervix/uterus": r"uterine|endometrial|trophoblastic|choriocarcinoma|"
                     r"hydatidiform|leiomyoma",
    "ovary": r"ovarian|granulosa cell tumor|thecoma|brenner tumor|"
             r"sertoli-leydig|luteoma|meigs",
    "bladder": r"urinary bladder neoplasms",
    "kidney": r"kidney neoplasms|renal cell|wilms|nephroblastoma|nephroma",
    "brain/CNS": r"brain neoplasms|glio|central nervous system neoplasms|"
                 r"astrocytoma|medulloblastoma|meningioma|ependymoma|"
                 r"oligodendroglioma|neurocytoma|pinealoma|craniopharyngioma|"
                 r"supratentorial|infratentorial|neuroectodermal tumors, primitive|"
                 r"cerebral ventricle neoplasms|skull base neoplasms|"
                 r"spinal cord neoplasms|pituitary neoplasms|hypothalamic neoplasms|"
                 r"meningeal (neoplasms|carcinomatosis)|gliosarcoma|"
                 r"neuroma, acoustic|nerve sheath",
    "leukaemia": r"leukemi|myelodysplastic|myeloproliferative|polycythemia vera|"
                 r"thrombocythemia|primary myelofibrosis|preleukemia|blast crisis|"
                 r"hematologic neoplasms",
    "lymphoma": r"lymphoma|hodgkin|multiple myeloma|plasmacytoma|myelomatosis|"
                r"waldenstrom|sezary|mycosis fungoides|lymphoproliferative|"
                r"immunoproliferative|paraproteinemias|monoclonal gammopathy",
    # `pharyn` alone reaches `Craniopharyngioma`, which is a brain tumour: a
    # substring trap of exactly the kind this repo has been caught by before.
    "head and neck": r"head and neck|mouth neoplasms|laryn|(?<!cranio)pharyn|"
                     r"tongue neoplasms|lip neoplasms|palatal neoplasms|"
                     r"gingival neoplasms|salivary gland neoplasms|parotid|"
                     r"submandibular gland neoplasms|sublingual gland|maxillary|"
                     r"jaw neoplasms|mandibular neoplasms|nose neoplasms|"
                     r"paranasal sinus|otorhinolaryngologic neoplasms|"
                     r"ear neoplasms|tonsillar|esthesioneuroblastoma|"
                     r"ameloblastoma|odontogenic",
    "skin/melanoma": r"melanoma|skin neoplasms|carcinoma, basal cell|"
                     r"carcinoma, merkel cell|keratoacanthoma|bowen|"
                     r"dermatofibrosarcoma|mycosis fungoides|nevus|"
                     r"sweat gland neoplasms|sebaceous gland neoplasms|"
                     r"hair follicle|paget disease, extramammary|"
                     r"keratosis, actinic",
    "thyroid": r"thyroid|parathyroid neoplasms",
}

# The STRICT variant drops two named groups, so the pair brackets the answer
# instead of one list being published as if it were the definition. A single
# number here would be a judgement wearing a measurement's clothes.
STRICT_EXCLUSIONS = {
    "animal or induced model": r"experimental|leukemia l\d|leukemia p\d|"
                              r"radiation-induced|avian|sarcoma 180|"
                              r"sarcoma, yoshida|carcinoma 256|"
                              r"carcinoma, ehrlich|krukenberg",
    "benign or precursor lesion": r"nevus|dysplasia|barrett|cyst$|cyst,|"
                                  r"nodule|intraepithelial neoplasia|"
                                  r"keratosis|leiomyoma|liver cell adenoma|"
                                  r"thyroid nodule",
}


def c04_labels() -> list:
    """The committed cancer definition, read rather than restated."""
    out = []
    for line in C04_TSV.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) >= 2:
            out.append(parts[1].strip())
    if not out:
        raise SystemExit(f"no descriptors read from {C04_TSV}")
    return out


def deep_map() -> dict:
    """site -> {'deep': set, 'strict': set, 'dropped': {reason: [names]}}"""
    labels = c04_labels()
    out = {}
    for site, pat in DEEP_RULES.items():
        deep = {x.lower() for x in labels if re.search(pat, x, re.I)}
        dropped = {}
        strict = set(deep)
        for reason, ex in STRICT_EXCLUSIONS.items():
            hit = {x for x in strict if re.search(ex, x, re.I)}
            if hit:
                dropped[reason] = sorted(hit)
                strict -= hit
        out[site] = {"deep": deep, "strict": strict, "dropped": dropped}
    return out


def write_map(dm: dict) -> None:
    """Commit the resolved map beside the report, so a placement is disputable."""
    L = ["# site\tvariant\tMeSH descriptor",
         "# Resolved from DEEP_RULES in scripts/atlas_site_coverage.py against",
         "# corpus/atlas/mesh/c04-descriptors.tsv. Regenerate with that script.",
         "# `strict` drops animal/induced models and benign or precursor lesions;",
         "# the report publishes both so neither list is taken as the definition."]
    for site in sorted(dm):
        for x in sorted(dm[site]["deep"]):
            v = "deep+strict" if x in dm[site]["strict"] else "deep"
            L.append(f"{site}\t{v}\t{x}")
    OUT_MAP.write_text("\n".join(L) + "\n", encoding="utf-8")


def _lut(mapping: dict) -> dict:
    lut = {}
    for site, descs in mapping.items():
        for d in descs:
            lut.setdefault(d, set()).add(site)
    return lut


GENERIC = "neoplasms"


def scan() -> dict:
    dm = deep_map()
    write_map(dm)
    luts = {
        "shallow": _lut(SITES),
        "deep": _lut({s: v["deep"] for s, v in dm.items()}),
        "strict": _lut({s: v["strict"] for s, v in dm.items()}),
    }
    per_site = {k: Counter() for k in luts}
    assigned = {k: 0 for k in luts}
    multi = {k: 0 for k in luts}
    total = 0
    adjacent = 0
    # the unassigned pile, decomposed rather than narrated
    un_deeper = 0
    un_generic = 0
    un_adjacent = 0
    un_desc = Counter()
    for f in sorted((ATLAS / "records").glob("*.jsonl.gz")):
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                total += 1
                adj = r.get("cancer_basis") != "C04"
                adjacent += adj
                mesh = {m.lower() for m in (r.get("mesh") or [])}
                hits = {}
                for k, lut in luts.items():
                    h = set()
                    for m in mesh & lut.keys():
                        h |= lut[m]
                    hits[k] = h
                    if h:
                        assigned[k] += 1
                        multi[k] += len(h) > 1
                        for s in h:
                            per_site[k][s] += 1
                if not hits["shallow"]:
                    if hits["deep"]:
                        un_deeper += 1
                    else:
                        if GENERIC in mesh:
                            un_generic += 1
                        if adj:
                            un_adjacent += 1
                        for m in mesh:
                            un_desc[m] += 1
    return {
        "census": total,
        "adjacent_basis": adjacent,
        "assigned": assigned["shallow"],
        "multi_site": multi["shallow"],
        "sites": dict(per_site["shallow"].most_common()),
        "n_sites": len(SITES),
        "variants": {
            k: {"assigned": assigned[k], "multi_site": multi[k],
                "sites": dict(per_site[k].most_common()),
                "n_descriptors": sum(len(v) for v in (
                    SITES if k == "shallow" else
                    {s: dm[s]["deep" if k == "deep" else "strict"]
                     for s in dm}).values())}
            for k in luts},
        "map_dropped": {s: dm[s]["dropped"] for s in sorted(dm)
                        if dm[s]["dropped"]},
        "unassigned": {
            "total": total - assigned["shallow"],
            "same_sites_deeper": un_deeper,
            "generic_neoplasms": un_generic,
            "no_c04_descriptor": un_adjacent,
            "top_descriptors": [[k, v] for k, v in un_desc.most_common(20)],
        },
        "excluded_streams": excluded_streams(),
    }


def excluded_streams() -> dict:
    """What this denominator leaves out, which the page never said.

    The `carry no MeSH at all` row it used to print could only ever be zero:
    `atlas_baseline.py` admits a record only when a DescriptorName matches, so
    no record in this stream can lack MeSH. The row measured the admission
    rule. What IS excluded is a second census stream and a sub-population of
    this one, and both change the denominator.
    """
    out = {}
    p = ATLAS / "unindexed-manifest.json"
    if p.exists():
        d = json.loads(p.read_text())
        out["text_matched_no_mesh"] = sum(
            v.get("cancer_text", 0) for v in dict(d.get("files") or {}).values())
    p = ATLAS / "manifest-c04only.json"
    if p.exists():
        d = json.loads(p.read_text())
        out["c04_core"] = sum(
            v.get("cancer", 0) for v in dict(d.get("files") or {}).values())
    return out


def render(d: dict) -> str:
    n, a = d["census"], d["assigned"]
    sites = d["sites"]
    var = d.get("variants") or {}
    ex = d.get("excluded_streams") or {}
    L = ["# Can the census assign articles to cancer sites?", ""]
    L += ["*Generated by `scripts/atlas_site_coverage.py`. This measures the "
          "denominator a burden-weighted analysis would need. It does NOT "
          "compute a burden ratio.*", ""]

    L += ["| | count | share of census |", "|---|--:|--:|"]
    L += [f"| census articles | {n:,} | |",
          f"| **assignable to a site** | **{a:,}** | **{100*a/n:.1f}%** |",
          f"| assigned to more than one | {d['multi_site']:,} | "
          f"{100*d['multi_site']/max(a,1):.1f}% of assigned |", ""]
    L += _denominator_section(d, n, a, ex)
    L += _depth_section(d, n, var)

    L += [f"Across {d['n_sites']} major sites, on the shallow list:", ""]
    L += ["| site | articles | share of census |", "|---|--:|--:|"]
    for s, c in sites.items():
        L.append(f"| {s} | {c:,} | {100*c/n:.2f}% |")
    L += [""]

    lo = min(sites.values()) if sites else 0
    hi = max(sites.values()) if sites else 0
    L += ["## What this says about the burden question", ""]
    L += [f"**{100*a/n:.1f}% of the census is assignable to one of these "
          f"sites** on the shallow list.", ""]
    L += _remainder_section(d, n, a)
    L += [f"The spread across sites is {lo:,} to {hi:,} articles, a factor of "
          f"{hi/max(lo,1):.0f}. That spread is the thing a burden ratio would "
          f"divide into mortality, so it carries directly into every "
          f"literature-per-death figure -- and the depth note above says how "
          f"much of it is the list rather than the literature.", ""]
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
          "* Which list to use. An earlier version of this bullet said a deep "
          "subtree walk would raise assignability \"at the cost of a mapping "
          "nobody can audit\", and NEITHER HALF HAD BEEN MEASURED. The gain is "
          "measured above. The cost is not auditability: this repo already "
          "commits its entire cancer definition as a 704-descriptor file, and "
          "the resolved site map is committed beside this report at "
          "`analysis/site-descriptor-map.tsv`. The real cost is that each "
          "placement is a judgement, so the map has to be reviewed alongside "
          "the counts -- a review burden, not an audit impossibility.",
          ""]
    return "\n".join(L) + "\n"


def _denominator_section(d, n, a, ex) -> list:
    """What the denominator leaves out. The row this replaces could only be 0."""
    L = ["## What is not in this denominator", ""]
    L += ["An earlier version of this table carried a row reading `carry no "
          "MeSH at all | 0 | 0.0%`. THAT ROW COULD NOT HAVE BEEN ANYTHING ELSE: "
          "`atlas_baseline.py` admits a record to this stream only when a MeSH "
          "DescriptorName matches, so no record here can lack MeSH. It measured "
          "the admission rule and read as a property of the literature. What is "
          "actually excluded is this:", ""]
    tm = ex.get("text_matched_no_mesh")
    core = ex.get("c04_core")
    adj = d.get("adjacent_basis") or 0
    if tm:
        both = n + tm
        L += [f"* **{tm:,} MeSH-less cancer articles** sit in a second census "
              f"stream, `corpus/atlas/records_unindexed/`, recovered by text "
              f"match and carrying no descriptors at all. They are excluded by "
              f"choice and cannot be assigned by any descriptor list. Over "
              f"both streams assignability is {a:,} / {both:,} = "
              f"**{100*a/both:.1f}%**, not {100*a/n:.1f}%."]
    if adj:
        L += [f"* **{adj:,} articles ({100*adj/n:.1f}% of this stream)** are "
              f"admitted only by the nine adjacent experimental-context "
              f"descriptors and carry NO C04 descriptor. Every site string is "
              f"a C04 descriptor, so these are unassignable by construction."
              + (f" Over the C04 core alone ({core:,} articles) assignability "
                 f"is **{100*a/core:.1f}%**." if core else "")]
    L += [""]
    return L


def _depth_section(d, n, var) -> list:
    """The shallow list is not uniformly shallow, and that is the real defect."""
    sh, dp, st = var.get("shallow"), var.get("deep"), var.get("strict")
    if not (sh and dp and st):
        return []
    L = ["## The list is shallow, but not uniformly shallow", ""]
    L += [f"The 18 sites are matched by {sh['n_descriptors']} descriptors "
          f"between them -- but not evenly. `stomach`, `ovary`, `bladder` and "
          f"`thyroid` get one each while `brain/CNS` and `head and neck` get "
          f"four, and the C04 tree holds far more for some sites than others. "
          f"So the per-site column is understated by a different factor for "
          f"every site, and it is the column a burden ratio divides into "
          f"mortality.", ""]
    L += [f"Measured against enumerated deeper lists of the SAME 18 sites "
          f"({dp['n_descriptors']} descriptors, and a strict variant of "
          f"{st['n_descriptors']} that drops animal or induced models and "
          f"benign or precursor lesions; both committed at "
          f"`analysis/site-descriptor-map.tsv`):", ""]
    L += ["| site | shallow | deep | strict | deep/shallow | rank shallow -> deep |",
          "|---|--:|--:|--:|--:|--:|"]
    r_sh = {s: i + 1 for i, s in enumerate(sh["sites"])}
    r_dp = {s: i + 1 for i, s in enumerate(dp["sites"])}
    rows = sorted(sh["sites"].items(), key=lambda kv: -kv[1])
    for s, c in rows:
        cd, cs = dp["sites"].get(s, 0), st["sites"].get(s, 0)
        mv = f"{r_sh[s]} -> {r_dp.get(s, 0)}"
        L.append(f"| {s} | {c:,} | {cd:,} | {cs:,} | {cd/max(c,1):.2f}x | "
                 + (f"**{mv}**" if r_sh[s] != r_dp.get(s) else mv) + " |")
    L += [""]
    ratios = {s: dp["sites"].get(s, 0) / max(c, 1) for s, c in sh["sites"].items()}
    worst = sorted(ratios.items(), key=lambda kv: -kv[1])[:4]
    flat = sorted(ratios.items(), key=lambda kv: kv[1])[:4]
    moved = [(s, r_sh[s], r_dp.get(s)) for s in r_sh
             if r_dp.get(s) and r_sh[s] != r_dp[s]]
    L += ["The understatement is not uniform and it is not small: "
          + ", ".join(f"`{s}` {v:.2f}x" for s, v in worst)
          + " against " + ", ".join(f"`{s}` {v:.2f}x" for s, v in flat)
          + ". " + (f"{len(moved)} of {len(r_sh)} sites change rank, "
                    + ", ".join(f"`{s}` {i} -> {j}" for s, i, j in
                                sorted(moved, key=lambda x: -abs(x[1] - x[2]))[:3])
                    + ". " if moved else "")
          + "So the per-site column is comparable within a list and not "
            "across sites, and any burden ratio built on it inherits that.", ""]
    L += [f"Assignability itself goes **{100*sh['assigned']/n:.1f}%** shallow "
          f"-> **{100*st['assigned']/n:.1f}%** strict -> "
          f"**{100*dp['assigned']/n:.1f}%** deep "
          f"(+{dp['assigned']-sh['assigned']:,} articles at the top of that "
          f"range). The shallow figure is the one this page leads with, "
          f"because it is the list whose every member can be read in one "
          f"screen -- but it is a floor, not the census's limit.", ""]
    return L


def _remainder_section(d, n, a) -> list:
    """The remainder was narrated. It is decomposed now."""
    u = d.get("unassigned") or {}
    if not u.get("total"):
        return []
    t = u["total"]
    L = [f"An earlier version of this page said the remainder \"is not a "
         f"failure of the census: much cancer literature is about biology, "
         f"methods or cancer in general rather than a site\". That was "
         f"narrated rather than measured, and it is wrong for a large share of "
         f"it. Of the {t:,} unassigned:", ""]
    L += [f"* **{u['same_sites_deeper']:,} ({100*u['same_sites_deeper']/t:.1f}%)** "
          f"are the SAME 18 sites, named by a deeper descriptor the shallow "
          f"list does not carry. Nothing about them is site-less.",
          f"* {u['generic_neoplasms']:,} ({100*u['generic_neoplasms']/t:.1f}%) "
          f"carry the generic `Neoplasms` descriptor, which is the reading the "
          f"sentence describes.",
          f"* {u['no_c04_descriptor']:,} "
          f"({100*u['no_c04_descriptor']/t:.1f}%) carry no C04 descriptor at "
          f"all and could not be assigned by any list of C04 strings."]
    top = [x for x in (u.get("top_descriptors") or [])
           if x[0] not in {"humans", "neoplasms", "animals", "female", "male",
                           "adult", "middle aged", "aged", "mice"}][:8]
    if top:
        L += [f"* the rest is not featureless. The commonest descriptors on "
              f"articles the list cannot place, excluding check-tags and the "
              f"generic term: "
              + ", ".join(f"`{k}` {v:,}" for k, v in top) + "."]
    L += ["", "So the remainder is substantially a limit of THIS 18-site list "
          "rather than of the census, and the honest version of the original "
          "sentence is much narrower.", ""]
    return L


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
