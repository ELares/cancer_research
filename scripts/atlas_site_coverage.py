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
# descriptor while `brain/CNS` and `head and neck` get four, while the tree
# holds far more under some of those nodes than others -- a figure the report
# DERIVES rather than states here, since a hand-written count of names inside
# the function whose thesis is that names are not membership repeats the
# mistake one level down (the count once written here, "~40 against 3", is 43
# against 2 on the tree). So the per-site
# column is understated by very different factors for different sites, and it
# is the column a burden ratio would divide into mortality.
#
# There are no rules and no patterns. The deep list is a walk of NLM's tree
# from the nodes the shallow list already occupies, and the resolved map is
# written out beside the report so a placement can be disputed the same way the
# partitions in #724 can.
# ---------------------------------------------------------------------------
# Experimental models, veterinary disease and named genetic syndromes. These
# sit legitimately under a site in MeSH and are not human disease at that site,
# so they are a scope cost the page names rather than a placement error.
# `syndrome`/`disease` ONCE APPEARED HERE and made the page tell a reader to
# strip `Hodgkin Disease` -- a descriptor the shallow SITES list itself uses to
# define lymphoma -- along with Sezary syndrome and Bowen's disease. 15 of the
# 25 placements it flagged were human cancers. A guard now refuses any shallow
# member here, which is the one-line invariant that would have caught it.
NON_HUMAN = re.compile(
    r"experimental|leukemia l\d|leukemia p\d|radiation-induced|avian|"
    r"bovine|feline|murine|, mouse", re.I)

C04_TSV = ATLAS / "mesh" / "c04-descriptors.tsv"
C04_TREE = ATLAS / "mesh" / "c04-tree-numbers.tsv"
OUT_MAP = PROJECT_ROOT / "analysis" / "site-descriptor-map.tsv"


def c04_labels() -> dict:
    """UI -> label, from the committed cancer definition. Read, not restated."""
    out = {}
    for line in C04_TSV.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) >= 2:
            out[parts[0].strip()] = parts[1].strip().lower()
    if not out:
        raise SystemExit(f"no descriptors read from {C04_TSV}")
    return out


def c04_tree() -> dict:
    """UI -> {tree numbers}. A descriptor may sit at several nodes."""
    out = {}
    for line in C04_TREE.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) >= 2:
            out.setdefault(parts[0].strip(), set()).add(parts[1].strip())
    if not out:
        raise SystemExit(f"no tree numbers read from {C04_TREE}")
    return out


def deep_map() -> dict:
    """site -> {'roots': [...], 'deep': {descriptors}}, DERIVED not written.

    THE FIRST VERSION OF THIS FUNCTION WAS A SUBSTRING MATCHER over descriptor
    labels, and it reproduced -- one function over -- the defect this page
    exists to correct: a hand-written rule per site, with the same class of
    trap. Measured, it put `Ganglion Cysts` and `Paraganglioma` under brain/CNS
    (`ganGLIOn`, `paraganGLIOma`), the benign salivary `Adenolymphoma` under
    lymphoma, a lung disease under cervix/uterus, and it merged 63,620
    plasma-cell records into lymphoma, which moved the headline rank.

    There is no rule here now. A site's deep list is every C04 descriptor
    sitting at or beneath the tree nodes THE SHALLOW LIST ALREADY OCCUPIES, so
    the only judgement is the shallow list, which was already published. NLM's
    tree decides the rest, and a substring accident cannot reach a node it does
    not sit under.
    """
    lab = c04_labels()
    tree = c04_tree()
    by_label = {v: k for k, v in lab.items()}
    out = {}
    for site, descs in SITES.items():
        roots = set()
        for d in descs:
            ui = by_label.get(d)
            if ui is None:
                raise SystemExit(
                    f"{site}: {d!r} is not a C04 descriptor, so the shallow "
                    "list is not a subset of the census definition")
            roots |= tree.get(ui, set())
        # a root beneath another root adds nothing and would double the prose
        roots = {r for r in roots
                 if not any(r != o and r.startswith(o + ".") for o in roots)}
        deep = {lab[u] for u, ts in tree.items()
                if any(t == r or t.startswith(r + ".") for t in ts for r in roots)}
        if not descs <= deep:
            raise SystemExit(
                f"{site}: the subtree walk does not contain its own roots "
                f"({sorted(descs - deep)}), which cannot happen unless the "
                "tree file and the descriptor file describe different builds")
        out[site] = {"roots": sorted(roots), "deep": deep}
    return out


def write_map(dm: dict) -> None:
    """Commit the resolved map, so a placement is disputable."""
    L = ["# site\ttree root\tMeSH descriptor",
         "# DERIVED: every C04 descriptor at or beneath the tree nodes the",
         "# shallow SITES list already occupies. No rule beyond that list --",
         "# NLM's tree decides membership, from corpus/atlas/mesh/",
         "# c04-tree-numbers.tsv. Regenerate with scripts/atlas_site_coverage.py."]
    lab = c04_labels()
    tree = c04_tree()
    by_label = {v: k for k, v in lab.items()}
    for site in sorted(dm):
        for x in sorted(dm[site]["deep"]):
            ts = sorted(tree.get(by_label.get(x, ""), {""}))
            root = next((r for r in dm[site]["roots"]
                         for t in ts if t == r or t.startswith(r + ".")), "?")
            L.append(f"{site}\t{root}\t{x}")
    OUT_MAP.write_text("\n".join(L) + "\n", encoding="utf-8")


def _pairs(v):
    """Ordered pairs, tolerating a dict from an older artifact.

    `json.dumps(sort_keys=True)` reorders a dict, so a ranking stored as one
    comes back alphabetical and `--render-only` renders a DIFFERENT report --
    the defect #730 already fixed in this repo and this file did not inherit.
    """
    if isinstance(v, dict):
        return sorted(v.items(), key=lambda kv: -kv[1])
    return [tuple(x) for x in v]


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
    res_desc = Counter()   # the RESIDUE alone -- see the note in _remainder_section
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
                        # THE PROSE'S "the descriptors above" pointed at
                        # `un_desc`, which accumulates over every unassigned
                        # record -- including the generic-`Neoplasms` and
                        # no-C04 ones the same sentence has just excluded. The
                        # residue gets its own counter.
                        if GENERIC not in mesh and not adj:
                            for m in mesh:
                                res_desc[m] += 1
    return {
        "census": total,
        "adjacent_basis": adjacent,
        "assigned": assigned["shallow"],
        "multi_site": multi["shallow"],
        "sites": [[k, v] for k, v in per_site["shallow"].most_common()],
        "n_sites": len(SITES),
        "variants": {
            k: {"assigned": assigned[k], "multi_site": multi[k],
                "sites": [[x, y] for x, y in per_site[k].most_common()],
                "n_descriptors": sum(
                    len(SITES[s] if k == "shallow" else dm[s]["deep"])
                    for s in SITES)}
            for k in luts},
        "site_tree_roots": {s: dm[s]["roots"] for s in sorted(dm)},
        # NLM's subtree carries more than benign and precursor entities, which
        # was the page's only stated scope cost. Derived and named.
        "non_human_disease_placements": {
            s: sorted(x for x in dm[s]["deep"] if NON_HUMAN.search(x))
            for s in sorted(dm)
            if any(NON_HUMAN.search(x) for x in dm[s]["deep"])},
        "n_placements": sum(len(dm[s]["deep"]) for s in dm),
        "n_distinct_descriptors": len(set().union(*(dm[s]["deep"] for s in dm))),
        "descriptors_in_more_than_one_site": sorted(
            x for x in set().union(*(dm[s]["deep"] for s in dm))
            if sum(1 for s in dm if x in dm[s]["deep"]) > 1),
        # NLM's tree does not agree with this page's 18 sites about where the
        # boundaries are. `Head and Neck Neoplasms` SUBSUMES oesophagus and
        # thyroid, which are listed here separately -- so the deep column
        # double-counts across the page's own list, and its rank order is
        # partly a statement about MeSH rather than about the literature.
        "deep_site_overlaps": {
            a: sorted(b for b in SITES
                      if b != a and SITES[b] <= dm[a]["deep"])
            for a in sorted(SITES)
            if any(b != a and SITES[b] <= dm[a]["deep"] for b in SITES)},
        "unassigned": {
            "total": total - assigned["shallow"],
            "same_sites_deeper": un_deeper,
            "generic_neoplasms": un_generic,
            "no_c04_descriptor": un_adjacent,
            "top_descriptors": [[k, v] for k, v in un_desc.most_common(20)],
            "residue_top_descriptors": [[k, v]
                                        for k, v in res_desc.most_common(24)],
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


def _roundtrip(d: dict) -> dict:
    """Render from what the artifact WILL contain, not from the live dict.

    `OUT_JSON` is written with `sort_keys=True`, so rendering the in-memory
    dict produces a document a `--render-only` run cannot reproduce. Both paths
    now render the same value.

    Checked rather than assumed: the committed report already matches the
    round-tripped render, so no published ordering changes and no table here
    turned out to depend on the declared SITES sequence. Where an order DOES
    carry meaning it must be re-established inside the renderer -- sorting the
    input replaces a rank with an alphabet, which flipped a published verdict
    elsewhere in this repo.
    """
    return json.loads(json.dumps(d, sort_keys=True))


def render(d: dict) -> str:
    n, a = d["census"], d["assigned"]
    sites = _pairs(d["sites"])
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
    for s, c in sites:
        L.append(f"| {s} | {c:,} | {100*c/n:.2f}% |")
    L += [""]

    lo = min((c for _s, c in sites), default=0)
    hi = max((c for _s, c in sites), default=0)
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
          "`analysis/site-descriptor-map.tsv`. Nor is each placement a "
          "judgement any more: the deep list is every descriptor at or "
          "beneath the tree nodes the shallow list ALREADY occupies, so NLM "
          "decides membership and the only thing to dispute is the shallow "
          "list. What the deeper list does cost is scope -- a site's subtree "
          "carries benign and precursor entities alongside the malignancies, "
          "and a burden ratio has to say whether it wants them.",
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
    sh, dp = var.get("shallow"), var.get("deep")
    if not (sh and dp):
        return []
    L = ["## The list is shallow, but not uniformly shallow", ""]
    L += [f"The 18 sites are matched by {sh['n_descriptors']} descriptors "
          f"between them -- but not evenly. `stomach`, `ovary`, `bladder` and "
          f"`thyroid` get one each while `brain/CNS` and `head and neck` get "
          f"four, and the tree holds far more under some of those nodes than "
          f"others. "
          f"So the per-site column is understated by a different factor for "
          f"every site, and it is the column a burden ratio divides into "
          f"mortality.", ""]
    L += [f"Measured against the SAME 18 sites walked down NLM's own tree -- "
          f"every C04 descriptor at or beneath the nodes the shallow list "
          f"already occupies, {dp['n_descriptors']} placements over "
          f"{d.get('n_distinct_descriptors', dp['n_descriptors'])} distinct "
          f"descriptors ({len(d.get('descriptors_in_more_than_one_site') or [])} "
          f"sit in more than one of these sites), committed at "
          f"`analysis/site-descriptor-map.tsv`. There is no rule to dispute "
          f"beyond the shallow list itself: an earlier version of this "
          f"section matched descriptor NAMES and put `Ganglion Cysts` and "
          f"`Paraganglioma` under brain/CNS, a benign salivary tumour under "
          f"lymphoma, and merged plasma-cell myeloma into lymphoma, which "
          f"moved the headline rank.", ""]
    L += ["| site | shallow | deep | deep/shallow | rank shallow -> deep |",
          "|---|--:|--:|--:|--:|"]
    sh_pairs, dp_pairs = _pairs(sh["sites"]), _pairs(dp["sites"])
    dp_by = dict(dp_pairs)
    rows = sh_pairs
    # A DEEP RANK IS MEANINGLESS FOR A SITE WHOSE SUBTREE CONTAINS ANOTHER OF
    # THESE SITES. `head and neck` reaches rank 1 only because MeSH puts
    # oesophagus and thyroid under it -- 53% of its deep gain is descriptors
    # this same table counts as separate rows -- so the cell is suppressed
    # rather than printed as a finding.
    # A SITE WITH NO RANK MUST NOT OCCUPY ONE. Suppressing the cell while
    # leaving the site in the ranking shifted every other site by one, so the
    # bolded "changes" the page presented were partly artifacts of the rank it
    # had just declared meaningless -- and one real move was hidden.
    merged = set(d.get("deep_site_overlaps") or {})
    dp_pairs = [(x, c) for x, c in dp_pairs if x not in merged]
    sh_rank_pairs = [(x, c) for x, c in sh_pairs if x not in merged]
    dp_by_all = dict(_pairs(dp["sites"]))
    r_sh = {x: i + 1 for i, (x, _c) in enumerate(sh_rank_pairs)}
    r_dp = {x: i + 1 for i, (x, _c) in enumerate(dp_pairs)}
    for s, c in rows:
        cd = dp_by_all.get(s, 0)
        if s in merged:
            cell = "n/a (subsumes " + ", ".join(
                f"`{x}`" for x in d["deep_site_overlaps"][s]) + ")"
        else:
            mv = f"{r_sh[s]} -> {r_dp.get(s, 0)}"
            cell = f"**{mv}**" if r_sh[s] != r_dp.get(s) else mv
        L.append(f"| {s} | {c:,} | {cd:,} | {cd/max(c,1):.2f}x | {cell} |")
    L += [""]
    ratios = {s: dp_by_all.get(s, 0) / max(c, 1) for s, c in sh_pairs}
    worst = sorted(ratios.items(), key=lambda kv: -kv[1])[:4]
    flat = sorted(ratios.items(), key=lambda kv: kv[1])[:4]
    moved = [(s, r_sh[s], r_dp.get(s)) for s in r_sh
             if s not in merged and r_dp.get(s) and r_sh[s] != r_dp[s]]
    ov = d.get("deep_site_overlaps") or {}
    if ov:
        L += ["**Read the deep column with its overlaps.** NLM's tree does not "
              "draw this page's 18 boundaries: "
              + "; ".join(f"`{a}` subsumes " + ", ".join(f"`{x}`" for x in b)
                          for a, b in sorted(ov.items()))
              + ". Those sites are listed separately here, so the deep column "
                "double-counts across the page's own list and its rank order "
                "is partly a statement about MeSH rather than about the "
                "literature. That is a reason to read the ratio column rather "
                "than the deep ranks, and a reason a burden analysis has to "
                "pick its boundaries before it picks its depth.", ""]
    def _lab(x):
        return f"`{x}`" + ("*" if x in merged else "")
    L += ["The gap between the two lists is not uniform and it is not small "
          "(a `*` marks a site whose subtree contains another of these sites, "
          "so its figure counts descriptors this table lists separately): "
          + ", ".join(f"{_lab(s)} {v:.2f}x" for s, v in worst)
          + " against " + ", ".join(f"{_lab(s)} {v:.2f}x" for s, v in flat)
          + ". " + (f"{len(moved)} of {len(r_sh)} rankable sites change rank, "
                    + ", ".join(f"`{s}` {i} -> {j}" for s, i, j in
                                sorted(moved, key=lambda x: -abs(x[1] - x[2]))[:3])
                    + ". " if moved else "")
          + "So the per-site column is comparable within a list and not "
            "across sites, and any burden ratio built on it inherits that.", ""]
    L += [f"Assignability itself goes **{100*sh['assigned']/n:.1f}%** shallow "
          f"-> **{100*dp['assigned']/n:.1f}%** on the subtree walk "
          f"(+{dp['assigned']-sh['assigned']:,} articles). The shallow figure "
          f"is the one this page leads with, because it is the shorter and "
          f"more conservative list -- NOT because it is more auditable, which "
          f"the bullet below retracts. It is a floor, not the census's limit, "
          f"and a deeper list can also over-reach: membership here is NLM's "
          f"tree, so an accident of naming cannot cause that, but a site's "
          f"subtree still carries entities its shallow row does not.", ""]
    nh = d.get("non_human_disease_placements") or {}
    if nh:
        tot = sum(len(v) for v in nh.values())
        L += [f"**And not only benign or precursor ones, which was this page's "
              f"whole stated scope cost.** {tot} placements are experimental "
              f"models, veterinary disease or named genetic syndromes: "
              + "; ".join(
                  f"`{k}` gets " + ", ".join(f"`{x}`" for x in v[:4])
                  + (f" and {len(v)-4} more" if len(v) > 4 else "")
                  for k, v in sorted(nh.items(), key=lambda kv: -len(kv[1]))[:3])
              + ". They sit legitimately under those nodes in MeSH and are not "
                "human disease at that site, so a burden analysis wanting the "
                "deeper list has to strip them.", ""]
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
          f"are the SAME 18 sites, named by a descriptor beneath the shallow "
          f"list's own tree nodes. The great majority of them name a site.",
          f"* {u['generic_neoplasms']:,} ({100*u['generic_neoplasms']/t:.1f}%) "
          f"carry the generic `Neoplasms` descriptor, which is the reading the "
          f"sentence describes.",
          f"* {u['no_c04_descriptor']:,} "
          f"({100*u['no_c04_descriptor']/t:.1f}%) carry no C04 descriptor at "
          f"all and could not be assigned by any list of C04 strings."]
    CHECKTAG = {"humans", "neoplasms", "animals", "female", "male", "adult",
                "middle aged", "aged", "mice", "adolescent", "child", "infant",
                "young adult", "aged, 80 and over", "child, preschool",
                "retrospective studies", "treatment outcome", "rats",
                "prognosis", "cell line, tumor", "antineoplastic agents",
                "cell proliferation"}
    top = [x for x in (u.get("residue_top_descriptors") or [])
           if x[0] not in CHECKTAG][:8]
    if top:
        L += [f"* the rest is not featureless, and it is largely CANCER AT A "
              f"SITE THIS LIST DOES NOT COVER -- the same correction bullet "
              f"one makes, one level out. Its commonest descriptors, "
              f"excluding check-tags and study-design terms: "
              + ", ".join(f"`{k}` {v:,}" for k, v in top) + ". An earlier "
              f"version pointed at a descriptor list accumulated over ALL "
              f"unassigned records, which included the generic-`Neoplasms` "
              f"and no-C04 buckets the same sentence had just excluded."]
    rest = t - u["same_sites_deeper"] - u["generic_neoplasms"] - u["no_c04_descriptor"]
    L += [f"* the remaining {rest:,} ({100*rest/t:.1f}%) is none of those "
          f"three, and is the largest single bucket."]
    L += ["", f"So the honest version of the original sentence is much "
          f"narrower: {100*u['same_sites_deeper']/t:.1f}% of the remainder is "
          f"a limit of THIS 18-site list rather than of the census, a further "
          f"{100*u['generic_neoplasms']/t:.1f}% is the reading the original "
          f"sentence described, and {100*rest/t:.1f}% is neither and is not "
          f"characterised here beyond the descriptors above.", ""]
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
    OUT_MD.write_text(render(_roundtrip(d)), encoding="utf-8")
    print(f"wrote {OUT_MD}")
    print(f"wrote {OUT_JSON}")
    print(f"  assignable: {d['assigned']:,} of {d['census']:,} "
          f"({100*d['assigned']/d['census']:.1f}%)")
    for s, c in _pairs(d["sites"])[:6]:
        print(f"    {s:16s} {c:>8,}")


if __name__ == "__main__":
    main()
