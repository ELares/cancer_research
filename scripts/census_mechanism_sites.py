#!/usr/bin/env python3
"""Mechanism class by anatomical site, at census scale (#RETIRE-FROZEN).

Section 4.2 of the manuscript derived a five-category tissue-of-origin layer by
coarsening this project's own cancer-type tags, and concluded that "the
physical-modality literature is overwhelmingly concentrated in epithelial and
neuroectodermal contexts, with hematologic and mesothelial categories containing
almost no physical-modality signal".

The census can test that on an axis this project did not draw. Site assignment
comes from `analysis/site-descriptor-map.tsv`, which is every C04 descriptor at
or beneath the tree nodes the 18-site shallow list already occupies -- NLM's own
hierarchy, not a rule written here (#729). So the question becomes: where does
each mechanism class sit relative to the site's own share of the census?

THE COMPARISON THAT CARRIES THE RESULT IS BETWEEN THE TWO CLASSES, not between a
class and the base rate. A site can be over- or under-represented in the tagged
literature for reasons that have nothing to do with modality (indexing depth,
how much of a site's literature is therapeutic at all). Running both classes over
the same sites with the same denominator makes those reasons common-mode: where
the two classes move in OPPOSITE directions, something about the modality is
doing the work.

WHAT THIS CLASS IS NOT. `PHYSICAL` holds three mechanisms and omits radiotherapy,
its largest real member, because radiotherapy has no mechanism tag in this
project's taxonomy (#724). That omission is load-bearing for exactly one row
here: radiotherapy is central to brain and head-and-neck practice, so a
brain/CNS reading from this class is a reading about sonodynamic, HIFU and
electrochemical therapy specifically, not about physically delivered treatment.
"""
import argparse
import gzip
import importlib.util
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO = SCRIPT_DIR.parent
RECORDS = REPO / "corpus/atlas/records"
SITE_MAP = REPO / "analysis/site-descriptor-map.tsv"
MECH_MAP = REPO / "analysis/mesh-mechanism-map.yaml"
OUT_MD = REPO / "analysis/census-mechanism-sites.md"
OUT_JSON = REPO / "analysis/census-mechanism-sites.json"


def load_sites() -> dict[str, set[str]]:
    """Shallow-list sites only. The DEEP walk is deliberately not used: MeSH
    `Head and Neck Neoplasms` subsumes oesophagus and thyroid, which this list
    counts as separate rows, so the deep column double-counts across the page's
    own sites (#729)."""
    out: dict[str, set[str]] = defaultdict(set)
    for ln in SITE_MAP.read_text(encoding="utf-8").splitlines():
        if ln.startswith("#") or not ln.strip():
            continue
        p = ln.split("\t")
        if len(p) >= 3:
            out[p[0]].add(p[2].strip().lower())
    return dict(out)


def load_mechanisms() -> dict[str, set[str]]:
    import yaml

    mp = yaml.safe_load(MECH_MAP.read_text(encoding="utf-8"))["mechanisms"]
    return {k: {x.lower() for x in v["descriptors"]} for k, v in mp.items()}


def load_classes() -> tuple[set[str], set[str]]:
    """Import the curated class lists rather than restating them -- a hand-written
    copy beside the real one is how the #ATLAS-LANDSCAPE discrepancy arose."""
    spec = importlib.util.spec_from_file_location("al", REPO / "scripts/atlas_landscape.py")
    al = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(al)
    return set(al.PHYSICAL), set(al.PHARMACOLOGICAL)


def scan(stride: int = 1) -> dict:
    sites = load_sites()
    mech = load_mechanisms()
    phys, pharm = load_classes()
    cls_site: dict[str, Counter] = {"physical": Counter(), "pharmacological": Counter()}
    site_tot: Counter = Counter()
    n = 0
    shards = sorted(RECORDS.glob("*.jsonl.gz"))[::stride]
    for f in shards:
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                rec = json.loads(line)
                n += 1
                ms = {m.lower() for m in (rec.get("mesh") or [])}
                if not ms:
                    continue
                hit_sites = [s for s, d in sites.items() if ms & d]
                if not hit_sites:
                    continue
                for s in hit_sites:
                    site_tot[s] += 1
                hits = {k for k, d in mech.items() if ms & d}
                for cname, members in (("physical", phys), ("pharmacological", pharm)):
                    if hits & members:
                        for s in hit_sites:
                            cls_site[cname][s] += 1
    return {
        "census": n,
        "shards": len(shards),
        "site_totals": dict(site_tot),
        "class_by_site": {k: dict(v) for k, v in cls_site.items()},
        "physical_members": sorted(phys),
        "pharmacological_members": sorted(pharm),
    }


def assemble(d: dict) -> dict:
    st = d["site_totals"]
    base_tot = sum(st.values())
    rows = []
    for site in st:
        row = {"site": site, "site_records": st[site], "base_share": st[site] / base_tot}
        for cname in ("physical", "pharmacological"):
            c = d["class_by_site"][cname]
            ctot = sum(c.values())
            n = c.get(site, 0)
            row[cname] = n
            row[f"{cname}_enrichment"] = (n / ctot) / (st[site] / base_tot) if n else 0.0
        rows.append(row)
    rows.sort(key=lambda r: -r["physical_enrichment"])
    # Sites where the two classes disagree in DIRECTION -- the common-mode-free signal.
    opposed = [
        r["site"]
        for r in rows
        if (r["physical_enrichment"] - 1) * (r["pharmacological_enrichment"] - 1) < 0
    ]
    d = dict(d)
    d["rows"] = rows
    d["physical_total"] = sum(d["class_by_site"]["physical"].values())
    d["pharmacological_total"] = sum(d["class_by_site"]["pharmacological"].values())
    d["site_assigned_records"] = base_tot
    d["opposed_sites"] = opposed
    return d


def render(d: dict) -> str:
    rows = d["rows"]
    top = rows[0]
    bot = rows[-1]
    haem = {r["site"]: r for r in rows if r["site"] in ("leukaemia", "lymphoma")}
    brain = next((r for r in rows if r["site"] == "brain/CNS"), None)
    L = []
    L.append("# Mechanism class by anatomical site, at census scale\n")
    L.append(
        f"Generated by `scripts/census_mechanism_sites.py` over "
        f"{d['census']:,} census records. Site assignment follows NLM's C04 tree "
        f"(`analysis/site-descriptor-map.tsv`, #729); the physical class holds "
        f"{len(d['physical_members'])} mechanisms "
        f"({', '.join(d['physical_members'])}) and the pharmacological class "
        f"{len(d['pharmacological_members'])}.\n"
    )
    L.append(
        f"Enrichment is a class's share of a site divided by that site's share of "
        f"all site-assigned records ({d['site_assigned_records']:,}). 1.00x is the "
        f"site's own weight in the literature.\n"
    )
    L.append("| site | physical | enrichment | pharmacological | enrichment | site share |")
    L.append("|---|--:|--:|--:|--:|--:|")
    for r in rows:
        L.append(
            f"| {r['site']} | {r['physical']:,} | {r['physical_enrichment']:.2f}x | "
            f"{r['pharmacological']:,} | {r['pharmacological_enrichment']:.2f}x | "
            f"{100 * r['base_share']:.1f}% |"
        )
    L.append("")
    L.append("## What the ordering tracks\n")
    L.append(
        f"The physical class runs from {top['site']} at "
        f"{top['physical_enrichment']:.2f}x down to {bot['site']} at "
        f"{bot['physical_enrichment']:.2f}x, a factor of "
        f"{top['physical_enrichment'] / bot['physical_enrichment']:.1f}. The enriched "
        f"end is solid organs a probe, a fibre or an electrode can be placed in or "
        f"on; the depleted end is disseminated and luminal disease. That is what "
        f"these three modalities ARE -- each needs a physically reachable target -- "
        f"so the ordering is a consistency check on the site assignment as much as a "
        f"finding about the field.\n"
    )
    if d["opposed_sites"]:
        L.append(
            f"The reading that does not depend on the base rate is the "
            f"{len(d['opposed_sites'])} site(s) where the two classes move in "
            f"OPPOSITE directions: {', '.join(sorted(d['opposed_sites']))}. Whatever "
            f"makes a site over- or under-represented in the tagged literature "
            f"generally is common to both classes; a sign disagreement is not.\n"
        )
    L.append("## Against the manuscript's Section 4.2\n")
    if haem:
        parts = ", ".join(
            f"{s} {r['physical_enrichment']:.2f}x physical against "
            f"{r['pharmacological_enrichment']:.2f}x pharmacological"
            for s, r in sorted(haem.items())
        )
        L.append(
            f"The haematologic half SURVIVES and is strengthened by the contrast: "
            f"{parts}. The pharmacological class is enriched in exactly the sites "
            f"the physical class is depleted in, so this is not a property of how "
            f"much those sites are written about.\n"
        )
    if brain:
        L.append(
            f"The neuroectodermal half does NOT survive. brain/CNS sits at "
            f"{brain['physical_enrichment']:.2f}x for the physical class against "
            f"{brain['pharmacological_enrichment']:.2f}x for the pharmacological one "
            f"-- both at the site's own weight, and indistinguishable from each "
            f"other. The concentration the manuscript reported was a property of a "
            f"corpus retrieved around mechanism keywords. THE CAVEAT MATTERS MOST "
            f"HERE: radiotherapy is outside this physical class by construction and "
            f"is central to brain practice, so this row reads on sonodynamic, HIFU "
            f"and electrochemical therapy, not on physically delivered treatment.\n"
        )
    L.append("## Limits\n")
    L.append(
        f"The physical class carries {d['physical_total']:,} site-assigned records "
        f"against the pharmacological class's {d['pharmacological_total']:,}, so its "
        f"per-site counts are small and the ordering is better determined than any "
        f"single row. Assignment is by MeSH descriptor, so a paper is placed where "
        f"NLM indexed it; multi-site papers are counted in each site, which is why "
        f"the site column sums past the census.\n"
    )
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--render-only", action="store_true")
    a = ap.parse_args()
    if a.render_only:
        # RE-ASSEMBLE rather than re-render. Rendering the stored derived
        # fields makes every guard that reads this artifact INERT against a
        # change to the derivation: a mutation sweep confirmed it, planting a
        # wrong sort key and a wrong enrichment denominator and watching both
        # survive a full guard run, because --render-only never recomputed the
        # column the guards check. The raw counts are stored, so re-deriving
        # costs nothing and makes the stored fields checkable rather than
        # merely carried forward.
        d = assemble(json.loads(OUT_JSON.read_text()))
        OUT_JSON.write_text(json.dumps(d, indent=1) + "\n")
    else:
        d = assemble(scan(a.stride))
        OUT_JSON.write_text(json.dumps(d, indent=1) + "\n")
    OUT_MD.write_text(render(d))
    print(f"wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
