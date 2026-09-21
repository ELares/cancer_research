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
each mechanism class sit relative to the site's share of census site assignments?

Both classes use the same census assignment-share baseline. Opposite-direction
flags identify sites whose baseline share lies between the two class shares;
changing that baseline can change the flags even when class counts stay fixed.
The ratio of the two enrichments within a site cancels the shared baseline,
but neither comparison removes class-specific indexing or descriptor-coverage
differences. These are descriptive contrasts, not bias-adjusted effects.

WHAT THIS CLASS IS NOT. `PHYSICAL` holds three mechanisms and omits radiotherapy,
its largest real member, because radiotherapy has no mechanism tag in this
project's taxonomy (#724). That omission is load-bearing for exactly one row
here: radiotherapy is central to brain and head-and-neck practice, so a
brain/CNS reading from this class is a reading about sonodynamic, HIFU and
electrochemical therapy specifically, not about physically delivered treatment.
"""
import argparse
import importlib.util
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from census_input import atlas_root, census_shards, iter_census_shards  # noqa: E402

SCRIPT_DIR = Path(__file__).resolve().parent
REPO = SCRIPT_DIR.parent
RECORDS = atlas_root() / "records"
SITE_MAP = REPO / "analysis/site-descriptor-map.tsv"
MECH_MAP = REPO / "analysis/mesh-mechanism-map.yaml"
OUT_MD = REPO / "analysis/census-mechanism-sites.md"
OUT_JSON = REPO / "analysis/census-mechanism-sites.json"


def load_sites() -> dict[str, set[str]]:
    """Load the committed C04 descendant map; its site groups can overlap.

    For example, head-and-neck descendants include descriptors also assigned
    to oesophagus and thyroid. Counts are unique within each site, not across
    sites (#729).
    """
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
    shards = census_shards(RECORDS, stride)
    for rec in iter_census_shards(shards, RECORDS):
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
    # The two class shares straddle this census baseline. This flag depends on
    # the baseline; it is not a bias-adjusted or baseline-independent signal.
    opposed = [
        r["site"]
        for r in rows
        if (r["physical_enrichment"] - 1) * (r["pharmacological_enrichment"] - 1) < 0
    ]
    d = dict(d)
    d["rows"] = rows
    # Preserve the existing JSON names for consumers. All three totals count
    # site assignments, not unique articles; overlapping sites contribute
    # separately. A row's site_records is unique within that one site only.
    d["physical_total"] = sum(d["class_by_site"]["physical"].values())
    d["pharmacological_total"] = sum(d["class_by_site"]["pharmacological"].values())
    d["site_assigned_records"] = base_tot
    d["opposed_sites"] = opposed
    return d


def render(d: dict) -> str:
    rows = d["rows"]
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
        "Each article contributes once to each matching site and once per "
        "matching class within that site. Sites and classes can overlap, so "
        "summed site assignments are not unique article counts. The legacy JSON "
        "fields `site_assigned_records`, `physical_total` and "
        "`pharmacological_total` store assignment totals.\n"
    )
    if not rows or not d["physical_total"] or not d["pharmacological_total"]:
        L.append(
            f"The input contains {d['site_assigned_records']:,} site assignments, "
            f"with {d['physical_total']:,} physical and "
            f"{d['pharmacological_total']:,} pharmacological class assignments. "
            "At least one class has no site assignments, so the comparison "
            "of class enrichment is unavailable.\n"
        )
        return "\n".join(L)
    top = rows[0]
    bot = rows[-1]
    L.append(
        f"Enrichment divides a site's share of a class's site assignments by "
        f"its share of all census site assignments ({d['site_assigned_records']:,}). "
        "1.00x means equal assignment shares, not equal article prevalence.\n"
    )
    L.append("| site | physical articles | enrichment | pharmacological articles | "
             "enrichment | share of census site assignments |")
    L.append("|---|--:|--:|--:|--:|--:|")
    for r in rows:
        L.append(
            f"| {r['site']} | {r['physical']:,} | {r['physical_enrichment']:.2f}x | "
            f"{r['pharmacological']:,} | {r['pharmacological_enrichment']:.2f}x | "
            f"{100 * r['base_share']:.1f}% |"
        )
    L.append("")
    L.append("## What the ordering tracks\n")
    if bot["physical_enrichment"]:
        L.append(
            f"The physical class runs from {top['site']} at "
            f"{top['physical_enrichment']:.2f}x down to {bot['site']} at "
            f"{bot['physical_enrichment']:.2f}x, a factor of "
            f"{top['physical_enrichment'] / bot['physical_enrichment']:.1f}. "
            "The ordering describes indexed assignment shares, not treatment "
            "performance.\n"
        )
    else:
        L.append(
            "At least one site has no physical-class records, so the ratio of "
            "highest to lowest physical enrichment is undefined.\n"
        )
    if d["opposed_sites"]:
        L.append(
            f"At {len(d['opposed_sites'])} site(s) the two classes move in "
            f"OPPOSITE directions relative to the chosen census assignment-share "
            f"baseline: {', '.join(sorted(d['opposed_sites']))}.\n"
        )
    L.append(
        "Opposite-direction flags depend on the chosen baseline: changing its "
        "site shares can change the flags even when the class counts stay fixed. "
        "The within-site ratio of the two enrichments cancels that shared "
        "baseline, but neither comparison rules out class-specific indexing or "
        "descriptor-coverage differences. These are descriptive contrasts, "
        "not bias-adjusted effects.\n"
    )
    L.append("## Against the manuscript's Section 4.2\n")
    if haem:
        parts = ", ".join(
            f"{s} {r['physical_enrichment']:.2f}x physical against "
            f"{r['pharmacological_enrichment']:.2f}x pharmacological"
            for s, r in sorted(haem.items())
        )
        L.append(
            f"The haematologic rows report {parts}. Each enrichment compares "
            "the site's share of class assignments with its share of all census "
            "site assignments; these values do not establish a biological or "
            "clinical difference between classes.\n"
        )
    if brain:
        L.append(
            f"brain/CNS sits at "
            f"{brain['physical_enrichment']:.2f}x for the physical class against "
            f"{brain['pharmacological_enrichment']:.2f}x for the pharmacological one "
            f"relative to the census assignment-share baseline. The earlier "
            f"retrieved-corpus comparison uses a different population and "
            f"taxonomy, so it cannot isolate the effect of retrieval. "
            f"Radiotherapy is outside this physical class by construction and "
            f"is central to brain practice, so this row reads on sonodynamic, HIFU "
            f"and electrochemical therapy, not on physically delivered treatment.\n"
        )
    L.append("## Limits\n")
    L.append(
        f"The physical class contributes {d['physical_total']:,} site assignments "
        f"against the pharmacological class's {d['pharmacological_total']:,}. These "
        "totals do not measure unique article coverage or establish the precision "
        "of the ordering. Assignment is by MeSH descriptor; an article is counted "
        "once within each matching site, including overlapping sites. Summed "
        "site assignments can exceed the number of unique articles, but need "
        "not exceed the full census.\n"
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
    else:
        d = assemble(scan(a.stride))
    json_text = json.dumps(d, indent=1) + "\n"
    md_text = render(d)
    OUT_JSON.write_text(json_text)
    OUT_MD.write_text(md_text)
    print(f"wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
