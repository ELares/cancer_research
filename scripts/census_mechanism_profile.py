#!/usr/bin/env python3
"""One census profile per mechanism, joining what six manuscript sections need.

Sections 3.8 through 3.12 each narrate one mechanism, and each was written over
a retrieved corpus: a volume, a spread across cancer types, a set of convergence
partners. Rebuilding them at census scale from six separate ad-hoc scans is how
the same quantity ends up quoted three different ways in three sections, so this
computes all of it in ONE pass and writes ONE artifact the prose can be checked
against.

Per mechanism: census volume, clinical-trial count and share (NLM publication
types, assigned independently of this project), the anatomical sites it
concentrates in relative to each site's own weight, the mechanisms it most often
co-occurs with, and its 2015-to-2025 series.

TWO THINGS THIS DELIBERATELY DOES NOT DO. It does not rank mechanisms against
each other on volume, because descriptor breadth varies enormously -- 75% of
`epigenetic` comes from `DNA Methylation`, carried by any paper MEASURING
methylation -- so a cross-mechanism volume ranking is substantially a ranking of
how broad each descriptor is. And it does not report a co-occurrence RATE, for
the reason Section 3.13 sets out: that rate is a property of the labelling
instrument rather than of how often researchers combine mechanisms. Partner
ORDERINGS are stable under both instruments and are what is reported.
"""
import argparse
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RECORDS = REPO / "corpus/atlas/records"
SITE_MAP = REPO / "analysis/site-descriptor-map.tsv"
MECH_MAP = REPO / "analysis/mesh-mechanism-map.yaml"
OUT_MD = REPO / "analysis/census-mechanism-profile.md"
OUT_JSON = REPO / "analysis/census-mechanism-profile.json"
TRIAL_TYPES = {
    "Clinical Trial", "Randomized Controlled Trial", "Controlled Clinical Trial",
    "Clinical Trial, Phase I", "Clinical Trial, Phase II",
    "Clinical Trial, Phase III", "Clinical Trial, Phase IV",
    "Pragmatic Clinical Trial", "Adaptive Clinical Trial",
}
START, END = 2015, 2025
TOP_N = 6


def load_sites() -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    for ln in SITE_MAP.read_text(encoding="utf-8").splitlines():
        if ln.startswith("#") or not ln.strip():
            continue
        p = ln.split("\t")
        if len(p) >= 3:
            out[p[0]].add(p[2].strip().lower())
    return dict(out)


def scan(stride: int = 1) -> dict:
    import yaml

    mp = yaml.safe_load(MECH_MAP.read_text(encoding="utf-8"))["mechanisms"]
    mech = {k: {x.lower() for x in v["descriptors"]} for k, v in mp.items()}
    sites = load_sites()

    count: Counter = Counter()
    trials: Counter = Counter()
    by_year: dict[str, Counter] = defaultdict(Counter)
    by_site: dict[str, Counter] = defaultdict(Counter)
    partners: dict[str, Counter] = defaultdict(Counter)
    site_tot: Counter = Counter()
    n = 0
    for f in sorted(RECORDS.glob("*.jsonl.gz"))[::stride]:
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                n += 1
                ms = {m.lower() for m in (r.get("mesh") or [])}
                if not ms:
                    continue
                hit_sites = [s for s, d in sites.items() if ms & d]
                for s in hit_sites:
                    site_tot[s] += 1
                hits = [k for k, d in mech.items() if ms & d]
                if not hits:
                    continue
                is_trial = bool(set(r.get("pub_types") or []) & TRIAL_TYPES)
                y = r.get("year")
                for k in hits:
                    count[k] += 1
                    if is_trial:
                        trials[k] += 1
                    if isinstance(y, int):
                        by_year[k][y] += 1
                    for s in hit_sites:
                        by_site[k][s] += 1
                    for other in hits:
                        if other != k:
                            partners[k][other] += 1
    return {
        "census": n,
        "site_totals": dict(site_tot),
        "count": dict(count),
        "trials": dict(trials),
        "by_year": {k: {str(y): v[y] for y in sorted(v)} for k, v in by_year.items()},
        "by_site": {k: dict(v) for k, v in by_site.items()},
        "partners": {k: dict(v) for k, v in partners.items()},
    }


def assemble(d: dict) -> dict:
    st = d["site_totals"]
    base_tot = sum(st.values()) or 1
    rows = []
    for k, n in sorted(d["count"].items(), key=lambda x: -x[1]):
        sites = d["by_site"].get(k, {})
        assigned = sum(sites.values())
        # Enrichment, not raw rank: the raw ordering of sites within a mechanism
        # mostly reproduces the ordering of the sites themselves, which says
        # nothing about the mechanism.
        enr = sorted(
            (
                {
                    "site": s,
                    "n": v,
                    "enrichment": round((v / assigned) / (st[s] / base_tot), 2),
                }
                for s, v in sites.items()
                if v >= 20
            ),
            key=lambda r: -r["enrichment"],
        )
        yr = d["by_year"].get(k, {})
        a, b = yr.get(str(START), 0), yr.get(str(END), 0)
        rows.append({
            "mechanism": k,
            "census": n,
            "trials": d["trials"].get(k, 0),
            "trial_share": round(100 * d["trials"].get(k, 0) / n, 2) if n else None,
            "site_assigned": assigned,
            "top_sites": enr[:TOP_N],
            "top_partners": sorted(
                ({"mechanism": p, "n": v} for p, v in d["partners"].get(k, {}).items()),
                key=lambda r: -r["n"],
            )[:TOP_N],
            "start": a, "end": b,
            "growth": round(b / a, 2) if a >= 30 else None,
        })
    out = dict(d)
    out["rows"] = rows
    out["start_year"], out["end_year"] = START, END
    return out


def render(d: dict) -> str:
    L = ["# Per-mechanism census profile\n"]
    L.append(
        f"Generated by `scripts/census_mechanism_profile.py` over "
        f"{d['census']:,} census records. Mechanisms are labelled by MeSH "
        f"descriptor, sites by NLM's C04 tree, and trial status by NLM "
        f"publication type -- none of the three assigned by this project.\n"
    )
    L.append(
        "Volume is NOT comparable across mechanisms and no cross-mechanism "
        "ranking is drawn from it: descriptor breadth varies enormously, so a "
        "volume ordering is substantially an ordering of how broad each "
        "descriptor is. Trial share, site enrichment and partner ordering are "
        "the columns that survive that objection, because each is a ratio "
        "within a mechanism or a comparison against the mechanism's own base.\n"
    )
    L.append(f"| mechanism | census | trials | share | {d['start_year']} | "
             f"{d['end_year']} | growth |")
    L.append("|---|--:|--:|--:|--:|--:|--:|")
    for r in d["rows"]:
        g = f"x{r['growth']}" if r["growth"] is not None else "n/a"
        L.append(f"| {r['mechanism']} | {r['census']:,} | {r['trials']:,} | "
                 f"{r['trial_share']}% | {r['start']:,} | {r['end']:,} | {g} |")
    L.append("")
    L.append(
        "A growth ratio is reported only where the start year holds at least 30 "
        "articles; below that it measures the handful.\n"
    )
    for r in d["rows"]:
        L.append(f"## {r['mechanism']}\n")
        L.append(
            f"{r['census']:,} census articles, {r['trials']:,} carrying a "
            f"clinical-trial publication type ({r['trial_share']}%). "
            f"{r['site_assigned']:,} are assignable to a site.\n"
        )
        if r["top_sites"]:
            L.append("Concentrates in (enrichment against the site's own share "
                     "of site-assigned records): "
                     + ", ".join(f"{s['site']} {s['enrichment']}x ({s['n']:,})"
                                 for s in r["top_sites"]) + ".\n")
        else:
            L.append("No site holds 20 or more of its articles, so it has no "
                     "measurable anatomical concentration here.\n")
        if r["top_partners"]:
            L.append("Most frequent co-occurring mechanisms: "
                     + ", ".join(f"{p['mechanism']} ({p['n']:,})"
                                 for p in r["top_partners"]) + ".\n")
        else:
            L.append("No mechanism co-occurs with it in the census.\n")
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--render-only", action="store_true")
    a = ap.parse_args()
    if a.render_only:
        d = json.loads(OUT_JSON.read_text())
    else:
        d = assemble(scan(a.stride))
        OUT_JSON.write_text(json.dumps(d, indent=1) + "\n")
    OUT_MD.write_text(render(d))
    print(f"wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
