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
how broad each descriptor is. And it does not report a co-occurrence RATE as
evidence of tested combinations: co-tagging depends on the labelling instrument
and does not show that mechanisms were tested together. Partner orderings rank
raw co-occurrence counts in the descriptor-selected articles. Trial share,
site enrichment and partner orderings can all change with descriptor coverage;
normalization alone does not establish comparability. This profile does not
compare partner rankings from different labelling methods on the same articles.

The raw census is not committed. Set FERRO_ATLAS_ROOT when it lives outside
corpus/atlas/, or use --render-only to regenerate from the committed counts
without reading raw records. A missing or empty input never replaces a report.
"""
import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from census_input import atlas_root, census_shards, iter_census_shards  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
RECORDS = atlas_root() / "records"
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
    for line_number, ln in enumerate(
            SITE_MAP.read_text(encoding="utf-8").splitlines(), start=1):
        if ln.startswith("#") or not ln.strip():
            continue
        p = ln.split("\t")
        if len(p) < 3 or any(not field.strip() for field in p[:3]):
            raise SystemExit(
                f"Invalid site map row at {SITE_MAP}:{line_number}; expected "
                "site, tree root and descriptor. Existing reports were not changed."
            )
        out[p[0]].add(p[2].strip().lower())
    if not out:
        raise SystemExit(
            f"No site descriptors found at {SITE_MAP}; existing reports were not changed."
        )
    return dict(out)


def scan(stride: int = 1) -> dict:
    import yaml

    shards = census_shards(RECORDS, stride)

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
    for r in iter_census_shards(shards, RECORDS):
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
    for k, n in sorted(d["count"].items(), key=lambda x: (-x[1], x[0])):
        sites = d["by_site"].get(k, {})
        assigned = sum(sites.values())
        # Enrichment compares the mechanism's site-assignment shares with the
        # census baseline. This normalization does not remove selection effects
        # from the mechanism's descriptor coverage.
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
            # Same total-order requirement: equal enrichment otherwise falls
            # back to dict order, which decides both the ranking AND which
            # rows survive the `[:TOP_N]` cut.
            key=lambda r: (-r["enrichment"], r["site"]),
        )
        yr = d["by_year"].get(k, {})
        a, b = yr.get(str(START), 0), yr.get(str(END), 0)
        rows.append({
            "mechanism": k,
            "census": n,
            "trials": d["trials"].get(k, 0),
            "trial_share": round(100 * d["trials"].get(k, 0) / n, 2) if n else None,
            # Legacy field name: this is a sum of site assignments, not a
            # unique article count. One article may match several sites.
            "site_assigned": assigned,
            "top_sites": enr[:TOP_N],
            "top_partners": sorted(
                ({"mechanism": p, "n": v} for p, v in d["partners"].get(k, {}).items()),
                # Tie-break by name so the ranking is TOTAL. The first
                # attempt keyed on `r.get("name", r.get("site", ""))`, and
                # these rows carry neither field, so the second key was the
                # empty string for every row and the fix was inert while
                # carrying a comment saying it worked.
                key=lambda r: (-r["n"], r["mechanism"]),
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
        f"{d['census']:,} census records. Articles carry NLM-assigned MeSH "
        f"descriptors and publication types. Mechanisms use the project's "
        f"curated descriptor map; sites use its selected roots in NLM's C04 "
        f"tree; trial status uses the configured publication-type set.\n"
    )
    L.append(
        "Volume is NOT comparable across mechanisms and no cross-mechanism "
        "ranking is drawn from it: descriptor breadth varies enormously, so a "
        "volume ordering is substantially an ordering of how broad each "
        "descriptor is. Trial share and site enrichment describe articles "
        "selected by each mechanism's descriptors. Partner ordering ranks raw "
        "co-occurrence counts; it is not adjusted for partner prevalence. All "
        "three summaries can change with descriptor coverage and indexing. "
        "Normalization does not establish comparability across mechanism "
        "definitions, and this profile does not test rank stability across "
        "labeling methods. A separate [paired frozen-corpus comparison]"
        "(paired-partner-agreement.md) measures agreement on the same articles "
        "and common groups; it does not establish census-wide ranking stability.\n"
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
    L.append(
        "Site enrichment divides a site's share of a mechanism's site assignments "
        "by its share of all census site assignments. Each article contributes "
        "once to each matching site; sites can overlap, so assignment totals are "
        "not unique article counts. The legacy JSON field `site_assigned` stores "
        "this assignment total. Trial shares retain an article denominator.\n"
    )
    for r in d["rows"]:
        L.append(f"## {r['mechanism']}\n")
        L.append(
            f"{r['census']:,} census articles, {r['trials']:,} carrying a "
            f"clinical-trial publication type ({r['trial_share']}%). "
            f"These articles contribute {r['site_assigned']:,} site assignments.\n"
        )
        if r["top_sites"]:
            L.append("Highest site enrichments (relative to the site's share "
                     "of all census site assignments): "
                     + ", ".join(f"{s['site']} {s['enrichment']}x ({s['n']:,})"
                                 for s in r["top_sites"]) + ".\n")
        else:
            L.append("No site meets the 20-article reporting threshold, so no "
                     "site-enrichment ranking is reported for these records.\n")
        if r["top_partners"]:
            L.append("Most frequent co-occurring mechanisms: "
                     + ", ".join(f"{p['mechanism']} ({p['n']:,})"
                                 for p in r["top_partners"]) + ".\n")
        else:
            L.append("No co-occurrences with other mapped mechanisms were "
                     "observed in these records.\n")
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stride", type=int, default=1,
                    help="read every Nth sorted shard (default: all shards)")
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
    # Complete both representations before replacing either published report.
    json_text = json.dumps(d, indent=1) + "\n"
    md_text = render(d)
    OUT_JSON.write_text(json_text)
    OUT_MD.write_text(md_text)
    print(f"wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
