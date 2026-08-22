#!/usr/bin/env python3
"""How large are the literatures MeSH cannot name?

This book reports nine modalities as *not measurable* rather than as zero,
which is the right call -- a zero invites the reading that nobody works on
something, and TTFields has FDA approval in two indications. But "not
measurable" is a shrug unless somebody says how big the thing is, and a reader
meeting the phrase nine times is entitled to ask.

They are measurable, on one instrument. This project's keyword vocabulary reads
title and abstract, so it can count a modality whose name authors write even
when NLM has no descriptor for it. That is what the phrase should have meant
all along: not measurable BY DESCRIPTOR, which is a statement about MeSH.

WHAT MAKES THESE ROWS WEAKER THAN THE REST, and it is not a small caveat. Every
other modality in this project can be measured twice, and the agreement between
the two arms is a descriptor-validity test that caught a 15-fold artifact in
sonodynamic therapy's trial share. For these nine that test is unavailable BY
CONSTRUCTION -- there is no second arm, because the absence of a second arm is
what puts them on this list. So their numbers rest on a single instrument whose
mechanism precision is 82.5%, with no independent check available, and they
should be read as less certain than the descriptor-backed figures beside them,
not more.

TWO KNOWN FAILURE MODES OF THAT INSTRUMENT are visible here rather than
hypothetical, and both are reported per row: a term that is a common English
word or an abbreviation collides (`echt` inside "liechtenstein"), and a term
that predates the modality it now names counts old work (`copper ionophore`
applied to disulfiram in 1982, four decades before cuproptosis was named).
"""
import argparse
import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RECORDS = REPO / "corpus/atlas/records"
MECH_MAP = REPO / "analysis/mesh-mechanism-map.yaml"
OUT_MD = REPO / "analysis/census-unnamed-modalities.md"
OUT_JSON = REPO / "analysis/census-unnamed-modalities.json"

TRIAL_TYPES = {
    "Clinical Trial", "Randomized Controlled Trial", "Controlled Clinical Trial",
    "Clinical Trial, Phase I", "Clinical Trial, Phase II",
    "Clinical Trial, Phase III", "Clinical Trial, Phase IV",
    "Pragmatic Clinical Trial", "Adaptive Clinical Trial",
}
LATE = {"Clinical Trial, Phase III", "Clinical Trial, Phase IV"}
# A modality whose count is carried by pre-1990 records is almost certainly
# matching an older use of its own words, since none of these modalities
# existed then.
VINTAGE_YEAR = 1990


def _arms():
    import sys
    import yaml

    sys.path.insert(0, str(REPO / "scripts"))
    import config

    mp = yaml.safe_load(MECH_MAP.read_text(encoding="utf-8"))["mechanisms"]
    mesh = {k.lower() for k, v in mp.items() if v["descriptors"]}
    text = {k.lower(): v for k, v in config.MECHANISM_KEYWORDS.items() if v}
    invisible = {k: v for k, v in text.items() if k not in mesh}
    return invisible


def scan(stride: int = 1) -> dict:
    invisible = _arms()
    pats = {k: re.compile("|".join(rf"\b{re.escape(t.lower())}\b" for t in v))
            for k, v in invisible.items()}
    pre = re.compile("|".join(p.pattern for p in pats.values()))

    stats = defaultdict(Counter)
    years = defaultdict(Counter)
    n = 0
    for f in sorted(RECORDS.glob("*.jsonl.gz"))[::stride]:
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                n += 1
                blob = f"{r.get('title') or ''} {r.get('abstract') or ''}".lower()
                if not blob or not pre.search(blob):
                    continue
                pts = set(r.get("pub_types") or [])
                y = r.get("year")
                for k, p in pats.items():
                    if not p.search(blob):
                        continue
                    stats[k]["articles"] += 1
                    if pts & TRIAL_TYPES:
                        stats[k]["trials"] += 1
                    if pts & LATE:
                        stats[k]["late_phase"] += 1
                    if isinstance(y, int):
                        years[k][y] += 1
    return {
        "census": n,
        "vintage_year": VINTAGE_YEAR,
        "terms": {k: sorted(v) for k, v in invisible.items()},
        "stats": {k: dict(v) for k, v in stats.items()},
        "years": {k: {str(y): c for y, c in sorted(v.items())}
                  for k, v in years.items()},
    }


def assemble(d: dict) -> dict:
    rows = []
    for k in sorted(d["stats"]):
        s = d["stats"][k]
        ys = {int(y): c for y, c in d["years"].get(k, {}).items()}
        total = sum(ys.values())
        pre = sum(c for y, c in ys.items() if y < d["vintage_year"])
        rows.append({
            "modality": k,
            "articles": s.get("articles", 0),
            "trials": s.get("trials", 0),
            "late_phase": s.get("late_phase", 0),
            "trial_share": (round(100 * s.get("trials", 0) / s["articles"], 2)
                            if s.get("articles") else None),
            "first_year": min(ys) if ys else None,
            "median_year": (sorted(y for y in ys for _ in range(ys[y]))[total // 2]
                            if total else None),
            "pre_vintage": pre,
            "pre_vintage_share": round(100 * pre / total, 1) if total else None,
            "n_terms": len(d["terms"].get(k, [])),
        })
    rows.sort(key=lambda r: -r["articles"])
    # A modality whose literature substantially predates its own existence is
    # matching an older use of its words, not counting itself.
    for r in rows:
        r["vintage_suspect"] = bool(r["pre_vintage_share"]
                                    and r["pre_vintage_share"] >= 10)
    out = dict(d)
    out["rows"] = rows
    out["suspect"] = sorted(r["modality"] for r in rows if r["vintage_suspect"])
    out["total_articles"] = sum(r["articles"] for r in rows)
    out["total_trials"] = sum(r["trials"] for r in rows)
    return out


def render(d: dict) -> str:
    L = ["# The literatures MeSH cannot name\n"]
    L.append(
        f"Generated by `scripts/census_unnamed_modalities.py` over "
        f"{d['census']:,} census records. Nine modalities are reported "
        f"throughout this project as *not measurable* rather than as zero, "
        f"because MeSH has no descriptor for them. That is the right call and "
        f"an incomplete one: a reader meeting the phrase nine times is "
        f"entitled to ask how big the thing is.\n"
    )
    L.append(
        f"They are measurable on one instrument. This project's keyword "
        f"vocabulary reads title and abstract, so it counts a modality whose "
        f"name authors write even where NLM has no term for it. Together these "
        f"nine carry **{d['total_articles']:,} articles and "
        f"{d['total_trials']:,} clinical trials** -- which is what *not "
        f"measurable by descriptor* actually amounts to.\n"
    )
    L.append("| modality | articles | trials | trial share | Phase III/IV | "
             "first | median | pre-1990 |")
    L.append("|---|--:|--:|--:|--:|--:|--:|--:|")
    for r in d["rows"]:
        mark = " *" if r["vintage_suspect"] else ""
        L.append(
            f"| {r['modality']}{mark} | {r['articles']:,} | {r['trials']:,} | "
            f"{r['trial_share']}% | {r['late_phase']:,} | {r['first_year']} | "
            f"{r['median_year']} | {r['pre_vintage_share']}% |")
    L.append("")
    L.append("## Why these rows are weaker than the ones beside them\n")
    L.append(
        "Every other modality in this project can be measured twice, and the "
        "agreement between a descriptor arm and a keyword arm is a validity "
        "test -- it caught a 15-fold artifact in sonodynamic therapy's trial "
        "share, where a descriptor was counting ultrasound hyperthermia. "
        "**For these nine that test is unavailable BY CONSTRUCTION**: there is "
        "no second arm, and the absence of one is exactly what puts a modality "
        "on this list.\n"
    )
    L.append(
        "So these numbers rest on a single instrument whose measured mechanism "
        "precision is 82.5%, with no independent check available. They should "
        "be read as LESS certain than the descriptor-backed figures elsewhere "
        "in this project, not more -- which is the opposite of how a "
        "precise-looking count in a table tends to read.\n"
    )
    if d["suspect"]:
        L.append("## The failure mode is visible, not hypothetical\n")
        L.append(
            f"{len(d['suspect'])} modalities have {d['vintage_year']}-and-"
            f"earlier literature above a tenth of their total, marked * above: "
            + ", ".join(f"`{m}`" for m in d["suspect"])
            + ". None of these modalities existed then, so those records are "
              "the vocabulary matching an older use of its own words.\n"
        )
        L.append(
            "Two mechanisms, both already traced elsewhere in this project. A "
            "term can be a common word or a short abbreviation that collides "
            "-- `echt` matches inside \"liechtenstein\". And a term can predate "
            "the modality it now names: `copper ionophore` was applied to "
            "disulfiram in 1982, four decades before cuproptosis was named, so "
            "`cuproptosis` counts 1980s chemistry as its own literature.\n"
        )
        L.append(
            "The affected counts are upper bounds. They are published rather "
            "than filtered because a date cut-off would be a second judgement "
            "layered on the first, and the per-row pre-1990 share lets a "
            "reader apply their own.\n"
        )
    L.append("## What this changes about the phrase\n")
    L.append(
        "*Not measurable* should mean **not measurable by descriptor**, which "
        "is a statement about MeSH rather than about the field. Used without "
        "that qualifier it invites exactly the reading it was meant to "
        "prevent -- that there is nothing there. There is: radioligand therapy "
        "alone carries hundreds of indexed trials while having no MeSH "
        "descriptor of its own.\n"
    )
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--render-only", action="store_true")
    a = ap.parse_args()
    d = assemble(json.loads(OUT_JSON.read_text()) if a.render_only
                 else scan(a.stride))
    OUT_JSON.write_text(json.dumps(d, indent=1) + "\n")
    OUT_MD.write_text(render(d))
    print(f"wrote {OUT_MD}")
    print(f"  {len(d['rows'])} modalities, {d['total_articles']:,} articles, "
          f"{d['total_trials']:,} trials, {len(d['suspect'])} vintage-suspect")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
