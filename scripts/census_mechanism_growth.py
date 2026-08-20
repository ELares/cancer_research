#!/usr/bin/env python3
"""Per-mechanism growth at census scale, and the denominator a growth claim needs.

`manuscript_vs_census.py` tests the manuscript's Section 3.7 growth claim by
comparing the retrieved corpus's growth (x30.71 over 2015-2025) against the
CANCER LITERATURE AS A WHOLE (x1.10) and against how much of each year is open
access (x1.67), concluding that neither confound accounts for the corpus's rise
and therefore that "the mechanisms the corpus tracks did grow far faster than
cancer literature as a whole".

THAT COMPARISON IS AGAINST THE WRONG DENOMINATOR, and the direction of the error
flatters the conclusion. A corpus assembled from queries about emerging
therapeutic mechanisms will outgrow all of cancer research whether or not
anything unusual happened, because the whole field includes epidemiology,
surgery, supportive care and decades of established practice. The matched
denominator is the growth of THOSE MECHANISMS in the census -- the same concepts,
labelled by MeSH rather than by keyword, over articles the project did not
select.

Measured here, that denominator is x2.91. So the claim does not simply survive or
fail: it splits. The mechanisms genuinely outgrew the field by a factor of about
2.6, which is the part of the manuscript's story that holds. The remaining factor
of roughly ten between x2.91 and x30.71 belongs to the retrieval, not to the
literature -- and reporting the corpus figure against x1.10 attributed all of it
to the mechanisms.
"""
import argparse
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RECORDS = REPO / "corpus/atlas/records"
MECH_MAP = REPO / "analysis/mesh-mechanism-map.yaml"
OUT_MD = REPO / "analysis/census-mechanism-growth.md"
OUT_JSON = REPO / "analysis/census-mechanism-growth.json"
START, END = 2015, 2025
# A second window, because the manuscript quotes one. 2015-2025 spans a decade
# and 2020-2025 is where the "post-2020 surge" claims live; a mechanism can
# look flat on one and steep on the other (epigenetic: x0.91 and +2%), so both
# are rendered rather than left for a reader to derive from the raw series.
RECENT_START = 2020
# Below this, a ratio is a statement about a handful of articles.
MIN_BASE = 30


def scan(stride: int = 1) -> dict:
    import yaml

    mp = yaml.safe_load(MECH_MAP.read_text(encoding="utf-8"))["mechanisms"]
    mech = {k: {x.lower() for x in v["descriptors"]} for k, v in mp.items()}
    by: dict[str, Counter] = defaultdict(Counter)
    union: Counter = Counter()
    field: Counter = Counter()
    n = 0
    shards = sorted(RECORDS.glob("*.jsonl.gz"))[::stride]
    for f in shards:
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                n += 1
                y = r.get("year")
                if not isinstance(y, int):
                    continue
                field[y] += 1
                ms = {m.lower() for m in (r.get("mesh") or [])}
                if not ms:
                    continue
                hit = False
                for k, d in mech.items():
                    if ms & d:
                        by[k][y] += 1
                        hit = True
                if hit:
                    union[y] += 1
    return {
        "census": n,
        "shards": len(shards),
        "start_year": START,
        "end_year": END,
        "field_by_year": {str(y): field[y] for y in sorted(field)},
        "union_by_year": {str(y): union[y] for y in sorted(union)},
        "mechanism_by_year": {k: {str(y): v[y] for y in sorted(v)} for k, v in by.items()},
    }


def assemble(d: dict) -> dict:
    s, e = str(d["start_year"]), str(d["end_year"])
    field = d["field_by_year"]
    union = d["union_by_year"]

    def g(a, b):
        return round(b / a, 2) if a else None

    rows = []
    for k, v in d["mechanism_by_year"].items():
        a, b = v.get(s, 0), v.get(e, 0)
        rows.append({
            "mechanism": k, "start": a, "end": b,
            # WITHHELD below the base floor, not merely flagged. An earlier
            # version computed the ratio regardless and let `base_sufficient`
            # gate only the RENDER, so the .md withheld microbiome's x13.62
            # while the committed JSON carried it -- and a downstream reader
            # takes the JSON. A caveat that lives in one representation of an
            # artifact is not a caveat.
            "growth": g(a, b) if a >= MIN_BASE else None,
            "base_sufficient": a >= MIN_BASE,
        })
    rows.sort(key=lambda r: (-(r["growth"] or 0)) if r["base_sufficient"] else 1e9)
    rs = str(RECENT_START)
    for r in rows:
        v = d["mechanism_by_year"][r["mechanism"]]
        a2, b2 = v.get(rs, 0), v.get(e, 0)
        r["recent_start"] = a2
        r["recent_end"] = b2
        r["recent_pct"] = round(100 * (b2 - a2) / a2, 0) if a2 >= 20 else None
    field_growth = g(field.get(s, 0), field.get(e, 0))
    union_growth = g(union.get(s, 0), union.get(e, 0))
    out = dict(d)
    out["rows"] = rows
    out["field_growth"] = field_growth
    out["union_start"] = union.get(s, 0)
    out["union_end"] = union.get(e, 0)
    out["union_growth"] = union_growth
    out["field_start"] = field.get(s, 0)
    out["field_end"] = field.get(e, 0)
    out["mechanisms_over_field"] = (
        round(union_growth / field_growth, 2) if field_growth else None
    )
    out["min_base"] = MIN_BASE
    out["recent_start_year"] = RECENT_START
    # DERIVED, never asserted. A first draft of the manuscript called
    # electrochemical therapy "the one mechanism whose literature is
    # shrinking" -- an extremum claimed over a set nobody had enumerated, and
    # false: epigenetic declines too (x0.91). The set is computed here so the
    # prose beside it cannot name the wrong number of members.
    out["shrinking"] = [r["mechanism"] for r in rows
                        if r["growth"] is not None and r["growth"] < 1.0]
    out["recent_shrinking"] = [r["mechanism"] for r in rows
                               if r["recent_pct"] is not None and r["recent_pct"] < 0]
    return out


def render(d: dict) -> str:
    s, e = d["start_year"], d["end_year"]
    ok = [r for r in d["rows"] if r["base_sufficient"]]
    thin = [r for r in d["rows"] if not r["base_sufficient"]]
    L = []
    L.append("# Per-mechanism growth at census scale\n")
    L.append(
        f"Generated by `scripts/census_mechanism_growth.py` over {d['census']:,} "
        f"census records, {s} to {e}, with mechanisms labelled by MeSH descriptor "
        f"rather than by keyword.\n"
    )
    L.append("## The denominator a growth claim needs\n")
    L.append(f"| | {s} | {e} | growth |")
    L.append("|---|--:|--:|--:|")
    L.append(f"| cancer literature (census) | {d['field_start']:,} | "
             f"{d['field_end']:,} | x{d['field_growth']} |")
    L.append(f"| **articles carrying any mechanism descriptor** | "
             f"{d['union_start']:,} | {d['union_end']:,} | "
             f"**x{d['union_growth']}** |")
    L.append("")
    L.append(
        f"The mechanisms this project tracks grew "
        f"**x{d['mechanisms_over_field']}** faster than cancer literature as a "
        f"whole. That is a real result and it is the part of the manuscript's "
        f"growth story that survives: these are not simply riding a rising tide, "
        f"because there is no rising tide -- the field is close to flat.\n"
    )
    L.append(
        f"It is also the denominator `manuscript_vs_census.py` should have used "
        f"and did not. Comparing the retrieved corpus against the whole field "
        f"attributes ALL of the corpus's rise to the mechanisms, when the "
        f"mechanisms account for x{d['union_growth']} of it. A corpus built from "
        f"queries about emerging therapies outgrows all of cancer research "
        f"whether or not anything unusual happened; the question a growth claim "
        f"is asking is whether it outgrew the thing it is about.\n"
    )
    rs = d["recent_start_year"]
    L.append("## Per mechanism\n")
    L.append(f"| mechanism | {s} | {e} | growth | {rs} | {e} | change |")
    L.append("|---|--:|--:|--:|--:|--:|--:|")
    for r in ok:
        rp = (f"{r['recent_pct']:+.0f}%" if r["recent_pct"] is not None else "n/a")
        L.append(f"| {r['mechanism']} | {r['start']:,} | {r['end']:,} | "
                 f"x{r['growth']} | {r['recent_start']:,} | {r['recent_end']:,} | "
                 f"{rp} |")
    L.append("")
    shrink = d["shrinking"]
    if shrink:
        L.append(
            f"{len(shrink)} mechanism(s) are SMALLER in {e} than in {s}: "
            + ", ".join(f"`{m}`" for m in shrink)
            + ". The count is derived rather than described, because an "
            "extremum stated over a set nobody enumerated is how a second "
            "member goes unnoticed.\n"
        )
    else:
        L.append(f"No mechanism is smaller in {e} than in {s}.\n")
    if thin:
        L.append(
            f"{len(thin)} mechanism(s) start below the {d['min_base']}-article "
            f"base this table requires and are reported without a ratio, because "
            f"a ratio off a handful of articles measures the handful: "
            + ", ".join(
                f"{r['mechanism']} ({r['start']} to {r['end']})" for r in thin
            )
            + ".\n"
        )
    L.append("## What this does not say\n")
    L.append(
        "Growth in indexed articles is growth in attention, not in evidence, "
        "results or clinical progress. MeSH indexing also lags publication, so "
        f"the {e} column is under-counted and every ratio here is a lower bound; "
        "the lag applies to numerator and denominator alike, so the comparison "
        "between rows is less affected than any single row.\n"
    )
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--render-only", action="store_true")
    a = ap.parse_args()
    if a.render_only:
        # RE-ASSEMBLE rather than render the stored derived fields. The JSON
        # carries the raw per-year series, so every derived column can be
        # recomputed from it -- which means adding a column does not strand
        # --render-only against an artifact written before that column
        # existed, and means the stored derived fields are checkable against a
        # fresh derivation rather than merely trusted.
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
