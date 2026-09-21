#!/usr/bin/env python3
"""What full text could and could not fix in the census's design labels.

44.5% of the census carries no design-informative NLM label -- 1,958,401
records that are neither trial, nor patient study, nor animal, nor cell
culture, nor review. That is the largest single hole in the study-design
layer, and the obvious response is to read the papers: the census ships an
open-access full-text layer holding 1,116,481 articles.

THE OBVIOUS RESPONSE IS BOUNDED BEFORE IT IS ATTEMPTED, which is the point of
this script. Measuring the ceiling costs one pass and settles whether the work
is worth doing; building the classifier first and discovering the bound
afterwards is how a project spends months closing a fifth of a gap.

THREE THINGS IT MEASURES, and each one narrows the answer further:

1. THE CEILING. A record is reachable only if it has a PMC identifier. Whatever
   classifier is built, however good, it cannot label a paper it cannot read.

2. THE DESIGN SKEW. Open-access availability is not independent of study
   design -- the whole reason a design-label gap exists is that some kinds of
   work are published differently from others. So a full-text correction does
   not sample the hole; it samples the part of the hole that happens to be
   readable, and the reachable rate differs several-fold across the classes
   that ARE labelled.

3. THE ERA SKEW, which is larger than the design skew and runs the same way as
   every other coverage measure in this project. Open access is a recent
   arrangement; the undetermined pile is not.

The consequence is not that full text is useless. It is that a distribution
recovered from it describes the READABLE literature, and quoting it as a
correction to the census would silently swap one population for another --
this repo's recurring defect, in a new place.
"""
import argparse
import importlib.util
import json
import sys
import statistics
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from census_input import atlas_root, iter_census_records  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
RECORDS = atlas_root() / "records"
OUT_MD = REPO / "analysis/census-fulltext-ceiling.md"
OUT_JSON = REPO / "analysis/census-fulltext-ceiling.json"
ERA_SPLIT = 2000


def _classifier():
    spec = importlib.util.spec_from_file_location(
        "ced", REPO / "scripts/census_evidence_design.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m.classify


def scan(stride: int = 1) -> dict:
    classify = _classifier()
    total: Counter = Counter()
    reachable: Counter = Counter()
    und_year: Counter = Counter()
    und_year_reach: Counter = Counter()
    n = 0
    for r in iter_census_records(RECORDS, stride):
        n += 1
        c = classify(r.get("pub_types"), r.get("mesh"))
        has = bool(r.get("pmcid"))
        total[c] += 1
        if has:
            reachable[c] += 1
        if c == "undetermined":
            y = r.get("year")
            if isinstance(y, int):
                und_year[y] += 1
                if has:
                    und_year_reach[y] += 1
    return {
        "census": n,
        "by_class_total": dict(total),
        "by_class_reachable": dict(reachable),
        "undetermined_by_year": {str(y): [und_year[y], und_year_reach[y]]
                                 for y in sorted(und_year)},
    }


def _median_year(pairs, which: int) -> int | None:
    """Median year weighted by count. `which` selects total (0) or reachable (1)."""
    years = []
    for y, v in pairs.items():
        years.extend([int(y)] * v[which])
    return int(statistics.median(years)) if years else None


def assemble(d: dict) -> dict:
    tot, reach = d["by_class_total"], d["by_class_reachable"]
    rows = []
    for c in sorted(tot, key=lambda k: -tot[k]):
        rows.append({
            "class": c, "total": tot[c], "reachable": reach.get(c, 0),
            "rate": round(100 * reach.get(c, 0) / tot[c], 1) if tot[c] else None,
        })
    # Zero reachability is an observed rate, not a missing design class.
    rated = [r for r in rows
             if r["class"] != "undetermined" and r["rate"] is not None]
    hi = max(rated, key=lambda r: r["rate"]) if len(rated) >= 2 else None
    lo = min(rated, key=lambda r: r["rate"]) if len(rated) >= 2 else None
    uy = d["undetermined_by_year"]
    old = {y: v for y, v in uy.items() if int(y) < ERA_SPLIT}
    new = {y: v for y, v in uy.items() if int(y) >= ERA_SPLIT}

    def rate(sub):
        t = sum(v[0] for v in sub.values())
        p = sum(v[1] for v in sub.values())
        return (t, p, round(100 * p / t, 1) if t else None)

    old_t, old_p, old_r = rate(old)
    new_t, new_p, new_r = rate(new)
    u = next((r for r in rows if r["class"] == "undetermined"),
             {"total": 0, "reachable": 0, "rate": None})
    out = dict(d)
    out["rows"] = rows
    out["ceiling_records"] = u["reachable"]
    out["ceiling_share_of_undetermined"] = u["rate"]
    out["ceiling_share_of_census"] = (
        round(100 * u["reachable"] / d["census"], 1) if d["census"] else None)
    out["unreachable_records"] = u["total"] - u["reachable"]
    out["design_skew"] = {
        "highest": hi, "lowest": lo,
        "fold": round(hi["rate"] / lo["rate"], 1) if lo and lo["rate"] else None,
    }
    out["era_skew"] = {
        "split": ERA_SPLIT,
        "before": {"total": old_t, "reachable": old_p, "rate": old_r},
        "since": {"total": new_t, "reachable": new_p, "rate": new_r},
        "fold": round(new_r / old_r, 1) if old_r and new_r is not None else None,
        "before_share_of_undetermined": (
            round(100 * old_t / (old_t + new_t), 1) if old_t + new_t else None),
    }
    out["median_year_pile"] = _median_year(uy, 0)
    out["median_year_reachable"] = _median_year(uy, 1)
    return out


def render(d: dict) -> str:
    ds, es = d["design_skew"], d["era_skew"]
    L = ["# What full text could and could not fix\n"]
    L.append(
        f"Generated by `scripts/census_fulltext_ceiling.py` over "
        f"{d['census']:,} census records. A record is REACHABLE if it carries a "
        f"PMC identifier; no classifier, however good, can label a paper it "
        f"cannot read.\n"
    )
    L.append("## The ceiling\n")
    u = next((r for r in d["rows"] if r["class"] == "undetermined"), None)
    if u and u["total"]:
        L.append(
            f"**{d['ceiling_records']:,} of {u['total']:,} undetermined records are "
            f"reachable ({d['ceiling_share_of_undetermined']}%).** The other "
            f"{d['unreachable_records']:,} lack a PMC identifier. Even a perfect "
            f"classifier operating through PMC could reach at most "
            f"{d['ceiling_share_of_undetermined']}% of this design-label gap, "
            f"or {d['ceiling_share_of_census']}% of the census.\n")
    else:
        L.append(
            "No undetermined records were observed in this input. The number "
            "available for design-label recovery is zero; its share of an "
            "undetermined population is unavailable.\n")
    L.append("| class | records | reachable | rate |")
    L.append("|---|--:|--:|--:|")
    for r in d["rows"]:
        L.append(f"| {r['class']} | {r['total']:,} | {r['reachable']:,} | "
                 f"{str(r['rate']) + '%' if r['rate'] is not None else 'unavailable'} |")
    L.append("")
    L.append("## Reachability by design and era\n")
    if ds["highest"] is None:
        L.append("**Design.** A comparison is unavailable: fewer than two "
                 "labelled design classes were observed.\n")
    else:
        L.append(
            f"**Design.** Reachability runs from {ds['highest']['class']} at "
            f"{ds['highest']['rate']}% to {ds['lowest']['class']} at "
            f"{ds['lowest']['rate']}%"
            + (f", a factor of {ds['fold']}.\n" if ds["fold"] is not None else
               ". The fold ratio is unavailable because the lowest rate is zero.\n"))
    old_rate, new_rate = es["before"]["rate"], es["since"]["rate"]
    if old_rate is None or new_rate is None:
        L.append("**Era.** A comparison is unavailable: dated undetermined "
                 "records were not observed on both sides of "
                 f"{es['split']}.\n")
    else:
        L.append(
            f"**Era.** Among dated undetermined records the rate is "
            f"{old_rate}% before {es['split']} against {new_rate}% since"
            + (f", a factor of {es['fold']}" if es["fold"] is not None else
               "; the fold ratio is unavailable because the earlier rate is zero")
            + f". {es['before_share_of_undetermined']}% of the dated "
              "undetermined records belong to the older era.\n")
    pile_year, reach_year = d["median_year_pile"], d["median_year_reachable"]
    if pile_year is not None and reach_year is not None:
        delta = reach_year - pile_year
        shift = (f"{delta} years later" if delta >= 0 else f"{-delta} years earlier")
        L.append(f"The median year of the dated undetermined records is {pile_year}; "
                 f"the reachable subset has median year {reach_year}, {shift}.\n")
    else:
        L.append("The median-year shift is unavailable without dated "
                 "undetermined records in both the complete and reachable sets.\n")
    L.append("## What this licenses\n")
    L.append(
        "Full text is worth reading for what it can answer directly. What it "
        "cannot do by itself is CORRECT the census's design distribution: a share "
        "recovered from reachable records describes that subset. It cannot "
        "automatically be generalized to records outside it. If a full-text distribution is ever "
        "published here it belongs in its own column, against its own "
        "denominator, never merged into an NLM-labelled one.\n"
    )
    L.append(
        "The ceiling is also an UPPER bound in a second way. Carrying a PMC "
        "identifier means a record is in PMC, not that its full text is in the "
        "open-access subset this project holds, nor that the text parsed. Every "
        "figure here is therefore the most optimistic version of the answer.\n"
    )
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stride", type=int, default=1,
                    help="sample every Nth shard; shards are CHRONOLOGICAL, so "
                         "a prefix samples one era and would destroy the era "
                         "measurement this script exists for")
    ap.add_argument("--render-only", action="store_true")
    a = ap.parse_args()
    if a.render_only:
        d = assemble(json.loads(OUT_JSON.read_text()))
    else:
        d = assemble(scan(a.stride))
    json_text = json.dumps(d, indent=1) + "\n"
    md_text = render(d)
    OUT_JSON.write_text(json_text)
    OUT_MD.write_text(md_text)
    print(f"wrote {OUT_MD}")
    print(f"  ceiling {d['ceiling_records']:,} "
          f"({d['ceiling_share_of_undetermined']}% of undetermined)")
    print(f"  design skew {d['design_skew']['fold']}x, "
          f"era skew {d['era_skew']['fold']}x")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
