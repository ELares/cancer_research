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
import gzip
import importlib.util
import json
import statistics
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RECORDS = REPO / "corpus/atlas/records"
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
    for f in sorted(RECORDS.glob("*.jsonl.gz"))[::stride]:
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
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
    rated = [r for r in rows if r["class"] != "undetermined" and r["rate"]]
    hi = max(rated, key=lambda r: r["rate"])
    lo = min(rated, key=lambda r: r["rate"])
    uy = d["undetermined_by_year"]
    old = {y: v for y, v in uy.items() if int(y) < ERA_SPLIT}
    new = {y: v for y, v in uy.items() if int(y) >= ERA_SPLIT}

    def rate(sub):
        t = sum(v[0] for v in sub.values())
        p = sum(v[1] for v in sub.values())
        return (t, p, round(100 * p / t, 1) if t else None)

    old_t, old_p, old_r = rate(old)
    new_t, new_p, new_r = rate(new)
    u = next(r for r in rows if r["class"] == "undetermined")
    out = dict(d)
    out["rows"] = rows
    out["ceiling_records"] = u["reachable"]
    out["ceiling_share_of_undetermined"] = u["rate"]
    out["ceiling_share_of_census"] = round(100 * u["reachable"] / d["census"], 1)
    out["unreachable_records"] = u["total"] - u["reachable"]
    out["design_skew"] = {
        "highest": hi, "lowest": lo,
        "fold": round(hi["rate"] / lo["rate"], 1) if lo["rate"] else None,
    }
    out["era_skew"] = {
        "split": ERA_SPLIT,
        "before": {"total": old_t, "reachable": old_p, "rate": old_r},
        "since": {"total": new_t, "reachable": new_p, "rate": new_r},
        "fold": round(new_r / old_r, 1) if old_r else None,
        "before_share_of_undetermined": round(100 * old_t / (old_t + new_t), 1),
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
    u = next(r for r in d["rows"] if r["class"] == "undetermined")
    L.append(
        f"**{d['ceiling_records']:,} of {u['total']:,} undetermined records are "
        f"reachable ({d['ceiling_share_of_undetermined']}%).** The other "
        f"{d['unreachable_records']:,} are not, and no amount of classifier "
        f"work changes that. Reading every open-access paper in the census "
        f"perfectly would close a fifth of the design-label gap and leave "
        f"four-fifths exactly where it is -- worth "
        f"{d['ceiling_share_of_census']}% of the census.\n"
    )
    L.append("| class | records | reachable | rate |")
    L.append("|---|--:|--:|--:|")
    for r in d["rows"]:
        L.append(f"| {r['class']} | {r['total']:,} | {r['reachable']:,} | "
                 f"{r['rate']}% |")
    L.append("")
    L.append("## The two skews, and why they matter more than the ceiling\n")
    L.append(
        f"**Design.** Reachability runs from {ds['highest']['class']} at "
        f"{ds['highest']['rate']}% down to {ds['lowest']['class']} at "
        f"{ds['lowest']['rate']}%, a factor of {ds['fold']}. Open-access "
        f"availability is not independent of study design -- which is the same "
        f"reason a design-label gap exists at all -- so a full-text pass does "
        f"not sample the hole. It samples the readable part of it, and the "
        f"readable part is enriched for the kinds of work that are already "
        f"best represented.\n"
    )
    L.append(
        f"**Era, and this one is larger.** Among undetermined records the rate "
        f"is {es['before']['rate']}% before {es['split']} against "
        f"{es['since']['rate']}% since -- a factor of {es['fold']} -- while "
        f"{es['before_share_of_undetermined']}% of the pile sits in the older "
        f"era. The median year of the undetermined pile is "
        f"{d['median_year_pile']}; the median year of the part full text could "
        f"reach is {d['median_year_reachable']}, "
        f"{d['median_year_reachable'] - d['median_year_pile']} years later.\n"
    )
    L.append("## What this licenses\n")
    L.append(
        "Full text is worth reading for what it can answer directly. What it "
        "cannot do is CORRECT the census's design distribution: a share "
        "recovered from the reachable fifth describes the readable literature, "
        "and reporting it as a correction would swap one population for "
        "another without saying so. If a full-text distribution is ever "
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
    OUT_JSON.write_text(json.dumps(d, indent=1) + "\n")
    OUT_MD.write_text(render(d))
    print(f"wrote {OUT_MD}")
    print(f"  ceiling {d['ceiling_records']:,} "
          f"({d['ceiling_share_of_undetermined']}% of undetermined)")
    print(f"  design skew {d['design_skew']['fold']}x, "
          f"era skew {d['era_skew']['fold']}x")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
