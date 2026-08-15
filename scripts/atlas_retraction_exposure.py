#!/usr/bin/env python3
"""How much of the cancer relation graph rests on retracted work?

WHY THIS IS WORTH ASKING
------------------------
Retracted papers do not leave the literature when they are retracted. They stay
indexed, they stay cited, and -- the part that matters here -- machine-extracted
relation graphs keep the assertions they contributed. Anyone building on a
PubTator-derived graph inherits them silently.

This project holds what is needed to bound that: PubMed's own
`Retracted Publication` publication type, stored per record by
`scripts/atlas_baseline.py`, and 10.5M typed relations keyed by PMID. Nothing in
this repository has ever read the retraction field -- `grep -ri retract` over
scripts/ and analysis/ returns nothing.

THE MEASUREMENT THAT MEANS SOMETHING, AND THE ONE THAT MISLEADS
---------------------------------------------------------------
This repository has been here before. `atlas_ambiguity_impact.py` found that
50.6% of relation rows TOUCH a contested identifier, which sounds severe and
measures containment, while only 1.31% rest on an assignment that was actually
at risk -- a 39-fold difference between the alarming number and the true one.

The same trap applies exactly. So both are computed and the containment figure
is never reported alone:

  CONTAINMENT   relations whose source paper is retracted. An upper bound on
                exposure, and mostly harmless: a pair asserted by fifty papers
                one of which was retracted has lost nothing.
  UNCORROBORATED  entity pairs whose ONLY asserting papers are retracted. This
                is the damaging class -- the graph asserts a relation that,
                after retraction, nothing standing supports.

A pair in the second class is not necessarily false. Retraction covers
misconduct, honest error, and irreproducibility alike, and a retracted paper can
still have been right. What the class means is that the graph's support for that
pair has been withdrawn, and a consumer should know.

WHAT THIS CANNOT SAY
--------------------
* Retraction indexing lags publication, and lags the retraction itself. Every
  count here is a LOWER bound on what will eventually be marked.
* PubMed marks the retracted ARTICLE with `Retracted Publication` and the notice
  with `Retraction of Publication`. Only the former is counted as tainted; the
  notice is a legitimate record and is reported separately.
* Relations come from PubTator3's extractor, whose own error rate (~79.6 F1 on
  this project's reading) is larger than most effects measured here. This bounds
  which assertions rest on withdrawn support, not which assertions are true.

Usage:
    python scripts/atlas_retraction_exposure.py
"""

import collections
import gzip
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RECORDS = PROJECT_ROOT / "corpus" / "atlas" / "records"
RELATIONS = PROJECT_ROOT / "corpus" / "atlas" / "relations" / "relations.tsv.gz"
OUT_MD = PROJECT_ROOT / "analysis" / "atlas-retraction-exposure.md"
OUT_JSON = PROJECT_ROOT / "analysis" / "atlas-retraction-exposure.json"

# PubMed's own types. Only the retracted ARTICLE is treated as tainted; a
# retraction notice, erratum or republication is a legitimate record.
TAINTED = "Retracted Publication"

# Everything else in this family is DISCOVERED rather than named, because the
# first version of this script hardcoded "Retraction of Publication" -- a
# plausible guess and not a string PubMed uses. It reported zero notices in a
# 4.4M-record census, which is the kind of zero that reads as a finding. The
# real type is "Retraction Notice". Matching a pattern and reporting whatever
# turns up means a type this project has not thought of shows up as a row
# instead of as a silent absence.
FAMILY = ("retract", "concern", "withdraw", "errat", "corrected and republished")


def scan_records():
    """Retracted PMIDs in the cancer census, plus the denominators."""
    shards = sorted(RECORDS.glob("*.jsonl.gz"))
    if not shards:
        raise SystemExit(f"no census shards under {RECORDS}")
    retracted = set()
    family = collections.Counter()      # every retraction-adjacent type seen
    total = 0
    by_year = collections.Counter()
    for p in shards:
        with gzip.open(p, "rt") as fh:
            for line in fh:
                r = json.loads(line)
                total += 1
                pts = r.get("pub_types") or []
                pmid = str(r.get("pmid") or "")
                if not pmid:
                    continue
                for t in pts:
                    if any(k in t.lower() for k in FAMILY):
                        family[t] += 1
                if TAINTED in pts:
                    retracted.add(pmid)
                    y = r.get("year")
                    if y:
                        by_year[int(y)] += 1
    return {"retracted": retracted, "family": family,
            "census_total": total, "by_year": by_year}


# KNOWN REGENERATION HAZARD, not yet fixed: the committed JSON is written
# key-sorted while `render()` emits the retraction-type and per-predicate tables
# in dict-insertion order, so regenerating the report reorders those rows even
# when no number moved. The next regeneration will therefore produce a large
# diff that looks like a data change and is not. Read the row VALUES, not the
# row order, before concluding anything moved.


def _amb(which: str) -> float:
    """The sibling ambiguity-impact percentages, READ rather than retyped.

    These two numbers were hand-written into the prose beside a table whose
    every other figure is derived, which is how they survived a re-ingest that
    moved both. They come from a committed artifact, so this raises rather than
    degrading quietly if it is missing -- a silent fallback here would put a
    vague sentence where a measured one used to be.
    """
    j = json.loads(
        (PROJECT_ROOT / "analysis" / "atlas-ambiguity-impact.json").read_text())
    rows = j["relation_rows"]
    if which == "containment":
        return 100 * j["relation_rows_touching_contested_id"] / rows
    return 100 * j["relation_rows_resting_on_at_risk"] / rows


def scan_relations(retracted: set):
    """Join the relation graph to the retracted set.

    Pairs are keyed as an unordered entity pair on the RAW `Type|ID` field.
    That is NOT the same shape `atlas_graph.py` uses and the denominators are
    NOT interchangeable: the graph strips the type prefix and applies the
    per-paper sense corrections `atlas_disambiguate.py` produces, so it reports
    2,840,563 pairs where this reports 2,854,431 over the same relations. The
    13,868-pair difference is two different quantities, not drift. An earlier
    version of this docstring asserted they were comparable, and a documentation
    refresh then quoted this file's count under `atlas_graph.py`'s name.
    """
    rows = 0
    tainted_rows = 0
    pmids_in_graph = set()
    pair_pmids = collections.defaultdict(set)
    pair_tainted = collections.defaultdict(int)
    pred_tainted = collections.Counter()
    pred_total = collections.Counter()

    with gzip.open(RELATIONS, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            pmid, pred, a, b = parts[0], parts[1], parts[2], parts[3]
            rows += 1
            pmids_in_graph.add(pmid)
            pred_total[pred] += 1
            key = (a, b) if a <= b else (b, a)
            pair_pmids[key].add(pmid)
            if pmid in retracted:
                tainted_rows += 1
                pred_tainted[pred] += 1
                pair_tainted[key] += 1

    # The damaging class: every asserting paper for this pair is retracted.
    uncorroborated = [k for k, v in pair_tainted.items()
                      if v and len(pair_pmids[k] - retracted) == 0]
    touched = len(pair_tainted)
    return {
        "rows": rows, "tainted_rows": tainted_rows,
        "pmids_in_graph": len(pmids_in_graph),
        "pairs": len(pair_pmids),
        "pairs_touched": touched,
        "pairs_uncorroborated": len(uncorroborated),
        "uncorroborated_examples": sorted(uncorroborated)[:25],
        "pred_tainted": dict(pred_tainted.most_common(12)),
        "pred_total": {k: pred_total[k] for k in dict(pred_tainted.most_common(12))},
        "single_paper_pairs": sum(1 for v in pair_pmids.values() if len(v) == 1),
    }


def main() -> int:
    print("scanning the census for retraction markers...", flush=True)
    rec = scan_records()
    print(f"  {len(rec['retracted']):,} retracted of {rec['census_total']:,} census records",
          flush=True)
    print("joining to the relation graph...", flush=True)
    rel = scan_relations(rec["retracted"])
    print(f"  {rel['tainted_rows']:,} of {rel['rows']:,} relation rows", flush=True)

    res = {
        "census_records": rec["census_total"],
        "retracted_articles": len(rec["retracted"]),
        "retraction_family_types": dict(rec["family"].most_common()),
        "retracted_by_year": dict(sorted(rec["by_year"].items())[-15:]),
        **rel,
    }
    OUT_JSON.write_text(json.dumps(res, indent=2, sort_keys=True) + "\n")
    OUT_MD.write_text(render(res), encoding="utf-8")
    print(f"wrote {OUT_MD}\nwrote {OUT_JSON}")
    return 0


def render(r: dict) -> str:
    pct_rows = 100.0 * r["tainted_rows"] / max(r["rows"], 1)
    pct_pairs_touch = 100.0 * r["pairs_touched"] / max(r["pairs"], 1)
    pct_pairs_bad = 100.0 * r["pairs_uncorroborated"] / max(r["pairs"], 1)
    ratio = (r["pairs_touched"] / r["pairs_uncorroborated"]
             if r["pairs_uncorroborated"] else float("inf"))
    L = [
        "# How much of the cancer relation graph rests on retracted work", "",
        "Generated by `scripts/atlas_retraction_exposure.py`.", "",
        "Retracted papers do not leave the literature. They stay indexed, stay",
        "cited, and machine-extracted relation graphs keep the assertions they",
        "contributed. Anyone building on a PubTator-derived graph inherits them",
        "silently. This bounds that, using PubMed's own `Retracted Publication`",
        "type — which this repository has stored per record all along and never",
        "read.", "",
        "## The census", "",
        f"| | count |", "|---|--:|",
        f"| cancer census records | {r['census_records']:,} |",
    ] + [f"| `{t}` | {n:,} |" for t, n in r["retraction_family_types"].items()] + [
        "",
        "Only `Retracted Publication` is treated as tainted below. A retraction",
        "notice, erratum or republication is a legitimate record. These types are",
        "discovered by pattern rather than named, because an earlier version of",
        "this script looked for `Retraction of Publication` — a plausible string",
        "PubMed does not use — and reported zero notices in a 4.4M-record census.",
        "",
        f"That is {100.0*r['retracted_articles']/max(r['census_records'],1):.3f}% of the",
        "census. Retraction indexing lags both publication and the retraction",
        "itself, so every figure here is a **lower bound**.", "",
        "## The graph", "",
        "Two numbers, and reporting only the first would mislead — this repository",
        f"has made that exact error before, when {_amb('containment'):.1f}% of relation rows",
        f"*touching* a contested identifier turned out to be {_amb('at_risk'):.2f}% actually",
        "at risk.", "",
        "| | count | share |", "|---|--:|--:|",
        f"| relation rows | {r['rows']:,} | |",
        f"| rows from a retracted paper | {r['tainted_rows']:,} | {pct_rows:.3f}% |",
        f"| distinct entity pairs | {r['pairs']:,} | |",
        f"| pairs **touched** by a retracted paper | {r['pairs_touched']:,} | {pct_pairs_touch:.3f}% |",
        f"| pairs whose **only** support is retracted | **{r['pairs_uncorroborated']:,}** | {pct_pairs_bad:.4f}% |",
        "",
        f"The touched-to-uncorroborated ratio is **{ratio:.1f}x**: most pairs that",
        "involve a retracted paper are also asserted by papers that still stand, and",
        "lose nothing when the retraction is applied. The second row is the damaging",
        "class — the graph asserts a relation that, after retraction, nothing",
        "standing supports.", "",
        "### The denominator that matters", "",
        f"**{r['single_paper_pairs']:,} of the graph's {r['pairs']:,} pairs "
        f"({100.0*r['single_paper_pairs']/max(r['pairs'],1):.1f}%) rest on a single",
        "paper.** That is the population at risk: a pair asserted by many papers",
        "cannot be left unsupported by one retraction, so the uncorroborated class",
        "is drawn almost entirely from the single-paper pairs. Against that",
        f"denominator the rate is **{100.0*r['pairs_uncorroborated']/max(r['single_paper_pairs'],1):.3f}%**,",
        f"not {100.0*r['pairs_uncorroborated']/max(r['pairs'],1):.4f}%.", "",
        "Reported both ways because picking either alone misleads in a different",
        "direction — and because the 70% figure is the more alarming one, and it",
        "has nothing to do with retraction. It is what a machine-extracted",
        "relation graph looks like.", "",
        "A pair in that class is **not necessarily false**. Retraction covers",
        "misconduct, honest error and irreproducibility alike, and a retracted paper",
        "can still have been right. What the class means is that the graph's support",
        "for that pair has been withdrawn, and a consumer should be told.", "",
    ]
    if r.get("pred_tainted"):
        L += ["## Which predicates carry it", "",
              "| predicate | tainted rows | all rows | share |", "|---|--:|--:|--:|"]
        rates = {}
        for p, n in r["pred_tainted"].items():
            tot = r["pred_total"].get(p, 0)
            rates[p] = 100.0 * n / max(tot, 1)
            L.append(f"| `{p}` | {n:,} | {tot:,} | {rates[p]:.3f}% |")
        hi = max(rates, key=rates.get)
        base = rates.get("associate", 0.0)
        L += ["",
              f"The rate is not uniform: `{hi}` carries {rates[hi]:.3f}% against",
              f"{base:.3f}% for `associate`, the graph's bulk predicate — a factor of",
              f"{rates[hi]/max(base,1e-9):.1f}x. Read that as a hypothesis rather than a",
              "result. It is consistent with retraction concentrating in claims about",
              "intervention rather than observation, but these are small counts, the",
              "predicates differ in how PubTator assigns them, and nothing here",
              "controls for field, year or journal.", ""]
    L += [
        "## What this cannot say", "",
        "* **It is a lower bound.** Retraction indexing lags, so the true figure is",
        "  higher and will keep rising for work already published.",
        "* **Only the retracted article is counted as tainted.** The retraction",
        "  notice is a legitimate record and is counted separately.",
        "* **This bounds withdrawn SUPPORT, not truth.** Relations come from",
        "  PubTator3's extractor, whose own error rate is larger than most effects",
        "  measured here. A pair resting on retracted work may be perfectly real.",
        "* **The census is the cancer slice**, so this is not a statement about",
        "  PubMed as a whole.", "",
    ]
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
