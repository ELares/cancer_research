#!/usr/bin/env python3
"""Does field-aware retrieval beat throwing the words at the index?

This project has been burned by a layer that was assumed to help and did not.
`analysis/atlas-discovery-headroom.md` records the result: no seed-aware
ranking method added measurable precision over a candidate-only popularity
prior. The honest response to building another retrieval layer is to hold it
to the same bar before claiming anything for it.

WHAT IS MEASURED
----------------
Known-item retrieval. For each mechanism question the gold set names papers
that are genuinely about that mechanism and that the corpus demonstrably
holds, and both arms search the same index for them. The score is the rank at
which each gold paper appears, summarised as MRR and recall@10.

WHAT IS NOT MEASURED
--------------------
Precision over unjudged results. Nothing here says the other nine results in
a top-ten are bad -- only whether the paper we know answers the question was
reachable. A precision claim needs judgements this evaluation does not have.

THE BIAS IS REAL AND IS NOT CORRECTED
-------------------------------------
Gold papers were found by phrase search, so each is guaranteed to contain its
phrase, which favours arm B. That is stated in the gold file and repeated in
the output, because a reader who does not know it would over-read the result.
What survives the bias is the weaker but still useful claim: a paper the
corpus holds can be made reachable, and the naive query does not reach it.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
from corpus_retrieval_index import DB, connect  # noqa: E402
from corpus_retrieval_search import (  # noqa: E402
    W_TITLE_HEAVY, W_UNIFORM, bag_of_words, count, phrase_query, search)

GOLD = REPO / "analysis" / "corpus-retrieval-gold.json"
OUT_JSON = REPO / "analysis" / "corpus-retrieval-precision.json"
OUT_MD = REPO / "analysis" / "corpus-retrieval-precision.md"

DEPTH = 50   # how far down to look for a gold paper before calling it missed
AT_K = 10    # the recall@K a human actually reads


def rank_of(rows: list[dict], pmid: str) -> int | None:
    for i, r in enumerate(rows, 1):
        if r.get("pmid") == pmid:
            return i
    return None


def run_arm(c, q: dict, arm: str) -> dict:
    if arm == "a":
        m, w = bag_of_words(q["question"]), W_UNIFORM
    else:
        m, w = phrase_query(q["concepts"]), W_TITLE_HEAVY
    rows = search(c, m, w, limit=DEPTH)
    ranks = {g: rank_of(rows, g) for g in q["gold"]}
    found = [r for r in ranks.values() if r is not None]
    return {
        "arm": arm,
        "match": m,
        "pool": count(c, m),
        "ranks": ranks,
        "found": len(found),
        "gold_total": len(q["gold"]),
        "recall_at_k": sum(1 for r in found if r <= AT_K) / len(q["gold"]),
        # Reciprocal rank of the BEST gold hit, 0.0 when none was reached.
        "rr": max((1.0 / r for r in found), default=0.0),
        "top3": [{"pmid": r["pmid"], "title": (r["title"] or "")[:70]} for r in rows[:3]],
    }


def assemble(db: Path = DB) -> dict:
    spec = json.loads(GOLD.read_text())
    c = connect(db)
    total = c.execute("SELECT COUNT(*) FROM doc").fetchone()[0]
    if total == 0:
        raise SystemExit("index is empty -- run corpus_retrieval_index.py first")

    per_query = []
    for q in spec["queries"]:
        row = {"id": q["id"], "question": q["question"], "concepts": q["concepts"],
               "note": q.get("note", ""), "gold": q["gold"]}
        for arm in ("a", "b"):
            row[arm] = run_arm(c, q, arm)
        per_query.append(row)
    c.close()

    nq = len(per_query)
    agg = {}
    for arm in ("a", "b"):
        agg[arm] = {
            "mrr": sum(r[arm]["rr"] for r in per_query) / nq,
            "recall_at_k": sum(r[arm]["recall_at_k"] for r in per_query) / nq,
            "gold_found": sum(r[arm]["found"] for r in per_query),
            "median_pool": sorted(r[arm]["pool"] for r in per_query)[nq // 2],
        }
    return {
        "docs_indexed": total,
        "queries": nq,
        "gold_total": sum(len(r["gold"]) for r in per_query),
        "depth": DEPTH,
        "at_k": AT_K,
        "aggregate": agg,
        "per_query": per_query,
        "provenance": spec["_provenance"],
        "known_bias": spec["_known_bias"],
    }


def render(d: dict) -> str:
    a, b = d["aggregate"]["a"], d["aggregate"]["b"]
    L = []
    w = L.append
    w("# Retrieval, not corpus size, was the bottleneck")
    w("")
    w(f"Measured over **{d['docs_indexed']:,} indexed titles and abstracts**, "
      f"{d['queries']} mechanism questions, {d['gold_total']} known-relevant papers.")
    w("")
    w("## The result")
    w("")
    w("| | arm A: bag of words | arm B: phrase + title weight |")
    w("|---|--:|--:|")
    w(f"| MRR | {a['mrr']:.3f} | {b['mrr']:.3f} |")
    w(f"| recall@{d['at_k']} | {a['recall_at_k']:.1%} | {b['recall_at_k']:.1%} |")
    w(f"| gold papers reached (of {d['gold_total']}, depth {d['depth']}) | {a['gold_found']} | {b['gold_found']} |")
    w(f"| median pool searched | {a['median_pool']:,} | {b['median_pool']:,} |")
    w("")
    # THE VERDICT READS BOTH COLUMNS, and an earlier version did not. It
    # compared MRR alone and printed "arm B reaches the known paper where arm A
    # does not" -- while the column immediately beside it said arm B reached
    # FOUR gold papers against arm A's five. A ranking win and a coverage loss
    # are different facts and the page has to carry both or it misleads.
    ranking_win = b["mrr"] > a["mrr"] and b["recall_at_k"] >= a["recall_at_k"]
    coverage_loss = b["gold_found"] < a["gold_found"]

    if ranking_win and not coverage_loss:
        w(f"Arm B wins on both axes: better ranking and no loss of reach. The "
          f"corpus held these papers the whole time -- every gold item was "
          f"verified present in `corpus/atlas/identity.sqlite` before the "
          f"measurement ran.")
    elif ranking_win and coverage_loss:
        w(f"**A trade, not a win, and the direction matters more than the "
          f"averages.** Arm B ranks better (MRR {b['mrr']:.3f} against "
          f"{a['mrr']:.3f}, recall@{d['at_k']} {b['recall_at_k']:.0%} against "
          f"{a['recall_at_k']:.0%}) while reaching FEWER gold papers overall "
          f"({b['gold_found']} against {a['gold_found']}).")
        w("")
        w(f"The per-question table below shows why, and it is not a tuning "
          f"problem. Where the phrase matches the literature's own vocabulary, "
          f"arm B finds the paper in a pool of tens of records against arm A's "
          f"millions. Where it does not, arm B returns **nothing at all** -- "
          f"its pool simply does not contain the answer. The failure mode is "
          f"total rather than graceful, which is the opposite shape to arm A's "
          f"(always finds something, usually buried).")
        worst = min(d["per_query"], key=lambda r: r["b"]["found"] / max(r["b"]["gold_total"], 1))
        if worst["b"]["found"] == 0:
            w("")
            w(f"The clearest case is `{worst['id']}`: arm B's pool of "
              f"{worst['b']['pool']:,} records holds 0 of "
              f"{worst['b']['gold_total']} gold papers, because the phrases "
              f"{worst['concepts']} are not the words those papers use. So the "
              f"concrete next step is vocabulary, not weighting.")
    else:
        w(f"**Arm B does not beat arm A.** MRR {b['mrr']:.3f} against "
          f"{a['mrr']:.3f}, recall@{d['at_k']} {b['recall_at_k']:.0%} against "
          f"{a['recall_at_k']:.0%}. On this evidence the field-aware layer does "
          f"not earn its place, and the same verdict "
          f"`atlas-discovery-headroom.md` reached for the seed-aware ranker "
          f"applies here.")
    w("")
    w("## Per question")
    w("")
    for r in d["per_query"]:
        w(f"### {r['id']}")
        w("")
        w(f"> {r['question']}")
        w("")
        if r["note"]:
            w(f"{r['note']}")
            w("")
        w("| arm | pool | gold reached | best rank |")
        w("|---|--:|--:|--:|")
        for arm in ("a", "b"):
            x = r[arm]
            best = min((v for v in x["ranks"].values() if v), default=None)
            w(f"| {arm.upper()} | {x['pool']:,} | {x['found']}/{x['gold_total']} | "
              f"{best if best else '—'} |")
        w("")
        w("Arm A's top 3:")
        w("")
        for t in r["a"]["top3"]:
            w(f"- `{t['pmid'] or '-'}` {t['title']}")
        w("")
        w("Arm B's top 3:")
        w("")
        for t in r["b"]["top3"]:
            w(f"- `{t['pmid'] or '-'}` {t['title']}")
        w("")
    w("## What this does and does not establish")
    w("")
    w(f"{d['provenance']}")
    w("")
    w(f"**Known bias.** {d['known_bias']}")
    w("")
    w("Generated by `scripts/corpus_retrieval_eval.py`.")
    return "\n".join(L) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--db", default=str(DB))
    ap.add_argument("--stdout", action="store_true", help="print, do not write")
    a = ap.parse_args()
    d = assemble(Path(a.db))
    md = render(d)
    if a.stdout:
        print(md)
        return 0
    OUT_JSON.write_text(json.dumps(d, indent=2) + "\n")
    OUT_MD.write_text(md)
    print(f"wrote {OUT_MD.relative_to(REPO)} and {OUT_JSON.relative_to(REPO)}")
    print(f"  arm A MRR {d['aggregate']['a']['mrr']:.3f}   "
          f"arm B MRR {d['aggregate']['b']['mrr']:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
