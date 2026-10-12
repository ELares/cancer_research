#!/usr/bin/env python3
"""Two ways to ask the corpus a question, so the difference can be measured.

The project's working assumption was that more articles would produce better
answers. It did not, and the reason is visible in one comparison: asked for
"engineered bacteria tag a tumour so T cells recognise it", a broad
conjunction ranked by citation count returned global mortality statistics, an
exosome review, and "Cancer genome landscapes". The papers that answer the
question were already held -- 10 of 11 known-relevant papers are in the corpus,
6 with full text. They were never retrievable.

So this module implements the naive method and the field-aware one side by
side, under one interface, specifically so an evaluation can run both over the
same index and report which finds a paper we know is there.

ARM A -- BAG OF WORDS
    Every token OR'd together, both columns weighted equally. This is what
    "search the corpus for X" means if you do not think about it, and it is
    the local analogue of the broad conjunction that failed. A paper matching
    one common word competes with a paper matching the whole question.

ARM B -- PHRASE-AWARE, TITLE-WEIGHTED
    Multi-word concepts stay phrases and are AND'd, and the title is weighted
    an order of magnitude above the abstract. The hypothesis being tested is
    that a paper ABOUT a mechanism names it in its title, while a paper that
    merely mentions it does so once in an abstract.

NEITHER ARM RANKS BY CITATIONS, deliberately. Citation rank is what produced
the mortality tables: it sorts by fame, and fame is the thing a discovery
query is trying to see past. Europe PMC records in this corpus carry no
citation count anyway, so a popularity arm could only be run over the 18,938
records that carry one -- 0.2% of the index, too thin a slice to compare
against the rest. The index size is deliberately not quoted here: it is
reported by the generated artifacts, and a figure hand-written into a
docstring goes stale the next time the corpus grows.
"""
from __future__ import annotations

import argparse
import re
import sqlite3
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
from corpus_retrieval_index import DB, connect  # noqa: E402

# FTS5 reserves these as query syntax; a bare one in user text is a syntax
# error, not a search. Stripped rather than escaped because none of them
# carries meaning in a mechanism question.
_PUNCT = re.compile(r'[^\w\s-]')

# Weights are POSITIONAL over the FTS5 columns, in the order
# `corpus_retrieval_index.COLUMNS` declares them: title, abstract, mesh, then
# the UNINDEXED ones, which cannot match and so cannot contribute whatever
# weight they are given.
#
# MeSH IS WEIGHTED ZERO, AND THAT IS NOT ENOUGH ON ITS OWN. Adding the column
# silently gave it the default weight 1.0 in both arms. Zeroing it removes it
# from the SCORE, but FTS5 matching is separate from ranking, so an
# unrestricted query still MATCHES on MeSH and the pool still grows -- a
# two-row check confirmed it: `"ferroptosis"` hit a row whose only occurrence
# was its MeSH list. Both levers are needed, which is why every query below is
# column-restricted as well.
W_UNIFORM = (1.0, 1.0, 0.0)
W_TITLE_HEAVY = (10.0, 1.0, 0.0)
# For a query that WANTS the expert labels. Never used by a measurement that
# scores a keyword, because MeSH already encodes the answer such a measurement
# is trying to find -- the circularity `mechanism_recall.py` exists for.
W_WITH_MESH = (10.0, 1.0, 5.0)

# The fields a query reads unless told otherwise. Named here so a caller
# cannot reach MeSH by forgetting to restrict.
TEXT_FIELDS = ("title", "abstract")


def _tokens(text: str) -> list[str]:
    return [t for t in _PUNCT.sub(" ", text).split() if t]


def _restrict(body: str, fields=TEXT_FIELDS) -> str:
    """Wrap a query so it can only match the named columns."""
    if not fields:
        return body
    return "{" + " ".join(fields) + "}: (" + body + ")"


def bag_of_words(question: str, fields=TEXT_FIELDS) -> str:
    """Arm A. Every token, OR'd, over the text fields only."""
    toks = _tokens(question)
    if not toks:
        raise ValueError("empty question")
    return _restrict(" OR ".join(f'"{t}"' for t in toks), fields)


def phrase_query(concepts: list[str], fields=TEXT_FIELDS) -> str:
    """Arm B. Each concept as a phrase, all required, text fields only.

    A concept is given as a phrase by the caller rather than guessed from the
    question, because "T cell" and "CAR-T" are two words that mean one thing
    and no tokeniser can know that from the string alone.
    """
    if not concepts:
        raise ValueError("no concepts")
    body = " AND ".join(f'"{c.strip()}"' for c in concepts if c.strip())
    return _restrict(body, fields)


def search(c: sqlite3.Connection, match: str, weights=W_TITLE_HEAVY,
           limit: int = 10) -> list[dict]:
    """Run one FTS5 query and return ranked rows, best first.

    bm25() returns a NEGATIVE score where more negative is a better match, so
    the ordering is ASCENDING. Sorting descending would return the worst
    matches in the pool and look superficially plausible.
    """
    w = ", ".join(str(float(x)) for x in weights)
    sql = (f"SELECT pmid, pmcid, doi, title, journal, year, has_text, "
           f"bm25(doc, {w}) AS score FROM doc WHERE doc MATCH ? "
           f"ORDER BY score ASC LIMIT ?")
    rows = c.execute(sql, (match, limit)).fetchall()
    cols = ("pmid", "pmcid", "doi", "title", "journal", "year", "has_text", "score")
    return [dict(zip(cols, r)) for r in rows]


def count(c: sqlite3.Connection, match: str) -> int:
    """Pool size for a query -- the 18,766-versus-8 number."""
    return c.execute("SELECT COUNT(*) FROM doc WHERE doc MATCH ?", (match,)).fetchone()[0]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("question", help="free text, or --concept phrases for arm B")
    ap.add_argument("--concept", action="append", default=[],
                    help="a phrase that must appear; repeatable (selects arm B)")
    ap.add_argument("--arm", choices=("a", "b"), help="force an arm")
    ap.add_argument("--limit", type=int, default=10)
    ap.add_argument("--db", default=str(DB))
    a = ap.parse_args()

    c = connect(Path(a.db))
    if c.execute("SELECT COUNT(*) FROM doc").fetchone()[0] == 0:
        print("index is empty -- run corpus_retrieval_index.py first", file=sys.stderr)
        return 2

    arm = a.arm or ("b" if a.concept else "a")
    if arm == "b":
        m = phrase_query(a.concept or [a.question])
        w = W_TITLE_HEAVY
    else:
        m = bag_of_words(a.question)
        w = W_UNIFORM

    print(f"arm {arm.upper()}  match: {m}")
    print(f"pool: {count(c, m):,}\n")
    for i, r in enumerate(search(c, m, w, a.limit), 1):
        ft = "+text" if r["has_text"] == 1 else "     "
        print(f"{i:3}. {ft} {r['year']:>4}  pmid={r['pmid'] or '-':>9}  "
              f"{(r['title'] or '')[:72]}")
        print(f"              {(r['journal'] or '')[:60]}   score={r['score']:.2f}")
    c.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
