#!/usr/bin/env python3
"""Fetch cancer literature OpenAlex indexes and Europe PMC does not.

WHY, AND THE MEASUREMENT THAT DECIDED IT
----------------------------------------
Europe PMC covers PubMed, PMC, preprint servers and a handful of specialist
sources. OpenAlex indexes 5,709,238 works whose TITLE OR ABSTRACT is about
cancer, and 2,894,767 of them carry no PMID at all -- so they are outside
PubMed, and therefore outside the crawler that walks Europe PMC.

The first version of this measurement used OpenAlex's `default.search`, which
reads FULL TEXT, and reported 4.9M PMID-less works. Sampling those found hop
powdery mildew, AIDS activism, and a paper titled "PEOPLE": a full-text search
matches anything that mentions cancer once, so most of that pile was noise
rather than missed oncology. Restricting to `title_and_abstract.search` is what
makes this defensible, and a sample of THAT set is genuinely oncology --
tamoxifen 1978, tarextumab NOTCH2/3, nivolumab+ipilimumab -- alongside
conference abstracts, dissertations and book chapters.

CONFERENCE ABSTRACTS ARE THE POINT. In oncology an ASCO or AACR abstract often
reports a trial result years before any paper, and PubMed does not index them.
A census assembled only from indexed journals cannot see that literature exists,
which is a different and worse problem than not being able to read it.

WHAT THIS COLLECTS
------------------
Metadata and abstracts. OpenAlex stores an abstract as an inverted index --
word -> positions -- which reconstructs exactly, and it is CC0, so there is no
licence question about keeping it. It holds no full text to give, and this
script does not go looking for any elsewhere.

NOTHING IS FETCHED TWICE, and the durability contract is IMPORTED from
`corpus_expand_fetch` rather than reimplemented, for the same reason the trials
fetcher imports it: every one of those properties was got wrong at least once
and fixed under review, and a second copy would drift.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import urllib.parse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from corpus_expand_fetch import (  # noqa: E402
    _RateLimit, Shards, TransientFetchError, _get, norm_doi, norm_pmid,
)
from corpus_identity_index import connect as id_connect, is_held  # noqa: E402

API = "https://api.openalex.org/works"
OUT_ROOT = Path(os.getenv(
    "FERRO_OPENALEX_OUT", str(Path.home() / "nas" / "cancer-atlas" / "openalex")))
PAGE = 200                      # the API's documented maximum
SLEEP = float(os.getenv("FERRO_OPENALEX_SLEEP", "0.15"))
SHARD_RECORDS = 4000
# OpenAlex asks for a contact address in the User-Agent and gives the polite
# pool in return. Sending one is the price of the faster queue, not a trick.
MAILTO = os.getenv("FERRO_CONTACT_EMAIL", "research@example.org")

# The subject filter. TITLE AND ABSTRACT, not `default.search`: the latter
# reads full text and pulls in anything that mentions cancer once.
SUBJECT = "title_and_abstract.search:cancer"
# OpenAlex ids are their own namespace -- they collide with no PMID, PMC id or
# DOI -- so a work with none of those is still dedupable.
KIND = "openalex"


def abstract_text(inv: dict | None) -> str | None:
    """Reconstruct an abstract from OpenAlex's inverted index.

    Stored as word -> [positions]. Rebuilding is exact when positions are
    contiguous, which they are for an abstract OpenAlex derived itself; a gap
    means a word was dropped upstream, and the surrounding text is still worth
    more than nothing.
    """
    if not inv:
        return None
    pos = {}
    for word, idxs in inv.items():
        for i in idxs:
            pos[i] = word
    return " ".join(pos[i] for i in sorted(pos)) or None


def _page(cursor: str, extra: str = ""):
    q = {"filter": SUBJECT + extra, "per-page": PAGE, "cursor": cursor,
         "mailto": MAILTO}
    raw = _get(f"{API}?{urllib.parse.urlencode(q)}", tries=4, timeout=120)
    if raw is None:
        # A refusal on a PAGED search is not "no more results"; treating it as
        # one truncates the crawl silently.
        raise TransientFetchError(f"OpenAlex refused at cursor={cursor[:24]!r}")
    return json.loads(raw)


def _record(w: dict) -> dict | None:
    ids = w.get("ids") or {}
    oid = (ids.get("openalex") or w.get("id") or "").rsplit("/", 1)[-1]
    if not oid:
        return None
    loc = w.get("primary_location") or {}
    src = loc.get("source") or {}
    return {
        "source": "OPENALEX",
        "openalex_id": oid,
        "doi": norm_doi(w.get("doi")),
        "pmid": norm_pmid((ids.get("pmid") or "").rsplit("/", 1)[-1] or None),
        "title": w.get("title"),
        "year": w.get("publication_year"),
        "type": w.get("type"),
        "language": w.get("language"),
        "venue": src.get("display_name"),
        "is_oa": bool((w.get("open_access") or {}).get("is_oa")),
        "oa_url": (w.get("open_access") or {}).get("oa_url"),
        "cited_by": w.get("cited_by_count"),
        "topics": [t.get("display_name") for t in (w.get("topics") or [])[:5]],
        "abstract": abstract_text(w.get("abstract_inverted_index")),
        "fetched": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }


def run(limit: int | None = None, only_missing_pmid: bool = True,
        verbose: bool = True) -> dict:
    """Walk the subject filter, storing what this project does not already hold.

    `only_missing_pmid` is the default because a work WITH a PMID is reachable
    through Europe PMC, which is the source that can also give us its text.
    Fetching it here would spend requests to arrive at a worse copy.
    """
    ident = id_connect()
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    shards = Shards(OUT_ROOT, "openalex")
    limiter = _RateLimit(SLEEP)
    totals = {"seen": 0, "new": 0, "with_abstract": 0, "held_already": 0}
    extra = ",has_pmid:false" if only_missing_pmid else ""
    cursor, total = "*", None
    try:
        while cursor:
            limiter.wait()
            d = _page(cursor, extra)
            if total is None:
                total = (d.get("meta") or {}).get("count")
                if verbose:
                    print(f"  {total:,} works match {SUBJECT}{extra}", flush=True)
            results = d.get("results") or []
            if not results:
                break
            for w in results:
                totals["seen"] += 1
                rec = _record(w)
                if rec is None:
                    continue
                if is_held(ident, rec["pmid"], None, rec["doi"]):
                    totals["held_already"] += 1
                    continue
                if ident.execute("SELECT 1 FROM held WHERE kind=? AND key=?",
                                 (KIND, rec["openalex_id"])).fetchone():
                    totals["held_already"] += 1
                    continue
                # WRITE THEN INDEX, the ordering that keeps the index from
                # vouching for bytes that were never stored.
                shards.write(rec)
                for kind, key in (("pmid", rec["pmid"]), ("doi", rec["doi"]),
                                  (KIND, rec["openalex_id"])):
                    if key:
                        ident.execute(
                            "INSERT INTO held(kind,key,source,has_fulltext) "
                            "VALUES (?,?,?,0) ON CONFLICT(kind,key) DO NOTHING",
                            (kind, key, "openalex"))
                totals["new"] += 1
                totals["with_abstract"] += bool(rec["abstract"])
                if limit and totals["new"] >= limit:
                    ident.commit()
                    return totals
            ident.commit()
            if verbose and totals["seen"] % 5000 < PAGE:
                print(f"  seen {totals['seen']:,}  new {totals['new']:,}  "
                      f"already held {totals['held_already']:,}", flush=True)
            cursor = (d.get("meta") or {}).get("next_cursor")
    finally:
        shards.close()
        try:
            ident.commit()
            ident.close()
        except Exception:  # noqa: BLE001
            pass
    return totals


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--include-pmid", action="store_true",
                    help="also walk works Europe PMC can already reach")
    a = ap.parse_args()
    print(json.dumps(run(limit=a.limit, only_missing_pmid=not a.include_pmid),
                     indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
