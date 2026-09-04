#!/usr/bin/env python3
"""Fetch cancer literature this project does not already hold.

WHAT IT WILL AND WILL NOT TAKE
------------------------------
Only content served by an API built to serve it. Europe PMC's full-text
endpoint returns text for the open-access subset and refuses the rest, so the
paywall boundary is enforced by the source rather than by this script guessing.
Publisher pages are never fetched.

That line matters more than "is it free to read". This project's own survey
found that of census records with a DOI and no PMC id, ~22% are readable by
some route but only ~6% carry a licence permitting redistribution -- the
difference is mostly BRONZE, free on the publisher's site with no licence at
all. MISSION.md says neither scraping nor redistributing that is intended, and
a bronze article is exactly the case where "free" and "ours to keep" come
apart. Every record stores the licence Europe PMC reports, so what was taken
under what terms stays answerable.

NOTHING IS FETCHED TWICE
------------------------
Every candidate is checked against `corpus_identity_index` on all three of
PMID, PMC id and DOI before anything is downloaded, and every item written is
added to that index in the same transaction. A re-run therefore skips what it
already has, and an interrupted run resumes from its stored cursor rather than
from the beginning.

WHY THE SEARCH IS PER-SOURCE AND PER-YEAR
-----------------------------------------
Europe PMC's cursor paging is capped, so a single 2.2M-hit query cannot be
walked to the end. Splitting by source and publication year keeps every slice
under the cap and makes progress resumable at a granularity that survives
interruption. It also makes the ledger legible: a slice that yielded nothing
is visible as a slice, not lost in an average.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import sqlite3
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from corpus_identity_index import (  # noqa: E402
    connect as id_connect, is_held, norm_doi, norm_pmcid, norm_pmid,
)

REPO = Path(__file__).resolve().parent.parent
OUT_ROOT = Path(os.getenv(
    "FERRO_EXPAND_OUT",
    str(Path.home() / "nas" / "cancer-atlas" / "expanded"),
))
STATE_DB = REPO / "corpus" / "atlas" / "expand_state.sqlite"
SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
CONTACT = os.getenv("FERRO_CONTACT", "cancer-research-corpus")
UA = f"cancer-research-corpus/1.0 (+{CONTACT}) python-urllib"

# Sources measured to carry cancer content, with whether any of it is open
# access. The zero-OA sources are still worth their METADATA -- a patent or a
# thesis is a record this project has never had -- so they are fetched for
# abstracts and skipped for text.
# MEDLINE first, because that is where the FULL TEXT is: 2.1M of its cancer
# records are open access, against ~16k for preprints. Preprints are 99.9% new
# to this project and worth having, but Europe PMC holds text for almost none
# of them, so putting them first buys metadata and defers everything else.
# Every slice is resumable, so the order changes only what arrives soonest.
SOURCES = [
    ("MED", True), ("PMC", True), ("PPR", True),
    ("PAT", False), ("ETH", False), ("AGR", False),
    ("CBA", False), ("HIR", False), ("CTX", False),
]
PAGE = 1000
SLEEP = float(os.getenv("FERRO_EXPAND_SLEEP", "0.34"))
# Small enough that an in-progress run is visible and a hard kill costs
# little, large enough that shard count stays manageable over millions of
# records. A finished shard is a complete gzip file; a buffered one is not.
SHARD_RECORDS = 4000

STATE_SCHEMA = """
CREATE TABLE IF NOT EXISTS slice (
    src TEXT NOT NULL, year INTEGER NOT NULL,
    cursor TEXT NOT NULL DEFAULT '*',
    seen INTEGER NOT NULL DEFAULT 0,
    kept INTEGER NOT NULL DEFAULT 0,
    fulltext INTEGER NOT NULL DEFAULT 0,
    hits INTEGER NOT NULL DEFAULT -1,
    done INTEGER NOT NULL DEFAULT 0,
    updated REAL NOT NULL DEFAULT 0,
    PRIMARY KEY (src, year)
);
"""


def state() -> sqlite3.Connection:
    STATE_DB.parent.mkdir(parents=True, exist_ok=True)
    c = sqlite3.connect(STATE_DB, timeout=60)
    c.executescript(STATE_SCHEMA)
    c.execute("PRAGMA journal_mode=WAL")
    return c


def _get(url: str, tries: int = 4, timeout: int = 120):
    """Retry on transient failure, give up loudly on a persistent one.

    A silent give-up here would look exactly like a slice with no new records,
    which is the one thing this ledger must never confuse.
    """
    last = None
    for a in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            last = e
        except Exception as e:  # noqa: BLE001
            last = e
        time.sleep(1.5 * (a + 1))
    raise RuntimeError(f"giving up after {tries}: {url[:110]} :: {last}")


def search(src: str, year: int, cursor: str):
    q = f"SRC:{src} AND PUB_YEAR:{year} AND cancer"
    url = (f"{SEARCH}?format=json&pageSize={PAGE}&resultType=core"
           f"&cursorMark={urllib.parse.quote(cursor)}&query={urllib.parse.quote(q)}")
    raw = _get(url)
    return json.loads(raw) if raw else {}


OPEN_LICENCES = ("cc0", "cc by", "cc-by", "cc by-nc", "cc by-sa", "cc by-nd",
                 "cc by-nc-sa", "cc by-nc-nd", "public domain")


def _text_available(rec) -> bool:
    """Does Europe PMC actually HOLD the text?

    `inEPMC` says so directly. Preprints are the case that forces the
    distinction: they are overwhelmingly CC-BY, so a permission-only gate says
    yes to all 108,000 of them, and Europe PMC holds the text for barely a
    sixth. The rest answer 404 -- correctly, but only after a request and a
    sleep each, which is most of a day spent being told no.
    """
    return rec.get("inEPMC") == "Y" or rec.get("isOpenAccess") == "Y"


def _text_permitted(rec) -> bool:
    """Are we allowed to take it?

    Separate from whether it is there. `isOpenAccess` means "in Europe PMC's
    open-access subset", NOT "openly licensed" -- about 118,000 cancer records
    are CC-BY while flagged `OPEN_ACCESS:N`, and gating on the flag alone
    silently skipped every one. The endpoint still refuses anything it should
    not serve, so this widens what is ASKED and cannot widen what is taken.
    """
    if rec.get("isOpenAccess") == "Y":
        return True
    lic = (rec.get("license") or "").strip().lower()
    return any(lic.startswith(x) for x in OPEN_LICENCES)


def fetch_fulltext(rec) -> str | None:
    """Europe PMC serves full text only for what it is allowed to serve.

    Asking for a paywalled article returns 404, which is the refusal doing the
    work: this script never has to decide what is free, and never sees the
    inside of an article it should not have.
    """
    pmcid = norm_pmcid(rec.get("pmcid"))
    if pmcid:
        raw = _get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML",
                   tries=2, timeout=120)
        if raw:
            return raw.decode("utf-8", "replace")
    if rec.get("source") == "PPR" and rec.get("id"):
        raw = _get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/PPR/{rec['id']}/fullTextXML",
                   tries=2, timeout=120)
        if raw:
            return raw.decode("utf-8", "replace")
    return None


class Shards:
    """Rotating gzip shards, flushed and closed on every rotation so an
    interrupted run leaves complete files rather than one truncated one."""

    def __init__(self, root: Path, name: str):
        self.dir = root / name
        self.dir.mkdir(parents=True, exist_ok=True)
        self.name, self.n, self.fh = name, 0, None
        self.idx = self._resume()

    def _resume(self) -> int:
        existing = sorted(self.dir.glob(f"{self.name}-*.jsonl.gz"))
        return int(existing[-1].stem.split("-")[-1].split(".")[0]) + 1 if existing else 0

    def write(self, obj) -> None:
        if self.fh is None:
            self.fh = gzip.open(self.dir / f"{self.name}-{self.idx:05d}.jsonl.gz", "wt",
                                encoding="utf-8")
        self.fh.write(json.dumps(obj, ensure_ascii=False) + "\n")
        self.n += 1
        if self.n % SHARD_RECORDS == 0:
            self.close()
            self.idx += 1

    def close(self) -> None:
        if self.fh:
            self.fh.close()
            self.fh = None


def run(years, sources, limit_new=None, verbose=True) -> dict:
    st, ident = state(), id_connect()
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    shards = Shards(OUT_ROOT, "expanded")
    totals = {"seen": 0, "new": 0, "fulltext": 0, "slices": 0}
    try:
        for src, oa_possible in sources:
            for year in years:
                row = st.execute("SELECT cursor, done FROM slice WHERE src=? AND year=?",
                                 (src, year)).fetchone()
                if row and row[1]:
                    continue
                cursor = row[0] if row else "*"
                st.execute("INSERT OR IGNORE INTO slice(src,year) VALUES (?,?)", (src, year))
                seen = kept = ft = 0
                while True:
                    d = search(src, year, cursor)
                    hits = d.get("resultList", {}).get("result", [])
                    nxt = d.get("nextCursorMark")
                    st.execute("UPDATE slice SET hits=? WHERE src=? AND year=? AND hits<0",
                               (d.get("hitCount", 0), src, year))
                    for rec in hits:
                        seen += 1
                        pmid = norm_pmid(rec.get("pmid"))
                        pmcid = norm_pmcid(rec.get("pmcid"))
                        doi = norm_doi(rec.get("doi"))
                        if is_held(ident, pmid, pmcid, doi):
                            continue
                        if not (pmid or pmcid or doi or rec.get("id")):
                            continue  # nothing to dedup on later; refuse it
                        text = None
                        if oa_possible and _text_available(rec) and _text_permitted(rec):
                            text = fetch_fulltext(rec)
                            time.sleep(SLEEP)
                        out = {
                            "source": src, "epmc_id": rec.get("id"),
                            "pmid": pmid, "pmcid": pmcid, "doi": doi,
                            "title": rec.get("title"), "journal": (
                                rec.get("journalInfo", {}) or {}).get("journal", {}).get("title"),
                            "year": rec.get("pubYear"), "abstract": rec.get("abstractText"),
                            "is_open_access": rec.get("isOpenAccess"),
                            "licence": rec.get("license"),
                            "pub_types": rec.get("pubTypeList", {}).get("pubType"),
                            "has_text": bool(text), "text": text,
                            "fetched": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                        }
                        shards.write(out)
                        kept += 1
                        ft += bool(text)
                        # Indexed IMMEDIATELY: a crash after writing but before
                        # indexing would re-fetch this record on the next run.
                        for kind, key in (("pmid", pmid), ("pmcid", pmcid), ("doi", doi)):
                            if key:
                                ident.execute(
                                    "INSERT INTO held(kind,key,source,has_fulltext) VALUES (?,?,?,?) "
                                    "ON CONFLICT(kind,key) DO UPDATE SET "
                                    "has_fulltext=MAX(held.has_fulltext, excluded.has_fulltext)",
                                    (kind, key, f"expand-{src}", 1 if text else 0))
                        if limit_new and totals["new"] + kept >= limit_new:
                            break
                    st.execute("UPDATE slice SET cursor=?, seen=seen+?, kept=kept+?, "
                               "fulltext=fulltext+?, updated=? WHERE src=? AND year=?",
                               (nxt or cursor, len(hits), kept, ft, time.time(), src, year))
                    st.commit()
                    ident.commit()
                    totals["seen"] += len(hits)
                    totals["new"] += kept
                    totals["fulltext"] += ft
                    kept = ft = 0
                    if limit_new and totals["new"] >= limit_new:
                        return _finish(shards, st, ident, totals)
                    if not nxt or nxt == cursor or not hits:
                        break
                    cursor = nxt
                    time.sleep(SLEEP)
                st.execute("UPDATE slice SET done=1 WHERE src=? AND year=?", (src, year))
                st.commit()
                totals["slices"] += 1
                if verbose and seen:
                    print(f"  {src:4} {year}  seen {seen:>6,}  new {totals['new']:>7,}  "
                          f"text {totals['fulltext']:>7,}", flush=True)
    finally:
        pass
    return _finish(shards, st, ident, totals)


def _finish(shards, st, ident, totals):
    shards.close()
    st.commit(); ident.commit(); st.close(); ident.close()
    return totals


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--from-year", type=int, default=1900)
    ap.add_argument("--to-year", type=int, default=2026)
    ap.add_argument("--limit-new", type=int, default=None,
                    help="stop after this many new records (a dry run)")
    ap.add_argument("--sources", default=None, help="comma list, e.g. PPR,PAT")
    ap.add_argument("--status", action="store_true")
    a = ap.parse_args()
    if a.status:
        c = state()
        rows = c.execute("SELECT src, SUM(seen), SUM(kept), SUM(fulltext), SUM(done), COUNT(*) "
                         "FROM slice GROUP BY src ORDER BY SUM(kept) DESC").fetchall()
        print(f"  {'src':5} {'seen':>10} {'new':>9} {'fulltext':>9} {'slices done':>12}")
        for r in rows:
            print(f"  {r[0]:5} {r[1]:>10,} {r[2]:>9,} {r[3]:>9,} {r[4]:>6}/{r[5]:<5}")
        c.close()
        return 0
    srcs = ([(s, dict(SOURCES).get(s, False)) for s in a.sources.split(",")]
            if a.sources else SOURCES)
    years = list(range(a.to_year, a.from_year - 1, -1))
    t = run(years, srcs, limit_new=a.limit_new)
    print(json.dumps(t, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
