#!/usr/bin/env python3
"""A full-text index over the corpus, because holding a paper is not finding it.

WHY THIS EXISTS
---------------
This project holds 2.5 million cancer records and had no way to search them.
Every question was answered by asking Europe PMC again -- over the network,
under a rate limit, against THEIR index -- which means the corpus on disk was
costing storage and returning nothing.

The absence was measured, not assumed. Asked for the mechanism "engineered
bacteria tag a tumour so T cells recognise it", a broad full-text conjunction
ranked by citations returned, in order: global mortality statistics for 235
causes of death, a review of exosome biology, and "Cancer genome landscapes".
None is on the subject. The same subject asked as a FIELD-RESTRICTED PHRASE
query returned 8 records with the correct paper first (Probiotic-guided CAR-T
cells for solid tumor targeting, Science 2023, PMID 37824640).

18,766 records versus 8, and the difference was not the corpus. It was the
query. A broad conjunction over full text matches any paper that happens to
contain the words anywhere, and ranking that pool by citation count sorts by
FAME, so the winners are famous papers that mention the terms in passing. The
corpus was never the bottleneck; retrieval was.

WHAT THIS INDEXES AND WHY ONLY THAT
-----------------------------------
Title and abstract, not full text. This is the lever the measurement above
identified: the phrase has to appear where the paper STATES ITS SUBJECT, and a
full-text match is precisely the signal that failed. Full text is also 60 GB
against 4 GB of titles and abstracts, so indexing it would cost fifteen times
as much to reproduce the behaviour we are trying to get away from. If a later
measurement shows full text adds recall the abstracts cannot, it can be added
as a second table; it is not added on the assumption that more text is better.

DEDUPLICATION IS NOT OPTIONAL HERE
----------------------------------
The same paper arrives in an `expanded` shard and again in an `openalex` one.
Indexed twice it occupies two of the ten result slots a human will read, so a
duplicate does not merely waste a row -- it displaces a different paper from
the answer. Keyed on PMID, then PMC id, then DOI, using the same
normalisation as the identity index so the two agree about what one paper is.
"""
from __future__ import annotations

import argparse
import glob
import gzip
import html
import json
import os
import re
import sqlite3
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
from corpus_identity_index import norm_doi, norm_pmcid, norm_pmid  # noqa: E402

DB = REPO / "corpus" / "atlas" / "retrieval.sqlite"
ROOT = Path(os.getenv("FERRO_ATLAS_ROOT", str(Path.home() / "nas" / "cancer-atlas")))
CENSUS = REPO / "corpus" / "atlas"

# WHICH CORPUS, AND THE MISTAKE THAT MADE THIS EXPLICIT.
#
# The first version indexed `expanded` and `openalex` only -- the expansion
# crawl, 2.5 million records, the newest and most visible body of work. The
# evaluation then scored 0.000 on BOTH arms, and the reason was not retrieval:
# not one of the ten known-relevant papers was in the index. All ten are
# credited to `census-mesh` in the identity index. The crawl does not contain
# them.
#
# So the useful corpus for a mechanism question is the CENSUS -- the PubMed
# MeSH baseline that was already on local disk -- and indexing the crawl alone
# measured the wrong 2.5 million records. Both are indexed now, and the
# evaluation is what establishes which one actually answers a question.
#
# `fulltext` shards stay out: they hold {pmcid, pmid, text, licence} and no
# title, so there is nothing here for them to contribute. Listing them would
# scan 16 GB to store nothing.
SOURCES = (
    CENSUS / "records",             # census-mesh, 4.4M MeSH-indexed
    CENSUS / "records_unindexed",   # census-text, 783k text-recovered
    ROOT / "expanded",              # the expansion crawl
    ROOT / "openalex",
)

# The row shape `to_row` returns, named so nothing downstream counts
# positions by hand. Inserting `mesh` at index 2 silently moved `cited_by`
# from 8 to 9, and a test asserting `to_row(...)[8]` then checked a different
# field while still passing.
COLUMNS = ("title", "abstract", "mesh", "pmid", "pmcid", "doi", "journal",
           "year", "has_text", "cited_by")

BATCH = 2000          # rows per transaction
PROGRESS_EVERY = 25   # shards between progress lines

SCHEMA = """
CREATE VIRTUAL TABLE IF NOT EXISTS doc USING fts5(
    title,
    abstract,
    pmid UNINDEXED,
    pmcid UNINDEXED,
    doi UNINDEXED,
    mesh,
    journal UNINDEXED,
    year UNINDEXED,
    has_text UNINDEXED,
    cited_by UNINDEXED,
    tokenize = 'porter unicode61 remove_diacritics 2'
);
-- One row per paper already in `doc`, so a second source cannot index it
-- again and spend a result slot on a duplicate.
CREATE TABLE IF NOT EXISTS seen (
    key TEXT PRIMARY KEY
);
-- Keyed on mtime AND size, not name: a shard the crawler reopened and
-- extended must be re-read, and an untouched one must not be.
CREATE TABLE IF NOT EXISTS scanned (
    path TEXT PRIMARY KEY,
    mtime REAL NOT NULL,
    size INTEGER NOT NULL,
    records INTEGER NOT NULL DEFAULT 0,
    bad INTEGER NOT NULL DEFAULT 0
);
CREATE TABLE IF NOT EXISTS meta (k TEXT PRIMARY KEY, v TEXT);
"""


def connect(db: Path | None = None) -> sqlite3.Connection:
    """Opens DB unless told otherwise.

    `db=None` and a module-global lookup, NOT `db: Path = DB`: a default
    argument is evaluated once at definition time, so reassigning the module
    global had no effect and a smoke test aimed at a scratch file silently
    wrote 18,000 rows into the real index.
    """
    db = Path(db) if db is not None else DB
    db.parent.mkdir(parents=True, exist_ok=True)
    c = sqlite3.connect(str(db))
    c.executescript(SCHEMA)
    # WAL survives a crash mid-build without losing the whole transaction;
    # the build is hours long and must be resumable.
    c.execute("PRAGMA journal_mode=WAL")
    c.execute("PRAGMA synchronous=NORMAL")
    return c


def doc_key(rec: dict) -> str | None:
    """The identity of one paper, in the order the identifiers are trustworthy.

    PMID first because it is a single integer namespace with no spelling
    variants; DOI last because it has the most. Returns None for a record with
    no usable identifier -- such a record is indexed WITHOUT a dedup key rather
    than dropped, since an unidentified paper is still a findable paper.
    """
    p = norm_pmid(rec.get("pmid"))
    if p:
        return f"pmid:{p}"
    m = norm_pmcid(rec.get("pmcid"))
    if m:
        return f"pmcid:{m}"
    d = norm_doi(rec.get("doi"))
    if d:
        return f"doi:{d}"
    return None


_TAG = re.compile(r"<[^>]{1,40}>")
_WS = re.compile(r"\s+")


def clean_text(s: str) -> str:
    """Strip publisher markup so phrases match across it.

    MEASURED, not tidied: 5.9% of indexed titles and 53.7% of abstracts carry
    HTML, and it clusters on exactly the tokens a cancer query is made of --
    gene symbols, which publishers italicise by convention.

    The harm is specific to the tokeniser. FTS5's unicode61 treats `<`, `>`
    and `/` as separators, so `<i>KRAS</i> -mutant` becomes the token stream
    `i · kras · i · mutant`. A phrase search for "KRAS mutant" requires the two
    tokens to be ADJACENT, and the stray `i` from the markup sits between
    them, so the phrase cannot match a title that plainly contains it.

    Entities are unescaped before tags are stripped because the markup arrives
    both ways -- `<i>` from one source and `&lt;i&gt;` from another -- and
    twice, because some records are doubly escaped. Whitespace is collapsed for
    legibility only; it has no effect on matching, since the tokeniser already
    treats any run of separators as one boundary.
    """
    if not s:
        return ""
    for _ in range(2):
        u = html.unescape(s)
        if u == s:
            break
        s = u
    return _WS.sub(" ", _TAG.sub(" ", s)).strip()


def to_row(rec: dict) -> tuple | None:
    """Flatten one source record into the index's columns.

    Returns None when the record has neither a title nor an abstract: there is
    nothing to match on, so indexing it could only ever add noise.
    """
    title = clean_text(rec.get("title") or "")
    abstract = clean_text(rec.get("abstract") or "")
    if not title and not abstract:
        return None
    # MeSH descriptors, as their own searchable column and NEVER merged into
    # the title or abstract. They are NLM's expert subject labels, so they are
    # the strongest retrieval signal here -- and that is exactly why they are
    # kept apart: `mechanism_recall.py` exists because folding MeSH into the
    # matched text makes any keyword measurement circular, the keyword being
    # scored against a label that already encodes the answer. Separate columns
    # let a query use MeSH deliberately and let a measurement exclude it.
    mesh = rec.get("mesh") or []
    if isinstance(mesh, list):
        mesh_s = clean_text("; ".join(str(m) for m in mesh))
    else:
        mesh_s = clean_text(str(mesh))

    return (
        title,
        abstract,
        mesh_s,
        norm_pmid(rec.get("pmid")) or "",
        norm_pmcid(rec.get("pmcid")) or "",
        norm_doi(rec.get("doi")) or "",
        (rec.get("journal") or rec.get("venue") or "").strip(),
        str(rec.get("year") or ""),
        1 if rec.get("has_text") else 0,
        # OpenAlex carries a citation count; Europe PMC records here do not.
        # Stored as -1 rather than 0 so "no data" cannot be read as "never
        # cited" by anything that later ranks on it.
        int(rec["cited_by"]) if isinstance(rec.get("cited_by"), int) else -1,
    )


def shard_paths(root: Path | None = None, sources=SOURCES) -> list[Path]:
    """Every shard across every source.

    `sources` are ABSOLUTE directories, because they no longer share a root:
    the census sits on local disk inside the repo and the crawl sits on
    external storage. `root`, when given, relocates only the sources that live
    under the default external root, so a test or a drive-relative run still
    works.
    """
    out: list[Path] = []
    for d in sources:
        d = Path(d)
        if root is not None and str(d).startswith(str(ROOT)):
            d = Path(root) / d.relative_to(ROOT)
        out.extend(Path(p) for p in sorted(
            glob.glob(str(d / "**" / "*.jsonl.gz"), recursive=True)))
    return out


def index_shard(c: sqlite3.Connection, path: Path, seen: set[str],
                verbose=True) -> tuple[int, int]:
    """Index one shard. Returns (records added, unparseable lines).

    `seen` is the live set of dedup keys and is MUTATED here. It is passed in
    rather than queried per record because the alternative is 2.5 million
    single-row SELECTs; it is a set rather than a list because the first
    version scanned the pending batch linearly, which is 4 million string
    comparisons per 2000-row batch.

    A shard already scanned at this exact mtime and size is skipped and
    returns (0, 0).
    """
    try:
        st = path.stat()
    except OSError:
        return (0, 0)
    row = c.execute("SELECT mtime, size FROM scanned WHERE path=?", (str(path),)).fetchone()
    if row and abs(row[0] - st.st_mtime) < 0.001 and row[1] == st.st_size:
        return (0, 0)

    added = bad = 0
    batch: list[tuple] = []
    keys: list[tuple[str]] = []
    try:
        with gzip.open(path, "rt", errors="replace") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    rec = json.loads(line)
                except Exception:
                    # One corrupt line costs one record, not the shard. Counted
                    # so the loss is visible without being re-suffered.
                    bad += 1
                    continue
                key = doc_key(rec)
                # Covers both the already-indexed case and a duplicate within
                # this same batch, because an accepted key joins `seen`
                # immediately rather than at flush time.
                if key is not None and key in seen:
                    continue
                row = to_row(rec)
                if row is None:
                    continue
                # CLAIMED ONLY ONCE THE ROW IS ACCEPTED. The first version
                # added the key the moment it was read, one line ABOVE this
                # rejection, so a record with neither title nor abstract --
                # an ID-only stub, which the crawl emits -- burned the key for
                # the real paper. The stub was dropped, the key stayed in
                # `seen`, and the same paper arriving later WITH text was
                # skipped as a duplicate. It also could not heal: the stub's
                # key never reached the `seen` TABLE, so a resumed run would
                # accept the record, but the shard is marked scanned and never
                # re-read. Order matters because `seen` means "indexed", and
                # a dropped record was not indexed.
                if key is not None:
                    seen.add(key)
                    keys.append((key,))
                batch.append(row)
                if len(batch) >= BATCH:
                    _flush(c, batch, keys)
                    added += len(batch)
                    batch, keys = [], []
        if batch:
            _flush(c, batch, keys)
            added += len(batch)
    except (OSError, EOFError, gzip.BadGzipFile) as e:
        # A truncated shard is the crawler's live write target, so this path is
        # ORDINARY, not exceptional -- which is why losing data here was a real
        # loss and not a corner case. The pending batch MUST be flushed: the
        # first version dropped it, and shard 00000 holds 392 records, fewer
        # than one BATCH, so all 392 were parsed and then discarded.
        if batch:
            _flush(c, batch, keys)
            added += len(batch)
        if verbose:
            print(f"  ~ {path.name}: no gzip trailer, kept all {added} readable records")
        # MARKED SCANNED, at the size and mtime just observed. The first
        # version deliberately refused to mark it "so the rest is picked up
        # later", which was wrong twice over: `scanned` is keyed on mtime AND
        # size, so a shard the crawler later extends is re-read anyway, and
        # refusing to mark it guarantees a PERMANENT rescan of a 49 MB file
        # over the network on every single run.
        #
        # EOFError here does not mean data was lost. It means the stream has no
        # trailer, which is what Z_SYNC_FLUSH leaves behind when a writer is
        # killed before close(). Measured on expanded-00019: Python's reader
        # returns 3,518 complete lines and then raises; `gzcat` recovers 3,517.
        # The reader is not losing records -- it is reporting a missing CRC.
        c.execute(
            "INSERT OR REPLACE INTO scanned(path, mtime, size, records, bad) "
            "VALUES (?,?,?,?,?)", (str(path), st.st_mtime, st.st_size, added, bad))
        c.commit()
        return (added, bad)

    c.execute(
        "INSERT OR REPLACE INTO scanned(path, mtime, size, records, bad) VALUES (?,?,?,?,?)",
        (str(path), st.st_mtime, st.st_size, added, bad))
    c.commit()
    return (added, bad)


def _flush(c: sqlite3.Connection, batch: list[tuple], keys: list[tuple[str]]) -> None:
    c.executemany(
        "INSERT INTO doc(title, abstract, mesh, pmid, pmcid, doi, journal, "
        "year, has_text, cited_by) VALUES (?,?,?,?,?,?,?,?,?,?)", batch)
    if keys:
        c.executemany("INSERT OR IGNORE INTO seen(key) VALUES (?)", keys)
    c.commit()


def summary(c: sqlite3.Connection) -> dict:
    q = lambda s: c.execute(s).fetchone()[0]  # noqa: E731
    return {
        "docs": q("SELECT COUNT(*) FROM doc"),
        "deduped_keys": q("SELECT COUNT(*) FROM seen"),
        "shards_scanned": q("SELECT COUNT(*) FROM scanned"),
        "bad_lines": q("SELECT COALESCE(SUM(bad),0) FROM scanned"),
        # `has_text=1`, not `has_text='1'`: FTS5 keeps an UNINDEXED column's
        # storage class, and SQLite does not coerce INTEGER to TEXT, so the
        # string form counted zero against 2,161 real rows in the first six
        # shards alone. A summary that cannot go above zero measures nothing.
        #
        # What this counts is the flag the SHARD RECORD carries, which is NOT
        # the corpus's full-text total. Joined against identity.sqlite over the
        # whole 7,727,217-record index: the flag is never 1 where the ledger
        # says 0, but it is 0 for 1,102,481 records the ledger knows have text,
        # so it finds 1,524,514 of 2,626,995 -- 58%. The census writes full-text
        # availability to the ledger rather than back into the record it already
        # wrote, so a one-directional flag is the most this field can support.
        # Read it as a floor on full text, never as a measurement of it.
        "with_fulltext": q("SELECT COUNT(*) FROM doc WHERE has_text=1"),
        "with_citations": q("SELECT COUNT(*) FROM doc WHERE cited_by > -1"),
        "with_mesh": q("SELECT COUNT(*) FROM doc WHERE mesh != ''"),
    }


def build(root: Path = ROOT, limit: int | None = None, verbose=True) -> dict:
    c = connect()
    paths = shard_paths(root)
    if not paths:
        print(f"no shards under {root} -- is the corpus mounted?", file=sys.stderr)
        return summary(c)
    if limit:
        paths = paths[:limit]
    seen: set[str] = {r[0] for r in c.execute("SELECT key FROM seen")}
    if verbose:
        print(f"  resuming with {len(seen):,} known keys")
    t0 = time.time()
    total = badtotal = 0
    for i, p in enumerate(paths, 1):
        a, b = index_shard(c, p, seen, verbose=verbose)
        total += a
        badtotal += b
        if verbose and (i % PROGRESS_EVERY == 0 or i == len(paths)):
            el = time.time() - t0
            rate = total / el if el > 0 else 0
            print(f"  [{i}/{len(paths)}] +{total:,} docs  {rate:,.0f}/s  "
                  f"{el/60:.1f}m  bad={badtotal}", flush=True)
    c.execute("INSERT OR REPLACE INTO meta(k,v) VALUES ('built', ?)",
              (time.strftime("%Y-%m-%dT%H:%M:%S"),))
    c.commit()
    s = summary(c)
    c.close()
    return s


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--root", default=str(ROOT), help="corpus root holding the source dirs")
    ap.add_argument("--limit", type=int, help="index only the first N shards (smoke test)")
    ap.add_argument("--summary", action="store_true", help="report the index and exit")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()

    if a.summary:
        c = connect()
        for k, v in summary(c).items():
            print(f"{k:18} {v:,}")
        c.close()
        return 0

    s = build(Path(a.root), limit=a.limit, verbose=not a.quiet)
    print()
    for k, v in s.items():
        print(f"{k:18} {v:,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
