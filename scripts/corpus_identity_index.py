#!/usr/bin/env python3
"""Every identifier this project already holds, in one place.

WHY THIS EXISTS
---------------
An ingestion that cannot say what it already has will re-download it. This
project holds articles in five places that grew separately and key on different
things: the MeSH census and the text-recovered stream (PMID, DOI, PMC id), the
frozen full-text corpus and the abstract archive (PMID), and 1.1 million
open-access full texts on external storage (PMC id + PMID). Nothing joined
them, so "is this article new?" had no cheap answer.

THE ANSWER HAS TO BE CHEAP, because it is asked once per candidate article and
there are millions of candidates. A set in memory would work and would have to
be rebuilt on every run over ~6 GB of gzipped JSON; SQLite makes it durable, so
a resumed or repeated run pays nothing.

THREE KEYS, NOT ONE
-------------------
The same article arrives under different identifiers from different sources: a
census record has a PMID, a full-text shard has a PMC id, an OpenAlex or
Crossref result has only a DOI, and a preprint has a DOI and no PMID at all.
Indexing one key would let the same paper in through another door. All three
are stored, normalised, and any hit on any of them means "held".

WHAT A HIT DOES AND DOES NOT MEAN
---------------------------------
`held` means an identifier is already known to this project, NOT that its full
text is on disk. The two are different questions and the `has_fulltext` column
keeps them apart: the census knows about 4.4M articles and holds full text for
a fraction of them, so an ingestion looking for TEXT must ask the second
question, while one looking for new RECORDS asks the first.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import sqlite3
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DB = REPO / "corpus" / "atlas" / "identity.sqlite"
NAS_FULLTEXT = Path(
    os.getenv("FERRO_ATLAS_FULLTEXT", str(Path.home() / "nas" / "cancer-atlas" / "fulltext"))
)

_DOI_PREFIX = re.compile(r"^\s*(?:https?://(?:dx\.)?doi\.org/|doi:)\s*", re.I)


def norm_doi(v) -> str | None:
    """Lowercased, prefix-stripped. DOIs are case-insensitive by spec, and the
    same DOI arrives as a bare string from PubMed and as a URL from OpenAlex --
    indexing both forms would report a held article as new."""
    if not v:
        return None
    s = _DOI_PREFIX.sub("", str(v)).strip().rstrip(".").lower()
    # STRUCTURE, not length. The first version required more than six
    # characters, which is a guess at what a DOI looks like and rejects the
    # short-but-valid `10.1/x` shape. A DOI is `10.<registrant>/<suffix>`, so
    # that is what is checked: the prefix, a slash, and something either side
    # of it. A length threshold would keep drifting; the syntax does not.
    if not s.startswith("10.") or "/" not in s:
        return None
    registrant, _, suffix = s.partition("/")
    return s if len(registrant) > 3 and suffix else None


def norm_pmid(v) -> str | None:
    if not v:
        return None
    s = str(v).strip()
    return s if s.isdigit() and s != "0" else None


def norm_pmcid(v) -> str | None:
    """`PMC` + digits, uppercased. Sources emit `PMC123`, `pmc123` and a bare
    `123`; the bare form is REFUSED rather than guessed, because a bare number
    is indistinguishable from a PMID and guessing would collide two namespaces."""
    if not v:
        return None
    s = str(v).strip().upper()
    if s.startswith("PMC") and s[3:].isdigit():
        return s
    return None


SCHEMA = """
CREATE TABLE IF NOT EXISTS held (
    kind          TEXT NOT NULL,      -- pmid | pmcid | doi
    key           TEXT NOT NULL,
    source        TEXT NOT NULL,
    has_fulltext  INTEGER NOT NULL DEFAULT 0,
    PRIMARY KEY (kind, key)
) WITHOUT ROWID;
CREATE INDEX IF NOT EXISTS held_ft ON held(has_fulltext);
CREATE TABLE IF NOT EXISTS scanned (
    path      TEXT PRIMARY KEY,
    mtime     REAL NOT NULL,
    size      INTEGER NOT NULL,
    records   INTEGER NOT NULL,
    finished  REAL NOT NULL
);
"""


def connect() -> sqlite3.Connection:
    DB.parent.mkdir(parents=True, exist_ok=True)
    c = sqlite3.connect(DB, timeout=60)
    c.executescript(SCHEMA)
    c.execute("PRAGMA journal_mode=WAL")
    c.execute("PRAGMA synchronous=NORMAL")
    return c


def _already(c, path: Path) -> bool:
    """A file is re-scanned when its mtime or size moves, and skipped
    otherwise. Keying on the name alone would make a regenerated shard
    invisible; keying on a hash would mean reading every byte to decide whether
    to read every byte."""
    try:
        st = path.stat()
    except OSError:
        return True
    row = c.execute("SELECT mtime, size FROM scanned WHERE path=?", (str(path),)).fetchone()
    return bool(row and abs(row[0] - st.st_mtime) < 1e-6 and row[1] == st.st_size)


def _record(c, path: Path, n: int) -> None:
    st = path.stat()
    c.execute(
        "INSERT OR REPLACE INTO scanned(path, mtime, size, records, finished) VALUES (?,?,?,?,?)",
        (str(path), st.st_mtime, st.st_size, n, time.time()),
    )


def _add(rows, kind, key, source, ft):
    if key:
        rows.append((kind, key, source, 1 if ft else 0))


def _flush(c, rows):
    if not rows:
        return
    # `has_fulltext` is raised but never lowered: learning that a held id also
    # has text is new information, and forgetting it would be a regression.
    c.executemany(
        "INSERT INTO held(kind,key,source,has_fulltext) VALUES (?,?,?,?) "
        "ON CONFLICT(kind,key) DO UPDATE SET has_fulltext=MAX(held.has_fulltext, excluded.has_fulltext)",
        rows,
    )
    rows.clear()


def scan_jsonl_gz(c, path: Path, source: str, fulltext: bool, verbose=True) -> int:
    if _already(c, path):
        return 0
    rows, n = [], 0
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
            for ln in fh:
                try:
                    r = json.loads(ln)
                except Exception:
                    continue
                n += 1
                has_text = fulltext and bool(r.get("text"))
                _add(rows, "pmid", norm_pmid(r.get("pmid")), source, has_text)
                _add(rows, "pmcid", norm_pmcid(r.get("pmcid")), source, has_text)
                _add(rows, "doi", norm_doi(r.get("doi")), source, has_text)
                if len(rows) >= 20000:
                    _flush(c, rows)
    except (OSError, EOFError, gzip.BadGzipFile) as e:
        print(f"  ! unreadable {path.name}: {e}", file=sys.stderr)
        return 0
    _flush(c, rows)
    _record(c, path, n)
    c.commit()
    if verbose:
        print(f"  {source:22s} {path.name[:44]:44s} {n:>8,}")
    return n


def scan_markdown_dir(c, d: Path, source: str, fulltext: bool) -> int:
    """The frozen corpus and abstract archive are one file per PMID, named by
    it, with the DOI in front-matter."""
    if not d.is_dir():
        return 0
    rows, n = [], 0
    for p in d.glob("*.md"):
        n += 1
        _add(rows, "pmid", norm_pmid(p.stem), source, fulltext)
        try:
            head = p.read_text(encoding="utf-8", errors="replace")[:1200]
        except OSError:
            head = ""
        m = re.search(r"^doi:\s*(\S+)", head, re.M)
        if m:
            _add(rows, "doi", norm_doi(m.group(1)), source, fulltext)
        m = re.search(r"^pmcid:\s*(\S+)", head, re.M)
        if m:
            _add(rows, "pmcid", norm_pmcid(m.group(1)), source, fulltext)
        if len(rows) >= 20000:
            _flush(c, rows)
    _flush(c, rows)
    c.commit()
    print(f"  {source:22s} {str(d.relative_to(REPO)):44s} {n:>8,}")
    return n


def build(include_nas=True) -> dict:
    c = connect()
    t0 = time.time()
    total = 0
    for sub, src in (("records", "census-mesh"), ("records_unindexed", "census-text"),
                     ("records_updates", "census-updates")):
        d = REPO / "corpus" / "atlas" / sub
        if d.is_dir():
            for p in sorted(d.glob("*.jsonl.gz")):
                total += scan_jsonl_gz(c, p, src, fulltext=False, verbose=False)
            print(f"  {src:22s} {sub:44s} scanned")
    total += scan_markdown_dir(c, REPO / "corpus" / "by-pmid", "frozen-fulltext", True)
    total += scan_markdown_dir(c, REPO / "corpus" / "abstracts" / "by-pmid", "abstract-archive", False)
    for d in sorted((REPO / "corpus" / "living").glob("*")):
        if d.is_dir():
            total += scan_markdown_dir(c, d, "living-review", False)
    if include_nas:
        shards = NAS_FULLTEXT / "shards"
        if shards.is_dir():
            for p in sorted(shards.glob("*.jsonl.gz")):
                total += scan_jsonl_gz(c, p, "pmc-oa-fulltext", fulltext=True)
        else:
            print(f"  ! full-text store not present at {shards} -- SKIPPED, so "
                  "this index does not yet know what text is held", file=sys.stderr)
    stats = summary(c)
    stats["records_scanned"] = total
    stats["seconds"] = round(time.time() - t0, 1)
    c.close()
    return stats


def summary(c) -> dict:
    q = lambda s, *a: c.execute(s, a).fetchone()[0]
    return {
        "pmid": q("SELECT COUNT(*) FROM held WHERE kind='pmid'"),
        "pmcid": q("SELECT COUNT(*) FROM held WHERE kind='pmcid'"),
        "doi": q("SELECT COUNT(*) FROM held WHERE kind='doi'"),
        "with_fulltext": q("SELECT COUNT(*) FROM held WHERE has_fulltext=1"),
        "files_scanned": q("SELECT COUNT(*) FROM scanned"),
    }


def is_held(c, pmid=None, pmcid=None, doi=None) -> bool:
    """Any hit on any key. The point of three keys is that one source's only
    identifier may be another's missing one."""
    for kind, key in (("pmid", norm_pmid(pmid)), ("pmcid", norm_pmcid(pmcid)),
                      ("doi", norm_doi(doi))):
        if key and c.execute("SELECT 1 FROM held WHERE kind=? AND key=? LIMIT 1",
                             (kind, key)).fetchone():
            return True
    return False


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-nas", action="store_true",
                    help="skip the external full-text store")
    ap.add_argument("--stats", action="store_true", help="report and exit")
    a = ap.parse_args()
    if a.stats:
        c = connect()
        print(json.dumps(summary(c), indent=2))
        c.close()
        return 0
    s = build(include_nas=not a.no_nas)
    print(json.dumps(s, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
