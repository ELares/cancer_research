#!/usr/bin/env python3
"""Audit the expansion crawl for duplicates against pre-existing holdings.

The requirement this answers is "every article we download must be new". That
sounds like a one-line lookup, and the obvious version of it is WRONG in a way
that reports perfect success:

    is_held(idx, pmid, pmcid, doi)   # <-- always True, for every record

because `corpus_expand_fetch` indexes each record into the SAME identity index
in the pass that writes it. By the time an audit runs, every record it wants to
check has been inserted by the crawl itself, so the audit re-discovers its own
writes and reports that everything is a duplicate. Run naively on live data
this reported 6,807 duplicates out of 7,472 -- all of them phantom.

The fix is that the crawl tags its own rows `expand-<SRC>`, so provenance can
be excluded. `held_by_others` asks only about sources that predate the crawl.

Both directions are checked, because a predicate that always returns False
would also report zero duplicates:

  forward : no written record may be held by a pre-existing source
  control : keys known to be held MUST report held

Usage:  corpus_expand_verify.py [--shards DIR] [--index PATH] [--sample N]
"""
from __future__ import annotations

import argparse
import glob
import gzip
import json
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from corpus_identity_index import norm_doi, norm_pmcid, norm_pmid  # noqa: E402

DEFAULT_SHARDS = Path.home() / "nas/cancer-atlas/expanded/expanded"
DEFAULT_INDEX = Path(__file__).resolve().parent.parent / "corpus/atlas/identity.sqlite"

# Rows written by the crawl itself. Excluded from the "did we already have it?"
# question -- they ARE the answer being tested, not evidence about it.
OWN_SOURCE_PREFIX = "expand-"


def held_by_others(idx: sqlite3.Connection, pmid, pmcid, doi) -> bool:
    """Is this article held by something that predates this crawl?

    Any one of the three keys hitting is enough: the same article arrives under
    different identifiers from different sources, so a PMID miss says nothing
    about the DOI.
    """
    for kind, key in (("pmid", norm_pmid(pmid)),
                      ("pmcid", norm_pmcid(pmcid)),
                      ("doi", norm_doi(doi))):
        if key is None:
            continue
        if idx.execute(
            "SELECT 1 FROM held WHERE kind=? AND key=? AND source NOT LIKE ?",
            (kind, key, OWN_SOURCE_PREFIX + "%"),
        ).fetchone():
            return True
    return False


def iter_records(shard_dir: Path):
    """Yield records from every shard, tolerating a shard still being written.

    An in-progress gzip raises EOFError partway through; the records read
    before that point are still valid and still worth auditing.
    """
    for path in sorted(glob.glob(str(shard_dir / "*.jsonl.gz"))):
        try:
            with gzip.open(path, "rt") as fh:
                for line in fh:
                    try:
                        yield json.loads(line)
                    except Exception:
                        continue  # a torn final line, not a reason to stop
        except (EOFError, OSError):
            continue


def verify(shard_dir: Path, index: Path, sample: int = 500) -> dict:
    idx = sqlite3.connect(f"file:{index}?mode=ro", uri=True)
    try:
        written = 0
        dupes: list[str] = []
        for rec in iter_records(shard_dir):
            written += 1
            if held_by_others(idx, rec.get("pmid"), rec.get("pmcid"), rec.get("doi")):
                dupes.append(str(rec.get("pmid") or rec.get("doi") or rec.get("pmcid")))

        # Control. Without this, a predicate broken to always-False would
        # produce a flawless-looking zero above.
        ctl_keys = [r[0] for r in idx.execute(
            "SELECT key FROM held WHERE kind='pmid' AND source NOT LIKE ? LIMIT ?",
            (OWN_SOURCE_PREFIX + "%", sample))]
        ctl_hits = sum(1 for k in ctl_keys if held_by_others(idx, k, None, None))

        return {
            "written": written,
            "duplicates": len(dupes),
            "duplicate_examples": dupes[:10],
            "control_checked": len(ctl_keys),
            "control_hits": ctl_hits,
            "control_ok": bool(ctl_keys) and ctl_hits == len(ctl_keys),
        }
    finally:
        idx.close()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--shards", type=Path, default=DEFAULT_SHARDS)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--sample", type=int, default=500)
    a = ap.parse_args()

    if not a.index.exists():
        print(f"no identity index at {a.index}; build it with corpus_identity_index.py")
        return 2

    r = verify(a.shards, a.index, a.sample)
    pct = 100 * r["duplicates"] / r["written"] if r["written"] else 0.0
    print(f"records written : {r['written']:,}")
    print(f"true duplicates : {r['duplicates']:,}  ({pct:.3f}%)")
    if r["duplicate_examples"]:
        print(f"  e.g. {r['duplicate_examples']}")
    print(f"control         : {r['control_hits']}/{r['control_checked']} known-held report HELD")

    if not r["control_ok"]:
        print("FAIL: control did not hold -- the duplicate check cannot be trusted")
        return 1
    if r["duplicates"]:
        print("FAIL: the crawl wrote articles already held")
        return 1
    print("OK: every written record is new to this project")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
