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
import os
import sqlite3
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from corpus_identity_index import norm_doi, norm_pmcid, norm_pmid  # noqa: E402

# The SAME env var the fetcher and the reporter read. Hardcoding this path
# meant the one tool in the trio that can DELETE data was the only one that
# could not follow a relocation -- point FERRO_EXPAND_OUT elsewhere and
# --repair would find no shards and call every stored record an orphan.
DEFAULT_SHARDS = Path(os.getenv(
    "FERRO_EXPAND_OUT",
    str(Path.home() / "nas" / "cancer-atlas" / "expanded"))) / "expanded"
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


def iter_records(shard_dir: Path, truncated: list | None = None):
    """Yield records from every shard, tolerating one still being written.

    An in-progress gzip raises EOFError part-way through; the records read
    before that point are valid and worth auditing. A shard that raises is
    RECORDED in `truncated` rather than passed over in silence -- silently
    skipping is how three truncated shards sat unnoticed on disk while the
    identity index claimed every record in their lost tails was held.
    """
    for path in sorted(glob.glob(str(shard_dir / "*.jsonl.gz"))):
        n = 0
        try:
            with gzip.open(path, "rt") as fh:
                for line in fh:
                    try:
                        rec = json.loads(line)
                    except Exception:
                        continue  # a torn final line, not a reason to stop
                    n += 1
                    yield rec
        except (EOFError, OSError) as e:
            if truncated is not None:
                truncated.append((Path(path).name, n, type(e).__name__))


def stored_keys(shard_dir: Path) -> tuple[set, list]:
    """Every identifier actually recoverable from the shards, plus the damage.

    This is the set the identity index's own `expand-*` rows are checked
    against: a row with no matching stored record describes an article whose
    bytes are not there.
    """
    truncated: list = []
    keys = set()
    for rec in iter_records(shard_dir, truncated):
        for kind, key in (("pmid", norm_pmid(rec.get("pmid"))),
                          ("pmcid", norm_pmcid(rec.get("pmcid"))),
                          ("doi", norm_doi(rec.get("doi"))),
                          # The crawl indexes the Europe PMC id as a fourth
                          # key for records carrying none of the other three.
                          # Omitting it here would make every such row look
                          # orphaned, and --repair would delete rows that
                          # describe records sitting on disk.
                          ("epmc", rec.get("epmc_id"))):
            if key:
                keys.add((kind, key))
    return keys, truncated


def live_shard(shard_dir: Path) -> str | None:
    """The shard a RUNNING crawl is writing, or None if nothing is running.

    Asked of the OS rather than inferred from the filenames: with no crawl
    running there is no open shard, and every truncated file is real damage.
    That gate is the whole change -- an earlier docstring claimed the switch
    from `max(...)` to `sorted(...)[-1]` guarded against a second shard prefix,
    which is false twice over: the old code already compared name STRINGS, so
    the semantics are identical, and with two prefixes sorted-last picks the
    LATER prefix, not the newest shard. Returns None when nothing is running.
    """
    # Guarded, and narrower than a bare filename match. `pgrep -f
    # corpus_expand_fetch.py` matches ANY command line containing that string
    # -- an editor, a grep, or `pytest tests/test_corpus_expand_fetch.py` --
    # and a false positive here stamps "open, being written" on a shard that
    # is genuinely truncated, which is the reassurance this label exists to
    # avoid giving. A missing pgrep is worse still: it used to raise
    # FileNotFoundError and take the whole audit down, on an image without
    # procps, for a cosmetic label.
    try:
        out = subprocess.run(
            ["pgrep", "-f", r"python[^ ]* .*scripts/corpus_expand_fetch\.py"],
            capture_output=True, text=True).stdout.strip()
    except (OSError, ValueError):
        return None
    if not out:
        return None
    names = sorted(Path(x).name for x in glob.glob(str(shard_dir / "*.jsonl.gz")))
    return names[-1] if names else None


def orphaned_index_rows(index: Path, shard_dir: Path) -> tuple[list, list]:
    """Index rows this crawl wrote for records that are NOT on disk.

    THE FAILURE THIS REPAIRS. The fetcher used to index a record while its
    bytes sat unflushed in a gzip buffer, so a kill kept the row and lost the
    record -- and because the row says "held", every later run skips that
    article forever. It is silent, permanent, and the exact inverse of what
    the dedup is for. Two kills during development did it to real data.

    Only `expand-*` rows are considered: every other source is backed by
    storage this crawl does not own.
    """
    keys, truncated = stored_keys(shard_dir)
    c = sqlite3.connect(f"file:{index}?mode=ro", uri=True)
    try:
        rows = [(k, key) for k, key, in c.execute(
            "SELECT kind, key FROM held WHERE source LIKE ?",
            (OWN_SOURCE_PREFIX + "%",))]
    finally:
        c.close()
    return [r for r in rows if r not in keys], truncated


def verify(shard_dir: Path, index: Path, sample: int = 500) -> dict:
    idx = sqlite3.connect(f"file:{index}?mode=ro", uri=True)
    try:
        written = 0
        dupes: list[str] = []
        truncated: list = []
        for rec in iter_records(shard_dir, truncated):
            written += 1
            if held_by_others(idx, rec.get("pmid"), rec.get("pmcid"), rec.get("doi")):
                dupes.append(str(rec.get("pmid") or rec.get("doi") or rec.get("pmcid")))

        # Control. Without this, a predicate broken to always-False would
        # produce a flawless-looking zero above.
        ctl_keys = [r[0] for r in idx.execute(
            "SELECT key FROM held WHERE kind='pmid' AND source NOT LIKE ? LIMIT ?",
            (OWN_SOURCE_PREFIX + "%", sample))]
        ctl_hits = sum(1 for k in ctl_keys if held_by_others(idx, k, None, None))

        orphans, _ = orphaned_index_rows(index, shard_dir)
        return {
            "truncated_shards": truncated,
            "live_shard": live_shard(shard_dir),
            "orphaned_index_rows": len(orphans),
            "orphan_examples": orphans[:10],
            "written": written,
            "duplicates": len(dupes),
            "duplicate_examples": dupes[:10],
            "control_checked": len(ctl_keys),
            "control_hits": ctl_hits,
            "control_ok": bool(ctl_keys) and ctl_hits == len(ctl_keys),
        }
    finally:
        idx.close()


# A repair may delete at most this share of the crawl's index rows before it
# refuses and asks to be looked at. Chosen because the failure it guards is
# ALL-OR-NOTHING -- an absent or relocated shard directory orphans 100% of the
# rows, not 30% of them -- so any threshold well under 1.0 catches it, while a
# genuine repair after a kill touches a tiny fraction (374 rows of 40,814, or
# 0.9%, in the incident this tool was written for).
REPAIR_MAX_FRACTION = 0.10


def repair(index: Path, shard_dir: Path, force: bool = False,
           dry_run: bool = False) -> int:
    """Drop index rows describing records that are not on disk.

    THIS FUNCTION DELETES DATA AND ITS INPUT CAN BE SILENTLY EMPTY, which is
    the whole reason for the refusals below. `stored_keys` globs a directory;
    a directory that is missing, unmounted, or simply not the one the crawl
    wrote to yields zero records WITHOUT raising, and every stored row then
    looks like an orphan. The first version had no check at all: pointed at an
    empty directory it deleted the entire `expand-*` index -- 40,814 rows on
    the live index -- printed a cheerful success line and exited 0.

    That is worse than the defect it repairs. A dropped row means the article
    is fetched and written AGAIN, and because the duplicate audit deliberately
    ignores `expand-*` provenance, those duplicates are invisible to it. The
    tool built to fix silent, audit-proof loss would have caused a bigger one.
    """
    orphans, truncated = orphaned_index_rows(index, shard_dir)
    for name, n, err in truncated:
        print(f"  truncated {name}: {n:,} records recovered ({err})")

    stored, _ = stored_keys(shard_dir)
    c = sqlite3.connect(f"file:{index}?mode=ro", uri=True)
    try:
        total = c.execute("SELECT COUNT(*) FROM held WHERE source LIKE ?",
                          (OWN_SOURCE_PREFIX + "%",)).fetchone()[0]
    finally:
        c.close()

    print(f"shard dir       : {shard_dir}")
    print(f"stored keys     : {len(stored):,}")
    print(f"crawl index rows: {total:,}")
    print(f"orphans         : {len(orphans):,}")

    # NOT behind --force. The other two refusals describe a directory that
    # exists and is genuinely empty or genuinely mostly-orphaned, which an
    # operator can legitimately overrule. A path that does not exist cannot be
    # any of those: forcing it deletes the whole crawl index on a typo, which
    # is the round-2 hazard restored by a single flag. --force's own help text
    # presupposes a directory ("really is that empty"), so this check is
    # outside its scope by construction.
    if not shard_dir.is_dir():
        print(f"REFUSING: {shard_dir} is not a directory. If the crawl wrote "
              "elsewhere, set FERRO_EXPAND_OUT or pass --shards. This refusal "
              "is not overridable by --force.")
        return 2

    if not force:
        if not stored:
            print("REFUSING: no records were read from the shard directory, so "
                  "EVERY crawl row looks orphaned. That is what an unmounted or "
                  "wrong directory looks like, not what a repair looks like.")
            return 2
        if total and len(orphans) / total > REPAIR_MAX_FRACTION:
            print(f"REFUSING: {len(orphans)/total:.1%} of crawl rows look "
                  f"orphaned, over the {REPAIR_MAX_FRACTION:.0%} ceiling. A real "
                  "repair touches a small fraction; a wrong shard directory "
                  "touches nearly all of them. Check the path, then --force.")
            return 2

    if dry_run:
        print(f"dry run: would drop {len(orphans):,} rows")
        return 0
    if not orphans:
        print("nothing to repair")
        return 0

    c = sqlite3.connect(index, timeout=120)
    try:
        c.executemany("DELETE FROM held WHERE kind=? AND key=? AND source LIKE ?",
                      [(k, v, OWN_SOURCE_PREFIX + "%") for k, v in orphans])
        c.commit()
    finally:
        c.close()
    print(f"repaired: dropped {len(orphans):,} index rows with no stored record")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--shards", type=Path, default=DEFAULT_SHARDS)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--sample", type=int, default=500)
    ap.add_argument("--repair", action="store_true",
                    help="delete index rows this crawl wrote for records that "
                         "are not on disk, so the next run fetches them again")
    ap.add_argument("--dry-run", action="store_true",
                    help="with --repair: report what would be deleted, delete nothing")
    ap.add_argument("--force", action="store_true",
                    help="with --repair: proceed past the safety refusals. Only "
                         "when you have confirmed the shard directory is the "
                         "right one and really is that empty.")
    a = ap.parse_args()

    # Checked BEFORE --repair, not after: a missing index used to reach
    # sqlite3 directly and raise, instead of the message written for it.
    if not a.index.exists():
        print(f"no identity index at {a.index}; build it with corpus_identity_index.py")
        return 2

    if a.repair:
        return repair(a.index, a.shards, force=a.force, dry_run=a.dry_run)

    r = verify(a.shards, a.index, a.sample)
    pct = 100 * r["duplicates"] / r["written"] if r["written"] else 0.0
    print(f"records written : {r['written']:,}")
    print(f"true duplicates : {r['duplicates']:,}  ({pct:.3f}%)")
    if r["duplicate_examples"]:
        print(f"  e.g. {r['duplicate_examples']}")
    print(f"control         : {r['control_hits']}/{r['control_checked']} known-held report HELD")
    if r["truncated_shards"]:
        # An open gzip has no trailer, so the shard being written always reads
        # as truncated DURING a run. Labelling it keeps a real truncation from
        # being lost among the expected ones -- but only if the label is
        # conditional. The first version tagged the highest-numbered truncated
        # shard unconditionally, so with NO crawl running (the state after a
        # kill, which is exactly when you audit) it stamped "open, being
        # written" on the shard most likely to be genuinely damaged. A guard
        # that reassures you about the thing it was built to flag is worse
        # than none.
        live = r["live_shard"]
        print(f"truncated shards: {len(r['truncated_shards'])}")
        for name, n, err in r["truncated_shards"]:
            tag = "  <- open, being written" if name == live else ""
            print(f"  {name}  {n:,} records recovered  ({err}){tag}")
    print(f"orphaned rows   : {r['orphaned_index_rows']:,}"
          + ("  (index claims articles whose bytes are not on disk)"
             if r["orphaned_index_rows"] else ""))

    if r["orphaned_index_rows"]:
        print("FAIL: the index claims articles this crawl did not store; "
              "re-run with --repair to drop those rows so they are fetched again")
        return 1
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
