#!/usr/bin/env python3
"""What the expansion actually gathered, by source and by licence.

A bulk download reports itself in bytes, and bytes are the least interesting
thing about it. The questions worth answering are what came in, under what
terms it may be used, and how much of it carries full text rather than a
title -- because those decide what the corpus can be USED for, and because
"free to read" and "ours to redistribute" are not the same fact.

Reads only the shards on disk. Nothing here re-queries a source, so the report
can be regenerated offline and cannot disagree with what was stored.
"""
from __future__ import annotations

import argparse
import collections
import gzip
import json
import os
import re
import sqlite3
import subprocess
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
OUT_ROOT = Path(os.getenv("FERRO_EXPAND_OUT",
                          str(Path.home() / "nas" / "cancer-atlas" / "expanded")))
STATE_DB = REPO / "corpus" / "atlas" / "expand_state.sqlite"
OUT_MD = REPO / "analysis" / "corpus-expansion.md"
OUT_JSON = REPO / "analysis" / "corpus-expansion.json"

# Which licences permit REDISTRIBUTION is decided by `_redistributable` below,
# from decomposed clause tokens. Two module constants used to hold the answer
# as literal strings; they were left behind by the tokeniser rewrite, referenced
# by nothing, and had drifted out of agreement with the code (`"cc pd"` was
# listed as redistributable while the function refused it). A constant that no
# longer describes the behaviour beside it is worse than no constant, so the
# rule lives in one place now.
def _licence_tokens(lic: str) -> set:
    """A licence string reduced to comparable clause tokens.

    EVERY separator becomes whitespace, including `/`, so a licence URL is
    tokenised rather than truncated. An earlier version split on `/` and on
    version separators and DISCARDED everything after the match, which threw
    clauses away: `"cc by 4.0 nd"` became `"cc by"` and classified as
    redistributable, and every CC URL collapsed to `https:`.

    Tokens are then DECOMPOSED, because a clause can arrive glued to its
    neighbours. `cc by-nd4.0` and `CC BY-NCND` both reach an exact-token test
    as one opaque word (`nd4.0`, `ncnd`) and sail past a check for `nd`, which
    is how a NoDerivatives licence was still classified redistributable after
    the tokeniser was supposedly fixed.
    """
    s = (lic or "").strip().lower()
    for ch in "-_/,;()":
        s = s.replace(ch, " ")
    out = set()
    for t in s.split():
        t = t.strip(".:")
        if not t or t in _NOISE or _VERSION.match(t):
            continue
        # The Creative Commons host IS the "cc" in "cc by": dropping it as URL
        # noise left `.../licenses/by/4.0/` with the bare token `by`.
        if t in _CC_HOST:
            out.add("cc")
            continue
        out |= _decompose(t)
    return out


def _decompose(tok: str) -> set:
    """One token to the clause codes it actually contains.

    `nd4.0` -> {nd}; `ncnd` -> {nc, nd}; `cc0` -> {cc0} (NOT {cc, 0}, which is
    why the version strip is applied after the cc0 check). Anything that is not
    a recognised concatenation is returned unchanged, so an unknown publisher
    string is never invented into a clause.
    """
    if tok in _ATOMIC:
        return {tok}
    glued = _GLUED_VERSION.match(tok)
    if glued:
        tok = glued.group(1)
        if tok in _ATOMIC:
            return {tok}
    if _CLAUSE_RUN.fullmatch(tok):
        return set(_CLAUSE_CODE.findall(tok))
    return {tok}


_NOISE = {"http", "https", "www", "licenses", "license", "licence",
          "org", "int", "deed", "en"}
_CC_HOST = {"creativecommons.org", "creativecommons"}
_VERSION = re.compile(r"^v?\d+(\.\d+)*$")
# A clause with a version glued straight onto it: `nd4.0`, `by4`, `sa3.0`.
_GLUED_VERSION = re.compile(r"^([a-z]+?)v?\d+(?:\.\d+)*$")
# Two or more clause codes run together: `ncnd`, `bync`, `bysa`.
_CLAUSE_RUN = re.compile(r"(?:nc|nd|sa|by){2,}")
_CLAUSE_CODE = re.compile(r"nc|nd|sa|by")
_ATOMIC = {"cc0"}

# Words that mean the opposite of an open licence wherever they appear. Without
# these, "Zero rights granted, all rights reserved" matched on `zero` and
# "not in the public domain" matched on `public domain` -- both classified as
# redistributable by a rule looking only for a positive marker.
_DISQUALIFYING = {"not", "no", "none", "reserved", "copyright", "restricted",
                  "subscription", "proprietary", "except", "unless"}

_RESTRICTIVE = {"nc", "nd", "noncommercial", "noderivatives", "noderivs",
                "noderiv"}


def _redistributable(lic: str) -> bool:
    """May this article be REPUBLISHED, not merely read.

    Restrictive clauses and disqualifying words are checked FIRST and before
    any positive marker, because that is the direction where being wrong
    publishes an overclaim. Reordering these -- putting the `cc0`/`zero`
    shortcut above the restrictive check -- is a mutation the tests catch.
    """
    toks = _licence_tokens(lic)
    if not toks:
        return False
    if toks & _RESTRICTIVE or toks & _DISQUALIFYING:
        return False
    if "cc0" in toks:
        return True
    # `zero` and `publicdomain` only count as Creative Commons markers, not as
    # ordinary English words appearing in a sentence.
    cc = bool(toks & {"cc", "publicdomain"})
    if cc and "zero" in toks:
        return True
    if "publicdomain" in toks or {"public", "domain"} <= toks:
        return True
    # CC BY and CC BY-SA. `by` alone is the English word, not a licence.
    return "by" in toks and cc


def scan() -> dict:
    by_source = collections.Counter()
    with_text = collections.Counter()
    licences = collections.Counter()
    years = collections.Counter()
    text_bytes = 0
    n = 0
    redist = redist_text = 0
    for shard in sorted(OUT_ROOT.rglob("*.jsonl.gz")):
        try:
            with gzip.open(shard, "rt", encoding="utf-8", errors="replace") as fh:
                for ln in fh:
                    try:
                        r = json.loads(ln)
                    except Exception:
                        continue
                    n += 1
                    src = r.get("source") or "?"
                    by_source[src] += 1
                    lic = (r.get("licence") or "none").strip().lower()
                    licences[lic] += 1
                    if r.get("year"):
                        years[str(r["year"])] += 1
                    if r.get("has_text"):
                        with_text[src] += 1
                        text_bytes += len(r.get("text") or "")
                    if _redistributable(lic):
                        redist += 1
                        redist_text += bool(r.get("has_text"))
        except (OSError, EOFError, gzip.BadGzipFile):
            continue
    slices = []
    if STATE_DB.exists():
        c = sqlite3.connect(STATE_DB)
        try:
            slices = [dict(zip(("src", "seen", "kept", "fulltext", "done", "slices"), row))
                      for row in c.execute(
                          "SELECT src, SUM(seen), SUM(kept), SUM(fulltext), SUM(done), COUNT(*) "
                          "FROM slice GROUP BY src ORDER BY SUM(kept) DESC")]
        except sqlite3.Error:
            pass
        c.close()
    # Stamped because a crawl that runs for days makes every intermediate
    # report look like a final one. `crawl_running` is asked of the OS rather
    # than assumed, so a snapshot cannot silently claim to be complete.
    running = bool(subprocess.run(
        ["pgrep", "-f", "corpus_expand_fetch.py"],
        capture_output=True, text=True).stdout.strip())
    return {
        "generated": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        "crawl_running": running,
        "records": n,
        "with_fulltext": sum(with_text.values()),
        "fulltext_chars": text_bytes,
        "redistributable": redist,
        "redistributable_with_text": redist_text,
        "by_source": dict(by_source.most_common()),
        "fulltext_by_source": dict(with_text.most_common()),
        "licences": dict(licences.most_common()),
        "years": dict(sorted(years.items(), reverse=True)[:15]),
        "slices": slices,
        "shards": len(list(OUT_ROOT.rglob("*.jsonl.gz"))),
    }


def _n(v):
    """Rank key for a cell that may be a bare count or a {records, ...} dict.

    Ties break on the label so equal counts still have ONE published order;
    without that, two sources with the same count swap places between runs and
    the artifact is order-dependent again in the one case hardest to notice.
    """
    return v["records"] if isinstance(v, dict) else v


def render(d: dict) -> str:
    L = [
        "# Corpus expansion: what came in (#corpus-expand)",
        "",
        "*Generated by `scripts/corpus_expand_report.py`. Reads only the stored "
        "shards, so it cannot disagree with what is on disk.*",
        "",
        f"**Snapshot taken {d['generated']}.**" + (
            " The crawl was still running when this was written, so these are"
            " running totals and not a final result." if d["crawl_running"] else
            " No crawl was running when this was written."),
        "",
        f"**{d['records']:,} records new to this project**, of which "
        f"**{d['with_fulltext']:,} carry full text** "
        f"({d['fulltext_chars']/1e9:.2f} GB of it), across {d['shards']} shards.",
        "",
        "Every record was checked against `corpus_identity_index` on PMID, PMC "
        "id and DOI before it was fetched. That check is a claim about what the "
        "fetcher did, and this report cannot confirm it -- it reads shards, not "
        "the index. `scripts/corpus_expand_verify.py` is what confirms it, by "
        "re-testing every stored record against holdings that predate the crawl.",
        "",
        "## By source",
        "",
        "| source | records | with full text |",
        "|---|--:|--:|",
    ]
    names = {"MED": "MEDLINE / PubMed", "PMC": "PubMed Central",
             "PPR": "preprints", "PAT": "patents", "ETH": "theses (EThOS)",
             "AGR": "Agricola", "CBA": "Chinese Biological Abstracts",
             "HIR": "NHS evidence", "CTX": "CiteXplore"}
    for src, cnt in sorted(d["by_source"].items(),
                           key=lambda kv: (-_n(kv[1]), kv[0])):
        L.append(f"| {names.get(src, src)} | {cnt:,} | "
                 f"{d['fulltext_by_source'].get(src, 0):,} |")
    L += [
        "",
        "## By licence, and why that is the column that matters",
        "",
        "Free to read and free to redistribute are different facts, and this "
        "project's mission statement turns on the difference. The counts are "
        "kept apart rather than summed.",
        "",
        "| licence | records |",
        "|---|--:|",
    ]
    for lic, cnt in sorted(d["licences"].items(),
                           key=lambda kv: (-_n(kv[1]), kv[0])):
        L.append(f"| {lic} | {cnt:,} |")
    L += [
        "",
        f"**{d['redistributable']:,} records carry a licence permitting "
        f"redistribution** ({d['redistributable_with_text']:,} of them with full "
        "text). The remainder was still legitimately fetched -- Europe PMC "
        "serves only what it is allowed to serve, and refused the rest with a "
        "404 -- but it is not all onward-shareable, and a later release must "
        "read this column rather than the total.",
        "",
        "## Progress",
        "",
        "| source | pages seen | new | full text | slices done |",
        "|---|--:|--:|--:|--:|",
    ]
    # Ranked HERE, like the other two tables. The stored list arrives from
    # SQL whose ORDER BY leaves ties unspecified, so two equal sources could
    # swap places between runs -- an order-dependence a list cannot be
    # shuffle-tested for, since shuffling a list preserves it. Sorting on the
    # ranked quantity with a label tie-break is what makes the published order
    # a property of the renderer rather than of the query plan.
    for s in sorted(d["slices"], key=lambda r: (-r["kept"], -r["seen"], r["src"])):
        L.append(f"| {s['src']} | {s['seen']:,} | {s['kept']:,} | "
                 f"{s['fulltext']:,} | {s['done']}/{s['slices']} |")
    L += [
        "",
        "## What was not taken",
        "",
        "- **Paywalled full text.** Never requested from a publisher; the only "
        "text endpoint asked is Europe PMC's, which serves the open subset and "
        "returns 404 otherwise. The boundary is enforced by the source, not by "
        "a judgement made here.",
        "- **Bronze content.** Free to read on a publisher's own page, no "
        "licence attached, and not served by an open API. Its size is NOT "
        "measured here and the figure this line used to carry (\"roughly 12% "
        "of census records that have a DOI and no PMC id\") rested on no "
        "analysis in this repository. What IS measured, in "
        "`analysis/europepmc-access-ceiling.md`, is the gap bronze sits "
        "inside: 57.9% of Europe PMC's cancer records are readable there "
        "against 21.3% openly licensed.",
        "- **Anything already held.** Checked on three identifiers before "
        "fetching, not after.",
        "",
    ]
    return "\n".join(L) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()
    d = scan()
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n")
    OUT_MD.write_text(render(d))
    if not a.quiet:
        print(f"{d['records']:,} new records | {d['with_fulltext']:,} with full text "
              f"| {d['fulltext_chars']/1e9:.2f} GB | {d['redistributable']:,} redistributable")
        print(f"  by source: {d['by_source']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
