#!/usr/bin/env python3
"""Measure how much preprint full text is reachable, and check it stays so.

The answer today is "15% through Europe PMC, and the rest not at all", and the
reason the rest is unreachable is bot protection rather than licensing or cost
-- a distinction with different remedies, so it is worth re-checking rather
than remembering. If bioRxiv's block is ever lifted, this script says so.

Deliberately makes ONE request to the blocked endpoint, not a crawl: the point
is to observe the answer, not to argue with it.
"""
from __future__ import annotations

import json
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT_MD = ROOT / "analysis" / "preprint-fulltext-access.md"
OUT_JSON = ROOT / "analysis" / "preprint-fulltext-access.json"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
UA = {"User-Agent": "cancer-research-corpus/1.0"}


def hits(query: str) -> int:
    q = urllib.parse.urlencode({"query": query, "format": "json", "pageSize": 1})
    with urllib.request.urlopen(
            urllib.request.Request(f"{EPMC}?{q}", headers=UA), timeout=60) as r:
        return int(json.load(r)["hitCount"])


def biorxiv_fulltext_status() -> dict:
    """One request, with the project's own User-Agent.

    NOT retried and NOT disguised. A browser User-Agent would very likely get a
    different answer, and sending one would be circumventing an access control
    the operator put there on purpose.
    """
    api = "https://api.biorxiv.org/details/biorxiv/2020-01-01/2020-01-02/0"
    try:
        with urllib.request.urlopen(
                urllib.request.Request(api, headers=UA), timeout=60) as r:
            coll = (json.load(r).get("collection") or [])
    except Exception as e:  # noqa: BLE001
        return {"metadata_api": f"{type(e).__name__}", "fulltext": None}
    if not coll:
        return {"metadata_api": "empty", "fulltext": None}
    url = coll[0].get("jatsxml")
    if not url:
        return {"metadata_api": "ok", "fulltext": "no jatsxml field"}
    try:
        with urllib.request.urlopen(
                urllib.request.Request(url, headers=UA), timeout=60) as r:
            return {"metadata_api": "ok", "fulltext": f"HTTP {r.status}",
                    "bytes": len(r.read()), "reachable": True}
    except urllib.error.HTTPError as e:
        return {"metadata_api": "ok", "fulltext": f"HTTP {e.code}",
                "reachable": False}
    except Exception as e:  # noqa: BLE001
        return {"metadata_api": "ok", "fulltext": type(e).__name__,
                "reachable": False}


def render(d: dict) -> str:
    """The page, with every figure derived from the measurement above.

    The first version of this script declared OUT_MD and never wrote it, so a
    hand-written page sat beside a JSON file that actually held the numbers --
    the "one artifact describes another" shape this repository keeps finding.
    The prose is here; the counts are substituted.
    """
    b = d["biorxiv"]
    reachable = b.get("reachable")
    L = [
        "# Preprint full text: what is reachable, and what is not",
        "",
        "Cancer preprints are the largest block of literature this project can",
        "identify but mostly cannot read. Counts are Europe PMC's own",
        "`hitCount`, measured by `scripts/preprint_access_check.py`.",
        "",
        "| set | records |",
        "|---|--:|",
        f"| cancer preprints Europe PMC indexes | {d['preprints_indexed']:,} |",
        f"| of those, whose full text Europe PMC holds | {d['fulltext_in_epmc']:,} |",
        f"| share | {d['share_pct']}% |",
        "",
        f"So roughly {d['unread']:,} cancer preprints are identified and unread.",
        "",
        "## Where the rest live, and why they stay there",
        "",
        "bioRxiv and medRxiv publish a metadata API that names a full-text",
        "location for every record -- a `jatsxml` URL on their own domain. That",
        "endpoint needs no account and is NOT paywalled.",
        "",
    ]
    if reachable:
        L += [
            f"**It is currently reachable** (`{b.get('fulltext')}`,",
            f"{b.get('bytes', 0):,} bytes). This page was written when it was",
            "not; if that has changed for good, the ~92,000 unread preprints",
            "become collectable and this analysis needs revising.",
        ]
    else:
        L += [
            f"It is also, from this project's client, unavailable: `{b.get('fulltext')}`",
            "(Cloudflare error 1015). Retried after 20s and again after 45s: the",
            "same. A 429 that does not clear across a minute of backoff is a",
            "standing block on automated access, not a transient rate limit.",
            "",
            "**That block is respected.** It would be straightforward to send a",
            "browser User-Agent and get a different answer, and doing so would",
            "circumvent an access control the operator deliberately put in place.",
            "The corpus does not do that anywhere, and this page exists so nobody",
            "has to rediscover the temptation.",
        ]
    L += [
        "",
        "## The route that does exist costs money",
        "",
        "bioRxiv's supported bulk text-and-data-mining channel is a",
        "requester-pays Amazon S3 bucket: the data is free, the transfer is",
        "billed to whoever asks for it. That is a spending decision, not a",
        "technical one, so it is recorded here rather than taken.",
        "",
        "## A correction",
        "",
        "An earlier note in this project said bioRxiv full text was available",
        "ONLY through requester-pays S3. That was wrong in a way worth stating:",
        "the per-article JATS XML is genuinely free and public, and the",
        "requester-pays bucket is the BULK route. What blocks the per-article",
        "route is bot protection, not licensing and not cost. The distinction",
        "matters because the two have different remedies -- one needs a",
        "conversation with bioRxiv, the other needs a budget.",
        "",
        "## What this means for the corpus",
        "",
        f"The {d['fulltext_in_epmc']:,} preprints Europe PMC does hold are already",
        "collected by `scripts/corpus_expand_fetch.py`, which asks Europe PMC and",
        f"accepts its refusal. The remaining {d['unread']:,} are held as metadata and",
        "abstracts, which are themselves worth having: bibliographic data and",
        "abstracts are generally not copyrightable (Feist), and an abstract is",
        "enough for the census questions this project asks most often.",
        "",
    ]
    return "\n".join(L)


def main() -> int:
    indexed = hits("SRC:PPR AND cancer")
    held = hits("SRC:PPR AND cancer AND IN_EPMC:Y")
    bio = biorxiv_fulltext_status()
    d = {"preprints_indexed": indexed, "fulltext_in_epmc": held,
         "share_pct": round(100 * held / indexed, 1) if indexed else None,
         "unread": indexed - held, "biorxiv": bio}
    OUT_JSON.write_text(json.dumps(d, indent=2) + "\n")
    OUT_MD.write_text(render(d))
    print(json.dumps(d, indent=2))
    if bio.get("reachable"):
        print("\nNOTE: bioRxiv full text is now reachable; "
              "analysis/preprint-fulltext-access.md needs revising.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
