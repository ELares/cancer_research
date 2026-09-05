#!/usr/bin/env python3
"""Fetch cancer clinical trial records from ClinicalTrials.gov.

WHY THIS SOURCE, AND WHY IT IS DIFFERENT FROM THE OTHERS
--------------------------------------------------------
`config.py` has declared CLINICALTRIALS_API since the beginning with zero
callers, the same way EUROPEPMC_API sat unused. It is worth opening for a
reason no journal source can match: a registry entry is the study as
REGISTERED -- its planned arms, endpoints and status -- which is the record
that exists whether or not anyone published the result. Publication bias is
exactly the gap a literature census cannot see from the literature, and 123,069
cancer studies are the population it is invisible against.

It is also the one source in this project with NO redistribution question at
all. ClinicalTrials.gov is a work of the United States government, not subject
to domestic copyright, and its terms explicitly permit reuse. Everywhere else
this project has to separate "readable" from "republishable"; here they are the
same set.

WHAT IT DOES NOT DO
-------------------
It does not turn trials into evidence. A registration is a statement of intent:
many trials never report, some report elsewhere, and a status field is not an
outcome. The corpus stores what the registry says and nothing more.

NOTHING IS FETCHED TWICE, AND THE DURABILITY CONTRACT IS IMPORTED
-----------------------------------------------------------------
The write-then-index ordering, the rotating flushed shards and the transient
failure handling are IMPORTED from `corpus_expand_fetch` rather than
reimplemented. Every one of those was got wrong at least once and fixed under
review; a second copy would be a second chance to get them wrong, and the two
would drift. This module contributes the source, the identity key and the
record shape, and borrows everything that decides what is durable.
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
    _RateLimit, Shards, TransientFetchError, _get,
)
from corpus_identity_index import connect as id_connect  # noqa: E402

API = "https://clinicaltrials.gov/api/v2/studies"
OUT_ROOT = Path(os.getenv(
    "FERRO_TRIALS_OUT", str(Path.home() / "nas" / "cancer-atlas" / "trials")))
PAGE = 100          # the API's documented maximum
SLEEP = float(os.getenv("FERRO_TRIALS_SLEEP", "0.25"))
SHARD_RECORDS = 4000

# The identity index keys on PMID, PMC id, DOI and Europe PMC id -- all
# ARTICLE namespaces. A trial is not an article and its NCT id collides with
# none of them, so it gets its own kind rather than being forced into one.
KIND = "nct"


def _page(token: str | None, condition: str):
    q = {"query.cond": condition, "pageSize": PAGE, "format": "json"}
    if token:
        q["pageToken"] = token
    else:
        q["countTotal"] = "true"
    raw = _get(f"{API}?{urllib.parse.urlencode(q)}", tries=4, timeout=120)
    if raw is None:
        # The refusal set is for ONE item. A registry search refusing is not
        # "no more trials", and treating it as one would silently truncate.
        raise TransientFetchError(f"search refused at token={token!r}")
    return json.loads(raw)


def _record(study: dict) -> dict | None:
    """The fields worth keeping, flattened. Returns None without an NCT id."""
    proto = study.get("protocolSection") or {}
    ident = proto.get("identificationModule") or {}
    nct = ident.get("nctId")
    if not nct:
        return None
    status = proto.get("statusModule") or {}
    design = proto.get("designModule") or {}
    desc = proto.get("descriptionModule") or {}
    arms = proto.get("armsInterventionsModule") or {}
    outcomes = proto.get("outcomesModule") or {}
    refs = proto.get("referencesModule") or {}
    # A registry entry that cites its own publications is the join to the
    # article corpus -- the one field that makes a trial findable from a paper.
    pmids = [r["pmid"] for r in (refs.get("references") or []) if r.get("pmid")]
    return {
        "source": "CTGOV",
        "nct_id": nct,
        "title": ident.get("briefTitle"),
        "official_title": ident.get("officialTitle"),
        "status": status.get("overallStatus"),
        "why_stopped": status.get("whyStopped"),
        "start": (status.get("startDateStruct") or {}).get("date"),
        "completion": (status.get("completionDateStruct") or {}).get("date"),
        "phases": design.get("phases"),
        "study_type": design.get("studyType"),
        "enrollment": (design.get("enrollmentInfo") or {}).get("count"),
        "conditions": (proto.get("conditionsModule") or {}).get("conditions"),
        "interventions": [i.get("name") for i in (arms.get("interventions") or [])],
        "primary_outcomes": [o.get("measure")
                             for o in (outcomes.get("primaryOutcomes") or [])],
        "summary": desc.get("briefSummary"),
        "pmids": pmids,
        "has_results": "resultsSection" in study,
        "fetched": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }


def run(condition: str = "cancer", limit: int | None = None,
        verbose: bool = True) -> dict:
    ident = id_connect()
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    shards = Shards(OUT_ROOT, "trials")
    limiter = _RateLimit(SLEEP)
    totals = {"seen": 0, "new": 0, "with_results": 0, "linked_pmids": 0}
    token, total = None, None
    try:
        while True:
            limiter.wait()
            d = _page(token, condition)
            if total is None and d.get("totalCount"):
                total = d["totalCount"]
                if verbose:
                    print(f"  {total:,} studies match {condition!r}", flush=True)
            studies = d.get("studies") or []
            if not studies:
                break
            for study in studies:
                totals["seen"] += 1
                rec = _record(study)
                if rec is None:
                    continue
                held = ident.execute(
                    "SELECT 1 FROM held WHERE kind=? AND key=?",
                    (KIND, rec["nct_id"])).fetchone()
                if held:
                    continue
                # WRITE THEN INDEX, the same ordering the article crawler uses
                # and for the same reason: write() returns only once the bytes
                # are flushed, so the index never claims a record whose bytes
                # could still be lost.
                shards.write(rec)
                ident.execute(
                    "INSERT INTO held(kind,key,source,has_fulltext) VALUES (?,?,?,0) "
                    "ON CONFLICT(kind,key) DO NOTHING",
                    (KIND, rec["nct_id"], "trials-ctgov"))
                totals["new"] += 1
                totals["with_results"] += bool(rec["has_results"])
                totals["linked_pmids"] += len(rec["pmids"])
                if limit and totals["new"] >= limit:
                    ident.commit()
                    return totals
            ident.commit()
            if verbose and totals["seen"] % 2000 < PAGE:
                print(f"  seen {totals['seen']:,}  new {totals['new']:,}", flush=True)
            token = d.get("nextPageToken")
            if not token:
                break
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
    ap.add_argument("--condition", default="cancer")
    ap.add_argument("--limit", type=int, default=None,
                    help="stop after this many NEW records (a dry run)")
    a = ap.parse_args()
    t = run(a.condition, limit=a.limit)
    print(json.dumps(t, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
