#!/usr/bin/env python3
"""Fetch cancer literature this project does not already hold.

WHAT IT WILL AND WILL NOT TAKE
------------------------------
Only content served by an API built to serve it. Europe PMC's full-text
endpoint returns text for the open-access subset and refuses the rest, so the
paywall boundary is enforced by the source rather than by this script guessing.
Publisher pages are never fetched.

That line matters more than "is it free to read". Readable and redistributable
are different sets, and the gap between them is measured in
`analysis/europepmc-access-ceiling.md`: 57.9% of Europe PMC's cancer records
are readable there against 21.3% carrying a licence that permits
republication. The difference is largely BRONZE -- free on the publisher's
site with no licence at all -- which is exactly the case where "free" and
"ours to keep" come apart. MISSION.md says neither scraping nor redistributing
that is intended. Every record stores the licence Europe PMC reports, so what
was taken under what terms stays answerable.

(An earlier version of this docstring quoted ~22% readable and ~6%
redistributable "of census records with a DOI and no PMC id", attributed to a
survey by this project. No such analysis exists in this repository, and the
two figures were not consistent with each other. They are replaced above by
numbers a committed artifact derives.)

NOTHING IS FETCHED TWICE
------------------------
Every candidate is checked against `corpus_identity_index` on all three of
PMID, PMC id and DOI before anything is downloaded, and every item written is
added to that index. A re-run therefore skips what it already has, and an
interrupted run resumes from its stored cursor rather than from the beginning.

The ORDER of those two writes is load-bearing and is not a transaction. A
gzip buffer and a SQLite transaction cannot be made atomic with each other, so
the record's bytes are flushed to disk BEFORE the index is told the record is
held (see `Shards`). Getting this backwards -- as an earlier version did --
means a kill leaves the index claiming articles whose bytes never landed, and
those articles are then skipped forever. Ordered this way the residual failure
is a re-fetch, which `corpus_expand_verify.py` reports as a duplicate rather
than hiding as silence.

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
import threading
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
# Bumped when a migration in `state()` must run once over an existing ledger.
SCHEMA_VERSION = 1
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
# WHAT COUNTS AS A CANCER PAPER, for the purpose of asking Europe PMC.
#
# The crawler asked for the literal word `cancer` and nothing else, which
# sounds complete and is not: measured against this set, a one-word query sees
# 5,575,342 records where the set sees 9,872,234. FOUR MILLION of the
# difference say `carcinoma`, `leukemia`, `melanoma` or `glioblastoma` and
# never once say `cancer` -- 593,583 for carcinoma alone -- and 1,702,738 of
# the missed pile have full text available. The crawl was structurally blind to
# 44% of its own subject, and no amount of paging or concurrency reaches a
# record the query never asks for.
#
# GROUNDED IN THE CENSUS DEFINITION rather than invented beside it: these stems
# are checked against the 704 MeSH C04 descriptor names this project already
# commits as its definition of cancer, and cover 79.5% of them
# (tests/test_corpus_expand_fetch.py measures it, so the figure cannot drift).
# The residue is mostly named entities -- astrocytoma, apudoma, angiokeratoma --
# which is why the `-oma` stem is included despite catching glaucoma and
# trachoma too. Recall is the priority here and the cost of a false positive is
# one stored record, while the cost of a false negative is a paper this project
# can never see.
#
# MeSH was tried first and rejected on measurement, not preference: Europe PMC's
# `MESH:` field matches exact descriptors, and OR-ing 25 of them returned FEWER
# hits than OR-ing 5, so the field does not behave like the tree the census uses.
CANCER_TERMS = (
    "cancer", "neoplasm", "neoplasms", "neoplastic", "tumor", "tumour",
    "tumors", "tumours", "carcinoma", "carcinomas", "sarcoma", "leukemia",
    "leukaemia", "lymphoma", "melanoma", "myeloma", "glioma", "glioblastoma",
    "blastoma", "adenocarcinoma", "adenoma", "oncology", "oncological",
    "malignancy", "malignant", "metastasis", "metastatic", "mesothelioma",
    "teratoma", "thymoma", "chordoma", "papilloma",
)


def cancer_clause() -> str:
    """The subject half of every search, as one parenthesised OR."""
    return "(" + " OR ".join(CANCER_TERMS) + ")"


def query_fingerprint() -> str:
    """Identifies WHICH question a slice was walked with.

    A slice marked done under a narrower query is not done under a wider one --
    it was fully walked for a question nobody is asking any more. Without this,
    widening the subject clause would leave every completed slice permanently
    holding its old, smaller answer, and the 44% blind spot would be fixed for
    future slices and invisible in past ones.
    """
    import hashlib
    return hashlib.sha256("|".join(sorted(CANCER_TERMS)).encode()).hexdigest()[:12]


PAGE = 1000
# Consecutive per-article failures that mean "the service is down" rather than
# "this article is briefly unavailable". One 503 must never stop a crawl; a
# hundred in a row must, or the run quietly degrades into a metadata-only pass.
MAX_CONSECUTIVE_DEFERRALS = 100
# Ceiling on a backoff we chose ourselves. Reachable: at tries=4 the schedule
# is 2/4/8, and a longer `tries` walks into it -- the previous cap of 60.0 was
# unreachable at every value this code uses, so it bounded nothing.
BACKOFF_CAP = 30.0
# Full-text fetches running at once. The crawler is latency-bound -- median
# fetch 1.52s, measured, with a long tail (one sampled request took 16.9s) --
# so this is the only lever that materially changes how long a full pass takes.
# The rate never exceeds SLEEP's ceiling however high this goes, so raising it
# buys utilisation, not requests. FERRO_EXPAND_WORKERS=1 restores the strictly
# serial behaviour.
FETCH_WORKERS = max(1, int(os.getenv("FERRO_EXPAND_WORKERS", "8")))
# A server-supplied wait is obeyed within these bounds. The floor exists
# because `Retry-After: 0` is an instruction to hammer; the cap because a
# 300-second wait on one article is worse than deferring it.
RETRY_AFTER_FLOOR = 1.0
RETRY_AFTER_CAP = 60.0
# A server may legitimately say "come back tomorrow". OpenAlex's free tier is
# 1,000 requests a day and, when it is spent, answers 429 with
# `Retry-After: 60875` -- seventeen hours, until the budget resets at midnight
# UTC. Capping that at 60s and retrying four times is not politeness with a
# backoff; it is knocking on a door that has told you when it opens. A wait
# this long is honoured in full, up to a day, because the alternative is
# hammering a service that answered the question clearly.
LONG_WAIT_CAP = 86_400.0
# What counts as "not a retry any more" -- past this the wait is a scheduled
# return, and it is logged as one so a sleeping crawler is not read as hung.
LONG_WAIT_THRESHOLD = 120.0
# How many scheduled returns one fetch may keep before it gives up. A
# scheduled return does NOT spend a retry -- the server named a time and we
# kept it, which is not a failed attempt -- so without a separate bound a
# server that always answers with a long wait would hold one worker forever.
# Two is enough to cross one budget reset and notice if the next answer is the
# same, and bounds the worst case at two days rather than unbounded.
MAX_SCHEDULED_RETURNS = 2


# How long to keep trying one search page before giving up on the whole run.
# Generous, because the alternative is a stopped crawl that nobody notices:
# `_get` already spends about 14s on its own retries, so five attempts with
# this backoff covers roughly ten minutes of the service being unavailable.
SEARCH_PAGE_ATTEMPTS = 5
SEARCH_PAGE_BACKOFF = (15.0, 45.0, 90.0, 180.0)


def retry_page(fetch, *args, label: str = "page", verbose: bool = True,
               **kwargs):
    """Fetch ONE page of a paged search, surviving a briefly unavailable service.

    SHARED BY EVERY FETCHER IN THIS PROJECT, and it exists because the same
    defect appeared three times. A paged search that raises on a transient
    failure ends the whole crawl: it cost a two-day Europe PMC run at 302,973
    records (one 503), and then -- after that was fixed HERE and not carried
    across -- it killed the OpenAlex fetcher on its first 429. Writing the fix
    once in a place all three import is the only version of it that stays
    fixed.

    Why a page failure cannot simply be swallowed: an empty page reads as "no
    more results", so returning one would retire the slice, or end the walk,
    with the ledger recording an outage as a completed pass. Raising is right;
    raising IMMEDIATELY is what was wrong.

    Retrying the same cursor is safe because a cursor is only advanced after
    its page is fully processed, so a retry re-reads a page rather than
    skipping one.
    """
    last = None
    for attempt in range(SEARCH_PAGE_ATTEMPTS):
        try:
            return fetch(*args, **kwargs)
        except TransientFetchError as e:
            last = e
            if attempt + 1 >= SEARCH_PAGE_ATTEMPTS:
                break
            wait = SEARCH_PAGE_BACKOFF[min(attempt, len(SEARCH_PAGE_BACKOFF) - 1)]
            if verbose:
                print(f"  ! {label} unavailable "
                      f"(attempt {attempt + 1}/{SEARCH_PAGE_ATTEMPTS}), "
                      f"waiting {wait:.0f}s", flush=True)
            time.sleep(wait)
    raise RuntimeError(
        f"{label} failed {SEARCH_PAGE_ATTEMPTS} times over "
        f"~{sum(SEARCH_PAGE_BACKOFF):.0f}s of backoff; the service looks down "
        f"rather than busy: {last}") from last


def _search_with_retry(src: str, year: int, cursor: str, verbose: bool = True):
    """The Europe PMC search, through the shared page retry."""
    return retry_page(search, src, year, cursor,
                      label=f"search {src} {year}", verbose=verbose)


def _should_abort(consecutive: int) -> bool:
    """Is a run of per-article failures an outage rather than bad luck?

    A function, not an inline comparison, because a guard written inline can
    only be checked by grepping the source for it -- and a grep passes while
    the comparison sits behind `if False:`. That mutation survived until this
    was extracted.
    """
    return consecutive >= MAX_CONSECUTIVE_DEFERRALS
# Minimum interval between full-text request STARTS, shared across workers.
#
# THE EARLIER RATIONALE HERE WAS WRONG AND IS WITHDRAWN. It said lowering this
# to 0.10 was followed by a 503 and treated that as a reason to keep it high.
# Measured afterwards, the median full-text latency is 1.52s -- not the ~0.6s
# assumed -- so a SERIAL crawler at sleep=0.10 was issuing about 0.6 requests
# per second. No rate limiter fires at 0.6 req/s. The 503 was Europe PMC
# having a moment, and the crawl died because one 503 was fatal, which is
# fixed. I blamed my own change for someone else's blip and reasoned from it.
#
# WHAT IS MEASURED: median full-text latency 1.52s over a sample, with a long
# tail (one request took 16.9s). WHAT IS ARITHMETIC: the request rate is
# min(workers/latency, 1/interval), so 8 workers at 0.25s is capped at 4/s and
# reaches about 4/s. That is an order of magnitude under what a single client
# could extract from this API, and it cannot rise however high the worker count
# goes -- the cap is the point of the shared limiter.
#
# WHAT IS NOT MEASURED, and I am not going to imply otherwise: a clean
# throughput A/B between (4, 0.34) and (8, 0.25). Three attempts were each
# contaminated -- two harnesses died mid-run and the third had the crawler
# restarted underneath it by a stray child -- and the obvious metric,
# records-written per hour, depends on how dense new records are in whatever
# region the cursor happens to be in, which is why two early readings
# contradicted each other. The change rests on the latency measurement and the
# rate ceiling, not on a speedup number.
SLEEP = float(os.getenv("FERRO_EXPAND_SLEEP", "0.25"))
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
    -- Records this slice could not fetch text for. A slice with a non-zero
    -- count is NOT finished and must not be marked done: the deferral is the
    -- record of an article we walked past and did not store.
    deferred INTEGER NOT NULL DEFAULT 0,
    -- WHICH question this slice was walked with. A slice fully walked for a
    -- narrower subject clause is not finished for a wider one; without this
    -- column, widening the query fixes the future and leaves every completed
    -- slice holding its old, smaller answer forever.
    query_id TEXT NOT NULL DEFAULT '',
    updated REAL NOT NULL DEFAULT 0,
    PRIMARY KEY (src, year)
);
"""


def state() -> sqlite3.Connection:
    STATE_DB.parent.mkdir(parents=True, exist_ok=True)
    c = sqlite3.connect(STATE_DB, timeout=60)
    c.executescript(STATE_SCHEMA)
    c.execute("PRAGMA journal_mode=WAL")
    # CREATE TABLE IF NOT EXISTS does not add a column to an existing table.
    cols = {r[1] for r in c.execute("PRAGMA table_info(slice)")}
    if "deferred" not in cols:
        c.execute("ALTER TABLE slice ADD COLUMN deferred INTEGER NOT NULL DEFAULT 0")
    if "query_id" not in cols:
        c.execute("ALTER TABLE slice ADD COLUMN query_id TEXT NOT NULL DEFAULT ''")
    _migrate(c)
    c.commit()
    return c


def _migrate(c) -> None:
    """One-time migrations, gated on PRAGMA user_version.

    v1 -- RESTORE THE "NEVER ASKED" SENTINEL (#913). Before #912 an absent
    `hitCount` was written as 0, so `hits = 0` on an existing row means either
    "the service reported an empty slice" or "the response was degraded and the
    total was lost", and NOTHING IN THE LEDGER DISTINGUISHES THEM. The schema
    default is -1 for exactly that reason, and no row carried it: all 139
    zero-total rows predate the fix.

    That ambiguity is what hid the larger half of #913. A row reading
    `seen=0, hits=0, done=1` is indistinguishable from a genuinely empty year,
    and 24 consecutive MED slices (1991-2014) sat in it -- each marked complete
    having walked NOTHING, while the service reports 64,965 records for 1991
    and 274,868 for 2014. Restoring -1 makes them say "unknown" instead of
    "empty", which is what the detector reads.

    Bounded to run ONCE, because after #912 a `hits = 0` written by a real
    response is TRUSTWORTHY -- a genuinely empty slice must be allowed to
    settle at 0 and stop being flagged, so this must never become a standing
    rule.
    """
    if c.execute("PRAGMA user_version").fetchone()[0] >= SCHEMA_VERSION:
        return
    n = c.execute("UPDATE slice SET hits=-1 WHERE hits=0").rowcount
    c.execute(f"PRAGMA user_version={SCHEMA_VERSION}")
    c.commit()
    if n:
        print(f"ledger migration v{SCHEMA_VERSION}: {n} slice(s) carrying a "
              f"zero total written before #912 now read as UNKNOWN (-1); "
              f"run --incomplete to see which cannot be certified complete")


def reopen_stale_slices(st, verbose: bool = True) -> int:
    """Re-open every slice walked with a different subject clause.

    Called once per run. A slice carrying an older fingerprint was walked to
    completion for a question this crawl no longer asks, so `done` is a claim
    about the wrong thing -- it must be withdrawn and the slice re-walked. That
    is cheap: every record already stored is skipped without a fetch, so the
    re-walk costs paging and finds only what the wider query adds.
    """
    fp = query_fingerprint()
    rows = st.execute(
        "SELECT COUNT(*) FROM slice WHERE query_id != ? AND (done=1 OR seen>0)",
        (fp,)).fetchone()[0]
    if rows:
        st.execute(
            "UPDATE slice SET done=0, cursor='*', query_id=? WHERE query_id != ?",
            (fp, fp))
        st.commit()
        if verbose:
            print(f"  subject clause changed ({fp}): re-opening {rows:,} slice(s) "
                  f"walked with an older question", flush=True)
    return rows


# Statuses that mean "this item, no" rather than "try again later". Passed in
# per call site, because the right answer differs between them: for ONE
# article a refusal is the expected answer and the whole design leans on it,
# while for the SEARCH endpoint the same status means a page of results was
# not delivered, and treating that as "no results" silently discards a slice.
REFUSE_ITEM = (401, 403, 404, 410)
REFUSE_SEARCH = (404,)


class TransientFetchError(RuntimeError):
    """One item could not be fetched right now. NOT a reason to stop.

    A crawl measured in weeks died on a single 503, because every give-up
    raised the same RuntimeError and nothing distinguished "this article is
    briefly unavailable" from "the service is gone". The article served a
    154 KB document on the very next request. A crawl that cannot survive one
    bad response from one URL is not a crawl.
    """


def _get(url: str, tries: int = 4, timeout: int = 120,
         refuse: tuple = REFUSE_ITEM):
    """Fetch, retrying transient failures with exponential backoff.

    Returns bytes on success, None when the server REFUSES (see `refuse`), and
    raises TransientFetchError when it kept failing for a reason that may pass.

    A silent give-up would look exactly like a slice with no new records, which
    is the one thing this ledger must never confuse -- so giving up is loud.
    But loud must not mean fatal for the whole crawl: the caller decides
    whether one unavailable article is worth stopping for, and for a single
    full text it never is.

    `Retry-After` is honoured when the server sends it, because guessing a
    backoff while being told the answer is how a rate limiter turns into an
    outage.
    """
    last = None
    # A WHILE LOOP, NOT `for attempt in range(tries)`, and the difference is
    # the whole point of the scheduled-return path below. Under a `for` loop
    # `continue` advances the iterator, so honouring a seventeen-hour wait
    # would SPEND one of `tries` -- the opposite of the intent, and the first
    # version of this change shipped that, guarded only by a test that grepped
    # the source for the word `continue`. A text scan cannot see loop
    # semantics. `attempt` is incremented explicitly below, after a real
    # failure, and a scheduled return skips it.
    attempt = 0
    scheduled = 0
    while attempt < tries:
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            # A REFUSAL IS AN ANSWER, and 404 is not the only way to say no:
            # 401/403/410 are the server declining to serve this ITEM, which is
            # the boundary this crawl leans on rather than deciding for itself.
            # WHICH statuses count is the caller's decision: applied to the
            # SEARCH endpoint a 403 would make an empty page look like the end
            # of a slice, retiring a whole source-year with seen=0.
            if e.code in refuse:
                return None
            last = e
            wait = _retry_after(e)
            long_wait = _long_retry_after(e)
            if long_wait is not None and scheduled < MAX_SCHEDULED_RETURNS:
                scheduled += 1
                print(f"  server asked for {long_wait / 3600:.1f}h "
                      f"({int(long_wait)}s); waiting rather than retrying "
                      f"[{scheduled}/{MAX_SCHEDULED_RETURNS}]: {url[:70]}",
                      flush=True)
                time.sleep(long_wait)
                # `attempt` deliberately NOT incremented: the server named a
                # time and we kept it, so the retry budget is still intact for
                # the transient failures it exists for.
                continue
        except Exception as e:  # noqa: BLE001 - network, DNS, timeouts, resets
            last = e
            wait = None
        attempt += 1
        if attempt < tries:
            # MEASURED against what it replaced, because the first version of
            # this comment claimed exponential beat the old linear schedule and
            # the arithmetic said otherwise: 2**attempt gives 1+2+4 = 7.0s over
            # four tries, while the linear 1.5*(n+1) it replaced gave
            # 1.5+3+4.5 = 9.0s. Every individual wait was shorter too. Starting
            # at 2 makes it 2+4+8 = 14.0s, which is actually more room.
            time.sleep(wait if wait is not None
                       else min(BACKOFF_CAP, 2.0 ** attempt))
    raise TransientFetchError(f"giving up after {tries}: {url[:110]} :: {last}")


def _long_retry_after(e) -> float | None:
    """A server-supplied wait too long to be a retry: honour it as scheduled.

    Read from the RAW header rather than from `_retry_after`, which clamps to
    RETRY_AFTER_CAP -- clamping is right for a short wait and exactly wrong
    here, because it turns "come back in seventeen hours" into "come back in a
    minute" and then does so four times.
    """
    h = getattr(e, "headers", None)
    if h is None:
        return None
    try:
        v = h.get("Retry-After")
    except Exception:  # noqa: BLE001
        return None
    if v is None:
        return None
    if isinstance(v, bytes):
        v = v.decode("ascii", "replace")
    try:
        secs = float(v)
    except (TypeError, ValueError):
        secs = _http_date_delay(v)
        if secs is None:
            return None
    # `secs != secs` is the NaN test: float("nan") parses, and every
    # comparison against it is False, so an unguarded NaN would fall through
    # as a valid wait and `time.sleep(nan)` raises.
    if secs != secs or secs <= LONG_WAIT_THRESHOLD:
        return None
    return min(LONG_WAIT_CAP, secs)


def _retry_after(e) -> float | None:
    """The server's own answer to "how long should I wait?", when it gives one.

    RFC 7231 allows TWO forms and the first version handled one: `Retry-After`
    may be delay-seconds OR an HTTP-date, and the date form is what CDNs and
    WAFs in front of an origin commonly send on a 503. `float()` raised on it
    and the instruction was silently thrown away.

    FLOORED, not just clamped. `max(0.0, ...)` turned `Retry-After: 0` and any
    negative value into a zero-second wait, so a proxy emitting `Retry-After: 0`
    got four back-to-back requests -- strictly worse than the fixed backoff it
    was overriding. Zero is not a useful instruction to obey literally.
    """
    h = getattr(e, "headers", None)
    if h is None:
        return None
    try:
        v = h.get("Retry-After")
    except Exception:  # noqa: BLE001
        return None
    if v is None or (isinstance(v, str) and not v.strip()):
        return None
    if isinstance(v, bytes):
        v = v.decode("ascii", "replace")
    try:
        secs = float(v)
    except (TypeError, ValueError):
        secs = _http_date_delay(v)
        if secs is None:
            return None
    if secs != secs:  # NaN
        return None
    return min(RETRY_AFTER_CAP, max(RETRY_AFTER_FLOOR, secs))


def _http_date_delay(v) -> float | None:
    """Seconds until an RFC 7231 HTTP-date, or None if it is not one."""
    try:
        import email.utils
        parsed = email.utils.parsedate_to_datetime(str(v))
    except Exception:  # noqa: BLE001
        return None
    if parsed is None:
        return None
    try:
        import datetime as _dt
        now = _dt.datetime.now(_dt.timezone.utc)
        if parsed.tzinfo is None:
            parsed = parsed.replace(tzinfo=_dt.timezone.utc)
        return max(0.0, (parsed - now).total_seconds())
    except Exception:  # noqa: BLE001
        return None


def search(src: str, year: int, cursor: str):
    q = f"SRC:{src} AND PUB_YEAR:{year} AND {cancer_clause()}"
    url = (f"{SEARCH}?format=json&pageSize={PAGE}&resultType=core"
           f"&cursorMark={urllib.parse.quote(cursor)}&query={urllib.parse.quote(q)}")
    raw = _get(url, refuse=REFUSE_SEARCH)
    if raw is None:
        # Only a 404 reaches here, and the search endpoint does not 404 on an
        # empty result set -- it returns hitCount 0. So this is the endpoint
        # itself being wrong, not a slice being empty, and it must not be
        # recorded as one.
        raise RuntimeError(f"search endpoint refused: SRC:{src} PUB_YEAR:{year}")

    return json.loads(raw)


# Licences under which Europe PMC may be ASKED for the text. Deliberately
# wider than the set this project may REPUBLISH -- an NC or ND licence permits
# reading and analysis while forbidding redistribution, and this project's use
# is analysis. `corpus_expand_report._redistributable` decides the narrower
# question, and the two must not be conflated: reporting an NC record as
# redistributable overstates what may be passed on. Every record stores its
# own licence, so the narrower question stays answerable per record.
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
    open-access subset", NOT "openly licensed": 118,307 cancer records carry
    CC BY or CC0 while flagged `OPEN_ACCESS:N` (measured, see
    analysis/europepmc-access-ceiling.md), and gating on the flag alone skipped
    every one of them without a word.

    WHAT WIDENING THIS BUYS, stated honestly rather than both ways. An earlier
    version of this docstring said the change "cannot widen what is taken",
    which contradicts the reason for making it: either the endpoint serves
    those records, in which case more IS taken, or it refuses them, in which
    case the change buys only extra 404s. Both cannot be true. What is
    actually guaranteed is narrower and worth stating on its own: asking more
    often cannot make the endpoint serve something it would otherwise refuse,
    so widening the ASK cannot widen what we are PERMITTED to take. Whether it
    yields anything is an empirical question the ledger answers per source.

    The licence set here is the read-and-analyse set, not the republish set --
    see OPEN_LICENCES.
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
                   tries=4, timeout=120)
        if raw:
            return raw.decode("utf-8", "replace")
    if rec.get("source") == "PPR" and rec.get("id"):
        raw = _get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/PPR/{rec['id']}/fullTextXML",
                   tries=4, timeout=120)
        if raw:
            return raw.decode("utf-8", "replace")
    return None


class _RateLimit:
    """A minimum interval between request STARTS, shared by every worker.

    WITHOUT THIS, ADDING WORKERS DELETES THE RATE LIMIT. The serial path slept
    `SLEEP` after each fetch, so the pool -- which slept not at all -- was
    quietly a request-rate multiplier as well as a latency-hider. A measured
    6.5x on four workers was partly concurrency and partly the politeness pause
    going away, which is the same mistake as shortening it, made larger and
    less visibly. Lowering that pause preceded a 503 that ended a run.

    Spacing admissions rather than completions keeps request latency from
    multiplying the allowed rate. Network requests can overlap after admission;
    a delayed wakeup does not accumulate credit for a burst of admissions.
    """

    def __init__(self, interval: float):
        self.interval = interval
        self._lock = threading.Lock()
        self._next = 0.0

    def wait(self) -> None:
        if self.interval <= 0:
            return
        with self._lock:
            now = time.monotonic()
            while now < self._next:
                time.sleep(self._next - now)
                now = time.monotonic()
            # Space grants from the actual wakeup, not a reserved future slot.
            # A delayed worker must not release alongside workers whose slots
            # also expired while the scheduler was busy. Only the wait is
            # serialized; network requests still run concurrently.
            self._next = now + self.interval


def _fetch_texts(wanted, limiter=None) -> dict:
    """Fetch full text for a page's records, several at a time.

    THE ONLY CONCURRENT PART OF THIS CRAWLER, and deliberately the smallest one
    that helps. The process sits near-idle waiting on Europe PMC, so
    overlapping the waits is the lever; shortening the politeness pause is not,
    and trying that preceded a 503 that ended a run.

    Nothing here touches a shard, the identity index or the ledger. Every
    guarantee that makes this crawler trustworthy -- flush-before-index, the
    dedup lookup, the deferral bookkeeping -- stays on the single writer
    thread, in page order. This function returns data and raises nothing.

    Keyed by `id(rec)` rather than by any field: a record may carry only a
    Europe PMC id, and two records in a page can legitimately share a null.
    Object identity is the one key every record has, and the dict never
    outlives the page.

    Failures are RETURNED, not raised, so the caller decides per record. A
    worker that raised would abandon the rest of the page mid-flight.
    """
    if not wanted:
        return {}
    limiter = limiter or _RateLimit(SLEEP)
    out = {}

    def one(rec):
        limiter.wait()
        try:
            return fetch_fulltext(rec)
        except TransientFetchError as e:
            return e

    workers = max(1, min(FETCH_WORKERS, len(wanted)))
    if workers == 1:
        for rec, *_ in wanted:
            out[id(rec)] = one(rec)
        return out

    import concurrent.futures as _cf
    with _cf.ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(one, rec): rec for rec, *_ in wanted}
        for fut in _cf.as_completed(futures):
            rec = futures[fut]
            try:
                out[id(rec)] = fut.result()
            except Exception as e:  # noqa: BLE001 - one record never kills a page
                out[id(rec)] = TransientFetchError(str(e))
    return out


class Shards:
    """Rotating gzip shards whose bytes are readable as soon as write() returns.

    THE FLUSH IS THE WHOLE POINT, and it was missing. A record used to sit in
    the gzip object's compression buffer while the identity index was told the
    record was held. Kill the process and the index keeps the row while the
    bytes are gone -- so the article is permanently invisible to every later
    run, which is the exact failure the dedup exists to prevent, inverted. It
    is not hypothetical: two `pkill`s during development left shards of 392 and
    2,831 records against a 4,000-record rotation, and every record in the lost
    tails was already marked held.

    So write() ends with a Z_SYNC_FLUSH, which closes the deflate block and
    pushes the bytes to the OS. A killed process then loses nothing; only a
    machine crash can, and `os.fsync` is deliberately NOT called because this
    writes to a network mount where a per-record fsync costs more than the
    failure it would prevent. The residual window is one record wide, and it
    fails in the SAFE direction -- a re-fetch, which the audit reports as a
    duplicate rather than silence.
    """

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
        # Durable BEFORE the caller is allowed to index it. Ordering, not
        # speed, is what makes the index trustworthy.
        self.fh.flush()
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
    # A run this long meets transient failures. One is nothing; a run of them
    # means the service is gone and continuing would walk the whole corpus
    # storing metadata and calling it a crawl.
    deferred = consecutive_deferred = 0
    reopen_stale_slices(st, verbose=verbose)
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    shards = Shards(OUT_ROOT, "expanded")
    totals = {"seen": 0, "new": 0, "fulltext": 0, "slices": 0, "deferred": 0,
              # Records stored WITHOUT full text because a retry also failed.
              # Reported so "the slice completed" never hides how much text
              # was given up on.
              "no_fulltext": 0}
    try:
        for src, oa_possible in sources:
            for year in years:
                row = st.execute(
                    "SELECT cursor, done, deferred FROM slice WHERE src=? AND year=?",
                    (src, year)).fetchone()
                if row and row[1]:
                    continue
                cursor = row[0] if row else "*"
                # IS THIS A RETRY OF A SLICE THAT ALREADY DEFERRED? If so, a
                # full text that fails AGAIN is treated as unavailable rather
                # than deferred a second time, because deferring it again is
                # what made this non-terminating.
                #
                # MEASURED, not assumed. PMC 2023 went 88 deferred -> 7 -> 7
                # across two passes while `seen` went 10,611 -> 21,155 against
                # 10,551 hits: it re-walked the whole slice to retry records
                # that cannot succeed. The first one in both passes was
                # PMC10206698 ("E-Poster abstracts"), which EPMC advertises as
                # isOpenAccess=Y inEPMC=Y and whose fullTextXML returns HTTP
                # 500 on three consecutive direct requests. A permanently
                # broken document is not bad luck, and a slice that waits for
                # it never retires -- for the MED slices that is 2.9M search
                # requests per pass, forever.
                retry_pass = bool(row and row[2])
                st.execute("INSERT OR IGNORE INTO slice(src,year,query_id) "
                           "VALUES (?,?,?)", (src, year, query_fingerprint()))
                seen = kept = ft = 0
                slice_deferred = 0
                slice_abandoned = consecutive_abandoned = 0
                # Page-local bookkeeping: `slice_abandoned` accumulates across
                # the whole slice, so the run total adds only the delta rather
                # than the running count once per page.
                reported_abandoned = 0
                while True:
                    # A SEARCH FAILURE MUST NOT END A CRAWL MEASURED IN WEEKS.
                    #
                    # `search` raises rather than returning an empty page, and
                    # that is right: an empty page reads as "no more results"
                    # and would retire the slice with seen=0. But raising all
                    # the way out of run() is a different mistake, and it cost
                    # a two-day crawl at 302,973 records -- one 503 on one
                    # search page, after four retries, and the process was
                    # gone. The item path was made survivable and this one was
                    # left fatal, which is exactly the asymmetry a reviewer
                    # flagged and I did not act on.
                    #
                    # Retrying the SAME cursor is safe: the cursor is only
                    # advanced after a page is fully processed, so a retry
                    # re-reads a page rather than skipping one. Nothing is
                    # written and nothing is indexed in between.
                    d = _search_with_retry(src, year, cursor, verbose=verbose)
                    # `page`, not `hits`: the COLUMN named `hits` holds the
                    # answer's SIZE and this holds the records. One word for
                    # both is what made #910 hard to rule out -- the retirement
                    # test `not hits` reads this list while the UPDATE two lines
                    # down writes the count, and they are unrelated quantities.
                    page = d.get("resultList", {}).get("result", [])
                    nxt = d.get("nextCursorMark")
                    # Recorded on every page rather than only when unset: the
                    # count is the size of THIS question's answer, and a slice
                    # re-walked under a wider clause has a larger one. Keeping a
                    # stale value would make seen/hits read as complete while
                    # the slice was only complete for a question nobody asks.
                    #
                    # But ABSENT is not ZERO (#910). A response that parses as
                    # JSON without a `hitCount` -- a degraded or error document
                    # -- used to write 0, and a retired slice carrying hits=0
                    # satisfies `seen >= hits` trivially, so the ledger loses
                    # the only check it has on whether the walk was complete.
                    #
                    # AN EARLIER VERSION OF THIS COMMENT CALLED THAT "the
                    # ledger's readability, not the crawl's correctness". That
                    # was too generous and is retracted. `seen` is an exact
                    # multiple of PAGE for exactly 2 of 122 walked slices, and
                    # they are precisely the two that carried hits=0 AND
                    # done=1 -- so both walks stopped ON A PAGE BOUNDARY, i.e.
                    # on an empty or degraded page rather than at the end of
                    # the data, and were then retired. Against their healthy
                    # neighbours (every other slice has seen ~= hits) MED 2022
                    # holds 138,000 where 2021 and 2023 bracket ~486,000, and
                    # MED 1990 holds 10,000 where 1989 has 64,194: together
                    # roughly 402,000 records missing from slices the ledger
                    # reports as complete.
                    #
                    # The zeroing does not CAUSE that -- the `not page` break
                    # does, and it is legitimate for a genuinely empty page --
                    # but the zeroing is what HID it. Keeping the count honest
                    # turns a silent early retirement into a visible
                    # `seen < hits` on a done slice, which is the detection
                    # signal this ledger did not have. Refusing to retire on
                    # that signal is the prevention and is NOT done here: an
                    # unconditional rewind-on-short is the same shape as the
                    # non-terminating deferral loop fixed in #907, so it needs
                    # the same "optimistic once" bounding. Tracked in #910.
                    reported = d.get("hitCount")
                    if reported is not None:
                        st.execute("UPDATE slice SET hits=? WHERE src=? AND year=?",
                                   (reported, src, year))
                    stopped_early = False
                    examined = 0

                    # PHASE 1 -- SELECT, serially. The dedup lookup is a
                    # SQLite read and the decision to fetch depends on it, so
                    # this stays single-threaded and in order.
                    candidates = []
                    page_keys = set()
                    for rec in page:
                        seen += 1
                        examined += 1
                        pmid = norm_pmid(rec.get("pmid"))
                        pmcid = norm_pmcid(rec.get("pmcid"))
                        doi = norm_doi(rec.get("doi"))
                        if is_held(ident, pmid, pmcid, doi):
                            continue
                        if not (pmid or pmcid or doi or rec.get("id")):
                            continue  # nothing to dedup on later; refuse it
                        # A page can contain the same article twice. The
                        # serial version caught that because each record was
                        # indexed before the next was tested; batching removes
                        # that, so the batch dedups against itself.
                        keys = {("pmid", pmid), ("pmcid", pmcid), ("doi", doi)}
                        keys = {k for k in keys if k[1]}
                        if keys & page_keys:
                            continue
                        page_keys |= keys
                        candidates.append((rec, pmid, pmcid, doi))
                        if limit_new and totals["new"] + len(candidates) >= limit_new:
                            stopped_early = True
                            break

                    # PHASE 2 -- FETCH, concurrently. Only the HTTP GETs run in
                    # parallel; nothing here touches a shard, the index or the
                    # ledger. That is the whole reason this is safe: the
                    # write-then-index ordering that makes the index
                    # trustworthy lives entirely in phase 3, on one thread.
                    wanted = [c for c in candidates
                              if oa_possible and _text_available(c[0])
                              and _text_permitted(c[0])]
                    texts = _fetch_texts(wanted)

                    # PHASE 3 -- WRITE, serially, in page order.
                    for rec, pmid, pmcid, doi in candidates:
                        text = None
                        if id(rec) in texts:
                            got = texts[id(rec)]
                            if isinstance(got, TransientFetchError) and not retry_pass:
                                # SKIP THE RECORD, AND RECORD THAT WE DID.
                                # An earlier version said "skipping leaves it
                                # NEW, so a later pass fetches it properly".
                                # There was no later pass: the cursor advanced
                                # and the slice was marked done=1, so the
                                # record was skipped forever while `seen`
                                # certified the slice as fully walked. The
                                # deferral is persisted on the slice now, the
                                # slice is not retired, and its cursor is
                                # rewound so the next run re-walks it.
                                #
                                # FIRST PASS ONLY. Optimistic is right once:
                                # most of these are transient. Deferring the
                                # same record on every subsequent pass is what
                                # made a slice with one permanently broken
                                # document unretirable.
                                deferred += 1
                                slice_deferred += 1
                                consecutive_deferred += 1
                                if _should_abort(consecutive_deferred):
                                    raise RuntimeError(
                                        f"{consecutive_deferred} consecutive "
                                        f"full-text failures; the service looks "
                                        f"down rather than busy: {got}") from got
                                if verbose and deferred % 25 == 1:
                                    print(f"  ! deferred {deferred} "
                                          f"(latest: {got})", flush=True)
                                # NOT counted as examined: `seen == hits` is
                                # what certifies a slice as fully walked.
                                examined -= 1
                                continue
                            if isinstance(got, TransientFetchError):
                                # RETRY PASS, AND IT FAILED AGAIN. Store the
                                # record with its title, abstract and MeSH and
                                # without full text, counted as examined, so
                                # the slice can retire. The record is NOT lost
                                # -- only its full text is, and that full text
                                # has now failed four in-run retries twice,
                                # across separate runs.
                                #
                                # GUARDED AGAINST A DOWN SERVICE, because
                                # otherwise an outage during a retry pass would
                                # quietly retire a whole slice as abstract-only
                                # and there would be no second chance. A long
                                # run of consecutive failures here is an
                                # outage, not a set of broken documents, and it
                                # raises exactly as the first pass does.
                                slice_abandoned += 1
                                consecutive_abandoned += 1
                                consecutive_deferred = 0
                                if _should_abort(consecutive_abandoned):
                                    raise RuntimeError(
                                        f"{consecutive_abandoned} consecutive "
                                        f"full-text failures on a RETRY pass; "
                                        f"the service looks down rather than "
                                        f"holding broken documents, and "
                                        f"retiring this slice without text "
                                        f"would be permanent: {got}") from got
                                if verbose and slice_abandoned % 25 == 1:
                                    print(f"  ~ no full text after a second run "
                                          f"({slice_abandoned}); storing the "
                                          f"record without it: {got}", flush=True)
                            else:
                                text = got
                                consecutive_deferred = 0
                                consecutive_abandoned = 0
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
                        # WRITE THEN INDEX, never the reverse. write() returns
                        # only once the bytes are flushed, so the index never
                        # claims a record whose bytes could still be lost.
                        shards.write(out)
                        kept += 1
                        ft += bool(text)
                        # THE EPMC ID IS INDEXED TOO, and it has to be: the
                        # admission test above accepts a record on the strength
                        # of `rec["id"]` alone, so a record carrying only that
                        # was written and then indexed under nothing -- and
                        # would be re-downloaded on every future run, forever.
                        # Latent on MED and PPR, which always carry a DOI or a
                        # PMID; the patent, thesis and abstract sources are
                        # where records with none of the three live, and they
                        # have not been crawled yet.
                        for kind, key in (("pmid", pmid), ("pmcid", pmcid),
                                          ("doi", doi), ("epmc", rec.get("id"))):
                            if key:
                                ident.execute(
                                    "INSERT INTO held(kind,key,source,has_fulltext) VALUES (?,?,?,?) "
                                    "ON CONFLICT(kind,key) DO UPDATE SET "
                                    "has_fulltext=MAX(held.has_fulltext, excluded.has_fulltext)",
                                    (kind, key, f"expand-{src}", 1 if text else 0))
                    # A page abandoned part-way by --limit-new must NOT advance
                    # the cursor: the unexamined tail would be skipped forever
                    # on resume, silently, because `seen` had been credited the
                    # whole page. Staying on `cursor` re-reads the page next
                    # time, which costs one request and loses nothing.
                    st.execute("UPDATE slice SET cursor=?, seen=seen+?, kept=kept+?, "
                               "fulltext=fulltext+?, deferred=?, updated=? "
                               "WHERE src=? AND year=?",
                               (cursor if stopped_early else (nxt or cursor),
                                examined, kept, ft, slice_deferred,
                                time.time(), src, year))
                    st.commit()
                    ident.commit()
                    totals["seen"] += examined
                    totals["new"] += kept
                    totals["fulltext"] += ft
                    totals["deferred"] = deferred
                    totals["no_fulltext"] += slice_abandoned - reported_abandoned
                    reported_abandoned = slice_abandoned
                    kept = ft = 0
                    if limit_new and totals["new"] >= limit_new:
                        return _finish(shards, st, ident, totals)
                    if not nxt or nxt == cursor or not page:
                        break
                    cursor = nxt
                    time.sleep(SLEEP)
                if slice_deferred:
                    # NOT done, and rewound. Leaving the cursor where it
                    # stopped would mean the re-walk starts after the records
                    # it needs to recover.
                    st.execute("UPDATE slice SET done=0, cursor='*' "
                               "WHERE src=? AND year=?", (src, year))
                    if verbose:
                        print(f"  {src:4} {year}  NOT done: {slice_deferred:,} "
                              f"deferred, will re-walk", flush=True)
                else:
                    # #913: an empty result list is read as end-of-data by the
                    # break above, and a degraded response that parses with an
                    # empty `resultList` is indistinguishable from a genuinely
                    # exhausted cursor. Two slices retired that way before
                    # anyone noticed -- MED 2022 at 138,000 where its
                    # neighbours bracket ~486,000, and MED 1990 at 10,000
                    # against 1989's 64,194, about 402,000 records absent from
                    # slices the ledger called complete. The fingerprint was
                    # that `seen` landed on an exact multiple of PAGE for
                    # exactly those two of 122 walked slices.
                    #
                    # Now that an absent hitCount no longer zeroes the count,
                    # `seen < hits` IS the incompleteness signal, so SAY SO.
                    # Read back from the ledger rather than using the local
                    # `seen`, which counts only THIS run while the column
                    # accumulates across re-walks -- that also makes the
                    # comparison conservative, since a re-walked slice has an
                    # inflated count and cannot false-positive.
                    #
                    # Deliberately NOT a rewind. An unconditional
                    # rewind-on-short is the same shape as the non-terminating
                    # deferral loop fixed in #907, which burned ~2.9M search
                    # requests per pass on a permanently broken record. The
                    # explicit recovery pass is #913.
                    tally = st.execute(
                        "SELECT seen, hits FROM slice WHERE src=? AND year=?",
                        (src, year)).fetchone()
                    walked, total = (tally or (0, -1))
                    short = total > 0 and walked < total
                    # AND the case the fix above creates. `hits` defaults to -1
                    # and an absent hitCount now leaves it alone, so a slice
                    # whose responses NEVER carried one retires with total = -1
                    # and walked > 0 -- no recorded total to compare against.
                    # `short` is False there, so without this it would be
                    # exactly as silent as the thing being fixed.
                    #
                    # This comment used to read the proof off the live ledger:
                    # "the two damaged slices carry hits=0 rather than -1,
                    # which is itself the proof an UPDATE ran with a zero".
                    # That was true when written and `_migrate` v1 has since
                    # restored -1 on all 139 such rows, so it is recorded here
                    # as history rather than as something a reader can still
                    # check in the data. Both values behave identically from
                    # here on -- `total > 0` is false and `total <= 0` true for
                    # each -- which is why the migration could be safe.
                    no_total = total <= 0 and walked > 0
                    st.execute("UPDATE slice SET done=1, deferred=0 "
                               "WHERE src=? AND year=?", (src, year))
                    if short:
                        print(f"  {src:4} {year}  RETIRED SHORT: walked "
                              f"{walked:,} of {total:,} -- the page walk "
                              f"ended before the recorded total, so this slice "
                              f"is INCOMPLETE and marked done (#913)",
                              flush=True)
                    elif no_total:
                        print(f"  {src:4} {year}  RETIRED WITH NO TOTAL: walked "
                              f"{walked:,} and no response ever reported a "
                              f"hitCount, so completeness is UNKNOWN rather "
                              f"than verified (#913)", flush=True)
                    if verbose and slice_abandoned:
                        print(f"  {src:4} {year}  done, but {slice_abandoned:,} "
                              f"records stored WITHOUT full text after a second "
                              f"failed attempt", flush=True)
                st.commit()
                totals["slices"] += 1
                if verbose and seen:
                    print(f"  {src:4} {year}  seen {seen:>6,}  new {totals['new']:>7,}  "
                          f"text {totals['fulltext']:>7,}", flush=True)
    except BaseException:
        # `_get` raises deliberately on a persistent HTTP failure, and a run
        # this long will also meet Ctrl-C. Either way the shard must be closed
        # and both databases committed: the previous `finally: pass` was dead,
        # so any exception skipped _finish entirely and abandoned the open
        # shard mid-record. Re-raised, because a crawl that stopped early must
        # not look like one that finished.
        _finish(shards, st, ident, totals)
        raise
    return _finish(shards, st, ident, totals)


def _finish(shards, st, ident, totals):
    """Close the shard and commit both ledgers. Safe to call twice.

    The exception path calls this and then re-raises, and the limit-new path
    calls it before returning -- so if the FIRST call fails part-way, the
    handler calls it again on already-closed connections and the resulting
    ProgrammingError replaces the original exception, hiding why the run
    actually stopped. Each step is therefore independently guarded: this
    function's job is to salvage what it can on the way out, not to be the
    thing that reports a problem.
    """
    for step in (shards.close,
                 st.commit, ident.commit,
                 st.close, ident.close):
        try:
            step()
        except Exception:  # noqa: BLE001 - a close failure must not mask the cause
            pass
    return totals


def incomplete_slices(c) -> list:
    """Slices marked done that the ledger cannot certify as complete (#913).

    TWO detectors, because the historical cases are invisible to the general
    one:

      1. `seen < hits` -- the general case, and it only became expressible when
         #910 stopped an absent hitCount from zeroing the total. Retrospectively
         it finds one slice repo-wide (PPR 2022, short by 25 of 23,365), which
         is the right order of noise for a check anyone will actually read.

      2. `hits <= 0` with `seen > 0` -- the slices the zeroing BLINDED, which
         detector 1 cannot see because it needs `hits > 0`. The schema default
         is -1 and a healthy walk always records a real total, so this cannot
         arise on a sound path. Bounded: 3 slices match today (MED 1990, MED
         2015, MED 2022) and the clause can be dropped once they are recovered.

      3. `hits < 0` -- NO TOTAL WAS EVER OBTAINED, whatever was walked. The
         first two detectors BOTH REQUIRE A NON-ZERO COUNT (`hits > 0` or
         `seen > 0`) and so are blind to the worst row the ledger can hold:
         one marked done having walked NOTHING. That is not hypothetical and
         it is the larger half of #913 -- 24 consecutive MED slices (1991-2014)
         read `seen=0, hits=0, done=1`, while the service reports 64,965
         records for 1991 and 274,868 for 2014. #912's own analysis missed them
         because it reasoned over "2 of 122 WALKED slices", a denominator that
         excludes every slice with `seen = 0` by construction.

         ZERO SATISFIED BOTH CHECKS. A detector asking "did the walk fall short
         of the total?" cannot fire when neither number exists, so the damaged
         rows were the ones it could not see. The fix is a clause keyed on the
         ABSENCE of a total rather than on a shortfall between two totals.

         It is keyed on `hits < 0` and NOT on `seen = 0`, which is the whole
         difficulty: a genuinely empty year also walked nothing, and 112 PPR
         slices (1900-2011) are legitimately empty because preprints did not
         exist. Flagging `seen = 0` would reopen those forever. After #912 a
         zero total written by a real response is trustworthy, so an empty
         slice settles at `hits = 0` and stops being flagged, while a slice
         whose total was never obtained keeps the -1 the schema defaults to.
         The one-time `_migrate` v1 restores that sentinel on legacy rows.

         A slice the service PERSISTENTLY refuses a total for stays flagged
         every pass, and that is the honest outcome rather than a loop: this
         function only REPORTS, the rewind needs `--reopen-incomplete`, so an
         operator sees "still cannot be certified" instead of the crawl
         re-walking it forever. That distinction is why the bounding #907
         needed is not needed here.

    Deliberately NOT `seen % PAGE == 0`. That is how the two damaged slices were
    FOUND -- exactly 2 of 122 walked slices landed on a page boundary and they
    were precisely the two carrying hits=0 and done=1 -- but it is unsafe as a
    standing rule: a slice whose true total happens to be a round multiple of
    the page size would be reopened forever.
    """
    return list(c.execute(
        "SELECT src, year, seen, hits FROM slice "
        "WHERE done=1 AND ((hits > 0 AND seen < hits) "
        "               OR (hits <= 0 AND seen > 0) "
        "               OR hits < 0) "
        "ORDER BY src, year"))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--from-year", type=int, default=1900)
    ap.add_argument("--to-year", type=int, default=2026)
    ap.add_argument("--limit-new", type=int, default=None,
                    help="stop after this many new records (a dry run)")
    ap.add_argument("--sources", default=None, help="comma list, e.g. PPR,PAT")
    ap.add_argument("--status", action="store_true")
    ap.add_argument("--incomplete", action="store_true",
                    help="report slices marked done that the ledger cannot "
                         "certify as complete (#913); read-only")
    ap.add_argument("--reopen-incomplete", action="store_true",
                    help="reopen those slices (done=0, cursor='*') so the next "
                         "pass re-walks them; writes to the ledger")
    a = ap.parse_args()
    if a.incomplete or a.reopen_incomplete:
        c = state()
        rows = incomplete_slices(c)
        if not rows:
            print("no slice is marked done without being certifiable complete")
            c.close()
            return 0
        print(f"  {'src':5} {'year':>5} {'walked':>10} {'total':>10}  why")
        for src, year, seen, hits in rows:
            if hits > 0:
                why = "walked short of the recorded total"
            elif seen > 0:
                why = "no total was ever recorded, so completeness is UNKNOWN"
            else:
                why = "NOTHING was walked and no total recorded -- never ran"
            total = f"{hits:,}" if hits > 0 else "-"
            print(f"  {src:5} {year:>5} {seen:>10,} {total:>10}  {why}")
        if not a.reopen_incomplete:
            print(f"\n{len(rows)} slice(s). Re-run with --reopen-incomplete to "
                  f"reopen them for the next pass.")
            c.close()
            return 0
        # WRITES. Separate flag rather than a prompt: this is run from a script
        # as often as by hand, and a slice reopened by accident costs a full
        # re-walk of a source-year.
        for src, year, _seen, _hits in rows:
            c.execute("UPDATE slice SET done=0, cursor='*' WHERE src=? AND year=?",
                      (src, year))
        c.commit()
        print(f"\nreopened {len(rows)} slice(s); the next pass will re-walk them")
        c.close()
        return 0
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
