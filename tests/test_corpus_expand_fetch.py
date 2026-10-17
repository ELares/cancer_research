"""Guards for the expansion fetcher.

Two requirements drive every test here, and both were stated as one sentence
each: never download an article we already have, and do not take paywalled
content. The first is enforced by asking the identity index before fetching
and updating it in the same pass; the second by asking a source that refuses,
rather than by this code deciding what is free.
"""
import ast
import errno
import gzip
import json
import pathlib
import sqlite3
import sys
import urllib.parse
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
import corpus_expand_fetch as fx  # noqa: E402
import corpus_identity_index as ix  # noqa: E402


def test_text_is_asked_for_on_an_open_licence_not_only_the_oa_flag():
    """The flag and the licence are different facts, and the first version
    conflated them.

    `isOpenAccess` means "the full text is in Europe PMC's OA subset". It is
    `N` for about 118,000 cancer records that are nonetheless CC-BY, and every
    one of them was being skipped without being asked for.
    """
    assert fx._text_permitted({"isOpenAccess": "Y"})
    assert fx._text_permitted({"isOpenAccess": "N", "license": "cc by"})
    assert fx._text_permitted({"isOpenAccess": "N", "license": "CC BY-NC-ND"})
    assert fx._text_permitted({"license": "cc0"})
    # And a record with neither is not asked for at all.
    assert not fx._text_permitted({"isOpenAccess": "N", "license": None})
    assert not fx._text_permitted({"isOpenAccess": "N", "license": ""})
    assert not fx._text_permitted({})


def test_availability_and_permission_are_asked_separately():
    """Two different questions, and collapsing them cost a day of 404s.

    Preprints are overwhelmingly CC-BY, so a permission-only gate says yes to
    all 108,000 of them while Europe PMC holds text for barely a sixth. The
    rest answer 404 -- correctly, but only after a request and a sleep each.
    """
    here_and_allowed = {"inEPMC": "Y", "license": "cc by"}
    allowed_not_here = {"inEPMC": "N", "isOpenAccess": "N", "license": "cc by"}
    assert fx._text_available(here_and_allowed) and fx._text_permitted(here_and_allowed)
    assert fx._text_permitted(allowed_not_here), "a CC-BY record is permitted"
    assert not fx._text_available(allowed_not_here), (
        "a record Europe PMC does not hold is still being asked for, which is "
        "a guaranteed 404 per record")
    src = Path(fx.__file__).read_text()
    assert "_text_available(c[0])" in src and "_text_permitted(c[0])" in src, (
        "the fetch is gated on only one of the two questions")


def test_a_restrictive_licence_is_not_read_as_open():
    """The gate must not be a substring free-for-all: 'no licence' or a
    publisher-specific term is not CC."""
    for lic in ("all rights reserved", "copyright acs", "subscription",
                "unknown", "none"):
        assert not fx._text_permitted({"isOpenAccess": "N", "license": lic}), lic


def test_the_source_is_what_refuses_paywalled_text(monkeypatch):
    """This code never decides what is free. It asks, and a 404 is the answer.

    That is the whole reason no publisher page is ever fetched: the endpoint
    that serves open access is the only one asked, and it declines the rest.
    """
    calls = []

    def fake_get(url, tries=4, timeout=120):
        calls.append(url)
        return None  # what _get returns for a 404

    monkeypatch.setattr(fx, "_get", fake_get)
    assert fx.fetch_fulltext({"pmcid": "PMC1", "source": "MED"}) is None
    assert calls and calls[0].endswith("/PMC1/fullTextXML")

    # THE HOST IS PARSED, not searched for. `"ebi.ac.uk" in url` is satisfied by
    # https://evil.example/?ref=ebi.ac.uk and by https://ebi.ac.uk.attacker.test/
    # -- a substring standing in for a host check, which is the same shape as
    # the quote-splitting scan that let a planted publisher fetch through. Here
    # it would have passed a URL pointing anywhere at all, in the one test whose
    # stated job is that text comes only from Europe PMC. CodeQL flagged it as
    # py/incomplete-url-substring-sanitization and was right.
    host = urllib.parse.urlsplit(calls[0]).hostname
    assert host == "www.ebi.ac.uk", (
        f"text is fetched from {host!r}, not Europe PMC")


def _fetch_hosts(src: str) -> set:
    """Every host named by a string literal in the module, however quoted.

    Parsed with `ast` rather than split on `"`. The first version of this guard
    tokenised on the double-quote character alone, so a single-quoted URL was
    invisible to it -- planting

        return _get('https://www.sciencedirect.com/science/article/pii/' + doi)

    left the whole file GREEN. A guard against scraping that a scraper walks
    straight through is worse than no guard, because it is quoted as evidence.
    Walking the AST also sees f-string prefixes, concatenations and multi-line
    strings, none of which the character split could reach.
    """
    hosts = set()
    for node in ast.walk(ast.parse(src)):
        parts = []
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            parts = [node.value]
        elif isinstance(node, ast.JoinedStr):  # f-string: check its literal parts
            parts = [v.value for v in node.values
                     if isinstance(v, ast.Constant) and isinstance(v.value, str)]
        for s in parts:
            if s.startswith(("http://", "https://")):
                # urlsplit, not hand-rolled splitting: userinfo, ports and
                # query strings all move where the host appears to be, and a
                # hand-rolled parser is how a lookalike domain gets read as
                # the real one.
                host = urllib.parse.urlsplit(s).hostname
                if host:
                    hosts.add(host)
    return hosts


def test_only_europe_pmc_is_ever_fetched_from():
    """A publisher-page fetch would be scraping, which MISSION.md rules out.

    EQUALITY, not a subset. `hosts <= {"www.ebi.ac.uk"}` is satisfied by the
    EMPTY set, so any rewrite that assembled URLs from a variable would leave
    the guard green while proving nothing. Requiring the host to actually be
    present means the guard fails if it stops being able to see the URLs.
    """
    hosts = _fetch_hosts(Path(fx.__file__).read_text())
    assert hosts == {"www.ebi.ac.uk"}, (
        f"expected exactly the Europe PMC host, found {sorted(hosts)}")


def test_the_host_guard_actually_catches_a_publisher_fetch():
    """The guard is exercised on a scraping module rather than assumed to work.

    Both quote styles, because the defect this replaces was a quote-style blind
    spot and a fix for one style would look identical to a fix for both.
    """
    for url in ('"https://www.sciencedirect.com/science/article/pii/X"',
                "'https://www.sciencedirect.com/science/article/pii/X'",
                'f"https://link.springer.com/article/{doi}"'):
        found = _fetch_hosts(f"def scrape(doi):\n    return _get({url})\n")
        assert found and found != {"www.ebi.ac.uk"}, f"guard blind to {url}"


EPMC_HOST = "www.ebi.ac.uk"

# Every URL here CONTAINS the Europe PMC host and is served by someone else, so
# the substring form of the check passed all of them. They are built by
# interpolating the host rather than typed out, which makes that containment a
# property of the construction instead of something a second assertion has to
# state -- and that second assertion was itself a substring test against a URL,
# so CodeQL flagged the demonstration exactly as it flagged the original.
# Deriving it is better than suppressing it: the property is now impossible to
# get wrong rather than merely checked.
LOOKALIKE_URLS = (
    f"https://{EPMC_HOST}.attacker.test/x",   # host is a longer domain
    f"https://evil.example/?ref={EPMC_HOST}",  # host appears in the query
    f"https://{EPMC_HOST}@evil.example/x",     # host is userinfo, not the host
    f"https://not{EPMC_HOST}.co/x",            # host is a suffix of another
)


def test_a_lookalike_host_is_not_mistaken_for_europe_pmc():
    """A host test that a lookalike domain passes is not a host test."""
    for url in LOOKALIKE_URLS:
        # No `EPMC_HOST in url` assertion: the f-strings above interpolate the
        # host, so containment is guaranteed by construction. Asserting it
        # would be a substring test against a URL -- the very pattern under
        # test -- and would be flagged for being exactly what it demonstrates.
        assert _fetch_hosts(f'x = "{url}"') != {EPMC_HOST}, url


def test_a_record_with_no_identifier_is_refused():
    """Writing a record we could never recognise again would guarantee we
    re-download it, which is the one thing this fetcher must not do."""
    src = Path(fx.__file__).read_text()
    assert 'if not (pmid or pmcid or doi or rec.get("id")):' in src
    assert "continue  # nothing to dedup on later; refuse it" in src


def test_every_written_record_is_indexed(tmp_path, monkeypatch):
    """BEHAVIOURAL. This was a source grep for `INSERT INTO held` appearing
    between two string anchors, which broke the moment the loop was
    restructured -- and would have passed just as happily on a loop that
    indexed the wrong record. It now runs the loop and asks the index.
    """
    import corpus_identity_index as _ix
    import sqlite3 as _s3

    _drive(tmp_path, monkeypatch, pages=[_page(["A", "B"], None)])
    c = _s3.connect(tmp_path / "id.sqlite")
    held = {k for (k,) in c.execute("SELECT key FROM held WHERE kind='doi'")}
    c.close()
    assert held == {"10.1234/a", "10.1234/b"}, (
        f"records were written but the index holds {held}; an interrupted run "
        "would re-download them")


def test_a_record_is_never_indexed_without_being_written(tmp_path, monkeypatch):
    """The ordering that makes the index trustworthy. A deferred record must
    appear in NEITHER, or the index vouches for bytes that do not exist."""
    import sqlite3 as _s3

    totals, row, stored = _drive(
        tmp_path, monkeypatch, pages=[_page(["A", "B"], None)], fail_ids={"B"})
    assert stored == ["A"]
    c = _s3.connect(tmp_path / "id.sqlite")
    held = {k for (k,) in c.execute("SELECT key FROM held WHERE kind='doi'")}
    c.close()
    assert held == {"10.1234/a"}, (
        f"index holds {held} but only {stored} were written")


def _read_shards(d: Path) -> list:
    """Every record readable from the shards, tolerating an unclosed one.

    APPENDED ONE AT A TIME, not built as a comprehension. A comprehension over
    a stream that raises part-way through discards everything it had already
    produced, so `seen += [... for x in fh]` returned the EMPTY list for a
    shard with five perfectly readable records followed by a missing gzip
    trailer -- and this helper then reported total data loss where there was
    none. It cost a false failure here; in a reporting script it would have
    been a published zero.
    """
    seen = []
    for f in sorted(d.glob("*.jsonl.gz")):
        try:
            with gzip.open(f, "rt") as fh:
                for x in fh:
                    seen.append(json.loads(x)["i"])
        except EOFError:
            pass  # a shard the writer never closed; keep what landed
    return sorted(seen)


def test_shards_rotate_into_complete_files(tmp_path, monkeypatch):
    """Records must survive a rotation.

    `monkeypatch` rather than `fx.SHARD_RECORDS = 3`: the bare assignment
    mutated the shared module for every test that ran afterwards in the same
    session and was never restored.
    """
    monkeypatch.setattr(fx, "SHARD_RECORDS", 3)
    s = fx.Shards(tmp_path, "t")
    for i in range(7):
        s.write({"i": i})
    s.close()
    files = sorted(tmp_path.joinpath("t").glob("*.jsonl.gz"))
    assert len(files) >= 2
    assert _read_shards(tmp_path / "t") == list(range(7)), (
        "records were lost across a rotation")


def test_an_interrupted_run_loses_no_record_that_was_written(tmp_path, monkeypatch):
    """THE FAILURE THAT ACTUALLY HAPPENED, reproduced.

    The old version of this test claimed to check an interrupted run and only
    ever called `close()`, which is the clean path -- so it passed on an
    implementation with no crash-safety at all. Two `pkill`s during development
    then left shards of 392 and 2,831 records against a 4,000-record rotation,
    and every record in the lost tails had already been marked held in the
    identity index: permanently invisible to every later run.

    Here the writer is ABANDONED without close(), the way a killed process
    abandons it. Every record written must still be readable, which is what
    the Z_SYNC_FLUSH in Shards.write buys.
    """
    monkeypatch.setattr(fx, "SHARD_RECORDS", 1000)  # no rotation; one open shard
    s = fx.Shards(tmp_path, "t")
    for i in range(5):
        s.write({"i": i})

    # Read while `s` is still OPEN and unclosed. `del s` does not simulate a
    # kill: CPython finalises the object, GzipFile.__del__ closes the stream
    # and the buffer is flushed after all -- so the deleting version of this
    # test passed with the flush REMOVED, and was verified by mutation to be
    # measuring the finaliser rather than the write path. A killed process runs
    # no finaliser, and this is the state its file is left in.
    assert _read_shards(tmp_path / "t") == list(range(5)), (
        "a record was written and acknowledged but its bytes are not yet on "
        "disk; a kill here would leave the identity index claiming we hold an "
        "article whose bytes never landed")
    s.close()


def test_a_resumed_run_continues_rather_than_restarting(tmp_path, monkeypatch):
    """Shard numbering must resume past what is on disk, or a restart
    overwrites the previous run's output."""
    d = tmp_path / "t"
    d.mkdir()
    (d / "t-00000.jsonl.gz").write_bytes(b"")
    (d / "t-00003.jsonl.gz").write_bytes(b"")
    assert fx.Shards(tmp_path, "t").idx == 4


def test_a_finished_slice_is_not_walked_again():
    src = Path(fx.__file__).read_text()
    assert "if row and row[1]:" in src and "continue" in src
    assert "UPDATE slice SET done=1" in src


def test_a_persistent_http_failure_is_loud():
    """A silent give-up is indistinguishable from a slice with nothing new,
    which is the one confusion this ledger cannot afford."""
    src = Path(fx.__file__).read_text()
    assert "raise RuntimeError" in src, (
        "_get swallows a persistent failure, so a network outage would be "
        "recorded as 'no new records'")


def test_every_refusal_status_is_accepted_as_an_answer(monkeypatch):
    """MISSION.md says a refusal is accepted, and 404 is not the only refusal.

    401/403/410 are the server declining to serve an item -- the same boundary
    this crawl leans on rather than deciding for itself. They used to be
    treated as transient: retried twice, then raised, which stopped the entire
    run on a single article the source simply would not hand over.
    """
    import urllib.error

    for code in (401, 403, 404, 410):
        calls = []

        def boom(url, *a, **k):
            calls.append(url)
            raise urllib.error.HTTPError(url, code, "no", {}, None)

        monkeypatch.setattr(fx.urllib.request, "urlopen", boom)
        assert fx._get("https://www.ebi.ac.uk/x", tries=3, timeout=1) is None, code
        assert len(calls) == 1, f"status {code} was retried; it is not transient"


def test_a_search_refusal_is_not_read_as_an_empty_slice(monkeypatch):
    """THE HAZARD OF WIDENING THE REFUSAL SET, and it is not symmetric.

    For one ARTICLE a 403 is the expected answer. For the SEARCH endpoint it
    means a page of results was not delivered -- and the caller reads an empty
    result as "no more pages" and writes `done=1`, so the whole source-year
    slice is skipped forever with seen=0. That is precisely the confusion
    `_get`'s own docstring says this ledger must never make, and applying the
    item rule at the shared layer introduced it.
    """
    import urllib.error

    calls = []

    def boom(url, *a, **k):
        calls.append(url)
        raise urllib.error.HTTPError(url, 403, "rate limited", {}, None)

    monkeypatch.setattr(fx.urllib.request, "urlopen", boom)
    monkeypatch.setattr(fx.time, "sleep", lambda *_: None)
    try:
        fx.search("MED", 2026, "*")
    except RuntimeError:
        pass
    else:
        raise AssertionError(
            "a 403 on search returned quietly; the caller will record the "
            "slice as done and never walk it again")
    assert len(calls) > 1, "a 403 on search is transient and must be retried"


def test_the_item_and_search_refusal_sets_are_not_the_same():
    """If they are ever collapsed back into one, the slice-loss bug returns."""
    assert 403 in fx.REFUSE_ITEM and 403 not in fx.REFUSE_SEARCH
    assert 404 in fx.REFUSE_ITEM and 404 in fx.REFUSE_SEARCH


def test_a_transient_status_is_still_retried_and_then_loud(monkeypatch):
    """The counterpart: widening the refusal set must not swallow a real
    outage, which would be recorded as 'no new records'."""
    import urllib.error

    calls = []

    def boom(url, *a, **k):
        calls.append(url)
        raise urllib.error.HTTPError(url, 503, "later", {}, None)

    monkeypatch.setattr(fx.urllib.request, "urlopen", boom)
    monkeypatch.setattr(fx.time, "sleep", lambda *_: None)
    try:
        fx._get("https://www.ebi.ac.uk/x", tries=3, timeout=1)
    except RuntimeError:
        pass
    else:
        raise AssertionError("a persistent 503 must raise, not return None")
    assert len(calls) == 3, "a transient status must still be retried"


def test_a_partly_examined_page_does_not_advance_the_cursor():
    """`--limit-new` broke out of a page mid-way and then stored the NEXT
    page's cursor, so the unexamined tail was skipped forever on resume --
    silently, because `seen` had been credited the whole page."""
    src = Path(fx.__file__).read_text()
    # Bounded FORWARD from the start marker. `src.index(end)` searches from the
    # top of the file, so any earlier occurrence of the end marker -- one
    # appeared when the query migration added a commit -- makes the slice run
    # backwards and produce an empty string that satisfies every assertion
    # about it. Third time this shape has bitten in this file.
    i = src.index("stopped_early = False")
    body = src[i:src.index("st.commit()", i)]
    assert body.strip(), "the page-update block could not be located"
    assert "stopped_early = True" in body, "the early exit is not recorded"
    assert "cursor if stopped_early else" in body, (
        "the cursor advances even when the page was abandoned part-way")
    seg = body.split("UPDATE slice SET cursor")[1]
    # POSITIVE and rename-proof. This was spelled `"len(hits)" not in seg`,
    # which #910's rename of that local to `page` silently made
    # UNFALSIFIABLE: the string could no longer appear under any mutation, so
    # the guard passed for the wrong reason while protecting nothing. Assert
    # what must be TRUE -- the seen credit is the per-record counter -- plus
    # the absence of ANY length in the segment, which survives renaming the
    # page variable again.
    assert "examined," in seg, (
        "the seen credit is no longer `examined`, the per-record counter")
    assert "len(" not in seg, (
        "a length is credited to seen, so the whole page counts even after a "
        "partial pass; it must be `examined`")


def test_the_run_closes_its_shard_on_the_way_out():
    """`finally: pass` was dead, so any exception -- including the RuntimeError
    _get raises by design -- skipped _finish and abandoned an open shard."""
    src = Path(fx.__file__).read_text()
    assert "finally:\n        pass" not in src, "the dead finally block is back"
    assert "except BaseException:" in src and "_finish(shards, st, ident, totals)\n        raise" in src


def test_finish_actually_closes_and_commits():
    """Emptying `_finish` to `return totals` was a SURVIVING mutation: nothing
    anywhere asserted it does its job, only that calling it twice is safe. A
    salvage function that salvages nothing passes a safety test perfectly."""
    import sqlite3 as s3

    calls = []

    class _S:
        def close(self): calls.append("shard")

    class _DB:
        def commit(self): calls.append("commit")
        def close(self): calls.append("close")

    fx._finish(_S(), _DB(), _DB(), {"seen": 0})
    assert calls.count("shard") == 1, "the open shard was not closed"
    assert calls.count("commit") == 2, "both ledgers must be committed"
    assert calls.count("close") == 2, "both connections must be closed"


def test_finish_is_safe_to_call_twice():
    """The exception path calls it and re-raises; a second call on closed
    connections must not replace the original exception with its own."""
    import sqlite3 as s3

    class _S:
        def close(self): pass

    a, b = s3.connect(":memory:"), s3.connect(":memory:")
    fx._finish(_S(), a, b, {"seen": 0})
    fx._finish(_S(), a, b, {"seen": 0})  # must not raise


def test_one_unavailable_article_does_not_stop_the_crawl(monkeypatch):
    """THE FAILURE THAT ACTUALLY HAPPENED, and it cost a running crawl.

    Every give-up raised the same RuntimeError, so a single 503 on one
    full-text URL propagated out of `run()` and ended the whole thing. The
    article served a 154 KB document on the next request. Transient per-item
    failures are their own type now, and the record loop skips the record
    instead of dying.
    """
    import urllib.error

    def boom(url, *a, **k):
        raise urllib.error.HTTPError(url, 503, "busy", {}, None)

    monkeypatch.setattr(fx.urllib.request, "urlopen", boom)
    monkeypatch.setattr(fx.time, "sleep", lambda *_: None)
    try:
        fx._get("https://www.ebi.ac.uk/x", tries=2, timeout=1)
    except fx.TransientFetchError:
        pass
    else:
        raise AssertionError("a 503 must raise TransientFetchError")

    # DISTINCT from RuntimeError, not merely a subclass of it. Aliasing
    # `TransientFetchError = RuntimeError` passed the subclass check (every
    # class is a subclass of itself) while making every fatal error catchable
    # as transient -- a surviving mutation. The type must be its own.
    assert fx.TransientFetchError is not RuntimeError
    assert issubclass(fx.TransientFetchError, RuntimeError)
    assert not isinstance(RuntimeError("x"), fx.TransientFetchError)
    src = Path(fx.__file__).read_text()
    assert "except TransientFetchError as e:" in src, (
        "the record loop does not catch a per-article failure, so one "
        "unavailable article still ends the crawl")


def test_the_server_is_obeyed_when_it_says_how_long_to_wait(monkeypatch):
    """Guessing a backoff while being told the answer is how a rate limiter
    becomes an outage."""
    import urllib.error

    class H(dict):
        def get(self, k, d=None): return "7" if k == "Retry-After" else d

    e = urllib.error.HTTPError("u", 503, "busy", H(), None)
    assert fx._retry_after(e) == 7.0
    assert fx._retry_after(urllib.error.HTTPError("u", 503, "b", {}, None)) is None
    # Absurd values are clamped rather than trusted, and the cap is 60s not
    # 300s: a five-minute wait on ONE article is worse than deferring it and
    # moving on, since the deferral is now recorded and re-walked.
    class H2(dict):
        def get(self, k, d=None): return "99999" if k == "Retry-After" else d
    assert fx._retry_after(urllib.error.HTTPError("u", 503, "b", H2(), None)) == \
        fx.RETRY_AFTER_CAP == 60.0


def test_get_actually_waits_the_time_the_server_asked_for(monkeypatch):
    """`_retry_after` being correct is not the same as `_get` USING it --
    replacing the call with `wait = None` left every test green."""
    import urllib.error

    class H(dict):
        def get(self, k, d=None): return "9" if k == "Retry-After" else d

    slept = []
    monkeypatch.setattr(fx.time, "sleep", lambda s: slept.append(s))

    def boom(url, *a, **k):
        raise urllib.error.HTTPError(url, 503, "busy", H(), None)

    monkeypatch.setattr(fx.urllib.request, "urlopen", boom)
    try:
        fx._get("https://www.ebi.ac.uk/x", tries=3, timeout=1)
    except fx.TransientFetchError:
        pass
    assert slept and all(s == 9.0 for s in slept), (
        f"server asked for 9s between tries; slept {slept}")


def test_a_zero_or_negative_retry_after_is_floored_not_obeyed(monkeypatch):
    """`Retry-After: 0` is an instruction to hammer. `max(0.0, ...)` obeyed it
    literally and issued four back-to-back requests -- worse than the fixed
    backoff it was overriding."""
    import urllib.error

    for raw in ("0", "-5", "0.0"):
        class H(dict):
            def get(self, k, d=None): return raw if k == "Retry-After" else d
        got = fx._retry_after(urllib.error.HTTPError("u", 503, "b", H(), None))
        assert got >= 1.0, f"Retry-After {raw!r} produced a {got}s wait"


def test_an_http_date_retry_after_is_understood(monkeypatch):
    """RFC 7231 allows delay-seconds OR an HTTP-date, and the date form is what
    CDNs in front of an origin send. `float()` raised on it and the server's
    instruction was silently discarded."""
    import datetime
    import email.utils
    import urllib.error

    when = email.utils.format_datetime(
        datetime.datetime.now(datetime.timezone.utc) + datetime.timedelta(seconds=40))

    class H(dict):
        def get(self, k, d=None): return when if k == "Retry-After" else d

    got = fx._retry_after(urllib.error.HTTPError("u", 503, "b", H(), None))
    assert got is not None and 30 <= got <= 50, got


def test_backoff_gives_more_room_than_the_schedule_it_replaced(monkeypatch):
    """MEASURED against the old schedule, not asserted to be better.

    The first version claimed exponential beat the linear 1.5*(n+1) it
    replaced, and the arithmetic said otherwise: 1+2+4 = 7.0s against
    1.5+3+4.5 = 9.0s, with every individual wait shorter too. Starting the
    exponent at 1 gives 2+4+8 = 14.0s, which actually is more room.
    """
    import urllib.error

    slept = []
    monkeypatch.setattr(fx.time, "sleep", lambda s: slept.append(s))
    monkeypatch.setattr(fx.urllib.request, "urlopen",
                        lambda url, *a, **k: (_ for _ in ()).throw(
                            urllib.error.HTTPError(url, 503, "busy", {}, None)))
    try:
        fx._get("https://www.ebi.ac.uk/x", tries=4, timeout=1)
    except fx.TransientFetchError:
        pass
    assert slept == [2.0, 4.0, 8.0], slept
    old_linear = [1.5 * (n + 1) for n in range(3)]
    assert sum(slept) > sum(old_linear), (
        f"new schedule {slept} (total {sum(slept)}s) gives the server LESS "
        f"room than the linear one it replaced ({old_linear}, "
        f"total {sum(old_linear)}s)")
    assert all(a > b for a, b in zip(slept, old_linear)), (
        "some individual waits got shorter")


# --- End-to-end: what a deferral does to the LEDGER and the SHARDS ----------
#
# Every other test in this file about deferrals reads the source text. None
# drove `run()`, which is why all of them passed while a deferred record was
# being lost forever: the cursor advanced past its page and the slice was then
# marked done=1, so no later pass could ever reach it. Source greps cannot see
# a cursor. This drives the real loop against a fake endpoint.

def _drive(tmp_path, monkeypatch, pages, fail_ids=(), sources=(("MED", True),)):
    """Run the real `run()` over canned search pages, isolated from live data."""
    import corpus_identity_index as _ix

    monkeypatch.setattr(fx, "STATE_DB", tmp_path / "state.sqlite")
    monkeypatch.setattr(fx, "OUT_ROOT", tmp_path / "out")
    monkeypatch.setattr(_ix, "DB", tmp_path / "id.sqlite")
    monkeypatch.setattr(fx, "SLEEP", 0)
    monkeypatch.setattr(fx.time, "sleep", lambda *_: None)

    seq = list(pages)

    def fake_search(src, year, cursor):
        return seq.pop(0) if seq else {"resultList": {"result": []}}

    def fake_fulltext(rec):
        if rec.get("id") in fail_ids:
            raise fx.TransientFetchError(f"503 on {rec['id']}")
        return f"TEXT-{rec['id']}"

    monkeypatch.setattr(fx, "search", fake_search)
    monkeypatch.setattr(fx, "fetch_fulltext", fake_fulltext)
    totals = fx.run([2026], list(sources), verbose=False)

    import sqlite3 as _s3
    c = _s3.connect(tmp_path / "state.sqlite")
    row = c.execute("SELECT cursor, seen, kept, fulltext, deferred, done, "
                    "query_id, hits FROM slice").fetchone()
    c.close()

    import glob as _g
    import gzip as _gz
    import json as _j
    stored = []
    for f in sorted(_g.glob(str(tmp_path / "out" / "**" / "*.jsonl.gz"), recursive=True)):
        try:
            with _gz.open(f, "rt") as fh:
                for ln in fh:
                    stored.append(_j.loads(ln)["epmc_id"])
        except EOFError:
            pass
    return totals, dict(zip(("cursor", "seen", "kept", "fulltext",
                             "deferred", "done", "query_id", "hits"),
                            row)), stored


def _page(ids, nxt):
    # Preserve identity across pages without Python's randomized, truncated
    # hashes: a collision correctly deduplicates two otherwise distinct test
    # records. The nonzero prefix makes this UTF-8 encoding injective.
    return {"resultList": {"result": [
        {"id": i, "pmid": None,
         "pmcid": f"PMC{int.from_bytes(('record:' + i).encode('utf-8'), 'big')}",
         "doi": f"10.1234/{i}",
         "inEPMC": "Y", "isOpenAccess": "Y", "license": "cc by",
         "title": i, "pubYear": "2026"} for i in ids]},
        "nextCursorMark": nxt, "hitCount": len(ids)}


def test_a_slice_that_deferred_a_record_is_not_marked_done(tmp_path, monkeypatch):
    """THE DEFECT, end to end.

    Reproduced before the fix: pass 1 deferred one record, advanced the cursor
    and set done=1; pass 2 skipped the slice entirely and the record was gone
    permanently, while the ledger read seen=3 kept=2 done=1 -- a healthy-looking
    slice over lost content.
    """
    totals, slice_row, stored = _drive(
        tmp_path, monkeypatch,
        pages=[_page(["R0", "R1", "R2"], "PAGE2")],
        fail_ids={"R1"})

    assert stored == ["R0", "R2"], stored
    assert totals["deferred"] == 1
    assert slice_row["deferred"] == 1, "the loss is not recorded anywhere"
    assert slice_row["done"] == 0, (
        "a slice that walked past a record it could not store was retired; "
        "no later pass can ever reach that record")
    assert slice_row["cursor"] == "*", (
        "the cursor was left past the deferred record, so the re-walk starts "
        "after the thing it needs to recover")
    assert slice_row["seen"] == 2, (
        f"seen={slice_row['seen']} counts the deferred record as examined, so "
        "seen==hits would certify a slice that dropped content")


def test_a_clean_slice_is_still_retired(tmp_path, monkeypatch):
    """The refusal must not block the normal case: with nothing deferred the
    slice completes, or the crawl re-walks everything forever."""
    totals, slice_row, stored = _drive(
        tmp_path, monkeypatch,
        pages=[_page(["A", "B"], None)])
    assert stored == ["A", "B"]
    assert (slice_row["deferred"], slice_row["done"]) == (0, 1)


def test_a_partial_outage_does_not_silently_lose_half_a_slice(tmp_path, monkeypatch):
    """`_should_abort` only fires on CONSECUTIVE failures, so an interleaved
    50% failure rate never trips it. Before the fix that completed with exit 0,
    slice done=1, and half the records gone."""
    ids = [f"R{i}" for i in range(20)]
    totals, slice_row, stored = _drive(
        tmp_path, monkeypatch,
        pages=[_page(ids, None)],
        fail_ids={i for n, i in enumerate(ids) if n % 2})
    assert len(stored) == 10
    assert slice_row["deferred"] == 10
    assert slice_row["done"] == 0, (
        "half the slice was lost and the slice was still retired")


def test_a_failing_fetch_is_still_rate_limited(monkeypatch):
    """`continue` used to jump over the politeness pause, so the crawler
    removed its own rate limit exactly when the server asked for room.

    The pause lives in the shared limiter now, so the property is that a
    record whose fetch FAILS still passes through it -- an outage must not
    become a burst.
    """
    calls = []

    class Spy(fx._RateLimit):
        def wait(self):
            calls.append(1)
            super().wait()

    monkeypatch.setattr(fx, "FETCH_WORKERS", 1)
    monkeypatch.setattr(fx, "fetch_fulltext",
                        lambda rec: (_ for _ in ()).throw(
                            fx.TransientFetchError("503")))
    recs = [({"id": f"R{i}"}, None, None, None) for i in range(6)]
    got = fx._fetch_texts(recs, limiter=Spy(0))
    assert len(got) == 6 and all(isinstance(v, fx.TransientFetchError)
                                 for v in got.values())
    assert len(calls) == 6, (
        f"only {len(calls)} of 6 failing fetches were rate limited; an outage "
        "turns into a burst")


def test_a_sustained_outage_aborts_the_run_end_to_end(tmp_path, monkeypatch):
    """Drives `run()` rather than unit-testing the predicate beside it.

    Three surviving mutations hid here: `consecutive_deferred += 0` (the
    counter never moves, so the abort is dead), moving the reset out of the
    open-access block (it fires after every record, so a run never
    accumulates), and a threshold too loose to reach. None is visible to a
    test of `_should_abort` in isolation, because each breaks the WIRING
    rather than the predicate.
    """
    n = fx.MAX_CONSECUTIVE_DEFERRALS + 5
    ids = [f"R{i}" for i in range(n)]
    try:
        _drive(tmp_path, monkeypatch, pages=[_page(ids, None)], fail_ids=set(ids))
    except RuntimeError as e:
        assert "consecutive" in str(e).lower(), e
    else:
        raise AssertionError(
            f"{n} consecutive full-text failures did not stop the run; the "
            "crawl would walk the whole corpus storing nothing")


def test_a_run_of_failures_below_the_threshold_does_not_abort(tmp_path, monkeypatch):
    """The other side of the boundary: transient failures must be survivable,
    which is the entire point of the change."""
    n = fx.MAX_CONSECUTIVE_DEFERRALS - 1
    ids = [f"R{i}" for i in range(n)]
    totals, row, stored = _drive(
        tmp_path, monkeypatch, pages=[_page(ids, None)], fail_ids=set(ids))
    assert totals["deferred"] == n and row["done"] == 0


def test_the_abort_threshold_is_a_reachable_number():
    """A threshold nobody can reach is not a guard. 100 consecutive failures
    at ~1s each is under two minutes of a dead service; 500 would be eight."""
    assert 10 <= fx.MAX_CONSECUTIVE_DEFERRALS <= 200, fx.MAX_CONSECUTIVE_DEFERRALS


def test_a_network_error_is_never_a_silent_give_up(monkeypatch):
    """`_get` returning None means THE SERVER REFUSED. A timeout or DNS
    failure returning None instead would make the caller store the record with
    no text and index it as held -- permanently, and looking exactly like a
    paywalled article. Every existing give-up test used HTTPError, so this
    whole branch was untested."""
    import socket

    for exc in (socket.timeout("timed out"),
                OSError("dns"),
                ConnectionResetError("reset")):
        monkeypatch.setattr(fx.time, "sleep", lambda *_: None)
        monkeypatch.setattr(fx.urllib.request, "urlopen",
                            lambda url, *a, **k: (_ for _ in ()).throw(exc))
        try:
            got = fx._get("https://www.ebi.ac.uk/x", tries=2, timeout=1)
        except fx.TransientFetchError:
            continue
        raise AssertionError(
            f"{type(exc).__name__} returned {got!r} instead of raising; a "
            "network failure is being recorded as 'the server said no'")


def test_full_text_is_retried_more_than_twice():
    """The crawl died because one 503 exhausted two attempts. The retry count
    is part of the fix, so a silent revert to 2 must fail."""
    src = Path(fx.__file__).read_text()
    body = src[src.index("def fetch_fulltext("):src.index("class Shards")]
    tries = [int(t) for t in __import__("re").findall(r"tries=(\d+)", body)]
    assert tries, "fetch_fulltext no longer passes an explicit retry count"
    assert all(t >= 4 for t in tries), f"full text retried only {tries} times"


def test_the_failure_counter_survives_records_that_need_no_full_text(tmp_path, monkeypatch):
    """A run is a MIX, and the counter must only be reset by a SUCCESS.

    Moving the reset out of the open-access block survived every other test,
    because those runs are uniformly open-access and a deferred record
    `continue`s before reaching it. Interleave records that need no text --
    which is most of a real crawl -- and the reset fires on each one, zeroing
    the counter forever: the service can be down for the entire run and the
    abort never trips.
    """
    n = fx.MAX_CONSECUTIVE_DEFERRALS + 5
    recs = []
    for i in range(n):
        recs.append({"id": f"OA{i}", "pmid": None, "pmcid": f"PMC{i}",
                     "doi": None, "inEPMC": "Y", "isOpenAccess": "Y",
                     "license": "cc by", "title": "t", "pubYear": "2026"})
        # A record with no text to fetch at all: not open access, no licence.
        recs.append({"id": f"NO{i}", "pmid": None, "pmcid": None,
                     "doi": f"10.1/{i}", "inEPMC": "N", "isOpenAccess": "N",
                     "license": None, "title": "t", "pubYear": "2026"})

    page = {"resultList": {"result": recs}, "nextCursorMark": None,
            "hitCount": len(recs)}
    try:
        _drive(tmp_path, monkeypatch, pages=[page],
               fail_ids={f"OA{i}" for i in range(n)})
    except RuntimeError as e:
        assert "consecutive" in str(e).lower(), e
    else:
        raise AssertionError(
            "the service failed on every full text for the whole run and the "
            "crawl never aborted, because records needing no text reset the "
            "failure counter")


def test_a_record_with_only_a_europe_pmc_id_is_still_dedupable(tmp_path, monkeypatch):
    """The admission test accepts a record on `rec["id"]` alone, so the index
    must be able to hold it by that key -- otherwise it is written, indexed
    under nothing, and re-downloaded on every future run forever.

    Latent on MEDLINE and preprints, which always carry a DOI or a PMID. The
    patent, thesis and abstract sources are where records with none of the
    three live.
    """
    import sqlite3 as _s3

    page = {"resultList": {"result": [
        {"id": "PPR999", "pmid": None, "pmcid": None, "doi": None,
         "inEPMC": "N", "isOpenAccess": "N", "license": None,
         "title": "t", "pubYear": "2026"}]},
        "nextCursorMark": None, "hitCount": 1}
    totals, row, stored = _drive(tmp_path, monkeypatch, pages=[page])
    assert stored == ["PPR999"], stored

    c = _s3.connect(tmp_path / "id.sqlite")
    held = c.execute("SELECT kind, key FROM held").fetchall()
    c.close()
    assert ("epmc", "PPR999") in held, (
        f"the record was stored but indexed under {held}; the next run would "
        "download it again")


def test_a_page_containing_the_same_article_twice_stores_it_once(tmp_path, monkeypatch):
    """Batching the fetches removed the guarantee that made this work.

    The serial loop indexed each record before testing the next, so a repeat
    inside one page was caught by the normal dedup. With selection and writing
    split apart, the batch has to dedup against itself.
    """
    rec = {"id": "X", "pmid": None, "pmcid": None, "doi": "10.1/dup",
           "inEPMC": "N", "isOpenAccess": "N", "license": None,
           "title": "t", "pubYear": "2026"}
    page = {"resultList": {"result": [dict(rec), dict(rec)]},
            "nextCursorMark": None, "hitCount": 2}
    totals, row, stored = _drive(tmp_path, monkeypatch, pages=[page])
    assert stored == ["X"], f"the same article was stored twice: {stored}"
    assert totals["new"] == 1


def test_adding_workers_does_not_raise_the_request_rate(monkeypatch):
    """CONCURRENCY MUST HIDE LATENCY, NOT DELETE THE RATE LIMIT.

    The serial path slept `SLEEP` after each fetch and the first pooled version
    slept not at all, so adding workers quietly multiplied the request rate as
    well as overlapping the waits -- a measured 6.5x on four workers that was
    partly concurrency and partly the politeness pause going away. That is the
    same mistake as shortening the pause, made larger and harder to see, and
    shortening the pause preceded a 503 that ended a run.
    """
    import time as _t

    starts = []
    monkeypatch.setattr(fx, "SLEEP", 0.02)
    monkeypatch.setattr(fx, "fetch_fulltext",
                        lambda rec: (starts.append(_t.monotonic()), "x")[1])
    recs = [({"id": f"R{i}"}, None, None, None) for i in range(12)]

    rates = {}
    for w in (1, 4, 8):
        starts.clear()
        monkeypatch.setattr(fx, "FETCH_WORKERS", w)
        fx._fetch_texts(recs)
        span = max(starts) - min(starts)
        rates[w] = len(starts) / span if span > 0 else float("inf")

    ceiling = 1.0 / 0.02
    for w, r in rates.items():
        assert r <= ceiling * 1.6, (
            f"{w} workers issued {r:.0f} requests/s against a {ceiling:.0f}/s "
            f"ceiling; the pool is a rate multiplier, not a latency hider")


def test_the_rate_limiter_spaces_requests_from_every_worker(monkeypatch):
    """Admissions are spaced across workers; admitted requests may overlap."""
    import time as _t

    lim = fx._RateLimit(0.05)
    stamps = []

    def hit():
        lim.wait()
        stamps.append(_t.monotonic())

    threads = [__import__("threading").Thread(target=hit) for _ in range(6)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    stamps.sort()
    gaps = [b - a for a, b in zip(stamps, stamps[1:])]
    assert all(g >= 0.04 for g in gaps), f"requests bunched: {gaps}"


def test_a_zero_interval_limiter_does_not_deadlock():
    """SLEEP can legitimately be 0 in tests and dry runs."""
    lim = fx._RateLimit(0)
    for _ in range(5):
        lim.wait()


def test_a_delayed_wakeup_does_not_accumulate_request_credit(monkeypatch):
    """A scheduler pause must not turn expired reservations into a burst."""
    from types import SimpleNamespace

    clock = {"now": 0.0}

    def sleep(delay):
        clock["now"] += delay + 0.2  # the scheduler wakes us up late

    monkeypatch.setattr(fx, "time", SimpleNamespace(
        monotonic=lambda: clock["now"], sleep=sleep))
    limiter = fx._RateLimit(0.05)
    starts = []
    for _ in range(3):
        limiter.wait()
        starts.append(clock["now"])
    assert all(b - a >= 0.05 for a, b in zip(starts, starts[1:])), starts


def test_an_early_wakeup_does_not_release_a_request_early(monkeypatch):
    """The monotonic deadline, rather than sleep's return, grants a slot."""
    from types import SimpleNamespace

    clock = {"now": 0.0, "sleeps": 0}

    def sleep(delay):
        clock["sleeps"] += 1
        clock["now"] += delay / 2 if clock["sleeps"] == 1 else delay

    monkeypatch.setattr(fx, "time", SimpleNamespace(
        monotonic=lambda: clock["now"], sleep=sleep))
    limiter = fx._RateLimit(0.05)
    limiter.wait()
    limiter.wait()
    assert clock["now"] >= 0.05


def test_a_briefly_unavailable_search_page_is_retried_not_fatal(monkeypatch):
    """THE FAILURE THAT KILLED A 302,973-RECORD CRAWL.

    `search` raises rather than returning an empty page, which is right -- an
    empty page reads as "no more results" and would retire the slice with
    seen=0. But letting it raise all the way out of `run()` is a different
    mistake: one 503 on one search page ended a two-day run. The item path was
    made survivable and this one was left fatal, which is the asymmetry a
    reviewer flagged and I did not act on.
    """
    calls = {"n": 0}

    def flaky(src, year, cursor):
        calls["n"] += 1
        if calls["n"] < 3:
            raise fx.TransientFetchError("503")
        return {"resultList": {"result": []}, "hitCount": 0}

    monkeypatch.setattr(fx, "search", flaky)
    monkeypatch.setattr(fx.time, "sleep", lambda *_: None)
    d = fx._search_with_retry("MED", 2026, "*", verbose=False)
    assert d == {"resultList": {"result": []}, "hitCount": 0}
    assert calls["n"] == 3, "the page was not retried"


def test_a_sustained_search_outage_still_stops_the_run(monkeypatch):
    """The other direction: a service that is genuinely gone must stop the
    crawl, because the ledger must never record an outage as an empty slice."""
    monkeypatch.setattr(fx, "search",
                        lambda *a: (_ for _ in ()).throw(fx.TransientFetchError("503")))
    monkeypatch.setattr(fx.time, "sleep", lambda *_: None)
    try:
        fx._search_with_retry("MED", 2026, "*", verbose=False)
    except RuntimeError as e:
        assert "looks down" in str(e)
        assert not isinstance(e, fx.TransientFetchError), (
            "a sustained search outage must not be catchable as a transient "
            "per-item failure")
        return
    raise AssertionError("a permanently failing search returned quietly")


def test_a_retried_search_page_does_not_advance_the_cursor(tmp_path, monkeypatch):
    """Retrying the SAME cursor is what makes this safe: a retry re-reads a
    page rather than skipping one. If the retry advanced, a 503 would silently
    drop up to a thousand records."""
    src_text = Path(fx.__file__).read_text()
    i = src_text.index("d = _search_with_retry(")
    j = src_text.index("UPDATE slice SET cursor=?", i)
    between = src_text[i:j]
    assert "cursor = nxt" not in between, (
        "the cursor is advanced between the search and the page being "
        "processed, so a retry would skip records")

    calls = []

    def flaky(src, year, cursor):
        calls.append(cursor)
        if len(calls) < 3:
            raise fx.TransientFetchError("503")
        return {"resultList": {"result": []}, "hitCount": 0}

    monkeypatch.setattr(fx, "search", flaky)
    monkeypatch.setattr(fx.time, "sleep", lambda *_: None)
    fx._search_with_retry("MED", 2026, "CURSOR-A", verbose=False)
    assert calls == ["CURSOR-A"] * 3, f"retries moved the cursor: {calls}"


def test_the_search_retry_schedule_covers_a_real_outage():
    """A schedule too short to outlast a restart is decoration. `_get` already
    spends ~14s of its own retries per attempt."""
    assert fx.SEARCH_PAGE_ATTEMPTS >= 3
    assert sum(fx.SEARCH_PAGE_BACKOFF) >= 120, (
        f"total backoff {sum(fx.SEARCH_PAGE_BACKOFF)}s will not outlast a "
        "brief service interruption")


# --- What the query can ever see ------------------------------------------
#
# The crawler asked for the literal word `cancer`. Measured against the term
# set, that sees 5,575,342 records where the set sees 9,872,234 -- and 593,583
# records say `carcinoma` while never saying `cancer`. No amount of paging,
# concurrency or retrying reaches a record the query never asks for, which
# makes this the one defect in this crawler that all the others could not have
# compensated for.

def test_the_subject_clause_covers_the_projects_own_cancer_definition():
    """GROUNDED, not invented. The census defines cancer as the MeSH C04 tree,
    and this repository commits those 704 descriptor names. A term set that
    replaces a one-word query has to cover that definition, and the coverage is
    MEASURED here so it cannot drift downward unnoticed.
    """
    import re
    tsv = REPO / "corpus" / "atlas" / "mesh" / "c04-descriptors.tsv"
    if not tsv.exists():
        import pytest
        pytest.skip("C04 definition not present in this checkout")
    names = []
    for ln in tsv.read_text(encoding="utf-8").split("\n"):
        if ln.strip() and not ln.startswith("#"):
            parts = ln.split("\t")
            if len(parts) >= 2:
                names.append(parts[1].strip().lower())
    assert len(names) > 600, f"only {len(names)} descriptors read"

    stems = tuple(t.lower() for t in fx.CANCER_TERMS) + ("oma",)
    covered = sum(1 for n in names if any(s in n for s in stems))
    share = covered / len(names)
    assert share >= 0.75, (
        f"the subject clause covers only {share:.1%} of the {len(names)} C04 "
        "descriptor names; the crawl cannot see what it does not ask for")


def test_the_clause_is_wider_than_the_single_word_it_replaced():
    """The specific regression that matters: quietly reverting to `cancer`
    alone would restore a 44% blind spot and change nothing else observable."""
    assert len(fx.CANCER_TERMS) > 10
    for essential in ("carcinoma", "leukemia", "lymphoma", "melanoma",
                      "sarcoma", "neoplasm", "glioma"):
        assert essential in fx.CANCER_TERMS, (
            f"{essential!r} is absent; papers using that word and never "
            "'cancer' are invisible to this crawl")
    clause = fx.cancer_clause()
    assert clause.startswith("(") and clause.endswith(")")
    assert " OR " in clause


def test_the_search_query_uses_the_clause_not_a_bare_word(monkeypatch):
    seen = {}
    monkeypatch.setattr(fx, "_get", lambda url, **k: seen.setdefault("u", url) and None)
    try:
        fx.search("MED", 2026, "*")
    except Exception:
        pass
    import urllib.parse
    q = urllib.parse.unquote(seen.get("u", ""))
    assert "carcinoma" in q and "leukemia" in q, f"query is still narrow: {q[:200]}"


def test_a_slice_walked_with_an_older_question_is_reopened(tmp_path, monkeypatch):
    """A slice marked done under a NARROWER clause is not done under a wider
    one -- it was fully walked for a question nobody asks any more. Without
    this, widening the query fixes future slices and leaves every completed one
    holding its old, smaller answer forever.
    """
    monkeypatch.setattr(fx, "STATE_DB", tmp_path / "state.sqlite")
    st = fx.state()
    st.execute("INSERT INTO slice(src,year,cursor,seen,kept,done,query_id) "
               "VALUES ('MED',2026,'ABC',1000,900,1,'oldfingerprint')")
    st.commit()

    n = fx.reopen_stale_slices(st, verbose=False)
    assert n == 1, "the stale slice was not re-opened"
    row = st.execute("SELECT cursor, done, query_id FROM slice").fetchone()
    assert row[0] == "*", "the cursor was left past records the wider query adds"
    assert row[1] == 0, "a slice walked for a different question is still marked done"
    assert row[2] == fx.query_fingerprint()

    # And a second call must be a no-op, or every run re-walks everything.
    assert fx.reopen_stale_slices(st, verbose=False) == 0
    st.close()


def test_the_fingerprint_changes_with_the_terms(monkeypatch):
    """If it did not, a widened clause would look identical to the old one and
    no slice would ever be re-opened."""
    a = fx.query_fingerprint()
    monkeypatch.setattr(fx, "CANCER_TERMS", fx.CANCER_TERMS + ("chondrosarcoma",))
    assert fx.query_fingerprint() != a
    # ...and it must not depend on the ORDER the terms happen to be written in.
    monkeypatch.setattr(fx, "CANCER_TERMS", tuple(reversed(fx.CANCER_TERMS)))
    monkeypatch.setattr(fx, "CANCER_TERMS", tuple(sorted(fx.CANCER_TERMS)))
    b = fx.query_fingerprint()
    monkeypatch.setattr(fx, "CANCER_TERMS", tuple(reversed(sorted(fx.CANCER_TERMS))))
    assert fx.query_fingerprint() == b, "reordering the terms re-walks the corpus"


def test_run_reopens_stale_slices_before_walking(tmp_path, monkeypatch):
    """Calling `reopen_stale_slices` directly proves the function works, not
    that `run()` calls it -- removing the call left every test green while a
    completed slice kept its old, narrower answer forever."""
    monkeypatch.setattr(fx, "STATE_DB", tmp_path / "state.sqlite")
    st = fx.state()
    st.execute("INSERT INTO slice(src,year,cursor,seen,kept,done,query_id) "
               "VALUES ('MED',2026,'OLDCURSOR',1000,900,1,'oldfingerprint')")
    st.commit()
    st.close()

    totals, row, stored = _drive(
        tmp_path, monkeypatch, pages=[_page(["N1", "N2"], None)])

    assert stored == ["N1", "N2"], (
        "the slice was skipped as done, so what the wider query adds is never "
        "fetched")
    assert row["query_id"] == fx.query_fingerprint()


def test_a_reopened_slice_records_the_new_questions_total(tmp_path, monkeypatch):
    """`hits` is the size of THIS question's answer. Keeping the old value
    makes seen/hits read as complete while the slice is complete only for a
    question nobody asks."""
    monkeypatch.setattr(fx, "STATE_DB", tmp_path / "state.sqlite")
    st = fx.state()
    st.execute("INSERT INTO slice(src,year,cursor,seen,kept,hits,done,query_id) "
               "VALUES ('MED',2026,'*',0,0,111,0,'oldfingerprint')")
    st.commit()
    st.close()

    page = _page(["A"], None)
    page["hitCount"] = 999
    totals, row, stored = _drive(tmp_path, monkeypatch, pages=[page])
    assert row["hits"] == 999, (
        f"hits is {row['hits']}, still the old question's total; a wider "
        "query has a larger answer")


# --- A server that says "come back tomorrow" ------------------------------
#
# OpenAlex's free tier is 1,000 requests a day. When it is spent it answers
# 429 with `Retry-After: 60875` -- seventeen hours, until the budget resets at
# midnight UTC -- and the clamped `_retry_after` turned that into 60 seconds,
# retried four times, gave up, and let the supervisor restart and do it again.
# Politely, on a schedule, forever, against a service that had answered the
# question precisely.

def _hdr(value):
    """An HTTPError carrying one Retry-After header."""
    import urllib.error

    class H(dict):
        def get(self, k, d=None):
            return value if k == "Retry-After" else d

    return urllib.error.HTTPError("u", 429, "budget", H(), None)


def test_a_long_retry_after_is_read_from_the_raw_header():
    """`_retry_after` CLAMPS to 60s, which is right for a short wait and
    exactly wrong for a budget reset: it turns seventeen hours into a minute."""
    e = _hdr("60875")
    assert fx._retry_after(e) == fx.RETRY_AFTER_CAP == 60.0
    assert fx._long_retry_after(e) == 60875.0, (
        "the long wait is being clamped, so the crawler will knock again in a "
        "minute on a door that told it when it opens")


def test_a_short_retry_after_is_not_treated_as_a_scheduled_return():
    """The two paths must not merge: a 5s wait is a retry, and sleeping through
    it as a 'scheduled return' would log every blip as an outage."""
    for raw in ("5", "30", "120"):
        assert fx._long_retry_after(_hdr(raw)) is None, raw


def test_a_long_wait_is_bounded():
    """A server could say 'come back in a year'. Honouring that literally is
    indistinguishable from hanging."""
    assert fx._long_retry_after(_hdr("99999999")) == fx.LONG_WAIT_CAP <= 86_400.0


def test_a_nan_retry_after_is_not_a_wait():
    """`float("nan")` PARSES, and every comparison against it is False, so an
    unguarded NaN falls through as a valid wait and `time.sleep(nan)` raises."""
    assert fx._long_retry_after(_hdr("nan")) is None


def test_the_long_wait_does_not_consume_the_retry_budget(monkeypatch):
    """BEHAVIOURAL, because the first version of this guard was a text scan.

    It asserted that the source contained the word `continue` and called that
    proof the retry budget survived a scheduled return. But `_get` looped with
    `for attempt in range(tries)`, where `continue` ADVANCES the iterator and
    therefore spends an attempt -- the exact opposite of the claim, invisible
    to a grep, because a text scan cannot see loop semantics.

    So this measures the budget instead. With `tries=2`: one scheduled return,
    then one genuine transient failure, then success. That needs both retries
    to survive the wait. Under the old `for` loop the wait consumed the first
    attempt, the 503 consumed the second, and the call raised.
    """
    seq = [("long", "3600"), ("transient", None), ("ok", None)]
    calls = {"n": 0}
    slept = []
    monkeypatch.setattr(fx.time, "sleep", lambda s: slept.append(s))

    def fake(url, *a, **k):
        import urllib.error
        kind, val = seq[calls["n"]]
        calls["n"] += 1
        if kind == "long":
            raise _hdr(val)
        if kind == "transient":
            raise urllib.error.URLError("reset")

        class R:
            status = 200

            def read(self):
                return b"{}"

            def __enter__(self):
                return self

            def __exit__(self, *x):
                return False

        return R()

    monkeypatch.setattr(fx.urllib.request, "urlopen", fake)
    assert fx._get("https://www.ebi.ac.uk/x", tries=2, timeout=1) == b"{}", (
        "the scheduled return spent a retry, so a transient failure after a "
        "budget reset had no attempt left")
    assert 3600.0 in slept, f"waited {slept} instead of the hour asked for"
    assert calls["n"] == 3, (
        f"made {calls['n']} requests; the wait should pause rather than "
        "replace an attempt")


def test_a_server_that_always_asks_for_tomorrow_still_gives_up(monkeypatch):
    """The complement, so the fix above is not 'wait forever'.

    Not spending a retry on a scheduled return means a server that answers
    every request with a long wait would hold one worker indefinitely, so the
    scheduled returns have their own bound.
    """
    slept = []
    monkeypatch.setattr(fx.time, "sleep", lambda s: slept.append(s))

    def always_long(url, *a, **k):
        raise _hdr("60875")

    monkeypatch.setattr(fx.urllib.request, "urlopen", always_long)
    with pytest.raises(fx.TransientFetchError):
        fx._get("https://www.ebi.ac.uk/x", tries=3, timeout=1)
    long_sleeps = [s for s in slept if s > fx.LONG_WAIT_THRESHOLD]
    assert len(long_sleeps) == fx.MAX_SCHEDULED_RETURNS, (
        f"slept long {len(long_sleeps)} times against a bound of "
        f"{fx.MAX_SCHEDULED_RETURNS}; an unbounded scheduled return is a hang")


def test_the_backoff_schedule_survived_the_loop_rewrite(monkeypatch):
    """`attempt` moved from a `for` index to an explicit counter incremented
    BEFORE the backoff, so `2.0 ** (attempt + 1)` had to become `2.0 **
    attempt` to keep the measured 2+4+8 schedule. An off-by-one here doubles
    or halves every wait in the crawl, silently."""
    slept = []
    monkeypatch.setattr(fx.time, "sleep", lambda s: slept.append(s))

    def always_down(url, *a, **k):
        import urllib.error
        raise urllib.error.URLError("down")

    monkeypatch.setattr(fx.urllib.request, "urlopen", always_down)
    with pytest.raises(fx.TransientFetchError):
        fx._get("https://www.ebi.ac.uk/x", tries=4, timeout=1)
    assert slept == [2.0, 4.0, 8.0], (
        f"backoff was {slept}, not the measured 2+4+8 = 14s schedule")


# --- A permanently broken full text must not make a slice unretirable -----
#
# MEASURED on the live corpus, not hypothesised. PMC 2023 went 88 deferred ->
# 7 -> 7 across two passes while `seen` went 10,611 -> 21,155 against 10,551
# hits: the slice re-walked itself in full to retry records that cannot
# succeed. The first deferral in both passes was PMC10206698 ("E-Poster
# abstracts"), which Europe PMC advertises as isOpenAccess=Y inEPMC=Y and
# whose fullTextXML returns HTTP 500 on three consecutive direct requests.
# Deferring forever costs 2.9M search requests per pass on the MED slices.

def _rerun(tmp_path, monkeypatch, pages, fail_ids, outage=False):
    """A SECOND run over the same state, which is what `retry_pass` keys on.

    `outage=True` fails EVERY full-text request including the breaker's control
    probes, which is the only way to model a real outage now that the probe
    exists. Without it, failing a chosen set of ids models PERMANENTLY BROKEN
    DOCUMENTS -- a distinction the fixture could not previously express, which
    is why a test named for an outage was actually driving broken documents.
    """
    import corpus_identity_index as _ix
    monkeypatch.setattr(fx, "STATE_DB", tmp_path / "state.sqlite")
    monkeypatch.setattr(fx, "OUT_ROOT", tmp_path / "out")
    monkeypatch.setattr(_ix, "DB", tmp_path / "id.sqlite")
    monkeypatch.setattr(fx, "SLEEP", 0)
    monkeypatch.setattr(fx.time, "sleep", lambda *_: None)
    seq = list(pages)
    monkeypatch.setattr(fx, "search",
                        lambda s, y, c: seq.pop(0) if seq else {"resultList": {"result": []}})

    def ft(rec):
        if outage:
            raise fx.TransientFetchError("500 (the service is down)")
        # A control record from the identity index carries only a pmcid; it is
        # an article already proven fetchable, so it answers.
        if rec.get("id") is None:
            return "TEXT-control"
        if rec.get("id") in fail_ids:
            raise fx.TransientFetchError(f"500 on {rec['id']}")
        return f"TEXT-{rec['id']}"

    monkeypatch.setattr(fx, "fetch_fulltext", ft)
    totals = fx.run([2026], [("MED", True)], verbose=False)
    import sqlite3 as _s3
    c = _s3.connect(tmp_path / "state.sqlite")
    row = c.execute("SELECT seen, kept, deferred, done, hits FROM slice").fetchone()
    c.close()
    return totals, dict(zip(("seen", "kept", "deferred", "done", "hits"), row))


def test_a_second_failure_retires_the_slice_instead_of_deferring_forever(tmp_path, monkeypatch):
    """First pass defers; second pass stores the record without full text.

    Optimistic is right ONCE -- most full-text failures are transient. Being
    optimistic on every pass is what made a slice holding one permanently
    broken document impossible to retire.
    """
    page = _page(["a", "b", "c"], None)
    totals, row, stored = _drive(tmp_path, monkeypatch, [page], fail_ids={"b"})
    assert row["deferred"] == 1, "the first pass must defer, not give up"
    assert row["done"] == 0, "a slice with a deferral must not retire on pass one"

    totals2, row2 = _rerun(tmp_path, monkeypatch, [_page(["a", "b", "c"], None)], {"b"})
    assert row2["done"] == 1, (
        f"the slice still will not retire after a second failure "
        f"(deferred={row2['deferred']}); this is the non-terminating loop")
    assert row2["deferred"] == 0, "a retired slice must not keep a deferral count"
    assert totals2["no_fulltext"] == 1, (
        f"the record given up on is not reported: no_fulltext="
        f"{totals2.get('no_fulltext')}")


def test_the_record_survives_even_though_its_full_text_does_not(tmp_path, monkeypatch):
    """Retiring the slice must not drop the record.

    Only the full text is unavailable. The title, abstract and MeSH are in the
    search result already, and those are what the retrieval index matches on,
    so discarding the record would lose a paper to recover a byte count.
    """
    _drive(tmp_path, monkeypatch, [_page(["a", "b", "c"], None)], fail_ids={"b"})
    _rerun(tmp_path, monkeypatch, [_page(["a", "b", "c"], None)], {"b"})

    import glob as _g
    import gzip as _gz
    import json as _j
    ids = []
    for f in sorted(_g.glob(str(tmp_path / "out" / "**" / "*.jsonl.gz"), recursive=True)):
        try:
            with _gz.open(f, "rt") as fh:
                for ln in fh:
                    ids.append(_j.loads(ln)["epmc_id"])
        except EOFError:
            pass
    assert "b" in ids, (
        f"the record whose full text failed twice was dropped entirely: {ids}")


def test_an_outage_during_a_retry_pass_raises_instead_of_retiring_the_slice(tmp_path, monkeypatch):
    """The dangerous case, and the reason this is not just 'accept failures'.

    If a retry pass treats every failure as a permanently broken document,
    then an outage DURING that pass retires the whole slice as abstract-only
    and there is no third chance -- the loss is permanent and silent.

    WHY IT RAISES CHANGED WITH #964, and the behaviour it guards did not. A
    long consecutive run is no longer read AS an outage -- on a retry pass that
    premise is backwards, since the records being retried are the ones that
    already failed. The run now re-probes an article it fetched successfully
    and only stops if that fails too. This fixture is a retry-ONLY pass, so it
    So this fixture now fails EVERY request, control probes included, which is
    what an outage actually looks like -- and the run still refuses to retire
    the slice. Before `outage=True` existed the fixture failed only a chosen
    set of ids, which is permanently broken DOCUMENTS; the test was named for
    an outage while driving the other case entirely.
    """
    # Pass one must defer MANY records without tripping its OWN abort, so the
    # failures alternate with successes and no consecutive run reaches the
    # threshold. Pass two then sees only the deferred records -- the stored
    # ones are held and skipped -- so they fail back to back and the retry
    # guard fires.
    n = 2 * fx.MAX_CONSECUTIVE_DEFERRALS + 20
    ids = [f"r{i}" for i in range(n)]
    doomed = {i for i in ids if int(i[1:]) % 2 == 0}
    assert len(doomed) > fx.MAX_CONSECUTIVE_DEFERRALS, "fixture too small to trip the guard"
    _drive(tmp_path, monkeypatch, [_page(ids, None)], fail_ids=doomed)

    import sqlite3 as _s3
    c = _s3.connect(tmp_path / "state.sqlite")
    d1, done1 = c.execute("SELECT deferred, done FROM slice").fetchone()
    c.close()
    assert done1 == 0 and d1 == len(doomed), (
        f"pass one did not set up the retry state: deferred={d1} done={done1}")

    with pytest.raises(RuntimeError) as err:
        _rerun(tmp_path, monkeypatch, [_page(ids, None)], doomed, outage=True)
    assert "retry pass" in str(err.value), str(err.value)
    assert "really is down" in str(err.value), (
        "with every request failing -- the control probes included -- the "
        f"reason must be a MEASURED outage: {err.value}")

    import sqlite3 as _s3
    c = _s3.connect(tmp_path / "state.sqlite")
    done = c.execute("SELECT done FROM slice").fetchone()[0]
    c.close()
    assert done == 0, (
        "the slice was retired during an outage, so its full texts are gone "
        "for good")


# ---------------------------------------------------------------------------
# #910: an absent hitCount must not be read as a hitCount of zero.
#
# `hits` was written as `d.get("hitCount", 0)` on every page, so a response
# that parsed as JSON without the key wrote 0. That does not retire a slice
# early -- retirement is cursor exhaustion -- but a retired slice carrying
# hits=0 satisfies `seen >= hits` trivially, so the ledger loses the only
# check it has on whether the walk was complete. Measured on the live ledger
# when found: MED 2022 (seen 138,000) and MED 1990 (seen 10,000) were both
# done=1 with hits=0.
# ---------------------------------------------------------------------------

def _seed_slice(tmp_path, monkeypatch, hits):
    """A slice already carrying a known total, mid-walk on the CURRENT question."""
    monkeypatch.setattr(fx, "STATE_DB", tmp_path / "state.sqlite")
    st = fx.state()
    st.execute("INSERT INTO slice(src,year,cursor,seen,kept,hits,done,query_id) "
               "VALUES ('MED',2026,'*',0,0,?,0,?)",
               (hits, fx.query_fingerprint()))
    st.commit()
    st.close()


def test_a_response_without_a_hitcount_leaves_the_recorded_total_alone(
        tmp_path, monkeypatch):
    """ABSENT is not ZERO. The prior total survives a degraded response."""
    _seed_slice(tmp_path, monkeypatch, hits=777)
    page = _page(["A", "B"], None)
    del page["hitCount"]
    totals, row, stored = _drive(tmp_path, monkeypatch, pages=[page])
    assert stored == ["A", "B"], "the page was still walked"
    assert row["hits"] == 777, (
        f"hits is {row['hits']}: a response carrying no hitCount overwrote the "
        "recorded total, so a retired slice can read as complete against zero")


def test_a_genuine_hitcount_of_zero_is_still_recorded(tmp_path, monkeypatch):
    """The fix must distinguish absent from zero, not merely refuse to write 0.

    Skipping every zero would be the easy over-correction: an empty slice
    really does have a total of zero and must be able to say so, otherwise the
    stale value persists and reads as a slice that was never finished.
    """
    _seed_slice(tmp_path, monkeypatch, hits=777)
    empty = {"resultList": {"result": []}, "nextCursorMark": None, "hitCount": 0}
    totals, row, stored = _drive(tmp_path, monkeypatch, pages=[empty])
    assert stored == [], "nothing to store"
    assert row["hits"] == 0, (
        f"hits is {row['hits']}: an explicit hitCount of 0 must be recorded, or "
        "a genuinely empty slice keeps a total it no longer has")


def test_the_update_is_conditional_on_the_key_being_present(tmp_path):
    """Spelled against the source, because the two cases above share a path.

    A `d.get("hitCount", 0)` restores the defect while both behavioural tests
    could still pass if the fixtures happened to carry the key.
    """
    src = Path(fx.__file__).read_text()
    assert 'd.get("hitCount", 0)' not in src, (
        "an absent hitCount defaults to 0 again (#910)")
    assert 'reported = d.get("hitCount")' in src
    assert "if reported is not None:" in src, (
        "the UPDATE is unconditional again, so absent reads as zero")


def test_the_page_list_and_the_hits_column_no_longer_share_a_name(tmp_path):
    """One word for two quantities is what made #910 hard to rule out.

    The local held the RECORDS and the column holds the SIZE, and the
    retirement test reads the list two lines from an UPDATE that writes the
    count. Pinned so the confusing spelling cannot come back.
    """
    src = Path(fx.__file__).read_text()
    i = src.index("stopped_early = False")
    head = src[src.index("_search_with_retry(src, year, cursor", 0):i]
    assert 'page = d.get("resultList", {}).get("result", [])' in head
    assert 'hits = d.get("resultList"' not in src, (
        "the page list is called `hits` again, which is also the column name")


# ---------------------------------------------------------------------------
# #913 step 1 (folded in here): a walk that retires short of its recorded
# total must SAY SO. Two slices retired that way before anyone noticed --
# MED 2022 at 138,000 where its neighbours bracket ~486,000, and MED 1990 at
# 10,000 against 1989's 64,194, about 402,000 records absent from slices the
# ledger called complete. Deliberately a report and not a rewind: an
# unconditional rewind-on-short is the same shape as the non-terminating
# deferral loop fixed in #907.
# ---------------------------------------------------------------------------

def test_a_walk_that_retires_short_of_the_recorded_total_says_so(
        tmp_path, monkeypatch, capsys):
    """The alarm the ledger never had, and it must fire even when quiet.

    `_drive` runs with `verbose=False`; this print is deliberately NOT gated on
    it, because an incomplete slice marked complete is not progress chatter.
    """
    page = _page(["A", "B"], None)
    page["hitCount"] = 5000
    totals, row, stored = _drive(tmp_path, monkeypatch, pages=[page])
    assert row["done"] == 1, "the slice still retires; this is a report, not a rewind"
    assert row["seen"] == 2 and row["hits"] == 5000
    out = capsys.readouterr().out
    assert "RETIRED SHORT" in out, (
        "a slice that walked 2 of 5,000 retired silently, which is exactly how "
        "~402,000 records went missing without anyone noticing")
    assert "2 of 5,000" in out, "the alarm must name both numbers"


def test_a_complete_walk_raises_no_short_alarm(tmp_path, monkeypatch, capsys):
    """Non-vacuity: the alarm must be able to stay quiet.

    `_page` sets hitCount to the page length, so this walk is complete. An
    alarm that fires on every slice would be worthless and would pass the test
    above.
    """
    totals, row, stored = _drive(tmp_path, monkeypatch, pages=[_page(["A", "B"], None)])
    assert row["done"] == 1 and row["seen"] == row["hits"] == 2
    out = capsys.readouterr().out
    assert "RETIRED SHORT" not in out
    # BOTH alarms, because they sit in an if/elif: an always-true `no_total`
    # never fires while `short` is true, so a test that only exercises the
    # short branch cannot see it. That mutation survived the first sweep.
    assert "RETIRED WITH NO TOTAL" not in out, (
        "a complete walk with a known total reported its completeness as "
        "unknown")


def test_a_rewalked_slice_does_not_trip_the_short_alarm(tmp_path, monkeypatch, capsys):
    """`seen` accumulates across runs while `hits` is one question's total.

    A slice re-walked under the deferral path has an INFLATED count, so
    comparing the two could only ever false-positive in the safe direction --
    but that has to be true rather than assumed, because the obvious
    implementation compares the per-run local instead of the column.
    """
    # The fixture has to put the LOCAL and the COLUMN on opposite sides of the
    # total, or it cannot tell the two implementations apart. A first version
    # seeded hits=2 and drove `_page(["A"])`, which sets hitCount to the page
    # length -- so hits was overwritten to 1 and NEITHER the correct code nor
    # the per-run-local mutant fired. It passed while testing nothing.
    #
    # Here: the column accumulates to 101 while this run walks 1 record against
    # a recorded total of 50. The column says complete, the local says short.
    _seed_slice(tmp_path, monkeypatch, hits=50)
    st = fx.state()
    st.execute("UPDATE slice SET seen=100 WHERE src='MED' AND year=2026")
    st.commit()
    st.close()
    page = _page(["A"], None)
    page["hitCount"] = 50
    totals, row, stored = _drive(tmp_path, monkeypatch, pages=[page])
    assert row["seen"] == 101 and row["hits"] == 50, (
        f"fixture broken: seen={row['seen']} hits={row['hits']}; the column must "
        "exceed the total while this run's own count falls short of it")
    assert "RETIRED SHORT" not in capsys.readouterr().out, (
        "a resumed slice was called short, so the check is reading the per-run "
        "local instead of the accumulated column")


def test_the_short_check_reads_the_column_not_the_per_run_local(tmp_path):
    """Spelled against the source: the local `seen` counts only THIS run.

    Using it would call every resumed slice short. The check must read the
    accumulated column back from the ledger.
    """
    src = Path(fx.__file__).read_text()
    assert 'SELECT seen, hits FROM slice WHERE src=? AND year=?' in src, (
        "the short check no longer reads the accumulated column")
    # The DATA-FLOW CHAIN, not a character window. This was `src[i-900:i]`,
    # which broke the moment a comment was added above it -- a fixed offset
    # fails on unrelated prose and passes on a real defect that happens to fit.
    assert "walked, total = (tally or (0, -1))" in src, (
        "the alarm no longer unpacks the ledger read")
    assert "short = total > 0 and walked < total" in src, (
        "the short test is not computed from the ledger read, so it may be "
        "using the per-run local and calling every resumed slice short")
    assert src.index("tally = st.execute(") < src.index("RETIRED SHORT"), (
        "the ledger is read after the alarm, which cannot work")

# ---------------------------------------------------------------------------
# #913: finding and reopening slices the ledger cannot certify as complete.
# ---------------------------------------------------------------------------

def _ledger(tmp_path, monkeypatch, rows):
    """A ledger seeded with (src, year, seen, hits, done) rows and nothing else."""
    monkeypatch.setattr(fx, "STATE_DB", tmp_path / "state.sqlite")
    c = fx.state()
    for src, year, seen, hits, done in rows:
        c.execute("INSERT INTO slice(src,year,cursor,seen,kept,hits,done,query_id) "
                  "VALUES (?,?,'*',?,0,?,?,'q')", (src, year, seen, hits, done))
    c.commit()
    return c


def test_a_slice_that_walked_short_of_its_total_is_found(tmp_path, monkeypatch):
    c = _ledger(tmp_path, monkeypatch, [("PPR", 2022, 23340, 23365, 1)])
    assert [(r[0], r[1]) for r in fx.incomplete_slices(c)] == [("PPR", 2022)]
    c.close()


def test_a_slice_with_no_recorded_total_is_found(tmp_path, monkeypatch):
    """The historical shape detector 1 cannot see, because it needs hits > 0."""
    c = _ledger(tmp_path, monkeypatch,
                [("MED", 2022, 138000, 0, 1), ("MED", 1990, 10000, -1, 1)])
    found = [(r[0], r[1]) for r in fx.incomplete_slices(c)]
    assert set(found) == {("MED", 2022), ("MED", 1990)}, (
        "both the zeroed and the never-recorded totals must be found; the "
        "zeroing is exactly what hid ~402,000 records")
    c.close()


def test_a_complete_slice_is_not_found(tmp_path, monkeypatch):
    """Non-vacuity. A detector that flags everything would pass the tests above."""
    c = _ledger(tmp_path, monkeypatch,
                [("MED", 1989, 64206, 64194, 1), ("MED", 1973, 24323, 24323, 1)])
    assert fx.incomplete_slices(c) == [], (
        "a healthy slice (seen >= hits) was flagged; note 1989 has seen > hits, "
        "which is normal for a re-walked slice and must not read as damage")
    c.close()


def test_an_open_slice_is_not_found(tmp_path, monkeypatch):
    """`done=0` is already going to be re-walked; flagging it is noise."""
    c = _ledger(tmp_path, monkeypatch, [("MED", 2015, 101993, 0, 0)])
    assert fx.incomplete_slices(c) == []
    c.close()


def test_a_round_walk_total_is_not_treated_as_evidence(tmp_path, monkeypatch):
    """The page-boundary fingerprint must NOT be a detector.

    It is how the two damaged slices were found, but a slice whose true total
    happens to be an exact multiple of PAGE would be reopened forever.
    """
    c = _ledger(tmp_path, monkeypatch, [("MED", 2000, 10000, 10000, 1)])
    assert fx.incomplete_slices(c) == [], (
        "a complete slice was flagged for the shape of its number, which would "
        "reopen it on every pass")
    c.close()


def test_reopening_clears_done_and_rewinds_only_the_flagged_slices(
        tmp_path, monkeypatch, capsys):
    rows = [("MED", 2022, 138000, 0, 1),          # flagged
            ("MED", 1989, 64206, 64194, 1)]       # healthy, must be untouched
    c = _ledger(tmp_path, monkeypatch, rows)
    c.execute("UPDATE slice SET cursor='somewhere' WHERE src='MED'")
    c.commit()
    c.close()
    monkeypatch.setattr(sys, "argv", ["corpus_expand_fetch.py", "--reopen-incomplete"])
    assert fx.main() == 0
    c = fx.state()
    got = dict((f"{s}{y}", (cur, done)) for s, y, cur, done in
               c.execute("SELECT src, year, cursor, done FROM slice"))
    c.close()
    assert got["MED2022"] == ("*", 0), "the flagged slice was not rewound"
    assert got["MED1989"] == ("somewhere", 1), (
        "a healthy slice was reopened; only the flagged ones may move")
    out = capsys.readouterr().out
    assert "reopened 1 slice" in out


def test_the_report_alone_writes_nothing(tmp_path, monkeypatch, capsys):
    """`--incomplete` must be read-only; the write needs its own flag."""
    c = _ledger(tmp_path, monkeypatch, [("MED", 2022, 138000, 0, 1)])
    c.close()
    monkeypatch.setattr(sys, "argv", ["corpus_expand_fetch.py", "--incomplete"])
    assert fx.main() == 0
    c = fx.state()
    cur, done = c.execute("SELECT cursor, done FROM slice").fetchone()
    c.close()
    assert done == 1 and cur == "*", "the report mutated the ledger"
    assert "--reopen-incomplete" in capsys.readouterr().out, (
        "the report must say how to act on it")


def test_a_slice_marked_done_having_walked_nothing_is_found(tmp_path, monkeypatch):
    """The larger half of #913, invisible to both earlier detectors.

    `seen = 0, hits = -1, done = 1` -- nothing walked, no total, retired. The
    shortfall detector needs `hits > 0` and the zeroed detector needs
    `seen > 0`, so ZERO SATISFIED BOTH and the worst row the ledger can hold
    was the one it could not see. 24 consecutive MED slices (1991-2014) sat
    here holding 3,110,537 records by the service's own count.
    """
    c = _ledger(tmp_path, monkeypatch, [("MED", 2000, 0, -1, 1)])
    assert [(r[0], r[1]) for r in fx.incomplete_slices(c)] == [("MED", 2000)], (
        "a slice marked done having walked nothing was not flagged; this is "
        "the shape that hid 3.1M records behind a complete-looking ledger")
    c.close()


def test_a_genuinely_empty_slice_is_not_flagged(tmp_path, monkeypatch):
    """Why the clause is keyed on `hits < 0` and NOT on `seen = 0`.

    An empty year walks nothing too, and it is not damage: 112 PPR slices
    (1900-2011) are legitimately empty because preprints did not exist, and
    the service returns hitCount 0 for them (measured, not assumed). Keying
    the detector on `seen = 0` would reopen all of them on every pass forever.
    After #912 a zero total written by a real response is trustworthy, so an
    empty slice must be allowed to settle at 0.
    """
    c = _ledger(tmp_path, monkeypatch,
                [("PPR", 1950, 0, 0, 1), ("PPR", 2010, 0, 0, 1)])
    assert fx.incomplete_slices(c) == [], (
        "a legitimately empty slice was flagged; keyed on seen=0 rather than "
        "on the absence of a total, this reopens 112 PPR slices forever")
    c.close()


def test_the_legacy_zero_total_migration_runs_exactly_once(tmp_path, monkeypatch):
    """v1 restores the -1 sentinel on pre-#912 rows, and then stops.

    The idempotence is the load-bearing half. If the migration ran on every
    open it would keep de-certifying the genuine zeros #912 made trustworthy,
    so every empty PPR slice would be reopened on every pass -- the exact
    behaviour the detector is keyed to avoid.
    """
    monkeypatch.setattr(fx, "STATE_DB", tmp_path / "state.sqlite")

    # A ledger as it stood before #912: a zero total, and no version marker.
    c = fx.state()
    c.execute("INSERT INTO slice(src,year,cursor,seen,kept,hits,done,query_id) "
              "VALUES ('MED',2000,'*',0,0,0,1,'q')")
    c.execute("PRAGMA user_version=0")
    c.commit()
    c.close()

    c = fx.state()
    hits, = c.execute("SELECT hits FROM slice WHERE year=2000").fetchone()
    ver, = c.execute("PRAGMA user_version").fetchone()
    assert hits == -1, (
        f"legacy zero total is still {hits}; it reads as a certified-empty "
        "slice and the detector cannot distinguish it from a never-walked one")
    assert ver == fx.SCHEMA_VERSION, f"version marker not set: {ver}"
    assert [(r[0], r[1]) for r in fx.incomplete_slices(c)] == [("MED", 2000)], (
        "the migrated row must now be visible to the detector")

    # A genuine post-#912 empty walk, recorded AFTER the migration ran.
    c.execute("INSERT INTO slice(src,year,cursor,seen,kept,hits,done,query_id) "
              "VALUES ('PPR',1950,'*',0,0,0,1,'q')")
    c.commit()
    c.close()

    c = fx.state()
    hits, = c.execute("SELECT hits FROM slice WHERE src='PPR'").fetchone()
    c.close()
    assert hits == 0, (
        f"a genuine zero recorded after the migration was de-certified to "
        f"{hits}; the migration re-ran and will reopen every empty slice")


def test_retiring_with_no_recorded_total_at_all_also_says_so(
        tmp_path, monkeypatch, capsys):
    """The gap the #910 fix CREATED, on a fresh slice.

    `hits` defaults to -1 and an absent hitCount now leaves it alone, so a slice
    whose responses never carried one retires with total = -1. `short` is False
    there (it needs total > 0), so without a second alarm this would be exactly
    as silent as the defect being fixed.
    """
    page = _page(["A", "B"], None)
    del page["hitCount"]
    totals, row, stored = _drive(tmp_path, monkeypatch, pages=[page])
    assert row["done"] == 1 and row["seen"] == 2
    assert row["hits"] == -1, (
        f"hits is {row['hits']}, not the schema default; the fixture is not "
        "exercising a never-recorded total")
    out = capsys.readouterr().out
    assert "RETIRED WITH NO TOTAL" in out, (
        "a slice retired with no total ever reported and said nothing, which is "
        "the same silence #910 was about")
    assert "RETIRED SHORT" not in out, (
        "the two alarms must be distinguishable: there is no total to be short of")


def test_the_two_alarms_are_mutually_exclusive(tmp_path, monkeypatch, capsys):
    """A known total that was not reached is SHORT, not unknown."""
    page = _page(["A"], None)
    page["hitCount"] = 900
    totals, row, stored = _drive(tmp_path, monkeypatch, pages=[page])
    out = capsys.readouterr().out
    assert "RETIRED SHORT" in out
    assert "RETIRED WITH NO TOTAL" not in out, (
        "a slice with a recorded total of 900 reported its completeness as "
        "unknown; the total is known and was missed")


def test_an_empty_slice_is_not_called_unknown(tmp_path, monkeypatch, capsys):
    """Nothing walked and no total is an EMPTY slice, not an unverifiable one.

    Without the `walked > 0` clause the alarm fires on every empty source-year
    -- and this crawl has dozens, every PPR year before 1946 among them. An
    alarm that fires on a third of the ledger is one nobody reads, which is how
    the thing it exists to catch gets missed.
    """
    empty = {"resultList": {"result": []}, "nextCursorMark": None}
    totals, row, stored = _drive(tmp_path, monkeypatch, pages=[empty])
    assert row["seen"] == 0 and row["hits"] == -1, (
        f"fixture broken: seen={row['seen']} hits={row['hits']}")
    out = capsys.readouterr().out
    assert "RETIRED WITH NO TOTAL" not in out, (
        "an empty slice was reported as having unknown completeness")
    assert "RETIRED SHORT" not in out


# --- The output root must not be invented (#960) ---------------------------
#
# `OUT_ROOT.mkdir(parents=True, exist_ok=True)` is the obvious line and it is
# the dangerous one, because A MOUNT POINT IS JUST A DIRECTORY. With the NAS
# unmounted `~/nas` is an empty LOCAL directory, so `parents=True` creates
# `~/nas/cancer-atlas/expanded` on the boot disk and the crawl writes the
# expanded corpus -- 51 GB of it, against 16 GB free -- into the wrong
# filesystem, silently, while the ledger (which lives in the repo, not on the
# share) records every slice as walked. Those shards are then invisible to
# every later run: the same permanent-invisibility failure the write-then-
# index ordering exists to prevent, reached from the other side.
#
# No behavioural test on a machine with the share mounted can catch that, and
# the defect is ONE KEYWORD, which is why the structural guard below carries
# as much weight as the two behavioural ones.


def test_an_unmounted_output_root_is_refused_rather_than_created(tmp_path):
    """The live failure. `parent` absent is what an unmounted mount point
    looks like from inside this process."""
    root = tmp_path / "nas" / "cancer-atlas" / "expanded"
    with pytest.raises(SystemExit) as e:
        fx.require_output_root(root)
    assert not root.exists(), "the refusal must not leave a directory behind"
    assert not root.parent.exists(), "and must not create the path to it either"
    msg = str(e.value)
    assert str(root.parent) in msg, "the message must name what is missing"
    assert "FERRO_EXPAND_OUT" in msg, "and how to point it elsewhere"
    assert "resume" in msg, (
        "a refusal that reads as data loss would invite someone to force past "
        "it; the ledger is local, so nothing is lost by stopping here")


def test_an_existing_parent_lets_the_crawl_create_its_own_leaf(tmp_path):
    """The normal case, and the one every test in this file relies on: the
    parent exists, the leaf does not, and the crawl owns the leaf."""
    root = tmp_path / "expanded"
    fx.require_output_root(root)
    assert root.is_dir()
    fx.require_output_root(root)      # idempotent across resumed runs


def test_the_crawl_never_creates_its_output_path_with_parents():
    """STRUCTURAL. The behavioural guards above can only fail on a machine
    where the share happens to be absent; this one fails anywhere, and it is
    the reintroduction of `parents=True` that it exists to stop.

    Scoped to the shard-writing path: `Shards.__init__` legitimately creates
    its own subdirectory beneath a root this function has already vouched for.
    """
    src = (REPO / "scripts" / "corpus_expand_fetch.py").read_text()
    tree = ast.parse(src)
    offenders = []
    for fn in ast.walk(tree):
        if not isinstance(fn, ast.FunctionDef) or fn.name != "run":
            continue
        for node in ast.walk(fn):
            if (isinstance(node, ast.Call)
                    and isinstance(node.func, ast.Attribute)
                    and node.func.attr == "mkdir"):
                for kw in node.keywords:
                    if kw.arg == "parents" and getattr(kw.value, "value", False):
                        offenders.append(node.lineno)
    assert not offenders, (
        f"run() creates its output path with parents=True at line(s) "
        f"{offenders}: on an unmounted share that writes the corpus to the "
        "boot disk instead of refusing")
    # And the refusal must actually be reached, or the above is vacuous.
    run = next(n for n in ast.walk(tree)
               if isinstance(n, ast.FunctionDef) and n.name == "run")
    called = {n.func.id for n in ast.walk(run)
              if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)}
    assert "require_output_root" in called, (
        "run() must consult the refusal; a guard against parents=True is "
        "worthless if nothing checks the path at all")


# --- A closed shard is durable, not merely written (#960) ------------------


def test_closing_a_shard_syncs_it_to_storage():
    """STRUCTURAL, because an fsync has no observable effect to assert on.

    The per-record flush reaches the OS and deliberately stops there: an
    fsync per record across a network mount costs more than the one-record
    window it closes. That trade is right for a process kill, where the bytes
    are already the kernel's, and WRONG for a mount that disappears, because
    an SMB client holds written bytes of its own and loses them with the
    connection. Measured: a dropped mount left expanded-01009 with 830
    readable records and a torn tail whose records were already indexed.

    fsync here cannot save the OPEN shard's tail and does not claim to. It
    bounds the damage to one shard, so a later drop cannot reach back and
    truncate work finished minutes ago.
    """
    src = (REPO / "scripts" / "corpus_expand_fetch.py").read_text()
    cls = next(n for n in ast.walk(ast.parse(src))
               if isinstance(n, ast.ClassDef) and n.name == "Shards")
    close = next(n for n in cls.body
                 if isinstance(n, ast.FunctionDef) and n.name == "close")
    syncs = [n for n in ast.walk(close)
             if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)
             and n.func.attr == "fsync"]
    assert syncs, "Shards.close must fsync the finished shard"


def test_the_shard_filename_has_exactly_one_definition():
    """write() and close() must name the same file. They were two copies of
    one format string, so the fsync could have been sent to a path the writer
    had stopped using -- an fsync on the wrong file succeeds silently, which
    is the worst available outcome for a durability fix."""
    src = (REPO / "scripts" / "corpus_expand_fetch.py").read_text()
    # The INDEX-to-filename mapping, not every mention of the prefix: the
    # sibling `glob(f"{self.name}-*.jsonl.gz")` in _resume enumerates shards
    # rather than naming one, so it is a different operation and a first
    # version of this guard wrongly flagged it.
    assert src.count(":05d}.jsonl.gz") == 1, (
        "the shard filename is built from its index in more than one place; a "
        "durability call that disagrees with the writer about the path still "
        "returns 0")


def test_a_closed_shard_is_a_complete_readable_gzip(tmp_path):
    """Behavioural companion: whatever close() now does, it must still leave a
    shard a reader can finish without raising."""
    s = fx.Shards(tmp_path, "t")
    for i in range(5):
        s.write({"pmid": str(i)})
    s.close()
    p = tmp_path / "t" / "t-00000.jsonl.gz"
    with gzip.open(p, "rt", encoding="utf-8") as fh:
        got = [json.loads(line)["pmid"] for line in fh]
    assert got == ["0", "1", "2", "3", "4"]


def test_close_syncs_the_shard_it_actually_wrote(tmp_path, monkeypatch):
    """BEHAVIOURAL, and it exists because a structural "does it fsync?" guard
    SURVIVED the mutation that matters.

    Pointing the sync one index past the open shard still passed every other
    guard: `os.open` on a path that does not exist raises FileNotFoundError,
    an OSError, which the handler tolerates by design -- so the fsync became
    a silent no-op while the code still plainly called fsync. Identity is
    checked by inode, which is the only thing that distinguishes syncing the
    right file from syncing successfully.
    """
    import os as _os
    seen = []
    real = _os.fsync
    monkeypatch.setattr(_os, "fsync",
                        lambda fd: (seen.append(_os.fstat(fd).st_ino), real(fd))[1])
    s = fx.Shards(tmp_path, "t")
    s.write({"pmid": "1"})
    s.close()
    p = tmp_path / "t" / "t-00000.jsonl.gz"
    assert seen == [p.stat().st_ino], (
        "close must fsync the shard it wrote; got "
        f"{seen} against inode {p.stat().st_ino}")


# --- The breaker must measure the outage, not infer it (#964) --------------
#
# `run()` aborts after MAX_CONSECUTIVE_DEFERRALS consecutive full-text
# failures, and the retry-pass site used to justify that with "a long run of
# consecutive failures here is an outage, not a set of broken documents".
#
# ON A RETRY PASS THAT IS BACKWARDS. The records being retried are the ones
# that already failed -- a population selected for failing -- so a long
# consecutive run is the EXPECTED outcome of a set of permanently broken
# documents.
#
# Measured when it stopped the #913 recovery walk at 426 minutes: five PMCIDs
# the crawl had given up on returned HTTP 500 on a hand probe, five whose text
# it had successfully stored returned HTTP 200 with real content, and the
# search endpoint answered normally. Europe PMC serves 500 rather than 404 for
# articles with no retrievable full text, so the run aborted on a permanent
# condition it had diagnosed as transient.


def test_the_service_is_up_when_a_control_still_answers(monkeypatch):
    calls = []

    def fake_ft(rec):
        calls.append(rec["pmcid"])
        return "TEXT"

    monkeypatch.setattr(fx, "fetch_fulltext", fake_ft)
    assert fx._service_is_up([{"pmcid": "PMC1"}, {"pmcid": "PMC2"}], None) is True
    assert len(calls) == 1, "one success is enough; it must stop probing"
    assert calls[0] == "PMC2", "it must probe the NEWEST success first"


def test_the_service_is_down_when_no_control_answers(monkeypatch):
    calls = []

    def fake_ft(rec):
        calls.append(rec["pmcid"])
        raise fx.TransientFetchError("500")

    monkeypatch.setattr(fx, "fetch_fulltext", fake_ft)
    assert fx._service_is_up(
        [{"pmcid": f"PMC{i}"} for i in range(1, 5)], None) is False
    assert len(calls) == fx.ABORT_CONTROL_PROBES, (
        f"it must try {fx.ABORT_CONTROL_PROBES} controls before concluding the "
        f"service is down; tried {len(calls)}")


def test_a_refusal_still_proves_the_service_answered(monkeypatch):
    """`_get` returns None when the server REFUSES. A 404 is an answer, so a
    control that has since been withdrawn still shows the service is up --
    the question is reachability, not whether one article survives."""
    monkeypatch.setattr(fx, "fetch_fulltext", lambda rec: None)
    assert fx._service_is_up([{"pmcid": "PMC1"}], None) is True


def test_no_control_is_a_THIRD_state_and_aborts(monkeypatch):
    """A run that begins on a retry pass can hit the breaker having fetched no
    text at all. "Unknown" must not collapse into "up" -- that would retire a
    slice without its text during a real outage, which is permanent."""
    monkeypatch.setattr(fx, "fetch_fulltext",
                        lambda *a, **k: pytest.fail("must not probe with no control"))
    assert fx._service_is_up([], None) is None
    with pytest.raises(RuntimeError) as e:
        fx._abort_unless_documents_are_broken(100, [], None, Exception("x"),
                                              "a retry pass")
    assert "UNKNOWN" in str(e.value), (
        "the message must distinguish 'no control to check' from 'the controls "
        f"failed'; got {e.value}")


def test_broken_documents_do_not_abort_the_run(monkeypatch):
    """The regression that matters: this is what stopped the #913 walk."""
    monkeypatch.setattr(fx, "fetch_fulltext", lambda rec: "TEXT")
    # Returns rather than raising == the caller carries on.
    assert fx._abort_unless_documents_are_broken(
        100, [{"pmcid": "PMC1"}], None, Exception("500"), "a retry pass") is None


def test_a_real_outage_still_aborts(monkeypatch):
    """The guard the breaker exists for is NOT weakened: with the service
    genuinely unreachable it still stops, because retiring a slice without
    its text cannot be undone."""
    monkeypatch.setattr(fx, "fetch_fulltext",
                        lambda rec: (_ for _ in ()).throw(fx.TransientFetchError("500")))
    with pytest.raises(RuntimeError) as e:
        fx._abort_unless_documents_are_broken(
            100, [{"pmcid": "PMC1"}], None, Exception("500"), "a retry pass")
    msg = str(e.value)
    assert "really is down" in msg and "permanent" in msg, msg


def test_both_breaker_sites_consult_the_same_decision():
    """STRUCTURAL. The two sites make the same judgement call, and the retry
    site is the one that fired -- fixing it alone would leave the first-pass
    site inferring an outage from the same bad premise."""
    src = (REPO / "scripts" / "corpus_expand_fetch.py").read_text()
    run = next(n for n in ast.walk(ast.parse(src))
               if isinstance(n, ast.FunctionDef) and n.name == "run")
    calls = [n for n in ast.walk(run)
             if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)]
    aborts = [n for n in calls if n.func.id == "_should_abort"]
    guards = [n for n in calls if n.func.id == "_abort_unless_documents_are_broken"]
    assert len(aborts) == 2, f"expected both breaker sites; found {len(aborts)}"
    assert len(guards) == len(aborts), (
        f"{len(aborts)} breaker sites but {len(guards)} consult the probe -- a "
        "site that still infers an outage from the failures themselves")
    # And no site may raise an outage RuntimeError on its own again.
    raises = [n for n in ast.walk(run) if isinstance(n, ast.Raise)
              and isinstance(n.exc, ast.Call)
              and getattr(n.exc.func, "id", "") == "RuntimeError"]
    for r in raises:
        seg = ast.get_source_segment(src, r) or ""
        assert "looks down" not in seg, (
            "an outage verdict is being asserted inline again rather than measured")


def test_the_rest_base_url_has_one_definition():
    """The probe must ask the same host the fetcher does. Three call sites
    spelling out one URL is how a check ends up answering about somewhere
    else."""
    src = (REPO / "scripts" / "corpus_expand_fetch.py").read_text()
    assert src.count("europepmc/webservices/rest") == 1, (
        "the Europe PMC base URL is spelled out more than once")


def test_the_control_buffer_is_bounded_and_actually_filled():
    """A run stores millions of records; the controls are the newest few.

    The first version checked `"if pmcid:" in src`, which SURVIVED the mutation
    that stopped recording controls altogether -- that substring appears
    elsewhere in the file, so the guard passed by matching a different site.
    The same over-broad-text-scan shape as a `pgrep` pattern matching its own
    command line. AST now, over the statement that does the recording.
    """
    import ast as _ast
    src = (REPO / "scripts" / "corpus_expand_fetch.py").read_text()
    tree = _ast.parse(src)
    appends = [n for n in _ast.walk(tree)
               if isinstance(n, _ast.Call) and isinstance(n.func, _ast.Attribute)
               and n.func.attr == "append"
               and getattr(n.func.value, "id", "") == "recent_ok"]
    assert appends, "nothing records a control, so the probe can never run"
    # Not sitting behind a constant-false test, which is how the mutation hid.
    guards = [n for n in _ast.walk(tree) if isinstance(n, _ast.If)
              and any(a in _ast.walk(n) for a in appends)]
    for g in guards:
        seg = _ast.get_source_segment(src, g.test) or ""
        assert seg.strip() not in {"False", "0", "None"}, (
            f"the control recording is unreachable: `if {seg}`")
    assert "del recent_ok[0]" in src, "the control buffer must be capped"


def test_more_than_one_control_is_probed_before_calling_the_service_down():
    """The VALUE is bounded, because the mechanism test derives its expectation
    from this constant and therefore cannot see it change -- the
    guard-computes-its-own-expectation shape, confirmed by a mutation to 1
    surviving.

    Two is the floor for a reason: a single control that happens to blip would
    abort a run holding days of work on one request.
    """
    assert fx.ABORT_CONTROL_PROBES >= 2, (
        f"ABORT_CONTROL_PROBES={fx.ABORT_CONTROL_PROBES}: one control makes a "
        "transient blip indistinguishable from an outage")
    assert fx.ABORT_CONTROL_PROBES <= 5, (
        "and it must stay small: these are extra requests issued at the moment "
        "the service is already suspected of being in trouble")


def test_a_retry_pass_retires_the_slice_when_the_service_still_answers(
        tmp_path, monkeypatch):
    """THE #913 BLOCKER, as a regression test.

    The recovery walk stopped here at 426 minutes with "the service looks down
    rather than holding broken documents". It was not down: five articles the
    crawl had given up on returned HTTP 500 on a hand probe while five whose
    text it had stored returned HTTP 200, and the search endpoint answered
    normally. Europe PMC serves 500 rather than 404 for articles with no
    retrievable full text, so those documents fail forever.

    With the service answering, the slice must RETIRE -- records stored with
    title, abstract and MeSH and without full text, which is already the
    documented behaviour for a text that has failed four in-run retries twice
    across separate runs. The alternative is a slice that can never complete.
    """
    n = 2 * fx.MAX_CONSECUTIVE_DEFERRALS + 20
    ids = [f"r{i}" for i in range(n)]
    doomed = {i for i in ids if int(i[1:]) % 2 == 0}
    _drive(tmp_path, monkeypatch, [_page(ids, None)], fail_ids=doomed)

    # NOTHING IS PATCHED AT THE DECISION. The probe finds its own control --
    # pass one stored full text for the odd ids, so the identity index can name
    # an article with known-good text even though this pass has fetched none.
    # That index fallback is what makes the check usable on a retry-only pass,
    # which is exactly the state that blocked #913.
    _rerun(tmp_path, monkeypatch, [_page(ids, None)], doomed)

    import sqlite3 as _s3
    c = _s3.connect(tmp_path / "state.sqlite")
    done, seen, hits = c.execute("SELECT done, seen, hits FROM slice").fetchone()
    c.close()
    assert done == 1, (
        "with the service answering, a slice whose remaining texts are "
        "permanently broken must be able to retire; otherwise it is "
        "unretirable forever, which is what blocked #913")
    # `seen` is CUMULATIVE across passes and `hits` is the service's own
    # total, so a re-walked slice legitimately exceeds it -- the live MED 2014
    # sits at seen=368,562 against hits=274,899. The certification rule is
    # "not walked SHORT", which is what `incomplete_slices` tests.
    assert seen >= hits, (
        f"a retired slice must not be walked short: seen={seen} hits={hits}")
    c = _s3.connect(tmp_path / "state.sqlite")
    try:
        assert not fx.incomplete_slices(c), (
            "the slice retired without being certifiable, which is the #913 "
            "defect itself")
    finally:
        c.close()


def test_a_successful_fetch_in_the_same_pass_supplies_the_control(
        tmp_path, monkeypatch):
    """END TO END, with the real probe rather than a patched verdict.

    This is what the mutation "stop recording controls" needed in order to be
    caught: without it the buffer is empty, `_service_is_up` returns None and
    the run aborts, so a test that patches the verdict cannot tell the
    difference. Here the retry pass ALSO fetches fresh records successfully, so
    the buffer fills, the real probe finds a control and the slice retires.

    The fresh ids come FIRST in the page deliberately. Records are processed in
    page order, so a buffer filled only after the breaker trips would not help
    -- and that ordering dependence is real rather than an artifact: a pass
    that has not yet fetched anything genuinely has no control, and aborting is
    the conservative branch there.
    """
    n = 2 * fx.MAX_CONSECUTIVE_DEFERRALS + 20
    doomed_ids = [f"r{i}" for i in range(n)]
    doomed = {i for i in doomed_ids if int(i[1:]) % 2 == 0}
    _drive(tmp_path, monkeypatch, [_page(doomed_ids, None)], fail_ids=doomed)

    fresh = [f"ok{i}" for i in range(3)]
    _rerun(tmp_path, monkeypatch, [_page(fresh + doomed_ids, None)], doomed)

    import sqlite3 as _s3
    c = _s3.connect(tmp_path / "state.sqlite")
    try:
        done = c.execute("SELECT done FROM slice").fetchone()[0]
        assert done == 1, (
            "a pass that successfully fetched text has a control to probe, so "
            "the permanently broken documents must not stop it")
        assert not fx.incomplete_slices(c)
    finally:
        c.close()


def test_the_probe_does_not_build_its_own_request():
    """STRUCTURAL, and it exists because the first version DID.

    `_service_is_up` called `_get` with a URL it assembled itself. That
    duplicated construction the fetcher owns, so the probe could ask a
    different question than the fetcher asks -- and it bypassed the one seam
    the tests patch, so the guards for this function issued LIVE requests to
    Europe PMC with synthetic PMCIDs. The tell was the runtime: 22.2s for this
    file, against 1.5s once the probe went through `fetch_fulltext`.
    """
    src = (REPO / "scripts" / "corpus_expand_fetch.py").read_text()
    fn = next(n for n in ast.walk(ast.parse(src))
              if isinstance(n, ast.FunctionDef) and n.name == "_service_is_up")
    called = {n.func.id for n in ast.walk(fn)
              if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)}
    assert "fetch_fulltext" in called, (
        "the probe must re-run the same call that succeeded")
    assert "_get" not in called, (
        "the probe must not assemble its own request: that duplicates the "
        "fetcher's URL construction and escapes the test seam, which made "
        "these guards hit the live network")
    assert "urlopen" not in ast.dump(fn), "no direct HTTP from the probe"


def test_controls_come_from_this_run_first_then_the_index(tmp_path):
    """The fallback is what makes the probe usable, and the ORDER matters.

    A retry pass fetches only the records that already failed, so it can reach
    the breaker having succeeded at nothing. With no control the honest answer
    is "unknown" and the run aborts -- which is safe, and is also precisely the
    state that left #913 unable to finish: a slice whose remaining texts are
    permanently broken could never retire no matter what the service was doing.

    The identity index already records `has_fulltext` per key, so it can name
    articles whose text this crawl HAS held. This run's own successes come
    first because they are the freshest evidence.
    """
    import sqlite3
    import corpus_identity_index as ix
    db = tmp_path / "id.sqlite"
    c = sqlite3.connect(db)
    c.executescript(ix.SCHEMA)
    # Real-SHAPED ids. An earlier version of this fixture used readable names
    # like "PMC_WITH_TEXT", which `norm_pmcid` rejects -- so once controls were
    # filtered to those that can actually produce a request, the fixture's own
    # ids stopped qualifying. The hardening caught the test, not the code.
    c.executemany("INSERT INTO held(kind,key,source,has_fulltext) VALUES (?,?,?,?)",
                  [("pmcid", "PMC111", "expand-MED", 1),
                   ("pmcid", "PMC222", "expand-MED", 0),
                   ("pmcid", "PMC333", "frozen", 1)])
    c.commit()

    # No in-run success: every control comes from the index, and only the row
    # that actually HAS text qualifies.
    got = fx._control_records([], c, 3)
    assert [r["pmcid"] for r in got] == ["PMC111"], (
        f"only an expand-* row with has_fulltext=1 is a valid control; got {got}")

    # An in-run success takes precedence and is used alone when it suffices.
    got = fx._control_records([{"pmcid": "PMC999"}], c, 1)
    assert got == [{"pmcid": "PMC999"}], (
        f"this run's own success is the freshest evidence; got {got}")

    # And a run with no index at all degrades to "unknown" rather than raising.
    assert fx._control_records([], None, 3) == []
    c.close()


def test_a_broken_probe_never_takes_the_run_down(tmp_path):
    """A diagnostic must not be the thing that kills a run holding days of
    work. A closed or malformed index returns no controls, which the caller
    reads as "unknown" and handles, rather than propagating."""
    import sqlite3
    c = sqlite3.connect(tmp_path / "empty.sqlite")
    c.close()                      # closed connection: every query raises
    assert fx._control_records([], c, 3) == []


def test_a_control_that_would_cause_no_request_is_not_a_control(tmp_path):
    """A probe that asks nothing must not report success.

    `fetch_fulltext` returns None WITHOUT touching the network when a record
    has no usable PMCID and is not a PPR, and None is read as "the service
    answered" -- correctly, since a 404 is an answer. So a control with no
    usable PMCID would report the service up having issued no request at all:
    the degenerate input satisfying the check, which is the vacuity shape this
    repo keeps finding.
    """
    import sqlite3
    import corpus_identity_index as ix
    assert fx._control_records([{"pmcid": None}, {"pmcid": ""}], None, 3) == [], (
        "records that cannot produce a request must be filtered out")
    assert fx._control_records([{"doi": "10.1/x"}], None, 3) == [], (
        "a record with no pmcid at all is not a usable control")
    # And the index side is filtered on the same rule.
    db = tmp_path / "id.sqlite"
    c = sqlite3.connect(db)
    c.executescript(ix.SCHEMA)
    c.executemany("INSERT INTO held(kind,key,source,has_fulltext) VALUES (?,?,?,?)",
                  [("pmcid", "", "expand-MED", 1),
                   ("pmcid", "PMC77", "expand-MED", 1)])
    c.commit()
    got = [r["pmcid"] for r in fx._control_records([], c, 3)]
    assert got == ["PMC77"], f"the unusable key must be dropped; got {got}"
    c.close()
    # The whole point: with no usable control the verdict is UNKNOWN, which the
    # caller treats as unsafe -- NOT "up".
    assert fx._service_is_up([{"pmcid": None}], None) is None


# --- A vanished mount must not end the run (#966) --------------------------
#
# Two #913 passes died in Shards.write, 96 hours and 404 minutes in, both on a
# flush to a network mount that had gone away (Errno 5, then Errno 6). Neither
# failure said anything about the record: the bytes were fine and the
# destination was not there.


class _Flaky:
    """A gzip handle whose flush fails for the first N calls, like a mount
    that goes away and comes back."""

    def __init__(self, fails, errno_):
        self.fails, self.errno = fails, errno_
        self.written, self.closed, self.flushes = [], False, 0

    def write(self, s):
        self.written.append(s)

    def flush(self):
        self.flushes += 1
        if self.fails > 0:
            self.fails -= 1
            raise OSError(self.errno, "boom")

    def close(self):
        self.closed = True


def _shards_with(monkeypatch, tmp_path, handles):
    """A Shards whose opens hand back the given fakes in order."""
    s = fx.Shards(tmp_path, "t")
    seq = list(handles)
    opened = []

    def fake_open(path, *a, **k):
        opened.append(pathlib.Path(path).name)
        return seq.pop(0)

    monkeypatch.setattr(fx.gzip, "open", fake_open)
    monkeypatch.setattr(fx.time, "sleep", lambda *_: None)
    monkeypatch.setattr(fx, "STORAGE_POLL", 0.0)
    return s, opened


@pytest.mark.parametrize("err", [errno.EIO, errno.ENXIO, errno.ENODEV,
                                 errno.ESTALE, errno.ENOENT])
def test_a_vanished_mount_is_survived_and_the_record_still_lands(
        tmp_path, monkeypatch, err):
    """Every errno a disappearing network mount produces. EIO and ENXIO are
    the two that actually killed passes; the others are the same event with a
    different server or a different moment."""
    dead, live = _Flaky(1, err), _Flaky(0, err)
    s, opened = _shards_with(monkeypatch, tmp_path, [dead, live])
    s.write({"pmid": "1"})
    assert live.written, f"the record must land after the storage returns: {opened}"
    assert "1" in live.written[0]
    assert len(opened) == 2, f"a fresh shard must be opened; opened {opened}"
    assert opened[0] != opened[1], (
        f"the retry must not reuse the damaged shard's name: {opened}")


def test_the_dead_handle_is_discarded_rather_than_closed(tmp_path, monkeypatch):
    """Closing it would FLUSH, and flushing is what just failed -- so closing
    turns a recoverable drop into the same fatal error one line later."""
    dead, live = _Flaky(1, errno.EIO), _Flaky(0, errno.EIO)
    s, _ = _shards_with(monkeypatch, tmp_path, [dead, live])
    s.write({"pmid": "1"})
    assert not dead.closed, (
        "the handle whose flush failed must not be closed; its flush would "
        "raise again and re-kill the run")


def test_a_write_error_that_is_not_a_lost_mount_propagates(tmp_path, monkeypatch):
    """EACCES says the write is wrong, not that the destination left. Retrying
    it would loop on a real problem and hide it."""
    bad = _Flaky(1, errno.EACCES)
    s, _ = _shards_with(monkeypatch, tmp_path, [bad])
    with pytest.raises(OSError) as e:
        s.write({"pmid": "1"})
    assert e.value.errno == errno.EACCES


def test_storage_that_never_returns_still_ends_the_run(tmp_path, monkeypatch):
    """Bounded, not infinite. Giving up is safe rather than lossy -- the
    ledger is local and every slice resumes from its own cursor -- while
    waiting forever would hide a dead NAS behind a running process."""
    dead = _Flaky(99, errno.EIO)
    s, _ = _shards_with(monkeypatch, tmp_path, [dead])
    monkeypatch.setattr(fx, "STORAGE_RETURN_TIMEOUT", 0.0)
    with pytest.raises(OSError):
        s.write({"pmid": "1"})


def test_only_one_retry_is_attempted(tmp_path, monkeypatch):
    """A second failure after the storage came back is not a blip."""
    dead, also_dead = _Flaky(1, errno.EIO), _Flaky(1, errno.EIO)
    s, _ = _shards_with(monkeypatch, tmp_path, [dead, also_dead])
    with pytest.raises(OSError):
        s.write({"pmid": "1"})


def test_a_failed_write_never_reports_success_to_its_caller(
        tmp_path, monkeypatch):
    """THE INVARIANT THAT MAKES ALL OF THIS SAFE, asserted rather than assumed.

    The identity index is trustworthy because the caller indexes a record only
    after `write` RETURNS. So the property this change must preserve is that a
    record whose bytes never landed makes `write` raise -- on the first
    attempt and on the retry alike. If it ever returned after a failed write,
    the row would claim an article whose bytes are gone, which is the
    permanent-invisibility failure the ordering exists to prevent.
    """
    s, _ = _shards_with(monkeypatch, tmp_path, [_Flaky(99, errno.EIO)])
    monkeypatch.setattr(fx, "STORAGE_RETURN_TIMEOUT", 0.0)
    raised = False
    try:
        s.write({"pmid": "1"})
    except OSError:
        raised = True
    assert raised, (
        "write returned despite the record never reaching storage; the caller "
        "would now index an article whose bytes do not exist")


def test_the_classifier_is_a_function_so_it_can_be_mutated_and_caught():
    """STRUCTURAL, for the reason `_should_abort` is a function: a condition
    written inline can only be checked by grepping for it, and a grep passes
    while the comparison sits behind `if False:`."""
    assert fx._storage_vanished(OSError(errno.EIO, "x")) is True
    assert fx._storage_vanished(OSError(errno.EACCES, "x")) is False
    src = (REPO / "scripts" / "corpus_expand_fetch.py").read_text()
    cls = next(n for n in ast.walk(ast.parse(src))
               if isinstance(n, ast.ClassDef) and n.name == "Shards")
    called = {n.func.id for n in ast.walk(cls)
              if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)}
    assert "_storage_vanished" in called, (
        "Shards must consult the classifier rather than testing errnos inline")


def test_close_tolerates_a_lost_mount_but_not_a_real_error(tmp_path, monkeypatch):
    """`close` has its own classifier and it was UNGUARDED -- the mutation that
    made it swallow every error survived the whole sweep.

    The two cases are not the same. A mount that went away between the last
    write and the close leaves a damaged shard the audit reports by name, and
    raising on the way out would discard a run's totals for nothing. Any other
    error is a real problem and must not be hidden by a code path that exists
    for storage loss.
    """
    class _Closes:
        def __init__(self, err):
            self.err = err

        def write(self, s):
            pass

        def flush(self):
            pass

        def close(self):
            raise OSError(self.err, "on close")

    for err, should_raise in ((errno.EIO, False), (errno.EACCES, True)):
        s = fx.Shards(tmp_path, f"t{err}")
        monkeypatch.setattr(fx.gzip, "open", lambda *a, **k: _Closes(err))
        s.write({"pmid": "1"})
        if should_raise:
            with pytest.raises(OSError) as e:
                s.close()
            assert e.value.errno == errno.EACCES, (
                "a real close error must not be swallowed by the path that "
                "exists for a lost mount")
        else:
            s.close()          # tolerated
            assert s.fh is None
