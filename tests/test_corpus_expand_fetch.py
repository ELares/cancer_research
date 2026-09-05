"""Guards for the expansion fetcher.

Two requirements drive every test here, and both were stated as one sentence
each: never download an article we already have, and do not take paywalled
content. The first is enforced by asking the identity index before fetching
and updating it in the same pass; the second by asking a source that refuses,
rather than by this code deciding what is free.
"""
import ast
import gzip
import json
import sqlite3
import sys
import urllib.parse
from pathlib import Path

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
    body = src[src.index("stopped_early = False"):src.index("st.commit()")]
    assert "stopped_early = True" in body, "the early exit is not recorded"
    assert "cursor if stopped_early else" in body, (
        "the cursor advances even when the page was abandoned part-way")
    assert "len(hits)" not in body.split("UPDATE slice SET cursor")[1], (
        "the whole page is still credited as seen after a partial pass")


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
    row = c.execute("SELECT cursor, seen, kept, fulltext, deferred, done "
                    "FROM slice").fetchone()
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
                             "deferred", "done"), row)), stored


def _page(ids, nxt):
    return {"resultList": {"result": [
        {"id": i, "pmid": None, "pmcid": f"PMC{abs(hash(i)) % 10**7}",
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
    """Gated on request STARTS, so the rate is the same whatever the latency
    does -- workers overlap their waiting, not their requests."""
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
