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
    assert "_text_available(rec) and _text_permitted(rec)" in src, (
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


def test_a_lookalike_host_is_not_mistaken_for_europe_pmc():
    """The substring form of this check accepted every one of these. A host
    test that a lookalike domain passes is not a host test."""
    for url in ("https://www.ebi.ac.uk.attacker.test/x",
                "https://evil.example/?ref=www.ebi.ac.uk",
                "https://www.ebi.ac.uk@evil.example/x",
                "https://notwww.ebi.ac.uk.co/x"):
        assert _fetch_hosts(f'x = "{url}"') != {"www.ebi.ac.uk"}, url
        assert "www.ebi.ac.uk" in url, "the substring check would have passed this"


def test_a_record_with_no_identifier_is_refused():
    """Writing a record we could never recognise again would guarantee we
    re-download it, which is the one thing this fetcher must not do."""
    src = Path(fx.__file__).read_text()
    assert 'if not (pmid or pmcid or doi or rec.get("id")):' in src
    assert "continue  # nothing to dedup on later; refuse it" in src


def test_every_written_record_is_indexed_in_the_same_pass():
    """A crash between writing and indexing would re-fetch the record. The
    index write must sit inside the same loop as the shard write."""
    src = Path(fx.__file__).read_text()
    body = src[src.index("shards.write(out)"):src.index("if limit_new and totals")]
    assert "INSERT INTO held" in body, (
        "a record is written to a shard without being added to the index in "
        "the same pass, so an interrupted run will re-download it")
    assert "MAX(held.has_fulltext" in body


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
