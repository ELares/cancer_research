"""The page-retry contract, tested once for every fetcher that imports it.

THE DEFECT THIS EXISTS TO STOP RECURRING. A paged search that raises on a
transient failure ends the whole crawl. It happened three times:

  Europe PMC full text  -- one 503 ended a run.
  Europe PMC search     -- one 503 ended a two-day run at 302,973 records.
  OpenAlex search       -- one 429 killed the fetcher on its first live run,
                           AFTER the Europe PMC fix, because the fix was
                           written in one place and not carried across.

The third is the reason the retry is shared rather than copied, and the reason
this file tests the SHARED function plus the fact that each fetcher routes
through it -- a fetcher that stops calling it is the same bug returning.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
import corpus_expand_fetch as fx  # noqa: E402
import openalex_fetch as oaf  # noqa: E402
import trials_fetch as tf  # noqa: E402

FETCHERS = (fx, oaf, tf)


def test_a_briefly_unavailable_page_is_retried(monkeypatch):
    monkeypatch.setattr(fx.time, "sleep", lambda *_: None)
    calls = {"n": 0}

    def flaky(x):
        calls["n"] += 1
        if calls["n"] < 3:
            raise fx.TransientFetchError("503")
        return {"ok": x}

    assert fx.retry_page(flaky, "ARG", verbose=False) == {"ok": "ARG"}
    assert calls["n"] == 3


def test_the_same_arguments_are_replayed(monkeypatch):
    """Retrying a DIFFERENT cursor would skip a page. A cursor is only advanced
    after its page is fully processed, so the retry must replay exactly."""
    monkeypatch.setattr(fx.time, "sleep", lambda *_: None)
    seen = []

    def flaky(cursor, extra):
        seen.append((cursor, extra))
        if len(seen) < 3:
            raise fx.TransientFetchError("429")
        return {}

    fx.retry_page(flaky, "CURSOR-A", ",filter", verbose=False)
    assert seen == [("CURSOR-A", ",filter")] * 3, seen


def test_a_sustained_outage_still_stops_the_run(monkeypatch):
    """An empty page reads as 'no more results', so swallowing the failure
    would record an outage as a completed pass."""
    monkeypatch.setattr(fx.time, "sleep", lambda *_: None)

    def dead(*a, **k):
        raise fx.TransientFetchError("503")

    try:
        fx.retry_page(dead, label="x", verbose=False)
    except RuntimeError as e:
        assert "looks down" in str(e)
        assert not isinstance(e, fx.TransientFetchError), (
            "a sustained outage must not be catchable as a per-item failure")
        return
    raise AssertionError("a permanently failing page returned quietly")


def test_a_non_transient_error_is_not_retried(monkeypatch):
    """Only TransientFetchError means 'try again'. Retrying a KeyError four
    times just delays the traceback."""
    monkeypatch.setattr(fx.time, "sleep", lambda *_: None)
    calls = {"n": 0}

    def broken(*a, **k):
        calls["n"] += 1
        raise KeyError("schema changed")

    try:
        fx.retry_page(broken, verbose=False)
    except KeyError:
        assert calls["n"] == 1, f"a KeyError was retried {calls['n']} times"
        return
    raise AssertionError("a programming error was swallowed")


def test_every_paged_fetcher_routes_through_the_shared_retry():
    """THE POINT OF SHARING IT. The Europe PMC fix was written once and not
    carried across, and OpenAlex died on its first 429 as a result. A fetcher
    that stops calling this is that bug returning.
    """
    for mod in FETCHERS:
        src = Path(mod.__file__).read_text()
        assert "retry_page(" in src, (
            f"{Path(mod.__file__).name} pages a search without the shared "
            "retry; one transient failure will end its crawl")


def test_no_fetcher_reimplements_the_retry():
    """A local copy is how the three drift apart again."""
    for mod in (oaf, tf):
        src = Path(mod.__file__).read_text()
        assert "def retry_page(" not in src
        assert "SEARCH_PAGE_BACKOFF = " not in src


def test_the_backoff_outlasts_a_brief_interruption():
    assert fx.SEARCH_PAGE_ATTEMPTS >= 3
    assert sum(fx.SEARCH_PAGE_BACKOFF) >= 120
