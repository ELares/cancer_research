"""Guards for the crawl supervisor.

A supervisor is a small amount of code with two ways to be actively harmful:
restarting something that finished (a busy loop against a public API) and
restarting something that cannot start (the same, faster). Both look like
"keeping the crawl going" from the outside.

The decision logic is separated from the subprocess call precisely so these can
run without a crawl.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
import crawl_supervise as cs  # noqa: E402


def _runner(results):
    """Returns (exit_code, elapsed) from a canned list, then 0 forever."""
    seq = list(results)

    def run(env):
        return seq.pop(0) if seq else (0, 3600.0)
    return run


def test_a_finished_crawl_is_not_restarted():
    """`corpus_expand_fetch` exits 0 when every slice is done. Restarting that
    is a busy loop against a public API, dressed as diligence."""
    calls = []

    def run(env):
        calls.append(1)
        return (0, 10.0)

    s = cs.supervise(runner=run, log=lambda *_: None)
    assert len(calls) == 1
    assert s["reason"] == "crawl finished" and s["restarts"] == 0


def test_a_crash_after_real_work_is_restarted(monkeypatch):
    """The case this exists for: two multi-day runs were ended by single
    transient failures, one of them after 302,973 records."""
    monkeypatch.setattr(cs.time, "sleep", lambda *_: None)
    s = cs.supervise(runner=_runner([(1, 7200.0), (1, 3600.0)]),
                     log=lambda *_: None)
    assert s["restarts"] == 2 and s["reason"] == "crawl finished"


def test_a_crawler_that_cannot_start_is_not_restarted_forever(monkeypatch):
    """A syntax error restarted every 30 seconds is worse than one reported
    once: the failure is hidden by the retrying."""
    monkeypatch.setattr(cs.time, "sleep", lambda *_: None)
    calls = []

    def run(env):
        calls.append(1)
        return (1, 0.4)          # died immediately, every time

    s = cs.supervise(runner=run, log=lambda *_: None)
    assert len(calls) == cs.MAX_IMMEDIATE_FAILURES
    assert "cannot start" in s["reason"]


def test_a_long_run_resets_the_immediate_failure_count(monkeypatch):
    """A crash after six hours is a different event from one after six
    seconds. Counting them together would retire a healthy crawl."""
    monkeypatch.setattr(cs.time, "sleep", lambda *_: None)
    # two instant failures, then a real run, then two more instant failures:
    # never three in a row, so it must keep going and reach the finish.
    s = cs.supervise(
        runner=_runner([(1, 0.5), (1, 0.5), (1, 9000.0), (1, 0.5), (1, 0.5)]),
        log=lambda *_: None)
    assert s["reason"] == "crawl finished", s
    assert s["restarts"] == 5


def test_the_backoff_escalates_and_then_holds():
    """The same failure twice usually needs longer than the same failure once,
    and an unbounded escalation would stop supervising anything."""
    waits = [cs._backoff(i) for i in range(8)]
    assert waits == sorted(waits), waits
    assert waits[0] < waits[-1]
    assert waits[-1] == waits[len(cs.BACKOFF) - 1], "the last value must repeat"
    assert waits[-1] <= 3600, "a wait this long is not supervision"


def test_max_restarts_is_honoured(monkeypatch):
    monkeypatch.setattr(cs.time, "sleep", lambda *_: None)
    s = cs.supervise(max_restarts=2, runner=_runner([(1, 100.0)] * 10),
                     log=lambda *_: None)
    assert s["restarts"] == 2 and "max-restarts" in s["reason"]


def test_the_supervisor_does_not_reimplement_the_crawl():
    """Its ONLY responsibility is noticing the crawl stopped. Every safety
    property -- flush before index, per-page commits, rewinding a lossy slice --
    belongs to the crawler, which is why restarting is safe at all."""
    src = Path(cs.__file__).read_text()
    for forbidden in ("INSERT INTO held", "gzip.open", "cursor", "shards"):
        assert forbidden not in src, (
            f"the supervisor touches {forbidden!r}; it must only start and "
            "watch the crawler")
