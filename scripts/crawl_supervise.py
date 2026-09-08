#!/usr/bin/env python3
"""Keep the expansion crawl running across the failures that end it.

WHY A SUPERVISOR AND NOT JUST A LONGER RETRY
--------------------------------------------
The crawl is a job measured in weeks against a service that is occasionally
unavailable, and it has now been stopped twice by single transient failures:
once by a 503 on one full-text URL, and once -- after 302,973 records and two
days -- by a 503 on one search page. Both are fixed at the point they happened.
Neither fix generalises to the next thing, because the failure that ends the
next run will be one nobody has seen yet: a disk hiccup, an OOM, a laptop
sleeping, a Python-level bug in a branch reached after a million records.

The lesson from those two is not "add another retry". It is that an unattended
job needs something whose ONLY responsibility is noticing it stopped.

WHAT THIS DOES NOT DO
---------------------
It does not retry forever, and it does not restart a crawl that FINISHED.
`corpus_expand_fetch` exits 0 when every slice is done, and a supervisor that
treats success as a reason to start again is a busy loop against a public API.
Exit 0 ends the supervision, and so does a run of failures fast enough to mean
the crawler cannot start at all -- a syntax error restarted forever is worse
than one reported once.

THE RESTART IS SAFE BECAUSE THE CRAWLER MAKES IT SAFE, not because this script
is careful. Every record's bytes are flushed before its identity is indexed,
progress is committed per page, and a slice that lost a record is rewound
rather than retired. This can start the crawler as often as it likes; the worst
case is re-reading a page.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

# Every long-running fetcher in this project, by name. They share the same
# durability contract -- write before index, per-page commits, a resumable
# cursor -- which is exactly what makes restarting any of them safe, so they
# can share the supervisor too. Adding a fetcher here rather than writing a
# second supervisor is the same reasoning that made the page retry shared:
# three copies of a thing drift, and the drift is invisible until one of them
# is the copy that failed.
TARGETS = {
    "expand": REPO / "scripts" / "corpus_expand_fetch.py",
    "openalex": REPO / "scripts" / "openalex_fetch.py",
    "trials": REPO / "scripts" / "trials_fetch.py",
}
CRAWLER = TARGETS["expand"]

# Waits between restarts, in seconds. Escalating, because the same failure
# twice in a row usually means the service needs longer than the same failure
# once. The last value repeats for every subsequent restart.
BACKOFF = (30, 120, 300, 900)
# A run shorter than this did not do any work -- the crawler could not start,
# or died immediately. Restarting THAT is a busy loop, not supervision.
MIN_USEFUL_RUN = 60.0
MAX_IMMEDIATE_FAILURES = 3


def _backoff(n: int) -> int:
    return BACKOFF[min(n, len(BACKOFF) - 1)]


def supervise(max_restarts: int | None = None, env: dict | None = None,
              runner=None, log=print, target: str = "expand") -> dict:
    """Run the crawler until it finishes, or until restarting is pointless.

    `runner` exists so the decision logic can be tested without a crawl: it is
    called with the environment and returns (exit_code, seconds_elapsed).
    """
    runner = runner or _run_crawler
    stats = {"runs": 0, "restarts": 0, "immediate_failures": 0, "reason": None}
    while True:
        stats["runs"] += 1
        code, elapsed = runner(env or {}, target)
        if code == 0:
            stats["reason"] = "crawl finished"
            log(f"crawl finished after {stats['runs']} run(s)")
            return stats

        if elapsed < MIN_USEFUL_RUN:
            stats["immediate_failures"] += 1
        else:
            # A run that did real work resets the counter: a crash after six
            # hours is a different event from one after six seconds, and
            # counting them together would retire a healthy crawl.
            stats["immediate_failures"] = 0

        if stats["immediate_failures"] >= MAX_IMMEDIATE_FAILURES:
            stats["reason"] = (
                f"{stats['immediate_failures']} runs failed within "
                f"{MIN_USEFUL_RUN:.0f}s; the crawler cannot start")
            log(f"giving up: {stats['reason']}")
            return stats

        if max_restarts is not None and stats["restarts"] >= max_restarts:
            stats["reason"] = f"reached --max-restarts={max_restarts}"
            log(f"stopping: {stats['reason']}")
            return stats

        wait = _backoff(stats["restarts"])
        stats["restarts"] += 1
        log(f"crawl exited {code} after {elapsed:.0f}s; "
            f"restart {stats['restarts']} in {wait}s")
        time.sleep(wait)


def _run_crawler(extra_env: dict, target: str = "expand") -> tuple[int, float]:
    script = TARGETS[target]
    env = {**os.environ, **extra_env}
    t = time.time()
    p = subprocess.run([sys.executable, str(script)], env=env)
    return p.returncode, time.time() - t


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target", default="expand", choices=sorted(TARGETS),
                    help="which fetcher to supervise")
    ap.add_argument("--workers", default=None)
    ap.add_argument("--sleep", default=None)
    ap.add_argument("--max-restarts", type=int, default=None)
    a = ap.parse_args()
    env = {}
    if a.workers:
        env["FERRO_EXPAND_WORKERS"] = str(a.workers)
    if a.sleep:
        env["FERRO_EXPAND_SLEEP"] = str(a.sleep)

    def stamped(msg):
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)

    s = supervise(max_restarts=a.max_restarts, env=env, log=stamped,
                  target=a.target)
    stamped(f"supervisor done: {s}")
    return 0 if s["reason"] == "crawl finished" else 1


if __name__ == "__main__":
    raise SystemExit(main())
