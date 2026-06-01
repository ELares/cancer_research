"""Pin the documented Python test count against the live collected count (#260).

CLAUDE.md and README.md both advertise "NNN Python tests". That number was
re-derived by hand in #257 and has drifted before. This test fails — in the
normal `pytest tests/` CI run, no extra workflow step — whenever the collected
count diverges from the documented one, and the failure message says exactly
which constant + doc lines to update.

Breakdown: `EXPECTED_TESTS_DIR` tests under `tests/` (including this one) plus
`EXPECTED_BINDINGS` under `simulations/ferroptosis-python/test_bindings.py`. The
binding tests need the compiled `ferroptosis_core` extension, which the Python
CI does NOT build — so they are accounted for as a constant and only re-counted
when the extension happens to be importable (never flakes when it isn't).

To update after intentionally adding/removing Python tests:
  python -m pytest tests/ --collect-only -q                       # → EXPECTED_TESTS_DIR
  python -m pytest simulations/ferroptosis-python/test_bindings.py --collect-only -q  # → EXPECTED_BINDINGS
then set the constants below AND the "NNN Python tests" number in CLAUDE.md and
README.md (EXPECTED_TESTS_DIR + EXPECTED_BINDINGS).
"""

import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]

# Collected-test expectations. Bump these (and the doc strings) in lockstep with
# intentional test additions/removals — that lockstep is the whole point.
EXPECTED_TESTS_DIR = 94  # `tests/` total, including this meta-test
EXPECTED_BINDINGS = 12  # `simulations/ferroptosis-python/test_bindings.py`
EXPECTED_TOTAL = EXPECTED_TESTS_DIR + EXPECTED_BINDINGS  # the documented "NNN Python tests"

BINDINGS_FILE = REPO_ROOT / "simulations" / "ferroptosis-python" / "test_bindings.py"
# (doc file, human label) pairs whose "NNN Python tests" must equal EXPECTED_TOTAL.
DOC_FILES = [
    (REPO_ROOT / "CLAUDE.md", "CLAUDE.md ('NNN Python tests' in the simulation-suite line)"),
    (REPO_ROOT / "README.md", "README.md ('NNN Python tests' in the tests/ table row)"),
]

_COLLECTED_RE = re.compile(r"(\d+) tests? collected")
_DOC_COUNT_RE = re.compile(r"(\d+) Python tests")


def _collect_count(target: Path) -> "int | None":
    """Return the number of tests pytest collects from `target`, or None if
    collection failed (e.g. the compiled binding is unavailable) — so the
    binding count never flakes the suite."""
    res = subprocess.run(
        [sys.executable, "-m", "pytest", str(target), "--collect-only", "-q", "-p", "no:cacheprovider"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    match = _COLLECTED_RE.search(res.stdout)
    if res.returncode != 0 or match is None:
        return None
    return int(match.group(1))


def test_documented_python_test_count_matches_collection():
    # 1. tests/ must collect exactly the documented non-binding count.
    tests_n = _collect_count(REPO_ROOT / "tests")
    assert tests_n is not None, "could not collect tests/ — is pytest set up?"
    assert tests_n == EXPECTED_TESTS_DIR, (
        f"tests/ collects {tests_n}, expected {EXPECTED_TESTS_DIR}. If this was "
        f"intentional, set EXPECTED_TESTS_DIR = {tests_n} in tests/test_doc_test_count.py "
        f"and bump the 'NNN Python tests' number to {tests_n + EXPECTED_BINDINGS} in CLAUDE.md + README.md."
    )

    # 2. The binding tests only when the compiled extension is importable (CI
    #    does not build it); otherwise accounted for as a constant, no flake.
    bindings_n = _collect_count(BINDINGS_FILE)
    if bindings_n is not None:
        assert bindings_n == EXPECTED_BINDINGS, (
            f"{BINDINGS_FILE.name} collects {bindings_n}, expected {EXPECTED_BINDINGS}. "
            f"Set EXPECTED_BINDINGS = {bindings_n} and bump 'NNN Python tests' to "
            f"{EXPECTED_TESTS_DIR + bindings_n} in CLAUDE.md + README.md."
        )

    # 3. Both docs must advertise the total (tests/ + bindings).
    for doc_path, label in DOC_FILES:
        text = doc_path.read_text(encoding="utf-8")
        documented = _DOC_COUNT_RE.search(text)
        assert documented is not None, f"{label}: no 'NNN Python tests' string found"
        documented_n = int(documented.group(1))
        assert documented_n == EXPECTED_TOTAL, (
            f"{label} says {documented_n} Python tests, but the collected total is "
            f"{EXPECTED_TOTAL} ({EXPECTED_TESTS_DIR} in tests/ + {EXPECTED_BINDINGS} bindings). "
            f"Update that doc line to {EXPECTED_TOTAL}."
        )
