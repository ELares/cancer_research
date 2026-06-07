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

Those two locations are the *only* counted scope, and that scope is **enforced**:
a Python test file added anywhere else fails this test loudly (rather than being
silently uncounted by both the docs and this guard) — see the scope assertion.

To update after intentionally adding/removing Python tests:
  python -m pytest tests/ --collect-only -q                       # → EXPECTED_TESTS_DIR
  python -m pytest simulations/ferroptosis-python/test_bindings.py --collect-only -q  # → EXPECTED_BINDINGS
then set the constants below AND the "NNN Python tests" number in CLAUDE.md and
README.md (EXPECTED_TESTS_DIR + EXPECTED_BINDINGS).
"""

import os
import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]

# Collected-test expectations. Bump these (and the doc strings) in lockstep with
# intentional test additions/removals — that lockstep is the whole point.
EXPECTED_TESTS_DIR = 278  # `tests/` total, including this meta-test
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


# Build/VCS/cache dirs pruned from the repo-wide test-file scan (skip the big
# ones — `target/`, `.git/` — both for speed and to avoid vendored test files).
_PRUNE_DIRS = {
    "__pycache__",
    ".git",
    ".pytest_cache",
    "target",
    "build",
    "dist",
    "node_modules",
    ".venv",
    "venv",
    "site-packages",
}


def _all_python_test_files() -> "list[Path]":
    """Every `test_*.py` / `*_test.py` under the repo, pruning build/VCS/cache
    dirs. Used to enforce that the counted scope (`tests/` + the one bindings
    file) is exhaustive — a test added in a new location must not slip the count."""
    found: "list[Path]" = []
    for dirpath, dirnames, filenames in os.walk(REPO_ROOT):
        dirnames[:] = [d for d in dirnames if d not in _PRUNE_DIRS]
        for fn in filenames:
            if fn.endswith(".py") and (fn.startswith("test_") or fn.endswith("_test.py")):
                found.append(Path(dirpath) / fn)
    return found


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

    # 3. Both docs must advertise the total (tests/ + bindings) — at EVERY
    #    "NNN Python tests" occurrence, not just the first.
    for doc_path, label in DOC_FILES:
        text = doc_path.read_text(encoding="utf-8")
        found = _DOC_COUNT_RE.findall(text)
        assert found, f"{label}: no 'NNN Python tests' string found"
        for documented_n in map(int, found):
            assert documented_n == EXPECTED_TOTAL, (
                f"{label} says {documented_n} Python tests, but the collected total is "
                f"{EXPECTED_TOTAL} ({EXPECTED_TESTS_DIR} in tests/ + {EXPECTED_BINDINGS} bindings). "
                f"Update that doc line to {EXPECTED_TOTAL}."
            )

    # 4. The counted scope must be exhaustive: no Python test file outside
    #    `tests/` or the one bindings file (else it would be silently uncounted
    #    by both the docs and this guard). A new test location fails loudly here.
    tests_dir = REPO_ROOT / "tests"
    for path in _all_python_test_files():
        in_scope = path == BINDINGS_FILE or tests_dir in path.parents
        assert in_scope, (
            f"Python test file outside the counted scope: {path.relative_to(REPO_ROOT)}. "
            f"This guard counts only tests/ + {BINDINGS_FILE.relative_to(REPO_ROOT)} — "
            f"extend tests/test_doc_test_count.py (and the docs) to include the new location."
        )
