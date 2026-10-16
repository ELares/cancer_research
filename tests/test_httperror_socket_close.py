"""An HTTPError IS the response, so every handler must release it (#954).

`urlopen` RAISES `urllib.error.HTTPError` instead of returning it, so a
`with urllib.request.urlopen(...) as r:` never binds and its `__exit__` never
runs. `HTTPError` subclasses `addinfourl` and wraps the open socket, so an
error reply leaves that socket in CLOSE_WAIT until the garbage collector
happens to reclaim it.

Two guards, and the structural one is the load-bearing half: the behavioural
test can only cover the handlers that exist today, while a NEW unguarded
handler in a new fetcher is exactly how this class comes back.
"""

import ast
import io
import sys
import urllib.error
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "scripts"
sys.path.insert(0, str(SCRIPTS))

import corpus_expand_fetch as cef  # noqa: E402


# --------------------------------------------------------------------------
# structural: every handler in the repo releases the response
# --------------------------------------------------------------------------

def _httperror_handlers(tree: ast.AST):
    """Yield (lineno, bound_name, body) for each `except HTTPError as e:`."""
    for node in ast.walk(tree):
        if not isinstance(node, ast.ExceptHandler) or node.name is None:
            continue
        kind = node.type
        names = []
        if isinstance(kind, ast.Tuple):
            names = list(kind.elts)
        elif kind is not None:
            names = [kind]
        for n in names:
            if ast.unparse(n).endswith("HTTPError"):
                yield node.lineno, node.name, node.body
                break


def _releases(bound: str, body) -> bool:
    """Does this handler body call `<bound>.close()`?

    DELIBERATELY LITERAL, and that cuts both ways. A handler that releases the
    response some other way -- `with contextlib.closing(e):`, or a helper like
    `_release(e)` -- reads here as a leak, so this gate can produce a false
    positive on correct code. That is the accepted trade: the repo already
    learned from `audit_heading_only_assertions.py` that a gate at middling
    precision trains exemption-adding rather than fixing, so the rule is kept
    narrow and obvious enough that the fix for a false positive is to call
    `.close()` directly. If a second release idiom ever becomes common here,
    widen this rather than exempting the file.
    """
    for stmt in body:
        for node in ast.walk(stmt):
            if (isinstance(node, ast.Call)
                    and isinstance(node.func, ast.Attribute)
                    and node.func.attr == "close"
                    and isinstance(node.func.value, ast.Name)
                    and node.func.value.id == bound):
                return True
    return False


def _scan(path: Path):
    """Every unreleased handler in one file, as (lineno, name) pairs."""
    tree = ast.parse(path.read_text())
    return [(ln, nm) for ln, nm, body in _httperror_handlers(tree)
            if not _releases(nm, body)]


# The five sites this landed for, NAMED rather than counted: a count cannot
# tell a handler that was legitimately removed from a detector that stopped
# finding one, and those need different responses.
KNOWN_HANDLER_SITES = frozenset({
    "atlas_fulltext.py",
    "atlas_relations.py",
    "burden_data_feasibility.py",
    "corpus_expand_fetch.py",
    "preprint_access_check.py",
})


def test_every_httperror_handler_in_scripts_releases_the_response():
    """Inverse coverage: a NEW fetcher cannot land with an unclosed handler.

    The floor is load-bearing, not decoration. `assert not offenders` passes
    when the scan finds NOTHING, so a detector that silently stopped matching
    -- a renamed exception, an `ast` change, a bad edit to
    `_httperror_handlers` -- would read exactly like a clean repo. Measured:
    blinding the detector leaves `offenders` empty and that assertion green.
    """
    scanned = 0
    offenders = []
    for path in sorted(SCRIPTS.glob("*.py")):
        tree = ast.parse(path.read_text())
        handlers = list(_httperror_handlers(tree))
        scanned += len(handlers)
        offenders += [f"{path.name}:{ln} (`{nm}` never closed)"
                      for ln, nm, body in handlers if not _releases(nm, body)]
    assert scanned >= len(KNOWN_HANDLER_SITES), (
        f"the detector found {scanned} HTTPError handlers but at least "
        f"{len(KNOWN_HANDLER_SITES)} exist, so this scan is not reading what "
        "it claims to -- a broken detector looks identical to a clean repo")
    assert not offenders, (
        "an HTTPError IS the response; these handlers leak its socket:\n  "
        + "\n  ".join(offenders))


def test_the_detector_finds_every_site_this_landed_for():
    """Names them, so a vanished handler is distinguishable from a blind scan."""
    found = {p.name for p in SCRIPTS.glob("*.py")
             if any(_httperror_handlers(ast.parse(p.read_text())))}
    missing = KNOWN_HANDLER_SITES - found
    assert not missing, (
        f"no HTTPError handler detected in {sorted(missing)}. If one was "
        "deliberately removed, drop it from KNOWN_HANDLER_SITES; if not, the "
        "detector has stopped seeing it.")


def test_the_scanner_finds_a_planted_leak(tmp_path):
    """A scan reporting zero because it is BROKEN looks like a clean repo."""
    leaky = tmp_path / "leaky.py"
    leaky.write_text(
        "import urllib.error\n"
        "def f(u):\n"
        "    try:\n"
        "        pass\n"
        "    except urllib.error.HTTPError as e:\n"
        "        return e.code\n"
    )
    assert _scan(leaky) == [(5, "e")]


def test_the_scanner_accepts_a_handler_that_closes(tmp_path):
    """And it must not flag the fixed shape, or it would train exemptions."""
    clean = tmp_path / "clean.py"
    clean.write_text(
        "import urllib.error\n"
        "def f(u):\n"
        "    try:\n"
        "        pass\n"
        "    except urllib.error.HTTPError as e:\n"
        "        code = e.code\n"
        "        e.close()\n"
        "        return code\n"
    )
    assert _scan(clean) == []


def test_the_scanner_sees_a_handler_catching_a_tuple(tmp_path):
    """`except (URLError, HTTPError) as e` is the same handler."""
    leaky = tmp_path / "tup.py"
    leaky.write_text(
        "import urllib.error\n"
        "def f(u):\n"
        "    try:\n"
        "        pass\n"
        "    except (urllib.error.URLError, urllib.error.HTTPError) as e:\n"
        "        return e.code\n"
    )
    assert _scan(leaky) == [(5, "e")]


def test_the_scanner_reads_the_repository_and_not_a_fixture():
    """Point it at an empty tree and it must report nothing to find."""
    assert list((REPO / "nonexistent-tree").glob("*.py")) == []
    assert (SCRIPTS / "corpus_expand_fetch.py").exists()


# --------------------------------------------------------------------------
# behavioural: _get releases the response on all three exits
# --------------------------------------------------------------------------

class _Recording(urllib.error.HTTPError):
    """An HTTPError that records whether anyone released it."""

    def __init__(self, code, headers=None):
        self.closed = False
        super().__init__("http://example.invalid/x", code,
                         f"error {code}", headers or {}, io.BytesIO(b""))

    def close(self):
        self.closed = True
        super().close()


@pytest.fixture
def no_sleep(monkeypatch):
    """So a scheduled multi-hour wait does not actually wait."""
    slept = []
    monkeypatch.setattr(cef.time, "sleep", lambda s: slept.append(s))
    return slept


def _raise_each(monkeypatch, errors):
    """Make urlopen raise each error in turn, so retries see fresh ones."""
    it = iter(errors)

    def fake(req, timeout=None):
        raise next(it)

    monkeypatch.setattr(cef.urllib.request, "urlopen", fake)


def test_a_refused_status_releases_the_response(monkeypatch, no_sleep):
    """The `return None` path: refused means answered, not retried."""
    err = _Recording(404)
    _raise_each(monkeypatch, [err])
    assert cef._get("http://example.invalid/x", tries=4) is None
    assert err.closed, "a refused reply returned without closing the response"


def test_an_exhausted_retry_budget_releases_every_response(monkeypatch,
                                                           no_sleep):
    """The raising path: `last` must not carry a live socket out of the loop."""
    errs = [_Recording(500) for _ in range(3)]
    _raise_each(monkeypatch, errs)
    with pytest.raises(cef.TransientFetchError):
        cef._get("http://example.invalid/x", tries=3)
    unclosed = [i for i, e in enumerate(errs) if not e.closed]
    assert not unclosed, f"retry attempts {unclosed} left a response open"


def test_the_give_up_message_still_names_the_status(monkeypatch, no_sleep):
    """Releasing the response must not cost the diagnostic it carried."""
    _raise_each(monkeypatch, [_Recording(500) for _ in range(3)])
    with pytest.raises(cef.TransientFetchError) as exc:
        cef._get("http://example.invalid/x", tries=3)
    assert "500" in str(exc.value), (
        "keeping the text instead of the response lost the status")


def test_a_scheduled_long_wait_does_not_sleep_holding_the_socket(monkeypatch,
                                                                 no_sleep):
    """The worst path: it sleeps for HOURS, previously still holding it."""
    long_wait = {"Retry-After": "7200"}
    first = _Recording(503, long_wait)
    rest = [_Recording(503) for _ in range(3)]
    _raise_each(monkeypatch, [first] + rest)
    with pytest.raises(cef.TransientFetchError):
        cef._get("http://example.invalid/x", tries=3)
    assert first.closed, (
        "the response was still open across a scheduled multi-hour wait")
    assert 7200.0 in no_sleep, "the scheduled-wait path was not driven"


def test_the_scheduled_wait_path_is_the_one_being_driven(no_sleep):
    """Or the test above would pass through the ordinary retry branch."""
    assert cef._long_retry_after(_Recording(503, {"Retry-After": "7200"})) \
        == 7200.0
    assert cef._long_retry_after(_Recording(503, {"Retry-After": "5"})) is None
    assert cef.LONG_WAIT_THRESHOLD < 7200.0 <= cef.LONG_WAIT_CAP


def test_the_headers_survive_the_close(no_sleep):
    """The fix reads headers BEFORE closing; this pins that it is safe after."""
    err = _Recording(503, {"Retry-After": "7200"})
    err.close()
    assert cef._long_retry_after(err) == 7200.0, (
        "headers must outlive the socket, or the ordering becomes load-bearing")
