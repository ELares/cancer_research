"""No tool may invent the path to externally-configured storage (#960).

A MOUNT POINT IS JUST A DIRECTORY. With a share unmounted its mountpoint is an
ordinary empty local directory, so `root.mkdir(parents=True, exist_ok=True)`
creates the whole path on the boot disk and the caller writes its corpus into
the wrong filesystem -- silently, with the progress ledgers (which live in the
repository, not on the share) still reporting success, and the bytes then
invisible to every later run because the run looks for them on the share.

THIS GUARD IS AT CONCEPT LEVEL, BECAUSE FIXING THE SITE THAT BIT WAS NOT
FIXING THE CLASS. The crawl's shard root was the one that failed, and a sweep
for roots defaulting under the same mountpoint found the open-access full-text
root doing exactly the same thing. Both now go through
`storage_root.require_output_root`.
"""
from __future__ import annotations

import ast
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
SCRIPTS = REPO / "scripts"

# NAMED, not counted. A count cannot tell a root legitimately moved off the
# share from a detector that stopped recognising one -- the vacuity lesson
# from #954, where `assert not offenders` passed on an empty scan.
KNOWN_MOUNT_ROOTED = {
    "atlas_fulltext.py": {"DEFAULT_FT"},
    "corpus_expand_fetch.py": {"OUT_ROOT"},
    "corpus_expand_report.py": {"OUT_ROOT"},
    "corpus_identity_index.py": {"NAS_FULLTEXT"},
    "corpus_retrieval_index.py": {"ROOT"},
}

# The mountpoint component every one of these defaults through.
MOUNT_COMPONENT = '"nas"'


def _mount_rooted_names(tree: ast.AST, src: str) -> set[str]:
    """Module-level names whose DEFAULT value routes through the mountpoint."""
    names = set()
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        seg = ast.get_source_segment(src, node.value) or ""
        if MOUNT_COMPONENT in seg:
            names.update(t.id for t in node.targets if isinstance(t, ast.Name))
    return names


def _parents_true_mkdirs(tree: ast.AST, names: set[str]) -> list[tuple[int, str]]:
    """`X.mkdir(parents=True)` where X is one of those roots."""
    out = []
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr == "mkdir"):
            continue
        if not any(kw.arg == "parents" and getattr(kw.value, "value", False)
                   for kw in node.keywords):
            continue
        base = node.func.value
        # `ROOT.mkdir(...)` and `ROOT.parent.mkdir(...)` alike.
        while isinstance(base, ast.Attribute):
            base = base.value
        if isinstance(base, ast.Name) and base.id in names:
            out.append((node.lineno, base.id))
    return out


def _scan() -> tuple[dict[str, set[str]], list[str]]:
    found: dict[str, set[str]] = {}
    offenders: list[str] = []
    for path in sorted(SCRIPTS.glob("*.py")):
        src = path.read_text()
        tree = ast.parse(src)
        names = _mount_rooted_names(tree, src)
        if names:
            found[path.name] = names
        for lineno, name in _parents_true_mkdirs(tree, names):
            offenders.append(f"{path.name}:{lineno} ({name}.mkdir(parents=True))")
    return found, offenders


def test_every_mount_rooted_path_is_still_detected():
    """The floor. A scan that finds nothing is indistinguishable from a clean
    repo, so the roots are named and each must still be recognised."""
    found, _ = _scan()
    for fname, names in KNOWN_MOUNT_ROOTED.items():
        assert fname in found, (
            f"{fname} no longer presents any mount-rooted path; if its root "
            "genuinely moved off the share, remove it from KNOWN_MOUNT_ROOTED "
            "deliberately rather than letting the detector go quiet")
        missing = names - found[fname]
        assert not missing, f"{fname}: no longer detected: {sorted(missing)}"


def test_no_tool_creates_mount_rooted_storage_with_parents():
    _, offenders = _scan()
    assert not offenders, (
        "these create a path that defaults under a network mountpoint, so with "
        "the share unmounted they write to the boot disk instead of refusing:\n  "
        + "\n  ".join(offenders))


@pytest.mark.parametrize("writer", ["corpus_expand_fetch.py", "atlas_fulltext.py",
                                   "openalex_fetch.py", "trials_fetch.py"])
def test_the_two_writers_consult_the_shared_refusal(writer):
    """Absence of `parents=True` is not presence of a check: a tool that
    simply dropped the mkdir would satisfy the guard above and then fail
    later with a bare FileNotFoundError instead of the actionable refusal.

    All FOUR writers, which is the point of the sweep. Reading the code found
    the crawl and then the full-text root; the AST scan found two more
    (`openalex_fetch`, `trials_fetch`) whose roots are declared in a shape my
    grep did not match. Fixing the site that bit was not fixing the class, and
    enumerating the class by hand did not finish the job either.
    """
    src = (SCRIPTS / writer).read_text()
    tree = ast.parse(src)
    imported = {a.name for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)
                and n.module == "storage_root" for a in n.names}
    assert "require_output_root" in imported, (
        f"{writer} must import the shared refusal, not reimplement it -- two "
        "copies of a safety check is how they drift apart")
    called = {n.func.id for n in ast.walk(tree)
              if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)}
    assert "require_output_root" in called, f"{writer} imports it but never calls it"


def test_the_detector_actually_detects(tmp_path):
    """Self-validation. A scanner returning zero because it is BROKEN looks
    exactly like a clean repo (`test_no_vacuous_assertions.py` precedent)."""
    src = (
        'from pathlib import Path\n'
        'ROOT = Path.home() / "nas" / "share" / "leaf"\n'
        'LOCAL = Path("corpus") / "atlas"\n'
        'def go():\n'
        '    ROOT.mkdir(parents=True, exist_ok=True)\n'
        '    LOCAL.mkdir(parents=True, exist_ok=True)\n'
    )
    tree = ast.parse(src)
    names = _mount_rooted_names(tree, src)
    assert names == {"ROOT"}, f"the mount-rooted name must be found alone; got {names}"
    hits = _parents_true_mkdirs(tree, names)
    assert [n for _, n in hits] == ["ROOT"], (
        "the local path must NOT be flagged -- a gate that fires on repo-local "
        "mkdirs would train exemption-adding rather than fixing")
    # And the fixed shape must come back clean.
    fixed = src.replace("    ROOT.mkdir(parents=True, exist_ok=True)\n",
                        "    require_output_root(ROOT)\n")
    t2 = ast.parse(fixed)
    assert not _parents_true_mkdirs(t2, _mount_rooted_names(t2, fixed))
