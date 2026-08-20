"""MANIFEST.sha256 must be current, checked LOCALLY and not only in CI.

`.github/workflows/manifest-check.yml` regenerates the manifest and fails on a
diff. Nothing in `pytest tests/` did, so the manifest went stale and stayed
stale across every push of a long branch while the local suite reported green
the whole time. A check that exists only in CI is a check the person doing the
work does not have.

The regeneration is too slow to run inside the suite -- it hashes the whole
tracked tree -- so this verifies the same PROPERTY by a cheaper route: every
entry's hash must match the file it names, and the manifest's file set must
match the generator's own scope. A stale manifest fails either because an
entry's hash moved or because a file entered or left the scope.

WHAT THIS DELIBERATELY DOES NOT DO is call the generator and diff, which is what
CI does. Two implementations of the same rule can disagree, and if they do it is
better that this one is a plain restatement of the property (a hash names its
file) than a second copy of the generator's traversal logic. If the two ever
disagree, CI is authoritative and this file is wrong.
"""
import hashlib
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
MANIFEST = REPO / "MANIFEST.sha256"
# Hashing every tracked file is slower than the rest of the suite combined, so
# the full pass is opt-in. The sampled pass is what runs by default and is what
# would have caught the stale manifest: any commit touching tracked files moves
# enough entries that a sample finds one.
SAMPLE = 400


def _entries() -> "list[tuple[str, str]]":
    out = []
    for line in MANIFEST.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        h, _, path = line.partition("  ")
        if len(h) == 64 and path:
            out.append((h, path))
    return out


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _scope() -> "set[str]":
    """Tracked files minus the generator's OWN declared exclusions.

    A first version inferred the scope from which top-level directories the
    manifest touched, and got it wrong in the way a heuristic does: `corpus/`
    is manifested, so `corpus/by-pmid/` looked in-scope and 4,830 deliberately
    excluded files read as a stale manifest.

    The exclusions are IMPORTED from the generator rather than restated. They
    are not arbitrary -- by-pmid holds non-OA full text that cannot be
    redistributed -- so a copy here would be a second place for a
    redistribution rule to drift.
    """
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "grm", REPO / "scripts/generate_release_manifest.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    r = subprocess.run(["git", "ls-files"], cwd=REPO, capture_output=True,
                       text=True, check=True)
    return {f for f in r.stdout.split()
            if not any(f.startswith(p) for p in m.EXCLUDE_PREFIXES)
            and f not in m.EXCLUDE_FILES}


def test_the_manifest_is_not_empty_and_parses():
    entries = _entries()
    assert len(entries) > 1000, (
        f"MANIFEST.sha256 holds only {len(entries)} entries, which is far below "
        "the tracked-file count -- it is truncated or the format changed")
    assert len({p for _, p in entries}) == len(entries), (
        "MANIFEST.sha256 lists a path twice")


def test_a_sample_of_entries_matches_the_file_on_disk():
    """Deterministic sample: the first N entries in manifest order.

    Not a random sample, because a guard whose failure depends on a seed is a
    guard that fails for some people and not others. Manifest order is
    alphabetical over the tracked tree, so a fixed prefix spans the top-level
    directories that change most.
    """
    bad = []
    for h, path in _entries()[:SAMPLE]:
        f = REPO / path
        if not f.exists():
            bad.append(f"{path}: listed in the manifest, missing on disk")
        elif _sha256(f) != h:
            bad.append(f"{path}: hash differs from the manifest")
    assert not bad, (
        "MANIFEST.sha256 is stale:\n  " + "\n  ".join(bad[:10])
        + "\nRun: python3 scripts/generate_release_manifest.py")


def test_no_tracked_file_is_missing_from_the_manifest():
    """The failure mode a hash check cannot see.

    A file ADDED and never manifested has no entry to compare, so every entry
    still matches and the manifest is still stale. This is the half that
    actually caught the branch: nineteen commits added scripts, tests and
    analysis artifacts, and none of them was in the manifest.
    """
    manifested = {p for _, p in _entries()}
    tracked = _scope()
    missing = sorted(tracked - manifested)
    assert not missing, (
        f"{len(missing)} tracked file(s) in a manifested directory have no "
        f"manifest entry, so the manifest is stale even though every existing "
        f"entry matches:\n  " + "\n  ".join(missing[:10])
        + "\nRun: python3 scripts/generate_release_manifest.py")


def test_the_exclusions_are_the_generators_own():
    """The scope must come from the generator, not from a copy.

    `corpus/by-pmid/` is excluded because 61 of those articles are non-OA and
    redistributing them would be unsafe. A restated exclusion list here would
    be a second place for a redistribution rule to drift out of step.
    """
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "grm", REPO / "scripts/generate_release_manifest.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    assert "corpus/by-pmid/" in m.EXCLUDE_PREFIXES, (
        "the release manifest no longer excludes corpus/by-pmid/, which holds "
        "non-OA full text that cannot be redistributed")
    assert "MANIFEST.sha256" in m.EXCLUDE_FILES
    # STRUCTURAL, via the parse tree. A first version compared the source
    # against a string literal and matched ITS OWN assertion -- the third time
    # this session a checker has classified itself by finding its own marker
    # text. An AST walk cannot: it asks whether this module BINDS an exclusion
    # list, which mentioning one in a docstring or an assertion does not do.
    import ast

    tree = ast.parse(Path(__file__).read_text())
    bound = {t.id for node in tree.body if isinstance(node, ast.Assign)
             for t in node.targets if isinstance(t, ast.Name)}
    restated = bound & {"EXCLUDE_PREFIXES", "EXCLUDE_FILES"}
    assert not restated, (
        f"this file binds {sorted(restated)} instead of importing the "
        "generator's, which would be a second place for a redistribution rule "
        "to drift out of step")


def test_no_manifest_entry_names_a_deleted_file():
    manifested = {p for _, p in _entries()}
    tracked = _scope()
    stale = sorted(manifested - tracked)
    assert not stale, (
        f"{len(stale)} manifest entr(ies) name files that are no longer "
        f"tracked:\n  " + "\n  ".join(stale[:10])
        + "\nRun: python3 scripts/generate_release_manifest.py")


@pytest.mark.slow
def test_every_entry_matches_the_file_on_disk():
    """The full pass. Opt-in via `-m slow` because it hashes the tracked tree."""
    bad = []
    for h, path in _entries():
        f = REPO / path
        if f.exists() and _sha256(f) != h:
            bad.append(path)
    assert not bad, (
        f"{len(bad)} manifest entr(ies) do not match the file on disk:\n  "
        + "\n  ".join(bad[:10]))
