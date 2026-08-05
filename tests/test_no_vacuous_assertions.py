"""No test may assert something that cannot fail when the world changes.

This session found the same defect four times, in four different places: a claim
checked where it was ANNOUNCED rather than where it was ARGUED. Twice in
generated prose, once in a guard that searched a whole document for a substring
its own heading supplied, and once in a helper that was supposed to strip that
heading and silently did not. The common shape is an assertion that passes for
a reason unrelated to the thing it names.

One sub-class of that is mechanically decidable, and this file enumerates it: an
`assert` whose entire expression is built from literal module constants -- with
at most string/list methods called ON those constants -- and which therefore
compares the test file to itself. `assert EXPECTED.count("\\n\\n") == 1` was a
real instance here: it read a module-level literal and a literal, so blanking
every value in the artifact it claimed to guard would not have moved it.

WHAT THIS DOES NOT COVER, stated so the green is not over-read
---------------------------------------------------------------
Only the decidable sub-class. It cannot see a guard that reads a real file and
then asserts something trivially true of it, nor one whose substring is
satisfied by a different part of the document than the part it names -- the two
shapes that actually cost the most time this session. Those still need mutation
testing. This closes the one hole a machine can close.

THE DETECTOR VALIDATES ITSELF
------------------------------
A scan that returns zero because it is broken is indistinguishable from a scan
that returns zero because the suite is clean. So `test_the_detector_can_fire`
runs it against planted samples and requires it to flag them. Without that, this
whole file could rot into a no-op and look like a pass.
"""

import ast
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
TESTS = REPO_ROOT / "tests"

LITERAL = (ast.Constant, ast.Tuple, ast.List, ast.Set, ast.Dict, ast.JoinedStr)

# Assertions deliberately exempted, each with the reason it is not vacuous.
# Keyed by "<file>:<lineno>". Empty today; a future entry must say WHY, because
# "it was easier than fixing it" is the failure this file exists to prevent.
ALLOW: dict[str, str] = {}


def _literal_consts(tree: ast.Module) -> set:
    """Module-level names bound to a pure literal (no call anywhere inside)."""
    out = set()
    for node in tree.body:
        if isinstance(node, ast.Assign) and isinstance(node.value, LITERAL):
            if not any(isinstance(n, ast.Call) for n in ast.walk(node.value)):
                for tgt in node.targets:
                    if isinstance(tgt, ast.Name):
                        out.add(tgt.id)
    return out


def _root_name(node: ast.AST):
    while isinstance(node, (ast.Attribute, ast.Subscript)):
        node = node.value
    return node.id if isinstance(node, ast.Name) else None


def find_constant_only_assertions(source: str, filename: str) -> list:
    """(filename, lineno, test name, source) for each constant-only assertion."""
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return []
    consts = _literal_consts(tree)
    if not consts:
        return []
    hits = []
    for fn in ast.walk(tree):
        if not isinstance(fn, ast.FunctionDef) or not fn.name.startswith("test"):
            continue
        # Anything bound inside the function is evidence the assertion touches
        # something computed at run time.
        local = {a.arg for a in fn.args.args}
        for n in ast.walk(fn):
            if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Store):
                local.add(n.id)
            if isinstance(n, (ast.Import, ast.ImportFrom)):
                for al in n.names:
                    local.add((al.asname or al.name).split(".")[0])
        for node in ast.walk(fn):
            if not isinstance(node, ast.Assert):
                continue
            expr = node.test
            # A call is allowed only when it is a METHOD on one of the literal
            # constants (str.count, list.index). A call to a Name is a function
            # under test, which means the assertion exercises something.
            if any(not (isinstance(n.func, ast.Attribute)
                        and _root_name(n.func) in consts)
                   for n in ast.walk(expr) if isinstance(n, ast.Call)):
                continue
            names = {n.id for n in ast.walk(expr)
                     if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load)}
            if not names or (names & local) or not (names <= consts):
                continue
            hits.append((filename, node.lineno, fn.name, ast.unparse(expr)[:130]))
    return hits


def test_the_detector_can_fire():
    """Prove the scan works before trusting a zero from it.

    Both planted shapes must be caught, and the file-reading one must not be --
    a detector that flagged everything would be as useless as one that flagged
    nothing, just louder.
    """
    planted = (
        'EXPECTED = "a\\n\\nb"\n'
        "THRESH = 3\n"
        "\n"
        "def test_vacuous_method_on_constant():\n"
        '    assert EXPECTED.count("\\n\\n") == 1\n'
        "\n"
        "def test_vacuous_plain_constant():\n"
        "    assert THRESH > 1\n"
        "\n"
        "def test_legitimate_reads_a_file(tmp_path):\n"
        "    p = tmp_path / 'f'\n"
        "    p.write_text(EXPECTED)\n"
        '    assert p.read_text().count("\\n\\n") == 1\n'
    )
    hits = find_constant_only_assertions(planted, "planted.py")
    caught = {h[2] for h in hits}
    assert caught == {"test_vacuous_method_on_constant",
                      "test_vacuous_plain_constant"}, (
        f"the detector is broken; it caught {sorted(caught)}. A zero from it "
        "would mean nothing.")


def test_no_assertion_compares_the_test_file_to_itself():
    found = []
    for f in sorted(TESTS.glob("test_*.py")):
        for filename, lineno, fn, src in find_constant_only_assertions(
                f.read_text(), f.name):
            if f"{filename}:{lineno}" in ALLOW:
                continue
            found.append(f"  {filename}:{lineno} in {fn}\n      {src}")
    assert not found, (
        "these assertions are built only from literal constants in their own "
        "file, so they cannot fail however the artifacts they claim to guard "
        "change:\n" + "\n".join(found)
        + "\n\nEither assert against something read at run time, or add the "
          "line to ALLOW with the reason it is not vacuous.")


def test_the_allowlist_entries_still_exist():
    """A stale exemption silently widens the hole it was cut for."""
    for key in ALLOW:
        name, _, lineno = key.rpartition(":")
        path = TESTS / name
        assert path.exists(), f"ALLOW names {name}, which is gone"
        assert int(lineno) <= len(path.read_text().split("\n")), (
            f"ALLOW names {key}, past the end of that file")
