"""Every committed analysis artifact must still be what its generator produces.

THE DEFECT THIS EXISTS FOR, recorded in this repo three times before it got a
gate. A guard that reads only a generated artifact CANNOT FAIL when its
generator changes, because the artifact is committed and does not move until
somebody re-runs it -- so the `.md` and the `.json` go stale TOGETHER and every
guard comparing one to the other stays green. `scope_audit.py` shipped a stale
row on the repo's front door with seventeen assertions passing;
`engine_time_audit.py` had the same hole; `atlas_taxonomy_reach.py` regenerated
a DIFFERENT document from the documented command and nothing noticed.

Two ad-hoc gates were written for the two scripts where it was caught. This is
the general one.

WHAT IT CHECKS, in process and without writing a single file. First, that the
committed `.md` is byte-for-byte what the generator's own `--render-only`
branch produces -- replicating that branch rather than assuming its shape,
because two generators turn out not to re-assemble at all. Second, that the
committed `.json` is what the generator would write, which is how a
formatting drift was found. Third, and structurally rather than by calling
anything: that a generator owning an `assemble()` actually USES it on the
render-only path.

THE THIRD CHECK IS THE ONE WITH TEETH. This repo documents the contract as
"`--render-only` must RE-ASSEMBLE from stored raw counts, not re-render stored
derived fields (else guards are inert)", and two generators violate it: they
load the already-assembled JSON and render it. Every guard reading a derived
field of those two is comparing the artifact to itself.

WHAT IT CANNOT CHECK, and the limit is the whole reason the census is scanned
on a schedule rather than per-PR: whether the JSON is stale with respect to
the CENSUS. Re-assembling from stored counts cannot notice that the counts
themselves are old. That needs a full re-scan, which takes minutes per
generator over 1,334 shards. So this gate covers renderer and prose drift --
the drift that actually happens when somebody edits a docstring -- and says so
rather than implying more.
"""
import ast
import importlib
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
SCRIPTS = REPO / "scripts"
sys.path.insert(0, str(SCRIPTS))

# Generators that legitimately cannot round-trip in process, each with the
# reason. This list is checked for staleness by a test below: an entry that no
# longer names a real generator fails, so it cannot quietly become a dumping
# ground.
EXEMPT = {
    "atlas_thesis_position": "has no render() -- it prints and writes inline",
    "atlas_ambiguity": "writes no OUT_MD/OUT_JSON pair at module level",
}


def _generators():
    """Every script exposing the --render-only contract, found by parsing."""
    out = []
    for f in sorted(list(SCRIPTS.glob("census_*.py"))
                    + list(SCRIPTS.glob("atlas_*.py"))):
        src = f.read_text(encoding="utf-8")
        if "render-only" not in src:
            continue
        tree = ast.parse(src)
        fns = {n.name for n in tree.body if isinstance(n, ast.FunctionDef)}
        consts = {t.id for n in tree.body if isinstance(n, ast.Assign)
                  for t in n.targets if isinstance(t, ast.Name)}
        out.append((f.stem, "render" in fns, "assemble" in fns,
                    "OUT_MD" in consts and "OUT_JSON" in consts))
    return out


GENERATORS = _generators()
LIVE = [g[0] for g in GENERATORS if g[1] and g[3] and g[0] not in EXEMPT]


def test_the_generator_list_is_discovered_not_listed():
    """A hand-maintained list goes stale, which is this file's whole subject."""
    assert len(GENERATORS) >= 25, (
        f"only {len(GENERATORS)} generators discovered; the parse rule has "
        "probably stopped matching")
    assert len(LIVE) >= 23


def test_every_exemption_still_names_a_real_generator():
    """An exemption for a script that no longer exists, or that has since
    grown the interface, is an exemption nobody re-examined."""
    names = {g[0] for g in GENERATORS}
    for name, reason in EXEMPT.items():
        assert name in names, f"exemption {name!r} names no generator"
        assert reason.strip(), f"exemption {name!r} has no reason"
        g = next(g for g in GENERATORS if g[0] == name)
        assert not (g[1] and g[3]), (
            f"{name} now has render() and OUT_MD/OUT_JSON, so its exemption "
            f"({reason}) is stale and it should be gated")


def _render_only_reassembles(name: str) -> bool:
    """Does the generator's --render-only branch pass through assemble()?

    Parsed rather than called. The question is about the code path the
    documented command takes, and a call cannot distinguish "re-assembled" from
    "the stored fields happened to match".
    """
    tree = ast.parse((SCRIPTS / f"{name}.py").read_text(encoding="utf-8"))
    main = next((n for n in tree.body
                 if isinstance(n, ast.FunctionDef) and n.name == "main"), None)
    if main is None:
        return False
    for node in ast.walk(main):
        # The render-only branch is `if a.render_only:` / `if args.render_only:`
        if not (isinstance(node, ast.If)
                and isinstance(node.test, ast.Attribute)
                and node.test.attr == "render_only"):
            continue
        return any(isinstance(c, ast.Call) and getattr(c.func, "id", None) == "assemble"
                   for stmt in node.body for c in ast.walk(stmt))
    # Some generators inline it as `assemble(load(...) if render_only else scan())`
    for node in ast.walk(main):
        if isinstance(node, ast.Call) and getattr(node.func, "id", None) == "assemble":
            if any(isinstance(a, ast.IfExp) for a in ast.walk(node)):
                return True
    return False


def _reassemble(mod, stored):
    """Call `assemble` the way the generator's render-only branch does.

    An `assemble` taking several scan products needs those products back, so
    the generator supplies a `_split_stored` that recovers them from the merged
    artifact. Guessing the arity is what made the first version of this gate
    report two false failures.
    """
    if hasattr(mod, "_split_stored"):
        return mod.assemble(*mod._split_stored(stored))
    return mod.assemble(stored)


def _reproduce(mod, name):
    """Replicate the generator's own render-only path."""
    stored = json.loads(mod.OUT_JSON.read_text())
    if hasattr(mod, "assemble") and _render_only_reassembles(name):
        return _reassemble(mod, stored)
    return stored


@pytest.mark.parametrize("name", LIVE)
def test_the_committed_markdown_is_what_the_renderer_produces(name):
    mod = importlib.import_module(name)
    produced = mod.render(_reproduce(mod, name))
    committed = mod.OUT_MD.read_text()
    assert produced == committed, (
        f"analysis/{mod.OUT_MD.name} is not what {name} renders. Either the "
        "generator changed and the artifact was not regenerated, or the "
        f"artifact was hand-edited. Run `python scripts/{name}.py --render-only`.")


def _dump_kwargs(name: str) -> dict:
    """The json.dumps options the generator actually writes with.

    Assuming them is how the first version of this gate reported two false
    failures: two generators pass `sort_keys=True` and the check did not, so
    every key order looked wrong. A gate that guesses the format is testing its
    own guess.
    """
    tree = ast.parse((SCRIPTS / f"{name}.py").read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr == "write_text"
                and isinstance(node.func.value, ast.Name)
                and node.func.value.id == "OUT_JSON"):
            continue
        for arg in ast.walk(node):
            if (isinstance(arg, ast.Call)
                    and getattr(arg.func, "attr", None) == "dumps"):
                return {k.arg: ast.literal_eval(k.value) for k in arg.keywords}
    return {}


@pytest.mark.parametrize("name", LIVE)
def test_the_committed_json_is_what_the_generator_writes(name):
    """A committed artifact that differs from its own generator's output makes
    every regenerate-and-diff check permanently dirty, which trains people to
    ignore a dirty tree -- the state this gate exists to make meaningful."""
    mod = importlib.import_module(name)
    committed = mod.OUT_JSON.read_text()
    produced = json.dumps(_reproduce(mod, name), **_dump_kwargs(name)) + "\n"
    assert produced == committed, (
        f"analysis/{mod.OUT_JSON.name} does not round-trip through {name}. "
        f"Regenerate it with `python scripts/{name}.py --render-only`.")


ASSEMBLERS = [n for n in LIVE if any(g[0] == n and g[2] for g in GENERATORS)]


@pytest.mark.parametrize("name", ASSEMBLERS)
def test_render_only_reassembles_from_the_stored_raw_counts(name):
    """The contract that keeps every --render-only guard from being inert.

    A generator that owns an `assemble()` and then renders the stored,
    already-assembled JSON is re-rendering its own derived fields. A guard
    checking one of those fields is then comparing the artifact to itself and
    cannot fail, whatever the generator does.
    """
    assert _render_only_reassembles(name), (
        f"{name} has an assemble() and its --render-only branch does not call "
        "it, so the documented command re-renders stored derived fields "
        "instead of recomputing them. Every guard reading a derived field of "
        "this artifact is inert.")


def test_assemble_is_idempotent_where_it_is_used(name=None):
    """Re-assembling an assembled dict must not move anything.

    If it does, the raw counts are not sufficient to rebuild the derived
    fields, and the render-only contract cannot be satisfied.
    """
    checked = 0
    for n in ASSEMBLERS:
        if not _render_only_reassembles(n):
            continue
        mod = importlib.import_module(n)
        stored = json.loads(mod.OUT_JSON.read_text())
        once = _reassemble(mod, stored)
        twice = _reassemble(mod, json.loads(json.dumps(once)))
        assert once == twice, f"{n}.assemble() is not idempotent"
        checked += 1
    assert checked >= 10, f"only {checked} generators exercised"
