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

WHAT IT CANNOT CHECK, stated at full width because an adversarial review
found the first version of this paragraph understated it.

The obvious limit is census staleness: re-assembling from stored counts cannot
notice the counts are old, and a full re-scan is minutes per generator over
1,334 shards. But the real limit is broader. **This gate cannot detect any
wrong value that re-assembly reproduces**, and that includes every RAW COUNT
in every artifact -- which is most of their numeric content. Corrupt a stored
scan total and re-assembly propagates the corruption consistently into the
derived fields and the prose, and both checks pass.

IT ALSO CANNOT TELL A CORRECT RANK FROM A WRONG-BUT-EXPLICIT ONE. The
order-dependence gate below catches an order INHERITED from serialisation; a
renderer that sorts explicitly in the wrong direction is deterministic,
reproducible, and green. Re-introducing this campaign's own regression -- an
ascending sort where the sweep steps outward from a bound -- leaves the whole
suite passing with the bound on the last row under prose describing the first.
That is precisely the failure this campaign suffered, so the individual
orderings need their own regression guards; `tests/test_report_orderings.py`
pins them.

What it does catch, precisely: an artifact that no longer matches the code
that produces it. A renderer edited without regenerating, a `.md` hand-edited,
a `.json` whose format drifted, a derived field that stopped being recomputed.
That is renderer and prose drift -- the drift that actually happens when
somebody edits a docstring -- and it is genuinely less than "the artifacts are
correct".

A quantitative version of this limit was attempted and withdrawn: measuring
the fraction of perturbed values that re-assembly restores conflates a raw
count, which SHOULD move the output, with a derived field that is merely
copied through. Two attempts got the polarity wrong in opposite directions, so
the honest deliverable is the scoped claim above rather than a third guess at
a metric.
"""
import ast
import importlib
import inspect
import json
import sys
import textwrap
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
SCRIPTS = REPO / "scripts"
sys.path.insert(0, str(SCRIPTS))

# Generators that legitimately cannot round-trip in process, each with the
# reason. This list is checked for staleness by a test below: an entry that no
# longer names a real generator fails, so it cannot quietly become a dumping
# ground.
# The names a pure renderer goes by in this repo. Keying on ONE name was a
# dodge hatch: an earlier version excluded 13 scripts for "having no
# renderer", and 12 of them had a top-level `write_report(r)` taking the
# stored dict and gated with zero refactoring. A quarter of the corpus was
# excluded on a false premise, in a gate whose whole thesis is that the
# previous exclusions were unjustified.
RENDERER_NAMES = ("render", "write_report")

# Owns an artifact pair but cannot be exercised in process. Two distinct
# reasons, and the FIRST version of this list gave a false one -- it said these
# scripts "build their report inline while writing" and would need refactoring,
# when 12 of the 13 expose a top-level `write_report(r)`. A reviewer was right
# that the justification was untrue; the real obstacle is different and is
# CHECKED below rather than asserted:
#
#   `write_report(d)` WRITES the artifact and returns None. Calling it to
#   compare would write to the repository, which is the property this gate
#   promises not to have -- and writing during a check is exactly how a
#   strided sample scan once clobbered a committed sidecar. Making these
#   gateable means splitting each into `render(d) -> str` plus a writer, which
#   is a real improvement and is not folded in here.
#
# `rare_event_analysis` has no renderer under any recognised name at all.
#
# The check below tries EVERY name in RENDERER_NAMES and requires a
# write-only renderer to actually return None, so a script cannot be parked
# here by renaming its renderer or by claiming a reason nobody verified.
NO_RENDERER = {
    "abc_joint_posterior", "abc_posterior", "calibrate_erastin",
    "calibrate_kill_switch", "calibrate_pk", "embed_evidence_leg",
    "identifiability_report", "rare_event_analysis", "validate_penetration",
    "validate_rd_vs_biofvm", "validate_spheroid_kill",
    "validate_spheroid_structure", "validate_trigger_wave",
}
EXEMPT: dict = {}


def _generators():
    """Every script OWNING a committed artifact pair, found by parsing.

    Keyed on the interface, not on a filename. The first version globbed
    `census_*` and `atlas_*` and filtered on the substring "render-only", which
    missed FOURTEEN generators with the identical interface -- including
    `scope_audit` and `engine_time_audit`, the two this file's docstring cites
    as the motivating defects, and five that carry `--render-only` and were
    excluded purely by filename. A gate advertised as the general one was not
    general over either of the cases that motivated it.
    """
    out = []
    for f in sorted(SCRIPTS.glob("*.py")):
        try:
            tree = ast.parse(f.read_text(encoding="utf-8"))
        except SyntaxError:
            continue
        fns = {n.name for n in tree.body if isinstance(n, ast.FunctionDef)}
        consts = {t.id for n in tree.body if isinstance(n, ast.Assign)
                  for t in n.targets if isinstance(t, ast.Name)}
        if not ("OUT_MD" in consts and "OUT_JSON" in consts):
            continue
        out.append((f.stem, any(n in fns for n in RENDERER_NAMES),
                    "assemble" in fns, True))
    return out


GENERATORS = _generators()
LIVE = [g[0] for g in GENERATORS
        if g[1] and g[0] not in EXEMPT and g[0] not in NO_RENDERER]
# Pinned EXACTLY, not as a floor. A floor with slack lets a generator drop out
# of the gate silently: at `>= 25` against 26, deleting the marker from one
# script left the suite green with two parametrised cases quietly gone.
EXPECTED_GENERATORS = 61


def test_the_generator_list_is_discovered_not_listed():
    """A hand-maintained list goes stale, which is this file's whole subject --
    and so does a floor with slack."""
    assert len(GENERATORS) == EXPECTED_GENERATORS, (
        f"{len(GENERATORS)} generators discovered, expected exactly "
        f"{EXPECTED_GENERATORS}. If a generator was added or removed, update "
        "EXPECTED_GENERATORS in the same commit; if the parse rule stopped "
        "matching, fix the rule.")
    assert len(LIVE) == EXPECTED_GENERATORS - len(NO_RENDERER) - len(EXEMPT), (
        f"{len(LIVE)} gated of {EXPECTED_GENERATORS} discovered, with "
        f"{len(NO_RENDERER)} lacking a renderer and {len(EXEMPT)} exempt")


def test_every_known_generator_is_actually_gated():
    """Coverage by NAME, so a generator cannot leave the gate unnoticed.

    The count check above catches a generator disappearing; this catches one
    that is still discovered but silently stops being exercised.
    """
    names = {g[0] for g in GENERATORS}
    for name in sorted(names - set(EXEMPT) - NO_RENDERER):
        assert name in LIVE, (
            f"{name} owns a committed artifact pair and is neither gated, "
            "exempt, nor listed as having no renderer")


def test_every_ungated_generator_is_ungated_for_a_checked_reason():
    """An exemption nobody re-examines is how a gate hollows out.

    Both lists are verified against the source rather than trusted: a script in
    NO_RENDERER that has since grown a `render()` fails, and an EXEMPT entry
    must name a real generator and carry a reason.
    """
    names = {g[0] for g in GENERATORS}
    renders = {g[0] for g in GENERATORS if g[1]}
    for name in sorted(NO_RENDERER):
        assert name in names, f"NO_RENDERER lists {name!r}, which is no generator"
        if name not in renders:
            continue                      # no renderer at all: nothing to check
        # It HAS a renderer, so the exclusion rests on that renderer writing
        # rather than returning. Verify it, on a dict it cannot act on, so the
        # claim is measured rather than believed.
        mod = importlib.import_module(name)
        f = _renderer(mod)
        src = inspect.getsource(f)
        assert "write_text" in src, (
            f"{name} exposes {f.__name__}() which does not write, so it can be "
            "exercised in process and should be gated rather than excluded")
        # Only the function's OWN returns, not those of helpers nested in it.
        fn = ast.parse(textwrap.dedent(src)).body[0]
        own = []
        for node in ast.walk(fn):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node is not fn:
                continue
            if isinstance(node, ast.Return):
                own.append(node)
        nested = {id(n) for f2 in ast.walk(fn)
                  if isinstance(f2, (ast.FunctionDef, ast.AsyncFunctionDef)) and f2 is not fn
                  for n in ast.walk(f2) if isinstance(n, ast.Return)}
        returns_str = any(n.value is not None for n in own if id(n) not in nested)
        assert not returns_str, (
            f"{name}.{f.__name__}() now returns a value, so it can be compared "
            "without writing and should be gated")
    for name, reason in EXEMPT.items():
        assert name in names, f"exemption {name!r} names no generator"
        assert reason.strip(), f"exemption {name!r} has no reason"


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
    """Replicate the generator's own render-only path.

    ROUND-TRIPPED, always. Generators write with `sort_keys=True` and several
    then rendered the in-memory dict, producing a document that could not be
    reproduced from its own artifact -- four were found that way. Rendering the
    round-tripped value is what the artifact will actually contain, and it also
    normalises tuples, which JSON turns into lists and a renderer may one day
    distinguish.
    """
    stored = json.loads(mod.OUT_JSON.read_text())
    if hasattr(mod, "assemble") and _render_only_reassembles(name):
        # Round-trip with the generator's OWN dump options, read from its
        # source. Using a guessed set was the exact failure `_dump_kwargs`'s
        # docstring warns about: 16 of 18 assemblers change key order under
        # `sort_keys`, and for two of them the rendered markdown differs, so a
        # gate that omitted it would disagree with the generator the moment
        # anyone adopted the documented round-trip pattern. Measured: 18 of 18
        # assemblers change key order under `sort_keys`, and for one --
        # `burden_data_feasibility` -- the rendered markdown differs. (An
        # earlier hand-written "16 of 18 ... two of them" was wrong on both
        # counts, and this campaign's own ordering fixes are what moved the
        # second figure to one.)
        kw = {k: v for k, v in _dump_kwargs(name).items() if k == "sort_keys"}
        return json.loads(json.dumps(_reassemble(mod, stored), **kw))
    return stored


def _renderer(mod):
    """The generator's renderer, whatever it is called here."""
    for n in RENDERER_NAMES:
        f = getattr(mod, n, None)
        if callable(f):
            return f
    return None


def _render(mod, d):
    """Call the renderer the way the generator does.

    One renderer needs an identifier->name resolver built from a committed
    label table. The generator exposes a factory for it rather than closing
    over it inside `main`, which is what let that one be gated at all.
    """
    f = _renderer(mod)
    if hasattr(mod, "make_namer"):
        return f(d, mod.make_namer())
    return f(d)


def _remediation(name: str) -> str:
    """The command that actually regenerates THIS artifact.

    Nine gated generators have no `--render-only`, so a message prescribing it
    sends the reader to a flag that errors out. Read from the script's own
    argparse rather than assumed.
    """
    src = (SCRIPTS / f"{name}.py").read_text(encoding="utf-8")
    if "--render-only" in src:
        return f"python scripts/{name}.py --render-only"
    return (f"python scripts/{name}.py  (no --render-only: this one "
            "re-runs its full scan or simulation)")


@pytest.mark.parametrize("name", LIVE)
def test_the_committed_markdown_is_what_the_renderer_produces(name):
    mod = importlib.import_module(name)
    produced = _render(mod, _reproduce(mod, name))
    committed = mod.OUT_MD.read_text()
    assert produced == committed, (
        f"analysis/{mod.OUT_MD.name} is not what {name} renders. Either the "
        "generator changed and the artifact was not regenerated, or the "
        f"artifact was hand-edited. Run `{_remediation(name)}`.")


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
        f"Regenerate it with `{_remediation(name)}`.")


ASSEMBLERS = [n for n in LIVE if any(g[0] == n and g[2] for g in GENERATORS)]


@pytest.mark.parametrize("name", ASSEMBLERS)
def test_render_only_reassembles_from_the_stored_raw_counts(name):
    """The contract that keeps every --render-only guard from being inert.

    STRUCTURAL, and that is all it claims: it proves `assemble` is CALLED on
    the render-only path, not that the call re-derives any particular amount.
    An `assemble` that copied its input would satisfy this and the idempotency
    check both. What rules that out here is reading the generators, not this
    test.

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


def test_assemble_is_idempotent_where_it_is_used():
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
    # EXACT, not a floor -- this file argues at the top that a floor with
    # slack lets a generator drop out silently, and then shipped `>= 10`
    # against a true 18.
    assert checked == len(ASSEMBLERS), (
        f"{checked} of {len(ASSEMBLERS)} assemblers exercised")


def _shuffled(obj, rng):
    """The same data with every dict's key order randomised."""
    if isinstance(obj, dict):
        items = list(obj.items())
        rng.shuffle(items)
        return {k: _shuffled(v, rng) for k, v in items}
    if isinstance(obj, list):
        return [_shuffled(v, rng) for v in obj]
    return obj


# Generators whose rendered order is still inherited from dict iteration
# rather than established by the renderer. A RATCHET, not an allowlist: the
# test below fails if a listed generator has been fixed and not removed, so the
# set can only shrink, and it fails if any UNLISTED generator joins them, so it
# cannot grow silently.
#
# Left rather than fixed in one sitting deliberately. Eight ordering defects
# were repaired in this campaign and one of the repairs flipped a published
# verdict, because each fix needs the ranked quantity identified individually
# and checked against what the artifact published. Batching the last five at
# the end of a long change is how that regression happened.
ORDER_DEPENDENT = {
    "atlas_ingest_sensitivity",
    "atlas_recent_window",
    "atlas_site_coverage",
    "census_mechanism_growth",
}


def _order_dependent(name: str) -> bool:
    import random

    mod = importlib.import_module(name)
    stored = json.loads(mod.OUT_JSON.read_text())
    base = _render(mod, _reproduce(mod, name))
    reassembles = hasattr(mod, "assemble") and _render_only_reassembles(name)
    rng = random.Random(20260823)
    for _ in range(5):
        s = _shuffled(stored, rng)
        d = _reassemble(mod, s) if reassembles else s
        if _render(mod, d) != base:
            return True
    return False


@pytest.mark.parametrize("name", LIVE)
def test_the_report_does_not_depend_on_dict_order(name):
    """A rendered order must be established by the renderer, not inherited.

    THIS IS THE GATE FOR THE WORST THING THAT HAPPENED IN THIS CAMPAIGN. A
    report that renders a dict in iteration order is reproducible only until
    something serialises it -- and `json.dumps(sort_keys=True)` does. Four
    tables silently lost a `Counter.most_common()` ranking that way, putting
    `Expression of Concern` (172) above `Retracted Publication` (7,900); and a
    monotone test reading dict order flipped a published verdict from "3 of 3"
    to "2 of 3", failing the one admissible case.

    Shuffling every dict and requiring the output to be unchanged catches the
    whole class, including the orders nobody has thought about yet. When it
    fires the fix is an explicit sort INSIDE the renderer -- never a sort of
    the input, which replaces a rank with an alphabet.
    """
    assert ORDER_DEPENDENT <= set(LIVE), (
        f"ORDER_DEPENDENT holds names that are not gated generators: "
        f"{sorted(ORDER_DEPENDENT - set(LIVE))}. A renamed or removed "
        "generator leaves a dead ratchet entry, which is the unexamined "
        "exclusion list this file exists to argue against.")
    dependent = _order_dependent(name)
    if name in ORDER_DEPENDENT:
        assert dependent, (
            f"{name} is listed in ORDER_DEPENDENT and no longer depends on "
            "dict order. Remove it from the set -- the list is a ratchet and "
            "may only shrink.")
        pytest.skip(f"{name} is a known-order-dependent generator (ratcheted)")
    assert not dependent, (
        f"{name} renders differently when its stored dicts are reordered, so "
        "its published order is inherited from serialisation rather than "
        "established by the renderer. Sort inside render(), on the quantity "
        "the table is ranked by -- never sort the input, which replaces a rank "
        "with an alphabet.")
