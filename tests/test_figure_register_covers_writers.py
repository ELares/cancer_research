"""Every figure a generator writes must be in the register, and vice versa.

THE DEFECT THIS FILE EXISTS FOR: `scripts/generate_figures.py` drew seven
figures that FIGURES.yaml did not list, that the manuscript did not include,
and that were committed nowhere. A clean generator run therefore left fourteen
untracked files in the working tree, and nothing said whether that was a
backlog or a bug. The register was checked against the figures DIRECTORY, so a
stem that was never committed was invisible to it -- the denominator was the
committed files, not the writes.

This file keys the denominator on the WRITE instead: it parses each generator
for the stems it saves and requires the register to account for every one, in
both directions. A new figure that nobody registers fails here, and a register
entry naming a generator that does not write it fails here too.

The scan must find BOTH spellings matplotlib is called with in this repo --
`fig.savefig(FIG_DIR / "fig1_x.pdf")` and the f-string loop form
`fig.savefig(FIG_DIR / f"fig28_x.{ext}")` -- because a scan that sees only the
literal form silently drops fig28 and reports a false gap.
"""
import ast
import re
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parent.parent
FIGURES = REPO / "FIGURES.yaml"

# DERIVED FROM THE REGISTER, which is what the comment here used to claim
# while the list was a four-element literal nothing checked. A new generator
# script -- even one the register NAMES -- was then scanned by nothing: a
# figure it drew and did not register passed, which is the defect this file
# exists to close, reopened one script over.
def _python_generators():
    d = yaml.safe_load(FIGURES.read_text())
    entries = d["figures"] if isinstance(d, dict) and "figures" in d else d
    if isinstance(entries, dict):
        entries = list(entries.values())
    out = set()
    for e in entries or []:
        script = ((e or {}).get("generator") or {}).get("script", "")
        if script.endswith(".py"):
            out.add(script)
    return tuple(sorted(out))

_STEM = re.compile(r"(fig[0-9]+[a-z]?_[A-Za-z0-9_]+)")


def _written_stems(src: str) -> set[str]:
    """Stems saved by this source, from savefig call sites only.

    Walking the AST rather than grepping the whole file keeps a stem mentioned
    in a docstring or a print() out of the denominator -- both occur here, and
    either would make this guard report figures nothing writes.
    """
    stems: set[str] = set()
    tree = ast.parse(src)
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        fn = node.func
        if not (isinstance(fn, ast.Attribute) and fn.attr == "savefig"):
            continue
        # ARGS AND KEYWORDS. `fname` is savefig's first parameter, so
        # `savefig(fname=...)` is ordinary Python and bypassed a scan that
        # read `node.args` only -- an unregistered figure written that way
        # passed.
        for arg in list(node.args) + [k.value for k in node.keywords]:
            for piece in ast.walk(arg):
                if isinstance(piece, ast.Constant) and isinstance(piece.value, str):
                    m = _STEM.search(piece.value)
                    if m:
                        stems.add(m.group(1))
    return stems


@pytest.fixture(scope="module")
def register():
    d = yaml.safe_load(FIGURES.read_text())
    figs = d["figures"] if isinstance(d, dict) and "figures" in d else d
    return {e["filename"]: e for e in figs}


@pytest.fixture(scope="module")
def writes():
    """stem -> the scripts whose savefig calls name it as a literal.

    ONLY TWO OF THE GENERATORS CONTRIBUTE STEMS, and the two that do not are
    blind for DIFFERENT reasons, which an earlier version of this note got
    wrong by describing one of them twice:

      * `generate_conceptual_diagrams.py` has no `savefig` at all -- fig19 and
        fig22 are produced by `subprocess.run(["dot", ...])`, and its other
        figures save through an `f"{name}.{ext}"` loop helper;
      * `rare_event_analysis.py` passes a LITERAL stem into a helper
        (`figure(by_cond, FIG_DIR / "fig29_rare_event_resolution")`) and the
        `savefig` lives inside that helper, so the literal is recoverable from
        the source -- the scan simply does not follow calls.

    So "every write is registered" is vacuously true for those two: an
    unregistered figure added to either passes. That is the limit, stated
    plainly rather than described as a future risk, and it is why the reverse
    direction below is a definition check and not a write check.
    """
    out: dict[str, set[str]] = {}
    gens = _python_generators()
    assert gens, "FIGURES.yaml names no Python generator; the scan has no input"
    for rel in gens:
        p = REPO / rel
        assert p.exists(), f"{rel} is named as a generator in FIGURES.yaml but absent"
        for stem in _written_stems(p.read_text()):
            out.setdefault(stem, set()).add(rel)
    return out


def test_the_scan_still_finds_writes(writes):
    """A denominator that silently empties makes every guard below vacuous.

    `_written_stems` keys on `.savefig(...)`. If a refactor routes saves through
    a helper, the scan returns nothing and 'every write is registered' becomes
    trivially true. Pin the two scripts it can read, and pin fig28 by name --
    it is the only f-string save, so it is the case a literal-only scan drops.
    """
    # THE BY-NAME PINS BELOW ARE THE REAL CHECK. A count with a margin of one
    # (21 actual against 20) fails for the wrong reason when a figure is
    # legitimately retired, and passes at exactly 20 when one figure's saves
    # are routed through a helper. It is kept only as a floor that a
    # WHOLESALE change of save mechanism would trip.
    assert len(writes) >= 10, (
        f"the savefig scan found only {len(writes)} stems, against 21 today. "
        "Saves have stopped going through .savefig(), so every assertion in "
        "this file has gone vacuous rather than failing")
    assert "fig28_census_capture" in writes, (
        "fig28 is saved as f\"fig28_census_capture.{ext}\" -- losing it means "
        "the scan has stopped reading f-string saves, and any figure written "
        "that way is now outside the denominator")
    assert "fig1_ferroptosis_comparison" in writes, (
        "the plain-literal save form is no longer being read")


def test_every_figure_a_generator_writes_is_in_the_register(writes, register):
    """The direction that caught the real defect.

    Seven stems were drawn on every run, registered nowhere, included by
    nothing and committed nowhere, so a clean run left fourteen untracked files
    behind. Checking the register against the figures DIRECTORY could not see
    them: a figure that was never committed is not in that denominator.
    """
    unregistered = {s: sorted(v) for s, v in writes.items() if s not in register}
    assert not unregistered, (
        "these figures are drawn but have no FIGURES.yaml entry, so nothing "
        f"records whether they are wanted: {unregistered}. Either register the "
        "figure or stop drawing it -- a stem that is drawn and unregistered is "
        "committed nowhere and dirties the tree on every generator run")


def test_every_register_entry_names_a_generator_that_defines_it(register):
    """The reverse direction, as far as it can honestly be checked.

    Two generators build their output stem in a variable, so 'this script
    writes that stem' is not recoverable from the source. What IS recoverable
    is whether the declared function exists in the declared script, which is
    what a repointed or deleted generator breaks.
    """
    broken = []
    for name, e in sorted(register.items()):
        g = e.get("generator") or {}
        script, fn = g.get("script"), g.get("function")
        assert script, f"{name} has no generator script recorded"
        src_path = REPO / script
        if not src_path.exists():
            broken.append(f"{name}: generator script {script} does not exist")
            continue
        if fn and f"def {fn}(" not in src_path.read_text():
            broken.append(f"{name}: {script} defines no `{fn}`")
    assert not broken, (
        "FIGURES.yaml points at generators that cannot produce these "
        f"figures: {broken}")


# The corpus figures this repo stopped drawing, each mapped to the census entry
# whose `note` records why. The mapping is the evidence that retiring them was a
# decision already taken and written down, not one made by deleting the code.
RETIRED = {
    "fig2_mechanism_heatmap": "fig2c_census_volume",
    "fig3_literature_disconnect": "fig5c_mechanism_site_matrix",
    "fig5_publication_trends": "fig1c_ratio_straddle",
    "fig9_evidence_tiers": "fig9c_design_composition",
    "fig14_tissue_mechanism_heatmap": "fig14c_class_by_site",
    "fig15_designed_combinations": "fig15c_mechanism_pairs",
    # RETIRED THE SAME WAY AND OMITTED. `fig23_census_flow`'s note opens
    # `Replaces fig23_prisma_flow.`, identical in form to the seven, and it is
    # drawn nowhere and registered nowhere. The pin covered 7 of 8.
    "fig23_prisma_flow": "fig23_census_flow",
    "fig16_weighted_evidence": "fig16c_trial_share",
}


def test_each_retired_figure_is_named_by_the_entry_that_replaced_it(register):
    """Deleting a generator function deletes the only other record of it.

    Once the code is gone, `FIGURES.yaml` is the sole place saying these seven
    ever existed and what took over. If a later edit tidies those sentences
    away, the repository silently loses the reason seven figures vanished --
    and the next person to notice the manuscript has no Figure 2 has nothing
    to read. This pins the record to the successor that carries it.
    """
    for retired, successor in sorted(RETIRED.items()):
        e = register.get(successor)
        assert e is not None, (
            f"{successor} is gone from the register, and it held the note "
            f"recording that it replaced {retired}")
        assert retired in (e.get("note") or ""), (
            f"{successor}'s note no longer mentions {retired}. That sentence "
            "is the only surviving record of why that figure was retired")


def test_no_retired_figure_is_drawn_or_registered_again(writes, register):
    """The retired names must not come back without a decision.

    Coming back is not forbidden -- it needs someone to say so by editing
    RETIRED below. What this catches is the silent case: a generator function
    restored from history, drawing to a stem that nothing includes and nothing
    tracks, exactly as before.

    BOTH ROUTES, because the name says both. This took `register` as a
    parameter and never read it, so re-adding a retired figure to FIGURES.yaml
    passed a test called `..._or_registered_again`. An unused fixture is the
    same tell as an unused local: the check it was fetched for is missing.
    """
    redrawn = sorted(set(RETIRED) & set(writes))
    assert not redrawn, (
        f"{redrawn} are drawn again. These were retired because a census "
        "figure replaced each of them; if one is genuinely wanted back, "
        "register it in FIGURES.yaml and remove it from RETIRED here")
    reregistered = sorted(set(RETIRED) & set(register))
    assert not reregistered, (
        f"{reregistered} are back in FIGURES.yaml. Registering a retired "
        "figure is the legitimate way to bring one back, but it is a decision "
        "-- remove it from RETIRED here so the record says so deliberately")
