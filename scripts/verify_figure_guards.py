#!/usr/bin/env python3
"""Mutation-verify the figure guards in BOTH directions.

`tests/test_simulation_figures_draw_their_data.py` asserts what fig24, fig25
and fig26 actually DRAW. A guard like that has two ways to be worthless, and a
green suite shows neither:

  * it PASSES on a figure whose meaning is wrong -- the #790 defect class,
    where the extracted text is byte-identical and only the arrangement moved;
  * it FAILS on a figure that is correct -- a plain rcParam, a restyle, or a
    legitimate re-run of the simulation.

This runs both. Each case mutates the GENERATOR (or an input, or an rcParam),
regenerates the three figures into a scratch directory, and records whether the
SEMANTIC tests fail. The drawing fingerprint is excluded throughout: it fires on
any change at all, so it cannot distinguish the two directions and a case it
alone catches is not caught.

WHY THIS IS A LOCAL TOOL, NOT A CI TEST. It regenerates figures, and the
generator reads `simulations/output/`, which is gitignored -- the same reason
the guards themselves compare against committed fixtures rather than
regenerating. CI cannot run this; a person with the simulation outputs can.
That is the contract `scripts/validate_rd_vs_biofvm.py` and
`scripts/fetch_calibration_data.py` already use.

WHY IT IS COMMITTED. Thirty adversarial review rounds produced these cases, and
until now they lived in a scratch directory and were re-derived each round. Two
consequences, both measured: a reviewer could not reproduce the verification a
commit message claimed, and -- worse -- when a review showed a guard rejecting a
CORRECT figure, the fix was made without re-running the mutation that motivated
the original assertion. That produced a vacuous guard three times, including a
bound (`<= len(names)`) that ZERO satisfies, under which a bar chart drawn with
no bars at all passed every semantic assertion.

It deletes the sandbox's figures before the baseline run, so a generator that
silently skips (no `simulations/output/`) is reported as such rather than
passing a hash comparison against files it never touched.

Run it after touching the guards or the generator:

    python3 scripts/verify_figure_guards.py            # all cases
    python3 scripts/verify_figure_guards.py --list     # names only
    python3 scripts/verify_figure_guards.py -k arrow   # a subset

Exit status is non-zero if any case behaves the wrong way, so it can gate a
commit locally. A case that cannot run (its anchor text is gone from the
generator) is reported as SKIP and counted, never silently dropped -- an
anchor that stops matching is how a mutation harness quietly stops testing.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
TESTS = "tests/test_simulation_figures_draw_their_data.py"
FIGS = ("fig24_hypoxia_killcurve", "fig25_bliss_synergy",
        "fig26_vulnerability_window")

# Each case: (name, expect_fail, [(find, replace), ...], rcParams, input_edit)
#   expect_fail True  -- a WRONG figure: a semantic test must reject it
#   expect_fail False -- a CORRECT figure: every semantic test must accept it
CASES: list[tuple] = []


def case(name, expect_fail, subs=(), rc=None, edit=None):
    CASES.append((name, expect_fail, list(subs), rc or {}, edit))


# ---------------------------------------------------------------- inversions
case("fig24 bar offsets swapped", True,
     [("axA.bar(x - w / 2, [rsl3_norm, sdt_norm], w, label=\"Normoxic (uniform O$_2$)\", color=\"#4C72B0\")\n"
       "    axA.bar(x + w / 2, [rsl3_hyp, sdt_hyp], w, label=\"Hypoxic (O$_2$ gradient)\", color=\"#C44E52\")",
       "axA.bar(x + w / 2, [rsl3_norm, sdt_norm], w, label=\"Normoxic (uniform O$_2$)\", color=\"#4C72B0\")\n"
       "    axA.bar(x - w / 2, [rsl3_hyp, sdt_hyp], w, label=\"Hypoxic (O$_2$ gradient)\", color=\"#C44E52\")")])
case("fig24 collapse arrow reversed", True,
     [("{rsl3_norm:.1f}% $\\\\to$ {rsl3_hyp:.1f}%", "{rsl3_norm:.1f}% $\\\\leftarrow$ {rsl3_hyp:.1f}%")])
case("fig24 collapse note over the SDT group", True,
     [("xytext=(0.55, 55)", "xytext=(1.0, 55)")])
case("fig24 group tick labels reversed", True,
     [("axA.set_xticks(x)", "axA.set_xticks(x[::-1])")])
case("fig25 bars de-paired from their scores", True,
     [("axB.barh(names[::-1], scores[::-1], color=cols[::-1])",
       "axB.barh(names[::-1], scores, color=cols[::-1])")])
case("fig25 panel (b) highlight on the wrong pair", True,
     [("axB.barh(names[::-1], scores[::-1], color=cols[::-1])",
       "axB.barh(names[::-1], scores[::-1], color=cols)")])
case("fig25 panel (a) result colour on the prediction", True,
     [('colors = ["#4C72B0", "#55A868", "#999999", "#C44E52"]',
       'colors = ["#4C72B0", "#55A868", "#C44E52", "#999999"]')])
case("fig25 synergy arrow onto the monotherapy bar", True,
     [("xy=(3, vals[3]), xytext=(1.9, 93)", "xy=(0, vals[0]), xytext=(1.9, 93)")])
case("fig25 additive line moved to 1.9", True,
     [("axB.axvline(1.0", "axB.axvline(1.9")])
case("fig26 window caption names SDT", True,
     [('"RSL3 window\\nopen"', '"SDT window\\nopen"')])
case("fig26 window shade in the SDT colour", True,
     [('axA.axvspan(-0.3, win_end, color="#4C72B0"',
       'axA.axvspan(-0.3, win_end, color="#C44E52"')])
case("fig26 closes-arrow onto the last timepoint", True,
     [("xy=(win_end, rsl3[win_end])", "xy=(len(days) - 1, rsl3[-1])")])
case("fig26 closes-words at the far end", True,
     [("xytext=(win_end + 0.4, 24)", "xytext=(len(days) - 1.4, 24)")])

# ------------------------------------------------------- degenerate figures
case("fig25 every bar zero-length, scores still printed", True,
     [("axB.barh(names[::-1], scores[::-1], color=cols[::-1])",
       "axB.barh(names[::-1], [0.0] * len(scores), color=cols[::-1])")])
case("fig25 the flagship bar alone missing", True,
     [("axB.barh(names[::-1], scores[::-1], color=cols[::-1])",
       "axB.barh(names[::-1], [0.0 if n == \"RSL3+FSP1i\" else s "
       "for n, s in zip(names[::-1], scores[::-1])], color=cols[::-1])")])

# ------------------------ escapes found by review, each a one-token variant
# Round 30 showed the arrow bindings defeated four separate ways, each leaving
# the annotation pointing somewhere it does not belong. They are cases here so
# that a future relaxation of the arrow reading re-runs them.
case("fig25 headless arrow onto the monotherapy bar", True,
     [('arrowprops=dict(arrowstyle="->", color="#C44E52")',
       'arrowprops=dict(arrowstyle="-", color="#C44E52")'),
      ("xy=(3, vals[3]), xytext=(1.9, 93)", "xy=(0, vals[0]), xytext=(1.9, 93)")])
case("fig25 headless arrow onto the prediction bar", True,
     [('arrowprops=dict(arrowstyle="->", color="#C44E52")',
       'arrowprops=dict(arrowstyle="-", color="#C44E52")'),
      ("xy=(3, vals[3]), xytext=(1.9, 93)", "xy=(2, vals[2]), xytext=(1.9, 93)")])
case("fig26 headless arrow onto the last timepoint", True,
     [('arrowprops=dict(arrowstyle="->", color="#4C72B0")',
       'arrowprops=dict(arrowstyle="-", color="#4C72B0")'),
      ("xy=(win_end, rsl3[win_end])", "xy=(len(days) - 1, rsl3[-1])")])
case("fig26 triangle markers hide the arrow's tip", True,
     [('axA.plot(x, rsl3, "s-"', 'axA.plot(x, rsl3, "^-"'),
      ("xy=(win_end, rsl3[win_end])", "xy=(len(days) - 1, rsl3[-1])")])
case("fig26 closes-arrow onto day 0", True,
     [("xy=(win_end, rsl3[win_end])", "xy=(0, rsl3[0])")])
case("fig26 closes-arrow onto the SDT curve", True,
     [("xy=(win_end, rsl3[win_end])", "xy=(win_end, sdt[win_end])")])
case("fig24 collapse note onto a value-colliding group", True,
     [("xytext=(0.55, 55)", "xytext=(0.92, 55)")],
     edit=lambda root: _sdt_uniform_equals_rsl3(root))
case("fig25 colour swap with a bar shrunk below the size cut", True,
     [('colors = ["#4C72B0", "#55A868", "#999999", "#C44E52"]',
       'colors = ["#4C72B0", "#55A868", "#C44E52", "#999999"]')],
     edit=lambda root: _shrink_rate_b(root))

# --------------------------- correct figures whose ARROWS are drawn oddly
# `_annotation_arrows` identifies the arrow by PROVENANCE, not appearance: a
# path with an endpoint beside the annotation's words, not enclosing them, not
# running frame-to-frame, and last in paint order (matplotlib draws an
# annotation after the artists it annotates). These are the ways an arrow or
# its text legitimately gains geometry; each defeated one of the six earlier
# appearance-based rules.
case("fig25 arrow with connectionstyle arc3", False,
     [('arrowprops=dict(arrowstyle="->", color="#C44E52")',
       'arrowprops=dict(arrowstyle="->", color="#C44E52", '
       'connectionstyle="arc3,rad=0.4")')])
case("fig25 arrow with connectionstyle angle3", False,
     [('arrowprops=dict(arrowstyle="->", color="#C44E52")',
       'arrowprops=dict(arrowstyle="->", color="#C44E52", '
       'connectionstyle="angle3")')])
case("fig26 arrow with connectionstyle arc3", False,
     [('arrowprops=dict(arrowstyle="->", color="#4C72B0")',
       'arrowprops=dict(arrowstyle="->", color="#4C72B0", '
       'connectionstyle="arc3,rad=0.3")')])
case("fig25 annotation in a rounded bbox", False,
     [('fontsize=10, fontweight="bold", color="#C44E52"',
       'fontsize=10, fontweight="bold", color="#C44E52", '
       'bbox=dict(boxstyle="round", fc="w")')])
case("legend.fancybox", False, rc={"legend.fancybox": True})

# ---------------- arrow styles and connection styles, all correct figures
# Round 31 showed a shape-based arrow rule rejecting six of these. The rule is
# "it starts at the annotation's words" now, and these pin that it stays so.
for _style in ("simple", "fancy", "wedge"):
    case(f"fig25 arrowstyle={_style}", False,
         [('arrowprops=dict(arrowstyle="->", color="#C44E52")',
           f'arrowprops=dict(arrowstyle="{_style}", color="#C44E52")')])
for _conn in ("angle,angleA=0,angleB=90,rad=5", "bar,fraction=0.2"):
    case(f"fig26 connectionstyle={_conn.split(',')[0]}", False,
         [('arrowprops=dict(arrowstyle="->", color="#4C72B0")',
           f'arrowprops=dict(arrowstyle="->", color="#4C72B0", '
           f'connectionstyle="{_conn}")')])
case("fig24 bar width 0.45", False, [("w = 0.36", "w = 0.45")])
case("fig24 bar width 0.55 interleaves the groups", True,
     [("w = 0.36", "w = 0.55")])
case("fig24 wide bars hide the note over SDT", True,
     [("w = 0.36", "w = 0.55"), ("xytext=(0.55, 55)", "xytext=(1.0, 55)")])
case("fig26 arrow at the right day, on neither curve", True,
     [("xy=(win_end, rsl3[win_end])", "xy=(win_end, 40)")])

# ------------------- arrows pointing AT the frame, which must still be found
# `_annotation_arrows` drops a path with BOTH endpoints on the panel frame, to
# exclude a gridline drawn through the annotation. An arrow pointing at a value
# that happens to lie on the frame -- a zero, which sits on the x axis -- has
# only its target end there and must survive. Measured: both do.
case("fig26 arrow at a zero value, on the axis", False,
     [("xy=(win_end, rsl3[win_end])", "xy=(win_end, 0)")])
case("fig25 arrow at the observed bar's base", False,
     [("xy=(3, vals[3]), xytext=(1.9, 93)", "xy=(3, 0), xytext=(1.9, 93)")])

# ------------------ round 32: artists that out-reach or out-anchor the arrow
# A bar rectangle reaches the panel floor and a data curve reaches its own last
# timepoint, so both dwarf an arrow. On fig25 the winner was the very bar the
# assertion then compared against itself, which cannot fail.
case("fig25 arrow onto the monotherapy, taller combo bar", True,
     [("xy=(3, vals[3]), xytext=(1.9, 93)", "xy=(0, vals[0]), xytext=(1.9, 93)")],
     edit=lambda root: _set_rate_combo(root, 0.90))
case("fig25 taller combo bar alone", False, edit=lambda root: _set_rate_combo(root, 0.90))
case("fig25 arrow onto the monotherapy, annotation lowered", True,
     [("xy=(3, vals[3]), xytext=(1.9, 93)", "xy=(0, vals[0]), xytext=(1.9, 90)")])
case("fig26 annotation at day 0, level with the curve", False,
     [("xytext=(win_end + 0.4, 24)", "xytext=(0.0, 40)")])
case("fig26 annotation at day 0, above the curve", False,
     [("xytext=(win_end + 0.4, 24)", "xytext=(0.0, 47)")])

# ------------------------------- artists drawn AFTER the annotation
# The arrow is taken as the last candidate in paint order, so anything drawn
# later is the obvious attack. fig26's generator ALREADY draws its legend after
# the `closes` annotation.
#
# THESE THREE CASES DO NOT DECIDE THAT, and the commit that added them said
# they did. Measured: `lower left` puts the legend frame 69.2pt from the
# annotation's words and `lower center` 28.0pt, both outside the 12.65pt
# admission threshold -- so no legend path is ever a candidate and paint order
# is never consulted for it. They are kept because they pin the placements,
# and the cases that actually reach the legend are the round-33 group below:
# there the legend's marker SWATCH (a `re` a few points across, not the frame)
# was admitted and won.
case("fig26 legend over the annotation", False,
     [('axA.legend(fontsize=8, loc="center right")',
       'axA.legend(fontsize=8, loc="lower left")')])
case("fig26 legend over the annotation, arrow to the last day", True,
     [('axA.legend(fontsize=8, loc="center right")',
       'axA.legend(fontsize=8, loc="lower left")'),
      ("xy=(win_end, rsl3[win_end])", "xy=(len(days) - 1, rsl3[-1])")])
case("fig26 legend below the annotation, arrow to the last day", True,
     [('axA.legend(fontsize=8, loc="center right")',
       'axA.legend(fontsize=8, loc="lower center")'),
      ("xy=(win_end, rsl3[win_end])", "xy=(len(days) - 1, rsl3[-1])")])

# ------------------ round 33: the annotation's row, bounded to its own words
# A row scan with no horizontal bound swept the LEGEND'S LABELS into "the
# annotation", so `near` slid beside the legend's marker swatches -- painted
# after the annotation, so paint order promoted a 6pt swatch to "the arrow" and
# `closes ~day 3` pointing at day 28 passed. The mirror rejected correct
# figures. Both sites walk gaps now, as fig25's sibling has since round 22.
case("fig26 annotation raised clear of the legend", False,
     [("xytext=(win_end + 0.4, 24)", "xytext=(win_end + 0.4, 60)")])
case("fig26 annotation level with the legend top", False,
     [("xytext=(win_end + 0.4, 24)", "xytext=(win_end + 0.4, 55)")])
case("fig26 legend labels on the caption's row", False, rc={"legend.fontsize": 16})

# ------------------------------- the gap walk's own edges (round 34, pre-empted)
# Bounding the annotation's row by a gap walk introduces its own risks: the
# word it keys on can vanish or repeat, and the annotation can WRAP, which puts
# that word alone on its row. The wrapped case rejected a correct figure -- the
# one-word-box defect of round 31 returning by another route -- so the block
# now takes adjacent overlapping rows the way a wrapped title is taken.
case("fig26 annotation wrapped to two lines", False,
     [('"closes ~day 3"', '"closes\\n~day 3"')])
case("fig26 annotation renamed away from `closes`", True,
     [('"closes ~day 3"', '"shuts ~day 3"')])
case("fig26 caption renamed away from `window`", True,
     [('"RSL3 window\\nopen"', '"RSL3 period\\nopen"')])
case("fig26 font.size 14, wider word gaps", False, rc={"font.size": 14})
case("fig26 font.size 7, tighter word gaps", False, rc={"font.size": 7})

# ------------------------------------------------- correct figures, restyled
case("committed figures, unmodified", False)
case("xtick.direction inout", False, rc={"xtick.direction": "inout"})
case("tick sizes 20", False, rc={"xtick.major.size": 20, "ytick.major.size": 20})
case("titlepad 0", False, rc={"axes.titlepad": 0})
case("ytick.alignment center", False, rc={"ytick.alignment": "center"})
case("xtick.alignment left", False, rc={"xtick.alignment": "left"})
case("axes.xmargin 0", False, rc={"axes.xmargin": 0})
case("every spine off individually (bottom)", False, rc={"axes.spines.bottom": False})
case("axes.linewidth 0", False, rc={"axes.linewidth": 0})
case("dashed grid", False, rc={"axes.grid": True, "grid.linestyle": "--"})
case("mathtext stixsans", False, rc={"mathtext.fontset": "stixsans"})
case("labelsize 16", False, rc={"axes.labelsize": 16})
case("fig25 width_ratios 1.6:1", False,
     [("fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.5))",
       "fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.5), "
       "gridspec_kw={\"width_ratios\": [1.6, 1]})")])
case("fig25 annotation fontsize 18", False,
     [('fontsize=10, fontweight="bold", color="#C44E52"',
       'fontsize=18, fontweight="bold", color="#C44E52"')])
case("fig24 ylim starting below zero", False,
     [("axA.set_ylim(0, 105)", "axA.set_ylim(-5, 105)")])
case("fig25 xlim with no headroom", False,
     [("axB.set_xlim(0, max(scores) * 1.25)", "axB.set_xlim(0, max(scores))")])
case("fig24 title wrapped to two lines", False,
     [('axA.set_title("(a) Kill collapse under hypoxia")',
       'axA.set_title("(a) Kill collapse\\nunder hypoxia")')])
case("fig26 both titles wrapped", False,
     [('axA.set_title("(a) Treatment window: RSL3 closes, SDT stays open")',
       'axA.set_title("(a) Treatment window:\\nRSL3 closes, SDT stays open")'),
      ('axB.set_title("(b) Why: GPX4 re-expression closes the window")',
       'axB.set_title("(b) Why: GPX4 re-expression\\ncloses the window")')])
case("fig26 closes-words at day 0, inside the window", False,
     [("xytext=(win_end + 0.4, 24)", "xytext=(0.0, 24)")])
case("fig24 legend below the axes", False,
     [('axA.legend(fontsize=8, loc="upper left")',
       'axA.legend(fontsize=8, loc="upper center", ncol=2, '
       'bbox_to_anchor=(0.5, -0.22))')])
case("fig24 legend outside, to the right", False,
     [('axA.legend(fontsize=8, loc="upper left")',
       'axA.legend(fontsize=8, loc="center left", bbox_to_anchor=(1.02, 0.5))')])


# ------------------------------------------------------- legitimate re-runs
def _sdt_o2_independent(root: Path) -> None:
    """SDT's gradient rates equal to its uniform rate.

    fig24's own footnote says the model treats SDT as O2-independent, so this
    is a legitimate run, not a defect. It also makes two bar annotations read
    the same value, which is what collapsed a text-keyed lookup.
    """
    for rel in ("tests/fixtures/hypoxia_killcurve_rows.json",
                "simulations/output/tme/tme_summary.json"):
        p = root / rel
        if not p.exists():
            continue
        d = json.loads(p.read_text())
        conds = d["conditions"] if isinstance(d, dict) and "conditions" in d else d
        uni = next(c["overall_kill_rate"] for c in conds
                   if c["treatment"] == "SDT" and c["o2_condition"] == "uniform"
                   and c.get("immune_mode") == "off")
        for c in conds:
            if (c["treatment"] == "SDT" and c["o2_condition"] != "uniform"
                    and c.get("immune_mode") == "off"):
                c["overall_kill_rate"] = uni
        p.write_text(json.dumps(d, indent=1))


def _sdt_uniform_equals_rsl3(root: Path) -> None:
    """SDT's uniform rate set equal to RSL3's, so two bar annotations collide.

    A legitimate value, and it defeated a membership test keyed on annotation
    TEXT: SDT's leftmost label then read the same as one the collapse note
    quotes, dropped out of the "other group", and the boundary moved 53pt.
    """
    for rel in ("tests/fixtures/hypoxia_killcurve_rows.json",
                "simulations/output/tme/tme_summary.json"):
        p = root / rel
        if not p.exists():
            continue
        d = json.loads(p.read_text())
        conds = d["conditions"] if isinstance(d, dict) and "conditions" in d else d
        rsl3 = next(c["overall_kill_rate"] for c in conds
                    if c["treatment"] == "RSL3" and c["o2_condition"] == "uniform"
                    and c.get("immune_mode") == "off")
        for c in conds:
            if (c["treatment"] == "SDT" and c["o2_condition"] == "uniform"
                    and c.get("immune_mode") == "off"):
                c["overall_kill_rate"] = rsl3
        p.write_text(json.dumps(d, indent=1))


def _shrink_rate_b(root: Path) -> None:
    """FSP1i-alone below the height cut that `a_bars` applies.

    3.7% is 9.2pt today; at 1.5% the bar drops out of the rectangle scan, and
    a colour check gated on "exactly four bars" switched itself off.

    BOTH THE RUN AND THE FIXTURE, and an earlier version edited only the run.
    `tests/fixtures/bliss_synergy.json` is shaped `{_comment, rsl3_fsp1i, ...}`
    -- no `combinations` key and no top-level `rate_b` -- so the loop matched
    nothing there while still rewriting the file. The figure then disagreed
    with its own fixture, the DRIFT assertion fired first, and this case was
    reported as caught while the colour guard it exists for was never reached.
    A control (the same edit with the colour swap removed) failed identically,
    which is how that was established.

    `bliss_prediction` and `synergy_score` are recomputed, because a run with
    rate_a 40%, rate_b 1.5% and a stale 42.2% prediction is not the
    "legitimate re-run" this case claims to be.
    """
    def _fix(c):
        c["rate_b"] = 0.015
        a, b = c["rate_a"], c["rate_b"]
        c["bliss_prediction"] = a + b - a * b
        if c.get("rate_combo") is not None and c["bliss_prediction"]:
            c["synergy_score"] = c["rate_combo"] / c["bliss_prediction"]

    run = root / "simulations/output/combo-mech/combo_summary.json"
    if run.exists():
        d = json.loads(run.read_text())
        for c in d["combinations"]:
            if (c.get("drug_a"), c.get("drug_b")) == ("RSL3", "FSP1i"):
                _fix(c)
        run.write_text(json.dumps(d, indent=1))

    fix = root / "tests/fixtures/bliss_synergy.json"
    if fix.exists():
        d = json.loads(fix.read_text())
        _fix(d["rsl3_fsp1i"])
        fix.write_text(json.dumps(d, indent=1))


def _set_rate_combo(root: Path, value: float) -> None:
    """RSL3+FSP1i's observed combination rate, in the run AND the fixture.

    A legitimate re-run: a higher observed kill raises the bar's top, which is
    what brought a bar rectangle within reach of the annotation and let it be
    read as the arrow. `bliss_prediction` and `synergy_score` are recomputed so
    the figure and its fixture still agree and the drift guard is not what
    fires.
    """
    for rel, key in (("simulations/output/combo-mech/combo_summary.json",
                      "combinations"),
                     ("tests/fixtures/bliss_synergy.json", None)):
        p = root / rel
        if not p.exists():
            continue
        d = json.loads(p.read_text())
        rows = d[key] if key else [d["rsl3_fsp1i"]]
        for c in rows:
            if key and (c.get("drug_a"), c.get("drug_b")) != ("RSL3", "FSP1i"):
                continue
            c["rate_combo"] = value
            a, b = c["rate_a"], c["rate_b"]
            c["bliss_prediction"] = a + b - a * b
            c["synergy_score"] = value / c["bliss_prediction"]
        p.write_text(json.dumps(d, indent=1))


def _extra_pairs(n: int):
    def edit(root: Path) -> None:
        p = root / "simulations/output/combo-mech/combo_summary.json"
        if not p.exists():
            return
        d = json.loads(p.read_text())
        combos = d["combinations"]
        base, i = combos[0], 0
        while len(combos) < n:
            c = dict(base)
            c["drug_a"], c["drug_b"] = f"DrugX{i}", "HDACi"
            combos.append(c)
            i += 1
        p.write_text(json.dumps(d, indent=1))
    return edit


case("SDT modelled O2-independent (two bars share a height)", False,
     edit=_sdt_o2_independent)
case("fig25 with 15 pairs", False, edit=_extra_pairs(15))
case("fig25 with 29 pairs", False, edit=_extra_pairs(29))


# ------------------------------------------------------------------- runner
def _sandbox(tmp: Path) -> Path:
    root = tmp / "sb"
    root.mkdir()
    for d in ("scripts", "tests", "simulations"):
        src = REPO / d
        if src.exists():
            shutil.copytree(src, root / d, symlinks=True,
                            ignore=shutil.ignore_patterns("target", "__pycache__"))
    figs = root / "article" / "figures"
    figs.mkdir(parents=True)
    for stem in FIGS:
        shutil.copy(REPO / "article" / "figures" / f"{stem}.pdf", figs)
    (root / "regen.py").write_text(
        "import os, sys, json, io, contextlib\n"
        "R = os.path.dirname(os.path.abspath(__file__))\n"
        "sys.path.insert(0, os.path.join(R, 'scripts'))\n"
        "os.environ['FERRO_FIG_DIR'] = os.path.join(R, 'article', 'figures')\n"
        "import matplotlib; matplotlib.use('Agg')\n"
        "import matplotlib.pyplot as plt\n"
        "import generate_figures as gf\n"
        "rc = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}\n"
        "plt.rcParams.update(rc)\n"
        "with contextlib.redirect_stdout(io.StringIO()):\n"
        "    for n in %r:\n"
        "        getattr(gf, n)()\n" % (FIGS,))
    return root


def _bit_identical(root: Path):
    """Did regeneration actually WRITE figures matching the committed ones?

    THE FILES ARE DELETED BEFORE THE BASELINE RUN, and that is the point.
    Comparing hashes without deleting is satisfied by regeneration never
    happening: with `simulations/output/` absent -- a fresh clone, and CI --
    the generator prints "not found, skipping" and returns, the sandbox's
    copies of the committed PDFs sit there untouched, and they are of course
    identical to themselves. Measured: the harness then ran all 39 cases
    against unmutated figures and reported 15 wrong-figure cases as MISSED
    DEFECTS and one correct-figure case as a false rejection -- true
    statements about that run, and completely misleading about the guards.

    A check that ABSENCE satisfies is the defect class this harness exists to
    catch, and it appeared here first. Worth saying plainly.
    """
    import hashlib
    figs = root / "article" / "figures"
    missing = [st for st in FIGS if not (figs / f"{st}.pdf").exists()]
    if missing:
        return False, (
            "the generator wrote no figures (" + ", ".join(missing) + ").\n"
            "It reads simulations/output/, which is gitignored, so a fresh clone\n"
            "cannot run this harness -- see the module docstring.")
    for stem in FIGS:
        a = hashlib.sha256((REPO / "article" / "figures" / f"{stem}.pdf").read_bytes())
        b = hashlib.sha256((figs / f"{stem}.pdf").read_bytes())
        if a.hexdigest() != b.hexdigest():
            return False, (
                f"{stem}.pdf regenerates differently from the committed copy, so a\n"
                "mutation's effect cannot be told from drift already present.\n"
                "Regenerate and commit the figures first.")
    return True, ""


def run_case(root: Path, case_tuple) -> tuple[str, str]:
    name, expect_fail, subs, rc, edit = case_tuple
    gen = root / "scripts" / "generate_figures.py"
    saved = {p: p.read_text() for p in [gen]}
    # EVERY FILE AN `edit=` CALLBACK CAN TOUCH. `bliss_synergy.json` was
    # missing from this list while a callback rewrote it -- harmless only
    # because that rewrite happened to be value-preserving.
    for rel in ("tests/fixtures/hypoxia_killcurve_rows.json",
                "tests/fixtures/bliss_synergy.json",
                "tests/fixtures/vulnerability_window.json",
                "simulations/output/tme/tme_summary.json",
                "simulations/output/combo-mech/combo_summary.json"):
        p = root / rel
        if p.exists():
            saved[p] = p.read_text()
    try:
        src = gen.read_text()
        for find, repl in subs:
            if find not in src:
                return "SKIP", f"anchor absent: {find[:48]!r}"
            src = src.replace(find, repl, 1)
        gen.write_text(src)
        if edit:
            # COMPARE PARSED VALUES, NOT BYTES. Every callback re-serialises
            # with `json.dumps(indent=1)` while four of the five tracked
            # inputs are committed at `indent=2`, so a byte comparison was
            # non-empty unconditionally and this check could never fire --
            # verified by neutering a callback's value loop so it changed
            # nothing at all, which still reported OK. The guard against a
            # vacuous case was itself vacuous, which is the defect it exists
            # to report.
            before = {p: json.loads(v) for p, v in saved.items()
                      if p != gen and p.suffix == ".json"}
            edit(root)
            changed = [p for p, v in before.items()
                       if json.loads(p.read_text()) != v]
            if not changed:
                return "SKIP", ("the edit callback changed no tracked value; "
                                "this case tests nothing")
        r = subprocess.run([sys.executable, str(root / "regen.py"), json.dumps(rc)],
                           capture_output=True, text=True, cwd=root)
        if r.returncode:
            tail = (r.stderr.strip().splitlines() or [""])[-1]
            return "SKIP", f"regeneration failed: {tail[:60]}"
        t = subprocess.run(
            [sys.executable, "-m", "pytest", TESTS, "-q", "-p", "no:randomly",
             "-k", "not moved_unnoticed"],
            capture_output=True, text=True, cwd=root)
        failed = t.returncode != 0
        why = next((l.strip()[2:].strip() for l in t.stdout.splitlines()
                    if l.startswith("E ")), "")
        if failed == expect_fail:
            return "OK", why[:96] if failed else ""
        return "WRONG", (
            "a wrong figure PASSED every semantic test" if expect_fail
            else f"a correct figure was REJECTED: {why[:80]}")
    finally:
        for p, text in saved.items():
            p.write_text(text)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("-k", metavar="SUBSTR", help="only cases whose name contains this")
    ap.add_argument("--list", action="store_true", help="list case names and exit")
    args = ap.parse_args()

    selected = [c for c in CASES if not args.k or args.k.lower() in c[0].lower()]
    if args.list:
        for name, expect_fail, *_ in selected:
            print(f"  {'wrong-figure' if expect_fail else 'correct-figure'}  {name}")
        return 0
    if not selected:
        print(f"no case matches {args.k!r}")
        return 2

    with tempfile.TemporaryDirectory() as tmp:
        root = _sandbox(Path(tmp))
        # DELETE FIRST, so "the generator skipped" cannot look like "the
        # generator reproduced them exactly".
        for stem in FIGS:
            (root / "article" / "figures" / f"{stem}.pdf").unlink()
        r = subprocess.run([sys.executable, str(root / "regen.py"), "{}"],
                           capture_output=True, text=True, cwd=root)
        if r.returncode:
            print("REFUSING TO RUN: the generator errored.")
            print((r.stderr.strip().splitlines() or [""])[-1])
            return 2
        ok, why = _bit_identical(root)
        if not ok:
            print(f"REFUSING TO RUN: {why}")
            return 2
        print(f"baseline regeneration is bit-identical to the committed PDFs\n"
              f"running {len(selected)} cases\n")
        counts = {"OK": 0, "WRONG": 0, "SKIP": 0}
        for c in selected:
            verdict, detail = run_case(root, c)
            counts[verdict] += 1
            kind = "wrong " if c[1] else "correct"
            print(f"  [{verdict:5}] {kind}  {c[0]}")
            if detail and verdict != "OK":
                print(f"            {detail}")
            elif detail:
                print(f"            caught by: {detail}")
        print(f"\n{counts['OK']} as expected, {counts['WRONG']} wrong, "
              f"{counts['SKIP']} skipped")
        if counts["SKIP"]:
            print("a SKIP means a case could not run, so that mutation is no longer\n"
                  "being tested. Two causes: its anchor text is gone from the\n"
                  "generator, or its edit callback changed no tracked value. The\n"
                  "line above says which.")
        return 1 if counts["WRONG"] or counts["SKIP"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
