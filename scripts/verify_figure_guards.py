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


def _bit_identical(root: Path) -> bool:
    import hashlib
    for stem in FIGS:
        a = hashlib.sha256((REPO / "article" / "figures" / f"{stem}.pdf").read_bytes())
        b = hashlib.sha256((root / "article" / "figures" / f"{stem}.pdf").read_bytes())
        if a.hexdigest() != b.hexdigest():
            return False
    return True


def run_case(root: Path, case_tuple) -> tuple[str, str]:
    name, expect_fail, subs, rc, edit = case_tuple
    gen = root / "scripts" / "generate_figures.py"
    saved = {p: p.read_text() for p in [gen]}
    for rel in ("tests/fixtures/hypoxia_killcurve_rows.json",
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
            edit(root)
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
        r = subprocess.run([sys.executable, str(root / "regen.py"), "{}"],
                           capture_output=True, text=True, cwd=root)
        if r.returncode:
            print("cannot regenerate the figures at all -- is simulations/output/ present?")
            print((r.stderr.strip().splitlines() or [""])[-1])
            return 2
        if not _bit_identical(root):
            print("REFUSING TO RUN: regeneration is not bit-identical to the committed\n"
                  "PDFs, so a mutation's effect cannot be told from the baseline's.\n"
                  "Regenerate and commit the figures first.")
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
            print("a SKIP means a case could not run -- its anchor text is gone from\n"
                  "the generator, so that mutation is no longer being tested")
        return 1 if counts["WRONG"] or counts["SKIP"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
