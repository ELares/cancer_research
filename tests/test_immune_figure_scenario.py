"""The immune figure's panels and its annotations must be the same run.

`sim-tme` writes 33 conditions, three of which carry `immune_mode: immune_on` --
the baseline run, a stromal-shielded run, and a pH-gradient run. The DAMP
heatmap CSVs are written by the baseline section alone
(`simulations/sim-tme/src/main.rs:1265`, inside a loop over immune modes), but
`generate_figures.py` looped over every `immune_on` row and kept whichever the
file listed LAST, which is the pH run.

So fig17 rendered the baseline run's DAMP fields and labelled them with the pH
run's kill counts: `0 / 0 / 496` drawn over panels whose own counts are
`0 / 5 / 521`. The manuscript's `104:1` is 521/5 -- the panels' scenario -- so
the figure disagreed with the section it illustrates, and nothing said which
scenario either was.

These guards pin the selection, not the numbers: the counts move whenever the
engine changes, and the property that must hold is that the figure's numbers
come from the run its panels come from.
"""
import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
GEN = REPO / "scripts/generate_figures.py"
SUMMARY = REPO / "simulations/output/tme/tme_summary.json"
FIG = REPO / "article/figures/fig17_damp_heatmap.pdf"


def _immune_on_rows():
    rows = json.loads(SUMMARY.read_text())["conditions"]
    return [r for r in rows if r.get("immune_mode") == "immune_on"]


def _baseline_rows(rows):
    return [r for r in rows
            if r.get("stromal_mode") == "off" and "ph_mode" not in r]


def test_the_generator_selects_the_baseline_condition_by_name():
    """Read from the source, because the defect was an iteration-order accident
    that no output inspection distinguishes from a deliberate choice."""
    src = GEN.read_text()
    body = src[src.index("def fig17_damp_heatmap"):]
    body = body[:body.index("\ndef ")]
    flat = " ".join(body.split())
    assert 'r.get("stromal_mode") == "off"' in flat, (
        "fig17 no longer restricts to the baseline stromal condition; it will "
        "silently take whichever immune_on row the summary lists last")
    assert '"ph_mode" not in r' in flat, (
        "fig17 no longer excludes the pH condition")
    assert "imm_kills_map[r[\"treatment\"]] = r.get(\"immune_kills\", 0)" in flat
    # And it must REFUSE rather than fall back if the condition is absent: a
    # zero here would render as "(0 immune kills)" and read as a result.
    #
    # BY STRUCTURE, not by substring: `if False: assert matched` still contains
    # the text and refuses nothing, which a string check passed.
    import ast

    fn = next(n for n in ast.parse(src).body
              if isinstance(n, ast.FunctionDef) and n.name == "fig17_damp_heatmap")
    def _asserts_matched(body):
        return any(isinstance(st, ast.Assert)
                   and ast.unparse(st.test).strip() == "matched" for st in body)
    live = _asserts_matched(fn.body) or any(
        _asserts_matched(st.body) for st in fn.body
        if isinstance(st, (ast.If, ast.With, ast.For))
        and not (isinstance(st, ast.If)
                 and ast.unparse(st.test).strip() in ("False", "0")))
    assert live, (
        "fig17 does not refuse unconditionally when the baseline condition is "
        "missing; a guarded or removed assertion lets it label the panels with "
        "zeros, which render as a result")


@pytest.mark.skipif(not SUMMARY.exists(),
                    reason="simulations/output is gitignored; see #788")
def test_exactly_one_baseline_condition_exists_per_treatment():
    rows = _immune_on_rows()
    base = _baseline_rows(rows)
    assert len(rows) > len(base), (
        "there is no longer more than one immune_on scenario, so the selector "
        "this file guards is no longer load-bearing -- check whether the "
        "figure still needs it before deleting the guard")
    treatments = [r["treatment"] for r in base]
    assert sorted(treatments) == sorted(set(treatments)), (
        f"the baseline condition is not unique per treatment: {treatments}")


@pytest.mark.skipif(not (SUMMARY.exists() and FIG.exists()),
                    reason="simulations/output is gitignored; see #788")
def test_the_committed_figure_draws_the_baseline_counts():
    """The end-to-end property: what is printed on the figure is what the run
    behind its panels produced."""
    try:
        import pymupdf
    except ImportError:
        try:
            import fitz as pymupdf
        except ImportError:
            pytest.skip("no PDF reader available")

    base = {r["treatment"]: r.get("immune_kills", 0)
            for r in _baseline_rows(_immune_on_rows())}
    doc = pymupdf.open(FIG)
    try:
        text = " ".join(" ".join(p.get_text().split()) for p in doc)
    finally:
        doc.close()
    drawn = [int(n) for n in re.findall(r"\((\d+) immune kills\)", text)]
    assert drawn, f"fig17 no longer annotates immune kills: {text[:200]}"
    assert sorted(drawn) == sorted(base.values()), (
        f"the figure draws {sorted(drawn)} and the baseline condition produced "
        f"{sorted(base.values())}. If the engine moved, regenerate the figure; "
        "if the selector moved, the panels and the numbers are of different "
        "runs, which is the defect this file exists for")


@pytest.mark.skipif(not FIG.exists(), reason="figure absent")
def test_the_figure_names_which_scenario_it_is():
    """Three immune_on scenarios exist and a reader cannot tell them apart from
    a heatmap. The caption has to say."""
    try:
        import pymupdf
    except ImportError:
        try:
            import fitz as pymupdf
        except ImportError:
            pytest.skip("no PDF reader available")

    doc = pymupdf.open(FIG)
    try:
        text = " ".join(" ".join(p.get_text().split()) for p in doc)
    finally:
        doc.close()
    assert "baseline run" in text, (
        "fig17 does not say which of the three immune-coupling scenarios it "
        "shows")
    for absent in ("no stromal shielding", "no pH gradient"):
        assert absent in text, (
            f"fig17's caption no longer states {absent!r}, so a reader cannot "
            "tell it apart from the stromal or pH runs")
