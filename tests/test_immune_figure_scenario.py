"""The immune figure's panels and its annotations must be the same run.

`sim-tme` writes 33 conditions, nine of which carry `immune_mode: immune_on` --
three scenarios by three treatments:
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


FIXTURE = REPO / "tests/fixtures/flagship_tme_rows.json"


def _conditions():
    """The sim-tme conditions, from the COMMITTED fixture when the live output
    is absent.

    `simulations/output/` is gitignored, so a skip-if-absent guard does not run
    in CI at all -- which is where it matters, and this repository already
    treats a check that reports green while never running as the defect it is.
    `tests/fixtures/flagship_tme_rows.json` exists for exactly this reason
    (`test_flagship_figure_data.py` says so), and now carries the baseline
    Control row this file needs alongside the RSL3 and SDT ones it already had.

    When the live output IS present it is used, and a second guard asserts the
    two agree, so the fixture cannot drift away from the engine unnoticed.
    """
    src = SUMMARY if SUMMARY.exists() else FIXTURE
    return json.loads(src.read_text())["conditions"], src


def _immune_on_rows(rows=None):
    if rows is None:
        rows, _ = _conditions()
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
    # AT THE FUNCTION'S OWN STATEMENT LEVEL. The first version accepted an
    # assert nested one level inside any `If` that was not a literal `False`,
    # which certified `if summary_path.exists(): assert matched` as
    # unconditional -- and a missing summary file then skipped the refusal
    # entirely and drew "(0 immune kills)" on all three panels.
    assert any(isinstance(st, ast.Assert)
               and ast.unparse(st.test).strip() == "not missing"
               for st in fn.body), (
        "fig17 does not refuse at its own statement level when a treatment has "
        "no baseline condition. Nested inside `if summary_path.exists()` the "
        "refusal is skipped when the file is absent, and a per-list `assert "
        "matched` is satisfied by one surviving treatment while the others "
        "fall through to a default of zero -- which renders as a result")
    # And the per-panel lookup must not supply that default.
    assert "imm_kills_map.get(tx_label, 0)" not in flat, (
        "fig17 defaults a missing treatment to zero again")


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


def test_the_committed_figure_draws_the_baseline_counts():
    """The end-to-end property: what is printed on the figure is what the run
    behind its panels produced.

    PER TREATMENT. A first version compared sorted lists, so permuting the map
    -- labelling Control with RSL3's count, RSL3 with SDT's -- passed: the
    multiset was unchanged and every panel carried another treatment's number,
    which is this file's own defect one step further in.
    """
    try:
        import pymupdf
    except ImportError:
        try:
            import fitz as pymupdf
        except ImportError:
            pytest.skip("no PDF reader available")

    rows, src = _conditions()
    base = {r["treatment"]: r.get("immune_kills", 0)
            for r in _baseline_rows(_immune_on_rows(rows))}
    doc = pymupdf.open(FIG)
    try:
        text = " ".join(" ".join(p.get_text().split()) for p in doc)
    finally:
        doc.close()
    drawn = {t: int(n) for t, n in
             re.findall(r"(Control|RSL3|SDT) \((\d+) immune kills\)", text)}
    assert drawn, f"fig17 no longer annotates immune kills: {text[:200]}"
    assert drawn == base, (
        f"the figure draws {drawn} and the baseline condition produced {base} "
        f"(read from {src.name}). If the engine moved, regenerate the figure; "
        "if the mapping moved, a panel carries another run's or another "
        "treatment's number, which is the defect this file exists for")


def test_the_fixture_agrees_with_the_live_engine():
    """The fixture is what CI checks against, so it must not drift."""
    if not SUMMARY.exists():
        pytest.skip("no live sim output here; the fixture is the reference")
    live = _baseline_rows(_immune_on_rows(
        json.loads(SUMMARY.read_text())["conditions"]))
    fixed = _baseline_rows(_immune_on_rows(
        json.loads(FIXTURE.read_text())["conditions"]))
    as_map = lambda rs: {r["treatment"]: r.get("immune_kills") for r in rs}
    assert as_map(live) == as_map(fixed), (
        f"the committed fixture says {as_map(fixed)} and this machine's "
        f"sim-tme says {as_map(live)}; CI checks the figure against the "
        "fixture, so a stale fixture makes that check meaningless")


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


def test_only_the_baseline_section_writes_the_damp_heatmaps():
    """The premise the whole selector rests on, read from the Rust source.

    "The DAMP CSVs are written by the baseline section alone" is why the figure
    must be labelled with the baseline run's counts, and it existed only as
    prose in a docstring and a comment. Moving `write_heatmap_csv` into the
    stromal or pH loop would put the panels and the annotations back on
    different runs with every other guard in this file still green.

    Located by structure rather than by line number, which rots.
    """
    src = (REPO / "simulations/sim-tme/src/main.rs").read_text()
    writes = [i for i, line in enumerate(src.splitlines())
              if "damp_field_" in line and "format!" in line]
    assert len(writes) == 1, (
        f"{len(writes)} sites write damp_field_*.csv; the figure's annotation "
        "assumes exactly one, so which run the panels are of is no longer "
        "determined")
    # It must sit inside the dedicated immune section -- the one that loops
    # over `immune_modes` -- and not inside the stromal or pH sections.
    before = "\n".join(src.splitlines()[:writes[0]])
    last_loop = max(
        (before.rfind("for (immune_label, immune_cfg) in &immune_modes"),
         before.rfind("for (immune_label, immune_cfg) in &stromal_immune_modes")))
    assert last_loop != -1, "cannot locate the enclosing immune loop"
    assert before[last_loop:].startswith(
        "for (immune_label, immune_cfg) in &immune_modes"), (
        "the DAMP heatmap write is no longer inside the baseline immune "
        "section, so the panels are of a different run from the one the "
        "figure's counts are selected from")
    # AND UNDER THE `immune_on` ARM. That loop runs `immune_on` and
    # `immune_anti_pd1`; flipping the branch string alone puts the panels on
    # the anti-PD-1 run (SDT 1,477 immune kills) while the annotations stay on
    # `immune_on` (521) -- the original defect one arm over, and invisible to
    # the end-to-end check, because the drawn numbers would still equal the
    # selected ones.
    guard = before.rfind('if *immune_label == "immune_on"')
    assert guard > last_loop, (
        "the DAMP heatmap write is not guarded by `immune_label == "
        '"immune_on"` inside the baseline loop, so the panels may come from '
        "the anti-PD-1 arm while the annotations come from the immune-on one")


@pytest.mark.parametrize("doctor", ["drop_sdt", "no_summary"])
def test_the_generator_refuses_rather_than_labelling_a_panel_zero(tmp_path, doctor):
    """BEHAVIOURAL, because the structural check cannot see a vacuous refusal.

    An AST guard asserting `assert not missing` appears at statement level is
    satisfied by `missing = []` -- the shape survives and the check is empty.
    So this runs the real function against doctored inputs and requires it to
    raise: once with a summary that has no baseline SDT row, and once with no
    summary at all, which is how a missing file used to skip the refusal
    entirely and draw "(0 immune kills)" on every panel.
    """
    import shutil
    import sys as _sys

    import matplotlib
    matplotlib.use("Agg")
    _sys.path.insert(0, str(REPO / "scripts"))
    import generate_figures as g

    rows, _src = _conditions()
    tme = tmp_path / "tme"
    tme.mkdir()
    for name in ("control", "rsl3", "sdt"):
        src_csv = REPO / f"simulations/output/tme/damp_field_{name}.csv"
        if src_csv.exists():
            shutil.copy(src_csv, tme / src_csv.name)
        else:
            # The CSVs are gitignored; a 2x2 stand-in is enough, because the
            # refusal must fire before any field is read.
            (tme / f"damp_field_{name}.csv").write_text("0,1\n1,0\n")
    if doctor == "drop_sdt":
        kept = [r for r in rows
                if not (r.get("treatment") == "SDT"
                        and r.get("immune_mode") == "immune_on"
                        and r.get("stromal_mode") == "off"
                        and "ph_mode" not in r)]
        (tme / "tme_summary.json").write_text(
            json.dumps({"schema_version": 1, "conditions": kept}))

    out = tmp_path / "figs"
    out.mkdir()
    old_tme, old_fig = g.TME_DIR, g.FIG_DIR
    g.TME_DIR, g.FIG_DIR = tme, out
    try:
        with pytest.raises(AssertionError, match="baseline immune_on"):
            g.fig17_damp_heatmap()
    finally:
        g.TME_DIR, g.FIG_DIR = old_tme, old_fig
    assert not list(out.glob("*.pdf")), (
        "fig17 wrote a figure despite refusing; a partially written artifact "
        "is worse than none, because it looks regenerated")


def test_the_latex_caption_states_the_baseline_counts():
    """The caption is prose beside a measurement, so it is derived here.

    `generate_latex.py`'s fig 14 entry carried "139,641 kills, 539 immune
    kills / 163 kills, 2 immune kills" -- main's figure, one off on the
    ferroptosis count -- while this branch moved the figure to 521 and 5. A
    reviewer found it; nothing in the suite could, and reverting it plus
    regenerating the manifest (the routine next step) left every test green.

    `test_figure_captions_agree.py` compares the markdown and LaTeX captions on
    shared content WORDS, which strips digits, so it cannot see this either.
    """
    rows, src = _conditions()
    base = {r["treatment"]: r for r in _baseline_rows(_immune_on_rows(rows))}
    latex = (REPO / "scripts/generate_latex.py").read_text()
    entry = re.search(r"'14': \('fig17_damp_heatmap', '(.*?)'\),", latex, re.S)
    assert entry, "the fig17 caption entry is gone or was renamed"
    caption = entry.group(1)

    for treatment in ("SDT", "RSL3"):
        row = base[treatment]
        ferro = f"{row['ferroptosis_kills']:,}"
        immune = str(row["immune_kills"])
        assert f"({ferro} kills, {immune} immune kills)" in caption, (
            f"the LaTeX caption does not state {treatment}'s baseline counts "
            f"({ferro} kills, {immune} immune kills) from {src.name}. It is "
            "the caption that ships: release-pdf.yml builds from v1.tex")
    # And it must say which of the three scenarios, for the same reason the
    # figure's own suptitle does.
    assert "baseline run" in caption, (
        "the LaTeX caption does not name the scenario, so the shipped PDF "
        "gives the reader less than the figure does")
