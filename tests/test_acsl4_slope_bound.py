"""Guards for the ACSL4 slope bound (#934).

What needs protecting is not the bound's value but the reasoning that produces
it: two requirements that must BOTH hold, a grid that must CONTAIN the crossing
it reports, an inversion whose order flips, a bridge that must still agree with
the Rust it claims to reproduce, and a refusal to fit a slope that the page must
not quietly start doing. Every guard on a derived field drives `assemble` on
doctored measurements rather than reading the committed artifact.
"""

from __future__ import annotations

import ast
import copy
import importlib.util
import json
import re
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
SCRIPTS = REPO / "scripts"
CORE = REPO / "simulations" / "ferroptosis-core" / "src"


def _load(name):
    sys.path.insert(0, str(SCRIPTS))
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


M = _load("acsl4_slope_bound")
SCRIPT = SCRIPTS / "acsl4_slope_bound.py"
PARAMS_RS = (CORE / "params.rs").read_text()
ACSL4_RS = (CORE / "acsl4.rs").read_text()


@pytest.fixture(scope="module")
def committed():
    return json.loads(M.OUT_JSON.read_text())


@pytest.fixture(scope="module")
def raw(committed):
    return {"engine": copy.deepcopy(committed["engine"]),
            "cell_lines": M.measure_cell_lines()}


def _doctor(raw):
    return {"engine": {p: [dict(r) for r in rows] for p, rows in raw["engine"].items()},
            "cell_lines": copy.deepcopy(raw["cell_lines"])}


# ------------------------------------------------------------------ freshness

def test_the_committed_json_is_what_the_generator_produces(raw, committed):
    assert M.assemble(raw) == committed


def test_the_committed_markdown_is_what_the_renderer_produces(committed):
    assert M.OUT_MD.read_text() == M.render(committed)


def test_the_render_only_path_reassembles(raw):
    tree = ast.parse(SCRIPT.read_text())
    main = next(n for n in tree.body
                if isinstance(n, ast.FunctionDef) and n.name == "main")
    branch = next(n for n in ast.walk(main)
                  if isinstance(n, ast.If) and isinstance(n.test, ast.Attribute)
                  and n.test.attr == "render_only")
    assert any(isinstance(c, ast.Call) and getattr(c.func, "id", None) == "assemble"
               for stmt in branch.body for c in ast.walk(stmt))


def test_the_cell_line_leg_is_remeasured_on_every_render(committed):
    """It reads a committed CSV, so it must never be taken from the artifact.
    `_split_stored` recovers only the engine leg."""
    recovered = M._split_stored(committed)[0]
    assert recovered["engine"] == committed["engine"]
    assert recovered["cell_lines"] == M.measure_cell_lines()
    assert "cell_lines" not in json.loads(M.OUT_JSON.read_text()).get("engine", {})
    # STRUCTURAL, and it has to be: on clean data an echo of the stored bins and
    # a fresh measurement are EQUAL, so a value check cannot tell them apart. A
    # sweep rewrote `_split_stored` to rebuild the leg from the artifact and the
    # value assertions above all passed. The call must be in THAT function.
    tree = ast.parse(SCRIPT.read_text())
    fn = next(n for n in tree.body
              if isinstance(n, ast.FunctionDef) and n.name == "_split_stored")
    calls = [getattr(c.func, "id", None) for c in ast.walk(fn)
             if isinstance(c, ast.Call)]
    assert "measure_cell_lines" in calls, (
        "_split_stored does not re-measure the cell-line leg, so --render-only "
        "can echo a stored copy of a file that is in the repository")
    names = {n.id for n in ast.walk(fn) if isinstance(n, ast.Name)}
    assert "measured_bins" not in names, (
        "_split_stored reads the artifact's own derived bins")


# ------------------------------------------------------------ the constraint

def test_the_baseline_constraint_has_one_python_source():
    """Re-declaring the number is how two copies disagree."""
    assert "from headline_at_fitted import BASELINE_MAX" in SCRIPT.read_text()
    assert not re.search(r"^BASELINE_MAX\s*=", SCRIPT.read_text(), re.M)
    from headline_at_fitted import BASELINE_MAX
    assert M.BASELINE_MAX == BASELINE_MAX


def test_the_baseline_constraint_still_matches_the_manuscript():
    """The sibling's guard parses it; this one refuses to drift from the prose
    independently, because the bound's lower end rests on it."""
    text = (REPO / "article" / "drafts" / "v1.md").read_text()
    m = re.search(r"all phenotypes show less than (\d+(?:\.\d+)?)% death under Control",
                  text)
    assert m, "the manuscript no longer states the untreated-death constraint"
    assert M.BASELINE_MAX == pytest.approx(float(m.group(1)) / 100.0)


# ------------------------------------------------------------ the bridge itself

def test_the_bridge_reproduces_the_rust_constants_it_claims_to():
    """PARSED from `acsl4.rs`, not restated.

    The Python bridge exists only to sweep the Rust one. If the Rust constants
    move, this page is describing a model that no longer exists, and a
    hand-copied table would not notice.
    """
    constants = dict(re.findall(
        r"pub const ACSL4_(NEGATIVE|LOW|NORMAL|HIGH): f64 = ([0-9.]+);", ACSL4_RS))
    assert set(constants) == {"NEGATIVE", "LOW", "NORMAL", "HIGH"}
    documented = dict(re.findall(
        r"z\s*=\s*([+-]?\d+)\s*.{0,4}\s*\[`ACSL4_(NEGATIVE|LOW|NORMAL|HIGH)`\]",
        ACSL4_RS))
    assert documented, "acsl4.rs no longer documents its z-to-constant mapping"
    for z, name in documented.items():
        assert M.status_from_zscore(float(z)) == pytest.approx(float(constants[name])), (
            f"the Python bridge disagrees with acsl4.rs at z={z}")


def test_the_boost_is_the_status_minus_one_with_a_floor():
    assert M.boost_from_zscore(0.0) == 0.0, "the baseline must stay inert"
    assert M.boost_from_zscore(1.0) == pytest.approx(0.5)
    assert M.boost_from_zscore(-2.0) == pytest.approx(-1.0)
    # The floor binds below the collapse point rather than going more negative.
    assert M.boost_from_zscore(-10.0) == pytest.approx(-1.0)


def test_the_slope_enters_as_z_over_slope(committed):
    """The inversion the bound depends on: at the top of the bridge's range the
    boost is one over the slope."""
    top = committed["bridge"]["z_high"]
    for slope in (1.0, 2.0, 4.0):
        assert M.boost_from_zscore(top, slope) == pytest.approx(top / slope)


def test_the_parameter_is_reachable_from_the_override_hook():
    """The enabling change. Without the arm there is no path to the field, and
    this whole page would be unmeasurable."""
    assert f'"{M.BOOST_PARAM}" => params.{M.BOOST_PARAM} = val,' in PARAMS_RS
    assert f"pub {M.BOOST_PARAM}: f64," in PARAMS_RS


def test_adding_the_arm_cannot_move_an_undriven_run():
    """Byte-identity as a PROPERTY: the match consults only keys it is handed,
    so the Rust test asserting that must exist rather than the claim resting on
    a matrix re-run alone."""
    assert "acsl4_status_boost_is_drivable_and_adding_it_cannot_move_an_undriven_run" \
        in PARAMS_RS
    assert "the new arm touched a field other than its own" in PARAMS_RS


# ------------------------------------------------- both requirements, not one

def test_a_window_needs_responsiveness_and_admissibility_together(raw):
    """Make the untreated arm die everywhere and no window may survive, even
    though the treated arm still responds."""
    d = _doctor(raw)
    for rows in d["engine"].values():
        for r in rows:
            r["untreated_kill"] = 1.0
    after = M.assemble(d)
    assert after["phenotypes_usable"] == []
    for w in after["windows"].values():
        assert not w["window_exists"]
    assert after["verdict"] == "slope_not_bounded_by_this_route"


def test_a_flat_treated_arm_gives_no_window(raw):
    d = _doctor(raw)
    for rows in d["engine"].values():
        for r in rows:
            r["treated_kill"] = 0.0
    after = M.assemble(d)
    assert after["phenotypes_usable"] == []
    assert after["verdict"] == "slope_not_bounded_by_this_route"


def test_the_effect_floor_binds(raw):
    """A floor above every measured kill must remove every window."""
    original = M.EFFECT_FLOOR
    try:
        M.EFFECT_FLOOR = 1.01
        after = M.assemble(raw)
    finally:
        M.EFFECT_FLOOR = original
    assert after["phenotypes_usable"] == []


def test_the_effect_floor_clears_the_sampling_resolution(committed):
    """A floor below the noise would call a rounding an effect.

    Derived from the intervals the measurements carry, not asserted: `n` is a
    declared constant and lowering it must break this rather than silently
    making the floor meaningless.
    """
    res = committed["sampling_resolution"]
    assert res["usable"], "no interval was carried near the floor"
    assert res["n_points"] > 0
    assert res["floor_clears_resolution"] is True
    assert res["effect_floor"] > res["worst_half_width"], (
        f"the effect floor {res['effect_floor']} is inside the sampling "
        f"half-width {res['worst_half_width']}")


def test_the_measurement_carries_its_interval_at_the_source():
    """`measure_engine` needs the compiled extension, so it does not run in CI
    and an artifact-reading guard cannot see it change.

    A sweep replaced the carried interval with `None` and nothing failed: the
    committed artifact still held the old values, and `assemble` never calls the
    measurement. Pinned in SOURCE, which is the only place this is visible --
    the same reason the Rust property test is asserted by text above.
    """
    src = SCRIPT.read_text()
    tree = ast.parse(src)
    fn = next(n for n in tree.body
              if isinstance(n, ast.FunctionDef) and n.name == "measure_engine")
    body = ast.get_source_segment(src, fn)
    for key in ("treated_ci_low", "treated_ci_high"):
        assert f'"{key}": treated.get(' in body, (
            f"measure_engine no longer carries {key}, so the effect floor "
            "cannot be checked against the sampling resolution")


def test_the_resolution_verdict_is_derived(raw):
    """Widen the intervals past the floor and the flag must go false."""
    d = _doctor(raw)
    for rows in d["engine"].values():
        for r in rows:
            r["treated_ci_low"] = 0.0
            r["treated_ci_high"] = 1.0
    after = M.assemble(d)["sampling_resolution"]
    assert after["usable"]
    assert after["floor_clears_resolution"] is False


def test_the_page_states_the_z_its_bound_is_anchored_at(committed):
    """A bound quoted without its anchor is meaningless: at z = 0 the boost is
    zero for every slope and no bound exists at all."""
    md = M.render(committed)
    assert "anchored at one z" in md
    assert "no bound exists" in md
    top = committed["bridge"]["z_high"]
    assert f"z {top:+.0f}" in md


def test_the_effect_floor_is_a_declared_fraction():
    assert 0.0 < M.EFFECT_FLOOR < M.BASELINE_MAX * 10
    assert M.EFFECT_FLOOR > 0.0, "a floor of zero calls any rounding an effect"


# --------------------------------------------- the grid must contain the crossing

def test_an_unbracketed_crossing_is_reported_as_such(committed):
    """The defect this column exists for: a bound read off the last grid point
    is where the grid stopped, not where the model changed. At least one real
    phenotype must show it, or the column is decoration."""
    flags = {p: w["admissibility_crossing_bracketed"]
             for p, w in committed["windows"].items()}
    assert any(flags.values()), "no crossing is bracketed anywhere"
    assert not all(flags.values()), (
        "every crossing is bracketed, so this guard cannot tell whether the "
        "column is derived -- check the grid still extends past one phenotype")


def test_a_grid_that_stops_early_is_not_bracketed(raw):
    """Truncate the grid below every crossing and the flag must go false."""
    d = _doctor(raw)
    for p, rows in d["engine"].items():
        d["engine"][p] = [r for r in rows if r["untreated_kill"] < M.BASELINE_MAX]
    after = M.assemble(d)
    for p, w in after["windows"].items():
        assert not w["admissibility_crossing_bracketed"], p


def test_the_grid_extends_past_the_bridge(committed):
    """So "unresponsive in range" is distinguishable from "unresponsive anywhere"."""
    top = committed["bridge"]["boost_at_z_high"]
    boosts = [r["boost"] for r in next(iter(committed["engine"].values()))]
    assert max(boosts) > top * 2


# ---------------------------------------------------------- the bound's shape

def test_the_bound_inverts_the_boost_window_and_flips_its_order(raw):
    """Driven through `assemble`, not read off the artifact.

    A sweep swapped the two division operands and only the freshness guard
    fired, because this guard was comparing the committed file to itself.
    """
    for source in (json.loads(M.OUT_JSON.read_text()), M.assemble(raw)):
        top = source["bridge"]["z_high"]
        checked = 0
        for p, s in source["slope_bounds"].items():
            w = source["windows"][p]
            if not w["usable_for_calibration"]:
                assert s["slope_low"] is None and s["slope_high"] is None, p
                continue
            assert s["slope_low"] == pytest.approx(top / w["boost_admissible_to"]), p
            assert s["slope_high"] == pytest.approx(top / w["boost_responsive_from"]), p
            assert s["slope_low"] < s["slope_high"], (
                f"{p}: the inversion did not flip the order")
            checked += 1
        assert checked > 0, "no usable window, so the inversion is untested"


def test_a_phenotype_without_a_window_gets_no_bound_rather_than_a_wide_one(raw):
    """Driven, for the same reason: dropping the `usable` gate, or widening it
    to ignore the bridge range, left this guard green when it read the file."""
    for source in (json.loads(M.OUT_JSON.read_text()), M.assemble(raw)):
        unusable = [p for p, w in source["windows"].items()
                    if not w["usable_for_calibration"]]
        assert unusable, "every phenotype is usable, so this guard is inert"
        for p in unusable:
            assert source["slope_bounds"][p]["slope_low"] is None, p
            assert source["slope_bounds"][p]["basis"] == "no usable window", p
    # And `usable` must REQUIRE the bridge range, not merely a window: a
    # phenotype responsive only far outside it is not usable.
    d = _doctor(raw)
    top = json.loads(M.OUT_JSON.read_text())["bridge"]["boost_at_z_high"]
    for rows in d["engine"].values():
        for r in rows:
            r["treated_kill"] = 1.0 if r["boost"] > top * 2 else 0.0
            r["untreated_kill"] = 0.0
    after = M.assemble(d)
    for p, w in after["windows"].items():
        assert w["window_exists"], p
        assert not w["usable_for_calibration"], (
            f"{p}: responsive only outside the bridge is being called usable")
        assert after["slope_bounds"][p]["slope_low"] is None, p


def test_the_shipped_slope_is_reported_inside_or_outside_rather_than_assumed(raw):
    """Move the window so the shipped slope falls outside it and the flag must
    follow."""
    d = _doctor(raw)
    for rows in d["engine"].values():
        for r in rows:
            # Responsive only at the very top, so the slope bound sits low.
            r["treated_kill"] = 1.0 if r["boost"] >= max(
                x["boost"] for x in rows) else 0.0
    after = M.assemble(d)
    for p, s in after["slope_bounds"].items():
        if s["slope_low"] is None:
            continue
        assert s["shipped_inside"] == (s["slope_low"] <= M.SHIPPED_SLOPE
                                       <= s["slope_high"]), p


def test_the_verdict_requires_the_shipped_slope_to_be_bounded(raw, committed):
    assert committed["verdict"] == "slope_bounded_not_fitted"
    d = _doctor(raw)
    for rows in d["engine"].values():
        for r in rows:
            r["untreated_kill"] = 1.0
    assert M.assemble(d)["verdict"] == "slope_not_bounded_by_this_route"


# ------------------------------------------------------------ the measured leg

def test_the_measured_direction_is_derived_not_asserted(raw, committed):
    for c, m in committed["measured_shape"].items():
        if not m.get("usable"):
            continue
        assert m["direction_matches_engine"] == (m["absolute_span"] > 0), c
    d = _doctor(raw)
    for leg in d["cell_lines"].values():
        for b in leg.get("bins", []):
            b["mean_sensitivity"] = -b["mean_sensitivity"]
    after = M.assemble(d)
    assert after["compounds_agreeing_in_direction"] == []


def test_the_measured_shape_reports_monotonicity_as_a_count(committed):
    for c, m in committed["measured_shape"].items():
        if not m.get("usable"):
            continue
        assert 0 <= m["steps_rising"] <= m["steps_total"], c
        assert m["monotone"] == (m["steps_rising"] == m["steps_total"]), c


def test_a_thin_bin_is_dropped_rather_than_averaged(committed, tmp_path,
                                                   monkeypatch):
    """The filter is INERT on the real join -- every bin clears it -- so a
    sweep that removed it entirely survived both this guard and the freshness
    check. Exercised against a panel built to contain a thin bin.
    """
    for c, bins in committed["measured_bins"].items():
        for b in bins:
            assert b["n_lines"] >= M.MIN_LINES_PER_BIN, (c, b["z"])
    assert 1 < M.MIN_LINES_PER_BIN, "a bin of one line is not a mean"
    compound = M.COMPOUNDS[0]
    panel = tmp_path / "panel.csv"
    rows = [f"ModelID,OncotreeLineage,expr_{M.BRIDGE_GENE},viability_auc_{compound}"]
    # A dense bin at the centre and a single far-out line, which must be dropped.
    for i in range(M.MIN_LINES_PER_BIN * 3):
        rows.append(f"ACH-{i:06d},L,5.0,0.5")
    for i in range(M.MIN_LINES_PER_BIN * 3):
        rows.append(f"ACH-1{i:05d},L,6.0,0.4")
    rows.append("ACH-999999,L,99.0,0.1")
    panel.write_text("\n".join(rows) + "\n")
    monkeypatch.setattr(M, "PANEL", panel)
    legs = M.measure_cell_lines()
    kept = legs[compound]["bins"]
    assert kept, "the doctored panel produced no bins at all"
    for b in kept:
        assert b["n_lines"] >= M.MIN_LINES_PER_BIN, b
    assert sum(b["n_lines"] for b in kept) < M.MIN_LINES_PER_BIN * 6 + 1, (
        "the single far-out line was kept, so the thin-bin filter is gone")


def test_a_blank_field_is_absent_not_zero():
    assert M._f("") is None
    assert M._f(None) is None
    assert M._f("0") == 0.0


def test_the_measured_leg_reads_the_committed_join(committed):
    assert M.PANEL.exists()
    for c in M.COMPOUNDS:
        assert committed["cell_lines"][c]["n_lines"] > 100, c


# ------------------------------------------------------------ the refusals

def test_the_page_fits_no_slope(committed):
    """The whole point of the restraint: a reader must not find a fitted value."""
    md = M.render(committed)
    assert "bounded, not fitted" in md
    assert "cannot distinguish a wrong functional form" in md
    for forbidden in ("best fit", "fitted slope of", "calibrated slope ="):
        assert forbidden not in md, forbidden


def test_the_page_states_the_identifiability_count_does_not_move(committed):
    md = M.render(committed)
    assert "data-conditioned-headline count" in md
    report = json.loads((REPO / "analysis" / "identifiability-report.json").read_text())
    assert report.get("data_constrained_in_production") in (0, None), (
        "the identifiability report now claims a data-conditioned headline; this "
        "page's refusal has to be revisited rather than left standing")


def test_the_page_names_the_averaging_confound(committed):
    md = M.render(committed)
    assert "single-cell death rate" in md
    assert "population viability averaged over a dose range" in md


def test_the_page_says_how_many_phenotypes_carry_it(committed):
    md = M.render(committed)
    n = len(committed["phenotypes_usable"])
    assert f"{n} of {len(committed['windows'])} phenotypes" in md


# ------------------------------------------------------------- housekeeping

def test_ordered_tables_survive_the_sorted_dump(committed):
    """`sort_keys` reorders dicts, so a boost grid presented in order must be a
    LIST. This repo has shipped alphabetised tables presented as rankings."""
    for p, rows in committed["engine"].items():
        assert isinstance(rows, list), p
        assert [r["boost"] for r in rows] == sorted(r["boost"] for r in rows), p
    for c, bins in committed["measured_bins"].items():
        assert isinstance(bins, list), c
        assert [b["z"] for b in bins] == sorted(b["z"] for b in bins), c


def test_the_json_dump_options_are_stable():
    assert "json.dumps(data, indent=2, sort_keys=True)" in SCRIPT.read_text()


def test_the_cell_count_and_seed_are_declared(committed):
    cons = committed["constraints"]
    assert cons["n_cells"] >= 1000, "too few cells to read a percent-level kill"
    assert cons["seed"] == M.SEED
    assert cons["baseline_viability_max"] == M.BASELINE_MAX
    assert cons["effect_floor"] == M.EFFECT_FLOOR


def test_the_module_docstring_quotes_no_measurement():
    doc = re.sub(r"#\d+", "", M.__doc__)
    standalone = re.findall(r"(?<![A-Za-z0-9_])\d[\d.,%]*", doc)
    assert not standalone, f"the docstring carries numbers: {standalone}"


def test_the_generator_needs_no_network_and_no_census():
    src = SCRIPT.read_text()
    for forbidden in ("urllib", "requests", "corpus/atlas", "subprocess"):
        assert forbidden not in src, forbidden
