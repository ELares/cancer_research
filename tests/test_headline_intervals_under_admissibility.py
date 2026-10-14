"""Guards for the conditioned spatial and combination headline intervals (#938).

The load-bearing property is the same one `#936`'s guards protect and is worth
restating because it is easy to lose: the FULL-PRIOR half must reproduce the
three committed reports. Without that this page corrects a re-run rather than
what is published, and every shift it reports is measured against the wrong
baseline.

Three further properties are specific to this page and each is driven rather
than asserted:

  * the three families' draws are NESTED (the sim-tme family's are a prefix of
    the others'), which is what lets one admissibility pass cover all three;
  * the Bliss score is a RATIO, so a rise in it has two possible sources and the
    page must be able to say which;
  * the penetration ordering is robust BY CONSTRUCTION (#926), so it cannot
    break here -- registered in advance, and guarded by driving it to a value it
    could never take and requiring the VERDICT not to move. A page that let that
    fraction into its verdict would be reporting an entailment as evidence.
"""

from __future__ import annotations

import ast
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
SCRIPTS = REPO / "scripts"


def _load(name):
    sys.path.insert(0, str(SCRIPTS))
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


M = _load("headline_intervals_under_admissibility")
SCRIPT = SCRIPTS / "headline_intervals_under_admissibility.py"


def _src():
    return SCRIPT.read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def raw():
    return M.load()


@pytest.fixture(scope="module")
def built(raw):
    return M.assemble(raw)


def _doctor(raw):
    """A shallow per-row copy. `copy.deepcopy` of the whole table per doctored
    test was enough to get pytest killed beside a long-running crawl once, and
    the rows hold only strings."""
    return {"rows": [dict(r) for r in raw["rows"]],
            "rust": raw["rust"],
            "sibling": json.loads(json.dumps(raw["sibling"])),
            "identifiability": dict(raw["identifiability"]),
            "subadditive_citations": dict(raw["subadditive_citations"]),
            "published": dict(raw["published"])}


def _obs(built, column):
    return next(o for o in built["observables"] if o["column"] == column)


# --------------------------------------------------------------------------
# the committed artifacts are what this generator produces, now
# --------------------------------------------------------------------------
def test_committed_page_is_what_the_generator_renders(built):
    """The strong round trip, taken directly rather than through the stored JSON.

    The freshness gate separately checks that `assemble` is idempotent and that
    the page does not depend on the stored dicts' key order; this checks the
    whole chain from the committed INPUTS, which is what a reader of the page
    actually relies on.
    """
    assert M.OUT_MD.read_text() == M.render(built), (
        "analysis/headline-intervals-under-admissibility.md is not what the "
        "generator renders from the committed inputs; regenerate it")


def test_committed_json_is_what_the_generator_assembles(built):
    stored = json.loads(M.OUT_JSON.read_text())
    assert stored == json.loads(json.dumps(built, sort_keys=True)), (
        "the committed JSON disagrees with a fresh assemble of the same inputs")


# --------------------------------------------------------------------------
# ONE source for the published design -- the premise everything else rests on
# --------------------------------------------------------------------------
def test_the_design_is_imported_and_not_restated():
    """A reproduction whose draw count, seed or ranges are retyped is a
    reproduction of something else. Checked on the IMPORTS, not on values."""
    tree = ast.parse(_src())
    imported = {alias.name: node.module
                for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)
                for alias in node.names}
    for name, want in (("sample_prior", "headline_uncertainty"),
                       ("DEFAULT_SAMPLES", "headline_uncertainty"),
                       ("DEFAULT_TME_SAMPLES", "headline_uncertainty"),
                       ("DEFAULT_PEN_SAMPLES", "headline_uncertainty"),
                       ("SEED", "headline_uncertainty"),
                       ("BASELINE_MAX", "headline_at_fitted"),
                       ("N_CELLS", "uncertainty_intervals"),
                       ("SIM_SEED", "uncertainty_intervals"),
                       ("QUANTILES", "uncertainty_intervals"),
                       ("PARAM_NAMES", "headline_sensitivity"),
                       ("CONTROL_CONDITIONS", "intervals_under_admissibility"),
                       ("_selection", "intervals_under_admissibility")):
        assert imported.get(name) == want, (
            f"{name} must be imported from {want}, not restated here")


def test_the_constraint_is_not_a_literal_in_this_file():
    """Two copies of a threshold are free to disagree."""
    assert f"{M.BASELINE_MAX}" not in _src().replace("BASELINE_MAX", ""), (
        "the baseline-viability ceiling appears as a literal in this generator; "
        "it must only ever arrive through the import")


def test_the_three_families_draws_are_nested():
    """One admissibility pass covers all three families ONLY because the
    smaller samples are prefixes of the larger. If that stops holding, the
    sim-tme partition silently starts using another family's draws."""
    big = M.sample_prior(M.N_ROWS, seed=M.SEED)
    for n in sorted(set(M.FAMILY_SAMPLES.values())):
        small = M.sample_prior(n, seed=M.SEED)
        assert np.array_equal(big[:n], small), (
            f"the {n}-draw sample is not a prefix of the {M.N_ROWS}-draw one, so "
            "the admissibility mask does not line up with that family's draws")


def test_the_prior_is_the_same_box_as_the_published_single_cell_design():
    """A draw's admissibility is only well defined if both designs sweep the
    same parameters over the same ranges. Measured, not assumed."""
    import headline_sensitivity as H
    import uncertainty_intervals as U
    assert list(H.PARAM_NAMES) == [p[0] for p in U.PARAMS]
    assert np.array_equal(H.LOWS, np.array([p[1] for p in U.PARAMS]))
    assert np.array_equal(H.HIGHS, np.array([p[2] for p in U.PARAMS]))


def test_the_confidence_label_matches_the_quantiles_actually_computed():
    """The label is derived from `QUANTILES` while the numbers come from
    `_pctiles`. If those ever disagree the page would advertise one interval and
    print another."""
    low, mid, high = M.QUANTILES
    assert (low, mid, high) == (0.025, 0.5, 0.975)
    arr = np.arange(0.0, 1.0, 0.001)
    p = M._pctiles(arr)
    assert math.isclose(p["p2_5"], float(np.percentile(arr, 100 * low)))
    assert math.isclose(p["p97_5"], float(np.percentile(arr, 100 * high)))
    assert M._confidence_label() == "95%"


# --------------------------------------------------------------------------
# the reproduction gate
# --------------------------------------------------------------------------
def test_every_observable_reproduces_what_is_published(built):
    bad = [o["label"] for o in built["observables"] if not o["reproduces"]]
    assert not bad, f"these no longer reproduce their published figures: {bad}"
    assert built["all_reproduce"] is True
    # and the gate is not vacuous: each observable must have CHECKED something
    for o in built["observables"]:
        assert len(o["reproduction"]) >= 5, (
            f"{o['label']} reproduces only {len(o['reproduction'])} figures; a "
            "reproduction gate checking one number is not one")


def test_the_tolerance_is_derived_from_the_published_precision():
    """Both directions. Too loose accepts a real discrepancy while reporting a
    reproduction; too tight refuses a table the page does reproduce. #936's
    first attempt was 22x looser than the largest real deviation."""
    assert M._tolerance(M.PUBLISHED_DECIMALS) == pytest.approx(
        M.ROUNDING_SLACK * 0.5 * 10 ** (-M.PUBLISHED_DECIMALS))
    floor = 0.5 * 10 ** (-M.PUBLISHED_DECIMALS)
    tol = M._tolerance(M.PUBLISHED_DECIMALS)
    assert tol >= floor, "the tolerance is below the rounding floor it must absorb"
    assert tol <= 10 * floor, "the tolerance has drifted far above the rounding it absorbs"


def test_the_largest_real_deviation_is_well_inside_the_tolerance(built):
    """A tolerance is only informative beside the deviation it has to absorb.

    No lower bound is asserted on that deviation: these figures can and do match
    the published ones exactly, so requiring a non-zero gap would be demanding
    an error. (An earlier line here read `assert worst > 0.0 or True`, which is
    `True` -- the no-op shape this repo has a gate for.)
    """
    deviations = [c["deviation"]
                  for o in built["observables"] for c in o["reproduction"].values()]
    assert len(deviations) >= 5 * len(built["observables"]), (
        "fewer deviations than published figures, so some figure stopped being "
        "compared at all")
    worst = max(deviations)
    assert worst <= built["tolerance"], (
        f"largest deviation {worst} exceeds the tolerance "
        f"{built['tolerance']}, so something stopped reproducing")


def test_a_disagreeing_published_figure_withholds_the_verdict(raw):
    """Drives the refusal: doctor one published median and the page must refuse
    rather than compare against a re-run."""
    d = _doctor(raw)
    d["published"]["bliss"] = d["published"]["bliss"].replace(
        "| prior-predictive median | 1.348 |",
        "| prior-predictive median | 9.999 |")
    assert d["published"]["bliss"] != raw["published"]["bliss"], "the doctoring missed"
    out = M.assemble(d)
    o = _obs(out, M.BLISS_COLUMN)
    assert o["reproduces"] is False
    assert o["reportable"] is False
    assert out["verdict"] is False
    assert out["all_reproduce"] is False


def test_a_moved_measurement_also_fails_the_gate(raw):
    """The other side of the same gate: the published file is fine and the
    MEASUREMENT moved."""
    d = _doctor(raw)
    for r in d["rows"]:
        if r["pen_cns_bbb"] not in ("", M.FAILED):
            r["pen_cns_bbb"] = f"{float(r['pen_cns_bbb']) * 0.5 + 0.2:.6f}"
    out = M.assemble(d)
    assert _obs(out, "pen_cns_bbb")["reproduces"] is False
    assert out["verdict"] is False


# --------------------------------------------------------------------------
# admissibility
# --------------------------------------------------------------------------
def test_admissibility_is_measured_on_the_committed_controls(built):
    fam = built["families"]
    for name, f in fam.items():
        assert 0 < f["n_admissible"] <= f["draws"]
        assert f["fraction"] == pytest.approx(f["n_admissible"] / f["draws"])
    assert set(built["constraint"]["controls"]) == set(M.CONTROL_COLUMNS)
    assert built["constraint"]["baseline_max"] == M.BASELINE_MAX


def test_making_every_draw_admissible_removes_every_shift(raw):
    d = _doctor(raw)
    for r in d["rows"]:
        for c in M.CONTROL_COLUMNS:
            r[c] = "0.000000"
    out = M.assemble(d)
    for f in out["families"].values():
        assert f["fraction"] == 1.0
    for o in out["observables"]:
        for s in o["shifts"].values():
            assert s["delta"] == pytest.approx(0.0), (
                "with nothing removed the admissible interval must equal the full one")
        assert not any(s["material"] for s in o["shifts"].values())
    assert out["n_material"] == 0


def test_making_every_draw_inadmissible_withholds_every_interval(raw):
    d = _doctor(raw)
    for r in d["rows"]:
        r[M.CONTROL_COLUMNS[0]] = "1.000000"
    out = M.assemble(d)
    for f in out["families"].values():
        assert f["fraction"] == 0.0
        assert f["usable"] is False
    for o in out["observables"]:
        assert o["admissible_only"] is None
        assert o["reportable"] is False
    assert out["verdict"] is False


def test_the_constraint_is_strict_at_its_own_boundary(raw):
    """A draw whose untreated death is exactly the stated ceiling does NOT
    satisfy a constraint of "less than" that ceiling."""
    d = _doctor(raw)
    for r in d["rows"]:
        for c in M.CONTROL_COLUMNS:
            r[c] = f"{M.BASELINE_MAX:.6f}"
    out = M.assemble(d)
    assert all(f["n_admissible"] == 0 for f in out["families"].values()), (
        "a draw sitting exactly on the ceiling is being counted as admissible")


def test_the_minimum_admissible_floor_actually_binds(raw):
    """An inert floor is a floor that protects nothing -- the shape found inert
    twice in #930 and #932 and once in #934."""
    keep = M.MIN_ADMISSIBLE_DRAWS - 1
    d = _doctor(raw)
    for i, r in enumerate(d["rows"]):
        for c in M.CONTROL_COLUMNS:
            r[c] = "0.000000" if i < keep else "1.000000"
    out = M.assemble(d)
    assert out["families"]["bliss"]["n_admissible"] == keep
    assert out["families"]["bliss"]["usable"] is False
    o = _obs(out, M.BLISS_COLUMN)
    assert o["interval_quotable"] is False, "the floor did not withhold the tails"
    assert o["shifts"]["low"]["quotable"] is False
    assert o["shifts"]["high"]["quotable"] is False
    assert o["shifts"]["median"]["quotable"] is True, (
        "a median needs far fewer draws than a tail; withholding it too would "
        "throw away the one figure that survives")


# --------------------------------------------------------------------------
# the material-shift bar
# --------------------------------------------------------------------------
def test_the_material_bar_is_relative_to_the_published_width(raw):
    """The decisive non-vacuity check: the SAME delta must change its verdict
    when the published WIDTH changes. An absolute threshold would not move."""
    base = M.assemble(raw)
    o = _obs(base, "pen_well_vascularized")
    delta = abs(o["shifts"]["median"]["delta"])
    assert delta > 0, "this fixture needs an observable whose median actually moves"

    d = _doctor(raw)
    wide = delta / M.MATERIAL_FRACTION * 2.0            # bar above the delta
    d["published"]["penetration"] = d["published"]["penetration"].replace(
        "| [0.002, 0.926] |", f"| [0.000, {wide:.3f}] |")
    out = M.assemble(d)
    assert _obs(out, "pen_well_vascularized")["published_width"] == pytest.approx(wide)
    assert _obs(out, "pen_well_vascularized")["shifts"]["median"]["material"] is False, (
        "the material bar is not scaled by the published interval width")
    assert o["shifts"]["median"]["material"] is True


def test_material_shifts_exist_on_the_committed_data(built):
    assert built["n_material"] >= 1, (
        "no headline carries a material shift, so this page would be reporting "
        "a correction nobody needs -- check the measurement, not the bar")


# --------------------------------------------------------------------------
# the direction claim, which must be able to come back either way
# --------------------------------------------------------------------------
def test_the_two_directions_are_derived_and_not_asserted(raw):
    """The page's sharpest sentence is that one mechanism moves some medians
    down and others up. Doctor the data so they all move ONE way and the
    sentence must follow."""
    base = M.render(M.assemble(raw))
    assert "medians DOWN and" in base

    # Fix an admissible subset, then halve every observable INSIDE it, so each
    # admissible median must land below its full-prior one. Selecting draws on
    # one observable does not control the others -- an earlier version of this
    # fixture selected on the Bliss score and five of the six medians went UP.
    d = _doctor(raw)
    columns = [c for c, _, _, _ in M.OBSERVABLES]
    for i, r in enumerate(d["rows"]):
        keep = i % 2 == 0
        for c in M.CONTROL_COLUMNS:
            r[c] = "0.000000" if keep else "1.000000"
        if keep:
            for c in columns:
                if r[c] not in ("", M.FAILED):
                    r[c] = f"{float(r[c]) * 0.5:.6f}"
    out = M.assemble(d)
    ups = [o["label"] for o in out["observables"]
           if o["shifts"].get("median", {}).get("delta", 0) > 0]
    assert not ups, f"a doctored one-directional world still reports rises: {ups}"
    assert "0 UP" in M.render(out)


def test_the_point_crossing_claim_is_derived(raw):
    """Where the manuscript's own point estimate lands is a derived statement,
    so a world where nothing crosses must not print the crossing section."""
    built = M.assemble(raw)
    crossed = [o["label"] for o in built["observables"] if o["point"]["crossed"]]
    d = _doctor(raw)
    for r in d["rows"]:                      # nothing removed => nothing crosses
        for c in M.CONTROL_COLUMNS:
            r[c] = "0.000000"
    out = M.assemble(d)
    assert not any(o["point"]["crossed"] for o in out["observables"])
    assert "crosses from one side" not in M.render(out)
    if crossed:
        assert "crosses from one side" in M.render(built)


# --------------------------------------------------------------------------
# the Bliss ratio decomposition
# --------------------------------------------------------------------------
def test_the_stored_bliss_parts_satisfy_the_engines_own_identity(built):
    """A column mix-up is otherwise invisible: the score would still look like a
    plausible ratio."""
    assert built["bliss_decomposition"]["identity_max_error"] < 1e-6, (
        "the stored Bliss parts do not reproduce the engine's own "
        "rate_a + rate_b - rate_a*rate_b, so the columns are not what they say")


def test_a_mixed_up_part_column_breaks_the_identity(raw):
    d = _doctor(raw)
    for r in d["rows"]:
        r["bliss_rate_rsl3"], r["bliss_rate_combo"] = (
            r["bliss_rate_combo"], r["bliss_rate_rsl3"])
    out = M.assemble(d)
    assert out["bliss_decomposition"]["identity_max_error"] > 1e-6, (
        "the identity check cannot see the two columns swapped")


def test_which_half_of_the_ratio_moves_is_measured(raw):
    """`ratio_follows_denominator` must be able to come back either way, or the
    page's mechanism sentence is decoration."""
    base = M.assemble(raw)["bliss_decomposition"]
    assert base["numerator_ratio"] is not None and base["denominator_ratio"] is not None
    # Collapse the numerator in the ADMISSIBLE draws only. Scaling every row
    # leaves `admissible / full` unchanged, which is how the first version of
    # this fixture tested nothing at all.
    d = _doctor(raw)
    for r in d["rows"]:
        if max(float(r[c]) for c in M.CONTROL_COLUMNS) < M.BASELINE_MAX:
            r["bliss_rate_combo"] = f"{float(r['bliss_rate_combo']) * 0.01:.6f}"
    out = M.assemble(d)["bliss_decomposition"]
    assert out["ratio_follows_denominator"] != base["ratio_follows_denominator"], (
        "collapsing the numerator did not change which half the page blames")


def test_the_undefined_bliss_floor_is_parsed_from_the_engine(raw):
    """Repoint the parser at a doctored source and the reported floor must move.
    A constant here would be a second copy of an engine property."""
    assert M.undefined_bliss_floor("if bliss > 0.001 {") == 0.001
    assert M.undefined_bliss_floor("if bliss > 0.25 {") == 0.25
    d = _doctor(raw)
    d["rust"] = d["rust"].replace("if bliss > 0.001", "if bliss > 0.5")
    out = M.assemble(d)
    assert out["bliss_decomposition"]["floor"] == 0.5
    assert out["bliss_decomposition"]["denominator_clears_floor"] is False, (
        "raising the engine's own refusal threshold above the smallest "
        "admissible denominator must be visible on the page")
    with pytest.raises(SystemExit):
        M.undefined_bliss_floor("no floor here")


def test_the_denominator_clears_the_floor_on_the_committed_data(built):
    dec = built["bliss_decomposition"]
    assert dec["denominator_clears_floor"] is True, (
        "the conditioned Bliss interval is being computed beside the corner "
        "where the engine refuses the ratio")
    assert dec["admissible_min_denominator"] > dec["floor"]


def test_the_single_agent_rates_are_resolved_by_drug_name(raw):
    """The shared reader matches the two drugs as a SET, so reading `rate_a`
    positionally would swap the agents on any run that emits them the other way
    round. Driven on a swapped row, without the binary."""
    a = {"drug_a": "RSL3", "drug_b": "FSP1i", "rate_a": "0.4", "rate_b": "0.037",
         "rate_combo": "0.841", "bliss_prediction": "0.4222", "synergy_score": "1.99"}
    b = dict(a, drug_a="FSP1i", drug_b="RSL3", rate_a="0.037", rate_b="0.4")
    assert M._bliss_from_row(a) == M._bliss_from_row(b)
    assert M._bliss_from_row(b)["bliss_rate_rsl3"] == 0.4


# --------------------------------------------------------------------------
# the ordering, registered in advance as NOT evidence
# --------------------------------------------------------------------------
def test_the_ordering_reproduces_the_published_figure(built):
    o = built["ordering"]
    assert o["published"] is not None, "the published ordering line stopped parsing"
    assert o["reproduces_published"] is True
    assert o["structural"] is True


def test_the_ordering_fraction_cannot_reach_the_verdict():
    """#926 established this ordering is robust BY CONSTRUCTION, so letting it
    into the verdict would report an entailment as evidence.

    Checked STRUCTURALLY, and that is not a shortcut: the ordering is a
    within-row property of three columns, so any doctoring that breaks it also
    changes those columns' own distributions and therefore their reproduction
    against the published report. A behavioural test would move the verdict for
    that reason instead -- which is exactly what the first version of this guard
    did, and it read as the ordering feeding the verdict when it was not.
    """
    tree = ast.parse(_src())
    fn = next(n for n in ast.walk(tree)
              if isinstance(n, ast.FunctionDef) and n.name == "assemble")
    verdicts = [n for n in ast.walk(fn)
                if isinstance(n, ast.Dict)
                for k, v in zip(n.keys, n.values)
                if isinstance(k, ast.Constant) and k.value == "verdict"]
    assert verdicts, "the assembled verdict is no longer a dict entry named 'verdict'"
    for node in verdicts:
        expr = next(v for k, v in zip(node.keys, node.values)
                    if isinstance(k, ast.Constant) and k.value == "verdict")
        names = {n.id for n in ast.walk(expr) if isinstance(n, ast.Name)}
        assert "ordering" not in names, (
            "the ordering is reachable from the verdict expression, which turns "
            "an entailment into evidence")


def test_a_broken_ordering_is_still_reported(raw):
    """Structural above, so this is the other half: the page must still MEASURE
    the ordering and say when it breaks, rather than not looking."""
    d = _doctor(raw)
    for r in d["rows"]:                       # invert the ordering everywhere
        r["pen_well_vascularized"], r["pen_cns_bbb"] = (
            r["pen_cns_bbb"], r["pen_well_vascularized"])
    out = M.assemble(d)
    assert out["ordering"]["admissible_only"]["fraction"] < 1.0
    assert out["ordering"]["reproduces_published"] is False
    assert "does NOT reproduce" in M.render(out)


# --------------------------------------------------------------------------
# the circularity falsifier
# --------------------------------------------------------------------------
def test_the_drug_parameters_are_untouched_by_the_truncation(built):
    sel = built["selection"]
    assert sel["drug_parameters"], "the falsifier is measuring nothing"
    assert sel["drug_parameters_untouched"] is True, (
        "conditioning selected a drug parameter, so the kill reductions on this "
        "page are partly an artifact of keeping weaker-drug draws")
    assert sel["max_drug_shift"] < M.STRONG_SELECTION


def test_a_selected_drug_parameter_withholds_the_verdict(raw):
    """The case that does not occur naturally, so it has to be driven: admit
    only draws whose drug parameter is low and the falsifier must fire.

    The fixture admits exactly the SIBLING'S measured share, so the
    cross-artifact agreement check still passes and the falsifier is the only
    thing left that can withhold the verdict. An earlier version admitted a
    quarter of the draws, which tripped the agreement check instead -- and the
    guard then passed for that reason, which is how a mutation removing the
    falsifier from the verdict survived a sweep.
    """
    d = _doctor(raw)
    # Read the drug parameters off the DERIVED field rather than re-importing
    # the constant, so this guard breaks if the page stops naming them.
    name = M.assemble(raw)["selection"]["drug_parameters"][0]
    target = round(raw["sibling"]["admissibility"]["fraction"] * len(d["rows"]))
    order = sorted(range(len(d["rows"])), key=lambda i: float(d["rows"][i][name]))
    keep = set(order[:target])
    for i, r in enumerate(d["rows"]):
        for c in M.CONTROL_COLUMNS:
            r[c] = "0.000000" if i in keep else "1.000000"
    out = M.assemble(d)
    assert out["agreement"]["agrees"] is True, (
        "this fixture must not trip the agreement check, or it cannot isolate "
        "the falsifier")
    assert out["all_reproduce"] is True, (
        "this fixture must not disturb the full-prior half")
    assert out["selection"]["per_parameter"][name]["strongly_selected"] is True
    assert out["selection"]["drug_parameters_untouched"] is False
    assert out["verdict"] is False, (
        "a strongly selected drug parameter must withhold the verdict on its "
        "own: the kill reductions would otherwise be partly an artifact of "
        "keeping weaker-drug draws")


# --------------------------------------------------------------------------
# cross-artifact agreement about the prior
# --------------------------------------------------------------------------
def test_the_two_pages_agree_about_the_prior(built):
    ag = built["agreement"]
    assert ag["agrees"] is True, (
        "this page and #936 no longer describe the same prior, which means the "
        "two definitions of admissible have drifted apart")
    assert ag["max_abs_z"] < ag["bar"]
    assert ag["sibling_draws"] > max(M.FAMILY_SAMPLES.values())


def test_a_disagreeing_sibling_withholds_the_verdict(raw):
    d = _doctor(raw)
    for r in d["rows"]:                      # force an extreme fraction
        for c in M.CONTROL_COLUMNS:
            r[c] = "0.000000"
    out = M.assemble(d)
    assert out["agreement"]["agrees"] is False, (
        "an admissible fraction of 1.0 against the sibling's measured fraction "
        "must not count as agreement")
    assert out["verdict"] is False


# --------------------------------------------------------------------------
# the restraint
# --------------------------------------------------------------------------
def test_the_data_conditioned_headline_count_is_still_zero(built):
    assert built["data_conditioned_headlines"] == 0
    report = json.loads((REPO / "analysis" / "identifiability-report.json").read_text())
    assert report.get("data_constrained_in_production") in (0, None), (
        "the identifiability report now claims a data-conditioned headline, so "
        "this page's refusal has to be re-examined rather than left standing")


def test_the_page_states_the_refusal(built):
    md = M.render(built)
    assert "Conditioning is not inference" in md
    assert "not a posterior" in md
    assert "does NOT establish" in md


# --------------------------------------------------------------------------
# completeness and the three kinds of missing value
# --------------------------------------------------------------------------
def test_an_incomplete_leg_refuses_rather_than_publishing(raw):
    d = _doctor(raw)
    d["rows"][7]["tme_immune"] = ""
    with pytest.raises(SystemExit) as e:
        M.assemble(d)
    assert "tme_immune" in str(e.value)


def test_the_completeness_check_is_not_vacuous(raw):
    M._refuse_incomplete(raw["rows"])        # the committed file must pass it


def test_a_failed_run_is_not_an_observation(raw):
    """`fail`, empty and `nan` are three different things. A failure counted as
    a zero would publish a kill rate of zero as a measurement."""
    d = _doctor(raw)
    d["rows"][3][M.BLISS_COLUMN] = M.FAILED
    out = M.assemble(d)
    o = _obs(out, M.BLISS_COLUMN)
    assert o["n_failed"] == 1
    assert o["n_full"] == _obs(M.assemble(raw), M.BLISS_COLUMN)["n_full"] - 1


def test_an_undefined_observable_is_counted_and_not_averaged(raw):
    d = _doctor(raw)
    d["rows"][4][M.BLISS_COLUMN] = "nan"
    out = M.assemble(d)
    o = _obs(out, M.BLISS_COLUMN)
    assert o["n_undefined"] == 1
    assert math.isfinite(o["full"]["median"]), "one NaN poisoned the percentiles"


def test_the_writer_records_a_non_finite_value_as_undefined():
    assert M._fmt(float("nan")) == "nan"
    assert M._fmt(float("inf")) == "nan"
    assert M._fmt(None) == M.FAILED
    assert M._fmt(0.5) == "0.500000"


# --------------------------------------------------------------------------
# source pins: the measurement half never runs in CI, so only the source can
# say that this page measures the PUBLISHED observables
# --------------------------------------------------------------------------
def test_the_measurement_drives_the_published_observable_extractors():
    """These runners are the published generators' own. A local copy would let
    this page measure a different quantity from the intervals it corrects, and
    no artifact-reading guard could see it."""
    src = _src()
    for name in ("run_sim_tme_observables", "run_tissue_pk_observables",
                 "run_bliss_row"):
        assert name in src, f"{name} is no longer used to take the measurement"
    tree = ast.parse(src)
    measure = next(n for n in ast.walk(tree)
                   if isinstance(n, ast.FunctionDef) and n.name == "measure")
    called = {n.id for n in ast.walk(measure) if isinstance(n, ast.Name)}
    assert {"run_sim_tme_observables", "run_tissue_pk_observables",
            "_bliss_observables"} <= called, (
        "measure() no longer dispatches to the shared observable extractors")


def test_assemble_drives_every_load_bearing_derivation():
    """An artifact-reading guard cannot tell a derivation from an echo, so the
    calls are pinned in source. This is the defect this repo found seven times
    across the three preceding PRs."""
    tree = ast.parse(_src())
    fn = next(n for n in ast.walk(tree)
              if isinstance(n, ast.FunctionDef) and n.name == "assemble")
    called = {n.func.id for n in ast.walk(fn)
              if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)}
    for name in ("_refuse_incomplete", "parse_published", "undefined_bliss_floor",
                 "_selection", "_pctiles"):
        assert name in called, f"assemble() no longer calls {name}"


def test_the_admissibility_leg_uses_the_published_cell_count_and_seed():
    src = _src()
    assert "n=N_CELLS" in src and "seed=SIM_SEED" in src, (
        "the untreated arms must be evaluated at the published design's own "
        "per-draw cell count and simulation seed, or they are not comparable "
        "with the treated arms they gate")


def test_the_legs_cover_every_column_the_page_reads():
    """A column no leg fills could never be complete, and one no leg claims
    would silently escape the completeness refusal."""
    claimed = {c for columns, _, _ in M.LEGS.values() for c in columns}
    derived = set(M.CONTROL_COLUMNS) | {M.BLISS_COLUMN} | set(M.BLISS_PART_COLUMNS) \
        | set(M.TME_COLUMNS) | set(M.PEN_COLUMNS)
    assert claimed == derived
    assert all(c in M.COLUMNS for c in claimed)
    for column, fam, _, _ in M.OBSERVABLES:
        assert column in claimed
        assert fam in M.FAMILY_SAMPLES

# --------------------------------------------------------------------------
# the thinnest family cannot be given a conditioned interval at all
# --------------------------------------------------------------------------
def test_the_floor_is_unreachable_for_the_costliest_family_by_construction():
    """The page claims this is structural, so the structure is asserted: the
    imported floor exceeds the smallest published draw count, so no data could
    ever let that family quote a conditioned tail."""
    assert M.MIN_ADMISSIBLE_DRAWS >= min(M.FAMILY_SAMPLES.values()), (
        "the floor is now reachable for the smallest family without the "
        "constraint being a no-op, so the page's structural paragraph is no "
        "longer true and must be rewritten rather than left standing")
    # And the page must claim exactly that, not more: the floor IS met in the
    # degenerate world where the constraint removes nothing, and an earlier
    # draft of this page said "unreachable whatever the data says", which that
    # world falsifies.
    md = M.render(M.assemble(M.load()))
    assert "unreachable" not in md.lower()
    assert "removed\nNOTHING" in md or "removed NOTHING" in md


def test_the_thin_familys_tails_are_withheld_and_its_median_is_kept(built):
    thin = min(M.FAMILY_SAMPLES, key=lambda k: M.FAMILY_SAMPLES[k])
    thin_obs = [o for o in built["observables"] if o["family"] == thin]
    assert thin_obs, "the thinnest family contributes no observable"
    for o in thin_obs:
        assert o["interval_quotable"] is False
        assert o["shifts"]["median"]["quotable"] is True
        assert o["shifts"]["low"]["quotable"] is False
        assert o["label"] in built["withheld_tails"]
    md = M.render(built)
    assert "withheld, and the reason is structural" in md
    assert "cannot be given a conditioned interval" in md


def test_a_withheld_tail_cannot_count_as_a_material_correction(raw):
    """Otherwise the page would report a correction drawn from a quantile it has
    just said it will not quote.

    Driven on a world where a tail is material AND withheld while its median is
    not material: the observable must then contribute nothing to the count.
    """
    built = M.assemble(raw)
    quoted = sum(1 for o in built["observables"]
                 if any(sh["material"] for sh in o["shifts"].values()
                        if sh["quotable"]))
    assert built["n_material"] == quoted, (
        "n_material is not counting exactly the quotable material shifts")

    # Admit a handful of draws chosen so the admissible MEDIAN barely moves
    # while the admissible TAIL moves a long way -- the tail is then material
    # and withheld, and must not be counted.
    keep = M.MIN_ADMISSIBLE_DRAWS - 1
    d = _doctor(raw)
    vals = sorted((float(r["pen_cns_bbb"]), i) for i, r in enumerate(d["rows"]))
    mid = len(vals) // 2
    chosen = {i for _, i in vals[mid - keep // 2: mid - keep // 2 + keep]}
    for i, r in enumerate(d["rows"]):
        for c in M.CONTROL_COLUMNS:
            r[c] = "0.000000" if i in chosen else "1.000000"
    out = M.assemble(d)
    o = _obs(out, "pen_cns_bbb")
    assert o["interval_quotable"] is False, "this fixture needs a withheld tail"
    withheld_material = [q for q, sh in o["shifts"].items()
                         if sh["material"] and not sh["quotable"]]
    assert withheld_material, (
        "this fixture was meant to make a WITHHELD tail material; it did not, "
        "so it is not exercising the rule")
    counted = any(sh["material"] for sh in o["shifts"].values() if sh["quotable"])
    assert out["n_material"] == sum(
        1 for x in out["observables"]
        if any(sh["material"] for sh in x["shifts"].values() if sh["quotable"]))
    if not counted:
        assert o["label"] not in [
            x["label"] for x in out["observables"]
            if any(sh["material"] for sh in x["shifts"].values() if sh["quotable"])], (
            "an observable whose only material shift is a withheld tail is "
            "being counted as a material correction")


def test_a_withheld_tail_cannot_place_the_point_estimate(built):
    """Three states, not two. Saying the published point falls OUTSIDE an
    interval the page refuses to quote would manufacture a finding out of a
    missing bound."""
    for o in built["observables"]:
        if o["interval_quotable"]:
            assert o["point"]["inside_admissible_interval"] in (True, False)
        else:
            assert o["point"]["inside_admissible_interval"] is None
    md = M.render(built)
    if any(o["point"]["inside_admissible_interval"] is None
           for o in built["observables"]):
        assert "cannot be placed either way" in md


def test_the_bliss_decomposition_survives_an_empty_admissible_subset(raw):
    """It reads the admissible minimum of a denominator, which does not exist
    when nothing is admissible. That crashed the page before it was guarded."""
    d = _doctor(raw)
    for r in d["rows"]:
        r[M.CONTROL_COLUMNS[0]] = "1.000000"
    out = M.assemble(d)                      # must not raise
    dec = out["bliss_decomposition"]
    assert dec["admissible_min_denominator"] is None
    assert dec["denominator_clears_floor"] is None
    assert dec["ratio_follows_denominator"] is False
    M.render(out)                            # and must still render


def test_a_short_per_draw_file_is_refused(raw):
    d = _doctor(raw)
    d["rows"] = d["rows"][:-1]
    with pytest.raises(SystemExit) as e:
        M.assemble(d)
    assert "rows against" in str(e.value)


def test_the_identity_check_skips_a_failed_run_rather_than_crashing(raw):
    d = _doctor(raw)
    for col in M.BLISS_PART_COLUMNS:
        d["rows"][11][col] = M.FAILED
    out = M.assemble(d)                      # must not raise
    assert out["bliss_decomposition"]["identity_max_error"] < 1e-6

def test_a_report_claiming_a_data_conditioned_headline_withholds_the_verdict(raw):
    """The restraint has to be drivable, or it is a sentence rather than a
    check. Doctored here because no state of this page's own data can produce
    it."""
    d = _doctor(raw)
    d["identifiability"]["data_constrained_in_production"] = 1
    out = M.assemble(d)
    assert out["data_conditioned_headlines"] == 1
    assert out["verdict"] is False, (
        "the page would still publish a verdict while the identifiability "
        "report claims a data-conditioned headline this page cannot have created")


def test_whether_a_direction_flips_is_measured_not_asserted(built, raw):
    """The sibling identifiability report's verdicts are DIRECTION claims, so
    the page has to measure sign changes rather than reason about them. Driven
    by forcing a flip."""
    assert built["sign_flips"] == [], (
        "a conditioned median changed sign, so a per-headline verdict in "
        "identifiability-report.md is now contradicted and must be revisited")
    assert "no conditioned median changes sign" in M.render(built)

    # Force one: negate the only signed observable so its admissible median
    # lands on the other side of zero from its full-prior median.
    # Negate a MINORITY of the family, so the full-prior median stays positive,
    # and admit exactly those draws. Negating half dragged the full median
    # negative too, and then there was no sign to change.
    d = _doctor(raw)
    n_tme = M.FAMILY_SAMPLES["tme"]
    vals = sorted((float(r["tme_hypoxia_gap"]), i)
                  for i, r in enumerate(d["rows"][:n_tme]))
    few = {i for _, i in vals[: n_tme // 5]}
    for i, r in enumerate(d["rows"]):
        if i in few:
            r["tme_hypoxia_gap"] = f"{-abs(float(r['tme_hypoxia_gap'])) - 1.0:.6f}"
        for c in M.CONTROL_COLUMNS:
            r[c] = "0.000000" if i in few else "1.000000"
    out = M.assemble(d)
    o = next(x for x in out["observables"] if x["column"] == "tme_hypoxia_gap")
    assert o["full"]["median"] > 0 > o["admissible_only"]["median"], (
        "this fixture must leave the full median positive and the admissible "
        "one negative, or there is no sign change to detect")
    assert out["sign_flips"], "forcing the signed observable negative flipped nothing"
    assert "have to be re-examined" in M.render(out)


def test_the_page_says_the_family_rows_are_nested(built):
    """Three rows in a table read as three independent estimates. Two of these
    families draw the SAME rows and the third draws a prefix of them, so saying
    so is the difference between one measurement and three."""
    md = M.render(built)
    assert "NESTED, not independent" in md
    counts = {M.FAMILY_SAMPLES[f] for f in M.FAMILY_SAMPLES}
    assert len(counts) < len(M.FAMILY_SAMPLES), (
        "every family now draws a different count, so the nesting sentence is "
        "no longer true as written")
    same = [f for f in M.FAMILY_SAMPLES
            if M.FAMILY_SAMPLES[f] == max(M.FAMILY_SAMPLES.values())]
    assert len(same) > 1
    # and families sharing a draw count really do get identical admissibility
    fractions = {built["families"][f]["fraction"] for f in same}
    assert len(fractions) == 1, (
        "families drawing the same rows report different admissible shares, so "
        "they are not actually reading the same draws")


def test_the_no_synergy_line_is_the_definitional_value_of_a_ratio(raw):
    """1.0 is not a figure anybody measured -- it is what a ratio equals when
    numerator and denominator agree. That only holds while the engine really
    computes a ratio, so the engine is checked rather than trusted."""
    assert M.BLISS_NO_SYNERGY == 1.0
    assert "cr.death_rate / bliss" in raw["rust"], (
        "sim-combo-mech no longer computes the synergy as a ratio of the "
        "combination kill to the Bliss prediction, so 1.0 is no longer the "
        "no-synergy line and the page's lower-bound section is wrong")


def test_the_sub_additive_corner_is_wholly_inadmissible(built):
    dec = built["bliss_decomposition"]
    assert dec["n_at_or_below_line"] > 0, (
        "no draw reaches the no-synergy line, so the published report's "
        "lower-bound claim has nothing to rest on and this section is vacuous")
    assert dec["n_at_or_below_line_admissible"] == 0
    assert dec["sub_additive_corner_is_wholly_inadmissible"] is True
    o = _obs(built, M.BLISS_COLUMN)
    assert o["admissible_only"]["p2_5"] > M.BLISS_NO_SYNERGY, (
        "the conditioned lower end does not clear the no-synergy line, so the "
        "section claiming it does is wrong")


def test_an_admissible_sub_additive_draw_flips_that_claim(raw):
    """Driven, because on the committed data the corner happens to be empty of
    admissible draws and a claim that cannot come back the other way is not a
    measurement."""
    d = _doctor(raw)
    target = next(i for i, r in enumerate(d["rows"])
                  if float(r[M.BLISS_COLUMN]) <= M.BLISS_NO_SYNERGY)
    for c in M.CONTROL_COLUMNS:
        d["rows"][target][c] = "0.000000"
    out = M.assemble(d)
    dec = out["bliss_decomposition"]
    assert dec["n_at_or_below_line_admissible"] >= 1
    assert dec["sub_additive_corner_is_wholly_inadmissible"] is False
    assert "not one of them is admissible" not in M.render(out)


def test_the_on_line_saturation_figure_is_measured(raw):
    """Reported as the weakest single-arm kill among those draws rather than
    classified against a saturation threshold nobody measured, so lowering one
    of them must move the figure."""
    base = M.assemble(raw)["bliss_decomposition"]
    assert base["n_on_line"] >= 1
    assert base["min_single_arm_kill_on_line"] is not None
    d = _doctor(raw)
    target = next(i for i, r in enumerate(d["rows"])
                  if float(r[M.BLISS_COLUMN]) == M.BLISS_NO_SYNERGY)
    d["rows"][target]["bliss_rate_rsl3"] = "0.100000"
    out = M.assemble(d)["bliss_decomposition"]
    assert out["min_single_arm_kill_on_line"] == pytest.approx(0.1), (
        "the saturation figure is not read from the draws it describes")


def test_the_draws_are_regenerated_from_the_published_seed_unmodified():
    """Source-pinned, because nothing in CI re-measures.

    `draw_matrix` runs only under `--measure`, so perturbing the seed changes
    no committed figure until someone re-measures -- at which point the
    reproduction gate would catch it, far too late to be a guard. A sweep that
    changed `SEED` to `SEED + 1` survived every behavioural test for exactly
    this reason.
    """
    tree = ast.parse(_src())
    fn = next(n for n in ast.walk(tree)
              if isinstance(n, ast.FunctionDef) and n.name == "draw_matrix")
    calls = [n for n in ast.walk(fn)
             if isinstance(n, ast.Call) and getattr(n.func, "id", None) == "sample_prior"]
    assert len(calls) == 1, "draw_matrix no longer makes exactly one draw call"
    call = calls[0]
    assert [a.id for a in call.args if isinstance(a, ast.Name)] == ["N_ROWS"], (
        "the draw count is no longer the module's own N_ROWS")
    seed = next((k.value for k in call.keywords if k.arg == "seed"), None)
    assert isinstance(seed, ast.Name) and seed.id == "SEED", (
        "the seed passed to sample_prior is not the published design's SEED "
        "imported unmodified, so these are not the draws the published "
        "intervals were computed from")


# --------------------------------------------------------------------------
# the note this page puts on the three published reports
# --------------------------------------------------------------------------
def test_the_published_reports_carry_a_note_derived_from_this_artifact(tmp_path,
                                                                       monkeypatch):
    """A hand-written note in those reports would be a second copy of this
    page's conclusions, free to go stale. Driven by pointing the note at a
    doctored artifact and requiring its figures to follow."""
    pub = _load("headline_uncertainty")
    real = json.loads(M.OUT_JSON.read_text())
    before = "\n".join(pub.conditioning_note("bliss"))
    assert str(real["families"]["bliss"]["n_admissible"]) in before

    doctored = json.loads(json.dumps(real))
    doctored["families"]["bliss"]["n_admissible"] = 7
    doctored["families"]["bliss"]["fraction"] = 7 / doctored["families"]["bliss"]["draws"]
    path = tmp_path / "doctored.json"
    path.write_text(json.dumps(doctored))
    monkeypatch.setattr(pub, "ADMISSIBILITY_JSON", path)
    after = "\n".join(pub.conditioning_note("bliss"))
    assert " 7 of " in after, "the note's admissible count is not read from the artifact"
    assert after != before


def test_the_note_fails_open_when_this_page_has_not_been_generated(tmp_path,
                                                                   monkeypatch):
    """The published generator must still run standalone. A note that raised
    would make this page a hard dependency of three reports that predate it."""
    pub = _load("headline_uncertainty")
    monkeypatch.setattr(pub, "ADMISSIBILITY_JSON", tmp_path / "absent.json")
    assert pub.conditioning_note("bliss") == []


def test_the_note_names_only_quantiles_this_page_will_quote():
    """It lists materially moved quantiles, and a withheld tail is not one this
    page quotes -- so the note must not surface a bound its own source refuses."""
    pub = _load("headline_uncertainty")
    built = M.assemble(M.load())
    withheld = {o["label"] for o in built["observables"] if not o["interval_quotable"]}
    for family in M.FAMILY_SAMPLES:
        text = "\n".join(pub.conditioning_note(family))
        for label in withheld:
            if label in text:
                assert "median" in text.split(label)[1].split("\n")[0], (
                    f"the note surfaces a withheld tail for {label}")


def test_a_thin_family_is_not_told_its_bounds_survived():
    """"No quantile moves materially" is the OPPOSITE of the truth for a family
    whose tails are withheld: they move a great deal and simply cannot be
    quoted. The note has to say which."""
    pub = _load("headline_uncertainty")
    built = M.assemble(M.load())
    thin = min(M.FAMILY_SAMPLES, key=lambda k: M.FAMILY_SAMPLES[k])
    assert any(not o["interval_quotable"] for o in built["observables"]
               if o["family"] == thin), "this fixture needs a withheld family"
    text = "\n".join(pub.conditioning_note(thin))
    assert "cannot be restated at all" in text
    assert "not a finding that the bounds survive" in text
    # and a family whose tails ARE quotable must not get that caveat
    fat = max(M.FAMILY_SAMPLES, key=lambda k: M.FAMILY_SAMPLES[k])
    assert "cannot be restated at all" not in "\n".join(pub.conditioning_note(fat))


def test_the_uncorrected_citation_sites_are_checked_not_remembered(built, raw):
    """The page names two handwritten sites whose conclusion rests on the
    inadmissible corner. A remembered list would outlive those sentences, so the
    naming is checked against their current text -- and driven here by claiming
    the phrase is gone."""
    assert built["subadditive_citations"], "the page names no citation site"
    assert all(built["subadditive_citations"].values()), (
        "a named site no longer carries the sentence this page points at; "
        "re-read it and update the note rather than leaving the pointer")
    md = M.render(built)
    assert "NOT edited here" in md
    for path in built["subadditive_citations"]:
        assert path in md

    d = _doctor(raw)
    d["subadditive_citations"] = {k: False for k in raw["subadditive_citations"]}
    assert "NOT edited here" not in M.render(M.assemble(d)), (
        "the page still names sites whose sentence has gone, which is a pointer "
        "to something that no longer exists")


def test_the_citation_check_can_answer_both_ways():
    """The predicate behind that note, driven in both directions -- and pinned
    to `load`, because hard-coding the answer there is unreachable from any
    guard that doctors the loaded inputs."""
    path, phrase = M.SUBADDITIVE_CITATIONS[0]
    assert M._phrase_present(path, phrase) is True
    assert M._phrase_present(path, "a sentence this file does not contain") is False
    tree = ast.parse(_src())
    fn = next(n for n in ast.walk(tree)
              if isinstance(n, ast.FunctionDef) and n.name == "load")
    called = {n.func.id for n in ast.walk(fn)
              if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)}
    assert "_phrase_present" in called, (
        "load() no longer consults the citation sites, so the page's note about "
        "them is remembered rather than checked")


def test_load_reports_a_missing_citation_as_missing(monkeypatch):
    """Drives `load` itself, which is the only thing that can see the answer
    being ignored.

    Pinning the CALL is not enough: `True or _phrase_present(...)` keeps the
    call and discards its result, and that mutation survived both the predicate
    unit test and the AST check. Pointing the loader at a phrase that is not
    there, and requiring it to say so, is what closes it.
    """
    path = M.SUBADDITIVE_CITATIONS[0][0]
    monkeypatch.setattr(M, "SUBADDITIVE_CITATIONS",
                        ((path, "a sentence this file does not contain"),))
    assert M.load()["subadditive_citations"] == {path: False}, (
        "the loader reports a citation present when its phrase is absent, so "
        "the page's note about those sites is not checked at all")
