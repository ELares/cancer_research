"""Guards for the cyst(e)inase dose-equivalence page.

The page's value is mostly in what it REFUSES, so the refusal is what needs
guarding hardest. Three shapes would make it hollow:

  * a dose-equivalence map that collapsed to the identity, making the whole
    composition vacuous;
  * an `admissible` flag read back from the artifact rather than recomputed,
    which is the `test_generator_guards_pin_source.py` defect -- the stored
    value cannot disagree with itself;
  * a refusal pinned as prose, which would survive the day a cascade became
    admissible and the refusal stopped being true.

Every one is checked against a recomputed quantity. The cascade measurements
need the compiled extension, which CI does not build, so those guards skip --
the arithmetic ones do not.
"""

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from cysteine_depletion_dose_equivalence import (  # noqa: E402
    MissingInputError,
    cystine_factor,
    erastin_dose_for_inhibition,
    escape_delivery,
    inhibition_from_depletion,
    load_inputs,
    saturation,
)

CAL = REPO_ROOT / "analysis" / "calibration"
RAW = CAL / "cysteine-depletion-dose-equivalence.json"
DOC = CAL / "cysteine-depletion-dose-equivalence.md"


def _data():
    return json.loads(RAW.read_text())


def test_the_hill_inversion_round_trips_through_the_PRODUCTION_map():
    """Round-trip through `calibrate_erastin.dose_to_inhib`, not a restatement.

    THE FIRST VERSION OF THIS GUARD RESTATED THE FORM INLINE as
    `dose**h / (dose**h + k)` -- the same wrong form the inverse used -- so it
    inverted this module's function with this module's own assumption and
    passed while every published dose was low by `k^(1-1/h)` = 4.451x. That is
    the guard-computes-its-own-expectation shape, and importing the production
    function is the only version that can fail.
    """
    from calibrate_erastin import dose_to_inhib  # noqa: PLC0415

    inp = _data()["inputs"]
    k, h = inp["k_erastin_uM"], inp["hill"]
    for i in (0.05, 0.134, 0.317, 0.5, 0.582, 0.807, 0.898, 0.979):
        dose = erastin_dose_for_inhibition(i, k, h)
        assert abs(dose_to_inhib(dose, k, h) - i) < 1e-12, (
            f"inhibition {i} did not round-trip through the production map")


def test_the_committed_k_is_the_half_maximal_dose():
    """The property that would have caught the wrong form in one line.

    For `dose^h/(dose^h + k^h)`, `k` IS the half-maximal dose by construction.
    The earlier inverse put it at 1.35 uM against a committed k of 6.
    """
    from calibrate_erastin import dose_to_inhib  # noqa: PLC0415

    inp = _data()["inputs"]
    k, h = inp["k_erastin_uM"], inp["hill"]
    assert abs(dose_to_inhib(k, k, h) - 0.5) < 1e-12, (
        "the production map does not put half-maximal inhibition at k")
    assert abs(erastin_dose_for_inhibition(0.5, k, h) - k) < 1e-12, (
        f"the inverse puts the half-maximal dose at "
        f"{erastin_dose_for_inhibition(0.5, k, h):.3f} uM, not the committed "
        f"k of {k}")


def test_the_endpoints_are_refused_rather_than_numbered():
    """A zero block needs no dose and a total block an infinite one."""
    assert erastin_dose_for_inhibition(0.0, 6.0, 6.0) is None
    assert erastin_dose_for_inhibition(1.0, 6.0, 6.0) is None
    assert erastin_dose_for_inhibition(-0.1, 6.0, 6.0) is None


def test_the_depletion_map_is_not_the_identity():
    """The #909 finding, re-pinned here because this page composes on top of it."""
    inp = _data()["inputs"]
    km, s0 = inp["km_uM"], inp["reference_uM"]
    for d in (0.1, 0.25, 0.5, 0.75, 0.9):
        assert inhibition_from_depletion(d, s0, km) < d


def test_the_erastin_window_is_narrow_and_inside_the_calibrated_grid():
    """The headline. Recomputed from the committed constants, not read."""
    d = _data()
    inp = d["inputs"]
    km, s0, k, h = (inp["km_uM"], inp["reference_uM"],
                    inp["k_erastin_uM"], inp["hill"])
    doses = [erastin_dose_for_inhibition(inhibition_from_depletion(x, s0, km), k, h)
             for x in (0.25, 0.50, 0.75, 0.90, 0.95, 0.99)]
    lo, hi = min(doses), max(doses)
    assert abs(d["erastin_window_uM"][0] - lo) < 1e-9
    assert abs(d["erastin_window_uM"][1] - hi) < 1e-9
    grid = inp["erastin_dose_grid_uM"]
    assert grid[0] <= lo and hi <= grid[-1], (
        f"the window {lo:.2f}-{hi:.2f} leaves the calibrated grid "
        f"{grid[0]}-{grid[-1]}, so the equivalence extrapolates")
    assert hi / lo < 5.0, (
        "the window is no longer narrow, which is the claim the Hill "
        "coefficient supports")


def test_the_escape_is_an_exact_no_op_at_a_full_floor():
    """Exact, not approximate, and independent of the inhibition entirely."""
    for i in (0.0, 0.137, 0.25, 0.5, 0.75, 0.898, 1.0):
        assert cystine_factor(i, 1.0) == 1.0, (
            f"at a full floor the factor must be exactly 1; inhibition {i} gave "
            f"{cystine_factor(i, 1.0)}")
    assert _data()["full_floor_factor_is_exactly_one"] is True


def test_the_delivered_fraction_does_not_depend_on_the_inhibition():
    """`1 - floor`, for every floor, and it must FALL as the floor rises."""
    vals = [escape_delivery(f) for f in (0.0, 0.25, 0.5, 0.75, 1.0)]
    assert vals == [1.0, 0.75, 0.5, 0.25, 0.0]
    # And the identity it rests on: supply removed is inhib*(1-floor).
    for i in (0.3, 0.9):
        for f in (0.0, 0.5, 1.0):
            assert abs((1.0 - cystine_factor(i, f)) - i * escape_delivery(f)) < 1e-12


def test_saturation_is_monotone_and_zero_at_zero():
    km = 179.4
    assert saturation(0.0, km) == 0.0 and saturation(-1.0, km) == 0.0
    v = [saturation(s, km) for s in (1, 10, 100, 1000)]
    assert v == sorted(v) and v[-1] < 1.0


def test_the_viability_constraint_is_the_repos_own_not_this_pages():
    """A page that invented its own threshold could refuse or not at will."""
    src = json.loads((REPO_ROOT / "analysis" / "headline-at-fitted-cascade.json")
                     .read_text())
    expected = src["ctrpv2_point"]["admissibility"]["constraint"]
    assert _data()["inputs"]["viability_constraint"] == expected
    # And it must be loaded, not defaulted: an absent record is a refusal.
    assert "viability_constraint" in load_inputs()


@pytest.mark.parametrize("name", ["in_vitro_fitted", "in_vivo_default"])
def test_the_admissibility_verdict_is_recomputed_from_the_controls(name):
    """RECOMPUTED. Reading back the stored flag compares the artifact to itself.

    This is the shape `tests/test_generator_guards_pin_source.py` exists to
    block: a verdict that travels with the numbers it judges cannot disagree
    with them.
    """
    d = _data()
    c = d.get("cascades")
    if not c:
        pytest.skip("cascade measurements not committed")
    m = c[name]
    limit = d["inputs"]["viability_constraint"]
    recomputed = max(m["untreated_controls"].values()) <= limit
    assert m["admissible"] is recomputed, (
        f"{name}: stored admissible={m['admissible']} but its own controls "
        f"(worst {max(m['untreated_controls'].values())}) against the "
        f"{limit} constraint give {recomputed}")


def test_neither_cascade_can_support_a_kill_number():
    """The refusal, DERIVED. If this ever fails the page must be rewritten.

    Not pinned as prose: the day a cascade is both admissible and responsive,
    this assertion is the thing that should break, not a sentence.
    """
    d = _data()
    c = d.get("cascades")
    if not c:
        pytest.skip("cascade measurements not committed")
    limit = d["inputs"]["viability_constraint"]
    usable = [n for n, m in c.items()
              if max(m["untreated_controls"].values()) <= limit
              and m["complete_block_death"] > 0.0]
    assert not usable, (
        f"{usable} is both admissible and responsive, so the magnitude is now "
        "predictable and this page's refusal is obsolete -- rewrite it rather "
        "than relaxing this guard")


def test_the_page_states_the_refusal_and_names_what_would_lift_it():
    doc = DOC.read_text()
    assert "refuses to say" in doc.lower()
    assert "No absolute kill rate" in doc
    assert "would unblock" in doc.lower(), (
        "a refusal that does not say what would lift it is a dead end, not a "
        "finding")


def test_the_measured_no_op_agrees_with_the_algebra():
    """The identity is algebra; this checks the engine actually behaves that way."""
    c = _data().get("cascades")
    if not c:
        pytest.skip("cascade measurements not committed")
    for name, m in c.items():
        assert m["full_floor_is_exact_noop"] is True, (
            f"{name}: a complete block at a full floor did not reproduce the "
            "untreated rate exactly, so the factor identity does not hold in "
            "the engine")
        assert m["complete_block_at_full_floor_death"] == m["untreated_glycolytic"]


def test_it_refuses_rather_than_rendering_from_absent_inputs(tmp_path, monkeypatch):
    import cysteine_depletion_dose_equivalence as m
    monkeypatch.setattr(m, "CONST_SRC", tmp_path / "nope.csv")
    with pytest.raises(MissingInputError):
        m.load_inputs()


def test_a_calibration_without_the_hill_fit_is_refused(tmp_path, monkeypatch):
    import cysteine_depletion_dose_equivalence as m
    bad = tmp_path / "erastin.json"
    bad.write_text(json.dumps({"calibrated_params": {"lp_rate": 0.4}}))
    monkeypatch.setattr(m, "ERASTIN_SRC", bad)
    with pytest.raises(MissingInputError):
        m.load_inputs()


def test_render_only_reproduces_the_page():
    before = DOC.read_text()
    r = subprocess.run(
        [sys.executable,
         str(REPO_ROOT / "scripts" / "cysteine_depletion_dose_equivalence.py"),
         "--render-only"],
        capture_output=True, text=True, cwd=REPO_ROOT)
    assert r.returncode == 0, r.stderr
    assert DOC.read_text() == before


def test_the_committed_constants_are_the_same_ones_the_feasibility_page_uses():
    """Two pages, one Km. A second copy is a second thing to go stale."""
    rows = {r["constant"]: r for r in
            csv.DictReader((CAL / "cysteine_depletion_constants.csv").open())}
    inp = _data()["inputs"]
    assert inp["km_uM"] == float(rows["xct_cystine_km"]["value"])
    expected_ref = (float(rows["rpmi1640_cystine_2hcl_mass"]["value"])
                    / float(rows["cystine_2hcl_molar_mass"]["value"]) * 1000.0)
    assert abs(inp["reference_uM"] - expected_ref) < 1e-9


def test_the_constraint_follows_the_repos_record_rather_than_a_literal(tmp_path, monkeypatch):
    """BEHAVIOURAL, because the value-level check cannot see a hardcoded literal.

    A mutation planting the correct `0.02` directly survived the first sweep:
    comparing the stored value to the repo's record passes whether it was read
    or coincidentally typed. Pointing the loader at a doctored record and
    requiring the loaded value to MOVE is what distinguishes them.
    """
    import cysteine_depletion_dose_equivalence as m
    doctored = tmp_path / "admissibility.json"
    doctored.write_text(json.dumps(
        {"ctrpv2_point": {"admissibility": {"constraint": 0.07}}}))
    monkeypatch.setattr(m, "ADMISSIBILITY_SRC", doctored)
    assert m.load_inputs()["viability_constraint"] == 0.07, (
        "the constraint did not follow the record, so it is not being read")


def test_a_record_without_a_constraint_is_refused(tmp_path, monkeypatch):
    """Absent must refuse, not default: a default would silently set the bar."""
    import cysteine_depletion_dose_equivalence as m
    empty = tmp_path / "admissibility.json"
    empty.write_text(json.dumps({"ctrpv2_point": {"admissibility": {}}}))
    monkeypatch.setattr(m, "ADMISSIBILITY_SRC", empty)
    with pytest.raises(MissingInputError):
        m.load_inputs()


def test_the_page_does_not_inherit_its_row_order_from_serialisation():
    """`json.dumps(sort_keys=True)` reorders a stored dict.

    A table walking `.items()` then publishes whatever serialisation produced,
    and `--render-only` can differ from a fresh run. Third instance of this
    defect in this repo; `test_artifact_freshness.py` caught this one.
    """
    import random

    import cysteine_depletion_dose_equivalence as m

    def shuffled(o, rng):
        if isinstance(o, dict):
            items = list(o.items())
            rng.shuffle(items)
            return {k: shuffled(v, rng) for k, v in items}
        if isinstance(o, list):
            return [shuffled(v, rng) for v in o]
        return o

    stored = _data()
    base = m.render(m.assemble(stored))
    for seed in range(8):
        alt = m.render(m.assemble(
            shuffled(json.loads(json.dumps(stored)), random.Random(seed))))
        assert alt == base, f"render depends on dict order (seed {seed})"


def test_an_undeclared_cascade_is_refused_rather_than_ordered_arbitrarily():
    """A new cascade must be given a place in the argument, not appended by luck."""
    import cysteine_depletion_dose_equivalence as m

    stored = json.loads(json.dumps(_data()))
    if "cascades" not in stored:
        pytest.skip("cascade measurements not committed")
    stored["cascades"]["mystery"] = stored["cascades"]["in_vivo_default"]
    with pytest.raises(MissingInputError):
        m.render(m.assemble(stored))


def test_the_ceiling_is_cross_checked_against_the_committed_calibration():
    """The ceiling stops being pure topology, and the harness is validated.

    Two claims in one: a complete cystine block must reproduce the committed
    saturating-erastin increment, which confirms the ceiling AND confirms this
    page drives the engine the way `calibrate_erastin.py` did. A harness that
    disagreed would be measuring something else and every figure on the page
    would be suspect.
    """
    x = _data().get("ceiling_cross_check")
    if not x:
        pytest.skip("cascade measurements not committed")
    # Recomputed from the committed calibration, not read back.
    er = json.loads((CAL / "erastin-calibration.json").read_text())
    assert x["committed_top_dose_increment"] == er["erastin_increment_top_dose"]
    assert x["absolute_disagreement"] < 0.02, (
        f"a complete block differs from the committed saturating-erastin "
        f"increment by {x['absolute_disagreement']:.4f}; the harness is not "
        "reproducing the committed calibration")
    # And the premise: the committed top dose really is a near-complete block.
    assert x["inhibition_at_committed_top_dose"] > 0.99, (
        "the committed top dose is not a near-complete block, so its increment "
        "cannot stand in for one")
    # Non-vacuity: it must not be EXACTLY complete either, or the page would be
    # rounding a real asymptote to 100%.
    assert x["inhibition_at_committed_top_dose"] < 1.0
