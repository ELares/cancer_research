"""Guards for the per-dead-cell DAMP measurement (#928).

Every guard on a DERIVED field drives `assemble()` on doctored raw input rather
than reading the committed artifact, because a test that reads only a generated
file cannot fail when its generator changes -- the artifact does not move until
someone re-runs it.
"""
from __future__ import annotations

import importlib.util
import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "scripts" / "icd_damp_per_death.py"
REPORT = REPO / "analysis" / "calibration" / "icd-damp-per-death.md"
ARTIFACT = REPO / "analysis" / "calibration" / "icd-damp-per-death.json"

CORE_IMMUNE = REPO / "simulations" / "ferroptosis-core" / "src" / "immune.rs"
ICD_MAIN = REPO / "simulations" / "sim-icd" / "src" / "main.rs"
ICD_README = REPO / "simulations" / "sim-icd" / "README.md"
PANEL_MAIN = REPO / "simulations" / "sim-modality-panel" / "src" / "main.rs"


def _src() -> str:
    return SCRIPT.read_text(encoding="utf-8")


def _mod():
    spec = importlib.util.spec_from_file_location("icd_damp_per_death", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def mod():
    return _mod()


@pytest.fixture(scope="module")
def stored():
    assert ARTIFACT.is_file(), f"missing {ARTIFACT}; run the generator with --measure"
    return json.loads(ARTIFACT.read_text())


def _raw(stored: dict) -> dict:
    """The measured legs only, with every derived field dropped."""
    return {k: json.loads(json.dumps(stored[k])) for k in ("icd", "panel", "spatial")}


# --------------------------------------------------------------------------
# The quoted predictions must still be what their sources say
# --------------------------------------------------------------------------

def test_the_documented_predictions_are_still_the_quoted_ones(mod):
    """The numbers live in the generator; the QUOTES are pinned to their files.

    Without this the page could go on comparing measurements against a band
    nobody claims any more -- and deleting the predictions from those documents
    would erase the falsification along with the error, so they are kept and
    checked instead.
    """
    assert mod.PREDICTIONS, "no predictions declared"
    for p in mod.PREDICTIONS:
        text = (REPO / p["source"]).read_text(encoding="utf-8")
        assert p["quote"] in text, (
            f"{p['source']} no longer contains the quoted prediction "
            f"{p['quote']!r}; this page compares against it")
        # The declared band must be the one the quote states, or the generator
        # is measuring against a number of its own.
        assert str(p["low"]).rstrip("0").rstrip(".") in p["quote"]
        assert str(p["high"]).rstrip("0").rstrip(".") in p["quote"]


# --------------------------------------------------------------------------
# The headline: direction holds, magnitude does not
# --------------------------------------------------------------------------

def test_every_regime_falls_short_of_the_weakest_prediction(mod, stored):
    d = mod.assemble(_raw(stored))
    weakest = min(p["low"] for p in mod.PREDICTIONS)
    assert d["ratio_max"] < weakest, (
        f"measured max {d['ratio_max']} reaches the weakest documented "
        f"prediction {weakest}; the page's whole finding is that it does not")
    assert d["every_prediction_refuted"] is True
    assert d["verdict"] == "direction_holds_magnitude_refuted"


def test_a_ratio_reaching_the_band_would_stop_the_refutation(mod, stored):
    """The refutation is DERIVED. Hand it a regime that reaches 1.5x."""
    raw = _raw(stored)
    for row in raw["icd"]["rows"]:
        if row["phenotype"] == "Persister" and row["treatment"] == "SDT":
            row["damp_per_dead_cell"] = row["damp_per_dead_cell"] * 3.0
    d = mod.assemble(raw)
    assert d["ratio_max"] > 1.5
    assert d["every_prediction_refuted"] is False
    assert d["verdict"] != "direction_holds_magnitude_refuted"
    assert mod.render(d).count("**yes**") >= 1


def test_the_direction_is_reported_separately_from_the_magnitude(mod, stored):
    """Direction and magnitude are independent conclusions, so one must be
    able to fail without the other. A single inverted row flips only the
    direction."""
    raw = _raw(stored)
    raw["spatial"]["arms"]["SDT"]["mean_release_lp"] = (
        raw["spatial"]["arms"]["RSL3"]["mean_release_lp"] * 0.5)
    d = mod.assemble(raw)
    assert d["direction_holds_everywhere"] is False
    assert d["every_prediction_refuted"] is True, (
        "the magnitude verdict must not move when only the direction does")
    assert d["verdict"] == "direction_does_not_hold"
    assert "does NOT hold in every readable comparison" in mod.render(d)


# --------------------------------------------------------------------------
# The mechanism: at the moment of death there is no contrast to read
# --------------------------------------------------------------------------

def test_both_arms_are_near_the_threshold_at_death(mod, stored):
    d = mod.assemble(_raw(stored))
    sp = d["spatial"]["arms"]
    # NON-VACUITY FIRST. If death LP and release LP were the same number the
    # two-time-point comparison would be empty and every assertion below would
    # pass for the wrong reason.
    for arm in ("SDT", "RSL3"):
        assert sp[arm]["mean_death_lp"] != sp[arm]["mean_release_lp"]
        assert sp[arm]["mean_release_lp"] > sp[arm]["mean_death_lp"]
    assert d["both_arms_near_threshold_at_death"] is True
    assert d["at_death_contrast_exists"] is False
    assert abs(d["at_death_ratio"] - 1.0) < 0.01, (
        "the two arms must be within 1% of each other at death for the "
        "'no contrast at that moment' reading to hold")
    assert d["release_ratio"] > d["at_death_ratio"], (
        "the gap must be larger after death than at it, which is the claim")


def test_an_arm_far_from_threshold_at_death_would_restore_the_contrast(mod, stored):
    raw = _raw(stored)
    raw["spatial"]["arms"]["SDT"]["mean_death_lp"] = 30.0
    d = mod.assemble(raw)
    assert d["both_arms_near_threshold_at_death"] is False
    assert d["at_death_contrast_exists"] is True
    assert "IS readable at that moment" in mod.render(d)


def test_the_death_threshold_is_read_from_the_engine_not_hardcoded(mod, stored, tmp_path):
    """`near threshold` is a claim about a number that lives in the engine.

    A literal 10.0 here would stop tracking it the moment anyone retuned the
    switch, so the value is parsed -- and this points the parser at a different
    threshold and requires the derived flags to move.
    """
    assert mod.death_threshold() == 10.0
    fake = tmp_path / "params.rs"
    fake.write_text("            death_threshold: 40.0,\n", encoding="utf-8")
    original = mod.PARAMS_RS
    try:
        mod.PARAMS_RS = fake
        d = mod.assemble(_raw(stored))
        assert d["death_threshold"] == 40.0
        assert d["both_arms_near_threshold_at_death"] is False, (
            "at a threshold of 40 a death LP of ~10.25 is nowhere near it")
    finally:
        mod.PARAMS_RS = original


# --------------------------------------------------------------------------
# Absence is not zero
# --------------------------------------------------------------------------

def test_an_absent_cohort_is_reported_as_absent_and_not_as_zero(mod, stored):
    d = mod.assemble(_raw(stored))
    assert d["panel_rows_with_absent_cohort"] > 0, (
        "no absent cohorts at all would make this guard vacuous")
    assert d["absent_cohorts_all_have_zero_deaths"] is True
    # And the sentinel must be gone from the producing binary.
    panel = PANEL_MAIN.read_text(encoding="utf-8")
    assert "fn damp_per_death(lps: &[f64], immune: &ImmuneParams) -> Option<f64>" in panel
    assert '"deaths": {' in panel, (
        "the denominator must be emitted beside the value, or a null is "
        "unreadable")


def test_a_null_without_a_zero_denominator_is_refused(mod, stored):
    """A `null` must mean 'no cohort'. Give one a non-zero kill count."""
    raw = _raw(stored)
    for row in raw["panel"]["rows"]:
        if row["damp_per_death"].get("RSL3") is None:
            row["deaths"]["RSL3"] = 17
            break
    else:
        pytest.fail("no null cohort in the stored panel rows to doctor")
    d = mod.assemble(raw)
    assert d["absent_cohorts_all_have_zero_deaths"] is False
    assert "NOT all of them at a kill count of zero" in mod.render(d)


# --------------------------------------------------------------------------
# The two readings that fell out of the measurement
# --------------------------------------------------------------------------

def test_the_panel_sweep_is_counted_once(mod, stored):
    """`ratios` is a display table whose panel entry summarises many rows.

    Summing it with the per-row ratios double-counts the sweep, which the first
    version of this page did -- it reported 36 comparisons where there are 35.
    """
    d = mod.assemble(_raw(stored))
    summary_rows = sum(1 for r in d["ratios"] if "span" in r)
    assert summary_rows == 1
    expected = (len(d["ratios"]) - summary_rows) + d["panel_readable_rows"]
    assert d["readable_comparisons"] == expected


def test_the_adc_inversion_is_measured_not_asserted(mod, stored):
    d = mod.assemble(_raw(stored))
    assert d["adc_outranks_sdt_rows"] > 0
    assert d["adc_outranks_sdt_is_hypoxic"] is True
    raw = _raw(stored)
    for row in raw["panel"]["rows"]:
        q = row["damp_per_death"]
        if q.get("AntibodyDrugConjugate") and q.get("SDT"):
            q["AntibodyDrugConjugate"] = q["SDT"] * 0.5
    d2 = mod.assemble(raw)
    assert d2["adc_outranks_sdt_rows"] == 0
    assert "not top of the per-death ranking" not in mod.render(d2)


def test_the_saturation_inversion_is_derived_from_the_death_rates(mod, stored):
    d = mod.assemble(_raw(stored))
    assert d["immune_kill_ranking_inverts_under_saturation"] is True
    k = d["persister_immune_kills"]
    assert k["SDT"] < k["RSL3"]
    raw = _raw(stored)
    for row in raw["icd"]["rows"]:
        if row["phenotype"] == "Persister" and row["treatment"] == "SDT":
            row["death_rate"] = 0.5
    d2 = mod.assemble(raw)
    assert d2["immune_kill_ranking_inverts_under_saturation"] is False, (
        "the inversion is only reportable as structural while the arm "
        "saturates; without that it is just a lower number")
    assert "ranks survivors" not in mod.render(d2)


# --------------------------------------------------------------------------
# Generator hygiene
# --------------------------------------------------------------------------

def test_render_only_reassembles_rather_than_echoing(mod):
    """`--render-only` must recompute the derived block, not re-print it.

    Echoing stored derived fields would make every guard above inert on a
    re-render, which `tests/test_artifact_freshness.py` caught on a sibling
    generator.
    """
    import ast

    body = _src().split('"""', 2)[2]
    assert len(body) > 2000, "docstring split returned a near-empty body"

    tree = ast.parse(_src())
    main = next(n for n in tree.body
                if isinstance(n, ast.FunctionDef) and n.name == "main")
    branch = next(
        n for n in ast.walk(main)
        if isinstance(n, ast.If) and isinstance(n.test, ast.Attribute)
        and n.test.attr == "render_only")

    def calls(nodes, name):
        return [c for stmt in nodes for c in ast.walk(stmt)
                if isinstance(c, ast.Call) and getattr(c.func, "id", None) == name]

    assert calls(branch.body, "assemble"), (
        "the --render-only branch does not call assemble(), so the documented "
        "command re-renders stored derived fields instead of recomputing them")
    reads = {n.slice.value for stmt in branch.body for n in ast.walk(stmt)
             if isinstance(n, ast.Subscript) and isinstance(n.slice, ast.Constant)
             and isinstance(n.slice.value, str)}
    assert {"icd", "panel"} <= reads, (
        f"render-only must rebuild its input from the stored MEASURED legs; "
        f"it reads {sorted(reads)}")
    assert "spatial" not in reads, (
        "the spatial leg must be recomputed, never read back from the artifact")
    # Recomputed in BOTH arms, so it can never be echoed on either path.
    assert calls(branch.body, "measure_spatial"), "render-only echoes the spatial leg"
    assert calls(branch.orelse, "measure_spatial"), "the measure path skips the spatial leg"


def test_the_report_states_the_claim_it_does_not_make(mod, stored):
    md = mod.render(mod.assemble(_raw(stored)))
    for clause in (
        "uncalibrated",
        "not a measurement of immunogenicity",
        "Nothing is recalibrated",
        "the manuscript is not edited",
        "that is the kill count",
        # Self-review of #928 itself: "all 35 readable comparisons" invites
        # reading the 32 panel rows as 32 independent tests when they share one
        # biochemical path with sim-icd. The page has to say so.
        "not independent replicates",
    ):
        assert clause in md, f"the report no longer states: {clause}"


# --------------------------------------------------------------------------
# The engine must no longer state the refuted magnitude
# --------------------------------------------------------------------------

REFUTED_PHRASES = (
    "release MORE DAMPs per cell",
    "release MORE DAMPs per dead cell",
    "dies with far more LP",
    "dies carrying far more lipid peroxidation",
    "LP at death is MUCH higher than threshold",
    "LP at death is NEAR threshold",
)


def _comment_blocks(text: str) -> list[str]:
    """Runs of consecutive Rust doc-comment lines, as one string each."""
    blocks: list[list[str]] = []
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("///") or s.startswith("//!") or s.startswith("//"):
            if blocks and blocks[-1] is not None:
                blocks[-1].append(s)
            else:
                blocks.append([s])
        elif blocks and blocks[-1] is not None:
            blocks.append(None)  # terminator
    return ["\n".join(b) for b in blocks if b]


@pytest.mark.parametrize("path", [CORE_IMMUNE, ICD_MAIN, PANEL_MAIN])
def test_no_engine_source_still_asserts_the_refuted_magnitude(path):
    """The retraction has to reach every site, not the one that was edited.

    These are the three sources that carried the claim: the library function
    itself, the binary whose README quantified it, and the panel that restated
    it twice.

    The old wording is QUOTED in each correction rather than deleted, because a
    reader needs to know what the comment used to say. So the check is not that
    the phrase is absent -- that would force the corrections to paraphrase -- but
    that it never appears in a comment block which does not also carry the
    measurement that retracts it.
    """
    text = path.read_text(encoding="utf-8")
    blocks = _comment_blocks(text)
    assert blocks, f"no comment blocks parsed out of {path.name}"
    found = 0
    for block in blocks:
        for phrase in REFUTED_PHRASES:
            if phrase in block:
                found += 1
                assert "#928" in block, (
                    f"{path.name} states {phrase!r} in a comment block that does "
                    f"not carry the measurement retracting it:\n{block[:400]}")
    assert found > 0, (
        f"{path.name} quotes none of the refuted phrases, so this guard is "
        "vacuous there -- the correction should say what it is correcting")


def test_the_sim_icd_field_is_named_for_the_quantity_it_holds():
    """`sim_cell` returns `final_lp`; the field used to be `avg_lp_at_death`."""
    main = ICD_MAIN.read_text(encoding="utf-8")
    assert '"avg_final_lp": avg_final_lp' in main
    assert '"avg_lp_at_death"' not in main
    readme = ICD_README.read_text(encoding="utf-8")
    assert "| avg_final_lp |" in readme
    assert "| avg_lp_at_death |" not in readme, (
        "the README data dictionary still documents the old field name")


def test_the_frozen_core_source_is_registered_before_being_edited():
    """`immune.rs` is a declared numerical input of eight committed archives.

    Editing it -- even a comment -- makes their declared input unresolvable
    unless the superseded bytes are registered, which is the gate working as
    designed rather than an obstacle.
    """
    manifest = json.loads(
        (REPO / "analysis/calibration/frozen-sources/manifest.json").read_text())
    rel = "simulations/ferroptosis-core/src/immune.rs"
    assert rel in manifest["sources"], (
        f"{rel} was edited without registering its superseded bytes")
    for digest in manifest["sources"][rel]:
        snapshot = REPO / "analysis/calibration/frozen-sources" / f"{digest}.source"
        assert snapshot.is_file(), f"missing snapshot for {digest}"


def test_the_committed_report_matches_the_committed_measurements(mod, stored):
    """Regenerate-and-compare: the prose must be the one these numbers produce."""
    assert REPORT.read_text(encoding="utf-8") == mod.render(mod.assemble(_raw(stored)))


# --------------------------------------------------------------------------
# The registered prediction this measurement lands on
# --------------------------------------------------------------------------

PREREG = REPO / "PREREGISTRATION.md"


def test_the_p5_band_is_parsed_from_the_registration(mod, stored):
    """P5 registers a falsification threshold on exactly this quantity.

    The band is parsed, not typed: a literal 1.5 here would stop tracking the
    registered number, and the registered number is the whole point.
    """
    d = mod.assemble(_raw(stored))
    assert d["p5"]["band"] == 1.5
    assert mod.P5_QUOTE in PREREG.read_text(encoding="utf-8"), (
        "P5's registered falsification sentence has changed; re-read it before "
        "editing this guard, because this page measures the model against it")


def test_a_different_registered_band_would_move_the_verdict(mod, stored, tmp_path):
    """Asserting `band == 1.5` cannot tell a PARSE from a literal that happens
    to be right today, so the parser is repointed at a different registration
    and the derived flag must follow it."""
    fake = tmp_path / "prereg.md"
    fake.write_text(
        "- *Falsification threshold:* DAMP release and DC maturation per dead "
        "cell agree within 1.02x between dense (SDT or 3D) and sparse (RSL3) "
        "kill.\n", encoding="utf-8")
    original = mod.PREREG
    try:
        mod.PREREG = fake
        d = mod.assemble(_raw(stored))
        assert d["p5"]["band"] == 1.02
        assert d["p5"]["both_legs_inside_band"] is False, (
            "at a 1.02x band the measured 1.05-1.13x readouts sit OUTSIDE it")
    finally:
        mod.PREREG = original


def test_a_registration_without_the_band_is_refused(mod, stored, tmp_path):
    """A missing band must stop the page, not let it supply its own number."""
    fake = tmp_path / "prereg.md"
    fake.write_text("P5 has no falsification threshold any more.\n",
                    encoding="utf-8")
    original = mod.PREREG
    try:
        mod.PREREG = fake
        with pytest.raises(RuntimeError, match="registered falsification band"):
            mod.assemble(_raw(stored))
    finally:
        mod.PREREG = original


def test_both_registered_readouts_sit_inside_the_registered_band(mod, stored):
    """The finding: the MODEL satisfies P5's falsification condition.

    Both readouts are required, because the registration names both -- a guard
    on DAMP alone would miss the DC-maturation leg entirely.
    """
    d = mod.assemble(_raw(stored))
    legs = d["p5"]["legs"]
    assert legs, "no readable phenotype, so this guard would be vacuous"
    for lg in legs:
        assert 1.0 < lg["damp_per_death_ratio"] < d["p5"]["band"]
        assert 1.0 < lg["dc_maturation_per_death_ratio"] < d["p5"]["band"]
    assert d["p5"]["both_legs_inside_band"] is True
    md = mod.render(d)
    assert "satisfies P5's registered falsification condition" in md
    # The restraint is guarded as hard as the finding.
    assert "nothing here falsifies the biology" in md
    assert "THE THRESHOLD IS NOT AMENDED" in md.upper()


def test_a_readout_outside_the_band_would_be_reported_as_such(mod, stored):
    """Derived, so a model that DID separate the arms per cell must read
    differently -- including when only the DC leg separates, which a guard on
    DAMP alone could not tell apart."""
    raw = _raw(stored)
    for row in raw["icd"]["rows"]:
        if row["phenotype"] == "Persister" and row["treatment"] == "SDT":
            row["mature_dcs"] = row["mature_dcs"] * 4.0
    d = mod.assemble(raw)
    assert d["p5"]["both_legs_inside_band"] is False
    md = mod.render(d)
    assert "remains a model prediction" in md
    assert "satisfies P5's registered falsification condition" not in md


def test_the_registration_carries_the_measurement_without_amending_it():
    """The note must record the measurement and leave the threshold alone.

    Amending a registered threshold after measuring against it is what
    preregistration exists to prevent, so the registered sentence has to be
    present verbatim alongside the note.
    """
    text = PREREG.read_text(encoding="utf-8")
    assert ("DAMP release and DC maturation per dead cell agree within 1.5x"
            in text), "the registered threshold is no longer stated verbatim"
    assert "THE THRESHOLD ABOVE IS NOT AMENDED" in text
    assert "#928" in text
    assert "1.1332x" in text and "1.0950x" in text, (
        "the note must carry both measured readouts, not only the DAMP one")
    assert "does NOT falsify the biological hypothesis" in text, (
        "the note must keep the model/biology distinction the finding rests on")


def test_the_page_says_which_leg_carries_the_p5_conclusion(mod, stored):
    """Self-review of #928: P5's dense-versus-sparse contrast is SPATIAL, and
    two of the three measurement sites have no geometry at all. The page must
    name the leg whose mechanism matches the registration rather than letting
    35 comparisons stand in for it."""
    md = mod.render(mod.assemble(_raw(stored)))
    assert "P5's contrast is DENSE versus SPARSE kill" in md
    assert "`sim-icd` has no geometry at all" in md
    assert "in the SAME geometry" in md
    # And the loose reading has to stay out.
    assert "the null P5 would be falsified by" not in md
    assert "not a null, but not distinguishable from one" in md
    assert "SAME ratio to five decimals" in md, (
        "the two DC normalisations agreeing is the check that the readout is "
        "per-cell at all; the page must state it")
    text = PREREG.read_text(encoding="utf-8")
    assert "spatial property" in text and "1.0469x" in text, (
        "the registration note must carry the same scoping, not only the page")


def test_the_dc_leg_is_checked_through_two_normalisations(mod, stored):
    """`mature_dcs / kills` and `dc_activation_fraction` must give the same
    ratio, or the DC figure is not a per-cell readout and P5's DC leg is not
    being measured at all. Derived, and the disagreeing branch is driven."""
    d = mod.assemble(_raw(stored))
    assert d["p5"]["dc_normalisations_agree"] is True
    assert "SAME ratio to five decimals" in mod.render(d)
    raw = _raw(stored)
    for row in raw["icd"]["rows"]:
        if row["phenotype"] == "Persister" and row["treatment"] == "SDT":
            row["dc_activation_fraction"] = row["dc_activation_fraction"] * 2.0
    d2 = mod.assemble(raw)
    assert d2["p5"]["dc_normalisations_agree"] is False
    md = mod.render(d2)
    assert "They DISAGREE" in md
    assert "SAME ratio to five decimals" not in md
