"""Tests for the practical-identifiability synthesis (#503).

Pure-Python stdlib (no compiled extension), so these run in CI. They guard the
load-bearing facts (degrees of freedom, the non-identifiable count, that no
headline is point-estimable, that zero headlines are data-conditioned in the
production regime) and that the committed report is not stale.
"""

import json
import sys
import pytest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import identifiability_report as ir  # noqa: E402


def test_degrees_of_freedom_matches_prcc():
    """The DOF count must equal the actual PRCC swept-parameter set (11), so the
    accounting cannot drift away from the real parameter set."""
    r = ir.build()
    assert r["degrees_of_freedom"] == 11
    assert len(r["swept_parameters"]) == 11


def test_non_identifiable_params_are_real_and_count_holds():
    """6 of 11 parameters are non-identifiable from the kill rate, and each named
    one is a genuine PRCC parameter (build() asserts this; re-check here)."""
    r = ir.build()
    sc = r["single_cell_sobol"]
    assert sc["non_identifiable_count"] == 6
    for p in sc["non_identifiable"]:
        assert p in r["swept_parameters"], f"{p} not a swept parameter"
    # The dominant three are also real parameters and dominate.
    for p in sc["dominant_ST"]:
        assert p in r["swept_parameters"]
    assert sc["dominant_ST"]["lp_propagation"] > sc["dominant_ST"]["gpx4_rate"]


def test_zero_headlines_data_conditioned_in_production():
    """The load-bearing finding: no headline is conditioned on data in the
    regime that produces it (the only fit is the disjoint in-vitro switch)."""
    r = ir.build()
    assert r["data_constrained_in_production"] == 0


def test_no_headline_is_point_estimable():
    """Uncertain direction is permitted; no headline is point-estimable."""
    r = ir.build()
    allowed = {"directional_only", "direction_robust_magnitude_not",
               "direction_and_magnitude_uncertain",
               # #926: robust, but entailed by a monotone dose-response at
               # three fixed concentrations rather than measured. Still not
               # point-estimable, which is what this guard is about.
               "direction_robust_by_construction_magnitude_not",
               # #940: robust only once the model's own untreated-viability
               # constraint removes the sub-additive corner. Conditioning adds
               # no data, so this is further from point-estimable than the
               # unconditional label above it, not closer.
               "direction_robust_only_under_admissibility_magnitude_not"}
    for h in r["headlines"]:
        assert h["verdict"] in allowed, f"{h['key']} claims {h['verdict']}"
    # At least the single-cell kill rate and the immune ratio are directional-only.
    by_key = {h["key"]: h for h in r["headlines"]}
    assert by_key["single_cell_kill_rate"]["verdict"] == "directional_only"
    assert by_key["immune_amplification_ratio"]["verdict"] == "directional_only"


def _published(source, minimum, lower):
    source.write_text(
        "| default point estimate (no perturbation) | 1.992 |\n"
        "| prior-predictive median | 1.348 |\n"
        f"| 95% prior-predictive interval | [{lower:.3f}, 5.242] |\n"
        f"| full range (min, max) | [{minimum:.3f}, 7.800] |\n"
    )


def _admissibility(path, *, at_or_below, admissible_corner, reproduces=True,
                   admissible_min=1.027, drop=None):
    """An admissibility artifact with a chosen corner, for #940's second input.

    Written field by field rather than copied from the committed artifact, so a
    fixture cannot pass by echoing the real measurement.
    """
    doc = {
        "bliss_decomposition": {
            "n_at_or_below_line": at_or_below,
            "n_at_or_below_line_admissible": admissible_corner,
            "n_on_line": 2,
            "min_single_arm_kill_on_line": 1.0,
            "no_synergy_line": 1.0,
            "numerator_ratio": 0.740,
            "denominator_ratio": 0.374,
            "corner_min_worst_untreated": 0.2975,
        },
        "constraint": {"baseline_max": 0.02},
        "observables": [{
            "column": "bliss_synergy",
            "reproduces": reproduces,
            "reproduction": {"min": {"reproduces": reproduces}},
            "n_admissible": 125,
            "n_full": 300,
            "interval_quotable": True,
            "admissible_only": {"min": admissible_min},
        }],
    }
    if drop:
        section, key = drop
        target = (doc["observables"][0] if section == "observable"
                  else doc["bliss_decomposition"])
        assert key in target, f"the fixture cannot drop absent key {key!r}"
        target.pop(key)
    path.write_text(json.dumps(doc))
    return path


def test_bliss_verdict_tracks_the_sampled_range_not_a_rounded_null(tmp_path, monkeypatch):
    """Exercise both sides of the inference, including an exactly additive draw.

    #940 added a SECOND input, so each fixture now sets both. Doctoring only
    the published range would leave the two pages disagreeing about the same
    quantity, which the derivation refuses rather than silently resolving.
    """
    source = tmp_path / "uncertainty.md"
    adm = tmp_path / "admissibility.json"
    monkeypatch.setattr(ir, "BLISS_UNCERTAINTY", source)
    monkeypatch.setattr(ir, "ADMISSIBILITY", adm)
    for minimum, lower, at_or_below, expected in (
        (0.953, 1.000, 6, "direction_robust_only_under_admissibility_magnitude_not"),
        (1.000, 1.100, 6, "direction_robust_only_under_admissibility_magnitude_not"),
        (1.050, 1.100, 0, "direction_robust_magnitude_not"),
    ):
        _published(source, minimum, lower)
        _admissibility(adm, at_or_below=at_or_below, admissible_corner=0)
        report = ir.build()
        bliss = next(h for h in report["headlines"] if h["key"] == "bliss_synergy")
        assert bliss["verdict"] == expected
        assert f"[{minimum:.3f}, 7.800]" in bliss["prior_predictive"]
        assert ("not established" in report["overall"]) == (minimum <= 1.0)
        if minimum == 1.0:
            assert "precision does not establish" in bliss["rationale"]
            assert "is not uniform" not in bliss["rationale"]
        assert "not equivalent to exceeding Bliss independence" in bliss["rationale"]


def test_missing_bliss_evidence_cannot_reuse_a_hardcoded_verdict(tmp_path, monkeypatch):
    source = tmp_path / "uncertainty.md"
    source.write_text("| full range (min, max) | [0.953, 7.800] |\n")
    monkeypatch.setattr(ir, "BLISS_UNCERTAINTY", source)
    with pytest.raises(ValueError, match="Missing or malformed Bliss uncertainty row"):
        ir.build()


def _form_artifact(tmp_path, monkeypatch, zero, full, reproduces=True):
    """A hypoxia form-sensitivity artifact with a chosen pair of endpoints."""
    import json
    src = tmp_path / "form.json"
    src.write_text(json.dumps({
        "gap_at_zero": zero,
        "gap_at_full": full,
        "sign_changes_across_form_axis": (zero > 0.0) and (full <= 0.0),
        "full_dependence_outside_published_interval": not (0.163 <= full <= 0.999),
        "reproduces_published_point": reproduces,
        "published": {"interval_95": [0.163, 0.999], "median": 0.960,
                      "default_point": 0.866},
    }))
    monkeypatch.setattr(ir, "HYPOXIA_FORM", src)
    return src


def test_the_hypoxia_verdict_tracks_the_measured_form_axis(tmp_path, monkeypatch):
    """BOTH sides. The §7.1 interval samples the 11 rate constants and never
    passes `--sdt-o2-dependence`, so it is computed at the O2 setting most
    favourable to a positive gap. A gap that collapses THROUGH zero along the
    form axis has no robust sign, and the row may not claim one."""
    # Sign survives the form axis -> the direction claim stands.
    _form_artifact(tmp_path, monkeypatch, 0.866065, 0.400000)
    report = ir.build()
    row = next(h for h in report["headlines"]
               if h["key"] == "hypoxia_kill_collapse")
    assert row["verdict"] == "direction_robust_magnitude_not"
    assert "survives BOTH" in row["rationale"]
    # The SUMMARY must agree with the row. A static sentence here called the
    # asymmetry direction-robust while the row withdrew it.
    # The two clauses were SPLIT in #926: bundled, one static half could
    # contradict the other's row. So the hypoxia clause is checked on its own,
    # and the penetration clause is checked against its own artifact in
    # test_the_penetration_verdict_tracks_the_measured_lambda_axis.
    assert "The hypoxia asymmetry is direction-robust but " \
        "magnitude-uncalibrated." in report["overall"]

    # Sign flips -> the direction claim must be withdrawn.
    _form_artifact(tmp_path, monkeypatch, 0.866065, -0.000097)
    report = ir.build()
    row = next(h for h in report["headlines"]
               if h["key"] == "hypoxia_kill_collapse")
    assert row["verdict"] == "direction_and_magnitude_uncertain", (
        "a gap measured to change sign across the O2 form axis was still "
        "published as direction-robust")
    assert "changes" in row["rationale"] and "sign" in row["rationale"], (
        "the rationale must say the sign changes, not merely that the magnitude "
        "is uncertain -- the magnitude collapse was already on record")
    assert "withdrawn" in row["rationale"]
    # The measured endpoints must REACH the prose, so the row cannot carry a
    # verdict derived from numbers it does not show.
    assert "-0.000097" in row["rationale"]
    assert "+0.866065" in row["prior_predictive"]
    # ... and the summary must withdraw it too, not just the row.
    assert "not direction-robust uncond" in report["overall"], (
        "the overall summary still calls the hypoxia asymmetry direction-robust "
        "while the row it summarises withdraws that claim")


def test_missing_hypoxia_form_evidence_cannot_reuse_a_hardcoded_verdict(
        tmp_path, monkeypatch):
    """The failure mode this derivation exists to remove: with the measurement
    absent the row must NOT fall back to the old `direction_robust` label."""
    import json
    src = tmp_path / "form.json"
    src.write_text(json.dumps({"gap_at_zero": 0.866}))
    monkeypatch.setattr(ir, "HYPOXIA_FORM", src)
    with pytest.raises(ValueError, match="missing"):
        ir.build()


def test_a_form_run_that_does_not_reproduce_the_published_point_is_refused(
        tmp_path, monkeypatch):
    """A form sweep that cannot reproduce the published default point is not
    comparable with the interval it qualifies, so it must not be used."""
    _form_artifact(tmp_path, monkeypatch, 0.5, -0.1, reproduces=False)
    with pytest.raises(ValueError, match="not comparable"):
        ir.build()


def test_report_is_deterministic_and_committed_json_is_fresh():
    """build() is deterministic, and the committed JSON must match a fresh build
    (so the committed report cannot silently go stale)."""
    import json

    r1 = ir.build()
    r2 = ir.build()
    assert r1 == r2
    committed = json.loads((REPO_ROOT / "analysis" / "identifiability-report.json").read_text())
    assert committed == r1, "analysis/identifiability-report.json is stale; re-run scripts/identifiability_report.py"


def _axis_artifact(tmp_path, monkeypatch, *, cancels=True, entailed=True,
                   reproduces=True, drop=None):
    """A penetration access-vs-depth artifact with chosen derived flags.

    Written as the generator writes it -- raw measurements plus a `derived`
    block -- because that is the contract `_penetration_verdict` reads."""
    import json
    tissues = ["Epithelial (well-vascularized)",
               "Epithelial (poorly-vascularized)",
               "Neuroectodermal (CNS/BBB)"]
    derived = {
        "lambda_cancels_at_the_wall": cancels,
        "lambda_ratio": 2.0,
        "headline_is_one_dose_response": reproduces,
        "dose_response_monotone": entailed,
        "ordering_entailed_not_measured": entailed,
        "access_order": tissues,
        "distance_order": [tissues[0], tissues[2], tissues[1]],
        "access_and_distance_agree": False,
        "deep_crossover_lambda_um": 50.98,
        "drug_attributable": {
            t: {"headline": h, "control": 0.014,
                "drug_attributable": h - 0.014,
                "fold_over_control": h / 0.014}
            for t, h in zip(tissues, (0.121, 0.026, 0.018))
        },
    }
    if drop:
        derived.pop(drop)
    src = tmp_path / "axis.json"
    src.write_text(json.dumps({
        "prior_ordering": {"n_draws": 200, "ordering_preserved": 200,
                           "cns_above_own_control": 173},
        "derived": derived,
    }))
    monkeypatch.setattr(ir, "PENETRATION_AXIS", src)
    return src


def test_the_penetration_verdict_tracks_the_measured_lambda_axis(
        tmp_path, monkeypatch):
    """When the observable is measured to be lambda-free, the row may not call
    the ordering a parameter-robustness result."""
    _axis_artifact(tmp_path, monkeypatch, cancels=True)
    v = ir._penetration_verdict()
    assert v["verdict"] == "direction_robust_by_construction_magnitude_not"
    assert "by construction" in v["rationale"]
    assert "vascular_permeability" in v["rationale"]
    # and the summary paragraph must agree with the row it summarises
    built = ir.build()
    row = [h for h in built["headlines"] if h["key"] == "penetration_gap"][0]
    assert row["verdict"] == v["verdict"]
    assert "by CONSTRUCTION" in built["overall"]
    # INDEPENDENCE: the hypoxia artifact here is the real (sign-flipping) one,
    # so both clauses must be present with their own wording at once. Bundled,
    # one branch decided the other's sentence.
    assert "not direction-robust uncond" in built["overall"], (
        "the hypoxia clause stopped tracking its own row when the penetration "
        "clause was split out")


def test_a_lambda_responsive_observable_gets_the_plain_verdict(
        tmp_path, monkeypatch):
    """The inverse direction: if the wall kill DID move with the penetration
    length, reading it as a penetration result would be supported, and the row
    must say so rather than keeping the stronger caveat."""
    _axis_artifact(tmp_path, monkeypatch, cancels=False)
    v = ir._penetration_verdict()
    assert v["verdict"] == "direction_robust_magnitude_not"
    assert "by construction" not in v["rationale"]
    built = ir.build()
    assert "by CONSTRUCTION" not in built["overall"]


def test_missing_penetration_axis_evidence_cannot_reuse_a_hardcoded_verdict(
        tmp_path, monkeypatch):
    """A missing field must RAISE. Falling back would silently republish the
    hardcoded `penetration-gradient ordering is parameter-robust` row that this
    derivation exists to remove."""
    _axis_artifact(tmp_path, monkeypatch, drop="lambda_cancels_at_the_wall")
    with pytest.raises(ValueError, match="lambda_cancels_at_the_wall"):
        ir._penetration_verdict()
    missing = tmp_path / "gone.json"
    monkeypatch.setattr(ir, "PENETRATION_AXIS", missing)
    with pytest.raises(OSError):
        ir._penetration_verdict()


def test_an_axis_run_that_does_not_reproduce_the_wall_kills_is_refused(
        tmp_path, monkeypatch):
    """If the single-cell binding cannot reproduce the binary's own wall kills,
    the run is not describing the observable the published interval used."""
    _axis_artifact(tmp_path, monkeypatch, reproduces=False)
    with pytest.raises(ValueError, match="does not reproduce"):
        ir._penetration_verdict()


def test_the_immune_row_claims_counts_not_per_cell_immunogenicity():
    """The manuscript is careful here and the row must not be looser than it.

    v1.md states: "The approximately 104:1 ratio compares total immune-kill
    counts, not DAMP potency or immunogenicity per dead cell." The row used to
    retain the direction as "SDT >> RSL3 immune PRIMING", which is exactly the
    per-cell reading that sentence declines -- SDT leads on counts because it
    kills far more cells (139,640 against 163 in the committed 2D run), and per
    ferroptotic death the comparison does not run its way.

    Pinned against the manuscript so the two cannot drift apart: this is the
    `one-artifact-describes-another` shape, where a derived label outruns the
    source it summarises.
    """
    row = [h for h in ir.build()["headlines"]
           if h["key"] == "immune_amplification_ratio"][0]
    rationale = row["rationale"]
    assert "immune-kill counts" in rationale, (
        "the row must say what the ratio compares")
    assert "immune priming) is claimed" not in rationale, (
        "the per-cell 'priming' reading is back; the manuscript declines it")

    ms = (REPO_ROOT / "article" / "drafts" / "v1.md").read_text(encoding="utf-8")
    assert ("compares total immune-kill counts, not DAMP potency or "
            "immunogenicity per dead cell") in ms, (
        "the manuscript sentence this row is pinned to has changed; re-check "
        "which reading the row is allowed to claim before editing this guard")


def _per_death_artifact(tmp_path, monkeypatch, *, ratio_max=1.1475,
                        refuted=True, direction=True, at_death=1.0018,
                        near=True, drop=None):
    """A per-death DAMP artifact with chosen derived values."""
    import json
    payload = {
        "ratio_min": 1.0469,
        "ratio_max": ratio_max,
        "readable_comparisons": 35,
        "direction_holds_everywhere": direction,
        "every_prediction_refuted": refuted,
        "at_death_ratio": at_death,
        "both_arms_near_threshold_at_death": near,
        "predictions": [
            {"low": 1.5, "high": 2.0, "source": "simulations/sim-icd/README.md"},
            {"low": 2.0, "high": 5.0, "source": "plans/simulation-suite-v2-plan.md"},
        ],
    }
    if drop:
        payload.pop(drop)
    src = tmp_path / "per-death.json"
    src.write_text(json.dumps(payload))
    monkeypatch.setattr(ir, "IMMUNE_PER_DEATH", src)
    return src


def test_the_immune_row_carries_the_measured_per_cell_magnitude(tmp_path, monkeypatch):
    """#927 corrected this row by READING the manuscript; #928 measured it.

    The clause must carry the measured magnitude and the documented bands it
    falls short of, so the row's refusal of the per-cell reading is evidenced
    rather than only careful.
    """
    _per_death_artifact(tmp_path, monkeypatch)
    row = [h for h in ir.build()["headlines"]
           if h["key"] == "immune_amplification_ratio"][0]
    assert "1.0469x to 1.1475x" in row["rationale"]
    assert "1.5-2x" in row["rationale"] and "2-5x" in row["rationale"]
    assert "no documented band is reached" in row["rationale"]
    assert "icd-damp-per-death.md" in row["source"]
    # AND THE VERDICT IS DELIBERATELY UNTOUCHED. This row reads
    # `directional_only` because the 104:1 moves ~25x with geometry on the
    # COUNT axis, which is true whatever the per-cell magnitude is. Deriving
    # the verdict from #928 would manufacture a dependency that does not exist.
    assert row["verdict"] == "directional_only"


def test_a_per_cell_magnitude_reaching_the_band_changes_the_clause(
        tmp_path, monkeypatch):
    """The clause is derived, so a measurement that reached the documented band
    would have to be reported as such rather than as a refutation."""
    _per_death_artifact(tmp_path, monkeypatch, ratio_max=2.4, refuted=False)
    row = [h for h in ir.build()["headlines"]
           if h["key"] == "immune_amplification_ratio"][0]
    assert "at least one documented band IS reached" in row["rationale"]
    assert "no documented band is reached" not in row["rationale"]
    assert row["verdict"] == "directional_only", (
        "the count-axis verdict must not move with the per-cell magnitude")


def test_a_per_cell_direction_that_inverted_would_be_reported(tmp_path, monkeypatch):
    """The row's wording accounts for a small positive difference. If the
    direction itself failed, the row would be understating the problem and the
    clause has to say so."""
    _per_death_artifact(tmp_path, monkeypatch, direction=False)
    row = [h for h in ir.build()["headlines"]
           if h["key"] == "immune_amplification_ratio"][0]
    assert "direction does NOT hold in every comparison" in row["rationale"]
    assert "only the magnitude is withdrawn" not in row["rationale"]


def test_missing_per_death_evidence_cannot_reuse_a_hardcoded_clause(
        tmp_path, monkeypatch):
    """A static clause quoting a measured number goes stale silently, which is
    how this report's own summary came to contradict the row above it."""
    monkeypatch.setattr(ir, "IMMUNE_PER_DEATH", tmp_path / "absent.json")
    with pytest.raises(OSError):
        ir._immune_per_cell_clause()
    _per_death_artifact(tmp_path, monkeypatch, drop="ratio_max")
    with pytest.raises(ValueError, match="missing 'ratio_max'"):
        ir._immune_per_cell_clause()


# --- #940: the Bliss row's direction clause, and the draws it rests on -------
#
# The row used to branch on the published sampled minimum ALONE, so its clause
# rested on a corner of the prior the model's own untreated-viability constraint
# rejects entirely, and nothing in the derivation could see that. These guards
# drive each part of the condition separately, because a verdict built from an
# AND of several conditions needs a fixture per condition -- each leaving the
# others true -- or removing one of them changes nothing a test can observe.


def test_the_bliss_row_is_conditional_on_the_committed_measurement():
    """On the real artifacts, the row must name the condition and not drop it."""
    bliss = next(h for h in ir.build()["headlines"] if h["key"] == "bliss_synergy")
    assert bliss["verdict"] == "direction_robust_only_under_admissibility_magnitude_not"
    r = bliss["rationale"]
    # The three restraints, each of which is why this is not a robustness result.
    assert "not a posterior" in r, "the conditioning is being read as inference"
    # NOT a bare "saturate" substring: the clause separately calls a saturated
    # arm the easiest way to fail the constraint, and "saturated" contains it, so
    # a sweep showed the bare check passing with the mechanism clause removed.
    assert "both arms saturate" in r, "the corner's cause is missing"
    assert "1 by arithmetic" in r, "the cause is stated without its consequence"
    assert "not more synergy" in r, "the median rise is being read as synergy"
    # And the measured figures, from the artifact rather than from this file.
    adm = ir._bliss_admissibility()
    assert f"{adm['admissible_min']:.3f}" in r
    assert f"{adm['n_admissible']} of {adm['n_full']}" in r
    assert str(adm["n_at_or_below_line"]) in r


def test_a_surviving_subadditive_draw_leaves_the_verdict_unqualified(tmp_path, monkeypatch):
    """THE FALSIFIER. If one sub-additive draw is admissible, the counterexample
    survives the model's own constraint and the row may not qualify itself.

    None is admissible as measured, so this case exists only as a fixture --
    which is exactly why it has to be driven rather than reasoned about.
    """
    source = tmp_path / "uncertainty.md"
    adm = tmp_path / "admissibility.json"
    monkeypatch.setattr(ir, "BLISS_UNCERTAINTY", source)
    monkeypatch.setattr(ir, "ADMISSIBILITY", adm)
    _published(source, 0.953, 1.000)

    _admissibility(adm, at_or_below=6, admissible_corner=0)
    wholly = ir.build()
    _admissibility(adm, at_or_below=6, admissible_corner=1, admissible_min=0.953)
    survives = ir.build()

    def row(r):
        return next(h for h in r["headlines"] if h["key"] == "bliss_synergy")
    assert row(wholly)["verdict"] == "direction_robust_only_under_admissibility_magnitude_not"
    assert row(survives)["verdict"] == "direction_and_magnitude_uncertain"
    # The ONLY thing that changed is the corner's admissibility: the published
    # range is byte-identical across the two runs, so the clause cannot be
    # coming from the published numbers.
    assert row(wholly)["prior_predictive"] == row(survives)["prior_predictive"]
    assert "inadmissible under the model" not in row(survives)["rationale"]


def test_a_page_that_does_not_reproduce_cannot_qualify_what_is_published(tmp_path, monkeypatch):
    """Conditioning a re-run would describe a re-run. The clause is withheld,
    and the row says WHY rather than silently reverting to the old sentence."""
    source = tmp_path / "uncertainty.md"
    adm = tmp_path / "admissibility.json"
    monkeypatch.setattr(ir, "BLISS_UNCERTAINTY", source)
    monkeypatch.setattr(ir, "ADMISSIBILITY", adm)
    _published(source, 0.953, 1.000)
    _admissibility(adm, at_or_below=6, admissible_corner=0, reproduces=False)
    bliss = next(h for h in ir.build()["headlines"] if h["key"] == "bliss_synergy")
    assert bliss["verdict"] == "direction_and_magnitude_uncertain"
    assert "does not reproduce the published range" in bliss["rationale"]
    assert "NOT reported here" in bliss["rationale"]
    assert "admissible minimum" not in bliss["rationale"]


def test_the_two_pages_must_agree_about_reaching_the_additive_null(tmp_path, monkeypatch):
    """An agreement being relied on is checked, not assumed.

    A published minimum above the line beside a page claiming to reproduce that
    minimum while counting draws below it is a contradiction, and resolving it
    silently would pick one of two numbers for the reader.
    """
    source = tmp_path / "uncertainty.md"
    adm = tmp_path / "admissibility.json"
    monkeypatch.setattr(ir, "BLISS_UNCERTAINTY", source)
    monkeypatch.setattr(ir, "ADMISSIBILITY", adm)
    _published(source, 1.050, 1.100)
    _admissibility(adm, at_or_below=6, admissible_corner=0)
    with pytest.raises(ValueError, match="disagree about whether any draw reaches"):
        ir.build()
    # The same disagreement is tolerated when the page does NOT claim to
    # reproduce the minimum, because then it is entitled to differ.
    _admissibility(adm, at_or_below=6, admissible_corner=0, reproduces=False)
    assert ir.build()["headlines"]


def test_a_missing_admissibility_artifact_raises_rather_than_republishing(tmp_path, monkeypatch):
    """Failing open here would quietly restore the unqualified row, which is the
    defect the derivation exists to remove."""
    source = tmp_path / "uncertainty.md"
    _published(source, 0.953, 1.000)
    monkeypatch.setattr(ir, "BLISS_UNCERTAINTY", source)
    monkeypatch.setattr(ir, "ADMISSIBILITY", tmp_path / "absent.json")
    with pytest.raises(ValueError, match="is missing; the Bliss row"):
        ir.build()


@pytest.mark.parametrize("drop", [
    ("decomposition", "n_at_or_below_line_admissible"),
    ("decomposition", "min_single_arm_kill_on_line"),
    ("decomposition", "denominator_ratio"),
    ("observable", "admissible_only"),
    ("observable", "n_admissible"),
])
def test_a_partial_admissibility_artifact_raises(tmp_path, monkeypatch, drop):
    """A clause built from a half-present measurement is worse than none."""
    source = tmp_path / "uncertainty.md"
    adm = tmp_path / "admissibility.json"
    _published(source, 0.953, 1.000)
    monkeypatch.setattr(ir, "BLISS_UNCERTAINTY", source)
    monkeypatch.setattr(ir, "ADMISSIBILITY", adm)
    _admissibility(adm, at_or_below=6, admissible_corner=0, drop=drop)
    with pytest.raises(ValueError, match="missing"):
        ir.build()


def test_the_minimum_reproduction_leg_is_read_and_not_assumed(tmp_path, monkeypatch):
    """The endpoint the claim rests on is the one that must be reproduced.

    A page reproducing the median and missing the minimum reproduces everything
    except the number in question, so the overall flag alone is not enough.
    """
    source = tmp_path / "uncertainty.md"
    adm = tmp_path / "admissibility.json"
    _published(source, 0.953, 1.000)
    monkeypatch.setattr(ir, "BLISS_UNCERTAINTY", source)
    monkeypatch.setattr(ir, "ADMISSIBILITY", adm)
    doc = json.loads(_admissibility(adm, at_or_below=6,
                                    admissible_corner=0).read_text())
    doc["observables"][0]["reproduction"]["min"]["reproduces"] = False
    adm.write_text(json.dumps(doc))
    bliss = next(h for h in ir.build()["headlines"] if h["key"] == "bliss_synergy")
    assert bliss["verdict"] == "direction_and_magnitude_uncertain", (
        "the row qualified itself on a page that did not reproduce the endpoint "
        "its own clause rests on")
    doc["observables"][0].pop("reproduction")
    adm.write_text(json.dumps(doc))
    with pytest.raises(ValueError, match="published sampled MINIMUM"):
        ir.build()


def test_the_overall_summary_cannot_contradict_the_bliss_row(tmp_path, monkeypatch):
    """#925's defect, one row over: a static summary beside a derived verdict.

    Drives every state the row can take and requires the paragraph to move with
    it, so a two-way branch left behind a three-way verdict fails here.
    """
    source = tmp_path / "uncertainty.md"
    adm = tmp_path / "admissibility.json"
    monkeypatch.setattr(ir, "BLISS_UNCERTAINTY", source)
    monkeypatch.setattr(ir, "ADMISSIBILITY", adm)
    seen = {}
    for minimum, lower, at_or_below, corner in ((0.953, 1.000, 6, 0),
                                                (0.953, 1.000, 6, 1),
                                                (1.050, 1.100, 0, 0)):
        _published(source, minimum, lower)
        _admissibility(adm, at_or_below=at_or_below, admissible_corner=corner,
                       admissible_min=0.953 if corner else 1.027)
        r = ir.build()
        bliss = next(h for h in r["headlines"] if h["key"] == "bliss_synergy")
        seen[bliss["verdict"]] = r["overall"]
    assert len(seen) == 3, f"only {len(seen)} verdict states reachable: {sorted(seen)}"
    assert len(set(seen.values())) == 3, (
        "two different Bliss verdicts produced the same overall paragraph, so "
        "the summary is not tracking the row")
    conditional = seen["direction_robust_only_under_admissibility_magnitude_not"]
    assert "inadmissible" in conditional and "conditioning adds no data" in conditional
    assert "inadmissible" not in seen["direction_and_magnitude_uncertain"]
    assert "not established" not in seen["direction_robust_magnitude_not"]


def test_the_conditional_state_does_not_move_the_data_conditioned_count(tmp_path, monkeypatch):
    """Conditioning removes draws. It does not add data, and the headline count
    of data-conditioned outputs must stay where it is."""
    source = tmp_path / "uncertainty.md"
    adm = tmp_path / "admissibility.json"
    monkeypatch.setattr(ir, "BLISS_UNCERTAINTY", source)
    monkeypatch.setattr(ir, "ADMISSIBILITY", adm)
    _published(source, 0.953, 1.000)
    _admissibility(adm, at_or_below=6, admissible_corner=0)
    r = ir.build()
    bliss = next(h for h in r["headlines"] if h["key"] == "bliss_synergy")
    assert bliss["verdict"] == "direction_robust_only_under_admissibility_magnitude_not"
    assert r["data_constrained_in_production"] == 0
    assert "this is not a data-conditioned claim" in bliss["rationale"]


def test_the_row_never_claims_unconditional_robustness_below_the_null():
    """THE REGISTERED RESTRAINT (#940 P3), as a property of the source.

    The unconditional label may be produced only by the published range clearing
    the null. An edit routing the admissibility leg into that branch instead of
    its own would read as this work having found robustness, which conditioning
    cannot establish.
    """
    import ast
    src = (REPO_ROOT / "scripts" / "identifiability_report.py").read_text()
    fn = next(n for n in ast.walk(ast.parse(src))
              if isinstance(n, ast.FunctionDef) and n.name == "_bliss_verdict")
    ret = next(n for n in ast.walk(fn) if isinstance(n, ast.Return))
    verdict = next(v for k, v in zip(ret.value.keys, ret.value.values)
                   if getattr(k, "value", None) == "verdict")
    names = {n.id for n in ast.walk(verdict) if isinstance(n, ast.Name)}
    assert "uniformly_supra_additive" in names and "conditional" in names, (
        "the verdict no longer consults both the published range and the "
        "admissibility condition")
    # The unconditional label and the conditional one must sit on DIFFERENT
    # tests, which an inline conditional expression makes checkable: the
    # unconditional string may not be reachable from the `conditional` branch.
    strings = [c.value for c in ast.walk(verdict)
               if isinstance(c, ast.Constant) and isinstance(c.value, str)]
    assert strings.count("direction_robust_magnitude_not") == 1
    assert strings.count("direction_robust_only_under_admissibility_magnitude_not") == 1
    assert isinstance(verdict, ast.IfExp)
    assert {n.id for n in ast.walk(verdict.test) if isinstance(n, ast.Name)} == {
        "uniformly_supra_additive"}, (
        "the unconditional robustness label is no longer gated on the published "
        "range alone")


def test_the_rationale_figures_are_read_not_typed(tmp_path, monkeypatch):
    """A literal that happens to EQUAL the true value is invisible to equality.

    That shape is on record in this repo, so every figure in the conditional
    clause is driven at a DOCTORED artifact holding values the real one does not
    have. An echo of the committed measurement then fails instead of passing.
    """
    source = tmp_path / "uncertainty.md"
    adm = tmp_path / "admissibility.json"
    monkeypatch.setattr(ir, "BLISS_UNCERTAINTY", source)
    monkeypatch.setattr(ir, "ADMISSIBILITY", adm)
    _published(source, 0.953, 1.000)
    doc = {
        "bliss_decomposition": {
            "n_at_or_below_line": 17,
            "n_at_or_below_line_admissible": 0,
            "n_on_line": 9,
            "min_single_arm_kill_on_line": 0.831,
            "no_synergy_line": 1.0,
            "numerator_ratio": 0.611,
            "denominator_ratio": 0.222,
            "corner_min_worst_untreated": 0.345,
        },
        "constraint": {"baseline_max": 0.05},
        "observables": [{
            "column": "bliss_synergy",
            "reproduces": True,
            "reproduction": {"min": {"reproduces": True}},
            "n_admissible": 44,
            "n_full": 191,
            "interval_quotable": True,
            "admissible_only": {"min": 1.413},
        }],
    }
    adm.write_text(json.dumps(doc))
    r = next(h for h in ir.build()["headlines"] if h["key"] == "bliss_synergy")["rationale"]
    real = json.loads(
        (REPO_ROOT / "analysis"
         / "headline-intervals-under-admissibility.json").read_text())
    live = next(o for o in real["observables"] if o["column"] == "bliss_synergy")
    for label, shown in (("corner size", "17"), ("on-line count", "9"),
                         ("saturated arm", "0.831"), ("admissible minimum", "1.413"),
                         ("admissible of full", "44 of 191"),
                         ("denominator ratio", "0.222"),
                         ("numerator ratio", "0.611"),
                         # the corner margin: 0.345 untreated against a 0.05
                         # limit is 6.9x, and both the share and the multiple
                         # must be read from the artifact
                         ("corner untreated share", "34.50%"),
                         ("corner margin", "6.9x")):
        assert shown in r, f"the {label} is not read from the artifact"
    # And the committed values must be ABSENT, so the clause cannot be carrying
    # both a read and a leftover literal.
    for stale in (f"{live['admissible_only']['min']:.3f}",
                  f"{live['n_admissible']} of {live['n_full']}"):
        assert stale not in r, (
            f"the committed figure {stale!r} survives in a clause built from a "
            "doctored artifact, so it is typed rather than read")


def test_a_withheld_interval_cannot_supply_the_bound_this_row_rests_on(tmp_path, monkeypatch):
    """The sibling page refuses a family's interval when too few draws survive,
    and a MINIMUM is the most tail-sensitive statistic there is.

    INERT ON REAL DATA, which is why it is driven here: the Bliss family clears
    that page's admissible-draw floor, asserted below rather than assumed, so
    nothing in the committed artifacts exercises this branch. An inert gate that
    only a fixture can reach is a recorded defect class in this repository.
    """
    live = next(o for o in json.loads(
        (REPO_ROOT / "analysis"
         / "headline-intervals-under-admissibility.json").read_text())["observables"]
        if o["column"] == "bliss_synergy")
    assert live["interval_quotable"] is True, (
        "the Bliss family's interval is now withheld, so this branch is no "
        "longer inert and the row's own clause must be re-read")

    source = tmp_path / "uncertainty.md"
    adm = tmp_path / "admissibility.json"
    monkeypatch.setattr(ir, "BLISS_UNCERTAINTY", source)
    monkeypatch.setattr(ir, "ADMISSIBILITY", adm)
    _published(source, 0.953, 1.000)
    doc = json.loads(_admissibility(adm, at_or_below=6,
                                    admissible_corner=0).read_text())
    doc["observables"][0]["interval_quotable"] = False
    adm.write_text(json.dumps(doc))
    bliss = next(h for h in ir.build()["headlines"] if h["key"] == "bliss_synergy")
    assert bliss["verdict"] == "direction_and_magnitude_uncertain", (
        "a direction claim was rested on a bound the sibling page withholds")
    assert "WITHHOLDS this family's interval" in bliss["rationale"]
    assert f"{doc['observables'][0]['admissible_only']['min']:.3f}" not in (
        bliss["rationale"]), "the withheld minimum is quoted anyway"
    # The corner fact is a COUNT of draws rather than a tail statistic, so it
    # survives the withholding and is still reported.
    assert "none of those draws is admissible" in bliss["rationale"]
    doc["observables"][0].pop("interval_quotable")
    adm.write_text(json.dumps(doc))
    with pytest.raises(ValueError, match="interval_quotable"):
        ir.build()


def test_without_the_corner_margin_the_row_says_so_rather_than_implying_it(tmp_path, monkeypatch):
    """The saturation clause accounts for the draws ON the line only.

    If the sibling page stops reporting how far the rest of the corner fails,
    the row must say that only those draws are accounted for -- not keep a
    sentence that reads as covering the whole corner.
    """
    source = tmp_path / "uncertainty.md"
    adm = tmp_path / "admissibility.json"
    monkeypatch.setattr(ir, "BLISS_UNCERTAINTY", source)
    monkeypatch.setattr(ir, "ADMISSIBILITY", adm)
    _published(source, 0.953, 1.000)
    doc = json.loads(_admissibility(adm, at_or_below=6,
                                    admissible_corner=0).read_text())
    doc["bliss_decomposition"].pop("corner_min_worst_untreated")
    adm.write_text(json.dumps(doc))
    bliss = next(h for h in ir.build()["headlines"] if h["key"] == "bliss_synergy")
    r = bliss["rationale"]
    assert bliss["verdict"] == "direction_robust_only_under_admissibility_magnitude_not"
    assert "only the draws on the line are accounted for" in r
    assert "x the limit" not in r, "a margin is quoted with nothing to compute it from"
    # The same when the constraint itself is unavailable, since a multiple of a
    # missing limit is not a number.
    doc["bliss_decomposition"]["corner_min_worst_untreated"] = 0.2975
    doc.pop("constraint")
    adm.write_text(json.dumps(doc))
    r2 = next(h for h in ir.build()["headlines"]
              if h["key"] == "bliss_synergy")["rationale"]
    assert "only the draws on the line are accounted for" in r2
