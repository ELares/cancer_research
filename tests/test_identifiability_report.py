"""Tests for the practical-identifiability synthesis (#503).

Pure-Python stdlib (no compiled extension), so these run in CI. They guard the
load-bearing facts (degrees of freedom, the non-identifiable count, that no
headline is point-estimable, that zero headlines are data-conditioned in the
production regime) and that the committed report is not stale.
"""

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
               "direction_robust_by_construction_magnitude_not"}
    for h in r["headlines"]:
        assert h["verdict"] in allowed, f"{h['key']} claims {h['verdict']}"
    # At least the single-cell kill rate and the immune ratio are directional-only.
    by_key = {h["key"]: h for h in r["headlines"]}
    assert by_key["single_cell_kill_rate"]["verdict"] == "directional_only"
    assert by_key["immune_amplification_ratio"]["verdict"] == "directional_only"


def test_bliss_verdict_tracks_the_sampled_range_not_a_rounded_null(tmp_path, monkeypatch):
    """Exercise both sides of the inference, including an exactly additive draw."""
    source = tmp_path / "uncertainty.md"
    monkeypatch.setattr(ir, "BLISS_UNCERTAINTY", source)
    for minimum, lower, expected in (
        (0.953, 1.000, "direction_and_magnitude_uncertain"),
        (1.000, 1.100, "direction_and_magnitude_uncertain"),
        (1.050, 1.100, "direction_robust_magnitude_not"),
    ):
        source.write_text(
            "| default point estimate (no perturbation) | 1.992 |\n"
            "| prior-predictive median | 1.348 |\n"
            f"| 95% prior-predictive interval | [{lower:.3f}, 5.242] |\n"
            f"| full range (min, max) | [{minimum:.3f}, 7.800] |\n"
        )
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
