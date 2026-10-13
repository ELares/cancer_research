"""Practical-identifiability synthesis for the headline simulation outputs (#503).

The 2026 fresh-eyes review noted that the headline numbers are dimensionless and
many parameters are practically non-identifiable, so the point estimates carry
little information beyond "a plausible mechanism CAN produce a differential."
This script makes that accounting explicit and reproducible: for each HEADLINE
output it consolidates which parameters drive it, which are identifiable, how
wide the prior-predictive interval is, whether it is data-conditioned, and the
resulting verdict (point-estimable vs directional-only).

It does NOT re-run the sensitivity analyses (Sobol/Morris/ABC need the compiled
extension). It SYNTHESIZES the already-committed analyses, citing each source,
and cross-checks the structural facts it can against the machine-readable
`analysis/prcc-results.json` (the 11 swept parameters and their ranges, so the
degrees-of-freedom count and the named non-identifiable parameters are validated
against the real parameter set rather than hand-asserted).

Sources synthesized:
- analysis/prcc-results.json (PRCC #134): the swept parameter set + ranges.
- analysis/sobol-sensitivity-report.md (Sobol #331): total-effect ST per
  parameter; identifiable iff ST >= 0.05.
- analysis/uncertainty-intervals-report.md (#332): prior-predictive death-rate
  intervals.
- analysis/headline-sensitivity-report.md / headline-uncertainty-*.md (#331/#332
  spatial): per-headline Morris drivers + prior-predictive headline intervals.
- analysis/calibration/abc-posterior-report.md (#332): the in-vitro / in-vivo
  prior disjunction.

Run `python scripts/identifiability_report.py` to (re)write the committed report.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
PRCC = REPO_ROOT / "analysis" / "prcc-results.json"
BLISS_UNCERTAINTY = REPO_ROOT / "analysis" / "headline-uncertainty-report.md"
HYPOXIA_FORM = (REPO_ROOT / "analysis" / "calibration"
                / "hypoxia-gap-form-sensitivity.json")
PENETRATION_AXIS = (REPO_ROOT / "analysis" / "calibration"
                    / "penetration-access-vs-depth.json")
OUT_MD = REPO_ROOT / "analysis" / "identifiability-report.md"
OUT_JSON = REPO_ROOT / "analysis" / "identifiability-report.json"

# Parameters the single-cell Sobol screen found NON-identifiable from the kill
# rate (ST < 0.05). Verbatim from analysis/sobol-sensitivity-report.md (#331).
SOBOL_NON_IDENTIFIABLE = [
    "gsh_scav_efficiency",
    "nrf2_gsh_rate",
    "fsp1_rate",
    "fenton_rate",
    "death_threshold",
    "gpx4_degradation_by_ros",
]
# The three that dominate (ST), from the same report.
SOBOL_DOMINANT = {"lp_propagation": 0.504, "gpx4_rate": 0.285, "lp_rate": 0.177}

# Per-headline identifiability synthesis. Numbers are quoted from the committed
# reports (cited in `source`); the structural facts are cross-checked below.
HEADLINES = [
    {
        "key": "single_cell_kill_rate",
        "observable": "single-cell death rate (e.g. Persister x RSL3), the Figure 7 numbers",
        "drivers": ["lp_propagation", "gpx4_rate", "lp_rate"],
        "non_identifiable_params": list(SOBOL_NON_IDENTIFIABLE),
        "prior_predictive": "Persister x RSL3 point 42.5%, but 95% prior-predictive [1.6%, 99.7%] (width 98.1%); PersisterNrf2 x RSL3 point 0.0%, interval [0.0%, 37.8%]",
        "data_conditioned": "in-vitro only (single-inducer ABC posterior, #332); the in-vivo priors that produce the Figure 7 numbers are numerically DISJOINT from that in-vitro fit. Read that carefully: those priors are sensitivity ranges around the defaults themselves, with both LP ranges spanning 0.5–2 times default (scripts/run_prcc.py), so the disjunction says the in-vitro fit falls outside the defaults' own neighbourhood -- it restates the falsification rather than supplying independent grounds to discount it (analysis/calibration/in-vivo-prior-provenance.md). The current joint run is assessed separately below.",
        "verdict": "directional_only",
        "rationale": "the point estimate is essentially uninformative under the documented parameter uncertainty (the interval nearly spans [0,1]); the robust claim is that the differential between phenotypes exists, not its magnitude",
        "source": "sobol-sensitivity-report.md (#331); uncertainty-intervals-report.md (#332); abc-posterior-report.md (#332); joint-posterior.md (#500); abc-information-content.md; headline-at-fitted-cascade.md",
    },
    {
        "key": "bliss_synergy",
        "observable": "RSL3 + FSP1i Bliss synergy_score (the ~1.99x)",
        "drivers": ["lp_propagation", "gsh_scav_efficiency", "gpx4_rate"],
        "non_identifiable_params": ["sdt_ros (structural zero, no SDT in the pair)", "rsl3_gpx4_inhib (structural zero, fixed DrugEffect)"],
        "data_conditioned": "no (prior-predictive; the combo fit is not data-conditioned)",
        "source": "headline-sensitivity-report.md (#331); headline-uncertainty-report.md (#332)",
    },
    {
        "key": "hypoxia_kill_collapse",
        "observable": "SDT-minus-RSL3 hypoxic-zone kill gap (the kill-collapse asymmetry)",
        "drivers": ["sdt_ros", "lp_propagation", "lp_rate", "gsh_scav_efficiency"],
        "non_identifiable_params": [],
        "prior_predictive": "gap stays POSITIVE across its 95% interval ~[0.16, 1.00] (median 0.96, point 0.87) under the O2-independent assumption",
        "data_conditioned": "no (prior-predictive); additionally ASSUMPTION-bracketed: the magnitude collapses from ~87.8% to ~0% under full SDT O2-dependence (#336/#358), which the lead clinical agent SONALA-001 occupies",
        "verdict": "direction_robust_magnitude_not",
        "rationale": "the asymmetry SIGN is parameter-robust, but the magnitude is both wide (interval) and assumption-dependent (the contested SDT O2-dependence), so it is an assumption restated quantitatively, not a calibrated prediction",
        "source": "headline-sensitivity-report.md (#331); headline-uncertainty-tme-report.md (#332); manuscript Section 7.1",
    },
    {
        "key": "penetration_gap",
        "observable": "per-tissue vessel-wall RSL3 kill (the 40% -> 1.8% behind the BBB)",
        "drivers": ["lp_propagation", "lp_rate", "gpx4_rate"],
        "non_identifiable_params": [],
        "prior_predictive": "very wide per-tissue intervals (e.g. well-vascularized median 0.23 ~[0.00, 0.93]; CNS/BBB median 0.04 ~[0.00, 0.77]) because the bistable switch dominates; BUT the within-draw across-tissue ORDERING (well >= poorly >= CNS) held in 300/300 draws",
        "data_conditioned": "no (prior-predictive; transport params at fixed uncalibrated presets, their own ranges not swept)",
        "verdict": "direction_robust_magnitude_not",
        "rationale": "the penetration-gradient ordering is parameter-robust (per-draw, not inferred from the overlapping marginals); the absolute per-tissue magnitudes are not point-estimable",
        "source": "headline-uncertainty-penetration-report.md (#332)",
    },
    {
        "key": "immune_amplification_ratio",
        "observable": "SDT:RSL3 immune-kill ratio (the 104:1)",
        "drivers": ["lp_propagation", "lp_rate", "sdt_ros"],
        "non_identifiable_params": [],
        "prior_predictive": "the 104:1 (2D, near DAMP saturation) falls to ~4:1 in 3D geometry; the separate SDT immune kill fraction using final non-ferroptotic counts spans ~[0.009, 0.171], low-but-positive within this prior-predictive interval. The listed Morris drivers describe that final-population fraction, not the SDT:RSL3 count ratio. The fraction does not measure eligibility at earlier immune windows or isolate per-cell DAMP amplification",
        "data_conditioned": "no; geometry-dependent (the 2D-vs-3D shrink is a structural, not parametric, effect)",
        "verdict": "directional_only",
        "rationale": "the ratio is presented as a directional ceiling, not a number: it changes ~25x with geometry alone, which a parametric analysis cannot capture, so only the direction (SDT >> RSL3 immune priming) is claimed",
        "source": "headline-sensitivity-report.md (#331); CALIBRATION_STATUS.md immune row",
    },
]


def _bliss_verdict() -> dict:
    """Read the reported precision and sampled range, not a rounded headline.

    A ratio of one is Bliss independence, not synergy. A positive combination
    benefit over either single agent is also a different comparison.
    """
    source = BLISS_UNCERTAINTY.read_text(encoding="utf-8")
    rows = {}
    for line in source.splitlines():
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if line.startswith("|") and len(cells) == 2:
            rows[cells[0]] = cells[1]

    def numbers(label, count):
        value = rows.get(label, "")
        pattern = r"\d+(?:\.\d+)?"
        expected = pattern if count == 1 else rf"\[({pattern}), ({pattern})\]"
        if not re.fullmatch(expected, value):
            raise ValueError(f"Missing or malformed Bliss uncertainty row: {label}")
        return [float(v) for v in re.findall(pattern, value)]

    point, = numbers("default point estimate (no perturbation)", 1)
    median, = numbers("prior-predictive median", 1)
    lower, upper = numbers("95% prior-predictive interval", 2)
    minimum, maximum = numbers("full range (min, max)", 2)
    if not minimum <= lower <= median <= upper <= maximum:
        raise ValueError("Bliss uncertainty interval and sampled range are inconsistent")
    uniformly_supra_additive = minimum > 1.0
    if uniformly_supra_additive:
        direction = (
            "every sampled ratio exceeds the additive null of 1 at the reported precision; "
            "this supports supra-additivity in the sampled ensemble, not a guarantee for "
            "unsampled parameters or biological experiments. "
        )
    elif minimum < 1.0:
        direction = (
            "the reported sampled range includes sub-additive draws, so "
            "supra-additivity is not uniform across the sampled ensemble. "
        )
    else:
        direction = (
            "the sampled minimum is reported as the additive null of 1.000; "
            "that precision does not establish whether every sampled ratio exceeds the null. "
        )
    return {
        "prior_predictive": (
            f"point {point:.3f}x, 95% prior-predictive [{lower:.3f}, {upper:.3f}], "
            f"median {median:.3f}x; full sampled range [{minimum:.3f}, {maximum:.3f}]; "
            "nonlinearity and/or interactions (Morris sigma > mu* for every active parameter)"
        ),
        "verdict": ("direction_robust_magnitude_not" if uniformly_supra_additive
                    else "direction_and_magnitude_uncertain"),
        "rationale": (
            direction + "A bound rounded to 1.000 cannot establish strict exclusion of the null; "
            "outperforming a single agent is not equivalent to exceeding Bliss independence. "
            "The magnitude is not point-estimable or data-conditioned"
        ),
    }


def _hypoxia_verdict() -> dict:
    """Derive the hypoxia row from the MEASURED form axis, not from a label.

    The §7.1 interval (#332) samples the 11 rate constants and never passes
    `--sdt-o2-dependence`, so it is computed entirely at the O2 setting most
    favourable to a positive gap. This row used to be hardcoded
    `direction_robust_magnitude_not` with the rationale "the asymmetry SIGN is
    parameter-robust, but the magnitude is ... assumption-dependent".

    The magnitude collapse was already recorded. The SIGN at the collapsed end
    was not, and a difference that collapses THROUGH zero does not have a robust
    sign to report. So the verdict is computed from the committed measurement:
    if the gap changes sign across the form axis, the direction is not robust
    and the row may not say that it is.
    """
    d = json.loads(HYPOXIA_FORM.read_text(encoding="utf-8"))
    for key in ("gap_at_zero", "gap_at_full", "sign_changes_across_form_axis",
                "reproduces_published_point", "published"):
        if key not in d:
            raise ValueError(
                f"hypoxia form-sensitivity artifact is missing {key!r}; the "
                "verdict is derived from it and must not fall back to a label")
    if not d["reproduces_published_point"]:
        raise ValueError(
            "the form-sensitivity run does not reproduce the published default "
            "point, so it is not comparable with the interval it qualifies")
    lo, hi = d["published"]["interval_95"]
    zero, full = float(d["gap_at_zero"]), float(d["gap_at_full"])
    flips = bool(d["sign_changes_across_form_axis"])
    outside = bool(d.get("full_dependence_outside_published_interval"))
    spread = (
        f"gap stays POSITIVE across its 95% interval ~[{lo}, {hi}] (median "
        f"{d['published']['median']}, point {d['published']['default_point']}) "
        "under the O2-independent assumption -- an interval over the 11 RATE "
        "CONSTANTS at a FIXED O2 form, which the sweep never varies. Measured "
        f"along the form axis instead, the gap runs {zero:+.6f} at zero "
        f"dependence to {full:+.6f} at full"
        + f", {'outside' if outside else 'inside'} that interval."
    )
    if flips:
        return {
            "prior_predictive": spread,
            "verdict": "direction_and_magnitude_uncertain",
            "rationale": (
                "the asymmetry SIGN is robust across the 11 rate constants but "
                "NOT across the O2 response form: measured, the gap changes "
                f"sign ({zero:+.6f} -> {full:+.6f}), and full O2-dependence is "
                "the regime the lead clinical agent SONALA-001 occupies. A "
                "difference that collapses through zero has no robust sign to "
                "report, so the earlier `direction_robust` reading was scoped "
                "to the rate constants at a fixed form and is withdrawn as an "
                "unconditional claim. The magnitude is not point-estimable "
                "either"
            ),
        }
    return {
        "prior_predictive": spread,
        "verdict": "direction_robust_magnitude_not",
        "rationale": (
            "the asymmetry sign survives BOTH the rate constants and the "
            f"measured O2 form axis ({zero:+.6f} -> {full:+.6f}), while the "
            "magnitude is wide and uncalibrated"
        ),
    }


def _penetration_verdict() -> dict:
    """Derive the penetration row from what the observable actually contains.

    This row used to be hardcoded `direction_robust_magnitude_not` with the
    rationale "the penetration-gradient ordering is parameter-robust (per-draw,
    not inferred from the overlapping marginals)". That parenthetical defends
    against the wrong objection. The ordering is not in doubt because the
    marginals overlap; it is in doubt as EVIDENCE, for two measured reasons
    (#926, `penetration-access-vs-depth.md`):

      - The observable is `vessel_wall_death_rate`, read at `r = 0`, where
        `exp(-r/lambda)` is 1 for every lambda. Measured: the two drugs
        `sim-tissue-pk` ships differ two-fold in penetration length and return
        BIT-IDENTICAL wall kills. So the quantity is a `vascular_permeability`
        gradient, and lambda -- the one transport parameter #335 anchored to
        measured data -- does not enter it.
      - All three tissues run identical cells with identical RNG at bin 0, so
        the three numbers are ONE monotone dose-response at three fixed
        decreasing concentrations. A paired comparison of that cannot invert, so
        the draw count is entailed rather than informative.

    The claim is not false and this does not say it is. It is robust by
    construction, and it is about access rather than penetration.
    """
    d = json.loads(PENETRATION_AXIS.read_text(encoding="utf-8"))
    v = d.get("derived")
    if not isinstance(v, dict):
        raise ValueError(
            "penetration axis artifact has no `derived` block; the verdict is "
            "derived from it and must not fall back to a label")
    for key in ("lambda_cancels_at_the_wall", "lambda_ratio",
                "headline_is_one_dose_response", "dose_response_monotone",
                "ordering_entailed_not_measured", "access_and_distance_agree",
                "deep_crossover_lambda_um", "drug_attributable"):
        if key not in v:
            raise ValueError(
                f"penetration axis artifact is missing {key!r}; the verdict is "
                "derived from it and must not fall back to a label")
    if not v["headline_is_one_dose_response"]:
        raise ValueError(
            "the penetration run does not reproduce the binary's own wall "
            "kills from the single-cell binding, so it is not describing the "
            "observable the published interval used")

    cns = v["access_order"][-1]
    cns_add = float(v["drug_attributable"][cns]["drug_attributable"])
    ctrl = float(v["drug_attributable"][cns]["control"])
    prior = d["prior_ordering"]
    spread = (
        "very wide per-tissue intervals (e.g. well-vascularized median 0.23 "
        "~[0.00, 0.93]; CNS/BBB median 0.04 ~[0.00, 0.77]) because the bistable "
        "switch dominates; the within-draw across-tissue ORDERING "
        f"(well >= poorly >= CNS) held in {prior['ordering_preserved']}/"
        f"{prior['n_draws']} draws here and 300/300 as published -- but that is "
        "ENTAILED, not measured: one monotone dose-response at three fixed "
        "decreasing concentrations with a shared RNG stream cannot invert. "
        f"Measured against its own no-drug control ({ctrl:.3f}), the CNS figure "
        f"is {cns_add:+.3f} of drug-attributable kill"
    )

    if not v["lambda_cancels_at_the_wall"]:
        return {
            "prior_predictive": spread,
            "source": "headline-uncertainty-penetration-report.md (#332); calibration/penetration-access-vs-depth.md (#926)",
            "verdict": "direction_robust_magnitude_not",
            "rationale": (
                "the ordering is parameter-robust and the observable does "
                "respond to the penetration length, so reading it as a "
                "penetration result is supported; magnitudes are not "
                "point-estimable"
            ),
        }
    cross = v["deep_crossover_lambda_um"]
    disagree = not v["access_and_distance_agree"]
    return {
        "prior_predictive": spread,
        "source": "headline-uncertainty-penetration-report.md (#332); calibration/penetration-access-vs-depth.md (#926)",
        "verdict": "direction_robust_by_construction_magnitude_not",
        "rationale": (
            "the ordering IS robust to the 11 rate constants, but by "
            "construction rather than by measurement, and it is not a "
            "penetration result: the observable sits at r = 0 where "
            f"exp(-r/lambda) = 1, and two drugs differing {v['lambda_ratio']:.1f}"
            "-fold in penetration length return bit-identical wall kills, so "
            "what is ordered is the uncalibrated `vascular_permeability` and "
            "not the lambda #335 anchored to measured data"
            + (f". Where lambda does act -- each tissue's own deepest cell -- "
               f"the poorly/CNS ordering inverts below lambda ~ {cross:.0f} um, "
               "because access and distance rank these three tissues "
               "differently; so the direction is robust at the wall and not at "
               "depth" if disagree else "")
            + ". The absolute per-tissue magnitudes are not point-estimable"
        ),
    }


def build() -> dict:
    prcc = json.loads(PRCC.read_text(encoding="utf-8"))
    params = prcc["metadata"]["parameters"]
    ranges = prcc["metadata"]["parameter_ranges"]
    dof = len(params)
    # Cross-check: every named single-cell non-identifiable parameter is a real
    # swept PRCC parameter (so the list cannot drift away from the parameter set).
    for p in SOBOL_NON_IDENTIFIABLE:
        assert p in params, f"non-identifiable param {p!r} not in the PRCC parameter set"
    for p in SOBOL_DOMINANT:
        assert p in params, f"dominant param {p!r} not in the PRCC parameter set"
    bliss = _bliss_verdict()
    hypoxia = _hypoxia_verdict()
    penetration = _penetration_verdict()
    derived = {"bliss_synergy": bliss, "hypoxia_kill_collapse": hypoxia,
               "penetration_gap": penetration}
    # A key that matches NOTHING would silently leave the hardcoded row in
    # place, which is the whole defect this derivation exists to remove.
    known = {h["key"] for h in HEADLINES}
    missing = sorted(set(derived) - known)
    if missing:
        raise ValueError(
            f"derived verdict keys {missing} match no headline; the hardcoded "
            "row would have been published instead")
    headlines = [{**h, **derived[h["key"]]} if h["key"] in derived else dict(h)
                 for h in HEADLINES]
    bliss_summary = (
        "Bliss supra-additivity holds across the sampled ensemble, with an uncalibrated magnitude. "
        if bliss["verdict"] == "direction_robust_magnitude_not" else
        "Uniform Bliss supra-additivity is not established by the reported ensemble, and its magnitude is uncalibrated. "
    )
    # DERIVED, because a static sentence here contradicted the row it
    # summarises the moment the hypoxia verdict started tracking the measured
    # form axis: the row said `direction_and_magnitude_uncertain` while this
    # line still called the asymmetry direction-robust.
    # DERIVED, and split from the penetration clause it used to share a
    # sentence with. Bundling them meant one static half could contradict the
    # other's row: the sentence called the penetration ordering direction-robust
    # while the penetration verdict started qualifying that as true only by
    # construction. Each clause now tracks its own row.
    hypoxia_summary = (
        "The hypoxia asymmetry is direction-robust only at a FIXED O2 form -- "
        "measured along that axis its sign changes, so it is not direction-robust "
        "unconditionally. "
        if hypoxia["verdict"] != "direction_robust_magnitude_not" else
        "The hypoxia asymmetry is direction-robust but magnitude-uncalibrated. "
    )
    penetration_summary = (
        "The penetration ordering holds in every draw, but by CONSTRUCTION -- "
        "one monotone dose-response at three fixed decreasing concentrations -- "
        "and it orders vascular access, not penetration: the observable sits "
        "where the penetration length cancels. "
        if penetration["verdict"] == "direction_robust_by_construction_magnitude_not"
        else
        "The penetration ordering is direction-robust but "
        "magnitude-uncalibrated. "
    )
    joint_information = _joint_information()
    if joint_information.get("unassessable"):
        headlines[0]["data_conditioned"] += (
            " Joint inference unavailable: " + joint_information["unassessable"]
            + "; no joint parameter-information or interval claim is supported.")
    return {
        "degrees_of_freedom": dof,
        "swept_parameters": params,
        "parameter_ranges": ranges,
        "single_cell_sobol": {
            "dominant_ST": SOBOL_DOMINANT,
            "non_identifiable_count": len(SOBOL_NON_IDENTIFIABLE),
            "non_identifiable": SOBOL_NON_IDENTIFIABLE,
            "note": "6 of 11 swept parameters are practically non-identifiable from the kill rate (ST < 0.05)",
        },
        "data_constrained_in_production": 0,
        "data_constrained_note": (
            "The production simulation matrix uses fixed in-vivo defaults. The "
            "in-vitro single-cell fits do not condition the spatial headline outputs. "
            "The current joint run's inferential status is reported separately; "
            "historical joint intervals are not substituted for it."
        ),
        "joint_posterior_information": joint_information,
        "headlines": headlines,
        "overall": (
            "No headline output is fully point-estimable. The single-cell kill rate "
            "and the immune ratio are directional-only. " + bliss_summary +
            hypoxia_summary +
            penetration_summary +
            "With 11 free rate constants, 6 non-identifiable "
            "from the kill rate, and 0 of the headlines data-conditioned in the "
            "production regime, the honest reading of every reported magnitude is "
            "order-of-magnitude / directional, exactly as the manuscript labels them."
        ),
    }


def _joint_information():
    """Assess the current source, never a possibly stale derived information file."""
    from abc_posterior_information import assess

    path = REPO_ROOT / "analysis" / "calibration" / "joint-posterior.json"
    try:
        return assess(path)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        return {"artifact": "analysis/calibration/joint-posterior.json",
                "unassessable": f"current joint artifact unavailable or invalid ({type(exc).__name__})"}


def _fitted_cascade_facts():
    """Numbers for the closing section, READ FROM THE ARTIFACTS.

    These were hand-typed beside two files that compute them. A figure written
    into prose next to the artifact that produces it is the failure this
    repository keeps rediscovering: the artifact moves, the sentence does not,
    and the stale number reads as freshly checked. Missing artifacts produce an
    explicit "not available" rather than a number that might be from last month.
    """
    root = REPO_ROOT
    out = {}
    try:
        h = json.loads((root / "analysis" / "headline-at-fitted-cascade.json").read_text())
        worst = max(("ctrpv2_point", "posterior_median"),
                    key=lambda k: h[k]["admissibility"]["worst_rate"])
        a = h[worst]["admissibility"]
        out["baseline_default"] = h["default"]["admissibility"]["worst_rate"]
        out["baseline_fitted"] = a["worst_rate"]
        out["baseline_condition"] = a["worst_condition"]
        out["baseline_constraint"] = a["constraint"]
        out["bliss_fitted"] = h[worst]["bliss"]
    except (OSError, ValueError, KeyError):
        pass
    joint = _joint_information()
    out["joint_information"] = joint
    if joint.get("parameters"):
        params = joint["parameters"]
        out["informed"] = sum(1 for v in params.values() if v["informed"])
        out["total_params"] = len(params)
        out["uninformed"] = [k for k, v in params.items() if not v["informed"]]
    return out


def write_report(r: dict) -> None:
    def hsec(h: dict) -> str:
        return (
            f"### {h['observable']}\n\n"
            f"- **Drivers:** {', '.join(h['drivers'])}\n"
            f"- **Non-identifiable:** {', '.join(h['non_identifiable_params']) or 'none flagged for this headline'}\n"
            f"- **Prior-predictive spread:** {h['prior_predictive']}\n"
            f"- **Data-conditioned:** {h['data_conditioned']}\n"
            f"- **Verdict:** `{h['verdict']}` ({h['rationale']})\n"
            f"- **Source:** {h['source']}\n"
        )

    f = _fitted_cascade_facts()
    if "informed" in f:
        un = ", ".join(f"`{u}`" for u in f["uninformed"]) or "none"
        informed_clause = (
            f"{f['informed']} of its {f['total_params']} parameters are informed by "
            f"the marginal-width criterion. No marginal contraction was detected "
            f"for {un}; this does not establish equality with the prior")
    else:
        informed_clause = (
            "the current joint run is NOT ASSESSED: "
            + f["joint_information"].get("unassessable", "no parameters parsed")
            + ". No joint credible-interval or parameter-contraction claim is made")
    if "baseline_fitted" in f:
        baseline_clause = (
            f"untreated {f['baseline_condition']} death goes from "
            f"{f['baseline_default']*100:.2f}% to {f['baseline_fitted']*100:.2f}%, "
            f"against the model's own stated constraint of under "
            f"{f['baseline_constraint']*100:.0f}%.")
        bliss_clause = (
            f"the Bliss ratio is exactly {f['bliss_fitted']:.1f} because both single "
            "arms saturate")
    else:
        baseline_clause = ("the comparison is not available "
                           "(headline-at-fitted-cascade.json missing).")
        bliss_clause = "the headline table is not available"
    md = f"""# Practical-identifiability of the headline outputs (#503)

A consolidated, reproducible accounting of which headline simulation results are
point-estimable and which are directional-only, synthesizing the committed
sensitivity and uncertainty analyses (it does not re-run them; the Sobol/Morris/
ABC steps need the compiled extension). The structural facts (the {r['degrees_of_freedom']}
swept parameters and the named non-identifiable set) are cross-checked against
`analysis/prcc-results.json`.

## Headline accounting

The degrees of freedom: **{r['degrees_of_freedom']} free rate constants** are swept
(`{', '.join(r['swept_parameters'])}`). Of these, **{r['single_cell_sobol']['non_identifiable_count']}
are practically non-identifiable from the kill rate** (Sobol total-effect ST < 0.05:
{', '.join(r['single_cell_sobol']['non_identifiable'])}); three dominate
({', '.join(f'{k} ST={v}' for k, v in r['single_cell_sobol']['dominant_ST'].items())}).

**Data-constrained in the production regime: {r['data_constrained_in_production']}.**
{r['data_constrained_note']}

## Per-headline verdicts

{''.join(hsec(h) for h in r['headlines'])}

## Overall

{r['overall']}

## What would make a headline point-estimable

A headline becomes point-estimable when (1) its driving parameters are identified
(narrowed) by data in the regime that produces it, and (2) the prior-predictive
interval collapses to a usable width.

The current joint fit and the historical substitution experiment answer
different questions. Neither makes a production headline point-estimable.

- **Current joint inference:** {informed_clause}
  (`analysis/calibration/abc-information-content.md`).
- **Historical substitution experiment:** carrying the recorded historical
  in-vitro parameter vectors into the spatial model was INADMISSIBLE:
  {baseline_clause} Under those historical vectors, {bliss_clause}
  (`analysis/headline-at-fitted-cascade.md`).

That historical test is not an admissibility result for the corrected joint
run. An underpowered run requires enough accepted draws under its unchanged
criterion before its intervals or information content can be interpreted.
Independent measurements matched to the model's observables and experimental
regime remain necessary to condition the spatial headlines. No such matching
dataset has been established here; that is a limit of the evidence currently
held, not proof that no suitable dataset exists.
"""
    OUT_MD.write_text(md, encoding="utf-8")


def main() -> int:
    r = build()
    OUT_JSON.write_text(json.dumps(r, indent=2) + "\n", encoding="utf-8")
    write_report(r)
    print(
        f"identifiability report: DOF={r['degrees_of_freedom']}, "
        f"non-identifiable={r['single_cell_sobol']['non_identifiable_count']}, "
        f"data-conditioned headlines={r['data_constrained_in_production']}"
    )
    print(f"wrote {OUT_MD.relative_to(REPO_ROOT)} + {OUT_JSON.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
