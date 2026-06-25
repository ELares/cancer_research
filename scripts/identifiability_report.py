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
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
PRCC = REPO_ROOT / "analysis" / "prcc-results.json"
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
        "data_conditioned": "in-vitro only (ABC posterior); the in-vivo priors that produce the Figure 7 numbers are DISJOINT from the in-vitro data, so this headline cannot be conditioned on the data we hold",
        "verdict": "directional_only",
        "rationale": "the point estimate is essentially uninformative under the documented parameter uncertainty (the interval nearly spans [0,1]); the robust claim is that the differential between phenotypes exists, not its magnitude",
        "source": "sobol-sensitivity-report.md (#331); uncertainty-intervals-report.md (#332); abc-posterior-report.md (#332)",
    },
    {
        "key": "bliss_synergy",
        "observable": "RSL3 + FSP1i Bliss synergy_score (the ~1.99x)",
        "drivers": ["lp_propagation", "gsh_scav_efficiency", "gpx4_rate"],
        "non_identifiable_params": ["sdt_ros (structural zero, no SDT in the pair)", "rsl3_gpx4_inhib (structural zero, fixed DrugEffect)"],
        "prior_predictive": "point ~1.99x, 95% prior-predictive ~[1.0x, 5.2x], median ~1.35x; strongly interaction-laden (Morris sigma > mu* for every active parameter)",
        "data_conditioned": "no (prior-predictive; the combo fit is not data-conditioned)",
        "verdict": "direction_robust_magnitude_not",
        "rationale": "the supra-additive DIRECTION holds at the lower bound (interval stays >= 1.0x), so dual-pathway depletion beating single-pathway is defensible; the 1.99x magnitude is not",
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
        "prior_predictive": "the 104:1 (2D, near DAMP saturation) falls to ~4:1 in 3D geometry; SDT de-confounded immune rate ~[0.009, 0.171], robustly low-but-positive",
        "data_conditioned": "no; geometry-dependent (the 2D-vs-3D shrink is a structural, not parametric, effect)",
        "verdict": "directional_only",
        "rationale": "the ratio is presented as a directional ceiling, not a number: it changes ~25x with geometry alone, which a parametric analysis cannot capture, so only the direction (SDT >> RSL3 immune priming) is claimed",
        "source": "headline-sensitivity-report.md (#331); CALIBRATION_STATUS.md immune row",
    },
]


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
            "The production simulation matrix uses fixed in-vivo defaults; the only "
            "data-conditioned fit is the in-vitro single-cell switch (#330), whose "
            "posterior is DISJOINT from the in-vivo/spatial regime that carries the "
            "headlines. So zero of the headline outputs are conditioned on data."
        ),
        "headlines": HEADLINES,
        "overall": (
            "No headline output is fully point-estimable. The single-cell kill rate "
            "and the immune ratio are directional-only; the Bliss synergy, the "
            "hypoxia asymmetry, and the penetration gap are direction-robust but "
            "magnitude-uncalibrated. With 11 free rate constants, 6 non-identifiable "
            "from the kill rate, and 0 of the headlines data-conditioned in the "
            "production regime, the honest reading of every reported magnitude is "
            "order-of-magnitude / directional, exactly as the manuscript labels them."
        ),
    }


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
interval collapses to a usable width. Concretely: the multi-inducer joint fit
(#500) plus System Xc- in the core (#502) would condition the LP-cascade and
defense constants in a calibrated regime; until then, the manuscript's
order-of-magnitude / directional labeling is the correct one, and this report is
the standing evidence for it.
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
