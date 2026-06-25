# Practical-identifiability of the headline outputs (#503)

A consolidated, reproducible accounting of which headline simulation results are
point-estimable and which are directional-only, synthesizing the committed
sensitivity and uncertainty analyses (it does not re-run them; the Sobol/Morris/
ABC steps need the compiled extension). The structural facts (the 11
swept parameters and the named non-identifiable set) are cross-checked against
`analysis/prcc-results.json`.

## Headline accounting

The degrees of freedom: **11 free rate constants** are swept
(`fenton_rate, gsh_scav_efficiency, lp_rate, lp_propagation, gpx4_rate, fsp1_rate, nrf2_gsh_rate, gpx4_degradation_by_ros, death_threshold, sdt_ros, rsl3_gpx4_inhib`). Of these, **6
are practically non-identifiable from the kill rate** (Sobol total-effect ST < 0.05:
gsh_scav_efficiency, nrf2_gsh_rate, fsp1_rate, fenton_rate, death_threshold, gpx4_degradation_by_ros); three dominate
(lp_propagation ST=0.504, gpx4_rate ST=0.285, lp_rate ST=0.177).

**Data-constrained in the production regime: 0.**
The production simulation matrix uses fixed in-vivo defaults; the only data-conditioned fit is the in-vitro single-cell switch (#330), whose posterior is DISJOINT from the in-vivo/spatial regime that carries the headlines. So zero of the headline outputs are conditioned on data.

## Per-headline verdicts

### single-cell death rate (e.g. Persister x RSL3), the Figure 7 numbers

- **Drivers:** lp_propagation, gpx4_rate, lp_rate
- **Non-identifiable:** gsh_scav_efficiency, nrf2_gsh_rate, fsp1_rate, fenton_rate, death_threshold, gpx4_degradation_by_ros
- **Prior-predictive spread:** Persister x RSL3 point 42.5%, but 95% prior-predictive [1.6%, 99.7%] (width 98.1%); PersisterNrf2 x RSL3 point 0.0%, interval [0.0%, 37.8%]
- **Data-conditioned:** in-vitro only (ABC posterior); the in-vivo priors that produce the Figure 7 numbers are DISJOINT from the in-vitro data, so this headline cannot be conditioned on the data we hold
- **Verdict:** `directional_only` (the point estimate is essentially uninformative under the documented parameter uncertainty (the interval nearly spans [0,1]); the robust claim is that the differential between phenotypes exists, not its magnitude)
- **Source:** sobol-sensitivity-report.md (#331); uncertainty-intervals-report.md (#332); abc-posterior-report.md (#332)
### RSL3 + FSP1i Bliss synergy_score (the ~1.99x)

- **Drivers:** lp_propagation, gsh_scav_efficiency, gpx4_rate
- **Non-identifiable:** sdt_ros (structural zero, no SDT in the pair), rsl3_gpx4_inhib (structural zero, fixed DrugEffect)
- **Prior-predictive spread:** point ~1.99x, 95% prior-predictive ~[1.0x, 5.2x], median ~1.35x; strongly interaction-laden (Morris sigma > mu* for every active parameter)
- **Data-conditioned:** no (prior-predictive; the combo fit is not data-conditioned)
- **Verdict:** `direction_robust_magnitude_not` (the supra-additive DIRECTION holds at the lower bound (interval stays >= 1.0x), so dual-pathway depletion beating single-pathway is defensible; the 1.99x magnitude is not)
- **Source:** headline-sensitivity-report.md (#331); headline-uncertainty-report.md (#332)
### SDT-minus-RSL3 hypoxic-zone kill gap (the kill-collapse asymmetry)

- **Drivers:** sdt_ros, lp_propagation, lp_rate, gsh_scav_efficiency
- **Non-identifiable:** none flagged for this headline
- **Prior-predictive spread:** gap stays POSITIVE across its 95% interval ~[0.16, 1.00] (median 0.96, point 0.87) under the O2-independent assumption
- **Data-conditioned:** no (prior-predictive); additionally ASSUMPTION-bracketed: the magnitude collapses from ~87.8% to ~0% under full SDT O2-dependence (#336/#358), which the lead clinical agent SONALA-001 occupies
- **Verdict:** `direction_robust_magnitude_not` (the asymmetry SIGN is parameter-robust, but the magnitude is both wide (interval) and assumption-dependent (the contested SDT O2-dependence), so it is an assumption restated quantitatively, not a calibrated prediction)
- **Source:** headline-sensitivity-report.md (#331); headline-uncertainty-tme-report.md (#332); manuscript Section 7.1
### per-tissue vessel-wall RSL3 kill (the 40% -> 1.8% behind the BBB)

- **Drivers:** lp_propagation, lp_rate, gpx4_rate
- **Non-identifiable:** none flagged for this headline
- **Prior-predictive spread:** very wide per-tissue intervals (e.g. well-vascularized median 0.23 ~[0.00, 0.93]; CNS/BBB median 0.04 ~[0.00, 0.77]) because the bistable switch dominates; BUT the within-draw across-tissue ORDERING (well >= poorly >= CNS) held in 300/300 draws
- **Data-conditioned:** no (prior-predictive; transport params at fixed uncalibrated presets, their own ranges not swept)
- **Verdict:** `direction_robust_magnitude_not` (the penetration-gradient ordering is parameter-robust (per-draw, not inferred from the overlapping marginals); the absolute per-tissue magnitudes are not point-estimable)
- **Source:** headline-uncertainty-penetration-report.md (#332)
### SDT:RSL3 immune-kill ratio (the 104:1)

- **Drivers:** lp_propagation, lp_rate, sdt_ros
- **Non-identifiable:** none flagged for this headline
- **Prior-predictive spread:** the 104:1 (2D, near DAMP saturation) falls to ~4:1 in 3D geometry; SDT de-confounded immune rate ~[0.009, 0.171], robustly low-but-positive
- **Data-conditioned:** no; geometry-dependent (the 2D-vs-3D shrink is a structural, not parametric, effect)
- **Verdict:** `directional_only` (the ratio is presented as a directional ceiling, not a number: it changes ~25x with geometry alone, which a parametric analysis cannot capture, so only the direction (SDT >> RSL3 immune priming) is claimed)
- **Source:** headline-sensitivity-report.md (#331); CALIBRATION_STATUS.md immune row


## Overall

No headline output is fully point-estimable. The single-cell kill rate and the immune ratio are directional-only; the Bliss synergy, the hypoxia asymmetry, and the penetration gap are direction-robust but magnitude-uncalibrated. With 11 free rate constants, 6 non-identifiable from the kill rate, and 0 of the headlines data-conditioned in the production regime, the honest reading of every reported magnitude is order-of-magnitude / directional, exactly as the manuscript labels them.

## What would make a headline point-estimable

A headline becomes point-estimable when (1) its driving parameters are identified
(narrowed) by data in the regime that produces it, and (2) the prior-predictive
interval collapses to a usable width. Concretely: the multi-inducer joint fit
(#500) plus System Xc- in the core (#502) would condition the LP-cascade and
defense constants in a calibrated regime; until then, the manuscript's
order-of-magnitude / directional labeling is the correct one, and this report is
the standing evidence for it.
