# CTRPv2 ferroptosis-inducer dose-response: the #330 calibration target

This is the first **real-data** calibration target for the ferroptosis-core kill
switch, which until now has been honest-but-uncalibrated scaffolding (see
`simulations/calibration/CALIBRATION_STATUS.md`). It is the **data leg** of #330;
the model fit + held-out validation is the follow-up leg (scope at the bottom).

Generated artifacts (regenerate with `python3 scripts/fetch_calibration_data.py`):
- `analysis/calibration/ctrpv2_ferroptosis_curves.csv`: per-cell-line 4-parameter
  logistic fits (committed; 3,021 curves).
- `analysis/calibration/ctrpv2_ferroptosis_summary.json`: per-compound summary plus
  provenance (committed).

## Source and provenance

- **Dataset**: CTRPv2 (Cancer Therapeutics Response Portal v2), the standard
  small-molecule sensitivity screen that includes ferroptosis inducers. Primary
  refs: Seashore-Ludlow 2015 Cancer Discov (PMID 26181016), Rees 2016 Nat Chem
  Biol (PMID 26656090).
- **Why not GDSC**: GDSC2 does not screen the canonical ferroptosis inducers
  (its panel has only Cisplatin / Elesclomol / Sorafenib among ferroptosis-
  relevant compounds), so it cannot calibrate a GPX4-inhibitor kill switch.
- **Access**: the original NCI CTD2 portal URL is dead (301 to studycatalog.cancer.gov).
  We pull the reprocessed CTRPv2 that DepMap redistributes, release
  **"Harmonized CTD^2 25Q2"**, file `CTRPResponseCurves.csv`
  (md5 `36115a26e2cabab33906793ba04c47f8`), via the login-free DepMap download API.
  The raw 45 MB file is NOT committed; only the filtered ferroptosis-inducer
  derivative is. CI never downloads (offline contract preserved); the committed
  CSV is the reproducible artifact downstream.

## The compounds and what they map to

| compound | role | model analog |
|---|---|---|
| ML162, ML210 | direct GPX4 inhibitors | the RSL3 / `rsl3_gpx4_inhib` kill switch |
| erastin | system-xc- (cystine import) inhibitor | the GSH-depletion / `gsh` arm |
| CIL55, CIL56 | additional ferroptosis probes | secondary checks |

## Empirical dose-response (the target numbers)

Each cell line has a fitted 4-parameter logistic. The fitted `Slope` is
**negative** for a cytotoxic curve, so viability is reconstructed as
`v(d) = lower + (upper - lower) / (1 + (d/EC50)^(-slope))`, which is `upper` (~1) at
low dose and `lower` (the residual viability) at saturating dose. `kill_ceiling =
1 - residual_viability`. `AUC` is the mean viability across the screened log-dose
range (lower = more potent).

| compound | n cell lines | median EC50 (µM) | EC50 IQR (µM) | median kill ceiling | median AUC |
|---|---:|---:|---:|---:|---:|
| ML210 | 762 | 0.53 | 0.11 to 1.58 | 0.89 | 0.69 |
| ML162 | 795 | 0.75 | 0.12 to 1.84 | 0.98 | 0.64 |
| erastin | 795 | 4.65 | 2.72 to 11.61 | 0.78 | 0.85 |
| CIL56 | 380 | 0.66 | 0.23 to 1.28 | 0.98 | 0.90 |
| CIL55 | 289 | 7.23 | 4.30 to 10.1 | 0.93 | 0.97 |

What the GPX4-inhibitor target (ML162/ML210) tells the model:
1. **Potency**: median EC50 ~0.5 to 0.75 µM for direct GPX4 inhibition.
2. **Kill ceiling**: GPX4 inhibition is near-complete at saturating dose for the
   median line (ML162 ~98%, ML210 ~89%), i.e. most cell lines are killable by
   GPX4i but a residual-viable fraction remains, and that fraction varies.
3. **Heterogeneity**: EC50 spans roughly an order of magnitude across cell lines
   (IQR ~0.1 to 1.8 µM, with a non-responder tail of very high fitted EC50). This is
   the real-world analog of the model's phenotype/parameter heterogeneity, and is
   the most useful signal: the model should reproduce a *distribution* of
   sensitivities, not a single number.

erastin (system-xc-) is less potent (median EC50 ~4.7 µM) and has a higher
residual viability (~0.22), consistent with its more indirect, GSH-depletion
mechanism, a useful contrast to direct GPX4 inhibition.

## Caveats the fit leg must handle (do NOT skip these)

This file is the target, not a calibration. Turning it into a calibrated kill
switch requires two mappings that are genuinely non-trivial and must be stated:

1. **Concentration units**: CTRPv2 dose is µM of a specific compound; the model's
   drug-intensity parameter is dimensionless. The fit must define an explicit
   µM to model-intensity mapping (e.g. fit a scale so the model EC50 matches the
   empirical median EC50), and report it as a fitted nuisance parameter, not a
   physical constant.
2. **Cell line to phenotype**: CTRPv2 cell lines are not the model's abstract
   phenotypes (Glycolytic / OXPHOS / Persister / NRF2-rescued). The defensible
   target is the *distribution* of dose-response across cell lines, which the
   model's phenotype/parameter spread should reproduce, rather than a per-line
   identity match.

Because of (1) and (2), the calibration is a **distributional** one: fit the kill
rates so the model's dose-response (potency, kill ceiling, and spread) matches the
GPX4-inhibitor target, with **held-out validation** (e.g. fit on ML162, validate
on ML210; or fit on a random half of cell lines, validate on the other half).
Until that lands, the kill switch stays labelled uncalibrated in
`CALIBRATION_STATUS.md`.
