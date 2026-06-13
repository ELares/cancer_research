# Preregistration of Falsifiable Predictions

This document registers the model's specific, directional, falsifiable predictions
and their pre-stated falsification criteria **before** the calibration and wet-lab
work that would confirm or refute them. The point is to lock in what would count
as success and what would count as failure ahead of time, so that no prediction
can be quietly re-fit after the data arrives. This operationalizes the first
guiding principle of the project (let the evidence lead): the repository is here
to help people, not to defend one hypothesis, and a preregistration is the
clearest way to commit to that in public.

The predictions and experiment briefs below are reproduced from
`analysis/contribution-plan-2026.md` (Parts 2 and 3), which remains the working
source. This file is the registrable, time-stampable version.

## Registration status

> **OSF registration is a human step and is not yet complete.** Registering on the
> Open Science Framework (osf.io) requires a free account and a few clicks, and it
> mints an immutable, time-stamped DOI. Once that is done, paste the DOI and the
> public registration URL here and in `analysis/contribution-plan-2026.md`.
>
> - **OSF registration URL:** _pending (TODO: paste after registering)_
> - **DOI:** _pending (TODO: paste after registering)_
> - **Registered on:** _pending_
>
> Until the OSF DOI exists, the git commit history of this file is the time stamp:
> the predictions below are fixed as of the commit that added this document, and
> any later change is visible in the diff.

## Part 1: Falsifiable predictions (P1 to P8)

Each prediction is **directional**: it states the sign of an effect, not a
calibrated magnitude, because most of the simulation layers are uncalibrated
mechanistic scaffolding (see `simulations/calibration/CALIBRATION_STATUS.md` and
`MODEL_CARD.md`). The falsification column states the outcome that would kill the
prediction.

| # | Prediction (directional) | Pre-stated falsification |
|---|--------------------------|--------------------------|
| P1 | GPX4 plus FSP1 dual inhibition is synergistic (Bliss greater than 1, the model gives about 1.99x) in persister-enriched cells. | Combination index greater than 0.8 (additive or antagonistic), or CI greater than 1.0. |
| P2 | Physical-ROS modalities (PDT and SDT) are less penetration-limited with depth than systemic RSL3 in spheroids of at least 500 um radius. | All three modalities show similar depth-dependent kill profiles. |
| P3 | The post-withdrawal ferroptosis-vulnerability window closes on a days timescale, with defenses recovering sequentially (FSP1 and GSH first). | Defenses recover simultaneously, or recovery is hours not days (the window is too narrow to exploit). |
| P4 | SDT retains more efficacy than RSL3 under hypoxia (direction only; magnitude explicitly contested). | RSL3 and SDT show similar hypoxic efficacy loss (this kills the ROS-source-asymmetry assumption). This is the model's least-certain leg; the off-by-default O2-dependent SDT mode and dynamic-iron hypoxia coupling exist precisely so the reverse can be modeled. |
| P5 | Dense ferroptotic kill produces more immunogenic-cell-death signal per cell than sparse kill (direction; the 104:1 ratio is a 2D ceiling that falls to about 4:1 in 3D and can flip sign under the immunosuppressive-ferroptosis arm). | DAMP release is proportional to dead-cell count regardless of geometry or quality. |
| P6 | CAF coculture protects RSL3 more than SDT (stromal GSH and MUFA supply blunts pharmacologic but not exogenous ROS). | CAFs protect both equally, or neither. |
| P7 | RSL3 efficacy drops at acidic pH (6.5 versus 7.4) via ion trapping (most-uncertain parameter: RSL3 pKa). | RSL3 efficacy is pH-independent (this invalidates the pH-resistance leg). |
| P8 | A persister-targeting ferroptosis inducer (RSL3) has the OPPOSITE size-dependence to generic cytotoxics: small, all-proliferating spheroids resist it (no persister target) and vulnerability emerges as the persister core appears (around 280 um radius), whereas generic cytotoxics kill smaller spheroids better. There are two competing size-effects, the supply gradient (bigger resists more) versus the persister-targeting (bigger has more target), so the net direction depends on the inducer's phenotype specificity. (#333; `analysis/calibration/spheroid-kill-vs-size.md`.) | A ferroptosis inducer shows the SAME monotone bigger-resists-more size-dependence as generic cytotoxics (no persister-targeting inversion at small size), or its kill is size-independent. |

### Honesty clause

P4 (SDT oxygen-dependence) and P7 (RSL3 pKa) are flagged as the **least certain**
predictions: the SDT oxygen-dependence is genuinely contested in the field (the
lead clinical sonosensitizer is itself oxygen-dependent), and the RSL3 pKa that
the pH-resistance leg rests on is poorly characterized. This registration commits
to reporting failures of these two predictions as prominently as any success. If
the data refutes P4 or P7, that refutation will be stated plainly in the
manuscript and in the calibration status, not buried.

## Part 2: Collaborator-facing experiment briefs (E1 to E6)

These are the cheapest wet-lab experiments that would validate or falsify the
headline predictions, ordered by cost and accessibility. Each lists the model's
prediction, the measurement, and the falsifying outcome, so a collaborator can run
one in isolation. They map to the manuscript's "How to test this prediction" boxes
in Chapters 6 and 7. The brief covers P1 and P3 through P7; P2's depth-penetration
leg is folded into E1's confocal depth-sectioning, so it has no separate
experiment.

### E1. Spheroid RSL3 versus SDT kill at measured hypoxia (tests P4, the keystone)

- **Setup:** persister-enriched cells in spheroids of at least 500 um radius with
  verified hypoxic cores (pimonidazole), or a hypoxia chamber (21% versus 1% O2).
  Apply RSL3 and SDT.
- **Measure:** viability (and depth-resolved viability by confocal sectioning);
  pO2 to anchor the hypoxia axis.
- **Model says:** RSL3 kill collapses under hypoxia far more than SDT.
- **Falsifies if:** both collapse similarly (the model's central ROS-source
  asymmetry is wrong) OR SDT collapses more (it is O2-dependent, as the lead
  clinical agent suggests, and the optimistic 2D upper bound is unjustified).
- **Cost:** low to moderate (standard 3D culture plus hypoxia readout).

### E2. CAF-coculture IC50 shift (tests P6)

- **Setup:** tumor cells with or without patient-derived or established CAFs, 4
  arms (alone versus CAF) by (RSL3 versus SDT), dose-response.
- **Measure:** RSL3 and SDT IC50 shift with CAFs; C11-BODIPY lipid peroxidation; GSH.
- **Model says:** CAFs raise the RSL3 IC50 (stromal GSH and MUFA shielding) more
  than the SDT IC50.
- **Falsifies if:** CAFs shield both equally (CAF antioxidants neutralize even the
  exogenous ROS burst, or the basal-versus-exogenous ROS-dose asymmetry is
  overstated) or neither (the metabolic-cooperation model is weak).
- **Cost:** low to moderate (coculture plus viability/IC50).

### E3. Spheroid-supernatant DAMP plus DC-maturation assay (tests P5)

- **Setup:** SDT versus RSL3 on 2D monolayer (uniform kill) versus 3D spheroid
  (spatially structured kill).
- **Measure:** calreticulin surface exposure (flow), HMGB1 (ELISA), ATP
  (luminescence); then load supernatant onto dendritic cells and measure DC
  maturation and cross-presentation.
- **Model says:** dense ferroptotic kill yields more ICD signal per cell and more
  DC maturation; but the immunosuppressive arm (extracellular GPX4 and oxidized
  lipids) can flip the net at high death density.
- **Falsifies if:** DAMP and DC-maturation track total dead-cell count regardless
  of modality, geometry, or kill quality (the LP-overshoot DAMP-quality
  differential is wrong).
- **Cost:** moderate (adds the DC-maturation readout).

### E4. GPX4 plus FSP1 dose-matrix synergy (tests P1)

- **Setup:** RSL3 (or ML162) by iFSP1 (or brequinar) dose matrix in
  persister-enriched cells.
- **Measure:** Chou-Talalay combination index or Bliss; C11-BODIPY to confirm the
  lipid-peroxidation pathway.
- **Model says:** synergy (CI less than 1, about 1.99x Bliss).
- **Falsifies if:** CI greater than 0.8 (additive or antagonistic), revising the
  parallel-independent-repair assumption.
- **Cost:** low (standard combination assay).

### E5. Sequential defense recovery after drug withdrawal (tests P3)

- **Setup:** withdraw a persister-inducing therapy; serial timepoints (0h, 6h, 1d,
  3d, 1wk, 2wk, 4wk).
- **Measure:** GPX4, FSP1, NRF2, GSH; matched RSL3-sensitivity at each timepoint.
- **Model says:** a transient vulnerability window with sequential recovery.
- **Falsifies if:** simultaneous recovery, or recovery in hours not days.
- **Cost:** moderate (time-course).

### E6. RSL3 efficacy and intracellular concentration versus pH (tests P7)

- **Setup:** RSL3 at pH 7.4 versus 6.5.
- **Measure:** efficacy; intracellular RSL3 (HPLC or fluorescent analog) to test
  ion trapping directly.
- **Model says:** lower efficacy at acidic pH via ion trapping.
- **Falsifies if:** pH-independent efficacy (RSL3 chloroacetamide chemistry does
  not behave as a classic weak base; the pH-resistance leg is invalidated).
- **Cost:** low.

## Part 3: Calibration status at registration

Registering before calibration is the point. For transparency, here is what was
and was not anchored to independent data as of registration. The full per-layer
ledger is `simulations/calibration/CALIBRATION_STATUS.md`.

- **Calibrated (in-vitro, held-out):** the single-cell RSL3 kill switch, fit to
  CTRPv2 GPX4-inhibitor dose-response on ML162 and validated on ML210 (#330).
- **Partially anchored to published data:** spheroid radial zone geometry versus
  Browning 2021 (#333, geometry only); tumor-PK disposition versus imidazole
  ketone erastin plus a sorafenib forward check (#334); Krogh penetration form and
  reference length versus Primeau/Tannock (#335).
- **Prior-predictive only (parameter, not data-conditioned):** the spatial and
  combination headlines that P1, P2, P4, P5, P6, and P7 rest on. An ABC analysis
  (#332) shows the in-vivo priors and the in-vitro posterior are disjoint, so the
  in-vivo and spatial headlines cannot be conditioned on the available in-vitro
  data and stay prior-predictive until an in-vivo ferroptosis dataset exists.

As calibration data arrives, the plan is to report calibrated-versus-preregistered
for each prediction, failures included.

## Sequencing

1. **Now:** register P1 through P8 above on OSF and mint the DOI, then paste it
   into the registration-status block at the top of this file.
2. **Recruit:** circulate the Part 2 briefs to find a wet-lab collaborator for E1
   (the keystone hypoxia leg) and E4 (the cheapest, the synergy claim). These two
   issues are tracked as `help wanted` (#442 spheroid kill, #448 in-vivo PK).
3. **Calibrate honestly:** as data arrives, report calibrated-versus-preregistered
   for each prediction, failures stated as plainly as confirmations.

## Cross-references

- `analysis/contribution-plan-2026.md` (the working source for Parts 1 and 2)
- `simulations/calibration/CALIBRATION_STATUS.md` (per-layer calibration ledger)
- `MODEL_CARD.md` (intended use, out-of-scope, per-layer status)
- `analysis/research-roadmap-2026.md` (the gap analysis and issue backlog)
