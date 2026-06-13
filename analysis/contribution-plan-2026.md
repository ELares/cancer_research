# Contribution Plan (2026): preprint, preregistration, and a collaborator brief

A concrete plan to turn this repository from "interesting open work" into a
**citable, falsifiable contribution**. It has three parts, matching issue #353:

1. a **venue split** (the corpus synthesis and the simulation methods are
   different contributions and should be published separately);
2. a **preregistration draft** of the model's specific, falsifiable predictions,
   to be registered BEFORE any calibration so that calibrated-versus-predicted
   stays honest;
3. a **collaborator-facing experiment brief**: the cheapest wet-lab experiments
   that would validate or falsify the headline predictions.

This is a **plan, not a result**. Venue names below are candidate TARGETS with
rationale, not acceptances. The predictions are **directional / order-of-
magnitude** (see `MODEL_CARD.md` and `article/drafts/v1.md` Section 8.4); the
point of preregistering them is to be honest about which survive contact with
data.

---

## Part 1: Venue split

The repository contains two genuinely different contributions that should not be
forced into one paper:

### (a) Cross-literature corpus synthesis

**What it is:** a reproducible map of the cancer-therapy literature (4,830
full-text articles across 19 mechanisms and 22 cancer types), with honest
coverage caveats (96% tagger precision, 55% recall, taxonomy-dependent gap
counts, open-access skew), and the methodological point that several apparent
"gaps" are artifacts of search design rather than biology.

- **Preprint:** post to a general life-sciences / biomedical preprint server so
  the synthesis is citable immediately and open to comment.
- **Journal target (candidates, by fit):** an evidence-synthesis /
  research-on-research / scientometrics venue, or a cancer-informatics venue.
  The contribution is a *methods + landscape* paper (reproducible tagging
  pipeline + honest coverage analysis), NOT a clinical-claims paper, so target
  venues that publish reproducible evidence-mapping methodology.
- **Framing discipline:** lead with the reproducible method and the coverage
  caveats; the mechanism rankings are descriptive, and absence claims are
  reported as "not detected in the local keyword-derived analysis."

### (b) Simulation methods + honest validation status

**What it is:** an embeddable, open ferroptosis biochemistry engine
(`ferroptosis-core`) plus 2D/3D tumor-microenvironment layers, with a *model
card* and a *per-layer calibration ledger* that states plainly the suite is
broad but mostly uncalibrated.

- **Preprint:** post to a quantitative-biology / computational-biology preprint
  category.
- **Journal target (candidates, by fit):** a computational / systems-biology
  *methods* venue. The honest framing (a reusable engine + an explicit
  uncalibrated-status ledger + a falsification roadmap) is itself the
  contribution; a methods venue is the right home, NOT a venue expecting
  validated quantitative predictions.
- **Framing discipline:** the headline is the *tool and the honesty
  infrastructure* (model card, byte-identity discipline, calibration tiers),
  not the specific numbers. Every quantitative figure is labeled directional /
  order-of-magnitude.

### Why split

The two have different audiences, different review criteria, and different
failure modes. Bundling them would let a weakness in one (e.g. the uncalibrated
sim) undermine the other (the corpus method), and would invite reviewers to
judge the sim as if it claimed calibrated predictions, which it does not.

---

## Part 2: Preregistration draft (register BEFORE calibration)

Register these on a preregistration platform (e.g. OSF) as the model's standing,
falsifiable predictions, with their pre-stated falsification criteria, BEFORE
fitting any layer to external data. Each is a *directional* prediction; the
preregistration locks in the direction and the rough magnitude so that a later
calibration cannot quietly move the goalposts. Sources are the manuscript's
prediction/falsification pairs (the boxes live in Chapters 6 and 7; Chapter 9
organizes them into the falsification roadmap).

| # | Prediction (directional) | Pre-stated falsification |
|---|---|---|
| P1 | **GPX4 + FSP1 dual inhibition is synergistic** (Bliss > 1, model gives ~1.99x) in persister-enriched cells. | Combination index > 0.8 (additive/antagonistic), or CI > 1.0. |
| P2 | **Physical-ROS modalities (PDT/SDT) are less penetration-limited with depth than systemic RSL3** in >=500 um spheroids. | All three modalities show similar depth-dependent kill profiles. |
| P3 | **The post-withdrawal ferroptosis-vulnerability window closes on a days timescale, with defenses recovering sequentially** (FSP1/GSH first). | Defenses recover simultaneously, or recovery is hours not days (window too narrow to exploit). |
| P4 | **SDT retains more efficacy than RSL3 under hypoxia** (direction only; magnitude explicitly contested). | RSL3 and SDT show similar hypoxic efficacy loss (kills the ROS-source-asymmetry assumption). NOTE this is the model's least-certain leg; the off-by-default O2-dependent SDT mode and dynamic-iron hypoxia coupling exist precisely so the reverse can be modeled. |
| P5 | **Dense ferroptotic kill produces more immunogenic-cell-death signal per cell than sparse kill** (direction; the 104:1 ratio is a 2D ceiling that falls to ~4:1 in 3D and can flip sign under the immunosuppressive-ferroptosis arm). | DAMP release is proportional to dead-cell count regardless of geometry/quality. |
| P6 | **CAF coculture protects RSL3 more than SDT** (stromal GSH/MUFA supply blunts pharmacologic but not exogenous ROS). | CAFs protect both equally, or neither. |
| P7 | **RSL3 efficacy drops at acidic pH (6.5 vs 7.4) via ion trapping** (most-uncertain parameter: RSL3 pKa). | RSL3 efficacy is pH-independent (invalidates the pH-resistance leg). |
| P8 | **A persister-targeting ferroptosis inducer (RSL3) has the OPPOSITE size-dependence to generic cytotoxics**: small, all-proliferating spheroids resist it (no persister target) and vulnerability emerges as the persister core appears (~280 um radius), whereas generic cytotoxics kill smaller spheroids better. Two competing size-effects, the supply gradient (bigger resists more) vs the persister-targeting (bigger has more target), so the net direction depends on the inducer's phenotype specificity. (#333; `analysis/calibration/spheroid-kill-vs-size.md`.) | A ferroptosis inducer shows the SAME monotone bigger-resists-more size-dependence as generic cytotoxics (no persister-targeting inversion at small size), or its kill is size-independent. |

**Honesty clause to include in the registration:** P4 and P7 are flagged as the
least certain (contested SDT oxygen-dependence; poorly-characterized RSL3 pKa);
the registration commits to reporting failures of these as prominently as
successes.

---

## Part 3: Collaborator-facing experiment brief

The cheapest wet-lab experiments that would validate or falsify the headline
predictions. Ordered by cost/accessibility. Each lists the model's prediction,
the measurement, and the falsifying outcome, so a collaborator can run one in
isolation. (These map to the manuscript's "How to test this prediction" boxes in
Chapters 6 and 7.) The brief covers P1 and P3 through P7; P2's depth-penetration
leg is folded into E1's confocal depth-sectioning (and is the standalone spheroid
proposal in manuscript Section 6.1), so it has no separate experiment here.

### E1. Spheroid RSL3 vs SDT kill at measured hypoxia (tests P4, the keystone)

- **Setup:** persister-enriched cells in >=500 um spheroids with verified hypoxic
  cores (pimonidazole), or a hypoxia chamber (21% vs 1% O2). Apply RSL3 and SDT.
- **Measure:** viability (and depth-resolved viability by confocal sectioning);
  pO2 to anchor the hypoxia axis.
- **Model says:** RSL3 kill collapses under hypoxia far more than SDT.
- **Falsifies if:** both collapse similarly (the model's central ROS-source
  asymmetry is wrong) OR SDT collapses more (it is O2-dependent, as the lead
  clinical agent suggests, and the optimistic 2D upper bound is unjustified).
- **Cost:** low-moderate (standard 3D culture + hypoxia readout).

### E2. CAF-coculture IC50 shift (tests P6)

- **Setup:** tumor cells +/- patient-derived or established CAFs, 4 arms (alone vs
  CAF) x (RSL3 vs SDT), dose-response.
- **Measure:** RSL3 and SDT IC50 shift with CAFs; C11-BODIPY lipid peroxidation;
  GSH.
- **Model says:** CAFs raise the RSL3 IC50 (stromal GSH/MUFA shielding) more than
  the SDT IC50.
- **Falsifies if:** CAFs shield both equally (CAF antioxidants neutralize even the
  exogenous ROS burst, or the basal-vs-exogenous ROS-dose asymmetry is
  overstated) or neither (the metabolic-cooperation model is weak).
- **Cost:** low-moderate (coculture + viability/IC50).

### E3. Spheroid-supernatant DAMP + DC-maturation assay (tests P5)

- **Setup:** SDT vs RSL3 on 2D monolayer (uniform kill) vs 3D spheroid (spatially
  structured kill).
- **Measure:** calreticulin surface exposure (flow), HMGB1 (ELISA), ATP
  (luminescence); then load supernatant onto dendritic cells and measure DC
  maturation / cross-presentation.
- **Model says:** dense ferroptotic kill yields more ICD signal per cell and more
  DC maturation; but the immunosuppressive arm (extracellular GPX4 / oxidized
  lipids) can flip the net at high death density.
- **Falsifies if:** DAMP/DC-maturation tracks total dead-cell count regardless of
  modality, geometry, or kill quality (the LP-overshoot DAMP-quality differential
  is wrong).
- **Cost:** moderate (adds the DC-maturation readout).

### E4. GPX4 + FSP1 dose-matrix synergy (tests P1)

- **Setup:** RSL3 (or ML162) x iFSP1 (or brequinar) dose matrix in
  persister-enriched cells.
- **Measure:** Chou-Talalay combination index or Bliss; C11-BODIPY to confirm the
  lipid-peroxidation pathway.
- **Model says:** synergy (CI < 1, ~1.99x Bliss).
- **Falsifies if:** CI > 0.8 (additive/antagonistic), revising the
  parallel-independent-repair assumption.
- **Cost:** low (standard combination assay).

### E5. Sequential defense recovery after drug withdrawal (tests P3)

- **Setup:** withdraw a persister-inducing therapy; serial timepoints (0h, 6h,
  1d, 3d, 1wk, 2wk, 4wk).
- **Measure:** GPX4, FSP1, NRF2, GSH; matched RSL3-sensitivity at each timepoint.
- **Model says:** a transient vulnerability window with sequential recovery.
- **Falsifies if:** simultaneous recovery, or recovery in hours not days.
- **Cost:** moderate (time-course).

### E6. RSL3 efficacy and intracellular concentration vs pH (tests P7)

- **Setup:** RSL3 at pH 7.4 vs 6.5.
- **Measure:** efficacy; intracellular RSL3 (HPLC or fluorescent analog) to test
  ion trapping directly.
- **Model says:** lower efficacy at acidic pH via ion trapping.
- **Falsifies if:** pH-independent efficacy (RSL3 chloroacetamide chemistry does
  not behave as a classic weak base; the pH-resistance leg is invalidated).
- **Cost:** low.

---

## Sequencing

1. **Now:** post both preprints (Part 1) and register P1-P7 (Part 2) BEFORE any
   calibration work begins (calibration is tracked in the #330-#335 issues).
2. **Recruit:** circulate the Part 3 brief to recruit a wet-lab collaborator for
   E1 (the keystone hypoxia leg) and E4 (the cheapest, the synergy claim).
3. **Calibrate honestly:** as calibration data arrives (#330-#335), report
   calibrated-versus-preregistered for each prediction, failures included.

## Cross-references

- Falsification chapter: `article/drafts/v1.md` Chapter 9.
- Honest scope and validation status: `MODEL_CARD.md`,
  `simulations/calibration/CALIBRATION_STATUS.md`.
- Gap analysis and backlog: `analysis/research-roadmap-2026.md`.
- Scientific cross-check: `analysis/manuscript-scientific-review.md`.
