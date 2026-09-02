# Preregistration of Falsifiable Predictions

This document registers the model's specific, directional, falsifiable predictions
and their pre-stated falsification criteria **before any wet-lab data exists**, and
**independently of** — not prior to — the model-side calibration legs.

That distinction is deliberate and it corrects an earlier version of this sentence,
which claimed registration came before the calibration work. It did not. Four of the
calibration legs were committed before this file was: the CTRPv2 kill-switch fit and
the spheroid zone geometry on 2026-06-06, the tumor-PK and Krogh penetration anchors
on 2026-06-07, against this document's 2026-06-13. The git history is public, so the
precedence claim was checkable and wrong.

What is true, and is the claim that actually carries the weight, is DISJOINTNESS:
none of those datasets is used to set any prediction below. The predictions come
from the model's own outputs; the calibration legs anchor different quantities
against independent data. A reader wanting the per-leg accounting will find it in
Part 3, which has always been honest about what was already anchored at
registration. The point is to lock in what would count
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

## Part 1: Falsifiable predictions (P1 to P13)

Each prediction is **directional**: it states the sign of an effect, not a
calibrated magnitude, because most of the simulation layers are uncalibrated
mechanistic scaffolding (see `simulations/calibration/CALIBRATION_STATUS.md` and
`MODEL_CARD.md`). Each prediction below states (a) a quantitative model output and (b) a
pre-stated numeric falsification threshold, the outcome that would kill it. The
claims stay directional (the sign of an effect), but the model output and the
threshold are now numeric so a result can be scored against them without
re-interpretation.

**P1. GPX4 plus FSP1 dual inhibition is synergistic in FSP1-low persister-enriched cells.**
- *Quantitative model output:* Bliss excess about 1.99x (95% prior-predictive interval about 1.0x to 5.2x; the supra-additive direction is robust at the lower bound, the magnitude is not).
- *Falsification threshold:* Chou-Talalay combination index greater than 0.8 at the matched-effect dose ratio (additive or antagonistic), or measured combined kill at or below the Bliss-independence prediction within assay error.

**P2. Physical-ROS modalities (PDT and SDT) are less depth-limited than systemic RSL3 in spheroids of at least 500 um radius.**
- *Quantitative model output:* at the spheroid core RSL3 kill falls toward zero while SDT retains most of its rim kill, so the RSL3 core-to-rim kill ratio is far below the SDT core-to-rim ratio (the same penetration asymmetry behind the 40 percent to 1.8 percent RSL3-like penetration headline).
- *Falsification threshold:* the core-to-rim kill ratio for RSL3 is within 1.5x of the SDT core-to-rim ratio (no differential depth penetration), or all three modalities' depth-kill half-distances agree within 25 percent.

**P3. The post-withdrawal ferroptosis-vulnerability window closes on a days timescale, with defenses recovering sequentially (FSP1 and GSH first, GPX4 and NRF2 later).**
- *Quantitative model output:* matched RSL3-sensitivity returns toward baseline over roughly 3 to 7 days, with FSP1 and GSH recovering before GPX4 and NRF2.
- *Falsification threshold:* matched RSL3-sensitivity returns to baseline within 24 hours (the window is too narrow to exploit), or all four defenses recover within the same timepoint (no sequential order).

**P4. SDT retains more efficacy than RSL3 under hypoxia (direction only; magnitude explicitly contested).**
- *Quantitative model output:* in the hypoxic zone RSL3 kill collapses to near zero while SDT kill is 86.6 percent under the O2-independent upper bound and falls to near zero only under the O2-dependent lower bound, so the modeled SDT-to-RSL3 hypoxic-kill ratio is at least 1 and bracketed wide (Section 7.1).
- *Falsification threshold:* the SDT hypoxic-kill loss (normoxic minus hypoxic) is greater than or equal to the RSL3 hypoxic-kill loss, i.e. the SDT-to-RSL3 hypoxic-kill ratio is at or below 1.0 (SDT collapses as much as or more than RSL3, killing the ROS-source-asymmetry assumption). This is the model's least-certain leg; the off-by-default O2-dependent SDT mode and dynamic-iron hypoxia coupling exist precisely so the reverse can be modeled.

**P5. Dense ferroptotic kill produces more immunogenic-cell-death signal per dead cell than sparse kill.**
- *Quantitative model output:* the dense-to-sparse DAMP-per-cell (or DC-maturation-per-cell) ratio is about 4:1 in 3D (down from the 2D 104:1 ceiling; it can fall below 1 under the immunosuppressive-ferroptosis arm at high death density).
- *Falsification threshold:* DAMP release and DC maturation per dead cell agree within 1.5x between dense (SDT or 3D) and sparse (RSL3) kill (signal proportional to dead-cell count regardless of geometry or kill quality).

**P6. CAF coculture protects RSL3 more than SDT.**
- *Quantitative model output:* CAFs raise the RSL3 IC50 fold-shift above the SDT IC50 fold-shift (boundary-cell kill: RSL3 halves, 3.0 percent to 1.5 percent, while SDT barely moves, 96.1 percent to 91.2 percent).
- *Falsification threshold:* the RSL3 and SDT IC50 fold-shifts with CAFs agree within 1.5x (equal shielding), or neither shifts more than 1.2x (no metabolic-cooperation effect).

**P7. RSL3 efficacy drops at acidic pH (6.5 versus 7.4) via ion trapping.**
- *Quantitative model output:* ferroptosis kills fall about 53 percent at pH 6.5 versus 7.4 (163 to 77 in the immune-free counter) while SDT is unaffected, i.e. the RSL3 IC50 rises with acidity.
- *Falsification threshold:* RSL3 efficacy (or IC50) agrees within 1.2x between pH 6.5 and 7.4 (pH-independent; RSL3 chloroacetamide chemistry does not behave as a classic weak base, invalidating the pH-resistance leg).

**P8. A persister-targeting ferroptosis inducer (RSL3) has the OPPOSITE size-dependence to generic cytotoxics.**
- *Quantitative model output:* with size-aware zone thresholds, RSL3 kill is near zero below about 280 um radius and rises as the persister core appears (a non-monotone profile), whereas generic cytotoxics fall monotonically with size (the model's fixed-threshold supply-gradient leg gives RSL3 kill 3.6 percent to 0.19 percent over 144 to 540 um). The net direction depends on the inducer's phenotype specificity (#333; `analysis/calibration/spheroid-kill-vs-size.md`).
- *Falsification threshold:* a ferroptosis inducer shows the same monotone bigger-resists-more profile as generic cytotoxics (no persister-targeting inversion below about 280 um radius), or its kill is size-independent.

### Predictions for the modality arms (P9 to P13)

P1 to P8 all concern ferroptosis or the physical-ROS modalities. Before the
five below were registered `analysis/scope-audit.md` reported that as 8 of 8;
it now reports 8 of 13. That is the sharpest measure of
this project's narrowness -- sharper than code volume or figure count, because
a falsifiable prediction is the currency this repository treats as real. An
engine that can EXPRESS nine modalities while having committed to being wrong
about only one is still a ferroptosis project by its own chosen measure. The
five below are registered against the arms added in the modality campaign.
Model outputs are derived in `analysis/modality-predictions.md`; the
thresholds are stated here.

**P9. PARP inhibition enhances radiation MORE at low dose per fraction than at high.**
- *Quantitative model output:* an alpha-only boost raises the single-hit term while beta is untouched, so the sensitizer enhancement ratio decays with dose. One boost must hold the published 1.2 to 1.7 band at surviving fractions [0.5, 0.1, 0.01], i.e. 2.24 to 9.38 Gy; that admits boosts 0.544 to 0.948, and across the whole of that window SER(2 Gy)/SER(6 Gy) is 1.113 to 1.168. The decay is analytic in dose, not an artefact of the sampled points.
- *Falsification threshold:* EITHER measured clonogenic SER at 6 Gy at or above SER at 2 Gy in the same line (ratio at or below 1.0, against the 1.113 the model requires at its most conservative admissible boost), OR measured SER outside 1.2 to 1.7 anywhere in 2.24 to 9.38 Gy -- the second clause matters because the ratio alone can be satisfied while the band the window was derived from is violated at every dose.
- *Stated scope:* the window is a claim about 2.24 to 9.38 Gy and nothing wider. Requiring the same band at 1.8 Gy and 20 Gy as well empties it, which is a limit of the alpha-only form and is registered as such rather than discovered later.

**P10. The ADC bystander effect is STARVED by the antigen escape it is supposed to answer.**
- *Quantitative model output:* bystander payload comes from cells that took up the ADC, so it scales with the dying antigen-positive population, and it is apportioned across the surviving pool rather than aimed at antigen-negative cells. The share of the antigen-NEGATIVE pool it reaches falls from 77.1% at 0.9 antigen-positive to 2.6% at 0.1, while the relative advantage over a non-cleavable linker stays exactly flat at 1.3x.
- *Falsification threshold:* in a graded-antigen co-culture the cleavable arm's kill of antigen-NEGATIVE cells, as a share of the antigen-negative population, is flat or rising as the antigen-positive fraction falls; or the cleavable-to-non-cleavable total-kill ratio rises by more than 10 percent across the same range.
- *Note:* this has the OPPOSITE sign to the intuition the module was built on, and to a guard this repository shipped and has now retracted. An earlier version of this entry published the share as 216 percent, which no experiment can produce -- the bystander term is bounded by every surviving cell, antigen-positive ones included, so dividing it by the antigen-negative pool is not a share of anything.

**P11. The adoptive-cell barriers MULTIPLY, so opening one buys only its own reciprocal.**
- *Quantitative model output:* delivery efficiency is 0.06 as the product of three barriers, of which trafficking is 0.3. Intratumoural rather than intravenous administration sets trafficking to roughly 1 and leaves infiltration and activation untouched, predicting a 3.333x gain in kill -- not the 632.6x that separates the two diseases.
- *Falsification threshold:* the measured gain from intratumoural delivery exceeds 2x the predicted single-barrier factor (delivery gates something the product does not represent), or falls below half of it (the barriers are not independent). Either refutes the multiplicative FORM, which is the only part of this layer a measurement can reach -- every barrier VALUE is a placeholder.

**P12. Whether an oncolytic infection establishes is dose-INDEPENDENT.**
- *Quantitative model output:* simulated, not asserted. Across 5 orders of magnitude of initial infected fraction (1e-06 to 0.1) the verdict is identical and the cumulative lysed fraction spans only 0.0084. Effective replication is 0.63 against a removal rate of 0.35, and it is that comparison -- not the titre -- which decides the outcome. Note the crate's own `spread_threshold_ratio` compares replication against CLEARANCE alone (3.15), so a config above that ratio can still die out once lysis is counted.
- *Falsification threshold:* at fixed tumour permissiveness, raising the input titre flips a non-establishing infection to an establishing one, or moves the cumulative lysed fraction by more than 10 percentage points across the 5 orders of magnitude tested.

**P13. Recurrence after ablation tracks margin GEOMETRY and not delivered energy.**
- *Quantitative model output:* read from `ablation.rs`, not asserted: above threshold `margin_survival_fraction` returns one minus the covered fraction, its signature takes only the config and the coverage, and while its body DOES read temperature, duration and field strength, it reads them only to test whether the threshold is crossed -- above it the return value does not vary with them. An earlier version of this entry said the body read none of them, which was false. Doubling delivered energy at fixed coverage therefore changes predicted survival by exactly zero.
- *Falsification threshold:* recurrence correlates with delivered energy at matched coverage fraction (a dose-response above threshold), or coverage fails to predict recurrence at matched energy.

**P14. A radiotherapy trial's two schedules imply the tissue's fractionation sensitivity, and the implication holds for prostate and fails for breast.**
- *Quantitative model output:* `radiation::isoeffect_alpha_beta` inverts the linear-quadratic model over two schedules a randomised trial reported as not differing. For CHHiP (PMID 27339115), 37 x 2 Gy against 20 x 3 Gy implies alpha/beta = 2.29 Gy, which falls inside the 1.2 to 3.0 Gy band the radiobiology literature estimates for prostate from other data. For START-B (PMID 24055415), 25 x 2 Gy against 15 x 2.67 Gy implies 0.70 Gy against a published 2.3 to 4.5 Gy, and the model therefore FAILS to reproduce that equivalence.
- *Falsification threshold:* the prediction is that the inversion recovers a published site-specific alpha/beta for a trial whose arms were designed to be isoeffective, and does not for a trial whose arms differ in total dose and elapsed time. It fails if a new isoeffective pair from a designed-equivalent trial lands outside its site's published band by more than the band's own width, or if the breast pair is shown to be reproducible under the same single-alpha/beta arithmetic.
- *Note:* a non-inferiority result bounds a difference rather than establishing equality, so the implied ratio is an anchor and not an estimate with an interval. The full accounting, including the repopulation prediction (`ln2/(alpha*T_p)` = 0.77 Gy/day against Withers' published band, PMID 3390344) and its sensitivity to which tumour's alpha is used, is `analysis/calibration/fractionation-validation.md`.

**P15. Shortening a chemotherapy interval helps only inside a window of regrowth rates, and the window has two ends.**
- *Quantitative model output:* `chemo::dose_density_advantage` compares six cycles at 14-day and 21-day intervals at the SAME total dose, with Gompertzian regrowth between cycles. The advantage peaks at about 5.0x near a regrowth rate of 0.0125/day and is worth at least a tenth only between roughly 0.002 and 0.1/day. At zero regrowth it is 0.99x -- below one, because the longer gap lets survivors redistribute into sensitive phases -- and at 0.4/day it is 1.00x, because both schedules' tumours return to the Gompertz plateau between cycles.
- *Falsification threshold:* the prediction is the SHAPE, not the numbers. It fails if a dose-density benefit is demonstrated in a setting whose regrowth kinetics place it outside the window (either a tumour with no measurable regrowth between cycles, or one whose regrowth is fast relative to the interval), or if a within-model change produces a monotone rather than an interior-peaked dependence on regrowth rate.
- *Note:* this is NOT a reproduction of CALGB 9741 (PMID 12668651), which found a shorter interval better in early breast cancer and is the reason the question is asked. Whether breast micrometastatic disease sits inside the window is not determinable from anything in this repository, and no absolute magnitude here is defensible: the potency is a placeholder because the CTRPv2 dose-response route is access-blocked. The hypothesis is Norton and Simon's (PMID 3510732). Full accounting in `analysis/calibration/chemo-validation.md`.

**P16. A checkpoint response ratio between mutational-burden strata is reproducible where the absolute response rate is not.**
- *Quantitative model output:* `checkpoint::response_ratio` between a representative tTMB-high tumour (20 mut/Mb) and a non-high one (3 mut/Mb) returns 3.10x at the shipped shape constants, against KEYNOTE-158's measured 4.83x (PMID 32919526: 29% of 102 against 6% of 688) and a conservative 2.63-7.80x band from the interval endpoints. The ratio is independent of the brake, the occupancy and the unknown response-to-kill mapping, all of which cancel; the absolute index is not comparable to anything.
- *Falsification threshold:* the prediction fails if a stratified response pair from a single trial, drug and endpoint, stratified by a factor this model treats as moving antigenicity alone, gives a ratio outside the model's admissible shape region -- or if a stratification is shown to move the brake as well, in which case the cancellation argument does not hold and the comparison is invalid rather than wrong.
- *Structural companion, which is stronger than the numerical claim:* a tumour that has lost antigen presentation (B2M-null, PMID 27433843) has a response index of exactly zero at any mutational burden, and a ratio against it is undefined rather than infinite. That follows from the layer's structure and not from a fit.
- *Note:* this constrains a SHAPE and identifies nothing. One ratio is one equation and the brake remains unidentified. The representative-burden choice -- which the trial does not make for us -- moves the model's answer from 1.5x to 7.3x, a spread comparable to the width of the target band. Full accounting in `analysis/calibration/checkpoint-validation.md`.

**P17. A CAR-T failure caused by low antigen density cannot be rescued by escalating the dose, and one caused by poor delivery can.**
- *Quantitative model output:* `adoptive::dose_escalation_gain` at matched barriers (`AdoptiveBarriers::solid_tumour`), 1e5 infused effectors against 2e4 tumour cells. With antigen density five times the engagement threshold the gain at a tenfold dose is 10.0x; with density one fifth of the threshold it is 1.00x. The two cases differ in ONE variable and share their barriers, their tumour and their dose.
- *Falsification threshold:* the prediction fails if a dose-escalation series in cohorts matched for infiltration but stratified by target-antigen density shows comparable relative gains in both strata -- specifically, if the low-density stratum's gain exceeds half the high-density stratum's. It also fails if a low-density non-responder is rescued by dose escalation alone at unchanged antigen expression.
- *Note:* the density threshold is a PLACEHOLDER and varies by orders of magnitude with receptor affinity, costimulatory domain and target; the prediction is about the SHAPE of the two responses, not the position of the threshold. The threshold is also sharper here than in a patient, where density varies cell to cell, so a real experiment should see a blend of the two modes rather than either alone. This arm is deliberately NOT fitted to the ELIANA remission rate (PMID 29385370): a remission is not a kill fraction, and unlike the checkpoint arm no ratio is available that cancels the mapping, because blood and solid CAR-T are different trials with different endpoints. Full accounting in `analysis/calibration/adoptive-validation.md`.

**P18. Suppressing immunity to help an oncolytic virus spread is optimal only below a crossover in priming efficiency.**
- *Quantitative model output:* `oncolytic::optimal_immune_competence` over a scan of priming efficiency. Below about 2.6 (this model's units) the optimum sits at zero immune competence; above it the optimum is interior. `oncolytic::priming_efficiency_for_interior_optimum` returns the crossover, and both arms of the trade-off -- `virus_survival` falling with competence, `antitumour_priming` rising -- are individually monotonic, so the non-monotonic sum is emergent rather than imposed.
- *Falsification threshold:* the prediction fails if immunosuppression (or immune-competent versus immunodeficient hosts) changes durable outcome MONOTONICALLY across a range of viruses whose priming efficiency differs substantially -- that is, if no crossover behaviour is observable at all. It also fails if a virus with demonstrably weak anti-tumour priming does worse under immunosuppression, which is the regime this model says suppression should help.
- *Note:* the crossover's POSITION is not a prediction. Priming efficiency is a placeholder in units this model invented and nothing in this repository constrains it; what is claimed is that a crossover exists and what it depends on. OPTiM (PMID 26014293, 16.3% vs 2.1%) anchors the DIRECTION only: those are two different agents, so unlike the checkpoint arm no mapping cancels in a ratio between them. A third of the scanned range saturates and is excluded, because there the optimum's position no longer reflects the trade-off. Full accounting in `analysis/calibration/oncolytic-validation.md`.

**P19. The delivered payload of an antibody-drug conjugate has an interior optimum in drug loading, and it sits near four.**
- *Quantitative model output:* `adc::payload_delivered_per_dose` is loading divided by `adc::clearance_multiplier`, the latter anchored on Hamblett 2004 (PMID 15501986): DAR 8 clears 3x faster than DAR 4 and 5x faster than DAR 2. `adc::optimal_dar` returns 4.0. Both factors are monotonic in loading; the ratio is not, and the peak's position follows from the measured ratios rather than from any choice here. The model also gives delivered(8)/delivered(4) = 0.67 -- twice the payload per antibody delivering LESS -- against the same study's report that E4 and E8 had comparable in-vivo activity at equal antibody dose.
- *Falsification threshold:* the prediction fails if a conjugate series in which clearance is measured across at least three loadings shows delivered payload (or exposure-corrected efficacy at equal antibody dose) rising monotonically with loading. It also fails if the clearance penalty is shown NOT to accelerate -- that is, if `c(4)/c(2)` and `c(8)/c(4)` are equal within measurement error -- because the interior optimum depends entirely on that acceleration.
- *Note:* one conjugate, one payload, one antibody. The optimum at four is a property of this study's molecule and NOT a general claim about ADCs; newer conjugates with more stable linkers are deliberately built at higher loading. Two ratios do not determine a curve, and outside DAR 2 to 8 nothing is measured. Delivered payload is not efficacy. Full accounting in `analysis/calibration/adc-validation.md`.

**P20. Thermal ablation and irreversible electroporation fail in geometrically different places.**
- *Quantitative model output:* `ablation::perivascular_failure_radius_mm` returns a surviving sleeve around a vessel of ~4.4 mm at a 50 C applicator falling to ~0.5 mm at 90 C (5 minutes, CEM43 threshold 240, cooling length 2 mm), and `ablation::electroporation_failure_radius_mm` returns exactly 0 -- structurally, because electroporation is non-thermal and a heat sink removes nothing that matters to it.
- *Falsification threshold:* the prediction fails if local progression after thermal ablation is NOT spatially concentrated near large vessels at matched margin coverage, or if electroporation shows the same perivascular concentration as thermal ablation in a matched series. It is a claim about the GEOMETRY of recurrence, and it is testable on follow-up imaging that records where the recurrence sat rather than only whether one occurred.
- *Note:* the sleeve's SIZE is not predicted. The cooling length is a placeholder standing in for vessel calibre and flow rate, neither of which this layer represents, and the sleeve scales with it almost proportionally. The direction is anchored on the recognised problem of perivascular local progression (PMID 35114665), which supports the effect's existence and not its magnitude. This arm's fitted calibration remains UNCONSTRAINED and this prediction does not repair it: a threshold observable cannot identify a threshold parameter, and the ledger carries both facts. Full accounting in `analysis/calibration/ablation-validation.md`.

**P21. The optimal photodynamic fluence rate falls as the photosensitizer clears more slowly, by more than an order of magnitude across the clinical range of half-lives.**
- *Quantitative model output:* `photosensitizer_pk::optimal_fluence_rate`, holding total fluence at 150 J/cm2 and varying only the rate, returns 315 mW/cm2 for a 15-minute half-life, 66 at 4 hours and 18 at 48 -- a 17-fold span. Two opposing monotonic effects produce it: a faster delivery depletes the tissue oxygen a Type II sensitizer consumes, and a slower one runs further into the drug's own clearance. Every optimum is strictly interior to the scanned 5-400 mW/cm2 range, checked rather than assumed.
- *Falsification threshold:* the prediction fails if a fluence-rate dose-finding series run on two sensitizers whose tissue half-lives differ by at least tenfold shows the SAME optimal rate within a factor of two, or shows the faster-clearing drug preferring the SLOWER rate. It is a claim about the ORDERING across drugs, not about any single rate.
- *Note:* the POSITION of each optimum is not predicted and this is measured rather than conceded. It scales as `phi_crit^0.43` -- roughly the square root of the fluence rate at which photochemical consumption matches perfusive resupply -- and nothing in this repository measures `phi_crit`, so a 20-fold span in it moves the optimum 3.66-fold. The direction is anchored on the measured fluence-rate effect (PMID 16615136); the milliwatt figures are a restatement of an assumption. The oxygen model is quasi-steady and site-local, and the clinic bounds the slow end where the model does not. Full accounting in `analysis/calibration/pdt-fluence-rate-validation.md`.

**P22. The optimal sonodynamic frequency falls as one over the focal depth. THIS PREDICTION IS ALREADY CONTRADICTED, and it is registered anyway.**
- *Quantitative model output:* `sonodynamic::optimal_frequency_mhz` is `10 / (alpha * ln(10) * z)` in closed form, agreeing with a brute scan at every depth tested. At the engine's `sdt_alpha = 0.7` dB/cm/MHz it gives 1241 kHz at 50 mm, 620 at 100 and 414 at 150 -- a 3.0-fold fall.
- *Falsification threshold:* the prediction fails if the frequency maximising delivered cavitation at a fixed acoustic power varies by less than 1.5-fold across a threefold range of focal depth in a matched applicator. **It has already failed against the nearest available comparator:** Ellens & Hynynen 2015 (PMID 26233216), an independent full-wave study of the same three mechanisms, report the same 750 kHz optimum at 50, 100 and 150 mm.
- *Why it is registered.* Two adjacent claims from the same closed form SURVIVE that comparison -- an interior optimum exists, and it falls as attenuation rises (they see 750 kHz drop to 500 when their attenuation doubles; the model predicts an exact halving). Registering only the surviving two would let a reader take the third for granted, and the disagreement is the most informative thing this arm produced. The missing term is one the comparator names itself: near-field heat accumulation, which their conclusion calls the rate limiter for large-volume ablation and which this model has no representation of at all. The comparison is also NOT numeric on either side -- their observable is thermal ablation efficiency where this model's is cavitation likelihood, and absorption enters those with opposite sign; and their attenuation is quoted on pressure amplitude in Np/m/MHz where this engine's is a dB figure on intensity. Only the DIRECTIONS are compared. Full accounting in `analysis/calibration/sonodynamic-validation.md`.

**P23. The oxygen enhancement ratio a spatial run exhibits depends on the steepness of the oxygen gradient, and is largest at an intermediate steepness.**
- *Quantitative model output:* the dose-modifying factor -- the single-fraction dose ratio for equal kill between the oxygenated rim and the hypoxic core of a 60³ tumour -- measured across O₂ gradients, gives 2.27 at λ = 30 µm, **2.55** at 50, 2.36 at 80, 1.67 at 120 and 1.24 at 200. The peak is strictly interior to the swept range. Averaged over kill levels 0.3 to 0.9, across which the factor moves by at most 0.14.
- *Falsification threshold:* the prediction fails if the factor is monotonic in gradient steepness across a threefold range, or if it varies by less than 1.5-fold across it — either would mean the spatial measurement is the Alper–Howard-Flanders formula restated, which contains no gradient term and returns 2.86 regardless.
- *Why it is not circular.* The engine already contains that formula, and `radiation::dna_channel_dose_modifying_factor` is documented and tested as a RESTATEMENT of it. Scoring THAT against the published 2.5–3.0 would be a guard computing its own expectation. What is registered here is a different quantity: a population average over a nonlinear function, read from a grid, that moves with a spatial parameter the formula has no term for.
- *Note:* the mechanism is structural rather than biological. ONE λ sets both the rim and the core, so a gradient steep enough to make the core anoxic leaves the rim at 11.2 mmHg and hypoxic itself, narrowing the contrast from the wrong end. Consequently **every measured factor is a LOWER bound** on the oxic-versus-anoxic ratio the published band describes, and the model reaches that band at one gradient of five. At the engine's own zone reference (120 µm) it gives 1.67, and reporting that as a failure would be a category error: the deep zone there sits near 1.43 mmHg and is not anoxic. Only the DNA channel is wired; the ferroptotic channel's gray-to-ROS conversion is unmeasured and enabling it would make this a function of an unanchored knob. Full accounting in `analysis/calibration/radiation-oer-validation.md`.

### Honesty clause for P9 to P23

**P9 to P13 are DIRECTIONAL and every barrier value behind them is an
uncalibrated placeholder** (`CALIBRATION_STATUS.md` records each as feeding no
reported number). A magnitude in P10 to P13 is a property of a preset, not a
prediction. P9 is the deliberate exception: its band is published, which is
why it is the only one whose model output is quantitative, and it is
correspondingly the easiest of the five to kill.

**P10 is registered against this project's own prior belief.** The module was
built expecting the bystander advantage to GROW as antigen is lost, and shipped
a guard asserting it. The guard was vacuous, the belief was wrong, and the
prediction registered here has the opposite sign. Reporting that is the point
of preregistering at all.

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

### E1. Spheroid RSL3 versus SDT kill at measured hypoxia (tests P4, the most novel test)

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

### E4. GPX4 plus FSP1 dose-matrix synergy (tests P1, **the keystone**)

- **Setup:** RSL3 (or ML162) by iFSP1 (or brequinar) dose matrix in
  persister-enriched cells.
- **Measure:** Chou-Talalay combination index or Bliss; C11-BODIPY to confirm the
  lipid-peroxidation pathway.
- **Model says:** synergy (CI less than 1, about 1.99x Bliss).
- **Falsifies if:** CI greater than 0.8 (additive or antagonistic), revising the
  parallel-independent-repair assumption.
- **Cost:** low (standard combination assay).
- **Why this is the keystone:** because it is the most DECISIVE, not the most
  novel. Both of its outcomes are informative. It rests on GPX4/FSP1 parallel
  independence, which is established (Doll 2019, Bersuker 2019), so a negative
  result falsifies the parallel-independent-repair assumption cleanly rather
  than leaving the mechanism and the execution confounded. It also tests the
  single most quantitative claim the model makes, the ~1.99x Bliss synergy, and
  it is the cheapest experiment in the program to staff.

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

What matters is that no prediction above is set by the calibration data, not that
registration came first — it did not, for four of the legs. For transparency, here
is what was and was not anchored to independent data as of registration. The full per-layer
ledger is `simulations/calibration/CALIBRATION_STATUS.md`.

- **Calibrated (in-vitro, held-out):** the single-cell RSL3 kill switch, fit to
  CTRPv2 GPX4-inhibitor dose-response on ML162 and validated on ML210 (#330).
- **Partially anchored to published data:** spheroid radial zone geometry versus
  Browning 2021 (#333, geometry only); tumor-PK disposition versus imidazole
  ketone erastin plus a sorafenib forward check (#334); Krogh penetration form and
  reference length versus Primeau/Tannock (#335).
- **Prior-predictive only (parameter, not data-conditioned):** the spatial and
  combination headlines that P1, P2, P3, P4, P5, P6, and P7 rest on. An ABC analysis
  (#332) shows the in-vivo priors and the in-vitro posterior are disjoint, so the
  in-vivo and spatial headlines cannot be conditioned on the available in-vitro
  data and stay prior-predictive until an in-vivo ferroptosis dataset that maps onto
  these headline parameters exists (the in-vivo ferroptosis readouts that DO exist —
  e.g. IKE pharmacokinetics and in-vivo SCD1/MUFA data — measure different
  observables that do not condition the headlines).

- **Derived from a calibration artifact, and therefore NOT independent of it:** P8.
  Its numbers (kill near zero below about 280 um radius; 3.6 percent to 0.19 percent
  over 144 to 540 um) are reproduced from `analysis/calibration/spheroid-kill-vs-size.md`,
  which was committed the same day P8 was written. P8 is a statement of what the model
  already computed, not a prediction made in advance of computing it. It is still
  falsifiable by experiment — no ferroptosis-inducer size-kill data exists to test it
  against — but it must not be counted as an independent hit if it later agrees with
  the model.

An earlier version of this section listed six predictions and silently omitted P3 and
P8. P3 is prior-predictive like its neighbours; P8 needed the separate classification
above, which is presumably why it was easier to leave out.

As calibration data arrives, the plan is to report calibrated-versus-preregistered
for each prediction, failures included.

## Literature position of each prediction (measured 2026-08-03)

Added after the census made it measurable. Each prediction chains ferroptosis to
another concept, and MeSH indexes both sides, so the size of the literature each
one lands in is countable across 4,403,994 cancer articles
(`analysis/atlas-prediction-position.md`):

| prediction | ferroptosis-intersecting articles |
|---|---|
| P1, P3 (persister / drug resistance) | 479 |
| P4 (hypoxia) | 64 |
| P5 (immunogenic cell death) | 41 |
| P7 (pH / ion trapping) | 37 |
| P6 (CAF coculture) | 22 |
| P2, P8 (spheroid size) | 15 |

**This is not a quality ranking, and low support is not a reason to drop a
prediction.** A sparse leg is where novelty lives. But it bears on SEQUENCING,
and it settled a disagreement this document once had with
`analysis/p1-wetlab-protocol.md`: this file used to call E1 (testing P4) the
keystone while that protocol called P1 the highest-leverage prediction in the
program, and neither designation had been made with these numbers. Resolved
below (#619).

The relevant asymmetry is interpretability of a NEGATIVE result. P4 additionally
depends on sonodynamic therapy, which has 32 ferroptosis-intersecting articles in
the entire indexed literature (`analysis/atlas-thesis-position.md`), so a
negative P4 could mean the mechanism is wrong or that nobody has yet worked out
how to run the experiment well. P1 rests on GPX4/FSP1 parallel independence,
which is established, so a negative P1 is simply a negative result.

### The decision (#619)

**E4, testing P1, is the keystone, on the basis of DECISIVENESS.** E1 remains in
the program as the most novel test and is not demoted; it is sequenced second.
`analysis/p1-wetlab-protocol.md` defers to this designation rather than making
its own.

Four things decided it, and the first is the one that would reverse the call if
it were wrong:

1. **Interpretability of a negative**, as set out above. A negative E4 is a
   negative result. A negative E1 is ambiguous between "the mechanism is wrong"
   and "nobody has worked out how to run this yet".

2. **That ambiguity is not shrinking.** The sonodynamic leg has sat at roughly a
   quarter of one percent of the ferroptosis field throughout a 25-fold
   expansion of that field, statistically flat in every window pair tested
   (`analysis/atlas-thesis-position.md`). So the hope that the technique matures
   on its own is not supported by the literature's behaviour. The honest limit
   cuts both ways: that test could not have detected a rise below about 4.7x, so
   this is an absence of evidence of movement, not evidence of stasis.

3. **Cost and recruitability.** E4 is a standard combination assay. E1 needs an
   ultrasound rig, a sonosensitizer, spheroid culture, pimonidazole hypoxia
   verification and confocal depth-sectioning. This program's stated bottleneck
   is finding a wet-lab collaborator at all, and naming as keystone the
   experiment that is hardest to staff sequences the whole program behind its
   rarest resource.

4. **P4 tests a claim this project has already softened.** The scientific review
   found the "physical ROS bypasses hypoxia" framing over-extended, since SDT
   ROS is widely oxygen-dependent and the lead clinical agent is Type II, and
   Section 7.1 was rebalanced accordingly. A keystone is a poor use of an
   experiment whose claim has already been partly retracted.

**What this does NOT say.** P4 is not less important, and a sparse leg is not a
weak one. E1 is the most distinctive test in the program and the one whose
positive result would matter most. The two designations genuinely point at
different experiments, and this document is choosing decisiveness over novelty
and saying so, rather than leaving the choice implicit.

This was delegated to the assistant by the owner rather than decided by
measurement; the table above is retained so the reasoning stays auditable and
the call can be revisited.

## Sequencing

1. **Now:** register P1 through P8 above on OSF and mint the DOI, then paste it
   into the registration-status block at the top of this file.
2. **Recruit:** circulate the Part 2 briefs to find a wet-lab collaborator for E4
   FIRST (the keystone, and the cheapest to staff), then E1 (the most novel, and
   the one needing an ultrasound rig). These two issues are tracked as
   `help wanted` (#442 spheroid kill, #448 in-vivo PK).
3. **Calibrate honestly:** as data arrives, report calibrated-versus-preregistered
   for each prediction, failures stated as plainly as confirmations.

## Cross-references

- `analysis/contribution-plan-2026.md` (the working source for Parts 1 and 2)
- `simulations/calibration/CALIBRATION_STATUS.md` (per-layer calibration ledger)
- `MODEL_CARD.md` (intended use, out-of-scope, per-layer status)
- `analysis/research-roadmap-2026.md` (the gap analysis and issue backlog)
