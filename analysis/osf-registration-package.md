# OSF preregistration package — ready-to-file (#543)

> **Purpose.** This is a fill-in transformation of [`PREREGISTRATION.md`](../PREREGISTRATION.md)
> (#497) into the field structure of the **OSF Preregistration** template, so that
> filing reduces to *paste each section and click submit*. Filing on osf.io is the
> only remaining human step; once the DOI is minted, paste it into
> `PREREGISTRATION.md` "Registration status" and `analysis/contribution-plan-2026.md`.
>
> **What this registers.** A set of *directional, pre-calibration* predictions of an
> open ferroptosis simulation engine, plus the wet-lab experiments that would
> falsify them. It is registered BEFORE fitting any layer to the external data
> (#330/#502/#500/#333/#334/#335/#482) so that calibrated-versus-predicted stays
> honest. The predictions are directional (sign of an effect); each carries a
> numeric model output and a pre-stated numeric falsification threshold.

---

## 1. Title
Pre-registered falsifiable predictions of an open ferroptosis tumor-simulation engine, with falsification thresholds and wet-lab test designs.

## 2. Authors
Ezequiel Lares (and any co-author recruited per `analysis/collaborator-brief.md`).

## 3. Description (Study Information → Description)
Most layers of the simulation suite are uncalibrated mechanistic scaffolding (see `MODEL_CARD.md`, `simulations/calibration/CALIBRATION_STATUS.md`). This registration locks in eight directional predictions (P1–P8) of the engine — each with a quantitative model output and a numeric falsification threshold — prior to and independent of the data-anchored calibration legs, so that any later calibration cannot move the goalposts. The aim is falsification-first: the predictions are framed so that a single inexpensive wet-lab experiment (E1–E6) can refute each.

## 4. Hypotheses (Study Information → Hypotheses)
Paste P1–P8 verbatim from `PREREGISTRATION.md` Part 1. Summary (full text + thresholds in that file):
- **P1** GPX4+FSP1 dual inhibition is synergistic in FSP1-low persister-enriched cells (Bliss ~1.99×, 95% PPI ~1.0–5.2×). *Falsified if* Chou-Talalay CI > 0.8, or combined kill ≤ Bliss-independence within assay error.
- **P2** Physical-ROS (PDT/SDT) less depth-limited than systemic RSL3 in ≥500 µm spheroids. *Falsified if* RSL3 core/rim within 1.5× of SDT core/rim, or all modalities' depth half-distances agree within 25%.
- **P3** Post-withdrawal vulnerability window closes over days, defenses recover sequentially (FSP1/GSH before GPX4/NRF2). *Falsified if* sensitivity returns within 24 h, or all defenses recover at the same timepoint.
- **P4** SDT retains more efficacy than RSL3 under hypoxia (direction only; least-certain). *Falsified if* SDT hypoxic-kill loss ≥ RSL3's (ratio ≤ 1.0).
- **P5** Dense ferroptotic kill yields more ICD signal per dead cell than sparse (~4:1 in 3D; can invert under the immunosuppressive arm). *Falsified if* DAMP/DC-maturation per dead cell agree within 1.5× between dense and sparse.
- **P6** CAF coculture protects RSL3 more than SDT. *Falsified if* RSL3 and SDT IC50 fold-shifts agree within 1.5×, or neither shifts > 1.2×.
- **P7** RSL3 efficacy drops at acidic pH (6.5 vs 7.4) via ion trapping (least-certain). *Falsified if* efficacy/IC50 agree within 1.2× across pH.
- **P8** A persister-targeting inducer (RSL3) has the OPPOSITE size-dependence to generic cytotoxics (near-zero kill < ~280 µm, rising as the persister core appears). *Falsified if* monotone bigger-resists-more, or size-independent.

**Honesty clause (paste).** P4 and P7 are flagged least-certain; this registration commits to reporting their failures as prominently as any success.

## 5. Design Plan (Design Plan → Study type / Blinding / Study design)
- **Study type:** Registration of falsifiable computational-model predictions, to be tested by independent wet-lab experiments (E1–E6) described in §8. This is not a re-analysis of existing data; the model outputs are fixed at registration.
- **Blinding:** N/A for the model outputs (deterministic, seed-fixed, byte-identity-gated). For the wet-lab tests, analysts scoring viability / lipid-peroxidation / DC-maturation should be blinded to arm where feasible.
- **Study design:** Each prediction maps to a controlled in-vitro comparison (dose-matrix, ± modality, ± CAF, ± hypoxia, pH, spheroid size, or withdrawal time-course) with a pre-stated falsifying outcome.

## 6. Sampling Plan (Sampling Plan → Existing data / Data collection / Sample size / Stopping rule)
- **Existing data:** Registration prior to data collection for the wet-lab tests. The model side uses only already-public datasets for the separate calibration legs (CTRPv2, Browning 2021, IKE PK, Primeau/Tannock, Co 2024), which are NOT used to set these predictions.
- **Data collection / sample size:** Per the E-series (§8): standard dose-response/checkerboard designs with biological replicates ≥ 3 and the power to resolve the pre-stated fold thresholds (e.g. a 1.5× IC50 shift). Final n is set by the executing lab and recorded before unblinding.
- **Stopping rule:** Fixed replicate count per E-series design; no optional stopping.

## 7. Variables (Variables → Manipulated / Measured / Indices)
- **Manipulated:** drug identity (RSL3/ML162 vs iFSP1/brequinar vs SDT/PDT), dose, O₂ tension (21% vs 1%), pH (7.4 vs 6.5), ± CAF coculture, spheroid radius, time since drug withdrawal.
- **Measured:** viability / IC50, C11-BODIPY lipid peroxidation, GSH, surface calreticulin / HMGB1 / ATP (DAMPs), DC maturation / cross-presentation, defense-protein levels (FSP1/GSH/GPX4/NRF2), intracellular RSL3 (HPLC), pO₂ (pimonidazole).
- **Indices:** Chou-Talalay combination index / Bliss excess (P1); core-to-rim kill ratio (P2); IC50 fold-shift (P6); hypoxic-kill loss ratio (P4); DAMP/DC-maturation per dead cell (P5).

## 8. Analysis Plan (Analysis Plan → Statistical models / Inference criteria / Exclusions)
- **Inference / falsification criteria:** the numeric thresholds in §4 are the pre-stated decision rules. A prediction is *refuted* if its threshold is crossed at the matched-effect comparison; *supported* (directionally) otherwise. Magnitudes are NOT claimed — only the sign and the threshold.
- **Statistical models:** dose-response fits (4-parameter logistic) for IC50/EC50; Chou-Talalay / Bliss for synergy; per-zone kill fractions for spatial predictions; two-group comparisons at matched effect for the ratio thresholds.
- **Calibrated-vs-preregistered reporting:** as each calibration leg lands, the calibrated value is reported against the preregistered prediction, failures included.
- **The experiment briefs (E1–E6):** paste from `PREREGISTRATION.md` Part 2 (E1 keystone hypoxia spheroid → P4; E2 CAF IC50 → P6; E3 DAMP/DC → P5; E4 GPX4+FSP1 matrix → P1; E5 sequential recovery → P3; E6 pH → P7).

## 9. Other (Other → Prior work / Conflicts)
- **Prior work:** the manuscript (`article/drafts/v1.md` Ch.6–7, 9), `MODEL_CARD.md`, `simulations/calibration/CALIBRATION_STATUS.md`, and `analysis/identifiability-report.md` (which states no headline magnitude is point-estimable). The two-paper split outline is `analysis/manuscript-split-outline.md`.
- **Conflicts of interest:** none declared; the work is MIT-licensed and non-commercial.

---

### Filing checklist (the human step)
1. Create/log in at osf.io → New Registration → "OSF Preregistration".
2. Paste sections 1–9 above into the matching fields (P1–P8 and E1–E6 verbatim from `PREREGISTRATION.md`).
3. Submit to mint the immutable DOI.
4. Paste the DOI + public URL + date into `PREREGISTRATION.md` "Registration status" and `analysis/contribution-plan-2026.md`.
