# An embeddable ferroptosis simulation engine with a data-anchored validation core and byte-identity honesty discipline

**Authors:** Ezequiel Lares (lead); _open for a co-author with ferroptosis/redox or computational-oncology expertise — see [`analysis/collaborator-brief.md`](../../analysis/collaborator-brief.md)._

> **STATUS: SKELETON — not a submission.** This is the Paper B scaffold called for by #523, operationalizing the two-paper split outline ([`analysis/manuscript-split-outline.md`](../../analysis/manuscript-split-outline.md), #504). Section bodies are content notes + source pointers, not finished prose. **Venue and length targets below are provisional placeholders the maintainer can change.**

**Provisional target:** a computational / systems-biology *methods* venue that publishes reusable tools with explicit validation-status reporting. Preprint first (a quantitative-biology category).

**Length target:** methods article ~4,000–6,000 words main text + a software/figure supplement.

---

## Lead claim (the one sentence this paper defends)

A small, reusable, byte-identity-disciplined ferroptosis biochemistry engine reproduces **several independently measured** ferroptosis phenomena — an in-vitro GPX4-inhibitor dose-response, spheroid zone geometry, a paired plasma/tumor PK course, a drug-penetration length, and a propagating ferroptotic trigger-wave speed — and the **engineering discipline that keeps those anchored results honest** (a model card, a per-layer calibration ledger, a calibration-regression CI gate, and a byte-identity containment gate) is itself a contribution.

## What this paper is NOT

Not a paper of validated quantitative *clinical* predictions. The dozens of off-by-default tumor-microenvironment realism layers are **explicitly not the headline**: they live in a single supplementary catalogue, labeled directional / uncalibrated, and are used in **zero** reported numbers.

---

## Abstract (draft stub, ~210 words — replace)

> Mechanistic tumor simulations are most useful when their validation status is legible — when a reader can tell which numbers are anchored to data and which are illustrative. We present `ferroptosis-core`, a small embeddable ferroptosis biochemistry engine (Rust, with C-FFI and Python bindings) built around a bistable recover-or-collapse lipid-peroxidation switch, and we lead with **only** the legs tied to independent published data, each reported with a held-out or out-of-fit metric: an in-vitro GPX4-inhibitor dose-response (CTRPv2; ML162 fit RMSE 0.05, held-out *compound* ML210 RMSE 0.07), spheroid zone geometry (Browning 2021 confocal; size-aware thresholds, mean boundary error 0.025 vs 0.338), a paired plasma/tumor PK course (imidazole-ketone-erastin; partition Kp 0.90 and a multi-compartment finding the summary NCA proves), a Krogh drug-penetration length (within the measured 25–75 µm), and a propagating ferroptotic trigger-wave speed (Co 2024; 5.52 µm/min baseline, c∝√iron). We then state plainly that **no headline magnitude is point-estimable** under the documented parameter uncertainty, and we describe the honesty infrastructure — a model card, a per-layer calibration ledger, a calibrate-or-cut policy, a calibration-regression CI gate, and a byte-identity containment gate — as a reusable pattern for broad-but-honest mechanistic modeling. The dozens of uncalibrated microenvironment realism layers are demoted to a single supplementary catalogue used in no reported number.

---

## Section outline (skeleton)

### 1. Introduction
- A reusable engine is useful only if its validation status is legible.
- The contribution = the tool **plus** the honesty infrastructure, not a set of clinical numbers.

### 2. Engine design
- The single-cell biochemical core (bistable recover-or-collapse switch); System Xc-/GSH/GPX4 + FSP1/CoQ10 backup; the System Xc-/erastin inducer leg (#502).
- FFI / Python bindings (PyO3, cdylib); the off-by-default byte-identical discipline (every realism layer ships inert).
- _Source: `simulations/ferroptosis-core/`, v1.md Ch.5._

### 3. The data-anchored validation core — the headline
Lead with ONLY the legs tied to independent data, each with its metric:
- **#330** in-vitro kill switch vs CTRPv2 GPX4-inhibitor dose-response (ML162 fit RMSE 0.05, ML210 held-out-**compound** RMSE 0.07; same dataset/mechanism scope stated).
- **#333** spheroid zone geometry vs Browning 2021 confocal (size-aware thresholds, mean boundary error 0.025 vs 0.338 fixed).
- **#334** tumor PK vs IKE paired plasma/tumor course (Kp 0.90; the e·Tmax floor / multi-compartment finding, reached in the ka→ke limit).
- **#335** Krogh penetration form + reference length vs Primeau/Tannock (λ 50 µm within measured 25–75 µm).
- **#464** PDT/SDT source-independent singlet-oxygen kill threshold vs Zhu 2015 (~0.5 mM).
- **#482** ferroptotic trigger-wave speed vs Co 2024 (5.52 µm/min baseline; c∝√iron). **State the honest scope:** the baseline is a one-point D·k calibration and the iron-dose folds are back-solved (see `iron_dose_shape_independently_validated=false`), so the response *shape* is not independently validated.
- Cross-validation vs published ODE models (#344) + the Arbatskiy 73-ODE structural cross-check (#471).
- _Source: `analysis/calibration/*.md`, `tests/test_calibration_regression.py`._

### 4. Identifiability and uncertainty
- No headline magnitude is point-estimable under the documented parameter uncertainty (`analysis/identifiability-report.md`, #503): 11 swept parameters, 6 non-identifiable from the kill rate, prior-predictive intervals, the in-vitro/in-vivo ABC disjunction (#332/#500).

### 5. Honesty infrastructure — the second half of the contribution
- The model card (`MODEL_CARD.md`), `CALIBRATION_STATUS.md` per-layer ledger, the calibrate-or-cut policy (#501), the calibration-regression CI gate (#499), and the byte-identity containment gate (#253) as a reusable "broad but honest" pattern.

### 6. Discussion
- What an anchored-core + uncalibrated-periphery tool is good for (hypothesis generation, falsifiable predictions, preregistration — `PREREGISTRATION.md`) and what it cannot do (stand in for in-vivo dose-response).

### 7. Supplementary catalogue
- **ONE table** listing every off-by-default realism layer with its tier, named calibration target, and `used-in-any-reported-number = N` status. The layer catalogue lives here, demoted from the body.

## Figures (provisional)
- F1: Engine architecture (core → FFI/Python → realism layers, all off-by-default).
- F2: The six anchored validation legs, each measured-vs-model.
- F3: Identifiability summary (point-estimable vs directional-only per headline).
- F4: The honesty-infrastructure gate diagram (calibration-regression + byte-identity).

## What moves OUT of Paper B
All corpus/literature content (→ Paper A). The uncalibrated layer catalogue → one supplementary table. The book's narrative chapters → background.

## Source map (book → this paper)
| v1.md content | Section here |
|---|---|
| Ch.5 single-cell switch + identifiability | §3 + §4 |
| Ch.6–7 penetration, combinations, TME layers | anchored legs → §3; the rest → §7 catalogue |
| Ch.8.4 structural limitations | §6 discussion |
| 30+ off-by-default realism layers | §7 — ONE supplementary table |
