# Model Card: ferroptosis-core / sim-tme-3d

A concise intended-use and limitations summary for the simulation suite in
this repository. It consolidates the per-layer accounting in
[`simulations/calibration/CALIBRATION_STATUS.md`](simulations/calibration/CALIBRATION_STATUS.md)
and the manuscript honesty section (`article/drafts/v1.md` Section 8.4) into the
structured "model card" format used in computational biology, with an explicit
assumptions/scope checklist (ARRIVE-style for the in-silico experiments,
TRIPOD-style for any predictive framing, adapting the TRIPOD-AI checklist even
though this is a mechanistic model, not an AI predictor).

If you read only one thing: **this suite is broad but mostly uncalibrated. It is
mechanistic scaffolding that shows the shape of an effect (direction,
qualitative behavior, spatial structure), not a calibrated quantitative
predictor. Breadth of coverage is not evidence of depth.**

---

## 1. Model details

| Field | Value |
|---|---|
| Name | `ferroptosis-core` (embeddable engine) + `sim-tme-3d` (3D capstone binary) and 10 sibling binaries |
| Version | ferroptosis-core 0.66.0 (see `simulations/ferroptosis-core/Cargo.toml` for the authoritative current value) |
| Author | Ezequiel Lares |
| License | MIT |
| Language | Rust (core + binaries), Python bindings (PyO3), C FFI (PhysiCell-style ABI) |
| Type | Mechanistic, stochastic (Monte Carlo) single-cell ferroptosis ODE engine + 2D/3D spatial tumor-microenvironment layers. NOT a statistical / machine-learning predictor. |
| What it is NOT | Not a trained model, not fit to a patient dataset, not a clinical decision tool, not a validated biomarker. |

## 2. Intended use

**Primary intended use:** mechanistic hypothesis exploration. Take a specific,
mechanistically-stated claim about ferroptosis-based therapy (for example "dual
GPX4 + FSP1 inhibition is synergistic," or "physical-ROS modalities are less
sensitive than RSL3 to the hypoxia / stromal / pH barriers") and test whether it
is internally consistent and what direction and rough magnitude the model
produces, so the claim can be turned into a falsifiable wet-lab experiment.

**Intended users:** researchers, engineers, and students with enough domain
context to read the caveats. Every quantitative output should be read alongside
[`CALIBRATION_STATUS.md`](simulations/calibration/CALIBRATION_STATUS.md).

**In-scope outputs (read as directional / order-of-magnitude):**

- Relative kill-rate differences between treatments / phenotypes / microenvironment conditions.
- The direction and qualitative shape of an effect (e.g. a barrier collapses RSL3 more than SDT).
- Spatial structure (radial gradients, depth-kill curves, kill-density patterns).
- Whether a mechanistic combination is plausibly synergistic or antagonistic.

## 3. Out-of-scope / do-not-use

Do **NOT** use this model for, or cite it as evidence of, any of the following:

- **Patient-specific or clinical decisions.** No part of this is calibrated to or validated against patient data.
- **Precise quantitative predictions.** The specific numbers in the manuscript (1.99x Bliss synergy, 53% pH-driven reduction, 104:1 immune ratio, depth-kill percentages) are explicitly order-of-magnitude estimates from the 2D engine with estimated parameters; they are predictions to test, not measurements.
- **Dosing, scheduling, or safety guidance.**
- **Absence claims about biology.** The model not producing an effect is not evidence the effect does not exist; the model only includes the mechanisms that were coded.
- **Any uncalibrated 3D realism layer as a quantitative result.** The 3D `sim-tme-3d` layers (vasculature, slab, spheroid, clonal, persister + locking, suppressor, multi-checkpoint, contact, nutrient, senescence, dynamic iron, immunosuppressive ferroptosis, ...) are off-by-default, uncalibrated, and deliberately excluded from every quantitative claim in the manuscript.

## 4. Inputs and outputs

- **Inputs:** cell phenotype parameters (`Cell`: iron, GSH, GPX4, FSP1, basal ROS, lipid unsaturation, NRF2), biochemistry rate constants (`Params`), a treatment (Control / RSL3 / SDT / PDT), and optional spatial / microenvironment configs (oxygen, pH, stromal, immune, drug-transport, and the off-by-default realism layers). A fixed RNG seed makes every run deterministic and reproducible.
- **Outputs:** per-cell death / final lipid-peroxide / GSH / GPX4 state; aggregate kill rates; spatial fields (DAMP, O2, supply); and a `summary.json` for `sim-tme-3d`. The default 24-condition `sim-tme-3d` matrix is byte-identical run to run (SHA-guarded in CI), so off-by-default layers cannot silently perturb the headline output.

## 5. Architecture (brief)

A single-cell ferroptosis ODE engine (`biochem`): `total_ros = basal_ros + exogenous + Fenton`; antioxidant quench from GPX4/GSH + FSP1 (+ optional DHODH, GCH1/BH4); autocatalytic lipid-peroxidation propagation gated into a bistable recover-or-collapse switch; death when lipid peroxide crosses a threshold. The single-cell engine is composed into 2D and 3D spatial grids with microenvironment fields. Full module list and current unit-test count: [`simulations/ferroptosis-core/README.md`](simulations/ferroptosis-core/README.md).

## 6. Assumptions and scope checklist

### 6.1 Core modeling assumptions

- [x] Ferroptosis is driven by iron-dependent lipid peroxidation overwhelming GPX4/FSP1 (and optional DHODH/GCH1) repair. Upstream signaling (System xc-, ACSL4 remodeling, full iron import/export dynamics) is represented only partially (e.g. the `contact`, `nutrient`, dynamic-iron layers) or not at all.
- [x] A single treatment window (180 steps). **No adaptive evolution** across cycles in the 2D analyses that produce the manuscript numbers; the off-by-default `clonal` and `persister` layers begin to model evolution/tolerance but are uncalibrated.
- [x] Stochastic single-cell variation via per-cell parameter sampling; spatial coupling via diffusion fields.
- [x] Determinism: identical seed -> identical output (position-independent per-cell RNG).

### 6.2 In-silico experiment reporting (ARRIVE-style)

- [x] **Objective stated** per experiment (manuscript Chapters 6-7; each has a falsification criterion in Chapter 9).
- [x] **Model + parameters documented**: `parameter_provenance.md` (per-parameter, with `Grounded?` column) and source code (authoritative for defaults).
- [x] **Randomization / seeds**: fixed, documented seeds; runs reproducible.
- [x] **Replication / reproducibility**: pinned Rust toolchain (`simulations/rust-toolchain.toml`, 1.96.0) and Python lockfile; CI re-runs the suite and a production byte-identity SHA.
- [x] **No silent truncation**: off-by-default layers are identity (byte-identical) until explicitly enabled.

### 6.3 Predictive-claim reporting (TRIPOD-style, where the work is framed predictively)

- [x] **Intended use and population**: in-silico mechanistic exploration, NOT a patient-outcome model. No human-subjects predictors or outcomes.
- [x] **Predictors / outcome**: biochemistry parameters -> simulated cell death; not a fitted risk score.
- [x] **Validation status**: see Section 7. There is **no external clinical validation**; the model is not a clinical prediction model and should not be reported as one.
- [x] **Limitations stated**: Section 3 (out-of-scope) and the manuscript Section 8.4.

## 7. Calibration and validation status

### 7.1 Calibration tiers (per layer)

The authoritative, per-layer table lives in
[`CALIBRATION_STATUS.md`](simulations/calibration/CALIBRATION_STATUS.md). The
tiers are:

| Tier | Meaning |
|---|---|
| **Calibrated** | At least one key parameter fit/anchored to an independent published measurement. |
| **Partially anchored** | Structure or a parameter subset has a literature basis, but the headline magnitude is not fit to data. |
| **Self-consistency only** | A regression-guard that the model reproduces its own hard-coded physics; not an independent calibration. |
| **Uncalibrated (illustrative)** | Placeholder parameters chosen for plausible behavior; the mechanism/direction is the claim, the magnitude is not. |

**Distribution, stated plainly:** a small core (depth physics, parts of the
single-cell engine) is anchored or sensitivity-tested; the bulk of the 3D
realism layers are **uncalibrated (illustrative)**, and the three 3D "validation"
targets are **self-consistency** checks, not independent data. None of the
uncalibrated layers feed the manuscript's quantitative claims.

### 7.2 What HAS been checked

- **Sensitivity / robustness:** the headline qualitative result (Persister > Glycolytic vulnerability under SDT) held under +/-50% perturbation on 11 rate constants (22/22 conditions); weight-, taxonomy-, and PRCC global sensitivity analyses were pre-registered and run.
- **Physics self-consistency:** PDT optical attenuation and SDT acoustic penetration reproduce published tissue optics/acoustics constants.
- **Determinism / regression:** golden kill-count tests and a full-production `summary.json` SHA gate in CI.
- **In-vitro kill-switch calibration (#330):** the single-cell RSL3 kill switch is fit to CTRPv2 GPX4-inhibitor median dose-response (fit on ML162 RMSE 0.05, held-out validated on ML210 RMSE 0.07; the default in-vivo-tuned switch is ~11x worse, far too RSL3-resistant for in-vitro data). This anchors the single-cell switch to an independent in-vitro dataset with held-out comparison; it does NOT calibrate the in-vivo / spatial layers (see `analysis/calibration/kill-switch-calibration.md`).
- **Spheroid zone-geometry validation (#333):** the spheroid module's fixed radial zone thresholds were checked against the Browning 2021 size-resolved confocal dataset (994 spheroids); they match the data only for radius >= 300 µm. The size-aware refinement (`SizeAwareZones`, opt-in/byte-identical) now ramps the thresholds with size so small spheroids are correctly all-proliferating with no core (size-aware mean boundary error 0.025 vs the fixed 0.338 on the committed bins), reducing to the validated fixed structure at large radius. The per-zone biochem gradient strengths remain uncalibrated, and the size-resolved KILL leg is data-blocked (no ferroptosis-inducer size-kill data) (see `analysis/calibration/spheroid-structure-validation.md`).
- **Tumor-PK measured-data anchor (#334):** the `tumor_pk` plasma+tumor disposition is anchored to imidazole ketone erastin (the only public ferroptosis-specific in-vivo PK with a paired plasma+tumor course) plus a sorafenib human popPK forward check. Measured tissue:plasma partition Kp=0.90 (vs the presets' estimated 0.15 to 0.5), a ~2 h plasma to tumor delay, and a closed-form proof that the disposition is multi-compartment. The per-tumor presets are NOT recalibrated (no public per-tumor PK), and the RSL3/ML tool-compound PK gap is flagged not fabricated (see `analysis/calibration/pk-calibration.md`).
- **Krogh penetration validation (#335):** the `drug_transport` exponential penetration form is validated against measured in-vivo data (Primeau 2005 and Tannock 2002 both report exponential drug-vs-distance-from-vessel decline), and the doxorubicin transport reference (model half-distance 34.7 µm) sits within Tannock 2002's measured 25 to 75 µm and ~13% below Primeau 2005 / Minchinton 2006 (40 to 80 µm). NOT validated: the RSL3-like penetration length (no ferroptosis-inducer penetration data exists). The dose-dependent binding-site barrier is deliberately NOT added: a #335 data-availability review found it is an antibody phenomenon (Fujimori 1990, Saga 1995, Thurber/Wittrup 2008), physically weak for the small molecules this model targets (deep binding sink barely saturates, El-Kareh/Secomb), and unvalidatable from public small-molecule data, so the clause is resolved by evidence rather than by building unvalidatable complexity. See `analysis/calibration/penetration-validation.md`.
- **Corpus evidence tagger (separate from the simulation):** gold-set evaluation measured **46% exact / 96% precision / 55% recall**. The low recall is why corpus absence claims are reported as provisional ("not detected in the local keyword-derived analysis"), NOT as definitive.

### 7.3 What has NOT been validated

- Four layers are now anchored to or validated against independent published data: the single-cell RSL3 kill switch (in-vitro CTRPv2 GPX4 inhibitors with held-out comparison, #330), the spheroid zone geometry (Browning 2021 confocal structure, #333, geometry only), the tumor-PK partition + disposition structure (IKE/sorafenib, #334, structure + partition not per-tumor magnitudes), and the Krogh penetration form + reference-drug length (Primeau/Tannock, #335, form + doxorubicin λ only). The in-vivo kill magnitudes, the spatial immune coupling, the per-tumor PK presets, the ferroptosis-inducer penetration length, the dose-dependent binding-site barrier, and the 3D-realism layers are NOT validated against any independent ferroptosis dataset.
- No clinical / patient-outcome validation of any kind.
- The 3D self-consistency targets check the model against itself, not against measurements.

## 8. Known failure modes and limitations

- **Magnitudes are not trustworthy.** Read direction and order of magnitude only.
- **The hypoxia / SDT leg is contested.** SDT is modeled as oxygen-independent (optimistic upper bound); the lead clinical agent is oxygen-dependent. An off-by-default oxygen-dependent SDT mode and a dynamic-iron hypoxia coupling exist to test the reverse, but the magnitude of the SDT-vs-RSL3 hypoxia gap is unresolved.
- **The immune-coupling ratio is geometry-sensitive** (104:1 in 2D shrinks to roughly 4:1 in 3D), and the net sign can flip once the off-by-default immunosuppressive-ferroptosis arm is enabled.
- **Several directions are genuinely contested in the literature** and the model encodes them as configurable/bidirectional rather than single-signed (e.g. ether-lipid plasmalogen sub-step, nutrient stress, and senescence, which is a senolytic target under direct GPX4 inhibition but resistant to upstream triggers).
- **Corpus limitations:** open-access skew, missing landmark full text, taxonomy-dependent gap counts, 55% tagger recall.

## 9. Ethical considerations

This is open, MIT-licensed research released so others can validate, falsify, or
build on it. It makes **no clinical claims** and must not be used to guide
patient care. The repository's guiding principles (`CLAUDE.md`) require
evidence-led, honestly-caveated reporting and explicitly warn against treating
broad mechanism coverage as therapeutic depth.

## 10. How to cite and reproduce

- Citation metadata: `CITATION.cff`. Provenance/redistribution: `PROVENANCE.yaml`.
- Reproduce: pinned `simulations/rust-toolchain.toml` (1.96.0) + `requirements-lock.txt`; `cargo test --workspace` and the Python pipeline are re-run in CI.
- Figure traceability: `FIGURES.yaml`. Per-layer calibration: `CALIBRATION_STATUS.md`. Manuscript: `article/drafts/v1.{md,tex,pdf}` (the honesty section is 8.4).

## 11. Maintenance

This card is handwritten interpretation, kept in sync by hand. When a layer's
calibration tier changes, update both this card and `CALIBRATION_STATUS.md`. The
machine-checked authorities (source code for defaults, `parameter_provenance.md`,
`targets.yaml`) take precedence over prose if they ever disagree.
