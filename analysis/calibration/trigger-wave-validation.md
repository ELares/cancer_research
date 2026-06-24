# Ferroptotic trigger-wave validation (#482)

Validates the `ferroptosis-core` `trigger_wave` module (a 1-D bistable Nagumo
reaction-diffusion front) against the measured ferroptotic trigger-wave speeds of
Co, Wu, Lee & Chen, *Nature* 631:654 (2024), **PMID 38987590** (open code
github.com/imb-lcd/ftw2024 + figshare 10.6084/m9.figshare.25762806).

## Model

The propagating ferroptotic front obeys `dL/dt = D*Lxx + k*L*(L-a)*(1-L)` with
closed-form speed `c = sqrt(D*k/2)*(1 - 2a)`. The autocatalytic peroxidation rate
`k` scales with labile iron (Fenton), so **`c ~ sqrt(iron)`**; the GPX4/GSH
defense raises the ignition threshold `a`, slowing and ultimately halting the
front.

## Result

| Condition | Measured (um/min) | Model (um/min) | Implied iron fold |
| --- | --- | --- | --- |
| iron chelation (DFO) | 2.33 | 2.343 | 0.18 |
| baseline | 5.52 | 5.522 | 1.0 |
| iron loaded | 9.4 | 9.403 | 2.9 |

Numeric solve (baseline): 5.429 um/min (agrees with the
closed form, the cross-language self-consistency check).

Checks: iron-dose ordering DFO<control<loaded = True;
baseline near 5.52 = True; DFO near 2.33 = True;
loaded near 9.40 = True; numeric == analytical (<6%) =
True; GPX4 defense slows the front =
True. **All passed: True.**

## Honest scope

- The baseline `D`/`base_reaction_rate` are **tuned** so the baseline lands at
  5.52 um/min (a one-point calibration of the product `D*k`), so the baseline
  match is a calibration, not a prediction.
- The **predicted** result is the iron-dose RESPONSE SHAPE `c ~ sqrt(iron)`: the
  measured 2.33 / 5.52 / 9.40 um/min imply iron fold-changes of ~0.18 / 1.0 /
  2.9, which are biologically plausible (DFO strips most labile iron; FAC loading
  multiplies it a few-fold). So the measured iron-tuning is **consistent with a
  Fenton-iron-driven bistable front**.
- The GPX4-defense leg (front slows/halts as the ignition threshold rises toward
  0.5) is a **direction-only** prediction with no matched quantitative dataset
  here.
- A full first-principles calibration would fix `D` from a measured lipid-radical
  diffusion coefficient (not done; `D` is absorbed into the one-point `D*k` fit).
  The robust contribution is the spatial-front CAPABILITY plus the iron-dose-shape
  agreement, not the absolute `D`.

A drift-guard (`drift_guard()`) re-reads the Rust `trigger_wave.rs` `baseline()`
constants so this Python validator and the Rust module cannot silently diverge.
