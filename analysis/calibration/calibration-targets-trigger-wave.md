# Calibration target: ferroptotic trigger-wave speed (#482)

## Source

Co, H. K. C., Wu, C.-Y., Lee, Y.-C. & Chen, S.-h. (and colleagues), "Emergence
of large-scale cell death through ferroptotic trigger waves", *Nature*
631:654-662 (2024). **PMID 38987590**, DOI 10.1038/s41586-024-07623-6. Open code:
github.com/imb-lcd/ftw2024 (MATLAB simulation + analysis). Open microscopy source
data: Springer Nature figshare, DOI 10.6084/m9.figshare.25762806.

This is the SPATIAL counterpart of the single-cell bistable switch we already
cross-validate against the same group's work (`analysis/ode-cross-validation.md`,
#344).

## Measured quantity

A ferroptotic lipid-peroxide / ROS front that propagates across a cell monolayer
at a CONSTANT speed, each dying cell igniting its neighbours. The headline,
perturbation-resolved measurements:

| Condition | Front speed (um/min) |
| --- | --- |
| baseline (LPO front 5.83, ROS front 5.46) | **5.52 +/- 0.09** |
| iron chelation (deferoxamine, DFO) | **2.33** |
| iron supplementation (loaded) | **9.40** |

Additional reported behaviour (not all individually committed here): NADPH-oxidase
inhibitors slow the front dose-dependently; erastin drives a monostable ->
bistable ROS transition over 0.39-11.7 uM.

## Why this is a strong calibration target

1. **Quantitative + perturbation-resolved**: a measured speed AND its response to
   a clean iron manipulation (chelation vs loading), so it constrains both the
   absolute speed and the iron dependence.
2. **Open data + open code**: the offline contract holds (we commit the published
   summary values as the derived target; CI never downloads).
3. **Hits a previously-uncalibrated axis**: the `trigger_wave` module is a
   propagating Fenton-iron-driven bistable front, the spatial complement of the
   single-cell switch; its iron dependence is exactly what DFO/loading probe.

## What the model is held to

The `trigger_wave` module (a 1-D bistable Nagumo reaction-diffusion front,
`dL/dt = D*Lxx + k*L*(L-a)*(1-L)`, closed-form speed `c = sqrt(D*k/2)*(1-2a)`,
`k ~ iron`) is validated by `scripts/validate_trigger_wave.py`:

- **baseline speed** lands at the measured 5.52 um/min (a one-point calibration
  of the diffusion-rate product `D*k`);
- the **iron-dose response shape** `c ~ sqrt(iron)` reproduces the measured
  2.33 / 5.52 / 9.40 um/min at iron fold-changes ~0.18 / 1.0 / 2.9 (biologically
  plausible: DFO strips most labile iron, FAC loading multiplies it a few-fold),
  the predicted (not fit) result;
- the numerical solve agrees with the closed form (cross-language self-check);
- GPX4 defense slows/halts the front (direction-only, no matched dataset here).

Committed artifacts: `trigger_wave_measured_data.csv` (this target),
`trigger-wave-validation.{md,json}` (the result). The validator drift-guards the
Python constants against the Rust `trigger_wave.rs` `baseline()`.

## Honest limits

- `D` is absorbed into the one-point `D*k` fit; a first-principles calibration
  would fix `D` from a measured lipid-radical diffusion coefficient (not done).
- The robust contributions are the spatial-front CAPABILITY and the iron-dose
  SHAPE agreement, not the absolute diffusion coefficient.
- The figshare microscopy frames are not re-analysed here; we use the paper's
  published summary speeds as the target (consistent with the other calibration
  legs that anchor to published summary statistics).
