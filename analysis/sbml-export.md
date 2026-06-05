# SBML export of the ferroptosis-core biochemistry (#351)

Interoperability is how a model gets reused. The repo already ships a PhysiCell
C-FFI; this adds an **SBML** export of the core single-cell ferroptosis ODE
network so others can import it into COPASI, Tellurium / roadrunner,
PhysiCell-style frameworks, or the BioModels repository.

- Exported model: `simulations/sbml/ferroptosis_core.xml` (SBML Level 3 Version 2).
- Generator + round-trip validator: `simulations/sbml/export_ferroptosis_sbml.py`.
- Round-trip figure: `simulations/sbml/sbml_roundtrip.png`.

## What is exported

The **deterministic mean-field** of the core single-cell ferroptosis ODE network
in `ferroptosis-core/src/biochem.rs::sim_cell_step`, under the 2D default
parameters with **RSL3** (a direct GPX4 inhibitor) applied at t=0:

- **State variables:** `LP` (lipid peroxide), `GSH`, `GPX4`.
- **Rate laws** (verbatim from the engine, stochastic noise replaced by its mean):
  iron/Fenton ROS production; Michaelis-Menten GSH scavenging + NRF2-driven GSH
  resynthesis; the autocatalytic LP propagation gated by the GSH/GPX4
  antioxidant quench (`lp_propagation / (1 + quench*5)`, the bistable-switch
  term) plus GPX4/FSP1 repair; and GPX4 dynamic regulation (ROS degradation +
  NRF2 upregulation).
- **Parameters** carry the 2D default values; provenance is in
  `simulations/calibration/parameter_provenance.md`.

## What is NOT exported (and why)

SBML/COPASI users get the deterministic single-cell core. These parts of the
suite are deliberately out of scope for the SBML export:

| Not exported | Why |
|---|---|
| Per-step stochastic noise (`norm(...)` on Fenton/exo/LP) | SBML is the deterministic mean-field; the noise mean is 1.0/0.0, so the SBML is the expectation. The full stochastic engine stays in Rust. |
| The death threshold + post-death LP accumulation | A discrete event/cutoff, not part of the continuous repair-vs-generation ODE. |
| Per-cell parameter sampling (`gen_cell`) | The SBML is one representative (OXPHOS-mean) cell; population variation lives in the Monte Carlo engine. |
| MUFA / ether-lipid / dynamic-iron / persister / senescence / etc. realism layers | Off-by-default, uncalibrated; they would each add species/parameters and are tracked separately. The export is the calibrated-status core, not the scaffolding. |
| SDT/PDT time-varying exogenous-ROS envelope | RSL3 (exo = 0) gives an autonomous ODE; the SDT/PDT bolus + decay envelope is time-dependent and would need an SBML event/assignment per dose. |
| All 2D/3D spatial fields (oxygen, pH, immune, vasculature, ...) | The SBML is a single well-mixed cell; spatial coupling is a separate modeling layer. |

## Round-trip validation

The generator does a genuine round-trip with an independent SBML toolchain:

1. Builds the model in **Antimony**, converts to **SBML L3V2**, and validates it
   with **libSBML**: **0 errors, 0 fatal**.
2. Reloads the SBML and simulates it with **roadrunner** (Tellurium's CVODE
   integrator) over 180 steps.
3. Compares against a forward-Euler `dt=1` reference that mirrors the discrete
   `biochem.rs` engine, and against the **actual** `ferroptosis-core` engine's
   mean final state (via the Python bindings).

Result (RSL3 on an OXPHOS-mean cell), final state:

| | LP | GSH | GPX4 |
|---|---|---|---|
| exported SBML (roadrunner, continuous) | 0.211 | 0.523 | 0.935 |
| engine ODE (Euler dt=1 reference) | 0.210 | 0.525 | 0.936 |
| ferroptosis-core engine (stochastic mean, n=4000) | 0.462 | 0.592 | 0.935 |

- **The exported SBML reproduces the intended ODE.** All three state variables
  track the deterministic engine closely: GSH and GPX4 are essentially exact
  (the curves overlap in the figure), and the SBML-vs-Euler-reference final-LP
  difference is ~0.3%. This is the round-trip integrity check: a third-party
  SBML simulator reproduces the engine's deterministic dynamics. (The GPX4
  ROS-degradation term is gated behind `total_ros > 1` in the engine; the SBML
  reproduces that with `max(0, total_ros - 1)`, so the two GPX4 trajectories
  match.)
- **The deterministic mean-field under-estimates the stochastic mean LP**
  (0.21 vs 0.46). This is expected and honest: because the autocatalytic
  propagation is a nonlinear bistable switch, the per-cell parameter variation
  and per-step noise push a tail of cells to high LP, inflating the true mean
  above the central-parameter deterministic trajectory (a Jensen-inequality /
  nonlinear-averaging effect). GSH and GPX4, which are closer to linear over this
  window, match closely. The SBML is therefore a faithful export of the
  deterministic core, not a replacement for the stochastic engine.

## Use it

```bash
pip install tellurium python-libsbml     # round-trip toolchain (not a repo runtime dep)
python3 simulations/sbml/export_ferroptosis_sbml.py
```

Then `simulations/sbml/ferroptosis_core.xml` can be opened directly in COPASI or
Tellurium, or deposited in BioModels. The script is standalone (it is not part of
the pytest suite, and `tellurium`/`libsbml` are validation-only dependencies, not
added to `requirements-lock.txt`).
