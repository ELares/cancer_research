# ODE cross-validation: ferroptosis-core vs published ferroptosis models

`ferroptosis-core`'s biochemistry ODEs were built from first principles plus
estimated constants. This note cross-validates the engine's *qualitative
dynamics* against independent, published ferroptosis / redox dynamical-systems
models: same kind of system, do we get the same kind of behavior? It is a cheap,
high-value structural check that catches the worst failure mode (a model that
produces the wrong *shape* of dynamics), distinct from parameter calibration
(which is separate, uncalibrated work tracked in `CALIBRATION_STATUS.md` and
issues #330 to #335).

Figure: `simulations/calibration/ode_cross_validation.png`
(generator: `simulations/calibration/cross_validate_odes.py`, which runs the
ACTUAL engine via the `ferroptosis_core` Python bindings).

## The published comparators (all PMIDs verified via NCBI esummary)

| Model | PMID | Type | Core behavior |
|---|---|---|---|
| Co et al., *Nature* 2024 | 38987590 | Continuous ODE, Fenton + NADPH-oxidase positive feedback minus glutathione clearance | Single-cell ROS steady state **bifurcates monostable -> bistable** as antioxidant defense is suppressed; an unstable-steady-state threshold separates recover from runaway. The canonical bistability precedent. |
| Seidel et al., *Front Cell Dev Biol* 2026 | 41960191 | Minimal 2-ODE (population + lipid-ROS), logistic growth + Fenton-amplified ROS + threshold death | **Two stable steady states** (ferroptosis-insensitive vs sensitive) separated by a tipping point. The most directly re-implementable. (Peer-reviewed version of bioRxiv 2025.09.15.676259, which has no PMID.) |
| Konstorum et al., *J Theor Biol* 2020 | 32114023 | Discrete, stochastic, multistate **logical** model (11 species) | GPX4 is the critical brake on lipid-peroxide (LOOH) accumulation; "high ACSL4 is necessary but not sufficient." Structural, not a trajectory comparator (discrete logic). |
| Pannala et al., *Free Radic Res* 2014 | 24456207 | Enzyme-kinetic ODE of glutathione peroxidase (ping-pong, pH, GSSG) | Mechanistic rate law for the GPX/GSH step; a reference for the *functional form* of our GPX4/GSH repair term, not the system switch. |

The single behavior all three dynamical models report, and that a structurally
sound ferroptosis model MUST reproduce: **a bistable switch with a GSH/GPX4-set
threshold (a separatrix in ROS / lipid-peroxide space) below which lipid peroxide
is repaired (recover) and above which Fenton-driven positive feedback runs it
away to a high death state (collapse), with the system tipping monostable ->
bistable as antioxidant defense is suppressed.**

## What ferroptosis-core produces (actual engine)

Running the real engine under an exogenous-ROS (SDT) dose sweep on OXPHOS cells:

- **Panel A (the key result): the single-cell final lipid-peroxide distribution
  at the tipping dose is BIMODAL.** At `sdt_ros = 2.0` (near the population 50%
  point), 29% of cells finish at low lipid peroxide (recovered, LP < 2) and 68%
  at high (collapsed, LP > 8), with only ~3% in the middle. A near-empty middle
  is the separatrix signature of a bistable switch: a cell either stays in the
  low basin or runs away to the high one, almost never settles in between. This
  is the direct single-cell analog of Co 2024's unstable-threshold separatrix
  and Seidel 2026's two stable states.
- **Panel B: the population death rate vs exogenous-ROS dose is a sharp sigmoid
  threshold** (18% at the low dose to 100% at the high dose, steepest near
  `sdt_ros` 1.5 to 2.0). This is the population-level manifestation of the
  underlying single-cell bifurcation: as the drive crosses the bistable region,
  the fraction of cells tipped into the high basin rises steeply.
- **Panel C: the shared canonical structure.** A minimal one-variable bistable
  ODE (sigmoidal Fenton-positive-feedback ROS production minus linear antioxidant
  clearance, the structure shared by Co 2024 and Seidel 2026, NOT a verbatim
  reproduction of either paper's parameters) shows the monostable -> bistable
  fold: a high-clearance (well-defended) cell has one low stable state, while a
  suppressed-clearance cell has three intersections (low stable, unstable
  threshold, high stable). ferroptosis-core's bimodal recover/collapse is exactly
  the dynamics this structure predicts.

The mechanism in the code that produces this is the autocatalytic propagation
term in `biochem::sim_cell_step`:
`propagation_rate = lp_propagation / (1 + antioxidant_quench * 5)`, where
`antioxidant_quench = gpx4 * (gsh/(gsh+0.5)) + fsp1 (+ optional gch1)`. When
quench is high the propagation is damped (recover); when GPX4 is inhibited and
GSH depleted, the propagation runs away (collapse). That GSH/GPX4-gated
positive-feedback loop IS the bistable switch the published models describe.

## Agreements, divergences, and parameter implications

**Agreements (structural cross-validation PASSES):**
- Bistability / recover-or-collapse with a threshold separatrix: present
  (Panel A), matching Co 2024 and Seidel 2026.
- Fenton positive feedback as the engine of the high-ROS branch: present (the
  `fenton` term + autocatalytic propagation), matching both continuous models.
- GPX4 as the critical brake: present (GPX4 inhibition via RSL3 is what tips
  cells over the threshold), matching Konstorum 2020's central finding.
- The threshold behavior under antioxidant suppression (Panel B sigmoid):
  matches Co 2024's monostable -> bistable transition direction.

**Divergences (expected, and what they imply):**
- ferroptosis-core is a STOCHASTIC Monte-Carlo single-cell model; the comparators
  are deterministic. So our "bistability" appears as a *bimodal distribution over
  cells* rather than two fixed points of one trajectory. This is the correct
  stochastic analog, not a discrepancy.
- The single-cell engine has no explicit population or spatial axis, so it cannot
  reproduce Seidel 2026's density-dependent basin selection or Co 2024's spatial
  trigger waves directly (the 3D suite's spatial layers are where that would
  live, and remain uncalibrated).
- Konstorum 2020 is discrete-logic, so only a structural (which-species-gates-
  what) comparison is meaningful, not a trajectory overlay.
- The GPX/GSH repair term here is a lumped Michaelis-Menten form, simpler than
  Pannala 2014's full ping-pong enzyme kinetics; matching the *functional shape*
  (saturating in GSH) is the relevant check, and it does.

**No structural discrepancy was found.** The pre-stated red flag (a single
monotone decay, or a non-hysteretic graded unimodal response, instead of two
basins separated by a threshold) does NOT occur: the LP distribution is bimodal
with a near-empty middle. Had it been unimodal-graded, that would have indicated
the autocatalytic-propagation gating is too weak to produce a real switch.

**Parameter implication:** the location of the tipping dose (and the Panel B
sigmoid midpoint) is set by the GSH/GPX4 antioxidant setpoint and the
`lp_propagation` / `lp_rate` constants. This cross-validation confirms the
STRUCTURE is right; it does NOT calibrate those constants (the position of the
threshold on a real dose axis still requires the calibration data in #330 to
#335). Cross-model agreement on *shape* plus pending calibration of *position* is
the honest status.

## Reproduce

```
python3 simulations/calibration/cross_validate_odes.py
```

References for the comparator models are recorded in
`simulations/calibration/parameter_provenance.md`.
