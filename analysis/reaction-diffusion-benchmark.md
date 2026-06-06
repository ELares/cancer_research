# Reaction-diffusion vs exponential-proxy supply benchmark (#343)

What the explicit-vasculature O2/drug supply model gets wrong when it uses the monotonic `exp(-dist_to_nearest_vessel / λ)` proxy instead of solving the steady-state reaction-diffusion field, and where the proxy is nonetheless adequate.

Companion to the library module `ferroptosis-core::reaction_diffusion` (#343 PR 1) and its sim-tme-3d wiring (`Overrides.reaction_diffusion`, `--snapshot=reaction-diffusion`, #343 PR 2).

## Background and the gap

The TME supply layers (`oxygen::radial_o2_field`, `vasculature::vessel_supply_field`, `slab::slab_supply_field`) all model O2/drug availability as a function that is **monotonic in distance to the nearest source**: a cell's supply is `exp(-d / λ)` where `d` is the distance to the closest vessel (or the edge, or the +z slab face) and `λ` is a fixed penetration length. This is the exact 1-D analytical solution of the steady-state reaction-diffusion equation around a *single isolated planar source*, so it is correct in that limit. It is the wrong shape for a real, irregular vessel network, because:

1. **It ignores cumulative consumption.** The proxy length `λ = sqrt(D/k)` bakes consumption into a single source's decay, but a cell deep in tissue is screened by all the tissue between it and the vessel consuming O2 along the way. The proxy counts only the straight-line distance to the nearest vessel.
2. **It ignores 3-D geometric spreading.** A 3-D point source's true field falls off faster than `exp(-r/λ)` (roughly `exp(-r/λ)/r`), so the proxy overestimates supply away from a point vessel.
3. **It ignores superposition.** Diffusion sums the contributions of all nearby vessels, so a point between several vessels can be better supplied than its single nearest neighbour implies (the non-monotonic "pocket").

The `reaction_diffusion` module instead solves

```
D ∇²c − k·c = 0   in tumor tissue (vessel Dirichlet sources, c=1; no-flux boundary)
```

over the **same** vessel network at the **same** λ, by deterministic SOR. Because the single-source limit reproduces `exp(-d/λ)`, every difference below is a genuine multi-vessel / consumption / geometry effect, not a re-parameterization.

## Validation (analytical self-consistency, PR 1)

Before comparing to the proxy, the solver is checked against the closed-form 1-D reaction-diffusion slab solution `c(x) = cosh((W−x)/λ)/cosh(W/λ)` (the standard textbook result; the diffusion-limited-distance concept is the Thomlinson & Gray 1955 lineage). On a 1-D bar the converged SOR field matches the closed form to **max abs error < 0.01** (`reaction_diffusion::tests::solver_reproduces_1d_analytical_slab`). This confirms the discretization converges to the exact solution of its own continuous equation; it is numerical verification, not an external cross-check (see "external benchmark" below).

## Proxy vs RD field comparison

Generated on a 36³ generated tumor sphere, λ = 120 µm (the snapshot `ZONE_REF_LAMBDA`), comparing `vessel_supply_field` (proxy) to `reaction_diffusion_supply_field` (RD) over the 17,845 tumor cells, for both a sparse and a dense vessel network. Regenerate with:

```
cd simulations
cargo test -p sim-tme-3d --release reaction_diffusion_supply_differs -- --nocapture
```

| network | vessels | mean&#124;RD−proxy&#124; | max&#124;RD−proxy&#124; | cells RD<proxy | hypoxic frac (supply<0.1) proxy → RD |
|---------|--------:|--------:|--------:|---------------:|:---------------------------:|
| sparse (`poorly_vascularized`) | 2 | 0.104 | 0.677 | 100.0% | **0.623 → 0.988** |
| dense (`well_vascularized`) | 42 | 0.233 | 0.545 | 99.7% | 0.000 → 0.000 |

Two robust findings:

- **The proxy systematically overestimates supply.** The RD field is lower than the proxy at essentially every cell (100% sparse, 99.7% dense), because the proxy omits cumulative consumption and 3-D geometric spreading. The mean per-cell overestimate is 0.10 (sparse) to 0.23 (dense) on a [0,1] scale.
- **In a sparse network the error is qualitative.** The proxy calls the tumor 62% hypoxic; the RD field calls it 99% hypoxic. A study ranking "where does drug reach but O2 does not" off the proxy would badly under-count the hypoxic compartment.

The non-monotonic *enrichment* pocket (a point between several vessels supplied better than its nearest-vessel value) is real but small in these realistic 3-D spheres: only 0.3% of cells have RD > proxy in the dense network. It is demonstrated cleanly at the controlled two-source level in `reaction_diffusion::tests::two_sources_enrich_the_midpoint_above_the_nearest_vessel_proxy` (PR 1). In whole-tumor 3-D geometry the consumption-driven *depletion* dominates the superposition-driven *enrichment*.

## End-to-end (through the full sim)

A 30³ RSL3 run, 120 steps, sparse vasculature, proxy vs RD supply (`reaction_diffusion_supply_differs...`, end-to-end section):

| supply model | reported vascular hypoxic fraction | RSL3 total_dead |
|--------------|:----------------------------------:|:---------------:|
| proxy | 0.580 | 0 |
| RD | 0.992 | 0 |

The supply model carries all the way through to the emitted `vascular_hypoxic_fraction` (0.58 → 0.99). The downstream RSL3 **kill count** is threshold-damped (both ≈ 0 on a sparse, mostly-hypoxic tumor, because RSL3 is hypoxia-sensitive and fails either way), so the hypoxic fraction, not the kill count, is the sensitive end-to-end signal. In a well-vascularized 24³ RSL3 run the kill count does differ slightly (proxy 2 vs RD 0), in the expected direction (RD's lower supply → marginally fewer kills), but the absolute numbers are small because the RSL3 kill switch is bistable.

## Where the proxy is adequate vs misleading

- **Adequate:** near vessels and in densely vascularized tissue *for a binary hypoxic/normoxic call* (both models agree no cell crosses the hypoxia threshold in the dense network). The proxy is also exact in the genuine single-isolated-source 1-D limit it was derived for.
- **Misleading:** (1) for the *absolute* supply magnitude everywhere (the proxy overestimates by 0.1–0.2 on average), which matters for any dose-graded or sub-threshold effect; (2) for the *hypoxic fraction* in sparse / poorly-vascularized tumors, where the proxy can be off by a factor that flips the qualitative picture (62% vs 99% hypoxic here); (3) anywhere the *ranking* of where-drug-reaches-but-O2-does-not depends on the supply being a true diffusion field rather than a nearest-vessel distance.

Practical guidance: keep the cheap proxy for dense-vasculature, binary-hypoxia uses and as the byte-identical default; switch on `Overrides.reaction_diffusion` (or `--snapshot=reaction-diffusion`) when the question is the absolute supply level or the hypoxic burden of a sparse network.

## External-model benchmark (scope)

Issue #343's acceptance asks to benchmark "against PhysiCell (via the existing FFI) and/or a published tumor reaction-diffusion model." This report delivers:

- the comparison vs the exponential proxy (above), and
- the validation vs the closed-form 1-D analytical reaction-diffusion solution (PR 1), which is a published-model result for the single-source geometry.

A full **PhysiCell / BioFVM cross-check** of the 3-D multi-vessel field is *not* run here. It would require building PhysiCell and exporting the same vessel geometry through the existing C-FFI, which is a heavier, separately-scoped effort; it is left as a future extension. The honest status is therefore: solver verified against the analytical model, benchmarked against the in-repo proxy; an independent-simulator cross-check is outstanding. The qualitative conclusions above (the proxy is optimistic; the error is largest and qualitative in sparse networks) follow from the consumption + geometry physics and do not depend on the absolute calibration of `λ`, `D`, or `k`, all of which remain uncalibrated placeholders (see `simulations/calibration/CALIBRATION_STATUS.md`).

## Reproduce

```
cd simulations
# field + end-to-end comparison numbers (this report):
cargo test -p sim-tme-3d --release reaction_diffusion_supply_differs -- --nocapture
# gating / byte-identity proof:
cargo test -p sim-tme-3d reaction_diffusion_flag_is_inert_without_vasculature
# analytical-model validation (PR 1):
cargo test -p ferroptosis-core reaction_diffusion
# visualization:
cargo run --release -p sim-tme-3d -- --snapshot=reaction-diffusion
python ../scripts/render_tme_3d_trajectory.py output/...   # vessel_supply.npy panel shows the RD field
```
