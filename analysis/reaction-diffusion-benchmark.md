# Reaction-diffusion vs exponential-proxy supply benchmark (#343)

How the explicit-vasculature O2/drug supply model differs when it solves the steady-state reaction-diffusion field instead of the monotonic `exp(-dist_to_nearest_vessel / λ)` proxy, what drives the difference, and where each is adequate.

Companion to the library module `ferroptosis-core::reaction_diffusion` (#343 PR 1) and its sim-tme-3d wiring (`Overrides.reaction_diffusion`, `--snapshot=reaction-diffusion`, #343 PR 2).

## Background and the gap

The TME supply layers (`oxygen::radial_o2_field`, `vasculature::vessel_supply_field`, `slab::slab_supply_field`) model O2/drug availability as a function that is **monotonic in distance to the nearest source**: `exp(-d / λ)` where `d` is the distance to the closest vessel and `λ = sqrt(D/k)` is a fixed penetration length. This is the exact 1-D analytical solution of the steady-state reaction-diffusion equation around a *single isolated planar source*, so it is correct in that limit. It is the wrong shape for a real, irregular 3-D vessel network. Three effects separate the proxy from a true diffusion solve, in rough order of how much they move the numbers here:

1. **Source geometry (the dominant term).** The proxy is the *planar* single-source solution. A real 3-D *point* vessel's field is the Yukawa form `~ exp(-r/λ)/r`, which falls off faster than `exp(-r/λ)`. This is present already at a *single isolated vessel* with zero superposition and zero extra consumption: on an all-tumor cube (λ=120 µm, h=20 µm) the first-neighbour supply is **0.85 (proxy) vs 0.31 (RD)**, a 0.54 gap. Setting the RD diffusion length equal to the proxy's λ equalizes the decay *length* but not the source *geometry*, so this term is a modeling *choice* (planar vs point Green's function), not the proxy being "wrong".
2. **Cumulative consumption — already largely in λ.** A cell is screened by tissue consuming O2 along the path to the vessel. But `λ = sqrt(D/k)` already bakes straight-line-path consumption into the proxy's decay: in the RD solver a single isolated *planar* source reproduces `exp(-x/λ)` to < 5e-4 *despite* per-voxel consumption. So consumption is **not** a separate omitted driver for a single source; it contributes only the *extra* screening along the bent paths of a multi-vessel field (a second-order term here).
3. **Superposition.** Diffusion sums all nearby vessels, so a point between several can be better supplied than its single nearest neighbour implies (the non-monotonic "pocket"). Demonstrated cleanly in the library test `two_sources_enrich_the_midpoint_above_the_nearest_vessel_proxy`. In whole-tumor 3-D geometry at λ ≈ vessel spacing this is small (≈ 0.3% of cells); it dominates only when λ is large relative to the spacing (see "λ regime" below).

The `reaction_diffusion` module solves `D ∇²c − k·c = 0` (vessel Dirichlet sources `c=1`, no-flux boundary) over the **same** vessel network at the **same** λ, by deterministic SOR.

## Verification (analytical self-consistency, PR 1)

The solver is checked against the closed-form 1-D reaction-diffusion slab solution `c(x) = cosh((W−x)/λ)/cosh(W/λ)` (the standard textbook result; the diffusion-limited-distance concept is the Thomlinson & Gray 1955 lineage). On a 1-D bar the converged SOR field matches the closed form to **max abs error < 0.01** (`reaction_diffusion::tests::solver_reproduces_1d_analytical_slab`). This confirms the discretization converges to the exact solution of its own continuous equation; it is numerical verification, not an external cross-check (see "external benchmark" below).

## Proxy vs RD field comparison

Generated on a 36³ generated tumor sphere, λ = 120 µm (the snapshot `ZONE_REF_LAMBDA`, ≈ 6 cells at h = 20 µm), comparing `vessel_supply_field` (proxy) to `reaction_diffusion_supply_field` (RD) over the 17,845 tumor cells, plus the single-isolated-point-source control. Regenerate with:

```
cd simulations
cargo test -p sim-tme-3d --release reaction_diffusion_supply_differs -- --nocapture
```

| case | vessels | proxy → RD (representative) | mean&#124;RD−proxy&#124; | cells RD<proxy | hypoxic frac (supply<0.1) proxy → RD |
|------|--------:|:---------------------------:|--------:|---------------:|:-----------------------------------:|
| single point source (control) | 1 | 0.85 → 0.31 (first neighbour) | — | — | — |
| sparse (`poorly_vascularized`) | 2 | — | 0.104 | 100.0% | **0.623 → 0.988** |
| dense (`well_vascularized`) | 42 | — | 0.233 | 99.7% | 0.000 → 0.000 |

Findings, **in the λ ≈ vessel-spacing regime** (λ = 120 µm here):

- **The RD field is lower than the proxy almost everywhere** (100% sparse, 99.7% dense). The single-point-source control shows the bulk of this is *source geometry* (point vs planar), present before any multi-vessel effect; consumption and superposition are second-order. The mean per-cell gap is 0.10 (sparse) to 0.23 (dense). The dense gap is *larger* not because the network is "worse" but because its mean proxy supply is higher (≈ 0.50 vs 0.15), leaving more headroom to be lower.
- **In a sparse network the gap is qualitative for the hypoxic call.** The proxy calls the tumor 62% hypoxic; RD calls it 99% hypoxic. A study ranking "where does drug reach but O2 does not" off the proxy would badly under-count the hypoxic compartment *in this regime*.

### λ regime: the direction is not calibration-independent

The sign of the proxy-RD gap depends on λ relative to the vessel spacing. Consumption enters the solver as `γ = (h/λ)²`: as λ grows large relative to the cell spacing, `γ → 0`, the equation degenerates to Laplace, and a harmonic field between Dirichlet sources is *enriched* above the nearest-vessel exponential (the `two_sources_enrich...` limit) — so the proxy *under*estimates. As λ shrinks toward the spacing, point-source geometry + consumption dominate and the proxy *over*estimates (the table above). The crossover sits near λ ≈ inter-vessel spacing. The headline "proxy is optimistic / RD more hypoxic" therefore holds **for λ ≲ vessel spacing** (the λ = 120 µm ≈ 6-cell case shown); it is not a calibration-independent statement, and λ here is an uncalibrated placeholder.

## End-to-end (through the full sim)

A 30³ RSL3 run, 120 steps, sparse vasculature, proxy vs RD supply (`reaction_diffusion_supply_differs...`, end-to-end section):

| supply model | reported vascular hypoxic fraction | RSL3 total_dead |
|--------------|:----------------------------------:|:---------------:|
| proxy | 0.580 | 0 |
| RD | 0.992 | 0 |

The supply model carries all the way through to the emitted `vascular_hypoxic_fraction` (0.58 → 0.99). The downstream RSL3 **kill count** is threshold-damped (both = 0 on a sparse, mostly-hypoxic tumor, because RSL3 is hypoxia-sensitive and fails either way), so the hypoxic fraction, not the kill count, is the sensitive end-to-end signal in this regime.

## Where the proxy is adequate vs misleading

- **Adequate:** in the genuine single-isolated-*planar*-source 1-D limit it was derived for (exact); and for a *binary* hypoxic/normoxic call in densely vascularized tissue (both models agree no cell crosses the threshold in the dense case).
- **Misleading (at λ ≲ vessel spacing):** for the *absolute* supply magnitude essentially everywhere (the proxy reads ≈ 0.1–0.2 high on average, dominated by point-vs-planar geometry), which matters for any dose-graded or sub-threshold effect; and for the *hypoxic fraction* in sparse networks (62% vs 99% here).

Practical guidance: keep the cheap proxy for dense-vasculature binary-hypoxia uses and as the byte-identical default; switch on `Overrides.reaction_diffusion` (or `--snapshot=reaction-diffusion`) when the question is the absolute supply level or the hypoxic burden of a sparse network, and read the result at the λ you actually set.

## Caveats (what sets the magnitude)

The numbers above are **directional, not calibrated**:

- **λ regime** (above): the *sign* flips with λ relative to vessel spacing; the over-estimate headline is the λ ≲ spacing regime.
- **Vessel representation / grid resolution.** Each vessel is a single Dirichlet voxel (`vessel_mask`), so it behaves like a ≈ 5–10 µm-radius capillary at h = 20 µm. The near-source falloff — and therefore the 0.85→0.31 single-source gap and the headline fractions — depend on h and on the assumed vessel radius (a finer grid only deepens the depletion; the *direction* does not flip). The specific "62% vs 99%" is representation-dependent; the robust result is the direction.
- **λ, D, k** are uncalibrated placeholders (see `simulations/calibration/CALIBRATION_STATUS.md`). The qualitative conclusion (in the λ ≲ spacing regime: the proxy reads optimistic, worst in sparse networks, dominated by point-source geometry) follows from the physics; the magnitudes do not transfer without calibration.

## External-model benchmark (scope)

Issue #343's acceptance asks to benchmark "against PhysiCell (via the existing FFI) and/or a published tumor reaction-diffusion model." This report delivers the comparison vs the exponential proxy (above) and the validation vs the closed-form 1-D analytical reaction-diffusion solution (PR 1, a published-model result for the single-source geometry). A full **PhysiCell / BioFVM cross-check** of the 3-D multi-vessel field is *not* run here — it would require building PhysiCell and exporting the same vessel geometry through the existing C-FFI, a heavier separately-scoped effort, left as a future extension. Honest status: solver verified against the analytical model, benchmarked against the in-repo proxy; an independent-simulator cross-check is outstanding.

## Reproduce

```
cd simulations
# single-source control + field + end-to-end comparison numbers (this report):
cargo test -p sim-tme-3d --release reaction_diffusion_supply_differs -- --nocapture
# gating / byte-identity proof:
cargo test -p sim-tme-3d reaction_diffusion_flag_is_inert_without_vasculature
# snapshot preset wiring guard:
cargo test -p sim-tme-3d reaction_diffusion_snapshot_preset_is_wired
# analytical-model validation + two-source enrichment (PR 1):
cargo test -p ferroptosis-core reaction_diffusion
# visualization:
cargo run --release -p sim-tme-3d -- --snapshot=reaction-diffusion
python ../scripts/render_tme_3d_trajectory.py output/...   # vessel_supply.npy panel shows the RD field
```
