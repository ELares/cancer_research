//! 3D spatial immune coupling: DAMP diffusion + activation primitives.
//!
//! Ferroptotic cell death releases damage-associated molecular patterns
//! (DAMPs) — calreticulin, HMGB1, ATP — that diffuse through the
//! extracellular space, activate nearby dendritic cells via TLR/RAGE
//! signaling, and prime T-cell-mediated tumor killing. This module
//! provides the 3D-spatial primitives that downstream consumers
//! (#195 sim-tme-3d) compose into a full immune-step orchestrator.
//!
//! **Scope vs the existing [`crate::immune`] module.** That module is the
//! *dimensionless* single-event ICD cascade (one death → one DAMP burst →
//! one DC activation → one T-cell kill). This module is the *spatial*
//! complement: how DAMPs diffuse across a 3D spheroid grid and how local
//! DAMP concentration drives per-cell kill probability. The two compose;
//! `immune` answers "what does one death contribute?" and `immune_3d`
//! answers "where does it spread and who does it affect?"
//!
//! **The 104:1 question (issue #188).** Sim-tme's 2D model finds SDT
//! produces ~104× more immune kills than RSL3 because SDT's dense kill
//! field creates a high local DAMP concentration. The issue asks
//! whether this ratio holds in 3D. **Answering it requires a full
//! multi-step simulation** (sim-tme-3d, #195) — not a library unit
//! test. This module provides the diffusion primitive; the kill-ratio
//! comparison lands with #196 (3D validation).
//!
//! ## ⚠️ Stability requirement (critical for 3D)
//!
//! [`diffuse_damp_3d_step`] is mathematically stable only when
//! `diffusion_fraction × max_neighbor_count < 1`. In 3D with up to 26
//! Moore neighbors, that means **`diffusion_fraction < ≈ 0.038`**.
//!
//! Sim-tme's 2D default `0.08` is **unsafe in 3D** (0.08 × 26 = 2.08 > 1):
//! the source cell loses more DAMP per step than it has, the defensive
//! `.max(0.0)` clamp destroys mass, and the field silently produces
//! nonsense. `debug_assert!` rejects unsafe values in tests; release
//! builds silently mass-destroy. A consumer porting sim-tme's parameters
//! verbatim would hit this immediately — hence the prominent warning.
//!
//! **Suggested 3D-safe value: 0.025**, which gives the same per-step
//! total-diffusion fraction as 2D's `0.08 × 8 = 0.64` (compare to
//! `0.025 × 26 = 0.65`).
//!
//! ## API design — scratch-buffer pattern
//!
//! [`diffuse_damp_3d_step`] takes both `damp_field: &mut [f64]` AND
//! `scratch: &mut [f64]` (both length `grid.cells.len()`). The scratch
//! buffer ensures the spread step is order-independent (otherwise
//! source-ordering would bias the result). **Allocate scratch ONCE
//! before the simulation loop and reuse per step** — at 100³ × 180
//! steps, per-step allocation would be ~1.4 GB of churn.
//!
//! ## Quick example
//!
//! ```
//! use ferroptosis_core::grid::TumorGrid3D;
//! use ferroptosis_core::immune_3d::{diffuse_damp_3d_step, dc_activation, immune_kill_probability};
//!
//! let g = TumorGrid3D::generate(10, 10, 10, 20.0, 42);
//! let n = g.cells.len();
//!
//! // Allocate state once, reuse per step.
//! let mut damp_field = vec![0.0_f64; n];
//! let mut scratch = vec![0.0_f64; n];
//!
//! // Simulate a single death at the center: DAMP burst.
//! let center_idx = g.flat_index(5, 5, 5);
//! damp_field[center_idx] = 10.0;
//!
//! // 3D-safe diffusion_fraction (see stability requirement above).
//! diffuse_damp_3d_step(&mut damp_field, &mut scratch, &g, 0.025, 0.03);
//!
//! // Per-cell immune activation and kill probability.
//! for &local_damp in &damp_field {
//!     let activation = dc_activation(local_damp, 50.0);
//!     let kill_prob = immune_kill_probability(activation, 0.02, 0.21);
//!     // consumer rolls: if rng.gen() < kill_prob { ... }
//!     let _ = (activation, kill_prob);
//! }
//! ```

use crate::grid::TumorGrid3D;

/// Maximum Moore-neighbor count in 3D (3×3×3 cube − self).
const MAX_3D_NEIGHBORS: usize = 26;

/// One step of DAMP diffusion + exponential clearance on a 3D spheroid
/// grid. Mutates `damp_field` in place using `scratch` to avoid
/// source-ordering bias (each step, every source spreads `share = local
/// × diffusion_fraction` to each of its up-to-26 Moore neighbors, then
/// every cell decays by `(1 − clearance_rate)`).
///
/// **Stability requirement** (see module doc): `diffusion_fraction ×
/// MAX_3D_NEIGHBORS < 1.0`. `debug_assert!` catches the unsafe case in
/// tests; release silently mass-destroys via the defensive `.max(0.0)`
/// clamp. **Suggested 3D-safe value: `0.025`** (matches 2D's
/// per-step total diffusion of ~64%).
///
/// **Length contract**: `damp_field.len() == scratch.len() == grid.cells.len()`.
/// Validated with `assert!` (matches stromal pattern — programming-
/// contract bug deserves a clear release-mode panic).
///
/// **Cost**: O(N × 26) for N = `grid.cells.len()`. Same per-call
/// recompute concern as the rest of the 3D code; #194 hoisting applies.
pub fn diffuse_damp_3d_step(
    damp_field: &mut [f64],
    scratch: &mut [f64],
    grid: &TumorGrid3D,
    diffusion_fraction: f64,
    clearance_rate: f64,
) {
    let n = grid.cells.len();
    assert!(
        damp_field.len() == n,
        "diffuse_damp_3d_step: damp_field.len() {} must equal grid.cells.len() {}",
        damp_field.len(),
        n
    );
    assert!(
        scratch.len() == n,
        "diffuse_damp_3d_step: scratch.len() {} must equal grid.cells.len() {}",
        scratch.len(),
        n
    );
    debug_assert!(
        diffusion_fraction.is_finite() && diffusion_fraction >= 0.0,
        "diffuse_damp_3d_step: diffusion_fraction must be finite and ≥ 0, got {diffusion_fraction}"
    );
    // **Stability check is `assert!` (not `debug_assert!`)** because the
    // failure mode is silent: violation causes the source's `local + delta`
    // to go negative, then the defensive `.max(0.0)` clamp destroys mass
    // with no panic. A consumer porting sim-tme's 2D default (0.08) would
    // get wrong DAMP fields in release with no indication. The per-call
    // cost (one multiply + one compare) is negligible.
    assert!(
        diffusion_fraction * (MAX_3D_NEIGHBORS as f64) < 1.0,
        "diffuse_damp_3d_step: diffusion_fraction × 26 = {} must be < 1 for stability. \
         Got diffusion_fraction = {}. Sim-tme's 2D default (0.08) is UNSAFE in 3D \
         (0.08 × 26 = 2.08 > 1) — use ≤ 0.038, suggested 0.025 to match 2D's per-step total.",
        diffusion_fraction * (MAX_3D_NEIGHBORS as f64),
        diffusion_fraction
    );
    debug_assert!(
        clearance_rate.is_finite() && (0.0..=1.0).contains(&clearance_rate),
        "diffuse_damp_3d_step: clearance_rate must be in [0, 1], got {clearance_rate}"
    );

    // Compute the spread into scratch.
    scratch.fill(0.0);
    for r in 0..grid.rows {
        for c in 0..grid.cols {
            for l in 0..grid.layers {
                let idx = grid.flat_index(r, c, l);
                let local = damp_field[idx];
                if local < 0.001 {
                    continue;
                }
                let share = local * diffusion_fraction;
                let (neighbors, count) = grid.neighbors(r, c, l);
                for &(nr, nc, nl) in &neighbors[..count] {
                    let nidx = grid.flat_index(nr, nc, nl);
                    scratch[nidx] += share;
                }
                // Source loses exactly `share × count` (where count is the
                // ACTUAL neighbor count, not always 26 at boundaries).
                // This is what makes mass conservation hold across edges.
                scratch[idx] -= share * count as f64;
            }
        }
    }

    // Apply spread + clearance.
    for i in 0..n {
        damp_field[i] = (damp_field[i] + scratch[i]).max(0.0);
        damp_field[i] *= 1.0 - clearance_rate;
    }
}

/// Dendritic-cell activation as a function of local DAMP concentration.
///
/// Michaelis-Menten saturation: `activation = damp / (damp + kd)`, in
/// `[0, 1]` for non-negative inputs. At `damp = 0`: returns `0.0`
/// exactly (IEEE: `0/kd = 0`). At `damp = kd`: returns `0.5` exactly
/// (IEEE: `x/(2x) = 0.5` for finite x). For `damp >> kd`: asymptotes
/// to `1.0`.
///
/// **Pure scalar function.** No clamp; trust caller. Sim-tme uses
/// `kd = 50.0` (the `dc_activation_kd` field of `ImmuneConfig`).
///
/// **Release behavior for invalid inputs**:
///
/// | Bad input | Output |
/// |-----------|--------|
/// | `damp = NaN` | `NaN` (propagates) |
/// | `kd = NaN` | `NaN` |
/// | `damp + kd = 0` (both 0) | `NaN` (`0/0`) |
/// | `damp < 0`, `damp + kd < 0` | negative; consumer's problem |
/// | `kd < 0` | math still works but biologically meaningless |
///
/// Consumers passing finite `damp >= 0` and `kd > 0` (the normal regime)
/// don't hit any of these.
#[inline]
#[must_use = "the activation is the function's only output; ignoring it suggests a logic bug"]
pub fn dc_activation(local_damp: f64, kd: f64) -> f64 {
    local_damp / (local_damp + kd)
}

/// Per-cell immune kill probability per step.
///
/// `probability = (activation × kill_rate × (1 − effective_brake)).min(0.99)`
///
/// The `.min(0.99)` cap matches sim-tme: even at full activation with no
/// PD-1 brake, kills are never guaranteed (preserves stochasticity over
/// long simulations). At `activation = 0` or `kill_rate = 0` or
/// `effective_brake = 1`: returns `0.0` exactly.
///
/// **Pure scalar function.** Sim-tme uses `kill_rate = 0.02` and computes
/// `effective_brake = pd1_brake × (1 - anti_pd1_efficacy)` upstream.
///
/// **Lower bound**: NOT clamped at zero. For pathological inputs (e.g.,
/// `activation < 0`) the return value can be negative. Caller's
/// `rng.gen::<f64>() < probability` correctly produces zero kills in
/// that case (random in `[0, 1)` is never less than a negative), so
/// downstream semantics are safe — but the return-type contract is
/// `(-∞, 0.99]`, not `[0, 0.99]`.
///
/// **Release behavior for invalid inputs**:
///
/// | Bad input | Output |
/// |-----------|--------|
/// | any `NaN` argument | `NaN` (`f64::min` returns NaN if self is NaN) |
/// | `activation < 0` or `effective_brake > 1` | negative; clamped only at upper bound |
#[inline]
#[must_use = "the kill probability is the function's only output; ignoring it suggests a logic bug"]
pub fn immune_kill_probability(activation: f64, kill_rate: f64, effective_brake: f64) -> f64 {
    (activation * kill_rate * (1.0 - effective_brake)).min(0.99)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::TumorGrid;

    // =============================
    // diffuse_damp_3d_step tests
    // =============================

    /// **v2 addition**: stability `debug_assert` rejects sim-tme's 2D
    /// default (0.08), which is unsafe in 3D (0.08 × 26 = 2.08 > 1).
    /// Critical bug-class guard.
    #[test]
    #[should_panic(expected = "diffusion_fraction × 26")]
    fn diffusion_fraction_stability_assertion() {
        let g = TumorGrid3D::generate(5, 5, 5, 20.0, 42);
        let mut field = vec![0.0; g.cells.len()];
        let mut scratch = vec![0.0; g.cells.len()];
        diffuse_damp_3d_step(&mut field, &mut scratch, &g, 0.08, 0.03);
    }

    /// Field length mismatch panics in release (matches stromal pattern).
    #[test]
    #[should_panic(expected = "damp_field.len()")]
    fn field_length_validation() {
        let g = TumorGrid3D::generate(5, 5, 5, 20.0, 42);
        let mut field = vec![0.0; g.cells.len() / 2];
        let mut scratch = vec![0.0; g.cells.len()];
        diffuse_damp_3d_step(&mut field, &mut scratch, &g, 0.025, 0.03);
    }

    /// Scratch length mismatch panics in release.
    #[test]
    #[should_panic(expected = "scratch.len()")]
    fn scratch_length_validation() {
        let g = TumorGrid3D::generate(5, 5, 5, 20.0, 42);
        let mut field = vec![0.0; g.cells.len()];
        let mut scratch = vec![0.0; g.cells.len() / 2];
        diffuse_damp_3d_step(&mut field, &mut scratch, &g, 0.025, 0.03);
    }

    /// A single DAMP source at an INTERIOR cell spreads to all 26
    /// neighbors after one step. With no clearance, each neighbor
    /// receives exactly `share = local × diffusion_fraction`, and the
    /// source retains `local × (1 - 26 × diffusion_fraction)`.
    #[test]
    fn interior_source_spreads_to_26_neighbors() {
        let g = TumorGrid3D::generate(10, 10, 10, 20.0, 42);
        let n = g.cells.len();
        let mut field = vec![0.0; n];
        let mut scratch = vec![0.0; n];

        let center = g.flat_index(5, 5, 5); // interior cell
        let source_value = 100.0_f64;
        let fraction = 0.025_f64;
        field[center] = source_value;

        diffuse_damp_3d_step(&mut field, &mut scratch, &g, fraction, 0.0);

        // Source retains local × (1 - 26 × fraction) = 100 × (1 - 0.65) = 35
        let expected_source = source_value * (1.0 - 26.0 * fraction);
        assert!(
            (field[center] - expected_source).abs() < 1e-9,
            "source should retain {expected_source}, got {}",
            field[center]
        );

        // Each of the 26 neighbors should have received exactly share = 2.5
        let expected_share = source_value * fraction;
        let (neighbors, count) = g.neighbors(5, 5, 5);
        assert_eq!(count, 26);
        for &(nr, nc, nl) in &neighbors[..26] {
            let nidx = g.flat_index(nr, nc, nl);
            assert!(
                (field[nidx] - expected_share).abs() < 1e-9,
                "neighbor ({nr},{nc},{nl}) should have {expected_share}, got {}",
                field[nidx]
            );
        }
    }

    /// A source at the grid CORNER has only 7 neighbors. Source loses
    /// `share × 7` (not `share × 26`); each of the 7 neighbors gains
    /// `share`. Mass is still conserved.
    #[test]
    fn corner_source_spreads_to_7_neighbors() {
        let g = TumorGrid3D::generate(10, 10, 10, 20.0, 42);
        let n = g.cells.len();
        let mut field = vec![0.0; n];
        let mut scratch = vec![0.0; n];

        let corner = g.flat_index(0, 0, 0);
        let source_value = 100.0_f64;
        let fraction = 0.025_f64;
        field[corner] = source_value;

        diffuse_damp_3d_step(&mut field, &mut scratch, &g, fraction, 0.0);

        // Source retains 100 × (1 - 7 × 0.025) = 100 × 0.825 = 82.5
        let expected_source = source_value * (1.0 - 7.0 * fraction);
        assert!(
            (field[corner] - expected_source).abs() < 1e-9,
            "corner source should retain {expected_source}, got {}",
            field[corner]
        );

        let (_neighbors, count) = g.neighbors(0, 0, 0);
        assert_eq!(
            count, 7,
            "grid corner should have exactly 7 Moore neighbors"
        );
    }

    /// Total DAMP is conserved (modulo clearance) across a diffusion step.
    /// `sum_after = sum_before × (1 - clearance_rate)` to numerical
    /// tolerance.
    #[test]
    fn total_damp_conserved_modulo_clearance() {
        let g = TumorGrid3D::generate(10, 10, 10, 20.0, 42);
        let n = g.cells.len();
        let mut field = vec![0.0; n];
        let mut scratch = vec![0.0; n];

        // Sparse, varied initial conditions.
        field[g.flat_index(5, 5, 5)] = 100.0;
        field[g.flat_index(3, 7, 2)] = 50.0;
        field[g.flat_index(8, 1, 6)] = 25.0;

        let sum_before: f64 = field.iter().sum();
        let clearance = 0.03_f64;
        diffuse_damp_3d_step(&mut field, &mut scratch, &g, 0.025, clearance);
        let sum_after: f64 = field.iter().sum();

        let expected = sum_before * (1.0 - clearance);
        let rel_error = (sum_after - expected).abs() / expected;
        assert!(
            rel_error < 1e-12,
            "mass not conserved: before={sum_before}, after={sum_after}, expected≈{expected} (rel error {rel_error:.2e})"
        );
    }

    /// Zero diffusion → only clearance acts → every cell decays by
    /// `(1 - clearance)` exactly.
    #[test]
    fn zero_diffusion_just_applies_clearance() {
        let g = TumorGrid3D::generate(10, 10, 10, 20.0, 42);
        let n = g.cells.len();
        let mut field = vec![1.0_f64; n]; // uniform
        let mut scratch = vec![0.0; n];

        diffuse_damp_3d_step(&mut field, &mut scratch, &g, 0.0, 0.03);

        let expected = 0.97_f64;
        for (i, &v) in field.iter().enumerate() {
            assert!(
                (v - expected).abs() < 1e-12,
                "cell {i} should be {expected}, got {v}"
            );
        }
    }

    /// Determinism: same inputs → same outputs.
    #[test]
    fn deterministic_same_inputs_same_outputs() {
        let g = TumorGrid3D::generate(10, 10, 10, 20.0, 42);
        let n = g.cells.len();
        let init: Vec<f64> = (0..n).map(|i| (i as f64 % 7.0) * 0.5).collect();

        let mut field1 = init.clone();
        let mut scratch1 = vec![0.0; n];
        diffuse_damp_3d_step(&mut field1, &mut scratch1, &g, 0.025, 0.03);

        let mut field2 = init.clone();
        let mut scratch2 = vec![0.0; n];
        diffuse_damp_3d_step(&mut field2, &mut scratch2, &g, 0.025, 0.03);

        assert_eq!(field1, field2);
    }

    /// **Cross-geometry** (AC #4, scope-bounded): at matched
    /// `diffusion_fraction` and an interior source, a 3D source decays
    /// **faster** than a 2D source per step because 26 > 8 neighbors take
    /// share. This is the geometric component of the issue's kill-ratio
    /// question; the full kill-ratio answer (whether SDT's 104:1
    /// differential holds in 3D) requires multi-step simulation and lands
    /// with #195/#196.
    ///
    /// Uses 3D-safe `0.025` so the test exercises real diffusion math
    /// rather than mass-destroying via the stability violation.
    #[test]
    fn three_d_source_decays_faster_than_two_d_at_matched_diffusion() {
        let fraction = 0.025_f64;
        let source = 100.0_f64;

        // 2D: 10×10 interior cell, 8 neighbors. Inline sim-tme math
        // (source of truth: sim-tme/src/main.rs:701-716).
        let g2 = TumorGrid::generate(10, 10, 20.0, 42);
        let n2 = g2.cells.len();
        let mut field2 = vec![0.0; n2];
        let center_2d = 5 * g2.cols + 5;
        field2[center_2d] = source;
        let (_, count_2d) = g2.neighbors(5, 5);
        let source_after_2d = source * (1.0 - count_2d as f64 * fraction);

        // 3D: 10×10×10 interior cell, 26 neighbors.
        let g3 = TumorGrid3D::generate(10, 10, 10, 20.0, 42);
        let n3 = g3.cells.len();
        let mut field3 = vec![0.0; n3];
        let mut scratch3 = vec![0.0; n3];
        let center_3d = g3.flat_index(5, 5, 5);
        field3[center_3d] = source;
        diffuse_damp_3d_step(&mut field3, &mut scratch3, &g3, fraction, 0.0);

        assert!(
            field3[center_3d] < source_after_2d,
            "3D source ({}) should retain LESS than 2D source ({}) at matched diffusion — \
             geometry: 26 neighbors take share vs 8. Empirical: 3D loses 26·share = {}, \
             2D loses 8·share = {}.",
            field3[center_3d],
            source_after_2d,
            26.0 * source * fraction,
            8.0 * source * fraction
        );

        // Quantitative: 3D should lose ~26/8 = 3.25× more per step.
        let loss_3d = source - field3[center_3d];
        let loss_2d = source - source_after_2d;
        let ratio = loss_3d / loss_2d;
        let expected_ratio = 26.0 / 8.0; // = 3.25
        assert!(
            (ratio - expected_ratio).abs() < 0.01,
            "loss ratio 3D/2D should be ≈ {expected_ratio} (26/8), got {ratio}"
        );
    }

    // =============================
    // dc_activation tests
    // =============================

    /// `dc_activation(0, kd) = 0` exactly (IEEE: 0/kd = 0 for kd > 0).
    #[test]
    fn dc_activation_at_zero_damp_returns_zero() {
        for &kd in &[1.0, 50.0, 1000.0] {
            assert_eq!(dc_activation(0.0, kd), 0.0);
        }
    }

    /// `dc_activation(kd, kd) = 0.5` exactly (IEEE: x/(2x) = 0.5 for
    /// finite x; even though `kd` may not be IEEE-exact, the doubling
    /// `kd + kd = 2*kd` is exact, then x/(2x) when 2x is finite is
    /// IEEE-exactly 0.5).
    #[test]
    fn dc_activation_at_damp_eq_kd_returns_half() {
        for &kd in &[1.0_f64, 50.0, 1000.0, 7.4] {
            assert_eq!(dc_activation(kd, kd), 0.5);
        }
    }

    /// `dc_activation(damp >> kd, kd) ≈ 1.0` asymptote.
    #[test]
    fn dc_activation_large_damp_approaches_one() {
        let kd = 50.0;
        assert!(dc_activation(1e9, kd) > 0.999999);
        assert!(dc_activation(1e9, kd) <= 1.0);
    }

    /// Monotone non-decreasing: more DAMP → higher activation.
    #[test]
    fn dc_activation_monotone_non_decreasing() {
        let kd = 50.0;
        let probes = [0.0, 1.0, 10.0, 50.0, 100.0, 1000.0];
        let mut prev = dc_activation(probes[0], kd);
        for &d in &probes[1..] {
            let cur = dc_activation(d, kd);
            assert!(
                cur >= prev,
                "activation decreased from {prev} to {cur} at damp={d}"
            );
            prev = cur;
        }
    }

    // =============================
    // immune_kill_probability tests
    // =============================

    /// `min(0.99)` cap: even at activation=1, rate=1, brake=0 → result
    /// is exactly 0.99. Hardcoded cap matches sim-tme.
    #[test]
    fn immune_kill_clamps_at_0_99_for_extreme_inputs() {
        let result = immune_kill_probability(1.0, 1.0, 0.0);
        assert_eq!(result, 0.99);

        // Even more extreme: activation > 1 (shouldn't happen but test the
        // clamp behavior).
        let result2 = immune_kill_probability(2.0, 1.0, 0.0);
        assert_eq!(result2, 0.99);
    }

    /// `effective_brake = 1.0` → full PD-1 suppression → zero kill.
    /// IEEE-exact: `(activation × rate × 0.0) = 0.0` for any finite
    /// activation, rate.
    #[test]
    fn immune_kill_full_brake_returns_zero() {
        for &activation in &[0.0_f64, 0.5, 1.0] {
            for &rate in &[0.01_f64, 0.02, 0.05] {
                assert_eq!(immune_kill_probability(activation, rate, 1.0), 0.0);
            }
        }
    }

    /// Default sim-tme numerical example: activation=0.5, rate=0.02,
    /// brake=0.21 (= pd1_brake 0.7 × (1 - anti_pd1 0.7)) → ≈ 0.0079.
    /// Tight tolerance via libm.
    #[test]
    fn immune_kill_at_sim_tme_defaults_matches_expected() {
        let activation = 0.5;
        let rate = 0.02;
        let effective_brake = 0.7 * (1.0 - 0.7); // = 0.21
        let prob = immune_kill_probability(activation, rate, effective_brake);
        // Expected: 0.5 × 0.02 × (1 - 0.21) = 0.5 × 0.02 × 0.79 = 0.0079
        let expected = 0.5 * 0.02 * 0.79;
        assert!(
            (prob - expected).abs() < 1e-9,
            "expected ≈ {expected}, got {prob}"
        );
        assert!(prob < 0.99, "should not be capped at this activation");
    }
}
