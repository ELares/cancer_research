//! 3D radial pH gradient for spheroid tumors.
//!
//! Glycolytic tumors produce lactic acid (Warburg effect), lowering
//! extracellular pH from ~7.4 at the well-perfused periphery to ~6.5 in
//! the core. Two competing effects on ferroptosis:
//!
//! 1. **Iron release** (pro-ferroptosis): low pH destabilizes ferritin →
//!    more Fenton-available Fe²⁺ → more ROS.
//! 2. **Ion trapping** (anti-ferroptosis for weak-base drugs like RSL3):
//!    low extracellular pH → drug protonated → less intracellular delivery
//!    (Henderson-Hasselbalch, linearized over the narrow tumor pH range).
//!
//! This module provides the field-computation primitive and the two
//! per-cell modulation factors that downstream 3D consumers (#188 immune
//! coupling, #195 sim-tme-3d, #197 cell-level biochem) share. Parallels
//! the [`crate::oxygen`] module structurally.
//!
//! **Physical model.** `pH(d) = ph_edge - (ph_edge - ph_core) · (1 - exp(-d/λ))`,
//! clamped to `[ph_core, ph_edge]`. Same form as the 2D `sim-tme` binary.
//! Defaults (sim-tme `PhConfig::default`): `ph_edge = 7.4`, `ph_core = 6.5`,
//! `lambda_ph_um = 120` (matches O₂ reference λ — both are
//! perfusion-limited diffusion lengths). All parameters are ESTIMATES;
//! tumor pH ranges from primary literature (Stubbs 2000, Vaupel 1989).
//!
//! **Stromal convention.** Cells outside the spheroid (`is_tumor == false`)
//! return `ph_edge` (well-perfused bulk tissue). Same convention as 2D.
//!
//! **API design — pure functions, no mutation.** Sim-tme's 2D
//! `apply_ph_gradient` mutates `cell.iron *= iron_multiplier` in place;
//! these 3D analogs return values. The consumer chooses what to do with
//! the field — mutate `cell.iron`, apply ion-trapping correction to
//! drug delivery, snapshot for analysis, etc. Same rationale as
//! [`crate::oxygen`].
//!
//! ## Quick example
//!
//! ```
//! use ferroptosis_core::grid::TumorGrid3D;
//! use ferroptosis_core::ph::{radial_ph_field, iron_multiplier_from_ph, ion_trap_factor_from_ph};
//!
//! let g = TumorGrid3D::generate(40, 40, 40, 20.0, 42);
//! let (ph_edge, ph_core, lambda) = (7.4, 6.5, 120.0);
//! let ph_field = radial_ph_field(&g, ph_edge, ph_core, lambda);
//! assert_eq!(ph_field.len(), g.cells.len());
//!
//! // Per-cell modulation: apply iron release + ion trapping at each pH.
//! let (iron_sens, ion_sens) = (1.5, 0.4);
//! for (idx, &local_ph) in ph_field.iter().enumerate() {
//!     let iron_mult = iron_multiplier_from_ph(local_ph, ph_edge, iron_sens);
//!     let drug_factor = ion_trap_factor_from_ph(local_ph, ph_edge, ion_sens);
//!     // consumer applies: cell.iron *= iron_mult; effective_drug = base_drug * drug_factor;
//!     let _ = (idx, iron_mult, drug_factor);
//! }
//! ```

use crate::grid::TumorGrid3D;

/// Lower-bound clamp on the ion-trap drug-availability factor. Matches
/// the hardcoded `[0.3, 1.0]` floor in sim-tme's `apply_ph_gradient`
/// (rationale: even at extreme acidity, some weak-base drug crosses
/// membranes; 30% bioavailability is a conservative empirical floor).
const ION_TRAP_FLOOR: f64 = 0.3;

/// Per-cell pH on a 3D spheroidal grid.
///
/// Returns a `Vec<f64>` of length `grid.cells.len()` in the flat order
/// (`r·cols·layers + c·layers + l`):
/// - **Stromal cells** (`is_tumor == false`): `ph_edge` (well-perfused
///   bulk tissue convention)
/// - **Tumor cells**: `pH(d) = ph_edge - (ph_edge - ph_core) · (1 - exp(-d/λ))`,
///   clamped to `[ph_core, ph_edge]`. `d` is `radial_depth_um.max(0.0)`
///   (defensive clip against floating-point roundoff at the surface).
///
/// Note the **opposite sign convention from O₂**: pH *decreases* with
/// depth (more acidic core), while O₂ *decreases* with depth toward
/// hypoxia. Both formulas are exponential but the pH form uses the
/// `(1 - exp)` transformation to start at `ph_edge` and asymptote to
/// `ph_core`, whereas O₂ uses plain `exp` to start at 1 and decay to 0.
///
/// **Validation.** All parameters must be finite; `ph_edge > ph_core`
/// (delta > 0); `lambda_ph_um > 0`. Invalid values trigger `debug_assert!`
/// in tests. **Release behavior is more complex than `oxygen.rs`'s** —
/// the `clamp(ph_core, ph_edge)` step has an unconditional internal
/// `assert!(min <= max)` (in `f64::clamp`) that runs *even in release*,
/// so some classes of bad input panic instead of producing undefined
/// values:
///
/// | Bad input | Per-cell behavior (release) |
/// |-----------|------------------------------|
/// | `λ = 0`, depth > 0 | `ph_core` (`exp(-∞) = 0`, formula collapses) |
/// | `λ = 0`, depth = 0 (surface cells) | `NaN` (`0/0 → NaN`, propagates) |
/// | `λ < 0` | `ph_edge` (`exp(+positive)` makes raw > `ph_edge`, clamp routes to `ph_edge` — NOT `ph_core`) |
/// | `λ = +∞` | `ph_edge` (`exp(0) = 1`, raw = `ph_edge`) |
/// | `λ = NaN` | per-cell `NaN` (NaN propagates through arithmetic into clamp's *self*; `f64::clamp` preserves NaN-self) |
/// | `ph_edge = NaN` or `ph_core = NaN` | **PANIC** (`f64::clamp` panics when min or max is NaN, since `min <= max` returns false for NaN comparisons) |
/// | `ph_edge < ph_core` | **PANIC** (`f64::clamp` panics on `min > max`) |
/// | `ph_edge == ph_core` | constant `ph_edge` everywhere (delta = 0; no panic since `min <= max` holds) — but `debug_assert` rejects this in tests as a likely configuration bug |
///
/// In short: callers loading parameters from untrusted sources must
/// validate at the boundary, particularly the `ph_edge > ph_core`
/// invariant and finite-NaN-free `ph_edge`/`ph_core`. Otherwise a
/// release run will panic on the first tumor-cell clamp.
///
/// **Cost.** O(N). Per-cell calls `radial_depth_um` which recomputes
/// geometry constants; same perf TODO as `radial_o2_field` (#194).
pub fn radial_ph_field(
    grid: &TumorGrid3D,
    ph_edge: f64,
    ph_core: f64,
    lambda_ph_um: f64,
) -> Vec<f64> {
    debug_assert_ph_inputs(ph_edge, ph_core, lambda_ph_um, "radial_ph_field");

    let mut out = Vec::with_capacity(grid.cells.len());
    for r in 0..grid.rows {
        for c in 0..grid.cols {
            for l in 0..grid.layers {
                let ph = if !grid.get(r, c, l).is_tumor {
                    ph_edge
                } else {
                    let depth_um = grid.radial_depth_um(r, c, l).max(0.0);
                    ph_at_depth(depth_um, ph_edge, ph_core, lambda_ph_um)
                };
                out.push(ph);
            }
        }
    }
    out
}

/// 2D analog of [`radial_ph_field`] — disc-shaped pH field over a
/// [`TumorGrid`] using the same first-order radial-decay formula.
///
/// Lifted from the inline math in sim-tme's `apply_ph_gradient` and
/// the cross-geometry library test in `ph::tests::matched_lambda_2d_vs_3d_acidic_fraction`
/// (#224 item 1b). Returns `ph_edge` for non-tumor cells (same
/// convention as the 3D version).
///
/// **Caveats:** same first-order single-source radial-decay
/// approximation as the 3D variant — see `radial_ph_field`'s rustdoc
/// for biological validity, calibration status (#190 → #194), and
/// debug-assert validation of `ph_edge > ph_core`.
pub fn radial_ph_field_2d(
    grid: &crate::grid::TumorGrid,
    ph_edge: f64,
    ph_core: f64,
    lambda_ph_um: f64,
) -> Vec<f64> {
    debug_assert_ph_inputs(ph_edge, ph_core, lambda_ph_um, "radial_ph_field_2d");

    let mut out = Vec::with_capacity(grid.cells.len());
    for r in 0..grid.rows {
        for c in 0..grid.cols {
            let ph = if !grid.get(r, c).is_tumor {
                ph_edge
            } else {
                let depth_um = grid.radial_depth_um(r, c).max(0.0);
                ph_at_depth(depth_um, ph_edge, ph_core, lambda_ph_um)
            };
            out.push(ph);
        }
    }
    out
}

/// Scalar pH-at-depth formula (first-order radial decay + clamp).
///
/// `ph = ph_edge − (ph_edge − ph_core) · (1 − exp(−depth/λ))`,
/// clamped to `[ph_core, ph_edge]`.
///
/// Dimensionality-agnostic primitive — both [`radial_ph_field`] (3D)
/// and [`radial_ph_field_2d`] use it. Sim-tme's `apply_ph_gradient`
/// calls it via the 2D field function after #224's lift.
///
/// **No tumor-cell check.** The caller decides whether to call this
/// (tumor cells) or default to `ph_edge` (non-tumor cells). Matches
/// the inline pattern in both wrappers.
#[inline]
pub fn ph_at_depth(depth_um: f64, ph_edge: f64, ph_core: f64, lambda_ph_um: f64) -> f64 {
    let delta = ph_edge - ph_core;
    let raw = ph_edge - delta * (1.0 - (-depth_um / lambda_ph_um).exp());
    raw.clamp(ph_core, ph_edge)
}

#[inline]
fn debug_assert_ph_inputs(ph_edge: f64, ph_core: f64, lambda_ph_um: f64, ctx: &'static str) {
    debug_assert!(
        ph_edge.is_finite() && ph_core.is_finite() && lambda_ph_um.is_finite(),
        "{ctx}: ph_edge, ph_core, lambda_ph_um must all be finite; got ph_edge={ph_edge}, ph_core={ph_core}, lambda_ph_um={lambda_ph_um}"
    );
    debug_assert!(
        ph_edge > ph_core,
        "{ctx}: ph_edge ({ph_edge}) must be strictly greater than ph_core ({ph_core}); equal or inverted values are likely a configuration bug"
    );
    debug_assert!(
        lambda_ph_um > 0.0,
        "{ctx}: lambda_ph_um must be > 0, got {lambda_ph_um}"
    );
}

/// Iron-release multiplier from local pH.
///
/// `multiplier = 1.0 + sensitivity · (ph_edge - local_ph)`
///
/// Models ferritin destabilization at low pH releasing labile Fe²⁺ that
/// fuels Fenton ROS production. At default sim-tme parameters
/// (`sensitivity = 1.5`, `ph_edge = 7.4`):
/// - `local_ph = ph_edge` → `1.0` (no perturbation)
/// - `local_ph = 6.5` → `1.0 + 1.5 · 0.9 = 2.35×`
///
/// **Pure scalar function.** Not clamped — at `local_ph > ph_edge` the
/// multiplier drops below 1.0; at extreme acidity it can grow large.
/// Matches sim-tme's unclamped 2D implementation; consumer responsible
/// for valid `local_ph` (typically from [`radial_ph_field`] which is
/// already clamped to `[ph_core, ph_edge]`).
///
/// **Sensitivity convention.** Must be `>= 0` (negative would invert the
/// biology: ferritin would *stabilize* at low pH, contradicting Yu et al.
/// 2017 and sim-tme's parameterization). `debug_assert` rejects negative
/// values; release accepts any sign and silently inverts the multiplier
/// curve.
///
/// **Release behavior for invalid inputs** (no clamping anywhere in this
/// helper, so propagation is straightforward):
///
/// | Bad input | Output |
/// |-----------|--------|
/// | `local_ph = NaN` | `NaN` (propagates through arithmetic; no clamp to catch it) |
/// | `ph_edge = NaN` | `NaN` |
/// | `sensitivity = NaN` | `NaN` |
/// | `sensitivity = +∞`, `local_ph != ph_edge` | `±∞` (sign matches `ph_edge - local_ph`) |
/// | `sensitivity = +∞`, `local_ph == ph_edge` | `NaN` (`∞ · 0 = NaN` in IEEE) |
///
/// Consumers passing `local_ph` from [`radial_ph_field`] (which is
/// clamped to `[ph_core, ph_edge]`) and finite sane sensitivities won't
/// hit any of these. The helpers are deliberately panic-free for
/// composability inside hot loops.
#[inline]
pub fn iron_multiplier_from_ph(local_ph: f64, ph_edge: f64, sensitivity: f64) -> f64 {
    debug_assert!(
        sensitivity >= 0.0,
        "iron_multiplier_from_ph: sensitivity must be >= 0 (negative inverts biology), got {sensitivity}"
    );
    1.0 + sensitivity * (ph_edge - local_ph)
}

/// Ion-trapping drug-availability factor for weak-base drugs.
///
/// `factor = 1.0 - sensitivity · (ph_edge - local_ph)`, clamped to
/// `[ION_TRAP_FLOOR, 1.0]` where `ION_TRAP_FLOOR = 0.3`.
///
/// Linearized Henderson-Hasselbalch over the narrow tumor pH range
/// (6.5-7.4). At low extracellular pH, weak-base drugs (e.g., RSL3) are
/// protonated and trapped outside the cell, reducing intracellular
/// concentration. At default sim-tme parameters (`sensitivity = 0.4`,
/// `ph_edge = 7.4`):
/// - `local_ph = ph_edge` → `1.0` (full bioavailability)
/// - `local_ph = 6.5` → `1.0 - 0.4 · 0.9 = 0.64` (36% drug lost)
/// - extreme acidity → clamped to `0.3` (model floor; matches sim-tme)
///
/// The `0.3` floor is hardcoded to match sim-tme's behavior. If you need
/// to tune it, the formula is exposed enough that you can implement your
/// own version — but reconsider whether the linearized model is still
/// valid in that range.
///
/// **Sensitivity convention.** Must be `>= 0` (negative would imply the
/// drug is *more* bioavailable at low pH, contradicting Henderson-
/// Hasselbalch for weak bases). `debug_assert` rejects negative values;
/// release silently inverts the curve.
///
/// **Release behavior for invalid inputs**:
///
/// | Bad input | Output |
/// |-----------|--------|
/// | `local_ph = NaN` | `NaN` (clamp preserves NaN-self per IEEE) |
/// | `ph_edge = NaN` | `NaN` (propagates through arithmetic into self, then clamp preserves) |
/// | `sensitivity = NaN` | `NaN` |
/// | `sensitivity = +∞` | clamped to `[ION_TRAP_FLOOR, 1.0]` (±∞ clamps to one of the bounds) |
///
/// Note the **asymmetry** with [`iron_multiplier_from_ph`]: the latter
/// has no clamp, so `±∞`/`NaN` flow straight through. This helper's
/// clamp bounds are constants (`ION_TRAP_FLOOR` and `1.0`), so unlike
/// [`radial_ph_field`] there's no panic risk — but NaN still propagates
/// through clamp's self-passthrough rule.
#[inline]
pub fn ion_trap_factor_from_ph(local_ph: f64, ph_edge: f64, sensitivity: f64) -> f64 {
    debug_assert!(
        sensitivity >= 0.0,
        "ion_trap_factor_from_ph: sensitivity must be >= 0 (negative inverts biology), got {sensitivity}"
    );
    let raw = 1.0 - sensitivity * (ph_edge - local_ph);
    raw.clamp(ION_TRAP_FLOOR, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::TumorGrid;

    /// Output length matches grid cell count.
    #[test]
    fn radial_ph_field_length_matches_grid() {
        let g = TumorGrid3D::generate(7, 5, 11, 20.0, 42);
        let ph = radial_ph_field(&g, 7.4, 6.5, 120.0);
        assert_eq!(ph.len(), g.cells.len());
        assert_eq!(ph.len(), 7 * 5 * 11);
    }

    /// Stromal cells (outside spheroid) return ph_edge regardless of λ.
    #[test]
    fn stromal_cells_get_edge_ph() {
        let g = TumorGrid3D::generate(10, 10, 10, 20.0, 42);
        for &lambda in &[10.0_f64, 100.0, 1_000.0] {
            let ph = radial_ph_field(&g, 7.4, 6.5, lambda);
            let mut stromal_count = 0usize;
            for (i, gc) in g.cells.iter().enumerate() {
                if !gc.is_tumor {
                    assert_eq!(
                        ph[i], 7.4,
                        "stromal cell at flat idx {i} got pH={} for λ={lambda}",
                        ph[i]
                    );
                    stromal_count += 1;
                }
            }
            assert!(stromal_count > 0, "expected some stromal cells in 10³ grid");
        }
    }

    /// Surface tumor cell (depth = 0 exactly) returns ph_edge.
    ///
    /// IEEE-exact: at `depth = 0`, `(-0.0/λ).exp() = exp(0) = 1.0`
    /// (IEEE-required-correct), then `1.0 - 1.0 = 0.0` exactly, then
    /// `ph_edge - delta · 0 = ph_edge - 0 = ph_edge` (the multiply-by-zero
    /// kills any rounding error in `delta`), then the clamp passes through
    /// since `ph_edge ∈ [ph_core, ph_edge]`. Strict equality is safe by
    /// construction, not by happy accident.
    #[test]
    fn surface_tumor_cell_has_edge_ph() {
        let g = TumorGrid3D::generate(20, 20, 20, 20.0, 42);
        // Cell (10, 10, 19) has radial_depth_um = 0 exactly (per
        // grid::tests_3d::radial_depth_at_sphere_surface_is_zero).
        assert!(
            g.get(10, 10, 19).is_tumor,
            "test precondition: this should be a tumor cell"
        );
        assert_eq!(g.radial_depth_um(10, 10, 19), 0.0, "test precondition");

        let ph = radial_ph_field(&g, 7.4, 6.5, 120.0);
        let flat = 10 * g.cols * g.layers + 10 * g.layers + 19;
        assert_eq!(ph[flat], 7.4);
    }

    /// Very large λ → `exp(-d/λ) ≈ 1` for all reasonable d, so `(1 - exp) ≈ 0`
    /// and `pH ≈ ph_edge` for all tumor cells. Asymptotic sanity check.
    #[test]
    fn very_large_lambda_gives_uniform_edge_ph() {
        let g = TumorGrid3D::generate(20, 20, 20, 20.0, 42);
        // λ orders of magnitude > grid size. Every cell should be ≈ 7.4.
        // For deepest cell (depth = 180), pH ≈ 7.4 - 0.9 · (1 - exp(-180/1e9))
        //                              ≈ 7.4 - 0.9 · 1.8e-7 ≈ 7.4 - 1.6e-7
        let ph = radial_ph_field(&g, 7.4, 6.5, 1e9);
        for (i, &v) in ph.iter().enumerate() {
            assert!(
                v > 7.4 - 1e-5,
                "cell {i} has pH={v}, expected ≈ 7.4 at huge λ"
            );
        }
    }

    /// Deep tumor (depth >> λ) approaches ph_core. With depth=180 µm,
    /// λ=20 µm: exp(-9) ≈ 1.23e-4 ≈ 0; pH ≈ ph_core + delta · 1.23e-4
    /// ≈ ph_core. Tight tolerance via libm.
    #[test]
    fn deep_tumor_approaches_core_ph() {
        let g = TumorGrid3D::generate(20, 20, 20, 20.0, 42);
        let (ph_edge, ph_core, lambda) = (7.4, 6.5, 20.0_f64);
        let ph = radial_ph_field(&g, ph_edge, ph_core, lambda);

        // Deepest cell is the center (10, 10, 10) with depth ≈ 180 µm.
        let flat = 10 * g.cols * g.layers + 10 * g.layers + 10;
        let center_ph = ph[flat];
        // (1 - exp(-9)) ≈ 0.99988; pH ≈ 7.4 - 0.9 · 0.99988 ≈ 6.5001
        assert!(
            (center_ph - ph_core).abs() < 1e-3,
            "deepest cell pH = {center_ph}, expected ≈ ph_core = {ph_core}"
        );
    }

    /// pH monotonically *decreases* with depth (opposite sign from O₂).
    #[test]
    fn ph_decreases_monotonically_with_depth() {
        let g = TumorGrid3D::generate(20, 20, 20, 20.0, 42);
        let ph = radial_ph_field(&g, 7.4, 6.5, 100.0);

        // Walk l from 0 (outside) to 10 (center) along (r=10, c=10).
        // For tumor cells, pH should be non-increasing as l → 10.
        let mut prev: Option<f64> = None;
        let mut tumor_samples = 0usize;
        for l in 0..=10 {
            let flat = 10 * g.cols * g.layers + 10 * g.layers + l;
            if !g.cells[flat].is_tumor {
                continue;
            }
            let cur = ph[flat];
            if let Some(p) = prev {
                assert!(
                    cur <= p + 1e-12,
                    "pH not monotone-decreasing toward center at l={l}: prev={p}, cur={cur}"
                );
            }
            prev = Some(cur);
            tumor_samples += 1;
        }
        assert!(
            tumor_samples >= 5,
            "expected several tumor cells on the radial line"
        );
    }

    /// pH stays within `[ph_core, ph_edge]` for every cell (clamping
    /// works). With pathological-but-allowed sensitivities the formula
    /// CAN go out of range; the clamp guarantees the output doesn't.
    #[test]
    fn ph_stays_in_clamped_range() {
        let g = TumorGrid3D::generate(20, 20, 20, 20.0, 42);
        let (ph_edge, ph_core) = (7.4, 6.5);
        let ph = radial_ph_field(&g, ph_edge, ph_core, 100.0);
        for (i, &v) in ph.iter().enumerate() {
            assert!(
                v >= ph_core && v <= ph_edge,
                "cell {i} pH={v} out of [{ph_core}, {ph_edge}]"
            );
        }
    }

    /// `iron_multiplier_from_ph(ph_edge, ph_edge, _) = 1.0` exactly.
    /// IEEE-exact: `(ph_edge - ph_edge) = 0`, `sens · 0 = 0`, `1.0 + 0 = 1.0`.
    /// Multiply-by-zero kills any rounding in `sens` or `ph_edge`.
    #[test]
    fn iron_multiplier_at_edge_is_one() {
        for &edge in &[7.4_f64, 7.0, 8.0] {
            for &sens in &[0.0_f64, 1.5, 100.0] {
                assert_eq!(iron_multiplier_from_ph(edge, edge, sens), 1.0);
            }
        }
    }

    /// `iron_multiplier_from_ph(ph_core, ph_edge, sens)` at defaults gives
    /// 2.35× (sim-tme docstring claim). Tolerance: `7.4 - 6.5` is not
    /// IEEE-exact (7.4 has repeating binary expansion); also `1.5 · 0.9`
    /// and the final `1.0 + ...` accumulate roundoff. Use tight tolerance.
    #[test]
    fn iron_multiplier_at_core_matches_sim_tme_default() {
        let m = iron_multiplier_from_ph(6.5, 7.4, 1.5);
        // Expected: 1 + 1.5 · 0.9 = 2.35
        assert!(
            (m - 2.35).abs() < 1e-9,
            "iron multiplier at pH 6.5 = {m}, expected 2.35 ± 1e-9"
        );
    }

    /// `ion_trap_factor_from_ph(ph_edge, ph_edge, _) = 1.0` exactly,
    /// and the clamp passes through (1.0 is within [0.3, 1.0]).
    #[test]
    fn ion_trap_at_edge_is_one() {
        for &edge in &[7.4_f64, 7.0, 8.0] {
            for &sens in &[0.0_f64, 0.4, 1.0] {
                assert_eq!(ion_trap_factor_from_ph(edge, edge, sens), 1.0);
            }
        }
    }

    /// `ion_trap_factor_from_ph(6.5, 7.4, 0.4)` at defaults gives 0.64
    /// (sim-tme docstring claim). Tight tolerance.
    #[test]
    fn ion_trap_at_core_matches_sim_tme_default() {
        let f = ion_trap_factor_from_ph(6.5, 7.4, 0.4);
        // Expected: 1 - 0.4 · 0.9 = 0.64
        assert!(
            (f - 0.64).abs() < 1e-9,
            "ion-trap factor at pH 6.5 = {f}, expected 0.64 ± 1e-9"
        );
    }

    /// Ion-trap factor clamps to the floor `ION_TRAP_FLOOR = 0.3` at
    /// extreme sensitivities. With sens=2.0, delta=0.9: raw = -0.8,
    /// clamped to 0.3.
    #[test]
    fn ion_trap_clamps_to_floor() {
        let f = ion_trap_factor_from_ph(6.5, 7.4, 2.0);
        assert_eq!(f, ION_TRAP_FLOOR, "extreme sens should clamp at floor");
        // Lock the floor to sim-tme's 0.3 convention as a contract assertion
        // (defends against an accidental change to the private constant —
        // public callers can rely on the floor being exactly 0.3).
        assert_eq!(ION_TRAP_FLOOR, 0.3);
    }

    /// Ion-trap factor clamps to the upper bound `1.0` when the raw
    /// formula exceeds 1 (i.e., `local_ph > ph_edge` — alkalosis above
    /// the reference, unusual but representable). Validates the upper
    /// half of the `[0.3, 1.0]` clamp contract.
    #[test]
    fn ion_trap_clamps_to_upper_bound() {
        // local_ph above ph_edge → ph_edge - local_ph < 0 → raw > 1.
        // With sens=0.4, ph_edge=7.4, local_ph=8.4: raw = 1.0 - 0.4·(-1.0) = 1.4.
        let f = ion_trap_factor_from_ph(8.4, 7.4, 0.4);
        assert_eq!(f, 1.0, "raw > 1 should clamp to upper bound 1.0");
    }

    /// **Property test (reviewer ask)**: for every tumor cell in a
    /// `radial_ph_field` output, `iron_multiplier_from_ph` with a
    /// non-negative sensitivity must return `>= 1.0` (i.e., low pH
    /// either keeps iron at baseline or releases more — never reduces
    /// it). Locks the biology invariant that the docstring promises
    /// only at the endpoints.
    #[test]
    fn iron_multiplier_is_at_least_one_for_field_inputs() {
        let g = TumorGrid3D::generate(20, 20, 20, 20.0, 42);
        let (ph_edge, ph_core) = (7.4, 6.5);
        let ph_field = radial_ph_field(&g, ph_edge, ph_core, 120.0);
        let sensitivity = 1.5;

        for (i, gc) in g.cells.iter().enumerate() {
            if !gc.is_tumor {
                continue;
            }
            let m = iron_multiplier_from_ph(ph_field[i], ph_edge, sensitivity);
            assert!(
                m >= 1.0,
                "iron multiplier should be >= 1.0 for tumor cell at flat idx {i}: \
                 local_ph={}, multiplier={m}",
                ph_field[i]
            );
        }
    }

    /// **Acceptance criterion #3** (issue #190; loose interpretation):
    /// "Compare 2D vs 3D pH kill effects at matched parameters." We
    /// compare pH-field *shape* (acidic-fraction with `pH < 6.8` as the
    /// bucket threshold), not actual kill simulation — kill comparison
    /// would require full multi-step runs and is out of scope here.
    ///
    /// Construction: matched 40×40 / 40³ grids, cell_size=20 µm →
    /// `tumor_radius_um = 360`. Defaults: ph_edge=7.4, ph_core=6.5, λ=120.
    /// Threshold: pH < 6.8 — conventional cutoff for moderate tumor
    /// acidosis (between the typical mild-acidity boundary ~7.0 and the
    /// severe-core regime ~6.5). Chosen to give non-empty acidic regions
    /// in both geometries at default λ.
    ///
    /// Drift-proof: both the 2D and 3D arms call library helpers
    /// ([`radial_ph_field_2d`] and [`radial_ph_field`], #224 item 1b)
    /// that share the same `ph_at_depth` scalar primitive. sim-tme's
    /// `apply_ph_gradient` also routes through `radial_ph_field_2d`,
    /// so the test and the binary can't diverge.
    ///
    /// **Expected result** (same cubic-vs-quadratic geometry as #187):
    /// at matched R and λ, the 3D acidic-volume fraction is *smaller*
    /// than the 2D acidic-area fraction (raising a fraction in (0,1) to
    /// a higher exponent makes it smaller). The biological "3D more
    /// acidic" observation reflects vasculature effects (smaller
    /// effective λ), not pure geometry.
    #[test]
    fn matched_lambda_2d_vs_3d_acidic_fraction() {
        let cell_size_um = 20.0;
        let (ph_edge, ph_core, lambda) = (7.4, 6.5, 120.0);
        let acidic_threshold_ph = 6.8;

        // 2D: 40×40 grid via the lifted `radial_ph_field_2d`
        // (#224 item 1b) — same first-order radial-decay formula as
        // the 3D version, sharing `ph_at_depth` as the scalar primitive
        // so a future tuning of the formula touches one place.
        let g2 = TumorGrid::generate(40, 40, cell_size_um, 42);
        let ph_2d_field = radial_ph_field_2d(&g2, ph_edge, ph_core, lambda);
        let (mut acidic_2, mut total_2) = (0usize, 0usize);
        for (i, gc) in g2.cells.iter().enumerate() {
            if !gc.is_tumor {
                continue;
            }
            total_2 += 1;
            if ph_2d_field[i] < acidic_threshold_ph {
                acidic_2 += 1;
            }
        }
        let frac_2d = acidic_2 as f64 / total_2 as f64;

        // 3D: 40³ via radial_ph_field.
        let g3 = TumorGrid3D::generate(40, 40, 40, cell_size_um, 42);
        let ph_3d = radial_ph_field(&g3, ph_edge, ph_core, lambda);
        let (mut acidic_3, mut total_3) = (0usize, 0usize);
        for (i, gc) in g3.cells.iter().enumerate() {
            if !gc.is_tumor {
                continue;
            }
            total_3 += 1;
            if ph_3d[i] < acidic_threshold_ph {
                acidic_3 += 1;
            }
        }
        let frac_3d = acidic_3 as f64 / total_3 as f64;

        assert!(
            acidic_2 > 0 && acidic_3 > 0,
            "test precondition: expected non-empty acidic regions in both geometries; got 2D={acidic_2}, 3D={acidic_3}"
        );

        // 3D acidic-core fraction is SMALLER than 2D at matched parameters
        // (same cubic-vs-quadratic geometry as the #187 O₂ cross-geometry
        // test). Vaupel-1989-style "3D more acidic" reflects biology
        // (poor vasculature → smaller effective λ), not pure geometry.
        assert!(
            frac_3d < frac_2d,
            "3D acidic-volume fraction ({frac_3d:.4}) should be SMALLER than 2D acidic-area \
             fraction ({frac_2d:.4}) at matched λ — pure-geometry consequence (see test docstring)"
        );
    }

    // ============================================================
    // Scalar `ph_at_depth` + 2D field (#224 item 1b lift).
    // ============================================================

    /// `ph_at_depth` at `depth = 0` returns `ph_edge` exactly. Same
    /// IEEE-exact rationale as the surface case for `radial_ph_field`:
    /// `exp(-0/λ) = 1.0` exact, so `ph_edge - (ph_edge - ph_core) * 0 = ph_edge`.
    #[test]
    fn ph_at_depth_surface_is_edge() {
        assert_eq!(ph_at_depth(0.0, 7.4, 6.5, 120.0), 7.4);
    }

    /// `ph_at_depth` clamps at `ph_core` for arbitrarily deep cells —
    /// no overshoot below `ph_core` regardless of how far inside.
    #[test]
    fn ph_at_depth_deep_clamps_at_core() {
        // At depth = 100×λ, exp(-100) ≈ 0 → raw → ph_edge - delta = ph_core.
        let ph = ph_at_depth(100.0 * 120.0, 7.4, 6.5, 120.0);
        assert!((ph - 6.5).abs() < 1e-9, "got {ph}");
    }

    /// `radial_ph_field_2d` length matches the 2D grid cell count and
    /// non-tumor cells get `ph_edge`. Mirrors the 3D test contract.
    #[test]
    fn radial_ph_field_2d_length_and_non_tumor_default() {
        let g = TumorGrid::generate(40, 40, 20.0, 42);
        let ph = radial_ph_field_2d(&g, 7.4, 6.5, 120.0);
        assert_eq!(ph.len(), g.cells.len());
        for (i, gc) in g.cells.iter().enumerate() {
            if !gc.is_tumor {
                assert_eq!(ph[i], 7.4, "non-tumor cell {i} did not return ph_edge");
            }
        }
    }

    /// Bit-identical to sim-tme's binary-local `apply_ph_gradient` pH
    /// formula (pre-clamp). Locks the lift; a future refactor of either
    /// side can't diverge silently.
    #[test]
    fn radial_ph_field_2d_matches_sim_tme_apply_ph_gradient_formula() {
        let g = TumorGrid::generate(40, 40, 20.0, 42);
        let (ph_edge, ph_core, lambda) = (7.4, 6.5, 120.0);
        let delta = ph_edge - ph_core;
        let lib_field = radial_ph_field_2d(&g, ph_edge, ph_core, lambda);
        for (idx, gc) in g.cells.iter().enumerate() {
            let (r, c) = (idx / g.cols, idx % g.cols);
            let expected = if !gc.is_tumor {
                ph_edge
            } else {
                let depth_um = g.radial_depth_um(r, c).max(0.0);
                let raw = ph_edge - delta * (1.0 - (-depth_um / lambda).exp());
                raw.clamp(ph_core, ph_edge)
            };
            assert_eq!(lib_field[idx], expected, "cell ({r},{c})");
        }
    }
}
