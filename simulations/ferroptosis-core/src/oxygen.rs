//! 3D radial oxygen gradients for spheroid tumors.
//!
//! In a spheroid, oxygen diffuses inward from the well-perfused surface;
//! deeper cells experience progressively less O₂, with a hypoxic (sometimes
//! necrotic) core. This module provides the field-computation and
//! zone-census primitives that downstream 3D consumers (#188 immune
//! coupling, #191 vasculature, #197 cell-level biochem) share.
//!
//! **Physical model.** O₂ at radial depth `d` from the spheroid surface is
//! modelled as `O2(d) = exp(-d/λ)` clamped to `[0, 1]` — a first-order
//! exponential-decay approximation. The *exact* spheroidal diffusion
//! solution (Riley equation) is
//! `O2(r) ∝ sinh(r·√(k/D)) / (r·sinh(R·√(k/D)))`, which the exponential
//! approximates well near the surface but underestimates deep in the
//! core (more shoulder, less tail). Same approximation level as the 2D
//! `sim-tme` binary; consistent across geometries. (Note: Krogh's
//! original 1919 work derived the *cylinder* model for O₂ transport
//! around capillaries — the spheroid form here is the standard
//! exponential generalization, not Krogh's cylinder per se.) Ref:
//! Krogh A, *J Physiol* 1919 (cylinder model); Vaupel P, *Cancer Res*
//! 1989 (spheroid O₂ and hypoxia).
//!
//! **Stromal convention.** Cells outside the spheroid (negative
//! `radial_depth_um`, `is_tumor == false`) return `O2 = 1.0`. They
//! represent bulk normal tissue near vasculature — well-oxygenated by
//! convention, not by computed diffusion.
//!
//! **API design — pure functions, no mutation.** `sim-tme`'s 2D
//! `apply_o2_gradient` mutates `cell.basal_ros *= o2_factor` in place;
//! these 3D analogs return values instead. The pure form composes
//! cleanly with O₂ cycling (re-call per step with a different λ; no
//! `original_ros` snapshot needed) and decouples O₂ computation from
//! consumer state-management. The consumer chooses whether to mutate
//! `basal_ros`, snapshot for cycling, or feed the field into something
//! else entirely.
//!
//! ## Quick example
//!
//! ```
//! use ferroptosis_core::grid::TumorGrid3D;
//! use ferroptosis_core::oxygen::{radial_o2_field, radial_o2_zone_kill_rates};
//!
//! let g = TumorGrid3D::generate(40, 40, 40, 20.0, 42);
//! let lambda_um = 100.0; // O2 penetration length (Vaupel 1989, ~100-150 µm)
//!
//! // Per-cell O2 factor (length matches g.cells).
//! let o2 = radial_o2_field(&g, lambda_um);
//! assert_eq!(o2.len(), g.cells.len());
//!
//! // Zone-based dead-cell census (no cells dead in a fresh grid).
//! let (norm, trans, hyp) = radial_o2_zone_kill_rates(&g, lambda_um);
//! assert_eq!((norm, trans, hyp), (0.0, 0.0, 0.0));
//!
//! // Cycling = consumer pattern; re-call with alternating λ per step.
//! let o2_normoxic = radial_o2_field(&g, lambda_um * 1.5);
//! let o2_hypoxic = radial_o2_field(&g, lambda_um * 0.5);
//! assert!(o2_normoxic.iter().zip(o2_hypoxic.iter()).any(|(a, b)| a != b));
//! ```

use crate::grid::{RadialDepthGeom, TumorGrid3D};

/// Per-cell O₂ factor on a 3D spheroidal grid.
///
/// Returns a `Vec<f64>` of length `grid.cells.len()` in the same flat order
/// (`r·cols·layers + c·layers + l`). Each entry is in `[0, 1]`:
/// - **Stromal cells** (`is_tumor == false`): `1.0` (well-oxygenated bulk
///   tissue, by convention)
/// - **Tumor cells**: `exp(-d/λ)` clamped to `[0, 1]`, where `d` is the
///   radial depth from the spheroid surface in µm (per
///   [`TumorGrid3D::radial_depth_um`]). For tumor cells `d ≥ 0`; the
///   `.max(0.0)` clip on the depth is defensive against floating-point
///   roundoff right at the boundary.
///
/// **Validation.** `λ` must be finite and strictly positive. Invalid
/// values trigger `debug_assert!` in tests; release builds produce
/// undefined values per the table below — callers loading `λ` from
/// untrusted sources should validate at the boundary. Matches the
/// validation posture of [`crate::physics`].
///
/// | `λ` value | Per-cell output (release) |
/// |-----------|---------------------------|
/// | finite, `> 0` | valid `[0, 1]` |
/// | `0.0` | `0.0` (`exp(-∞)` then clamp) |
/// | `< 0` | `1.0` (`exp(+depth/\|λ\|)` saturates, then clamp) |
/// | `NaN` | `NaN` (`f64::clamp` propagates NaN, does not clip) |
/// | `+∞` | `1.0` (`exp(-depth/∞) = 1.0`) |
///
/// **Cost.** O(N) for N = `grid.cells.len()`. The dimension-only depth
/// geometry (`center_{r,c,l}`, `tumor_radius`) is hoisted once via
/// [`RadialDepthGeom`] (#289) instead of being recomputed per cell; the
/// per-cell `depth_um` is bit-identical to [`TumorGrid3D::radial_depth_um`].
pub fn radial_o2_field(grid: &TumorGrid3D, lambda_um: f64) -> Vec<f64> {
    debug_assert!(
        lambda_um.is_finite() && lambda_um > 0.0,
        "radial_o2_field: lambda_um must be finite and positive, got {lambda_um}"
    );

    // Hoist the dimension-only depth geometry once (#289); per-cell
    // `geom.depth_um` is bit-identical to `grid.radial_depth_um`.
    let geom = RadialDepthGeom::new(grid);
    let mut out = Vec::with_capacity(grid.cells.len());
    for r in 0..grid.rows {
        for c in 0..grid.cols {
            for l in 0..grid.layers {
                let factor = if !grid.get(r, c, l).is_tumor {
                    1.0
                } else {
                    let depth_um = geom.depth_um(r, c, l).max(0.0);
                    (-depth_um / lambda_um).exp().clamp(0.0, 1.0)
                };
                out.push(factor);
            }
        }
    }
    out
}

/// Oxygen-dependence scaling factor for exogenous (SDT/PDT) ROS yield (#336).
///
/// SDT and PDT generate ROS via a sono/photosensitizer. The dominant clinical
/// mechanism is Type II (singlet oxygen via energy transfer), which is
/// **oxygen-dependent**: the lead clinical sonodynamic agent SONALA-001 (5-ALA
/// derived protoporphyrin IX, activated by MR-guided focused ultrasound) is
/// Type II, and its first-in-human recurrent-high-grade-glioma trial reported
/// only modest cell death (Sanai et al., Science Translational Medicine 2025,
/// DOI 10.1126/scitranslmed.ads5813). The default model treats the exogenous
/// ROS burst as O2-independent, an optimistic upper bound (manuscript §7.1);
/// this helper lets a consumer make a configurable fraction of it O2-dependent
/// so the contested hypoxia leg can be re-examined.
///
/// `o2_supply` is the local relative O2 availability (the same factor that
/// scales `cell.basal_ros`). `dependence` is the O2-dependent fraction of the
/// exogenous ROS yield (the "Type II fraction"): `0.0` = fully O2-independent
/// (the historical default, returns `1.0`, byte-identical), `1.0` = fully
/// O2-dependent (Type II, scales linearly with O2).
///
/// Returns `(1 − dependence) + dependence·o2_supply`, clamped to `[0, 1]`. At
/// `dependence = 0` the exogenous ROS is unaffected by hypoxia (the SDT hypoxia
/// advantage is maximal); raising it shrinks that advantage toward the
/// pharmacologic case.
pub fn o2_dependent_exo_factor(o2_supply: f64, dependence: f64) -> f64 {
    let d = dependence.clamp(0.0, 1.0);
    let s = o2_supply.clamp(0.0, 1.0);
    (1.0 - d + d * s).clamp(0.0, 1.0)
}

/// Dead-cell rate for each of three O₂-defined concentric shells.
///
/// Zone semantics (matches `sim-tme`'s 2D
/// [`zone_kill_rates`](https://github.com/ELares/cancer_research/blob/main/simulations/sim-tme/src/main.rs)
/// exactly — same depth thresholds, geometry-agnostic):
/// - **Normoxic shell**: `depth_um ∈ [0, shell_depth_um)`
///   (cells within `shell_depth_um` of the surface; `O2 > 1/e ≈ 0.37`
///   if `shell_depth_um == λ`)
/// - **Transition zone**: `depth_um ∈ [shell_depth_um, 3·shell_depth_um)`
/// - **Hypoxic core**: `depth_um ≥ 3·shell_depth_um`
///   (`O2 < e⁻³ ≈ 0.05` at the threshold if `shell_depth_um == λ`)
///
/// Returns `(normoxic_rate, transition_rate, hypoxic_rate)`. Each rate is
/// `dead_in_zone / total_tumor_in_zone`, or `0.0` for an empty zone (no
/// division-by-zero panic). **Stromal cells are excluded** from all counts.
///
/// **Validation.** `shell_depth_um` must be `≥ 0`. Zero is allowed (it
/// collapses normoxic+transition into empty zones and routes all tumor
/// cells into hypoxic — biologically odd but mathematically defined).
pub fn radial_o2_zone_kill_rates(grid: &TumorGrid3D, shell_depth_um: f64) -> (f64, f64, f64) {
    debug_assert!(
        shell_depth_um >= 0.0 && shell_depth_um.is_finite(),
        "radial_o2_zone_kill_rates: shell_depth_um must be finite and non-negative, got {shell_depth_um}"
    );

    let deep_threshold_um = shell_depth_um * 3.0;

    let (mut norm_dead, mut norm_total) = (0usize, 0usize);
    let (mut trans_dead, mut trans_total) = (0usize, 0usize);
    let (mut hyp_dead, mut hyp_total) = (0usize, 0usize);

    for r in 0..grid.rows {
        for c in 0..grid.cols {
            for l in 0..grid.layers {
                let gc = grid.get(r, c, l);
                if !gc.is_tumor {
                    continue;
                }
                // For tumor cells, radial_depth_um >= 0 by the
                // generate/radial_depth_um geometry contract (locked down
                // by grid::tests_3d::radial_depth_sign_agrees_with_generated_is_tumor).
                let depth_um = grid.radial_depth_um(r, c, l);

                let (dead_count, total_count) = if depth_um < shell_depth_um {
                    (&mut norm_dead, &mut norm_total)
                } else if depth_um < deep_threshold_um {
                    (&mut trans_dead, &mut trans_total)
                } else {
                    (&mut hyp_dead, &mut hyp_total)
                };

                *total_count += 1;
                if gc.state.dead {
                    *dead_count += 1;
                }
            }
        }
    }

    let rate = |d: usize, t: usize| if t > 0 { d as f64 / t as f64 } else { 0.0 };
    (
        rate(norm_dead, norm_total),
        rate(trans_dead, trans_total),
        rate(hyp_dead, hyp_total),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::TumorGrid;

    /// #336: the exo-ROS O2-dependence factor is identity at dependence 0 (the
    /// historical, byte-identical default) and scales linearly with O2 at
    /// dependence 1 (fully Type II / O2-dependent), monotone decreasing in
    /// dependence under hypoxia, and clamped on out-of-range inputs.
    #[test]
    fn o2_dependent_exo_factor_identity_at_zero_and_scales_at_one() {
        // dependence = 0 (default): always 1.0 regardless of O2 (no hypoxia effect).
        assert_eq!(o2_dependent_exo_factor(0.05, 0.0), 1.0);
        assert_eq!(o2_dependent_exo_factor(1.0, 0.0), 1.0);
        // dependence = 1 (fully Type II): factor equals the O2 supply.
        assert_eq!(o2_dependent_exo_factor(0.05, 1.0), 0.05);
        assert_eq!(o2_dependent_exo_factor(1.0, 1.0), 1.0);
        // Monotone decreasing in dependence at fixed low O2.
        let hyp = 0.1;
        assert!(o2_dependent_exo_factor(hyp, 0.5) < o2_dependent_exo_factor(hyp, 0.0));
        assert!(o2_dependent_exo_factor(hyp, 1.0) < o2_dependent_exo_factor(hyp, 0.5));
        // Out-of-range inputs are clamped, never panicking.
        assert_eq!(o2_dependent_exo_factor(2.0, 1.0), 1.0);
        assert_eq!(o2_dependent_exo_factor(-1.0, 1.0), 0.0);
        assert_eq!(
            o2_dependent_exo_factor(0.5, 2.0),
            o2_dependent_exo_factor(0.5, 1.0)
        );
    }

    /// #289: the hoisted-geometry field must be bit-for-bit identical to a
    /// per-cell baseline that uses the canonical `TumorGrid3D::radial_depth_um`.
    /// Locks the contract that hoisting `RadialDepthGeom` out of the loop is a
    /// timing-only change (so the `summary.json` matrix is unaffected).
    #[test]
    fn radial_o2_field_matches_per_cell_depth_baseline() {
        let g = TumorGrid3D::generate(20, 20, 20, 20.0, 42);
        let lambda = 100.0;
        let hoisted = radial_o2_field(&g, lambda);
        let mut baseline = Vec::with_capacity(g.cells.len());
        for r in 0..g.rows {
            for c in 0..g.cols {
                for l in 0..g.layers {
                    let factor = if !g.get(r, c, l).is_tumor {
                        1.0
                    } else {
                        let depth_um = g.radial_depth_um(r, c, l).max(0.0);
                        (-depth_um / lambda).exp().clamp(0.0, 1.0)
                    };
                    baseline.push(factor);
                }
            }
        }
        assert_eq!(hoisted, baseline, "hoisted O2 field must be bit-identical");
    }

    /// Output length matches grid cell count. Guards against any future
    /// refactor that diverges the iteration order or skips cells.
    #[test]
    fn radial_o2_field_length_matches_grid() {
        let g = TumorGrid3D::generate(7, 5, 11, 20.0, 42);
        let o2 = radial_o2_field(&g, 100.0);
        assert_eq!(o2.len(), g.cells.len());
        assert_eq!(o2.len(), 7 * 5 * 11);
    }

    /// Stromal cells (outside spheroid) get 1.0 regardless of λ. Walks
    /// every cell and confirms is_tumor == false ⇒ o2 == 1.0 exactly.
    #[test]
    fn stromal_cells_get_full_oxygen() {
        let g = TumorGrid3D::generate(10, 10, 10, 20.0, 42);
        for &lambda in &[10.0_f64, 100.0, 1_000.0] {
            let o2 = radial_o2_field(&g, lambda);
            let mut stromal_count = 0usize;
            for (i, gc) in g.cells.iter().enumerate() {
                if !gc.is_tumor {
                    assert_eq!(
                        o2[i], 1.0,
                        "stromal cell at flat idx {i} got O2={} for λ={lambda}",
                        o2[i]
                    );
                    stromal_count += 1;
                }
            }
            assert!(stromal_count > 0, "expected some stromal cells in 10³ grid");
        }
    }

    /// Surface tumor cell (depth=0 exactly) returns O2 = 1.0. IEEE-exact:
    /// `exp(0) = 1.0` is required-correct, `1.0.clamp(0,1) = 1.0` is
    /// trivially exact. Reuses the same surface cell coordinates as
    /// `grid::tests_3d::radial_depth_at_sphere_surface_is_zero`.
    #[test]
    fn surface_tumor_cell_has_full_oxygen() {
        let g = TumorGrid3D::generate(20, 20, 20, 20.0, 42);
        // Cell (10, 10, 19): radial_depth_um = 0 exactly (per the
        // surface-boundary test in grid::tests_3d). is_tumor = true
        // because dist == tumor_radius makes it tumor (not strictly outside).
        assert!(
            g.get(10, 10, 19).is_tumor,
            "test precondition: this should be a tumor cell"
        );
        assert_eq!(g.radial_depth_um(10, 10, 19), 0.0, "test precondition");

        let o2 = radial_o2_field(&g, 100.0);
        let flat = 10 * g.cols * g.layers + 10 * g.layers + 19;
        assert_eq!(o2[flat], 1.0);
    }

    /// At depth = λ, O₂ = exp(-1) = 1/e. Tight tolerance since the libm
    /// implementation of exp is well-tested but not bit-portable across
    /// platforms.
    #[test]
    fn depth_equals_lambda_gives_one_over_e() {
        // 20³ grid, cell_size=20 → tumor_radius_lattice = 9.0,
        // tumor_radius_um = 180. Center is at (10,10,10), depth there = 180.
        // We want a cell whose depth ≈ λ; pick λ = 180 - dist_um, so dist
        // along an axis = (180 - λ)/cell_size lattice units from center.
        let g = TumorGrid3D::generate(20, 20, 20, 20.0, 42);

        // Pick cell (10, 10, 14): dr=dc=0, dl=4 → dist=4 lattice = 80 µm
        // from center. radial_depth_um = (9.0 - 4.0) * 20 = 100 µm.
        let depth_um = g.radial_depth_um(10, 10, 14);
        assert_eq!(
            depth_um, 100.0,
            "test precondition: depth should be 100 µm exactly"
        );

        let o2 = radial_o2_field(&g, 100.0);
        let flat = 10 * g.cols * g.layers + 10 * g.layers + 14;
        let one_over_e = (-1.0_f64).exp();
        assert!(
            (o2[flat] - one_over_e).abs() < 1e-12,
            "depth==λ should give exp(-1)≈{}, got {}",
            one_over_e,
            o2[flat]
        );
    }

    /// λ → very large → O₂ ≈ 1 everywhere (including the deepest tumor
    /// cell). Asymptotic sanity check on the exp formula.
    #[test]
    fn large_lambda_gives_uniform_oxygen() {
        let g = TumorGrid3D::generate(20, 20, 20, 20.0, 42);
        // λ orders of magnitude > grid size. Every cell should be > 0.999.
        // For the deepest cell at depth=180, o2 = exp(-180/1e9) ≈ 1 - 1.8e-7.
        let o2 = radial_o2_field(&g, 1e9);
        for (i, &v) in o2.iter().enumerate() {
            assert!(v > 0.999, "cell {i} has O2={v} which is < 0.999 at huge λ");
        }
    }

    /// O₂ decreases monotonically with depth along a radial line. Picks the
    /// l-axis through the center and asserts non-increasing O₂ as l
    /// approaches the center from one side.
    #[test]
    fn oxygen_decreases_monotonically_with_depth() {
        let g = TumorGrid3D::generate(20, 20, 20, 20.0, 42);
        let o2 = radial_o2_field(&g, 100.0);

        // Walk l from 0 (outside) to 10 (center) along (r=10, c=10). For
        // tumor cells, O2 should be non-increasing as l → 10.
        let mut prev: Option<f64> = None;
        let mut tumor_samples = 0usize;
        for l in 0..=10 {
            let flat = 10 * g.cols * g.layers + 10 * g.layers + l;
            if !g.cells[flat].is_tumor {
                continue;
            }
            let cur = o2[flat];
            if let Some(p) = prev {
                assert!(
                    cur <= p + 1e-12, // tiny epsilon for floating-point safety
                    "O2 not monotone-decreasing toward center at l={l}: prev={p}, cur={cur}"
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

    /// Two calls to `radial_o2_field` with distinct `λ` produce distinct
    /// fields — sanity check that `λ` actually drives the math.
    ///
    /// **Not a cycling test** (the function is deterministic in λ, so
    /// this trivially holds); see
    /// [`o2_cycling_pattern_demonstrates_consumer_loop`] below for the
    /// actual square-wave-cycling demonstration that exercises AC #2.
    #[test]
    fn field_varies_with_lambda() {
        let g = TumorGrid3D::generate(20, 20, 20, 20.0, 42);
        let normoxic = radial_o2_field(&g, 150.0);
        let hypoxic = radial_o2_field(&g, 50.0);
        let n_diff = normoxic
            .iter()
            .zip(hypoxic.iter())
            .filter(|(a, b)| a != b)
            .count();
        assert!(
            n_diff > 0,
            "expected at least one tumor cell to differ between λ=150 and λ=50"
        );
        // Stromal cells are 1.0 in both. Tumor cells should differ — there
        // should be many of them.
        let n_tumor = g.cells.iter().filter(|gc| gc.is_tumor).count();
        // Some interior tumor cells right at the surface have depth=0 and
        // therefore O2=1 in both fields — they wouldn't differ. The
        // expected count is "less than all tumor cells".
        assert!(n_diff <= n_tumor);
        // But for a 20³ grid, the vast majority of tumor cells have
        // depth > 0 and should differ.
        assert!(
            n_diff > n_tumor / 2,
            "expected more than half of tumor cells to have different O2 at distinct λ: got {n_diff} differing of {n_tumor} tumor"
        );
    }

    /// **Acceptance criterion #2** ("O₂ cycling works in 3D"):
    /// demonstrates the consumer-side square-wave cycling pattern that
    /// `sim-tme`'s 2D loop uses (`cycling_lambda(step, period, λ_low,
    /// λ_high)`). Runs a short step loop, alternating `λ` between
    /// normoxic and hypoxic per the square-wave schedule, and asserts:
    /// - the field at a normoxic-phase step differs from the field at
    ///   a hypoxic-phase step (cycling actually changes the field)
    /// - the field at two normoxic-phase steps is identical (consumer
    ///   can rely on the function's determinism in λ — no hidden state)
    ///
    /// This is the canonical pattern downstream consumers (#188, #191,
    /// #197) should follow. The library function itself is dim-agnostic;
    /// "cycling" is wholly a caller responsibility.
    #[test]
    fn o2_cycling_pattern_demonstrates_consumer_loop() {
        let g = TumorGrid3D::generate(20, 20, 20, 20.0, 42);
        let lambda_high = 150.0_f64; // normoxic phase
        let lambda_low = 50.0_f64; // hypoxic phase
        let period = 4u32; // first half normoxic, second half hypoxic

        // Mirrors sim-tme's `cycling_lambda` helper (binary-local, not
        // imported here; same square-wave O2-cycling definition).
        let cycling_lambda = |step: u32| {
            if (step % period) < period / 2 {
                lambda_high
            } else {
                lambda_low
            }
        };

        let mut fields = Vec::with_capacity(period as usize);
        for step in 0..period {
            fields.push(radial_o2_field(&g, cycling_lambda(step)));
        }

        // step 0 = normoxic, step 2 = hypoxic — must differ.
        assert!(
            fields[0] != fields[2],
            "cycling: normoxic-phase field must differ from hypoxic-phase field"
        );
        // step 0 and step 1 both fall in the normoxic half — identical.
        assert_eq!(
            fields[0], fields[1],
            "cycling: two normoxic-phase steps should produce identical fields (function is deterministic in λ)"
        );
        // step 2 and step 3 both fall in the hypoxic half — identical.
        assert_eq!(
            fields[2], fields[3],
            "cycling: two hypoxic-phase steps should produce identical fields"
        );
    }

    /// All-dead tumor cells in a hand-built grid → rates are 1.0 for
    /// populated zones, 0.0 for empty zones. Stromal exclusion verified.
    /// Hand-crafted by toggling `state.dead` directly.
    #[test]
    fn zone_kill_rates_all_dead_gives_one() {
        let mut g = TumorGrid3D::generate(20, 20, 20, 20.0, 42);
        for gc in g.cells.iter_mut() {
            if gc.is_tumor {
                gc.state.dead = true;
            }
        }
        // shell_depth_um = 60µm → normoxic = [0, 60), transition = [60, 180),
        // hypoxic = [180, ∞). For tumor_radius_um = 180, deepest cell has
        // depth = 180 → goes into hypoxic. Most cells will be in the
        // 60-180 range.
        let (norm, trans, hyp) = radial_o2_zone_kill_rates(&g, 60.0);
        // Each populated zone is all-dead → rate = 1.0.
        // Possible an empty zone returns 0.0 — guard against that with .max.
        // For 20³ at shell=60, all three zones SHOULD be populated.
        assert_eq!(norm, 1.0, "normoxic should be 1.0 (all dead)");
        assert_eq!(trans, 1.0, "transition should be 1.0 (all dead)");
        assert_eq!(hyp, 1.0, "hypoxic should be 1.0 (all dead)");
    }

    /// Fresh grid (no dead cells) → all zone rates = 0.0. Empty-zone path
    /// returns 0.0 instead of NaN. Stromal-only grid would also return
    /// (0,0,0), but here we use a normal grid that does have tumor cells.
    #[test]
    fn zone_kill_rates_no_dead_gives_zero() {
        let g = TumorGrid3D::generate(20, 20, 20, 20.0, 42);
        let (norm, trans, hyp) = radial_o2_zone_kill_rates(&g, 60.0);
        assert_eq!((norm, trans, hyp), (0.0, 0.0, 0.0));
    }

    /// **Acceptance criterion #4** (issue #187): "Comparison: 2D vs 3D
    /// hypoxic volume fractions at matched λ."
    ///
    /// Construction: 40×40 (2D) and 40×40×40 (3D) at cell_size=20 µm →
    /// both have `tumor_radius_um = 40 × 0.45 × 20 = 360 µm`. With
    /// `shell_depth_um = 100 µm`, the hypoxic threshold is 300 µm; the
    /// deepest cell (center) is at depth 360 µm — so both geometries have
    /// non-empty hypoxic regions.
    ///
    /// **Finding — issue text was inverted.** The issue context claimed
    /// "the hypoxic VOLUME fraction is larger than the 2D hypoxic AREA
    /// fraction (cubic vs quadratic scaling)." Pure geometry says the
    /// opposite: for matched R and λ, the hypoxic fraction is
    /// `((R - 3λ) / R)^d` where `d = 2` (area) or `d = 3` (volume). Since
    /// `(R-3λ)/R ∈ (0, 1)`, raising to a higher power gives a *smaller*
    /// number — so **3D hypoxic fraction < 2D hypoxic fraction**.
    ///
    /// Where cubic-vs-quadratic scaling DOES dominate is the *normoxic
    /// shell*: 3D shell volume scales as the surface area (~R²) times
    /// shell thickness, while 2D shell area scales as the perimeter
    /// (~R) times shell thickness — so 3D has *more* near-surface cells
    /// relative to total, *less* deep-core cells. This test asserts both
    /// halves of that relationship.
    ///
    /// The Vaupel 1989 observation that "3D spheroids are more hypoxic
    /// than 2D cultures" reflects a *biological* effect (poor 3D
    /// vasculature reduces effective O₂ supply) that this pure-geometry
    /// model does not capture — a real-world spheroid's effective λ is
    /// smaller than a 2D culture's, which is a parameter choice, not a
    /// geometry consequence.
    #[test]
    fn matched_lambda_2d_vs_3d_zone_fractions() {
        let cell_size_um = 20.0;
        let shell_depth_um = 100.0;
        let hypoxic_threshold_um = 3.0 * shell_depth_um;

        // 2D: 40×40 grid. Uses the lifted `TumorGrid::radial_depth_um`
        // (#224 item 1a) — same depth-from-edge math sim-tme's
        // `compute_depth_map` calls, sharing `TUMOR_RADIUS_FRACTION`,
        // so this test can't drift from the binary.
        let g2 = TumorGrid::generate(40, 40, cell_size_um, 42);
        let (mut norm_2, mut hyp_2, mut total_2) = (0usize, 0usize, 0usize);
        for r in 0..g2.rows {
            for c in 0..g2.cols {
                let gc = g2.get(r, c);
                if !gc.is_tumor {
                    continue;
                }
                let depth_um = g2.radial_depth_um(r, c).max(0.0);
                total_2 += 1;
                if depth_um < shell_depth_um {
                    norm_2 += 1;
                }
                if depth_um >= hypoxic_threshold_um {
                    hyp_2 += 1;
                }
            }
        }
        let frac_norm_2d = norm_2 as f64 / total_2 as f64;
        let frac_hyp_2d = hyp_2 as f64 / total_2 as f64;

        // 3D: 40³ grid using my new infrastructure.
        let g3 = TumorGrid3D::generate(40, 40, 40, cell_size_um, 42);
        let (mut norm_3, mut hyp_3, mut total_3) = (0usize, 0usize, 0usize);
        for r in 0..g3.rows {
            for c in 0..g3.cols {
                for l in 0..g3.layers {
                    let gc = g3.get(r, c, l);
                    if !gc.is_tumor {
                        continue;
                    }
                    let depth_um = g3.radial_depth_um(r, c, l);
                    total_3 += 1;
                    if depth_um < shell_depth_um {
                        norm_3 += 1;
                    }
                    if depth_um >= hypoxic_threshold_um {
                        hyp_3 += 1;
                    }
                }
            }
        }
        let frac_norm_3d = norm_3 as f64 / total_3 as f64;
        let frac_hyp_3d = hyp_3 as f64 / total_3 as f64;

        assert!(
            hyp_2 > 0 && hyp_3 > 0 && norm_2 > 0 && norm_3 > 0,
            "test precondition: expected all four counts non-empty; got 2D(norm={norm_2}, hyp={hyp_2}), 3D(norm={norm_3}, hyp={hyp_3})"
        );

        // Cubic-vs-quadratic scaling acts on the *normoxic shell* (more
        // surface relative to total in 3D), not the hypoxic core.
        assert!(
            frac_norm_3d > frac_norm_2d,
            "3D normoxic shell fraction ({frac_norm_3d:.4}) should exceed 2D \
             ({frac_norm_2d:.4}) — surface area ~ R² (3D) vs perimeter ~ R (2D)"
        );

        // Corollary: 3D hypoxic core fraction is *smaller* than 2D at matched R, λ.
        // ((R-3λ)/R)^3 < ((R-3λ)/R)^2 because the ratio is in (0, 1).
        assert!(
            frac_hyp_3d < frac_hyp_2d,
            "3D hypoxic core fraction ({frac_hyp_3d:.4}) should be SMALLER than 2D \
             ({frac_hyp_2d:.4}) — pure-geometry consequence; biological effects \
             (poor 3D vasculature) are not modelled here (see test docstring)"
        );
    }
}
