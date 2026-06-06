//! Steady-state reaction-diffusion supply field (#343).
//!
//! The vasculature/oxygen layers ship two supply proxies: an edge-distance
//! gradient and an explicit-vessel `exp(-dist_to_nearest_vessel / λ)` field
//! ([`crate::vasculature::vessel_supply_field`]). Both are *monotonic in
//! distance to the nearest source*. Real tumor O2/drug fields are not: with an
//! irregular vessel network, diffusion superposes the contributions of *all*
//! nearby vessels while tissue consumption draws the field back down, so the
//! steady state has non-monotonic pockets (a point can be far from any single
//! vessel yet well-supplied because it sits between several, or close to one
//! vessel yet starved because consumption around it is high). An exponential
//! nearest-vessel proxy averages these away.
//!
//! This module solves the linear steady-state reaction-diffusion equation for a
//! normalized supply concentration `c ∈ [0, 1]` over the same vessel network the
//! exponential proxy uses, so the two can be compared on identical geometry:
//!
//! ```text
//!   D ∇²c − k·c = 0     in consuming (tumor) tissue
//!   D ∇²c       = 0     in non-consuming (stromal) tissue
//!   c = 1               at vessel voxels        (Dirichlet sources)
//!   ∂c/∂n = 0           at the domain boundary  (no-flux / reflective)
//! ```
//!
//! The single physical length is the **diffusion-consumption length**
//! `λ = sqrt(D / k)`. Setting `λ` equal to the exponential proxy's `lambda_um`
//! makes the two apples-to-apples: a *single isolated* planar source reproduces
//! `exp(-dist / λ)`, so any divergence between the RD field and the proxy is
//! exactly the multi-vessel / consumption effect the proxy cannot represent.
//!
//! **Discretization.** Uniform 6-point (von Neumann) stencil with spacing
//! `h = cell_size_um`. The interior balance `D ∇²c − k c = 0` discretizes to
//!
//! ```text
//!   c[i] = (Σ_{nb ∈ grid} c[nb]) / (N + γ),   γ = k·h²/D = (h/λ)²
//! ```
//!
//! where `N` is the count of in-grid neighbors (stromal voxels use `γ = 0`).
//! Excluding out-of-grid neighbors and reducing `N` is the finite-volume no-flux
//! boundary (a zero-flux outer face contributes nothing to the cell balance),
//! the same reflective convention the slab/oxygen layers use. Solved by
//! Gauss-Seidel with successive over-relaxation (SOR) in fixed linear-index
//! sweep order, so the result is fully deterministic.
//!
//! **Validation.** [`analytical_1d_slab`] is the closed-form 1-D solution
//! (Thomlinson & Gray 1955; the Krogh/diffusion-limited slab), used by the unit
//! tests as the published-model benchmark that the solver must reproduce.
//!
//! This is an *opt-in alternative* supply field: nothing in the production
//! matrix calls it, so adding the module is byte-identical.

use crate::grid::TumorGrid3D;

/// Configuration for the steady-state reaction-diffusion solver.
#[derive(Debug, Clone, PartialEq)]
pub struct ReactionDiffusionConfig {
    /// Diffusion-consumption length `λ = sqrt(D/k)` in µm. Set equal to the
    /// exponential proxy's `lambda_um` for an apples-to-apples comparison.
    pub diffusion_length_um: f64,
    /// Maximum SOR sweeps before giving up (the solution is returned either
    /// way; [`RdSolution::converged`] reports whether `tol` was met).
    pub max_iters: usize,
    /// Convergence tolerance on the max per-sweep change `max|c_new − c_old|`.
    pub tol: f64,
    /// SOR over-relaxation factor in `(0, 2)`. `1.0` is plain Gauss-Seidel;
    /// `> 1` accelerates convergence on the Laplace-like operator. `1.8` is a
    /// good general default for these grids.
    pub omega: f64,
}

impl ReactionDiffusionConfig {
    /// Default solver tuned to a given diffusion length: SOR `ω = 1.8`,
    /// `max_iters = 5000`, `tol = 1e-6`.
    pub fn new(diffusion_length_um: f64) -> Self {
        debug_assert!(
            diffusion_length_um.is_finite() && diffusion_length_um > 0.0,
            "diffusion_length_um must be finite and positive, got {diffusion_length_um}"
        );
        Self {
            diffusion_length_um,
            max_iters: 5000,
            tol: 1e-6,
            omega: 1.8,
        }
    }
}

/// Result of a reaction-diffusion solve: the raw field over *every* voxel plus
/// convergence diagnostics. (The drop-in supply field that overrides non-tumor
/// voxels to `1.0` is [`reaction_diffusion_supply_field`].)
#[derive(Debug, Clone)]
pub struct RdSolution {
    /// Solved concentration `c ∈ [0, 1]` for every voxel, in grid index order.
    pub field: Vec<f64>,
    /// Number of SOR sweeps actually performed.
    pub iters: usize,
    /// Final max per-sweep change `max|c_new − c_old|`.
    pub residual: f64,
    /// Whether `residual < cfg.tol` was reached within `cfg.max_iters`.
    pub converged: bool,
}

/// Map a continuous vessel coordinate to its nearest voxel and flag it as a
/// Dirichlet source. Coordinates rounding outside the grid are clamped to the
/// boundary (a vessel just past the edge still perfuses the edge voxel).
fn vessel_mask(grid: &TumorGrid3D, vessels: &[(f64, f64, f64)]) -> Vec<bool> {
    let mut mask = vec![false; grid.cells.len()];
    let clamp = |v: f64, hi: usize| -> usize {
        if v < 0.0 {
            0
        } else {
            (v.round() as usize).min(hi.saturating_sub(1))
        }
    };
    for &(r, c, l) in vessels {
        let ri = clamp(r, grid.rows);
        let ci = clamp(c, grid.cols);
        let li = clamp(l, grid.layers);
        mask[grid.flat_index(ri, ci, li)] = true;
    }
    mask
}

/// Solve the steady-state reaction-diffusion field over `grid` with `vessels`
/// as Dirichlet sources. Returns the raw field everywhere plus diagnostics.
///
/// Panics if `vessels` is empty (with no source the field relaxes to zero,
/// which is never the intent), matching [`crate::vasculature::vessel_supply_field`].
pub fn reaction_diffusion_solve(
    grid: &TumorGrid3D,
    vessels: &[(f64, f64, f64)],
    cfg: &ReactionDiffusionConfig,
) -> RdSolution {
    assert!(
        !vessels.is_empty(),
        "reaction_diffusion_solve needs ≥1 vessel source"
    );
    debug_assert!(
        cfg.diffusion_length_um.is_finite() && cfg.diffusion_length_um > 0.0,
        "diffusion_length_um must be finite and positive"
    );
    debug_assert!(
        cfg.omega > 0.0 && cfg.omega < 2.0,
        "SOR omega must be in (0, 2), got {}",
        cfg.omega
    );

    let (rows, cols, layers) = (grid.rows, grid.cols, grid.layers);
    let n = grid.cells.len();
    let h = grid.cell_size_um;
    // γ = k·h²/D = (h/λ)². Tumor voxels consume; stroma do not.
    let gamma = (h / cfg.diffusion_length_um).powi(2);

    let source = vessel_mask(grid, vessels);

    // Initialize: sources at 1.0, everything else at 0.0.
    let mut field = vec![0.0_f64; n];
    for (idx, &is_src) in source.iter().enumerate() {
        if is_src {
            field[idx] = 1.0;
        }
    }

    let mut iters = 0;
    let mut residual = f64::INFINITY;
    let mut converged = false;

    // SOR: Gauss-Seidel sweep in linear index order with over-relaxation.
    // Reads updated neighbor values in-place (deterministic given fixed order).
    while iters < cfg.max_iters {
        iters += 1;
        let mut max_change = 0.0_f64;
        for idx in 0..n {
            if source[idx] {
                continue; // Dirichlet, fixed at 1.0
            }
            let (r, c, l) = grid.coords(idx);
            let mut sum = 0.0_f64;
            let mut count = 0.0_f64;
            // 6-neighbor von Neumann stencil; out-of-grid neighbors are
            // dropped (finite-volume no-flux boundary).
            if r > 0 {
                sum += field[grid.flat_index(r - 1, c, l)];
                count += 1.0;
            }
            if r + 1 < rows {
                sum += field[grid.flat_index(r + 1, c, l)];
                count += 1.0;
            }
            if c > 0 {
                sum += field[grid.flat_index(r, c - 1, l)];
                count += 1.0;
            }
            if c + 1 < cols {
                sum += field[grid.flat_index(r, c + 1, l)];
                count += 1.0;
            }
            if l > 0 {
                sum += field[grid.flat_index(r, c, l - 1)];
                count += 1.0;
            }
            if l + 1 < layers {
                sum += field[grid.flat_index(r, c, l + 1)];
                count += 1.0;
            }
            // Consumption only in tumor voxels.
            let g = if grid.cells[idx].is_tumor { gamma } else { 0.0 };
            let gs = sum / (count + g);
            let old = field[idx];
            let new = ((1.0 - cfg.omega) * old + cfg.omega * gs).clamp(0.0, 1.0);
            field[idx] = new;
            let change = (new - old).abs();
            if change > max_change {
                max_change = change;
            }
        }
        residual = max_change;
        if residual < cfg.tol {
            converged = true;
            break;
        }
    }

    RdSolution {
        field,
        iters,
        residual,
        converged,
    }
}

/// Drop-in alternative to [`crate::vasculature::vessel_supply_field`]: the
/// steady-state reaction-diffusion supply factor for every voxel, with
/// non-tumor voxels overridden to `1.0` (matching the proxy's contract that
/// only tumor cells carry a sub-unity supply). Use [`reaction_diffusion_solve`]
/// when you need the raw stromal field or the convergence diagnostics.
pub fn reaction_diffusion_supply_field(
    grid: &TumorGrid3D,
    vessels: &[(f64, f64, f64)],
    cfg: &ReactionDiffusionConfig,
) -> Vec<f64> {
    let sol = reaction_diffusion_solve(grid, vessels, cfg);
    grid.cells
        .iter()
        .zip(sol.field)
        .map(|(cell, v)| if cell.is_tumor { v } else { 1.0 })
        .collect()
}

/// Closed-form 1-D steady-state reaction-diffusion concentration (the
/// published-model benchmark; Thomlinson & Gray 1955, the diffusion-limited
/// slab / Krogh model). Source `c = 1` at `x = 0`, no-flux at `x = width_um`,
/// uniform consumption with diffusion length `λ`:
///
/// ```text
///   c(x) = cosh((W − x) / λ) / cosh(W / λ)
/// ```
///
/// For `W ≫ λ` this reduces to `exp(-x/λ)` (the exponential proxy), confirming
/// the proxy is exactly the single-isolated-source limit of the RD field.
pub fn analytical_1d_slab(x_um: f64, width_um: f64, diffusion_length_um: f64) -> f64 {
    let lam = diffusion_length_um;
    (((width_um - x_um) / lam).cosh() / (width_um / lam).cosh()).clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vasculature::vessel_supply_field;

    /// A thin bar (cols = layers = 1) with a Dirichlet source plane at r = 0,
    /// all-tumor (uniform consumption), is the discrete analog of the 1-D slab.
    fn bar(rows: usize, h: f64) -> TumorGrid3D {
        let mut grid = TumorGrid3D::generate(rows, 1, 1, h, 7);
        for cell in grid.cells.iter_mut() {
            cell.is_tumor = true;
        }
        grid
    }

    #[test]
    fn analytical_1d_slab_reduces_to_exponential_for_wide_slab() {
        // W ≫ λ: cosh form collapses to exp(-x/λ).
        let lam = 10.0;
        let w = 2000.0;
        for &x in &[0.0, 5.0, 10.0, 20.0, 50.0] {
            let cosh = analytical_1d_slab(x, w, lam);
            let exp = (-x / lam).exp();
            assert!((cosh - exp).abs() < 1e-3, "x={x}: cosh={cosh} vs exp={exp}");
        }
    }

    #[test]
    fn solver_reproduces_1d_analytical_slab() {
        // Published-model benchmark: the converged solver must match the
        // closed-form cosh profile. Finite-volume no-flux puts the reflective
        // face half a cell beyond the last node, so W = (rows-1)·h + 0.5·h.
        let rows = 40;
        let h = 10.0;
        let lam = 80.0; // 8 cells
        let grid = bar(rows, h);
        let vessels: Vec<(f64, f64, f64)> = vec![(0.0, 0.0, 0.0)];
        let mut cfg = ReactionDiffusionConfig::new(lam);
        cfg.max_iters = 20000;
        cfg.tol = 1e-9;
        let sol = reaction_diffusion_solve(&grid, &vessels, &cfg);
        assert!(sol.converged, "solver did not converge: {sol:?}");
        let w = (rows - 1) as f64 * h + 0.5 * h;
        let mut max_err = 0.0_f64;
        for r in 0..rows {
            let x = r as f64 * h;
            let got = sol.field[grid.flat_index(r, 0, 0)];
            let want = analytical_1d_slab(x, w, lam);
            max_err = max_err.max((got - want).abs());
        }
        assert!(
            max_err < 0.01,
            "solver vs analytical cosh max abs error {max_err} too large"
        );
    }

    #[test]
    fn single_isolated_source_matches_exponential_proxy() {
        // The core apples-to-apples claim: with one vessel, the RD field near
        // the source tracks the exponential proxy's e-folding length λ.
        let rows = 60;
        let h = 10.0;
        let lam = 60.0; // 6 cells
        let grid = bar(rows, h);
        let vessels: Vec<(f64, f64, f64)> = vec![(0.0, 0.0, 0.0)];
        let mut cfg = ReactionDiffusionConfig::new(lam);
        cfg.max_iters = 20000;
        cfg.tol = 1e-9;
        let sol = reaction_diffusion_solve(&grid, &vessels, &cfg);
        // e-folding: c(x+λ)/c(x) ≈ e^-1 in the interior, away from boundaries.
        let r0 = 12; // x = 2λ
        let r1 = 18; // x = 3λ
        let ratio = sol.field[grid.flat_index(r1, 0, 0)] / sol.field[grid.flat_index(r0, 0, 0)];
        assert!(
            (ratio - (-1.0_f64).exp()).abs() < 0.02,
            "e-folding ratio {ratio} not ≈ e^-1"
        );
        // Source pinned, field monotonic-decreasing, in range.
        assert_eq!(sol.field[grid.flat_index(0, 0, 0)], 1.0);
        for r in 1..rows {
            let prev = sol.field[grid.flat_index(r - 1, 0, 0)];
            let cur = sol.field[grid.flat_index(r, 0, 0)];
            assert!(cur <= prev + 1e-9 && (0.0..=1.0).contains(&cur));
        }
    }

    #[test]
    fn two_sources_enrich_the_midpoint_above_the_nearest_vessel_proxy() {
        // Non-monotonicity the proxy cannot represent: between two vessels the
        // RD field superposes both contributions, so the midpoint is better
        // supplied than the nearest-vessel exponential proxy (which sees only
        // the closer of the two) predicts.
        let rows = 41;
        let h = 10.0;
        let lam = 80.0;
        let grid = bar(rows, h);
        // Sources at both ends; the geometric midpoint (r = 20) is equidistant.
        let vessels: Vec<(f64, f64, f64)> = vec![(0.0, 0.0, 0.0), (40.0, 0.0, 0.0)];
        let mut cfg = ReactionDiffusionConfig::new(lam);
        cfg.max_iters = 20000;
        cfg.tol = 1e-9;
        let sol = reaction_diffusion_solve(&grid, &vessels, &cfg);
        let rd_mid = sol.field[grid.flat_index(20, 0, 0)];
        // Exponential proxy at the same midpoint (nearest vessel = 20 cells).
        let proxy = vessel_supply_field(&grid, &vessels, lam);
        let proxy_mid = proxy[grid.flat_index(20, 0, 0)];
        assert!(
            rd_mid > proxy_mid + 1e-6,
            "RD midpoint {rd_mid} should exceed nearest-vessel proxy {proxy_mid}"
        );
    }

    #[test]
    fn supply_field_overrides_non_tumor_to_unity() {
        // generate() leaves a stromal shell outside the tumor sphere; the
        // drop-in supply field must report 1.0 there (proxy contract).
        let grid = TumorGrid3D::generate(20, 20, 20, 15.0, 3);
        let center = (10.0, 10.0, 10.0);
        let cfg = ReactionDiffusionConfig::new(150.0);
        let field = reaction_diffusion_supply_field(&grid, &[center], &cfg);
        let mut saw_stroma = false;
        for (cell, &v) in grid.cells.iter().zip(field.iter()) {
            assert!((0.0..=1.0).contains(&v));
            if !cell.is_tumor {
                assert_eq!(v, 1.0, "non-tumor voxel must be 1.0");
                saw_stroma = true;
            }
        }
        assert!(saw_stroma, "expected a stromal shell in a 20³ grid");
    }

    #[test]
    fn solver_is_deterministic() {
        let grid = TumorGrid3D::generate(16, 16, 16, 15.0, 5);
        let vessels: Vec<(f64, f64, f64)> = vec![(4.0, 4.0, 4.0), (11.0, 11.0, 11.0)];
        let cfg = ReactionDiffusionConfig::new(120.0);
        let a = reaction_diffusion_solve(&grid, &vessels, &cfg);
        let b = reaction_diffusion_solve(&grid, &vessels, &cfg);
        assert_eq!(a.field, b.field);
        assert_eq!(a.iters, b.iters);
    }

    #[test]
    #[should_panic(expected = "needs ≥1 vessel")]
    fn empty_vessels_panics() {
        let grid = bar(8, 10.0);
        let cfg = ReactionDiffusionConfig::new(50.0);
        let _ = reaction_diffusion_solve(&grid, &[], &cfg);
    }
}
