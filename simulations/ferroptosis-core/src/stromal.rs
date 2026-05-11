//! 3D stromal-shielding boundary detection for spheroid tumors.
//!
//! Cancer-associated fibroblasts (CAFs) in the peritumoral stroma supply
//! cysteine (via GGT1-mediated GSH cleavage) and oleic acid (via
//! ACSL3-mediated uptake) to adjacent tumor cells, boosting their GSH
//! antioxidant capacity and MUFA membrane protection respectively. The
//! effect is **single-cell-deep** — only tumor cells with at least one
//! stromal Moore neighbor are CAF-shielded.
//!
//! In 2D, sim-tme uses 8-Moore neighbors (`sim-tme/main.rs:431-452`); the
//! shielded shell is an annular boundary one cell thick. In 3D, we use
//! 26-Moore neighbors via [`TumorGrid3D::neighbors`]; the shielded shell
//! is a spherical surface one cell thick.
//!
//! **Surface-to-volume scaling.** The boundary fraction is roughly
//! `2 × t/R` in 2D (perimeter/area ~ 2/R for thickness t=1) and `3 × t/R`
//! in 3D (surface/volume ~ 3/R) at matched tumor radius `R`. So at
//! matched dimensions, 3D shielding affects ~1.5× more cells relative
//! to total than 2D — which is the biological point of issue #189.
//! Larger boundary fraction means stromal shielding has a *bigger*
//! impact in 3D.
//!
//! **API design — pure functions, no mutation.** Same posture as
//! [`crate::oxygen`] and [`crate::ph`]: the library identifies which
//! cells are shielded and reports a kill-rate metric; per-cell
//! GSH/MUFA boost application stays consumer-side (sim-tme's 2D code
//! applies boosts inline in `run_spatial_with_immune`, not as a library
//! helper). Consumers get the mask once and reuse it for boost
//! application, kill-rate reporting, and visualization.
//!
//! **Canonical boost magnitudes** live in sim-tme's `StromalConfig`
//! (`simulations/sim-tme/src/main.rs:406`): `gsh_boost_per_step = 0.06`,
//! `gsh_boost_cap = 18.0`, `mufa_boost_per_step = 0.003`, `mufa_boost_cap
//! = 0.25`. A 3D consumer (#195 sim-tme-3d, #197 cell-level biochem)
//! should either re-export these values from sim-tme or lift `StromalConfig`
//! into the library — both are explicit follow-up work, not in this PR's
//! scope.
//!
//! Refs: PMID 34373744 (CAF metabolic reprogramming),
//! PMID 31813804 (ACSL3-mediated oleic acid),
//! PMID 30842648 (MUFA ferroptosis).
//!
//! ## Quick example
//!
//! ```
//! use ferroptosis_core::grid::TumorGrid3D;
//! use ferroptosis_core::stromal::{stromal_adjacency_mask, stromal_adjacent_kill_rate};
//!
//! // 10³ grid is enough to demonstrate the API — 1k cells exercises
//! // the same code paths as larger sims at 1/64th the doctest cost.
//! let g = TumorGrid3D::generate(10, 10, 10, 20.0, 42);
//! let mask = stromal_adjacency_mask(&g);
//! assert_eq!(mask.len(), g.cells.len());
//!
//! // Fresh grid → no dead cells → kill rate is 0.0 even on shielded cells.
//! let rate = stromal_adjacent_kill_rate(&g, &mask);
//! assert_eq!(rate, 0.0);
//! ```

use crate::grid::TumorGrid3D;

/// Boolean mask identifying CAF-shielded boundary tumor cells.
///
/// Returns a `Vec<bool>` of length `grid.cells.len()` in the same flat
/// order (`r·cols·layers + c·layers + l`). An entry is `true` iff the
/// corresponding cell is a tumor cell (`is_tumor == true`) with at least
/// one stromal (`is_tumor == false`) **26-Moore neighbor** (the 3D analog
/// of sim-tme's 2D 8-Moore boundary detection). Stromal cells themselves
/// are always `false`.
///
/// **Geometry note.** With `TumorGrid3D::generate`'s default
/// `tumor_radius = 0.45 × min_dim`, the tumor never touches the grid
/// edges, so this function correctly identifies the spheroid surface.
/// For user-constructed grids where the tumor *does* touch the bounding
/// box, the **failure mode is a false negative, not a false positive**:
/// a tumor cell at a grid face/edge/corner has only its in-grid
/// neighbors checked (missing out-of-grid positions don't count as
/// "stromal"), so a grid-edge tumor cell whose in-grid neighbors are
/// all tumor will NOT be flagged, even though it geometrically sits at
/// the spheroid surface. This is the same semantics as sim-tme's 2D
/// implementation. Mitigation: keep `tumor_radius < min_dim / 2` so the
/// spheroid stays interior to the bounding box.
///
/// **Cost.** O(N × 26) for N = `grid.cells.len()`. Each cell costs a
/// 26-neighbor sweep (stack-allocated via [`TumorGrid3D::neighbors`]).
/// Negligible for hundreds of thousands of cells.
///
/// **Compute once and reuse** (perf hint for #194): the mask depends
/// only on `is_tumor` topology, which is **invariant under standard
/// simulation flow** — cells die but no cell ever changes from tumor
/// to stromal or vice versa. So this function should be called *once*
/// after `TumorGrid3D::generate` and the mask cached for the whole
/// simulation, NOT recomputed per timestep. At 200³ the per-call cost
/// is ~200M index ops; called per step over 180 steps that's 36B
/// wasted ops. `#[must_use]` enforces the return-value-is-meaningful
/// contract.
#[must_use = "the boundary mask must be observed; calling for side-effects suggests a logic bug — compute once after `generate` and cache for the simulation"]
pub fn stromal_adjacency_mask(grid: &TumorGrid3D) -> Vec<bool> {
    let n = grid.cells.len();
    let mut mask = vec![false; n];
    for r in 0..grid.rows {
        for c in 0..grid.cols {
            for l in 0..grid.layers {
                // Use `grid.get` for reads and `grid.flat_index` for the
                // mask-write index — both delegate the row-major layout
                // formula to TumorGrid3D so a future stride change can't
                // silently desync this module from `grid.cells`.
                if !grid.get(r, c, l).is_tumor {
                    continue;
                }
                let idx = grid.flat_index(r, c, l);
                let (neighbors, count) = grid.neighbors(r, c, l);
                for &(nr, nc, nl) in &neighbors[..count] {
                    if !grid.get(nr, nc, nl).is_tumor {
                        mask[idx] = true;
                        break;
                    }
                }
            }
        }
    }
    mask
}

/// Dead-cell rate among CAF-shielded boundary tumor cells.
///
/// Returns `dead_in_mask / total_in_mask`, or `0.0` if the mask flags no
/// tumor cells (no division-by-zero panic). **Stromal cells are excluded
/// from both numerator and denominator** even if `mask[stromal_idx]` were
/// somehow `true` (the function checks `is_tumor` per cell, so the mask
/// + is_tumor double-gate ensures only tumor cells contribute).
///
/// **Mask contract.** Must be in the same flat-index order as
/// `grid.cells`, with `mask.len() == grid.cells.len()`. Typically
/// produced by [`stromal_adjacency_mask`]. A length mismatch is a
/// programming error and panics in **both debug and release** via
/// `assert!` (one usize comparison per call — negligible perf cost,
/// clearer error message than the out-of-bounds index panic that would
/// otherwise occur on the first iteration).
#[must_use = "the kill rate is the function's only output; ignoring it suggests a logic bug"]
pub fn stromal_adjacent_kill_rate(grid: &TumorGrid3D, mask: &[bool]) -> f64 {
    assert!(
        mask.len() == grid.cells.len(),
        "stromal_adjacent_kill_rate: mask length {} must equal grid.cells.len() {}",
        mask.len(),
        grid.cells.len()
    );

    let mut total = 0usize;
    let mut dead = 0usize;
    for (idx, gc) in grid.cells.iter().enumerate() {
        if gc.is_tumor && mask[idx] {
            total += 1;
            if gc.state.dead {
                dead += 1;
            }
        }
    }
    if total > 0 {
        dead as f64 / total as f64
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::TumorGrid;

    /// Mask length always matches grid cell count, regardless of grid
    /// dimensions or seed.
    #[test]
    fn mask_length_matches_grid() {
        let g = TumorGrid3D::generate(7, 5, 11, 20.0, 42);
        let mask = stromal_adjacency_mask(&g);
        assert_eq!(mask.len(), g.cells.len());
        assert_eq!(mask.len(), 7 * 5 * 11);
    }

    /// Stromal cells (is_tumor = false) are never flagged in the mask,
    /// regardless of their neighbors. Walks every cell and asserts.
    #[test]
    fn stromal_cells_never_in_mask() {
        let g = TumorGrid3D::generate(10, 10, 10, 20.0, 42);
        let mask = stromal_adjacency_mask(&g);

        let mut stromal_count = 0usize;
        for (i, gc) in g.cells.iter().enumerate() {
            if !gc.is_tumor {
                assert!(
                    !mask[i],
                    "stromal cell at flat idx {i} should never be in mask"
                );
                stromal_count += 1;
            }
        }
        assert!(
            stromal_count > 0,
            "expected some stromal cells in a 10³ grid"
        );
    }

    /// A deep-interior tumor cell — far enough from the spheroid surface
    /// that all 26 Moore neighbors are also tumor cells — has `mask[i] = false`.
    ///
    /// For a 20³ grid: center=(10,10,10), tumor_radius_lattice = 9.0. The
    /// center cell (10,10,10) is 0 lattice from center, depth = 9 lattice.
    /// All 26 of its (10±1, 10±1, 10±1) neighbors are at most √3 ≈ 1.73
    /// lattice from center (well inside the 9.0-lattice sphere), so all
    /// are tumor. Mask should be false.
    #[test]
    fn deep_interior_tumor_cell_not_in_mask() {
        let g = TumorGrid3D::generate(20, 20, 20, 20.0, 42);
        let mask = stromal_adjacency_mask(&g);
        let flat_center = 10 * g.cols * g.layers + 10 * g.layers + 10;
        assert!(
            g.get(10, 10, 10).is_tumor,
            "test precondition: center is tumor"
        );
        assert!(
            !mask[flat_center],
            "deep-interior cell at flat idx {flat_center} should not be in the boundary mask"
        );
    }

    /// A surface tumor cell — within ~1 lattice unit of the spheroid edge —
    /// has at least one stromal neighbor and is flagged in the mask.
    ///
    /// For a 20³ grid: tumor_radius_lattice = 9.0. Cell (10, 10, 19) is
    /// 9 lattice from center (sits exactly on the spheroid surface per
    /// the grid::tests_3d invariant). Some of its 26 neighbors at
    /// (10±1, 10±1, 18..20) are at distance > 9 from center → stromal.
    #[test]
    fn surface_tumor_cell_in_mask() {
        let g = TumorGrid3D::generate(20, 20, 20, 20.0, 42);
        let mask = stromal_adjacency_mask(&g);
        let flat_surface = 10 * g.cols * g.layers + 10 * g.layers + 19;
        assert!(
            g.get(10, 10, 19).is_tumor,
            "test precondition: this is a tumor cell"
        );
        assert!(
            mask[flat_surface],
            "surface tumor cell at flat idx {flat_surface} should be in the boundary mask"
        );
    }

    /// Mask construction is deterministic — same seed/dims → identical
    /// mask. No RNG, just iteration order.
    #[test]
    fn mask_is_deterministic() {
        let g1 = TumorGrid3D::generate(15, 15, 15, 20.0, 42);
        let g2 = TumorGrid3D::generate(15, 15, 15, 20.0, 42);
        let mask1 = stromal_adjacency_mask(&g1);
        let mask2 = stromal_adjacency_mask(&g2);
        assert_eq!(mask1, mask2);
    }

    /// When all masked tumor cells are dead, the kill rate is exactly 1.0
    /// (every cell in the boundary shell counts as dead).
    #[test]
    fn kill_rate_all_dead_returns_one() {
        let mut g = TumorGrid3D::generate(20, 20, 20, 20.0, 42);
        let mask = stromal_adjacency_mask(&g);
        for (i, gc) in g.cells.iter_mut().enumerate() {
            if gc.is_tumor && mask[i] {
                gc.state.dead = true;
            }
        }
        let rate = stromal_adjacent_kill_rate(&g, &mask);
        assert_eq!(rate, 1.0);
    }

    /// Fresh grid (no dead cells) → kill rate is 0.0.
    #[test]
    fn kill_rate_no_dead_returns_zero() {
        let g = TumorGrid3D::generate(20, 20, 20, 20.0, 42);
        let mask = stromal_adjacency_mask(&g);
        let rate = stromal_adjacent_kill_rate(&g, &mask);
        assert_eq!(rate, 0.0);
    }

    /// All-false mask → kill rate is 0.0 (no division by zero).
    #[test]
    fn kill_rate_empty_mask_returns_zero() {
        let g = TumorGrid3D::generate(20, 20, 20, 20.0, 42);
        let mask = vec![false; g.cells.len()];
        let rate = stromal_adjacent_kill_rate(&g, &mask);
        assert_eq!(rate, 0.0);
    }

    /// **v2 addition**: mismatched mask length panics in **both debug and
    /// release** via a regular `assert!` (one usize comparison per call;
    /// negligible perf cost). Catches the programming-contract violation
    /// with a clear message rather than the cryptic out-of-bounds index
    /// panic that would otherwise hit on the first iteration.
    #[test]
    #[should_panic(expected = "mask length")]
    fn kill_rate_validates_mask_length() {
        let g = TumorGrid3D::generate(10, 10, 10, 20.0, 42);
        // Wrong length: half the expected size.
        let bad_mask = vec![false; g.cells.len() / 2];
        let _ = stromal_adjacent_kill_rate(&g, &bad_mask);
    }

    /// **v2 addition**: the `mask[i] && is_tumor` double-gate ensures
    /// stromal cells never contribute to the rate, even if a caller
    /// passes a hand-built mask that incorrectly flags stromal cells.
    /// Locks down assumption I from the planning rigor pass.
    #[test]
    fn kill_rate_excludes_stromal_cells_from_count() {
        let mut g = TumorGrid3D::generate(20, 20, 20, 20.0, 42);
        // Find a stromal cell and mark it dead.
        let stromal_idx = g
            .cells
            .iter()
            .position(|gc| !gc.is_tumor)
            .expect("expected some stromal cells in a 20³ grid");
        g.cells[stromal_idx].state.dead = true;

        // Hand-build a mask that flags ALL cells including the stromal one.
        let bad_mask = vec![true; g.cells.len()];
        let rate = stromal_adjacent_kill_rate(&g, &bad_mask);

        // The stromal cell's death must NOT contribute. The rate is
        // (tumor dead) / (tumor count) — both zero on a fresh grid.
        // Even with the stromal cell dead, rate stays 0.0 because the
        // function filters by is_tumor.
        assert_eq!(rate, 0.0, "stromal cell should not contribute to kill rate");
    }

    /// **Acceptance criterion #3** (issue #189): "Comparison: 2D vs 3D
    /// boundary fraction and stromal kill rates."
    ///
    /// At matched tumor radius R, geometry predicts 3D boundary fraction
    /// ≈ 3/R (surface-to-volume) and 2D boundary fraction ≈ 2/R
    /// (perimeter-to-area). So 3D boundary fraction is ~1.5× the 2D
    /// fraction at matched R.
    ///
    /// **Why not assert a specific ratio**: issue text quotes "3% vs 1.1%"
    /// but those numbers were measured at DIFFERENT tumor radii (2D was
    /// from sim-tme's 400×400 default with R=180; 3D was hypothetical
    /// "100-cell-radius sphere" with R=100). At matched R, the ratio is
    /// closer to 1.5×, not 2.7×. So we assert only the directional
    /// invariant: 3D fraction > 2D fraction.
    ///
    /// **Source of truth — keep in sync with `sim-tme/src/main.rs:433`
    /// (`stromal_adjacency_mask`).** This test inlines the 2D mask
    /// construction; same DRY caveat as the O₂/pH cross-geometry tests.
    /// No automated guard against drift.
    #[test]
    fn matched_dims_3d_boundary_fraction_exceeds_2d() {
        let cell_size_um = 20.0;

        // 2D: 40×40 grid → R = 18 lattice units. Replicate sim-tme's
        // 8-Moore boundary detection inline.
        let g2 = TumorGrid::generate(40, 40, cell_size_um, 42);
        let mut mask_2d = vec![false; g2.rows * g2.cols];
        for r in 0..g2.rows {
            for c in 0..g2.cols {
                let idx = r * g2.cols + c;
                if !g2.cells[idx].is_tumor {
                    continue;
                }
                let (neighbors, count) = g2.neighbors(r, c);
                for &(nr, nc) in &neighbors[..count] {
                    if !g2.cells[nr * g2.cols + nc].is_tumor {
                        mask_2d[idx] = true;
                        break;
                    }
                }
            }
        }
        let total_tumor_2d = g2.cells.iter().filter(|gc| gc.is_tumor).count();
        let boundary_2d = mask_2d.iter().filter(|&&b| b).count();
        let frac_2d = boundary_2d as f64 / total_tumor_2d as f64;

        // 3D: 40³ grid → R = 18 lattice units (matched).
        let g3 = TumorGrid3D::generate(40, 40, 40, cell_size_um, 42);
        let mask_3d = stromal_adjacency_mask(&g3);
        let total_tumor_3d = g3.cells.iter().filter(|gc| gc.is_tumor).count();
        let boundary_3d = mask_3d.iter().filter(|&&b| b).count();
        let frac_3d = boundary_3d as f64 / total_tumor_3d as f64;

        assert!(
            boundary_2d > 0 && boundary_3d > 0,
            "test precondition: expected non-empty boundary regions; got 2D={boundary_2d}, 3D={boundary_3d}"
        );

        // Theoretical ratio: 3D ~ 3/R, 2D ~ 2/R → 1.5× at matched R.
        // Empirically measured at this grid (40³, seed=42): 1.66× (lattice
        // quantization gives slightly more than continuum theory). Asserting
        // **`1.3 < ratio < 1.8`** (bounded interval) is stronger than a
        // directional `frac_3d > frac_2d` would be:
        // - lower bound catches a regression that cuts neighbor coverage
        //   (would give ratio ≈ 0.5–0.8× — clearly < 1.3 but possibly > 1.0)
        // - upper bound catches a regression that double-counts cells or
        //   misses the stromal filter (would inflate the 3D numerator)
        // Both bounds give comfortable margin from the empirical 1.66×.
        let ratio = frac_3d / frac_2d;
        assert!(
            ratio > 1.3 && ratio < 1.8,
            "3D/2D boundary-fraction ratio ({ratio:.3}) should be in (1.3, 1.8) — \
             surface-to-volume scaling (3/R vs 2/R) predicts ≈ 1.5×, empirical at \
             this grid is ≈ 1.66×. Outside this interval likely means the 26-neighbor \
             sweep is mis-coverage (under: < 1.3; over: > 1.8). \
             2D frac={frac_2d:.4}, 3D frac={frac_3d:.4}"
        );
    }
}
