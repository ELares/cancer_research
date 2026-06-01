//! Explicit 3D tumor vasculature (#191).
//!
//! The 2D simulations (and the 3D `oxygen::radial_o2_field`) use "distance
//! from the tumor edge" as a proxy for vasculature: the spheroid surface is
//! the only O2/drug source, so supply decays monotonically inward. Real 3D
//! tumors instead carry **internal** vessels that penetrate the volume, so
//! oxygenation is patchy: well-supplied near a vessel, hypoxic in the gaps
//! between them. Inter-vessel distance is ~100–200 µm in well-vascularized
//! tumors and ~300–500 µm in poorly-vascularized ones (Vaupel; Krogh-cylinder
//! supply geometry).
//!
//! This module implements **Option A** from #191: random vessel seed points in
//! the tumor volume, with each cell's supply set by the Krogh-style decay
//! `exp(-distance_to_nearest_vessel / λ)`. It is a drop-in alternative to
//! `oxygen::radial_o2_field` — same `Vec<f64>` per-cell-factor output, same
//! "non-tumor cells = 1.0 (well-perfused bulk)" convention — so the consumer
//! multiplies it into `cell.basal_ros` (O2) or a drug-availability field
//! identically. The same proximity factor supplies both O2 and drug.
//!
//! ## Design: independent RNG (byte-identity)
//!
//! [`place_vessels_3d`] draws vessel positions from its **own** `StdRng`, so it
//! never perturbs [`TumorGrid3D::generate`](crate::grid::TumorGrid3D::generate)'s
//! stream — the cell grid is bit-for-bit unchanged whether or not vessels are
//! placed. A consumer that doesn't opt into the vessel model keeps using
//! `radial_o2_field` and stays byte-identical.

use crate::grid::{TumorGrid3D, TUMOR_RADIUS_FRACTION};
use rand::prelude::*;

/// Vessel-network configuration. `inter_vessel_um` is the target mean spacing
/// between vessel seed points; the vessel count is derived from it and the
/// tumor volume in [`place_vessels_3d`]. The Krogh decay length λ is supplied
/// per call (it matches the condition's O2 reference λ, like `radial_o2_field`).
#[derive(Clone, Copy, Debug)]
pub struct VasculatureConfig {
    /// Target mean inter-vessel spacing (µm). Smaller ⇒ denser ⇒ better-oxygenated.
    pub inter_vessel_um: f64,
}

impl VasculatureConfig {
    /// Well-vascularized tumor (~150 µm inter-vessel spacing).
    pub fn well_vascularized() -> Self {
        VasculatureConfig {
            inter_vessel_um: 150.0,
        }
    }

    /// Poorly-vascularized tumor (~400 µm inter-vessel spacing).
    pub fn poorly_vascularized() -> Self {
        VasculatureConfig {
            inter_vessel_um: 400.0,
        }
    }
}

/// Place vessel seed points uniformly in the tumor sphere. The count is an
/// approximation from the tumor volume and target inter-vessel spacing
/// (`n ≈ tumor_volume / inter_vessel³`, assuming cubic packing, floored at 1).
/// Returns positions in **lattice (cell) coordinates**. Deterministic given
/// `(grid dims, cfg, seed)`.
///
/// Uses an **independent** `StdRng(seed)` so it never advances the RNG used by
/// [`TumorGrid3D::generate`], preserving byte-identity of the cell grid.
pub fn place_vessels_3d(
    grid: &TumorGrid3D,
    cfg: &VasculatureConfig,
    seed: u64,
) -> Vec<(f64, f64, f64)> {
    debug_assert!(
        cfg.inter_vessel_um.is_finite() && cfg.inter_vessel_um > 0.0,
        "inter_vessel_um must be finite and positive, got {}",
        cfg.inter_vessel_um
    );
    let cell_um3 = grid.cell_size_um.powi(3);
    let n_tumor = grid.cells.iter().filter(|gc| gc.is_tumor).count();
    let tumor_volume_um3 = n_tumor as f64 * cell_um3;
    let n_vessels = (tumor_volume_um3 / cfg.inter_vessel_um.powi(3))
        .round()
        .max(1.0) as usize;

    let mut rng = StdRng::seed_from_u64(seed);
    let center = (
        grid.rows as f64 / 2.0,
        grid.cols as f64 / 2.0,
        grid.layers as f64 / 2.0,
    );
    let tumor_radius = (grid.rows.min(grid.cols).min(grid.layers) as f64) * TUMOR_RADIUS_FRACTION;

    // Uniform-in-sphere sampling (cbrt radial avoids center bias — same
    // convention as generate's persister clusters / clonal seeds).
    (0..n_vessels)
        .map(|_| {
            let dist = rng.gen::<f64>().cbrt() * tumor_radius * 0.95;
            let theta = rng.gen::<f64>() * std::f64::consts::TAU;
            let cos_phi = 2.0 * rng.gen::<f64>() - 1.0;
            let sin_phi = (1.0 - cos_phi * cos_phi).sqrt();
            (
                center.0 + dist * cos_phi,
                center.1 + dist * sin_phi * theta.cos(),
                center.2 + dist * sin_phi * theta.sin(),
            )
        })
        .collect()
}

/// Per-cell supply factor from the explicit vessel network: `exp(-d/λ)` where
/// `d` is the distance (µm) to the **nearest** vessel. Drop-in replacement for
/// [`crate::oxygen::radial_o2_field`]: returns a `Vec<f64>` of length
/// `grid.cells.len()`, non-tumor cells = `1.0` (well-perfused bulk), tumor
/// cells clamped to `[0, 1]`. Supplies both O2 (× `basal_ros`) and drug.
///
/// Distances are computed in lattice units and scaled by `grid.cell_size_um`.
///
/// **Cost**: brute-force nearest-vessel, `O(tumor_cells × vessels)`. Cheap at
/// the 60³ matrix scale (~16M evals, one-time setup) but grows with tumor
/// volume (vessel count ∝ volume), so at patient scale (#240, e.g. a
/// well-vascularized 200³ ≈ 34B evals) it needs a spatial index (uniform grid
/// / kd-tree) for nearest-vessel. Deferred until #240 makes it bite.
///
/// # Panics
/// If `vessels` is empty (no source ⇒ undefined supply); callers pass the
/// output of [`place_vessels_3d`], which is floored at 1 vessel.
pub fn vessel_supply_field(
    grid: &TumorGrid3D,
    vessels: &[(f64, f64, f64)],
    lambda_um: f64,
) -> Vec<f64> {
    assert!(!vessels.is_empty(), "vessel_supply_field needs ≥1 vessel");
    debug_assert!(
        lambda_um.is_finite() && lambda_um > 0.0,
        "vessel_supply_field: lambda_um must be finite and positive, got {lambda_um}"
    );
    let cell_size = grid.cell_size_um;
    (0..grid.cells.len())
        .map(|idx| {
            if !grid.cells[idx].is_tumor {
                return 1.0;
            }
            let (r, c, l) = grid.coords(idx);
            let (rf, cf, lf) = (r as f64, c as f64, l as f64);
            let mut nearest_d2 = f64::INFINITY;
            for &(vr, vc, vl) in vessels {
                let d2 = (rf - vr).powi(2) + (cf - vc).powi(2) + (lf - vl).powi(2);
                if d2 < nearest_d2 {
                    nearest_d2 = d2;
                }
            }
            let dist_um = nearest_d2.sqrt() * cell_size;
            (-dist_um / lambda_um).exp().clamp(0.0, 1.0)
        })
        .collect()
}

/// Fraction of tumor cells whose supply factor is below `threshold` (the
/// hypoxic fraction). Used to compare the explicit-vessel field against the
/// edge-distance proxy (#191 AC: irregular vasculature shifts the hypoxic
/// fraction relative to a smooth radial gradient).
pub fn hypoxic_fraction(grid: &TumorGrid3D, field: &[f64], threshold: f64) -> f64 {
    let (hyp, tot) = grid
        .cells
        .iter()
        .zip(field)
        .fold((0usize, 0usize), |(h, t), (gc, &f)| {
            if gc.is_tumor {
                (h + usize::from(f < threshold), t + 1)
            } else {
                (h, t)
            }
        });
    if tot > 0 {
        hyp as f64 / tot as f64
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::oxygen::radial_o2_field;

    fn grid() -> TumorGrid3D {
        TumorGrid3D::generate(40, 40, 40, 20.0, 42)
    }

    #[test]
    fn placement_is_deterministic_and_nonempty() {
        let g = grid();
        let a = place_vessels_3d(&g, &VasculatureConfig::well_vascularized(), 7);
        let b = place_vessels_3d(&g, &VasculatureConfig::well_vascularized(), 7);
        assert_eq!(a, b);
        assert!(!a.is_empty());
        // Denser config places more vessels than a sparse one.
        let sparse = place_vessels_3d(&g, &VasculatureConfig::poorly_vascularized(), 7);
        assert!(
            a.len() > sparse.len(),
            "well={}, poor={}",
            a.len(),
            sparse.len()
        );
    }

    #[test]
    fn supply_is_one_for_stroma_and_in_range_for_tumor() {
        let g = grid();
        let v = place_vessels_3d(&g, &VasculatureConfig::well_vascularized(), 7);
        let field = vessel_supply_field(&g, &v, 100.0);
        assert_eq!(field.len(), g.cells.len());
        for (idx, &f) in field.iter().enumerate() {
            if g.cells[idx].is_tumor {
                assert!((0.0..=1.0).contains(&f), "tumor supply {f} out of [0,1]");
            } else {
                assert_eq!(f, 1.0, "stroma must be well-perfused (1.0)");
            }
        }
    }

    #[test]
    fn vessel_field_oxygenates_the_core_unlike_the_edge_proxy() {
        // The key #191 difference: the edge-distance proxy makes ALL deep-core
        // cells hypoxic (supply decays monotonically from the surface), whereas
        // explicit internal vessels reach the core, so some deep cells are
        // well-supplied. So among deep-core tumor cells, the mean vessel supply
        // exceeds the mean edge-proxy supply at matched λ — the irregular,
        // non-radial oxygenation the model is meant to capture.
        let g = grid();
        let lambda = 100.0;
        let edge = radial_o2_field(&g, lambda);
        let v = place_vessels_3d(&g, &VasculatureConfig::well_vascularized(), 7);
        let vessel = vessel_supply_field(&g, &v, lambda);

        let tumor_radius_um =
            (g.rows.min(g.cols).min(g.layers) as f64) * TUMOR_RADIUS_FRACTION * g.cell_size_um;
        let deep_threshold = 0.5 * tumor_radius_um;
        // Mean supply over deep-core tumor cells (depth > half the radius).
        let deep_mean = |field: &[f64]| -> f64 {
            let (sum, n) = (0..g.cells.len()).fold((0.0_f64, 0usize), |(s, n), idx| {
                let (r, c, l) = g.coords(idx);
                if g.cells[idx].is_tumor && g.radial_depth_um(r, c, l) > deep_threshold {
                    (s + field[idx], n + 1)
                } else {
                    (s, n)
                }
            });
            if n > 0 {
                sum / n as f64
            } else {
                0.0
            }
        };
        let (vessel_core, edge_core) = (deep_mean(&vessel), deep_mean(&edge));
        assert!(
            vessel_core > edge_core,
            "explicit vessels must oxygenate the deep core more than the edge proxy: \
             vessel_core_mean={vessel_core:.4}, edge_core_mean={edge_core:.4}"
        );
    }

    #[test]
    #[should_panic(expected = "needs ≥1 vessel")]
    fn empty_vessels_panics() {
        let g = grid();
        let _ = vessel_supply_field(&g, &[], 100.0);
    }
}
