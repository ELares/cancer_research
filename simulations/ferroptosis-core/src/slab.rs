//! Patient-scale slab geometry (#240, Option A).
//!
//! sim-tme-3d's spheroid is ~540 µm radius — in-vitro scale. Real patient
//! tumors are 5–50 mm, where drug/O2 penetration (Krogh ~150 µm) drops
//! catastrophically with depth, so a deep tumor is mostly drug-deprived. This
//! module models a **slab**: an all-tumor block (see
//! [`TumorGrid3D::generate_slab`](crate::grid::TumorGrid3D::generate_slab))
//! representing a chunk at a configurable **depth** in a larger virtual tumor.
//!
//! Geometry: the **+z face** (top layer, `l = layers-1`) is vessel-proximal at
//! `depth_offset`; supply decays going **−z** (toward `l = 0`, deeper into the
//! tumor) as `exp(-depth/λ)`. Uniform across `(row, col)` within a layer. This
//! is the planar (1-D) analog of `oxygen::radial_o2_field`, offset by the
//! slab's depth, and like [`crate::vasculature`] it is a unified supply factor
//! the consumer applies to BOTH O2 (× `basal_ros`) and drug.
//!
//! Boundary conditions (#240): +z = vessel (the supply maximum, immune-
//! accessible via the existing DAMP cascade); −z = continuing tumor, a no-flux
//! / reflective face — already satisfied because the iron/DAMP diffusion uses
//! bounded (no-wrap) neighbors, i.e. no flux across any face.
//!
//! Opt-in ⇒ a consumer that doesn't request slab mode keeps the spheroid and
//! stays byte-identical.

use crate::grid::TumorGrid3D;

/// Krogh-style default O2/drug penetration length (µm) for the slab when the
/// condition doesn't specify one. ~150 µm is the canonical inter-capillary
/// half-distance / O2 diffusion length in tumor tissue.
pub const KROGH_LAMBDA_UM: f64 = 150.0;

/// Slab placement within a virtual large tumor.
#[derive(Clone, Copy, Debug)]
pub struct SlabConfig {
    /// Depth (mm) of the slab's +z (vessel-proximal) face from the supply
    /// source. The slab spans `[depth_offset_mm, depth_offset_mm + slab_mm]`.
    pub depth_offset_mm: f64,
    /// Diameter (mm) of the virtual tumor this slab is embedded in — used only
    /// for the human-readable scale-interpretation string.
    pub virtual_tumor_mm: f64,
}

impl SlabConfig {
    /// A deep slab in a 5 mm virtual tumor: +z face at 4 mm, so the slab spans
    /// ~4.0–5.2 mm — essentially fully drug/O2-deprived (the patient-scale
    /// penetration collapse the spheroid scale misses).
    pub fn patient_5mm() -> Self {
        SlabConfig {
            depth_offset_mm: 4.0,
            virtual_tumor_mm: 5.0,
        }
    }

    /// A shallow slab: +z face at the surface (0 mm), so the +z face is a
    /// well-perfused vessel and supply decays across the slab — the
    /// in-vitro-spheroid-equivalent control for the depth comparison.
    pub fn surface() -> Self {
        SlabConfig {
            depth_offset_mm: 0.0,
            virtual_tumor_mm: 5.0,
        }
    }
}

/// Per-cell planar depth-graded supply for a slab: `exp(-depth/λ)` where the
/// depth of layer `l` is `depth_offset_um + (layers-1 - l)·cell_size_um`
/// (the +z face `l = layers-1` is shallowest, `l = 0` deepest). Uniform across
/// `(row, col)`. Returns a `Vec<f64>` of length `grid.cells.len()`, clamped to
/// `[0, 1]`. Drop-in for `oxygen::radial_o2_field`; supplies O2 and drug.
pub fn slab_supply_field(grid: &TumorGrid3D, depth_offset_um: f64, lambda_um: f64) -> Vec<f64> {
    debug_assert!(
        depth_offset_um >= 0.0 && lambda_um.is_finite() && lambda_um > 0.0,
        "slab_supply_field: depth_offset_um >= 0 and lambda_um finite > 0; got {depth_offset_um}, {lambda_um}"
    );
    let cell_size = grid.cell_size_um;
    let top = grid.layers.saturating_sub(1);
    // Precompute per-layer supply (it varies only with l), then broadcast.
    let per_layer: Vec<f64> = (0..grid.layers)
        .map(|l| {
            let depth_um = depth_offset_um + (top - l) as f64 * cell_size;
            (-depth_um / lambda_um).exp().clamp(0.0, 1.0)
        })
        .collect();
    (0..grid.cells.len())
        .map(|idx| {
            let (_, _, l) = grid.coords(idx);
            per_layer[l]
        })
        .collect()
}

/// Human-readable interpretation of what depth/scale a slab run represents,
/// for the output JSON (#240 AC). e.g. "slab spanning depth 4.0–5.2 mm of a
/// 5 mm virtual tumor (1.2 mm thick)".
pub fn scale_interpretation(grid: &TumorGrid3D, cfg: &SlabConfig) -> String {
    let slab_mm = grid.layers as f64 * grid.cell_size_um / 1000.0;
    format!(
        "slab spanning depth {:.1}–{:.1} mm of a {:.0} mm virtual tumor ({:.1} mm thick)",
        cfg.depth_offset_mm,
        cfg.depth_offset_mm + slab_mm,
        cfg.virtual_tumor_mm,
        slab_mm,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn grid() -> TumorGrid3D {
        TumorGrid3D::generate_slab(20, 20, 20, 20.0, 42)
    }

    #[test]
    fn slab_grid_is_all_tumor() {
        let g = grid();
        assert_eq!(g.cells.len(), 8000);
        assert!(
            g.cells.iter().all(|c| c.is_tumor),
            "every slab cell is tumor"
        );
    }

    #[test]
    fn supply_decreases_with_depth() {
        let g = grid();
        // Surface slab (depth_offset 0): +z (top layer) is well-supplied,
        // −z (bottom) is deprived.
        let f = slab_supply_field(&g, 0.0, KROGH_LAMBDA_UM);
        let top = g.flat_index(10, 10, g.layers - 1);
        let bottom = g.flat_index(10, 10, 0);
        assert!(
            f[top] > f[bottom],
            "+z {} should exceed −z {}",
            f[top],
            f[bottom]
        );
        assert!((f[top] - 1.0).abs() < 1e-9, "+z at depth 0 is saturated");
        // Uniform across (row, col) within a layer.
        assert_eq!(f[g.flat_index(0, 0, 5)], f[g.flat_index(19, 19, 5)]);
    }

    #[test]
    fn deep_slab_is_essentially_deprived() {
        let g = grid();
        // 4 mm offset: even the +z face is ~exp(-4000/150) ≈ 3e-12.
        let f = slab_supply_field(&g, 4000.0, KROGH_LAMBDA_UM);
        assert!(
            f.iter().all(|&s| s < 1e-6),
            "a 4 mm-deep slab is fully deprived"
        );
    }

    #[test]
    fn scale_interpretation_reports_depth_span() {
        let g = grid(); // 20 layers × 20 µm = 0.4 mm thick
        let s = scale_interpretation(&g, &SlabConfig::patient_5mm());
        assert!(s.contains("4.0"), "mentions the offset depth: {s}");
        assert!(
            s.contains("5 mm virtual tumor"),
            "mentions the virtual size: {s}"
        );
    }
}
