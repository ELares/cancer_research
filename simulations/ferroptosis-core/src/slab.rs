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

use crate::cell::{gen_cell, Phenotype};
use crate::grid::TumorGrid3D;
use rand::prelude::*;

/// Krogh-style default O2/drug penetration length (µm) for the slab when the
/// condition doesn't specify one. ~150 µm is the canonical inter-capillary
/// half-distance / O2 diffusion length in tumor tissue.
pub const KROGH_LAMBDA_UM: f64 = 150.0;

/// Diameter (mm) of the virtual patient-scale tumor the slab is embedded in.
/// 10 mm is a representative tumor at presentation (real tumors span 5–50 mm)
/// and comfortably contains the deep slab at production grid sizes, so the
/// `scale_interpretation` string never reports a slab extending past its own
/// tumor (see the `debug_assert` in [`scale_interpretation`]).
pub const VIRTUAL_TUMOR_MM: f64 = 10.0;

/// Slab placement within a virtual large tumor.
#[derive(Clone, Copy, Debug)]
pub struct SlabConfig {
    /// Depth (mm) of the slab's +z (vessel-proximal) face from the supply
    /// source. The slab spans `[depth_offset_mm, depth_offset_mm + slab_mm]`,
    /// which must stay within `virtual_tumor_mm`.
    pub depth_offset_mm: f64,
    /// Diameter (mm) of the virtual tumor this slab is embedded in — used for
    /// the human-readable scale-interpretation string and as the upper bound
    /// the slab span is checked against.
    pub virtual_tumor_mm: f64,
}

impl SlabConfig {
    /// A deep slab in a 10 mm patient-scale tumor: +z face at 4 mm, so the
    /// slab spans ~4.0–5.2 mm at the production grid size — essentially fully
    /// drug/O2-deprived (the patient-scale penetration collapse the spheroid
    /// scale misses). The 10 mm virtual tumor contains the slab with margin.
    pub fn patient_deep() -> Self {
        SlabConfig {
            depth_offset_mm: 4.0,
            virtual_tumor_mm: VIRTUAL_TUMOR_MM,
        }
    }

    /// A shallow slab: +z face at the surface (0 mm), so the +z face is a
    /// well-perfused vessel and supply decays across the slab — the
    /// in-vitro-spheroid-equivalent control for the depth comparison. Same
    /// virtual tumor as [`patient_deep`](Self::patient_deep) so the two are an
    /// apples-to-apples depth comparison.
    pub fn surface() -> Self {
        SlabConfig {
            depth_offset_mm: 0.0,
            virtual_tumor_mm: VIRTUAL_TUMOR_MM,
        }
    }
}

/// Per-**layer** planar depth-graded supply: a `Vec<f64>` of length
/// `grid.layers` where entry `l` is `exp(-depth/λ)`, depth =
/// `depth_offset_um + (layers-1 - l)·cell_size_um` (the +z face `l = layers-1`
/// is shallowest, `l = 0` deepest), clamped to `[0, 1]`. Supply varies only
/// with depth, so this is the compact form: [`slab_supply_field`] broadcasts it
/// across `(row, col)` and [`apply_depth_graded_cells_3d`] thresholds on it, so
/// both share one source of truth for the depth formula.
pub fn layer_supply(grid: &TumorGrid3D, depth_offset_um: f64, lambda_um: f64) -> Vec<f64> {
    debug_assert!(
        depth_offset_um >= 0.0 && lambda_um.is_finite() && lambda_um > 0.0,
        "layer_supply: depth_offset_um >= 0 and lambda_um finite > 0; got {depth_offset_um}, {lambda_um}"
    );
    let cell_size = grid.cell_size_um;
    let top = grid.layers.saturating_sub(1);
    (0..grid.layers)
        .map(|l| {
            let depth_um = depth_offset_um + (top - l) as f64 * cell_size;
            (-depth_um / lambda_um).exp().clamp(0.0, 1.0)
        })
        .collect()
}

/// Per-cell planar depth-graded supply for a slab: `exp(-depth/λ)` where the
/// depth of layer `l` is `depth_offset_um + (layers-1 - l)·cell_size_um`
/// (the +z face `l = layers-1` is shallowest, `l = 0` deepest). Uniform across
/// `(row, col)`. Returns a `Vec<f64>` of length `grid.cells.len()`, clamped to
/// `[0, 1]`. Drop-in for `oxygen::radial_o2_field`; supplies O2 and drug.
pub fn slab_supply_field(grid: &TumorGrid3D, depth_offset_um: f64, lambda_um: f64) -> Vec<f64> {
    let per_layer = layer_supply(grid, depth_offset_um, lambda_um);
    (0..grid.cells.len())
        .map(|idx| {
            let (_, _, l) = grid.coords(idx);
            per_layer[l]
        })
        .collect()
}

/// Depth-graded phenotype zones for a slab (#272). `generate_slab` assigns a
/// flat bulk mix (no spatial structure); a real chunk of tumor at depth is
/// layered: vessel-proximal cells are proliferating, chronically supply-starved
/// deep cells are quiescent/persister-like. Unlike the spheroid's geometric
/// (volume-fraction) zones, the slab models an *absolute* depth, so its
/// phenotype tracks the **planar supply** `exp(-depth/λ)` that already shapes
/// its O2/drug field. Thresholds are on that supply value (∈ [0, 1]).
#[derive(Clone, Copy, Debug)]
pub struct SlabPhenotypeConfig {
    /// Supply (∈ [0, 1]) at/above which a cell is Glycolytic (well-perfused,
    /// proliferating), the vessel-proximal +z layers.
    pub glycolytic_supply: f64,
    /// Supply at/above which (and below `glycolytic_supply`) a cell is OXPHOS
    /// (quiescent intermediate); below it, Persister-like (chronically
    /// supply-deprived, drug-tolerant), the deep (−z) layers.
    pub oxphos_supply: f64,
}

impl SlabPhenotypeConfig {
    /// Heuristic placeholder thresholds: a proliferating glycolytic zone where
    /// relative supply ≥ 0.5, a quiescent OXPHOS intermediate down to 0.15, and
    /// a persister-like core below that. UNLIKE the spheroid's
    /// literature-grounded zone *volumes* (Browning 2021), these supply
    /// cut-points are uncalibrated (see the CALIBRATION_STATUS slab row): the
    /// result is the DIRECTION (deep, supply-starved tissue is enriched for
    /// tolerant phenotypes), not the exact layer counts. The realized zone
    /// thicknesses depend on λ and the slab's `depth_offset_mm`: a
    /// [`SlabConfig::patient_deep`] slab whose every layer sits below 0.15
    /// supply is uniformly persister-like, which is the intended behavior for a
    /// 4 mm-deep chunk (drug/O2 essentially never reach it).
    pub fn literature() -> Self {
        SlabPhenotypeConfig {
            glycolytic_supply: 0.5,
            oxphos_supply: 0.15,
        }
    }
}

/// Phenotype for a slab cell at planar supply `supply` (∈ [0, 1], the
/// `slab_supply_field` / [`layer_supply`] value): glycolytic (well-perfused) →
/// OXPHOS (intermediate) → persister-like (chronically deprived). Monotone in
/// supply, so a deeper (lower-supply) cell is never assigned a *more*
/// proliferative phenotype than a shallower one.
pub fn depth_phenotype(supply: f64, cfg: &SlabPhenotypeConfig) -> Phenotype {
    let s = supply.clamp(0.0, 1.0);
    if s >= cfg.glycolytic_supply {
        Phenotype::Glycolytic
    } else if s >= cfg.oxphos_supply {
        Phenotype::OXPHOS
    } else {
        Phenotype::Persister
    }
}

/// Re-assign every (tumor) slab cell's phenotype by its layer's planar supply
/// `exp(-depth/λ)` (#272): the vessel-proximal +z layers become
/// proliferating/glycolytic, the chronically supply-deprived deep (−z) layers
/// become persister-like, the depth-axis analog of the spheroid's rim→core
/// structure ([`crate::spheroid::apply_radial_cells_3d`]).
///
/// Like the spheroid re-assignment, each cell is re-generated from its OWN
/// per-cell `StdRng` seeded from `seed`, so `generate_slab`'s RNG stream is
/// untouched, so a consumer that doesn't opt in keeps the flat bulk mix and stays
/// byte-identical. `depth_offset_um` / `lambda_um` MUST match the values the
/// consumer passes to [`slab_supply_field`], so the phenotype gradient and the
/// O2/drug supply gradient are coherent. Deterministic given
/// `(grid dims, depth_offset, lambda, cfg, seed)`.
///
/// Two scoping notes (documented in CALIBRATION_STATUS): (1) the supply used
/// here is the **planar depth** gradient only; internal vessels (#272
/// coupling) raise the *delivered* supply dynamically downstream but do not
/// reshape the chronic phenotype here, a future refinement; (2) only the
/// phenotype is depth-graded; the per-cell biochemical draw is `gen_cell`'s
/// phenotype default, since the supply field already deprives deep cells of
/// O2/drug dynamically (no separate static GSH/iron gradient as in the
/// spheroid).
pub fn apply_depth_graded_cells_3d(
    grid: &mut TumorGrid3D,
    depth_offset_um: f64,
    lambda_um: f64,
    cfg: &SlabPhenotypeConfig,
    seed: u64,
) {
    let per_layer = layer_supply(grid, depth_offset_um, lambda_um);
    for idx in 0..grid.cells.len() {
        if !grid.cells[idx].is_tumor {
            continue;
        }
        let (_, _, l) = grid.coords(idx);
        let pheno = depth_phenotype(per_layer[l], cfg);
        let mut rng = StdRng::seed_from_u64(seed.wrapping_add(idx as u64));
        let cell = gen_cell(pheno, &mut rng);
        grid.cells[idx].cell = cell;
        grid.cells[idx].phenotype = pheno;
    }
}

/// Human-readable interpretation of what depth/scale a slab run represents,
/// for the output JSON (#240 AC). e.g. "slab spanning depth 4.0–5.2 mm of a
/// 10 mm virtual tumor (1.2 mm thick)".
pub fn scale_interpretation(grid: &TumorGrid3D, cfg: &SlabConfig) -> String {
    let slab_mm = grid.layers as f64 * grid.cell_size_um / 1000.0;
    debug_assert!(
        cfg.depth_offset_mm + slab_mm <= cfg.virtual_tumor_mm,
        "slab span {:.1}–{:.1} mm exceeds its {:.1} mm virtual tumor — bump virtual_tumor_mm or reduce depth_offset_mm",
        cfg.depth_offset_mm,
        cfg.depth_offset_mm + slab_mm,
        cfg.virtual_tumor_mm,
    );
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
    fn layer_supply_matches_broadcast_field() {
        // The refactor invariant: slab_supply_field broadcasts layer_supply, so
        // every cell's field value equals its layer's per-layer supply. Guards
        // against the two formulas silently diverging.
        let g = grid();
        let per_layer = layer_supply(&g, 0.0, KROGH_LAMBDA_UM);
        let field = slab_supply_field(&g, 0.0, KROGH_LAMBDA_UM);
        assert_eq!(per_layer.len(), g.layers);
        for idx in 0..g.cells.len() {
            let (_, _, l) = g.coords(idx);
            assert_eq!(field[idx], per_layer[l]);
        }
    }

    #[test]
    fn depth_phenotype_thresholds_and_monotone() {
        let c = SlabPhenotypeConfig::literature();
        // Well-perfused → glycolytic; intermediate → OXPHOS; starved → persister.
        assert_eq!(depth_phenotype(1.0, &c), Phenotype::Glycolytic);
        assert_eq!(depth_phenotype(0.5, &c), Phenotype::Glycolytic); // boundary is inclusive
        assert_eq!(depth_phenotype(0.3, &c), Phenotype::OXPHOS);
        assert_eq!(depth_phenotype(0.15, &c), Phenotype::OXPHOS); // boundary is inclusive
        assert_eq!(depth_phenotype(0.05, &c), Phenotype::Persister);
        assert_eq!(depth_phenotype(0.0, &c), Phenotype::Persister);
        // Monotone in supply: a more-deprived (lower-supply) cell is never more
        // proliferative than a better-supplied one. Rank Glycolytic > OXPHOS >
        // Persister and check the assigned rank is non-decreasing in supply.
        let rank = |p: Phenotype| match p {
            Phenotype::Glycolytic => 2,
            Phenotype::OXPHOS => 1,
            _ => 0,
        };
        let mut prev = -1i32;
        for i in 0..=100 {
            let s = i as f64 / 100.0;
            let r = rank(depth_phenotype(s, &c));
            assert!(r >= prev, "rank dropped as supply rose at s={s}");
            prev = r;
        }
        // Out-of-range supply is clamped, not panicking.
        assert_eq!(depth_phenotype(2.0, &c), Phenotype::Glycolytic);
        assert_eq!(depth_phenotype(-1.0, &c), Phenotype::Persister);
    }

    #[test]
    fn depth_grading_is_deterministic_and_layered() {
        let cfg = SlabPhenotypeConfig::literature();
        let mut a = grid();
        let mut b = grid();
        // Surface slab (offset 0): top (+z) well-perfused, bottom (−z) deprived.
        apply_depth_graded_cells_3d(&mut a, 0.0, KROGH_LAMBDA_UM, &cfg, 7);
        apply_depth_graded_cells_3d(&mut b, 0.0, KROGH_LAMBDA_UM, &cfg, 7);
        // Deterministic: same (grid, depth, λ, cfg, seed) → identical phenotypes
        // AND identical cell draws.
        let ph_a: Vec<_> = a.cells.iter().map(|gc| gc.phenotype).collect();
        let ph_b: Vec<_> = b.cells.iter().map(|gc| gc.phenotype).collect();
        assert_eq!(ph_a, ph_b);
        assert_eq!(a.cells[0].cell.gpx4, b.cells[0].cell.gpx4);
        // Still all-tumor, every cell re-assigned (slab has no stroma).
        assert!(a.cells.iter().all(|c| c.is_tumor));
        // +z face is glycolytic (supply 1.0); deep −z face is persister.
        let top = a.flat_index(10, 10, a.layers - 1);
        let bottom = a.flat_index(10, 10, 0);
        assert_eq!(a.cells[top].phenotype, Phenotype::Glycolytic);
        assert_eq!(a.cells[bottom].phenotype, Phenotype::Persister);
        // A whole +z layer is uniform (supply is constant within a layer).
        let l = a.layers - 1;
        let p0 = a.cells[a.flat_index(0, 0, l)].phenotype;
        assert!(a
            .cells
            .iter()
            .enumerate()
            .filter(|(idx, _)| a.coords(*idx).2 == l)
            .all(|(_, gc)| gc.phenotype == p0));
    }

    #[test]
    fn deep_slab_is_uniformly_persister() {
        // A 4 mm-deep slab: even the +z face supply is exp(-4000/150) ≈ 3e-12,
        // far below oxphos_supply (0.15), so every layer is persister-like, the
        // intended behavior for a chunk drug/O2 essentially never reach.
        let cfg = SlabPhenotypeConfig::literature();
        let mut g = grid();
        apply_depth_graded_cells_3d(&mut g, 4000.0, KROGH_LAMBDA_UM, &cfg, 7);
        assert!(
            g.cells.iter().all(|c| c.phenotype == Phenotype::Persister),
            "a 4 mm-deep slab is uniformly persister-like"
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
        let s = scale_interpretation(&g, &SlabConfig::patient_deep());
        assert!(s.contains("4.0"), "mentions the offset depth: {s}");
        assert!(
            s.contains("10 mm virtual tumor"),
            "mentions the virtual size: {s}"
        );
    }

    /// At the production grid size (60 layers × 20 µm = 1.2 mm) the deep slab
    /// spans 4.0–5.2 mm, which must still fit inside its 10 mm virtual tumor —
    /// the `scale_interpretation` debug_assert (review #1) holds, so the
    /// reported string never describes a slab poking out of its own tumor.
    #[test]
    fn deep_slab_fits_within_its_virtual_tumor_at_production_size() {
        let g = TumorGrid3D::generate_slab(60, 60, 60, 20.0, 42);
        let cfg = SlabConfig::patient_deep();
        let slab_mm = g.layers as f64 * g.cell_size_um / 1000.0;
        assert!(
            cfg.depth_offset_mm + slab_mm <= cfg.virtual_tumor_mm,
            "span {:.1} mm must fit in {:.1} mm tumor",
            cfg.depth_offset_mm + slab_mm,
            cfg.virtual_tumor_mm
        );
        // The string must not contradict itself; this would panic on the
        // debug_assert if the span exceeded the tumor.
        let _ = scale_interpretation(&g, &cfg);
    }

    /// Reflective / no-flux −z boundary (#240 scope): the slab's deep face
    /// (`l = 0`) must not exchange iron/DAMP with the +z vessel face
    /// (`l = layers-1`). The 26-Moore neighborhood is bounded (no wrap), so a
    /// cell at `l = 0` has NO neighbor at `l = layers-1` — i.e. there is no
    /// flux across the −z face. Verified directly here rather than only argued.
    #[test]
    fn minus_z_face_is_reflective_no_wraparound() {
        let g = grid(); // 20³
        let top = g.layers - 1;
        // A −z face cell and a +z face cell.
        for &(r, c) in &[(0usize, 0usize), (10, 10), (19, 19)] {
            let (neigh, count) = g.neighbors(r, c, 0);
            assert!(
                neigh[..count].iter().all(|&(_, _, nl)| nl != top),
                "(-z) cell at l=0 must have no neighbor at the +z face l={top}"
            );
            assert!(
                neigh[..count].iter().all(|&(_, _, nl)| nl <= 1),
                "(-z) cell neighbors stay within l ∈ {{0, 1}} (no wrap)"
            );
        }
    }
}
