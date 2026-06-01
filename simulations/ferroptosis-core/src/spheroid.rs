//! 3D spheroid radial cell biology (#197).
//!
//! `TumorGrid3D::generate` assigns phenotypes by a coarse core/periphery split
//! with random rolls. A real spheroid is more structured: a proliferating,
//! well-nourished **glycolytic rim**, a quiescent **OXPHOS** intermediate
//! zone, and a hypoxic, nutrient-deprived **persister-like core**. And within a
//! phenotype the cascade itself is position-dependent — peripheral cells
//! accumulate more MUFA, core cells start GSH-poor (cysteine-limited) and
//! iron-rich (HIF-driven import). This module makes those gradients explicit,
//! to be run with [`Params::spheroid`](crate::params::Params::spheroid).
//!
//! ## Design: opt-in re-assignment via an independent RNG (byte-identity)
//!
//! [`apply_radial_cells_3d`] **re-generates** each tumor cell with its
//! radial-position phenotype, drawing from its own per-cell `StdRng` seeded
//! from `seed` — it never touches `TumorGrid3D::generate`'s stream, so a
//! consumer that doesn't opt in keeps the random grid and stays byte-identical.
//! The phenotype change rewrites the cell's biochemical draw (unlike the
//! parameter-only `clonal`/`vasculature` layers), so this is a deliberate
//! "different tumor model", not a perturbation of the default one. Position-
//! dependent MUFA is a `CellState` value the consumer applies after init via
//! [`radial_mufa_protection`] (it needs the freshly-initialized state, and
//! relies on `Params::spheroid`'s partially-active SCD1 to persist).

use crate::cell::{gen_cell, Phenotype};
use crate::grid::{TumorGrid3D, TUMOR_RADIUS_FRACTION};
use rand::prelude::*;

/// Radial structure + gradient strengths for a spheroid. Fractions are of the
/// tumor radius (0 = center/core, 1 = surface).
#[derive(Clone, Copy, Debug)]
pub struct SpheroidConfig {
    /// At/above this radial fraction, cells are Glycolytic (the rim).
    pub glycolytic_frac: f64,
    /// At/above this (and below `glycolytic_frac`), cells are OXPHOS; below it,
    /// the Persister-like core.
    pub oxphos_frac: f64,
    /// Per-cell MUFA carrying **cap** at the surface (high) and core (low) —
    /// the `cell.mufa_cap` `update_mufa_protection` saturates toward, so the
    /// position-dependent MUFA is durable (#270), not a transient initial
    /// condition that converges to the uniform M_ss. The same radial value also
    /// seeds the initial `state.mufa_protection`.
    pub mufa_surface: f64,
    pub mufa_core: f64,
    /// Initial-GSH multiplier at the core (cysteine-limited; < 1). Surface = 1.
    pub gsh_core_factor: f64,
    /// Labile-iron multiplier at the core (HIF-driven import; > 1). Surface = 1.
    pub iron_core_factor: f64,
}

impl SpheroidConfig {
    /// Literature-informed defaults (placeholders pending calibration):
    /// outer third glycolytic, middle third OXPHOS, inner third persister-like;
    /// MUFA 0.25→0.05 surface→core; core GSH ×0.5; core iron ×1.6.
    ///
    /// MUFA is now a **durable** position-dependent axis (#270): `mufa_surface`
    /// / `mufa_core` are the per-cell MUFA carrying caps (`cell.mufa_cap`), so
    /// `update_mufa_protection` relaxes each cell toward a steady state that
    /// scales with its cap instead of every cell converging to the global
    /// uniform M_ss. (Previously these set only the *initial* `state.mufa_protection`,
    /// which relaxed back — the #197-review transience caveat, resolved here.)
    /// `iron_core_factor` (a static `cell.iron` scale) is likewise durable;
    /// `gsh_core_factor` sets the *initial* GSH, which then evolves under NRF2
    /// resynthesis. Values remain placeholders pending calibration.
    ///
    /// **Zone geometry caveat**: the thresholds are *radial* fractions, so by
    /// volume (∝ r³) the persister core (`frac < 0.33`) is only ~4% and the
    /// glycolytic rim (`frac ≥ 0.66`) ~71%. Real spheroids have a thin
    /// proliferating rim and a larger quiescent/hypoxic core; tilting the core
    /// larger (raising `oxphos_frac`) is a calibration follow-up.
    pub fn literature() -> Self {
        SpheroidConfig {
            glycolytic_frac: 0.66,
            oxphos_frac: 0.33,
            // Per-cell MUFA caps (#270): rim saturates near the global cap,
            // core saturates lower → durable rim-vs-core MUFA spread.
            mufa_surface: 0.25,
            mufa_core: 0.05,
            gsh_core_factor: 0.5,
            iron_core_factor: 1.6,
        }
    }
}

/// Radial fraction of a cell: `0.0` at the spheroid center, `1.0` at the
/// surface (clamped). Geometry matches `TumorGrid3D::generate` (center =
/// dims/2, radius = min(dims) × `TUMOR_RADIUS_FRACTION`).
pub fn radial_fraction_3d(grid: &TumorGrid3D, idx: usize) -> f64 {
    let (r, c, l) = grid.coords(idx);
    let center = (
        grid.rows as f64 / 2.0,
        grid.cols as f64 / 2.0,
        grid.layers as f64 / 2.0,
    );
    let tumor_radius = (grid.rows.min(grid.cols).min(grid.layers) as f64) * TUMOR_RADIUS_FRACTION;
    let dist = ((r as f64 - center.0).powi(2)
        + (c as f64 - center.1).powi(2)
        + (l as f64 - center.2).powi(2))
    .sqrt();
    (dist / tumor_radius).clamp(0.0, 1.0)
}

/// Phenotype for a cell at radial fraction `frac`: glycolytic rim → OXPHOS mid
/// → persister-like core.
pub fn radial_phenotype(frac: f64, cfg: &SpheroidConfig) -> Phenotype {
    if frac >= cfg.glycolytic_frac {
        Phenotype::Glycolytic
    } else if frac >= cfg.oxphos_frac {
        Phenotype::OXPHOS
    } else {
        Phenotype::Persister
    }
}

#[inline]
fn lerp_core_surface(core: f64, surface: f64, frac: f64) -> f64 {
    core + (surface - core) * frac.clamp(0.0, 1.0)
}

/// Position-dependent MUFA level: high at the surface, low at the core. Used
/// for **both** the per-cell MUFA cap (`cell.mufa_cap`, set in
/// [`apply_radial_cells_3d`] so the value is durable, #270) and the consumer's
/// initial `state.mufa_protection`.
pub fn radial_mufa_protection(frac: f64, cfg: &SpheroidConfig) -> f64 {
    lerp_core_surface(cfg.mufa_core, cfg.mufa_surface, frac)
}

/// Initial-GSH multiplier: `gsh_core_factor` (< 1) at the core, 1.0 at surface.
pub fn radial_gsh_factor(frac: f64, cfg: &SpheroidConfig) -> f64 {
    lerp_core_surface(cfg.gsh_core_factor, 1.0, frac)
}

/// Labile-iron multiplier: `iron_core_factor` (> 1) at the core, 1.0 at surface.
pub fn radial_iron_factor(frac: f64, cfg: &SpheroidConfig) -> f64 {
    lerp_core_surface(cfg.iron_core_factor, 1.0, frac)
}

/// Re-assign every tumor cell radially: re-generate it with its
/// radial-position phenotype (per-cell independent RNG), then scale `cell.gsh`
/// (core-low) and `cell.iron` (core-high) by the radial gradients and set the
/// per-cell `cell.mufa_cap` (rim-high, core-low; #270) so the position-dependent
/// MUFA is **durable** — the MUFA steady state scales with this cap rather than
/// every cell converging to the global uniform M_ss. Non-tumor (stromal) cells
/// are left untouched. Deterministic given `(grid dims, cfg, seed)`. The
/// consumer still sets the *initial* `state.mufa_protection` after state init
/// via [`radial_mufa_protection`] (the same radial value the cap uses), so the
/// cell both starts at and relaxes toward its position-dependent level.
pub fn apply_radial_cells_3d(grid: &mut TumorGrid3D, cfg: &SpheroidConfig, seed: u64) {
    for idx in 0..grid.cells.len() {
        if !grid.cells[idx].is_tumor {
            continue;
        }
        let frac = radial_fraction_3d(grid, idx); // immutable borrow ends here
        let pheno = radial_phenotype(frac, cfg);
        let mut rng = StdRng::seed_from_u64(seed.wrapping_add(idx as u64));
        let mut cell = gen_cell(pheno, &mut rng);
        cell.gsh *= radial_gsh_factor(frac, cfg);
        cell.iron *= radial_iron_factor(frac, cfg);
        cell.mufa_cap = Some(radial_mufa_protection(frac, cfg));
        grid.cells[idx].cell = cell;
        grid.cells[idx].phenotype = pheno;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cfg() -> SpheroidConfig {
        SpheroidConfig::literature()
    }

    #[test]
    fn phenotype_zones_run_rim_to_core() {
        let c = cfg();
        assert_eq!(radial_phenotype(0.95, &c), Phenotype::Glycolytic); // rim
        assert_eq!(radial_phenotype(0.5, &c), Phenotype::OXPHOS); // mid
        assert_eq!(radial_phenotype(0.1, &c), Phenotype::Persister); // core
    }

    #[test]
    fn gradients_have_correct_direction_and_endpoints() {
        let c = cfg();
        // Surface endpoints are the neutral values; core endpoints are the cfg.
        assert!((radial_gsh_factor(1.0, &c) - 1.0).abs() < 1e-12);
        assert!((radial_iron_factor(1.0, &c) - 1.0).abs() < 1e-12);
        assert!((radial_gsh_factor(0.0, &c) - c.gsh_core_factor).abs() < 1e-12);
        assert!((radial_iron_factor(0.0, &c) - c.iron_core_factor).abs() < 1e-12);
        // Core is GSH-poor (<1) and iron-rich (>1); MUFA higher at surface.
        assert!(radial_gsh_factor(0.0, &c) < radial_gsh_factor(1.0, &c));
        assert!(radial_iron_factor(0.0, &c) > radial_iron_factor(1.0, &c));
        assert!(radial_mufa_protection(1.0, &c) > radial_mufa_protection(0.0, &c));
    }

    #[test]
    fn radial_fraction_is_zero_at_center_and_high_at_surface() {
        let g = TumorGrid3D::generate(40, 40, 40, 20.0, 42);
        let center_idx = g.flat_index(20, 20, 20);
        // A near-surface tumor cell: scan a ray out from center along +row.
        let mut surface_idx = center_idx;
        for r in 20..40 {
            let i = g.flat_index(r, 20, 20);
            if g.cells[i].is_tumor {
                surface_idx = i;
            }
        }
        assert!(radial_fraction_3d(&g, center_idx) < 0.1);
        assert!(radial_fraction_3d(&g, surface_idx) > radial_fraction_3d(&g, center_idx));
    }

    #[test]
    fn re_assignment_is_deterministic_and_radial() {
        let mut a = TumorGrid3D::generate(30, 30, 30, 20.0, 42);
        let mut b = TumorGrid3D::generate(30, 30, 30, 20.0, 42);
        apply_radial_cells_3d(&mut a, &cfg(), 7);
        apply_radial_cells_3d(&mut b, &cfg(), 7);
        // Deterministic: same seed → identical phenotypes + cell draws.
        let phenos_a: Vec<_> = a.cells.iter().map(|gc| gc.phenotype).collect();
        let phenos_b: Vec<_> = b.cells.iter().map(|gc| gc.phenotype).collect();
        assert_eq!(phenos_a, phenos_b);
        // The center tumor cell is Persister-like; a surface tumor cell is not.
        let center = a.flat_index(15, 15, 15);
        assert_eq!(a.cells[center].phenotype, Phenotype::Persister);
        // Stromal (non-tumor) cells are untouched.
        for gc in &a.cells {
            if !gc.is_tumor {
                assert_eq!(gc.phenotype, Phenotype::Stromal);
            }
        }
    }

    /// #270 wiring: `apply_radial_cells_3d` sets the per-cell MUFA cap radially
    /// (core-low, rim-high) and leaves stromal cells uncapped (`None`), so the
    /// durable-MUFA mechanism gets a per-cell cap to act on.
    #[test]
    fn apply_sets_radial_mufa_cap_core_below_rim() {
        let mut g = TumorGrid3D::generate(40, 40, 40, 20.0, 42);
        apply_radial_cells_3d(&mut g, &cfg(), 7);

        // A core tumor cell and a near-surface tumor cell along +row from center.
        let center = g.flat_index(20, 20, 20);
        let mut surface = center;
        for r in 20..40 {
            let i = g.flat_index(r, 20, 20);
            if g.cells[i].is_tumor {
                surface = i;
            }
        }
        let core_cap = g.cells[center]
            .cell
            .mufa_cap
            .expect("core tumor cell is capped");
        let rim_cap = g.cells[surface]
            .cell
            .mufa_cap
            .expect("rim tumor cell is capped");
        assert!(
            core_cap < rim_cap,
            "core MUFA cap should be below rim: core={core_cap}, rim={rim_cap}"
        );
        // Caps lie within the configured [core, surface] band.
        let c = cfg();
        assert!(core_cap >= c.mufa_core - 1e-9 && rim_cap <= c.mufa_surface + 1e-9);
        // Stromal (non-tumor) cells are left uncapped.
        for gc in &g.cells {
            if !gc.is_tumor {
                assert!(gc.cell.mufa_cap.is_none(), "stroma stays uncapped");
            }
        }
    }
}
