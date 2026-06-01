//! Clonal heterogeneity (#242).
//!
//! Real tumors are genetic mosaics: 4–8+ subclones with measurably different
//! ferroptosis vulnerabilities coexist, spatially organized into patches.
//! High-mesenchymal (ZEB1+) subclones are highly ferroptosis-vulnerable;
//! epithelial subclones are resistant, and the between-subclone variance often
//! exceeds the between-treatment variance in real drug screens (Conrad et al.,
//! Nat Rev Mol Cell Biol 2018; Marusyk et al., Cancer Cell 2014; Heindl et al.,
//! Nat Methods 2019 for spatial clonal mosaics).
//!
//! ## Design: separate field + independent RNG (byte-identity)
//!
//! [`assign_subclones_3d`] returns a consumer-owned `Vec<u8>` of per-cell
//! subclone ids (`0` = non-tumor, `1..=k` = subclone) computed by Voronoi
//! assignment from `k` seed points. It draws from its **own** `StdRng`, so it
//! never perturbs [`TumorGrid3D::generate`](crate::grid::TumorGrid3D::generate)'s
//! RNG stream — the generated cell grid (phenotypes, persister clusters) is
//! bit-for-bit unchanged whether or not subclones are assigned. The consumer
//! applies [`SubclonePerturbation`]s as RNG-neutral setup mutations (like the
//! `oxygen`/`ph` gradients). A `k = 1` identity config is therefore a no-op
//! and byte-identical to having no clonal model.

use crate::grid::{TumorGrid3D, TUMOR_RADIUS_FRACTION};
use rand::prelude::*;

/// Per-subclone perturbation of the baseline ferroptosis parameters, applied
/// by the consumer as one-time setup mutations (RNG-neutral). All-identity
/// (the [`SubclonePerturbation::identity`] value) leaves the cell unchanged.
#[derive(Clone, Copy, Debug)]
pub struct SubclonePerturbation {
    /// Multiplier on labile iron (>1 ⇒ more Fenton ROS ⇒ more vulnerable).
    pub iron_mul: f64,
    /// Multiplier on the GPX4 reserve (<1 ⇒ less repair ⇒ more vulnerable).
    pub gpx4_mul: f64,
    /// Additive MUFA membrane protection (>0 ⇒ less oxidizable PUFA ⇒ more
    /// resistant). Consumer clamps the resulting protection to a sane range.
    pub mufa_add: f64,
}

impl SubclonePerturbation {
    /// The no-op perturbation: iron ×1, GPX4 ×1, MUFA +0.
    pub fn identity() -> Self {
        SubclonePerturbation {
            iron_mul: 1.0,
            gpx4_mul: 1.0,
            mufa_add: 0.0,
        }
    }

    /// True when applying this perturbation changes nothing.
    pub fn is_identity(&self) -> bool {
        self.iron_mul == 1.0 && self.gpx4_mul == 1.0 && self.mufa_add == 0.0
    }
}

/// Clonal-heterogeneity configuration: one [`SubclonePerturbation`] per
/// subclone (`perturbations.len() == k`; subclone ids are `1..=k`).
#[derive(Clone, Debug)]
pub struct ClonalConfig {
    pub perturbations: Vec<SubclonePerturbation>,
}

impl ClonalConfig {
    /// Number of subclones (`k`).
    pub fn k(&self) -> usize {
        self.perturbations.len()
    }

    /// A single identity subclone (`k = 1`, no perturbation). A run with this
    /// config is byte-identical to one with no clonal model — the K=1
    /// byte-identity guarantee.
    pub fn single_identity() -> Self {
        ClonalConfig {
            perturbations: vec![SubclonePerturbation::identity()],
        }
    }

    /// Literature-informed 4-subclone table (placeholders pending calibration).
    /// Spans the mesenchymal⇄epithelial ferroptosis-vulnerability axis:
    /// - 1 high-mesenchymal (ZEB1+): iron-loaded, GPX4-low ⇒ most vulnerable.
    /// - 2 intermediate-mesenchymal: mildly vulnerable.
    /// - 3 intermediate-epithelial: mildly resistant (MUFA-enriched).
    /// - 4 epithelial: GPX4-high, MUFA-enriched ⇒ most resistant.
    pub fn literature_4() -> Self {
        ClonalConfig {
            perturbations: vec![
                SubclonePerturbation {
                    iron_mul: 1.5,
                    gpx4_mul: 0.6,
                    mufa_add: 0.0,
                },
                SubclonePerturbation {
                    iron_mul: 1.2,
                    gpx4_mul: 0.85,
                    mufa_add: 0.0,
                },
                SubclonePerturbation {
                    iron_mul: 0.9,
                    gpx4_mul: 1.1,
                    mufa_add: 0.1,
                },
                SubclonePerturbation {
                    iron_mul: 0.7,
                    gpx4_mul: 1.3,
                    mufa_add: 0.2,
                },
            ],
        }
    }

    /// True when every subclone perturbation is the identity (no effect).
    pub fn is_identity(&self) -> bool {
        self.perturbations.iter().all(|p| p.is_identity())
    }
}

/// Assign each grid cell to a subclone via Voronoi tessellation from `k` seed
/// points placed uniformly in the tumor sphere. Returns a `Vec<u8>` of length
/// `grid.cells.len()`: `0` for non-tumor cells, `1..=k` for tumor cells (the
/// id of the nearest seed). Deterministic given `(grid dims, k, seed)`.
///
/// Uses an **independent** `StdRng(seed)` so it never advances the RNG used by
/// [`TumorGrid3D::generate`], preserving byte-identity of the cell grid.
///
/// # Panics
/// If `k == 0` or `k > 255` (ids must fit in `u8`, `0` reserved for stroma).
pub fn assign_subclones_3d(grid: &TumorGrid3D, k: usize, seed: u64) -> Vec<u8> {
    assert!(
        k >= 1 && k <= u8::MAX as usize,
        "k must be in 1..=255, got {k}"
    );
    let mut rng = StdRng::seed_from_u64(seed);
    let center = (
        grid.rows as f64 / 2.0,
        grid.cols as f64 / 2.0,
        grid.layers as f64 / 2.0,
    );
    let tumor_radius = (grid.rows.min(grid.cols).min(grid.layers) as f64) * TUMOR_RADIUS_FRACTION;

    // Seed points uniformly distributed in the tumor sphere (cbrt radial
    // sampling to avoid center bias — same convention as generate's clusters).
    let seeds: Vec<(f64, f64, f64)> = (0..k)
        .map(|_| {
            let dist = rng.gen::<f64>().cbrt() * tumor_radius * 0.9;
            let theta = rng.gen::<f64>() * std::f64::consts::TAU;
            let cos_phi = 2.0 * rng.gen::<f64>() - 1.0;
            let sin_phi = (1.0 - cos_phi * cos_phi).sqrt();
            (
                center.0 + dist * cos_phi,
                center.1 + dist * sin_phi * theta.cos(),
                center.2 + dist * sin_phi * theta.sin(),
            )
        })
        .collect();

    (0..grid.cells.len())
        .map(|idx| {
            if !grid.cells[idx].is_tumor {
                return 0u8;
            }
            let (r, c, l) = grid.coords(idx);
            let (rf, cf, lf) = (r as f64, c as f64, l as f64);
            let mut best = 0usize;
            let mut best_d2 = f64::INFINITY;
            for (si, &(sr, sc, sl)) in seeds.iter().enumerate() {
                let d2 = (rf - sr).powi(2) + (cf - sc).powi(2) + (lf - sl).powi(2);
                if d2 < best_d2 {
                    best_d2 = d2;
                    best = si;
                }
            }
            (best + 1) as u8
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identity_perturbation_is_noop() {
        assert!(SubclonePerturbation::identity().is_identity());
        assert!(ClonalConfig::single_identity().is_identity());
    }

    #[test]
    fn literature_4_is_not_identity_and_spans_vulnerability() {
        let c = ClonalConfig::literature_4();
        assert_eq!(c.k(), 4);
        assert!(!c.is_identity());
        // Subclone 1 is more vulnerable than subclone 4 on every axis.
        let v = c.perturbations[0];
        let r = c.perturbations[3];
        assert!(v.iron_mul > r.iron_mul && v.gpx4_mul < r.gpx4_mul && v.mufa_add < r.mufa_add);
    }

    #[test]
    fn assignment_is_deterministic_given_seed() {
        let g = TumorGrid3D::generate(20, 20, 20, 20.0, 42);
        let a = assign_subclones_3d(&g, 4, 7);
        let b = assign_subclones_3d(&g, 4, 7);
        assert_eq!(a, b);
        // A different seed generally moves at least one boundary cell.
        let c = assign_subclones_3d(&g, 4, 8);
        assert_ne!(a, c);
    }

    #[test]
    fn ids_are_zero_for_stroma_and_in_range_for_tumor() {
        let g = TumorGrid3D::generate(20, 20, 20, 20.0, 42);
        let ids = assign_subclones_3d(&g, 4, 7);
        assert_eq!(ids.len(), g.cells.len());
        for (idx, &id) in ids.iter().enumerate() {
            if g.cells[idx].is_tumor {
                assert!((1..=4).contains(&id), "tumor cell id {id} out of 1..=4");
            } else {
                assert_eq!(id, 0, "non-tumor cell must be subclone 0");
            }
        }
    }

    #[test]
    fn k1_assigns_every_tumor_cell_to_subclone_one() {
        let g = TumorGrid3D::generate(20, 20, 20, 20.0, 42);
        let ids = assign_subclones_3d(&g, 1, 7);
        for (idx, &id) in ids.iter().enumerate() {
            assert_eq!(id, u8::from(g.cells[idx].is_tumor));
        }
    }

    #[test]
    #[should_panic(expected = "k must be in 1..=255")]
    fn k_zero_panics() {
        let g = TumorGrid3D::generate(10, 10, 10, 20.0, 42);
        let _ = assign_subclones_3d(&g, 0, 7);
    }
}
