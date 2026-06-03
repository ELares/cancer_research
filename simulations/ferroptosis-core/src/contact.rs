//! Cell-cell contact-mediated ferroptosis resistance (#270).
//!
//! Densely-packed, highly-contacting tumor cells resist ferroptosis: E-cadherin
//! junctions activate the Merlin (NF2)/Hippo pathway, which inhibits YAP and so
//! suppresses the YAP target genes ACSL4 (PUFA incorporation into membranes) and
//! TFRC (transferrin-receptor iron import) that fuel lipid peroxidation. Sparse
//! or surface cells (few contacts) lose this brake and are MORE
//! ferroptosis-sensitive. The effect is stronger in 3D, where interior cells
//! have up to 26 neighbours vs 8 in 2D.
//!
//! Ref: Wu J, Minikes AM, Gao M, et al. "Intercellular interaction dictates
//! cancer cell ferroptosis via NF2-YAP signalling." Nature 2019 (PMID 31341276).
//!
//! ## Model
//!
//! For each tumor cell, `contact = (tumor 26-Moore neighbours) / 26 ∈ [0, 1]`.
//! Resistance scales the two durable YAP-target axes DOWN with contact:
//!
//! ```text
//! cell.lipid_unsat *= max(0, 1 - lipid_strength * contact)   // ACSL4 / PUFA
//! cell.iron        *= max(0, 1 - iron_strength  * contact)    // TFRC / labile iron
//! ```
//!
//! Both are static [`Cell`](crate::cell::Cell) fields, so the resistance is
//! **durable** for the whole run (the autocatalytic LP cascade reads them every
//! step), not just an initial condition. This is effectively a per-cell
//! ferroptosis-threshold modulation: less peroxidizable substrate + less Fenton
//! iron ⇒ a denser cell needs more ROS to tip into death.
//!
//! ## Discipline
//!
//! Purely **geometric** — it reads `is_tumor` (set by grid generation) and uses
//! no RNG, so it never perturbs [`TumorGrid3D::generate`]'s stream. Off-by-default
//! ([`ContactConfig::default`] / [`is_identity`](ContactConfig::is_identity) ⇒
//! not applied) keeps the default matrix byte-identical. Composes
//! multiplicatively with the clonal `lipid_unsat_mul` / `iron_mul` axes.

use crate::grid::TumorGrid3D;

/// Maximum 3D Moore neighbourhood (a fully-interior cell), the denominator for
/// the contact fraction. A surface cell has fewer tumor neighbours ⇒ lower
/// contact ⇒ less resistance.
pub const MAX_NEIGHBORS_3D: f64 = 26.0;

/// Contact-resistance configuration. Both strengths default to 0 (identity ⇒
/// byte-identical). A strength `s` reduces the corresponding axis by up to a
/// factor `s` for a fully-contacted (contact = 1) interior cell.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactConfig {
    /// Fraction by which a fully-contacted cell's `lipid_unsat` (ACSL4/PUFA
    /// axis, the headline YAP target) is reduced.
    pub lipid_strength: f64,
    /// Fraction by which a fully-contacted cell's `iron` (TFRC axis) is reduced.
    pub iron_strength: f64,
}

impl Default for ContactConfig {
    fn default() -> Self {
        // Identity: no contact modulation ⇒ byte-identical.
        ContactConfig {
            lipid_strength: 0.0,
            iron_strength: 0.0,
        }
    }
}

impl ContactConfig {
    /// Literature-motivated (Wu 2019) qualitative strengths: a dense interior
    /// cell gets a marked PUFA-substrate reduction (ACSL4 is the headline YAP
    /// target) and a smaller iron reduction. **Magnitudes are UNCALIBRATED
    /// placeholders** encoding the documented direction (contact ⇒ resistance),
    /// not fit to data; calibrate vs density-resolved ferroptosis-sensitivity
    /// assays (sparse vs confluent culture; NF2/YAP knockdown).
    pub fn literature() -> Self {
        ContactConfig {
            lipid_strength: 0.4,
            iron_strength: 0.2,
        }
    }

    /// True when the config applies no modulation (both strengths 0).
    pub fn is_identity(&self) -> bool {
        self.lipid_strength == 0.0 && self.iron_strength == 0.0
    }
}

/// Contact fraction for one cell: `(tumor 26-Moore neighbours) / 26 ∈ [0, 1]`.
/// Non-tumor neighbours and out-of-bounds positions count as "no contact", so
/// surface cells (stroma on one side) score lower than interior cells.
pub fn contact_fraction_3d(grid: &TumorGrid3D, idx: usize) -> f64 {
    let (r, c, l) = grid.coords(idx);
    let (nbrs, n) = grid.neighbors(r, c, l);
    let tumor = nbrs[..n]
        .iter()
        .filter(|&&(nr, nc, nl)| grid.cells[grid.flat_index(nr, nc, nl)].is_tumor)
        .count();
    tumor as f64 / MAX_NEIGHBORS_3D
}

/// Apply contact-mediated resistance to every tumor cell (#270): scale the
/// durable `lipid_unsat` and `iron` axes down in proportion to contact fraction.
/// No-op when `cfg.is_identity()`. RNG-free and geometric, so byte-identity of
/// the cell grid is preserved when the layer is off.
pub fn apply_contact_resistance_3d(grid: &mut TumorGrid3D, cfg: &ContactConfig) {
    if cfg.is_identity() {
        return;
    }
    // `contact_fraction_3d` reads `is_tumor`, which this function never mutates
    // (only `iron`/`lipid_unsat`), so computing each fraction inside the same
    // pass is order-independent.
    let n = grid.cells.len();
    for idx in 0..n {
        if !grid.cells[idx].is_tumor {
            continue;
        }
        let frac = contact_fraction_3d(grid, idx);
        let cell = &mut grid.cells[idx].cell;
        cell.lipid_unsat *= (1.0 - cfg.lipid_strength * frac).max(0.0);
        cell.iron *= (1.0 - cfg.iron_strength * frac).max(0.0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn grid() -> TumorGrid3D {
        TumorGrid3D::generate(40, 40, 40, 20.0, 42)
    }

    #[test]
    fn identity_is_a_no_op() {
        let g = grid();
        let mut g2 = grid();
        apply_contact_resistance_3d(&mut g2, &ContactConfig::default());
        // Byte-identical iron/lipid_unsat for every cell under the identity config.
        for (a, b) in g.cells.iter().zip(g2.cells.iter()) {
            assert_eq!(a.cell.iron, b.cell.iron);
            assert_eq!(a.cell.lipid_unsat, b.cell.lipid_unsat);
        }
        assert!(ContactConfig::default().is_identity());
        assert!(!ContactConfig::literature().is_identity());
    }

    #[test]
    fn interior_cells_are_more_contacted_than_surface() {
        let g = grid();
        // Deep-interior cell (grid centre) is fully surrounded by tumor.
        let center = g.flat_index(20, 20, 20);
        assert!(g.cells[center].is_tumor);
        let c_center = contact_fraction_3d(&g, center);
        assert!(
            c_center > 0.95,
            "interior contact {c_center} should be ~1.0"
        );
        // A tumor cell near the spheroid surface has at least one non-tumor
        // neighbour ⇒ strictly lower contact than the saturated interior.
        let surface = g
            .cells
            .iter()
            .enumerate()
            .filter(|(_, gc)| gc.is_tumor)
            .map(|(i, _)| (i, contact_fraction_3d(&g, i)))
            .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
            .unwrap();
        assert!(
            surface.1 < c_center,
            "surface contact {} should be < interior {c_center}",
            surface.1
        );
        assert!((0.0..=1.0).contains(&surface.1));
    }

    #[test]
    fn contact_reduces_lipid_and_iron_monotonically() {
        let baseline = grid();
        let mut g = grid();
        let cfg = ContactConfig::literature();
        apply_contact_resistance_3d(&mut g, &cfg);
        // A fully-contacted interior cell is reduced by ~the configured strength;
        // every tumor cell's axes are reduced (or unchanged at contact 0), never
        // increased; non-tumor cells are untouched.
        let center = g.flat_index(20, 20, 20);
        let frac = contact_fraction_3d(&baseline, center);
        let exp_lipid = baseline.cells[center].cell.lipid_unsat * (1.0 - cfg.lipid_strength * frac);
        assert!((g.cells[center].cell.lipid_unsat - exp_lipid).abs() < 1e-12);
        for (b, a) in baseline.cells.iter().zip(g.cells.iter()) {
            if b.is_tumor {
                assert!(a.cell.lipid_unsat <= b.cell.lipid_unsat + 1e-12);
                assert!(a.cell.iron <= b.cell.iron + 1e-12);
                assert!(a.cell.lipid_unsat > 0.0 && a.cell.iron >= 0.0);
            } else {
                assert_eq!(a.cell.lipid_unsat, b.cell.lipid_unsat);
                assert_eq!(a.cell.iron, b.cell.iron);
            }
        }
    }
}
