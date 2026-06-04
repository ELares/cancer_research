//! Radial nutrient gradient (#270 item 3b, follow-up to #197).
//!
//! Beyond the O2 (#187) and pH (#190) radial fields, a spheroid develops a
//! radial NUTRIENT gradient: glucose and glutamine are abundant at the
//! well-perfused rim and consumed toward the core. Glucose metabolism feeds the
//! pentose-phosphate pathway that supplies the NADPH used to regenerate reduced
//! glutathione for the GPX4 antioxidant defense, so a nutrient-starved core has
//! LESS antioxidant capacity and is more ferroptosis-sensitive on that axis
//! (Dixon et al. 2012, PMID 22632970, the foundational GSH-dependent
//! lipid-peroxidation mechanism; glucose metabolic reprogramming regulates
//! ferroptosis, PMID 42190602).
//!
//! ## Direction caveat (read the DIRECTION, not the number)
//!
//! This models ONE documented direction: nutrient (glucose) deprivation lowers
//! antioxidant capacity (the issue #270 framing). It is NOT the whole story:
//! energy stress also activates AMPK, which can INHIBIT ferroptosis by cutting
//! PUFA synthesis, and glutaminolysis is REQUIRED for some ferroptosis (Gao et
//! al., "Role of Mitochondria in Ferroptosis," PMID 30581146). So the NET effect
//! of a nutrient-starved core is genuinely context-dependent. We model the
//! antioxidant axis with an UNCALIBRATED placeholder strength; compose with the
//! spheroid / O2 / pH layers to explore the net. The magnitude is not a result;
//! the direction (less glucose ⇒ less NADPH ⇒ less GSH/GPX4 regeneration) is.
//!
//! ## Model
//!
//! Nutrient availability shares the O2 field's radial-supply form (a different
//! penetration length `λ`, a different coupling):
//!
//! ```text
//! availability(depth) = exp(-depth / λ)      ∈ (0, 1]   (1 at the rim, →0 core)
//! deprivation         = 1 - availability
//! cell.nrf2 *= max(0, 1 - antioxidant_strength · deprivation)
//! ```
//!
//! `cell.nrf2` is the durable antioxidant setpoint (it scales both the GPX4
//! target and GSH resynthesis, #266), so reducing it is a durable, whole-run
//! antioxidant-capacity cut, not just an initial condition. Geometric (reads
//! `is_tumor` + radial depth via [`RadialDepthGeom`], no RNG), so off-by-default
//! identity ([`NutrientConfig::default`]) keeps the matrix byte-identical.
//! Composes multiplicatively with the clonal `gpx4_mul` antioxidant axis.

use crate::grid::{RadialDepthGeom, TumorGrid3D};

/// Nutrient-gradient configuration. `antioxidant_strength` defaults to 0
/// (identity ⇒ byte-identical). `lambda_um` is the nutrient penetration length
/// (same units/role as the O2 `λ`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NutrientConfig {
    /// Nutrient penetration length (µm). Larger ⇒ nutrients reach deeper.
    pub lambda_um: f64,
    /// Fraction by which a fully nutrient-deprived (availability → 0) cell's
    /// durable antioxidant setpoint `cell.nrf2` is reduced. `0.0` ⇒ identity.
    pub antioxidant_strength: f64,
}

impl Default for NutrientConfig {
    fn default() -> Self {
        // Identity: no antioxidant modulation ⇒ byte-identical.
        NutrientConfig {
            lambda_um: 120.0,
            antioxidant_strength: 0.0,
        }
    }
}

impl NutrientConfig {
    /// Literature-motivated qualitative strength: a deep, nutrient-starved core
    /// loses a meaningful slice of its NADPH-fed antioxidant capacity. The
    /// `λ = 120 µm` matches the O2 zone-reference length; `antioxidant_strength`
    /// is an **UNCALIBRATED placeholder** encoding the documented direction
    /// (less glucose ⇒ less GSH/GPX4 regeneration), not fit to data. Calibrate
    /// vs depth-resolved glucose/NADPH/GSH measurements in spheroids.
    pub fn literature() -> Self {
        NutrientConfig {
            lambda_um: 120.0,
            antioxidant_strength: 0.3,
        }
    }

    /// True when the config applies no modulation (`antioxidant_strength == 0`).
    pub fn is_identity(&self) -> bool {
        self.antioxidant_strength == 0.0
    }
}

/// Radial nutrient availability at `depth_um` (µm from the spheroid surface,
/// negative outside ⇒ clamped to 0): `exp(-depth / λ) ∈ (0, 1]`. Shares the O2
/// field's form. Non-finite/`λ ≤ 0` is a caller bug (`debug_assert`); release
/// returns a clamped value.
pub fn nutrient_availability(depth_um: f64, lambda_um: f64) -> f64 {
    debug_assert!(
        lambda_um.is_finite() && lambda_um > 0.0,
        "nutrient_availability: lambda_um must be finite and positive, got {lambda_um}"
    );
    (-depth_um.max(0.0) / lambda_um).exp().clamp(0.0, 1.0)
}

/// Apply nutrient-gradient antioxidant stress to every tumor cell (#270 item 3b):
/// scale the durable antioxidant setpoint `cell.nrf2` down in proportion to
/// radial nutrient deprivation. No-op when `cfg.is_identity()`. RNG-free and
/// geometric, so the cell grid stays byte-identical when the layer is off.
pub fn apply_nutrient_stress_3d(grid: &mut TumorGrid3D, cfg: &NutrientConfig) {
    if cfg.is_identity() {
        return;
    }
    // The depth geometry is dimension-only; hoist it once (#289).
    let geom = RadialDepthGeom::new(grid);
    let n = grid.cells.len();
    for idx in 0..n {
        if !grid.cells[idx].is_tumor {
            continue;
        }
        let (r, c, l) = grid.coords(idx);
        let depth_um = geom.depth_um(r, c, l).max(0.0);
        let availability = nutrient_availability(depth_um, cfg.lambda_um);
        let deprivation = 1.0 - availability;
        let cell = &mut grid.cells[idx].cell;
        cell.nrf2 *= (1.0 - cfg.antioxidant_strength * deprivation).max(0.0);
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
        apply_nutrient_stress_3d(&mut g2, &NutrientConfig::default());
        for (a, b) in g.cells.iter().zip(g2.cells.iter()) {
            assert_eq!(a.cell.nrf2, b.cell.nrf2);
        }
        assert!(NutrientConfig::default().is_identity());
        assert!(!NutrientConfig::literature().is_identity());
    }

    #[test]
    fn availability_is_one_at_surface_and_decays_inward() {
        // At the surface (depth 0) availability is 1 (no deprivation).
        assert!((nutrient_availability(0.0, 120.0) - 1.0).abs() < 1e-12);
        // Negative depth (outside the tumor) clamps to surface (1.0).
        assert!((nutrient_availability(-50.0, 120.0) - 1.0).abs() < 1e-12);
        // Deeper ⇒ strictly less available, bounded in (0, 1].
        let shallow = nutrient_availability(60.0, 120.0);
        let deep = nutrient_availability(300.0, 120.0);
        assert!(deep < shallow && shallow < 1.0);
        assert!(deep > 0.0);
    }

    #[test]
    fn core_loses_more_antioxidant_capacity_than_rim() {
        let baseline = grid();
        let mut g = grid();
        let cfg = NutrientConfig::literature();
        apply_nutrient_stress_3d(&mut g, &cfg);
        let geom = RadialDepthGeom::new(&baseline);

        // A deep-core tumor cell (grid centre) is more nutrient-deprived than a
        // near-surface tumor cell, so its nrf2 is reduced more (lower ratio).
        let center = baseline.flat_index(20, 20, 20);
        assert!(baseline.cells[center].is_tumor);

        // Every tumor cell's nrf2 is reduced (or unchanged at the rim), never
        // increased; non-tumor cells are untouched; the reduction equals the
        // geometric value.
        let mut any_reduced = false;
        for (i, (b, a)) in baseline.cells.iter().zip(g.cells.iter()).enumerate() {
            if b.is_tumor {
                assert!(a.cell.nrf2 <= b.cell.nrf2 + 1e-12);
                let (r, c, l) = baseline.coords(i);
                let depth = geom.depth_um(r, c, l).max(0.0);
                let dep = 1.0 - nutrient_availability(depth, cfg.lambda_um);
                let expected = b.cell.nrf2 * (1.0 - cfg.antioxidant_strength * dep);
                assert!((a.cell.nrf2 - expected).abs() < 1e-12);
                if dep > 0.0 {
                    any_reduced = true;
                }
            } else {
                assert_eq!(a.cell.nrf2, b.cell.nrf2);
            }
        }
        assert!(
            any_reduced,
            "interior tumor cells must lose antioxidant capacity"
        );

        // The centre (deepest) cell is reduced at least as much as a rim cell.
        let rim = baseline
            .cells
            .iter()
            .enumerate()
            .filter(|(_, gc)| gc.is_tumor)
            .map(|(i, _)| {
                let (r, c, l) = baseline.coords(i);
                (i, geom.depth_um(r, c, l).max(0.0))
            })
            .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
            .unwrap()
            .0;
        let center_ratio = g.cells[center].cell.nrf2 / baseline.cells[center].cell.nrf2;
        let rim_ratio = g.cells[rim].cell.nrf2 / baseline.cells[rim].cell.nrf2;
        assert!(
            center_ratio <= rim_ratio + 1e-12,
            "core nrf2 ratio {center_ratio} should be <= rim ratio {rim_ratio}"
        );
    }
}
