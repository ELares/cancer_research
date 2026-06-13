//! IFN-gamma -> System Xc- ferroptosis-sensitization coupling (#443).
//!
//! Activated CD8+ T cells secrete IFN-gamma, which transcriptionally downregulates
//! the cystine/glutamate antiporter System Xc- (SLC7A11/SLC3A2) on nearby tumor
//! cells. Less cystine import means less glutathione (GSH), so GPX4 loses its
//! reducing substrate and the cells become MORE sensitive to ferroptosis. This is
//! the molecular feedback loop that couples immune activation to ferroptosis, the
//! mechanistic link under the manuscript's immune-amplification thesis (Wang et al.,
//! Nature 2019, PMID 31043744, "CD8+ T cells regulate tumour ferroptosis ... IFN-gamma
//! downregulates SLC3A2 and SLC7A11"). IFN-gamma also raises ACSL4 (more oxidizable
//! PUFA), a second, parallel sensitizing arm.
//!
//! ## What this module provides
//!
//! Pure, off-by-default coupling functions and a config; the spatial wiring (seed an
//! IFN-gamma field at immune-active positions, diffuse it, apply it per cell) lives
//! in the consumer (sim-tme-3d). Two arms are provided, but only the dominant one is
//! currently wired into the spatial consumer:
//!   - [`system_xc_retention`]: the GSH-retention multiplier `ic50/(ifn+ic50)` in
//!     `(0, 1]`, 1 at `ifn == 0` (no IFN-gamma, full cystine uptake) and falling
//!     toward 0 as IFN-gamma rises (System Xc- shut down). A consumer multiplies the
//!     cell GSH pool by this, the DOMINANT Wang-2019 arm. **This is the arm wired
//!     into sim-tme-3d.**
//!   - [`acsl4_upregulation`]: the lipid-unsaturation boost multiplier `>= 1`, 1 at
//!     `ifn == 0`, rising with IFN-gamma (more ACSL4 -> more PUFA). The parallel arm.
//!     **Provided and unit-tested in isolation, but NOT yet wired into the per-cell
//!     static-lipid axis in sim-tme-3d (a deferred follow-up); see
//!     `simulations/calibration/CALIBRATION_STATUS.md`.**
//!
//! ## Honesty / calibration
//!
//! The DIRECTION (IFN-gamma sensitizes to ferroptosis via System Xc- suppression) is
//! a verified landmark (Wang 2019). The magnitudes here (`system_xc_ic50`,
//! `acsl4_strength`, the field seeding/diffusion constants) are UNCALIBRATED
//! placeholders pending IFN-gamma + ferroptosis co-culture data; the result is the
//! coupling direction, not the numbers. [`IFNGammaConfig::disabled`] (the default) is
//! identity (`retention == 1`, `acsl4 == 1`) so a consumer that does not opt in is
//! byte-identical.

/// IFN-gamma coupling configuration. [`disabled`](IFNGammaConfig::disabled) is the
/// identity (no coupling, byte-identical); [`literature`](IFNGammaConfig::literature)
/// is the uncalibrated placeholder that turns the coupling on.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IFNGammaConfig {
    /// IFN-gamma concentration (model units) at which System Xc- is half-inhibited
    /// (the GSH-retention multiplier hits 0.5). Larger ⇒ less sensitive to IFN-gamma.
    pub system_xc_ic50: f64,
    /// ACSL4-upregulation strength: the lipid-unsaturation boost saturates toward
    /// `1 + acsl4_strength` as IFN-gamma grows.
    pub acsl4_strength: f64,
    /// Seeding stoichiometry: a consumer seeds the IFN-gamma field as
    /// `per_damp * local_damp` at immune-active (DAMP-positive) positions, coupling
    /// IFN-gamma secretion to the local immune response.
    pub per_damp: f64,
    /// Diffusion fraction for the IFN-gamma field (consumer passes this to the shared
    /// `diffuse_damp_3d_step` operator).
    pub diffusion_fraction: f64,
    /// Per-step clearance rate for the IFN-gamma field.
    pub clearance_rate: f64,
}

impl IFNGammaConfig {
    /// Identity: no coupling. `system_xc_ic50 = INFINITY` so retention is always 1
    /// and the ACSL4 boost is 1; a consumer using this is byte-identical.
    pub fn disabled() -> Self {
        IFNGammaConfig {
            system_xc_ic50: f64::INFINITY,
            acsl4_strength: 0.0,
            per_damp: 0.0,
            diffusion_fraction: 0.0,
            clearance_rate: 0.0,
        }
    }

    /// Uncalibrated placeholder that enables the coupling. The IC50 and ACSL4
    /// strength set the sensitization magnitude; the field constants mirror the DAMP
    /// field's transport (a slightly slower-diffusing, similarly-cleared cytokine).
    /// All magnitudes are placeholders; the direction is the result.
    pub fn literature() -> Self {
        IFNGammaConfig {
            system_xc_ic50: 1.0,
            acsl4_strength: 0.3,
            per_damp: 0.3,
            diffusion_fraction: 0.02,
            clearance_rate: 0.03,
        }
    }

    /// True when the coupling applies no effect (the `disabled()` identity): a
    /// consumer can skip the whole IFN-gamma path and stay byte-identical.
    pub fn is_disabled(&self) -> bool {
        self.system_xc_ic50.is_infinite() && self.acsl4_strength == 0.0 && self.per_damp == 0.0
    }
}

/// GSH-retention multiplier under local IFN-gamma: `ic50 / (ifn + ic50)`, in `(0, 1]`.
/// `1.0` at `ifn == 0` (full cystine uptake) and falling toward `0` as IFN-gamma
/// rises (System Xc- shut down ⇒ cystine starvation ⇒ GSH drop). A consumer
/// multiplies the cell GSH pool by this. Negative `ifn` is clamped to 0.
pub fn system_xc_retention(local_ifn: f64, ic50: f64) -> f64 {
    if !ic50.is_finite() {
        return 1.0;
    }
    let ifn = local_ifn.max(0.0);
    (ic50 / (ifn + ic50)).clamp(0.0, 1.0)
}

/// ACSL4-upregulation multiplier for lipid unsaturation under local IFN-gamma:
/// `1 + strength * ifn/(ifn + 1)`, `>= 1`, `1.0` at `ifn == 0` and saturating toward
/// `1 + strength`. A consumer multiplies the cell PUFA / lipid-unsaturation by this.
pub fn acsl4_upregulation(local_ifn: f64, strength: f64) -> f64 {
    let ifn = local_ifn.max(0.0);
    1.0 + strength.max(0.0) * (ifn / (ifn + 1.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn disabled_is_identity() {
        let c = IFNGammaConfig::disabled();
        assert!(c.is_disabled());
        // Any IFN-gamma value leaves GSH retention at 1 and ACSL4 boost at 1.
        for ifn in [0.0, 0.5, 2.0, 100.0] {
            assert_eq!(system_xc_retention(ifn, c.system_xc_ic50), 1.0);
            assert_eq!(acsl4_upregulation(ifn, c.acsl4_strength), 1.0);
        }
        assert!(!IFNGammaConfig::literature().is_disabled());
    }

    #[test]
    fn system_xc_retention_falls_monotonically_with_ifn() {
        let ic50 = 1.0;
        assert_eq!(system_xc_retention(0.0, ic50), 1.0); // no IFN -> full retention
        let r1 = system_xc_retention(0.5, ic50);
        let r2 = system_xc_retention(1.0, ic50);
        let r3 = system_xc_retention(5.0, ic50);
        assert!(
            1.0 > r1 && r1 > r2 && r2 > r3,
            "monotone decreasing: 1 {r1} {r2} {r3}"
        );
        assert!(
            (r2 - 0.5).abs() < 1e-9,
            "retention is 0.5 at ifn == ic50: {r2}"
        );
        assert!(r3 > 0.0 && r3 < 0.2, "retention -> 0 at high IFN: {r3}");
        // negative IFN clamps to 0 (full retention)
        assert_eq!(system_xc_retention(-3.0, ic50), 1.0);
    }

    #[test]
    fn acsl4_upregulation_rises_monotonically_and_is_bounded() {
        let s = 0.3;
        assert_eq!(acsl4_upregulation(0.0, s), 1.0); // no IFN -> no boost
        let a1 = acsl4_upregulation(0.5, s);
        let a2 = acsl4_upregulation(2.0, s);
        let a3 = acsl4_upregulation(1000.0, s);
        assert!(
            1.0 < a1 && a1 < a2 && a2 < a3,
            "monotone increasing: 1 {a1} {a2} {a3}"
        );
        assert!(
            a3 <= 1.0 + s + 1e-9 && a3 > 1.0 + 0.9 * s,
            "saturates toward 1+strength: {a3}"
        );
    }

    #[test]
    fn arms_sensitize_to_ferroptosis_in_the_right_direction() {
        // The coupling should LOWER GSH (retention < 1) and RAISE PUFA (boost > 1)
        // under IFN-gamma, both ferroptosis-sensitizing, the Wang 2019 direction.
        let c = IFNGammaConfig::literature();
        let ifn = 2.0;
        assert!(
            system_xc_retention(ifn, c.system_xc_ic50) < 1.0,
            "GSH retention drops"
        );
        assert!(
            acsl4_upregulation(ifn, c.acsl4_strength) > 1.0,
            "PUFA rises"
        );
    }
}
