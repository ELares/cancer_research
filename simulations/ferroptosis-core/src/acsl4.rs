//! ACSL4-status biomarker stratification (#444).
//!
//! ACSL4 (acyl-CoA synthetase long-chain family member 4) is the enzyme that
//! ligates the polyunsaturated fatty acids (arachidonic / adrenic acid) ferroptosis
//! requires into membrane phospholipids. It is the single most discriminating
//! pro-ferroptotic lipid-metabolism gene: ACSL4-high tumors load their membranes
//! with the oxidizable PUFA-PE substrate and are ferroptosis-SENSITIVE, while
//! ACSL4-low / ACSL4-negative tumors incorporate little PUFA and are intrinsically
//! ferroptosis-REFRACTORY through a mechanism entirely distinct from the
//! GPX4/GSH/FSP1 antioxidant defenses the rest of the model carries (Doll et al.,
//! Nat Chem Biol 2017, PMID 27842070, ACSL4 dictates ferroptosis sensitivity;
//! Yang et al., PNAS 2016, PMID 27506793). Several tumor contexts are
//! constitutively ACSL4-low/negative (e.g. some hepatocellular carcinoma and AML
//! subtypes), which is why a fixed, uniform PUFA assumption mis-predicts their
//! ferroptosis response.
//!
//! ## What this module provides
//!
//! A single biomarker axis: an ACSL4 *status* (a relative expression scalar,
//! `1.0` = wild-type / model baseline) is mapped to an additive PUFA-incorporation
//! boost a consumer writes onto [`crate::params::Params::acsl4_status_boost`], which
//! the biochem engine folds into the oxidizable-PUFA augmentation
//! (`biochem::ether_augmented_pufa`). The boost is `status - 1`, clamped at `-1`
//! (the ACSL4-null floor: no PUFA incorporation):
//!   - [`ACSL4_NEGATIVE`] (0.0) → boost `-1` → the PUFA substrate collapses ⇒
//!     ferroptosis-refractory (the distinct, defense-independent escape).
//!   - [`ACSL4_LOW`] (0.5) → boost `-0.5` ⇒ partially resistant.
//!   - [`ACSL4_NORMAL`] (1.0) → boost `0` ⇒ the model baseline (byte-identical).
//!   - [`ACSL4_HIGH`] (1.5) → boost `+0.5` ⇒ more oxidizable PUFA ⇒ sensitive.
//!
//! This composes additively with the MCFA→ACSL4 boost (#446) and the ether-lipid
//! pool (#339) in the same augmentation: ACSL4 STATUS is the tumor-intrinsic
//! expression baseline, MCFA is a pharmacological perturbation that acts through
//! ACSL4, and ether-PUFA is a separate oxidizable pool.
//!
//! ## Honesty / calibration
//!
//! The DIRECTION (ACSL4-high ⇒ ferroptosis-sensitive; ACSL4-negative ⇒ refractory)
//! is a verified landmark (Doll 2017). The numeric status→boost mapping is an
//! UNCALIBRATED linear placeholder; the per-cancer-type ACSL4 prevalence and the
//! cell-line ACSL4-status-vs-dose-response validation are DATA-GATED (they need
//! TCGA/cBioPortal expression + a cell-line meta-analysis the repo does not carry)
//! and are NOT fabricated here. `status == 1.0` (the [`ACSL4_NORMAL`] baseline)
//! yields boost `0` exactly, so a consumer that does not opt in is byte-identical.

/// ACSL4-negative status (absent/deleted ACSL4): boost `-1` ⇒ PUFA substrate
/// collapses ⇒ intrinsically ferroptosis-refractory. The escape mechanism
/// distinct from GPX4/GSH/FSP1 (e.g. some HCC / AML subtypes).
pub const ACSL4_NEGATIVE: f64 = 0.0;
/// ACSL4-low status: boost `-0.5` ⇒ partially ferroptosis-resistant.
pub const ACSL4_LOW: f64 = 0.5;
/// ACSL4-normal / wild-type status: boost `0` ⇒ the model baseline
/// (byte-identical).
pub const ACSL4_NORMAL: f64 = 1.0;
/// ACSL4-high status (placeholder): boost `+0.5` ⇒ more oxidizable PUFA ⇒
/// ferroptosis-sensitive (e.g. lung, ER+ breast, cervical contexts).
pub const ACSL4_HIGH: f64 = 1.5;

/// Map an ACSL4 expression status (`>= 0`, `1.0` = wild-type) to the additive
/// PUFA-incorporation boost written to `Params::acsl4_status_boost`:
/// `status - 1`, clamped at `-1` (the ACSL4-null floor). `1.0` ⇒ `0` exactly
/// (byte-identical baseline); `< 1` ⇒ negative (less PUFA, resistant); `> 1` ⇒
/// positive (more PUFA, sensitive). Negative `status` is clamped to `0`.
pub fn pufa_boost_from_status(status: f64) -> f64 {
    (status.max(0.0) - 1.0).max(-1.0)
}

/// Map a within-cohort ACSL4 mRNA z-score to the ACSL4 [`status`](pufa_boost_from_status)
/// scalar: `max(0, 1 + z/2)`. This is the calibrated bridge from real expression
/// data to the model input (#462), anchored to cBioPortal TCGA PanCancer Atlas
/// per-cancer-type ACSL4 z-score distributions
/// (`analysis/calibration/acsl4-prevalence-calibration.md`). The shipped status
/// constants are exactly the integer-z points of this bridge, which is why they
/// were well chosen as placeholders and are now interpretable as z-scores:
///   - `z = +1` ⇒ [`ACSL4_HIGH`] (1.5)
///   - `z =  0` ⇒ [`ACSL4_NORMAL`] (1.0, the wild-type baseline)
///   - `z = -1` ⇒ [`ACSL4_LOW`] (0.5)
///   - `z = -2` ⇒ [`ACSL4_NEGATIVE`] (0.0, the PUFA-collapse floor)
///
/// A consumer that has a patient's ACSL4 mRNA z-score can feed
/// `pufa_boost_from_status(status_from_zscore(z))` to stratify that tumor. The
/// per-cancer-type fraction of tumors with `z < -1` (about 11 to 19% in TCGA,
/// median 14%) is the committed population prior for how many tumors fall in the
/// low-ACSL4 (refractory-leaning) tail. The slope (`/2`) is the placeholder that
/// reproduces the existing constants; only the within-cohort z interpretation is
/// data-anchored, not the absolute status→ferroptosis magnitude.
pub fn status_from_zscore(z: f64) -> f64 {
    (1.0 + z / 2.0).max(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normal_status_is_zero_boost() {
        // Wild-type ACSL4 is the model baseline ⇒ boost exactly 0 ⇒ byte-identical.
        assert_eq!(pufa_boost_from_status(ACSL4_NORMAL), 0.0);
    }

    #[test]
    fn status_orders_ferroptosis_sensitivity() {
        // Boost is monotone increasing in ACSL4 status: negative < low < normal <
        // high, the Doll-2017 direction (more ACSL4 ⇒ more oxidizable PUFA ⇒ more
        // ferroptosis).
        let neg = pufa_boost_from_status(ACSL4_NEGATIVE);
        let low = pufa_boost_from_status(ACSL4_LOW);
        let normal = pufa_boost_from_status(ACSL4_NORMAL);
        let high = pufa_boost_from_status(ACSL4_HIGH);
        assert!(
            neg < low && low < normal && normal < high,
            "{neg} {low} {normal} {high}"
        );
        // ACSL4-negative collapses the PUFA term (boost -1, the null floor).
        assert_eq!(neg, -1.0);
        assert_eq!(high, 0.5);
    }

    #[test]
    fn null_floor_and_negative_clamp() {
        // Even an over-deleted / negative status cannot drive the boost below -1
        // (the PUFA term floors at 0 in the engine, not negative).
        assert_eq!(pufa_boost_from_status(-3.0), -1.0);
        assert_eq!(pufa_boost_from_status(0.0), -1.0);
        // A very high ACSL4 keeps rising (no upper clamp; the engine consumes it).
        assert!(pufa_boost_from_status(3.0) > pufa_boost_from_status(2.0));
    }

    #[test]
    fn zscore_bridge_hits_the_status_constants() {
        // The calibrated z-score bridge (#462) reproduces the shipped status
        // constants exactly at integer z, so a real cBioPortal ACSL4 z-score maps
        // onto the existing NEGATIVE/LOW/NORMAL/HIGH scale.
        assert_eq!(status_from_zscore(1.0), ACSL4_HIGH);
        assert_eq!(status_from_zscore(0.0), ACSL4_NORMAL);
        assert_eq!(status_from_zscore(-1.0), ACSL4_LOW);
        assert_eq!(status_from_zscore(-2.0), ACSL4_NEGATIVE);
        // Floors at 0 below z = -2 (cannot go below the ACSL4-null collapse).
        assert_eq!(status_from_zscore(-5.0), 0.0);
        // Monotone increasing in z.
        assert!(status_from_zscore(2.0) > status_from_zscore(1.0));
        // Composing with the boost: a z = -2 tumor collapses the PUFA term (boost -1).
        assert_eq!(pufa_boost_from_status(status_from_zscore(-2.0)), -1.0);
    }
}
