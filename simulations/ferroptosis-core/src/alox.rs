//! ALOX (lipoxygenase) isoform-specific lipid-peroxidation activity + MCFA
//! sensitization (#446).
//!
//! The core engine drives lipid peroxidation through one generic propagation
//! rate (`Params::lp_propagation`), implicitly assuming a fixed, "average"
//! enzymatic-oxidation capacity. In reality the rate at which membrane PUFA is
//! oxidized to the ferroptotic lipid hydroperoxides depends on which
//! arachidonate-lipoxygenase isoforms a tumor expresses and how active they are.
//! ALOX15 (15-LOX) is the canonical ferroptosis-driving isoform (it oxidizes
//! PUFA-PE directly), with ALOX12 and ALOX5 contributing at isoform-specific
//! rates that vary several-fold (Yang & Stockwell; lipoxygenase-driven
//! ferroptosis reviewed in PNAS 2016, PMID 27506793). So an ALOX15/12/5-high
//! tumor peroxidizes faster (more ferroptosis-sensitive) and an ALOX-poor tumor
//! peroxidizes slower (more resistant), independent of the GPX4/GSH/FSP1 defense
//! axes the model already carries.
//!
//! A second, parallel sensitization axis is medium-chain fatty acids (MCFA):
//! MCFA exposure upregulates ACSL4 / CD36, raising the incorporation of
//! oxidizable PUFA into membranes and thereby raising ferroptosis susceptibility
//! (an emerging pharmacological strategy; Sci Rep 2024 s41598-024-55050-4;
//! medium-chain-fatty-acid ferroptosis sensitization PMC11901882).
//!
//! ## What this module provides
//!
//! [`AloxConfig`] collapses an isoform activity/expression mix into two
//! off-by-default, additive *boosts* a consumer writes onto [`crate::params::Params`]:
//!   - [`AloxConfig::lp_propagation_boost`]: added to the propagation multiplier
//!     as `1 + boost` (so `0` ⇒ unchanged; `>0` ALOX-high ⇒ faster propagation;
//!     down to `-1` ⇒ the ALOX-null/knockout limit with no enzymatic
//!     propagation). Written to `Params::alox_propagation_boost`.
//!   - [`AloxConfig::mcfa_pufa_boost`]: added to the oxidizable-PUFA augmentation
//!     (alongside the ether-lipid pool), `0` ⇒ unchanged, rising with MCFA.
//!     Written to `Params::mcfa_pufa_boost`.
//!
//! ## Honesty / calibration
//!
//! The DIRECTIONS (ALOX-high ⇒ more ferroptosis; MCFA ⇒ more ferroptosis) are
//! literature-anchored, but the per-isoform activity weights and the MCFA
//! kinetics are UNCALIBRATED placeholders (isoform Kcat/Km vary ~10-fold and are
//! not fit here; absolute MCFA kinetics are deferred to the experimental
//! E-series). [`AloxConfig::identity`] (the implicit default of the existing
//! model: balanced isoforms, no MCFA) yields both boosts `== 0` exactly, so a
//! consumer that does not opt in is byte-identical. Per-cell stochastic ALOX
//! heterogeneity (per-phenotype sampling like iron/gsh) is a deferred refinement;
//! this models a per-condition ALOX phenotype.

/// ALOX isoform mix + MCFA level. [`identity`](AloxConfig::identity) is the
/// balanced, no-MCFA baseline (both boosts `0`, byte-identical);
/// [`literature`](AloxConfig::literature) is the uncalibrated ALOX15-high + MCFA
/// placeholder that turns the sensitization on.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AloxConfig {
    /// Relative PUFA-oxidation activity of ALOX15 (15-LOX), the canonical
    /// ferroptosis-driving isoform. `1.0` is the model's implicit average.
    pub alox15_activity: f64,
    /// Relative PUFA-oxidation activity of ALOX12 (12-LOX).
    pub alox12_activity: f64,
    /// Relative PUFA-oxidation activity of ALOX5 (5-LOX).
    pub alox5_activity: f64,
    /// Expression fraction of ALOX15 (need not sum to 1 with the others; the
    /// effective activity is the fraction-weighted mean of the three).
    pub alox15_frac: f64,
    /// Expression fraction of ALOX12.
    pub alox12_frac: f64,
    /// Expression fraction of ALOX5.
    pub alox5_frac: f64,
    /// The fraction-weighted activity that maps to NO change (propagation boost
    /// `0`). The balanced [`identity`](AloxConfig::identity) baseline equals this
    /// exactly, so its boost is `0` and the run is byte-identical.
    pub reference_activity: f64,
    /// Medium-chain fatty acid level (arbitrary units, `0` = none). MCFA
    /// upregulates ACSL4/CD36 ⇒ more oxidizable PUFA incorporation.
    pub mcfa_level: f64,
    /// MCFA → PUFA sensitization strength: the PUFA boost saturates toward this
    /// as `mcfa_level` grows.
    pub mcfa_strength: f64,
}

impl AloxConfig {
    /// Identity: a balanced isoform mix (all activities at the reference) and no
    /// MCFA. Both boosts are exactly `0`, so a consumer using this is
    /// byte-identical to the model without an ALOX layer.
    pub fn identity() -> Self {
        AloxConfig {
            alox15_activity: 1.0,
            alox12_activity: 1.0,
            alox5_activity: 1.0,
            alox15_frac: 1.0,
            alox12_frac: 1.0,
            alox5_frac: 1.0,
            reference_activity: 1.0,
            mcfa_level: 0.0,
            mcfa_strength: 0.0,
        }
    }

    /// Uncalibrated placeholder: an ALOX15-high tumor (ALOX15 the most active and
    /// most expressed) plus a moderate MCFA exposure, so BOTH the propagation
    /// boost and the PUFA boost are positive (more ferroptosis). All magnitudes
    /// are placeholders; the directions are the result.
    pub fn literature() -> Self {
        AloxConfig {
            alox15_activity: 1.6,
            alox12_activity: 1.1,
            alox5_activity: 0.8,
            alox15_frac: 0.60,
            alox12_frac: 0.25,
            alox5_frac: 0.15,
            reference_activity: 1.0,
            mcfa_level: 1.0,
            mcfa_strength: 0.5,
        }
    }

    /// Fraction-weighted mean isoform activity. Falls back to the reference when
    /// no isoform is expressed (all fractions `0`), so the boost is `0` there.
    pub fn weighted_activity(&self) -> f64 {
        let f15 = self.alox15_frac.max(0.0);
        let f12 = self.alox12_frac.max(0.0);
        let f5 = self.alox5_frac.max(0.0);
        let total = f15 + f12 + f5;
        if total <= 0.0 {
            return self.reference_activity;
        }
        (self.alox15_activity * f15 + self.alox12_activity * f12 + self.alox5_activity * f5) / total
    }

    /// Additive boost for the propagation multiplier: `weighted_activity /
    /// reference_activity - 1`. `0` at the balanced baseline; `>0` ALOX-high
    /// (faster propagation); `< 0` ALOX-poor, floored at `-1` (the ALOX-null
    /// limit ⇒ propagation multiplier `0`). A consumer writes this to
    /// `Params::alox_propagation_boost`.
    pub fn lp_propagation_boost(&self) -> f64 {
        let reference = self.reference_activity.max(1e-9);
        (self.weighted_activity() / reference - 1.0).max(-1.0)
    }

    /// Additive boost for the oxidizable-PUFA augmentation under MCFA:
    /// `mcfa_strength * mcfa_level / (mcfa_level + 1)`, `>= 0`, `0` at
    /// `mcfa_level == 0` or `mcfa_strength == 0`, saturating toward
    /// `mcfa_strength`. A consumer writes this to `Params::mcfa_pufa_boost`.
    pub fn mcfa_pufa_boost(&self) -> f64 {
        let level = self.mcfa_level.max(0.0);
        self.mcfa_strength.max(0.0) * (level / (level + 1.0))
    }

    /// True when neither boost changes anything (the [`identity`](AloxConfig::identity)
    /// case): a consumer can skip the ALOX path and stay byte-identical.
    pub fn is_identity(&self) -> bool {
        self.lp_propagation_boost() == 0.0 && self.mcfa_pufa_boost() == 0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identity_is_zero_boost() {
        let c = AloxConfig::identity();
        assert_eq!(c.lp_propagation_boost(), 0.0);
        assert_eq!(c.mcfa_pufa_boost(), 0.0);
        assert!(c.is_identity());
        assert!(!AloxConfig::literature().is_identity());
    }

    #[test]
    fn alox15_high_raises_propagation_boost_low_lowers_it() {
        // ALOX15-high (most active isoform dominant) ⇒ positive boost.
        let high = AloxConfig::literature();
        assert!(
            high.lp_propagation_boost() > 0.0,
            "ALOX15-high boost: {}",
            high.lp_propagation_boost()
        );
        // ALOX-poor (shift expression to the least-active isoform, low activities)
        // ⇒ negative boost (more ferroptosis-resistant).
        let low = AloxConfig {
            alox15_activity: 0.3,
            alox12_activity: 0.4,
            alox5_activity: 0.5,
            alox15_frac: 0.1,
            alox12_frac: 0.2,
            alox5_frac: 0.7,
            ..AloxConfig::identity()
        };
        assert!(
            low.lp_propagation_boost() < 0.0,
            "ALOX-poor boost: {}",
            low.lp_propagation_boost()
        );
        // Floored at -1 (ALOX-null limit): all activities 0 ⇒ boost exactly -1.
        let null = AloxConfig {
            alox15_activity: 0.0,
            alox12_activity: 0.0,
            alox5_activity: 0.0,
            ..AloxConfig::identity()
        };
        assert_eq!(null.lp_propagation_boost(), -1.0);
    }

    #[test]
    fn weighted_activity_rises_with_alox15_fraction() {
        // Shifting expression toward the most-active isoform (ALOX15) raises the
        // weighted activity monotonically.
        let mostly_5 = AloxConfig {
            alox15_frac: 0.1,
            alox12_frac: 0.1,
            alox5_frac: 0.8,
            ..AloxConfig::literature()
        };
        let mostly_15 = AloxConfig {
            alox15_frac: 0.8,
            alox12_frac: 0.1,
            alox5_frac: 0.1,
            ..AloxConfig::literature()
        };
        assert!(mostly_15.weighted_activity() > mostly_5.weighted_activity());
        // No expression at all ⇒ falls back to reference (boost 0).
        let none = AloxConfig {
            alox15_frac: 0.0,
            alox12_frac: 0.0,
            alox5_frac: 0.0,
            ..AloxConfig::identity()
        };
        assert_eq!(none.weighted_activity(), none.reference_activity);
        assert_eq!(none.lp_propagation_boost(), 0.0);
    }

    #[test]
    fn mcfa_boost_rises_monotonically_and_is_bounded() {
        let s = 0.5;
        let cfg = |level: f64| AloxConfig {
            mcfa_level: level,
            mcfa_strength: s,
            ..AloxConfig::identity()
        };
        assert_eq!(cfg(0.0).mcfa_pufa_boost(), 0.0); // no MCFA ⇒ no boost
        let b1 = cfg(0.5).mcfa_pufa_boost();
        let b2 = cfg(2.0).mcfa_pufa_boost();
        let b3 = cfg(1000.0).mcfa_pufa_boost();
        assert!(0.0 < b1 && b1 < b2 && b2 < b3, "monotone: 0 {b1} {b2} {b3}");
        assert!(
            b3 <= s + 1e-9 && b3 > 0.9 * s,
            "saturates toward strength: {b3}"
        );
    }
}
