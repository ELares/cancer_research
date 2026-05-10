//! Photosensitizer plasma pharmacokinetics for PDT.
//!
//! Models the time-decaying concentration of a systemically administered
//! photosensitizer between dosing and light delivery (the drug-light
//! interval, DLI). Used by `physics::pdt_intensity_at_depth` to scale
//! depth-attenuated light intensity by the fraction of drug still present
//! when illumination occurs.
//!
//! v1 captures the *temporal* PK only — intra-drug C(t) decay between
//! administration and illumination. It does NOT capture inter-drug
//! ROS-yield differences (singlet-O2 quantum yield phi_so2), 5-ALA → PpIX
//! intracellular accumulation, photobleaching during illumination, or
//! formulation-dependent Vd / protein binding. Those are intentional
//! follow-ups; see issue #200.
//!
//! `t_half_h` represents *plasma* terminal half-life. Cellular
//! concentration is assumed to track plasma proportionally — a reasonable
//! approximation for porfimer (slow-distributing, weeks-scale t½) but
//! explicitly wrong for 5-ALA/PpIX, which accumulates intracellularly via
//! ferrochelatase deficiency rather than decaying. ALA kinetics need a
//! different variant.
//!
//! The model's `t = 0` represents *post-distribution peak*, not the moment
//! of IV bolus. Porfimer distribution takes 1–2 days; users setting
//! `t_drug_light_interval_h = 0` are modeling "light at peak," not
//! "light at injection."
//!
//! Reference: Bellnier DA et al., Lasers Surg Med 38(5):439-444, 2006 —
//! clinical PK of porfimer sodium, Photochlor, and 5-ALA-induced PpIX
//! in humans (PMID 16634075).

use std::fmt;

use serde::{Deserialize, Serialize};

/// Photosensitizer pharmacokinetic model used to scale PDT light dose by
/// the fraction of drug present at the moment of illumination.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Photosensitizer {
    /// Drug at constant fraction `c` of peak for all `t`. The default
    /// `Uniform(1.0)` is the identity case (no scaling, byte-exact pre-PK
    /// behavior). `Uniform(0.5)` asserts "drug constant at half-peak forever"
    /// — itself a PK assertion, not "no PK model."
    Uniform(f64),
    /// Porfimer sodium (Photofrin). Single-exponential plasma decay.
    /// Bellnier 2006 reports terminal t½ ~21 d (504 h) in humans, with
    /// substantial infusion-protocol-dependent variability (~250–500+ h).
    Porfimer {
        /// Plasma terminal half-life in hours. Must be strictly positive.
        t_half_h: f64,
    },
}

impl Default for Photosensitizer {
    fn default() -> Self {
        Self::Uniform(1.0)
    }
}

impl fmt::Display for Photosensitizer {
    /// Human-readable form that round-trips through the CLI spec parser
    /// (`name[=value]`, lowercase variant name). A user can copy a
    /// `Photosensitizer: ...` line from stderr and pass it back to
    /// `--photosensitizer` verbatim.
    ///
    /// Note: relies on `f64`'s `Display` impl rendering whole numbers
    /// without a decimal point (`504.0_f64` → `"504"`). This is
    /// long-standing stable behavior but technically implementation-
    /// defined; if it ever changes, both this output and the round-trip
    /// tests will need updating.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Uniform(c) => write!(f, "uniform={c}"),
            Self::Porfimer { t_half_h } => write!(f, "porfimer={t_half_h}"),
        }
    }
}

impl Photosensitizer {
    /// Validate that the photosensitizer's parameters are physically
    /// reasonable. Returns `Ok(())` for valid configurations and a
    /// descriptive `Err` otherwise.
    ///
    /// Default `Uniform(1.0)` is always valid. Use this when loading a
    /// `Photosensitizer` from a config file or external source where you
    /// don't trust the values.
    ///
    /// `concentration_at` itself uses `debug_assert!` against this contract
    /// — invalid configs panic in test/debug builds and produce non-physical
    /// outputs in release builds.
    pub fn validate(&self) -> Result<(), String> {
        match self {
            Self::Uniform(c) => {
                if !c.is_finite() {
                    Err(format!(
                        "Photosensitizer::Uniform multiplier must be finite, got {c}"
                    ))
                } else if *c < 0.0 {
                    Err(format!(
                        "Photosensitizer::Uniform multiplier must be >= 0, got {c}"
                    ))
                } else {
                    Ok(())
                }
            }
            Self::Porfimer { t_half_h } => {
                if !t_half_h.is_finite() {
                    Err(format!(
                        "Photosensitizer::Porfimer.t_half_h must be finite, got {t_half_h}"
                    ))
                } else if *t_half_h <= 0.0 {
                    Err(format!(
                        "Photosensitizer::Porfimer.t_half_h must be > 0, got {t_half_h}"
                    ))
                } else {
                    Ok(())
                }
            }
        }
    }

    /// Fractional drug presence at time `t_h` after administration,
    /// relative to peak (post-distribution).
    ///
    /// - `Uniform(c)` returns `c` for any `t_h`.
    /// - `Porfimer { t_half_h }` returns `exp(-ln(2) * t_h / t_half_h)`
    ///   for non-negative `t_h`. Negative `t_h` (illumination "before"
    ///   the post-distribution peak) saturates to peak (returns 1.0 for
    ///   the porfimer kinetics).
    ///
    /// Invalid configurations (NaN, non-positive `t_half_h`, negative
    /// `Uniform` multiplier) trigger `debug_assert!` failures in tests
    /// and CI. Use [`Self::validate`] to check ahead of time when loading
    /// from untrusted sources.
    pub fn concentration_at(&self, t_h: f64) -> f64 {
        debug_assert!(
            self.validate().is_ok(),
            "invalid Photosensitizer: {:?}",
            self.validate().err()
        );
        // f64::max treats NaN as smaller than any value, so NaN inputs
        // saturate to 0.0 and produce a deterministic result.
        let t_h = t_h.max(0.0);
        match self {
            Self::Uniform(c) => *c,
            Self::Porfimer { t_half_h } => {
                let k = 2.0_f64.ln() / *t_half_h;
                (-k * t_h).exp()
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_is_uniform_one() {
        assert_eq!(Photosensitizer::default(), Photosensitizer::Uniform(1.0));
        // Default photosensitizer must return 1.0 at every t for byte-exact
        // backwards-compat with pre-PK PDT physics.
        assert_eq!(Photosensitizer::default().concentration_at(0.0), 1.0);
        assert_eq!(Photosensitizer::default().concentration_at(24.0), 1.0);
        assert_eq!(Photosensitizer::default().concentration_at(1e6), 1.0);
    }

    #[test]
    fn porfimer_at_one_halflife_is_half() {
        let p = Photosensitizer::Porfimer { t_half_h: 504.0 };
        // `exp(-ln(2))` lands on 0.5 on most libm implementations but is
        // not guaranteed bit-exact across platforms (intermediate rounding
        // in `ln(2)/504 * 504` may not round-trip). Compare with a tight
        // tolerance instead of strict equality.
        let c = p.concentration_at(504.0);
        assert!(
            (c - 0.5).abs() < 1e-12,
            "concentration_at(t_half_h) = {c}, expected ~0.5"
        );
    }

    #[test]
    fn porfimer_at_zero_is_one() {
        // `exp(0.0) == 1.0` is required by IEEE 754 — strict equality is safe here.
        let p = Photosensitizer::Porfimer { t_half_h: 504.0 };
        assert_eq!(p.concentration_at(0.0), 1.0);
    }

    #[test]
    fn porfimer_decays_monotonically() {
        let p = Photosensitizer::Porfimer { t_half_h: 504.0 };
        let c0 = p.concentration_at(0.0);
        let c1 = p.concentration_at(24.0);
        let c2 = p.concentration_at(168.0);
        let c3 = p.concentration_at(504.0);
        assert!(c0 > c1 && c1 > c2 && c2 > c3);
        // exp(0.0) == 1.0 exactly per IEEE 754; one-half-life value is
        // libm-dependent so use epsilon comparison.
        assert_eq!(c0, 1.0);
        assert!((c3 - 0.5).abs() < 1e-12);
    }

    #[test]
    fn uniform_constant_at_all_times() {
        // Uniform(c) returns c regardless of t_h, including c != 1.0.
        let p = Photosensitizer::Uniform(0.5);
        assert_eq!(p.concentration_at(0.0), 0.5);
        assert_eq!(p.concentration_at(100.0), 0.5);
        assert_eq!(p.concentration_at(1e9), 0.5);
        // Non-default value must not be silently treated as 1.0.
        assert_ne!(p.concentration_at(50.0), 1.0);
    }

    // -- Display --

    // The Display assertions below rely on `f64`'s Display impl
    // formatting whole numbers without a decimal point (e.g. `1.0_f64`
    // renders as `"1"`). This is long-standing stable Rust behavior.

    #[test]
    fn display_uniform() {
        assert_eq!(format!("{}", Photosensitizer::Uniform(1.0)), "uniform=1");
        assert_eq!(format!("{}", Photosensitizer::Uniform(0.5)), "uniform=0.5");
    }

    #[test]
    fn display_porfimer() {
        assert_eq!(
            format!("{}", Photosensitizer::Porfimer { t_half_h: 504.0 }),
            "porfimer=504"
        );
        assert_eq!(
            format!("{}", Photosensitizer::Porfimer { t_half_h: 336.5 }),
            "porfimer=336.5"
        );
    }

    #[test]
    fn serde_roundtrip_uniform() {
        let p = Photosensitizer::Uniform(0.7);
        let json = serde_json::to_string(&p).unwrap();
        let q: Photosensitizer = serde_json::from_str(&json).unwrap();
        assert_eq!(p, q);
    }

    #[test]
    fn serde_roundtrip_porfimer() {
        let p = Photosensitizer::Porfimer { t_half_h: 504.0 };
        let json = serde_json::to_string(&p).unwrap();
        let q: Photosensitizer = serde_json::from_str(&json).unwrap();
        assert_eq!(p, q);
    }

    // -- input validation --

    #[test]
    fn validate_uniform_accepts_valid() {
        assert!(Photosensitizer::Uniform(0.0).validate().is_ok());
        assert!(Photosensitizer::Uniform(0.5).validate().is_ok());
        assert!(Photosensitizer::Uniform(1.0).validate().is_ok());
        // Values > 1.0 are intentionally permissive (forward-compat hook
        // for enrichment factors). Tighten in a follow-up if needed.
        assert!(Photosensitizer::Uniform(2.0).validate().is_ok());
    }

    #[test]
    fn validate_uniform_rejects_invalid() {
        assert!(Photosensitizer::Uniform(-0.5).validate().is_err());
        assert!(Photosensitizer::Uniform(f64::NAN).validate().is_err());
        assert!(Photosensitizer::Uniform(f64::INFINITY).validate().is_err());
        assert!(Photosensitizer::Uniform(f64::NEG_INFINITY)
            .validate()
            .is_err());
    }

    #[test]
    fn validate_porfimer_accepts_valid() {
        assert!(Photosensitizer::Porfimer { t_half_h: 504.0 }
            .validate()
            .is_ok());
        assert!(Photosensitizer::Porfimer { t_half_h: 0.001 }
            .validate()
            .is_ok());
    }

    #[test]
    fn validate_porfimer_rejects_invalid() {
        assert!(Photosensitizer::Porfimer { t_half_h: 0.0 }
            .validate()
            .is_err());
        assert!(Photosensitizer::Porfimer { t_half_h: -100.0 }
            .validate()
            .is_err());
        assert!(Photosensitizer::Porfimer { t_half_h: f64::NAN }
            .validate()
            .is_err());
        assert!(Photosensitizer::Porfimer {
            t_half_h: f64::INFINITY
        }
        .validate()
        .is_err());
    }

    #[test]
    fn negative_dli_saturates_to_peak() {
        // Negative DLI is non-physical (illumination before drug peak);
        // saturate to t=0 to keep outputs in [0, 1] for valid params.
        let p = Photosensitizer::Porfimer { t_half_h: 504.0 };
        assert_eq!(p.concentration_at(-1.0), 1.0);
        assert_eq!(p.concentration_at(-1e6), 1.0);
        assert_eq!(p.concentration_at(f64::NEG_INFINITY), 1.0);
    }

    #[test]
    fn nan_dli_saturates_to_peak() {
        // f64::max treats NaN as smaller than any value, so NaN.max(0.0)
        // returns 0.0 and concentration_at(NaN) is deterministic.
        let p = Photosensitizer::Porfimer { t_half_h: 504.0 };
        assert_eq!(p.concentration_at(f64::NAN), 1.0);
    }
}
