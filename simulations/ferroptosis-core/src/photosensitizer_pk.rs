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
    /// Parse a CLI `--photosensitizer` SPEC string into a `Photosensitizer`.
    ///
    /// Accepted forms (case-insensitive on the variant name):
    /// - `uniform` → `Uniform(1.0)`
    /// - `uniform=N` → `Uniform(N)`
    /// - `porfimer` → `Porfimer { t_half_h: 504.0 }` (Bellnier 2006 t½ in hours)
    /// - `porfimer=N` → `Porfimer { t_half_h: N }`
    ///
    /// Errors on unknown variant, unparseable number, or any value that
    /// fails [`Photosensitizer::validate`]. Round-trips with the
    /// `Display` impl: `Photosensitizer::from_cli_spec(&format!("{ps}"))
    /// == Ok(ps)` for every valid `Photosensitizer`.
    pub fn from_cli_spec(s: &str) -> Result<Self, String> {
        let s = s.trim();
        let (name, value) = match s.split_once('=') {
            Some((n, v)) => (n.trim(), Some(v.trim())),
            None => (s, None),
        };
        // `eq_ignore_ascii_case` avoids allocating a lowercased String
        // for every CLI parse. Match-on-lowercase would be tidier syntax
        // but allocates per call.
        let ps = if name.eq_ignore_ascii_case("uniform") {
            let c = match value {
                Some(v) => v
                    .parse::<f64>()
                    .map_err(|e| format!("uniform=N: cannot parse N={v:?}: {e}"))?,
                None => 1.0,
            };
            Self::Uniform(c)
        } else if name.eq_ignore_ascii_case("porfimer") {
            let t_half_h = match value {
                Some(v) => v
                    .parse::<f64>()
                    .map_err(|e| format!("porfimer=N: cannot parse N={v:?}: {e}"))?,
                None => 504.0, // Bellnier 2006 terminal t½ in humans, hours.
            };
            Self::Porfimer { t_half_h }
        } else {
            return Err(format!(
                "unknown photosensitizer {name:?}; expected one of: uniform, uniform=N, porfimer, porfimer=N (case-insensitive)"
            ));
        };
        ps.validate()?;
        Ok(ps)
    }

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

/// Validate a drug-light-interval value (hours): reject NaN, negative,
/// and infinite inputs so they cannot reach `concentration_at` and
/// produce silent non-physical PDT output. Pair with
/// [`Photosensitizer::from_cli_spec`] at the same parse-time gate.
///
/// Error messages refer to the field as `dli_h` (matching the parameter
/// name) rather than embedding any specific CLI flag name; callers are
/// free to prefix with their flag spelling — e.g. sim-spatial wraps the
/// error as `--dli-h: <message>`.
pub fn validate_dli_h(dli_h: f64) -> Result<(), String> {
    if !dli_h.is_finite() {
        return Err(format!("dli_h must be finite, got {dli_h}"));
    }
    if dli_h < 0.0 {
        return Err(format!("dli_h must be >= 0, got {dli_h}"));
    }
    Ok(())
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

    // -- from_cli_spec --

    #[test]
    fn from_cli_spec_uniform_default() {
        assert_eq!(
            Photosensitizer::from_cli_spec("uniform").unwrap(),
            Photosensitizer::Uniform(1.0)
        );
    }

    #[test]
    fn from_cli_spec_uniform_with_value() {
        assert_eq!(
            Photosensitizer::from_cli_spec("uniform=0.5").unwrap(),
            Photosensitizer::Uniform(0.5)
        );
    }

    #[test]
    fn from_cli_spec_porfimer_default_is_bellnier() {
        assert_eq!(
            Photosensitizer::from_cli_spec("porfimer").unwrap(),
            Photosensitizer::Porfimer { t_half_h: 504.0 }
        );
    }

    #[test]
    fn from_cli_spec_porfimer_with_value() {
        assert_eq!(
            Photosensitizer::from_cli_spec("porfimer=336").unwrap(),
            Photosensitizer::Porfimer { t_half_h: 336.0 }
        );
    }

    #[test]
    fn from_cli_spec_trims_whitespace() {
        assert_eq!(
            Photosensitizer::from_cli_spec("  porfimer = 504  ").unwrap(),
            Photosensitizer::Porfimer { t_half_h: 504.0 }
        );
    }

    #[test]
    fn from_cli_spec_case_insensitive() {
        assert_eq!(
            Photosensitizer::from_cli_spec("Uniform").unwrap(),
            Photosensitizer::Uniform(1.0)
        );
        assert_eq!(
            Photosensitizer::from_cli_spec("PORFIMER=504").unwrap(),
            Photosensitizer::Porfimer { t_half_h: 504.0 }
        );
        assert_eq!(
            Photosensitizer::from_cli_spec("PorFimEr").unwrap(),
            Photosensitizer::Porfimer { t_half_h: 504.0 }
        );
    }

    #[test]
    fn from_cli_spec_unknown_variant_errors() {
        let err = Photosensitizer::from_cli_spec("photochlor").unwrap_err();
        assert!(err.contains("photochlor"));
        assert!(err.contains("uniform"));
        assert!(err.contains("porfimer"));
    }

    #[test]
    fn from_cli_spec_unparseable_number_errors() {
        let err = Photosensitizer::from_cli_spec("porfimer=abc").unwrap_err();
        assert!(err.contains("porfimer=N"));
        assert!(err.contains("abc"));
    }

    #[test]
    fn from_cli_spec_negative_t_half_h_rejected() {
        let err = Photosensitizer::from_cli_spec("porfimer=-1").unwrap_err();
        assert!(err.contains("t_half_h"));
    }

    #[test]
    fn from_cli_spec_zero_t_half_h_rejected() {
        let err = Photosensitizer::from_cli_spec("porfimer=0").unwrap_err();
        assert!(err.contains("t_half_h"));
    }

    #[test]
    fn from_cli_spec_negative_uniform_rejected() {
        let err = Photosensitizer::from_cli_spec("uniform=-0.5").unwrap_err();
        assert!(err.contains("must be >= 0"));
    }

    #[test]
    fn from_cli_spec_nan_rejected() {
        let err = Photosensitizer::from_cli_spec("uniform=NaN").unwrap_err();
        assert!(err.contains("must be finite"));
    }

    #[test]
    fn display_round_trips_through_from_cli_spec() {
        // Display's contract is round-trip parseability via from_cli_spec.
        for ps in [
            Photosensitizer::Uniform(1.0),
            Photosensitizer::Uniform(0.5),
            Photosensitizer::Porfimer { t_half_h: 504.0 },
            Photosensitizer::Porfimer { t_half_h: 336.5 },
        ] {
            let rendered = format!("{ps}");
            let reparsed = Photosensitizer::from_cli_spec(&rendered)
                .unwrap_or_else(|e| panic!("round-trip failed for {ps:?}: {e}"));
            assert_eq!(reparsed, ps, "round-trip mismatch via {rendered:?}");
        }
    }

    // -- validate_dli_h --

    #[test]
    fn validate_dli_h_accepts_zero_and_positive() {
        assert!(validate_dli_h(0.0).is_ok());
        assert!(validate_dli_h(24.0).is_ok());
        assert!(validate_dli_h(504.0).is_ok());
        assert!(validate_dli_h(1e9).is_ok());
    }

    #[test]
    fn validate_dli_h_rejects_negative() {
        let err = validate_dli_h(-1.0).unwrap_err();
        // Error message is library-friendly (`dli_h`), not CLI-coupled
        // (`--dli-h`); callers are free to prefix with their flag name.
        assert!(err.contains("dli_h"));
        assert!(err.contains(">= 0"));
    }

    #[test]
    fn validate_dli_h_rejects_nan() {
        let err = validate_dli_h(f64::NAN).unwrap_err();
        assert!(err.contains("finite"));
    }

    #[test]
    fn validate_dli_h_rejects_infinity() {
        assert!(validate_dli_h(f64::INFINITY).unwrap_err().contains("finite"));
        assert!(validate_dli_h(f64::NEG_INFINITY)
            .unwrap_err()
            .contains("finite"));
    }
}
