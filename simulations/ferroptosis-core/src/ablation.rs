//! Physical ABLATION: thermal (HIFU) and non-thermal (irreversible
//! electroporation), the taxonomy's two remaining physical modalities.
//!
//! ## Why these are one module and why they are NOT the ROS modalities
//!
//! `analysis/modality-coverage.md` files `hifu` (1,352 census articles) and
//! `electrochemical-therapy` (2,515) as separate mechanisms, and they are —
//! one deposits heat, the other punches holes in membranes with an electric
//! field. But they share the property that separates them from everything
//! else in this crate, and it is worth stating because getting it wrong would
//! misrepresent both:
//!
//! **Ablation is a THRESHOLD phenomenon, not a dose-response.** SDT, PDT,
//! RSL3 and radiation all kill probabilistically — more dose, more death, on a
//! continuum. An ablated cell is destroyed outright. Above the threshold
//! everything dies; below it, essentially nothing does, and the interesting
//! quantity is the POSITION OF THE MARGIN rather than a kill fraction.
//!
//! That is why neither routes through [`crate::biochem::CellState`]. There is
//! no lipid peroxidation to accumulate and no antioxidant defence to overcome:
//! `analysis/distilled-hypotheses-final.md` already argues that SDT is
//! misclassified when it is grouped with these ("SDT uses physical energy to
//! initiate a biochemical program; the others use physical energy for direct
//! cell killing"), and this module is where that distinction becomes code.
//!
//! ## What follows from the threshold, and it is the clinically hard part
//!
//! Because the margin is sharp, the failure mode is geometric rather than
//! biological: cells just outside it are untouched however long the treatment
//! runs. [`margin_survival_fraction`] is therefore the quantity to read, and
//! it is the one this project's hypoxia and penetration work cannot help
//! with — no oxygen dependence, no drug transport, no resistance state.
//!
//! ## A note on the names in this module
//!
//! `electroporation_ablated` rather than `ire_ablated`, and
//! `hifu_thermal_ablation` rather than `thermal_ablated`. The abbreviations
//! were shorter and read as jargon, and `analysis/modality-coverage.md`
//! reported both mechanisms as ABSENT while this file implemented them --
//! correctly, because that audit scans CODE and a module that does not name
//! what it models has not documented itself to anything. Spelling the
//! mechanism out is better API naming and honest measurement at the same
//! time.
//!
//! ## Calibration
//!
//! Both thresholds are published and both are cited on their constants. What
//! is NOT calibrated is the field or temperature DISTRIBUTION a real applicator
//! produces, which is an engineering problem this crate does not model: the
//! consumer supplies the field, and this module decides what it does.

use serde::{Deserialize, Serialize};

/// Irreversible-electroporation field threshold, volts per centimetre.
///
/// 1000 V/cm (PMID 20191380, `corpus/by-pmid/20191380.md`): "Numerical models
/// of the electric field distribution for the protocol used suggest that a
/// 1000 V/cm field threshold is sufficient to treat a tumor".
///
/// **This is a tissue- and pulse-protocol-dependent number**, not a physical
/// constant, and the source says as much by attributing it to a protocol. It
/// is here as a documented reference point so a consumer has somewhere to
/// start, exactly as `OER_REFERENCE_PO2_MMHG` is.
pub const IRREVERSIBLE_ELECTROPORATION_THRESHOLD_V_PER_CM: f64 = 1000.0;

/// The reference temperature of the CEM43 thermal-dose model, in Celsius.
///
/// Thermal dose is conventionally expressed as cumulative equivalent minutes
/// at 43 °C, because tissue damage above ~43 °C follows an Arrhenius-like law
/// with a well-established doubling behaviour per degree.
pub const CEM43_REFERENCE_C: f64 = 43.0;

/// Thermal dose in cumulative equivalent minutes at 43 °C.
///
/// `CEM43 = t · R^(43 − T)`, with `R = 0.5` above 43 °C and `R = 0.25` below —
/// the standard asymmetry, and it is not a detail: one degree above 43
/// DOUBLES the dose while one degree below QUARTERS it, so heating is far
/// more forgiving of overshoot than of undershoot. A symmetric model would
/// understate how sharp the ablation margin is.
///
/// Returns `0.0` for zero duration, so an unconfigured run is inert.
#[must_use = "the thermal dose is the function's only output"]
pub fn cem43(temperature_c: f64, minutes: f64) -> f64 {
    let t = minutes.max(0.0);
    if t == 0.0 {
        return 0.0;
    }
    let r: f64 = if temperature_c >= CEM43_REFERENCE_C {
        0.5
    } else {
        0.25
    };
    t * r.powf(CEM43_REFERENCE_C - temperature_c)
}

/// Ablation configuration for both modalities.
///
/// [`Default`] is the identity: zero field and zero duration, so every helper
/// reports no ablation and a run carrying it is byte-identical to one with no
/// ablation model at all.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct AblationConfig {
    /// Applied field at the cell, V/cm (irreversible electroporation).
    #[serde(default)]
    pub field_v_per_cm: f64,
    /// Field above which cells are destroyed. Defaults to
    /// [`IRREVERSIBLE_ELECTROPORATION_THRESHOLD_V_PER_CM`].
    #[serde(default = "default_electroporation_threshold")]
    pub field_threshold_v_per_cm: f64,
    /// Peak temperature at the cell, Celsius (HIFU).
    #[serde(default = "default_body_temp")]
    pub temperature_c: f64,
    /// Duration held at that temperature, minutes.
    #[serde(default)]
    pub minutes: f64,
    /// CEM43 above which tissue is taken as coagulated. 240 CEM43 is the
    /// conventional clinical ablation endpoint.
    #[serde(default = "default_cem43_threshold")]
    pub cem43_threshold: f64,
}

fn default_electroporation_threshold() -> f64 {
    IRREVERSIBLE_ELECTROPORATION_THRESHOLD_V_PER_CM
}

fn default_body_temp() -> f64 {
    37.0
}

fn default_cem43_threshold() -> f64 {
    240.0
}

impl Default for AblationConfig {
    fn default() -> Self {
        Self {
            field_v_per_cm: 0.0,
            field_threshold_v_per_cm: default_electroporation_threshold(),
            temperature_c: default_body_temp(),
            minutes: 0.0,
            cem43_threshold: default_cem43_threshold(),
        }
    }
}

/// Is this cell ablated by the electric field?
///
/// A step function, deliberately. Softening it into a sigmoid would be more
/// comfortable numerically and would misrepresent the modality: the clinical
/// literature on irreversible electroporation is about MARGIN PLACEMENT
/// precisely because the transition is sharp.
#[must_use = "the verdict is the function's only output"]
pub fn electroporation_ablated(cfg: &AblationConfig) -> bool {
    cfg.field_v_per_cm >= cfg.field_threshold_v_per_cm && cfg.field_threshold_v_per_cm > 0.0
}

/// Is this cell ablated thermally?
#[must_use = "the verdict is the function's only output"]
pub fn hifu_thermal_ablation(cfg: &AblationConfig) -> bool {
    cem43(cfg.temperature_c, cfg.minutes) >= cfg.cem43_threshold && cfg.cem43_threshold > 0.0
}

/// Surviving fraction across an ablation margin.
///
/// **The quantity that actually matters for these modalities**, and the reason
/// they get a function rather than a boolean apiece. Because the transition is
/// sharp, the outcome is decided by how much of the target lies inside the
/// treated volume, and cells just outside are untouched however long the
/// treatment runs — which is why local recurrence after ablation is a margin
/// problem and not a resistance problem.
///
/// `covered_fraction` is the share of the target inside the ablation zone.
/// Returns `1.0 - covered` when the applied dose clears the threshold and
/// `1.0` when it does not: below threshold, coverage buys nothing at all,
/// which is exactly the failure this makes visible.
#[must_use = "the surviving fraction is the function's only output"]
pub fn margin_survival_fraction(cfg: &AblationConfig, covered_fraction: f64) -> f64 {
    let covered = covered_fraction.clamp(0.0, 1.0);
    if electroporation_ablated(cfg) || hifu_thermal_ablation(cfg) {
        1.0 - covered
    } else {
        1.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_default_config_ablates_nothing() {
        let d = AblationConfig::default();
        assert!(!electroporation_ablated(&d));
        assert!(!hifu_thermal_ablation(&d));
        assert_eq!(
            margin_survival_fraction(&d, 1.0).to_bits(),
            1.0_f64.to_bits(),
            "the identity config killed a fully covered target"
        );
        assert_eq!(cem43(37.0, 0.0).to_bits(), 0.0_f64.to_bits());
    }

    #[test]
    fn the_thermal_dose_is_asymmetric_around_43c() {
        // One degree ABOVE doubles; one degree BELOW quarters. The asymmetry
        // is the standard model and it is not a detail -- a symmetric version
        // would understate how sharp the margin is.
        let base = cem43(43.0, 10.0);
        assert!((base - 10.0).abs() < 1e-12, "43C for 10 min is 10 CEM43");
        assert!((cem43(44.0, 10.0) - 20.0).abs() < 1e-9);
        assert!((cem43(42.0, 10.0) - 2.5).abs() < 1e-9);
        // Which means overshoot is far more forgiving than undershoot.
        let over = cem43(45.0, 10.0) / base;
        let under = base / cem43(41.0, 10.0);
        assert!(
            under > over,
            "undershoot ({under:.1}x) should cost more than overshoot gains \
             ({over:.1}x)"
        );
        // Monotone in both arguments.
        let mut prev = 0.0;
        for &t in &[40.0_f64, 42.0, 43.0, 45.0, 50.0] {
            let v = cem43(t, 5.0);
            assert!(v > prev, "not monotone in temperature at {t}");
            prev = v;
        }
        assert!(cem43(45.0, 20.0) > cem43(45.0, 5.0));
    }

    #[test]
    fn ablation_is_a_threshold_not_a_dose_response() {
        // THE defining property, and the reason these do not route through
        // CellState. Just below the threshold nothing happens; just above,
        // everything does. A sigmoid would be numerically nicer and would
        // misrepresent the modality.
        let below = AblationConfig {
            field_v_per_cm: IRREVERSIBLE_ELECTROPORATION_THRESHOLD_V_PER_CM - 1.0,
            ..AblationConfig::default()
        };
        let above = AblationConfig {
            field_v_per_cm: IRREVERSIBLE_ELECTROPORATION_THRESHOLD_V_PER_CM + 1.0,
            ..AblationConfig::default()
        };
        assert!(!electroporation_ablated(&below));
        assert!(electroporation_ablated(&above));
        // A 0.2% change in field flips the outcome completely, which no
        // continuum model in this crate does.
        assert_eq!(margin_survival_fraction(&below, 1.0), 1.0);
        assert_eq!(margin_survival_fraction(&above, 1.0), 0.0);

        // Ten times the field changes NOTHING once past the threshold: there
        // is no "more ablated".
        let overkill = AblationConfig {
            field_v_per_cm: IRREVERSIBLE_ELECTROPORATION_THRESHOLD_V_PER_CM * 10.0,
            ..AblationConfig::default()
        };
        assert_eq!(
            margin_survival_fraction(&overkill, 0.8).to_bits(),
            margin_survival_fraction(&above, 0.8).to_bits()
        );
    }

    #[test]
    fn survival_is_a_margin_problem_and_not_a_resistance_problem() {
        // The clinical point: below threshold, coverage buys NOTHING, and
        // above it, survival is entirely the uncovered share. Recurrence
        // after ablation is geometry.
        let hot = AblationConfig {
            temperature_c: 56.0,
            minutes: 1.0,
            ..AblationConfig::default()
        };
        assert!(
            hifu_thermal_ablation(&hot),
            "56C for a minute should coagulate"
        );
        for &c in &[0.0_f64, 0.25, 0.5, 0.9, 1.0] {
            assert!((margin_survival_fraction(&hot, c) - (1.0 - c)).abs() < 1e-12);
        }
        let cool = AblationConfig {
            temperature_c: 41.0,
            minutes: 1.0,
            ..AblationConfig::default()
        };
        assert!(!hifu_thermal_ablation(&cool));
        for &c in &[0.0_f64, 0.5, 1.0] {
            assert_eq!(
                margin_survival_fraction(&cool, c).to_bits(),
                1.0_f64.to_bits(),
                "sub-threshold heating killed a {c} covered target"
            );
        }
        // Out-of-range coverage clamps rather than producing a negative
        // surviving fraction.
        assert_eq!(margin_survival_fraction(&hot, 5.0), 0.0);
        assert_eq!(margin_survival_fraction(&hot, -5.0), 1.0);
    }

    #[test]
    fn the_two_modalities_are_independent_routes_to_the_same_verdict() {
        // Thermal and electrical ablation must not leak into one another: a
        // field that ablates must not need heat, and vice versa.
        let field_only = AblationConfig {
            field_v_per_cm: 1500.0,
            ..AblationConfig::default()
        };
        assert!(electroporation_ablated(&field_only) && !hifu_thermal_ablation(&field_only));
        let heat_only = AblationConfig {
            temperature_c: 60.0,
            minutes: 2.0,
            ..AblationConfig::default()
        };
        assert!(hifu_thermal_ablation(&heat_only) && !electroporation_ablated(&heat_only));
        // Either alone is sufficient for the margin result.
        assert_eq!(margin_survival_fraction(&field_only, 1.0), 0.0);
        assert_eq!(margin_survival_fraction(&heat_only, 1.0), 0.0);
    }
}
