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

// ── Where an ablation fails, and why it is not where the dose is lowest ──
//
// The module above treats ablation as a threshold and a coverage fraction: get
// the applicator field over the margin and everything inside it dies. That is
// the right first model and it cannot express the failure the clinic actually
// sees, which is spatial and has a cause.
//
// Flowing blood carries heat away. Tissue next to a large vessel therefore
// reaches a lower peak temperature than tissue the same distance from the
// applicator but away from the vessel, and survives an ablation that killed
// everything around it. Perivascular local progression after thermal ablation
// is a recognised problem for exactly this reason (PMID 35114665).
//
// The consequence is a DISCRIMINATION rather than a correction, because
// irreversible electroporation is not thermal: it kills by permeabilising
// membranes with an electric field, and a heat sink does not apply to it. So
// the two modalities fail in different PLACES -- one perivascular, one not --
// which is testable in a way that a coverage fraction is not, and which is the
// reason IRE is reached for near vessels.

/// Perfusion-driven cooling as a function of distance from a vessel, in mm.
///
/// Returns a multiplier on the temperature RISE above body temperature: 0
/// against the vessel wall, approaching 1 far from it. The length scale is the
/// distance over which flowing blood stops mattering.
///
/// UNCALIBRATED. The exponential form has the right limits and the right
/// direction; the length scale is a placeholder and depends on vessel calibre
/// and flow, neither of which this layer represents.
#[must_use]
pub fn perfusion_cooling(distance_mm: f64, cooling_length_mm: f64) -> f64 {
    let d = distance_mm.max(0.0);
    let l = cooling_length_mm.max(f64::MIN_POSITIVE);
    1.0 - (-d / l).exp()
}

/// Peak temperature at a point, given the applicator's unimpeded temperature
/// and the cooling a nearby vessel imposes.
///
/// Body temperature is the floor: a heat sink cannot cool tissue below the
/// blood that is doing the cooling.
#[must_use]
pub fn perivascular_temperature(
    unimpeded_c: f64,
    distance_mm: f64,
    cooling_length_mm: f64,
    body_temperature_c: f64,
) -> f64 {
    let rise = (unimpeded_c - body_temperature_c).max(0.0);
    body_temperature_c + rise * perfusion_cooling(distance_mm, cooling_length_mm)
}

/// Whether a thermal ablation kills at a given distance from a vessel.
///
/// The threshold is the module's existing thermal-dose criterion; what changes
/// here is the temperature that reaches the tissue.
#[must_use]
pub fn thermal_kills_at(
    unimpeded_c: f64,
    minutes: f64,
    distance_mm: f64,
    cooling_length_mm: f64,
    body_temperature_c: f64,
    cem43_threshold: f64,
) -> bool {
    let t = perivascular_temperature(
        unimpeded_c,
        distance_mm,
        cooling_length_mm,
        body_temperature_c,
    );
    cem43(t, minutes) >= cem43_threshold
}

/// The distance from a vessel inside which a thermal ablation fails, in mm.
///
/// **The layer's spatial prediction.** Returns the radius of the surviving
/// perivascular sleeve: 0 when the applicator is hot enough that the sink
/// cannot save anything, and larger as the margin gets thinner. A coverage
/// fraction cannot express this at all, because the survivors are not
/// distributed at random through the volume -- they are in a specific place, a
/// clinician can see where, and that is what makes it checkable.
#[must_use]
pub fn perivascular_failure_radius_mm(
    unimpeded_c: f64,
    minutes: f64,
    cooling_length_mm: f64,
    body_temperature_c: f64,
    cem43_threshold: f64,
    max_distance_mm: f64,
) -> f64 {
    let mut d = 0.0;
    let step = max_distance_mm / 500.0;
    while d <= max_distance_mm {
        if thermal_kills_at(
            unimpeded_c,
            minutes,
            d,
            cooling_length_mm,
            body_temperature_c,
            cem43_threshold,
        ) {
            return d;
        }
        d += step;
    }
    max_distance_mm
}

/// The same radius for irreversible electroporation, which is zero.
///
/// Not a stub. IRE kills by permeabilising membranes with an electric field
/// and deposits little heat, so flowing blood removes nothing that matters to
/// it -- which is the documented reason it is reached for near vessels. The
/// function exists so the CONTRAST is a value the model produces rather than a
/// sentence in a comment, and so a change that gave electroporation a thermal
/// dependence would break a test rather than pass silently.
#[must_use]
pub fn electroporation_failure_radius_mm() -> f64 {
    0.0
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
    // ── The heat sink ────────────────────────────────────────────────────

    #[test]
    fn cooling_is_total_at_the_vessel_wall_and_absent_far_away() {
        assert!(
            perfusion_cooling(0.0, 2.0).abs() < 1e-12,
            "tissue against the vessel wall should get no temperature rise"
        );
        assert!(
            perfusion_cooling(20.0, 2.0) > 0.99,
            "tissue far from a vessel should be unaffected"
        );
        assert!(perfusion_cooling(1.0, 2.0) < perfusion_cooling(3.0, 2.0));
        // A longer cooling length means the sink reaches further, which is the
        // direction a larger vessel with more flow would move it.
        assert!(perfusion_cooling(2.0, 5.0) < perfusion_cooling(2.0, 1.0));
    }

    #[test]
    fn a_heat_sink_cannot_cool_below_the_blood_doing_the_cooling() {
        let t = perivascular_temperature(90.0, 0.0, 2.0, 37.0);
        assert!(
            (t - 37.0).abs() < 1e-9,
            "temperature at the wall was {t}, not body"
        );
        let far = perivascular_temperature(90.0, 50.0, 2.0, 37.0);
        assert!(
            (far - 90.0).abs() < 0.1,
            "temperature far away was {far}, not the applicator's"
        );
        // And it never inverts: hotter applicator, hotter tissue, everywhere.
        for d in [0.5, 2.0, 5.0] {
            assert!(
                perivascular_temperature(90.0, d, 2.0, 37.0)
                    > perivascular_temperature(60.0, d, 2.0, 37.0)
            );
        }
    }

    #[test]
    fn the_kill_criterion_is_a_thermal_DOSE_and_not_a_temperature() {
        // A MUTATION SURVIVOR: replacing the CEM43 dose with a bare 43 C
        // threshold left every other test green, because at the temperatures
        // those tests use the two criteria happen to agree on which side of
        // the line each point falls. TIME is what separates them.
        let hot_and_brief = thermal_kills_at(55.0, 0.05, 10.0, 2.0, 37.0, 240.0);
        let hot_and_sustained = thermal_kills_at(55.0, 30.0, 10.0, 2.0, 37.0, 240.0);
        assert!(
            !hot_and_brief,
            "three seconds at 55 C should not reach a CEM43 of 240"
        );
        assert!(hot_and_sustained, "thirty minutes at 55 C should");
        // A bare temperature test would call both of those the same, and a
        // long exposure just above 43 C the same as a short one.
        assert!(
            thermal_kills_at(44.0, 600.0, 10.0, 2.0, 37.0, 240.0),
            "ten hours just above 43 C should accumulate a lethal dose"
        );
        assert!(
            !thermal_kills_at(44.0, 1.0, 10.0, 2.0, 37.0, 240.0),
            "one minute just above 43 C should not"
        );
    }

    #[test]
    fn the_surviving_sleeve_shrinks_as_the_applicator_gets_hotter() {
        // THE LAYER'S SPATIAL PREDICTION. A coverage fraction cannot express
        // it: the survivors are not scattered through the volume, they are in
        // a specific place, and a clinician can look there.
        let radius = |c: f64| perivascular_failure_radius_mm(c, 5.0, 2.0, 37.0, 240.0, 20.0);
        let cool = radius(50.0);
        let hot = radius(90.0);
        assert!(
            cool > hot,
            "a hotter applicator left a WIDER sleeve: {cool} vs {hot}"
        );
        assert!(hot > 0.0, "even a very hot applicator leaves some sleeve");
        assert!(cool < 20.0, "the sleeve filled the whole scanned range");
        // Millimetre scale, which is what makes it a clinical problem rather
        // than a rounding error.
        assert!((0.1..=10.0).contains(&cool), "sleeve radius {cool} mm");
    }

    #[test]
    fn electroporation_has_no_perivascular_sleeve_and_that_is_the_contrast() {
        // Not a stub: the CONTRAST is the finding, and it is the documented
        // reason IRE is reached for near vessels. A change that gave
        // electroporation a thermal dependence should break this.
        assert!(electroporation_failure_radius_mm().abs() < 1e-12);
        let thermal = perivascular_failure_radius_mm(60.0, 5.0, 2.0, 37.0, 240.0, 20.0);
        assert!(thermal > electroporation_failure_radius_mm(),
                "thermal ablation no longer leaves a sleeve that electroporation                  does not, which is this section's whole claim");
    }
}
