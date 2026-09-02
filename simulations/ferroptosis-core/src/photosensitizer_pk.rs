//! Photosensitizer plasma pharmacokinetics for PDT.
//!
//! Models the time-resolved concentration and per-photon ROS yield of a
//! systemically administered photosensitizer between dosing and light
//! delivery (the drug-light interval, DLI). `physics::pdt_intensity_at_depth`
//! scales depth-attenuated light intensity by `Photosensitizer::yield_at`.
//!
//! Closed via #200 / #203:
//! - Intra-drug temporal PK — single-exponential plasma decay
//!   (`Photosensitizer::Porfimer.t_half_h`).
//! - Saturating distribution-phase hold (`Porfimer.t_distribution_h`):
//!   drug is held at peak for the first `t_distribution_h` hours then
//!   begins exponential decay. Default `0.0` recovers the pre-#203
//!   "light at peak" model bit-exactly. With `t_distribution_h > 0`,
//!   `t_drug_light_interval_h` can be the **clinical DLI from injection**.
//! - Inter-drug ROS-yield normalization (`Porfimer.phi_so2_relative`):
//!   `yield_at(t)` returns `concentration_at(t) × phi_so2_relative`.
//!   `Params::pdt_ros = 5.0` is calibrated to porfimer at peak (yield = 1.0);
//!   future drug variants set this field to `absolute_phi_so2 / 0.65`
//!   so the calibration carries through.
//!
//! Still out of scope (intentional, separate-issue follow-ups):
//! - 5-ALA → PpIX intracellular accumulation (ferrochelatase-deficiency
//!   biology, not just PK)
//! - Photobleaching during illumination (separate temporal axis)
//! - Two-phase rising-curve distribution (current is saturating-step
//!   approximation per #203's design choice)
//! - Formulation-dependent Vd / protein binding
//!
//! `t_half_h` represents *plasma* terminal half-life. Cellular
//! concentration is assumed to track plasma proportionally — a reasonable
//! approximation for porfimer (slow-distributing, weeks-scale t½) but
//! explicitly wrong for 5-ALA/PpIX, which accumulates intracellularly via
//! ferrochelatase deficiency rather than decaying. ALA kinetics need a
//! different variant.
//!
//! With `t_distribution_h = 0` (default), `t = 0` represents
//! *post-distribution peak*. With `t_distribution_h > 0`, `t = 0` is the
//! moment of administration and the drug is held at peak for the first
//! `t_distribution_h` hours of the model.
//!
//! Reference: Bellnier DA et al., Lasers Surg Med 38(5):439-444, 2006 —
//! clinical PK of porfimer sodium, Photochlor, and 5-ALA-induced PpIX
//! in humans (PMID 16634075).
//! Reference: Wilson BC, Patterson MS, Phys Med Biol 53(9):R61-109, 2008 —
//! porfimer absolute phi_so2 ≈ 0.65 in solution (community-anchored
//! across PDT literature).

use std::fmt;
use std::str::FromStr;

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
    /// Porfimer sodium (Photofrin). Single-exponential plasma decay,
    /// with optional pre-decay distribution phase and a relative
    /// singlet-O₂ quantum-yield scalar.
    ///
    /// `t_distribution_h` *approximates* the absorption / redistribution
    /// rise as instantaneous saturation: drug is held at peak for the
    /// first `t_distribution_h` hours after administration, then begins
    /// single-exponential decay with `t_half_h`. The real PK has a
    /// rising curve from 0 to peak (bi- or tri-exponential, formulation-
    /// dependent); this saturating-step model is a deliberate
    /// simplification per #203's design choice (issue body recommended
    /// the simpler subtract-before-decay path). For drugs where the
    /// rising-curve shape matters more than its area-under-the-curve
    /// approximation, a two-phase variant would be a separate issue.
    ///
    /// Default `0.0` recovers the previous "light at peak" behavior
    /// bit-exactly. Bellnier 2006 reports porfimer redistribution over
    /// ~24–48 h; default 0 is preserved for backwards compatibility.
    ///
    /// `phi_so2_relative` scales the per-photon ROS yield relative to a
    /// porfimer-equivalent baseline (`1.0` = porfimer). `Params::pdt_ros`
    /// is calibrated to porfimer at peak; new drug variants set this
    /// field to their absolute phi_so2 divided by porfimer's (~0.65).
    /// Default `1.0` recovers previous physics.
    ///
    /// Caveat: tissue phi_so2 can diverge from solution phi_so2 due to
    /// aggregation and microenvironment effects. See
    /// `simulations/calibration/parameter_provenance.md` (the
    /// "Photosensitizer pharmacokinetics" section) for the full
    /// citation chain and the tissue-vs-solution discussion.
    Porfimer {
        /// Plasma terminal half-life in hours. Must be strictly positive.
        t_half_h: f64,
        /// Hours from administration during which drug is held at peak
        /// before exponential decay begins. Models the absorption /
        /// redistribution rise; default 0 (instant peak, no rise).
        /// Must be finite and ≥ 0.
        #[serde(default)]
        t_distribution_h: f64,
        /// Singlet-O₂ quantum yield relative to porfimer (1.0). Scales
        /// `concentration_at` to give the per-photon ROS yield via
        /// `Photosensitizer::yield_at`. Default 1.0 (porfimer baseline).
        /// Must be finite and ≥ 0.
        #[serde(default = "porfimer_phi_default")]
        phi_so2_relative: f64,
    },
}

/// Serde default for `Porfimer::phi_so2_relative`. Returns 1.0 (the
/// porfimer-equivalent baseline) so legacy JSON without this field
/// deserializes to identity-preserving physics.
fn porfimer_phi_default() -> f64 {
    1.0
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
            Self::Porfimer {
                t_half_h,
                t_distribution_h,
                phi_so2_relative,
            } => {
                // Always emit all three positional fields so Display
                // round-trips unambiguously through `FromStr`. f64's
                // whole-number rendering keeps it tidy in the common
                // case (`porfimer=504,0,1`).
                write!(
                    f,
                    "porfimer={t_half_h},{t_distribution_h},{phi_so2_relative}"
                )
            }
        }
    }
}

impl FromStr for Photosensitizer {
    type Err = String;

    /// Parse a SPEC string into a `Photosensitizer`. Accepted forms
    /// (case-insensitive on the variant name):
    /// - `uniform` → `Uniform(1.0)`
    /// - `uniform=N` → `Uniform(N)`
    /// - `porfimer` → `Porfimer { t_half_h: 504.0, t_distribution_h: 0.0, phi_so2_relative: 1.0 }`
    /// - `porfimer=A` → `Porfimer { t_half_h: A, t_distribution_h: 0.0, phi_so2_relative: 1.0 }`
    /// - `porfimer=A,B` → `Porfimer { t_half_h: A, t_distribution_h: B, phi_so2_relative: 1.0 }`
    /// - `porfimer=A,B,C` → `Porfimer { t_half_h: A, t_distribution_h: B, phi_so2_relative: C }`
    ///
    /// Errors on unknown variant, unparseable number, empty positional
    /// value (`porfimer=,B`), or any value that fails
    /// [`Photosensitizer::validate`]. Round-trips with the `Display`
    /// impl: `format!("{ps}").parse::<Photosensitizer>() == Ok(ps)`
    /// for every valid `Photosensitizer`.
    ///
    /// Implementing `FromStr` (rather than a bespoke `from_cli_spec`)
    /// gives free clap integration via `#[arg(value_parser =
    /// clap::value_parser!(Photosensitizer))]` and ergonomic
    /// `"porfimer=504".parse()?` syntax for non-CLI consumers.
    fn from_str(s: &str) -> Result<Self, String> {
        let s = s.trim();
        let (name, value) = match s.split_once('=') {
            Some((n, v)) => (n.trim(), Some(v.trim())),
            None => (s, None),
        };
        // `eq_ignore_ascii_case` avoids allocating a lowercased String
        // for every parse. Match-on-lowercase would be tidier syntax
        // but allocates per call.
        let ps = if name.eq_ignore_ascii_case("uniform") {
            let c = match value {
                Some(v) => parse_f64_field("uniform=N", "N", v)?,
                None => 1.0,
            };
            Self::Uniform(c)
        } else if name.eq_ignore_ascii_case("porfimer") {
            let mut t_half_h = 504.0; // Bellnier 2006 terminal t½ in humans, hours.
            let mut t_distribution_h = 0.0;
            let mut phi_so2_relative = 1.0;
            if let Some(v) = value {
                let parts: Vec<&str> = v.split(',').map(str::trim).collect();
                if parts.len() > 3 {
                    return Err(format!(
                        "porfimer=t_half[,t_dist[,phi]]: expected 1-3 comma-separated values, got {} in {v:?}",
                        parts.len()
                    ));
                }
                if let Some(p) = parts.first() {
                    t_half_h = parse_f64_field("porfimer=t_half[,t_dist[,phi]]", "t_half", p)?;
                }
                if let Some(p) = parts.get(1) {
                    t_distribution_h =
                        parse_f64_field("porfimer=t_half[,t_dist[,phi]]", "t_dist", p)?;
                }
                if let Some(p) = parts.get(2) {
                    phi_so2_relative = parse_f64_field("porfimer=t_half[,t_dist[,phi]]", "phi", p)?;
                }
            }
            Self::Porfimer {
                t_half_h,
                t_distribution_h,
                phi_so2_relative,
            }
        } else {
            return Err(format!(
                "unknown photosensitizer {name:?}; expected one of: uniform, uniform=N, porfimer, porfimer=t_half[,t_dist[,phi]] (case-insensitive)"
            ));
        };
        ps.validate()?;
        Ok(ps)
    }
}

/// Parse a single positional f64 field with a uniform error format,
/// rejecting empty values explicitly so `porfimer=,36` gives a clear
/// message rather than the cryptic `parse::<f64>("")` error.
fn parse_f64_field(spec_form: &str, field_name: &str, raw: &str) -> Result<f64, String> {
    if raw.is_empty() {
        return Err(format!("{spec_form}: empty value for {field_name}"));
    }
    raw.parse::<f64>()
        .map_err(|e| format!("{spec_form}: cannot parse {field_name}={raw:?}: {e}"))
}

/// Validation helper: field must be finite and strictly > 0. Used for
/// `t_half_h` (a 0-half-life would mean infinite decay rate; `concentration_at`
/// would divide by zero).
fn check_finite_strict_positive(field: &str, value: f64) -> Result<(), String> {
    if !value.is_finite() {
        return Err(format!(
            "Photosensitizer::{field} must be finite, got {value}"
        ));
    }
    if value <= 0.0 {
        return Err(format!("Photosensitizer::{field} must be > 0, got {value}"));
    }
    Ok(())
}

/// Validation helper: field must be finite and ≥ 0. Used for the optional
/// `t_distribution_h` (zero is a sensible "no distribution phase" default)
/// and `phi_so2_relative` (zero means "drug present but emits no ROS",
/// a defensible thought-experiment use).
fn check_finite_nonneg(field: &str, value: f64) -> Result<(), String> {
    if !value.is_finite() {
        return Err(format!(
            "Photosensitizer::{field} must be finite, got {value}"
        ));
    }
    if value < 0.0 {
        return Err(format!(
            "Photosensitizer::{field} must be >= 0, got {value}"
        ));
    }
    Ok(())
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
            Self::Porfimer {
                t_half_h,
                t_distribution_h,
                phi_so2_relative,
            } => {
                check_finite_strict_positive("Porfimer.t_half_h", *t_half_h)?;
                check_finite_nonneg("Porfimer.t_distribution_h", *t_distribution_h)?;
                check_finite_nonneg("Porfimer.phi_so2_relative", *phi_so2_relative)?;
                Ok(())
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
            Self::Porfimer {
                t_half_h,
                t_distribution_h,
                ..
            } => {
                // Hold drug at peak during the absorption / redistribution
                // phase, then begin exponential decay. With default
                // `t_distribution_h = 0.0`, `(t - 0).max(0) == t.max(0)`
                // and the math is bit-exact identical to the pre-#203
                // single-exponential model.
                let effective_t = (t_h - t_distribution_h).max(0.0);
                let k = 2.0_f64.ln() / *t_half_h;
                (-k * effective_t).exp()
            }
        }
    }

    /// Per-photon ROS yield at time `t_h` post-administration: drug
    /// concentration scaled by the variant's singlet-O₂ quantum yield
    /// relative to porfimer.
    ///
    /// `Photosensitizer::Uniform` carries no quantum-yield model;
    /// `yield_at == concentration_at` (preserves the no-PK identity).
    ///
    /// `Photosensitizer::Porfimer { phi_so2_relative, .. }` returns
    /// `concentration_at(t_h) * phi_so2_relative`. The default
    /// `phi_so2_relative = 1.0` (porfimer baseline) keeps physics
    /// bit-exact identical to the pre-#203 model. New drug variants
    /// (Photochlor, 5-ALA-PpIX) would set this to their absolute
    /// phi_so2 divided by porfimer's (~0.65) so `Params::pdt_ros`'s
    /// porfimer-equivalent calibration carries through correctly.
    ///
    /// `pdt_intensity_at_depth` calls this rather than
    /// `concentration_at` so light × drug × yield composes correctly
    /// for inter-drug comparisons.
    pub fn yield_at(&self, t_h: f64) -> f64 {
        match self {
            Self::Uniform(_) => self.concentration_at(t_h),
            Self::Porfimer {
                phi_so2_relative, ..
            } => self.concentration_at(t_h) * phi_so2_relative,
        }
    }
}

/// Tissue pO₂ at the illuminated site under photochemical consumption,
/// as a fraction of its unilluminated value.
///
/// PDT CONSUMES THE OXYGEN IT NEEDS. A Type II photosensitizer works by
/// handing energy to ground-state O₂, so the illumination depletes its own
/// substrate at the site it is treating, and perfusion has to resupply it.
/// At quasi-steady state the two balance as `1 / (1 + phi/phi_crit)`, where
/// `phi_crit` is the fluence rate at which consumption equals resupply.
///
/// The consequence is the one PDT's own literature reports and no other arm
/// in this chapter has: delivering the SAME total light more slowly kills
/// MORE, because a lower fluence rate holds the tissue oxygenated (Henderson
/// & Busch, "Fluence rate as a modulator of PDT mechanisms", Lasers Surg Med
/// 2006, PMID 16615136, which reports O₂ depletion within seconds of
/// illumination at 75 mW/cm² in murine tumours).
///
/// `phi_crit <= 0` means "perfusion cannot resupply at all" and returns 0.0
/// for any illumination; an unilluminated site returns exactly 1.0.
pub fn oxygen_fraction_under_illumination(fluence_rate_mw_cm2: f64, phi_crit_mw_cm2: f64) -> f64 {
    let phi = fluence_rate_mw_cm2.max(0.0);
    if phi == 0.0 {
        return 1.0;
    }
    if !(phi_crit_mw_cm2 > 0.0) {
        return 0.0;
    }
    1.0 / (1.0 + phi / phi_crit_mw_cm2)
}

/// Relative singlet-oxygen yield per photon at a given fluence rate.
///
/// REUSES THE ENGINE'S MEASURED OXYGEN SHAPE rather than picking a second
/// one: [`crate::oxygen::oer_relative_efficacy`] is the Alper–Howard-Flanders
/// hyperbola #726 adopted after `oxygen_form_check.py` established the linear
/// form disagrees with eighty years of measurement below 10 mmHg. Radiation
/// and PDT are both oxygen-dependent modalities in this engine, and giving
/// them two different O₂ curves would let them drift apart for no reason
/// other than which module was written first.
pub fn fluence_rate_yield_factor(fluence_rate_mw_cm2: f64, phi_crit_mw_cm2: f64) -> f64 {
    let o2 = oxygen_fraction_under_illumination(fluence_rate_mw_cm2, phi_crit_mw_cm2);
    crate::oxygen::oer_relative_efficacy(o2, crate::oxygen::OER_REFERENCE_PO2_MMHG)
}

/// Reacted singlet oxygen delivered by one illumination, in arbitrary units.
///
/// Holds the TOTAL light fixed and varies only the rate at which it is
/// delivered, which is the comparison a clinician actually faces: the same
/// `total_fluence_j_cm2` over `total/rate` seconds. Two effects oppose across
/// that choice, and this is the fourth interior optimum in the chapter built
/// the same way:
///
///   * going FASTER depletes oxygen at the site, cutting the per-photon yield
///   * going SLOWER runs the illumination further into the sensitizer's own
///     pharmacokinetics, and drug that has cleared cannot be excited
///
/// The second term is not an assumption added for this function — it is
/// [`Photosensitizer::yield_at`], the module's existing PK, integrated over
/// the illumination window instead of sampled once at its start. Sampling
/// once is what the engine did before, and it is exactly the approximation
/// that makes a long illumination look free.
pub fn delivered_singlet_oxygen(
    ps: &Photosensitizer,
    dli_h: f64,
    total_fluence_j_cm2: f64,
    fluence_rate_mw_cm2: f64,
    phi_crit_mw_cm2: f64,
    steps: usize,
) -> f64 {
    let rate = fluence_rate_mw_cm2;
    if !(rate > 0.0) || !(total_fluence_j_cm2 > 0.0) || steps == 0 {
        return 0.0;
    }
    // J/cm2 divided by mW/cm2 is 1000 s; expressed in hours to match the PK.
    let duration_h = total_fluence_j_cm2 * 1000.0 / rate / 3600.0;
    let g = fluence_rate_yield_factor(rate, phi_crit_mw_cm2);
    // Midpoint rule over the illumination window. The integrand is the drug
    // present at each moment; the fluence rate is constant, so it factors out
    // and what is integrated is drug-hours.
    let dt = duration_h / steps as f64;
    let mut drug_hours = 0.0;
    for i in 0..steps {
        let t = dli_h + (i as f64 + 0.5) * dt;
        drug_hours += ps.yield_at(t) * dt;
    }
    // rate * drug_hours has units of (mW/cm2)*h; the constant 3600/1000
    // returns it to J/cm2 so the result is comparable across rates.
    g * rate * drug_hours * 3.6
}

/// The fluence rate maximising [`delivered_singlet_oxygen`], scanned.
///
/// Returns `(rate, delivered)`. NO CLOSED FORM IS OFFERED, because the drug
/// term is whatever [`Photosensitizer::yield_at`] is and a closed form would
/// silently assume the single-exponential variant.
///
/// A SCAN THAT RETURNS ITS OWN LOWER BOUND IS NOT AN OPTIMUM, and this
/// function does not hide that: when the sensitizer clears slowly compared
/// with the illumination, the answer is genuinely "as slow as you can go" and
/// the scan floor is returned. [`fluence_rate_optimum_is_interior`] is the
/// predicate a caller must consult before reporting a number, and the
/// degenerate case is the one the ablation and oncolytic sections both
/// shipped a wrong figure for before marking it.
pub fn optimal_fluence_rate(
    ps: &Photosensitizer,
    dli_h: f64,
    total_fluence_j_cm2: f64,
    phi_crit_mw_cm2: f64,
    lo_mw_cm2: f64,
    hi_mw_cm2: f64,
    scan_points: usize,
) -> (f64, f64) {
    if scan_points < 2 || !(hi_mw_cm2 > lo_mw_cm2) || !(lo_mw_cm2 > 0.0) {
        return (0.0, 0.0);
    }
    let mut best = (lo_mw_cm2, f64::NEG_INFINITY);
    for i in 0..scan_points {
        let frac = i as f64 / (scan_points - 1) as f64;
        // Logarithmic, because clinical fluence rates span two decades and a
        // linear scan spends almost all its points at the top.
        let rate = lo_mw_cm2 * (hi_mw_cm2 / lo_mw_cm2).powf(frac);
        let d = delivered_singlet_oxygen(ps, dli_h, total_fluence_j_cm2, rate, phi_crit_mw_cm2, 64);
        if d > best.1 {
            best = (rate, d);
        }
    }
    best
}

/// Whether the scanned optimum is genuinely interior, or sits on a scan edge.
///
/// An edge answer means the model is monotonic over the scanned range and the
/// real limit is something the model does not carry — the clinic's tolerance
/// for a two-hour illumination, at the slow end. Reporting an edge as an
/// optimum states a scan parameter as a finding.
pub fn fluence_rate_optimum_is_interior(rate: f64, lo_mw_cm2: f64, hi_mw_cm2: f64) -> bool {
    rate > lo_mw_cm2 * 1.001 && rate < hi_mw_cm2 * 0.999
}

/// Validate a drug-light-interval value (hours): reject NaN, negative,
/// and infinite inputs so they cannot reach `concentration_at` and
/// produce silent non-physical PDT output. Pair with the `FromStr` impl
/// for [`Photosensitizer`] at the same parse-time gate.
///
/// Error messages name only the constraint, not the field or CLI flag
/// (`"must be finite, got NaN"`). Callers are expected to prefix with
/// their own field/flag context — e.g. sim-spatial wraps as
/// `error: --dli-h: <message>` so the final user-facing line reads:
/// `error: --dli-h: must be finite, got NaN` (no stutter).
pub fn validate_dli_h(dli_h: f64) -> Result<(), String> {
    if !dli_h.is_finite() {
        return Err(format!("must be finite, got {dli_h}"));
    }
    if dli_h < 0.0 {
        return Err(format!("must be >= 0, got {dli_h}"));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── The fluence-rate layer (#831) ─────────────────────────────────

    #[test]
    fn an_unilluminated_site_keeps_all_its_oxygen() {
        assert_eq!(oxygen_fraction_under_illumination(0.0, 50.0), 1.0);
        // And a site perfusion cannot reach loses all of it under any light.
        assert_eq!(oxygen_fraction_under_illumination(10.0, 0.0), 0.0);
    }

    #[test]
    fn faster_illumination_depletes_more_oxygen() {
        // The whole premise. A mutation making the fraction rate-blind --
        // returning 1.0 whatever the rate -- fails here, and so does one
        // that inverts the sign.
        let a = oxygen_fraction_under_illumination(25.0, 50.0);
        let b = oxygen_fraction_under_illumination(200.0, 50.0);
        assert!(b < a, "expected {b} < {a}");
        // At the critical rate consumption exactly matches resupply, so half
        // the oxygen remains. That is what phi_crit MEANS, and pinning it
        // stops the parameter drifting into being an arbitrary scale.
        assert!((oxygen_fraction_under_illumination(50.0, 50.0) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn the_yield_uses_the_engines_measured_oxygen_shape() {
        // Not a second O2 curve. A mutation swapping in the linear
        // `o2_dependent_exo_factor` fails, because the hyperbola and the line
        // disagree by more than this tolerance away from the endpoints.
        let rate = 200.0;
        let o2 = oxygen_fraction_under_illumination(rate, 50.0);
        let expect =
            crate::oxygen::oer_relative_efficacy(o2, crate::oxygen::OER_REFERENCE_PO2_MMHG);
        assert!((fluence_rate_yield_factor(rate, 50.0) - expect).abs() < 1e-12);
        // No illumination means no depletion means no penalty, exactly.
        assert_eq!(fluence_rate_yield_factor(0.0, 50.0), 1.0);
    }

    #[test]
    fn the_delivered_dose_integrates_the_drug_and_does_not_sample_it_once() {
        // The defect this layer exists to fix: sampling `yield_at` at the
        // start of illumination makes a long treatment free. Under a
        // decaying sensitizer the integral must be STRICTLY BELOW the
        // start-sample, and the gap must widen as the illumination lengthens.
        let ps = Photosensitizer::Porfimer {
            t_half_h: 1.0,
            t_distribution_h: 0.0,
            phi_so2_relative: 1.0,
        };
        for &rate in &[20.0_f64, 200.0] {
            let d = delivered_singlet_oxygen(&ps, 0.0, 150.0, rate, 50.0, 256);
            let hours = 150.0 * 1000.0 / rate / 3600.0;
            // What a start-sampling model would report: drug held at t=0.
            let naive =
                fluence_rate_yield_factor(rate, 50.0) * rate * ps.yield_at(0.0) * hours * 3.6;
            assert!(
                d < naive,
                "rate {rate}: integral {d} not below start-sample {naive}"
            );
        }
        // The gap is larger for the slower, longer illumination.
        let slow = delivered_singlet_oxygen(&ps, 0.0, 150.0, 20.0, 50.0, 256);
        let slow_naive = fluence_rate_yield_factor(20.0, 50.0)
            * 20.0
            * 1.0
            * (150.0 * 1000.0 / 20.0 / 3600.0)
            * 3.6;
        let fast = delivered_singlet_oxygen(&ps, 0.0, 150.0, 200.0, 50.0, 256);
        let fast_naive = fluence_rate_yield_factor(200.0, 50.0)
            * 200.0
            * 1.0
            * (150.0 * 1000.0 / 200.0 / 3600.0)
            * 3.6;
        assert!(slow / slow_naive < fast / fast_naive);
    }

    #[test]
    fn the_optimum_is_interior_and_falls_as_the_drug_clears_more_slowly() {
        // The model's own prediction (P21), pinned so it cannot be tuned
        // away: a slow-clearing sensitizer wants a slower illumination.
        let (lo, hi) = (5.0, 400.0);
        let mut prev = f64::INFINITY;
        for &th in &[0.25_f64, 1.0, 4.0, 16.0, 48.0] {
            let ps = Photosensitizer::Porfimer {
                t_half_h: th,
                t_distribution_h: 0.0,
                phi_so2_relative: 1.0,
            };
            let (rate, dose) = optimal_fluence_rate(&ps, 0.0, 150.0, 50.0, lo, hi, 600);
            assert!(dose > 0.0);
            assert!(
                fluence_rate_optimum_is_interior(rate, lo, hi),
                "t_half {th}: optimum {rate} sits on a scan edge"
            );
            assert!(
                rate < prev,
                "t_half {th}: optimum {rate} did not fall below {prev}"
            );
            prev = rate;
        }
    }

    #[test]
    fn an_edge_answer_is_not_reported_as_an_optimum() {
        // The degenerate case, refused rather than rendered. A mutation
        // making the predicate unconditionally true fails here.
        assert!(!fluence_rate_optimum_is_interior(5.0, 5.0, 400.0));
        assert!(!fluence_rate_optimum_is_interior(400.0, 5.0, 400.0));
        assert!(fluence_rate_optimum_is_interior(100.0, 5.0, 400.0));
    }

    #[test]
    fn a_degenerate_illumination_delivers_nothing() {
        let ps = Photosensitizer::default();
        assert_eq!(
            delivered_singlet_oxygen(&ps, 0.0, 150.0, 0.0, 50.0, 64),
            0.0
        );
        assert_eq!(
            delivered_singlet_oxygen(&ps, 0.0, 0.0, 100.0, 50.0, 64),
            0.0
        );
        assert_eq!(
            delivered_singlet_oxygen(&ps, 0.0, 150.0, 100.0, 50.0, 0),
            0.0
        );
        assert_eq!(
            optimal_fluence_rate(&ps, 0.0, 150.0, 50.0, 5.0, 5.0, 10),
            (0.0, 0.0)
        );
    }

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
        let p = Photosensitizer::Porfimer {
            t_half_h: 504.0,
            t_distribution_h: 0.0,
            phi_so2_relative: 1.0,
        };
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
        let p = Photosensitizer::Porfimer {
            t_half_h: 504.0,
            t_distribution_h: 0.0,
            phi_so2_relative: 1.0,
        };
        assert_eq!(p.concentration_at(0.0), 1.0);
    }

    #[test]
    fn porfimer_decays_monotonically() {
        let p = Photosensitizer::Porfimer {
            t_half_h: 504.0,
            t_distribution_h: 0.0,
            phi_so2_relative: 1.0,
        };
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
        // Display always emits all three positional fields so the form
        // round-trips unambiguously through FromStr.
        assert_eq!(
            format!(
                "{}",
                Photosensitizer::Porfimer {
                    t_half_h: 504.0,
                    t_distribution_h: 0.0,
                    phi_so2_relative: 1.0
                }
            ),
            "porfimer=504,0,1"
        );
        assert_eq!(
            format!(
                "{}",
                Photosensitizer::Porfimer {
                    t_half_h: 336.5,
                    t_distribution_h: 36.0,
                    phi_so2_relative: 0.65
                }
            ),
            "porfimer=336.5,36,0.65"
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
        let p = Photosensitizer::Porfimer {
            t_half_h: 504.0,
            t_distribution_h: 0.0,
            phi_so2_relative: 1.0,
        };
        let json = serde_json::to_string(&p).unwrap();
        let q: Photosensitizer = serde_json::from_str(&json).unwrap();
        assert_eq!(p, q);
    }

    #[test]
    fn legacy_porfimer_json_deserializes_with_serde_defaults() {
        // JSON written before the t_distribution_h / phi_so2_relative
        // fields existed must still deserialize, with new fields filled
        // by `serde(default)` to identity-preserving values
        // (t_distribution_h=0, phi_so2_relative=1). This is the
        // backwards-compat guarantee the PR claims.
        let legacy_json = r#"{"Porfimer":{"t_half_h":504.0}}"#;
        let p: Photosensitizer = serde_json::from_str(legacy_json).unwrap();
        assert_eq!(
            p,
            Photosensitizer::Porfimer {
                t_half_h: 504.0,
                t_distribution_h: 0.0,
                phi_so2_relative: 1.0,
            }
        );
    }

    #[test]
    fn legacy_porfimer_json_with_partial_new_fields_deserializes() {
        // Mid-migration JSON that has t_distribution_h but not yet
        // phi_so2_relative — phi defaults to 1.0.
        let mid_json = r#"{"Porfimer":{"t_half_h":504.0,"t_distribution_h":36.0}}"#;
        let p: Photosensitizer = serde_json::from_str(mid_json).unwrap();
        assert_eq!(
            p,
            Photosensitizer::Porfimer {
                t_half_h: 504.0,
                t_distribution_h: 36.0,
                phi_so2_relative: 1.0,
            }
        );
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
        assert!(Photosensitizer::Porfimer {
            t_half_h: 504.0,
            t_distribution_h: 0.0,
            phi_so2_relative: 1.0
        }
        .validate()
        .is_ok());
        assert!(Photosensitizer::Porfimer {
            t_half_h: 0.001,
            t_distribution_h: 0.0,
            phi_so2_relative: 1.0
        }
        .validate()
        .is_ok());
    }

    #[test]
    fn validate_porfimer_rejects_invalid() {
        assert!(Photosensitizer::Porfimer {
            t_half_h: 0.0,
            t_distribution_h: 0.0,
            phi_so2_relative: 1.0
        }
        .validate()
        .is_err());
        assert!(Photosensitizer::Porfimer {
            t_half_h: -100.0,
            t_distribution_h: 0.0,
            phi_so2_relative: 1.0
        }
        .validate()
        .is_err());
        assert!(Photosensitizer::Porfimer {
            t_half_h: f64::NAN,
            t_distribution_h: 0.0,
            phi_so2_relative: 1.0
        }
        .validate()
        .is_err());
        assert!(Photosensitizer::Porfimer {
            t_half_h: f64::INFINITY,
            t_distribution_h: 0.0,
            phi_so2_relative: 1.0
        }
        .validate()
        .is_err());
    }

    #[test]
    fn negative_dli_saturates_to_peak() {
        // Negative DLI is non-physical (illumination before drug peak);
        // saturate to t=0 to keep outputs in [0, 1] for valid params.
        let p = Photosensitizer::Porfimer {
            t_half_h: 504.0,
            t_distribution_h: 0.0,
            phi_so2_relative: 1.0,
        };
        assert_eq!(p.concentration_at(-1.0), 1.0);
        assert_eq!(p.concentration_at(-1e6), 1.0);
        assert_eq!(p.concentration_at(f64::NEG_INFINITY), 1.0);
    }

    #[test]
    fn nan_dli_saturates_to_peak() {
        // f64::max treats NaN as smaller than any value, so NaN.max(0.0)
        // returns 0.0 and concentration_at(NaN) is deterministic.
        let p = Photosensitizer::Porfimer {
            t_half_h: 504.0,
            t_distribution_h: 0.0,
            phi_so2_relative: 1.0,
        };
        assert_eq!(p.concentration_at(f64::NAN), 1.0);
    }

    // -- FromStr --

    #[test]
    fn parse_uniform_default() {
        assert_eq!(
            "uniform".parse::<Photosensitizer>().unwrap(),
            Photosensitizer::Uniform(1.0)
        );
    }

    #[test]
    fn parse_uniform_with_value() {
        assert_eq!(
            "uniform=0.5".parse::<Photosensitizer>().unwrap(),
            Photosensitizer::Uniform(0.5)
        );
    }

    #[test]
    fn parse_porfimer_default_is_bellnier() {
        assert_eq!(
            "porfimer".parse::<Photosensitizer>().unwrap(),
            Photosensitizer::Porfimer {
                t_half_h: 504.0,
                t_distribution_h: 0.0,
                phi_so2_relative: 1.0
            }
        );
    }

    #[test]
    fn parse_porfimer_with_value() {
        assert_eq!(
            "porfimer=336".parse::<Photosensitizer>().unwrap(),
            Photosensitizer::Porfimer {
                t_half_h: 336.0,
                t_distribution_h: 0.0,
                phi_so2_relative: 1.0
            }
        );
    }

    #[test]
    fn parse_trims_whitespace() {
        assert_eq!(
            "  porfimer = 504  ".parse::<Photosensitizer>().unwrap(),
            Photosensitizer::Porfimer {
                t_half_h: 504.0,
                t_distribution_h: 0.0,
                phi_so2_relative: 1.0
            }
        );
    }

    #[test]
    fn parse_case_insensitive() {
        assert_eq!(
            "Uniform".parse::<Photosensitizer>().unwrap(),
            Photosensitizer::Uniform(1.0)
        );
        assert_eq!(
            "PORFIMER=504".parse::<Photosensitizer>().unwrap(),
            Photosensitizer::Porfimer {
                t_half_h: 504.0,
                t_distribution_h: 0.0,
                phi_so2_relative: 1.0
            }
        );
        assert_eq!(
            "PorFimEr".parse::<Photosensitizer>().unwrap(),
            Photosensitizer::Porfimer {
                t_half_h: 504.0,
                t_distribution_h: 0.0,
                phi_so2_relative: 1.0
            }
        );
    }

    #[test]
    fn parse_unknown_variant_errors() {
        let err = "photochlor".parse::<Photosensitizer>().unwrap_err();
        assert!(err.contains("photochlor"));
        assert!(err.contains("uniform"));
        assert!(err.contains("porfimer"));
    }

    #[test]
    fn parse_unparseable_number_errors() {
        let err = "porfimer=abc".parse::<Photosensitizer>().unwrap_err();
        assert!(err.contains("porfimer=t_half"));
        assert!(err.contains("abc"));
    }

    #[test]
    fn parse_negative_t_half_h_rejected() {
        let err = "porfimer=-1".parse::<Photosensitizer>().unwrap_err();
        assert!(err.contains("t_half_h"));
    }

    #[test]
    fn parse_zero_t_half_h_rejected() {
        let err = "porfimer=0".parse::<Photosensitizer>().unwrap_err();
        assert!(err.contains("t_half_h"));
    }

    #[test]
    fn parse_negative_uniform_rejected() {
        let err = "uniform=-0.5".parse::<Photosensitizer>().unwrap_err();
        assert!(err.contains("must be >= 0"));
    }

    #[test]
    fn parse_nan_rejected() {
        let err = "uniform=NaN".parse::<Photosensitizer>().unwrap_err();
        assert!(err.contains("must be finite"));
    }

    #[test]
    fn display_round_trips_through_parse() {
        // Display's contract is round-trip parseability via FromStr.
        // Includes Uniform(0.0) as the boundary case where validation
        // accepts zero but a regression in either direction (parse or
        // Display) would silently break the contract.
        for ps in [
            Photosensitizer::Uniform(0.0),
            Photosensitizer::Uniform(1.0),
            Photosensitizer::Uniform(0.5),
            Photosensitizer::Porfimer {
                t_half_h: 504.0,
                t_distribution_h: 0.0,
                phi_so2_relative: 1.0,
            },
            Photosensitizer::Porfimer {
                t_half_h: 336.5,
                t_distribution_h: 0.0,
                phi_so2_relative: 1.0,
            },
        ] {
            let rendered = format!("{ps}");
            let reparsed: Photosensitizer = rendered
                .parse()
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
        // Library returns only the constraint ("must be >= 0, got -1");
        // no field/flag name embedded. Callers prefix with their own
        // context, e.g. sim-spatial wraps as `error: --dli-h: <msg>`.
        assert!(err.contains("must be >= 0"));
        assert!(err.contains("-1"));
    }

    #[test]
    fn validate_dli_h_rejects_nan() {
        let err = validate_dli_h(f64::NAN).unwrap_err();
        assert!(err.contains("finite"));
    }

    #[test]
    fn validate_dli_h_rejects_infinity() {
        assert!(validate_dli_h(f64::INFINITY)
            .unwrap_err()
            .contains("finite"));
        assert!(validate_dli_h(f64::NEG_INFINITY)
            .unwrap_err()
            .contains("finite"));
    }

    // -- distribution-phase (#203 A) --

    #[test]
    fn porfimer_holds_at_peak_during_distribution() {
        // Drug is held at peak for `t_distribution_h` hours, then begins
        // exponential decay. So at any time t < t_distribution_h, the
        // concentration is exactly 1.0.
        let p = Photosensitizer::Porfimer {
            t_half_h: 504.0,
            t_distribution_h: 36.0,
            phi_so2_relative: 1.0,
        };
        assert_eq!(p.concentration_at(0.0), 1.0);
        assert_eq!(p.concentration_at(12.0), 1.0);
        assert_eq!(p.concentration_at(35.999), 1.0);
        assert_eq!(p.concentration_at(36.0), 1.0); // boundary: effective_t = 0
    }

    #[test]
    fn porfimer_decays_after_distribution_phase() {
        // After t_distribution_h, decay measures from t = t_distribution_h,
        // not from t = 0. So at t = t_distribution_h + t_half_h, drug is
        // at half peak — same as the no-distribution case at t = t_half_h.
        let with_dist = Photosensitizer::Porfimer {
            t_half_h: 504.0,
            t_distribution_h: 36.0,
            phi_so2_relative: 1.0,
        };
        let without_dist = Photosensitizer::Porfimer {
            t_half_h: 504.0,
            t_distribution_h: 0.0,
            phi_so2_relative: 1.0,
        };
        // exp(-ln(2)) is libm-dependent; use epsilon comparison.
        let c_with = with_dist.concentration_at(36.0 + 504.0);
        let c_without = without_dist.concentration_at(504.0);
        assert!((c_with - 0.5).abs() < 1e-12);
        assert!((c_without - 0.5).abs() < 1e-12);
        assert_eq!(c_with, c_without); // both pinned to one half-life past their respective t=0
    }

    #[test]
    fn porfimer_default_t_distribution_zero_preserves_legacy_math() {
        // The whole point of `t_distribution_h: 0.0` as the default is
        // bit-identical reproduction of the pre-#203 single-exponential
        // model. Verify by replicating the EXACT operation order
        // concentration_at uses: `k = ln(2) / t_half`, then `(-k * t).exp()`.
        // (Computing `(-ln(2) * t / t_half).exp()` differs by 1 ULP — float
        // arithmetic is non-associative; this is what got us a flaky test
        // the first time.)
        let p = Photosensitizer::Porfimer {
            t_half_h: 504.0,
            t_distribution_h: 0.0,
            phi_so2_relative: 1.0,
        };
        let k = 2.0_f64.ln() / 504.0;
        for t in [0.0, 1.0, 24.0, 168.0, 504.0, 1008.0] {
            let expected = (-k * t).exp();
            assert_eq!(
                p.concentration_at(t),
                expected,
                "concentration_at({t}) drift vs legacy math"
            );
        }
    }

    // -- phi_so2 / yield_at (#203 B) --

    #[test]
    fn yield_at_scales_concentration_by_phi_for_porfimer() {
        // yield_at(t) = concentration_at(t) * phi_so2_relative.
        let half_phi = Photosensitizer::Porfimer {
            t_half_h: 504.0,
            t_distribution_h: 0.0,
            phi_so2_relative: 0.5,
        };
        let unit_phi = Photosensitizer::Porfimer {
            t_half_h: 504.0,
            t_distribution_h: 0.0,
            phi_so2_relative: 1.0,
        };
        for t in [0.0, 24.0, 168.0, 504.0] {
            assert_eq!(
                half_phi.yield_at(t),
                unit_phi.yield_at(t) * 0.5,
                "phi=0.5 should halve yield at t={t}"
            );
        }
    }

    #[test]
    fn yield_at_uniform_matches_concentration_at() {
        // Uniform carries no phi (it's the no-PK identity); yield_at
        // delegates to concentration_at. This is what preserves the
        // default-invocation byte-identical snapshot.
        for c in [0.0, 0.5, 1.0, 2.0] {
            let u = Photosensitizer::Uniform(c);
            for t in [0.0, 24.0, 1e9] {
                assert_eq!(u.yield_at(t), u.concentration_at(t));
                assert_eq!(u.yield_at(t), c);
            }
        }
    }

    #[test]
    fn default_porfimer_yield_equals_concentration() {
        // Default Porfimer phi=1.0 means yield_at == concentration_at
        // exactly — the IEEE invariant `x * 1.0 = x` for finite x.
        // This is the byte-identical guarantee: anywhere physics calls
        // yield_at instead of concentration_at, default-photosensitizer
        // results don't change.
        let p = Photosensitizer::Porfimer {
            t_half_h: 504.0,
            t_distribution_h: 0.0,
            phi_so2_relative: 1.0,
        };
        for t in [0.0, 24.0, 168.0, 504.0, 1008.0] {
            assert_eq!(p.yield_at(t), p.concentration_at(t));
        }
    }

    // -- validate (new fields) --

    #[test]
    fn validate_porfimer_rejects_negative_t_distribution() {
        let err = Photosensitizer::Porfimer {
            t_half_h: 504.0,
            t_distribution_h: -1.0,
            phi_so2_relative: 1.0,
        }
        .validate()
        .unwrap_err();
        assert!(err.contains("t_distribution_h"));
        assert!(err.contains(">= 0"));
    }

    #[test]
    fn validate_porfimer_rejects_nan_t_distribution() {
        let err = Photosensitizer::Porfimer {
            t_half_h: 504.0,
            t_distribution_h: f64::NAN,
            phi_so2_relative: 1.0,
        }
        .validate()
        .unwrap_err();
        assert!(err.contains("t_distribution_h"));
        assert!(err.contains("finite"));
    }

    #[test]
    fn validate_porfimer_rejects_negative_phi() {
        let err = Photosensitizer::Porfimer {
            t_half_h: 504.0,
            t_distribution_h: 0.0,
            phi_so2_relative: -0.5,
        }
        .validate()
        .unwrap_err();
        assert!(err.contains("phi_so2_relative"));
        assert!(err.contains(">= 0"));
    }

    #[test]
    fn validate_porfimer_rejects_nan_phi() {
        let err = Photosensitizer::Porfimer {
            t_half_h: 504.0,
            t_distribution_h: 0.0,
            phi_so2_relative: f64::NAN,
        }
        .validate()
        .unwrap_err();
        assert!(err.contains("phi_so2_relative"));
        assert!(err.contains("finite"));
    }

    #[test]
    fn validate_porfimer_accepts_zero_phi_and_zero_distribution() {
        // Edge: phi=0 means drug emits no ROS yield (thought-experiment
        // / quenched-fluorophore use). t_dist=0 is the default. Both
        // are valid even at the boundary.
        assert!(Photosensitizer::Porfimer {
            t_half_h: 504.0,
            t_distribution_h: 0.0,
            phi_so2_relative: 0.0,
        }
        .validate()
        .is_ok());
    }

    // -- Parser: triple positional + legacy --

    #[test]
    fn parse_porfimer_triple_full() {
        assert_eq!(
            "porfimer=504,36,0.65".parse::<Photosensitizer>().unwrap(),
            Photosensitizer::Porfimer {
                t_half_h: 504.0,
                t_distribution_h: 36.0,
                phi_so2_relative: 0.65,
            }
        );
    }

    #[test]
    fn parse_porfimer_double_defaults_phi() {
        assert_eq!(
            "porfimer=504,36".parse::<Photosensitizer>().unwrap(),
            Photosensitizer::Porfimer {
                t_half_h: 504.0,
                t_distribution_h: 36.0,
                phi_so2_relative: 1.0,
            }
        );
    }

    #[test]
    fn parse_porfimer_single_legacy_form_defaults_new_fields() {
        // Backwards compat: existing `porfimer=504` parses to defaults
        // for the new fields, so legacy CLI invocations and serialized
        // configs continue to work without modification.
        assert_eq!(
            "porfimer=504".parse::<Photosensitizer>().unwrap(),
            Photosensitizer::Porfimer {
                t_half_h: 504.0,
                t_distribution_h: 0.0,
                phi_so2_relative: 1.0,
            }
        );
    }

    #[test]
    fn parse_porfimer_rejects_too_many_positional() {
        let err = "porfimer=504,36,0.65,extra"
            .parse::<Photosensitizer>()
            .unwrap_err();
        assert!(err.contains("expected 1-3"));
    }

    #[test]
    fn parse_porfimer_rejects_empty_positional() {
        let err = "porfimer=,36".parse::<Photosensitizer>().unwrap_err();
        assert!(err.contains("empty value"));
        assert!(err.contains("t_half"));
    }

    #[test]
    fn parse_porfimer_rejects_trailing_equals_with_no_value() {
        // `porfimer=` is a different shape from `porfimer=,36` — there's
        // no comma at all, so the value side is the single empty string.
        // Should give a clear empty-value error rather than the cryptic
        // f64::parse("") error or silently defaulting.
        let err = "porfimer=".parse::<Photosensitizer>().unwrap_err();
        assert!(err.contains("empty value"));
        assert!(err.contains("t_half"));
    }

    #[test]
    fn display_round_trips_with_new_fields() {
        // Round-trip property must hold for non-default new fields too,
        // including extreme f64 values where Display switches to
        // scientific notation (e.g. `0.0000001` → `0.0000001`,
        // `1e-10` → `0.0000000001` or `1e-10` depending on rustc).
        // FromStr accepts both decimal and scientific forms, so the
        // contract holds either way — but pinning these cases catches
        // regressions in either Display or parser.
        for ps in [
            Photosensitizer::Porfimer {
                t_half_h: 504.0,
                t_distribution_h: 36.0,
                phi_so2_relative: 0.65,
            },
            Photosensitizer::Porfimer {
                t_half_h: 250.0,
                t_distribution_h: 24.0,
                phi_so2_relative: 1.5, // > 1 enrichment hook
            },
            Photosensitizer::Porfimer {
                t_half_h: 504.0,
                t_distribution_h: 0.0,
                phi_so2_relative: 0.0, // boundary: zero-yield
            },
            // Extreme small / large values that may render in scientific
            // notation. FromStr handles both `1e10` and `10000000000`.
            Photosensitizer::Porfimer {
                t_half_h: 1e10,
                t_distribution_h: 1e-10,
                phi_so2_relative: 1e-15,
            },
        ] {
            let rendered = format!("{ps}");
            let reparsed: Photosensitizer = rendered
                .parse()
                .unwrap_or_else(|e| panic!("round-trip failed for {ps:?}: {e}"));
            assert_eq!(reparsed, ps, "round-trip mismatch via {rendered:?}");
        }
    }
}
