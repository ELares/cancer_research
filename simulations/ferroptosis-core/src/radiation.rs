//! Ionizing radiation as a treatment arm (#726, the half the title asks for).
//!
//! The oxygen half of #726 landed first: [`crate::oxygen::oxygen_enhancement_ratio`]
//! gave the engine's exogenous-ROS O₂ dependence the measured
//! Alper–Howard-Flanders shape. This module is the modality that consumes it.
//!
//! ## Radiation is TWO channels, and collapsing them would be a lie
//!
//! The obvious implementation — feed dose in as exogenous ROS and let the
//! existing lipid-peroxidation engine do the rest — models radiation as a
//! photodynamic agent. Radiation's dominant lethal lesion is the DNA
//! double-strand break, not membrane lipid peroxidation, and a model that
//! routes all of it through the ferroptosis machinery would predict that a
//! GPX4-competent cell is radioresistant. It is not.
//!
//! So this module supplies two channels and the consumer composes them:
//!
//! 1. **DNA-damage lethality** — [`lq_survival`], the linear-quadratic model
//!    `SF = exp(−αD − βD²)`. It does not touch [`crate::biochem::CellState`] at
//!    all: it is a per-cell survival probability the consumer rolls against.
//!    This is the channel with eighty years of clonogenic-survival data behind
//!    it, and it is the calibration target this layer arrives with.
//! 2. **Radiation-induced ferroptosis** — [`exo_ros_from_dose`], which enters
//!    through the same `exo_ros_peak` path SDT and PDT already use. Radiation
//!    does induce lipid peroxidation and ferroptosis, and lipid-metabolism
//!    manipulation does radiosensitize. The direction is anchored inside this
//!    project's own frozen corpus rather than by a citation carried in from
//!    outside -- PMID 35859734 is titled on exactly this coupling, SCD
//!    targeting enhancing radiation-induced ferroptosis and immunogenic cell
//!    death. No COUNT is quoted here: a number in a doc comment is a number
//!    nothing can check, so the size of that literature belongs in a derived
//!    analysis artifact, not in Rust prose. What the literature does NOT give
//!    is the ROS-per-gray conversion, so `ros_per_gy` is an explicitly
//!    uncalibrated knob and its CALIBRATION_STATUS row says so.
//!
//! ## Oxygen enters both channels through the same hyperbola
//!
//! Radiobiology's oxygen effect is a DOSE-MODIFYING factor: a hypoxic cell
//! needs OER-times more dose for the same survival. [`oer_scaled_dose`]
//! expresses exactly that, reusing [`crate::oxygen::oer_relative_efficacy`], so
//! the O₂ shape in this module is the same measured hyperbola the SDT path
//! already validated rather than a second, independent assumption.
//!
//! That is also what makes this arm testable rather than merely present. Fit
//! `ros_per_gy` at full O₂ and the OER is then PREDICTED, not fitted: the
//! iso-survival dose ratio between anoxia and full supply should land in the
//! published 2.5–3.0 band. If it does not, this layer is wrong and the report
//! must say so.
//!
//! ## Radiation is the physical modality with no depth limit
//!
//! [`crate::physics`] gives SDT and PDT exponential attenuation with depth, and
//! that attenuation IS this project's tissue-access argument: light dies in
//! millimetres, ultrasound reaches centimetres. Megavoltage photons attenuate
//! at roughly 3 % per centimetre of soft tissue, so over a whole tumour the
//! delivered dose is nearly flat. [`intensity_at_depth`] supplies that, which
//! gives the manuscript's depth comparison a control it has never had — the
//! one physical modality whose reach is not the limiting factor.
//!
//! ## Identity default
//!
//! [`RadiationConfig::default`] has `dose_gy = 0.0`, at which [`lq_survival`]
//! returns exactly `1.0`, [`dna_lethality`] exactly `0.0` and
//! [`exo_ros_from_dose`] exactly `0.0`. A run carrying the default config is
//! byte-identical to one with no radiation model at all.

use crate::oxygen::oer_relative_efficacy;
use crate::params::RadiationConfig;

/// Linear attenuation coefficient of 6 MV photons in soft tissue, per cm.
///
/// Roughly 3 % per centimetre. Quoted as a CONVENTIONAL value for
/// megavoltage photons in water-equivalent tissue, not as a fit performed
/// here; the point it carries is the ORDER, which is what separates radiation
/// from light and ultrasound rather than any third digit.
///
/// For contrast, `crate::physics::pdt_intensity_at_depth` falls by orders of
/// magnitude over the same distance. That contrast is the reason this constant
/// exists.
pub const MU_6MV_SOFT_TISSUE_PER_CM: f64 = 0.03;

/// α/β ratio for early-responding tumour tissue, in Gy.
///
/// The textbook split is ~10 Gy for tumours and early-responding normal
/// tissue against ~3 Gy for late-responding tissue. This project's frozen
/// corpus contains a worked use of it: PMID 32307022
/// (`corpus/by-pmid/32307022.md`) characterises its fractionation sensitivity
/// with `α/β = 10 Gy` alongside `SF2 = 0.5445`.
///
/// **What that paper is, stated because the first version of this comment
/// oversold it.** It is a seven-patient TTFields + IMRT dosimetry study, not a
/// clonogenic assay. Its SF2 is an EUD-model INPUT chosen to give a 30 %
/// tumour control probability (with the clonogen count from Suit 1992), and
/// it says so: "SF2 **was assumed** to be 0.5445". So this is a published
/// PARAMETERISATION that the engine can be checked against, not a measurement
/// the engine can be fitted to.
pub const ALPHA_BETA_TUMOUR_GY: f64 = 10.0;

/// α implied by PMID 32307022's `SF2 = 0.5445` at `α/β = 10`, per Gy.
///
/// `α = −ln(SF2) / (D + D²/(α/β))` at `D = 2`; β follows from the ratio.
pub const ALPHA_GBM_PARAMETERISATION_PER_GY: f64 = 0.25329;

/// The third number the same paper reports, `SF = 0.4244` at `2.67 Gy`.
///
/// **NOT an independent prediction, and calling it one was wrong.** The paper
/// COMPUTED it from SF2 and α/β through the LQ model — "the validity of the
/// linear quadratic cell survival model was assumed and the respective
/// surviving fraction of 0.4244 at 2.67 Gy was obtained" — which is the same
/// operation this module performs. Under LQ, `(SF2, α/β)` determines it.
///
/// What it is good for is a ROUND TRIP: reproducing it to 1e-4 shows
/// [`lq_survival`] implements the published parameterisation exactly, and a
/// sign error, a dropped β term or the late-responding 3 Gy convention all
/// move it off. That is an implementation check, not a test of the biology.
pub const SF_AT_2_67_GY_PUBLISHED: f64 = 0.4244;

/// Linear-quadratic surviving fraction, `exp(−αD − βD²)`.
///
/// The standard dose-response of radiobiology, and the reason radiation has an
/// external calibration literature that nothing else in this engine has:
/// clonogenic survival curves are measured per cell line and published with
/// their α and β.
///
/// Returns exactly `1.0` at `dose_gy == 0.0`, which is what makes the whole
/// layer inert by default.
///
/// ## Why not fold this into the ROS path
///
/// Because the two channels answer different questions and the interesting
/// result is how they separate. See the module docs.
#[must_use = "the surviving fraction is the function's only output"]
pub fn lq_survival(dose_gy: f64, alpha_per_gy: f64, beta_per_gy2: f64) -> f64 {
    debug_assert!(
        dose_gy >= 0.0 && dose_gy.is_finite(),
        "lq_survival: dose_gy must be finite and non-negative, got {dose_gy}"
    );
    debug_assert!(
        alpha_per_gy >= 0.0 && beta_per_gy2 >= 0.0,
        "lq_survival: alpha and beta must be non-negative, got alpha={alpha_per_gy}, beta={beta_per_gy2}"
    );
    let d = dose_gy.max(0.0);
    (-(alpha_per_gy.max(0.0) * d + beta_per_gy2.max(0.0) * d * d))
        .exp()
        .clamp(0.0, 1.0)
}

/// Dose scaled by the oxygen effect: what `dose_gy` is WORTH at `o2_supply`.
///
/// The oxygen effect is dose-modifying, so a hypoxic cell behaves as though it
/// received less dose. This multiplies by
/// [`crate::oxygen::oer_relative_efficacy`], the same normalised hyperbola the
/// SDT path uses, so the two modalities cannot drift apart in their O₂ shape.
///
/// `dependence ∈ [0, 1]` scales how much of the effect applies, exactly as
/// `crate::oxygen::oer_exo_factor`'s does: at `0.0` this returns `dose_gy`
/// unchanged, which is the O₂-independent limit and is the identity path.
#[must_use = "the scaled dose is the function's only output"]
pub fn oer_scaled_dose(dose_gy: f64, o2_supply: f64, dependence: f64, p_full_mmhg: f64) -> f64 {
    let dep = dependence.clamp(0.0, 1.0);
    let eff = oer_relative_efficacy(o2_supply, p_full_mmhg);
    dose_gy.max(0.0) * (1.0 - dep + dep * eff)
}

/// Relative delivered dose at depth `z_um`, `exp(−μ z)`.
///
/// Near-flat over a tumour: at `MU_6MV_SOFT_TISSUE_PER_CM` a 10 cm depth still
/// receives ~74 % of the surface dose. Compare `crate::physics`'s PDT and SDT
/// attenuation over the same distance.
///
/// Negative depths (outside the tumour) clamp to `0`, matching the sign
/// convention of `crate::grid::TumorGrid3D::radial_depth_um`.
#[must_use = "the intensity fraction is the function's only output"]
pub fn intensity_at_depth(z_um: f64, mu_per_cm: f64) -> f64 {
    debug_assert!(
        mu_per_cm >= 0.0 && mu_per_cm.is_finite(),
        "intensity_at_depth: mu_per_cm must be finite and non-negative, got {mu_per_cm}"
    );
    let z_cm = (z_um.max(0.0)) / 10_000.0;
    (-mu_per_cm.max(0.0) * z_cm).exp().clamp(0.0, 1.0)
}

/// Probability this cell is killed by DNA damage alone, `1 − SF`.
///
/// The channel that does NOT pass through the ferroptosis state. A consumer
/// rolls against it once, at the moment of delivery.
#[must_use = "the lethality is the function's only output"]
pub fn dna_lethality(cfg: &RadiationConfig, o2_supply: f64, z_um: f64) -> f64 {
    let d = delivered_dose(cfg, o2_supply, z_um);
    1.0 - lq_survival(d, cfg.alpha_per_gy, cfg.beta_per_gy2)
}

/// Exogenous-ROS peak this cell receives, for the ferroptosis channel.
///
/// `ros_per_gy` is UNCALIBRATED and the CALIBRATION_STATUS row says so: the
/// radiation-ferroptosis literature establishes the direction (radiation
/// induces lipid peroxidation; GPX4 inhibition radiosensitizes) and reports no
/// conversion from gray to an intracellular ROS peak. At the default
/// `ros_per_gy = 0.0` this returns `0.0`, so the ferroptosis channel is off
/// unless a caller opts in.
#[must_use = "the ROS peak is the function's only output"]
pub fn exo_ros_from_dose(cfg: &RadiationConfig, o2_supply: f64, z_um: f64) -> f64 {
    cfg.ros_per_gy.max(0.0) * delivered_dose(cfg, o2_supply, z_um)
}

/// Dose actually delivered to a cell: prescribed, attenuated, O₂-scaled.
#[must_use = "the delivered dose is the function's only output"]
pub fn delivered_dose(cfg: &RadiationConfig, o2_supply: f64, z_um: f64) -> f64 {
    let attenuated = cfg.dose_gy.max(0.0) * intensity_at_depth(z_um, cfg.mu_per_cm);
    oer_scaled_dose(attenuated, o2_supply, cfg.o2_dependence, cfg.p_full_mmhg)
}

/// Published sensitizer-enhancement-ratio band for PARP inhibition in glioma.
///
/// 1.2–1.7 (PMID 35205750, `corpus/by-pmid/35205750.md`): "The base excision
/// repair pathway is essential to repair damaged bases caused by ionizing
/// radiation. The pathway is blocked by PARP inhibitors, which result in
/// sensitizer enhancement ratios of 1.2–1.7 in glioma cells."
///
/// **This is a real calibration target and a real failable prediction**, which
/// is what the linear-quadratic α/β round trip is NOT. `parp_alpha_boost` is
/// one number; the SER it produces is DOSE-DEPENDENT under LQ, because raising
/// α alone shifts the linear term while β is unchanged. So a boost fitted at
/// one dose predicts the SER at every other dose, and the model is free to put
/// that outside the published band.
pub const PARP_SER_BAND: (f64, f64) = (1.2, 1.7);

/// α after PARP inhibition: `α · (1 + boost)`, β unchanged.
///
/// The single-hit term rises because unrepaired single-strand breaks convert
/// to double-strand breaks at replication — one-track damage, not pairwise —
/// so β, which models pairwise lesion interaction, is left alone. A model that
/// scaled both would be asserting a mechanism nobody measured.
#[must_use = "the sensitized alpha is the function's only output"]
pub fn sensitized_alpha(cfg: &RadiationConfig) -> f64 {
    cfg.alpha_per_gy * (1.0 + cfg.parp_alpha_boost.max(0.0))
}

/// Sensitizer enhancement ratio at a given surviving fraction.
///
/// The dose WITHOUT the drug divided by the dose WITH it, at equal survival —
/// the quantity the radiobiology literature actually reports, so the model can
/// be compared to [`PARP_SER_BAND`] rather than to itself.
///
/// Solved in closed form: `SF = exp(-αD - βD²)` inverts to
/// `D = (-α + sqrt(α² + 4βL)) / (2β)` with `L = -ln(SF)`; at `β = 0` it
/// degenerates to `L/α`, and there the SER is exactly `1 + boost` at every
/// dose, which is the tell that the DOSE-DEPENDENCE comes from β.
///
/// Returns `1.0` when no boost is configured, so an unsensitized run reports
/// no enhancement rather than a near-one artefact.
#[must_use = "the ratio is the function's only output"]
pub fn sensitizer_enhancement_ratio(cfg: &RadiationConfig, surviving_fraction: f64) -> f64 {
    debug_assert!(
        surviving_fraction > 0.0 && surviving_fraction < 1.0,
        "sensitizer_enhancement_ratio: surviving_fraction must be in (0, 1), got {surviving_fraction}"
    );
    let boost = cfg.parp_alpha_boost.max(0.0);
    if boost == 0.0 {
        return 1.0;
    }
    let l = -surviving_fraction.ln();
    let dose = |a: f64, b: f64| -> f64 {
        if b <= 0.0 {
            l / a
        } else {
            (-a + (a * a + 4.0 * b * l).sqrt()) / (2.0 * b)
        }
    };
    let plain = dose(cfg.alpha_per_gy, cfg.beta_per_gy2);
    let sens = dose(sensitized_alpha(cfg), cfg.beta_per_gy2);
    if sens <= 0.0 {
        return f64::INFINITY;
    }
    plain / sens
}

/// PARP-inhibitor MONOTHERAPY lethality — synthetic lethality proper.
///
/// The defining property of a synthetic-lethal pair is that neither hit alone
/// is lethal: PARP inhibition kills HR-DEFICIENT cells and spares
/// HR-proficient ones, which is why these drugs ship with a biomarker. So this
/// is the PRODUCT of the two conditions and returns exactly `0.0` when either
/// is absent, rather than a small number that would make the biomarker look
/// optional.
///
/// `hr_deficiency` and `parp_alpha_boost` both default to `0.0`, so an
/// unconfigured run gets no monotherapy kill at all.
#[must_use = "the lethality is the function's only output"]
pub fn parp_monotherapy_lethality(cfg: &RadiationConfig) -> f64 {
    let hrd = cfg.hr_deficiency.clamp(0.0, 1.0);
    let inhibition = (cfg.parp_alpha_boost.max(0.0)).min(1.0);
    (hrd * inhibition).clamp(0.0, 1.0)
}

/// The dose-modifying factor the O₂ scaling applies — a RESTATEMENT of the
/// OER, not a prediction of it.
///
/// Named and documented this way deliberately. [`oer_scaled_dose`] multiplies
/// dose by the normalised hyperbola, so the ratio between two oxygenations
/// comes back out as the hyperbola again: at full dependence and anoxia it is
/// `m(p_full)/m(0)`, about 2.86 for `p_full = 40` mmHg. Reporting that number
/// against the published 2.5–3.0 band would be a guard computing its own
/// expectation, and this module refuses to do it.
///
/// **The layer's actual independent prediction lives one level up.** The
/// ferroptosis channel carries its own O₂ dependence through
/// `crate::oxygen::oer_exo_factor`, and the DNA channel carries this one, so
/// the OER a full run EXHIBITS — the dose ratio at equal simulated survival
/// with both channels active — is an emergent quantity that neither input
/// fixes. Fit `ros_per_gy` at full supply, then measure that ratio: it is free
/// to land outside 2.5–3.0, and if it does, this layer is wrong.
///
/// This function exists so the analysis layer can report the DNA channel's
/// contribution separately from the combined result, which is the comparison
/// that makes the emergent number interpretable.
#[must_use = "the factor is the function's only output"]
pub fn dna_channel_dose_modifying_factor(cfg: &RadiationConfig, o2_supply: f64) -> f64 {
    let hypoxic = oer_scaled_dose(1.0, o2_supply, cfg.o2_dependence, cfg.p_full_mmhg);
    let full = oer_scaled_dose(1.0, 1.0, cfg.o2_dependence, cfg.p_full_mmhg);
    if hypoxic <= 0.0 {
        return f64::INFINITY;
    }
    full / hypoxic
}

// ── The schedule ─────────────────────────────────────────────────────────
//
// EVERYTHING ABOVE MODELS ONE FRACTION, AND NOBODY RECEIVES ONE FRACTION.
// A course of radiotherapy is a schedule, and the schedule is where the
// clinical reasoning lives: 60 Gy in 30 fractions and 55 Gy in 20 are
// prescribed as equivalent, and the arithmetic that makes them equivalent is
// the most-used calculation in the speciality. None of it was expressible
// here, which meant the arm modelled the physics of radiation and none of
// radiotherapy.
//
// The four Rs are the framework — repair, repopulation, redistribution,
// reoxygenation — and this section adds the first, the second and the fourth.
// Redistribution needs a cell cycle, which this engine does not have, and
// saying so is more useful than approximating it.

/// α/β for late-responding normal tissue, in Gy.
///
/// The other half of the ratio the tumour value ([`ALPHA_BETA_TUMOUR_GY`],
/// 10 Gy) belongs to. The whole therapeutic argument for fractionation is that
/// these two numbers differ: late-responding tissue is more sensitive to the
/// size of each fraction, so dividing a dose spares it more than it spares the
/// tumour. A model carrying one α/β cannot express the trade-off that
/// fractionation exists to exploit.
///
/// 3 Gy is the conventional value (Fowler 1989, PMID 2670032). It is a
/// convention rather than a measurement of any particular tissue, and the
/// spread across real late endpoints is wide.
pub const ALPHA_BETA_LATE_GY: f64 = 3.0;

/// Half-time of sublethal-damage repair, in hours.
///
/// Between two fractions a cell repairs the damage that would otherwise have
/// interacted to kill it. If the interval is short relative to this half-time
/// the repair is incomplete and the second fraction lands on a partly damaged
/// cell, which is why a six-hour gap is the standard minimum for
/// hyperfractionation.
///
/// ~1.5 h is the conventional value for human tissues; published estimates
/// span roughly 0.5 to 4 h and differ between early- and late-responding
/// tissue, so this is a placeholder in the sense the layer-freeze policy
/// means: the DIRECTION (shorter interval, less repair, more effect) is the
/// result, and a report resting on the precise number has to say so.
pub const SUBLETHAL_REPAIR_HALF_TIME_H: f64 = 1.5;

/// Days before accelerated clonogen repopulation begins, and the effective
/// clonogen doubling time once it has, from Withers 1988 (PMID 3390344).
///
/// The finding that changed how treatment gaps are regarded: after a lag of
/// about four weeks, surviving tumour clonogens proliferate fast enough that
/// each further day of treatment costs a substantial dose. It is measured for
/// head and neck cancer and is not a general constant.
pub const REPOP_KICKOFF_DAYS: f64 = 28.0;
/// See [`REPOP_KICKOFF_DAYS`].
pub const REPOP_DOUBLING_DAYS: f64 = 3.0;

/// α for head and neck squamous carcinoma, per Gy, the value Withers' own
/// repopulation analysis is built on.
///
/// Kept separate from [`ALPHA_GBM_PARAMETERISATION_PER_GY`] deliberately. The
/// two are different tumours and the prediction below is sensitive to which
/// one is used, which is a caveat the module states rather than hides.
pub const ALPHA_HEAD_NECK_PER_GY: f64 = 0.3;

/// The dose per day lost to repopulation, as Withers 1988 reports it, in Gy.
///
/// **This is the calibration target, and the model is free to miss it.**
/// [`repopulation_dose_per_day`] computes it from α and the doubling time —
/// two constants taken from the literature for other reasons — and the answer
/// is not fitted to this band.
pub const D_PROLIF_PUBLISHED_GY_PER_DAY: (f64, f64) = (0.5, 0.9);

/// A fractionation schedule: what is actually prescribed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Schedule {
    /// Number of fractions.
    pub n_fractions: u32,
    /// Dose per fraction, Gy.
    pub dose_per_fraction_gy: f64,
    /// Hours between consecutive fractions.
    pub interval_hours: f64,
    /// Elapsed days from the first fraction to the last, which is what
    /// repopulation consumes. NOT derivable from the other three: a schedule
    /// has weekends.
    pub overall_days: f64,
}

impl Schedule {
    /// A conventional once-daily, five-days-a-week course.
    #[must_use]
    pub fn conventional(n_fractions: u32, dose_per_fraction_gy: f64) -> Self {
        // Five fractions a week, so the elapsed span is the fraction count
        // stretched by 7/5, less the weekend that does not follow the last
        // fraction. An approximation, and the field it feeds says so.
        let weeks = f64::from(n_fractions) / 5.0;
        Self {
            n_fractions,
            dose_per_fraction_gy,
            interval_hours: 24.0,
            overall_days: (weeks * 7.0 - 2.0).max(0.0),
        }
    }

    /// Total prescribed dose, Gy.
    #[must_use]
    pub fn total_dose_gy(&self) -> f64 {
        f64::from(self.n_fractions) * self.dose_per_fraction_gy
    }
}

/// Biologically effective dose, `BED = nd(1 + d/(α/β))`.
///
/// The quantity that makes two schedules comparable at all. It is the total
/// dose scaled by a factor that grows with the size of each fraction, so a
/// hypofractionated course carries more biological effect per Gy — more so in
/// a tissue with a low α/β, which is the entire content of the fractionation
/// argument.
#[must_use]
pub fn bed(schedule: &Schedule, alpha_beta_gy: f64) -> f64 {
    debug_assert!(alpha_beta_gy > 0.0, "alpha/beta must be positive");
    schedule.total_dose_gy() * (1.0 + schedule.dose_per_fraction_gy / alpha_beta_gy)
}

/// The same dose expressed as an equivalent course given in 2 Gy fractions.
///
/// `EQD2 = BED / (1 + 2/(α/β))`. This is the form a clinic uses, because it
/// is in the units the conventional schedule is written in.
#[must_use]
pub fn eqd2(schedule: &Schedule, alpha_beta_gy: f64) -> f64 {
    bed(schedule, alpha_beta_gy) / (1.0 + 2.0 / alpha_beta_gy)
}

/// The Thames incomplete-repair factor `Hm` for `n` equally spaced fractions.
///
/// Standard BED assumes damage is fully repaired between fractions. When it
/// is not — fractions closer together than a few repair half-times — the
/// quadratic term is amplified by `1 + Hm`, and this returns `Hm`.
///
/// With `φ = exp(−μΔT)` and `μ = ln2 / T½`:
///
/// ```text
/// Hm = (2/n) · φ/(1−φ) · ( n − (1−φⁿ)/(1−φ) )
/// ```
///
/// (Thames 1985; the form used here follows Fowler 1989, PMID 2670032.)
///
/// **It reduces to zero, not approximately zero, as the interval grows.** At
/// the 24-hour interval of a conventional schedule and a 1.5-hour half-time,
/// `φ ≈ 1.5e−5` and the correction is below one part in ten thousand — which
/// is what makes this layer inert for every schedule the engine ran before it
/// existed.
#[must_use]
pub fn incomplete_repair_factor(
    n_fractions: u32,
    interval_hours: f64,
    repair_half_time_h: f64,
) -> f64 {
    if n_fractions <= 1 || repair_half_time_h <= 0.0 || !interval_hours.is_finite() {
        return 0.0;
    }
    let mu = std::f64::consts::LN_2 / repair_half_time_h;
    let phi = (-mu * interval_hours).exp();
    if phi <= 0.0 || phi >= 1.0 {
        return 0.0;
    }
    let n = f64::from(n_fractions);
    let geometric = (1.0 - phi.powf(n)) / (1.0 - phi);
    (2.0 / n) * (phi / (1.0 - phi)) * (n - geometric)
}

/// BED with the repair between fractions treated as incomplete.
///
/// `BED = nd · (1 + (1 + Hm)·d/(α/β))`. Equal to [`bed`] whenever `Hm` is 0,
/// which is every schedule with a day between fractions.
#[must_use]
pub fn bed_with_incomplete_repair(
    schedule: &Schedule,
    alpha_beta_gy: f64,
    repair_half_time_h: f64,
) -> f64 {
    let hm = incomplete_repair_factor(
        schedule.n_fractions,
        schedule.interval_hours,
        repair_half_time_h,
    );
    schedule.total_dose_gy() * (1.0 + (1.0 + hm) * schedule.dose_per_fraction_gy / alpha_beta_gy)
}

/// The dose each further treatment day costs once repopulation has begun,
/// `D_prolif = ln2 / (α · T_p)`, in Gy per day.
///
/// **The layer's failable prediction.** α comes from a survival
/// parameterisation and `T_p` from a repopulation study; neither was chosen to
/// make this number land anywhere. Whether it falls inside
/// [`D_PROLIF_PUBLISHED_GY_PER_DAY`] is therefore a test the model can lose,
/// and it is sensitive to which tumour's α is used — the glioblastoma value
/// this module also carries puts it outside the band.
#[must_use]
pub fn repopulation_dose_per_day(alpha_per_gy: f64, doubling_days: f64) -> f64 {
    if alpha_per_gy <= 0.0 || doubling_days <= 0.0 {
        return 0.0;
    }
    std::f64::consts::LN_2 / (alpha_per_gy * doubling_days)
}

/// BED reduced by the dose repopulation consumes over the treatment course.
///
/// Nothing is subtracted before `kickoff_days`, which is the shape of Withers'
/// finding: the loss is not proportional to treatment time from the start, it
/// begins after a lag. A course that finishes inside the lag pays nothing,
/// which is why the layer is inert for short schedules rather than merely
/// small.
#[must_use]
pub fn bed_with_repopulation(
    bed_gy: f64,
    overall_days: f64,
    kickoff_days: f64,
    alpha_per_gy: f64,
    doubling_days: f64,
) -> f64 {
    let over = overall_days - kickoff_days;
    if over <= 0.0 {
        return bed_gy;
    }
    (bed_gy - over * repopulation_dose_per_day(alpha_per_gy, doubling_days)).max(0.0)
}

/// Surviving fraction after a whole schedule, `exp(−n(αd + (1+Hm)βd²))`.
///
/// The fractionated form of [`lq_survival`]. Reduces to it exactly at one
/// fraction.
#[must_use]
pub fn schedule_survival(
    schedule: &Schedule,
    alpha_per_gy: f64,
    beta_per_gy2: f64,
    repair_half_time_h: f64,
) -> f64 {
    let hm = incomplete_repair_factor(
        schedule.n_fractions,
        schedule.interval_hours,
        repair_half_time_h,
    );
    let d = schedule.dose_per_fraction_gy;
    let per_fraction = alpha_per_gy * d + (1.0 + hm) * beta_per_gy2 * d * d;
    (-f64::from(schedule.n_fractions) * per_fraction).exp()
}

/// The α/β at which two schedules deliver the same EQD2, in Gy.
///
/// Solving `D₁(1 + d₁/x) = D₂(1 + d₂/x)` gives
/// `x = (D₂d₂ − D₁d₁) / (D₁ − D₂)`.
///
/// **This is the layer's second external check, and it runs backwards.**
/// Randomised fractionation trials compare two schedules and report whether
/// the outcomes differ; if they do not, the LQ model implies the α/β that
/// makes them isoeffective. That implied value can then be compared against
/// α/β estimates the radiobiology literature derives independently — so two
/// numbers from a trial protocol predict a third quantity nobody put in.
///
/// Returns `None` when the two schedules deliver the same total dose, where no
/// α/β is implied and the expression is a division by zero.
#[must_use]
pub fn isoeffect_alpha_beta(a: &Schedule, b: &Schedule) -> Option<f64> {
    let (d1, d2) = (a.total_dose_gy(), b.total_dose_gy());
    if (d1 - d2).abs() < 1e-9 {
        return None;
    }
    let x = (d2 * b.dose_per_fraction_gy - d1 * a.dose_per_fraction_gy) / (d1 - d2);
    Some(x)
}

/// Tumour effect per unit of late-normal-tissue effect, for a schedule.
///
/// `EQD2(tumour α/β) / EQD2(late α/β)`, and the quantitative form of the
/// reason radiotherapy is fractionated at all.
///
/// A single large dose kills more tumour than the same physical dose split up
/// — the quadratic term grows with the size of each fraction — and a model
/// carrying only the tumour would conclude that one fraction is better. It is
/// not, and the missing half is here: late-responding normal tissue has a
/// LOWER α/β, so it is hurt *more* by the same increase in fraction size. The
/// ratio falls as fractions grow.
///
/// Returns 1.0 for a schedule given in 2 Gy fractions, by construction: EQD2
/// is defined against that schedule, so the number is a comparison to
/// convention rather than an absolute.
#[must_use]
pub fn therapeutic_ratio(schedule: &Schedule, alpha_beta_tumour: f64, alpha_beta_late: f64) -> f64 {
    let late = eqd2(schedule, alpha_beta_late);
    if late <= 0.0 {
        return f64::INFINITY;
    }
    eqd2(schedule, alpha_beta_tumour) / late
}

/// The hypoxic fraction remaining after `k` fractions, given a reoxygenation
/// half-life measured IN FRACTIONS.
///
/// The fourth R, and the one that interacts with what this engine already has:
/// [`crate::oxygen::oxygen_enhancement_ratio`] makes a hypoxic cell about
/// three times harder to kill, and a single large dose leaves that population
/// behind. Between fractions the tumour reoxygenates — vasculature reaches
/// cells that were beyond it once the cells in front have died — so the
/// hypoxic compartment is refilled from a shrinking pool rather than being a
/// fixed shield.
///
/// `f(k) = f₀ · 2^(−k/h)`, the simplest form that expresses the direction.
/// **Uncalibrated.** The half-life is a free parameter here; published
/// reoxygenation kinetics vary by orders of magnitude across models, and no
/// dataset in this repository constrains it. The DIRECTION — that
/// fractionation converts a persistent hypoxic shield into a decaying one — is
/// the result.
#[must_use]
pub fn hypoxic_fraction_after(
    initial_hypoxic_fraction: f64,
    fractions_delivered: u32,
    reoxygenation_half_life_fractions: f64,
) -> f64 {
    if reoxygenation_half_life_fractions <= 0.0 {
        return initial_hypoxic_fraction;
    }
    let k = f64::from(fractions_delivered) / reoxygenation_half_life_fractions;
    (initial_hypoxic_fraction * (-k * std::f64::consts::LN_2).exp()).clamp(0.0, 1.0)
}

/// Surviving fraction over a schedule with a reoxygenating hypoxic
/// compartment.
///
/// Two populations. The oxic one takes the full dose; the hypoxic one takes
/// the dose scaled by [`oer_scaled_dose`], the engine's own oxygen-effect
/// function, so this cannot drift away from the OER the rest of the module
/// uses. Between fractions a share of the hypoxic survivors becomes oxic:
/// cells that were beyond the reach of a vessel come within it once the cells
/// in front of them have died.
///
/// **The comparison is the point.** With `reoxygenation_half_life_fractions`
/// at infinity nothing moves, the hypoxic survivors accumulate — they die
/// more slowly, so their SHARE of the surviving population climbs fraction by
/// fraction — and the course stalls against a shield it cannot remove. Turn
/// reoxygenation on and the same total dose in the same number of fractions
/// reaches further. That difference is the fourth R, and it is the reason a
/// single large dose is not equivalent to the fractionated course that
/// delivers the same physical dose.
///
/// **Uncalibrated.** The half-life is a free parameter; see
/// [`hypoxic_fraction_after`]. The direction is the result.
#[must_use]
pub fn reoxygenating_schedule_survival(
    schedule: &Schedule,
    alpha_per_gy: f64,
    beta_per_gy2: f64,
    initial_hypoxic_fraction: f64,
    reoxygenation_half_life_fractions: f64,
    hypoxic_o2_supply: f64,
    p_full_mmhg: f64,
) -> f64 {
    let d = schedule.dose_per_fraction_gy;
    // The engine's own oxygen effect, at full dependence: a hypoxic cell
    // behaves as though it received a smaller dose. An earlier draft of this
    // function divided two calls to `oxygen_enhancement_ratio` by hand and got
    // a factor that moved the dose the WRONG WAY -- reusing the module's
    // existing helper is not tidiness, it is the thing that makes a second
    // oxygen model impossible.
    let hypoxic_dose = oer_scaled_dose(d, hypoxic_o2_supply, 1.0, p_full_mmhg);

    // Per-fraction transfer out of the hypoxic compartment. At an infinite
    // half-life this is exactly zero, so the two arms differ in this term and
    // in nothing else.
    let moved_share = if reoxygenation_half_life_fractions.is_finite()
        && reoxygenation_half_life_fractions > 0.0
    {
        1.0 - (-std::f64::consts::LN_2 / reoxygenation_half_life_fractions).exp()
    } else {
        0.0
    };

    let f0 = initial_hypoxic_fraction.clamp(0.0, 1.0);
    let mut oxic = 1.0 - f0;
    let mut hypoxic = f0;
    let s_ox = (-(alpha_per_gy * d + beta_per_gy2 * d * d)).exp();
    let s_hy = (-(alpha_per_gy * hypoxic_dose + beta_per_gy2 * hypoxic_dose * hypoxic_dose)).exp();
    for _ in 0..schedule.n_fractions {
        oxic *= s_ox;
        hypoxic *= s_hy;
        let moved = hypoxic * moved_share;
        hypoxic -= moved;
        oxic += moved;
    }
    (oxic + hypoxic).clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::oxygen::oxygen_enhancement_ratio;
    use crate::params::SpatialParams;
    use crate::physics::{pdt_intensity_at_depth, sdt_intensity_at_depth};

    fn cfg() -> RadiationConfig {
        RadiationConfig {
            dose_gy: 2.0,
            alpha_per_gy: 0.3,
            beta_per_gy2: 0.03,
            ros_per_gy: 0.0,
            o2_dependence: 1.0,
            p_full_mmhg: crate::oxygen::OER_REFERENCE_PO2_MMHG,
            mu_per_cm: MU_6MV_SOFT_TISSUE_PER_CM,
            parp_alpha_boost: 0.0,
            hr_deficiency: 0.0,
        }
    }

    #[test]
    fn the_default_config_is_the_identity_element() {
        // The byte-identity contract: a run carrying the default must be
        // indistinguishable from a run with no radiation model. Asserted on
        // BIT PATTERNS, not with a tolerance -- `1.0 - 1e-17` would pass an
        // approximate check and still move a committed hash.
        let d = RadiationConfig::default();
        assert_eq!(d.dose_gy.to_bits(), 0.0_f64.to_bits());
        assert_eq!(
            lq_survival(d.dose_gy, d.alpha_per_gy, d.beta_per_gy2).to_bits(),
            1.0_f64.to_bits()
        );
        assert_eq!(dna_lethality(&d, 1.0, 0.0).to_bits(), 0.0_f64.to_bits());
        assert_eq!(exo_ros_from_dose(&d, 1.0, 0.0).to_bits(), 0.0_f64.to_bits());
        // And at every oxygenation and depth, not just the convenient one.
        for &o2 in &[0.0, 0.05, 0.5, 1.0] {
            for &z in &[0.0, 500.0, 10_000.0] {
                assert_eq!(dna_lethality(&d, o2, z).to_bits(), 0.0_f64.to_bits());
                assert_eq!(exo_ros_from_dose(&d, o2, z).to_bits(), 0.0_f64.to_bits());
            }
        }
    }

    #[test]
    fn the_lq_implementation_round_trips_the_published_parameterisation() {
        // A ROUND TRIP, not a prediction, and the earlier name claimed
        // otherwise. PMID 32307022 states three numbers, but it DERIVED the
        // third from the first two through the LQ model, and under LQ
        // `(SF2, alpha/beta)` determines `SF(2.67)` anyway. So this shows the
        // implementation matches the published parameterisation; it does not
        // test the biology.
        let beta = ALPHA_GBM_PARAMETERISATION_PER_GY / ALPHA_BETA_TUMOUR_GY;

        // The fit reproduces the point it was fitted to.
        let sf2 = lq_survival(2.0, ALPHA_GBM_PARAMETERISATION_PER_GY, beta);
        assert!(
            (sf2 - 0.5445).abs() < 1e-4,
            "the fitted alpha does not reproduce SF2 = 0.5445: {sf2}"
        );

        // And PREDICTS the one it was not.
        let sf267 = lq_survival(2.67, ALPHA_GBM_PARAMETERISATION_PER_GY, beta);
        assert!(
            (sf267 - SF_AT_2_67_GY_PUBLISHED).abs() < 1e-3,
            "predicted SF(2.67) = {sf267}, published {SF_AT_2_67_GY_PUBLISHED}"
        );

        // The control that MATTERS is the ratio, not alpha. Scaling both
        // alpha and beta by 1.2 also breaks SF2, so it only shows that SF
        // responds to alpha -- which nobody doubted. What this check has to
        // discriminate is the alpha/beta CONVENTION, since a reader could
        // reasonably reach for the late-responding 3 Gy value. So: re-derive
        // alpha at a different ratio, keeping SF2 EXACTLY right, and require
        // the 2.67 Gy point to move off the published number.
        for ratio in [3.0_f64, 5.0, 20.0] {
            let a = -0.5445_f64.ln() / (2.0 + 4.0 / ratio);
            let b = a / ratio;
            assert!(
                (lq_survival(2.0, a, b) - 0.5445).abs() < 1e-9,
                "the control must reproduce SF2 exactly, or it is testing SF2"
            );
            assert!(
                (lq_survival(2.67, a, b) - SF_AT_2_67_GY_PUBLISHED).abs() > 1e-3,
                "alpha/beta = {ratio} still hits the published SF(2.67), so \
                 this check cannot tell the conventions apart"
            );
        }
    }

    #[test]
    fn survival_falls_monotonically_with_dose() {
        let c = cfg();
        let mut prev = 1.0;
        for step in 0..40 {
            let d = step as f64 * 0.5;
            let sf = lq_survival(d, c.alpha_per_gy, c.beta_per_gy2);
            assert!(sf <= prev, "survival rose at {d} Gy: {sf} > {prev}");
            prev = sf;
        }
        assert!(prev < 1e-6, "20 Gy should be near-sterilising, got {prev}");
    }

    #[test]
    fn the_quadratic_term_makes_one_dose_deadlier_than_two_halves() {
        // The signature of the LQ model and the whole reason fractionation
        // works. Pinned in BOTH directions so the test cannot pass with beta
        // ignored: with beta = 0 the model is exponential and splitting the
        // dose changes nothing AT ALL (bit-equal), and with beta > 0 the
        // single dose must strictly win.
        let (a, b) = (0.3, 0.03);
        let d = 6.0;
        let one = lq_survival(d, a, b);
        let two = lq_survival(d / 2.0, a, b).powi(2);
        assert!(
            one < two,
            "beta>0 must make one dose deadlier: {one} vs {two}"
        );

        let one_lin = lq_survival(d, a, 0.0);
        let two_lin = lq_survival(d / 2.0, a, 0.0).powi(2);
        assert!(
            (one_lin - two_lin).abs() < 1e-12,
            "beta=0 must be exponential and fractionation-free: {one_lin} vs {two_lin}"
        );
    }

    #[test]
    fn the_oxygen_scaling_is_off_at_zero_dependence_and_is_the_hyperbola_at_one() {
        // Off must be BIT-identical, because this is the flag that keeps the
        // production matrix stable.
        for &o2 in &[0.0, 0.1, 0.5, 1.0] {
            assert_eq!(
                oer_scaled_dose(2.0, o2, 0.0, 40.0).to_bits(),
                2.0_f64.to_bits(),
                "dependence 0 moved the dose at o2={o2}"
            );
        }
        // On, it must be the SAME hyperbola the SDT path uses -- recomputed
        // here from `oxygen_enhancement_ratio` rather than restated, so a
        // change to that function fails this test instead of silently letting
        // the two modalities drift apart in their O2 shape.
        let m0 = oxygen_enhancement_ratio(0.0);
        let m40 = oxygen_enhancement_ratio(40.0);
        assert!(
            (oer_scaled_dose(2.0, 0.0, 1.0, 40.0) - 2.0 * m0 / m40).abs() < 1e-12,
            "the anoxic scaling is not m(0)/m(40)"
        );
    }

    #[test]
    fn the_dna_channel_factor_is_the_hyperbola_and_says_so() {
        // Asserted for what it IS -- a restatement of the OER -- so nobody
        // later reads it as the layer's prediction. The doc comment on
        // `dna_channel_dose_modifying_factor` makes the same point in prose;
        // this pins it.
        let c = cfg();
        let m0 = oxygen_enhancement_ratio(0.0);
        let m40 = oxygen_enhancement_ratio(c.p_full_mmhg);
        let f = dna_channel_dose_modifying_factor(&c, 0.0);
        assert!(
            (f - m40 / m0).abs() < 1e-12,
            "the factor should be exactly m(p_full)/m(0), got {f}"
        );
        // At full supply it is exactly 1: no modification.
        assert!((dna_channel_dose_modifying_factor(&c, 1.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn radiation_reaches_depth_where_light_and_ultrasound_do_not() {
        // The claim the module docs make about the manuscript's tissue-access
        // argument, measured against the two modalities it is compared with
        // rather than asserted. 1 cm is the scale that argument uses.
        let p = SpatialParams::default();
        let z = 10_000.0_f64; // 1 cm
        let rad = intensity_at_depth(z, MU_6MV_SOFT_TISSUE_PER_CM);
        let pdt = pdt_intensity_at_depth(z, &p);
        let sdt = sdt_intensity_at_depth(z, &p);
        assert!(
            rad > 0.95,
            "6 MV should be nearly flat over 1 cm, got {rad}"
        );
        assert!(
            rad > sdt,
            "radiation must outreach SDT at 1 cm: {rad} vs {sdt}"
        );
        assert!(
            rad > pdt,
            "radiation must outreach PDT at 1 cm: {rad} vs {pdt}"
        );
        // And the gap must be a GAP, not a rounding difference. Stated as
        // absolute bounds rather than a ratio: `pdt` can legitimately be zero
        // (the default drug-light interval yields no sensitizer), and a
        // division would then report `inf` as a pass.
        assert!(
            pdt < 0.1,
            "PDT should be nearly extinguished at 1 cm, got {pdt}"
        );
        assert!(sdt < 0.9, "SDT should be attenuated at 1 cm, got {sdt}");
    }

    #[test]
    fn attenuation_is_monotone_bounded_and_clamps_outside_the_tumour() {
        // STRICTLY decreasing, and pinned to a computed value. Every bound in
        // the depth tests was one-sided -- `rad > 0.95`, `rad > sdt`,
        // `rad > pdt` all bound mu from ABOVE, the direction the "radiation is
        // flat" claim wants -- so replacing this whole function with `1.0`
        // passed all thirteen tests.
        assert!(
            (intensity_at_depth(100_000.0, MU_6MV_SOFT_TISSUE_PER_CM) - 0.740_818).abs() < 1e-5,
            "10 cm should receive exp(-0.3) = 0.7408 of the surface dose, got {}",
            intensity_at_depth(100_000.0, MU_6MV_SOFT_TISSUE_PER_CM)
        );
        let mut prev = 1.0;
        for cm in 1..30 {
            let v = intensity_at_depth(cm as f64 * 10_000.0, MU_6MV_SOFT_TISSUE_PER_CM);
            assert!(
                v < prev,
                "attenuation is not strictly decreasing at {cm} cm: {v}"
            );
            assert!(v >= 0.0 && v <= 1.0, "at {cm} cm: {v}");
            prev = v;
        }
        // Negative depth is outside the tumour; the grid's sign convention
        // makes that reachable, and it must read as the surface rather than
        // amplifying the dose.
        assert_eq!(
            intensity_at_depth(-5_000.0, MU_6MV_SOFT_TISSUE_PER_CM).to_bits(),
            1.0_f64.to_bits()
        );
        // mu = 0 is the no-attenuation limit, exactly.
        assert_eq!(
            intensity_at_depth(10_000.0, 0.0).to_bits(),
            1.0_f64.to_bits()
        );
    }

    #[test]
    fn the_two_channels_are_independent() {
        // The design claim: DNA lethality does not depend on `ros_per_gy`, and
        // the ferroptosis channel does not depend on alpha or beta. A single
        // implementation that folded them together would fail this.
        let base = cfg();
        let mut ros_heavy = base;
        ros_heavy.ros_per_gy = 5.0;
        assert_eq!(
            dna_lethality(&base, 0.5, 0.0).to_bits(),
            dna_lethality(&ros_heavy, 0.5, 0.0).to_bits(),
            "ros_per_gy moved the DNA channel"
        );

        let mut lq_heavy = base;
        lq_heavy.alpha_per_gy = 0.9;
        lq_heavy.beta_per_gy2 = 0.5;
        let mut a = base;
        a.ros_per_gy = 2.0;
        let mut b = lq_heavy;
        b.ros_per_gy = 2.0;
        assert_eq!(
            exo_ros_from_dose(&a, 0.5, 0.0).to_bits(),
            exo_ros_from_dose(&b, 0.5, 0.0).to_bits(),
            "alpha/beta moved the ferroptosis channel"
        );
    }

    #[test]
    fn hypoxia_spares_a_cell_in_both_channels() {
        // Both channels must respond to O2 in the same DIRECTION, since both
        // read the same hyperbola. A sign error in either would show here.
        let mut c = cfg();
        c.ros_per_gy = 1.0;
        let (well, poor) = (1.0, 0.02);
        assert!(dna_lethality(&c, well, 0.0) > dna_lethality(&c, poor, 0.0));
        assert!(exo_ros_from_dose(&c, well, 0.0) > exo_ros_from_dose(&c, poor, 0.0));
    }

    #[test]
    fn radiation_through_from_cell_has_no_exogenous_ros() {
        // The CONTRACT `biochem.rs` documents, asserted rather than left as a
        // silent trap: `CellState::from_cell` cannot see dose, depth or local
        // O2, so a Radiation cell built through it carries zero exogenous ROS
        // and the ferroptosis channel is off. A caller wanting that channel
        // goes through `from_cell_with_ros` fed by `exo_ros_from_dose`.
        use crate::biochem::CellState;
        use crate::cell::{gen_cell, Phenotype, Treatment};
        use crate::params::Params;
        use rand::{rngs::StdRng, SeedableRng};

        let params = Params::default();
        let mut rng = StdRng::seed_from_u64(7);
        let cell = gen_cell(Phenotype::Glycolytic, &mut rng);
        let state = CellState::from_cell(&cell, Treatment::Radiation, &params, &mut rng);
        assert_eq!(state.exo_ros_peak.to_bits(), 0.0_f64.to_bits());

        // And the explicit path DOES deliver it, so the contract is a routing
        // rule rather than the channel being broken.
        let mut cfg = cfg();
        cfg.ros_per_gy = 2.0;
        let ros = exo_ros_from_dose(&cfg, 1.0, 0.0);
        assert!(ros > 0.0, "the dose-driven path yields nothing: {ros}");
        let state = CellState::from_cell_with_ros(&cell, Treatment::Radiation, &params, ros);
        assert_eq!(state.exo_ros_peak, ros);
    }

    #[test]
    fn the_depth_arm_is_reachable_through_the_engine_not_only_this_module() {
        // `physics::local_ros_multiplier` matches on the treatment, and the
        // depth claim is about THAT function rather than about
        // `intensity_at_depth` in isolation. Measured through the engine's
        // own entry point, at the 1 cm scale the tissue-access argument uses.
        use crate::cell::Treatment;
        use crate::params::SpatialParams;
        use crate::physics::local_ros_multiplier_3d;

        let p = SpatialParams::default();
        let z = 10_000.0_f64;
        let rad = local_ros_multiplier_3d(z, Treatment::Radiation, &p);
        let sdt = local_ros_multiplier_3d(z, Treatment::SDT, &p);
        let pdt = local_ros_multiplier_3d(z, Treatment::PDT, &p);
        assert!(rad > 0.95, "radiation should be near-flat at 1 cm: {rad}");
        assert!(rad > sdt && rad > pdt, "{rad} vs sdt {sdt} pdt {pdt}");
    }

    #[test]
    fn the_published_ser_band_constrains_the_one_free_parameter() {
        // WHY THIS IS A REAL CALIBRATION and the alpha/beta round trip is not.
        //
        // `parp_alpha_boost` is one number. The SER it produces is
        // DOSE-DEPENDENT under LQ, because raising alpha shifts the linear
        // term while beta is unchanged -- so a single boost has to satisfy
        // the published 1.2-1.7 band across the whole clinically relevant
        // survival range at once, and it is free to fail.
        //
        // It does not fail, and the band it survives in is NARROW: scanning
        // boosts from 0 to 2, only [0.544, 0.948] keeps the SER inside
        // 1.2-1.7 for every SF from 0.01 to 0.5. That is 20% of the scanned
        // range, which is what "the data constrains the parameter" means.
        let (lo_band, hi_band) = PARP_SER_BAND;
        let base = RadiationConfig {
            alpha_per_gy: ALPHA_GBM_PARAMETERISATION_PER_GY,
            beta_per_gy2: ALPHA_GBM_PARAMETERISATION_PER_GY / ALPHA_BETA_TUMOUR_GY,
            ..RadiationConfig::default()
        };
        let ser_at = |boost: f64, sf: f64| {
            sensitizer_enhancement_ratio(
                &RadiationConfig {
                    parp_alpha_boost: boost,
                    ..base
                },
                sf,
            )
        };

        let (sf_lo, sf_hi) = (0.01_f64, 0.5_f64);
        let mut admissible: Vec<f64> = Vec::new();
        let mut b = 0.0_f64;
        while b <= 2.0 {
            if ser_at(b, sf_lo) >= lo_band && ser_at(b, sf_hi) <= hi_band {
                admissible.push(b);
            }
            b += 0.001;
        }
        assert!(
            !admissible.is_empty(),
            "NO boost satisfies the published SER band across the survival \
             range -- the sensitization FORM is wrong, not the parameter"
        );
        let (first, last) = (admissible[0], admissible[admissible.len() - 1]);
        assert!(
            (first - 0.544).abs() < 0.01 && (last - 0.948).abs() < 0.01,
            "the admissible boost window moved to [{first:.3}, {last:.3}] \
             from [0.544, 0.948]; the published band now implies a different \
             parameter and the CALIBRATION_STATUS row must say so"
        );
        // NON-VACUOUS: the band must actually EXCLUDE most of the range, or
        // "constrained" means nothing.
        assert!(
            last - first < 0.6,
            "the admissible window spans {:.3}, which is most of the scanned \
             range -- the published band is not constraining anything",
            last - first
        );

        // And the DOSE-DEPENDENCE is what makes it a joint constraint rather
        // than a single equation. At beta = 0 the SER is exactly 1 + boost at
        // every dose, so any boost in the band would pass at every SF and the
        // test above would be measuring one number, not a range.
        let flat = RadiationConfig {
            beta_per_gy2: 0.0,
            parp_alpha_boost: 0.4,
            ..base
        };
        for &sf in &[0.5_f64, 0.1, 0.01] {
            assert!(
                (sensitizer_enhancement_ratio(&flat, sf) - 1.4).abs() < 1e-12,
                "with beta = 0 the SER must be exactly 1 + boost at every dose"
            );
        }
        // With beta > 0 it must MOVE, and downward as dose rises.
        let curved = RadiationConfig {
            parp_alpha_boost: 0.4,
            ..base
        };
        assert!(
            sensitizer_enhancement_ratio(&curved, 0.5)
                > sensitizer_enhancement_ratio(&curved, 0.01),
            "the SER is not dose-dependent, so one dose would fix the fit"
        );
    }

    #[test]
    fn no_sensitization_reports_no_enhancement_and_alpha_is_untouched() {
        let c = cfg();
        assert_eq!(c.parp_alpha_boost, 0.0);
        assert_eq!(sensitized_alpha(&c).to_bits(), c.alpha_per_gy.to_bits());
        for &sf in &[0.9_f64, 0.5, 0.1, 0.001] {
            assert_eq!(
                sensitizer_enhancement_ratio(&c, sf).to_bits(),
                1.0_f64.to_bits(),
                "an unsensitized run must report exactly 1.0, not a near-one \
                 artefact that a reader would take for a small real effect"
            );
        }
        // Beta is NOT scaled: the added lethality is one-track, not pairwise,
        // and a model scaling both would assert a mechanism nobody measured.
        let sens = RadiationConfig {
            parp_alpha_boost: 0.5,
            ..c
        };
        assert_eq!(sens.beta_per_gy2.to_bits(), c.beta_per_gy2.to_bits());
        assert!((sensitized_alpha(&sens) - c.alpha_per_gy * 1.5).abs() < 1e-12);
    }

    #[test]
    fn synthetic_lethality_needs_both_hits_and_says_so() {
        // The defining property: neither hit alone is lethal. That is why
        // these drugs ship with a biomarker, and a model returning a small
        // nonzero kill for the HR-proficient case would make the biomarker
        // look optional.
        let c = cfg();
        for (hrd, boost) in [(0.0, 0.0), (0.0, 0.8), (1.0, 0.0)] {
            let x = RadiationConfig {
                hr_deficiency: hrd,
                parp_alpha_boost: boost,
                ..c
            };
            assert_eq!(
                parp_monotherapy_lethality(&x).to_bits(),
                0.0_f64.to_bits(),
                "one hit alone (hrd={hrd}, boost={boost}) produced kill; synthetic lethality requires BOTH"
            );
        }
        // Both present: nonzero, monotone in each, bounded.
        let both = RadiationConfig {
            hr_deficiency: 0.8,
            parp_alpha_boost: 0.7,
            ..c
        };
        let l = parp_monotherapy_lethality(&both);
        assert!(l > 0.0 && l <= 1.0, "{l}");
        let more_hrd = RadiationConfig {
            hr_deficiency: 1.0,
            ..both
        };
        assert!(parp_monotherapy_lethality(&more_hrd) > l);
        let more_drug = RadiationConfig {
            parp_alpha_boost: 1.0,
            ..both
        };
        assert!(parp_monotherapy_lethality(&more_drug) > l);
        // Out-of-range inputs clamp rather than escaping [0, 1].
        let wild = RadiationConfig {
            hr_deficiency: 5.0,
            parp_alpha_boost: 9.0,
            ..c
        };
        assert!((parp_monotherapy_lethality(&wild) - 1.0).abs() < 1e-12);
    }

    #[test]
    #[should_panic(expected = "dose_gy must be finite and non-negative")]
    fn a_negative_dose_is_rejected() {
        lq_survival(-1.0, 0.3, 0.03);
    }
    // ── The schedule ─────────────────────────────────────────────────────

    #[test]
    fn eqd2_of_a_two_gray_schedule_is_its_own_total_dose() {
        // EQD2 is defined against 2 Gy fractions, so this is an identity and a
        // good place for an arithmetic slip to show: it holds for EVERY
        // alpha/beta, so a formula with the ratio on the wrong side of a
        // division breaks it at some ratio even if it survives one.
        for ab in [1.5, 3.0, 10.0, 20.0] {
            let s = Schedule::conventional(30, 2.0);
            assert!((eqd2(&s, ab) - 60.0).abs() < 1e-9, "alpha/beta {ab}");
        }
    }

    #[test]
    fn bed_grows_with_fraction_size_at_fixed_total_dose() {
        let total = 60.0;
        let mut previous = 0.0;
        for (n, d) in [(30u32, 2.0f64), (20, 3.0), (10, 6.0), (3, 20.0)] {
            let s = Schedule::conventional(n, d);
            assert!((s.total_dose_gy() - total).abs() < 1e-9);
            let b = bed(&s, ALPHA_BETA_TUMOUR_GY);
            assert!(b > previous, "{n}x{d} did not raise BED");
            previous = b;
        }
    }

    #[test]
    fn the_repair_correction_is_inert_at_a_day_between_fractions() {
        // THE INERTNESS CONTRACT for this layer. Every schedule the engine
        // could express before this section existed is once-daily, so the
        // incomplete-repair term has to be indistinguishable from zero there
        // or the layer changes results it was not supposed to touch.
        let hm = incomplete_repair_factor(30, 24.0, SUBLETHAL_REPAIR_HALF_TIME_H);
        assert!(hm < 1e-4, "Hm at 24 h was {hm}, not negligible");
        let s = Schedule::conventional(30, 2.0);
        let plain = bed(&s, ALPHA_BETA_TUMOUR_GY);
        let with_repair =
            bed_with_incomplete_repair(&s, ALPHA_BETA_TUMOUR_GY, SUBLETHAL_REPAIR_HALF_TIME_H);
        assert!((plain - with_repair).abs() < 1e-3);
    }

    #[test]
    fn bed_with_incomplete_repair_actually_uses_the_correction() {
        // A MUTATION SURVIVOR, and the reason to run the sweep. The inertness
        // test above passes whether or not `bed_with_incomplete_repair` reads
        // Hm at all, because at 24 hours the correction is zero either way --
        // so deleting the term entirely left every test green. This is the
        // case where the term has to bite: twice-daily treatment.
        let twice_daily = Schedule {
            n_fractions: 60,
            dose_per_fraction_gy: 1.2,
            interval_hours: 6.0,
            overall_days: 42.0,
        };
        let plain = bed(&twice_daily, ALPHA_BETA_LATE_GY);
        let corrected = bed_with_incomplete_repair(
            &twice_daily,
            ALPHA_BETA_LATE_GY,
            SUBLETHAL_REPAIR_HALF_TIME_H,
        );
        assert!(
            corrected > plain,
            "unrepaired damage between six-hourly fractions did not raise BED \
             ({corrected} vs {plain})"
        );
        // And it lands on the LATE-responding tissue hardest, which is the
        // clinical point of the six-hour rule: the correction scales the
        // quadratic term, and late tissue has the lower alpha/beta.
        let late_penalty = corrected / plain;
        let tumour_penalty = bed_with_incomplete_repair(
            &twice_daily,
            ALPHA_BETA_TUMOUR_GY,
            SUBLETHAL_REPAIR_HALF_TIME_H,
        ) / bed(&twice_daily, ALPHA_BETA_TUMOUR_GY);
        assert!(
            late_penalty > tumour_penalty,
            "{late_penalty} vs {tumour_penalty}"
        );
    }

    #[test]
    fn a_conventional_course_spans_the_calendar_it_should() {
        // ALSO A MUTATION SURVIVOR: nothing pinned the fractions-to-days
        // mapping, so a schedule could claim any elapsed time -- and elapsed
        // time is what repopulation charges for, so this is not cosmetic.
        // Five fractions a week: thirty fractions is six weeks.
        let six_weeks = Schedule::conventional(30, 2.0);
        assert!(
            (six_weeks.overall_days - 40.0).abs() < 0.5,
            "thirty daily fractions spanned {} days",
            six_weeks.overall_days
        );
        let one_week = Schedule::conventional(5, 2.0);
        assert!(
            (one_week.overall_days - 5.0).abs() < 0.5,
            "five fractions spanned {} days",
            one_week.overall_days
        );
        // The consequence, asserted here so the mapping cannot drift without
        // something failing: a six-week course is past the repopulation lag
        // and a one-week course is not.
        assert!(six_weeks.overall_days > REPOP_KICKOFF_DAYS);
        assert!(one_week.overall_days < REPOP_KICKOFF_DAYS);
    }

    #[test]
    fn shorter_gaps_leave_more_damage_unrepaired() {
        // And the ordering is the claim: six hours is the standard minimum
        // interval for treating twice a day, and the model has to say why.
        let six = incomplete_repair_factor(60, 6.0, SUBLETHAL_REPAIR_HALF_TIME_H);
        let two = incomplete_repair_factor(60, 2.0, SUBLETHAL_REPAIR_HALF_TIME_H);
        let day = incomplete_repair_factor(60, 24.0, SUBLETHAL_REPAIR_HALF_TIME_H);
        assert!(day < six && six < two, "{day} {six} {two}");
        assert!(six > 0.05, "a six-hour gap should be measurably incomplete");
    }

    #[test]
    fn the_prostate_isoeffect_recovers_a_low_alpha_beta() {
        // THE LAYER'S EXTERNAL CHECK, and it runs backwards. CHHiP (PMID
        // 27339115) randomised 74 Gy in 37 fractions against 60 Gy in 20 and
        // found the shorter schedule non-inferior. Two numbers from a trial
        // PROTOCOL then imply the alpha/beta that makes them isoeffective --
        // and the value that falls out is the low one the prostate
        // radiobiology literature independently estimates, which is the
        // observation that motivated hypofractionating prostate cancer at all.
        //
        // Nothing here is fitted: the schedules are the trial's, the formula
        // is the LQ model, and the band is somebody else's estimate.
        let conventional = Schedule::conventional(37, 2.0);
        let hypofractionated = Schedule::conventional(20, 3.0);
        let ab = isoeffect_alpha_beta(&conventional, &hypofractionated)
            .expect("two schedules of different total dose imply a ratio");
        assert!(
            (1.2..=3.0).contains(&ab),
            "the implied alpha/beta was {ab:.2} Gy, outside the published \
             prostate band -- either the arithmetic moved or the schedules did"
        );
        // A prostate-like ratio is the FINDING, so the control matters: the
        // same computation on a generic tumour ratio must NOT land there.
        assert!(ab < ALPHA_BETA_TUMOUR_GY / 2.0);
    }

    #[test]
    fn two_schedules_of_the_same_total_dose_imply_nothing() {
        // The degenerate case, which is a division by zero rather than a
        // number: equal total doses given in different fraction sizes are
        // never isoeffective under LQ, so there is no ratio to report.
        let a = Schedule::conventional(30, 2.0);
        let b = Schedule::conventional(20, 3.0);
        assert!((a.total_dose_gy() - b.total_dose_gy()).abs() < 1e-9);
        assert!(isoeffect_alpha_beta(&a, &b).is_none());
    }

    #[test]
    fn the_repopulation_dose_per_day_lands_in_the_published_band() {
        // THE SECOND PREDICTION, and the one that is free to fail. alpha comes
        // from a survival parameterisation and the doubling time from a
        // repopulation study; neither was chosen to make this land anywhere,
        // and their combination is compared against a third quantity Withers
        // 1988 (PMID 3390344) reports directly.
        let d = repopulation_dose_per_day(ALPHA_HEAD_NECK_PER_GY, REPOP_DOUBLING_DAYS);
        let (lo, hi) = D_PROLIF_PUBLISHED_GY_PER_DAY;
        assert!(
            (lo..=hi).contains(&d),
            "D_prolif came out at {d:.2} Gy/day against a published {lo}-{hi}"
        );
    }

    #[test]
    fn the_prediction_is_sensitive_to_which_tumours_alpha_is_used() {
        // Asserted rather than mentioned, because it is the honest limit of
        // the test above: the glioblastoma alpha this module also carries puts
        // the same computation OUTSIDE the band. The prediction is about head
        // and neck cancer, and a reader who transplants it has been warned by
        // a failing test rather than by a sentence.
        let gbm = repopulation_dose_per_day(ALPHA_GBM_PARAMETERISATION_PER_GY, REPOP_DOUBLING_DAYS);
        let (_, hi) = D_PROLIF_PUBLISHED_GY_PER_DAY;
        assert!(
            gbm > hi,
            "the GBM alpha now lands inside the head-and-neck band ({gbm:.2}); \
             the sensitivity caveat this test pins has stopped being true"
        );
    }

    #[test]
    fn repopulation_costs_nothing_until_the_lag_has_passed() {
        let b = 72.0;
        let inside = bed_with_repopulation(
            b,
            20.0,
            REPOP_KICKOFF_DAYS,
            ALPHA_HEAD_NECK_PER_GY,
            REPOP_DOUBLING_DAYS,
        );
        assert!(
            (inside - b).abs() < 1e-12,
            "a short course paid a repopulation cost"
        );
        let beyond = bed_with_repopulation(
            b,
            48.0,
            REPOP_KICKOFF_DAYS,
            ALPHA_HEAD_NECK_PER_GY,
            REPOP_DOUBLING_DAYS,
        );
        assert!(
            beyond < b - 10.0,
            "a long course paid almost nothing: {beyond}"
        );
    }

    #[test]
    fn the_therapeutic_ratio_falls_as_fractions_grow() {
        // The whole reason radiotherapy is fractionated, as a monotone series
        // rather than a sentence. A model carrying only the tumour would
        // prefer one large dose; the late-responding tissue is what makes that
        // wrong, and the ratio has to fall for the layer to express it.
        let mut previous = f64::INFINITY;
        for (n, d) in [(30u32, 2.0f64), (20, 3.0), (5, 7.25), (1, 24.0)] {
            let r = therapeutic_ratio(
                &Schedule::conventional(n, d),
                ALPHA_BETA_TUMOUR_GY,
                ALPHA_BETA_LATE_GY,
            );
            assert!(r < previous, "{n}x{d} did not lower the ratio");
            previous = r;
        }
        // And it is exactly 1 at the schedule EQD2 is defined against, which
        // is what makes the others readable as a comparison.
        let conventional = therapeutic_ratio(
            &Schedule::conventional(30, 2.0),
            ALPHA_BETA_TUMOUR_GY,
            ALPHA_BETA_LATE_GY,
        );
        assert!((conventional - 1.0).abs() < 1e-9);
    }

    #[test]
    fn schedule_survival_reduces_to_the_single_dose_form() {
        let one = Schedule {
            n_fractions: 1,
            dose_per_fraction_gy: 2.0,
            interval_hours: 24.0,
            overall_days: 1.0,
        };
        let a = schedule_survival(&one, 0.3, 0.03, SUBLETHAL_REPAIR_HALF_TIME_H);
        let b = lq_survival(2.0, 0.3, 0.03);
        assert!((a - b).abs() < 1e-12, "{a} vs {b}");
    }

    #[test]
    fn the_reported_hypoxic_trajectory_is_the_one_the_model_follows() {
        // A MUTATION SURVIVOR of the third kind: `hypoxic_fraction_after`
        // reports the trajectory while `reoxygenating_schedule_survival`
        // advances the compartments by a per-fraction transfer rate, and
        // NOTHING tied the two together -- so the reporting function could be
        // broken outright without a test noticing. It is a public function
        // with no production caller, which is exactly the shape this
        // repository has been caught shipping before.
        //
        // With no killing at all, the compartment model reduces to pure
        // reoxygenation, and its hypoxic share must then BE the reported
        // trajectory. Any other decay law fails here.
        let f0 = 0.4;
        let half_life = 5.0;
        for k in [0u32, 1, 3, 10, 25] {
            let s = Schedule {
                n_fractions: k,
                dose_per_fraction_gy: 0.0,
                interval_hours: 24.0,
                overall_days: 1.0,
            };
            // Survival is 1 by construction at zero dose; what is being
            // compared is the SHARE, recovered by running the same transfer.
            let survived = reoxygenating_schedule_survival(
                &s,
                0.0,
                0.0,
                f0,
                half_life,
                0.05,
                crate::oxygen::OER_REFERENCE_PO2_MMHG,
            );
            assert!((survived - 1.0).abs() < 1e-12, "zero dose killed something");
            let reported = hypoxic_fraction_after(f0, k, half_life);
            let modelled = f0 * (-(f64::from(k)) * std::f64::consts::LN_2 / half_life).exp();
            assert!(
                (reported - modelled).abs() < 1e-12,
                "after {k} fractions the reported share was {reported} and \
                     the model's own transfer gives {modelled}"
            );
        }
        // And the half-life means what it says: half the compartment gone
        // after that many fractions.
        assert!((hypoxic_fraction_after(f0, 5, 5.0) - f0 / 2.0).abs() < 1e-12);
    }

    #[test]
    fn reoxygenation_cannot_matter_in_a_single_fraction() {
        // The control that makes the next test mean something: with one
        // fraction there is no interval to reoxygenate in, so the two arms
        // must agree EXACTLY rather than closely.
        let one = Schedule {
            n_fractions: 1,
            dose_per_fraction_gy: 24.0,
            interval_hours: 24.0,
            overall_days: 1.0,
        };
        let on = reoxygenating_schedule_survival(
            &one,
            0.3,
            0.03,
            0.3,
            5.0,
            0.05,
            crate::oxygen::OER_REFERENCE_PO2_MMHG,
        );
        let frozen = reoxygenating_schedule_survival(
            &one,
            0.3,
            0.03,
            0.3,
            f64::INFINITY,
            0.05,
            crate::oxygen::OER_REFERENCE_PO2_MMHG,
        );
        assert!((on - frozen).abs() < 1e-15, "{on} vs {frozen}");
    }

    #[test]
    fn a_reoxygenating_tumour_loses_more_cells_than_a_shielded_one() {
        let s = Schedule::conventional(30, 2.0);
        let p = crate::oxygen::OER_REFERENCE_PO2_MMHG;
        let reoxygenating = reoxygenating_schedule_survival(&s, 0.3, 0.03, 0.3, 5.0, 0.05, p);
        let frozen = reoxygenating_schedule_survival(&s, 0.3, 0.03, 0.3, f64::INFINITY, 0.05, p);
        assert!(reoxygenating < frozen, "{reoxygenating} vs {frozen}");
        assert!(
            frozen / reoxygenating > 5.0,
            "reoxygenation moved the result by less than five-fold, which is \
                 not the effect the fourth R is supposed to be"
        );
        // And the hypoxic compartment must cost something in the first place,
        // or the comparison above is between two identical models.
        let no_hypoxia = reoxygenating_schedule_survival(&s, 0.3, 0.03, 0.0, 5.0, 0.05, p);
        assert!(no_hypoxia < reoxygenating);
    }

    #[test]
    fn the_hypoxic_dose_is_scaled_down_and_not_up() {
        // The direction the first draft of this layer got WRONG, by dividing
        // two calls to the OER hyperbola by hand instead of using the module's
        // own scaler. A hypoxic cell is HARDER to kill, so the dose it behaves
        // as though it received is SMALLER.
        let p = crate::oxygen::OER_REFERENCE_PO2_MMHG;
        let hypoxic = oer_scaled_dose(2.0, 0.05, 1.0, p);
        assert!(
            hypoxic < 2.0,
            "hypoxic effective dose {hypoxic} was not reduced"
        );
        assert!(hypoxic > 0.0);
        let oxic = oer_scaled_dose(2.0, 1.0, 1.0, p);
        assert!(
            (oxic - 2.0).abs() < 1e-12,
            "a fully oxygenated cell lost dose"
        );
    }
}
