//! Simulation parameters for all models.

use serde::{Deserialize, Serialize};

use crate::photosensitizer_pk::Photosensitizer;

/// Default ramp time constant (steps) for NCOA4-ferritinophagy iron release
/// (#340). Inert while `ferritinophagy_release == 0.0`.
fn default_ferritinophagy_tau() -> f64 {
    30.0
}

/// Core biochemistry parameters. Identical to v3 simulation defaults.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Params {
    // === ROS Generation ===
    pub fenton_rate: f64,

    // === GSH Dynamics ===
    pub gsh_scav_efficiency: f64,
    pub gsh_km: f64,
    pub nrf2_gsh_rate: f64,

    // === Lipid Peroxidation ===
    pub lp_rate: f64,
    pub lp_propagation: f64,

    // === Repair ===
    pub gpx4_rate: f64,
    pub fsp1_rate: f64,
    /// DHODH (dihydroorotate dehydrogenase) GPX4-independent ferroptosis
    /// suppressor (#338): a mitochondrial CoQ10/ubiquinol axis that reduces
    /// lipid radicals in parallel to FSP1 (Mao et al., Nature 2021). Modeled as
    /// an additional GPX4-independent repair rate. `0.0` (default) means off,
    /// keeping the matrix byte-identical; uncalibrated, direction-anchored.
    /// Inhibiting it (rate back to 0) is the DHODHi combination strategy.
    #[serde(default)]
    pub dhodh_rate: f64,
    /// GCH1/BH4 (tetrahydrobiopterin) GPX4-independent radical-trapping
    /// antioxidant (#338): a lipid-radical quench capacity that gates the
    /// autocatalytic propagation switch in parallel to GPX4/FSP1 (Kraft 2020;
    /// Soula 2020). `0.0` (default) means off, byte-identical; uncalibrated,
    /// direction-anchored.
    #[serde(default)]
    pub gch1_rate: f64,
    /// SCD1-driven MUFA lipid-remodeling rate. In 3D culture and in vivo,
    /// SCD1 (regulated by SREBP1/mTORC1, not NRF2) converts SFA→MUFA,
    /// displacing PUFAs from membranes and reducing ferroptosis susceptibility.
    /// Zero in 2D culture (default); non-zero in in-vivo contexts.
    /// (Dixon/Park, Cancer Res 2025; Tesfay et al., Cancer Res 2019)
    pub scd_mufa_rate: f64,
    /// Maximum fraction of PUFA vulnerability suppressed by MUFA enrichment.
    /// Literature range: 0.40–0.60 (40–60% PUFA displacement in 3D/in vivo).
    pub scd_mufa_max: f64,
    /// Starting MUFA protection level. In established 3D/in-vivo tumors,
    /// SCD1-driven remodeling has already reached steady state, so cells
    /// begin with pre-accumulated membrane MUFA. Zero in 2D culture.
    pub initial_mufa_protection: f64,
    /// Acute-dosing MUFA-kinetics override (#339). `invivo()` starts cells at
    /// the established steady state (`initial_mufa_protection` = M_ss), which
    /// assumes the ~48-72h SCD1/MUFA enrichment is already complete. A freshly
    /// dosed ("acute") tumor has not built that up yet, so this overstates
    /// early ferroptosis resistance. When `Some(x)`, the cell instead STARTS at
    /// the naive value `x` (e.g. 0.0) and accumulates MUFA over the run via the
    /// existing logistic `update_mufa_protection`, breaking the steady-state
    /// assumption for acute treatment. `None` (default) keeps the existing
    /// `initial_mufa_protection` start, so the matrix is byte-identical. The
    /// accumulation RATE is the existing `scd_mufa_rate` (uncalibrated against
    /// the literature 48-72h timescale; the acute-vs-established DIRECTION is
    /// the result, not a precise hours figure). The starting value is not
    /// clamped here, but `update_mufa_protection` clamps to `[0, mufa_max]` on
    /// step 1 before any consumer reads it, so an out-of-range `x` self-corrects.
    #[serde(default)]
    pub mufa_acute_start: Option<f64>,
    /// MUFA decay rate from natural phospholipid turnover. When SCD1 is
    /// active, accumulation outpaces decay and protection reaches steady
    /// state. When SCD1 is inhibited (rate=0), decay gradually depletes
    /// existing membrane MUFA. Membrane lipid half-life ~24-48h.
    pub scd_mufa_decay: f64,
    /// Ether-linked PUFA pool as a fraction of the base PUFA substrate (#339).
    /// Polyunsaturated ether phospholipids (ether-PUFA-PE, made via FAR1/AGPS +
    /// the peroxisomal ether-lipid pathway) are an extra pool of oxidizable
    /// membrane lipid that PROMOTES ferroptosis (Zou 2020 Nature PMID 32939090;
    /// Cui 2021 Cell Death Differ PMID 33731874). The peroxidizable PUFA term is
    /// scaled by `1 + ether_pufa_fraction`, so enabling it raises lipid-peroxide
    /// accumulation; the `0` limit is the FAR1/AGPS-null escape (ether-lipid
    /// loss confers resistance). `0.0` (default) ⇒ ×1.0 ⇒ byte-identical;
    /// uncalibrated, direction-anchored. The plasmalogen/TMEM189 vinyl-ether
    /// sub-step is deliberately not modeled (contested sign; see
    /// `ether_augmented_pufa`).
    #[serde(default)]
    pub ether_pufa_fraction: f64,
    /// MBOAT1/2 hormone-regulated MUFA-enrichment boost (#339). MBOAT1
    /// (estrogen-receptor-regulated) and MBOAT2 (androgen-receptor-regulated)
    /// remodel phospholipids toward MUFA-PE and suppress ferroptosis
    /// independently of GPX4 (Liang et al., Cell 2023, PMID 37267948). Modeled
    /// as a constant additive MUFA protection layered onto the dynamic SCD1
    /// `mufa_protection` at the peroxidation sites (see `total_mufa_protection`),
    /// so a higher value (e.g. an AR-driven MBOAT2-high tumor) lowers
    /// ferroptosis. Floored at `0` (enrichment is protective-only). `0.0`
    /// (default) ⇒ byte-identical; uncalibrated, direction-anchored.
    #[serde(default)]
    pub mboat_mufa_boost: f64,
    /// NCOA4-ferritinophagy labile-iron release (#340). The static-iron model
    /// holds `cell.iron` fixed; in reality NCOA4-mediated autophagy of ferritin
    /// releases stored iron into the labile pool over time, feeding Fenton
    /// chemistry (Mancias et al., Nature 2014, PMID 24695223; Hou et al.,
    /// Autophagy 2016, PMID 27245739). This is the asymptotic fractional
    /// increase in labile iron the release builds toward over the run; the
    /// Fenton iron is scaled by a time-dependent factor (see
    /// `biochem::ferritinophagy_iron_factor`) that ramps from `1.0` toward
    /// `1 + ferritinophagy_release`. `0.0` (default) ⇒ factor is exactly `1.0`
    /// for every step ⇒ byte-identical; uncalibrated, direction-anchored
    /// (more ferritinophagy ⇒ more labile iron ⇒ more ferroptosis).
    #[serde(default)]
    pub ferritinophagy_release: f64,
    /// Ramp time constant (in steps) for the NCOA4-ferritinophagy iron release
    /// (#340). Inert when `ferritinophagy_release == 0.0` (the factor is then a
    /// constant `1.0`), so the default keeps the matrix byte-identical. A
    /// smaller value releases iron faster.
    #[serde(default = "default_ferritinophagy_tau")]
    pub ferritinophagy_tau: f64,

    // === GPX4 Dynamic Regulation ===
    pub gpx4_degradation_by_ros: f64,
    pub gpx4_nrf2_upregulation: f64,

    // === Treatment ===
    pub sdt_ros: f64,
    pub pdt_ros: f64,
    pub rsl3_gpx4_inhib: f64,

    // === GSH Homeostasis ===
    /// Maximum intracellular GSH (mM). Ref: ~10-12 mM in healthy cells.
    pub gsh_max: f64,

    // === GPX4 Target ===
    /// Multiplier for NRF2-driven GPX4 target level. GPX4_target = nrf2 * this value.
    pub gpx4_nrf2_target_multiplier: f64,

    // === Death ===
    pub death_threshold: f64,

    /// Number of simulation steps to continue after LP crosses the death
    /// threshold. During post-death steps, LP continues to accumulate via
    /// the autocatalytic chain reaction with zero repair (defenses have
    /// failed). This makes LP at death treatment-dependent: high-ROS
    /// treatments (SDT/PDT) drive LP to ~14 (5 steps), while slow-
    /// accumulation treatments (RSL3) barely exceed ~10.5.
    pub post_death_steps: u32,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            fenton_rate: 0.02,
            gsh_scav_efficiency: 0.5,
            gsh_km: 2.0,
            nrf2_gsh_rate: 0.025,
            lp_rate: 0.06,
            lp_propagation: 0.10,
            gpx4_rate: 0.30,
            fsp1_rate: 0.08,
            dhodh_rate: 0.0,
            gch1_rate: 0.0,
            scd_mufa_rate: 0.0,
            scd_mufa_max: 0.0,
            initial_mufa_protection: 0.0,
            mufa_acute_start: None,
            scd_mufa_decay: 0.0,
            ether_pufa_fraction: 0.0,
            mboat_mufa_boost: 0.0,
            ferritinophagy_release: 0.0,
            ferritinophagy_tau: default_ferritinophagy_tau(),
            gpx4_degradation_by_ros: 0.002,
            gpx4_nrf2_upregulation: 0.008,
            sdt_ros: 5.0,
            pdt_ros: 5.0,
            rsl3_gpx4_inhib: 0.92,
            gsh_max: 12.0,
            gpx4_nrf2_target_multiplier: 1.0,
            death_threshold: 10.0,
            post_death_steps: 5,
        }
    }
}

impl Params {
    /// In-vivo / 3D culture parameters with SCD1-driven MUFA protection enabled.
    ///
    /// Cells start at the accumulation–decay steady state (M_ss ≈ 0.40),
    /// representing established in-vivo lipid remodeling. The `scd_mufa_rate`
    /// maintains protection while `scd_mufa_decay` models natural phospholipid
    /// turnover. When SCD1 is inhibited (rate=0), existing MUFA decays.
    ///
    /// `scd_mufa_max: 0.50` caps PUFA displacement, consistent with
    /// Dixon/Park 2025 lipidomics (40–60% range) and Tesfay 2019 showing
    /// ~3–5× ferroptosis resensitization upon SCD1 inhibition.
    pub fn invivo() -> Self {
        // Steady state with decay: rate*(1-M/max) = decay*M
        // → M_ss = rate*max / (rate + decay*max) = 0.01*0.5 / (0.01 + 0.005*0.5) = 0.40
        // Cells start at this steady state.
        let rate = 0.01;
        let max = 0.50;
        let decay = 0.005;
        let steady_state = rate * max / (rate + decay * max);
        Params {
            scd_mufa_rate: rate,
            scd_mufa_max: max,
            scd_mufa_decay: decay,
            initial_mufa_protection: steady_state,
            ..Params::default()
        }
    }

    /// 3D spheroid context (#197): intermediate between `default` (2D culture,
    /// no MUFA) and `invivo` (full SCD1-driven MUFA, M_ss = 0.40). Spheroid
    /// cells have *partially* active MUFA remodeling, so the consumer's
    /// position-dependent per-cell MUFA (peripheral high, core low) survives
    /// the per-step `update_mufa_protection` homeostasis instead of being reset
    /// to 0 (the #265 footgun — possible only because `scd_mufa_max > 0` here).
    /// Half the in-vivo rate/cap ⇒ M_ss ≈ 0.20.
    pub fn spheroid() -> Self {
        let rate = 0.005;
        let max = 0.25;
        let decay = 0.005;
        let steady_state = rate * max / (rate + decay * max);
        Params {
            scd_mufa_rate: rate,
            scd_mufa_max: max,
            scd_mufa_decay: decay,
            initial_mufa_protection: steady_state,
            ..Params::default()
        }
    }
}

/// Spatial model parameters for energy deposition and diffusion.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SpatialParams {
    /// Cell diameter in micrometers.
    pub cell_size_um: f64,
    /// Iron diffusion coefficient in tissue (µm²/s). Estimate: free aqueous
    /// Fe²⁺ ≈ 700 µm²/s, scaled ~2.5× down for tissue tortuosity ⇒ ≈ 281 µm²/s.
    /// (Not a directly cited measurement; see parameter_provenance.md.)
    pub iron_diffusion_coeff: f64,
    /// Labile iron released per dead cell (µM equivalent).
    pub iron_release_per_death: f64,
    /// PDT effective attenuation coefficient (1/mm).
    /// At 630nm: µ_eff ≈ sqrt(3 × 0.3 × (0.3 + 10)) ≈ 3.1 /cm ≈ 0.31 /mm → δ ≈ 3.2mm
    /// Ref: Jacques SL, Phys Med Biol 2013
    pub pdt_mu_eff: f64,
    /// PDT incident fluence (relative units, 1.0 = standard dose).
    pub pdt_i0: f64,
    /// SDT ultrasound attenuation coefficient (dB/cm/MHz).
    /// Ref: soft tissue ≈ 0.7, muscle ≈ 1.3, fat ≈ 0.6
    pub sdt_alpha: f64,
    /// SDT ultrasound frequency (MHz). Typical: 1.0-3.0 MHz.
    pub sdt_freq_mhz: f64,
    /// SDT incident intensity (relative units, 1.0 = standard dose).
    pub sdt_i0: f64,
    /// Fraction of released iron reaching each neighbor cell.
    pub neighbor_iron_fraction: f64,
    /// Photosensitizer PK model. `Uniform(1.0)` (default) reproduces
    /// pre-PK PDT physics exactly. `Porfimer { t_half_h, t_distribution_h,
    /// phi_so2_relative }` enables drug-light-interval-aware scaling of
    /// the PDT dose, with optional saturating distribution phase and
    /// relative singlet-O₂ yield. `physics::pdt_intensity_at_depth`
    /// composes via `Photosensitizer::yield_at`. See the
    /// `photosensitizer_pk` module for the full kinetics + parser.
    #[serde(default)]
    pub photosensitizer: Photosensitizer,
    /// Hours from administration to light delivery, passed to
    /// `photosensitizer.yield_at(t_h)`.
    ///
    /// Interpretation depends on `Porfimer.t_distribution_h`:
    /// - With `t_distribution_h = 0` (default), this is interpreted as
    ///   time from *post-distribution peak* — same semantics as the
    ///   pre-#203 model.
    /// - With `t_distribution_h > 0`, the model holds drug at peak for
    ///   the first `t_distribution_h` hours after administration, then
    ///   begins exponential decay. So this field can be the **clinical
    ///   DLI from injection** directly (Bellnier 2006 reports porfimer
    ///   redistribution over ~24–48 h; setting `t_distribution_h ≈ 36`
    ///   approximates the absorption phase as a saturating step).
    ///
    /// Default 0.0 combined with the default `Photosensitizer::Uniform(1.0)`
    /// reproduces pre-PK PDT physics exactly.
    #[serde(default)]
    pub t_drug_light_interval_h: f64,
}

impl Default for SpatialParams {
    fn default() -> Self {
        SpatialParams {
            cell_size_um: 20.0,
            iron_diffusion_coeff: 281.0,
            iron_release_per_death: 2.0,
            pdt_mu_eff: 0.31, // 1/mm, δ ≈ 3.2mm at 630nm
            pdt_i0: 1.0,
            sdt_alpha: 0.7, // dB/cm/MHz (soft tissue)
            sdt_freq_mhz: 1.0,
            sdt_i0: 1.0,
            neighbor_iron_fraction: 0.1,
            photosensitizer: Photosensitizer::default(),
            t_drug_light_interval_h: 0.0,
        }
    }
}

/// Immune cascade parameters for ICD modeling. Used by sim-icd + sim-combo.
///
/// Models the full DC→T-cell cascade with separate `dc_maturation_rate`,
/// `tcell_priming_rate`, `tcell_kill_rate` steps.
///
/// **See also [`SpatialImmuneConfig`]** for the spatial-DAMP-field variant
/// used by sim-tme + sim-tme-3d (single absorbed `immune_kill_rate`
/// plus `damp_diffusion_fraction` / `damp_clearance_rate`). Same biology,
/// different math; both valid. If you add a field here, check whether
/// the spatial variant needs the same change.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ImmuneParams {
    /// DAMP release proportional to LP at death.
    pub damp_per_lp: f64,
    /// DAMP concentration for half-maximal DC activation (Michaelis-Menten Kd).
    pub dc_activation_kd: f64,
    /// Fraction of activated DCs that successfully mature.
    pub dc_maturation_rate: f64,
    /// T cells primed per mature DC.
    pub tcell_priming_rate: f64,
    /// Tumor cells killed per primed T cell per cycle.
    pub tcell_kill_rate: f64,
    /// PD-1 brake: fraction of T-cell kill suppressed (0.0 = no brake, 1.0 = full suppression).
    pub pd1_brake: f64,
    /// Anti-PD-1 efficacy: fraction of PD-1 brake removed (0.0 = no drug, 1.0 = complete blockade).
    pub anti_pd1_efficacy: f64,
}

impl Default for ImmuneParams {
    fn default() -> Self {
        ImmuneParams {
            damp_per_lp: 1.0,
            dc_activation_kd: 50.0,
            dc_maturation_rate: 0.6,
            tcell_priming_rate: 10.0,
            tcell_kill_rate: 3.0,
            pd1_brake: 0.7,
            anti_pd1_efficacy: 0.8,
        }
    }
}

/// Spatial immune-coupling params used by sim-tme (2D) and sim-tme-3d (3D).
///
/// Distinct from [`ImmuneParams`] — that struct models the full DC→T-cell
/// cascade with separate maturation/priming/kill rates. This struct uses
/// a single absorbed `immune_kill_rate` plus spatial fields
/// (`damp_diffusion_fraction`, `damp_clearance_rate`) for the spatial
/// DAMP-field model. Same biology, different math; both are valid.
///
/// **2D vs 3D split**: `damp_diffusion_fraction` differs by geometry.
/// 2D uses 0.08 (8-Moore neighbors via `TumorGrid::neighbors`, stability
/// bound `< 1 / 8 = 0.125`); 3D uses 0.025 (26-Moore neighbors, stability
/// bound `< 1 / 26 ≈ 0.0385` enforced by `assert!` in
/// `immune_spatial::diffuse_damp_3d_step`). Use [`for_2d()`] or [`for_3d()`]
/// to pick the right default.
///
/// **No `Default` impl on purpose**: callers MUST pick a geometry. Using
/// the wrong default in 3D would panic at runtime via the diffuse-step
/// stability `assert!`; making that a compile error via explicit
/// constructors is preferable. Compare with [`StromalConfig`] /
/// [`PhConfig`], which do impl `Default` because they're
/// geometry-independent.
#[derive(Clone, Copy, Debug)]
pub struct SpatialImmuneConfig {
    /// DAMP released per unit LP at death.
    pub damp_per_lp: f64,
    /// Fraction of DAMP shared with each Moore neighbor per step. The
    /// total redistributed mass per source cell is `fraction * neighbor_count`
    /// and must stay below 1 to keep the diffusion update mass-conserving.
    /// 2D: 0.08 (×8 Moore = up to 0.64 redistributed; bound `< 1/8`).
    /// 3D: 0.025 (×26 Moore = up to 0.65 redistributed; bound `< 1/26`).
    pub damp_diffusion_fraction: f64,
    /// Exponential decay per step (models immune clearance of DAMPs).
    pub damp_clearance_rate: f64,
    /// Michaelis-Menten Kd for DC activation by DAMP concentration.
    pub dc_activation_kd: f64,
    /// Per-step immune kill rate (absorbs DC maturation + T cell priming + kill).
    pub immune_kill_rate: f64,
    /// PD-1 brake: fraction of T-cell kill suppressed.
    pub pd1_brake: f64,
    /// Anti-PD-1 efficacy: fraction of PD-1 brake removed.
    pub anti_pd1_efficacy: f64,
    /// T-cell exhaustion rate (#243, Phase 1): cumulative immune kills in a
    /// cell's Moore neighborhood suppress its further kill probability by
    /// `1/(1 + exhaustion_rate · cumulative)` (see
    /// [`crate::immune_spatial::exhaustion_factor`]). `0.0` (the default for
    /// both `for_2d`/`for_3d`) disables exhaustion, keeping output
    /// byte-identical to the pre-#243 single-PD-1-brake model.
    ///
    /// **3D-only in Phase 1**: currently consumed only by `sim-tme-3d`. A 2D
    /// (`sim-tme`) caller that sets this > 0 would silently get no effect
    /// until the 2D immune loop is wired up in a later phase.
    pub exhaustion_rate: f64,
}

impl SpatialImmuneConfig {
    /// 2D default (sim-tme): `damp_diffusion_fraction = 0.08`, safely
    /// below the `1/8` stability bound for the 8-Moore neighborhood used
    /// by `TumorGrid::neighbors`.
    pub fn for_2d() -> Self {
        SpatialImmuneConfig {
            damp_per_lp: 1.0,
            damp_diffusion_fraction: 0.08,
            damp_clearance_rate: 0.03,
            dc_activation_kd: 50.0,
            immune_kill_rate: 0.02,
            pd1_brake: 0.7,
            anti_pd1_efficacy: 0.0,
            exhaustion_rate: 0.0,
        }
    }

    /// 3D default (sim-tme-3d): `damp_diffusion_fraction = 0.025` for
    /// 26-Moore (matches `immune_spatial::diffuse_damp_3d_step`'s
    /// `assert!(0.025 * 26.0 < 1.0)` stability invariant).
    pub fn for_3d() -> Self {
        SpatialImmuneConfig {
            damp_per_lp: 1.0,
            damp_diffusion_fraction: 0.025,
            damp_clearance_rate: 0.03,
            dc_activation_kd: 50.0,
            immune_kill_rate: 0.02,
            pd1_brake: 0.7,
            anti_pd1_efficacy: 0.0,
            exhaustion_rate: 0.0,
        }
    }

    /// Return a copy with anti-PD-1 blockade applied (efficacy 0.8).
    pub fn with_anti_pd1(&self) -> Self {
        SpatialImmuneConfig {
            anti_pd1_efficacy: 0.8,
            ..*self
        }
    }

    /// Net PD-1 brake after anti-PD-1 modulation: `pd1_brake * (1 - anti_pd1_efficacy)`.
    pub fn effective_brake(&self) -> f64 {
        self.pd1_brake * (1.0 - self.anti_pd1_efficacy)
    }
}

/// CAF-mediated stromal protection params used by sim-tme + sim-tme-3d.
///
/// All ESTIMATED. Refs: PMID 34373744 (CAF metabolic reprogramming),
/// PMID 31813804 (ACSL3-mediated oleic acid), PMID 30842648 (MUFA
/// ferroptosis). Geometry-independent — no 2D/3D split (per-cell
/// shielding measured equal at 50.0% (2D) and 51.5% (3D) in the
/// PR #221 validation report).
#[derive(Clone, Copy, Debug)]
pub struct StromalConfig {
    /// Per-step GSH boost for stromal-adjacent tumor cells.
    pub gsh_boost_per_step: f64,
    /// Maximum GSH from CAF supply (1.5× normal gsh_max of 12.0).
    pub gsh_boost_cap: f64,
    /// Per-step MUFA boost from ACSL3-mediated oleic acid uptake.
    pub mufa_boost_per_step: f64,
    /// Maximum MUFA from CAF supply.
    pub mufa_boost_cap: f64,
}

impl Default for StromalConfig {
    fn default() -> Self {
        StromalConfig {
            gsh_boost_per_step: 0.06,
            gsh_boost_cap: 18.0,
            mufa_boost_per_step: 0.003,
            mufa_boost_cap: 0.25,
        }
    }
}

/// Tumor acidic pH gradient params (Warburg effect lactic acid + ion-trapping).
///
/// Used by sim-tme + sim-tme-3d. Linearized Henderson-Hasselbalch over
/// pH 6.5-7.4; valid only in that narrow window. Refs: Stubbs 2000,
/// Gatenby & Gillies 2004 (tumor pH); Harrison & Arosio 1996 (ferritin
/// iron release). Geometry-independent (same chemistry in 2D and 3D).
#[derive(Clone, Copy, Debug)]
pub struct PhConfig {
    /// pH at tumor edge (well-perfused).
    pub ph_edge: f64,
    /// pH at deep tumor core (Warburg lactic acid).
    pub ph_core: f64,
    /// pH penetration length (μm); matches O2 reference λ.
    pub lambda_ph_um: f64,
    /// Iron-pH sensitivity: ferritin releases Fe²⁺ at low pH.
    pub iron_ph_sensitivity: f64,
    /// Ion trapping sensitivity for weak-base drugs (RSL3).
    pub ion_trap_sensitivity: f64,
}

impl Default for PhConfig {
    fn default() -> Self {
        PhConfig {
            ph_edge: 7.4,
            ph_core: 6.5,
            lambda_ph_um: 120.0,
            iron_ph_sensitivity: 1.5,
            ion_trap_sensitivity: 0.4,
        }
    }
}

/// Drug-tolerant persister-cell parameters (#241).
///
/// Cells acquire an epigenetic ferroptosis-tolerant state under drug exposure
/// and revert once the drug clears. Consumed by [`crate::persister`] (pure
/// helpers; the consumer owns `CellState::persister_fraction` and applies the
/// effects). [`Default`] is the **identity element** (all rates zero ⇒ every
/// helper is a no-op ⇒ byte-identical to having no persister model), matching
/// the off-by-default discipline of `DoseSchedule` / `Photosensitizer`.
///
/// Values in [`PersisterConfig::enabled`] are plausible placeholders pending
/// calibration (the literature gives qualitative direction, not step-level
/// rates): Hangauer 2017 (persister ⇄ GPX4 dependence), Tsoi 2018 (MUFA lipid
/// rewiring), Viswanathan 2017 (mesenchymal ⇄ ferroptosis axis / reversion).
#[derive(Clone, Copy, Debug)]
pub struct PersisterConfig {
    /// Logistic acquisition rate per step under full drug exposure.
    pub acquisition_rate: f64,
    /// Exponential reversion rate per step when the drug is absent.
    pub reversion_rate: f64,
    /// Ceiling on `persister_fraction` (no cell becomes fully invulnerable).
    pub max_fraction: f64,
    /// Fraction of per-step GPX4 inactivation a full persister resists (0..1).
    pub gpx4_resistance: f64,
    /// Per-step additive MUFA protection at full persistence (same shape as
    /// `StromalConfig::mufa_boost_per_step`).
    pub mufa_boost_per_step: f64,
    /// Cap on persister-driven `CellState::mufa_protection`.
    pub mufa_boost_cap: f64,
}

impl Default for PersisterConfig {
    /// Identity: no acquisition, no reversion, no effect. A run with this
    /// config is byte-identical to one with no persister model.
    fn default() -> Self {
        PersisterConfig {
            acquisition_rate: 0.0,
            reversion_rate: 0.0,
            max_fraction: 0.0,
            gpx4_resistance: 0.0,
            mufa_boost_per_step: 0.0,
            mufa_boost_cap: 0.0,
        }
    }
}

impl PersisterConfig {
    /// Plausible (placeholder, pending calibration) values that produce an
    /// observable persister effect when the model is switched on.
    pub fn enabled() -> Self {
        PersisterConfig {
            acquisition_rate: 0.02,
            reversion_rate: 0.01,
            max_fraction: 0.8,
            gpx4_resistance: 0.5,
            mufa_boost_per_step: 0.004,
            mufa_boost_cap: 0.6,
        }
    }

    /// True when this config has no effect (every rate zero). Available for a
    /// consumer that holds a `PersisterConfig` directly and wants to skip the
    /// persister path; note `sim-tme-3d` instead gates on
    /// `Option<PersisterConfig>` (`None` vs `Some`) and `persister::acquire`
    /// short-circuits on `acquisition_rate == 0.0`. Omits `max_fraction` /
    /// `mufa_boost_cap` (caps are inert when all rates are zero). Since the
    /// fields are `pub` and the rate/cap pairing is not enforced, a hand-built
    /// `{ acquisition_rate: 0.0, max_fraction: 0.8 }` silently disables
    /// acquisition (the short-circuit fires) — a safe direction, but use
    /// `default()` / `enabled()` rather than partial literals.
    pub fn is_identity(&self) -> bool {
        self.acquisition_rate == 0.0
            && self.reversion_rate == 0.0
            && self.gpx4_resistance == 0.0
            && self.mufa_boost_per_step == 0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// SpatialParams JSON written before the photosensitizer/DLI fields
    /// existed must still deserialize, with new fields filled by serde
    /// defaults that reproduce pre-PK behavior.
    #[test]
    fn spatial_params_legacy_json_deserializes_with_defaults() {
        let legacy = r#"{
            "cell_size_um": 20.0,
            "iron_diffusion_coeff": 281.0,
            "iron_release_per_death": 2.0,
            "pdt_mu_eff": 0.31,
            "pdt_i0": 1.0,
            "sdt_alpha": 0.7,
            "sdt_freq_mhz": 1.0,
            "sdt_i0": 1.0,
            "neighbor_iron_fraction": 0.1
        }"#;
        let p: SpatialParams = serde_json::from_str(legacy).unwrap();
        assert_eq!(p.photosensitizer, Photosensitizer::Uniform(1.0));
        assert_eq!(p.t_drug_light_interval_h, 0.0);
        assert_eq!(
            p.photosensitizer
                .concentration_at(p.t_drug_light_interval_h),
            1.0
        );
    }

    /// `SpatialImmuneConfig::for_2d()` must reproduce the literal values
    /// sim-tme used before the lift (#220). Regression-guards a silent
    /// drift if anyone edits the constructor.
    #[test]
    fn spatial_immune_for_2d_matches_sim_tme_legacy_defaults() {
        let c = SpatialImmuneConfig::for_2d();
        assert_eq!(c.damp_per_lp, 1.0);
        assert_eq!(c.damp_diffusion_fraction, 0.08);
        assert_eq!(c.damp_clearance_rate, 0.03);
        assert_eq!(c.dc_activation_kd, 50.0);
        assert_eq!(c.immune_kill_rate, 0.02);
        assert_eq!(c.pd1_brake, 0.7);
        assert_eq!(c.anti_pd1_efficacy, 0.0);
        // 2D stability invariant: sim-tme's DAMP diffusion calls
        // `TumorGrid::neighbors` which returns up to 8 Moore neighbors.
        // Each step subtracts `fraction * count` from the source cell, so
        // `fraction * 8 < 1` is required to keep the update mass-conserving.
        assert!(c.damp_diffusion_fraction * 8.0 < 1.0);
    }

    /// `SpatialImmuneConfig::for_3d()` must reproduce the literal values
    /// sim-tme-3d used before the lift, in particular the 3D-safe
    /// `damp_diffusion_fraction = 0.025` (×26 Moore < 1 stability bound).
    #[test]
    fn spatial_immune_for_3d_matches_sim_tme_3d_legacy_defaults() {
        let c = SpatialImmuneConfig::for_3d();
        assert_eq!(c.damp_per_lp, 1.0);
        assert_eq!(c.damp_diffusion_fraction, 0.025);
        assert_eq!(c.damp_clearance_rate, 0.03);
        assert_eq!(c.dc_activation_kd, 50.0);
        assert_eq!(c.immune_kill_rate, 0.02);
        assert_eq!(c.pd1_brake, 0.7);
        assert_eq!(c.anti_pd1_efficacy, 0.0);
        // Stability invariant matched by immune_spatial::diffuse_damp_3d_step's
        // `assert!`. If this ever fails, that diffusion step panics in
        // release mode — the test exists to catch the drift first.
        assert!(c.damp_diffusion_fraction * 26.0 < 1.0);
    }

    /// `with_anti_pd1()` sets `anti_pd1_efficacy = 0.8` and leaves the
    /// rest unchanged. Matches sim-tme's `ImmuneConfig::with_anti_pd1`.
    #[test]
    fn spatial_immune_with_anti_pd1_changes_only_efficacy() {
        let base = SpatialImmuneConfig::for_2d();
        let blocked = base.with_anti_pd1();
        assert_eq!(blocked.anti_pd1_efficacy, 0.8);
        // Everything else must match.
        assert_eq!(blocked.damp_per_lp, base.damp_per_lp);
        assert_eq!(
            blocked.damp_diffusion_fraction,
            base.damp_diffusion_fraction
        );
        assert_eq!(blocked.damp_clearance_rate, base.damp_clearance_rate);
        assert_eq!(blocked.dc_activation_kd, base.dc_activation_kd);
        assert_eq!(blocked.immune_kill_rate, base.immune_kill_rate);
        assert_eq!(blocked.pd1_brake, base.pd1_brake);
    }

    /// `effective_brake` = `pd1_brake * (1 - anti_pd1_efficacy)`. With
    /// no anti-PD-1: full 0.7 brake. With anti-PD-1 (0.8 efficacy):
    /// 0.7 × 0.2 = 0.14.
    #[test]
    fn spatial_immune_effective_brake_math() {
        let base = SpatialImmuneConfig::for_2d();
        assert!((base.effective_brake() - 0.7).abs() < 1e-12);
        let blocked = base.with_anti_pd1();
        assert!((blocked.effective_brake() - 0.14).abs() < 1e-12);
    }

    /// `StromalConfig::default()` reproduces the literal values that
    /// sim-tme + sim-tme-3d used before the lift.
    #[test]
    fn stromal_default_matches_legacy() {
        let c = StromalConfig::default();
        assert_eq!(c.gsh_boost_per_step, 0.06);
        assert_eq!(c.gsh_boost_cap, 18.0);
        assert_eq!(c.mufa_boost_per_step, 0.003);
        assert_eq!(c.mufa_boost_cap, 0.25);
    }

    /// `PhConfig::default()` reproduces the literal values that
    /// sim-tme + sim-tme-3d used before the lift.
    #[test]
    fn ph_default_matches_legacy() {
        let c = PhConfig::default();
        assert_eq!(c.ph_edge, 7.4);
        assert_eq!(c.ph_core, 6.5);
        assert_eq!(c.lambda_ph_um, 120.0);
        assert_eq!(c.iron_ph_sensitivity, 1.5);
        assert_eq!(c.ion_trap_sensitivity, 0.4);
    }
}
