//! Simulation parameters for all models.

use serde::{Deserialize, Serialize};

use crate::photosensitizer_pk::Photosensitizer;

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
    /// MUFA decay rate from natural phospholipid turnover. When SCD1 is
    /// active, accumulation outpaces decay and protection reaches steady
    /// state. When SCD1 is inhibited (rate=0), decay gradually depletes
    /// existing membrane MUFA. Membrane lipid half-life ~24-48h.
    pub scd_mufa_decay: f64,

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
            scd_mufa_rate: 0.0,
            scd_mufa_max: 0.0,
            initial_mufa_protection: 0.0,
            scd_mufa_decay: 0.0,
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
}

/// Spatial model parameters for energy deposition and diffusion.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SpatialParams {
    /// Cell diameter in micrometers.
    pub cell_size_um: f64,
    /// Iron diffusion coefficient in tissue (µm²/s). Ref: D_eff ≈ 281 µm²/s
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
            pdt_mu_eff: 0.31,       // 1/mm, δ ≈ 3.2mm at 630nm
            pdt_i0: 1.0,
            sdt_alpha: 0.7,          // dB/cm/MHz (soft tissue)
            sdt_freq_mhz: 1.0,
            sdt_i0: 1.0,
            neighbor_iron_fraction: 0.1,
            photosensitizer: Photosensitizer::default(),
            t_drug_light_interval_h: 0.0,
        }
    }
}

/// Immune cascade parameters for ICD modeling.
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
        assert_eq!(p.photosensitizer.concentration_at(p.t_drug_light_interval_h), 1.0);
    }
}
