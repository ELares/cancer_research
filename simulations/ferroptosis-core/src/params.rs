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
    /// PROM2 / MVB-exosome labile-iron EFFLUX (#484), in `[0, 1]`. Pro-ferroptotic
    /// stress induces Prominin-2, which packages ferritin-bound iron into
    /// multivesicular bodies secreted as exosomes, DEPLETING the labile iron pool
    /// and starving the Fenton reaction (Brown et al., Dev Cell 2019, PMID
    /// 31761539). The OPPOSITE sign to `ferritinophagy_release`: the Fenton iron
    /// is scaled by a time-dependent factor (see `biochem::prom2_iron_factor`)
    /// that ramps from `1.0` toward `1 - prom2_iron_efflux` with the shared
    /// `ferritinophagy_tau` time constant, so a PROM2-high cell exports iron over
    /// the run and RESISTS ferroptosis. `0.0` (default) ⇒ factor is exactly `1.0`
    /// for every step ⇒ byte-identical; uncalibrated, direction-anchored (more
    /// PROM2 efflux ⇒ less labile iron ⇒ less ferroptosis). FFI defaults it to
    /// 0.0 so the C ABI is unchanged.
    #[serde(default)]
    pub prom2_iron_efflux: f64,
    /// ALOX (lipoxygenase) isoform-specific propagation boost (#446). The
    /// autocatalytic lipid-peroxidation propagation rate is multiplied by
    /// `1 + alox_propagation_boost` (clamped `>= 0`), so an ALOX15/12/5-high
    /// tumor peroxidizes faster (`> 0` ⇒ more ferroptosis) and an ALOX-poor
    /// tumor peroxidizes slower (`< 0`, down to `-1` ⇒ the ALOX-null limit, no
    /// enzymatic propagation), independent of the GPX4/GSH/FSP1 defenses
    /// (lipoxygenase-driven ferroptosis, PNAS 2016 PMID 27506793). A consumer
    /// computes this from an isoform mix via [`crate::alox::AloxConfig`].
    /// `0.0` (default) ⇒ ×1.0 ⇒ byte-identical; uncalibrated, direction-anchored.
    #[serde(default)]
    pub alox_propagation_boost: f64,
    /// MCFA (medium-chain fatty acid) → ACSL4/CD36 PUFA-incorporation boost
    /// (#446). Added to the oxidizable-PUFA augmentation alongside the
    /// ether-lipid pool (see `biochem::ether_augmented_pufa`), so MCFA exposure
    /// raises the peroxidizable substrate and thus ferroptosis susceptibility
    /// (Sci Rep 2024 s41598-024-55050-4; MCFA ferroptosis sensitization
    /// PMC11901882). Computed from [`crate::alox::AloxConfig`]. `0.0` (default)
    /// ⇒ ×1.0 ⇒ byte-identical; uncalibrated, direction-anchored.
    #[serde(default)]
    pub mcfa_pufa_boost: f64,
    /// ACSL4-status biomarker PUFA-incorporation boost (#444). ACSL4 ligates the
    /// PUFA that ferroptosis requires into membranes, so a tumor's ACSL4
    /// expression status sets its oxidizable-PUFA baseline: this additive boost
    /// (computed from a status scalar via `crate::acsl4::pufa_boost_from_status`,
    /// `status - 1` clamped at `-1`) is folded into the oxidizable-PUFA
    /// augmentation (`biochem::ether_augmented_pufa`). `0.0` (default, wild-type
    /// ACSL4) ⇒ ×1.0 ⇒ byte-identical; `-1` is the ACSL4-negative null floor
    /// (PUFA substrate collapses ⇒ ferroptosis-refractory through a mechanism
    /// distinct from GPX4/GSH/FSP1, e.g. some HCC/AML subtypes); `> 0` is
    /// ACSL4-high (more PUFA ⇒ sensitive). Doll et al., Nat Chem Biol 2017
    /// (PMID 27842070). Uncalibrated linear placeholder; direction-anchored.
    #[serde(default)]
    pub acsl4_status_boost: f64,
    /// Exogenous dietary-PUFA supply (#486): a double-bond-weighted oxidizable
    /// substrate added on top of the membrane PUFA. Exogenous n-3/n-6 PUFAs
    /// sensitize to ferroptosis in proportion to double-bond count, BUT only
    /// once the saturable lipid-droplet (triglyceride) storage sink below is
    /// exceeded (Dierge et al., Cell Metab 2021, PMID 34118189). The effective
    /// contribution folded into `biochem::ether_augmented_pufa` is
    /// `(dietary_pufa_supply - lipid_droplet_buffer).max(0)`. `0.0` (default) ⇒
    /// no contribution ⇒ byte-identical. FFI defaults it to 0.0 so the C ABI is
    /// unchanged. Uncalibrated placeholder; the direction (more dietary PUFA ⇒
    /// more ferroptosis, above the buffer) is the result.
    #[serde(default)]
    pub dietary_pufa_supply: f64,
    /// Lipid-droplet / DGAT triglyceride-storage buffer (#486): the saturable
    /// sink that must fill before `dietary_pufa_supply` raises peroxidation
    /// (esterifying exogenous PUFA into stored triglycerides is protective until
    /// storage saturates). DGAT inhibition LOWERS this buffer, so dietary-PUFA
    /// cytotoxicity emerges sooner (the DGATi synergy, Dierge 2021 PMID
    /// 34118189). Floored at 0. `0.0` (default) ⇒ the dietary-PUFA excess is the
    /// full supply (no buffering), but with `dietary_pufa_supply = 0.0` the net
    /// is still 0 ⇒ byte-identical. FFI defaults it to 0.0; not in the C ABI.
    #[serde(default)]
    pub lipid_droplet_buffer: f64,
    /// ESCRT-III membrane-repair rescue rate (#465): the per-step probability that
    /// a cell whose lipid peroxide has crossed `death_threshold` is resealed by
    /// ESCRT-III and survives that step (instead of dying), as long as repair
    /// budget remains. Acts on the death-EXECUTION step, not the lipid substrate
    /// (Dai et al., BBRC 2020, PMID 31761326: CHMP5/CHMP6 membrane repair blocks
    /// ferroptosis; knockdown sensitizes). The RNG roll is drawn only when this is
    /// `> 0`, so `0.0` (default) ⇒ the brake never fires ⇒ byte-identical.
    #[serde(default)]
    pub escrt_repair_rate: f64,
    /// ESCRT-III finite per-cell repair CAPACITY (#465): the number of rescue
    /// events available before the machinery is exhausted (tracked per cell on
    /// `CellState::escrt_budget_used`). A larger budget means a longer death delay
    /// / more resistance. `0.0` (default) ⇒ no rescues possible ⇒ byte-identical.
    #[serde(default)]
    pub escrt_repair_budget: f64,
    /// POR/CYB5R1 enzymatic NAD(P)H/O2-driven H2O2 source (#466). POR (cytochrome
    /// P450 reductase) and CYB5R1 transfer electrons from NAD(P)H to O2 to generate
    /// H2O2, the Fenton substrate that drives lipid peroxidation, an enzymatic,
    /// O2- and NADPH-dependent oxidant source distinct from the static `basal_ros`
    /// input and parallel to the ALOX enzymatic-propagation leg (Yan et al., 2021,
    /// PMID 33860083: POR/CYB5R1 catalyze lipid peroxidation to execute ferroptosis;
    /// Zou et al., Nat Chem Biol 2020). This rate is added to `total_ros` (more POR
    /// ⇒ more H2O2 ⇒ more Fenton-driven ROS ⇒ more ferroptosis). The single-cell
    /// term is uniform; the spatial consumer (sim-tme-3d) O2-couples it per cell via
    /// `oxygen::por_o2_factor`, which ties the H2O2 yield to local O2 and so helps
    /// correct the deep-core artifact (POR makes less H2O2 where O2 is low). `0.0`
    /// (default) ⇒ no added oxidant ⇒ byte-identical; FFI defaults it to 0.0 so the
    /// C ABI is unchanged.
    #[serde(default)]
    pub por_h2o2_rate: f64,
    /// 7-DHC sterol radical-trapping antioxidant pool (#467). 7-dehydrocholesterol
    /// (7-DHC), a distal cholesterol-biosynthesis intermediate, is a potent
    /// endogenous membrane-resident radical-trapping antioxidant that shields
    /// PUFA from peroxyl-radical autoxidation (like membrane vitamin E), gating the
    /// autocatalytic peroxidation chain (Freitas et al. / Li et al., Nature 2024,
    /// "7-dehydrocholesterol dictates ferroptosis sensitivity," PMID 38297130,
    /// DOI 10.1038/s41586-023-06983-9). It is added to the GPX4-independent
    /// radical-trapping quench term (alongside `gch1_rate`/FSP1), so a higher pool
    /// raises the quench and LOWERS the propagation rate ⇒ ferroptosis RESISTANCE.
    /// DHCR7 consumes 7-DHC, so DHCR7-loss raises this pool (resistance, the modeled
    /// escape) and EBP/SC5D-low lowers it (sensitization). This is a sterol-pathway
    /// radical-trapping defense distinct from GPX4/GSH/FSP1/DHODH and from the
    /// iron/MUFA axes. `0.0` (default) ⇒ no added quench ⇒ byte-identical; FFI
    /// defaults it to 0.0 so the C ABI is unchanged.
    #[serde(default)]
    pub dhc7_radical_trap: f64,

    /// Vitamin K / VKORC1L1 radical-trapping defense (#483), a SIXTH
    /// GPX4-independent ferroptosis-suppressor axis alongside `gch1_rate` (BH4),
    /// `dhc7_radical_trap` (7-DHC), and FSP1 (CoQ). VKORC1L1 reduces vitamin K to
    /// a radical-trapping antioxidant (vitamin K hydroquinone) that quenches
    /// phospholipid peroxyl radicals INDEPENDENT of GSH/GPX4; VKORC1L1 is a p53
    /// transcriptional target (Yang et al., Cell Metab 2023, PMID 37467745;
    /// mechanism origin Mishima et al., Nature 2022, PMID 35922516). Added to the
    /// `antioxidant_quench` pool, so a higher trap LOWERS the propagation rate ⇒
    /// ferroptosis RESISTANCE. The effective trap is reduced by
    /// `warfarin_vkor_inhibition` (the druggable knob below), via
    /// [`Params::effective_vitk_radical_trap`]. `0.0` (default) ⇒ no added quench
    /// ⇒ byte-identical; FFI defaults it to 0.0 so the C ABI is unchanged.
    #[serde(default)]
    pub vitk_radical_trap: f64,
    /// Warfarin inhibition of VKORC1L1 (#483), in `[0, 1]`. The FDA anticoagulant
    /// warfarin inhibits VKOR/VKORC1L1, collapsing the vitamin-K radical-trap and
    /// thereby DRIVING ferroptosis (a repurposable, p53-status-gated axis; Yang
    /// et al. 2023 PMID 37467745 showed warfarin suppresses tumors via this
    /// mechanism in immunocompetent mice). Scales the effective trap DOWN:
    /// `effective = vitk_radical_trap * (1 - warfarin_vkor_inhibition)`, so `1.0`
    /// fully removes the VKORC1L1 protection. Inert when `vitk_radical_trap == 0`.
    /// `0.0` (default) ⇒ no inhibition ⇒ byte-identical; not in the C ABI.
    #[serde(default)]
    pub warfarin_vkor_inhibition: f64,

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
            prom2_iron_efflux: 0.0,
            alox_propagation_boost: 0.0,
            mcfa_pufa_boost: 0.0,
            acsl4_status_boost: 0.0,
            dietary_pufa_supply: 0.0,
            lipid_droplet_buffer: 0.0,
            escrt_repair_rate: 0.0,
            escrt_repair_budget: 0.0,
            por_h2o2_rate: 0.0,
            dhc7_radical_trap: 0.0,
            vitk_radical_trap: 0.0,
            warfarin_vkor_inhibition: 0.0,
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
    /// Effective vitamin-K / VKORC1L1 radical-trapping quench (#483) after
    /// warfarin inhibition: `vitk_radical_trap * (1 - warfarin_vkor_inhibition)`,
    /// with `warfarin_vkor_inhibition` clamped to `[0, 1]` and the result floored
    /// at `0`. `0.0` when `vitk_radical_trap == 0` (the default) ⇒ the
    /// `antioxidant_quench` is unchanged ⇒ byte-identical.
    #[must_use]
    pub fn effective_vitk_radical_trap(&self) -> f64 {
        (self.vitk_radical_trap * (1.0 - self.warfarin_vkor_inhibition.clamp(0.0, 1.0))).max(0.0)
    }

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
    /// Immunosuppressive-ferroptosis strength (#337). The model treats
    /// ferroptotic ICD as net pro-immune (DAMP -> DC -> CD8 kill), but in vivo
    /// ferroptosis is frequently net IMMUNOSUPPRESSIVE: dying cells co-release
    /// DC-suppressing factors (extracellular GPX4 binding DC ZP3, Liu et al.
    /// Cell 2026 PMID 41494530; oxidized lipids / PGE2, Kim et al. Nature 2022
    /// PMID 36385526) that blunt DC maturation and CD8 priming (Wiernicki et al.
    /// Nat Commun 2022 PMID 35760796). Modeled as a per-cell kill multiplier
    /// `1/(1 + strength · local_damp)` keyed on the SAME local ferroptotic-death
    /// /DAMP signal that drives pro-immune `dc_activation`, so as ferroptotic-
    /// death density rises the suppressive arm grows and the NET immune effect
    /// can flip from pro- to anti-tumor (see
    /// [`crate::immune_spatial::ferroptotic_immunosuppression`]). `0.0` (the
    /// default for `for_2d`/`for_3d`) disables it, byte-identical. The direction
    /// is timing-dependent: a small early-ferroptotic fraction can be
    /// immunogenic (Efimova 2020 PMID 33188036), so this term should dominate
    /// only at sustained/high death density; magnitude is an uncalibrated
    /// placeholder, the sign is the result. 3D-only (consumed by `sim-tme-3d`).
    pub ferro_immunosuppression_strength: f64,
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
            ferro_immunosuppression_strength: 0.0,
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
            ferro_immunosuppression_strength: 0.0,
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
    /// Reversible-to-irreversible (epigenetic) locking rate (#342). Beyond a
    /// threshold of SUSTAINED drug exposure, persistence becomes epigenetically
    /// locked / effectively irreversible: the FSP1/HDAC-mediated suppression of
    /// alternative defenses that makes persistence durable is documented by PMID
    /// 41481741. The sustained-exposure timescale (on the order of days to weeks
    /// of continuous dosing in some reports) is illustrative context only, NOT
    /// from that source and NOT mapped to the dimensionless EMA `lock_threshold`
    /// used here. Per step, once `cumulative_exposure >= lock_threshold`,
    /// a fraction `lock_rate` of the reversible persister pool is moved into the
    /// locked pool, which does NOT revert (see
    /// [`crate::persister::step_with_locking`]). `0.0` (default) disables
    /// locking, so [`step_with_locking`](crate::persister::step_with_locking)
    /// reduces exactly to [`step`](crate::persister::step) and a consumer stays
    /// byte-identical. Uncalibrated placeholder.
    pub lock_rate: f64,
    /// Sustained-exposure threshold that triggers locking (#342). Compared to
    /// the `cumulative_exposure` EMA whose steady state is `avg_drug /
    /// exposure_decay`, so a threshold between the intermittent- and continuous-
    /// dosing steady states makes only CONTINUOUS dosing lock. Inert when
    /// `lock_rate == 0`.
    pub lock_threshold: f64,
    /// Decay of the sustained-exposure tracker per step (#342): the tracker is
    /// `cumulative_exposure = cumulative_exposure·(1 - exposure_decay) + drug`,
    /// an exponential moving average of drug exposure. A positive value makes
    /// drug-off windows decay the tracker, so intermittent dosing never reaches
    /// `lock_threshold` while continuous dosing does. Should be `> 0` when
    /// `lock_rate > 0` (else the tracker grows unbounded and any nonzero average
    /// dose eventually locks). Inert when `lock_rate == 0`.
    pub exposure_decay: f64,
    /// Non-drug stress-niche entry rate (#377). The classic drug-tolerant-
    /// persister biology has a second, NON-DRUG entry route: hypoxic /
    /// nutrient-poor drug-sanctuary microenvironments drive cells into a
    /// slow-cycling, drug-tolerant persister state independent of drug exposure
    /// (hypoxia-induced drug tolerance, e.g. the HIF1α-driven slow-cycling
    /// chemoresistant phenotype, Cuesta-Borràs et al. Cell Rep 2023 PMID 37537841;
    /// quiescent perivascular/hypoxic niches).
    /// This term is DECOUPLED from `drug_intensity`: a consumer applies
    /// [`crate::persister::stress_entry`] with a local stress signal (e.g.
    /// `1 - o2_supply`), which raises the REVERSIBLE pool only (a stress-niche
    /// persister reverts when the niche resolves) and does NOT feed the locking
    /// EMA or the drug-driven resistance. `0.0` (default) ⇒ no stress entry ⇒
    /// byte-identical. Uncalibrated placeholder; direction is the result.
    pub stress_entry_rate: f64,
    /// OXPHOS-ROS suppression at full persistence (#470). A distinct, MITOCHONDRIAL
    /// escape axis: drug-tolerant persisters survive GPX4 inhibition partly by
    /// DOWNREGULATING oxidative phosphorylation, a main source of the
    /// mitochondrial ROS / peroxidizable-substrate flux GPX4 inhibitors act on,
    /// so RSL3 has less to work with and kills persisters less ("FSP1 and histone
    /// deacetylases suppress cancer persister cell ferroptosis", PMID 40909720;
    /// the OXPHOS-suppression / mitochondrial-ROS-supply leg of persister
    /// ferroptosis tolerance, distinct from the existing `gpx4_resistance` and
    /// MUFA axes). A consumer scales a persister cell's basal / mitochondrial ROS
    /// DOWN by [`crate::persister::oxphos_ros_multiplier`], which interpolates
    /// from `1.0` at no persistence to `1 - oxphos_ros_suppression` at full
    /// persistence. `0.0` (default, AND in `enabled()`) ⇒ multiplier `1.0` ⇒
    /// byte-identical (the existing persister snapshot/tests are unaffected).
    /// Uncalibrated placeholder; only the direction (OXPHOS-low persister ⇒ less
    /// RSL3 kill) is claimed.
    pub oxphos_ros_suppression: f64,
    /// HDAC-inhibitor rescue of OXPHOS-ROS (#470), in `[0, 1]`. HDAC inhibitors
    /// re-raise OXPHOS / mitochondrial ROS in persisters and synergize with GPX4
    /// inhibition to kill them (PMID 40909720), so this knob REVERSES
    /// `oxphos_ros_suppression`: the effective suppression a persister keeps is
    /// `oxphos_ros_suppression · (1 - hdac_inhibitor)`, so `hdac_inhibitor = 1`
    /// fully restores the ROS (and the kill) while `0` leaves the suppression
    /// intact. Inert when `oxphos_ros_suppression == 0`. Uncalibrated placeholder.
    pub hdac_inhibitor: f64,
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
            lock_rate: 0.0,
            lock_threshold: 0.0,
            exposure_decay: 0.0,
            stress_entry_rate: 0.0,
            oxphos_ros_suppression: 0.0,
            hdac_inhibitor: 0.0,
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
            // Locking off even in `enabled()`: it is a distinct opt-in (#342),
            // so the existing persister snapshot/tests stay byte-identical.
            lock_rate: 0.0,
            lock_threshold: 0.0,
            exposure_decay: 0.0,
            // Stress-niche entry off even in `enabled()`: a distinct opt-in
            // (#377), so the existing persister snapshot/tests stay byte-identical.
            stress_entry_rate: 0.0,
            // OXPHOS-ROS suppression + HDAC rescue off even in `enabled()`: a
            // distinct opt-in (#470), so the existing persister snapshot/tests
            // stay byte-identical.
            oxphos_ros_suppression: 0.0,
            hdac_inhibitor: 0.0,
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
    ///
    /// Intentionally ignores the locking fields (#342): a locking-only config
    /// (`lock_rate > 0` with all other rates zero) reports `is_identity() ==
    /// true`. That is safe because with acquisition off there is no reversible
    /// pool for locking to ratchet, so the locking path is a no-op anyway; but a
    /// consumer that gates the locking path on `is_identity()` should gate on
    /// `lock_rate` directly instead. `step_with_locking` does NOT use
    /// `is_identity()`; it gates on `lock_rate == 0.0`.
    pub fn is_identity(&self) -> bool {
        self.acquisition_rate == 0.0
            && self.reversion_rate == 0.0
            && self.gpx4_resistance == 0.0
            && self.mufa_boost_per_step == 0.0
            // #377: a stress-entry-only config (all drug rates zero,
            // stress_entry_rate > 0) DOES raise the reversible pool from a stress
            // niche, so it is NOT identity. Include it so a consumer that gates
            // the persister path on is_identity() still runs the stress entry.
            && self.stress_entry_rate == 0.0
            // #470: a config with OXPHOS-ROS suppression set is NOT identity
            // (it scales a persister's basal ROS), mirroring the #377 precedent.
            // hdac_inhibitor alone is inert (it only modulates the suppression),
            // so only oxphos_ros_suppression participates in the contract.
            && self.oxphos_ros_suppression == 0.0
    }
}

// ============================================================
// Biochem parameter overrides (#331)
// ============================================================
//
// One name->field mapping, shared by the Python binding (`sim_batch`/`sim_cell`)
// AND the simulation binaries (sim-tme / sim-combo-mech), so global-sensitivity
// / calibration drivers can perturb the SAME biochemical rate constants on every
// code path without a per-consumer copy of the match arm. An EMPTY override set
// is a no-op, which is what keeps the default (un-driven) runs byte-identical.

/// Apply `(name, value)` biochem-parameter overrides onto `params`, in place.
/// The single source of truth for the override-name -> field mapping. Returns
/// `Err(unknown_name)` on the first unrecognized parameter so the caller can
/// report it. An empty iterator is a no-op.
///
/// The covered set is the biochemical rate constants the Python binding exposes
/// and the PRCC / Sobol sensitivity work (#134/#331) screens; structural toggles
/// and the spatial-layer configs are deliberately NOT here (they are separate
/// `*Config` structs).
pub fn apply_param_overrides<I, S>(params: &mut Params, overrides: I) -> Result<(), String>
where
    I: IntoIterator<Item = (S, f64)>,
    S: AsRef<str>,
{
    for (key, val) in overrides {
        match key.as_ref() {
            "fenton_rate" => params.fenton_rate = val,
            "gsh_scav_efficiency" => params.gsh_scav_efficiency = val,
            "gsh_km" => params.gsh_km = val,
            "nrf2_gsh_rate" => params.nrf2_gsh_rate = val,
            "lp_rate" => params.lp_rate = val,
            "lp_propagation" => params.lp_propagation = val,
            "gpx4_rate" => params.gpx4_rate = val,
            "fsp1_rate" => params.fsp1_rate = val,
            "scd_mufa_rate" => params.scd_mufa_rate = val,
            "scd_mufa_max" => params.scd_mufa_max = val,
            "initial_mufa_protection" => params.initial_mufa_protection = val,
            "scd_mufa_decay" => params.scd_mufa_decay = val,
            "gpx4_degradation_by_ros" => params.gpx4_degradation_by_ros = val,
            "gpx4_nrf2_upregulation" => params.gpx4_nrf2_upregulation = val,
            "sdt_ros" => params.sdt_ros = val,
            "pdt_ros" => params.pdt_ros = val,
            "rsl3_gpx4_inhib" => params.rsl3_gpx4_inhib = val,
            "gsh_max" => params.gsh_max = val,
            "gpx4_nrf2_target_multiplier" => params.gpx4_nrf2_target_multiplier = val,
            "death_threshold" => params.death_threshold = val,
            other => return Err(other.to_string()),
        }
    }
    Ok(())
}

/// Parse a JSON object `{"param_name": value, ...}` into override pairs (sorted
/// by name for determinism; a JSON object cannot carry duplicate keys, and
/// `apply_param_overrides` applies each once, so order does not affect the
/// result). Used by the binaries to read a driver-supplied override set.
pub fn parse_param_overrides_json(s: &str) -> Result<Vec<(String, f64)>, String> {
    let map: std::collections::BTreeMap<String, f64> =
        serde_json::from_str(s).map_err(|e| format!("invalid override JSON: {e}"))?;
    Ok(map.into_iter().collect())
}

/// Pure core of [`param_overrides_from_env`] (no process-env access, so it is
/// unit-testable without env-var races). `None`/blank ⇒ no overrides. A value
/// starting with `{` is treated as an inline JSON object; otherwise it is a path
/// to a JSON file.
pub fn overrides_from_env_value(value: Option<&str>) -> Result<Vec<(String, f64)>, String> {
    match value {
        None => Ok(Vec::new()),
        Some(s) if s.trim().is_empty() => Ok(Vec::new()),
        Some(s) => {
            let t = s.trim();
            let json = if t.starts_with('{') {
                t.to_string()
            } else {
                std::fs::read_to_string(t)
                    .map_err(|e| format!("FERRO_PARAM_OVERRIDES path '{t}': {e}"))?
            };
            parse_param_overrides_json(&json)
        }
    }
}

/// Read biochem overrides from the `FERRO_PARAM_OVERRIDES` env var, which holds
/// EITHER an inline JSON object (`{"lp_rate":0.05}`) OR a path to such a JSON
/// file. Unset/blank ⇒ `Ok(vec![])` (the no-op, byte-identical default path a
/// binary takes when no driver is perturbing it).
pub fn param_overrides_from_env() -> Result<Vec<(String, f64)>, String> {
    overrides_from_env_value(std::env::var("FERRO_PARAM_OVERRIDES").ok().as_deref())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_overrides_are_a_no_op() {
        // The byte-identity guarantee for the default (un-driven) binary path:
        // applying no overrides leaves Params bit-for-bit unchanged.
        let mut p = Params::default();
        let before = serde_json::to_string(&p).unwrap();
        apply_param_overrides(&mut p, Vec::<(String, f64)>::new()).unwrap();
        // Compare the serialized form (Params has no PartialEq; the JSON form is
        // also the surface the binaries' byte-identity is judged on).
        assert_eq!(serde_json::to_string(&p).unwrap(), before);
        // And the env reader returns nothing when the var is unset/blank.
        assert!(overrides_from_env_value(None).unwrap().is_empty());
        assert!(overrides_from_env_value(Some("   ")).unwrap().is_empty());
    }

    #[test]
    fn apply_param_overrides_sets_named_fields_and_rejects_unknown() {
        let mut p = Params::default();
        apply_param_overrides(
            &mut p,
            [
                ("lp_propagation".to_string(), 0.123),
                ("gpx4_rate".to_string(), 0.456),
                ("death_threshold".to_string(), 9.0),
            ],
        )
        .unwrap();
        assert_eq!(p.lp_propagation, 0.123);
        assert_eq!(p.gpx4_rate, 0.456);
        assert_eq!(p.death_threshold, 9.0);
        // Every PRCC-screened name resolves (no typo in the mapping).
        let prcc = [
            "fenton_rate",
            "gsh_scav_efficiency",
            "lp_rate",
            "lp_propagation",
            "gpx4_rate",
            "fsp1_rate",
            "nrf2_gsh_rate",
            "gpx4_degradation_by_ros",
            "death_threshold",
            "sdt_ros",
            "rsl3_gpx4_inhib",
        ];
        for name in prcc {
            apply_param_overrides(&mut Params::default(), [(name, 1.0)]).unwrap();
        }
        // Unknown name is reported, not silently ignored.
        let err =
            apply_param_overrides(&mut Params::default(), [("not_a_param", 1.0)]).unwrap_err();
        assert_eq!(err, "not_a_param");
    }

    #[test]
    fn overrides_from_inline_json_value_round_trip() {
        let pairs = overrides_from_env_value(Some(r#"{"lp_rate": 0.05, "sdt_ros": 7.5}"#)).unwrap();
        // BTreeMap-sorted, so deterministic order: lp_rate before sdt_ros.
        assert_eq!(
            pairs,
            vec![("lp_rate".to_string(), 0.05), ("sdt_ros".to_string(), 7.5)]
        );
        // Applying the parsed pairs sets the fields.
        let mut p = Params::default();
        apply_param_overrides(&mut p, pairs).unwrap();
        assert_eq!(p.lp_rate, 0.05);
        assert_eq!(p.sdt_ros, 7.5);
        // Malformed JSON is an error, not a panic.
        assert!(overrides_from_env_value(Some("{not json")).is_err());
    }

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
