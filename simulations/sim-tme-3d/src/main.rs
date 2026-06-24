//! sim-tme-3d: 3D spheroid tumor microenvironment simulation.
//!
//! Capstone 3D binary for the spheroid-validation series (#185–#197).
//! Integrates all five library primitives landed in v0.7.0–v0.11.0:
//! - 3D energy physics (#186) via `physics::local_ros_multiplier_3d`
//! - 3D radial O₂ gradient (#187) via `oxygen::radial_o2_field`
//! - 3D radial pH gradient (#190) via `ph::radial_ph_field` + helpers
//! - 3D CAF-shielded boundary detection (#189) via `stromal::stromal_adjacency_mask_3d`
//! - 3D spatial DAMP diffusion + activation (#188) via `immune_spatial::*`
//!
//! Produces a per-condition matrix of kill rates that the Python
//! comparison script (`scripts/generate_3d_comparison_table.py`) pairs
//! with sim-tme's existing 2D output to answer the four key questions
//! in issue #195:
//! 1. Does the hypoxia RSL3 collapse hold in 3D?
//! 2. Does the immune 104:1 ratio hold in 3D?
//! 3. Does stromal shielding have MORE impact in 3D?
//! 4. Does pH ion trapping produce similar RSL3 reduction in 3D?
//!
//! ## ⚠️ Scale-mismatch with sim-tme
//!
//! sim-tme uses a 500×500 grid (tumor radius ≈ 4500 µm — large in-vivo
//! tumor). sim-tme-3d uses a **60³ grid** (tumor radius ≈ 540 µm —
//! upper end of in-vitro spheroids) as the default for the 24-condition
//! matrix. Larger single-condition grids are now feasible after the #192
//! perf work: measured 200³ × 180 ≈ 41 s at ~1.29 GB on 10 cores (see
//! `--bench` + the "Performance & scalability" section of the README).
//! Only 500³+ (≈ 18 GB dense) remains out of reach. The 60³ default is a
//! deliberate matrix-throughput choice, not a ceiling.
//!
//! **The Python comparison script reports RATIOS** (e.g., RSL3 hypoxic
//! kill / RSL3 normoxic kill) which are dimensionally meaningful at
//! different scales. Absolute kill counts are NOT directly comparable
//! between the two binaries.
//!
//! ## Why λ sweep skips 150 µm
//!
//! 3D hypoxic-zone threshold = 3λ. At λ=150 µm, threshold = 450 µm,
//! but the 60³ grid's tumor radius is only 540 µm. The hypoxic zone
//! would be ≤1 cell — statistically meaningless. We sweep [80, 100,
//! 120] which give meaningful hypoxic shells at this scale.
//!
//! ## Stability requirement (immune diffusion)
//!
//! Sim-tme's 2D `damp_diffusion_fraction = 0.08` is **unsafe in 3D**
//! (0.08 × 26 = 2.08 > 1, would mass-destroy). We use 0.025 (matches
//! 2D's per-step total diffusion of ~64% — see `immune_spatial` rustdoc).
//! `immune_spatial::diffuse_damp_3d_step` enforces the stability invariant
//! with `assert!` (release-mode panic).

use std::fs;
use std::path::Path;
use std::time::Instant;

mod npy;
mod snapshot;

use ferroptosis_core::acsl4::{pufa_boost_from_status, ACSL4_NEGATIVE};
use ferroptosis_core::alox::AloxConfig;
use ferroptosis_core::biochem::{exo_decay_factor, sim_cell_step, CellState};
use ferroptosis_core::cell::{Phenotype, Treatment};
use ferroptosis_core::clonal::{assign_subclones_3d, repopulate_dead_sites_3d, ClonalConfig};
use ferroptosis_core::contact::{
    apply_contact_resistance_3d, apply_contact_resistance_at_3d, ContactConfig,
};
use ferroptosis_core::copper::CopperConfig;
use ferroptosis_core::dose_schedule::DoseSchedule;
use ferroptosis_core::grid::{TumorGrid3D, TUMOR_RADIUS_FRACTION};
use ferroptosis_core::ifngamma::{acsl4_upregulation, system_xc_retention, IFNGammaConfig};
use ferroptosis_core::immune_spatial::{
    dc_activation, dc_ferroptosis_survival, diffuse_damp_3d_step, exhaustion_factor,
    ferroptotic_immunosuppression, immune_kill_probability, sasp_field_kill_multiplier,
    suppressor_kill_multiplier, suppressor_source_mask_3d, CheckpointPanel, DcFerroptosisConfig,
    DcSubsetConfig, SuppressorConfig, DAMP_KILL_THRESHOLD,
};
use ferroptosis_core::nutrient::{apply_nutrient_stress_3d, NutrientConfig};
use ferroptosis_core::oxygen::{
    fenton_o2_factor, hypoxia_iron_factor, o2_dependent_exo_factor, por_o2_factor, radial_o2_field,
};
use ferroptosis_core::params::{
    Params, PersisterConfig, PhConfig, SpatialImmuneConfig, SpatialParams, StromalConfig,
};
use ferroptosis_core::persister;
use ferroptosis_core::ph::{ion_trap_factor_from_ph, iron_multiplier_from_ph, radial_ph_field};
use ferroptosis_core::phenotype_mufa::{
    apply_phenotype_mufa_3d, apply_phenotype_mufa_at_3d, PhenotypeMufaConfig,
};
use ferroptosis_core::physics::local_ros_multiplier_3d;
use ferroptosis_core::reaction_diffusion::{
    reaction_diffusion_supply_field, ReactionDiffusionConfig,
};
use ferroptosis_core::senescence::{
    apply_senescence_program_3d, sasp_immune_multiplier, SenescenceConfig,
};
use ferroptosis_core::slab::{
    apply_depth_graded_cells_3d, scale_interpretation, slab_supply_field, SlabConfig,
    SlabPhenotypeConfig, KROGH_LAMBDA_UM,
};
use ferroptosis_core::spheroid::{
    apply_radial_cells_3d, apply_radial_cells_sized_3d, radial_fraction_3d, radial_mufa_protection,
    SizeAwareZones, SpheroidConfig,
};
use ferroptosis_core::stromal::stromal_adjacency_mask_3d;
use ferroptosis_core::tumor_pk::RSL3_INACTIVATION_RATE;
use ferroptosis_core::vasculature::{
    hypoxic_fraction, place_vessels_3d, place_vessels_fractal_3d, place_vessels_in_slab_3d,
    vessel_supply_field, VasculatureConfig, VesselTopology,
};
use rand::distributions::Distribution;
use rand::prelude::*;
use rand::rngs::StdRng;
use rand_distr::Normal;
use rayon::prelude::*;
use serde::Serialize;

// ============================================================
// Constants
// ============================================================

const GRID_DIM: usize = 60;
const CELL_SIZE_UM: f64 = 20.0;
const N_STEPS: u32 = 180;
const SEED: u64 = 42;
/// Independent seed for clonal subclone assignment (#242). Distinct from
/// `SEED` so Voronoi seed-point sampling never advances the grid-generation
/// RNG stream — the cell grid stays byte-identical whether or not subclones
/// are assigned.
const SUBCLONE_SEED: u64 = 0x5c10_4e42;
/// Independent seed for vessel placement (#191). Distinct from `SEED` /
/// `SUBCLONE_SEED` so vessel sampling never advances the grid-generation or
/// subclone RNG streams.
const VESSEL_SEED: u64 = 0x7e55_e142;
/// Independent seed for radial spheroid cell re-generation (#197). Distinct
/// from the others so radial re-gen never touches the matrix RNG streams.
const SPHEROID_SEED: u64 = 0x5ade_0142;
/// Independent seed for heuristic Treg/MDSC suppressor-source placement (#264).
/// Distinct from the others so suppressor seeding never touches the matrix or
/// other realism-layer RNG streams. A FIXED constant (not per-condition), so
/// every suppressor-on condition in a run shares the same niche layout — the
/// intent for an A/B (e.g. Treg-present vs Treg-depleted = the same patient,
/// same niches, differing only in whether the field is applied).
const SUPPRESSOR_SEED: u64 = 0x5099_2e64;
/// Independent seed for the depth-graded slab phenotype re-generation (#272).
/// Distinct from the others so slab depth re-gen never touches the matrix RNG
/// streams. Slab and spheroid are mutually-exclusive geometries, so this layer
/// and `SPHEROID_SEED` never coexist in one run; a dedicated constant keeps the
/// invariant explicit and future-proof if that exclusion is ever relaxed.
const SLAB_PHENOTYPE_SEED: u64 = 0x51ab_0142;
/// Independent seed for therapy-induced-senescence cell marking (#341).
/// Distinct from the others so senescence marking never touches the matrix or
/// other realism-layer RNG streams.
const SENESCENCE_SEED: u64 = 0x5e4e_0341;
/// O2-supply factor below which a cell counts as hypoxic for
/// `vascular_hypoxic_fraction` reporting (#191). `exp(-d/λ) < 0.1` is ~2.3 λ
/// from the nearest vessel.
const VASCULAR_HYPOXIC_THRESHOLD: f64 = 0.1;

/// O₂ penetration sweep — 3λ must comfortably fit inside the
/// tumor_radius (60·0.45·20 = 540 µm). Skips λ=150 (3·150=450 →
/// only ~1 cell hypoxic zone, statistically meaningless).
const O2_LAMBDAS: &[f64] = &[80.0, 100.0, 120.0];

/// Zone-analysis reference λ — matches sim-tme.
const ZONE_REF_LAMBDA: f64 = 120.0;

/// DAMP diffusion fraction for 3D, retained as a const for legacy code
/// paths that emit the value into output metadata (see `RunMetadata` /
/// `summary.json`). Identical to `SpatialImmuneConfig::for_3d()
/// .damp_diffusion_fraction`. **Not** the source of truth for the
/// runtime value — that's the library const propagated through
/// `SpatialImmuneConfig`.
const DAMP_DIFFUSION_FRACTION_3D: f64 = 0.025;

/// Step at which immune activity begins (matches sim-tme).
const IMMUNE_START_STEP: u32 = 60;

/// Iron diffusion neighbor fraction for 3D — scaled from sim-tme's 2D
/// default `0.1` (which is per-Moore-8) to the per-Moore-26 equivalent:
/// `0.1 × 8/26 ≈ 0.0308`. Without scaling, 2D's `0.1` would over-spread
/// in 3D (10% × 26 = 260% per source-cell loss, vs 2D's 80%).
const IRON_DIFFUSE_FRACTION_3D: f64 = 0.1 * 8.0 / 26.0;

/// Diffusing SASP field (#376) transport constants, mirroring the Treg/MDSC
/// suppressor field's enabled() defaults (`SuppressorConfig::enabled`): per-step
/// replenishment at each senescent source cell, a 3D-safe diffusion fraction
/// (< 1/26 for stability), and a clearance rate. These set the field's
/// quasi-steady spatial scale; the single tunable knob is the signed
/// `SenescenceConfig::sasp_field_strength` (which absorbs the absolute
/// magnitude), so these are fixed UNCALIBRATED placeholders, not config surface.
const SASP_FIELD_REPLENISH_RATE: f64 = 0.15;
const SASP_FIELD_DIFFUSION_FRACTION: f64 = 0.025;
const SASP_FIELD_CLEARANCE_RATE: f64 = 0.03;

// `DAMP_KILL_THRESHOLD`, `SpatialImmuneConfig`, `StromalConfig`,
// `PhConfig` are imported from `ferroptosis-core` (#220). Previously
// these were binary-local copies duplicated with sim-tme.

// ============================================================
// Output records
// ============================================================

/// Per-condition kill-rate result. Mirrors `sim-tme::ConditionResult`
/// structure so the comparison script can pair fields across the two
/// JSONs.
#[derive(Clone, Debug, Serialize)]
struct ConditionResult {
    treatment: String,
    o2_condition: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    o2_lambda_um: Option<f64>,
    immune_mode: String,
    total_tumor: usize,
    total_dead: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    ferroptosis_kills: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    immune_kills: Option<usize>,
    overall_kill_rate: f64,
    normoxic_kill_rate: f64,
    transition_kill_rate: f64,
    hypoxic_kill_rate: f64,
    /// Peak per-cell DAMP value at end of simulation (max over all cells).
    /// Always populated (the binary computes the mask + DAMP regardless of
    /// immune-on toggle so cross-condition comparison is apples-to-apples).
    /// PR body's "DAMP fields" claim is vindicated by this + total_damp.
    peak_damp: f64,
    /// Sum of DAMP across all cells at end of simulation. Mass-balance
    /// indicator for the diffusion-clearance dynamics.
    total_damp: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    stromal_mode: Option<String>,
    /// Kill rate among cells flagged by `stromal_adjacency_mask_3d` (boundary
    /// cells with ≥1 stromal Moore-26 neighbor). **Always populated** —
    /// the binary computes the mask regardless of `stromal_on` toggle so
    /// Q3 in the comparison script can pair stromal-on adjacent rates
    /// with no-stromal baseline adjacent rates (reviewer-flagged fix).
    /// `None` only when grid produces zero boundary cells.
    #[serde(skip_serializing_if = "Option::is_none")]
    stromal_adjacent_kill_rate: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    stromal_adjacent_count: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    ph_mode: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    ph_edge: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    ph_core: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    ph_lambda_um: Option<f64>,
    /// Mean drug-tolerant persister fraction across surviving tumor cells at
    /// end of simulation (#241). `None` (field omitted) unless the persister
    /// model is enabled, so the default 24-condition matrix summary.json
    /// stays byte-identical to pre-#241 (guarded by #253).
    #[serde(skip_serializing_if = "Option::is_none")]
    persister_mean: Option<f64>,
    /// Mean LOCKED (irreversible) persister fraction across surviving tumor
    /// cells at end of simulation (#342). The epigenetically-locked sub-pool of
    /// `persister_mean` that does NOT revert after drug withdrawal: `> 0` only
    /// under sustained (continuous) exposure with a lock-enabled config, `0`
    /// under intermittent dosing or the default `lock_rate == 0`. `None` (field
    /// omitted) unless the persister model is enabled, so the default matrix
    /// summary.json stays byte-identical.
    #[serde(skip_serializing_if = "Option::is_none")]
    persister_locked_mean: Option<f64>,
    /// Per-subclone kill breakdown (#242). `None` (field omitted) unless the
    /// clonal model is enabled, so the default matrix summary.json stays
    /// byte-identical. One entry per subclone, ordered by id.
    #[serde(skip_serializing_if = "Option::is_none")]
    subclone_kills: Option<Vec<SubcloneKillStat>>,
    /// Hypoxic fraction (tumor cells with O2 supply < 0.1) under the explicit
    /// vessel field (#191). `None` (field omitted) unless vasculature is
    /// enabled; lets the analysis compare irregular-vessel oxygenation against
    /// the smooth edge-distance baseline. Byte-identical default path.
    #[serde(skip_serializing_if = "Option::is_none")]
    vascular_hypoxic_fraction: Option<f64>,
    /// Patient-scale slab interpretation (#240): which depth/scale of a virtual
    /// large tumor this run represents (e.g. "slab spanning depth 4.0–5.2 mm of
    /// a 10 mm virtual tumor (1.2 mm thick)" at the production grid size).
    /// `None` (field omitted) unless slab mode is on — byte-identical default
    /// path.
    #[serde(skip_serializing_if = "Option::is_none")]
    scale_interpretation: Option<String>,
    /// Number of Treg/MDSC suppressor-source (niche) cells (#264 Phase 2).
    /// `None` (field omitted) unless the suppressor field is **active** (config
    /// supplied AND `immune_on`) — byte-identical default path.
    #[serde(skip_serializing_if = "Option::is_none")]
    suppressor_source_count: Option<usize>,
    /// End-of-run suppressor-field max (#264) — representative of the run-time
    /// peak since the field replenishes every step to a quasi-steady state.
    /// `None` (field omitted) unless the suppressor field is active.
    #[serde(skip_serializing_if = "Option::is_none")]
    suppressor_peak: Option<f64>,
    /// Combined multi-checkpoint brake (#264 Phase 3): `1 − Π(1 − residualᵢ)`
    /// over the PD-1/CTLA-4/LAG-3/TIM-3 panel. `None` (field omitted) unless a
    /// checkpoint panel override is supplied — the single-PD-1 path leaves it off.
    #[serde(skip_serializing_if = "Option::is_none")]
    checkpoint_brake: Option<f64>,
    /// Fraction of tumor cells still in the therapy-induced-senescence program
    /// (#341) at end-of-run (clonal repopulation, if on, clears revived sites
    /// from the senescent set). `None` (field omitted) unless the senescence
    /// layer is enabled, so the default path stays byte-identical.
    #[serde(skip_serializing_if = "Option::is_none")]
    senescent_fraction: Option<f64>,
    /// Immune kills among NON-senescent tumor cells (#376). Populated only when
    /// the diffusing SASP field is active (a senescence mask AND `immune_on` AND
    /// non-zero `sasp_field_strength`); `None` (field omitted) otherwise, so the
    /// default matrix summary.json stays byte-identical. Lets a strength>0 vs
    /// strength<0 A/B attribute an immune-kill shift to the BYSTANDER (neighbor)
    /// population, discharging #376's neighbor-coupling acceptance criterion.
    #[serde(skip_serializing_if = "Option::is_none")]
    nonsenescent_immune_kills: Option<usize>,
}

/// Per-subclone kill statistics for one condition (#242). Lets the analysis
/// see how subclonal heterogeneity splits efficacy — the between-subclone
/// kill-rate spread that often exceeds the between-treatment spread.
#[derive(Clone, Debug, Serialize)]
struct SubcloneKillStat {
    subclone_id: u8,
    /// Tumor cells assigned to this subclone at the start (before any clonal
    /// expansion). Equals `total_tumor` when repopulation is off (#266 item 3);
    /// the gap `total_tumor − initial_tumor` is the subclone's net expansion.
    initial_tumor: usize,
    total_tumor: usize,
    /// End-of-run dead count. **With repopulation off this is the cells killed;
    /// with spatial expansion on (#266 item 3) it is the *currently*-dead count**
    /// — a site that died and was repopulated counts as alive here (its death is
    /// still tallied in the run-level `ferroptosis_kills`/`immune_kills`, which
    /// therefore count death *events* and can exceed the population). Read
    /// `total_tumor − initial_tumor` for expansion, not `kill_rate`, under
    /// repopulation.
    total_dead: usize,
    /// `total_dead / total_tumor` — see [`total_dead`](Self::total_dead): this
    /// is "fraction ever killed" only with repopulation off; with expansion on
    /// it is the end-state "fraction currently dead" and reads low.
    kill_rate: f64,
}

/// Schema version for `summary.json`. Bump when the output shape
/// changes; `scripts/generate_3d_comparison_table.py` cross-checks
/// this against `sim-tme/tme_summary.json`'s `schema_version` to catch
/// schema drift across the two binaries (#224 item 2).
const TME_3D_SCHEMA_VERSION: u32 = 1;

#[derive(Clone, Debug, Serialize)]
struct SimulationSummary {
    /// Bumped when the shape of this JSON changes. Currently `1`.
    schema_version: u32,
    grid_dim: usize,
    cell_size_um: f64,
    tumor_radius_um: f64,
    n_steps: u32,
    o2_lambdas: Vec<f64>,
    damp_diffusion_fraction_3d: f64,
    conditions: Vec<ConditionResult>,
    /// Scale-mismatch caveat embedded in the output so consumers see it.
    note: String,
}

// ============================================================
// Condition spec — every entry in the matrix is one of these.
// ============================================================

#[derive(Clone, Debug)]
struct Condition {
    name: String,
    treatment: Treatment,
    treatment_name: String,
    /// `None` = uniform O₂ (baseline); `Some(λ)` = radial O₂ gradient.
    o2_lambda: Option<f64>,
    immune_on: bool,
    stromal_on: bool,
    ph_on: bool,
    /// Drug-administration schedule over time (#239). `Constant` (the
    /// default) reproduces the historical steady-state behavior exactly;
    /// non-constant schedules drive per-step drug modulation.
    dose_schedule: DoseSchedule,
}

// ============================================================
// Zone-kill-rate analog of sim-tme's `zone_kill_rates`. Uses
// TumorGrid3D::radial_depth_um directly — same depth thresholds,
// geometry-agnostic.
// ============================================================

fn zone_kill_rates_3d(grid: &TumorGrid3D, shell_depth_um: f64) -> (f64, f64, f64) {
    let deep_threshold_um = shell_depth_um * 3.0;
    let (mut nd, mut nt, mut td, mut tt, mut hd, mut ht) = (0, 0, 0, 0, 0, 0);
    for r in 0..grid.rows {
        for c in 0..grid.cols {
            for l in 0..grid.layers {
                let gc = grid.get(r, c, l);
                if !gc.is_tumor {
                    continue;
                }
                let depth = grid.radial_depth_um(r, c, l);
                let (dead, total) = if depth < shell_depth_um {
                    (&mut nd, &mut nt)
                } else if depth < deep_threshold_um {
                    (&mut td, &mut tt)
                } else {
                    (&mut hd, &mut ht)
                };
                *total += 1;
                if gc.state.dead {
                    *dead += 1;
                }
            }
        }
    }
    let rate = |d: usize, t: usize| if t > 0 { d as f64 / t as f64 } else { 0.0 };
    (rate(nd, nt), rate(td, tt), rate(hd, ht))
}

// ============================================================
// Per-condition runner — the 3D analog of sim-tme's
// `run_spatial_with_immune`. Uses library primitives throughout.
// ============================================================

fn norm(rng: &mut StdRng, mean: f64, std: f64) -> f64 {
    if std <= 0.0 {
        return mean;
    }
    let dist = Normal::new(mean, std).expect("valid normal");
    dist.sample(rng)
}

/// Grid/step config for a run. Production uses `RunConfig::production()`;
/// tests use `RunConfig::for_test()` to avoid running the full 60³ × 180-step
/// simulation in CI debug builds (reviewer-flagged perf concern).
#[derive(Clone, Copy, Debug)]
struct RunConfig {
    grid_dim: usize,
    n_steps: u32,
}

impl RunConfig {
    fn production() -> Self {
        RunConfig {
            grid_dim: GRID_DIM,
            n_steps: N_STEPS,
        }
    }

    /// Tiny config for unit tests — exercises every code path but ~1000× cheaper.
    /// 10³ = 1000 cells, 20 steps. Runs in <100ms in debug mode.
    #[cfg(test)]
    fn for_test() -> Self {
        RunConfig {
            grid_dim: 10,
            n_steps: 20,
        }
    }
}

fn run_one_condition(condition: &Condition) -> ConditionResult {
    run_one_condition_with_config(condition, RunConfig::production(), None)
}

/// Optional realism-layer overrides for a single condition run. `Default` is
/// all-`None` — the byte-identical matrix path — and each field opts one layer
/// in: persister cells (#241), an immune-config override incl. T-cell
/// exhaustion (#243), and clonal heterogeneity (#242). Bundled into a struct
/// so the run signature stays readable as layers accrue, and so a call site
/// names exactly the layer it enables (`Overrides { clonal: Some(..),
/// ..Default::default() }`).
#[derive(Default)]
struct Overrides {
    persister: Option<PersisterConfig>,
    immune: Option<SpatialImmuneConfig>,
    clonal: Option<ClonalConfig>,
    vasculature: Option<VasculatureConfig>,
    spheroid: Option<SpheroidConfig>,
    /// Size-aware spheroid zone thresholds (#333). `None` ⇒ the fixed
    /// `SpheroidConfig` thresholds are used (the validated large-spheroid limiting
    /// structure), so this is byte-identical to the pre-#333 spheroid path. `Some`
    /// ramps the zone thresholds with the grid's spheroid radius so small spheroids
    /// are correctly mostly-proliferating with no necrotic core. Only meaningful
    /// when `spheroid` is also `Some` (it adjusts that config); ignored otherwise.
    spheroid_size_aware: Option<SizeAwareZones>,
    /// IFN-gamma -> System Xc- ferroptosis-sensitization coupling (#443). `None` /
    /// `disabled()` ⇒ no coupling ⇒ byte-identical. When set AND `immune_on`, an
    /// IFN-gamma field is seeded from the local DAMP (immune-active proxy), diffused,
    /// and used to scale each tumor cell's GSH DOWN (System Xc- downregulation),
    /// sensitizing to ferroptosis (Wang 2019 PMID 31043744). Gated on immune_on (off
    /// the matrix), so the production matrix is byte-identical.
    ifngamma: Option<IFNGammaConfig>,
    /// Copper-ionophore / cuproptosis crosstalk (#485): `Some(cfg)` overloads
    /// intracellular copper, depleting the GSH and GPX4 pools each step (so
    /// RSL3/SDT kills more), with `atp7b_efflux` protecting. `None` / disabled ⇒
    /// retention 1.0 ⇒ byte-identical. Off the production matrix.
    copper: Option<CopperConfig>,
    /// ALOX isoform-specific peroxidation rate + MCFA→ACSL4 PUFA sensitization
    /// (#446). `None` / `identity()` ⇒ both boosts `0` ⇒ byte-identical. When set,
    /// the config's `lp_propagation_boost()` and `mcfa_pufa_boost()` are written
    /// onto `params` before the run (ALOX-high ⇒ faster propagation; MCFA ⇒ more
    /// oxidizable PUFA; both raise ferroptosis). Cell-intrinsic, so not gated on
    /// immune; off the production matrix so the matrix is byte-identical.
    alox: Option<AloxConfig>,
    /// ACSL4-status biomarker (#444): the tumor's relative ACSL4 expression
    /// (`1.0` = wild-type baseline). `None` ⇒ no change ⇒ byte-identical. When set,
    /// `acsl4::pufa_boost_from_status(status)` is written onto
    /// `params.acsl4_status_boost`: `< 1` (ACSL4-low/negative) collapses the
    /// oxidizable-PUFA substrate ⇒ ferroptosis-refractory, `> 1` (ACSL4-high) ⇒
    /// sensitive. Off the production matrix so the matrix is byte-identical.
    acsl4_status: Option<f64>,
    /// ESCRT-III membrane-repair brake on death execution (#465): `Some((rate,
    /// budget))` sets `params.escrt_repair_rate`/`escrt_repair_budget` so a cell
    /// whose LP crosses the death threshold can be resealed for a finite per-cell
    /// budget (more repair ⇒ slower execution ⇒ more resistance). `None` ⇒ both
    /// `0.0` ⇒ the brake never fires ⇒ byte-identical. Off the production matrix.
    escrt: Option<(f64, f64)>,
    /// POR/CYB5R1 enzymatic O2-coupled H2O2 source (#466): `Some((rate, o2_dep))`
    /// injects an enzymatic H2O2 oxidant into each tumor cell's `basal_ros`,
    /// scaled per cell by `por_o2_factor(local_o2, o2_dep)` so POR makes less H2O2
    /// in the hypoxic core (tying the Fenton-feeding oxidant to O2). `None` ⇒ no
    /// injection ⇒ byte-identical. Off the production matrix.
    por: Option<(f64, f64)>,
    /// 7-DHC sterol radical-trapping defense (#467): `Some(pool)` sets
    /// `params.dhc7_radical_trap`, a GPX4-independent radical-trapping quench that
    /// LOWERS the propagation rate ⇒ ferroptosis resistance (a DHCR7-low tumor with
    /// a high 7-DHC pool). `None` ⇒ `0.0` ⇒ byte-identical. Off the production matrix.
    dhc7: Option<f64>,
    /// Vitamin K / VKORC1L1 radical-trapping defense + warfarin (#483):
    /// `Some((trap, warfarin))` sets `params.vitk_radical_trap` and
    /// `params.warfarin_vkor_inhibition`. The VKORC1L1 trap is a sixth
    /// GPX4-independent quench that LOWERS the propagation rate ⇒ ferroptosis
    /// resistance; warfarin (in `[0,1]`) inhibits VKORC1L1, collapsing the trap
    /// and DRIVING ferroptosis. `None` ⇒ both `0.0` ⇒ byte-identical. Off the
    /// production matrix.
    vitk: Option<(f64, f64)>,
    /// PROM2 / MVB-exosome labile-iron efflux (#484): `Some(efflux)` sets
    /// `params.prom2_iron_efflux`, draining the Fenton iron pool over the run
    /// (the OPPOSITE sign to ferritinophagy #340), so a PROM2-high tumor exports
    /// iron under RSL3 and RESISTS. `None` ⇒ `0.0` ⇒ byte-identical. Off the
    /// production matrix.
    prom2: Option<f64>,
    /// Dietary-PUFA supply + lipid-droplet/DGAT buffer (#486): `Some((supply,
    /// buffer))` sets `params.dietary_pufa_supply` and `params.lipid_droplet_buffer`.
    /// Exogenous PUFA above the saturable DGAT buffer raises ferroptosis; DGAT
    /// inhibition (a smaller buffer) makes it emerge sooner. `None` ⇒ both 0.0 ⇒
    /// byte-identical. Off the production matrix.
    dietary_pufa: Option<(f64, f64)>,
    slab: Option<SlabConfig>,
    /// Depth-graded slab phenotype (#272). `None` ⇒ the slab keeps `generate_slab`'s
    /// flat bulk phenotype mix. Only applied when `slab` is also `Some` (it needs the
    /// slab grid + depth offset); ignored otherwise. Off in the matrix ⇒ byte-identical.
    slab_phenotype: Option<SlabPhenotypeConfig>,
    /// Treg/MDSC immunosuppressor field (#264 Phase 2). `None` / disabled ⇒
    /// the suppressor multiplier is identity ⇒ byte-identical.
    suppressor: Option<SuppressorConfig>,
    /// Multi-checkpoint immune brake panel (#264 Phase 3). `None` ⇒ the single
    /// PD-1 `effective_brake` is used (byte-identical); `Some` replaces it with
    /// the combined PD-1/CTLA-4/LAG-3/TIM-3 brake.
    checkpoints: Option<CheckpointPanel>,
    /// Cell-cell contact-mediated ferroptosis resistance (#270). `None` /
    /// identity ⇒ no per-cell lipid/iron modulation ⇒ byte-identical.
    contact: Option<ContactConfig>,
    /// Radial nutrient gradient (#270 item 3b). `None` / identity ⇒ no
    /// antioxidant-setpoint modulation ⇒ byte-identical.
    nutrient: Option<NutrientConfig>,
    /// Phenotype-specific SCD1/MUFA accumulation rates (#363). `None` / identity
    /// ⇒ no per-cell `mufa_rate` is set (cells use the global `scd_mufa_rate`) ⇒
    /// byte-identical. Only meaningful in a MUFA-active context (e.g. the
    /// spheroid preset / an invivo Params); the matrix uses `scd_mufa_rate = 0`,
    /// so even a non-identity config is inert there.
    phenotype_mufa: Option<PhenotypeMufaConfig>,
    /// Dendritic-cell subset mix (#264 Phase 4). `None` / balanced ⇒ priming
    /// efficiency 1.0 ⇒ no immune-kill modulation ⇒ byte-identical.
    dc_subsets: Option<DcSubsetConfig>,
    /// DC / effector-cell ferroptosis susceptibility (#469). `None` / identity
    /// (`susceptibility == 0`) ⇒ DC-survival multiplier 1.0 ⇒ byte-identical.
    /// A ferroptotic TME (high local lipid-ROS / DAMP) kills the effector DCs
    /// themselves, lowering priming, with PD-L1 / checkpoint state protecting
    /// them (PMID 39423128).
    dc_ferroptosis: Option<DcFerroptosisConfig>,
    /// Therapy-induced senescence program (#341). `None` / identity (`fraction
    /// == 0`) ⇒ no senescent cells, no per-cell perturbation, no SASP→immune
    /// coupling ⇒ byte-identical.
    senescence: Option<SenescenceConfig>,
    /// NCOA4-ferritinophagy labile-iron release (#365): the time-dependent
    /// `biochem::ferritinophagy_iron_factor` ramp, consumed inside `sim_cell_step`
    /// via `Params.ferritinophagy_release`. `0.0` (default) ⇒ factor exactly 1.0
    /// every step ⇒ byte-identical.
    ferritinophagy_release: f64,
    /// Hypoxia-driven iron-import sensitivity (#365): `oxygen::hypoxia_iron_factor`
    /// scales each tumor cell's static `iron` UP where local O2 is low (HIF/TfR1),
    /// so the Fenton substrate rises in hypoxic zones even as the O2-dependent SDT
    /// yield (#336) falls. `0.0` (default) ⇒ factor 1.0 ⇒ byte-identical.
    hypoxia_iron_sensitivity: f64,
    /// Oxygen-dependent SDT/PDT exo-ROS yield (#336): the "Type II fraction" of
    /// the exogenous ROS that scales with local O2. `0.0` (default) ⇒ fully
    /// O2-independent (the historical optimistic upper bound) ⇒ the exo-ROS
    /// factor is exactly 1.0 ⇒ byte-identical. `1.0` ⇒ fully Type II /
    /// O2-dependent, so SDT loses efficacy in hypoxic zones like the clinical
    /// SONALA-001 agent (manuscript §7.1).
    sdt_o2_dependence: f64,
    /// O2-dependent Fenton H₂O₂ substrate (#383): `oxygen::fenton_o2_factor`
    /// scales each tumor cell's `iron` (its only consumer is the Fenton term)
    /// DOWN where local O2 is low, since the Fenton reaction needs O2-derived
    /// H₂O₂ (superoxide → SOD → H₂O₂). This is the counterweight to
    /// `hypoxia_iron_sensitivity` (#365): hypoxia raises the iron but lowers the
    /// H₂O₂ substrate, so the NET deep-core Fenton can fall instead of rise,
    /// correcting the §7.1 model artifact where the O2-independent Fenton let
    /// hypoxia-iron "rescue" the anoxic core. `0.0` (default) ⇒ factor exactly
    /// 1.0 ⇒ `iron` unchanged ⇒ byte-identical; `1.0` ⇒ fully O2-gated Fenton.
    fenton_o2_dependence: f64,
    /// Reaction-diffusion supply field (#343 PR 2): when `true` AND explicit
    /// vasculature is on, the per-cell O2/drug supply is solved as the
    /// steady-state reaction-diffusion field (`reaction_diffusion_supply_field`,
    /// vessel Dirichlet sources + tumor consumption + diffusion) instead of the
    /// monotonic `exp(-dist_to_nearest_vessel/λ)` proxy. Same λ, so the two are
    /// apples-to-apples; the RD field has the non-monotonic inter-vessel pockets
    /// the proxy averages away. `false` (default) ⇒ the proxy runs ⇒ the
    /// vasculature path is unchanged; and the whole branch is gated on
    /// vasculature being on (off the matrix), so the production matrix is
    /// byte-identical either way.
    reaction_diffusion: bool,
}

/// Thin wrapper running a condition with NO overrides (all realism layers
/// off). The matrix, dose-sweep, and every test use this 3-arg form, so they
/// stay byte-identical. `run_snapshot`/tests call [`run_one_condition_full`]
/// directly with a populated [`Overrides`] to opt a layer in.
fn run_one_condition_with_config(
    condition: &Condition,
    run_cfg: RunConfig,
    snapshot: Option<&mut snapshot::SnapshotBuffers>,
) -> ConditionResult {
    run_one_condition_full(condition, run_cfg, snapshot, Overrides::default())
}

/// Combine two per-cell supply fields by **element-wise max** (#272 slab +
/// vasculature coupling): each cell draws O2/drug from whichever source — the
/// planar depth gradient (#240) or the nearest internal vessel (#191) — is
/// stronger. Both inputs are in `[0,1]`, so the result is too. Used by both the
/// run path and the `--snapshot` overlay so they stay in lockstep.
fn combine_supply_max(planar: &[f64], vessel: &[f64]) -> Vec<f64> {
    debug_assert_eq!(
        planar.len(),
        vessel.len(),
        "combine_supply_max: field lengths differ ({} vs {})",
        planar.len(),
        vessel.len()
    );
    planar.iter().zip(vessel).map(|(&p, &q)| p.max(q)).collect()
}

/// Per-cell vessel supply factor, either the monotonic nearest-vessel proxy
/// `exp(-dist/λ)` (#191) or — when `reaction_diffusion` is on (#343 PR 2) — the
/// steady-state reaction-diffusion field solved over the same vessels at the
/// same λ. Single source of truth for the run path and the `--snapshot` overlay
/// so they stay in lockstep. `reaction_diffusion == false` (default) ⇒ the proxy
/// ⇒ the vasculature path is unchanged.
fn vessel_or_rd_supply(
    grid: &TumorGrid3D,
    vessels: &[(f64, f64, f64)],
    lambda: f64,
    reaction_diffusion: bool,
) -> Vec<f64> {
    if reaction_diffusion {
        reaction_diffusion_supply_field(grid, vessels, &ReactionDiffusionConfig::new(lambda))
    } else {
        vessel_supply_field(grid, vessels, lambda)
    }
}

fn run_one_condition_full(
    condition: &Condition,
    run_cfg: RunConfig,
    mut snapshot: Option<&mut snapshot::SnapshotBuffers>,
    overrides: Overrides,
) -> ConditionResult {
    // Destructure the optional realism layers. All-`None` (the matrix path)
    // keeps every layer inert → summary.json byte-identical (guarded by #253).
    let persister_cfg = overrides.persister;
    let immune_override = overrides.immune;
    let clonal_cfg = overrides.clonal;
    let vasculature_cfg = overrides.vasculature;
    let spheroid_cfg = overrides.spheroid;
    // Size-aware spheroid zone thresholds (#333): only meaningful with a spheroid.
    let spheroid_size_aware = overrides
        .spheroid_size_aware
        .filter(|_| spheroid_cfg.is_some());
    // IFN-gamma -> System Xc- coupling (#443): only meaningful with immune activity
    // (the field is seeded from DAMP). `disabled()` ⇒ skipped ⇒ byte-identical.
    let ifngamma_cfg = overrides
        .ifngamma
        .filter(|c| !c.is_disabled() && condition.immune_on);
    let ifngamma_ic50 = ifngamma_cfg.map_or(f64::INFINITY, |c| c.system_xc_ic50);
    // Copper-ionophore / cuproptosis crosstalk (#485): a non-identity config
    // depletes GSH/GPX4 each step. `None`/disabled ⇒ retention 1.0 ⇒
    // byte-identical. The retention multipliers are constant per condition, so
    // hoist them out of the per-cell hot loop.
    let copper_cfg = overrides.copper.filter(|c| !c.is_disabled());
    let copper_on = copper_cfg.is_some();
    let copper_gsh_ret = copper_cfg.map_or(1.0, |c| ferroptosis_core::copper::gsh_retention(&c));
    let copper_gpx4_ret = copper_cfg.map_or(1.0, |c| ferroptosis_core::copper::gpx4_retention(&c));
    // ACSL4 arm (#443 follow-up): only active when the strength is > 0 (the GSH /
    // System Xc- arm can run alone). `0.0` ⇒ the boost factor is always 1.0 ⇒ no
    // lipid_unsat change ⇒ byte-identical even when the field is non-zero.
    let ifngamma_acsl4_strength = ifngamma_cfg.map_or(0.0, |c| c.acsl4_strength);
    let ifngamma_acsl4_on = ifngamma_acsl4_strength > 0.0;
    let slab_cfg = overrides.slab;
    // Depth-graded slab phenotype (#272): only meaningful with a slab grid.
    let slab_phenotype_cfg = overrides.slab_phenotype.filter(|_| slab_cfg.is_some());
    let contact_cfg = overrides.contact;
    let nutrient_cfg = overrides.nutrient;
    // Phenotype-specific MUFA rates (#363). Skip identity entirely so the default
    // path sets no per-cell `mufa_rate` ⇒ byte-identical.
    let phenotype_mufa_cfg = overrides.phenotype_mufa.filter(|c| !c.is_identity());
    let dc_subsets_cfg = overrides.dc_subsets;
    // DC ferroptosis susceptibility (#469). `None`/identity (`susceptibility ==
    // 0`) ⇒ DC-survival multiplier 1.0 ⇒ byte-identical.
    let dc_ferroptosis_cfg = overrides.dc_ferroptosis.filter(|c| !c.is_identity());
    // Therapy-induced senescence (#341). `None`/identity ⇒ no cells marked, no
    // perturbation, no SASP coupling ⇒ byte-identical.
    let senescence_cfg = overrides.senescence.filter(|c| !c.is_identity());
    // Dynamic-iron knobs (#365). Both `0.0` (matrix path) ⇒ inert ⇒ byte-identical:
    // `ferritinophagy_release` ⇒ `ferritinophagy_iron_factor` returns 1.0 every
    // step; `hypoxia_iron_sensitivity` ⇒ `hypoxia_iron_factor` returns 1.0.
    let ferritinophagy_release = overrides.ferritinophagy_release;
    let hypoxia_iron_sensitivity = overrides.hypoxia_iron_sensitivity;
    // Oxygen-dependent SDT/PDT exo-ROS (#336). `0.0` (default/matrix) ⇒ the
    // exo-ROS O2 factor is exactly 1.0 ⇒ byte-identical.
    let sdt_o2_dependence = overrides.sdt_o2_dependence;
    // Oxygen-dependent Fenton H2O2 substrate (#383). `0.0` (default/matrix) ⇒
    // `fenton_o2_factor` returns 1.0 ⇒ `cell.iron` unchanged ⇒ byte-identical.
    let fenton_o2_dependence = overrides.fenton_o2_dependence;
    // Reaction-diffusion supply (#343 PR 2). `false` (default/matrix) ⇒ the
    // monotonic `vessel_supply_field` proxy runs; only meaningful with explicit
    // vasculature on (off the matrix), so byte-identical either way.
    let reaction_diffusion = overrides.reaction_diffusion;
    // Treg/MDSC suppressor (#264 Phase 2). `None`/disabled ⇒ no source mask,
    // no field, identity multiplier ⇒ byte-identical.
    let suppressor_cfg = overrides.suppressor.filter(|c| !c.is_disabled());
    // Multi-checkpoint brake (#264 Phase 3). `None` ⇒ the single PD-1
    // `effective_brake` (byte-identical); `Some` ⇒ the combined panel brake.
    let checkpoint_panel = overrides.checkpoints;
    // Slab and spheroid are mutually-exclusive geometries (#240 review): slab
    // builds an all-tumor block and slab supply wins, but spheroid would run
    // `Params::spheroid()` + radial phenotype re-assignment keyed on a spheroid
    // center that an all-tumor cube does not have — an incoherent mix. No
    // preset/path sets both; guard it so a future caller can't combine them by
    // accident. (Slab + vasculature, by contrast, COMPOSE since #272: the planar
    // depth gradient and the vessel proximity field combine by element-wise max
    // — see `supply_field` below.)
    debug_assert!(
        !(slab_cfg.is_some() && spheroid_cfg.is_some()),
        "slab and spheroid overrides are mutually exclusive (incompatible geometries)"
    );
    // Contact resistance (#270) assumes a centred sphere/spheroid: its fixed-26
    // contact denominator treats domain-boundary cells as low-contact (true for
    // a spheroid, whose tumor never touches the box face; WRONG for an all-tumor
    // slab, whose outer shell IS tumor at the box face and would be mis-scored as
    // a "surface"). Guard the combination rather than silently mis-score a slab.
    debug_assert!(
        !(slab_cfg.is_some() && contact_cfg.is_some()),
        "slab and contact overrides are mutually exclusive (the fixed-26 contact \
         denominator mis-scores a slab's domain-boundary shell as low-contact)"
    );
    // 3D spheroid context (#197): partially-active MUFA so position-dependent
    // per-cell MUFA persists. `None` (matrix path) ⇒ `Params::default()` ⇒
    // byte-identical.
    let mut params = if spheroid_cfg.is_some() {
        Params::spheroid()
    } else {
        Params::default()
    };
    // #365: NCOA4-ferritinophagy labile-iron release ramp. `0.0` (matrix) ⇒
    // `ferritinophagy_iron_factor` is exactly 1.0 for every step ⇒ byte-identical;
    // `ferritinophagy_tau` keeps its default. Consumed inside `sim_cell_step`.
    params.ferritinophagy_release = ferritinophagy_release;
    // ALOX isoform-specific propagation + MCFA→ACSL4 PUFA sensitization (#446).
    // `None` / `identity()` ⇒ both boosts `0.0` (the Params defaults) ⇒
    // byte-identical. Cell-intrinsic, applied uniformly per condition (per-cell
    // stochastic ALOX heterogeneity is a deferred refinement).
    if let Some(cfg) = overrides.alox.filter(|c| !c.is_identity()) {
        params.alox_propagation_boost = cfg.lp_propagation_boost();
        params.mcfa_pufa_boost = cfg.mcfa_pufa_boost();
    }
    // ACSL4-status biomarker stratification (#444): the tumor's ACSL4 expression
    // status sets its oxidizable-PUFA baseline. `None` (matrix path) ⇒ boost stays
    // 0.0 ⇒ byte-identical. Cell-intrinsic, applied uniformly per condition.
    if let Some(status) = overrides.acsl4_status {
        params.acsl4_status_boost = pufa_boost_from_status(status);
    }
    // ESCRT-III membrane-repair brake (#465): `None` (matrix path) ⇒ both 0.0 ⇒
    // byte-identical. When set, a cell crossing the death threshold can be resealed
    // for a finite budget, delaying execution (more repair ⇒ more resistance).
    if let Some((rate, budget)) = overrides.escrt {
        params.escrt_repair_rate = rate;
        params.escrt_repair_budget = budget;
    }
    // 7-DHC sterol radical-trapping defense (#467): `None` (matrix path) ⇒ 0.0 ⇒
    // byte-identical. A high pool (DHCR7-low tumor) raises the GPX4-independent
    // quench and lowers propagation ⇒ ferroptosis resistance.
    if let Some(pool) = overrides.dhc7 {
        params.dhc7_radical_trap = pool;
    }
    // Vitamin K / VKORC1L1 radical-trapping defense + warfarin (#483): `None`
    // (matrix path) ⇒ both 0.0 ⇒ byte-identical. A high trap (a VKORC1L1-high,
    // p53-competent tumor) resists ferroptosis; warfarin inhibition collapses it.
    if let Some((trap, warfarin)) = overrides.vitk {
        params.vitk_radical_trap = trap;
        params.warfarin_vkor_inhibition = warfarin;
    }
    // PROM2 / MVB-exosome labile-iron efflux (#484): `None` (matrix path) ⇒ 0.0 ⇒
    // byte-identical. A high efflux (a PROM2-high tumor) drains the Fenton iron
    // pool over the run, so RSL3 kills LESS (ferroptosis resistance).
    if let Some(efflux) = overrides.prom2 {
        params.prom2_iron_efflux = efflux;
    }
    // Dietary-PUFA / DGAT buffer (#486): `None` (matrix path) ⇒ both 0.0 ⇒
    // byte-identical. Exogenous PUFA above the saturable lipid-droplet buffer
    // adds oxidizable substrate (more ferroptosis); DGAT inhibition lowers the
    // buffer so it emerges sooner.
    if let Some((supply, buffer)) = overrides.dietary_pufa {
        params.dietary_pufa_supply = supply;
        params.lipid_droplet_buffer = buffer;
    }
    let spatial_params = SpatialParams {
        cell_size_um: CELL_SIZE_UM,
        // SpatialParams::default()'s `neighbor_iron_fraction = 0.1` is
        // 2D-tuned (per-Moore-8). Override here with the 3D-scaled value
        // so the field is correct in the struct, not just at the inline
        // call site (reviewer-flagged future-refactor footgun: if a
        // future `local_ros_multiplier_3d` ever reads this field, the
        // 2D value would silently propagate).
        neighbor_iron_fraction: IRON_DIFFUSE_FRACTION_3D,
        ..Default::default()
    };

    // Per-condition deterministic seed for the *runtime RNG* (per-cell ROS
    // noise, immune-kill rolls, etc.). NOT used for grid generation — sim-tme
    // uses the fixed `SEED` so every condition operates on the SAME tumor
    // geometry (same cell phenotypes, same Persister cluster placements).
    // Using a per-condition seed for the grid would mean each treatment is
    // evaluated on a DIFFERENT tumor — combining treatment effects with
    // grid-geometry noise. At ~82.5k tumor cells (60³ sphere × 0.45 radius
    // fraction), with ~268 cells/cluster × 10-20 random clusters, that
    // noise is material relative to subtle effects like the immune ratio.
    let cond_seed = SEED.wrapping_add(hash_condition_name(&condition.name));

    // Grid uses fixed SEED — matches sim-tme/src/main.rs (line 921 and
    // every other TumorGrid::generate call uses bare `SEED`). Same tumor
    // geometry across all conditions ⇒ valid treatment-effect comparison.
    // Slab mode (#240) uses an all-tumor block instead of the spheroid; gated
    // ⇒ the default matrix uses `generate` ⇒ byte-identical.
    let mut grid = if slab_cfg.is_some() {
        TumorGrid3D::generate_slab(
            run_cfg.grid_dim,
            run_cfg.grid_dim,
            run_cfg.grid_dim,
            CELL_SIZE_UM,
            SEED,
        )
    } else {
        TumorGrid3D::generate(
            run_cfg.grid_dim,
            run_cfg.grid_dim,
            run_cfg.grid_dim,
            CELL_SIZE_UM,
            SEED,
        )
    };
    let n_cells = grid.cells.len();

    // 3D spheroid radial cell biology (#197): re-assign tumor cells radially
    // (glycolytic rim / OXPHOS mid / persister core) + GSH(core-low) /
    // iron(core-high) gradients, via an INDEPENDENT RNG so generate's stream is
    // untouched. `None` (matrix path) ⇒ the random grid is left as-is ⇒
    // byte-identical. Runs FIRST (it rewrites cells); later layers (O2/pH/
    // subclone/vessel) then perturb the radially-assigned cells. Position-
    // dependent MUFA is applied after CellState init below.
    if let Some(cfg) = &spheroid_cfg {
        // #333: size-aware thresholds (ramp the zones with the spheroid radius) when
        // opted in, else the fixed limiting-structure thresholds (byte-identical).
        if let Some(size_aware) = &spheroid_size_aware {
            apply_radial_cells_sized_3d(&mut grid, cfg, size_aware, CELL_SIZE_UM, SPHEROID_SEED);
        } else {
            apply_radial_cells_3d(&mut grid, cfg, SPHEROID_SEED);
        }
    }

    // Depth-graded slab phenotype (#272): re-assign the slab's flat bulk mix so
    // vessel-proximal +z layers are proliferating/glycolytic and chronically
    // supply-starved deep (−z) layers are persister-like, the depth-axis analog
    // of the spheroid's rim-to-core structure. Uses the SAME lambda and depth
    // offset the supply field uses below, so phenotype and O2/drug gradients are
    // coherent. Independent per-cell RNG, so generate_slab's stream is untouched;
    // gated on slab mode plus an explicit `slab_phenotype` override, so the
    // matrix (no slab) stays byte-identical. Runs after spheroid (mutually
    // exclusive geometries, guarded above) and before the vessel/subclone
    // layers, which then perturb the depth-assigned cells.
    if let (Some(slab), Some(pheno)) = (&slab_cfg, &slab_phenotype_cfg) {
        let lambda = condition.o2_lambda.unwrap_or(KROGH_LAMBDA_UM);
        apply_depth_graded_cells_3d(
            &mut grid,
            slab.depth_offset_mm * 1000.0,
            lambda,
            pheno,
            SLAB_PHENOTYPE_SEED,
        );
    }

    // Clonal heterogeneity (#242): assign each cell to a subclone via Voronoi,
    // using an INDEPENDENT seed so generate's RNG stream (and the cell grid)
    // is unchanged. `None` on the matrix path → no assignment, no perturbation
    // → byte-identical.
    // `mut` because spatial clonal expansion (#266 item 3) rewrites a revived
    // dead site's subclone id over the run; static (repopulation off) leaves it
    // untouched.
    let mut subclone_ids: Option<Vec<u8>> = clonal_cfg
        .as_ref()
        .map(|c| assign_subclones_3d(&grid, c.k(), SUBCLONE_SEED));
    // Initial per-subclone tumor-cell census (at assignment, all living), so
    // the summary can show the composition shift under expansion. `None` unless
    // clonal is on. Index 0 = stroma (unused); 1..=k are the subclones.
    let initial_subclone_totals: Option<Vec<usize>> = subclone_ids
        .as_ref()
        .zip(clonal_cfg.as_ref())
        .map(|(ids, c)| {
            let mut totals = vec![0usize; c.k() + 1];
            for (idx, gc) in grid.cells.iter().enumerate() {
                if gc.is_tumor {
                    totals[ids[idx] as usize] += 1;
                }
            }
            totals
        });

    // Explicit vasculature (#191): place internal vessel seed points (own seed
    // → generate's RNG untouched), then build the per-cell supply factor
    // `exp(-dist_to_nearest_vessel / λ)`. This REPLACES the edge-distance O2
    // proxy and also supplies drug. `None` on the matrix path (no vasculature
    // or no λ) → the edge-distance `radial_o2_field` path runs → byte-identical.
    let vessels: Option<Vec<(f64, f64, f64)>> = vasculature_cfg.map(|cfg| {
        if slab_cfg.is_some() {
            // Slab (#240) is an all-tumor block, not a central sphere — scatter
            // vessels across the whole block (#272 coupling) so deep tissue
            // throughout, not just a central pocket, gets perfused.
            place_vessels_in_slab_3d(&grid, &cfg, VESSEL_SEED)
        } else {
            // Sphere/spheroid grid: random points (#191) or a fractal-branching
            // tree (#268), per the config's topology.
            match cfg.topology {
                VesselTopology::Fractal => place_vessels_fractal_3d(&grid, &cfg, VESSEL_SEED),
                VesselTopology::Random => place_vessels_3d(&grid, &cfg, VESSEL_SEED),
            }
        }
    });
    // Treg/MDSC suppressor source mask (#264 Phase 2): perivascular niches when
    // vessels are present, else heuristic patches (independent SUPPRESSOR_SEED).
    // `None` ⇒ suppressor off ⇒ byte-identical.
    let suppressor_sources: Option<Vec<bool>> = suppressor_cfg
        .map(|cfg| suppressor_source_mask_3d(&grid, &cfg, vessels.as_deref(), SUPPRESSOR_SEED));
    // Unified per-cell O2/drug **supply field**: the slab depth gradient (#240),
    // the vessel proximity field (#191), or — when both are set — their
    // element-wise combination (#272 coupling). All replace the edge-distance O2
    // proxy and scale drug delivery downstream. `None` on the matrix path ⇒ the
    // edge-distance `radial_o2_field` path runs ⇒ byte-identical.
    let supply_field: Option<Vec<f64>> = match (&slab_cfg, vessels.as_ref()) {
        (Some(cfg), Some(v)) => {
            // Slab + internal vessels (#272): combine the planar depth gradient
            // (#240) with per-vessel proximity (#191) by element-wise MAX — each
            // cell draws O2 from whichever source is stronger (the +z face or
            // the nearest internal vessel), so a deep slab cell next to a vessel
            // gets a focal well-perfused pocket instead of monotonic depth
            // collapse. Bounded in [0,1] (both factors are). λ as below.
            let lambda = condition.o2_lambda.unwrap_or(KROGH_LAMBDA_UM);
            let planar = slab_supply_field(&grid, cfg.depth_offset_mm * 1000.0, lambda);
            let vessel = vessel_or_rd_supply(&grid, v, lambda, reaction_diffusion);
            Some(combine_supply_max(&planar, &vessel))
        }
        (Some(cfg), None) => {
            // Slab only: planar depth gradient from the +z vessel face. Uses the
            // condition's λ when set (e.g. the `--snapshot=slab` preset's
            // ZONE_REF_LAMBDA = 120 µm), and only falls back to the Krogh
            // default (150 µm) for a condition that left λ unset. The gradient
            // is intrinsic to depth, independent of the radial-O2 condition flag.
            let lambda = condition.o2_lambda.unwrap_or(KROGH_LAMBDA_UM);
            Some(slab_supply_field(
                &grid,
                cfg.depth_offset_mm * 1000.0,
                lambda,
            ))
        }
        (None, _) => match (vessels.as_ref(), condition.o2_lambda) {
            (Some(v), Some(lambda)) => {
                Some(vessel_or_rd_supply(&grid, v, lambda, reaction_diffusion))
            }
            _ => None,
        },
    };

    // --- Apply O₂ gradient (mutates cell.basal_ros) ---
    // Per-cell O2 supply factor captured for the O2-dependent exo-ROS scaling
    // (#336): whichever O2 representation scales basal_ros (the supply field or
    // the edge-distance radial O2) is the same availability that gates a Type II
    // sonosensitizer's ROS yield. 1.0 (full O2) for stroma and when no O2 model
    // is active. Read ONLY through `o2_dependent_exo_factor`, which returns 1.0
    // when `sdt_o2_dependence == 0` ⇒ the matrix path is byte-identical.
    let mut o2_supply_for_exo: Vec<f64> = vec![1.0; n_cells];
    if let Some(supply) = &supply_field {
        // Slab depth gradient (#240) or explicit vessels (#191). Applies
        // regardless of the radial-O2 condition flag (the supply IS the O2
        // model in these modes).
        for (idx, &factor) in supply.iter().enumerate() {
            if grid.cells[idx].is_tumor {
                grid.cells[idx].cell.basal_ros *= factor;
                o2_supply_for_exo[idx] = factor;
            }
        }
    } else if let Some(lambda) = condition.o2_lambda {
        // Edge-distance proxy (default).
        let o2_factors = radial_o2_field(&grid, lambda);
        for (idx, &factor) in o2_factors.iter().enumerate() {
            if grid.cells[idx].is_tumor {
                grid.cells[idx].cell.basal_ros *= factor;
                o2_supply_for_exo[idx] = factor;
            }
        }
    }

    // --- Hypoxia-driven iron import (#365) ---
    // HIF/TfR1 raise the labile-iron pool where O2 is low, so scale each tumor
    // cell's static `iron` UP by `hypoxia_iron_factor(local_o2, sensitivity)`
    // using the same per-cell O2 supply that gates the SDT exo-ROS (#336). This
    // is the spatial counterweight to the O2-dependent SDT collapse: hypoxia can
    // RAISE the Fenton iron substrate even as it lowers the Type II exo-ROS yield,
    // qualifying the §7.1 "hypoxia uniformly protects RSL3" framing. `sensitivity
    // == 0` (matrix path) ⇒ factor 1.0 ⇒ `cell.iron` unchanged ⇒ byte-identical.
    // Composes multiplicatively with the pH iron-release below (order immaterial).
    if hypoxia_iron_sensitivity > 0.0 {
        for (idx, gc) in grid.cells.iter_mut().enumerate() {
            if gc.is_tumor {
                gc.cell.iron *=
                    hypoxia_iron_factor(o2_supply_for_exo[idx], hypoxia_iron_sensitivity);
            }
        }
    }

    // --- O2-dependent Fenton H2O2 substrate (#383) ---
    // The Fenton reaction (Fe2+ + H2O2 -> Fe3+ + OH·) needs H2O2, whose source
    // is O2 (superoxide -> SOD -> H2O2, the Haber-Weiss chain). The model's
    // Fenton term is O2-INDEPENDENT, so the #365 hypoxia-iron rise lets the model
    // wrongly "rescue" the anoxic core (§7.1 artifact). Couple the effective
    // Fenton flux to local O2 by scaling `cell.iron` (its ONLY consumer is the
    // Fenton term) DOWN where O2 is low. Composes multiplicatively with the #365
    // hypoxia-iron scaling above (order immaterial): hypoxia raises the iron but
    // gates the H2O2 substrate, so the NET deep-core Fenton can fall. `dependence
    // == 0` (matrix path) ⇒ factor 1.0 ⇒ `cell.iron` unchanged ⇒ byte-identical.
    // NOTE: this gates the STATIC `cell.iron` pool once at setup; the per-step
    // neighbor-death-diffused pool (`extra_iron`) is gated by the SAME factor
    // inside the per-step loop below, so the WHOLE Fenton substrate is O2-coupled.
    if fenton_o2_dependence > 0.0 {
        for (idx, gc) in grid.cells.iter_mut().enumerate() {
            if gc.is_tumor {
                gc.cell.iron *= fenton_o2_factor(o2_supply_for_exo[idx], fenton_o2_dependence);
            }
        }
    }

    // --- POR/CYB5R1 enzymatic O2-coupled H2O2 source (#466) ---
    // POR/CYB5R1 transfer electrons from NAD(P)H to O2 to make H2O2 (the Fenton
    // substrate), so inject an enzymatic oxidant into each tumor cell's `basal_ros`
    // (which feeds `total_ros`), scaled per cell by `por_o2_factor(local_o2, dep)`
    // using the same per-cell O2 supply as the SDT exo-ROS / hypoxia-iron legs. The
    // O2-coupling means POR makes less H2O2 in the hypoxic core, so unlike a uniform
    // oxidant it does NOT amplify the deep-core artifact. `None` (matrix path) ⇒ no
    // injection ⇒ byte-identical.
    if let Some((por_rate, por_o2_dep)) = overrides.por {
        for (idx, gc) in grid.cells.iter_mut().enumerate() {
            if gc.is_tumor {
                gc.cell.basal_ros += por_rate * por_o2_factor(o2_supply_for_exo[idx], por_o2_dep);
            }
        }
    }

    // --- Apply pH gradient if requested (mutates cell.iron via library helper) ---
    let ph_field = if condition.ph_on {
        let cfg = PhConfig::default();
        let field = radial_ph_field(&grid, cfg.ph_edge, cfg.ph_core, cfg.lambda_ph_um);
        for (idx, &local_ph) in field.iter().enumerate() {
            if grid.cells[idx].is_tumor {
                let iron_mult =
                    iron_multiplier_from_ph(local_ph, cfg.ph_edge, cfg.iron_ph_sensitivity);
                grid.cells[idx].cell.iron *= iron_mult;
            }
        }
        Some((field, cfg))
    } else {
        None
    };

    // --- Compute stromal adjacency mask ALWAYS ---
    // The mask is grid-geometry-dependent (not stromal_on-toggle-dependent),
    // so it's identical across all conditions on the same grid. Computing it
    // for every condition (cheap: O(N·26) one-time) lets the comparison
    // script's Q3 pair stromal-on adjacent rates against the matching
    // no-stromal baseline adjacent rates — measuring CAF shielding at the
    // boundary directly, not as a diluted whole-tumor delta.
    let adjacency_mask = stromal_adjacency_mask_3d(&grid);

    // Apply CAF GSH/MUFA boosts only when stromal_on is true.
    let stromal_cfg = if condition.stromal_on {
        Some(StromalConfig::default())
    } else {
        None
    };

    // --- Treatment-specific ROS multiplier (3D version) ---
    // PDT arm is intentionally kept for sim-tme parity (sim-tme has the same
    // treatment-ROS-multiplier pattern); `generate_conditions()` does not include PDT
    // in the v1 matrix (deferred to follow-up — manuscript focuses on the
    // RSL3-vs-SDT comparison, with PDT in a separate sim-spatial track).
    // The arm is dead code in this binary today but lets a future caller
    // pass `Treatment::PDT` without a match-exhaustiveness error.
    let base_ros = match condition.treatment {
        Treatment::SDT => params.sdt_ros,
        Treatment::PDT => params.pdt_ros,
        Treatment::RSL3 | Treatment::Control => 0.0,
    };

    // Time-varying dose schedule (#239). `dosed == false` for the default
    // `Constant` schedule, in which case EVERY dose-related branch below is
    // skipped and the run is byte-identical to the pre-#239 behavior.
    let schedule = &condition.dose_schedule;
    let dosed = !schedule.is_constant();
    // `base_exo` is read ONLY by the SDT/PDT per-step rescale arm. Allocate
    // (and init-write) it only when that arm can run — i.e. a non-constant
    // schedule on an exogenous-ROS treatment. Empty (never indexed) on the
    // Constant path AND on the RSL3 dosed path, which keeps the
    // "empty ⇒ never read" invariant tight (review #7).
    let dose_modulates_exo =
        dosed && matches!(condition.treatment, Treatment::SDT | Treatment::PDT);
    let mut base_exo: Vec<f64> = if dose_modulates_exo {
        vec![0.0; n_cells]
    } else {
        Vec::new()
    };

    // Initialize cell states with per-cell ROS peak
    for r in 0..grid.rows {
        for c in 0..grid.cols {
            for l in 0..grid.layers {
                let idx = grid.flat_index(r, c, l);
                let exo_ros_peak = {
                    let raw = if matches!(condition.treatment, Treatment::Control | Treatment::RSL3)
                    {
                        0.0
                    } else {
                        let depth_um = grid.radial_depth_um(r, c, l);
                        let ros_multiplier =
                            local_ros_multiplier_3d(depth_um, condition.treatment, &spatial_params);
                        let mut rng = StdRng::seed_from_u64(cond_seed.wrapping_add(idx as u64));
                        let peak = base_ros * ros_multiplier;
                        norm(&mut rng, peak, peak * 0.2).max(0.0)
                    };
                    // Vessel-delivered sonosensitizer/photosensitizer (#191):
                    // the exo-ROS dose scales by vessel proximity on EVERY path
                    // (constant + dosed). ×1.0 when off → byte-identical.
                    // Oxygen-dependent ROS yield (#336): a Type II sonosensitizer
                    // generates O2-dependent singlet oxygen, so the exo-ROS is
                    // scaled by local O2 availability. ×1.0 when
                    // sdt_o2_dependence == 0 → byte-identical.
                    raw * supply_field.as_ref().map_or(1.0, |s| s[idx])
                        * o2_dependent_exo_factor(o2_supply_for_exo[idx], sdt_o2_dependence)
                };
                if dose_modulates_exo {
                    base_exo[idx] = exo_ros_peak;
                }
                let gc = &mut grid.cells[idx];
                // For RSL3 under a non-constant schedule, skip the one-shot
                // init GPX4 knockdown — the schedule drives per-step covalent
                // inactivation instead (the `tumor_pk::sim_cell_with_pk`
                // model). Every other case keeps the historical one-shot
                // init via `from_cell_with_ros` (= `..._opts(.., true)`),
                // preserving byte-identical output on the default path.
                gc.state = if dosed && condition.treatment == Treatment::RSL3 {
                    CellState::from_cell_with_ros_opts(
                        &gc.cell,
                        condition.treatment,
                        &params,
                        exo_ros_peak,
                        false,
                    )
                } else {
                    CellState::from_cell_with_ros(
                        &gc.cell,
                        condition.treatment,
                        &params,
                        exo_ros_peak,
                    )
                };
                gc.extra_iron = 0.0;
                gc.newly_dead = false;
                // Init to NaN so any code path that reads `lp_at_grace_end`
                // before writing it (grace-end write or end-of-sim
                // catch-all) produces NaN downstream — calibration trips
                // instead of silently using a stale value (#225 review).
                gc.lp_at_grace_end = f64::NAN;
            }
        }
    }

    // Position-dependent MUFA (#197): peripheral cells reach higher membrane
    // MUFA than the nutrient-deprived core. Set on the freshly-initialized
    // state. Gated on spheroid ⇒ the default path never runs this ⇒
    // byte-identical.
    //
    // Durability: #270 made this axis durable. `apply_radial_cells_3d` (above)
    // sets a per-cell `cell.mufa_cap` (rim-high, core-low), and
    // `update_mufa_protection` saturates each cell toward ITS OWN cap rather
    // than the uniform `Params::spheroid` M_ss, so the rim-vs-core MUFA spread
    // persists for the whole run instead of relaxing away. This block seeds the
    // matching INITIAL `state.mufa_protection`; the durable carrying capacity is
    // the per-cell cap (`scd_mufa_max` = 0.25 is only the fallback when no
    // per-cell cap is set). `iron` (static `cell.iron`) is likewise fully
    // durable; `gsh` remains an initial condition (evolves under NRF2 resynthesis).
    if let Some(cfg) = &spheroid_cfg {
        for idx in 0..grid.cells.len() {
            if grid.cells[idx].is_tumor {
                let frac = radial_fraction_3d(&grid, idx);
                grid.cells[idx].state.mufa_protection = radial_mufa_protection(frac, cfg);
            }
        }
    }

    // Per-subclone parameter perturbations (#242), gated on clonal. RNG-neutral
    // one-time setup mutations (like the O2/pH gradients). Identity
    // perturbations (k=1 / single_identity) are no-ops → byte-identical.
    //
    // Durability of each axis (all read every step in sim_cell_step):
    //  - iron_mul → cell.iron (static): fully durable.
    //  - lipid_unsat_mul → cell.lipid_unsat (static): fully durable. This is
    //    the MUFA axis; perturbing state.mufa_protection instead would be
    //    overwritten on step 1 by update_mufa_protection under default params
    //    (scd_mufa_max == 0), so it must scale the static PUFA field (#265 rev).
    //  - gpx4_mul → state.gpx4 (initial) AND cell.nrf2 (static setpoint): now
    //    DURABLE (#266). Scaling only the initial state.gpx4 shaped the early
    //    autocatalytic window but relaxed back toward the shared NRF2 setpoint
    //    (gpx4_target = nrf2·gpx4_nrf2_target_multiplier, ~0.008/step), so a
    //    "GPX4-low" identity was transient over 180 steps. Also scaling the
    //    static cell.nrf2 setpoint keeps the axis differentiated for the whole
    //    run. cell.nrf2 also drives GSH resynthesis, so this axis is "general
    //    antioxidant capacity" (deliberate; NRF2 is the master regulator).
    //    Composes multiplicatively with the pH ion-trap correction below.
    if let (Some(ids), Some(cfg)) = (&subclone_ids, &clonal_cfg) {
        for (idx, gc) in grid.cells.iter_mut().enumerate() {
            if !gc.is_tumor {
                continue;
            }
            // ids[idx] is 1..=k for tumor cells (0 only for stroma, skipped).
            // `apply` scales iron + lipid_unsat (static), the initial gpx4
            // reserve, AND the durable cell.nrf2 setpoint (#266); unit-tested
            // in ferroptosis-core::clonal.
            let p = &cfg.perturbations[(ids[idx] - 1) as usize];
            p.apply(&mut gc.cell, &mut gc.state);
        }
    }

    // Cell-cell contact resistance (#270): dense / highly-contacting tumor cells
    // resist ferroptosis (E-cadherin → Merlin/NF2 → YAP inhibition → ACSL4/TFRC
    // down; Wu 2019, PMID 31341276). Scales the durable lipid_unsat (PUFA) and
    // iron axes down with each cell's tumor-neighbour fraction. Geometric (no
    // RNG), off-by-default identity ⇒ byte-identical. Runs after clonal so the
    // two multiplicative per-cell axes compose.
    if let Some(cfg) = &contact_cfg {
        apply_contact_resistance_3d(&mut grid, cfg);
    }

    // Radial nutrient gradient (#270 item 3b): the nutrient-starved core has
    // less glucose-derived NADPH to regenerate GSH, so its durable antioxidant
    // setpoint (cell.nrf2) is scaled down toward the centre. Geometric (radial
    // depth, no RNG), off-by-default identity ⇒ byte-identical. Runs after
    // contact so the nrf2 axis composes independently of the lipid/iron axes.
    if let Some(cfg) = &nutrient_cfg {
        apply_nutrient_stress_3d(&mut grid, cfg);
    }

    // Phenotype-specific SCD1/MUFA dynamics (#363 rate + #390 cap): scale each
    // tumor cell's MUFA build-up RATE and carrying CAP by per-phenotype
    // multipliers (a drug-tolerant persister remodels lipids toward MUFA at a
    // different rate, and to a different steady state, than a glycolytic rim
    // cell). Reads gc.phenotype (geometric, no RNG); base rate/cap are the run's
    // global scd_mufa_rate / scd_mufa_max. Off-by-default identity is filtered out
    // above, so the default path sets no per-cell mufa_rate/mufa_cap ⇒
    // byte-identical. Runs AFTER the phenotype-reassigning + radial-cap layers
    // (spheroid) so each cell's CURRENT phenotype drives its rate and the cap
    // multiplier COMPOSES with the spheroid's radial cap. Inert in the matrix
    // (scd_mufa_rate = 0 there).
    if let Some(cfg) = &phenotype_mufa_cfg {
        apply_phenotype_mufa_3d(&mut grid, params.scd_mufa_rate, params.scd_mufa_max, cfg);
    }

    // Therapy-induced senescence (#341): mark a `fraction` of tumor cells
    // senescent (independent RNG, matrix-untouched) and apply the multi-axis
    // biochem perturbation (iron + gpx4 + nrf2 + fsp1). The returned per-cell
    // mask drives both the SASP→immune coupling in the kill loop and the
    // `senescent_fraction` report. Runs after nutrient so its per-cell axes
    // compose with the earlier durable layers. Off-by-default identity ⇒
    // byte-identical (the mask is all-false / unused and `apply_*` is a no-op).
    let mut senescence_mask: Option<Vec<bool>> = senescence_cfg
        .as_ref()
        .map(|cfg| apply_senescence_program_3d(&mut grid, cfg, SENESCENCE_SEED));

    // pH-dependent RSL3 ion trapping correction (consumer-side, same pattern
    // as sim-tme). This is the CONSTANT-schedule path: it corrects the
    // one-shot init knockdown for spatial drug availability. For a non-constant
    // schedule the init knockdown was skipped, so this correction would have
    // nothing to correct — instead the per-cell availability is folded into
    // the per-step inactivation via `rsl3_drug_avail` below.
    if !dosed {
        if let (Some((ph_map, cfg)), Treatment::RSL3) = (&ph_field, condition.treatment) {
            for (idx, &local_ph) in ph_map.iter().enumerate() {
                if !grid.cells[idx].is_tumor {
                    continue;
                }
                let drug_factor =
                    ion_trap_factor_from_ph(local_ph, cfg.ph_edge, cfg.ion_trap_sensitivity);
                // Correct GPX4: from (1-inhib) to (1-inhib*drug_factor) — matches sim-tme:614
                let correction =
                    (1.0 - params.rsl3_gpx4_inhib * drug_factor) / (1.0 - params.rsl3_gpx4_inhib);
                grid.cells[idx].state.gpx4 *= correction;
            }
        }
    }

    // Vessel-delivered RSL3 correction (#191), CONSTANT-schedule path. Mirrors
    // the pH ion-trap correction above: cells far from a vessel see less drug,
    // so their one-shot GPX4 knockdown is weaker. Corrects (1-inhib) →
    // (1-inhib*supply). The dosed path instead folds vessel supply into
    // `rsl3_drug_avail` below. `None` ⇒ skipped ⇒ byte-identical. Composes
    // multiplicatively with the pH correction (order immaterial).
    if !dosed {
        if let (Some(supply), Treatment::RSL3) = (&supply_field, condition.treatment) {
            let inhib = params.rsl3_gpx4_inhib;
            for (idx, gc) in grid.cells.iter_mut().enumerate() {
                if gc.is_tumor {
                    let correction = (1.0 - inhib * supply[idx]) / (1.0 - inhib);
                    gc.state.gpx4 *= correction;
                }
            }
        }
    }

    // Per-cell spatial drug-availability multiplier for the DOSED RSL3 path
    // (#239): pH ion-trapping reduces effective drug in acidic regions. The
    // value is composed multiplicatively with the schedule's `factor_at(step)`
    // in the per-step inactivation below. `1.0` everywhere when pH is off;
    // empty (never indexed) on the Constant path or for non-RSL3 treatments.
    let rsl3_drug_avail: Vec<f64> = if dosed && condition.treatment == Treatment::RSL3 {
        let mut avail = vec![1.0_f64; n_cells];
        if let Some((ph_map, cfg)) = &ph_field {
            for (idx, &local_ph) in ph_map.iter().enumerate() {
                if grid.cells[idx].is_tumor {
                    avail[idx] =
                        ion_trap_factor_from_ph(local_ph, cfg.ph_edge, cfg.ion_trap_sensitivity);
                }
            }
        }
        // Supply-delivered drug (#191 vessels / #240 slab depth): scale
        // availability by the unified supply field. `None` ⇒ unchanged.
        if let Some(supply) = &supply_field {
            for (idx, gc) in grid.cells.iter().enumerate() {
                if gc.is_tumor {
                    avail[idx] *= supply[idx];
                }
            }
        }
        avail
    } else {
        Vec::new()
    };

    // --- Main 180-step loop ---
    let immune_cfg = immune_override.unwrap_or_else(SpatialImmuneConfig::for_3d);
    let mut damp_field = vec![0.0_f64; n_cells];
    let mut damp_scratch = vec![0.0_f64; n_cells];
    // IFN-gamma field (#443): allocated only when the coupling is on, else empty so
    // the per-cell application is skipped and the run stays byte-identical.
    let ifngamma_on = ifngamma_cfg.is_some();
    let mut ifngamma_field: Vec<f64> = if ifngamma_on {
        vec![0.0_f64; n_cells]
    } else {
        Vec::new()
    };
    let mut ifngamma_scratch: Vec<f64> = if ifngamma_on {
        vec![0.0_f64; n_cells]
    } else {
        Vec::new()
    };
    let mut ferroptosis_kills = 0usize;
    let mut immune_kills = 0usize;
    // Immune kills among NON-senescent tumor cells specifically (#376). Tracked
    // only when the diffusing SASP field is active, so a strength>0 (suppressive)
    // vs strength<0 (surveillance) A/B can attribute a shift to the BYSTANDER
    // population — proving the paracrine field reaches non-senescent neighbors,
    // not just self-acting on the senescent sources. `None`/0 otherwise.
    let mut nonsenescent_immune_kills = 0usize;

    // T-cell exhaustion (#243). `cumulative_kills[idx]` counts immune kills
    // accumulated in idx's Moore-26 neighborhood; it suppresses that cell's
    // future kill probability via `exhaustion_factor`. Allocated only when
    // exhaustion is enabled (rate > 0) — the default path never touches it,
    // and `exhaustion_factor(_, 0.0) == 1.0` keeps output byte-identical.
    let exhaustion_on = immune_cfg.exhaustion_rate > 0.0;
    let mut cumulative_kills: Vec<u32> = if exhaustion_on {
        vec![0u32; n_cells]
    } else {
        Vec::new()
    };

    // Treg/MDSC suppressor field (#264 Phase 2). A second diffusing field,
    // replenished at the source cells each step, that scales immune kill DOWN
    // locally. Allocated only when suppressor is on; `None` ⇒ never touched and
    // `suppressor_kill_multiplier(_, 0)` would be identity, so byte-identical.
    // Gated on `immune_on`: the field only evolves and only matters inside the
    // immune kill loop (#264 review #2), so an immune-off run neither allocates
    // it nor reports its metrics — the suppressor is meaningless without an
    // immune response to suppress.
    let suppressor_on = suppressor_sources.is_some() && condition.immune_on;
    let suppression_strength = suppressor_cfg.map_or(0.0, |c| c.suppression_strength);
    let mut suppressor_field: Vec<f64> = if suppressor_on {
        vec![0.0_f64; n_cells]
    } else {
        Vec::new()
    };
    let mut suppressor_scratch: Vec<f64> = if suppressor_on {
        vec![0.0_f64; n_cells]
    } else {
        Vec::new()
    };
    // DC subset mix (#264 Phase 4): a cDC1-poor tumor primes anti-tumor killing
    // less efficiently, so DAMP -> kill is scaled by a uniform priming-efficiency
    // scalar. Gated on `immune_on` (it only modulates the immune kill loop) and
    // on a non-identity mix; `None`/balanced ⇒ `dc_priming = 1.0` ⇒ byte-identical.
    let dc_subsets_on = condition.immune_on && dc_subsets_cfg.is_some_and(|c| !c.is_identity());
    // DC ferroptosis susceptibility (#469): gated on `immune_on` (it modulates
    // only the immune kill loop) and on a non-identity config. `dc_ferroptosis_cfg`
    // is already identity-filtered above, so `is_some()` is the enable test.
    let dc_ferroptosis_on = condition.immune_on && dc_ferroptosis_cfg.is_some();
    let (dc_ferro_susceptibility, dc_ferro_pdl1) =
        dc_ferroptosis_cfg.map_or((0.0, 1.0), |c| (c.susceptibility, c.pdl1_protection));
    let dc_priming = if dc_subsets_on {
        dc_subsets_cfg.map_or(1.0, |c| c.priming_efficiency())
    } else {
        1.0
    };

    // Immunosuppressive ferroptosis (#337): a per-cell kill multiplier keyed on
    // the local ferroptotic-death/DAMP signal that scales immune kill DOWN as
    // death density rises (extracellular GPX4 / oxidized-lipid DC suppression),
    // so the net immune effect of dense ferroptotic kill can flip sign. Gated on
    // `immune_on` (it only modulates the immune kill loop); strength 0.0 ⇒
    // `ferroptotic_immunosuppression(_, 0) == 1.0` ⇒ byte-identical.
    let ferro_immuno_strength = immune_cfg.ferro_immunosuppression_strength;
    let ferro_immuno_on = condition.immune_on && ferro_immuno_strength > 0.0;

    // Therapy-induced-senescence SASP to immune coupling (#341): a signed
    // per-cell multiplier on a senescent cell's immune-kill probability. `> 1`
    // is anti-tumor senescence immune surveillance (Kang 2011 PMID 22080947),
    // `< 1` is immunosuppressive SASP that blunts clearance (Di Mitri 2014 PMID
    // 25156255). Gated on `immune_on` (it only modulates the immune kill loop),
    // a populated senescence mask, AND a non-identity multiplier; senescence-off
    // or `mult == 1.0` means no coupling, byte-identical. The mask slice is
    // re-derived per step inside the kill block (not hoisted here) so clonal
    // repopulation can clear revived sites from the mask between steps without a
    // long-lived immutable borrow conflicting with that mutation.
    let sasp_immune_mult = senescence_cfg.map_or(1.0, |c| c.sasp_immune_mult);
    let sasp_on = condition.immune_on && senescence_mask.is_some() && sasp_immune_mult != 1.0;

    // Diffusing SASP field (#376): the paracrine/bystander extension of the
    // cell-autonomous SASP coupling above. Senescent cells are sources; the field
    // diffuses with `diffuse_damp_3d_step` (same operator as the DAMP and
    // suppressor fields) and then scales EVERY exposed cell's immune-kill
    // probability via `sasp_field_kill_multiplier` — including adjacent
    // NON-senescent tumor cells, the neighbor coupling #341 could not express.
    // `strength` is signed: `> 0` immunosuppressive (lowers neighbor kill), `< 0`
    // surveillance (raises it). Gated on `immune_on`, a populated senescence mask,
    // AND a non-zero strength; otherwise the field is never allocated and the
    // multiplier is identity, so the run stays byte-identical.
    let sasp_field_strength = senescence_cfg.map_or(0.0, |c| c.sasp_field_strength);
    let sasp_field_on =
        condition.immune_on && senescence_mask.is_some() && sasp_field_strength != 0.0;
    let mut sasp_field: Vec<f64> = if sasp_field_on {
        vec![0.0_f64; n_cells]
    } else {
        Vec::new()
    };
    let mut sasp_field_scratch: Vec<f64> = if sasp_field_on {
        vec![0.0_f64; n_cells]
    } else {
        Vec::new()
    };

    // The rich kill path handles any realism layer (each factor is identity when
    // its layer is off); the default allocation-free path runs only when ALL are
    // off, staying byte-identical to pre-#243.
    let realism_kill_path = exhaustion_on
        || suppressor_on
        || dc_subsets_on
        || dc_ferroptosis_on
        || ferro_immuno_on
        || sasp_on
        || sasp_field_on;

    // Hoisted out of the per-cell hot loop (these invariants don't vary by
    // cell or step): on the dosed path the relevant availability vec must be
    // populated. Compiled out in release.
    if dosed {
        match condition.treatment {
            Treatment::RSL3 => debug_assert!(
                !rsl3_drug_avail.is_empty(),
                "rsl3_drug_avail must be populated on the dosed RSL3 path"
            ),
            Treatment::SDT | Treatment::PDT => debug_assert!(
                !base_exo.is_empty(),
                "base_exo must be populated on the dosed SDT/PDT path"
            ),
            Treatment::Control => {}
        }
    }

    for step in 0..run_cfg.n_steps {
        // Drug availability for this step (#239). `1.0` on the Constant
        // default path (and the `dosed` guard below skips all modulation
        // there anyway, so it's never actually read for Constant).
        let dose_factor = if dosed { schedule.factor_at(step) } else { 1.0 };

        // Ferroptosis biochem + stromal protection.
        //
        // Parallelized over cells with rayon (#192). This is safe and
        // BYTE-IDENTICAL to the old serial r/c/l triple loop because:
        //   - `enumerate()` over the flat `cells` Vec yields the same `idx`
        //     as `flat_index(r,c,l)` (row-major), so the per-cell RNG seed
        //     `(cond_seed, idx, step)` is unchanged and position-independent;
        //   - each cell reads/writes ONLY its own `GridCell` and its own
        //     `damp_field[idx]` slot (paired via the zipped `par_iter_mut`);
        //     `extra_iron` is read-then-zeroed per own cell (cross-cell iron
        //     spread happens later in the serial `diffuse_iron` phase);
        //   - `ferroptosis_kills` is an integer sum, associative+commutative,
        //     so the reduction is order-independent.
        // Reads of `base_exo` / `rsl3_drug_avail` / `adjacency_mask` are
        // shared-immutable by own index. Iron + DAMP diffusion stay serial
        // (cross-cell deps) and run after the rayon join (a per-step barrier).
        //
        // The `zip` below relies on `cells` and `damp_field` having equal
        // length (both allocated `n_cells`); a shorter `damp_field` would
        // silently TRUNCATE the loop and drop cells. They're allocated full
        // size unconditionally above, but assert it so a future conditional
        // allocation can't introduce a silent correctness bug (review #255).
        debug_assert_eq!(grid.cells.len(), damp_field.len());
        let died_this_step: usize = grid
            .cells
            .par_iter_mut()
            .zip(damp_field.par_iter_mut())
            .enumerate()
            .map(|(idx, (gc, damp_slot))| {
                if !gc.is_tumor {
                    return 0;
                }
                if gc.state.dead {
                    if let Some(ds) = gc.state.death_step {
                        let grace_end = ds + params.post_death_steps;
                        if step == grace_end {
                            // `lp_at_grace_end` captures LP at the end of the
                            // post-death grace period (renamed from the
                            // misleading `lp_at_death` in #314).
                            gc.lp_at_grace_end = gc.state.lp;
                            // DAMP release gated on immune_on (else damp_field
                            // is never read/aggregated — PR #219 third-pass).
                            if condition.immune_on {
                                *damp_slot += gc.lp_at_grace_end * immune_cfg.damp_per_lp;
                            }
                        }
                        if step >= grace_end {
                            return 0;
                        }
                        // else: grace period — fall through to sim_cell_step
                        // for post-death LP accumulation.
                    } else {
                        return 0;
                    }
                }

                let mut rng = StdRng::seed_from_u64(
                    cond_seed
                        .wrapping_add(500_000)
                        .wrapping_add(idx as u64)
                        .wrapping_add(step as u64 * 1_000_000),
                );

                let extra_iron = gc.extra_iron;
                gc.extra_iron = 0.0;
                // O2-gate the neighbor-death-diffused iron too (#383). The Fenton
                // flux is `(cell.iron + extra_iron) · …`; the static `cell.iron`
                // pool was already gated once at setup, but `extra_iron` (iron
                // released from neighbor deaths via `diffuse_iron`) is a per-step
                // pool, so gate it by the SAME local-O2 factor each step to keep
                // the WHOLE Fenton substrate O2-coupled in hypoxia (otherwise the
                // diffused pool would be a residual O2-independent Fenton source
                // in the anoxic core). `fenton_o2_dependence == 0` (default) skips
                // this branch entirely ⇒ `extra_iron` untouched ⇒ byte-identical.
                let extra_iron = if fenton_o2_dependence > 0.0 {
                    extra_iron * fenton_o2_factor(o2_supply_for_exo[idx], fenton_o2_dependence)
                } else {
                    extra_iron
                };

                // Time-varying drug modulation (#239). Skipped on the Constant
                // default path (`dosed == false`) → byte-identical there.
                // Live cells only; a dead cell's stored drug state freezes at
                // death (see the longer note removed from here — preserved in
                // git history; the SDT post-death effective-exo nuance is
                // documented on `biochem::exo_decay_factor`).
                if dosed && !gc.state.dead {
                    match condition.treatment {
                        Treatment::RSL3 => {
                            // Covalent GPX4 inactivation ∝ availability
                            // (schedule × pH ion-trap). Mirrors sim_cell_with_pk.
                            // NOTE: the persister block below recomputes this
                            // same `dose_factor * rsl3_drug_avail[idx]` as its
                            // `drug_intensity`; keep the two in sync if edited.
                            let conc = (dose_factor * rsl3_drug_avail[idx]).clamp(0.0, 1.0);
                            // Persisters resist the covalent knockdown (#241).
                            // Reads `persister_fraction` as of the PREVIOUS step
                            // (the persister block writes it *after* this); this
                            // one-step explicit-Euler lag is intentional —
                            // reordering it would change the dosed-path output.
                            // `1.0` when persister off → byte-identical.
                            let presist = persister_cfg.as_ref().map_or(1.0, |c| {
                                persister::gpx4_inactivation_multiplier(
                                    gc.state.persister_fraction,
                                    c,
                                )
                            });
                            gc.state.gpx4 -=
                                RSL3_INACTIVATION_RATE * conc * presist * gc.state.gpx4;
                            gc.state.gpx4 = gc.state.gpx4.max(0.0);
                        }
                        Treatment::SDT | Treatment::PDT => {
                            // Divide out sim_cell_step's intrinsic single-bolus
                            // envelope so the schedule is the sole envelope
                            // (#239). `.max(1e-9)` guards 0/0; the ~1000×
                            // amplification at late steps is finite and never
                            // escapes sim_cell_step's immediate re-multiply.
                            let decay = exo_decay_factor(step).max(1e-9);
                            gc.state.exo_ros_peak = base_exo[idx] * dose_factor / decay;
                        }
                        Treatment::Control => {}
                    }
                }

                // IFN-gamma -> ACSL4 arm (#443 follow-up): transiently raise THIS
                // cell's PUFA / lipid unsaturation for THIS step by the local
                // IFN-gamma boost (more ACSL4 -> more oxidizable substrate -> more
                // ferroptosis; Wang 2019 PMID 31043744). Applied as a
                // save/multiply/restore around sim_cell_step, NOT a durable mutation:
                // lipid_unsat is a static per-cell parameter, so a persistent
                // per-step multiply would compound geometrically over the run
                // (unlike the regenerating GSH pool the System Xc- arm scales). The
                // transient form also composes cleanly with the contact/clonal
                // layers that set lipid_unsat. Off (acsl4_strength == 0, ifngamma
                // disabled, or dead cell) ⇒ boost is exactly 1.0 ⇒ no mutation ⇒
                // byte-identical. acsl4_upregulation(0, s) == 1.0 exactly, so cells
                // outside the IFN-gamma field are also untouched.
                let acsl4_boost = if ifngamma_acsl4_on && !gc.state.dead {
                    acsl4_upregulation(ifngamma_field[idx], ifngamma_acsl4_strength)
                } else {
                    1.0
                };
                let saved_unsat = gc.cell.lipid_unsat;
                if acsl4_boost != 1.0 {
                    gc.cell.lipid_unsat *= acsl4_boost;
                }
                // Persister OXPHOS-ROS suppression (#470): a drug-tolerant
                // persister downregulates OXPHOS (a main mitochondrial-ROS
                // source), so a GPX4 inhibitor has less peroxidizable-substrate
                // flux to act on and kills it less; an HDAC inhibitor re-raises
                // the ROS and restores the kill (PMID 40909720). Applied as a
                // transient save/multiply/restore of THIS cell's basal ROS for
                // THIS step (NOT durable: basal_ros is a static per-cell
                // parameter while persister_fraction evolves per step, so a
                // persistent multiply would compound). Reads persister_fraction
                // as of the PREVIOUS step, like the GPX4-resistance coupling
                // above. 1.0 when persister off OR oxphos_ros_suppression == 0
                // (incl. the default `enabled()` persister) ⇒ no mutation ⇒
                // byte-identical.
                let oxphos_mult = persister_cfg
                    .as_ref()
                    .filter(|_| !gc.state.dead)
                    .map_or(1.0, |c| {
                        persister::oxphos_ros_multiplier(gc.state.persister_fraction, c)
                    });
                let saved_basal_ros = gc.cell.basal_ros;
                if oxphos_mult != 1.0 {
                    gc.cell.basal_ros *= oxphos_mult;
                }
                let died =
                    sim_cell_step(&mut gc.state, &gc.cell, &params, step, extra_iron, &mut rng);
                if oxphos_mult != 1.0 {
                    gc.cell.basal_ros = saved_basal_ros;
                }
                if acsl4_boost != 1.0 {
                    gc.cell.lipid_unsat = saved_unsat;
                }
                if died {
                    // `newly_dead` is consumed by the later serial
                    // `diffuse_iron` to spread released iron to live neighbors.
                    gc.newly_dead = true;
                }

                // Stromal CAF protection for alive cells (stromal_on only).
                if !died && !gc.state.dead {
                    if let Some(cfg) = &stromal_cfg {
                        if adjacency_mask[idx] {
                            gc.state.gsh =
                                (gc.state.gsh + cfg.gsh_boost_per_step).min(cfg.gsh_boost_cap);
                            gc.state.mufa_protection = (gc.state.mufa_protection
                                + cfg.mufa_boost_per_step)
                                .min(cfg.mufa_boost_cap);
                        }
                    }
                }

                // IFN-gamma -> System Xc- downregulation (#443): scale the GSH pool
                // DOWN by the local IFN-gamma retention factor (cystine starvation),
                // sensitizing to ferroptosis (Wang 2019 PMID 31043744). The IFN-gamma
                // field reflects the previous step's diffused immune signal. Gated on
                // `ifngamma_on`; the field is empty (never indexed) when off, so the
                // matrix is byte-identical.
                if ifngamma_on && !died && !gc.state.dead {
                    gc.state.gsh *= system_xc_retention(ifngamma_field[idx], ifngamma_ic50);
                }

                // Copper-ionophore / cuproptosis crosstalk (#485): a copper
                // ionophore (elesclomol/disulfiram) overloads intracellular
                // copper, depleting the GSH and GPX4 pools each step (copper binds
                // GSH + drives GPX4 degradation), so RSL3/SDT kills MORE; ATP7B
                // efflux exports copper and protects. Uniform (no field). Gated on
                // `copper_on` (a non-identity config); 1.0 retention when off ⇒
                // byte-identical.
                if copper_on && !died && !gc.state.dead {
                    gc.state.gsh *= copper_gsh_ret;
                    gc.state.gpx4 *= copper_gpx4_ret;
                }

                // Persister-cell dynamics (#241). Gated on the model being
                // enabled (`Some`); `None` on the default matrix path so this
                // whole block is skipped → byte-identical. Applies to cells
                // still alive after this step's biochem.
                if let Some(pcfg) = persister_cfg.as_ref() {
                    if !gc.state.dead {
                        let frac = gc.state.persister_fraction;
                        // MUFA membrane remodeling (lipid-rewiring axis):
                        // additive per step like CAF supply, scaled by the
                        // current persister fraction, capped. `.max(m)` so it
                        // never pulls existing protection down.
                        let inc = persister::mufa_boost_increment(frac, pcfg);
                        let m = gc.state.mufa_protection;
                        gc.state.mufa_protection = (m + inc).min(pcfg.mufa_boost_cap.max(m));
                        // Per-step drug intensity drives the competing-rate
                        // persister update below (the `persister::step` call,
                        // #262); see its comment for the acquisition/reversion
                        // equilibrium semantics.
                        //
                        // `rsl3_drug_avail[idx]` is indexed only on the
                        // `dosed && RSL3` path; it is always populated there by
                        // the RSL3-treatment branch that allocates it (the same
                        // invariant the GPX4 inactivation above relies on — the
                        // pre-loop `debug_assert` only documents it in debug).
                        let drug_intensity = match condition.treatment {
                            Treatment::Control => 0.0,
                            Treatment::RSL3 => {
                                if dosed {
                                    (dose_factor * rsl3_drug_avail[idx]).clamp(0.0, 1.0)
                                } else {
                                    1.0
                                }
                            }
                            Treatment::SDT | Treatment::PDT => {
                                if dosed {
                                    dose_factor
                                } else {
                                    1.0
                                }
                            }
                        };
                        // Competing-rate update (#262) with reversible-to-
                        // irreversible epigenetic locking (#342). Reconstruct the
                        // persister state stored on the cell (reversible + locked
                        // pools + the sustained-exposure tracker), advance it one
                        // step, and cache its total() back into
                        // `persister_fraction` so the GPX4/MUFA couplings + the
                        // mean/snapshot keep reading one scalar. With
                        // `lock_rate == 0` (the default, including enabled())
                        // `step_with_locking` reduces EXACTLY to `step` on the
                        // reversible pool and total() == the old reversible value,
                        // so persister_fraction tracks the pre-#342 sequence and
                        // the matrix/golden runs stay byte-identical.
                        let pstate = persister::PersisterState {
                            reversible: gc.state.persister_reversible,
                            locked: gc.state.persister_locked,
                            cumulative_exposure: gc.state.persister_cum_exposure,
                        };
                        let pstate = persister::step_with_locking(pstate, drug_intensity, pcfg);
                        // #377: non-drug stress-niche entry. Hypoxic / nutrient-poor
                        // drug-sanctuary zones drive persister entry INDEPENDENT of
                        // drug; the stress signal is the local hypoxia deficit
                        // (1 - o2_supply), the same per-cell O2 field that gates the
                        // SDT exo-ROS (#336) and the Fenton substrate (#383). Applied
                        // AFTER the drug step, it raises only the reversible pool (not
                        // the locking EMA / resistance). `stress_entry_rate == 0`
                        // (default) ⇒ no-op ⇒ byte-identical.
                        let stress = (1.0 - o2_supply_for_exo[idx]).clamp(0.0, 1.0);
                        let pstate = persister::stress_entry(pstate, stress, pcfg);
                        gc.state.persister_reversible = pstate.reversible;
                        gc.state.persister_locked = pstate.locked;
                        gc.state.persister_cum_exposure = pstate.cumulative_exposure;
                        gc.state.persister_fraction = pstate.total(pcfg);
                    }
                }

                usize::from(died)
            })
            .sum();
        ferroptosis_kills += died_this_step;

        // Iron diffusion via TumorGrid3D. Now uses the value from
        // `spatial_params.neighbor_iron_fraction` (single source of truth);
        // the inline `0.031` is gone.
        grid.diffuse_iron(
            spatial_params.iron_release_per_death,
            spatial_params.neighbor_iron_fraction,
        );

        // DAMP diffusion if immune is enabled
        if condition.immune_on {
            diffuse_damp_3d_step(
                &mut damp_field,
                &mut damp_scratch,
                &grid,
                immune_cfg.damp_diffusion_fraction,
                immune_cfg.damp_clearance_rate,
            );

            // Treg/MDSC suppressor field (#264 Phase 2): replenish at the source
            // (niche) cells, then diffuse + clear (same scratch-buffer step as
            // DAMP). Off ⇒ skipped ⇒ byte-identical.
            if let (Some(scfg), Some(sources)) = (&suppressor_cfg, &suppressor_sources) {
                for (idx, &is_src) in sources.iter().enumerate() {
                    if is_src {
                        suppressor_field[idx] =
                            (suppressor_field[idx] + scfg.replenish_rate).min(1.0);
                    }
                }
                diffuse_damp_3d_step(
                    &mut suppressor_field,
                    &mut suppressor_scratch,
                    &grid,
                    scfg.diffusion_fraction,
                    scfg.clearance_rate,
                );
            }

            // Diffusing SASP field (#376): replenish at the senescent source
            // cells (the senescence mask), then diffuse + clear with the same
            // operator. The field then modulates EVERY exposed cell's kill in the
            // loop below (neighbor/bystander coupling). Off ⇒ skipped ⇒
            // byte-identical. Uses the current mask: a cell cleared from the mask
            // by clonal repopulation simply stops seeding new SASP (already-
            // diffused field clears at `SASP_FIELD_CLEARANCE_RATE`).
            if sasp_field_on {
                if let Some(mask) = &senescence_mask {
                    for (idx, &is_sen) in mask.iter().enumerate() {
                        if is_sen {
                            sasp_field[idx] =
                                (sasp_field[idx] + SASP_FIELD_REPLENISH_RATE).min(1.0);
                        }
                    }
                    diffuse_damp_3d_step(
                        &mut sasp_field,
                        &mut sasp_field_scratch,
                        &grid,
                        SASP_FIELD_DIFFUSION_FRACTION,
                        SASP_FIELD_CLEARANCE_RATE,
                    );
                }
            }

            // IFN-gamma field (#443): seed from the local DAMP at immune-active
            // (DAMP-positive) positions, coupling IFN-gamma secretion to where the
            // T-cell response is concentrated, then diffuse + clear with the same
            // operator. The field suppresses GSH (System Xc- downregulation) in the
            // next step's biochem loop. Off ⇒ skipped ⇒ byte-identical.
            if let Some(cfg) = &ifngamma_cfg {
                for idx in 0..n_cells {
                    if damp_field[idx] > DAMP_KILL_THRESHOLD {
                        ifngamma_field[idx] =
                            (ifngamma_field[idx] + cfg.per_damp * damp_field[idx]).min(1.0);
                    }
                }
                diffuse_damp_3d_step(
                    &mut ifngamma_field,
                    &mut ifngamma_scratch,
                    &grid,
                    cfg.diffusion_fraction,
                    cfg.clearance_rate,
                );
            }

            // Immune kill (after delay). Parallelized over cells with rayon
            // (#192) — byte-identical to the old serial triple loop: each cell
            // reads its own `damp_field[idx]` (immutable here; DAMP diffusion
            // already done) and writes only its own `state.dead`; the per-cell
            // RNG seed `(cond_seed, idx, step)` is position-independent; and
            // `immune_kills` is an order-independent integer sum.
            if step >= IMMUNE_START_STEP {
                // Multi-checkpoint panel (#264 Phase 3) replaces the single PD-1
                // brake when set; `None` ⇒ the single-PD-1 `effective_brake`
                // (byte-identical). A uniform scalar, hoisted out of the hot loop.
                let effective_brake =
                    checkpoint_panel.map_or(immune_cfg.effective_brake(), |p| p.combined_brake());
                // Two paths so the default matrix stays exactly the pre-#243
                // loop (allocation-free `map().sum()`, no exhaustion term) —
                // not just byte-identical output but the same hot-path code,
                // honoring the #192/#253 zero-cost discipline. The exhaustion
                // path instead `filter_map`-collects the killed flat indices so
                // their neighborhoods can be scattered into `cumulative_kills`.
                if realism_kill_path {
                    // SASP mask slice, re-derived each step so a borrow does not
                    // outlive this block (clonal repopulation clears revived bits
                    // from `senescence_mask` later this step). Empty (never
                    // indexed) unless `sasp_on`; `sasp_on` guarantees length
                    // `n_cells`, so the `[idx]` below is in bounds.
                    let senescence_mask_slice: &[bool] = senescence_mask.as_deref().unwrap_or(&[]);
                    debug_assert!(
                        !sasp_on || senescence_mask_slice.len() == n_cells,
                        "sasp_on implies a full-length senescence mask"
                    );
                    let killed: Vec<usize> = grid
                        .cells
                        .par_iter_mut()
                        .enumerate()
                        .filter_map(|(idx, gc)| {
                            if gc.state.dead || !gc.is_tumor {
                                return None;
                            }
                            let local_damp = damp_field[idx];
                            if local_damp < DAMP_KILL_THRESHOLD {
                                return None;
                            }
                            let activation = dc_activation(local_damp, immune_cfg.dc_activation_kd);
                            // Two opposing local modulators, each identity when
                            // its layer is off (its backing vec is empty then, so
                            // the `_on` guard avoids an out-of-bounds index):
                            // T-cell exhaustion (#243) and Treg/MDSC suppression
                            // (#264) both scale the kill probability DOWN.
                            let exh = if exhaustion_on {
                                exhaustion_factor(cumulative_kills[idx], immune_cfg.exhaustion_rate)
                            } else {
                                1.0
                            };
                            let supp = if suppressor_on {
                                suppressor_kill_multiplier(
                                    suppressor_field[idx],
                                    suppression_strength,
                                )
                            } else {
                                1.0
                            };
                            // Immunosuppressive ferroptosis (#337): the same
                            // local DAMP that drives `activation` also scales
                            // kill DOWN via the suppressive arm; 1.0 when off.
                            let ferro_supp = if ferro_immuno_on {
                                ferroptotic_immunosuppression(local_damp, ferro_immuno_strength)
                            } else {
                                1.0
                            };
                            // SASP to immune coupling (#341): for a senescent
                            // cell, a signed multiplier on its immune-kill
                            // probability (surveillance > 1, immunosuppression
                            // < 1); 1.0 for non-senescent cells and layer-off.
                            let sasp = if sasp_on {
                                sasp_immune_multiplier(senescence_mask_slice[idx], sasp_immune_mult)
                            } else {
                                1.0
                            };
                            // Diffusing SASP field (#376): the PARACRINE arm,
                            // signed and reaching EVERY exposed cell (not just
                            // senescent ones) — `> 0` strength lowers a neighbor's
                            // kill (immunosuppressive), `< 0` raises it
                            // (surveillance); 1.0 for field 0 and layer-off. The
                            // `sasp_field_on` guard keeps the empty backing vec
                            // out of the index on the default path.
                            let sasp_field_mult = if sasp_field_on {
                                sasp_field_kill_multiplier(sasp_field[idx], sasp_field_strength)
                            } else {
                                1.0
                            };
                            // DC ferroptosis susceptibility (#469): the SAME
                            // local DAMP/lipid-ROS that drives `activation` also
                            // kills the effector DCs themselves, scaling the kill
                            // probability DOWN (lower priming); PD-L1 protection
                            // restores it. 1.0 when off (identity).
                            let dc_ferro = if dc_ferroptosis_on {
                                dc_ferroptosis_survival(
                                    local_damp,
                                    dc_ferro_susceptibility,
                                    dc_ferro_pdl1,
                                )
                            } else {
                                1.0
                            };
                            // `dc_priming` is a uniform scalar (the cDC1/cDC2
                            // mix, #264 Phase 4): 1.0 when off, < 1.0 for a
                            // cDC1-poor tumor that primes killing less efficiently.
                            let kill_prob = immune_kill_probability(
                                activation,
                                immune_cfg.immune_kill_rate,
                                effective_brake,
                            ) * exh
                                * supp
                                * dc_priming
                                * dc_ferro
                                * ferro_supp
                                * sasp
                                * sasp_field_mult;
                            let mut rng = StdRng::seed_from_u64(
                                cond_seed
                                    .wrapping_add(900_000_000)
                                    .wrapping_add(idx as u64)
                                    .wrapping_add(step as u64 * 2_000_000),
                            );
                            if rng.gen::<f64>() < kill_prob {
                                gc.state.dead = true;
                                Some(idx)
                            } else {
                                None
                            }
                        })
                        .collect();
                    immune_kills += killed.len();

                    // Bystander accounting (#376): count this step's immune kills
                    // that landed on NON-senescent tumor cells. Only when the SASP
                    // field is on; `sasp_field_on` ⇒ `senescence_mask.is_some()` ⇒
                    // `senescence_mask_slice` is full-length, so the index is safe.
                    if sasp_field_on {
                        nonsenescent_immune_kills += killed
                            .iter()
                            .filter(|&&k| !senescence_mask_slice[k])
                            .count();
                    }

                    // Scatter exhaustion into the Moore-26 neighborhood of each
                    // cell killed this step (serial — runs after the par_iter
                    // join; integer adds commute, so order-independent and
                    // deterministic regardless of `collect` order). Only when
                    // exhaustion is on (suppressor-only runs leave it empty).
                    if exhaustion_on {
                        for &k in &killed {
                            let (r, c, l) = grid.coords(k);
                            let (nbrs, n) = grid.neighbors(r, c, l);
                            for &(nr, nc, nl) in &nbrs[..n] {
                                cumulative_kills[grid.flat_index(nr, nc, nl)] += 1;
                            }
                        }
                    }
                } else {
                    // DEFAULT path — identical to pre-#243: allocation-free,
                    // no exhaustion term, byte-identical output.
                    let killed_this_step: usize = grid
                        .cells
                        .par_iter_mut()
                        .enumerate()
                        .map(|(idx, gc)| {
                            if gc.state.dead || !gc.is_tumor {
                                return 0;
                            }
                            let local_damp = damp_field[idx];
                            if local_damp < DAMP_KILL_THRESHOLD {
                                return 0;
                            }
                            let activation = dc_activation(local_damp, immune_cfg.dc_activation_kd);
                            let kill_prob = immune_kill_probability(
                                activation,
                                immune_cfg.immune_kill_rate,
                                effective_brake,
                            );
                            let mut rng = StdRng::seed_from_u64(
                                cond_seed
                                    .wrapping_add(900_000_000)
                                    .wrapping_add(idx as u64)
                                    .wrapping_add(step as u64 * 2_000_000),
                            );
                            if rng.gen::<f64>() < kill_prob {
                                // Immune kills set `state.dead` but deliberately
                                // NOT `death_step`/`newly_dead`: immune-killed
                                // cells are apoptotic (no post-death LP grace,
                                // no DAMP burst, no iron release) — a modeling
                                // choice consistent with sim-tme.
                                gc.state.dead = true;
                                1
                            } else {
                                0
                            }
                        })
                        .sum();
                    immune_kills += killed_this_step;
                }
            }
        }

        // Spatial clonal expansion (#266 item 3): after all deaths this step,
        // repopulate dead tumor sites from living Moore-neighbors so resistant
        // subclones (more survivors ⇒ more donors) grow their territory. Gated
        // on clonal + repopulation_rate > 0 ⇒ off-by-default byte-identical. The
        // per-site RNG derives from `cond_seed` (per-step), distinct from the
        // setup seeds, so it never perturbs the assignment/placement streams.
        if let (Some(ids), Some(ccfg)) = (subclone_ids.as_mut(), &clonal_cfg) {
            if ccfg.repopulation_rate > 0.0 {
                let revived = repopulate_dead_sites_3d(
                    &mut grid,
                    ids,
                    ccfg,
                    &params,
                    cond_seed.wrapping_add(700_000_000),
                    step,
                );
                // #302: a revived dead site is a fresh, full-strength cell. If
                // the contact layer is on, re-apply its geometric resistance to
                // each revived cell so a dense interior site resists like its
                // neighbours (reduced at setup), instead of coming back as the
                // MOST ferroptosis-sensitive cell in the cluster — the opposite
                // of the modeled biology. Contact fractions are is_tumor-
                // geometric (death/repopulation never change them), so this
                // reproduces the setup-time reduction exactly. (clonal+contact
                // is otherwise unguarded, unlike slab+contact, because the two
                // are meant to compose.)
                if let Some(contact) = &contact_cfg {
                    for &idx in &revived {
                        apply_contact_resistance_at_3d(&mut grid, idx, contact);
                    }
                }
                // #363: a revived cell is fresh (gen_cell sets mufa_rate=None ⇒
                // the global rate). If the phenotype-MUFA layer is on, re-derive
                // the revived cell's per-phenotype rate from its CURRENT phenotype
                // so clonal repopulation stays coherent with #363 (same rationale
                // as the #302 contact re-application above). No-op unless the
                // layer is on.
                if let Some(pm) = &phenotype_mufa_cfg {
                    for &idx in &revived {
                        apply_phenotype_mufa_at_3d(
                            &mut grid,
                            idx,
                            params.scd_mufa_rate,
                            params.scd_mufa_max,
                            pm,
                        );
                    }
                }
                // #341: a revived dead site is a NEW cell grown from a living
                // neighbour, not the resurrection of the senescent cell that died
                // there, and its biochem was reset to a fresh cell. Clear its
                // senescence-mask bit so the SASP immune coupling and the
                // `senescent_fraction` report stop counting it as senescent (the
                // analogue of the #302 contact re-application: keep the senescence
                // layer coherent with clonal repopulation instead of leaving a
                // stale mask). No-op unless senescence is also on.
                if let Some(mask) = senescence_mask.as_mut() {
                    for &idx in &revived {
                        mask[idx] = false;
                    }
                }
            }
        }

        // Per-step trajectory capture for `--snapshot` runs. No-op for the
        // default 24-condition matrix path (snapshot is None there).
        // Captured *after* all per-step work (biochem + DAMP + immune kill +
        // any clonal repopulation) so the snapshot reflects end-of-step state.
        if let Some(buf) = snapshot.as_mut() {
            buf.capture_step(&grid, &damp_field);
        }
    }

    // Late DAMP release for cells still in their post-death grace period at
    // the end of the simulation. Iterate paired with damp_field by index to
    // satisfy clippy::needless_range_loop while keeping the dual-Vec access.
    // Run only when immune was on — otherwise damp_field is never read or
    // aggregated and this loop is dead writes (reviewer-flagged).
    if condition.immune_on {
        for (idx, gc) in grid.cells.iter_mut().enumerate() {
            if !(gc.is_tumor && gc.state.dead) {
                continue;
            }
            if let Some(ds) = gc.state.death_step {
                let grace_end = ds + params.post_death_steps;
                if grace_end >= run_cfg.n_steps {
                    gc.lp_at_grace_end = gc.state.lp;
                    damp_field[idx] += gc.lp_at_grace_end * immune_cfg.damp_per_lp;
                }
            }
        }
    }

    // --- Aggregate results ---
    let census = grid.census();
    let overall = census.total_dead as f64 / census.total_tumor.max(1) as f64;
    let (norm_r, trans_r, hyp_r) = zone_kill_rates_3d(&grid, ZONE_REF_LAMBDA);

    // DAMP aggregations — vindicate the PR-body's "DAMP fields" claim and
    // give the comparison script real numbers to work with. When immune
    // was off these are both 0.0 (damp_field was zero-initialized and
    // never written by the diffusion path).
    let peak_damp = damp_field.iter().cloned().fold(0.0_f64, f64::max);
    let total_damp: f64 = damp_field.iter().sum();

    // Stromal adjacency rates — ALWAYS computed (mask is independent of
    // the stromal_on toggle), so the comparison script's Q3 can pair
    // stromal-on adjacent rates with no-stromal baseline adjacent rates.
    // `stromal_mode` still only set when stromal_on is true (matches
    // sim-tme JSON convention).
    let mut adj_dead = 0usize;
    let mut adj_total = 0usize;
    for (idx, gc) in grid.cells.iter().enumerate() {
        if gc.is_tumor && adjacency_mask[idx] {
            adj_total += 1;
            if gc.state.dead {
                adj_dead += 1;
            }
        }
    }
    let stromal_adjacent_kill_rate = if adj_total > 0 {
        Some(adj_dead as f64 / adj_total as f64)
    } else {
        None
    };
    let stromal_adjacent_count = if adj_total > 0 { Some(adj_total) } else { None };
    let stromal_mode = if condition.stromal_on {
        Some("stromal_on".to_string())
    } else {
        None
    };

    let (ph_mode, ph_edge_v, ph_core_v, ph_lambda_v) = if let Some((_, cfg)) = &ph_field {
        (
            Some("ph_on".to_string()),
            Some(cfg.ph_edge),
            Some(cfg.ph_core),
            Some(cfg.lambda_ph_um),
        )
    } else {
        (None, None, None, None)
    };

    // Mean persister fraction over surviving tumor cells (#241). `Some` only
    // when the model is enabled — `None` omits the field (skip_serializing_if)
    // so the default-matrix summary.json is byte-identical to pre-#241.
    let persister_mean = persister_cfg.map(|_| {
        let (sum, n) = grid.cells.iter().fold((0.0_f64, 0usize), |(s, n), gc| {
            if gc.is_tumor && !gc.state.dead {
                (s + gc.state.persister_fraction, n + 1)
            } else {
                (s, n)
            }
        });
        if n > 0 {
            sum / n as f64
        } else {
            0.0
        }
    });
    // Mean LOCKED (irreversible) persister fraction (#342). Same surviving-tumor
    // denominator as `persister_mean`; reads the locked sub-pool stored on each
    // cell. `0` whenever `lock_rate == 0` (the default / `enabled()`), so a
    // lock-off persister run reports `Some(0.0)` and the matrix (persister off)
    // omits it entirely (byte-identical).
    let persister_locked_mean = persister_cfg.map(|_| {
        let (sum, n) = grid.cells.iter().fold((0.0_f64, 0usize), |(s, n), gc| {
            if gc.is_tumor && !gc.state.dead {
                (s + gc.state.persister_locked, n + 1)
            } else {
                (s, n)
            }
        });
        if n > 0 {
            sum / n as f64
        } else {
            0.0
        }
    });

    // Per-subclone kill breakdown (#242). `Some` only when clonal is enabled
    // (subclone_ids present), so `None` omits the field → byte-identical.
    let subclone_kills = subclone_ids.as_ref().map(|ids| {
        let k = clonal_cfg.as_ref().map_or(0, |c| c.k());
        // Index 0 unused (stroma); 1..=k are the subclones.
        let mut totals = vec![0usize; k + 1];
        let mut deads = vec![0usize; k + 1];
        for (idx, gc) in grid.cells.iter().enumerate() {
            if gc.is_tumor {
                let s = ids[idx] as usize;
                totals[s] += 1;
                if gc.state.dead {
                    deads[s] += 1;
                }
            }
        }
        (1..=k)
            .map(|s| SubcloneKillStat {
                subclone_id: s as u8,
                initial_tumor: initial_subclone_totals.as_ref().map_or(totals[s], |t| t[s]),
                total_tumor: totals[s],
                total_dead: deads[s],
                kill_rate: if totals[s] > 0 {
                    deads[s] as f64 / totals[s] as f64
                } else {
                    0.0
                },
            })
            .collect::<Vec<_>>()
    });

    // Supply-field hypoxic fraction (#191 vessels / #240 slab depth). `Some`
    // only when a supply field is active → `None` omits the field on the
    // default/spheroid path → byte-identical.
    let vascular_hypoxic_fraction = supply_field
        .as_ref()
        .map(|s| hypoxic_fraction(&grid, s, VASCULAR_HYPOXIC_THRESHOLD));

    // Patient-scale slab interpretation (#240). `Some` only in slab mode.
    let scale_interpretation_str = slab_cfg.map(|cfg| scale_interpretation(&grid, &cfg));

    // Treg/MDSC suppressor reporting (#264). `Some` only when the suppressor was
    // active (sources present AND immune_on, so the field actually evolved).
    // The field replenishes every step and reaches a quasi-steady state, so the
    // end-of-run field max is representative of the run-time peak.
    let suppressor_source_count = suppressor_on.then(|| {
        suppressor_sources
            .as_ref()
            .map_or(0, |s| s.iter().filter(|&&v| v).count())
    });
    let suppressor_peak =
        suppressor_on.then(|| suppressor_field.iter().cloned().fold(0.0_f64, f64::max));
    // Multi-checkpoint combined brake (#264 Phase 3). `Some` only when a panel
    // override is supplied; the single-PD-1 path omits it.
    let checkpoint_brake = checkpoint_panel.map(|p| p.combined_brake());
    // Senescent fraction (#341): fraction of tumor cells still flagged
    // senescent at end-of-run. `Some` only when senescence is enabled (mask
    // present), so `None` omits the field and the default path stays
    // byte-identical. Marking is a one-time setup assignment, but clonal
    // repopulation (if on) clears revived sites from the mask, so this reflects
    // the end-of-run senescent set rather than the initial one.
    let senescent_fraction = senescence_mask.as_ref().map(|mask| {
        let (count, total) = grid
            .cells
            .iter()
            .zip(mask.iter())
            .filter(|(gc, _)| gc.is_tumor)
            .fold((0usize, 0usize), |(c, t), (_, &is_sen)| {
                (c + usize::from(is_sen), t + 1)
            });
        if total > 0 {
            count as f64 / total as f64
        } else {
            0.0
        }
    });

    ConditionResult {
        treatment: condition.treatment_name.clone(),
        o2_condition: match condition.o2_lambda {
            None => "uniform".to_string(),
            Some(_) => "gradient".to_string(),
        },
        o2_lambda_um: condition.o2_lambda,
        // Match sim-tme's JSON convention exactly (sim-tme uses "immune_on" /
        // "immune_anti_pd1" / "off"). My earlier "on"/"off" prevented the
        // comparison script from finding 2D immune conditions — caught in
        // adversarial review.
        immune_mode: if condition.immune_on {
            "immune_on".to_string()
        } else {
            "off".to_string()
        },
        total_tumor: census.total_tumor,
        total_dead: census.total_dead,
        ferroptosis_kills: if condition.immune_on {
            Some(ferroptosis_kills)
        } else {
            None
        },
        immune_kills: if condition.immune_on {
            Some(immune_kills)
        } else {
            None
        },
        overall_kill_rate: overall,
        normoxic_kill_rate: norm_r,
        transition_kill_rate: trans_r,
        hypoxic_kill_rate: hyp_r,
        peak_damp,
        total_damp,
        stromal_mode,
        stromal_adjacent_kill_rate,
        stromal_adjacent_count,
        ph_mode,
        ph_edge: ph_edge_v,
        ph_core: ph_core_v,
        ph_lambda_um: ph_lambda_v,
        persister_mean,
        persister_locked_mean,
        subclone_kills,
        vascular_hypoxic_fraction,
        scale_interpretation: scale_interpretation_str,
        suppressor_source_count,
        suppressor_peak,
        checkpoint_brake,
        senescent_fraction,
        nonsenescent_immune_kills: if sasp_field_on {
            Some(nonsenescent_immune_kills)
        } else {
            None
        },
    }
}

// Cheap deterministic hash of condition name → seed offset.
fn hash_condition_name(name: &str) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325; // FNV-1a offset
    for b in name.bytes() {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

// ============================================================
// Condition matrix
// ============================================================

fn generate_conditions() -> Vec<Condition> {
    let mut conditions = Vec::new();
    let treatments = [
        (Treatment::Control, "Control"),
        (Treatment::RSL3, "RSL3"),
        (Treatment::SDT, "SDT"),
    ];

    // Baseline: uniform O₂, no immune/stromal/pH
    for (tx, name) in &treatments {
        conditions.push(Condition {
            name: format!("baseline_{}", name),
            treatment: *tx,
            treatment_name: name.to_string(),
            o2_lambda: None,
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        });
    }

    // O₂ gradient sweep (no immune/stromal/pH yet)
    for &lambda in O2_LAMBDAS {
        for (tx, name) in &treatments {
            conditions.push(Condition {
                name: format!("o2_{}_{}", lambda as i32, name),
                treatment: *tx,
                treatment_name: name.to_string(),
                o2_lambda: Some(lambda),
                immune_on: false,
                stromal_on: false,
                ph_on: false,
                dose_schedule: DoseSchedule::Constant,
            });
        }
    }

    // Immune on (at ZONE_REF_LAMBDA O₂)
    for (tx, name) in &treatments {
        conditions.push(Condition {
            name: format!("immune_{}", name),
            treatment: *tx,
            treatment_name: name.to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        });
    }

    // Stromal on + immune on (matches sim-tme's coupling: 2D stromal-on rows
    // all have `immune_mode == "immune_on"`. Standalone stromal-off-immune-off
    // wouldn't exist in 2D, so the comparison script couldn't pair against
    // anything matching — reviewer-flagged confounding in PR #219 review).
    for (tx, name) in &treatments {
        conditions.push(Condition {
            name: format!("stromal_{}", name),
            treatment: *tx,
            treatment_name: name.to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: true,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        });
    }

    // pH on + immune on (same coupling as stromal — sim-tme's pH-on rows
    // all have `immune_mode == "immune_on"`).
    for (tx, name) in &treatments {
        conditions.push(Condition {
            name: format!("ph_{}", name),
            treatment: *tx,
            treatment_name: name.to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: false,
            ph_on: true,
            dose_schedule: DoseSchedule::Constant,
        });
    }

    // Combined: immune + stromal + pH (at ZONE_REF_LAMBDA O₂)
    for (tx, name) in &treatments {
        conditions.push(Condition {
            name: format!("combined_{}", name),
            treatment: *tx,
            treatment_name: name.to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: true,
            ph_on: true,
            dose_schedule: DoseSchedule::Constant,
        });
    }

    conditions
}

// ============================================================
// Main
// ============================================================

/// Parse a `--snapshot[=NAME]` argument from the CLI args iterator.
/// Returns `Some(name)` if the flag is present (defaults to
/// `"combined"` if no `=NAME` suffix), `None` otherwise. Unknown
/// preset names are validated downstream by [`resolve_snapshot`].
fn parse_snapshot_arg<I: IntoIterator<Item = String>>(args: I) -> Option<String> {
    for a in args {
        if a == "--snapshot" {
            return Some("combined".to_string());
        }
        if let Some(name) = a.strip_prefix("--snapshot=") {
            return Some(name.to_string());
        }
    }
    None
}

/// One entry in the `--snapshot=NAME` registry. Each entry pins a
/// reproducible condition for visualization. Names are stable; adding
/// a new variant doesn't change existing ones. See [`SNAPSHOTS`].
struct SnapshotPreset {
    /// CLI name (e.g. `combined`, `bare`). Used as `--snapshot=NAME`.
    name: &'static str,
    /// Short human-readable description for the eprintln banner.
    desc: &'static str,
    /// Treatment modality for this preset.
    treatment: Treatment,
    /// Display name for `treatment` (avoids a Debug-format dependency).
    treatment_name: &'static str,
    /// True if immune kills + DAMP coupling are on (immune_mode = "immune_on").
    immune_on: bool,
    /// True if stromal protection (CAF GSH/MUFA shielding) is on.
    stromal_on: bool,
    /// True if the pH gradient + iron release + ion-trap are on.
    ph_on: bool,
    /// True if this preset uses the multi-dose schedule (#239) instead of
    /// the steady-state `Constant` default. The concrete `DoseSchedule` is
    /// built at runtime in `run_snapshot` (a `const` can't hold the
    /// `MultiDose` Vec).
    multidose: bool,
    /// True if the drug-tolerant persister model (#241) is enabled. Adds a
    /// `trajectory_persister.npy` overlay and emits `persister_mean`.
    persister: bool,
    /// True if the persister OXPHOS-ROS suppression axis (#470) is enabled on
    /// top of the persister model: a persister's basal/mitochondrial ROS is
    /// scaled down so a GPX4 inhibitor kills it less. Requires `persister` to
    /// also be true (it modulates the persister fraction). No extra overlay; the
    /// effect shows in the reduced RSL3 kill of the persister population.
    persister_oxphos: bool,
    /// True if clonal heterogeneity (#242) is enabled. Writes a static
    /// `subclone.npy` (u8, no time axis) for the renderer's subclone panel and
    /// emits per-subclone kill stats.
    clonal: bool,
    /// True if the explicit vessel model (#191) is enabled. Writes a static
    /// `vessel_supply.npy` (f32, no time axis) for the renderer's O2-supply
    /// panel and emits `vascular_hypoxic_fraction`.
    vasculature: bool,
    /// True if the 3D spheroid radial biology (#197) is enabled (radial
    /// phenotypes + GSH/iron/MUFA gradients, `Params::spheroid()`). Writes a
    /// static `phenotype.npy` (u8) for the renderer's phenotype panel.
    spheroid: bool,
    /// True if patient-scale slab mode (#240) is enabled. No extra overlay —
    /// the depth gradient is visible directly in the z-axis (layer) mid-slice
    /// of the dead/DAMP/LP panels.
    slab: bool,
    /// True if the Treg/MDSC suppressor field (#264) is enabled. Writes a
    /// static `suppressor.npy` (u8) source-niche mask for the renderer panel.
    suppressor: bool,
    /// True if the multi-checkpoint dual blockade (#264 Phase 3) is enabled
    /// (anti-PD-1 + anti-CTLA-4). No new overlay — the enhanced immune killing
    /// shows directly in the dead/DAMP panels.
    checkpoints: bool,
    /// True if cell-cell contact resistance (#270) is enabled. No extra static
    /// overlay — the contact-driven radial survival gradient (dense interior
    /// resists, sparse surface dies) shows directly in the dead/LP panels.
    contact: bool,
    /// True if the radial nutrient gradient (#270 item 3b) is enabled. No extra
    /// static overlay; the nutrient-starved core's reduced antioxidant capacity
    /// shows as more core killing in the dead/LP panels.
    nutrient: bool,
    /// True if the cDC1/cDC2 dendritic-cell subset mix (#264 Phase 4) is enabled.
    /// A cDC1-poor tumor primes anti-tumor killing less efficiently; no extra
    /// overlay, the reduced immune killing shows in the dead/DAMP panels.
    dc_subsets: bool,
    /// True if DC ferroptosis susceptibility (#469) is enabled: a ferroptotic
    /// TME kills the effector DCs themselves (PD-L1-gated), lowering immune
    /// amplification. No extra static overlay; the effect shows in the immune
    /// kill counts (run with the immune layer on).
    dc_ferroptosis: bool,
    /// True if the therapy-induced-senescence program (#341) is enabled. No
    /// extra static overlay; the senescent cells' shifted ferroptosis response
    /// (resistant or senolytic per the axis mix) plus the SASP immune coupling
    /// show directly in the dead/LP panels vs the immune-on baseline.
    senescence: bool,
    /// True if phenotype-specific SCD1/MUFA accumulation rates (#363) are enabled
    /// (`PhenotypeMufaConfig::literature()`). Runs in the spheroid context (the
    /// only MUFA-active path), so the spheroid `phenotype.npy` panel already shows
    /// which radial phenotype gets which rate; no extra static overlay. NOTE: the
    /// rate's effect is on the MUFA timecourse, not the kill count (the spheroid
    /// is cap-limited), so this preset's value is making the layer CLI-reachable.
    phenotype_mufa: bool,
    /// O2-dependent (Type II) SDT exo-ROS fraction (#358/#380). `0.0` in every
    /// matrix-adjacent preset (O2-independent, byte-identical); `1.0` in the
    /// `sdt-o2dep` preset makes SDT fully O2-gated, so the deep hypoxic core
    /// SURVIVES — the dead/LP panels show the hypoxic-core SDT survival front.
    sdt_o2_dependence: f64,
    /// True if the NCOA4-ferritinophagy + hypoxia-iron coupling (#365/#381) is
    /// enabled (`ferritinophagy_release` + `hypoxia_iron_sensitivity`). `false`
    /// in every existing preset (byte-identical); `true` in the `ferritinophagy`
    /// preset, so RSL3 on a hypoxic sphere shows the iron-amplified ferroptosis
    /// front (more kill where the boosted Fenton iron meets residual ROS).
    ferritinophagy: bool,
    /// True if the IFN-γ → System Xc⁻ + ACSL4 ferroptosis-sensitization coupling
    /// (#443) is enabled (`IFNGammaConfig::literature()`, both arms). No extra
    /// static overlay: the IFN-γ field seeds from the run's evolving DAMP signal
    /// (deaths during the run, not a static mask), so it is not reconstructable
    /// the way the `sasp-field` overlay is — the immune-primed sensitization
    /// shows directly as more kill in the high-immune-activity zones of the
    /// dead/LP panels vs the immune-on baseline.
    ifngamma: bool,
    /// True if the ALOX isoform-specific peroxidation + MCFA sensitization (#446)
    /// is enabled (`AloxConfig::literature()`: ALOX15-high + MCFA). No extra
    /// static overlay; the faster enzymatic propagation + extra oxidizable PUFA
    /// show as more kill in the dead/LP panels vs the baseline.
    alox: bool,
    /// True if the ACSL4-negative biomarker stratification (#444) is enabled
    /// (acsl4_status = ACSL4_NEGATIVE). The collapsed PUFA substrate makes the
    /// tumor ferroptosis-REFRACTORY, so RSL3 kills LESS than the (ACSL4-normal)
    /// baseline — the dead/LP panels show the survival of an ACSL4-negative tumor
    /// (HCC/AML-like) under the same RSL3 dose. No extra static overlay.
    acsl4_negative: bool,
    /// True if the ESCRT-III membrane-repair brake (#465) is enabled (a high
    /// repair rate + ample budget). A cell whose LP crosses the death threshold is
    /// resealed for a finite budget, so RSL3 kills LESS / later than the baseline
    /// (membrane repair delays death execution). No extra static overlay; the
    /// reduced/slower death front shows in the dead/LP panels vs the baseline.
    escrt: bool,
    /// True if the POR/CYB5R1 enzymatic O2-coupled H2O2 source (#466) is enabled.
    /// Injects an O2-scaled enzymatic oxidant into each tumor cell's basal ROS, so
    /// RSL3 kills MORE in the oxygenated rim and LESS in the hypoxic core (the
    /// O2-coupling keeps it from amplifying the deep-core artifact). No overlay; the
    /// rim-weighted extra death shows in the dead/LP panels.
    por: bool,
    /// True if the 7-DHC sterol radical-trapping defense (#467) is enabled (a high
    /// 7-DHC pool, i.e. a DHCR7-low tumor). The extra GPX4-independent quench lowers
    /// the propagation rate, so RSL3 kills LESS than the baseline (ferroptosis
    /// resistance). No overlay; the reduced death front shows in the dead/LP panels.
    dhc7: bool,
    /// True if the vitamin-K / VKORC1L1 radical-trapping defense (#483) is enabled
    /// (a VKORC1L1-high, p53-competent tumor). The extra GPX4-independent quench
    /// lowers the propagation rate, so RSL3 kills LESS than the baseline
    /// (ferroptosis resistance); a warfarin inhibitor reverses it (shown in the
    /// A/B test). No overlay; the reduced death front shows in the dead/LP panels.
    vitk: bool,
    /// True if PROM2 / MVB-exosome labile-iron efflux (#484) is enabled (a
    /// PROM2-high tumor that exports labile iron, starving the Fenton reaction).
    /// RSL3 kills LESS than the baseline (ferroptosis resistance). No overlay; the
    /// reduced death front shows in the dead/LP panels.
    prom2: bool,
    /// True if copper-ionophore / cuproptosis crosstalk (#485) is enabled (an
    /// elesclomol-like ionophore in an ATP7B-low tumor): intracellular copper
    /// overload depletes GSH + GPX4 each step, so RSL3 kills MORE than the
    /// baseline. No overlay; the increased death front shows in the dead/LP panels.
    copper: bool,
    /// True if dietary-PUFA / DGAT lipid-droplet buffering (#486) is enabled
    /// (exogenous PUFA above the saturable buffer, with DGAT inhibition). Adds
    /// oxidizable substrate, so RSL3 kills MORE than the baseline. No overlay; the
    /// increased death front shows in the dead/LP panels.
    dietary_pufa: bool,
}

/// Visualization presets for `--snapshot=NAME`. Keep this list small —
/// each entry costs ~333 MB on disk when its trajectory is generated.
const SNAPSHOTS: &[SnapshotPreset] = &[
    SnapshotPreset {
        name: "combined",
        desc: "RSL3 + immune + stromal + pH (all TME protections active)",
        treatment: Treatment::RSL3,
        treatment_name: "RSL3",
        immune_on: true,
        stromal_on: true,
        ph_on: true,
        multidose: false,
        persister: false,
        persister_oxphos: false,
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: false,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.0,
        ferritinophagy: false,
        ifngamma: false,
        alox: false,
        acsl4_negative: false,
        escrt: false,
        por: false,
        dhc7: false,
        vitk: false,
        prom2: false,
        copper: false,
        dietary_pufa: false,
    },
    SnapshotPreset {
        name: "bare",
        desc: "RSL3, no immune / stromal / pH (the unprotected baseline)",
        treatment: Treatment::RSL3,
        treatment_name: "RSL3",
        immune_on: false,
        stromal_on: false,
        ph_on: false,
        multidose: false,
        persister: false,
        persister_oxphos: false,
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: false,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.0,
        ferritinophagy: false,
        ifngamma: false,
        alox: false,
        acsl4_negative: false,
        escrt: false,
        por: false,
        dhc7: false,
        vitk: false,
        prom2: false,
        copper: false,
        dietary_pufa: false,
    },
    SnapshotPreset {
        name: "multidose",
        desc: "SDT multi-dose (4 pulses) + immune — death waves sync to each dose (#239)",
        treatment: Treatment::SDT,
        treatment_name: "SDT",
        immune_on: true,
        stromal_on: false,
        ph_on: false,
        multidose: true,
        persister: false,
        persister_oxphos: false,
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: false,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.0,
        ferritinophagy: false,
        ifngamma: false,
        alox: false,
        acsl4_negative: false,
        escrt: false,
        por: false,
        dhc7: false,
        vitk: false,
        prom2: false,
        copper: false,
        dietary_pufa: false,
    },
    SnapshotPreset {
        // SDT here visualizes the persister-fraction OVERLAY (the MUFA axis +
        // acquire/revert dynamics). SDT has no covalent GPX4 step, so the
        // GPX4-resistance kill-reduction does NOT apply here — that effect is
        // RSL3-specific (demonstrated by `persister_reduces_multidose_kills`).
        name: "persister",
        desc: "SDT multi-dose + immune — persister-fraction OVERLAY (accumulation/reversion; the RSL3 kill-drop is in the test suite) (#241)",
        treatment: Treatment::SDT,
        treatment_name: "SDT",
        immune_on: true,
        stromal_on: false,
        ph_on: false,
        multidose: true,
        persister: true,
        persister_oxphos: false,
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: false,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.0,
        ferritinophagy: false,
        ifngamma: false,
        alox: false,
        acsl4_negative: false,
        escrt: false,
        por: false,
        dhc7: false,
        vitk: false,
        prom2: false,
        copper: false,
        dietary_pufa: false,
    },
    SnapshotPreset {
        // RSL3 (covalent GPX4 inhibitor) on a persister population WITH the
        // OXPHOS-ROS suppression axis (#470): the persisters downregulate OXPHOS,
        // lowering the basal/mitochondrial ROS the GPX4 inhibitor needs, so RSL3
        // kills them LESS than the persister-without-OXPHOS baseline (an HDAC
        // inhibitor would re-raise the ROS, PMID 40909720). RSL3 (not SDT) so the
        // ROS-supply kill-reduction is on the covalent-knockdown path it actually
        // governs. No extra overlay; the reduced RSL3 death front IS the result.
        name: "persister-oxphos",
        desc: "RSL3 + persister with OXPHOS-ROS suppression (#470): OXPHOS-low persisters resist RSL3, HDACi reverses",
        treatment: Treatment::RSL3,
        treatment_name: "RSL3",
        immune_on: false,
        stromal_on: false,
        ph_on: false,
        multidose: true,
        persister: true,
        persister_oxphos: true,
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: false,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.0,
        ferritinophagy: false,
        ifngamma: false,
        alox: false,
        acsl4_negative: false,
        escrt: false,
        por: false,
        dhc7: false,
        vitk: false,
        prom2: false,
        copper: false,
        dietary_pufa: false,
    },
    SnapshotPreset {
        name: "clonal",
        desc: "SDT multi-dose + immune + clonal (4 subclones) — static subclone-id overlay (#242)",
        treatment: Treatment::SDT,
        treatment_name: "SDT",
        immune_on: true,
        stromal_on: false,
        ph_on: false,
        multidose: true,
        persister: false,
        persister_oxphos: false,
        clonal: true,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: false,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.0,
        ferritinophagy: false,
        ifngamma: false,
        alox: false,
        acsl4_negative: false,
        escrt: false,
        por: false,
        dhc7: false,
        vitk: false,
        prom2: false,
        copper: false,
        dietary_pufa: false,
    },
    SnapshotPreset {
        // RSL3 (hypoxia-sensitive) + explicit internal vessels: near-vessel
        // cells stay oxygenated and die; cells in the inter-vessel gaps go
        // hypoxic and survive — the irregular, non-radial kill pattern the
        // edge-distance proxy can't produce. immune/stromal/pH off to isolate
        // the O2→kill effect.
        name: "vasculature",
        desc: "RSL3 + explicit vessels (#191) — patchy O2 from internal vessels drives heterogeneous kill",
        treatment: Treatment::RSL3,
        treatment_name: "RSL3",
        immune_on: false,
        stromal_on: false,
        ph_on: false,
        multidose: false,
        persister: false,
        persister_oxphos: false,
        clonal: false,
        vasculature: true,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: false,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.0,
        ferritinophagy: false,
        ifngamma: false,
        alox: false,
        acsl4_negative: false,
        escrt: false,
        por: false,
        dhc7: false,
        vitk: false,
        prom2: false,
        copper: false,
        dietary_pufa: false,
    },
    SnapshotPreset {
        // RSL3 + explicit vessels, but the per-cell O2/drug supply is the
        // steady-state REACTION-DIFFUSION field (#343 PR 2) instead of the
        // monotonic nearest-vessel proxy. Same vessels, same λ — the
        // vessel_supply.npy panel shows the non-monotonic inter-vessel pockets
        // (well-supplied between several vessels, starved in avascular gaps)
        // that the proxy averages away. The RD path is name-keyed in
        // run_snapshot (like sasp-field), so vasculature stays `true` here.
        name: "reaction-diffusion",
        desc: "RSL3 + explicit vessels with the reaction-diffusion supply field (#343) — non-monotonic O2 the exp proxy misses",
        treatment: Treatment::RSL3,
        treatment_name: "RSL3",
        immune_on: false,
        stromal_on: false,
        ph_on: false,
        multidose: false,
        persister: false,
        persister_oxphos: false,
        clonal: false,
        vasculature: true,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: false,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.0,
        ferritinophagy: false,
        ifngamma: false,
        alox: false,
        acsl4_negative: false,
        escrt: false,
        por: false,
        dhc7: false,
        vitk: false,
        prom2: false,
        copper: false,
        dietary_pufa: false,
    },
    SnapshotPreset {
        // SDT + radial spheroid biology: the phenotype panel shows the
        // glycolytic rim / OXPHOS mid / persister core structure, and the
        // GSH/iron/MUFA gradients shape where SDT kills.
        name: "spheroid",
        desc: "SDT + 3D spheroid radial biochemistry (#197) — radial phenotype + GSH/iron/MUFA gradients",
        treatment: Treatment::SDT,
        treatment_name: "SDT",
        immune_on: true,
        stromal_on: false,
        ph_on: false,
        multidose: false,
        persister: false,
        persister_oxphos: false,
        clonal: false,
        vasculature: false,
        spheroid: true,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: false,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.0,
        ferritinophagy: false,
        ifngamma: false,
        alox: false,
        acsl4_negative: false,
        escrt: false,
        por: false,
        dhc7: false,
        vitk: false,
        prom2: false,
        copper: false,
        dietary_pufa: false,
    },
    SnapshotPreset {
        // SDT on a patient-scale slab at the SURFACE (+z face = vessel, depth
        // offset 0) so the supply gradient is fully visible across the block:
        // the top (+z) layers are well-perfused and die, deeper (−z) layers go
        // drug/O2-deprived and survive — the depth-dependent penetration the
        // in-vitro spheroid scale misses (#240). The death front in the z-axis
        // mid-slice IS the visualization (no extra static overlay needed). A
        // deep `patient_deep()` slab would kill ~nothing, so the surface slab is
        // the illustrative choice; the depth comparison lives in the tests.
        name: "slab",
        desc: "SDT on a patient-scale surface slab (#240) — depth-graded supply; +z dies, −z survives",
        treatment: Treatment::SDT,
        treatment_name: "SDT",
        immune_on: true,
        stromal_on: false,
        ph_on: false,
        multidose: false,
        persister: false,
        persister_oxphos: false,
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: true,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: false,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.0,
        ferritinophagy: false,
        ifngamma: false,
        alox: false,
        acsl4_negative: false,
        escrt: false,
        por: false,
        dhc7: false,
        vitk: false,
        prom2: false,
        copper: false,
        dietary_pufa: false,
    },
    SnapshotPreset {
        // Slab + internal vessels (#272 coupling). vessel_supply.npy (on a slab
        // grid) shows the combined planar-MAX-vessel field: focal well-perfused
        // pockets at depth around internal vessels, against the otherwise
        // monotonic depth collapse. Because supply scales drug/O2 DELIVERY, the
        // dead/LP panels show extra killing localized to those pockets — the
        // therapy reaches deep tissue near a vessel that the drug-starved bulk
        // escapes. So internal vasculature makes a patient-scale slab LESS
        // therapy-resistant at depth than the planar-only model (#240) implies.
        name: "slab-vessels",
        desc: "SDT on a patient-scale slab with internal vessels (#272) — drug reaches focal deep pockets",
        treatment: Treatment::SDT,
        treatment_name: "SDT",
        immune_on: true,
        stromal_on: false,
        ph_on: false,
        multidose: false,
        persister: false,
        persister_oxphos: false,
        clonal: false,
        vasculature: true,
        spheroid: false,
        slab: true,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: false,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.0,
        ferritinophagy: false,
        ifngamma: false,
        alox: false,
        acsl4_negative: false,
        escrt: false,
        por: false,
        dhc7: false,
        vitk: false,
        prom2: false,
        copper: false,
        dietary_pufa: false,
    },
    SnapshotPreset {
        // SDT + immune + Treg/MDSC suppressor (#264 Phase 2). Heuristic niche
        // patches (no vasculature in this preset) locally dampen immune kill:
        // the static `suppressor.npy` mask shows the niches, and the dead panel
        // shows immune killing suppressed in their neighborhoods.
        name: "suppressor",
        desc: "SDT + immune + Treg/MDSC suppressor field (#264) — niches locally dampen immune kill",
        treatment: Treatment::SDT,
        treatment_name: "SDT",
        immune_on: true,
        stromal_on: false,
        ph_on: false,
        multidose: false,
        persister: false,
        persister_oxphos: false,
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: true,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: false,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.0,
        ferritinophagy: false,
        ifngamma: false,
        alox: false,
        acsl4_negative: false,
        escrt: false,
        por: false,
        dhc7: false,
        vitk: false,
        prom2: false,
        copper: false,
        dietary_pufa: false,
    },
    SnapshotPreset {
        // SDT + immune + dual checkpoint blockade (#264 Phase 3): a PD-1 +
        // CTLA-4 tumor with BOTH inhibitors applied, so the combined brake is
        // low and immune killing is aggressive. No new overlay — the enhanced
        // death front shows directly in the dead/DAMP panels (contrast with the
        // `multidose`/`combined` presets' single-PD-1 brake).
        name: "checkpoint",
        desc: "SDT + immune + dual checkpoint blockade (#264) — anti-PD-1 + anti-CTLA-4 lifts both brakes",
        treatment: Treatment::SDT,
        treatment_name: "SDT",
        immune_on: true,
        stromal_on: false,
        ph_on: false,
        multidose: false,
        persister: false,
        persister_oxphos: false,
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: true,
        contact: false,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: false,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.0,
        ferritinophagy: false,
        ifngamma: false,
        alox: false,
        acsl4_negative: false,
        escrt: false,
        por: false,
        dhc7: false,
        vitk: false,
        prom2: false,
        copper: false,
        dietary_pufa: false,
    },
    SnapshotPreset {
        // Kitchen-sink composition (#278): several realism layers at once —
        // persister (#241) + clonal subclones (#242) + Treg/MDSC suppressor
        // (#264 P2) + dual checkpoint blockade (#264 P3), SDT multi-dose +
        // immune. Writes the persister + subclone + suppressor overlays
        // together (all geometry/seed-only, so they match the run). Excludes
        // vasculature/spheroid/slab: those re-grid or change the O2 source,
        // which would desync the static overlays from the actual run.
        name: "combined-realism",
        desc: "SDT multidose + immune + persister + clonal + suppressor + checkpoints (#278) — composed realism layers",
        treatment: Treatment::SDT,
        treatment_name: "SDT",
        immune_on: true,
        stromal_on: false,
        ph_on: false,
        multidose: true,
        persister: true,
        persister_oxphos: false,
        clonal: true,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: true,
        checkpoints: true,
        contact: false,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: false,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.0,
        ferritinophagy: false,
        ifngamma: false,
        alox: false,
        acsl4_negative: false,
        escrt: false,
        por: false,
        dhc7: false,
        vitk: false,
        prom2: false,
        copper: false,
        dietary_pufa: false,
    },
    SnapshotPreset {
        // RSL3 + cell-cell contact resistance (#270): dense interior cells
        // (E-cadherin/NF2-YAP → ACSL4/TFRC down) resist PUFA/iron-dependent
        // RSL3, while the sparse surface shell stays sensitive. No extra static
        // overlay — the radial survival gradient in the dead/LP panels IS the
        // visualization. Runs on the centred sphere (no spheroid/slab) where the
        // fixed-26 contact denominator is correct.
        name: "contact",
        desc: "RSL3 + cell-cell contact resistance (#270) — dense interior resists, sparse surface dies",
        treatment: Treatment::RSL3,
        treatment_name: "RSL3",
        immune_on: false,
        stromal_on: false,
        ph_on: false,
        multidose: false,
        persister: false,
        persister_oxphos: false,
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: true,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: false,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.0,
        ferritinophagy: false,
        ifngamma: false,
        alox: false,
        acsl4_negative: false,
        escrt: false,
        por: false,
        dhc7: false,
        vitk: false,
        prom2: false,
        copper: false,
        dietary_pufa: false,
    },
    SnapshotPreset {
        // SDT + radial nutrient gradient (#270 item 3b): the nutrient-starved
        // core loses glucose-derived NADPH for GSH/GPX4 regeneration, so its
        // durable antioxidant setpoint drops and the core kills MORE under SDT.
        // No extra static overlay; the shifted death front in the dead/LP
        // panels IS the visualization. Runs on the centred sphere.
        name: "nutrient",
        desc: "SDT + radial nutrient gradient (#270): starved core loses antioxidant capacity",
        treatment: Treatment::SDT,
        treatment_name: "SDT",
        immune_on: true,
        stromal_on: false,
        ph_on: false,
        multidose: false,
        persister: false,
        persister_oxphos: false,
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: true,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: false,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.0,
        ferritinophagy: false,
        ifngamma: false,
        alox: false,
        acsl4_negative: false,
        escrt: false,
        por: false,
        dhc7: false,
        vitk: false,
        prom2: false,
        copper: false,
        dietary_pufa: false,
    },
    SnapshotPreset {
        // SDT + cDC1/cDC2 dendritic-cell subset mix (#264 Phase 4): a cDC1-poor
        // tumor (the literature default) primes anti-tumor CD8 killing less
        // efficiently, so the same DAMP signal yields LESS immune killing. No
        // extra static overlay; the reduced death front in the dead/DAMP panels
        // (vs the immune-on baseline) IS the visualization.
        name: "dc-subsets",
        desc: "SDT + cDC1/cDC2 DC subset mix (#264): cDC1-poor tumor primes killing less efficiently",
        treatment: Treatment::SDT,
        treatment_name: "SDT",
        immune_on: true,
        stromal_on: false,
        ph_on: false,
        multidose: false,
        persister: false,
        persister_oxphos: false,
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: true,
        dc_ferroptosis: false,
        senescence: false,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.0,
        ferritinophagy: false,
        ifngamma: false,
        alox: false,
        acsl4_negative: false,
        escrt: false,
        por: false,
        dhc7: false,
        vitk: false,
        prom2: false,
        copper: false,
        dietary_pufa: false,
    },
    SnapshotPreset {
        // SDT + DC ferroptosis susceptibility (#469): the strong ferroptotic
        // DAMP/lipid-ROS signal that drives DC activation ALSO kills the
        // effector DCs themselves (PD-L1-low, the literature() config), lowering
        // priming, so a ferroptosis-inducing TME suppresses its own immune
        // amplification (PMID 39423128). No extra static overlay; the reduced
        // immune death front (vs the immune-on baseline) IS the visualization.
        name: "dc-ferroptosis",
        desc: "SDT + DC ferroptosis susceptibility (#469): a ferroptotic TME kills PD-L1-low effector DCs, lowering amplification",
        treatment: Treatment::SDT,
        treatment_name: "SDT",
        immune_on: true,
        stromal_on: false,
        ph_on: false,
        multidose: false,
        persister: false,
        persister_oxphos: false,
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: true,
        senescence: false,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.0,
        ferritinophagy: false,
        ifngamma: false,
        alox: false,
        acsl4_negative: false,
        escrt: false,
        por: false,
        dhc7: false,
        vitk: false,
        prom2: false,
        copper: false,
        dietary_pufa: false,
    },
    SnapshotPreset {
        // SDT + therapy-induced senescence (#341): a fraction of tumor cells
        // enter the senescence program (raised iron + antioxidant/GPX4 defenses)
        // AND couple to the immune layer via SASP. With the literature() config
        // the defense axis dominates the raised iron (net ferroptosis-resistant)
        // while the SASP multiplier < 1 (immunosuppressive, the established-tumor
        // arm, Di Mitri 2014 / Eggert 2016) ALSO lowers their immune clearance,
        // so these cells both resist cell-intrinsic ferroptosis and evade the
        // immune layer, the durable escape route. The net sign is axis/therapy/
        // stage-dependent (the module test drives both); no extra static overlay,
        // the shifted dead/LP front vs the immune-on baseline IS the
        // visualization. Runs on the centred sphere.
        name: "senescence",
        desc: "SDT + therapy-induced senescence (#341): ferroptosis resist/senolytic + SASP immune coupling",
        treatment: Treatment::SDT,
        treatment_name: "SDT",
        immune_on: true,
        stromal_on: false,
        ph_on: false,
        multidose: false,
        persister: false,
        persister_oxphos: false,
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: true,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.0,
        ferritinophagy: false,
        ifngamma: false,
        alox: false,
        acsl4_negative: false,
        escrt: false,
        por: false,
        dhc7: false,
        vitk: false,
        prom2: false,
        copper: false,
        dietary_pufa: false,
    },
    SnapshotPreset {
        // Like `senescence`, but adds the diffusing-SASP-field overlay (#376/#398):
        // writes a static `sasp_field.npy` (f32) of the quasi-steady SASP field
        // (seeded at the senescent cells, diffused with the same operator/constants
        // the run uses), so the renderer can show the paracrine source-to-neighbour
        // gradient that reaches non-senescent cells. The overlay is keyed off the
        // preset NAME in run_snapshot (no new struct field), and the run is the
        // identical senescence/SASP config, so the production matrix is untouched.
        name: "sasp-field",
        desc: "SDT + senescence with the diffusing SASP field overlay (#376/#398): paracrine immune coupling",
        treatment: Treatment::SDT,
        treatment_name: "SDT",
        immune_on: true,
        stromal_on: false,
        ph_on: false,
        multidose: false,
        persister: false,
        persister_oxphos: false,
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: true,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.0,
        ferritinophagy: false,
        ifngamma: false,
        alox: false,
        acsl4_negative: false,
        escrt: false,
        por: false,
        dhc7: false,
        vitk: false,
        prom2: false,
        copper: false,
        dietary_pufa: false,
    },
    SnapshotPreset {
        // RSL3 + 3D spheroid (#197) + phenotype-specific SCD1/MUFA rates (#363):
        // each radial phenotype (glycolytic rim / OXPHOS mid / persister core)
        // builds MUFA protection at its own rate. Runs in the spheroid context
        // (the only MUFA-active path); the spheroid phenotype.npy panel shows
        // which phenotype gets which rate. NOTE: the rate's effect is on the MUFA
        // timecourse, not the kill count (the spheroid is cap-limited), so this
        // preset exists to make the layer CLI-reachable, not to show a kill shift.
        name: "phenotype-mufa",
        desc: "RSL3 + spheroid + phenotype-specific SCD1/MUFA rates (#363) — per-phenotype MUFA build-up",
        treatment: Treatment::RSL3,
        treatment_name: "RSL3",
        immune_on: false,
        stromal_on: false,
        ph_on: false,
        multidose: false,
        persister: false,
        persister_oxphos: false,
        clonal: false,
        vasculature: false,
        spheroid: true,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: false,
        phenotype_mufa: true,
        sdt_o2_dependence: 0.0,
        ferritinophagy: false,
        ifngamma: false,
        alox: false,
        acsl4_negative: false,
        escrt: false,
        por: false,
        dhc7: false,
        vitk: false,
        prom2: false,
        copper: false,
        dietary_pufa: false,
    },
    SnapshotPreset {
        // SDT on a hypoxic sphere (the base edge-distance radial-O2 gradient, not
        // the spheroid biology layer) with the exo-ROS made fully O2-dependent
        // (#358/#380, sdt_o2_dependence=1.0): the Type II singlet-oxygen yield
        // scales with local O2, so SDT kills the oxygenated rim but the deep
        // hypoxic core SURVIVES. The dead/LP panels show that hypoxic-core SDT
        // survival front (the §7.1 contested-leg lower bound made visible).
        name: "sdt-o2dep",
        desc: "SDT + O2-dependent exo-ROS (#358): hypoxic-core SDT survival front",
        treatment: Treatment::SDT,
        treatment_name: "SDT",
        immune_on: false,
        stromal_on: false,
        ph_on: false,
        multidose: false,
        persister: false,
        persister_oxphos: false,
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: false,
        phenotype_mufa: false,
        sdt_o2_dependence: 1.0,
        ferritinophagy: false,
        ifngamma: false,
        alox: false,
        acsl4_negative: false,
        escrt: false,
        por: false,
        dhc7: false,
        vitk: false,
        prom2: false,
        copper: false,
        dietary_pufa: false,
    },
    SnapshotPreset {
        // RSL3 on a hypoxic sphere with the NCOA4-ferritinophagy + hypoxia-iron
        // coupling on (#365/#381): HIF/TfR1 raise the labile-iron pool where O2
        // is low while ferritinophagy releases stored iron over the run, so the
        // boosted Fenton iron amplifies RSL3 ferroptosis. The dead/LP panels show
        // the iron-amplified ferroptosis front (cf. the §7.1 hypoxia-iron leg; the
        // deep-core rise is the flagged O2-independent-Fenton artifact, corrected
        // by the off-by-default #383 fenton_o2_dependence not enabled here).
        name: "ferritinophagy",
        desc: "RSL3 + NCOA4-ferritinophagy + hypoxia-iron (#365): iron-amplified ferroptosis front",
        treatment: Treatment::RSL3,
        treatment_name: "RSL3",
        immune_on: false,
        stromal_on: false,
        ph_on: false,
        multidose: false,
        persister: false,
        persister_oxphos: false,
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: false,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.0,
        ferritinophagy: true,
        ifngamma: false,
        alox: false,
        acsl4_negative: false,
        escrt: false,
        por: false,
        dhc7: false,
        vitk: false,
        prom2: false,
        copper: false,
        dietary_pufa: false,
    },
    SnapshotPreset {
        // RSL3 + immune with the IFN-γ → System Xc⁻ + ACSL4 coupling on (#443):
        // ferroptotic deaths release DAMPs → CD8 priming → an IFN-γ field that
        // diffuses back onto nearby tumor cells, suppressing System Xc⁻ (GSH down)
        // and raising ACSL4 (PUFA up), so the cells re-sensitize to ferroptosis.
        // The dead/LP panels show MORE kill in the high-immune-activity zones than
        // the immune-on baseline (the molecular return arm of the immune-
        // amplification loop, Wang 2019 PMID 31043744). RSL3 (not SDT) so the
        // sensitization is read against a pure-ferroptosis death front.
        name: "ifngamma",
        desc: "RSL3 + immune + IFN-γ→System Xc⁻/ACSL4 coupling (#443): immune-primed ferroptosis sensitization",
        treatment: Treatment::RSL3,
        treatment_name: "RSL3",
        immune_on: true,
        stromal_on: false,
        ph_on: false,
        multidose: false,
        persister: false,
        persister_oxphos: false,
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: false,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.0,
        ferritinophagy: false,
        ifngamma: true,
        alox: false,
        acsl4_negative: false,
        escrt: false,
        por: false,
        dhc7: false,
        vitk: false,
        prom2: false,
        copper: false,
        dietary_pufa: false,
    },
    SnapshotPreset {
        // RSL3 on a sphere with an ALOX15-high, MCFA-exposed phenotype (#446):
        // the lipoxygenase isoform mix raises the enzymatic peroxidation rate and
        // MCFA→ACSL4/CD36 raises the oxidizable-PUFA pool, so RSL3 kills MORE than
        // the balanced-ALOX baseline — a lipid-machinery sensitization axis
        // distinct from the GPX4/GSH/FSP1 defenses. No extra overlay; the faster
        // death front shows in the dead/LP panels.
        name: "alox",
        desc: "RSL3 + ALOX15-high/MCFA (#446): lipoxygenase + PUFA-incorporation ferroptosis sensitization",
        treatment: Treatment::RSL3,
        treatment_name: "RSL3",
        immune_on: false,
        stromal_on: false,
        ph_on: false,
        multidose: false,
        persister: false,
        persister_oxphos: false,
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: false,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.0,
        ferritinophagy: false,
        ifngamma: false,
        alox: true,
        acsl4_negative: false,
        escrt: false,
        por: false,
        dhc7: false,
        vitk: false,
        prom2: false,
        copper: false,
        dietary_pufa: false,
    },
    SnapshotPreset {
        // RSL3 on an ACSL4-NEGATIVE tumor (#444): with ACSL4 absent the membrane
        // cannot incorporate the PUFA ferroptosis requires, so the oxidizable
        // substrate collapses and the tumor is intrinsically REFRACTORY to RSL3
        // (an escape distinct from the GPX4/GSH/FSP1 defenses; e.g. some HCC/AML
        // subtypes, Doll 2017 PMID 27842070). The dead/LP panels show MUCH LESS
        // kill than an ACSL4-normal tumor under the same RSL3 dose — the
        // patient-stratification headline (ACSL4-high respond, ACSL4-negative do
        // not). No extra overlay.
        name: "acsl4-negative",
        desc: "RSL3 on an ACSL4-negative tumor (#444): collapsed PUFA substrate ⇒ ferroptosis-refractory",
        treatment: Treatment::RSL3,
        treatment_name: "RSL3",
        immune_on: false,
        stromal_on: false,
        ph_on: false,
        multidose: false,
        persister: false,
        persister_oxphos: false,
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: false,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.0,
        ferritinophagy: false,
        ifngamma: false,
        alox: false,
        acsl4_negative: true,
        escrt: false,
        por: false,
        dhc7: false,
        vitk: false,
        prom2: false,
        copper: false,
        dietary_pufa: false,
    },
    SnapshotPreset {
        // RSL3 on a tumor with a high ESCRT-III membrane-repair capacity (#465).
        // Cells whose lipid peroxide crosses the death threshold are resealed for a
        // finite per-cell budget, so death execution is delayed and RSL3 kills LESS
        // / later than the no-repair baseline (membrane repair, not a redox/lipid
        // defense, is the resistance mechanism here; Dai 2020 PMID 31761326). The
        // dead/LP panels show the slower, reduced death front. No extra overlay.
        name: "escrt",
        desc: "RSL3 + ESCRT-III membrane repair (#465): death-execution brake delays/blocks ferroptosis",
        treatment: Treatment::RSL3,
        treatment_name: "RSL3",
        immune_on: false,
        stromal_on: false,
        ph_on: false,
        multidose: false,
        persister: false,
        persister_oxphos: false,
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: false,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.0,
        ferritinophagy: false,
        ifngamma: false,
        alox: false,
        acsl4_negative: false,
        escrt: true,
        por: false,
        dhc7: false,
        vitk: false,
        prom2: false,
        copper: false,
        dietary_pufa: false,
    },
    SnapshotPreset {
        // RSL3 on a tumor with a high POR/CYB5R1 enzymatic O2-coupled H2O2 source
        // (#466). POR transfers electrons from NAD(P)H to O2 to make H2O2 (the
        // Fenton substrate), so the extra oxidant raises RSL3 ferroptosis kill,
        // concentrated in the OXYGENATED rim because the O2-coupling makes POR make
        // little H2O2 in the hypoxic core (Yan 2021 PMID 33860083). The dead/LP
        // panels show the rim-weighted extra death front. No overlay.
        name: "por",
        desc: "RSL3 + POR/CYB5R1 O2-coupled H2O2 source (#466): rim-weighted enzymatic ferroptosis boost",
        treatment: Treatment::RSL3,
        treatment_name: "RSL3",
        immune_on: false,
        stromal_on: false,
        ph_on: false,
        multidose: false,
        persister: false,
        persister_oxphos: false,
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: false,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.0,
        ferritinophagy: false,
        ifngamma: false,
        alox: false,
        acsl4_negative: false,
        escrt: false,
        por: true,
        dhc7: false,
        vitk: false,
        prom2: false,
        copper: false,
        dietary_pufa: false,
    },
    SnapshotPreset {
        // RSL3 on a DHCR7-low tumor with a high 7-DHC sterol radical-trapping pool
        // (#467). 7-DHC is a membrane radical-trapping antioxidant that gates the
        // autocatalytic peroxidation chain, so the extra GPX4-independent quench
        // makes RSL3 kill LESS than the baseline (ferroptosis resistance; the
        // DHCR7-loss escape, Freitas/Li Nature 2024 PMID 38297130). The dead/LP
        // panels show the reduced death front. No overlay.
        name: "dhc7",
        desc: "RSL3 + 7-DHC sterol radical trap (#467): DHCR7-low ferroptosis resistance",
        treatment: Treatment::RSL3,
        treatment_name: "RSL3",
        immune_on: false,
        stromal_on: false,
        ph_on: false,
        multidose: false,
        persister: false,
        persister_oxphos: false,
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: false,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.0,
        ferritinophagy: false,
        ifngamma: false,
        alox: false,
        acsl4_negative: false,
        escrt: false,
        por: false,
        dhc7: true,
        vitk: false,
        prom2: false,
        copper: false,
        dietary_pufa: false,
    },
    SnapshotPreset {
        // RSL3 on a VKORC1L1-high, p53-competent tumor with a vitamin-K
        // radical-trapping defense (#483). VKORC1L1 reduces vitamin K to a
        // radical-trapping antioxidant that quenches the autocatalytic
        // peroxidation chain INDEPENDENT of GPX4/GSH, so the extra quench makes
        // RSL3 kill LESS than the baseline (ferroptosis resistance); the FDA
        // anticoagulant warfarin inhibits VKORC1L1 and reverses it, restoring the
        // kill (Yang et al. Cell Metab 2023 PMID 37467745; the warfarin reversal
        // is exercised by the warfarin_reverses_vkorc1l1_resistance A/B test).
        // This preset shows the DEFENDED (no-warfarin) state; the dead/LP panels
        // show the reduced death front. No overlay.
        name: "vkorc1l1",
        desc: "RSL3 + VKORC1L1 vitamin-K radical trap (#483): p53-competent ferroptosis resistance (warfarin reverses)",
        treatment: Treatment::RSL3,
        treatment_name: "RSL3",
        immune_on: false,
        stromal_on: false,
        ph_on: false,
        multidose: false,
        persister: false,
        persister_oxphos: false,
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: false,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.0,
        ferritinophagy: false,
        ifngamma: false,
        alox: false,
        acsl4_negative: false,
        escrt: false,
        por: false,
        dhc7: false,
        vitk: true,
        prom2: false,
        copper: false,
        dietary_pufa: false,
    },
    SnapshotPreset {
        // RSL3 on a PROM2-high tumor with MVB-exosome labile-iron efflux (#484).
        // Pro-ferroptotic stress induces Prominin-2, which exports ferritin-bound
        // iron in secreted exosomes, DEPLETING the labile iron pool and starving
        // the Fenton reaction (the OPPOSITE sign to ferritinophagy #340), so RSL3
        // kills LESS than the baseline (ferroptosis resistance; the EMT/metastatic
        // escape, Brown et al. Dev Cell 2019 PMID 31761539). The dead/LP panels
        // show the reduced death front. No overlay.
        name: "prom2",
        desc: "RSL3 + PROM2 iron efflux (#484): MVB-exosome iron export starves Fenton, ferroptosis resistance",
        treatment: Treatment::RSL3,
        treatment_name: "RSL3",
        immune_on: false,
        stromal_on: false,
        ph_on: false,
        multidose: false,
        persister: false,
        persister_oxphos: false,
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: false,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.0,
        ferritinophagy: false,
        ifngamma: false,
        alox: false,
        acsl4_negative: false,
        escrt: false,
        por: false,
        dhc7: false,
        vitk: false,
        prom2: true,
        copper: false,
        dietary_pufa: false,
    },
    SnapshotPreset {
        // RSL3 + a copper ionophore (elesclomol/disulfiram) on an ATP7B-low tumor
        // (#485). Copper overload depletes GSH (copper binds glutathione) and
        // drives GPX4 degradation each step, simultaneously enabling ferroptosis
        // and cuproptosis, so RSL3 kills MORE than the baseline (the FIN +
        // ionophore synergy, Gao et al. Mol Oncol 2021 PMID 34390123, elesclomol degrades ATP7A); an
        // ATP7B-efflux-competent tumor would resist (the A/B test). The dead/LP
        // panels show the increased death front. No overlay.
        name: "copper",
        desc: "RSL3 + copper ionophore (#485): elesclomol depletes GSH/GPX4 (cuproptosis crosstalk), raises kill",
        treatment: Treatment::RSL3,
        treatment_name: "RSL3",
        immune_on: false,
        stromal_on: false,
        ph_on: false,
        multidose: false,
        persister: false,
        persister_oxphos: false,
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: false,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.0,
        ferritinophagy: false,
        ifngamma: false,
        alox: false,
        acsl4_negative: false,
        escrt: false,
        por: false,
        dhc7: false,
        vitk: false,
        prom2: false,
        copper: true,
        dietary_pufa: false,
    },
    SnapshotPreset {
        // RSL3 + dietary-PUFA load with DGAT inhibition on an acidic tumor (#486).
        // Exogenous polyunsaturated fatty acids add oxidizable substrate once the
        // saturable lipid-droplet (DGAT) storage sink is exceeded, so RSL3 kills
        // MORE than the baseline; the effect is potentiated by tumor acidosis
        // (ph_on, so the existing pH layer lowers defenses in acidic zones and
        // composes with the added substrate, the Dierge et al. Cell Metab 2021
        // PMID 34118189 mechanism). The dead/LP panels show the increased death
        // front. No overlay.
        name: "dietary-pufa",
        desc: "RSL3 + dietary PUFA over the DGAT buffer on an acidic tumor (#486): more oxidizable substrate, more kill",
        treatment: Treatment::RSL3,
        treatment_name: "RSL3",
        immune_on: false,
        stromal_on: false,
        ph_on: true,
        multidose: false,
        persister: false,
        persister_oxphos: false,
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: false,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.0,
        ferritinophagy: false,
        ifngamma: false,
        alox: false,
        acsl4_negative: false,
        escrt: false,
        por: false,
        dhc7: false,
        vitk: false,
        prom2: false,
        copper: false,
        dietary_pufa: true,
    },
    SnapshotPreset {
        // SDT with a Type-I-heavy sonosensitizer (#468): sdt_o2_dependence = 0.3
        // means 70% Type I (oxygen-INDEPENDENT, hydroxyl/superoxide radicals) + 30%
        // Type II (oxygen-dependent singlet oxygen), so the exo-ROS yield retains
        // ~70% in the anoxic core where the pure Type II `sdt-o2dep` preset (dep =
        // 1.0) collapses. This is the hypoxia-tolerant SDT radical arm: the dead/LP
        // panels keep killing into the hypoxic core, the complement of `sdt-o2dep`.
        name: "sdt-typei",
        desc: "SDT + Type-I-heavy sonosensitizer (#468): O2-independent radical arm retains hypoxic-core kill",
        treatment: Treatment::SDT,
        treatment_name: "SDT",
        immune_on: false,
        stromal_on: false,
        ph_on: false,
        multidose: false,
        persister: false,
        persister_oxphos: false,
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
        dc_ferroptosis: false,
        senescence: false,
        phenotype_mufa: false,
        sdt_o2_dependence: 0.3,
        ferritinophagy: false,
        ifngamma: false,
        alox: false,
        acsl4_negative: false,
        escrt: false,
        por: false,
        dhc7: false,
        vitk: false,
        prom2: false,
        copper: false,
        dietary_pufa: false,
    },
];

/// The multi-dose schedule used by the `multidose` snapshot preset:
/// four ROS pulses across the 180-step run. Sharp half-life (8 steps) so
/// each pulse rises and fades distinctly, producing visible death waves.
///
/// These parameters are **illustrative** — chosen for visual clarity of
/// the death-wave dynamics, NOT calibrated to a clinical SDT protocol.
fn multidose_snapshot_schedule() -> DoseSchedule {
    DoseSchedule::MultiDose {
        dose_steps: vec![10, 55, 100, 145],
        peak: 1.0,
        half_life_steps: 8.0,
    }
}

/// Look up a snapshot preset by name, or print available choices and exit.
fn resolve_snapshot(name: &str) -> &'static SnapshotPreset {
    SNAPSHOTS
        .iter()
        .find(|s| s.name == name)
        .unwrap_or_else(|| {
            eprintln!("ERROR: unknown --snapshot name `{name}`. Available presets:");
            for s in SNAPSHOTS {
                eprintln!("  {:<10} — {}", s.name, s.desc);
            }
            std::process::exit(2);
        })
}

/// Run a single condition with per-step trajectory capture, then write
/// `trajectory_{dead,damp,lp}.npy` + `trajectory_meta.json` to
/// `output_dir`. Driven by the `--snapshot[=NAME]` CLI flag (#193, #239).
///
/// Map a `Phenotype` to a small int for the static phenotype-map overlay
/// (#197): `0` = stroma, `1..=4` = tumor phenotypes (rim→core ordering).
fn phenotype_to_u8(p: Phenotype) -> u8 {
    match p {
        Phenotype::Stromal => 0,
        Phenotype::Glycolytic => 1,
        Phenotype::OXPHOS => 2,
        Phenotype::Persister => 3,
        Phenotype::PersisterNrf2 => 4,
    }
}

/// Treatment, TME toggles, and dose schedule vary by [`SnapshotPreset`].
/// Files are written under the same names regardless of preset — running
/// a second time overwrites; rename manually to keep several side by side.
fn run_snapshot(output_dir: &Path, tumor_radius_um: f64, name: &str) {
    let preset = resolve_snapshot(name);
    let dose_schedule = if preset.multidose {
        multidose_snapshot_schedule()
    } else {
        DoseSchedule::Constant
    };
    // Administration steps for metadata / renderer dose markers (empty for
    // the steady-state Constant presets).
    let dose_steps = dose_schedule.dose_steps();

    let condition = Condition {
        name: format!("snapshot_{}_{}", preset.name, preset.treatment_name),
        treatment: preset.treatment,
        treatment_name: preset.treatment_name.to_string(),
        o2_lambda: Some(ZONE_REF_LAMBDA),
        immune_on: preset.immune_on,
        stromal_on: preset.stromal_on,
        ph_on: preset.ph_on,
        dose_schedule,
    };

    eprintln!(
        "=== --snapshot={}: {} ({}³ × {} steps) ===",
        preset.name, preset.desc, GRID_DIM, N_STEPS,
    );

    let run_cfg = RunConfig::production();
    // Enable the realism layers a preset opts into (#241 persister, #242
    // clonal, #191 vasculature).
    let clonal_cfg = preset.clonal.then(ClonalConfig::literature_4);
    let vasculature_cfg = preset
        .vasculature
        .then(VasculatureConfig::well_vascularized);
    // Reaction-diffusion supply (#343 PR 2): name-keyed like `sasp-field`, so no
    // per-preset struct field. The preset also sets `vasculature: true`, so the
    // vessels exist as the RD sources. `false` for every other preset ⇒ the
    // proxy runs ⇒ unchanged.
    let reaction_diffusion = preset.name == "reaction-diffusion";
    let spheroid_cfg = preset.spheroid.then(SpheroidConfig::literature);
    // Patient-scale slab (#240): surface slab so the supply gradient is visible
    // across the block in the z-axis mid-slice (no extra static overlay; the
    // death front in the dead/DAMP/LP panels IS the visualization).
    let slab_cfg = preset.slab.then(SlabConfig::surface);
    // Depth-graded slab phenotype (#272): the surface slab gets a layered
    // rim→core phenotype gradient (proliferating +z, persister-like deep −z)
    // matching its depth-graded supply, so the dead/LP panels show both the
    // supply-driven death front AND the intrinsic deep-cell tolerance.
    let slab_phenotype_cfg = preset.slab.then(SlabPhenotypeConfig::literature);
    let suppressor_cfg = preset.suppressor.then(SuppressorConfig::enabled);
    // Dual checkpoint blockade (#264 Phase 3): a PD-1 + CTLA-4 tumor with both
    // inhibitors applied (combined brake low ⇒ aggressive immune killing).
    let checkpoint_cfg = preset.checkpoints.then(|| {
        CheckpointPanel::pd1_ctla4_tumor()
            .with_anti_pd1(0.8)
            .with_anti_ctla4(0.8)
    });
    // Cell-cell contact resistance (#270): on the centred snapshot sphere the
    // fixed-26 contact denominator is correct (tumor never touches the box
    // face), so dense interior cells resist and the sparse surface dies — the
    // dead/LP panels show the radial survival gradient.
    let contact_cfg = preset.contact.then(ContactConfig::literature);
    let nutrient_cfg = preset.nutrient.then(NutrientConfig::literature);
    let dc_subsets_cfg = preset.dc_subsets.then(DcSubsetConfig::literature);
    // DC ferroptosis susceptibility (#469): the literature() config is a
    // PD-L1-low, ferroptosis-vulnerable DC population; `None` when the preset
    // does not request it ⇒ byte-identical.
    let dc_ferroptosis_cfg = preset.dc_ferroptosis.then(DcFerroptosisConfig::literature);
    // Therapy-induced senescence (#341): the literature() config (defense-
    // dominant ⇒ net ferroptosis-resistant under the trigger, plus an
    // immunosuppressive SASP multiplier < 1 ⇒ lowered immune clearance of the
    // same cells, the established-tumor escape arm).
    let senescence_cfg = preset.senescence.then(SenescenceConfig::literature);
    // Static viz overlays, recomputed from the same SEED + per-layer seed the
    // run uses internally, so they match the perturbations actually applied.
    // `None` unless the matching preset is active.
    let snapshot_grid = TumorGrid3D::generate(
        run_cfg.grid_dim,
        run_cfg.grid_dim,
        run_cfg.grid_dim,
        CELL_SIZE_UM,
        SEED,
    );
    let subclone_ids = clonal_cfg
        .as_ref()
        .map(|c| assign_subclones_3d(&snapshot_grid, c.k(), SUBCLONE_SEED));
    let vessel_supply = vasculature_cfg.map(|cfg| {
        if let Some(scfg) = slab_cfg {
            // Slab + vessels (#272): the overlay must match the run, which uses
            // a slab grid, slab-uniform vessels, and the combined planar-MAX-
            // vessel field. Regenerate on a slab grid so the panel shows the
            // focal well-perfused pockets the run actually applies.
            let g = TumorGrid3D::generate_slab(
                run_cfg.grid_dim,
                run_cfg.grid_dim,
                run_cfg.grid_dim,
                CELL_SIZE_UM,
                SEED,
            );
            let vessels = place_vessels_in_slab_3d(&g, &cfg, VESSEL_SEED);
            // Match the run path: RD field when reaction_diffusion is on (#343).
            let vessel = vessel_or_rd_supply(&g, &vessels, ZONE_REF_LAMBDA, reaction_diffusion);
            let planar = slab_supply_field(&g, scfg.depth_offset_mm * 1000.0, ZONE_REF_LAMBDA);
            combine_supply_max(&planar, &vessel)
        } else {
            let vessels = place_vessels_3d(&snapshot_grid, &cfg, VESSEL_SEED);
            vessel_or_rd_supply(
                &snapshot_grid,
                &vessels,
                ZONE_REF_LAMBDA,
                reaction_diffusion,
            )
        }
    });
    // Radial phenotype map (#197): re-run the same radial assignment, then map
    // each cell's Phenotype to a small int (0 = stroma, 1..=4 = tumor types).
    let phenotype_map = spheroid_cfg.map(|cfg| {
        let mut g = TumorGrid3D::generate(
            run_cfg.grid_dim,
            run_cfg.grid_dim,
            run_cfg.grid_dim,
            CELL_SIZE_UM,
            SEED,
        );
        apply_radial_cells_3d(&mut g, &cfg, SPHEROID_SEED);
        g.cells
            .iter()
            .map(|gc| phenotype_to_u8(gc.phenotype))
            .collect::<Vec<u8>>()
    });
    // Suppressor-source niche mask (#264): recompute the same mask the run uses
    // (heuristic here — this preset has no vasculature). 1 = Treg/MDSC source.
    let suppressor_mask = suppressor_cfg.map(|cfg| {
        suppressor_source_mask_3d(&snapshot_grid, &cfg, None, SUPPRESSOR_SEED)
            .iter()
            .map(|&s| u8::from(s))
            .collect::<Vec<u8>>()
    });
    // Diffusing SASP-field overlay (#376/#398), only for the `sasp-field` preset.
    // Reconstructs the quasi-steady SASP field deterministically from the same
    // senescence mask, diffusion operator, and constants the run uses (seed at the
    // senescent cells, diffuse over the immune window), so the renderer can show
    // the paracrine source-to-neighbour gradient. Computed on an identical
    // freshly-generated throwaway grid (apply_senescence_program_3d mutates and
    // TumorGrid3D is not Clone), so the run grid is untouched.
    let sasp_field_overlay: Option<Vec<f32>> = (preset.name == "sasp-field").then(|| {
        // Regenerate an identical throwaway grid (same args + SENESCENCE_SEED as
        // the run) so the mask matches the run's exactly; apply_senescence_program_3d
        // mutates, and TumorGrid3D is not Clone, so we do NOT touch snapshot_grid.
        let mut grid_clone = TumorGrid3D::generate(
            run_cfg.grid_dim,
            run_cfg.grid_dim,
            run_cfg.grid_dim,
            CELL_SIZE_UM,
            SEED,
        );
        let mask = apply_senescence_program_3d(
            &mut grid_clone,
            &SenescenceConfig::literature(),
            SENESCENCE_SEED,
        );
        let n = snapshot_grid.cells.len();
        let mut field = vec![0.0_f64; n];
        let mut scratch = vec![0.0_f64; n];
        // Match the run EXACTLY: in run_one_condition_full the SASP seed+diffuse
        // lives inside `if condition.immune_on`, which runs every step (it is the
        // immune KILL loop that is gated on `step >= IMMUNE_START_STEP`, not the
        // field update), so the field diffuses for the full run. Iterate the same
        // `0..n_steps` so the overlay reproduces the end-of-run field.
        for _ in 0..run_cfg.n_steps {
            for (idx, &is_sen) in mask.iter().enumerate() {
                if is_sen {
                    field[idx] = (field[idx] + SASP_FIELD_REPLENISH_RATE).min(1.0);
                }
            }
            diffuse_damp_3d_step(
                &mut field,
                &mut scratch,
                &snapshot_grid,
                SASP_FIELD_DIFFUSION_FRACTION,
                SASP_FIELD_CLEARANCE_RATE,
            );
        }
        field.iter().map(|&v| v as f32).collect()
    });
    let mut buffers =
        snapshot::SnapshotBuffers::new(run_cfg.grid_dim, run_cfg.n_steps, preset.persister);
    // Persister model (#241); the OXPHOS-ROS suppression axis (#470) is layered
    // on only when the preset requests it (`persister-oxphos`), at a placeholder
    // suppression so OXPHOS-low persisters resist RSL3. `persister == false` ⇒
    // `None` ⇒ no persister path ⇒ byte-identical.
    let persister_cfg = preset.persister.then(|| {
        let mut c = PersisterConfig::enabled();
        if preset.persister_oxphos {
            c.oxphos_ros_suppression = 0.7;
        }
        c
    });
    let result = run_one_condition_full(
        &condition,
        run_cfg,
        Some(&mut buffers),
        Overrides {
            persister: persister_cfg,
            clonal: clonal_cfg,
            vasculature: vasculature_cfg,
            spheroid: spheroid_cfg,
            slab: slab_cfg,
            slab_phenotype: slab_phenotype_cfg,
            suppressor: suppressor_cfg,
            checkpoints: checkpoint_cfg,
            contact: contact_cfg,
            nutrient: nutrient_cfg,
            dc_subsets: dc_subsets_cfg,
            dc_ferroptosis: dc_ferroptosis_cfg,
            senescence: senescence_cfg,
            phenotype_mufa: preset.phenotype_mufa.then(PhenotypeMufaConfig::literature),
            // #380: O2-dependent SDT exo-ROS (0.0 for every other preset ⇒ the
            // exo-ROS is O2-independent, unchanged).
            sdt_o2_dependence: preset.sdt_o2_dependence,
            // #381: NCOA4-ferritinophagy + hypoxia-iron coupling. The §7.1 headline
            // knob values (both 2.0); 0.0 ⇒ inert for every other preset.
            ferritinophagy_release: if preset.ferritinophagy { 2.0 } else { 0.0 },
            hypoxia_iron_sensitivity: if preset.ferritinophagy { 2.0 } else { 0.0 },
            // #343 PR 2: reaction-diffusion supply for the `reaction-diffusion`
            // preset (name-keyed above); `false` for every other preset.
            reaction_diffusion,
            // #443: IFN-γ → System Xc⁻ + ACSL4 coupling (both arms via literature()).
            // Gated on immune_on inside run_one_condition_full; `None` for every
            // other preset ⇒ no coupling.
            ifngamma: preset.ifngamma.then(IFNGammaConfig::literature),
            // #446: ALOX isoform-specific peroxidation + MCFA sensitization
            // (ALOX15-high + MCFA via literature()). `None` for every other
            // preset ⇒ both boosts 0 ⇒ unchanged.
            alox: preset.alox.then(AloxConfig::literature),
            // #444: ACSL4-negative biomarker stratification. `None` for every other
            // preset ⇒ no change; the acsl4-negative preset sets the null status so
            // the PUFA substrate collapses and the tumor resists RSL3.
            acsl4_status: preset.acsl4_negative.then_some(ACSL4_NEGATIVE),
            // #465: ESCRT-III membrane repair (high rate + ample budget) for the
            // `escrt` preset; `None` for every other preset ⇒ no brake.
            escrt: preset.escrt.then_some((0.6, 8.0)),
            // #466: POR/CYB5R1 enzymatic O2-coupled H2O2 source (rate 0.4, fully
            // O2-dependent so it is rim-weighted) for the `por` preset.
            por: preset.por.then_some((0.4, 1.0)),
            // #467: 7-DHC sterol radical-trapping pool for the `dhc7` preset (a
            // DHCR7-low resistant tumor); `None` for every other preset.
            dhc7: preset.dhc7.then_some(0.5),
            vitk: preset.vitk.then_some((1.0, 0.0)),
            prom2: preset.prom2.then_some(0.8),
            copper: preset.copper.then(CopperConfig::literature),
            dietary_pufa: preset.dietary_pufa.then_some((0.6, 0.2)),
            ..Default::default()
        },
    );
    eprintln!(
        "  done — total_kill={:.1}%, ferroptosis_kills={}, immune_kills={}, persister_mean={:.3}, subclones={}, captured {} steps",
        result.overall_kill_rate * 100.0,
        result.ferroptosis_kills.unwrap_or(0),
        result.immune_kills.unwrap_or(0),
        result.persister_mean.unwrap_or(0.0),
        result.subclone_kills.as_ref().map_or(0, |v| v.len()),
        buffers.steps_captured(),
    );

    buffers
        .write(output_dir)
        .expect("Failed to write trajectory .npy files");

    // Static subclone-id map (#242). Single 3D frame (no time axis), shape
    // (rows, cols, layers) to match the trajectory arrays' spatial axes so the
    // renderer can take the same mid-slice. Only when the preset is clonal.
    if let Some(ids) = &subclone_ids {
        let shape = [run_cfg.grid_dim, run_cfg.grid_dim, run_cfg.grid_dim];
        npy::write_u8_array(output_dir.join("subclone.npy"), &shape, ids)
            .expect("Failed to write subclone.npy");
    }

    // Static vessel O2-supply map (#191), same static 3D layout. Only when the
    // preset is vasculature.
    if let Some(supply) = &vessel_supply {
        let shape = [run_cfg.grid_dim, run_cfg.grid_dim, run_cfg.grid_dim];
        let supply_f32: Vec<f32> = supply.iter().map(|&v| v as f32).collect();
        npy::write_f32_array(output_dir.join("vessel_supply.npy"), &shape, &supply_f32)
            .expect("Failed to write vessel_supply.npy");
    }

    // Static radial phenotype map (#197), same static 3D layout. Only when the
    // preset is spheroid. 0 = stroma; 1..=4 = tumor phenotypes.
    if let Some(ph) = &phenotype_map {
        let shape = [run_cfg.grid_dim, run_cfg.grid_dim, run_cfg.grid_dim];
        npy::write_u8_array(output_dir.join("phenotype.npy"), &shape, ph)
            .expect("Failed to write phenotype.npy");
    }

    // Static Treg/MDSC suppressor-source mask (#264), same static 3D layout.
    // Only when the preset is suppressor. 1 = niche source cell, 0 = not.
    if let Some(mask) = &suppressor_mask {
        let shape = [run_cfg.grid_dim, run_cfg.grid_dim, run_cfg.grid_dim];
        npy::write_u8_array(output_dir.join("suppressor.npy"), &shape, mask)
            .expect("Failed to write suppressor.npy");
    }

    // Static diffusing-SASP-field overlay (#376/#398), f32 quasi-steady field.
    // Only when the preset is sasp-field.
    if let Some(field) = &sasp_field_overlay {
        let shape = [run_cfg.grid_dim, run_cfg.grid_dim, run_cfg.grid_dim];
        npy::write_f32_array(output_dir.join("sasp_field.npy"), &shape, field)
            .expect("Failed to write sasp_field.npy");
    }

    let meta = snapshot::TrajectoryMeta {
        schema_version: snapshot::TRAJECTORY_SCHEMA_VERSION,
        grid_dim: run_cfg.grid_dim,
        cell_size_um: CELL_SIZE_UM,
        tumor_radius_um,
        n_steps: run_cfg.n_steps,
        dose_steps,
        condition: snapshot::TrajectoryCondition {
            treatment: preset.treatment_name.to_string(),
            o2_condition: "gradient".to_string(),
            o2_lambda_um: Some(ZONE_REF_LAMBDA),
            immune_mode: if preset.immune_on { "immune_on" } else { "off" }.to_string(),
            stromal_mode: preset.stromal_on.then(|| "stromal_on".to_string()),
            ph_mode: preset.ph_on.then(|| "ph_on".to_string()),
        },
    };
    let meta_path = output_dir.join("trajectory_meta.json");
    fs::write(&meta_path, serde_json::to_string_pretty(&meta).unwrap())
        .expect("Failed to write trajectory_meta.json");

    eprintln!(
        "Wrote {} + trajectory_{{dead,damp,lp{}}}.npy{}{}{}{}",
        meta_path.display(),
        if preset.persister { ",persister" } else { "" },
        if preset.clonal { " + subclone.npy" } else { "" },
        if preset.vasculature {
            " + vessel_supply.npy"
        } else {
            ""
        },
        if preset.spheroid {
            " + phenotype.npy"
        } else {
            ""
        },
        if preset.suppressor {
            " + suppressor.npy"
        } else {
            ""
        },
    );
    if sasp_field_overlay.is_some() {
        eprintln!("  + sasp_field.npy (diffusing SASP field, #376/#398)");
    }
    eprintln!("Render with: python3 scripts/render_tme_3d_trajectory.py");
}

/// Schema version for `dose_comparison.json` (#239). Independent of
/// `summary.json`'s schema.
const DOSE_SWEEP_SCHEMA_VERSION: u32 = 1;

/// One row of the dose-sweep comparison: RSL3 kill outcome under a single
/// dosing protocol, all else equal.
#[derive(Serialize)]
struct DoseSweepEntry {
    /// Short protocol label (constant / bolus / multidose / infusion / frompk).
    schedule: String,
    /// Human-readable protocol description.
    description: String,
    total_tumor: usize,
    total_dead: usize,
    overall_kill_rate: f64,
    ferroptosis_kills: usize,
    immune_kills: usize,
}

#[derive(Serialize)]
struct DoseSweepResult {
    schema_version: u32,
    treatment: String,
    /// The fixed biological context shared by every protocol.
    context: String,
    grid_dim: usize,
    n_steps: u32,
    /// One entry per dosing protocol, all sharing the same grid + RNG seed
    /// so differences reflect the protocol alone, not stochastic noise.
    schedules: Vec<DoseSweepEntry>,
}

/// Normalized RSL3 interstitial-concentration factor series from the
/// two-compartment tumor PK ODE (#239 `DoseSchedule::FromPk` bridge).
///
/// Solves `tumor_pk::solve_tumor_pk` for an IV bolus into a breast-tumor
/// compartment, then normalizes the interstitial timecourse to peak 1.0 so
/// it reads as a drug-availability factor. This is what finally wires the
/// (previously orphaned) PK ODE into the spatial grid, without coupling the
/// grid to the solver — the grid just consumes the resulting `&[f64]`.
fn rsl3_pk_factor_series(n_steps: u32) -> Vec<f64> {
    use ferroptosis_core::tumor_pk::{breast_tumor, solve_tumor_pk, PlasmaModel};
    // IV bolus, ~35-step plasma half-life (k_el = ln2 / 35 ≈ 0.0198 /min).
    let plasma = PlasmaModel::IvBolus {
        c0: 1.0,
        k_el: 0.0198,
    };
    let tumor = breast_tumor();
    let res = solve_tumor_pk(&plasma, &tumor, n_steps as usize, 10);
    let max = res
        .c_interstitial
        .iter()
        .cloned()
        .fold(0.0_f64, f64::max)
        .max(1e-9);
    res.c_interstitial
        .iter()
        .map(|c| (c / max).clamp(0.0, 1.0))
        .collect()
}

/// `--dose-sweep`: run RSL3 across all five dosing protocols at the
/// combined-TME context (immune + stromal + pH, λ=120) and write
/// `dose_comparison.json` (#239).
///
/// **Controlled comparison**: every protocol uses the SAME tumor grid
/// (fixed `SEED`) and the SAME runtime RNG seed (all conditions share one
/// name), so the only thing that varies between rows is the dosing
/// schedule. Differences in kill rate therefore reflect the protocol, not
/// stochastic noise.
fn run_dose_sweep(output_dir: &Path) {
    eprintln!(
        "=== --dose-sweep: RSL3 across dosing protocols ({}³ × {} steps) ===",
        GRID_DIM, N_STEPS
    );

    // (label, description, schedule). All run RSL3 + immune + stromal + pH.
    //
    // The half-life / peak / level numbers below are ILLUSTRATIVE v1 values,
    // not calibrated to clinical protocols (RSL3_INACTIVATION_RATE itself was
    // tuned for sustained conc=1.0). The informative output is the
    // cross-protocol ORDERING of kill rates, not the absolute magnitudes —
    // dose_comparison.json's `context` field and the README record this.
    let protocols: Vec<(&str, &str, DoseSchedule)> = vec![
        (
            "constant",
            "Steady-state full availability (one-shot GPX4 knockdown)",
            DoseSchedule::Constant,
        ),
        (
            "bolus",
            "Single bolus at step 10, 20-step half-life",
            DoseSchedule::Bolus {
                dose_step: 10,
                peak: 1.0,
                half_life_steps: 20.0,
            },
        ),
        (
            "multidose",
            "4 doses at steps 10/55/100/145, 18-step half-life",
            DoseSchedule::MultiDose {
                dose_steps: vec![10, 55, 100, 145],
                peak: 1.0,
                half_life_steps: 18.0,
            },
        ),
        (
            "infusion",
            "Continuous infusion, level 0.5, steps 10-170",
            DoseSchedule::Infusion {
                start: 10,
                end: 170,
                level: 0.5,
            },
        ),
        (
            "frompk",
            "tumor_pk two-compartment IV-bolus interstitial curve (breast tumor)",
            DoseSchedule::FromPk {
                series: rsl3_pk_factor_series(N_STEPS),
            },
        ),
    ];

    let entries: Vec<DoseSweepEntry> = protocols
        .par_iter()
        .map(|(label, desc, sched)| {
            let cond = Condition {
                // Shared name → shared cond_seed → matched RNG across protocols.
                name: "dosesweep_RSL3".to_string(),
                treatment: Treatment::RSL3,
                treatment_name: "RSL3".to_string(),
                o2_lambda: Some(ZONE_REF_LAMBDA),
                immune_on: true,
                stromal_on: true,
                ph_on: true,
                dose_schedule: sched.clone(),
            };
            let r = run_one_condition(&cond);
            eprintln!(
                "  {label:<10} → kill={:.2}% (ferro={}, immune={})",
                r.overall_kill_rate * 100.0,
                r.ferroptosis_kills.unwrap_or(0),
                r.immune_kills.unwrap_or(0),
            );
            DoseSweepEntry {
                schedule: (*label).to_string(),
                description: (*desc).to_string(),
                total_tumor: r.total_tumor,
                total_dead: r.total_dead,
                overall_kill_rate: r.overall_kill_rate,
                ferroptosis_kills: r.ferroptosis_kills.unwrap_or(0),
                immune_kills: r.immune_kills.unwrap_or(0),
            }
        })
        .collect();

    let result = DoseSweepResult {
        schema_version: DOSE_SWEEP_SCHEMA_VERSION,
        treatment: "RSL3".to_string(),
        context: "immune + stromal + pH at λ=120 µm; shared grid + RNG seed across protocols"
            .to_string(),
        grid_dim: GRID_DIM,
        n_steps: N_STEPS,
        schedules: entries,
    };

    let path = output_dir.join("dose_comparison.json");
    fs::write(&path, serde_json::to_string_pretty(&result).unwrap())
        .expect("Failed to write dose_comparison.json");
    eprintln!("Wrote {}", path.display());
}

/// Read a positive `usize` env var, falling back to `default` if unset or
/// unparseable. Used by `--bench` to override grid size without touching the
/// byte-identical-locked `GRID_DIM` constant.
fn bench_env_usize(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&n| n > 0)
        .unwrap_or(default)
}

/// Read a positive `u32` env var (parsed directly into `u32` — values above
/// `u32::MAX` are rejected, not silently truncated). For `BENCH_N_STEPS`.
fn bench_env_u32(name: &str, default: u32) -> u32 {
    std::env::var(name)
        .ok()
        .and_then(|v| v.parse::<u32>().ok())
        .filter(|&n| n > 0)
        .unwrap_or(default)
}

/// `--bench`: run ONE representative condition (the combined-TME RSL3 cell —
/// immune + stromal + pH, the heaviest per-cell path) at a configurable grid
/// size and print wall-clock timing. The performance/scalability harness for
/// issue #192.
///
/// Grid size and step count come from env (`BENCH_GRID_DIM`, `BENCH_N_STEPS`;
/// defaults 60/180) so a sweep is scriptable, e.g.
/// `BENCH_GRID_DIM=200 BENCH_N_STEPS=180 cargo run --release -p sim-tme-3d -- --bench`.
/// Capture peak RSS externally: `/usr/bin/time -v <that command>`.
///
/// Runs a SINGLE condition with no condition-level rayon, so it measures the
/// within-condition cost — exactly the work the within-condition parallelism
/// targets. Writes no output files; the default matrix path is untouched.
fn run_bench() {
    let grid_dim = bench_env_usize("BENCH_GRID_DIM", GRID_DIM);
    let n_steps = bench_env_u32("BENCH_N_STEPS", N_STEPS);
    // saturating_pow: a huge BENCH_GRID_DIM saturates rather than overflow-
    // panicking in debug (the grid allocation would fail first anyway).
    let n_cells = grid_dim.saturating_pow(3);

    eprintln!(
        "=== --bench: single combined-TME RSL3 condition, {grid_dim}³ = {n_cells} cells × {n_steps} steps ==="
    );
    eprintln!(
        "(override via BENCH_GRID_DIM / BENCH_N_STEPS; wrap with `/usr/bin/time -v` for peak RSS)"
    );

    let condition = Condition {
        name: "bench_combined_RSL3".to_string(),
        treatment: Treatment::RSL3,
        treatment_name: "RSL3".to_string(),
        o2_lambda: Some(ZONE_REF_LAMBDA),
        immune_on: true,
        stromal_on: true,
        ph_on: true,
        dose_schedule: DoseSchedule::Constant,
    };
    let cfg = RunConfig { grid_dim, n_steps };

    let t0 = Instant::now();
    let r = run_one_condition_with_config(&condition, cfg, None);
    let elapsed = t0.elapsed();

    let secs = elapsed.as_secs_f64();
    let cell_steps = (n_cells as f64) * (n_steps as f64);
    eprintln!(
        "  done in {secs:.2}s — total_kill={:.2}% (ferro={}, immune={})",
        r.overall_kill_rate * 100.0,
        r.ferroptosis_kills.unwrap_or(0),
        r.immune_kills.unwrap_or(0),
    );
    // Machine-readable line for sweep collection (grep 'BENCH_RESULT').
    eprintln!(
        "BENCH_RESULT grid_dim={grid_dim} n_cells={n_cells} n_steps={n_steps} \
         wall_s={secs:.3} cell_steps_per_s={:.3e}",
        cell_steps / secs.max(1e-9)
    );
}

/// `--spheroid-size-sweep`: RSL3 kill vs spheroid SIZE, fixed vs size-aware zone
/// thresholds (#333 kill leg). Runs the spheroid context with the pharmacologic
/// RSL3 modality (the penetration/O2/core-limited modality that matches the
/// measured cytotoxic size-resistance direction; immune/stromal/pH off to isolate
/// it) across a range of grid sizes, twice each: fixed thresholds and the
/// size-aware thresholds (#333). Emits one machine-readable `SPHEROID_SIZE_SWEEP`
/// line per (size, mode) to stdout for `scripts/validate_spheroid_kill.py`, which
/// validates the bigger-spheroids-resist-more DIRECTION against published cytotoxic
/// size-resistance fold-ratios. Does NOT write summary.json (matrix untouched,
/// byte-identical). The O₂ λ is fixed (ZONE_REF_LAMBDA) so a larger spheroid has a
/// proportionally more hypoxic, less-penetrated core — the physical size-dependence.
fn run_spheroid_size_sweep() {
    let n_steps = bench_env_u32("BENCH_N_STEPS", N_STEPS);
    let dims: [usize; 5] = [16, 24, 32, 48, 60];
    eprintln!(
        "=== --spheroid-size-sweep: RSL3 kill vs spheroid size, fixed vs size-aware zones (#333), {n_steps} steps ==="
    );
    for &grid_dim in &dims {
        let r_um = (grid_dim as f64) * TUMOR_RADIUS_FRACTION * CELL_SIZE_UM;
        let cfg = RunConfig { grid_dim, n_steps };
        for size_aware in [false, true] {
            let condition = Condition {
                name: format!(
                    "spheroid_RSL3_d{grid_dim}_{}",
                    if size_aware { "sizeaware" } else { "fixed" }
                ),
                treatment: Treatment::RSL3,
                treatment_name: "RSL3".to_string(),
                o2_lambda: Some(ZONE_REF_LAMBDA),
                immune_on: false,
                stromal_on: false,
                ph_on: false,
                dose_schedule: DoseSchedule::Constant,
            };
            let overrides = Overrides {
                spheroid: Some(SpheroidConfig::literature()),
                spheroid_size_aware: if size_aware {
                    Some(SizeAwareZones::literature())
                } else {
                    None
                },
                ..Default::default()
            };
            let r = run_one_condition_full(&condition, cfg, None, overrides);
            println!(
                "SPHEROID_SIZE_SWEEP grid_dim={grid_dim} radius_um={r_um:.1} size_aware={size_aware} \
                 kill_rate={:.6} total_tumor={} total_dead={}",
                r.overall_kill_rate, r.total_tumor, r.total_dead
            );
        }
    }
}

fn main() {
    // Guard against silent drift between this binary's metadata const and
    // the library's runtime value: if a future PR tunes
    // `SpatialImmuneConfig::for_3d().damp_diffusion_fraction` without also
    // touching `DAMP_DIFFUSION_FRACTION_3D`, `summary.json` would report a
    // stale value while the simulation actually used the new one. This
    // debug_assert catches that in tests/dev runs without affecting release
    // bit-identical output.
    debug_assert_eq!(
        DAMP_DIFFUSION_FRACTION_3D,
        SpatialImmuneConfig::for_3d().damp_diffusion_fraction,
        "DAMP_DIFFUSION_FRACTION_3D metadata const out of sync with SpatialImmuneConfig::for_3d()",
    );

    // Use the library's TUMOR_RADIUS_FRACTION constant instead of the bare
    // 0.45 literal — keeps `tumor_radius_um` in lockstep with whatever
    // value `TumorGrid3D::generate` actually uses (reviewer-flagged drift
    // hazard if the library ever tunes the fraction).
    let tumor_radius_um = (GRID_DIM as f64) * TUMOR_RADIUS_FRACTION * CELL_SIZE_UM;
    eprintln!("=== sim-tme-3d: 3D Spheroid TME Simulation ===");
    eprintln!(
        "Grid: {0}³ ({1:.1} mm × {1:.1} mm × {1:.1} mm)",
        GRID_DIM,
        GRID_DIM as f64 * CELL_SIZE_UM / 1000.0
    );
    eprintln!("Tumor radius: {:.0} µm", tumor_radius_um);
    eprintln!(
        "O₂ λ sweep: {:?} µm (λ=150 skipped — 3λ > tumor radius)",
        O2_LAMBDAS
    );
    eprintln!(
        "DAMP diffusion fraction: {} (3D-safe; 2D's 0.08 would trigger immune_spatial stability assert!)",
        DAMP_DIFFUSION_FRACTION_3D
    );
    eprintln!();
    eprintln!(
        "⚠️  Scale mismatch with sim-tme (2D): tumor radius {:.0} µm vs sim-tme's ~4500 µm.",
        tumor_radius_um
    );
    eprintln!("    Compare via RATIOS, not absolute counts.");
    eprintln!();

    let output_dir = Path::new("output/tme-3d");
    fs::create_dir_all(output_dir).expect("Failed to create output/tme-3d");

    // `--snapshot[=NAME]` runs ONE visualization-focused condition with
    // per-step state capture for the Python animation. Default path (no
    // flag) is unchanged: runs the full 24-condition matrix →
    // summary.json. Bit-identical when the flag is absent (snapshot path
    // doesn't touch the matrix). Names listed in `SNAPSHOTS`; bare
    // `--snapshot` defaults to `combined`.
    if let Some(name) = parse_snapshot_arg(std::env::args()) {
        run_snapshot(output_dir, tumor_radius_um, &name);
        return;
    }

    // `--dose-sweep` runs RSL3 across all five dosing protocols (Constant /
    // Bolus / MultiDose / Infusion / FromPk) at the combined-TME context and
    // writes `dose_comparison.json` — the quantitative "does dosing protocol
    // change efficacy?" answer (#239). Separate file; does NOT touch
    // summary.json (which stays byte-identical to the 24-condition matrix).
    if std::env::args().any(|a| a == "--dose-sweep") {
        run_dose_sweep(output_dir);
        return;
    }

    // `--bench` runs ONE representative condition at a configurable grid size
    // (env BENCH_GRID_DIM / BENCH_N_STEPS, defaults 60/180) and prints
    // wall-clock timing — the performance/scalability harness for #192. It
    // measures the WITHIN-condition cost (no condition-level rayon), which is
    // exactly what the within-condition parallelism targets. Does NOT write
    // summary.json; the default matrix path below is untouched (byte-identical).
    if std::env::args().any(|a| a == "--bench") {
        run_bench();
        return;
    }

    // `--spheroid-size-sweep` runs RSL3 across a range of spheroid sizes, fixed vs
    // size-aware zone thresholds (#333 kill leg), emitting SPHEROID_SIZE_SWEEP lines.
    // Does NOT write summary.json; the default matrix path below is byte-identical.
    if std::env::args().any(|a| a == "--spheroid-size-sweep") {
        run_spheroid_size_sweep();
        return;
    }

    let conditions = generate_conditions();
    eprintln!(
        "Running {} conditions in parallel via rayon...",
        conditions.len()
    );

    // Condition-level parallelism. Each run_one_condition allocates its own
    // grid + DAMP/scratch buffers; no shared mutable state.
    let results: Vec<ConditionResult> = conditions
        .par_iter()
        .map(|cond| {
            eprintln!("  starting: {}", cond.name);
            let r = run_one_condition(cond);
            eprintln!(
                "  done:     {} — total_kill={:.1}% (norm={:.1}%, trans={:.1}%, hyp={:.1}%)",
                cond.name,
                r.overall_kill_rate * 100.0,
                r.normoxic_kill_rate * 100.0,
                r.transition_kill_rate * 100.0,
                r.hypoxic_kill_rate * 100.0
            );
            r
        })
        .collect();

    let summary = SimulationSummary {
        schema_version: TME_3D_SCHEMA_VERSION,
        grid_dim: GRID_DIM,
        cell_size_um: CELL_SIZE_UM,
        tumor_radius_um,
        n_steps: N_STEPS,
        o2_lambdas: O2_LAMBDAS.to_vec(),
        damp_diffusion_fraction_3d: DAMP_DIFFUSION_FRACTION_3D,
        conditions: results,
        note: format!(
            "sim-tme-3d uses {}³ grid (tumor radius {:.0} µm) vs sim-tme's 500×500 \
             (tumor radius ~4500 µm). Compare conditions via RATIOS, not absolute counts. \
             See generate_3d_comparison_table.py.",
            GRID_DIM, tumor_radius_um
        ),
    };

    let json_path = output_dir.join("summary.json");
    let json = serde_json::to_string_pretty(&summary).expect("serialize summary");
    fs::write(&json_path, json).expect("write summary.json");
    eprintln!("\nWrote {}", json_path.display());
    eprintln!("Done.");
}

// ============================================================
// Tests (smoke-level — primitives are already tested in
// ferroptosis-core)
// ============================================================

#[cfg(test)]
mod tests {
    use super::*;
    // Test-only: used by the clonal lipid-axis test.
    use ferroptosis_core::clonal::SubclonePerturbation;

    // Golden integer kill counts for the Constant-path regression guard
    // (`constant_path_golden_kill_counts`), captured from SDT + immune +
    // stromal + pH at grid_dim=20, n_steps=80. Deterministic on a fixed
    // platform/toolchain (seeded RNG + IEEE f64); may differ by an ULP-edge
    // count on other libm/CPU builds — see the test's doc. Update ONLY after
    // confirming a default-path change is intentional.
    const GOLDEN_TOTAL_TUMOR: usize = 3071;
    const GOLDEN_TOTAL_DEAD: usize = 2992;
    const GOLDEN_FERRO_KILLS: usize = 2990;
    const GOLDEN_IMMUNE_KILLS: usize = 2;
    // Dosed-path golden (MultiDose SDT + immune + stromal + pH, 20³×80).
    const GOLDEN_DOSED_TOTAL_DEAD: usize = 2660;
    const GOLDEN_DOSED_FERRO_KILLS: usize = 2658;
    const GOLDEN_DOSED_IMMUNE_KILLS: usize = 2;

    /// `--snapshot` (no `=NAME`) defaults to `combined`.
    #[test]
    fn parse_snapshot_arg_bare_defaults_to_combined() {
        let args = vec!["sim-tme-3d".to_string(), "--snapshot".to_string()];
        assert_eq!(parse_snapshot_arg(args), Some("combined".to_string()));
    }

    /// `--snapshot=bare` extracts the name after `=`.
    #[test]
    fn parse_snapshot_arg_equals_form_extracts_name() {
        let args = vec!["sim-tme-3d".to_string(), "--snapshot=bare".to_string()];
        assert_eq!(parse_snapshot_arg(args), Some("bare".to_string()));
    }

    /// Absent flag returns None (the default 24-condition matrix path).
    #[test]
    fn parse_snapshot_arg_absent_returns_none() {
        let args = vec!["sim-tme-3d".to_string()];
        assert_eq!(parse_snapshot_arg(args), None);
    }

    /// All `SNAPSHOTS` entries have unique names (used as the `=NAME`
    /// match key downstream). A typo or copy-paste duplicate would
    /// make `resolve_snapshot` ambiguous (first-match-wins).
    #[test]
    fn snapshot_preset_names_are_unique() {
        let mut names: Vec<&str> = SNAPSHOTS.iter().map(|s| s.name).collect();
        names.sort();
        let original = names.len();
        names.dedup();
        assert_eq!(names.len(), original, "duplicate snapshot preset name");
    }

    /// #302: lock the `--snapshot=contact` preset → `Overrides` wiring without
    /// running the full (~333 MB) render. Resolves the preset and asserts it
    /// turns contact ON and leaves the geometry layers OFF — contact is
    /// mutually-exclusive with slab, and must run on the centred sphere (no
    /// spheroid re-grid) where the fixed-26 contact denominator is correct.
    #[test]
    fn contact_snapshot_preset_is_wired() {
        let p = resolve_snapshot("contact");
        assert_eq!(p.name, "contact");
        assert!(
            p.contact,
            "the contact preset must enable the contact layer"
        );
        assert!(
            !p.slab && !p.spheroid && !p.clonal && !p.vasculature,
            "contact runs on the plain centred sphere (no conflicting geometry layer)"
        );
        assert!(
            matches!(p.treatment, Treatment::RSL3),
            "contact is shown under RSL3"
        );
    }

    /// Smoke test: condition matrix is non-empty and well-formed.
    #[test]
    fn condition_matrix_is_non_empty() {
        let conditions = generate_conditions();
        assert!(!conditions.is_empty(), "expected at least some conditions");
        // Names should be unique (used as seed-hash inputs).
        let mut names: Vec<&str> = conditions.iter().map(|c| c.name.as_str()).collect();
        names.sort();
        let original_count = names.len();
        names.dedup();
        assert_eq!(
            names.len(),
            original_count,
            "condition names must be unique"
        );
    }

    /// Smoke test: one baseline condition (Control, no toggles) runs
    /// end-to-end and produces sensible output.
    ///
    /// Uses `RunConfig::for_test()` (10³ × 20 steps) instead of production
    /// (60³ × 180 steps) — reviewer flagged that `cargo test --workspace`
    /// runs in DEBUG mode, where the full 60³ × 180 sim was costing minutes
    /// per test. The smaller config exercises every code path at ~1000×
    /// less cost.
    #[test]
    fn single_condition_runs_end_to_end() {
        let cond = Condition {
            name: "test_baseline_Control".to_string(),
            treatment: Treatment::Control,
            treatment_name: "Control".to_string(),
            o2_lambda: None,
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let r = run_one_condition_with_config(&cond, RunConfig::for_test(), None);
        assert_eq!(r.treatment, "Control");
        assert_eq!(r.o2_condition, "uniform");
        assert_eq!(r.immune_mode, "off");
        assert!(r.total_tumor > 0, "expected some tumor cells in 10³ grid");
        // Control has zero exo-ROS, no immune pressure, no pH stress. Production
        // 60³ run gives ~0.15% baseline kill rate. Threshold tightened from
        // <5% (too loose — would let a 10× baseline regression sail through)
        // to <2% (still loose for safety margin given the small test grid;
        // production-rate × 5 reasonable margin).
        assert!(
            r.overall_kill_rate < 0.02,
            "Control should have <2% kill rate, got {:.1}%",
            r.overall_kill_rate * 100.0
        );
    }

    /// Determinism: same condition seed → identical output.
    #[test]
    fn same_seed_same_output() {
        let cond = Condition {
            name: "deterministic_test".to_string(),
            treatment: Treatment::Control,
            treatment_name: "Control".to_string(),
            o2_lambda: None,
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let r1 = run_one_condition_with_config(&cond, RunConfig::for_test(), None);
        let r2 = run_one_condition_with_config(&cond, RunConfig::for_test(), None);
        assert_eq!(r1.total_dead, r2.total_dead);
        assert_eq!(r1.total_tumor, r2.total_tumor);
        assert_eq!(r1.overall_kill_rate, r2.overall_kill_rate);
    }

    /// **Reviewer-flagged invariant guard**: immune-on conditions with a
    /// treatment that produces ferroptotic kills (RSL3 or SDT) must
    /// produce at least one immune kill. If the immune block silently
    /// no-ops (e.g., DAMP threshold never crossed, or a bug in the
    /// activation chain), `immune_kills` would still be 0 and the JSON
    /// would still serialize — bug only surfaces in the manuscript.
    ///
    /// Uses RSL3 + O₂ gradient at λ=120 (the canonical "immune on"
    /// condition) and asserts `immune_kills > 0` over 20 steps. SDT
    /// also kills heavily so picking RSL3 makes the assertion stricter
    /// (RSL3 produces fewer ferroptotic deaths → fewer DAMP sources →
    /// harder for any immune kills to fire spuriously).
    #[test]
    fn immune_on_actually_fires() {
        let cond = Condition {
            name: "test_immune_RSL3".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        // Larger config than the smoke tests: need IMMUNE_START_STEP=60 +
        // buffer for DAMP to accumulate above the 0.01 kill threshold.
        // 25³ × 130 steps: ~1500 tumor cells, 70 steps post-immune-start
        // for DAMP to spread and reach kill threshold. Empirically reliable
        // (production 60³×180 produces 29 RSL3 immune kills; volume scales
        // 60³/25³ ≈ 14×, so expect ~2 kills at this config — comfortable
        // margin over the >0 invariant).
        let cfg = RunConfig {
            grid_dim: 25,
            n_steps: 130,
        };
        let r = run_one_condition_with_config(&cond, cfg, None);
        assert_eq!(r.immune_mode, "immune_on");
        let im_kills = r
            .immune_kills
            .expect("immune_on must populate immune_kills");
        // Strict invariant: at least one immune kill in this RSL3 + immune-on
        // run. If this fails the immune block has gone no-op.
        assert!(
            im_kills > 0,
            "immune-on RSL3 should produce ≥1 immune kill in {} steps on {grid}³; got {im_kills}. \
             Likely cause: DAMP threshold never crossed, or activation chain broken.",
            cfg.n_steps,
            grid = cfg.grid_dim
        );
    }

    /// **Reviewer-flagged invariant guard**: for immune-on conditions,
    /// `total_dead` must equal `ferroptosis_kills + immune_kills`.
    /// Catches double-counting drift (e.g., a cell counted in both
    /// branches) and missed-kill drift (e.g., a kill that mutates
    /// `state.dead` but doesn't increment either counter).
    #[test]
    fn total_dead_equals_ferro_plus_immune() {
        let cond = Condition {
            name: "test_invariant_RSL3".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let cfg = RunConfig {
            grid_dim: 15,
            n_steps: 80,
        };
        let r = run_one_condition_with_config(&cond, cfg, None);
        let ferro = r.ferroptosis_kills.expect("immune_on populates this");
        let im = r.immune_kills.expect("immune_on populates this");
        assert_eq!(
            r.total_dead,
            ferro + im,
            "total_dead {} != ferroptosis_kills {} + immune_kills {} — kill accounting drifted",
            r.total_dead,
            ferro,
            im
        );
    }

    // ============================================================
    // Time-varying dose schedule (#239)
    // ============================================================

    /// A zero-availability schedule (a bolus scheduled past the end of the
    /// run, so `factor_at ≡ 0`) must suppress SDT kills relative to the
    /// Constant full-availability default — proving the SDT exo-ROS rescale
    /// path is actually wired in.
    #[test]
    fn dosed_zero_availability_schedule_suppresses_sdt_kills() {
        let cfg = RunConfig {
            grid_dim: 20,
            n_steps: 60,
        };
        let constant = run_one_condition_with_config(
            &Condition {
                name: "test_sdt_constant".to_string(),
                treatment: Treatment::SDT,
                treatment_name: "SDT".to_string(),
                o2_lambda: Some(ZONE_REF_LAMBDA),
                immune_on: false,
                stromal_on: false,
                ph_on: false,
                dose_schedule: DoseSchedule::Constant,
            },
            cfg,
            None,
        );
        let zero = run_one_condition_with_config(
            &Condition {
                name: "test_sdt_zero_drug".to_string(),
                treatment: Treatment::SDT,
                treatment_name: "SDT".to_string(),
                o2_lambda: Some(ZONE_REF_LAMBDA),
                immune_on: false,
                stromal_on: false,
                ph_on: false,
                // Dose never arrives within the run → factor_at ≡ 0 → no SDT ROS.
                dose_schedule: DoseSchedule::Bolus {
                    dose_step: 9999,
                    peak: 1.0,
                    half_life_steps: 10.0,
                },
            },
            cfg,
            None,
        );
        assert!(
            constant.total_dead > 0,
            "Constant SDT should kill cells at 20³×60; got {}",
            constant.total_dead
        );
        assert!(
            zero.total_dead < constant.total_dead,
            "zero-availability schedule must suppress SDT kills: zero={}, constant={}",
            zero.total_dead,
            constant.total_dead
        );
    }

    /// A gentle single-bolus RSL3 schedule (continuous covalent inactivation)
    /// must kill fewer cells than the Constant one-shot 92% GPX4 knockdown —
    /// proving the RSL3 dosed path is live and the no-init-knockdown +
    /// per-step inactivation mechanism differs from the steady-state default.
    #[test]
    fn dosed_rsl3_bolus_kills_less_than_constant_oneshot() {
        let cfg = RunConfig {
            grid_dim: 20,
            n_steps: 60,
        };
        let constant = run_one_condition_with_config(
            &Condition {
                name: "test_rsl3_constant".to_string(),
                treatment: Treatment::RSL3,
                treatment_name: "RSL3".to_string(),
                o2_lambda: Some(ZONE_REF_LAMBDA),
                immune_on: false,
                stromal_on: false,
                ph_on: false,
                dose_schedule: DoseSchedule::Constant,
            },
            cfg,
            None,
        );
        let bolus = run_one_condition_with_config(
            &Condition {
                name: "test_rsl3_bolus".to_string(),
                treatment: Treatment::RSL3,
                treatment_name: "RSL3".to_string(),
                o2_lambda: Some(ZONE_REF_LAMBDA),
                immune_on: false,
                stromal_on: false,
                ph_on: false,
                dose_schedule: DoseSchedule::Bolus {
                    dose_step: 2,
                    peak: 1.0,
                    half_life_steps: 8.0,
                },
            },
            cfg,
            None,
        );
        assert!(
            constant.total_dead > 0,
            "Constant RSL3 (one-shot 92% knockdown) should kill cells; got {}",
            constant.total_dead
        );
        assert!(
            bolus.total_dead < constant.total_dead,
            "a single decaying RSL3 bolus (gentle continuous inactivation) must kill \
             fewer than the Constant one-shot: bolus={}, constant={}",
            bolus.total_dead,
            constant.total_dead
        );
    }

    /// **Byte-identity regression guard for the Constant default path.**
    /// Pins the exact integer kill counts of a representative Constant
    /// condition (SDT + immune + stromal + pH) at a fixed small config,
    /// catching ANY drift in the default (non-dosed) path — the load-bearing
    /// "summary.json byte-identical" property — at unit-test speed, without
    /// the full 60³×180 production run.
    ///
    /// **Scope: deterministic on a fixed platform + toolchain.** The kill
    /// counts are integers, but the death decision is a strict
    /// `lp > threshold` on values flowing through `powf` and the Ziggurat
    /// normal sampler, which are not bit-identical across libm/CPU
    /// implementations — so a cell within an ULP of the threshold could flip
    /// a count on a different platform. (This is a pre-existing
    /// whole-simulation property, not introduced by #239.) This test is the
    /// same-platform CI guard against accidental drift; the production
    /// SHA-256 check (done manually in #239's PR) is the cross-build
    /// authority. If these numbers change, the byte-identical claim is broken
    /// — investigate before updating.
    #[test]
    fn constant_path_golden_kill_counts() {
        let cond = Condition {
            name: "golden_constant_SDT".to_string(),
            treatment: Treatment::SDT,
            treatment_name: "SDT".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: true,
            ph_on: true,
            dose_schedule: DoseSchedule::Constant,
        };
        let cfg = RunConfig {
            grid_dim: 20,
            n_steps: 80,
        };
        let r = run_one_condition_with_config(&cond, cfg, None);
        // Golden values captured from the Constant path. A change here means
        // the default-path numerics drifted — investigate before updating.
        assert_eq!(
            r.total_tumor, GOLDEN_TOTAL_TUMOR,
            "tumor-cell count drifted"
        );
        assert_eq!(
            r.total_dead, GOLDEN_TOTAL_DEAD,
            "Constant-path total_dead drifted"
        );
        assert_eq!(
            r.ferroptosis_kills,
            Some(GOLDEN_FERRO_KILLS),
            "Constant-path ferroptosis_kills drifted"
        );
        assert_eq!(
            r.immune_kills,
            Some(GOLDEN_IMMUNE_KILLS),
            "Constant-path immune_kills drifted"
        );
    }

    /// **Golden guard for the DOSED parallel path** (review #255 substantive
    /// point 2). `constant_path_golden_kill_counts` pins the Constant path;
    /// this pins a MultiDose SDT + immune + stromal + pH condition, which
    /// exercises the machinery the Constant golden does NOT: the per-step
    /// SDT exo-envelope divide-out, the grace-end DAMP write, and the
    /// immune-kill loop — all through the rayon-parallelized loops. Same
    /// fixed-platform/toolchain caveat as the Constant golden (integer
    /// counts; `powf`/Ziggurat aren't bit-identical across libm).
    #[test]
    fn dosed_path_golden_kill_counts() {
        let cond = Condition {
            name: "golden_multidose_SDT".to_string(),
            treatment: Treatment::SDT,
            treatment_name: "SDT".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: true,
            ph_on: true,
            dose_schedule: DoseSchedule::MultiDose {
                dose_steps: vec![5, 30, 55],
                peak: 1.0,
                half_life_steps: 8.0,
            },
        };
        let cfg = RunConfig {
            grid_dim: 20,
            n_steps: 80,
        };
        let r = run_one_condition_with_config(&cond, cfg, None);
        assert_eq!(
            r.total_tumor, GOLDEN_TOTAL_TUMOR,
            "tumor-cell count drifted"
        );
        assert_eq!(
            r.total_dead, GOLDEN_DOSED_TOTAL_DEAD,
            "dosed-path total_dead drifted"
        );
        assert_eq!(
            r.ferroptosis_kills,
            Some(GOLDEN_DOSED_FERRO_KILLS),
            "dosed-path ferroptosis_kills drifted"
        );
        assert_eq!(
            r.immune_kills,
            Some(GOLDEN_DOSED_IMMUNE_KILLS),
            "dosed-path immune_kills drifted"
        );
    }

    /// The full dosed stack (RSL3 + immune + stromal + pH + MultiDose) must
    /// run end-to-end without panicking and be deterministic across runs.
    /// Exercises `rsl3_drug_avail` indexing (pH on), the no-knockdown init,
    /// per-step inactivation, and the immune/stromal interactions together.
    #[test]
    fn dosed_rsl3_full_stack_runs_and_is_deterministic() {
        let cond = Condition {
            name: "test_rsl3_multidose_full".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: true,
            ph_on: true,
            dose_schedule: DoseSchedule::MultiDose {
                dose_steps: vec![2, 10, 18],
                peak: 1.0,
                half_life_steps: 6.0,
            },
        };
        let cfg = RunConfig {
            grid_dim: 20,
            n_steps: 30,
        };
        let r1 = run_one_condition_with_config(&cond, cfg, None);
        let r2 = run_one_condition_with_config(&cond, cfg, None);
        assert_eq!(
            r1.total_dead, r2.total_dead,
            "dosed full-stack run must be deterministic"
        );
        assert!(r1.total_tumor > 0, "expected tumor cells");
        assert!(
            (0.0..=1.0).contains(&r1.overall_kill_rate),
            "kill rate must be a valid fraction, got {}",
            r1.overall_kill_rate
        );
    }

    /// End-to-end coverage for the `Infusion` variant (review #10): runs
    /// RSL3 under a continuous infusion through the full per-step dosed
    /// path and asserts it's deterministic and produces a valid fraction.
    #[test]
    fn dosed_rsl3_infusion_runs_and_is_deterministic() {
        let cond = Condition {
            name: "test_rsl3_infusion".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: false,
            ph_on: true,
            dose_schedule: DoseSchedule::Infusion {
                start: 2,
                end: 28,
                level: 1.0,
            },
        };
        let cfg = RunConfig {
            grid_dim: 20,
            n_steps: 30,
        };
        let r1 = run_one_condition_with_config(&cond, cfg, None);
        let r2 = run_one_condition_with_config(&cond, cfg, None);
        assert_eq!(
            r1.total_dead, r2.total_dead,
            "Infusion run must be deterministic"
        );
        assert!(r1.total_tumor > 0);
        assert!((0.0..=1.0).contains(&r1.overall_kill_rate));
    }

    /// End-to-end coverage for the `FromPk` variant (review #10): runs RSL3
    /// driven by the tumor_pk-derived availability series through the full
    /// per-step dosed path. Exercises the ODE-bridge → grid wiring beyond
    /// the `factor_at`-level unit test.
    #[test]
    fn dosed_rsl3_frompk_runs_and_is_deterministic() {
        let cfg = RunConfig {
            grid_dim: 20,
            n_steps: 30,
        };
        let cond = Condition {
            name: "test_rsl3_frompk".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: false,
            ph_on: true,
            dose_schedule: DoseSchedule::FromPk {
                series: rsl3_pk_factor_series(cfg.n_steps),
            },
        };
        let r1 = run_one_condition_with_config(&cond, cfg, None);
        let r2 = run_one_condition_with_config(&cond, cfg, None);
        assert_eq!(
            r1.total_dead, r2.total_dead,
            "FromPk run must be deterministic"
        );
        assert!(r1.total_tumor > 0);
        assert!((0.0..=1.0).contains(&r1.overall_kill_rate));
    }

    /// The FromPk bridge series must be a valid normalized availability
    /// factor: correct length, all values in [0, 1], and peak ≈ 1.0
    /// (normalization target). Guards the `tumor_pk` → spatial-grid wiring.
    #[test]
    fn rsl3_pk_factor_series_is_normalized() {
        let series = rsl3_pk_factor_series(180);
        assert_eq!(series.len(), 180, "series length must match n_steps");
        assert!(
            series.iter().all(|&f| (0.0..=1.0).contains(&f)),
            "all factors must be in [0, 1]"
        );
        let peak = series.iter().cloned().fold(0.0_f64, f64::max);
        assert!(
            (peak - 1.0).abs() < 1e-9,
            "series must be normalized to peak 1.0, got {peak}"
        );
    }

    /// Empirical NaN check on `norm()` per reviewer ask: confirms that a
    /// NaN `std` propagates to a panic (via `Normal::new(...).expect(...)`),
    /// not a silent NaN value. Documents the failure mode rather than
    /// making it "safe".
    #[test]
    #[should_panic(expected = "valid normal")]
    fn norm_panics_on_nan_std() {
        let mut rng = StdRng::seed_from_u64(0);
        norm(&mut rng, 1.0, f64::NAN);
    }

    /// **Reviewer-flagged invariant**: par_iter must produce identical
    /// results across runs. Per-condition closures are independent (each
    /// allocates its own grid + DAMP/scratch buffers), so determinism
    /// should hold — but an unverified invariant is one rayon refactor
    /// away from silent drift. This test exercises the actual
    /// `par_iter().map(run_one_condition)` path (smaller config for CI
    /// debug speed) and compares the resulting JSON.
    #[test]
    fn rayon_run_is_deterministic() {
        // Subset of 6 conditions (small enough for debug mode, large
        // enough to exercise inter-thread scheduling).
        let small_conditions: Vec<_> = generate_conditions().into_iter().take(6).collect();
        let cfg = RunConfig::for_test();

        let run = |conds: &[Condition]| -> Vec<ConditionResult> {
            conds
                .par_iter()
                .map(|c| run_one_condition_with_config(c, cfg, None))
                .collect()
        };

        let r1 = run(&small_conditions);
        let r2 = run(&small_conditions);

        assert_eq!(r1.len(), r2.len());
        for (a, b) in r1.iter().zip(r2.iter()) {
            assert_eq!(a.treatment, b.treatment);
            assert_eq!(a.total_dead, b.total_dead);
            assert_eq!(a.overall_kill_rate, b.overall_kill_rate);
            assert_eq!(a.peak_damp, b.peak_damp);
            assert_eq!(a.total_damp, b.total_damp);
        }
    }

    /// **Reviewer-flagged invariant**: the library functions
    /// `radial_o2_field`, `radial_ph_field`, and `stromal_adjacency_mask_3d`
    /// all return `Vec<T>` of length `grid.cells.len()`, and this binary
    /// indexes those vecs by `grid.flat_index(r, c, l)`. The implicit
    /// contract is that the library iterates in the SAME row-major order
    /// as the grid's flat-index formula. If a library refactor ever
    /// changes the iteration order (e.g., column-major instead of
    /// row-major), every cell would silently get the wrong O₂/pH/stromal
    /// value and nothing in the smoke tests would catch it.
    ///
    /// This test pins the contract with two spot-checks:
    /// 1. Vec lengths match `grid.cells.len()`.
    /// 2. A known cell has expected qualitative properties:
    ///    - The center cell has higher O₂ depth (less O₂) than a corner cell
    ///    - The center cell has lower pH than a corner cell
    ///    - The center cell is NOT in the stromal-adjacency mask (interior),
    ///      while a near-surface tumor cell IS.
    ///
    /// A column-major refactor would scramble these properties and fail.
    #[test]
    fn library_field_order_matches_flat_index() {
        use ferroptosis_core::oxygen::radial_o2_field;
        use ferroptosis_core::ph::radial_ph_field;

        let grid = TumorGrid3D::generate(20, 20, 20, CELL_SIZE_UM, SEED);
        let n = grid.cells.len();

        let o2 = radial_o2_field(&grid, 100.0);
        let ph = radial_ph_field(&grid, 7.4, 6.5, 120.0);
        let mask = stromal_adjacency_mask_3d(&grid);
        assert_eq!(o2.len(), n, "radial_o2_field length contract");
        assert_eq!(ph.len(), n, "radial_ph_field length contract");
        assert_eq!(mask.len(), n, "stromal_adjacency_mask_3d length contract");

        // Index the same cells via the binary's `flat_index` access path:
        let center_idx = grid.flat_index(10, 10, 10);
        let surface_idx = grid.flat_index(10, 10, 19); // depth=0 surface (per grid::tests_3d)

        // Center cell is interior tumor → high pH gradient depth → lower pH
        // than surface; lower O₂ than surface; mask=false (no stromal neighbors).
        assert!(
            grid.cells[center_idx].is_tumor,
            "test precondition: center is tumor"
        );
        assert!(
            grid.cells[surface_idx].is_tumor,
            "test precondition: surface is tumor"
        );
        assert!(
            ph[center_idx] < ph[surface_idx],
            "center pH {} should be lower than surface pH {} — \
             if this fails, library field order is desync'd from flat_index",
            ph[center_idx],
            ph[surface_idx]
        );
        assert!(
            o2[center_idx] < o2[surface_idx],
            "center O₂ {} should be lower than surface O₂ {} (deeper = more hypoxic)",
            o2[center_idx],
            o2[surface_idx]
        );
        // Surface cell at (10,10,19) is at the spheroid edge → has stromal
        // neighbors at (10,10,20)? No — (10,10,20) is OUTSIDE the grid bounds
        // for a 20³ grid. Surface cell's neighbors include some that are
        // outside the sphere → mask = true.
        assert!(
            mask[surface_idx],
            "surface tumor cell should be in adjacency mask (has stromal neighbors)"
        );
        assert!(
            !mask[center_idx],
            "interior center cell should NOT be in adjacency mask"
        );
    }

    // ===== Persister-cell model (#241) =====

    /// Persister model OFF (and the identity config) must be inert: the same
    /// kills as the un-modeled path, and `persister_mean` omitted when off.
    /// This is the per-PR side of the byte-identity guarantee (the production
    /// SHA guard is #253).
    #[test]
    fn persister_off_is_inert_and_unreported() {
        let cfg = RunConfig {
            grid_dim: 20,
            n_steps: 60,
        };
        let cond = Condition {
            name: "persister_off".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let wrapped = run_one_condition_with_config(&cond, cfg, None);
        let explicit_none = run_one_condition_full(&cond, cfg, None, Overrides::default());
        assert_eq!(wrapped.total_dead, explicit_none.total_dead);
        assert_eq!(
            wrapped.persister_mean, None,
            "off path must omit persister_mean"
        );
        assert_eq!(
            explicit_none.persister_mean, None,
            "explicit-None path must also omit persister_mean"
        );
        // The #342 locked-pool metric is likewise omitted when persister is off
        // (pins the skip_serializing_if guard behind the matrix byte-identity).
        assert_eq!(
            wrapped.persister_locked_mean, None,
            "off path must omit persister_locked_mean"
        );
        assert_eq!(
            explicit_none.persister_locked_mean, None,
            "explicit-None path must also omit persister_locked_mean"
        );
        // An identity PersisterConfig runs the (gated) block but every helper
        // is a no-op, so kills are unchanged.
        let identity = run_one_condition_full(
            &cond,
            cfg,
            None,
            Overrides {
                persister: Some(PersisterConfig::default()),
                ..Default::default()
            },
        );
        assert_eq!(
            identity.total_dead, wrapped.total_dead,
            "identity PersisterConfig must not change kills"
        );
    }

    /// Headline #241 result: under a repeated-dose RSL3 schedule, enabling the
    /// persister model leaves materially MORE cells alive than the
    /// no-persister baseline. Acquired tolerance resists each dose's covalent
    /// GPX4 knockdown (gpx4_inactivation_multiplier), so kill efficiency
    /// declines across cycles — the Hangauer et al. 2017 persister effect.
    /// RSL3 (not SDT) is used because the GPX4-resistance axis is RSL3's
    /// covalent mechanism; the other TME protections are off so the delta is
    /// purely the persister effect. Observed at the time of writing:
    /// off=79, on=27 (≈66% reduction).
    #[test]
    fn persister_reduces_multidose_kills() {
        // 120 steps: the dosed RSL3 per-step inactivation is gradual (no
        // one-shot init knockdown), so cumulative GPX4 loss only crosses the
        // lethal threshold after ~100 steps at this grid. Shorter runs kill
        // zero, leaving nothing to attenuate.
        let cfg = RunConfig {
            grid_dim: 20,
            n_steps: 120,
        };
        let cond = Condition {
            name: "persister_multidose".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            // Uniform oxygenation: the hypoxia gradient collapses RSL3 kill
            // (manuscript finding 3.7%→0.1%), which would leave no baseline
            // headroom to demonstrate the persister reduction.
            o2_lambda: None,
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            // Repeated RSL3 dosing across the run, dense enough (every 4 steps,
            // overlapping half-lives) that effective availability stays near
            // saturation — sparse pulses barely dent GPX4. As persistence
            // accrues it attenuates each dose's covalent knockdown
            // (gpx4_inactivation_multiplier), so fewer cells die than baseline.
            dose_schedule: DoseSchedule::MultiDose {
                dose_steps: (0..120).step_by(4).collect(),
                peak: 1.0,
                half_life_steps: 10.0,
            },
        };
        let off = run_one_condition_with_config(&cond, cfg, None);
        let on = run_one_condition_full(
            &cond,
            cfg,
            None,
            Overrides {
                persister: Some(PersisterConfig::enabled()),
                ..Default::default()
            },
        );
        // Strong precondition: if a grid/timing change drops the baseline kill
        // count near zero, there is no headroom to demonstrate the effect and
        // the test would pass/fail for the wrong reason. Fail loudly instead.
        assert!(
            off.total_dead >= 20,
            "baseline RSL3 must kill a meaningful number of cells for this test \
             to be informative; got {}. Adjust the schedule/grid/steps.",
            off.total_dead
        );
        // Material reduction, not a 1-cell rounding artifact. Observed ≈55%
        // (off=80, on=36) under the #262 competing-rate model — smaller than the
        // pre-#262 acquire-only ≈66% because reversion now also operates.
        let reduction = (off.total_dead - on.total_dead) as f64 / off.total_dead as f64;
        assert!(
            on.total_dead < off.total_dead && reduction > 0.2,
            "persister tolerance must materially reduce kills (>20%): \
             on={}, off={}, reduction={:.1}%",
            on.total_dead,
            off.total_dead,
            reduction * 100.0
        );
        let pm = on
            .persister_mean
            .expect("persister_mean must be reported when the model is enabled");
        let max_frac = PersisterConfig::enabled().max_fraction;
        assert!(
            pm > 0.0 && pm <= max_frac,
            "persister_mean must be in (0, max_fraction={max_frac}]; got {pm}"
        );
    }

    /// #470: the persister OXPHOS-ROS suppression axis. Drug-tolerant persisters
    /// downregulate OXPHOS, lowering the basal/mitochondrial ROS a GPX4 inhibitor
    /// needs, so RSL3 kills them LESS than the persister-without-OXPHOS baseline;
    /// an HDAC inhibitor re-raises the ROS and restores the kill (PMID 40909720).
    /// A/B with the OXPHOS-suppression config as the only difference. Both the
    /// HDAC-fully-rescued case and the suppression=0 case reproduce the persister
    /// baseline EXACTLY (the byte-identity invariant: `enabled()` keeps the axis
    /// off, so existing persister runs are unaffected). Same RSL3 multi-dose
    /// regime as `persister_reduces_multidose_kills`.
    #[test]
    fn persister_oxphos_suppression_reduces_rsl3_kill_and_hdac_restores() {
        let cfg = RunConfig {
            grid_dim: 20,
            n_steps: 120,
        };
        let cond = Condition {
            name: "persister_oxphos".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: None,
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::MultiDose {
                dose_steps: (0..120).step_by(4).collect(),
                peak: 1.0,
                half_life_steps: 10.0,
            },
        };
        let run = |p: PersisterConfig| {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    persister: Some(p),
                    ..Default::default()
                },
            )
            .total_dead
        };
        // Persister baseline (OXPHOS axis off, the existing #241 behavior).
        let baseline = run(PersisterConfig::enabled());
        assert!(
            baseline >= 20,
            "baseline RSL3 must kill a meaningful number of cells for this test \
             to be informative; got {baseline}"
        );
        // OXPHOS-low persisters: lower basal ROS ⇒ fewer RSL3 kills.
        let suppressed = run(PersisterConfig {
            oxphos_ros_suppression: 0.7,
            ..PersisterConfig::enabled()
        });
        assert!(
            suppressed < baseline,
            "OXPHOS-ROS suppression must reduce RSL3 kills of persisters: \
             baseline={baseline}, suppressed={suppressed}"
        );
        // Full HDAC inhibitor reverses the suppression (multiplier back to 1.0)
        // ⇒ reproduces the persister baseline EXACTLY.
        let hdac_rescued = run(PersisterConfig {
            oxphos_ros_suppression: 0.7,
            hdac_inhibitor: 1.0,
            ..PersisterConfig::enabled()
        });
        assert_eq!(
            hdac_rescued, baseline,
            "a full HDAC inhibitor must restore RSL3 kill to the persister baseline"
        );
        // suppression = 0 is the layer-off identity behind the matrix byte-identity.
        let zero = run(PersisterConfig {
            oxphos_ros_suppression: 0.0,
            ..PersisterConfig::enabled()
        });
        assert_eq!(
            zero, baseline,
            "oxphos_ros_suppression=0 must reproduce the persister baseline kills"
        );
    }

    /// #470: lock the `--snapshot=persister-oxphos` preset -> Overrides wiring.
    #[test]
    fn persister_oxphos_snapshot_preset_is_wired() {
        let p = resolve_snapshot("persister-oxphos");
        assert_eq!(p.name, "persister-oxphos");
        assert!(
            p.persister_oxphos,
            "the persister-oxphos preset must enable the OXPHOS-ROS axis"
        );
        assert!(
            p.persister,
            "OXPHOS-ROS suppression requires the persister model to also be on"
        );
        assert_eq!(
            p.treatment,
            Treatment::RSL3,
            "persister-oxphos uses RSL3 (the covalent-knockdown path the ROS supply governs)"
        );
    }

    /// Reversion (#241): once dosing truly stops, the epigenetic persister
    /// mark decays. This exercises the `revert` branch in the per-step loop
    /// (not just `acquire`), which requires `drug_intensity == 0.0`. A
    /// `MultiDose` schedule cannot do that — its `factor_at` is a sum of
    /// exponentially-decaying bolus tails that is strictly positive at every
    /// step after the first dose, so `acquire` would fire every step and
    /// `revert` would never run (review finding). `Infusion` returns a literal
    /// `0.0` outside its window, so the `early_window` run takes `revert` for
    /// the entire drug-free back half and ends well below the `sustained` run.
    #[test]
    fn persister_reverts_after_dosing_stops() {
        let cfg = RunConfig {
            grid_dim: 20,
            n_steps: 120,
        };
        let mk = |name: &str, end: u32| Condition {
            name: name.to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: None,
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            // Infusion: factor_at == level inside [start, end), literal 0.0
            // outside. The 0.0 forces the `revert` branch after `end`.
            dose_schedule: DoseSchedule::Infusion {
                start: 0,
                end,
                level: 1.0,
            },
        };
        let sustained = mk("sustained", 120); // drug present every step → only acquire
        let early_window = mk("early_window", 60); // drug-free steps 60..119 → revert fires
        let pm_sustained = run_one_condition_full(
            &sustained,
            cfg,
            None,
            Overrides {
                persister: Some(PersisterConfig::enabled()),
                ..Default::default()
            },
        )
        .persister_mean
        .unwrap();
        let pm_early = run_one_condition_full(
            &early_window,
            cfg,
            None,
            Overrides {
                persister: Some(PersisterConfig::enabled()),
                ..Default::default()
            },
        )
        .persister_mean
        .unwrap();
        // 60 steps of exponential reversion (rate 0.01/step → ×0.55 over 60)
        // leave the early-window mean materially below the sustained mean.
        assert!(
            pm_sustained > 0.0 && pm_early >= 0.0,
            "both runs should acquire some persistence: sustained={pm_sustained:.3}, early={pm_early:.3}"
        );
        assert!(
            pm_early < pm_sustained * 0.85,
            "persister fraction must materially revert after dosing stops: \
             early-window mean={pm_early:.3} should be well below sustained mean={pm_sustained:.3}"
        );
    }

    /// #342: the reversible-to-irreversible epigenetic-locking transition is
    /// dose-schedule dependent, and this exercises its spatial sim-tme-3d wiring
    /// (`step_with_locking` per cell; the library dynamics are unit-tested in
    /// ferroptosis-core::persister). Under SUSTAINED (continuous) drug the
    /// #377: non-drug stress-niche persister entry. Under Control (ZERO drug) on
    /// an O2-gradient grid (which has a hypoxic core), enabling the stress-entry
    /// term raises the OVERALL `persister_mean` above the stress-off baseline,
    /// which stays exactly 0 (no drug ⇒ no drug-driven entry). The rise is
    /// stress-driven and concentrated where `1 - o2_supply` is large (the hypoxic
    /// core), but this test asserts on the overall mean, not a zone-resolved
    /// split; the explicit hypoxic-vs-normoxic comparison (stress 0.9 vs 0.1) is
    /// the library test `stress_entry_is_noop_by_default_and_raises_reversible_under_stress`.
    /// This isolates the NON-DRUG entry route the issue adds. `stress_entry_rate
    /// = 0` is the byte-identical default (the production matrix never enters the
    /// persister path anyway, persister = None).
    #[test]
    fn stress_niche_drives_nondrug_persister_entry() {
        let rcfg = RunConfig {
            grid_dim: 24,
            n_steps: 120,
        };
        let cond = Condition {
            name: "stress_persister".to_string(),
            treatment: Treatment::Control, // ZERO drug: isolates non-drug stress entry
            treatment_name: "Control".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA), // O2 gradient ⇒ a hypoxic core
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let run = |pcfg: PersisterConfig| {
            run_one_condition_full(
                &cond,
                rcfg,
                None,
                Overrides {
                    persister: Some(pcfg),
                    ..Default::default()
                },
            )
            .persister_mean
            .unwrap_or(0.0)
        };
        let off = PersisterConfig::enabled(); // stress_entry_rate = 0
        let on = PersisterConfig {
            stress_entry_rate: 0.05,
            ..PersisterConfig::enabled()
        };
        let pm_off = run(off);
        let pm_on = run(on);
        // No drug + stress off ⇒ persister fraction stays exactly 0.
        assert_eq!(
            pm_off, 0.0,
            "no-drug, stress-off persister fraction must be 0; got {pm_off}"
        );
        // Stress entry drives non-drug persister entry in the hypoxic zones.
        assert!(
            pm_on > pm_off,
            "stress entry must raise the non-drug persister fraction: on={pm_on}, off={pm_off}"
        );
        assert!(
            pm_on > 0.001,
            "stress-driven persister fraction should be observable: {pm_on}"
        );
        // Deterministic (the stress signal is the geometric O2 field).
        assert_eq!(pm_on, run(on));
    }

    /// per-cell sustained-exposure EMA crosses `lock_threshold` and a fraction of
    /// the persister pool ratchets into the non-reverting `locked` sub-pool;
    /// under INTERMITTENT dosing the EMA decays in the drug-off gaps and never
    /// crosses, so nothing locks.
    #[test]
    fn continuous_dosing_locks_persisters_but_intermittent_does_not() {
        let cfg = RunConfig {
            grid_dim: 20,
            n_steps: 160,
        };
        // Lock-enabled config: under continuous drug=1.0 the EMA steady state is
        // 1.0/exposure_decay = 10 (> lock_threshold 5); under the ~25%-duty
        // intermittent schedule the EMA settles well below 5.
        let lock_cfg = PersisterConfig {
            lock_rate: 0.1,
            lock_threshold: 5.0,
            exposure_decay: 0.1,
            ..PersisterConfig::enabled()
        };
        let base = |name: &str, schedule: DoseSchedule| Condition {
            name: name.to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: None,
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: schedule,
        };
        let run = |cond: &Condition, pcfg: PersisterConfig| {
            run_one_condition_full(
                cond,
                cfg,
                None,
                Overrides {
                    persister: Some(pcfg),
                    ..Default::default()
                },
            )
        };
        // Continuous: drug = 1.0 every step (Infusion spanning the whole run).
        let continuous = base(
            "continuous_lock",
            DoseSchedule::Infusion {
                start: 0,
                end: cfg.n_steps,
                level: 1.0,
            },
        );
        // Intermittent: short sharp pulses every 8 steps (low duty, EMA stays low).
        let intermittent = base(
            "intermittent_lock",
            DoseSchedule::MultiDose {
                dose_steps: (0..cfg.n_steps).step_by(8).collect(),
                peak: 1.0,
                half_life_steps: 1.0,
            },
        );
        let cont_locked = run(&continuous, lock_cfg)
            .persister_locked_mean
            .expect("persister on reports a locked mean");
        let interm_locked = run(&intermittent, lock_cfg)
            .persister_locked_mean
            .expect("persister on reports a locked mean");
        // Continuous dosing locks a non-negligible irreversible pool.
        assert!(
            cont_locked > 0.01,
            "sustained dosing must lock a fraction of persisters: locked_mean={cont_locked:.4}"
        );
        // Intermittent dosing keeps the exposure EMA below threshold ⇒ no locking.
        assert!(
            interm_locked < 1e-4,
            "intermittent dosing must NOT lock (EMA stays below threshold): \
             locked_mean={interm_locked:.4}"
        );
        assert!(
            cont_locked > interm_locked,
            "continuous locks more than intermittent: {cont_locked:.4} vs {interm_locked:.4}"
        );
        // Off-by-default: with lock_rate=0 (enabled()), even continuous dosing
        // locks nothing, so the layer stays byte-identical to the pre-#342
        // persister (the locked sub-pool is exactly 0).
        assert_eq!(
            run(&continuous, PersisterConfig::enabled()).persister_locked_mean,
            Some(0.0),
            "lock_rate=0 must lock nothing even under continuous dosing"
        );
        // Deterministic.
        assert_eq!(
            cont_locked,
            run(&continuous, lock_cfg).persister_locked_mean.unwrap(),
            "locking is deterministic"
        );
    }

    /// #365: wiring the NCOA4-ferritinophagy ramp (#338/#340) + hypoxia-driven
    /// iron import (#340) into the spatial model RAISES the Fenton-iron
    /// contribution to RSL3 ferroptosis, qualifying the §7.1 "hypoxia uniformly
    /// protects RSL3" framing. Both knobs `0.0` ⇒ byte-identical to the
    /// no-override run. With them on, RSL3 kill rises several-fold, and the rise
    /// reaches the DEEP hypoxic core too, not only the oxygenated rim. That
    /// deep-core rescue is a MODEL ARTIFACT (flagged in §7.1): the model's Fenton
    /// ROS term is O2-independent (`iron × fenton_rate`, added to total ROS
    /// independently of basal ROS; the model modulates only `basal_ros`, NOT the
    /// Fenton rate), so the iron boost, largest where O2 is lowest, raises core
    /// Fenton ROS regardless of anoxia. A real anoxic core would also lose its
    /// superoxide/SOD-derived H2O2 substrate, which the model does not couple to
    /// O2, so the model OVERSTATES the deep-core rescue; the robust result is the
    /// direction, not the deep-core magnitude. This pins the §7.1 headline config
    /// (grid 60, λ=50 µm, 180 steps, both knobs at 2.0).
    #[test]
    fn hypoxia_iron_raises_rsl3_kill_including_the_deep_core() {
        let cfg = RunConfig {
            grid_dim: 60,
            n_steps: 180,
        };
        let cond = Condition {
            name: "rsl3_dyn_iron".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(50.0),
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let run = |fer: f64, hyp: f64| {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    ferritinophagy_release: fer,
                    hypoxia_iron_sensitivity: hyp,
                    ..Default::default()
                },
            )
        };
        let base = run(0.0, 0.0);
        // Knob-off identity: both 0.0 reproduces the no-override run exactly.
        let plain = run_one_condition_with_config(&cond, cfg, None);
        assert_eq!(
            base.overall_kill_rate, plain.overall_kill_rate,
            "ferritinophagy_release=0 + hypoxia_iron_sensitivity=0 must be byte-identical"
        );
        let both = run(2.0, 2.0); // the §7.1 headline knob values
        assert!(
            base.overall_kill_rate > 0.0,
            "RSL3 must produce some baseline kill to amplify; got {}",
            base.overall_kill_rate
        );
        // Overall kill rises several-fold (the §7.1 headline ~0.5% → ~2.9%).
        assert!(
            both.overall_kill_rate > base.overall_kill_rate * 3.0,
            "dynamic iron must materially raise overall RSL3 kill: base={:.4}, both={:.4}",
            base.overall_kill_rate,
            both.overall_kill_rate
        );
        // The DEEP hypoxic-zone kill rises from ~0 (the §7.1 headline 0% → ~1.3%):
        // the model's O2-independent Fenton means the iron boost rescues the core
        // (a model artifact, since a real anoxic core would lose its H2O2 source).
        assert!(
            both.hypoxic_kill_rate > base.hypoxic_kill_rate + 0.005,
            "the deep hypoxic core is rescued in the model (O2-independent Fenton): \
             base={:.4}, both={:.4}",
            base.hypoxic_kill_rate,
            both.hypoxic_kill_rate
        );
        // Deterministic (the iron scaling is a one-time geometric setup mutation).
        assert_eq!(both.overall_kill_rate, run(2.0, 2.0).overall_kill_rate);
    }

    /// #383: O2-coupling the Fenton H2O2 substrate REVERSES the #365 deep-core
    /// "rescue" artifact. With hypoxia-iron ON (the §7.1 artifact regime, where
    /// the O2-independent Fenton lets the iron boost kill the anoxic core),
    /// turning on `fenton_o2_dependence` scales the effective Fenton flux DOWN
    /// where O2 is low — so the deep hypoxic core is protected again, the
    /// biologically correct behavior (a real anoxic core loses its
    /// superoxide/SOD-derived H2O2). `fenton_o2_dependence = 0` reproduces the
    /// #365 run exactly (byte-identical). Same §7.1 headline config as the #365
    /// test (grid 60, λ=50 µm, 180 steps, iron knobs at 2.0) for direct
    /// comparability.
    #[test]
    fn o2_dependent_fenton_protects_the_anoxic_core_under_hypoxia_iron() {
        let cfg = RunConfig {
            grid_dim: 60,
            n_steps: 180,
        };
        let cond = Condition {
            name: "rsl3_o2_fenton".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(50.0),
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        // All runs share the #365 hypoxia-iron regime (the artifact); only the
        // new #383 Fenton-O2 dependence knob varies.
        let run = |fenton_o2: f64| {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    ferritinophagy_release: 2.0,
                    hypoxia_iron_sensitivity: 2.0,
                    fenton_o2_dependence: fenton_o2,
                    ..Default::default()
                },
            )
        };
        // Byte-identity: fenton_o2_dependence=0 reproduces the #365 hypoxia-iron
        // run exactly (the factor is 1.0 ⇒ `cell.iron` untouched).
        let artifact = run(0.0);
        let artifact_ref = run_one_condition_full(
            &cond,
            cfg,
            None,
            Overrides {
                ferritinophagy_release: 2.0,
                hypoxia_iron_sensitivity: 2.0,
                ..Default::default()
            },
        );
        assert_eq!(
            artifact.overall_kill_rate, artifact_ref.overall_kill_rate,
            "fenton_o2_dependence=0 must be byte-identical to the #365 run"
        );
        assert_eq!(
            artifact.hypoxic_kill_rate, artifact_ref.hypoxic_kill_rate,
            "fenton_o2_dependence=0 must be byte-identical in the hypoxic zone too"
        );
        // The artifact regime DOES rescue the deep core (precondition: the thing
        // we are correcting actually happens here).
        assert!(
            artifact.hypoxic_kill_rate > 0.005,
            "precondition: the #365 hypoxia-iron regime should rescue the deep core \
             (hypoxic_kill_rate={:.4}) so there is an artifact to reverse",
            artifact.hypoxic_kill_rate
        );
        // Fully O2-gated Fenton (dependence=1.0): the anoxic core loses its H2O2
        // substrate, so the deep-core rescue is REVERSED — hypoxic kill drops back
        // down toward the un-rescued baseline.
        let o2_gated = run(1.0);
        assert!(
            o2_gated.hypoxic_kill_rate < artifact.hypoxic_kill_rate,
            "O2-gated Fenton must protect the anoxic core relative to the artifact: \
             artifact={:.4}, o2_gated={:.4}",
            artifact.hypoxic_kill_rate,
            o2_gated.hypoxic_kill_rate
        );
        // The protection should be substantial in the deep core (the H2O2
        // substrate is nearly gone where O2 ≈ 0), recovering most of the artifact.
        assert!(
            o2_gated.hypoxic_kill_rate < artifact.hypoxic_kill_rate * 0.5,
            "O2-gated Fenton should recover most of the deep-core rescue: \
             artifact={:.4}, o2_gated={:.4}",
            artifact.hypoxic_kill_rate,
            o2_gated.hypoxic_kill_rate
        );
        // RIM-vs-CORE DIFFERENTIAL — the actual scientific content of the gate
        // (the manuscript §7.1 claim that the iron rise can still raise kill in
        // the oxygenated rim/mid-zone where some H2O2 remains, while the anoxic
        // core is starved). The O2 gate must suppress the hypoxic core's kill
        // FRACTIONALLY more than the oxygenated rim's, because the rim keeps its
        // O2-derived H2O2 (o2_supply ≈ 1 ⇒ factor ≈ 1) while the core loses it
        // (o2_supply ≈ 0 ⇒ factor ≈ 0). `frac` is the fractional reduction;
        // a zero baseline (no kill to reduce) maps to 0 so the comparison is safe.
        let frac = |before: f64, after: f64| {
            if before > 0.0 {
                (before - after) / before
            } else {
                0.0
            }
        };
        let core_frac = frac(artifact.hypoxic_kill_rate, o2_gated.hypoxic_kill_rate);
        let rim_frac = frac(artifact.normoxic_kill_rate, o2_gated.normoxic_kill_rate);
        assert!(
            core_frac > rim_frac,
            "the O2 gate must suppress the anoxic core MORE than the oxygenated rim: \
             core_frac={core_frac:.4}, rim_frac={rim_frac:.4} \
             (rim {:.4}->{:.4}, core {:.4}->{:.4})",
            artifact.normoxic_kill_rate,
            o2_gated.normoxic_kill_rate,
            artifact.hypoxic_kill_rate,
            o2_gated.hypoxic_kill_rate
        );
        // Deterministic (the iron scaling is a one-time geometric setup mutation
        // for the static pool, plus a per-step factor on the diffused pool).
        assert_eq!(o2_gated.hypoxic_kill_rate, run(1.0).hypoxic_kill_rate);
    }

    /// #383 (review): with the #365 hypoxia-iron regime OFF, the O2-dependent
    /// Fenton gate is well-behaved on its own — it can only REDUCE or hold the
    /// kill (it scales the Fenton iron substrate DOWN in hypoxic zones, never up),
    /// and `fenton_o2_dependence = 0` is byte-identical to the no-override run.
    /// A fast 30³ × 90 RSL3 config (no dynamic-iron knobs) isolates the gate.
    #[test]
    fn fenton_o2_gate_only_reduces_kill_and_is_identity_at_zero() {
        let cfg = RunConfig {
            grid_dim: 30,
            n_steps: 90,
        };
        let cond = Condition {
            name: "rsl3_fenton_o2_only".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let run = |fenton_o2: f64| {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    fenton_o2_dependence: fenton_o2,
                    ..Default::default()
                },
            )
        };
        // Identity: gate off reproduces the no-override run exactly.
        let off = run(0.0);
        let plain = run_one_condition_with_config(&cond, cfg, None);
        assert_eq!(
            off.overall_kill_rate, plain.overall_kill_rate,
            "fenton_o2_dependence=0 must be byte-identical with no other knobs"
        );
        // Monotone: gating the Fenton substrate down can only reduce or hold kill,
        // never raise it (even with no compensating hypoxia-iron rise).
        let gated = run(1.0);
        assert!(
            gated.overall_kill_rate <= off.overall_kill_rate + 1e-12,
            "the O2-Fenton gate must not raise kill on its own: off={:.4}, gated={:.4}",
            off.overall_kill_rate,
            gated.overall_kill_rate
        );
    }

    // ===== T-cell exhaustion (#243, Phase 1) =====

    /// Headline #243 result: enabling T-cell exhaustion suppresses the TOTAL
    /// immune kills relative to the no-exhaustion baseline. Sustained killing
    /// in a region drives local T cells toward dysfunction
    /// (`1/(1+rate·cumulative_neighborhood_kills)`), so later kills there get
    /// progressively rarer — the "cold tumor" emergence (Wherry 2011; Snell
    /// 2018). Uses SDT + immune (SDT's dense ferroptotic death builds the
    /// DAMP field that drives immune killing); the off run is the byte-identical
    /// `for_3d()` config, the on run overrides only `exhaustion_rate`.
    #[test]
    fn exhaustion_reduces_immune_kills() {
        // 30³ × 130: immune kills only fire after IMMUNE_START_STEP=60 once
        // DAMP crosses threshold, so a large grid + long run is needed for a
        // baseline kill count with headroom to show suppression.
        let cfg = RunConfig {
            grid_dim: 30,
            n_steps: 130,
        };
        let cond = Condition {
            name: "exhaustion_demo".to_string(),
            treatment: Treatment::SDT,
            treatment_name: "SDT".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        // Dense-killing regime: exhaustion only bites when kills CLUSTER (so
        // a cell's neighborhood has prior kills before its own check). At the
        // default immune_kill_rate (0.02) kills are too sparse to overlap, so
        // both runs share a boosted rate and differ ONLY in exhaustion_rate —
        // a fair A/B isolating the exhaustion effect. (Defaults stay
        // byte-identical; that path is covered by the golden + #253 tests.)
        let base = SpatialImmuneConfig {
            immune_kill_rate: 0.5,
            ..SpatialImmuneConfig::for_3d()
        };
        let exhausted = SpatialImmuneConfig {
            exhaustion_rate: 2.0,
            ..base
        };
        let off = run_one_condition_full(
            &cond,
            cfg,
            None,
            Overrides {
                immune: Some(base),
                ..Default::default()
            },
        );
        let on = run_one_condition_full(
            &cond,
            cfg,
            None,
            Overrides {
                immune: Some(exhausted),
                ..Default::default()
            },
        );
        let off_im = off.immune_kills.expect("immune_on populates immune_kills");
        let on_im = on.immune_kills.expect("immune_on populates immune_kills");
        assert!(
            off_im >= 50,
            "baseline must produce enough clustered immune kills for the test \
             to be informative; got {off_im} (adjust grid/steps/rate if the model changed)"
        );
        // Material reduction, not an off-by-one fluke. Observed ≈20%
        // (off=174, on=139) at the time of writing.
        let reduction = (off_im - on_im) as f64 / off_im as f64;
        assert!(
            on_im < off_im && reduction > 0.1,
            "T-cell exhaustion must materially suppress immune kills (>10%): \
             on={on_im}, off={off_im}, reduction={:.1}%",
            reduction * 100.0
        );
        // NOTE: ferroptosis kills are NOT held fixed. Exhaustion only gates the
        // immune-kill loop directly, but a cell spared an (apoptotic) immune
        // kill can instead die ferroptotically later — and that death releases
        // iron that couples to neighbors. So sparing immune kills shifts a few
        // deaths into the ferroptosis tally; that cross-coupling is expected,
        // not a bug, which is why this asserts on immune kills specifically.
    }

    // ===== Treg/MDSC suppressor field (#264 Phase 2) =====

    /// Headline #264 Phase 2 result: under anti-PD-1, **depleting Tregs**
    /// (turning the suppressor field off) recovers immune kills that the
    /// suppressor was damping — i.e. anti-PD-1 alone (Tregs present) is less
    /// effective than anti-PD-1 + Treg depletion (Tauriello 2018). Both arms
    /// share the SAME anti-PD-1 immune config + dense kill regime and differ
    /// ONLY in whether the suppressor field is present. (Defaults stay
    /// byte-identical — the disabled/None path is the matrix path.)
    #[test]
    fn treg_depletion_improves_anti_pd1_kills() {
        let cfg = RunConfig {
            grid_dim: 30,
            n_steps: 130,
        };
        let cond = Condition {
            name: "suppressor_demo".to_string(),
            treatment: Treatment::SDT,
            treatment_name: "SDT".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        // Anti-PD-1 immune config, dense regime (so there are enough immune
        // kills for the suppressor to measurably damp). Both arms share it.
        let immune = SpatialImmuneConfig {
            immune_kill_rate: 0.5,
            anti_pd1_efficacy: 0.5,
            ..SpatialImmuneConfig::for_3d()
        };
        // Tregs present: anti-PD-1 + suppressor on.
        let with_tregs = run_one_condition_full(
            &cond,
            cfg,
            None,
            Overrides {
                immune: Some(immune),
                suppressor: Some(SuppressorConfig::enabled()),
                ..Default::default()
            },
        );
        // Treg-depleted: anti-PD-1 alone (suppressor off).
        let depleted = run_one_condition_full(
            &cond,
            cfg,
            None,
            Overrides {
                immune: Some(immune),
                ..Default::default()
            },
        );
        let with_im = with_tregs
            .immune_kills
            .expect("immune_on populates immune_kills");
        let dep_im = depleted
            .immune_kills
            .expect("immune_on populates immune_kills");
        assert!(
            dep_im >= 50,
            "Treg-depleted baseline must produce enough immune kills to be \
             informative; got {dep_im}"
        );
        // Treg depletion materially recovers immune kills the suppressor damped.
        let recovery = (dep_im - with_im) as f64 / dep_im as f64;
        assert!(
            with_im < dep_im && recovery > 0.1,
            "Treg depletion must materially raise anti-PD-1 immune kills (>10%): \
             with_tregs={with_im}, depleted={dep_im}, recovery={:.1}%",
            recovery * 100.0
        );
        // The suppressor run reports its niche census.
        assert!(
            with_tregs.suppressor_source_count.unwrap_or(0) > 0,
            "suppressor run reports a non-empty niche source count"
        );
        assert!(
            depleted.suppressor_source_count.is_none(),
            "Treg-depleted run omits the suppressor census"
        );
    }

    /// #264 review #1: exercise the **perivascular** seeding branch end-to-end
    /// (suppressor + vasculature together — `vessels` is `Some`, so the niches
    /// are placed at perivascular positions rather than heuristic patches). The
    /// validation test above uses heuristic seeding (no vessels); this confirms
    /// the perivascular path runs through a full simulation and is deterministic.
    #[test]
    fn suppressor_perivascular_runs_with_vasculature() {
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 90,
        };
        let cond = Condition {
            name: "supp_perivascular".to_string(),
            treatment: Treatment::SDT,
            treatment_name: "SDT".to_string(),
            // Vasculature needs o2_lambda to build the supply field + vessels.
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let immune = SpatialImmuneConfig {
            immune_kill_rate: 0.5,
            ..SpatialImmuneConfig::for_3d()
        };
        let run = || {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    immune: Some(immune),
                    vasculature: Some(VasculatureConfig::well_vascularized()),
                    suppressor: Some(SuppressorConfig::enabled()),
                    ..Default::default()
                },
            )
        };
        let a = run();
        let b = run();
        // Perivascular niches were placed (the vessels-as-seed-points branch).
        assert!(
            a.suppressor_source_count.unwrap_or(0) > 0,
            "perivascular seeding marks niche cells near vessels"
        );
        // Deterministic across runs (independent fixed seeds throughout).
        assert_eq!(
            a.total_dead, b.total_dead,
            "perivascular run is deterministic"
        );
        assert_eq!(a.suppressor_source_count, b.suppressor_source_count);
        assert!(
            a.immune_kills.unwrap_or(0) > 0,
            "immune kills still fire under perivascular suppression"
        );
    }

    // ===== Multi-checkpoint brake (#264 Phase 3) =====

    /// Headline #264 Phase 3 result: on a PD-1 + CTLA-4 tumor, **dual blockade**
    /// (anti-PD-1 + anti-CTLA-4) out-kills anti-PD-1 **alone** — anti-PD-1 leaves
    /// CTLA-4 braking, so lifting both raises immune killing (the combination-
    /// immunotherapy result, Sharma & Allison 2015). Both arms share the SAME
    /// dense immune config + checkpoint panel and differ ONLY in whether
    /// anti-CTLA-4 is applied. (Defaults stay byte-identical — no panel ⇒ the
    /// single-PD-1 `effective_brake`.)
    #[test]
    fn dual_checkpoint_blockade_outkills_anti_pd1_alone() {
        let cfg = RunConfig {
            grid_dim: 30,
            n_steps: 130,
        };
        let cond = Condition {
            name: "checkpoint_demo".to_string(),
            treatment: Treatment::SDT,
            treatment_name: "SDT".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        // Dense regime so there are enough immune kills for the brake difference
        // to register. The panel REPLACES the single PD-1 brake, so the immune
        // config's pd1_brake/anti_pd1 are unused here.
        let immune = SpatialImmuneConfig {
            immune_kill_rate: 0.5,
            ..SpatialImmuneConfig::for_3d()
        };
        let tumor = CheckpointPanel::pd1_ctla4_tumor();
        let run = |panel: CheckpointPanel| {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    immune: Some(immune),
                    checkpoints: Some(panel),
                    ..Default::default()
                },
            )
        };
        let mono = run(tumor.with_anti_pd1(0.8)); // CTLA-4 still braking
        let combo = run(tumor.with_anti_pd1(0.8).with_anti_ctla4(0.8)); // both lifted
        let mono_im = mono.immune_kills.expect("immune_on populates immune_kills");
        let combo_im = combo
            .immune_kills
            .expect("immune_on populates immune_kills");
        assert!(
            mono_im >= 50,
            "anti-PD-1 monotherapy must produce enough immune kills to be \
             informative; got {mono_im}"
        );
        // Dual blockade materially out-kills the monotherapy.
        let gain = (combo_im - mono_im) as f64 / mono_im as f64;
        assert!(
            combo_im > mono_im && gain > 0.1,
            "dual checkpoint blockade must materially out-kill anti-PD-1 alone (>10%): \
             mono={mono_im}, combo={combo_im}, gain={:.1}%",
            gain * 100.0
        );
        // The combo's combined brake is the lower of the two (reported).
        assert!(
            combo.checkpoint_brake.unwrap() < mono.checkpoint_brake.unwrap(),
            "combo brake {:?} should be below mono brake {:?}",
            combo.checkpoint_brake,
            mono.checkpoint_brake
        );
    }

    /// #264 Phase 4: a cDC1-poor tumor (the literature DC-subset mix) primes
    /// anti-tumor CD8 killing less efficiently than a balanced/cDC1-rich one, so
    /// it produces FEWER immune kills under the same DAMP signal. A/B with the DC
    /// subset mix as the only difference; deterministic.
    #[test]
    fn cdc1_poor_dc_subsets_reduce_immune_kills() {
        let cfg = RunConfig {
            grid_dim: 30,
            n_steps: 130,
        };
        let cond = Condition {
            name: "dc_subset_demo".to_string(),
            treatment: Treatment::SDT,
            treatment_name: "SDT".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        // Dense regime (boosted kill rate, PD-1 brake lifted) so there are enough
        // immune kills for the priming-efficiency difference to register.
        let immune = SpatialImmuneConfig {
            immune_kill_rate: 0.5,
            anti_pd1_efficacy: 1.0,
            ..SpatialImmuneConfig::for_3d()
        };
        // Baseline: no DC-subset layer (full priming). cDC1-poor: the literature
        // mix (priming efficiency 0.37).
        let baseline = run_one_condition_full(
            &cond,
            cfg,
            None,
            Overrides {
                immune: Some(immune),
                ..Default::default()
            },
        );
        let run_cdc1_poor = || {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    immune: Some(immune),
                    dc_subsets: Some(DcSubsetConfig::literature()),
                    ..Default::default()
                },
            )
        };
        let cdc1_poor = run_cdc1_poor();
        let base_im = baseline
            .immune_kills
            .expect("immune_on populates immune_kills");
        let poor_im = cdc1_poor
            .immune_kills
            .expect("immune_on populates immune_kills");
        assert!(
            base_im >= 50,
            "baseline must produce enough immune kills to be informative; got {base_im}"
        );
        assert!(
            poor_im < base_im,
            "a cDC1-poor tumor must prime fewer immune kills than the balanced \
             baseline: baseline={base_im}, cDC1-poor={poor_im}"
        );
        // Deterministic (the priming scalar is uniform, no extra RNG).
        assert_eq!(poor_im, run_cdc1_poor().immune_kills.unwrap());
        // NB: ferroptosis kills are NOT held fixed. DC subsets only gate the
        // immune-kill loop directly (immune kills are apoptotic: no DAMP/iron
        // release), but a cell spared an immune kill can instead die
        // ferroptotically later, and THAT death releases iron that couples to
        // neighbors. So reducing immune kills shifts a few deaths into the
        // ferroptosis tally; that cross-coupling is expected, which is why this
        // asserts on immune kills specifically.
    }

    /// #337: the immunosuppressive arm of ferroptosis (extracellular GPX4 /
    /// oxidized-lipid DC suppression) scales immune kill DOWN as the local
    /// ferroptotic-death/DAMP signal rises, so enabling it produces FEWER net
    /// immune kills than the pro-immune-only baseline. A/B with the suppression
    /// strength as the only difference; deterministic. strength=0 reproduces the
    /// baseline exactly (byte-identical layer-off).
    #[test]
    fn immunosuppressive_ferroptosis_reduces_immune_kills() {
        let cfg = RunConfig {
            grid_dim: 30,
            n_steps: 130,
        };
        let cond = Condition {
            name: "ferro_immunosuppression_demo".to_string(),
            treatment: Treatment::SDT,
            treatment_name: "SDT".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        // Dense regime (boosted kill rate, PD-1 brake lifted) so there are enough
        // immune kills for the suppression to register.
        let base_immune = SpatialImmuneConfig {
            immune_kill_rate: 0.5,
            anti_pd1_efficacy: 1.0,
            ..SpatialImmuneConfig::for_3d()
        };
        let run = |im: SpatialImmuneConfig| {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    immune: Some(im),
                    ..Default::default()
                },
            )
        };
        let baseline = run(base_immune);
        let suppressed = run(SpatialImmuneConfig {
            ferro_immunosuppression_strength: 0.05,
            ..base_immune
        });
        let base_im = baseline
            .immune_kills
            .expect("immune_on populates immune_kills");
        let supp_im = suppressed
            .immune_kills
            .expect("immune_on populates immune_kills");
        assert!(
            base_im >= 50,
            "baseline must produce enough immune kills to be informative; got {base_im}"
        );
        assert!(
            supp_im < base_im,
            "immunosuppressive ferroptosis must reduce net immune kills: \
             baseline={base_im}, suppressed={supp_im}"
        );
        // Deterministic (the multiplier is a pure function of local DAMP).
        assert_eq!(
            supp_im,
            run(SpatialImmuneConfig {
                ferro_immunosuppression_strength: 0.05,
                ..base_immune
            })
            .immune_kills
            .unwrap()
        );
        // strength = 0 reproduces the baseline exactly (the layer-off invariant
        // behind the production-matrix byte-identity).
        let zero = run(SpatialImmuneConfig {
            ferro_immunosuppression_strength: 0.0,
            ..base_immune
        });
        assert_eq!(
            zero.immune_kills.unwrap(),
            base_im,
            "ferro_immunosuppression_strength=0 must reproduce the baseline immune kills"
        );
    }

    /// #264: lock the `--snapshot=dc-subsets` preset -> Overrides wiring.
    #[test]
    fn dc_subsets_snapshot_preset_is_wired() {
        let p = resolve_snapshot("dc-subsets");
        assert_eq!(p.name, "dc-subsets");
        assert!(p.dc_subsets, "the dc-subsets preset must enable the layer");
        assert!(
            p.immune_on,
            "DC subsets only matter with the immune response on"
        );
    }

    /// #469: DC / effector-cell ferroptosis susceptibility. A ferroptotic TME
    /// (the SAME local DAMP/lipid-ROS that drives DC activation) kills the
    /// effector DCs themselves, lowering immune amplification, and PD-L1 /
    /// checkpoint protection blocks that, restoring the baseline. A/B with the
    /// DC-ferroptosis config as the only difference; deterministic. The fully
    /// PD-L1-protected case and the susceptibility=0 case both reproduce the
    /// baseline exactly (the byte-identity invariant behind the production
    /// matrix). Anchor: Yao et al., Cell Reports 2024, PMID 39423128.
    #[test]
    fn dc_ferroptosis_reduces_immune_kills_and_pdl1_protects() {
        let cfg = RunConfig {
            grid_dim: 30,
            n_steps: 130,
        };
        let cond = Condition {
            name: "dc_ferroptosis_demo".to_string(),
            treatment: Treatment::SDT,
            treatment_name: "SDT".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        // Dense regime (boosted kill rate, PD-1 brake lifted) so there are
        // enough immune kills for the DC-ferroptosis reduction to register.
        let base_immune = SpatialImmuneConfig {
            immune_kill_rate: 0.5,
            anti_pd1_efficacy: 1.0,
            ..SpatialImmuneConfig::for_3d()
        };
        let run = |dc: Option<DcFerroptosisConfig>| {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    immune: Some(base_immune),
                    dc_ferroptosis: dc,
                    ..Default::default()
                },
            )
            .immune_kills
            .expect("immune_on populates immune_kills")
        };
        let baseline = run(None);
        // PD-L1-low, ferroptosis-vulnerable DCs ⇒ fewer net immune kills.
        let vulnerable = run(Some(DcFerroptosisConfig {
            susceptibility: 1.0,
            pdl1_protection: 0.0,
        }));
        assert!(
            baseline >= 50,
            "baseline must produce enough immune kills to be informative; got {baseline}"
        );
        assert!(
            vulnerable < baseline,
            "DC ferroptosis must reduce net immune kills: baseline={baseline}, vulnerable={vulnerable}"
        );
        // Deterministic (the survival multiplier is a pure function of local DAMP).
        assert_eq!(
            vulnerable,
            run(Some(DcFerroptosisConfig {
                susceptibility: 1.0,
                pdl1_protection: 0.0,
            }))
        );
        // Full PD-L1 protection ⇒ survival 1.0 ⇒ reproduces the baseline exactly
        // (DCs immune to ferroptosis regardless of the ferroptotic TME).
        assert_eq!(
            run(Some(DcFerroptosisConfig {
                susceptibility: 1.0,
                pdl1_protection: 1.0,
            })),
            baseline,
            "full PD-L1 protection must restore the baseline immune kills"
        );
        // susceptibility = 0 is the layer-off identity behind the matrix
        // byte-identity.
        assert_eq!(
            run(Some(DcFerroptosisConfig {
                susceptibility: 0.0,
                pdl1_protection: 0.0,
            })),
            baseline,
            "susceptibility=0 must reproduce the baseline immune kills"
        );
    }

    /// #469: lock the `--snapshot=dc-ferroptosis` preset -> Overrides wiring.
    #[test]
    fn dc_ferroptosis_snapshot_preset_is_wired() {
        let p = resolve_snapshot("dc-ferroptosis");
        assert_eq!(p.name, "dc-ferroptosis");
        assert!(
            p.dc_ferroptosis,
            "the dc-ferroptosis preset must enable the layer"
        );
        assert!(
            p.immune_on,
            "DC ferroptosis only matters with the immune response on"
        );
    }

    /// #341 (SASP→immune coupling): the senescence-associated secretory phenotype
    /// couples senescent cells to the immune layer with a SIGNED multiplier, so
    /// the layer is genuinely bidirectional. Using a SASP-ONLY config (all four
    /// biochem axes `1.0`, only `sasp_immune_mult` varies) isolates the immune
    /// coupling: the grid biochem and the non-senescent cells are byte-identical
    /// across arms, so the same senescent subset is killed MORE under surveillance
    /// (`> 1`, Kang 2011) and LESS under immunosuppressive SASP (`< 1`, Di Mitri
    /// 2014). Deterministic; `sasp_immune_mult=1.0` is identity ⇒ reproduces the
    /// no-senescence baseline exactly (the byte-identity invariant).
    #[test]
    fn sasp_immune_coupling_is_bidirectional_and_identity_is_baseline() {
        let cfg = RunConfig {
            grid_dim: 30,
            n_steps: 130,
        };
        let cond = Condition {
            name: "sasp_demo".to_string(),
            treatment: Treatment::SDT,
            treatment_name: "SDT".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        // Dense regime (boosted kill rate, PD-1 brake lifted) so there are enough
        // immune kills for the SASP multiplier to register.
        let immune = SpatialImmuneConfig {
            immune_kill_rate: 0.5,
            anti_pd1_efficacy: 1.0,
            ..SpatialImmuneConfig::for_3d()
        };
        // A SASP-only senescence config: 40% of tumor cells marked, NO biochem
        // perturbation (all muls 1.0), only the immune coupling varies. This keeps
        // the grid byte-identical across arms so the immune-kill delta is purely
        // the SASP multiplier.
        let sasp_only = |mult: f64| SenescenceConfig {
            fraction: 0.4,
            iron_mul: 1.0,
            gpx4_mul: 1.0,
            nrf2_mul: 1.0,
            fsp1_mul: 1.0,
            sasp_immune_mult: mult,
            sasp_field_strength: 0.0,
        };
        let run = |sen: Option<SenescenceConfig>| {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    immune: Some(immune),
                    senescence: sen,
                    ..Default::default()
                },
            )
            .immune_kills
            .expect("immune_on populates immune_kills")
        };
        let base_im = run(None);
        let surveillance_im = run(Some(sasp_only(1.6))); // anti-tumor surveillance
        let immunosuppress_im = run(Some(sasp_only(0.25))); // immunosuppressive SASP
        assert!(
            base_im >= 50,
            "baseline must produce enough immune kills to be informative; got {base_im}"
        );
        // Immunosuppressive SASP reduces net immune kills below baseline.
        assert!(
            immunosuppress_im < base_im,
            "immunosuppressive SASP (mult<1) must reduce immune kills: \
             baseline={base_im}, immunosuppressed={immunosuppress_im}"
        );
        // Surveillance does not reduce kills (raises or saturates) ...
        assert!(
            surveillance_im >= base_im,
            "surveillance SASP (mult>1) must not reduce immune kills: \
             baseline={base_im}, surveillance={surveillance_im}"
        );
        // ... and the bidirectional spread is real: surveillance > immunosuppression.
        assert!(
            surveillance_im > immunosuppress_im,
            "SASP is bidirectional: surveillance={surveillance_im} must outkill \
             immunosuppression={immunosuppress_im}"
        );
        // Deterministic (the multiplier is a pure function of the fixed mask).
        assert_eq!(immunosuppress_im, run(Some(sasp_only(0.25))));
        // Identity: mult=1.0 with identity biochem is filtered to None ⇒ the
        // no-senescence baseline EXACTLY (the byte-identity invariant in miniature).
        assert_eq!(
            run(Some(sasp_only(1.0))),
            base_im,
            "sasp_immune_mult=1.0 (identity) must reproduce the no-senescence baseline"
        );
    }

    /// #376 (diffusing SASP FIELD): unlike #341's cell-autonomous multiplier, the
    /// field is PARACRINE — it reaches non-senescent NEIGHBORS. Using a
    /// FIELD-only config (all four biochem axes `1.0` AND `sasp_immune_mult` `1.0`,
    /// only `sasp_field_strength` varies) isolates the field from the biochem and
    /// the cell-autonomous coupling, which stay byte-identical across arms.
    ///
    /// Two levels of assertion:
    /// 1. AGGREGATE: the immunosuppressive arm (`> 0`) LOWERS net kills below
    ///    baseline, the surveillance arm (`< 0`) RAISES them, and the two straddle
    ///    the baseline.
    /// 2. BYSTANDER (the #376 acceptance criterion): the per-run
    ///    `nonsenescent_immune_kills` metric — immune kills landing on
    ///    NON-senescent tumor cells specifically — shifts with the field sign
    ///    (surveillance > immunosuppressive). Because the biochem and mask are
    ///    identical across the two field arms, that shift is provably attributable
    ///    to the diffusing field reaching the non-senescent bystander population,
    ///    not just self-acting on the senescent sources.
    ///
    /// Deterministic; `sasp_field_strength=0.0` is identity ⇒ reproduces the
    /// no-senescence baseline EXACTLY (the production-matrix byte-identity).
    #[test]
    fn sasp_field_is_bidirectional_paracrine_and_identity_is_baseline() {
        let cfg = RunConfig {
            grid_dim: 30,
            n_steps: 130,
        };
        let cond = Condition {
            name: "sasp_field_demo".to_string(),
            treatment: Treatment::SDT,
            treatment_name: "SDT".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let immune = SpatialImmuneConfig {
            immune_kill_rate: 0.5,
            anti_pd1_efficacy: 1.0,
            ..SpatialImmuneConfig::for_3d()
        };
        // FIELD-only: 40% senescent, NO biochem perturbation AND no cell-autonomous
        // SASP coupling (`sasp_immune_mult` 1.0); only the diffusing field varies.
        let field_only = |strength: f64| SenescenceConfig {
            fraction: 0.4,
            iron_mul: 1.0,
            gpx4_mul: 1.0,
            nrf2_mul: 1.0,
            fsp1_mul: 1.0,
            sasp_immune_mult: 1.0,
            sasp_field_strength: strength,
        };
        let run_full = |sen: Option<SenescenceConfig>| {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    immune: Some(immune),
                    senescence: sen,
                    ..Default::default()
                },
            )
        };
        let im = |sen: Option<SenescenceConfig>| {
            run_full(sen)
                .immune_kills
                .expect("immune_on populates immune_kills")
        };
        let base_im = im(None);
        let surveillance = run_full(Some(field_only(-0.6))); // surveillance: raises kill
        let immunosuppress = run_full(Some(field_only(0.6))); // immunosuppressive: lowers
        let surveillance_im = surveillance.immune_kills.unwrap();
        let immunosuppress_im = immunosuppress.immune_kills.unwrap();
        assert!(
            base_im >= 50,
            "baseline must produce enough immune kills to be informative; got {base_im}"
        );
        // ---- Level 1: aggregate immune-kill shift, bidirectional. ----
        assert!(
            immunosuppress_im < base_im,
            "immunosuppressive SASP field (strength>0) must reduce net immune kills: \
             baseline={base_im}, suppressed={immunosuppress_im}"
        );
        assert!(
            surveillance_im > base_im,
            "surveillance SASP field (strength<0) must raise net immune kills: \
             baseline={base_im}, surveillance={surveillance_im}"
        );
        assert!(
            surveillance_im > immunosuppress_im,
            "the SASP field is bidirectional: surveillance={surveillance_im} must outkill \
             immunosuppression={immunosuppress_im}"
        );
        // ---- Level 2: the BYSTANDER (non-senescent neighbor) effect (#376). ----
        // The field arms populate `nonsenescent_immune_kills`; the baseline (no
        // field) leaves it `None` (so the production matrix summary.json omits it).
        assert!(
            run_full(None).nonsenescent_immune_kills.is_none(),
            "the no-field baseline must NOT populate nonsenescent_immune_kills"
        );
        let surveillance_bystander = surveillance
            .nonsenescent_immune_kills
            .expect("the field arm populates nonsenescent_immune_kills");
        let immunosuppress_bystander = immunosuppress
            .nonsenescent_immune_kills
            .expect("the field arm populates nonsenescent_immune_kills");
        // Kills on NON-senescent cells shift with the field sign — only possible
        // if the diffusing field reaches the bystander population (the biochem and
        // mask are identical across the two field arms).
        assert!(
            surveillance_bystander > immunosuppress_bystander,
            "the SASP field must reach NON-senescent neighbors: surveillance \
             non-senescent kills={surveillance_bystander} must exceed immunosuppressive \
             non-senescent kills={immunosuppress_bystander}"
        );
        // ---- Determinism + identity. ----
        // Deterministic (the field + multiplier are pure functions of the fixed mask).
        assert_eq!(immunosuppress_im, im(Some(field_only(0.6))));
        // Identity: strength=0.0 with identity biochem + cell-autonomous coupling is
        // filtered to None ⇒ the no-senescence baseline EXACTLY (byte-identity).
        assert_eq!(
            im(Some(field_only(0.0))),
            base_im,
            "sasp_field_strength=0.0 (identity) must reproduce the no-senescence baseline"
        );
    }

    /// #341: lock the `--snapshot=senescence` preset -> Overrides wiring, and
    /// that enabling it emits a `senescent_fraction` metric while the immune
    /// response (needed for the SASP coupling) is on.
    #[test]
    fn senescence_snapshot_preset_is_wired() {
        let p = resolve_snapshot("senescence");
        assert_eq!(p.name, "senescence");
        assert!(p.senescence, "the senescence preset must enable the layer");
        assert!(
            p.immune_on,
            "senescence SASP coupling only matters with the immune response on"
        );
        // The literature config is non-identity (it marks cells + couples SASP).
        assert!(!SenescenceConfig::literature().is_identity());
        // Enabling the layer emits the senescent_fraction metric (~0.2 marked).
        let cfg = RunConfig {
            grid_dim: 20,
            n_steps: 60,
        };
        let cond = Condition {
            name: "sen_metric".to_string(),
            treatment: Treatment::SDT,
            treatment_name: "SDT".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let baseline = run_one_condition_with_config(&cond, cfg, None);
        assert!(
            baseline.senescent_fraction.is_none(),
            "senescence off ⇒ the metric field is omitted"
        );
        let on = run_one_condition_full(
            &cond,
            cfg,
            None,
            Overrides {
                senescence: Some(SenescenceConfig::literature()),
                ..Default::default()
            },
        );
        let frac = on
            .senescent_fraction
            .expect("senescence on ⇒ senescent_fraction reported");
        assert!(
            (0.1..0.35).contains(&frac),
            "literature fraction ~0.2 of tumor cells marked; got {frac}"
        );
    }

    /// #398: lock the `--snapshot=sasp-field` preset wiring. It is the senescence
    /// preset plus the diffusing-SASP-field overlay (a static `sasp_field.npy`):
    /// `resolve_snapshot` must find it, it enables the senescence + immune layers,
    /// and `literature()` carries a non-zero `sasp_field_strength` so the run (and
    /// the overlay) actually exercise the diffusing field rather than just the mask.
    #[test]
    fn sasp_field_snapshot_preset_is_wired() {
        let p = resolve_snapshot("sasp-field");
        assert_eq!(p.name, "sasp-field");
        assert!(
            p.senescence,
            "the sasp-field preset must enable the senescence layer (the SASP source)"
        );
        assert!(
            p.immune_on,
            "the SASP field only couples with the immune response on"
        );
        // The preset's run config is SenescenceConfig::literature(); it must carry a
        // non-zero SASP field strength, else the field is inert and the overlay flat.
        assert!(
            SenescenceConfig::literature().sasp_field_strength != 0.0,
            "the sasp-field preset's literature() config must have a non-zero \
             SASP field strength, otherwise the overlay would be uniformly zero"
        );
    }

    /// #341 review: senescence must compose coherently with clonal repopulation.
    /// A repopulated dead site is a NEW cell grown from a living neighbour, so it
    /// is no longer the senescent cell that died there; its mask bit is cleared
    /// (the analogue of the #302 contact re-application). Therefore enabling
    /// repopulation can only LOWER the end-of-run senescent fraction (clearing
    /// revived sites never adds senescent cells), and with a senolytic config
    /// (senescent cells die under RSL3 and get repopulated) it strictly lowers
    /// it. immune_on=false here so the SASP arm is inert and this isolates the
    /// mask-clear coherence.
    #[test]
    fn senescence_composes_with_clonal_repopulation() {
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 120,
        };
        let cond = Condition {
            name: "sen_clonal".to_string(),
            treatment: Treatment::RSL3, // strong kill to drive turnover + repopulation
            treatment_name: "RSL3".to_string(),
            o2_lambda: None,
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        // Iron-dominant (senolytic) senescence: senescent cells peroxidize MORE
        // under RSL3, so they die and get repopulated, exercising the mask-clear.
        let sen = SenescenceConfig {
            fraction: 0.5,
            iron_mul: 3.0,
            gpx4_mul: 1.0,
            nrf2_mul: 1.0,
            fsp1_mul: 1.0,
            sasp_immune_mult: 1.0,
            sasp_field_strength: 0.0,
        };
        let run = |repop: f64| {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    senescence: Some(sen),
                    clonal: Some(ClonalConfig::literature_4().with_repopulation(repop)),
                    ..Default::default()
                },
            )
            .senescent_fraction
            .expect("senescence on ⇒ senescent_fraction reported")
        };
        let f_repop = run(0.3);
        let f_norepop = run(0.0);
        // Invariant: clearing revived sites can only lower (or hold) the fraction.
        assert!(
            f_repop <= f_norepop + 1e-12,
            "repopulation clears revived senescent sites, so the end-of-run \
             senescent fraction must not exceed the no-repopulation run: \
             repop={f_repop}, norepop={f_norepop}"
        );
        // Observable: the senolytic config kills+repopulates some senescent sites,
        // so the composition is not a silent no-op.
        assert!(
            f_repop < f_norepop,
            "senolytic senescence + repopulation should strictly lower the \
             senescent fraction: repop={f_repop}, norepop={f_norepop}"
        );
        // Deterministic.
        assert_eq!(
            f_repop,
            run(0.3),
            "senescence x clonal-repopulation composition is deterministic"
        );
    }

    // ===== Clonal heterogeneity (#242) =====

    /// K=1 with the identity perturbation must reproduce the no-clonal run
    /// exactly: subclone assignment uses an independent RNG (grid unchanged)
    /// and the perturbation is a no-op. Per-PR complement to the #253
    /// production-SHA guard.
    #[test]
    fn clonal_k1_identity_is_byte_identical() {
        let cfg = RunConfig {
            grid_dim: 20,
            n_steps: 80,
        };
        let cond = Condition {
            name: "clonal_k1".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: None,
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let baseline = run_one_condition_with_config(&cond, cfg, None);
        let k1 = run_one_condition_full(
            &cond,
            cfg,
            None,
            Overrides {
                clonal: Some(ClonalConfig::single_identity()),
                ..Default::default()
            },
        );
        assert_eq!(
            k1.total_dead, baseline.total_dead,
            "K=1 identity must not change kills"
        );
        assert!(
            baseline.subclone_kills.is_none(),
            "no-clonal run omits subclone_kills"
        );
        let sk = k1
            .subclone_kills
            .expect("clonal run reports subclone_kills");
        assert_eq!(sk.len(), 1, "K=1 ⇒ one subclone entry");
        assert_eq!(
            sk[0].total_dead, k1.total_dead,
            "the single subclone holds all kills"
        );
    }

    /// With a 4-subclone literature table, the most-vulnerable subclone (1:
    /// iron-loaded, GPX4-low) must die at a higher rate than the most-resistant
    /// (4: GPX4-high, MUFA-enriched) — the intratumoral-heterogeneity effect.
    /// Uses RSL3 + uniform O₂ so subclonal GPX4/iron differences drive the kill.
    #[test]
    fn clonal_subclones_differ_in_kill_rate() {
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 80,
        };
        let cond = Condition {
            name: "clonal_4".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: None,
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let r = run_one_condition_full(
            &cond,
            cfg,
            None,
            Overrides {
                clonal: Some(ClonalConfig::literature_4()),
                ..Default::default()
            },
        );
        let sk = r.subclone_kills.expect("clonal run reports subclone_kills");
        assert_eq!(sk.len(), 4, "literature_4 ⇒ four subclone entries");
        // Entries are ordered by id; ids 1..=4 span vulnerable→resistant.
        assert_eq!(sk[0].subclone_id, 1);
        assert_eq!(sk[3].subclone_id, 4);
        assert!(
            sk[0].kill_rate > sk[3].kill_rate,
            "vulnerable subclone 1 must out-die resistant subclone 4: \
             s1={:.3}, s4={:.3}",
            sk[0].kill_rate,
            sk[3].kill_rate
        );
        // Per-subclone tallies partition the tumor.
        let total_dead: usize = sk.iter().map(|s| s.total_dead).sum();
        let total_tumor: usize = sk.iter().map(|s| s.total_tumor).sum();
        assert_eq!(
            total_dead, r.total_dead,
            "subclone dead counts must sum to total_dead"
        );
        assert_eq!(
            total_tumor, r.total_tumor,
            "subclone tumor counts must sum to total_tumor"
        );
        // Without repopulation, the composition is static: final == initial.
        for s in &sk {
            assert_eq!(
                s.total_tumor, s.initial_tumor,
                "static clonal: subclone {} territory unchanged",
                s.subclone_id
            );
        }
    }

    /// #266 item 3: with spatial clonal **expansion** on, the resistant subclone
    /// gains territory relative to the vulnerable one — its share of living
    /// tumor cells rises from the initial assignment to the end of the run, as
    /// dead vulnerable sites are repopulated by surviving (mostly resistant)
    /// neighbors. Off-by-default leaves the composition static (asserted above).
    #[test]
    fn clonal_expansion_shifts_territory_toward_resistant() {
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 120,
        };
        let cond = Condition {
            name: "clonal_evolve".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: None, // uniform O2 ⇒ a strong baseline kill to drive turnover
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let r = run_one_condition_full(
            &cond,
            cfg,
            None,
            Overrides {
                clonal: Some(ClonalConfig::literature_4().with_repopulation(0.3)),
                ..Default::default()
            },
        );
        let sk = r.subclone_kills.expect("clonal run reports subclone_kills");
        let vuln = &sk[0]; // id 1: most vulnerable
        let resist = &sk[3]; // id 4: most resistant
                             // Territory shift: the resistant subclone's site count grows over the
                             // run while the vulnerable one's shrinks, as dead vulnerable sites are
                             // repopulated by surviving (mostly resistant) boundary neighbors.
        assert!(
            resist.total_tumor > resist.initial_tumor,
            "resistant subclone must expand: {} → {}",
            resist.initial_tumor,
            resist.total_tumor
        );
        assert!(
            vuln.total_tumor < vuln.initial_tumor,
            "vulnerable subclone must shrink: {} → {}",
            vuln.initial_tumor,
            vuln.total_tumor
        );
        // Deterministic.
        let r2 = run_one_condition_full(
            &cond,
            cfg,
            None,
            Overrides {
                clonal: Some(ClonalConfig::literature_4().with_repopulation(0.3)),
                ..Default::default()
            },
        );
        assert_eq!(
            r.total_dead, r2.total_dead,
            "expansion run is deterministic"
        );
    }

    /// #281 review: clonal expansion (repopulation on) composes with another
    /// realism layer (Treg/MDSC suppressor) — both run together, the run is
    /// deterministic, each layer reports its metric, and the resistant subclone
    /// still expands. Confirms repopulation (which only revives grid sites +
    /// rewrites subclone ids) doesn't break a layer it co-runs with.
    #[test]
    fn clonal_expansion_composes_with_suppressor() {
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 120,
        };
        let cond = Condition {
            name: "evolve_suppress".to_string(),
            treatment: Treatment::SDT,
            treatment_name: "SDT".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let immune = SpatialImmuneConfig {
            immune_kill_rate: 0.5,
            ..SpatialImmuneConfig::for_3d()
        };
        let run = || {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    immune: Some(immune),
                    clonal: Some(ClonalConfig::literature_4().with_repopulation(0.3)),
                    suppressor: Some(SuppressorConfig::enabled()),
                    ..Default::default()
                },
            )
        };
        let a = run();
        let b = run();
        assert_eq!(
            a.total_dead, b.total_dead,
            "composed evolution run is deterministic"
        );
        // Both layers report coherent metrics under composition.
        let sk = a.subclone_kills.expect("clonal reports subclone_kills");
        assert!(
            sk[3].total_tumor > sk[3].initial_tumor,
            "resistant subclone still expands under suppressor composition: {} → {}",
            sk[3].initial_tumor,
            sk[3].total_tumor
        );
        assert!(
            a.suppressor_source_count.unwrap_or(0) > 0,
            "suppressor still reports its niche census alongside expansion"
        );
    }

    /// Locks the MUFA/`lipid_unsat` axis in CI (#265 review): two K=1 configs
    /// differing ONLY in `lipid_unsat_mul` must produce different kill counts.
    /// The MUFA-enriched clone (lower oxidizable PUFA) dies less. Guards against
    /// the axis silently going inert again (e.g. if it were moved back onto the
    /// homeostatically-reset `state.mufa_protection`).
    #[test]
    fn clonal_lipid_unsat_axis_reduces_kills() {
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 80,
        };
        let cond = Condition {
            name: "lipid_axis".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: None,
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        // K=1, perturbing ONLY lipid_unsat (iron/gpx4 held at identity).
        let only_lipid = |lipid_unsat_mul: f64| ClonalConfig {
            perturbations: vec![SubclonePerturbation {
                iron_mul: 1.0,
                gpx4_mul: 1.0,
                lipid_unsat_mul,
            }],
            repopulation_rate: 0.0,
        };
        let baseline = run_one_condition_full(
            &cond,
            cfg,
            None,
            Overrides {
                clonal: Some(only_lipid(1.0)),
                ..Default::default()
            },
        );
        let mufa_enriched = run_one_condition_full(
            &cond,
            cfg,
            None,
            Overrides {
                clonal: Some(only_lipid(0.5)),
                ..Default::default()
            },
        );
        assert!(
            baseline.total_dead > 0,
            "baseline must kill some cells; got {}",
            baseline.total_dead
        );
        assert!(
            mufa_enriched.total_dead < baseline.total_dead,
            "MUFA-enriched (lower lipid_unsat) must reduce kills: enriched={}, baseline={}",
            mufa_enriched.total_dead,
            baseline.total_dead
        );
    }

    // ===== Explicit vasculature (#191) =====

    /// The #191 comparison at the simulation level: explicit internal vessels
    /// produce a different, irregular (non-radial) O₂ field than the
    /// edge-distance proxy, which materially changes the RSL3 kill outcome.
    /// (Direction is config-dependent: the edge proxy oxygenates the entire
    /// surface shell uniformly, whereas a sparse internal vessel set covers it
    /// irregularly — here the well-vascularized preset kills fewer, not more.)
    /// vasculature must report `vascular_hypoxic_fraction`; the edge default
    /// omits it.
    #[test]
    fn vasculature_oxygenates_core_and_changes_rsl3_kills() {
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 80,
        };
        let cond = Condition {
            name: "vasc_rsl3".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let edge = run_one_condition_with_config(&cond, cfg, None);
        let vasc = run_one_condition_full(
            &cond,
            cfg,
            None,
            Overrides {
                vasculature: Some(VasculatureConfig::well_vascularized()),
                ..Default::default()
            },
        );
        assert!(
            edge.vascular_hypoxic_fraction.is_none(),
            "edge-distance default must omit vascular_hypoxic_fraction"
        );
        assert!(
            vasc.vascular_hypoxic_fraction.is_some(),
            "vasculature run must report vascular_hypoxic_fraction"
        );
        // The vessel field changes the oxygenation pattern, so the kill
        // outcome must shift materially (>20%) vs the edge-distance proxy.
        assert!(edge.total_dead > 0, "edge baseline must kill some cells");
        let rel = (edge.total_dead as f64 - vasc.total_dead as f64).abs() / edge.total_dead as f64;
        assert!(
            rel > 0.2,
            "explicit vessels must materially change RSL3 kills vs the edge proxy \
             (>20%): vasc={}, edge={}, rel={:.2}",
            vasc.total_dead,
            edge.total_dead,
            rel
        );
    }

    /// #343 PR 2: the reaction-diffusion supply field differs from the
    /// monotonic nearest-vessel proxy on the SAME vessel network at the SAME λ
    /// (the genuine multi-vessel/consumption effect), is deterministic, and
    /// shifts the RSL3 kill outcome — while staying coherent (both fields in
    /// [0,1]). Also prints the proxy-vs-RD field divergence stats used by
    /// analysis/reaction-diffusion-benchmark.md (regenerate with
    /// `cargo test -p sim-tme-3d reaction_diffusion_supply_differs -- --nocapture`).
    #[test]
    fn reaction_diffusion_supply_differs_from_the_proxy_but_is_coherent() {
        // Field-level proxy-vs-RD comparison over the tumor cells of a generated
        // sphere, on the SAME vessels at the SAME λ. Returns
        // (n_vessels, mean_abs, max_abs, enriched%, depleted%, hyp_proxy, hyp_rd)
        // and prints the line consumed by analysis/reaction-diffusion-benchmark.md.
        let dim = 36;
        let lambda = ZONE_REF_LAMBDA;
        let grid = TumorGrid3D::generate(dim, dim, dim, CELL_SIZE_UM, SEED);

        // Single isolated 3-D point source on an all-tumor cube: isolates the
        // SOURCE-GEOMETRY term (planar proxy exp(-d/λ) vs a 3-D point source's
        // Yukawa exp(-r/λ)/r) with zero superposition and zero extra
        // consumption. The near-source gap here is most of the whole-tumor
        // overestimate, i.e. the dominant driver is geometry, not multi-vessel
        // effects. (Used by analysis/reaction-diffusion-benchmark.md.)
        {
            let mut cube = TumorGrid3D::generate(21, 21, 21, CELL_SIZE_UM, SEED);
            for c in cube.cells.iter_mut() {
                c.is_tumor = true;
            }
            let center = vec![(10.0, 10.0, 10.0)];
            let p = vessel_supply_field(&cube, &center, lambda);
            let r = reaction_diffusion_supply_field(
                &cube,
                &center,
                &ReactionDiffusionConfig::new(lambda),
            );
            let nb = cube.flat_index(11, 10, 10); // first neighbor, 1 cell from source
            eprintln!(
                "[RD-benchmark] single point source (all-tumor cube): first-neighbour proxy={:.3} RD={:.3} (Δ={:.3})",
                p[nb],
                r[nb],
                r[nb] - p[nb],
            );
            // Geometry alone (one vessel) already drops the near-source field
            // well below the proxy — the bulk of the whole-tumor overestimate.
            assert!(
                p[nb] - r[nb] > 0.3,
                "single-source geometry gap should be large: proxy={:.3} RD={:.3}",
                p[nb],
                r[nb]
            );
        }

        let compare = |vcfg: &VasculatureConfig, label: &str| -> (f64, f64, f64) {
            let vessels = place_vessels_3d(&grid, vcfg, VESSEL_SEED);
            let proxy = vessel_supply_field(&grid, &vessels, lambda);
            let rd = reaction_diffusion_supply_field(
                &grid,
                &vessels,
                &ReactionDiffusionConfig::new(lambda),
            );
            // Determinism (no RNG, fixed SOR sweep order).
            let rd2 = reaction_diffusion_supply_field(
                &grid,
                &vessels,
                &ReactionDiffusionConfig::new(lambda),
            );
            assert_eq!(rd, rd2, "RD supply must be deterministic");
            let eps = 1e-3;
            let (mut sum_abs, mut max_abs, mut n, mut enriched, mut depleted) =
                (0.0, 0.0_f64, 0usize, 0usize, 0usize);
            let (mut hyp_proxy, mut hyp_rd) = (0usize, 0usize);
            for (idx, cell) in grid.cells.iter().enumerate() {
                if !cell.is_tumor {
                    continue;
                }
                n += 1;
                let d = rd[idx] - proxy[idx];
                sum_abs += d.abs();
                max_abs = max_abs.max(d.abs());
                if d > eps {
                    enriched += 1;
                } else if d < -eps {
                    depleted += 1;
                }
                assert!((0.0..=1.0).contains(&proxy[idx]) && (0.0..=1.0).contains(&rd[idx]));
                if proxy[idx] < 0.1 {
                    hyp_proxy += 1;
                }
                if rd[idx] < 0.1 {
                    hyp_rd += 1;
                }
            }
            let mean_abs = sum_abs / n as f64;
            let (hp, hr) = (hyp_proxy as f64 / n as f64, hyp_rd as f64 / n as f64);
            eprintln!(
                "[RD-benchmark] {label}: dim={dim} vessels={} λ={lambda} tumor_cells={n}\n  \
                 mean|RD−proxy|={mean_abs:.4} max|RD−proxy|={max_abs:.4}\n  \
                 enriched(RD>proxy)={:.1}% depleted(RD<proxy)={:.1}%\n  \
                 hypoxic_fraction(<0.1): proxy={hp:.3} RD={hr:.3}",
                vessels.len(),
                100.0 * enriched as f64 / n as f64,
                100.0 * depleted as f64 / n as f64,
            );
            (mean_abs, hp, hr)
        };
        let (mean_sparse, hp_sparse, hr_sparse) =
            compare(&VasculatureConfig::poorly_vascularized(), "sparse");
        let (mean_dense, _, _) = compare(&VasculatureConfig::well_vascularized(), "dense");
        // The two fields genuinely differ (not the same model under a rename).
        assert!(
            mean_sparse > 0.01 && mean_dense > 0.01,
            "RD and proxy should differ materially: sparse={mean_sparse:.4}, dense={mean_dense:.4}"
        );
        // In a sparse 3D network the proxy is OPTIMISTIC: it ignores cumulative
        // consumption + 3D geometric spreading, so it under-reports hypoxia. RD
        // predicts a strictly higher hypoxic fraction (the "where the proxy is
        // misleading" result).
        assert!(
            hr_sparse > hp_sparse,
            "RD should report more hypoxia than the proxy in a sparse 3D network: \
             proxy={hp_sparse:.3}, RD={hr_sparse:.3}"
        );

        // End-to-end: through the full sim the supply difference shifts kills.
        // Sparse vessels, so the RD field's extra hypoxia (the proxy is
        // optimistic) translates into materially fewer hypoxia-sensitive RSL3
        // kills than the proxy predicts.
        let cfg = RunConfig {
            grid_dim: 30,
            n_steps: 120,
        };
        let cond = Condition {
            name: "rd_rsl3".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let run = |rd_on| {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    vasculature: Some(VasculatureConfig::poorly_vascularized()),
                    reaction_diffusion: rd_on,
                    ..Default::default()
                },
            )
        };
        let proxy_run = run(false);
        let rd_run = run(true);
        assert!(
            proxy_run.vascular_hypoxic_fraction.is_some()
                && rd_run.vascular_hypoxic_fraction.is_some(),
            "both supply models report a vascular hypoxic fraction"
        );
        let (hyp_proxy, hyp_rd) = (
            proxy_run.vascular_hypoxic_fraction.unwrap(),
            rd_run.vascular_hypoxic_fraction.unwrap(),
        );
        eprintln!(
            "[RD-benchmark] end-to-end RSL3 (dim=30, sparse): proxy total_dead={} hyp={hyp_proxy:.3} | RD total_dead={} hyp={hyp_rd:.3}",
            proxy_run.total_dead, rd_run.total_dead,
        );
        // The supply model carries all the way through to the emitted vascular
        // hypoxic fraction: RD reports materially more hypoxia than the
        // optimistic proxy. (The downstream RSL3 kill COUNT is a threshold-damped
        // readout — both are ~0 on a sparse, mostly-hypoxic tumor — so the
        // hypoxic fraction, not the kill count, is the sensitive end-to-end
        // signal.)
        assert!(
            hyp_rd > hyp_proxy + 0.05,
            "RD must report materially more hypoxia end-to-end than the proxy: \
             proxy={hyp_proxy:.3}, rd={hyp_rd:.3}"
        );
    }

    /// The reaction-diffusion flag is gated on explicit vasculature: with no
    /// vessels it is inert, so a run with `reaction_diffusion: true` and no
    /// vasculature is identical to the default. (This is why enabling the flag
    /// keeps the production matrix — which has no vasculature — byte-identical.)
    #[test]
    fn reaction_diffusion_flag_is_inert_without_vasculature() {
        let cfg = RunConfig {
            grid_dim: 16,
            n_steps: 40,
        };
        let cond = Condition {
            name: "rd_gate".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let base = run_one_condition_with_config(&cond, cfg, None);
        let rd_no_vessels = run_one_condition_full(
            &cond,
            cfg,
            None,
            Overrides {
                reaction_diffusion: true,
                ..Default::default()
            },
        );
        assert_eq!(
            base.total_dead, rd_no_vessels.total_dead,
            "reaction_diffusion must be inert without vasculature (gating ⇒ byte-identical)"
        );
        assert!(
            rd_no_vessels.vascular_hypoxic_fraction.is_none(),
            "no vasculature ⇒ no vascular hypoxic fraction"
        );
    }

    /// The `reaction-diffusion` snapshot preset is name-keyed in `run_snapshot`
    /// (`preset.name == "reaction-diffusion"`) and the RD field only has vessel
    /// Dirichlet sources because the preset also sets `vasculature: true`. Lock
    /// both invariants so a rename or a `vasculature` flip can't silently fall
    /// back to the exp proxy and emit the wrong field (same guard the other
    /// name-keyed presets carry).
    #[test]
    fn reaction_diffusion_snapshot_preset_is_wired() {
        let p = resolve_snapshot("reaction-diffusion");
        assert_eq!(
            p.name, "reaction-diffusion",
            "the name-key in run_snapshot matches on this exact string"
        );
        assert!(
            p.vasculature,
            "the reaction-diffusion preset must enable vasculature so the RD solver has vessel sources"
        );
        assert!(
            matches!(p.treatment, Treatment::RSL3),
            "the preset demonstrates RSL3 (hypoxia-sensitive) on the RD supply field"
        );
    }

    #[test]
    fn fractal_vasculature_is_hypoxier_than_random() {
        // #268: a fractal-branching vessel tree clusters along branches with
        // avascular gaps, so at near-equal POINT COUNT it leaves a HIGHER
        // vascular hypoxic fraction than uniform-random placement. Validates the
        // topology dispatch end-to-end through sim-tme-3d.
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 80,
        };
        // COUNT-PARITY GUARD: reconstruct the exact grid the run builds and check
        // both placers yield ~the same number of vessel points. Without this, a
        // higher hypoxic fraction could be a "fractal is just sparser" artifact
        // rather than a clustering (topology) effect.
        {
            let g =
                TumorGrid3D::generate(cfg.grid_dim, cfg.grid_dim, cfg.grid_dim, CELL_SIZE_UM, SEED);
            let vcfg = VasculatureConfig::well_vascularized();
            let n_random = place_vessels_3d(&g, &vcfg, VESSEL_SEED).len();
            let n_fractal = place_vessels_fractal_3d(&g, &vcfg.with_fractal(), VESSEL_SEED).len();
            assert!(
                n_fractal as f64 >= 0.9 * n_random as f64 && n_fractal <= n_random + 1,
                "fractal count {n_fractal} must be within 10% of random count {n_random} \
                 so the hypoxic comparison is count-controlled, not a sparsity artifact"
            );
        }
        let cond = Condition {
            name: "fractal_rsl3".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let run = |topo| {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    vasculature: Some(VasculatureConfig {
                        topology: topo,
                        ..VasculatureConfig::well_vascularized()
                    }),
                    ..Default::default()
                },
            )
        };
        let random = run(VesselTopology::Random);
        let fractal = run(VesselTopology::Fractal);
        let fractal_again = run(VesselTopology::Fractal);
        assert_eq!(
            fractal.total_dead, fractal_again.total_dead,
            "fractal vasculature must be deterministic"
        );
        let hyp_random = random
            .vascular_hypoxic_fraction
            .expect("random reports hypoxic fraction");
        let hyp_fractal = fractal
            .vascular_hypoxic_fraction
            .expect("fractal reports hypoxic fraction");
        assert!(
            hyp_fractal > hyp_random,
            "fractal network should be hypoxier (more avascular gaps) than random: \
             fractal={hyp_fractal:.3}, random={hyp_random:.3}"
        );
    }

    // ===== 3D spheroid radial biochemistry (#197) =====

    /// AC comparison: a radial spheroid (glycolytic rim / OXPHOS mid /
    /// persister core + GSH/iron/MUFA gradients, run under Params::spheroid())
    /// produces a materially different kill outcome than the default
    /// random-phenotype grid — answering "does radial structure change the
    /// kill rate, or just redistribute where cells die?" with: it changes it.
    #[test]
    fn contact_resistance_reduces_rsl3_kills() {
        // #270: cell-cell contact resistance (E-cadherin/Merlin/NF2-YAP, Wu
        // 2019) scales down the durable PUFA + iron axes for densely-contacted
        // cells, so a pharmacologic ferroptosis inducer (RSL3, which depends on
        // endogenous PUFA/iron) kills far fewer cells with the layer on. The
        // baseline run (no overrides) is the byte-identical matrix path.
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 80,
        };
        let cond = Condition {
            name: "contact_rsl3".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let baseline = run_one_condition_with_config(&cond, cfg, None);
        let run_contact = || {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    contact: Some(ContactConfig::literature()),
                    ..Default::default()
                },
            )
        };
        let contact = run_contact();
        assert!(
            baseline.total_dead > 0,
            "RSL3 baseline must kill some cells"
        );
        assert!(
            contact.total_dead < baseline.total_dead,
            "contact resistance should reduce RSL3 kills: baseline={}, contact={}",
            baseline.total_dead,
            contact.total_dead
        );
        // The effect is large (RSL3 depends entirely on endogenous PUFA/iron):
        // well under half the baseline kills survive the contact brake.
        assert!(
            (contact.total_dead as f64) < 0.5 * baseline.total_dead as f64,
            "contact should at least halve RSL3 kills: baseline={}, contact={}",
            baseline.total_dead,
            contact.total_dead
        );
        // Deterministic (geometric, no RNG).
        assert_eq!(contact.total_dead, run_contact().total_dead);
    }

    /// #363: phenotype-specific MUFA accumulation rates are wired into sim-tme-3d
    /// OFF-BY-DEFAULT. In the spheroid context (`Params::spheroid()`, MUFA active):
    /// an IDENTITY config is filtered out (is_identity) and reproduces the
    /// spheroid baseline byte-for-byte (the production-meaningful guarantee), and
    /// `None` likewise. A non-identity config runs deterministically. The RATE's
    /// effect on the MUFA timecourse is proven at the library level
    /// (`ferroptosis_core::biochem` `sim_cell_step_reads_per_cell_mufa_rate`,
    /// which directly measures `mufa_protection` divergence — the right readout).
    /// We deliberately do NOT assert a kill-COUNT delta here: the spheroid's
    /// bistable switch + cysteine-limited core cap make the kill count insensitive
    /// to the MUFA rate in this regime (the cells that die are not MUFA-rescuable
    /// and the survivors are not MUFA-dependent), so a count assertion would be
    /// fragile and misleading. Byte-identity + determinism is what matters for the
    /// production matrix.
    #[test]
    fn phenotype_mufa_off_by_default_byte_identical_in_spheroid() {
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 80,
        };
        let cond = Condition {
            name: "phenotype_mufa_rsl3".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let spheroid_only = || Overrides {
            spheroid: Some(SpheroidConfig::literature()),
            ..Default::default()
        };
        let run = |pm: Option<PhenotypeMufaConfig>| {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    phenotype_mufa: pm,
                    ..spheroid_only()
                },
            )
            .total_dead
        };
        let baseline = run_one_condition_full(&cond, cfg, None, spheroid_only()).total_dead;
        assert!(baseline > 0, "spheroid RSL3 baseline must kill some cells");
        // Off-by-default: None AND an identity config (filtered out) both
        // reproduce the spheroid baseline byte-for-byte.
        assert_eq!(
            run(None),
            baseline,
            "None phenotype_mufa must equal baseline"
        );
        assert_eq!(
            run(Some(PhenotypeMufaConfig::identity())),
            baseline,
            "identity phenotype_mufa must be byte-identical to the spheroid baseline (gate filters it)"
        );
        // A non-identity config runs and is deterministic (geometric per-cell
        // setup, no RNG), so the layer is reachable without breaking determinism.
        let nonident = PhenotypeMufaConfig {
            oxphos: 1.5,
            persister: 2.0,
            persister_nrf2: 2.0,
            persister_cap: 1.5,
            ..PhenotypeMufaConfig::identity()
        };
        assert_eq!(
            run(Some(nonident)),
            run(Some(nonident)),
            "a non-identity phenotype_mufa run must be deterministic"
        );
    }

    /// #270 item 3b: the radial nutrient gradient lowers the durable antioxidant
    /// setpoint (cell.nrf2) toward the nutrient-starved core, so the core has
    /// less GSH/GPX4 regeneration and an antioxidant-sensitive inducer (RSL3)
    /// kills MORE cells with the layer on. Deterministic (geometric, no RNG).
    #[test]
    fn nutrient_gradient_increases_rsl3_kills() {
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 80,
        };
        let cond = Condition {
            name: "nutrient_rsl3".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let baseline = run_one_condition_with_config(&cond, cfg, None);
        let run_nutrient = || {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    nutrient: Some(NutrientConfig::literature()),
                    ..Default::default()
                },
            )
        };
        let nutrient = run_nutrient();
        assert!(
            baseline.total_dead > 0,
            "RSL3 baseline must kill some cells"
        );
        assert!(
            nutrient.total_dead > baseline.total_dead,
            "nutrient deprivation (lower core antioxidant capacity) must raise RSL3 \
             kills: baseline={}, nutrient={}",
            baseline.total_dead,
            nutrient.total_dead
        );
        assert_eq!(
            nutrient.total_dead,
            run_nutrient().total_dead,
            "nutrient layer is deterministic"
        );
    }

    /// #270: lock the `--snapshot=nutrient` preset -> Overrides wiring without
    /// the full render. The preset enables nutrient and leaves the geometry
    /// layers off (it runs on the centred sphere).
    #[test]
    fn nutrient_snapshot_preset_is_wired() {
        let p = resolve_snapshot("nutrient");
        assert_eq!(p.name, "nutrient");
        assert!(
            p.nutrient,
            "the nutrient preset must enable the nutrient layer"
        );
        assert!(
            !p.slab && !p.spheroid && !p.contact,
            "nutrient runs on the plain centred sphere"
        );
    }

    /// #363: lock the `--snapshot=phenotype-mufa` preset -> Overrides wiring. The
    /// preset enables the phenotype-MUFA layer AND the spheroid context (the only
    /// MUFA-active path, so the rate is non-inert), and resolves to a non-identity
    /// `PhenotypeMufaConfig` that the run will actually apply.
    #[test]
    fn phenotype_mufa_snapshot_preset_is_wired() {
        let p = resolve_snapshot("phenotype-mufa");
        assert_eq!(p.name, "phenotype-mufa");
        assert!(
            p.phenotype_mufa,
            "the phenotype-mufa preset must enable the phenotype-MUFA layer"
        );
        assert!(
            p.spheroid,
            "phenotype-mufa must run in the spheroid context (the only MUFA-active path)"
        );
        // The preset must map to a non-identity config (else the gate filters it
        // out and the layer is silently inert).
        assert!(
            p.phenotype_mufa && !PhenotypeMufaConfig::literature().is_identity(),
            "the preset's literature() config must be non-identity so it is applied"
        );
    }

    /// #380: lock the `--snapshot=sdt-o2dep` preset -> Overrides wiring. The
    /// preset runs SDT with the exo-ROS made fully O2-dependent (`sdt_o2_dependence
    /// = 1.0`), so the hypoxic core survives. Every OTHER preset must keep
    /// `sdt_o2_dependence = 0.0` (the byte-identical O2-independent default).
    #[test]
    fn sdt_o2dep_snapshot_preset_is_wired() {
        let p = resolve_snapshot("sdt-o2dep");
        assert_eq!(p.name, "sdt-o2dep");
        assert_eq!(p.treatment, Treatment::SDT, "sdt-o2dep must run SDT");
        assert_eq!(
            p.sdt_o2_dependence, 1.0,
            "sdt-o2dep must make the exo-ROS fully O2-dependent"
        );
        // Off-by-default invariant: only the two intentional SDT-O2 presets
        // (sdt-o2dep = Type II / dep 1.0, and sdt-typei = Type-I-heavy / dep 0.3,
        // #468) set a non-zero dependence; every other preset keeps 0.0 so the
        // matrix stays byte-identical.
        for q in SNAPSHOTS {
            if q.name != "sdt-o2dep" && q.name != "sdt-typei" {
                assert_eq!(
                    q.sdt_o2_dependence, 0.0,
                    "preset {} must keep sdt_o2_dependence=0.0 (byte-identical)",
                    q.name
                );
            }
        }
    }

    /// #381: lock the `--snapshot=ferritinophagy` preset -> Overrides wiring. The
    /// preset runs RSL3 with the NCOA4-ferritinophagy + hypoxia-iron coupling on.
    /// Every OTHER preset must keep `ferritinophagy = false` (byte-identical).
    #[test]
    fn ferritinophagy_snapshot_preset_is_wired() {
        let p = resolve_snapshot("ferritinophagy");
        assert_eq!(p.name, "ferritinophagy");
        assert_eq!(p.treatment, Treatment::RSL3, "ferritinophagy must run RSL3");
        assert!(
            p.ferritinophagy,
            "the ferritinophagy preset must enable the dynamic-iron coupling"
        );
        for q in SNAPSHOTS {
            if q.name != "ferritinophagy" {
                assert!(
                    !q.ferritinophagy,
                    "preset {} must keep ferritinophagy=false (byte-identical)",
                    q.name
                );
            }
        }
    }

    #[test]
    fn radial_spheroid_changes_kills_vs_random() {
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 80,
        };
        let cond = Condition {
            name: "spheroid_cmp".to_string(),
            treatment: Treatment::SDT,
            treatment_name: "SDT".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let random = run_one_condition_with_config(&cond, cfg, None);
        let radial = run_one_condition_full(
            &cond,
            cfg,
            None,
            Overrides {
                spheroid: Some(SpheroidConfig::literature()),
                ..Default::default()
            },
        );
        assert!(
            random.total_dead > 0,
            "random baseline must kill some cells"
        );
        // The radial spheroid is net more SDT-resistant than the default
        // (random) phenotype split: its larger, volume-correct hypoxic/persister
        // core (#270 — core is now ~0.39 of volume, was ~0.04) kills fewer cells.
        // Deterministic (seeded), so the ~9.7% gap is stable; assert a material
        // change with margin. (Threshold was 10%, tuned to the pre-#270 inverted
        // config; the volume-correct config gives radial≈4601 vs random≈5093.)
        // NOTE: the < direction couples to the current net balance of the core's
        // MIXED mechanism (persister phenotype = resistant; low GSH + high iron =
        // vulnerable; hypoxic = resistant to O2-dependent ferroptosis). A future
        // gradient-strength calibration could flip the net, in which case revisit
        // this direction assertion (the magnitude/material-change check stands).
        assert!(
            radial.total_dead < random.total_dead,
            "volume-correct spheroid (larger hypoxic core) should kill fewer than random: \
             radial={}, random={}",
            radial.total_dead,
            random.total_dead
        );
        let rel =
            (random.total_dead as f64 - radial.total_dead as f64).abs() / random.total_dead as f64;
        assert!(
            rel > 0.05,
            "radial spheroid structure must materially change kills (>5%) vs random: \
             radial={}, random={}, rel={:.3}",
            radial.total_dead,
            random.total_dead,
            rel
        );
    }

    /// Spheroid mode is byte-identical-safe when off: the default (no spheroid
    /// override) path is unchanged. Also confirms `Params::spheroid` differs
    /// from `Params::default` (so the mode actually does something).
    #[test]
    fn spheroid_params_differ_from_default() {
        use ferroptosis_core::params::Params;
        let d = Params::default();
        let s = Params::spheroid();
        assert_eq!(d.scd_mufa_max, 0.0, "2D default has no MUFA cap");
        assert!(s.scd_mufa_max > 0.0, "spheroid has partial MUFA");
        assert!(
            s.initial_mufa_protection > 0.0 && s.initial_mufa_protection < 0.4,
            "spheroid M_ss is between 2D (0) and in-vivo (0.40); got {}",
            s.initial_mufa_protection
        );
    }

    // ===== Patient-scale slab (#240) =====

    /// Headline #240 result: a patient-scale slab (≥5 mm deep) has dramatically
    /// lower drug efficacy than the in-vitro spheroid, because drug/O2
    /// penetration (Krogh ~150 µm) leaves a deep slab essentially deprived.
    /// Compares the default 60³ spheroid against the deep slab under the same
    /// treatment, and checks the scale-interpretation reporting.
    #[test]
    fn patient_scale_slab_kills_far_less_than_spheroid() {
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 80,
        };
        let cond = Condition {
            name: "scale_cmp".to_string(),
            treatment: Treatment::SDT,
            treatment_name: "SDT".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        // In-vitro spheroid (default 60³-style grid, radial O2).
        let spheroid = run_one_condition_with_config(&cond, cfg, None);
        // Patient-scale slab: +z face at 4 mm ⇒ the whole slab is ~drug/O2-deprived.
        let slab = run_one_condition_full(
            &cond,
            cfg,
            None,
            Overrides {
                slab: Some(SlabConfig::patient_deep()),
                ..Default::default()
            },
        );
        assert!(
            spheroid.total_dead > 0,
            "spheroid baseline must kill some cells; got {}",
            spheroid.total_dead
        );
        assert!(
            spheroid.scale_interpretation.is_none(),
            "non-slab run omits scale_interpretation"
        );
        let interp = slab
            .scale_interpretation
            .as_ref()
            .expect("slab run reports scale_interpretation");
        assert!(
            interp.contains("mm"),
            "interpretation names the depth scale: {interp}"
        );
        // Dramatically lower efficacy at patient scale: a 4 mm-deep slab kills
        // a small fraction of what the in-vitro spheroid does.
        assert!(
            (slab.total_dead as f64) < 0.2 * spheroid.total_dead as f64,
            "patient-scale slab must kill far less than the spheroid: \
             slab={}, spheroid={}",
            slab.total_dead,
            spheroid.total_dead
        );
    }

    /// Within slab mode, a surface slab (vessel at +z) kills more than a deep
    /// slab — the depth-dependence the patient-scale model is built to show.
    #[test]
    fn shallow_slab_outkills_deep_slab() {
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 80,
        };
        let cond = Condition {
            name: "slab_depth".to_string(),
            treatment: Treatment::SDT,
            treatment_name: "SDT".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let run = |slab: SlabConfig| {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    slab: Some(slab),
                    ..Default::default()
                },
            )
            .total_dead
        };
        let shallow = run(SlabConfig::surface());
        let deep = run(SlabConfig::patient_deep());
        assert!(
            shallow > deep,
            "surface slab should out-kill the deep slab: shallow={shallow}, deep={deep}"
        );
    }

    /// Slab drug-coupling under RSL3 (#240 review #5). The SDT depth tests
    /// exercise the exo-ROS drug branch; RSL3 routes through the *other* two
    /// supply-scaled drug paths — the constant-knockdown GPX4 correction
    /// `(1 - inhib·supply)/(1 - inhib)` and the `rsl3_drug_avail *= supply`
    /// scaling. A surface vs deep RSL3 comparison pins that those RSL3-specific
    /// branches honor the slab supply field independently of the vessel path.
    #[test]
    fn slab_supply_scales_rsl3_kills() {
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 80,
        };
        let cond = Condition {
            name: "slab_rsl3".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let run = |slab: SlabConfig| {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    slab: Some(slab),
                    ..Default::default()
                },
            )
            .total_dead
        };
        let shallow = run(SlabConfig::surface());
        let deep = run(SlabConfig::patient_deep());
        assert!(
            shallow > deep,
            "RSL3 surface slab should out-kill the deep slab (drug/O2 supply \
             scaling on the RSL3 path): shallow={shallow}, deep={deep}"
        );
    }

    #[test]
    fn slab_internal_vessels_increase_deep_killing() {
        // #272 coupling: a DEEP slab is drug/O2-starved, so the planar-only
        // model kills little (depth collapse, #240). Adding an internal vessel
        // network (combined planar-MAX-vessel supply) delivers drug/O2 to focal
        // deep pockets, so total kills RISE vs the planar-only slab. (Supply
        // scales DELIVERY here, so well-perfused ⇒ more drug ⇒ more death — the
        // vessels make deep tissue LESS therapy-resistant, not "rescued".)
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 80,
        };
        let cond = Condition {
            name: "slab_vessels_rsl3".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let run = |vasc: Option<VasculatureConfig>| {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    slab: Some(SlabConfig::patient_deep()),
                    vasculature: vasc,
                    ..Default::default()
                },
            )
            .total_dead
        };
        let planar_only = run(None);
        let with_vessels = run(Some(VasculatureConfig::well_vascularized()));
        let with_vessels_again = run(Some(VasculatureConfig::well_vascularized()));
        assert_eq!(
            with_vessels, with_vessels_again,
            "slab+vasculature must be deterministic"
        );
        assert!(
            with_vessels > planar_only,
            "internal vessels should deliver drug to deep pockets and raise kills \
             on a deep slab: planar_only={planar_only}, with_vessels={with_vessels}"
        );
    }

    #[test]
    fn combine_supply_max_takes_elementwise_max() {
        // #272: directly pin the field invariant the slab+vessel run/overlay use
        // — the combined supply is the element-wise max of the planar and vessel
        // fields, never below either source, and bounded in [0,1] when both
        // inputs are. (Fast; no full sim run needed.)
        let planar = [0.0, 0.2, 0.9, 0.5, 1.0];
        let vessel = [0.1, 0.8, 0.3, 0.5, 0.0];
        let combined = combine_supply_max(&planar, &vessel);
        assert_eq!(combined, vec![0.1, 0.8, 0.9, 0.5, 1.0]);
        for i in 0..planar.len() {
            assert!(
                combined[i] >= planar[i] && combined[i] >= vessel[i],
                "combined[{i}]={} below a source (planar={}, vessel={})",
                combined[i],
                planar[i],
                vessel[i]
            );
            assert!(
                (0.0..=1.0).contains(&combined[i]),
                "combined[{i}]={} out of [0,1]",
                combined[i]
            );
        }
    }

    /// #336 oxygen-dependent SDT: under an edge-distance O2 gradient (a hypoxic
    /// spheroid core), an O2-dependent (Type II) SDT (`sdt_o2_dependence = 1.0`)
    /// kills fewer cells than the default O2-independent SDT
    /// (`sdt_o2_dependence = 0.0`), because the exogenous ROS yield collapses in
    /// the hypoxic core like a clinical Type II sonosensitizer. The default
    /// (0.0) reproduces the historical behavior (the matrix byte-identity is
    /// guarded separately by the golden tests); this asserts the new knob has
    /// the predicted DIRECTION and is deterministic.
    #[test]
    fn o2_dependent_sdt_reduces_hypoxic_kill() {
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 80,
        };
        let cond = Condition {
            name: "sdt_o2dep".to_string(),
            treatment: Treatment::SDT,
            treatment_name: "SDT".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA), // edge-distance hypoxia gradient on
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let run = |dep: f64| {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    sdt_o2_dependence: dep,
                    ..Default::default()
                },
            )
            .total_dead
        };
        let o2_independent = run(0.0); // default / historical (optimistic SDT)
        let o2_dependent = run(1.0); // fully Type II / O2-dependent
        let o2_dependent_again = run(1.0);
        assert_eq!(
            o2_dependent, o2_dependent_again,
            "O2-dependent SDT must be deterministic"
        );
        assert!(
            o2_dependent < o2_independent,
            "Type II (O2-dependent) SDT should kill fewer under hypoxia than the \
             O2-independent default: dep1={o2_dependent}, dep0={o2_independent}"
        );
    }

    /// #272 depth-graded slab phenotype wiring: (1) it is gated on slab mode, so
    /// a `slab_phenotype` override with NO slab grid is inert (it can never
    /// perturb the spheroid matrix); (2) turning it on actually changes the run
    /// (the depth-layered rim→core phenotype mix is a different tumor than the
    /// flat bulk mix); (3) it stays deterministic. The per-layer phenotype
    /// assignment itself is pinned in the `slab` module's unit tests.
    #[test]
    fn depth_graded_slab_phenotype_is_gated_and_deterministic() {
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 80,
        };
        let cond = Condition {
            name: "slab_pheno".to_string(),
            treatment: Treatment::SDT,
            treatment_name: "SDT".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        // (1) slab_phenotype WITHOUT a slab grid is inert: identical to a bare
        // run, so it can never change the (no-slab) production matrix.
        let bare = run_one_condition_full(&cond, cfg, None, Overrides::default());
        let pheno_no_slab = run_one_condition_full(
            &cond,
            cfg,
            None,
            Overrides {
                slab_phenotype: Some(SlabPhenotypeConfig::literature()),
                ..Default::default()
            },
        );
        assert_eq!(
            bare.total_dead, pheno_no_slab.total_dead,
            "slab_phenotype with no slab grid must be inert"
        );
        assert_eq!(bare.overall_kill_rate, pheno_no_slab.overall_kill_rate);

        // (2) + (3): on a surface slab, depth grading changes the outcome vs the
        // flat bulk mix, and both forms are deterministic.
        let run = |pheno: Option<SlabPhenotypeConfig>| {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    slab: Some(SlabConfig::surface()),
                    slab_phenotype: pheno,
                    ..Default::default()
                },
            )
            .total_dead
        };
        let flat = run(None);
        let graded = run(Some(SlabPhenotypeConfig::literature()));
        let graded_again = run(Some(SlabPhenotypeConfig::literature()));
        assert_eq!(
            graded, graded_again,
            "depth-graded slab must be deterministic"
        );
        assert_ne!(
            flat, graded,
            "depth grading should change the slab outcome vs the flat bulk mix: \
             flat={flat}, graded={graded}"
        );
    }

    // ===== Cross-layer composition (#278) =====

    /// The realism layers each seed an INDEPENDENT RNG with a distinct constant
    /// so they never correlate or perturb `generate`'s stream. A future layer
    /// reusing a constant would silently couple two layers' stochastic
    /// structure — assert the constants are pairwise distinct so that fails
    /// loudly here instead.
    #[test]
    fn realism_layer_seeds_are_pairwise_distinct() {
        // NOTE: when a new realism layer adds its own seed constant, add it
        // here too — this list is the (manual) registry of layer RNG seeds.
        let seeds = [
            ("SEED", SEED),
            ("SUBCLONE_SEED", SUBCLONE_SEED),
            ("VESSEL_SEED", VESSEL_SEED),
            ("SPHEROID_SEED", SPHEROID_SEED),
            ("SUPPRESSOR_SEED", SUPPRESSOR_SEED),
            ("SLAB_PHENOTYPE_SEED", SLAB_PHENOTYPE_SEED),
        ];
        for i in 0..seeds.len() {
            for j in (i + 1)..seeds.len() {
                assert_ne!(
                    seeds[i].1, seeds[j].1,
                    "realism-layer RNG seeds must be pairwise distinct: {} and {} collide",
                    seeds[i].0, seeds[j].0
                );
            }
        }
    }

    /// #443 A/B: the IFN-gamma -> System Xc- coupling sensitizes to ferroptosis.
    /// With immune on, an RSL3 grid run with the coupling ON kills at least as many
    /// tumor cells as the coupling-OFF baseline. IFN-gamma (seeded from the local
    /// DAMP of the immune response) suppresses GSH, so it can only RAISE ferroptosis,
    /// never protect; this confirms the immune-activation -> ferroptosis-sensitization
    /// link is wired and flows the right direction. Magnitude uncalibrated.
    #[test]
    fn ifngamma_couples_immune_activation_to_ferroptosis_sensitization() {
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 120,
        };
        let cond = Condition {
            name: "ifngamma_ab".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let immune = SpatialImmuneConfig {
            immune_kill_rate: 0.5,
            ..SpatialImmuneConfig::for_3d()
        };
        let run = |ifn: Option<IFNGammaConfig>| {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    immune: Some(immune),
                    ifngamma: ifn,
                    ..Default::default()
                },
            )
        };
        let baseline = run(None);
        let coupled = run(Some(IFNGammaConfig::literature()));
        // IFN-gamma seeded from the immune DAMP field suppresses GSH, so it can only
        // RAISE ferroptosis, never protect. Both runs are deterministic (identical
        // seeds; only the IFN-gamma field differs), so the delta is reproducible: at
        // this config the coupling roughly doubles kill (~0.076 -> ~0.154, delta
        // ~0.078). The 0.05 floor is well below that measured delta yet far above 0,
        // so it stays green on the real effect while still failing if the field
        // accumulation is attenuated (e.g. a halved per_damp / diffusion) or not
        // wired at all — a trivially-satisfied `>=` would catch neither.
        assert!(
            coupled.overall_kill_rate > baseline.overall_kill_rate + 0.05,
            "IFN-gamma coupling must measurably raise kill (expected delta ~0.078): coupled={} baseline={}",
            coupled.overall_kill_rate,
            baseline.overall_kill_rate
        );
    }

    /// #443 follow-up: the ACSL4 arm ALONE (System Xc⁻ arm switched off via an
    /// infinite IC50, so `system_xc_retention == 1.0`) still raises ferroptosis.
    /// This isolates the lipid arm — the IFN-gamma field seeds (per_damp > 0),
    /// diffuses, and transiently boosts each exposed cell's PUFA / lipid_unsat,
    /// which is the only active effect here — confirming the ACSL4 wiring flows
    /// the right direction independent of the GSH arm. Magnitude uncalibrated.
    #[test]
    fn ifngamma_acsl4_arm_alone_raises_ferroptosis() {
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 120,
        };
        let cond = Condition {
            name: "ifngamma_acsl4_ab".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let immune = SpatialImmuneConfig {
            immune_kill_rate: 0.5,
            ..SpatialImmuneConfig::for_3d()
        };
        let run = |ifn: Option<IFNGammaConfig>| {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    immune: Some(immune),
                    ifngamma: ifn,
                    ..Default::default()
                },
            )
        };
        let baseline = run(None);
        // ACSL4 ONLY: infinite IC50 ⇒ GSH retention is exactly 1.0 (System Xc⁻ arm
        // off), so the lipid_unsat boost is the sole sensitizer. per_damp > 0 so the
        // field still seeds from the immune DAMP.
        let acsl4_only = IFNGammaConfig {
            system_xc_ic50: f64::INFINITY,
            ..IFNGammaConfig::literature()
        };
        assert!(!acsl4_only.is_disabled(), "ACSL4-only config is active");
        let coupled = run(Some(acsl4_only));
        assert!(
            coupled.overall_kill_rate > baseline.overall_kill_rate,
            "ACSL4 arm alone must raise kill: coupled={} baseline={}",
            coupled.overall_kill_rate,
            baseline.overall_kill_rate
        );
    }

    /// #443 follow-up: the `--snapshot=ifngamma` preset is wired (enables both
    /// coupling arms under an immune RSL3 run).
    #[test]
    fn ifngamma_snapshot_preset_is_wired() {
        let p = resolve_snapshot("ifngamma");
        assert_eq!(p.name, "ifngamma");
        assert!(p.ifngamma, "the ifngamma preset must enable the coupling");
        assert!(
            p.immune_on,
            "the IFN-gamma field seeds from the immune DAMP signal, so immune must be on"
        );
        assert!(matches!(p.treatment, Treatment::RSL3));
        // literature() turns both arms on (non-disabled).
        assert!(!IFNGammaConfig::literature().is_disabled());
    }

    /// #446 A/B: the ALOX isoform-specific peroxidation + MCFA sensitization
    /// raises ferroptosis. An ALOX15-high, MCFA-exposed tumor (the `literature()`
    /// config: positive propagation boost + positive PUFA boost) kills MORE under
    /// RSL3 than the balanced-ALOX baseline, confirming the lipid-machinery
    /// sensitization axis is wired and flows the right direction (distinct from
    /// the GPX4/GSH/FSP1 defenses). Magnitude uncalibrated.
    #[test]
    fn alox_high_mcfa_raises_ferroptosis() {
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 120,
        };
        let cond = Condition {
            name: "alox_ab".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let run = |alox: Option<AloxConfig>| {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    alox,
                    ..Default::default()
                },
            )
        };
        let baseline = run(None);
        let alox_high = run(Some(AloxConfig::literature()));
        // ALOX15-high + MCFA can only RAISE peroxidation ⇒ more kill, never less.
        assert!(
            alox_high.overall_kill_rate > baseline.overall_kill_rate,
            "ALOX-high + MCFA must raise kill: alox_high={} baseline={}",
            alox_high.overall_kill_rate,
            baseline.overall_kill_rate
        );
        // Identity config is a no-op (byte-identical to baseline).
        let identity = run(Some(AloxConfig::identity()));
        assert_eq!(
            identity.overall_kill_rate, baseline.overall_kill_rate,
            "identity AloxConfig must be a no-op"
        );
    }

    /// #446: the `--snapshot=alox` preset is wired (RSL3 with the ALOX15-high +
    /// MCFA sensitization config).
    #[test]
    fn alox_snapshot_preset_is_wired() {
        let p = resolve_snapshot("alox");
        assert_eq!(p.name, "alox");
        assert!(p.alox, "the alox preset must enable the layer");
        assert!(matches!(p.treatment, Treatment::RSL3));
        // literature() is non-identity (both boosts positive).
        let lit = AloxConfig::literature();
        assert!(!lit.is_identity());
        assert!(lit.lp_propagation_boost() > 0.0 && lit.mcfa_pufa_boost() > 0.0);
    }

    /// #444 A/B: ACSL4 status stratifies ferroptosis sensitivity. Under the SAME
    /// RSL3 dose, an ACSL4-HIGH tumor kills MORE and an ACSL4-NEGATIVE tumor kills
    /// LESS than the ACSL4-normal baseline — the patient-stratification headline
    /// (ACSL4-high respond, ACSL4-negative are refractory via collapsed PUFA, a
    /// mechanism distinct from GPX4/GSH/FSP1; Doll 2017 PMID 27842070). Magnitude
    /// uncalibrated; the ordering is the result.
    #[test]
    fn acsl4_status_stratifies_ferroptosis_sensitivity() {
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 120,
        };
        let cond = Condition {
            name: "acsl4_ab".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let run = |status: Option<f64>| {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    acsl4_status: status,
                    ..Default::default()
                },
            )
        };
        let normal = run(None); // ACSL4-normal = the byte-identical baseline
        let high = run(Some(ferroptosis_core::acsl4::ACSL4_HIGH));
        let negative = run(Some(ACSL4_NEGATIVE));
        // ACSL4-high out-kills normal; ACSL4-negative is refractory (kills less).
        assert!(
            high.overall_kill_rate > normal.overall_kill_rate,
            "ACSL4-high must raise kill: high={} normal={}",
            high.overall_kill_rate,
            normal.overall_kill_rate
        );
        assert!(
            negative.overall_kill_rate < normal.overall_kill_rate,
            "ACSL4-negative must lower kill (refractory): negative={} normal={}",
            negative.overall_kill_rate,
            normal.overall_kill_rate
        );
    }

    /// #444: the `--snapshot=acsl4-negative` preset is wired (RSL3 on an
    /// ACSL4-negative, ferroptosis-refractory tumor).
    #[test]
    fn acsl4_negative_snapshot_preset_is_wired() {
        let p = resolve_snapshot("acsl4-negative");
        assert_eq!(p.name, "acsl4-negative");
        assert!(
            p.acsl4_negative,
            "the preset must enable ACSL4-negative status"
        );
        assert!(matches!(p.treatment, Treatment::RSL3));
        // ACSL4-negative maps to the null-floor boost (-1, collapsed PUFA).
        assert_eq!(pufa_boost_from_status(ACSL4_NEGATIVE), -1.0);
    }

    /// #465 A/B: the ESCRT-III membrane-repair brake reduces RSL3 kill. A cell
    /// whose lipid peroxide crosses the death threshold can be resealed for a
    /// finite per-cell budget, delaying/blocking death execution, so an ESCRT-high
    /// tumor kills LESS than the no-repair baseline under the same RSL3 dose. The
    /// resistance mechanism is membrane repair, not a redox/lipid defense (Dai 2020
    /// PMID 31761326). Magnitude uncalibrated; the direction is the result.
    #[test]
    fn escrt_repair_reduces_rsl3_kill() {
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 120,
        };
        let cond = Condition {
            name: "escrt_ab".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let run = |escrt: Option<(f64, f64)>| {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    escrt,
                    ..Default::default()
                },
            )
        };
        let baseline = run(None);
        // High repair rate + ample budget ⇒ many cells are resealed past the death
        // threshold for the run, so fewer die.
        let escrt_high = run(Some((0.8, 50.0)));
        assert!(
            escrt_high.overall_kill_rate < baseline.overall_kill_rate,
            "ESCRT repair must lower kill: escrt_high={} baseline={}",
            escrt_high.overall_kill_rate,
            baseline.overall_kill_rate
        );
    }

    /// #465: the `--snapshot=escrt` preset is wired (RSL3 with the membrane-repair
    /// brake enabled at a non-zero rate + budget).
    #[test]
    fn escrt_snapshot_preset_is_wired() {
        let p = resolve_snapshot("escrt");
        assert_eq!(p.name, "escrt");
        assert!(p.escrt, "the escrt preset must enable the repair brake");
        assert!(matches!(p.treatment, Treatment::RSL3));
    }

    /// #466 A/B: the POR/CYB5R1 enzymatic H2O2 source raises RSL3 kill. Injecting
    /// an enzymatic oxidant into each tumor cell's basal ROS feeds the Fenton/ROS
    /// pool, so an RSL3 run with POR on kills MORE than the no-POR baseline (more
    /// H2O2 ⇒ more ferroptosis; Yan 2021 PMID 33860083). O2-independent here
    /// (o2_dep = 0) to isolate the oxidant magnitude. Direction is the result.
    #[test]
    fn por_h2o2_raises_rsl3_kill() {
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 120,
        };
        let cond = Condition {
            name: "por_ab".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let run = |por: Option<(f64, f64)>| {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    por,
                    ..Default::default()
                },
            )
        };
        let baseline = run(None);
        let por_high = run(Some((0.6, 0.0))); // O2-independent to isolate magnitude
        assert!(
            por_high.overall_kill_rate > baseline.overall_kill_rate,
            "POR H2O2 must raise kill: por_high={} baseline={}",
            por_high.overall_kill_rate,
            baseline.overall_kill_rate
        );
    }

    /// #466: the `--snapshot=por` preset is wired (RSL3 with the O2-coupled POR
    /// H2O2 source enabled).
    #[test]
    fn por_snapshot_preset_is_wired() {
        let p = resolve_snapshot("por");
        assert_eq!(p.name, "por");
        assert!(p.por, "the por preset must enable the POR H2O2 source");
        assert!(matches!(p.treatment, Treatment::RSL3));
    }

    /// #467 A/B: the 7-DHC sterol radical-trapping defense reduces RSL3 kill. A high
    /// 7-DHC pool (DHCR7-low tumor) adds GPX4-independent radical-trapping quench,
    /// lowering the propagation rate, so an RSL3 run with 7-DHC on kills LESS than
    /// the baseline (ferroptosis resistance; Freitas/Li Nature 2024 PMID 38297130).
    /// Magnitude uncalibrated; the direction is the result.
    #[test]
    fn dhc7_radical_trap_reduces_rsl3_kill() {
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 120,
        };
        let cond = Condition {
            name: "dhc7_ab".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let run = |dhc7: Option<f64>| {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    dhc7,
                    ..Default::default()
                },
            )
        };
        let baseline = run(None);
        let dhc7_high = run(Some(1.0));
        assert!(
            dhc7_high.overall_kill_rate < baseline.overall_kill_rate,
            "7-DHC radical trap must lower kill: dhc7_high={} baseline={}",
            dhc7_high.overall_kill_rate,
            baseline.overall_kill_rate
        );
    }

    /// #467: the `--snapshot=dhc7` preset is wired (RSL3 on a DHCR7-low resistant tumor).
    #[test]
    fn dhc7_snapshot_preset_is_wired() {
        let p = resolve_snapshot("dhc7");
        assert_eq!(p.name, "dhc7");
        assert!(p.dhc7, "the dhc7 preset must enable the 7-DHC radical trap");
        assert!(matches!(p.treatment, Treatment::RSL3));
    }

    /// #483 A/B: the vitamin-K / VKORC1L1 radical-trapping defense reduces RSL3
    /// kill (a sixth GPX4-independent quench, so a VKORC1L1-high p53-competent
    /// tumor resists), and warfarin (which inhibits VKORC1L1) REVERSES it,
    /// restoring the kill toward the unprotected baseline (Yang et al. Cell Metab
    /// 2023 PMID 37467745). Magnitude uncalibrated; the direction is the result.
    #[test]
    fn warfarin_reverses_vkorc1l1_resistance() {
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 120,
        };
        let cond = Condition {
            name: "vitk_ab".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let run = |vitk: Option<(f64, f64)>| {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    vitk,
                    ..Default::default()
                },
            )
            .overall_kill_rate
        };
        let baseline = run(None);
        // VKORC1L1-defended, no warfarin: resists RSL3 (fewer kills).
        let defended = run(Some((1.0, 0.0)));
        assert!(
            defended < baseline,
            "VKORC1L1 vitamin-K trap must lower RSL3 kill: defended={defended} baseline={baseline}"
        );
        // Full warfarin inhibition collapses the trap, restoring the baseline kill
        // EXACTLY (effective trap 0). This is the druggable, repurposing-relevant
        // result: warfarin re-sensitizes a VKORC1L1-defended tumor to ferroptosis.
        let warfarin = run(Some((1.0, 1.0)));
        assert_eq!(
            warfarin, baseline,
            "full warfarin inhibition must restore the unprotected RSL3 kill"
        );
    }

    /// #483: the `--snapshot=vkorc1l1` preset is wired (RSL3 on a VKORC1L1-defended
    /// tumor, the no-warfarin resistant state).
    #[test]
    fn vkorc1l1_snapshot_preset_is_wired() {
        let p = resolve_snapshot("vkorc1l1");
        assert_eq!(p.name, "vkorc1l1");
        assert!(
            p.vitk,
            "the vkorc1l1 preset must enable the VKORC1L1 radical trap"
        );
        assert!(matches!(p.treatment, Treatment::RSL3));
    }

    /// #484 A/B: PROM2 / MVB-exosome labile-iron efflux reduces RSL3 kill. A
    /// PROM2-high tumor exports ferritin-bound iron in secreted exosomes,
    /// depleting the labile iron pool and starving the Fenton reaction (the
    /// OPPOSITE sign to ferritinophagy), so an RSL3 run with PROM2 efflux on kills
    /// LESS than the baseline (ferroptosis resistance; Brown et al. Dev Cell 2019
    /// PMID 31761539). Magnitude uncalibrated; the direction is the result.
    #[test]
    fn prom2_iron_efflux_reduces_rsl3_kill() {
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 120,
        };
        let cond = Condition {
            name: "prom2_ab".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let run = |prom2: Option<f64>| {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    prom2,
                    ..Default::default()
                },
            )
            .overall_kill_rate
        };
        let baseline = run(None);
        let prom2_high = run(Some(0.8));
        assert!(
            prom2_high < baseline,
            "PROM2 iron efflux must lower RSL3 kill: prom2_high={prom2_high} baseline={baseline}"
        );
        // efflux = 0 reproduces the baseline exactly (the byte-identity invariant).
        assert_eq!(run(Some(0.0)), baseline);
    }

    /// #484: the `--snapshot=prom2` preset is wired (RSL3 on a PROM2-high
    /// iron-exporting resistant tumor).
    #[test]
    fn prom2_snapshot_preset_is_wired() {
        let p = resolve_snapshot("prom2");
        assert_eq!(p.name, "prom2");
        assert!(p.prom2, "the prom2 preset must enable PROM2 iron efflux");
        assert!(matches!(p.treatment, Treatment::RSL3));
    }

    /// #485 A/B: a copper ionophore (elesclomol/disulfiram) depletes GSH + GPX4
    /// each step (copper overload), so an RSL3 run with copper on kills MORE than
    /// the baseline (ferroptosis-cuproptosis crosstalk; Gao et al. Cell Death Dis
    /// 2021 PMID 34390123), and an ATP7B-efflux-competent tumor RESISTS (efflux
    /// exports copper, restoring toward the baseline). Magnitude uncalibrated; the
    /// direction is the result.
    #[test]
    fn copper_ionophore_raises_rsl3_kill_and_atp7b_efflux_protects() {
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 120,
        };
        let cond = Condition {
            name: "copper_ab".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let run = |copper: Option<CopperConfig>| {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    copper,
                    ..Default::default()
                },
            )
            .overall_kill_rate
        };
        let baseline = run(None);
        // Copper ionophore in an ATP7B-low tumor: GSH/GPX4 depleted ⇒ more kill.
        let copper = run(Some(CopperConfig::literature()));
        assert!(
            copper > baseline,
            "copper ionophore must raise RSL3 kill: copper={copper} baseline={baseline}"
        );
        // ATP7B-efflux-competent tumor: copper exported ⇒ resists, fewer kills
        // than the efflux-low ionophore arm.
        let with_efflux = run(Some(CopperConfig {
            atp7b_efflux: 1.0,
            ..CopperConfig::literature()
        }));
        assert!(
            with_efflux < copper,
            "ATP7B efflux must protect: with_efflux={with_efflux} copper={copper}"
        );
        // Full efflux reduces to the baseline exactly (retention 1.0 on both pools).
        assert_eq!(with_efflux, baseline);
        // disabled() reproduces the baseline exactly (the byte-identity invariant).
        assert_eq!(run(Some(CopperConfig::disabled())), baseline);
    }

    /// #485: the `--snapshot=copper` preset is wired (RSL3 + a copper ionophore).
    #[test]
    fn copper_snapshot_preset_is_wired() {
        let p = resolve_snapshot("copper");
        assert_eq!(p.name, "copper");
        assert!(
            p.copper,
            "the copper preset must enable copper-ionophore crosstalk"
        );
        assert!(matches!(p.treatment, Treatment::RSL3));
    }

    /// #486 A/B: dietary PUFA above the saturable lipid-droplet (DGAT) buffer adds
    /// oxidizable substrate, so an RSL3 run with dietary PUFA over the buffer kills
    /// MORE than the baseline; an equal-or-larger buffer (no DGAT inhibition) stores
    /// the PUFA and reproduces the baseline (Dierge et al. Cell Metab 2021 PMID
    /// 34118189). Magnitude uncalibrated; the direction is the result.
    #[test]
    fn dietary_pufa_above_dgat_buffer_raises_rsl3_kill() {
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 120,
        };
        let cond = Condition {
            name: "dietary_pufa_ab".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let run = |dietary: Option<(f64, f64)>| {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    dietary_pufa: dietary,
                    ..Default::default()
                },
            )
            .overall_kill_rate
        };
        let baseline = run(None);
        // Dietary PUFA fully buffered (supply == buffer) => no excess => baseline.
        assert_eq!(
            run(Some((0.5, 0.5))),
            baseline,
            "dietary PUFA at or below the DGAT buffer must reproduce the baseline"
        );
        // Dietary PUFA over the buffer (DGAT-inhibited) => more kill.
        let excess = run(Some((0.6, 0.2)));
        assert!(
            excess > baseline,
            "dietary PUFA above the DGAT buffer must raise RSL3 kill: excess={excess} baseline={baseline}"
        );
        // disabled (both 0.0) reproduces the baseline exactly (byte-identity).
        assert_eq!(run(Some((0.0, 0.0))), baseline);
    }

    /// #486: the `--snapshot=dietary-pufa` preset is wired (RSL3 + dietary PUFA over
    /// the DGAT buffer on an acidic tumor).
    #[test]
    fn dietary_pufa_snapshot_preset_is_wired() {
        let p = resolve_snapshot("dietary-pufa");
        assert_eq!(p.name, "dietary-pufa");
        assert!(
            p.dietary_pufa,
            "the dietary-pufa preset must enable the dietary-PUFA layer"
        );
        assert!(
            p.ph_on,
            "dietary-pufa runs on an acidic tumor (the pH layer potentiates it)"
        );
        assert!(matches!(p.treatment, Treatment::RSL3));
    }

    /// #468: the model already expresses a Type-I (oxygen-INDEPENDENT) SDT radical
    /// arm through the #336 `sdt_o2_dependence` knob (the Type II fraction), no
    /// separate parameter needed. A Type-I-heavy sonosensitizer (dependence 0.3, so
    /// 70% Type I) retains MORE hypoxic-core SDT kill than a pure Type II agent
    /// (dependence 1.0, which collapses in the anoxic core). Demonstrated in a
    /// hypoxic config (small o2_lambda ⇒ a real hypoxic shell).
    #[test]
    fn type_i_sdt_retains_hypoxic_kill_where_type_ii_collapses() {
        let cfg = RunConfig {
            grid_dim: 60,
            n_steps: 180,
        };
        let cond = Condition {
            name: "sdt_typei".to_string(),
            treatment: Treatment::SDT,
            treatment_name: "SDT".to_string(),
            o2_lambda: Some(50.0), // steep gradient ⇒ a hypoxic core (§7.1 config)
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let run = |dep: f64| {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    sdt_o2_dependence: dep,
                    ..Default::default()
                },
            )
        };
        let type_ii = run(1.0); // fully O2-dependent ⇒ collapses in the hypoxic core
        let type_i_heavy = run(0.3); // 70% Type I (O2-independent) ⇒ retains hypoxic kill
                                     // At this config Type-I-heavy retains ~43% hypoxic kill while pure Type II
                                     // collapses to exactly 0 in the anoxic core.
        assert!(
            type_i_heavy.hypoxic_kill_rate > type_ii.hypoxic_kill_rate,
            "Type-I-heavy SDT must retain more hypoxic kill than pure Type II: \
             type_i={} type_ii={}",
            type_i_heavy.hypoxic_kill_rate,
            type_ii.hypoxic_kill_rate
        );
    }

    /// #468: the `--snapshot=sdt-typei` preset is wired (SDT, Type-I-heavy
    /// dependence 0.3), the hypoxia-tolerant complement of `sdt-o2dep` (dep 1.0).
    #[test]
    fn sdt_typei_snapshot_preset_is_wired() {
        let p = resolve_snapshot("sdt-typei");
        assert_eq!(p.name, "sdt-typei");
        assert!(matches!(p.treatment, Treatment::SDT));
        assert_eq!(p.sdt_o2_dependence, 0.3);
        // The complement: sdt-o2dep is the pure Type II (dep 1.0) collapse preset.
        assert_eq!(resolve_snapshot("sdt-o2dep").sdt_o2_dependence, 1.0);
    }

    /// Three independent realism layers enabled together — clonal subclones
    /// (#242) × Treg/MDSC suppressor (#264 P2) × multi-checkpoint brake (#264
    /// P3), all under immune — must compose: the run is deterministic, each
    /// layer reports a coherent metric (no field/seed collision silently
    /// dropping one), and — the directional cross-check — the checkpoint layer's
    /// modulation still *flows through* under composition (dual blockade
    /// out-kills anti-PD-1 alone even with clonal + suppressor also active).
    #[test]
    fn clonal_suppressor_checkpoints_compose() {
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 100,
        };
        let cond = Condition {
            name: "compose_csc".to_string(),
            treatment: Treatment::SDT,
            treatment_name: "SDT".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let immune = SpatialImmuneConfig {
            immune_kill_rate: 0.5,
            ..SpatialImmuneConfig::for_3d()
        };
        // Same clonal + suppressor composition; the checkpoint panel varies.
        let run = |panel: CheckpointPanel| {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    immune: Some(immune),
                    clonal: Some(ClonalConfig::literature_4()),
                    suppressor: Some(SuppressorConfig::enabled()),
                    checkpoints: Some(panel),
                    ..Default::default()
                },
            )
        };
        let tumor = CheckpointPanel::pd1_ctla4_tumor();
        let a = run(tumor.with_anti_pd1(0.8)); // anti-PD-1 alone (CTLA-4 braking)
        let b = run(tumor.with_anti_pd1(0.8)); // repeat for determinism
                                               // Deterministic across runs (no nondeterministic cross-layer coupling).
        assert_eq!(a.total_dead, b.total_dead, "composed run is deterministic");
        assert_eq!(a.ferroptosis_kills, b.ferroptosis_kills);
        assert_eq!(a.immune_kills, b.immune_kills);
        // Each layer reports a coherent metric (none silently dropped).
        assert_eq!(
            a.subclone_kills.as_ref().map_or(0, |v| v.len()),
            4,
            "clonal reports all 4 subclones"
        );
        assert!(
            a.suppressor_source_count.unwrap_or(0) > 0,
            "suppressor niches present under composition"
        );
        assert!(a.checkpoint_brake.is_some(), "checkpoint brake reported");
        let a_immune = a.immune_kills.expect("immune_on populates immune_kills");
        assert!(
            a_immune > 0,
            "immune killing actually fires under composition (not just reported); got {a_immune}"
        );
        assert!(a.total_dead > 0, "composition still kills some cells");
        // Directional cross-check: the checkpoint layer composes correctly —
        // dual blockade lifts CTLA-4 too, so it out-kills anti-PD-1 alone EVEN
        // with clonal + suppressor also active (the #264 P3 effect survives
        // composition rather than being masked by the other layers).
        let combo = run(tumor.with_anti_pd1(0.8).with_anti_ctla4(0.8));
        let combo_immune = combo
            .immune_kills
            .expect("immune_on populates immune_kills");
        assert!(
            combo_immune > a_immune,
            "dual blockade must out-kill anti-PD-1 alone under composition: \
             combo={combo_immune}, mono={a_immune}"
        );
    }

    /// Two grid-level layers enabled together — spheroid radial re-gen (#197,
    /// its own RNG + `Params::spheroid`) × explicit vasculature (#191, vessel
    /// supply replacing the edge-distance O2) — must compose deterministically
    /// (the supply field is built from the radially-regenerated grid).
    #[test]
    fn spheroid_vasculature_compose() {
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 90,
        };
        let cond = Condition {
            name: "compose_sv".to_string(),
            treatment: Treatment::SDT,
            treatment_name: "SDT".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let run = || {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    spheroid: Some(SpheroidConfig::literature()),
                    vasculature: Some(VasculatureConfig::well_vascularized()),
                    ..Default::default()
                },
            )
        };
        let a = run();
        let b = run();
        assert_eq!(
            a.total_dead, b.total_dead,
            "spheroid+vasculature composition is deterministic"
        );
        assert!(
            a.vascular_hypoxic_fraction.is_some(),
            "vasculature reports its hypoxic fraction under composition"
        );
        assert!(a.total_dead > 0, "composition still kills some cells");
    }

    // ===== Contact-layer composition + guards (#302) =====

    /// #302: contact + spheroid is the *intended, correct* pairing — the
    /// fixed-26 contact denominator is right for a centred spheroid whose tumor
    /// never touches the box face. Adding contact on top of the spheroid reduces
    /// RSL3 kills further, and the composition is deterministic.
    #[test]
    fn contact_composes_with_spheroid() {
        let cfg = RunConfig {
            grid_dim: 24,
            n_steps: 80,
        };
        let cond = Condition {
            name: "compose_cs".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let spheroid_only = run_one_condition_full(
            &cond,
            cfg,
            None,
            Overrides {
                spheroid: Some(SpheroidConfig::literature()),
                ..Default::default()
            },
        );
        let run_both = || {
            run_one_condition_full(
                &cond,
                cfg,
                None,
                Overrides {
                    spheroid: Some(SpheroidConfig::literature()),
                    contact: Some(ContactConfig::literature()),
                    ..Default::default()
                },
            )
        };
        let both = run_both();
        assert!(
            spheroid_only.total_dead > 0,
            "spheroid baseline must kill some cells"
        );
        assert!(
            both.total_dead < spheroid_only.total_dead,
            "contact must reduce RSL3 kills under a spheroid: spheroid={}, both={}",
            spheroid_only.total_dead,
            both.total_dead
        );
        assert_eq!(
            both.total_dead,
            run_both().total_dead,
            "spheroid+contact composition is deterministic"
        );
    }

    /// #302: slab + contact is mutually exclusive — the fixed-26 contact
    /// denominator mis-scores a slab's domain-boundary shell as low-contact. The
    /// guard is a `debug_assert!`; this `#[should_panic]` locks it (gated to
    /// debug builds, where the assert lives, so a `--release` run does not
    /// spuriously fail).
    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "mutually exclusive")]
    fn slab_plus_contact_panics() {
        let cfg = RunConfig {
            grid_dim: 12,
            n_steps: 1,
        };
        let cond = Condition {
            name: "slab_contact_guard".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let _ = run_one_condition_full(
            &cond,
            cfg,
            None,
            Overrides {
                slab: Some(SlabConfig::surface()),
                contact: Some(ContactConfig::literature()),
                ..Default::default()
            },
        );
    }

    /// #302: slab + spheroid is mutually exclusive (incompatible geometries —
    /// the spheroid's radial phenotype keys on a center an all-tumor block lacks).
    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "mutually exclusive")]
    fn slab_plus_spheroid_panics() {
        let cfg = RunConfig {
            grid_dim: 12,
            n_steps: 1,
        };
        let cond = Condition {
            name: "slab_spheroid_guard".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let _ = run_one_condition_full(
            &cond,
            cfg,
            None,
            Overrides {
                slab: Some(SlabConfig::surface()),
                spheroid: Some(SpheroidConfig::literature()),
                ..Default::default()
            },
        );
    }
}
