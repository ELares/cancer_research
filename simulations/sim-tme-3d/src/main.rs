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

use ferroptosis_core::biochem::{exo_decay_factor, sim_cell_step, CellState};
use ferroptosis_core::cell::{Phenotype, Treatment};
use ferroptosis_core::clonal::{assign_subclones_3d, repopulate_dead_sites_3d, ClonalConfig};
use ferroptosis_core::contact::{
    apply_contact_resistance_3d, apply_contact_resistance_at_3d, ContactConfig,
};
use ferroptosis_core::dose_schedule::DoseSchedule;
use ferroptosis_core::grid::{TumorGrid3D, TUMOR_RADIUS_FRACTION};
use ferroptosis_core::immune_spatial::{
    dc_activation, diffuse_damp_3d_step, exhaustion_factor, immune_kill_probability,
    suppressor_kill_multiplier, suppressor_source_mask_3d, CheckpointPanel, DcSubsetConfig,
    SuppressorConfig, DAMP_KILL_THRESHOLD,
};
use ferroptosis_core::nutrient::{apply_nutrient_stress_3d, NutrientConfig};
use ferroptosis_core::oxygen::{o2_dependent_exo_factor, radial_o2_field};
use ferroptosis_core::params::{
    Params, PersisterConfig, PhConfig, SpatialImmuneConfig, SpatialParams, StromalConfig,
};
use ferroptosis_core::persister;
use ferroptosis_core::ph::{ion_trap_factor_from_ph, iron_multiplier_from_ph, radial_ph_field};
use ferroptosis_core::physics::local_ros_multiplier_3d;
use ferroptosis_core::slab::{
    apply_depth_graded_cells_3d, scale_interpretation, slab_supply_field, SlabConfig,
    SlabPhenotypeConfig, KROGH_LAMBDA_UM,
};
use ferroptosis_core::spheroid::{
    apply_radial_cells_3d, radial_fraction_3d, radial_mufa_protection, SpheroidConfig,
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
    /// Dendritic-cell subset mix (#264 Phase 4). `None` / balanced ⇒ priming
    /// efficiency 1.0 ⇒ no immune-kill modulation ⇒ byte-identical.
    dc_subsets: Option<DcSubsetConfig>,
    /// Oxygen-dependent SDT/PDT exo-ROS yield (#336): the "Type II fraction" of
    /// the exogenous ROS that scales with local O2. `0.0` (default) ⇒ fully
    /// O2-independent (the historical optimistic upper bound) ⇒ the exo-ROS
    /// factor is exactly 1.0 ⇒ byte-identical. `1.0` ⇒ fully Type II /
    /// O2-dependent, so SDT loses efficacy in hypoxic zones like the clinical
    /// SONALA-001 agent (manuscript §7.1).
    sdt_o2_dependence: f64,
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
    let slab_cfg = overrides.slab;
    // Depth-graded slab phenotype (#272): only meaningful with a slab grid.
    let slab_phenotype_cfg = overrides.slab_phenotype.filter(|_| slab_cfg.is_some());
    let contact_cfg = overrides.contact;
    let nutrient_cfg = overrides.nutrient;
    let dc_subsets_cfg = overrides.dc_subsets;
    // Oxygen-dependent SDT/PDT exo-ROS (#336). `0.0` (default/matrix) ⇒ the
    // exo-ROS O2 factor is exactly 1.0 ⇒ byte-identical.
    let sdt_o2_dependence = overrides.sdt_o2_dependence;
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
    let params = if spheroid_cfg.is_some() {
        Params::spheroid()
    } else {
        Params::default()
    };
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
        apply_radial_cells_3d(&mut grid, cfg, SPHEROID_SEED);
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
            let vessel = vessel_supply_field(&grid, v, lambda);
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
            (Some(v), Some(lambda)) => Some(vessel_supply_field(&grid, v, lambda)),
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
    let mut ferroptosis_kills = 0usize;
    let mut immune_kills = 0usize;

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
    let dc_priming = if dc_subsets_on {
        dc_subsets_cfg.map_or(1.0, |c| c.priming_efficiency())
    } else {
        1.0
    };

    // The rich kill path handles any realism layer (each factor is identity when
    // its layer is off); the default allocation-free path runs only when ALL are
    // off, staying byte-identical to pre-#243.
    let realism_kill_path = exhaustion_on || suppressor_on || dc_subsets_on;

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

                let died =
                    sim_cell_step(&mut gc.state, &gc.cell, &params, step, extra_iron, &mut rng);
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
                        // Competing-rate update (#262): acquisition and
                        // reversion both act each step (not an either-or keyed
                        // on drug == 0), so sustained sub-saturating drug
                        // reaches a sub-cap equilibrium rather than ratcheting
                        // monotonically to the cap.
                        gc.state.persister_fraction = persister::step(frac, drug_intensity, pcfg);
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
                            // `dc_priming` is a uniform scalar (the cDC1/cDC2
                            // mix, #264 Phase 4): 1.0 when off, < 1.0 for a
                            // cDC1-poor tumor that primes killing less efficiently.
                            let kill_prob = immune_kill_probability(
                                activation,
                                immune_cfg.immune_kill_rate,
                                effective_brake,
                            ) * exh
                                * supp
                                * dc_priming;
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
                    for idx in revived {
                        apply_contact_resistance_at_3d(&mut grid, idx, contact);
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
        subclone_kills,
        vascular_hypoxic_fraction,
        scale_interpretation: scale_interpretation_str,
        suppressor_source_count,
        suppressor_peak,
        checkpoint_brake,
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
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
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
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
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
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
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
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
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
        clonal: true,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
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
        clonal: false,
        vasculature: true,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
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
        clonal: false,
        vasculature: false,
        spheroid: true,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
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
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: true,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
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
        clonal: false,
        vasculature: true,
        spheroid: false,
        slab: true,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
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
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: true,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: false,
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
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: true,
        contact: false,
        nutrient: false,
        dc_subsets: false,
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
        clonal: true,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: true,
        checkpoints: true,
        contact: false,
        nutrient: false,
        dc_subsets: false,
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
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: true,
        nutrient: false,
        dc_subsets: false,
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
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: true,
        dc_subsets: false,
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
        clonal: false,
        vasculature: false,
        spheroid: false,
        slab: false,
        suppressor: false,
        checkpoints: false,
        contact: false,
        nutrient: false,
        dc_subsets: true,
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
            let vessel = vessel_supply_field(&g, &vessels, ZONE_REF_LAMBDA);
            let planar = slab_supply_field(&g, scfg.depth_offset_mm * 1000.0, ZONE_REF_LAMBDA);
            combine_supply_max(&planar, &vessel)
        } else {
            let vessels = place_vessels_3d(&snapshot_grid, &cfg, VESSEL_SEED);
            vessel_supply_field(&snapshot_grid, &vessels, ZONE_REF_LAMBDA)
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
    let mut buffers =
        snapshot::SnapshotBuffers::new(run_cfg.grid_dim, run_cfg.n_steps, preset.persister);
    let result = run_one_condition_full(
        &condition,
        run_cfg,
        Some(&mut buffers),
        Overrides {
            persister: preset.persister.then(PersisterConfig::enabled),
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
