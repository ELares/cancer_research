//! Core ferroptosis biochemistry engine.
//!
//! Provides both the single-step function (for spatial model interleaving)
//! and the full 180-step loop (for single-cell simulations).

use rand::prelude::*;
use serde::{Deserialize, Serialize};

use crate::cell::{norm, Cell, Treatment};
use crate::params::Params;

/// Mutable state carried between timesteps.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CellState {
    pub gsh: f64,
    pub gpx4: f64,
    pub fsp1: f64,
    pub mufa_protection: f64,
    pub lp: f64,
    pub dead: bool,
    pub death_step: Option<u32>,
    pub exo_ros_peak: f64,
    /// Drug-tolerant persister fraction ∈ [0, 1] (#241). Consumer-owned:
    /// `sim_cell_step` never reads or writes it (the core engine stays
    /// byte-identical), so it is `0.0` for every code path that does not
    /// opt into the persister model. A consumer mutates it via
    /// [`crate::persister`] helpers around the step call.
    ///
    /// `#[serde(default)]` so older `CellState` JSON (written before #241)
    /// still deserializes, defaulting to `0.0` (the inert value).
    #[serde(default)]
    pub persister_fraction: f64,
    /// Reversible/locked persister sub-pools + the sustained-exposure tracker
    /// for the epigenetic-locking model (#342). Consumer-owned exactly like
    /// `persister_fraction`: `sim_cell_step` never touches them, and a consumer
    /// evolves them via [`crate::persister::step_with_locking`], caching their
    /// [`crate::persister::PersisterState::total`] back into `persister_fraction`
    /// so the existing biochem couplings keep reading one scalar. All `0.0` (the
    /// inert value) for every path that does not opt into locking, so the
    /// production matrix stays byte-identical. `#[serde(default)]` so pre-#342
    /// `CellState` JSON still deserializes.
    #[serde(default)]
    pub persister_reversible: f64,
    #[serde(default)]
    pub persister_locked: f64,
    #[serde(default)]
    pub persister_cum_exposure: f64,
    /// ESCRT-III membrane-repair budget consumed so far (#465). Read/written ONLY
    /// inside the death-threshold-crossing brake, and only when
    /// `Params::escrt_repair_rate > 0.0`; on the default path (`escrt_repair_rate
    /// == 0.0`) the brake is never entered, so this field is never touched and the
    /// engine stays byte-identical. `#[serde(default)]` (`0.0`) so pre-#465
    /// `CellState` JSON still deserializes.
    #[serde(default)]
    pub escrt_budget_used: f64,
}

impl CellState {
    /// Initialize cell state from a Cell and Treatment.
    /// Applies RSL3 GPX4 inhibition and samples exogenous ROS peak.
    pub fn from_cell(cell: &Cell, tx: Treatment, params: &Params, rng: &mut StdRng) -> Self {
        let mut gpx4 = cell.gpx4;
        let exo_ros_peak: f64 = match tx {
            Treatment::Control | Treatment::RSL3 => 0.0,
            Treatment::SDT => norm(rng, params.sdt_ros, 1.0).max(0.0),
            Treatment::PDT => norm(rng, params.pdt_ros, 1.0).max(0.0),
        };
        if let Treatment::RSL3 = tx {
            gpx4 *= 1.0 - params.rsl3_gpx4_inhib;
        }
        CellState {
            gsh: cell.gsh,
            gpx4,
            fsp1: cell.fsp1,
            // #339: an acute/naive start (Some) overrides the established
            // steady-state start; None ⇒ initial_mufa_protection ⇒ byte-identical.
            mufa_protection: params
                .mufa_acute_start
                .unwrap_or(params.initial_mufa_protection),
            lp: 0.0,
            dead: false,
            death_step: None,
            exo_ros_peak,
            persister_fraction: 0.0,
            persister_reversible: 0.0,
            persister_locked: 0.0,
            persister_cum_exposure: 0.0,
            escrt_budget_used: 0.0,
        }
    }

    /// Initialize with a custom exogenous ROS peak (for spatial model where
    /// ROS dose depends on depth/position).
    ///
    /// Applies RSL3's one-shot GPX4 knockdown at init (the steady-state
    /// model: drug is present from t=0 and immediately inhibits GPX4). This
    /// is the historical behavior; callers using a time-varying
    /// [`crate::dose_schedule::DoseSchedule`] for RSL3 should instead use
    /// [`from_cell_with_ros_opts`](Self::from_cell_with_ros_opts) with
    /// `apply_rsl3_init = false` and apply per-step inactivation themselves.
    pub fn from_cell_with_ros(
        cell: &Cell,
        tx: Treatment,
        params: &Params,
        exo_ros_peak: f64,
    ) -> Self {
        Self::from_cell_with_ros_opts(cell, tx, params, exo_ros_peak, true)
    }

    /// Like [`from_cell_with_ros`](Self::from_cell_with_ros) but with
    /// explicit control over the RSL3 one-shot init knockdown.
    ///
    /// `apply_rsl3_init = true` reproduces `from_cell_with_ros` exactly
    /// (byte-identical). `apply_rsl3_init = false` skips the
    /// `gpx4 *= 1 - rsl3_gpx4_inhib` step so a time-varying dose schedule
    /// can drive GPX4 inactivation per step instead of one-shot at t=0
    /// (the `tumor_pk::sim_cell_with_pk` model, #239). Has no effect for
    /// non-RSL3 treatments.
    pub fn from_cell_with_ros_opts(
        cell: &Cell,
        tx: Treatment,
        params: &Params,
        exo_ros_peak: f64,
        apply_rsl3_init: bool,
    ) -> Self {
        let mut gpx4 = cell.gpx4;
        if apply_rsl3_init {
            if let Treatment::RSL3 = tx {
                gpx4 *= 1.0 - params.rsl3_gpx4_inhib;
            }
        }
        CellState {
            gsh: cell.gsh,
            gpx4,
            fsp1: cell.fsp1,
            // #339: an acute/naive start (Some) overrides the established
            // steady-state start; None ⇒ initial_mufa_protection ⇒ byte-identical.
            mufa_protection: params
                .mufa_acute_start
                .unwrap_or(params.initial_mufa_protection),
            lp: 0.0,
            dead: false,
            death_step: None,
            exo_ros_peak,
            persister_fraction: 0.0,
            persister_reversible: 0.0,
            persister_locked: 0.0,
            persister_cum_exposure: 0.0,
            escrt_budget_used: 0.0,
        }
    }
}

/// Accumulate MUFA-style lipid-remodeling protection against peroxidation.
///
/// SCD1 (the enzyme converting SFA→MUFA) is regulated by SREBP1/mTORC1,
/// NOT by NRF2. In 3D culture and in vivo, SCD1 is constitutively active
/// and enriches membranes with MUFAs that displace PUFAs, reducing
/// ferroptosis susceptibility (Dixon/Park, Cancer Research 2025;
/// Magtanong et al., Cell Chem Biol 2019; Tesfay et al., Cancer Res 2019).
///
/// The rate is context-dependent: zero in 2D culture (Params::default),
/// non-zero in in-vivo-like conditions (Params::invivo).
/// `mufa_max` is the MUFA carrying capacity: the per-cell `cell.mufa_cap` when
/// set (the spheroid's radial cap, #270), else the global `params.scd_mufa_max`.
/// The logistic growth saturates at — and the value is clamped to — `mufa_max`,
/// so a per-cell cap yields a per-cell steady state (durable position-dependent
/// MUFA). `cell.mufa_cap = None` ⇒ `params.scd_mufa_max` ⇒ byte-identical.
///
/// `rate` is the per-cell SCD1/MUFA accumulation rate: the per-cell
/// `cell.mufa_rate` when set (a phenotype-specific rate, #363), else the global
/// `params.scd_mufa_rate`. It scales how fast MUFA protection accumulates toward
/// `mufa_max` (the cap sets the steady state, the rate sets the time constant).
/// `cell.mufa_rate = None` ⇒ `params.scd_mufa_rate` ⇒ byte-identical.
#[inline]
fn update_mufa_protection(current: f64, mufa_max: f64, rate: f64, params: &Params) -> f64 {
    let growth = rate * (1.0 - current / (mufa_max + 1e-9));
    let decay = params.scd_mufa_decay * current;
    (current + growth - decay).clamp(0.0, mufa_max.max(0.0))
}

/// Ether-linked PUFA augmentation of the peroxidizable lipid substrate (#339).
///
/// Polyunsaturated ether phospholipids (ether-PUFA-PE) are an additional pool
/// of peroxidation-vulnerable membrane lipid whose synthesis (FAR1/AGPS, then
/// the peroxisomal ether-lipid pathway) PROMOTES ferroptosis: more ether-PUFA
/// substrate means more lipid-peroxide accumulation (Zou et al., Nature 2020,
/// PMID 32939090; Cui et al., Cell Death Differ 2021, PMID 33731874).
/// Conversely, FAR1/AGPS loss (or the in-vivo downregulation Zou 2020
/// documents) shrinks this pool and confers ferroptosis RESISTANCE, which is
/// the modeled escape route. NOTE this is the opposite sign to the loose
/// "shift to ether-PE = escape" framing: the robust literature direction is
/// that ether-PUFA PROMOTES ferroptosis and its LOSS is the escape.
///
/// `ether_pufa_fraction` is the ether-PUFA pool as a fraction of the base PUFA
/// substrate (domain `>= 0`), so the peroxidizable PUFA term is scaled by
/// `1 + max(fraction, 0)`. `0.0` (default) ⇒ ×1.0 ⇒ byte-identical (FAR1/AGPS-
/// null is the `0` limit). The `max(_, 0)` floor keeps an out-of-contract
/// negative fraction from shrinking the substrate (the two post-death
/// `effective_unsat` sites have no `.max(0.05)` floor of their own); it is a
/// no-op for the documented `>= 0` domain, so the default stays byte-identical.
///
/// The plasmalogen / TMEM189 (PEDS1) vinyl-ether sub-step is deliberately NOT
/// folded in: its sign is genuinely contested (protective via a FAR1
/// degradation feedback in Cui 2021 and as a radical sink in Zoeller 1999,
/// PMID 10051451; dispensable for the pro-ferroptotic effect in Zou 2020;
/// sensitizing-on-loss in C. elegans, Perez 2022, PMID 36178986), so baking a
/// fixed direction would overstate certainty.
#[inline]
fn ether_augmented_pufa(lipid_unsat: f64, params: &Params) -> f64 {
    // The ether-lipid pool (#339), the MCFA→ACSL4/CD36 boost (#446), and the
    // tumor-intrinsic ACSL4-status boost (#444) all set the oxidizable-PUFA level,
    // so they augment the substrate additively. The ether/MCFA boosts are
    // protective-only floored at 0; the ACSL4-status boost is SIGNED (it reaches
    // `-1` for ACSL4-negative tumors, collapsing the PUFA substrate ⇒
    // ferroptosis-refractory), so it is added raw and the whole augmentation is
    // floored at 0. All boosts `0.0` by default ⇒ ×1.0 ⇒ byte-identical.
    lipid_unsat
        * (1.0
            + params.ether_pufa_fraction.max(0.0)
            + params.mcfa_pufa_boost.max(0.0)
            + params.acsl4_status_boost)
            .max(0.0)
}

/// Total MUFA-style ferroptosis protection at a peroxidation site: the dynamic
/// SCD1-driven `mufa_protection` plus the constant MBOAT1/2 hormone-regulated
/// MUFA enrichment (#339). MBOAT1 (ER-regulated) and MBOAT2 (AR-regulated)
/// remodel phospholipids toward MUFA-PE and suppress ferroptosis independently
/// of GPX4 (Liang et al., Cell 2023, PMID 37267948), so they act as a second,
/// constitutive MUFA-enrichment source layered onto the SCD1 dynamics. The
/// boost is floored at `0` because MBOAT is an enrichment (protective-only)
/// source. `mboat_mufa_boost = 0.0` (default) ⇒ unchanged ⇒ byte-identical.
/// The result can exceed the SCD1 cap; the consuming `(1 - protection)` term is
/// floored by the existing `.max(0.05)` substrate minimum, so full protection
/// leaves a small residual peroxidizable pool rather than going negative.
#[inline]
fn total_mufa_protection(dynamic: f64, params: &Params) -> f64 {
    dynamic + params.mboat_mufa_boost.max(0.0)
}

/// NCOA4-ferritinophagy labile-iron release factor (#340).
///
/// The static-iron model holds `cell.iron` fixed. In reality NCOA4-mediated
/// autophagy of ferritin (ferritinophagy) releases stored iron into the labile
/// pool over time, feeding Fenton chemistry (Mancias et al., Nature 2014, PMID
/// 24695223; Hou et al., Autophagy 2016, PMID 27245739). This scales the Fenton
/// iron by a time-dependent factor that ramps from `1.0` (step 0) toward
/// `1 + ferritinophagy_release` with time constant `ferritinophagy_tau` steps,
/// representing the gradual rise in labile iron as ferritinophagy proceeds.
///
/// `ferritinophagy_release = 0.0` (default) returns exactly `1.0` for every step
/// (a fast path, so `iron * factor == iron` bit-for-bit), keeping the production
/// matrix byte-identical. The release is floored at `0` (ferritinophagy adds
/// labile iron, it does not remove it).
#[inline]
pub fn ferritinophagy_iron_factor(step: u32, params: &Params) -> f64 {
    let release = params.ferritinophagy_release.max(0.0);
    if release == 0.0 {
        return 1.0;
    }
    let tau = params.ferritinophagy_tau.max(1e-9);
    1.0 + release * (1.0 - (-(step as f64) / tau).exp())
}

/// Deterministic exogenous-ROS decay envelope for the post-bolus phase.
///
/// Models singlet-oxygen / exogenous-ROS decay after a single treatment
/// bolus: `1.0` for `step < 30` (the pre-decay plateau, where
/// [`sim_cell_step`] instead applies multiplicative noise), then
/// `0.5^((step-30)/15)` — a 15-step half-life decay.
///
/// Exposed so a time-varying [`crate::dose_schedule::DoseSchedule`]
/// consumer can **divide it out**: for multi-dose SDT/PDT the schedule's
/// own per-dose rise+decay is the availability envelope, and this
/// single-bolus envelope (keyed to run start) would otherwise
/// double-count decay for later doses (#239).
///
/// **Contract (load-bearing):** [`sim_cell_step`] applies *exactly* this
/// factor to `exo_ros_peak` for `step >= 30`, and the dosed SDT/PDT path
/// in `sim-tme-3d` divides *exactly* this factor back out. Both the
/// producer (here, via `sim_cell_step`) and the consumer (the binary)
/// must call this one function so they cannot drift. **Do not inline the
/// `0.5^((step-30)/15)` formula at either site** — if the envelope shape
/// ever changes, both ends must change together, which only happens if
/// they share this function.
#[inline]
pub fn exo_decay_factor(step: u32) -> f64 {
    if step < 30 {
        1.0
    } else {
        0.5_f64.powf((step - 30) as f64 / 15.0)
    }
}

/// Execute a single timestep of the ferroptosis biochemistry.
///
/// Returns `true` if the cell died this step.
///
/// `extra_iron` is additional iron from neighbor deaths (spatial model).
/// For single-cell simulations, pass 0.0.
pub fn sim_cell_step(
    state: &mut CellState,
    cell: &Cell,
    params: &Params,
    step: u32,
    extra_iron: f64,
    rng: &mut StdRng,
) -> bool {
    if state.dead {
        // Post-death LP accumulation: the autocatalytic chain reaction
        // continues for post_death_steps after the threshold is crossed.
        // No repair (cell defenses have failed), only ROS → LP.
        if let Some(ds) = state.death_step {
            if step < ds + params.post_death_steps {
                let effective_iron =
                    (cell.iron + extra_iron) * ferritinophagy_iron_factor(step, params);
                let fenton = effective_iron * params.fenton_rate * norm(rng, 1.0, 0.08).max(0.0);
                let exo = if step < 30 {
                    state.exo_ros_peak * norm(rng, 1.0, 0.1).max(0.0)
                } else {
                    state.exo_ros_peak * exo_decay_factor(step)
                };
                // POR/CYB5R1 enzymatic NAD(P)H/O2-driven H2O2 (#466): an extra oxidant
                // source feeding the ROS pool. `0.0` default ⇒ unchanged ⇒ byte-identical.
                let total_ros = cell.basal_ros + exo + fenton + params.por_h2o2_rate.max(0.0);
                let effective_unsat = ether_augmented_pufa(cell.lipid_unsat, params); // no MUFA protection
                let lp_direct = total_ros * effective_unsat * params.lp_rate;
                // ALOX enzymatic capacity (#446) gates propagation in death too:
                // the lipoxygenase machinery is what oxidizes the membrane PUFA,
                // independent of the (already-failed) antioxidant defenses, so the
                // same `1 + alox_propagation_boost` multiplier applies here as in
                // life. `0.0` default ⇒ ×1.0 ⇒ byte-identical.
                let alox_mul = (1.0 + params.alox_propagation_boost).max(0.0);
                let lp_propagation = state.lp * effective_unsat * params.lp_propagation * alox_mul;
                state.lp += lp_direct + lp_propagation;
            }
        }
        return false;
    }

    // === ROS SOURCES ===
    let effective_iron = (cell.iron + extra_iron) * ferritinophagy_iron_factor(step, params);
    let fenton = effective_iron * params.fenton_rate * norm(rng, 1.0, 0.08).max(0.0);
    let exo = if step < 30 {
        state.exo_ros_peak * norm(rng, 1.0, 0.1).max(0.0)
    } else {
        state.exo_ros_peak * exo_decay_factor(step)
    };
    // POR/CYB5R1 enzymatic NAD(P)H/O2-driven H2O2 (#466): an extra oxidant source
    // feeding the ROS pool. `0.0` default ⇒ unchanged ⇒ byte-identical.
    let total_ros = cell.basal_ros + exo + fenton + params.por_h2o2_rate.max(0.0);

    // === GSH SCAVENGING (Michaelis-Menten, NO artificial cap) ===
    let gsh_fraction = state.gsh / (state.gsh + params.gsh_km);
    let scavenged = total_ros * params.gsh_scav_efficiency * gsh_fraction;
    state.gsh -= scavenged * 0.5;
    state.gsh = state.gsh.max(0.0);

    // === NRF2-DRIVEN GSH RESYNTHESIS ===
    let deficit_fraction = ((params.gsh_max - state.gsh) / params.gsh_max).max(0.0);
    state.gsh += cell.nrf2 * params.nrf2_gsh_rate * deficit_fraction;

    // === LIPID PEROXIDATION ===
    let unscav = (total_ros - scavenged).max(0.0);
    state.mufa_protection = update_mufa_protection(
        state.mufa_protection,
        cell.mufa_cap.unwrap_or(params.scd_mufa_max),
        cell.mufa_rate.unwrap_or(params.scd_mufa_rate),
        params,
    );

    let effective_unsat = (ether_augmented_pufa(cell.lipid_unsat, params)
        * (1.0 - total_mufa_protection(state.mufa_protection, params)))
    .max(0.05);
    let lp_direct = unscav * effective_unsat * params.lp_rate;
    // AUTOCATALYTIC PROPAGATION — GSH-gated bistable switch.
    // GCH1/BH4 (#338) adds GPX4-independent radical-trapping quench capacity
    // (`gch1_rate`, 0.0 by default ⇒ byte-identical).
    let antioxidant_quench =
        state.gpx4 * (state.gsh / (state.gsh + 0.5)) + state.fsp1 + params.gch1_rate;
    let propagation_rate = params.lp_propagation / (1.0 + antioxidant_quench * 5.0);
    // ALOX isoform-specific enzymatic-oxidation capacity (#446): scale the
    // propagation rate by `1 + alox_propagation_boost` (clamped >= 0), so an
    // ALOX-high tumor propagates faster and an ALOX-null one (boost -1) has no
    // enzymatic propagation. `0.0` default ⇒ ×1.0 ⇒ byte-identical.
    let alox_mul = (1.0 + params.alox_propagation_boost).max(0.0);
    let lp_propagation = state.lp * effective_unsat * propagation_rate * alox_mul;
    let lp_generation = lp_direct + lp_propagation;

    // === REPAIR ===
    let gpx4_repair = state.gpx4
        * (state.gsh / (state.gsh + 1.0))
        * params.gpx4_rate
        * (state.lp / (state.lp + 0.5));
    let fsp1_repair = state.fsp1 * params.fsp1_rate * (state.lp / (state.lp + 0.5));
    // DHODH (#338): GPX4-independent CoQ10 reduction, an extra repair term in
    // parallel to FSP1 (`dhodh_rate`, 0.0 by default ⇒ byte-identical).
    let dhodh_repair = params.dhodh_rate * (state.lp / (state.lp + 0.5));
    let total_repair = gpx4_repair + fsp1_repair + dhodh_repair;

    state.lp = (state.lp + lp_generation - total_repair).max(0.0);

    // === GPX4 DYNAMIC REGULATION ===
    if total_ros > 1.0 {
        state.gpx4 -= params.gpx4_degradation_by_ros * (total_ros - 1.0);
    }
    let gpx4_target = cell.nrf2 * params.gpx4_nrf2_target_multiplier;
    state.gpx4 += params.gpx4_nrf2_upregulation * (gpx4_target - state.gpx4);
    state.gpx4 = state.gpx4.max(0.0);

    // Small noise
    state.lp += norm(rng, 0.0, 0.003);
    state.lp = state.lp.max(0.0);

    // Death check. ESCRT-III membrane repair (#465) can rescue a cell that has
    // crossed the threshold, for a finite per-cell budget. The RNG roll is drawn
    // ONLY when `escrt_repair_rate > 0` and budget remains (the `escrt_can_attempt`
    // short-circuit), so the default path (`escrt_repair_rate == 0.0`) makes no
    // extra draw and stays byte-identical.
    if state.lp > params.death_threshold {
        if crate::repair::escrt_can_attempt(
            params.escrt_repair_rate,
            state.escrt_budget_used,
            params.escrt_repair_budget,
        ) && crate::repair::escrt_rescue(params.escrt_repair_rate, rng.gen::<f64>())
        {
            // Resealed: consume one repair event and survive this step. The cell's
            // defenses still run next step, so ESCRT buys time for GSH/GPX4 to
            // recover under transient stress while sustained GPX4 inhibition still
            // kills once the budget is spent (a finite delay).
            state.escrt_budget_used += 1.0;
            return false;
        }
        state.dead = true;
        state.death_step = Some(step);
        return true;
    }

    false
}

/// Full 180-step simulation for a single cell.
/// Returns (is_dead, final_lp, final_gsh, final_gpx4).
///
/// This retains the original structure of the v3 engine but now optionally
/// includes a generic in vivo-like MUFA protection term when the corresponding
/// params are non-zero.
pub fn sim_cell(
    cell: &Cell,
    tx: Treatment,
    params: &Params,
    rng: &mut StdRng,
) -> (bool, f64, f64, f64) {
    let mut gsh = cell.gsh;
    let mut gpx4 = cell.gpx4;
    let fsp1 = cell.fsp1;
    let mut death_step: Option<u32> = None;
    // ESCRT-III repair budget consumed so far (#465); only touched when
    // `escrt_repair_rate > 0`, so the default path is byte-identical.
    let mut escrt_budget_used = 0.0_f64;
    // #339: acute/naive MUFA start override; None ⇒ byte-identical.
    let mut mufa_protection = params
        .mufa_acute_start
        .unwrap_or(params.initial_mufa_protection);
    let mut lp: f64 = 0.0;

    // Treatment: exogenous ROS
    let exo_ros_peak: f64 = match tx {
        Treatment::Control | Treatment::RSL3 => 0.0,
        Treatment::SDT => norm(rng, params.sdt_ros, 1.0).max(0.0),
        Treatment::PDT => norm(rng, params.pdt_ros, 1.0).max(0.0),
    };

    // Treatment: GPX4 inhibition
    if let Treatment::RSL3 = tx {
        gpx4 *= 1.0 - params.rsl3_gpx4_inhib;
    }

    for step in 0..180_u32 {
        // === ROS SOURCES (used by both alive and post-death paths) ===
        let fenton = cell.iron
            * ferritinophagy_iron_factor(step, params)
            * params.fenton_rate
            * norm(rng, 1.0, 0.08).max(0.0);
        let exo = if step < 30 {
            exo_ros_peak * norm(rng, 1.0, 0.1).max(0.0)
        } else {
            exo_ros_peak * exo_decay_factor(step)
        };
        // POR/CYB5R1 enzymatic NAD(P)H/O2-driven H2O2 (#466): an extra oxidant
        // source feeding the ROS pool. `0.0` default ⇒ unchanged ⇒ byte-identical.
        let total_ros = cell.basal_ros + exo + fenton + params.por_h2o2_rate.max(0.0);

        if death_step.is_some() {
            // Post-death: LP-only accumulation (no GSH, no repair, no GPX4).
            // Break check BEFORE accumulation to match sim_cell_step, which
            // only accumulates when step < death_step + post_death_steps.
            if step >= death_step.unwrap() + params.post_death_steps {
                break;
            }
            let effective_unsat = ether_augmented_pufa(cell.lipid_unsat, params);
            let lp_direct = total_ros * effective_unsat * params.lp_rate;
            // ALOX enzymatic capacity gates post-death propagation too (#446);
            // `0.0` default ⇒ ×1.0 ⇒ byte-identical.
            let alox_mul = (1.0 + params.alox_propagation_boost).max(0.0);
            let lp_prop = lp * effective_unsat * params.lp_propagation * alox_mul;
            lp += lp_direct + lp_prop;
            continue;
        }

        // === GSH SCAVENGING ===
        let gsh_fraction = gsh / (gsh + params.gsh_km);
        let scavenged = total_ros * params.gsh_scav_efficiency * gsh_fraction;
        gsh -= scavenged * 0.5;
        gsh = gsh.max(0.0);

        // === NRF2-DRIVEN GSH RESYNTHESIS ===
        let deficit_fraction = ((params.gsh_max - gsh) / params.gsh_max).max(0.0);
        gsh += cell.nrf2 * params.nrf2_gsh_rate * deficit_fraction;

        // === LIPID PEROXIDATION ===
        let unscav = (total_ros - scavenged).max(0.0);
        mufa_protection = update_mufa_protection(
            mufa_protection,
            cell.mufa_cap.unwrap_or(params.scd_mufa_max),
            cell.mufa_rate.unwrap_or(params.scd_mufa_rate),
            params,
        );
        let effective_unsat = (ether_augmented_pufa(cell.lipid_unsat, params)
            * (1.0 - total_mufa_protection(mufa_protection, params)))
        .max(0.05);
        let lp_direct = unscav * effective_unsat * params.lp_rate;
        // GCH1/BH4 (#338) adds GPX4-independent quench (`gch1_rate`, 0.0 default).
        let antioxidant_quench = gpx4 * (gsh / (gsh + 0.5)) + fsp1 + params.gch1_rate;
        let propagation_rate = params.lp_propagation / (1.0 + antioxidant_quench * 5.0);
        // ALOX isoform enzymatic-oxidation capacity (#446): same `1 + boost`
        // multiplier as the spatial `sim_cell_step` path. `0.0` default ⇒ ×1.0 ⇒
        // byte-identical.
        let alox_mul = (1.0 + params.alox_propagation_boost).max(0.0);
        let lp_propagation = lp * effective_unsat * propagation_rate * alox_mul;
        let lp_generation = lp_direct + lp_propagation;

        // === REPAIR ===
        let gpx4_repair = gpx4 * (gsh / (gsh + 1.0)) * params.gpx4_rate * (lp / (lp + 0.5));
        let fsp1_repair = fsp1 * params.fsp1_rate * (lp / (lp + 0.5));
        // DHODH (#338): GPX4-independent repair in parallel to FSP1 (0.0 default).
        let dhodh_repair = params.dhodh_rate * (lp / (lp + 0.5));
        let total_repair = gpx4_repair + fsp1_repair + dhodh_repair;

        lp = (lp + lp_generation - total_repair).max(0.0);

        // === GPX4 DYNAMIC REGULATION ===
        if total_ros > 1.0 {
            gpx4 -= params.gpx4_degradation_by_ros * (total_ros - 1.0);
        }
        let gpx4_target = cell.nrf2 * params.gpx4_nrf2_target_multiplier;
        gpx4 += params.gpx4_nrf2_upregulation * (gpx4_target - gpx4);
        gpx4 = gpx4.max(0.0);

        // Small noise
        lp += norm(rng, 0.0, 0.003);
        lp = lp.max(0.0);

        // Death check. ESCRT-III repair (#465) can rescue across a finite budget,
        // only while the cell is still alive (death_step not yet set). The RNG roll
        // is drawn ONLY when `escrt_repair_rate > 0` (the `escrt_can_attempt`
        // short-circuit), so the default path is byte-identical.
        if lp > params.death_threshold {
            if death_step.is_none()
                && crate::repair::escrt_can_attempt(
                    params.escrt_repair_rate,
                    escrt_budget_used,
                    params.escrt_repair_budget,
                )
                && crate::repair::escrt_rescue(params.escrt_repair_rate, rng.gen::<f64>())
            {
                escrt_budget_used += 1.0;
            } else {
                death_step = Some(step);
            }
        }
    }

    (death_step.is_some(), lp, gsh, gpx4)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cell::{gen_cell, Phenotype};
    use rand::SeedableRng;

    /// #338: DHODH and GCH1/BH4 are GPX4-independent ferroptosis suppressors. Under
    /// RSL3 (GPX4 inhibited), activating either backup must reduce lipid
    /// peroxidation accumulation (the model previously had only the FSP1 backup,
    /// so it overstated RSL3 monotherapy kill). Inhibiting the backup (rate back
    /// to 0) restores kill, which is the GPX4i+DHODHi combination logic. The
    /// default rates (0.0) reproduce the historical behaviour (matrix
    /// byte-identity is guarded by the golden tests); this asserts the DIRECTION.
    #[test]
    fn dhodh_and_gch1_backups_reduce_rsl3_lipid_peroxidation() {
        // Run one Glycolytic cell under RSL3 to a fixed horizon and report the
        // peak LP reached (deterministic given the seed). A higher backup rate
        // must lower the LP the cell reaches.
        fn peak_lp(dhodh: f64, gch1: f64) -> f64 {
            let p = Params {
                dhodh_rate: dhodh,
                gch1_rate: gch1,
                ..Params::default()
            };
            // A modest basal ROS so RSL3 (GPX4 inhibition) can drive LP without
            // a hypoxic collapse.
            let mut gen_rng = StdRng::seed_from_u64(42);
            let cell = gen_cell(Phenotype::Glycolytic, &mut gen_rng);
            let mut rng = StdRng::seed_from_u64(7);
            let mut state = CellState::from_cell(&cell, Treatment::RSL3, &p, &mut rng);
            let mut peak = 0.0_f64;
            for step in 0..180 {
                let died = sim_cell_step(&mut state, &cell, &p, step, 0.0, &mut rng);
                peak = peak.max(state.lp);
                if died {
                    break;
                }
            }
            peak
        }
        let baseline = peak_lp(0.0, 0.0);
        let with_dhodh = peak_lp(2.0, 0.0);
        let with_gch1 = peak_lp(0.0, 5.0);
        assert!(
            with_dhodh < baseline,
            "DHODH backup should lower peak LP under RSL3: dhodh={with_dhodh} vs baseline={baseline}"
        );
        assert!(
            with_gch1 < baseline,
            "GCH1/BH4 backup should lower peak LP under RSL3: gch1={with_gch1} vs baseline={baseline}"
        );
        // Determinism.
        assert_eq!(with_dhodh, peak_lp(2.0, 0.0));
    }

    /// #339: break the MUFA steady-state assumption for ACUTE dosing.
    /// `invivo()` starts cells at M_ss (the established-tumor assumption that
    /// the ~48-72h SCD1/MUFA enrichment is already complete). An acute/naive
    /// start (`mufa_acute_start = Some(0.0)`) begins below M_ss and accumulates
    /// over the run, so the cell is LESS protected early and peroxidizes MORE
    /// under RSL3. `None` reproduces `invivo()` byte-for-byte.
    #[test]
    fn acute_naive_mufa_start_breaks_steady_state_for_acute_dosing() {
        // Peak LP and the MUFA level at step 20 over a 180-step RSL3 run.
        fn peak_lp_and_mufa20(acute_start: Option<f64>) -> (f64, f64) {
            let p = Params {
                mufa_acute_start: acute_start,
                ..Params::invivo()
            };
            let mut gen_rng = StdRng::seed_from_u64(42);
            let cell = gen_cell(Phenotype::OXPHOS, &mut gen_rng);
            let mut rng = StdRng::seed_from_u64(7);
            let mut state = CellState::from_cell(&cell, Treatment::RSL3, &p, &mut rng);
            let mut peak = 0.0_f64;
            let mut mufa20 = state.mufa_protection;
            for step in 0..180 {
                sim_cell_step(&mut state, &cell, &p, step, 0.0, &mut rng);
                if step <= 20 {
                    mufa20 = state.mufa_protection;
                }
                peak = peak.max(state.lp);
            }
            (peak, mufa20)
        }
        let (est_peak, est_mufa20) = peak_lp_and_mufa20(None); // established (M_ss start)
        let (acu_peak, acu_mufa20) = peak_lp_and_mufa20(Some(0.0)); // acute (naive start)

        // The acute cell is less MUFA-protected early...
        assert!(
            acu_mufa20 < est_mufa20,
            "acute MUFA at step 20 should be below established: acute={acu_mufa20} vs est={est_mufa20}"
        );
        // ...and therefore peroxidizes more under RSL3 (and it actually does
        // peroxidize, so the comparison is meaningful, not 0-vs-0).
        assert!(acu_peak > 0.0, "RSL3 should drive lipid peroxidation");
        assert!(
            acu_peak > est_peak,
            "acute (naive-MUFA) cell should peroxidize more under RSL3: acute={acu_peak} vs est={est_peak}"
        );

        // Byte-identity: a None override must reproduce invivo() exactly.
        let run_sim = |p: &Params| {
            let mut gen_rng = StdRng::seed_from_u64(42);
            let cell = gen_cell(Phenotype::OXPHOS, &mut gen_rng);
            let mut rng = StdRng::seed_from_u64(7);
            sim_cell(&cell, Treatment::RSL3, p, &mut rng)
        };
        let base = run_sim(&Params::invivo());
        let none_override = run_sim(&Params {
            mufa_acute_start: None,
            ..Params::invivo()
        });
        assert_eq!(
            base, none_override,
            "None override must reproduce invivo() byte-for-byte"
        );
    }

    /// #339 PR 2: the ether-PUFA pool is extra peroxidizable substrate, so
    /// enabling it INCREASES lipid peroxidation under RSL3 (Zou 2020 / Cui
    /// 2021). The `0.0` limit (FAR1/AGPS-null escape) is the base model.
    #[test]
    fn ether_pufa_pool_increases_peroxidation_under_rsl3() {
        fn peak_lp(ether_fraction: f64) -> f64 {
            let p = Params {
                ether_pufa_fraction: ether_fraction,
                ..Params::default()
            };
            let mut gen_rng = StdRng::seed_from_u64(42);
            let cell = gen_cell(Phenotype::OXPHOS, &mut gen_rng);
            let mut rng = StdRng::seed_from_u64(7);
            let mut state = CellState::from_cell(&cell, Treatment::RSL3, &p, &mut rng);
            let mut peak = 0.0_f64;
            for step in 0..180 {
                sim_cell_step(&mut state, &cell, &p, step, 0.0, &mut rng);
                peak = peak.max(state.lp);
            }
            peak
        }
        let base = peak_lp(0.0); // FAR1/AGPS-null escape (no ether pool)
        let with_ether = peak_lp(0.5); // +50% ether-PUFA substrate
        assert!(base > 0.0, "RSL3 should drive lipid peroxidation");
        assert!(
            with_ether > base,
            "ether-PUFA pool should raise peak LP under RSL3: ether={with_ether} vs base={base}"
        );
        // Determinism.
        assert_eq!(base, peak_lp(0.0));
    }

    /// #339 PR 3: the MBOAT1/2 hormone-regulated MUFA enrichment is a
    /// constitutive, GPX4-independent protective term, so enabling it REDUCES
    /// lipid peroxidation under RSL3 (Liang 2023). `0.0` is byte-identical.
    #[test]
    fn mboat_mufa_boost_reduces_peroxidation_under_rsl3() {
        fn peak_lp(mboat: f64) -> f64 {
            let p = Params {
                mboat_mufa_boost: mboat,
                ..Params::default()
            };
            let mut gen_rng = StdRng::seed_from_u64(42);
            let cell = gen_cell(Phenotype::OXPHOS, &mut gen_rng);
            let mut rng = StdRng::seed_from_u64(7);
            let mut state = CellState::from_cell(&cell, Treatment::RSL3, &p, &mut rng);
            let mut peak = 0.0_f64;
            for step in 0..180 {
                sim_cell_step(&mut state, &cell, &p, step, 0.0, &mut rng);
                peak = peak.max(state.lp);
            }
            peak
        }
        let base = peak_lp(0.0); // no MBOAT enrichment
        let with_mboat = peak_lp(0.3); // +0.3 constant MUFA protection
        assert!(base > 0.0, "RSL3 should drive lipid peroxidation");
        assert!(
            with_mboat < base,
            "MBOAT1/2 MUFA enrichment should lower peak LP under RSL3: mboat={with_mboat} vs base={base}"
        );
        // Determinism.
        assert_eq!(base, peak_lp(0.0));
    }

    /// #340: NCOA4-ferritinophagy releases stored iron into the labile pool over
    /// the run, so enabling it INCREASES peak lipid peroxidation under RSL3 (more
    /// Fenton). `release = 0.0` returns factor 1.0 at every step (byte-identical).
    #[test]
    fn ferritinophagy_iron_release_increases_peroxidation_under_rsl3() {
        fn peak_lp(release: f64) -> f64 {
            let p = Params {
                ferritinophagy_release: release,
                ..Params::default()
            };
            let mut gen_rng = StdRng::seed_from_u64(42);
            let cell = gen_cell(Phenotype::OXPHOS, &mut gen_rng);
            let mut rng = StdRng::seed_from_u64(7);
            let mut state = CellState::from_cell(&cell, Treatment::RSL3, &p, &mut rng);
            let mut peak = 0.0_f64;
            for step in 0..180 {
                sim_cell_step(&mut state, &cell, &p, step, 0.0, &mut rng);
                peak = peak.max(state.lp);
            }
            peak
        }
        let base = peak_lp(0.0); // static iron
        let with_ferritinophagy = peak_lp(1.0); // labile iron up to ~2x late in run
        assert!(base > 0.0, "RSL3 should drive lipid peroxidation");
        assert!(
            with_ferritinophagy > base,
            "ferritinophagy iron release should raise peak LP under RSL3: ferr={with_ferritinophagy} vs base={base}"
        );
        // Determinism.
        assert_eq!(base, peak_lp(0.0));

        // release = 0 ⇒ factor is exactly 1.0 at every step (byte-identical).
        let p0 = Params::default();
        for step in [0u32, 1, 30, 90, 179] {
            assert_eq!(ferritinophagy_iron_factor(step, &p0), 1.0);
        }
        // release > 0 ⇒ starts at exactly 1.0 (step 0) and rises monotonically.
        let pr = Params {
            ferritinophagy_release: 1.0,
            ..Params::default()
        };
        assert_eq!(ferritinophagy_iron_factor(0, &pr), 1.0);
        assert!(ferritinophagy_iron_factor(30, &pr) > ferritinophagy_iron_factor(0, &pr));
        assert!(ferritinophagy_iron_factor(179, &pr) > ferritinophagy_iron_factor(30, &pr));
    }

    #[test]
    fn exo_decay_factor_matches_envelope_formula() {
        // Plateau: 1.0 for every step < 30.
        for step in [0u32, 1, 15, 29] {
            assert_eq!(exo_decay_factor(step), 1.0, "plateau must be 1.0");
        }
        // At the decay onset (step 30): exactly 1.0 (0.5^0).
        assert_eq!(exo_decay_factor(30), 1.0);
        // One half-life (15 steps) later: exactly 0.5.
        assert!((exo_decay_factor(45) - 0.5).abs() < 1e-12);
        // Two half-lives later: 0.25.
        assert!((exo_decay_factor(60) - 0.25).abs() < 1e-12);
        // Matches the raw formula it replaced, for an arbitrary later step.
        let raw = 0.5_f64.powf((123 - 30) as f64 / 15.0);
        assert_eq!(exo_decay_factor(123), raw);
    }

    #[test]
    fn control_does_not_kill_glycolytic() {
        let params = Params::default();
        let mut rng = StdRng::seed_from_u64(0);
        let cell = gen_cell(Phenotype::Glycolytic, &mut rng);
        let mut sim_rng = StdRng::seed_from_u64(1);
        let (dead, lp, _, _) = sim_cell(&cell, Treatment::Control, &params, &mut sim_rng);
        assert!(!dead, "Glycolytic cell should survive Control");
        assert!(
            lp < params.death_threshold,
            "LP should stay below threshold"
        );
    }

    #[test]
    fn sdt_kills_persister() {
        let params = Params::default();
        let mut rng = StdRng::seed_from_u64(0);
        let cell = gen_cell(Phenotype::Persister, &mut rng);
        let mut sim_rng = StdRng::seed_from_u64(1);
        let (dead, _, _, _) = sim_cell(&cell, Treatment::SDT, &params, &mut sim_rng);
        assert!(dead, "Persister cell should die under SDT");
    }

    #[test]
    fn rsl3_inhibits_gpx4() {
        let params = Params::default();
        let mut rng = StdRng::seed_from_u64(0);
        let cell = gen_cell(Phenotype::Glycolytic, &mut rng);
        let mut sim_rng = StdRng::seed_from_u64(1);
        let state = CellState::from_cell(&cell, Treatment::RSL3, &params, &mut sim_rng);
        let expected = cell.gpx4 * (1.0 - params.rsl3_gpx4_inhib);
        assert!(
            (state.gpx4 - expected).abs() < 1e-10,
            "RSL3 should reduce GPX4 by {}%",
            params.rsl3_gpx4_inhib * 100.0
        );
    }

    #[test]
    fn mufa_protection_reduces_death_rate() {
        // In-vivo MUFA should protect persisters from RSL3 relative to 2D
        let params_2d = Params::default();
        let params_vivo = Params::invivo();
        let n = 1000;
        let mut deaths_2d = 0;
        let mut deaths_vivo = 0;
        for i in 0..n {
            let mut rng = StdRng::seed_from_u64(i * 2);
            let cell = gen_cell(Phenotype::Persister, &mut rng);
            let mut sr = StdRng::seed_from_u64(i * 2 + 1);
            if sim_cell(&cell, Treatment::RSL3, &params_2d, &mut sr).0 {
                deaths_2d += 1;
            }
            let mut sr = StdRng::seed_from_u64(i * 2 + 1);
            if sim_cell(&cell, Treatment::RSL3, &params_vivo, &mut sr).0 {
                deaths_vivo += 1;
            }
        }
        assert!(
            deaths_vivo < deaths_2d,
            "In-vivo MUFA should reduce RSL3 deaths: 2D={deaths_2d}, vivo={deaths_vivo}"
        );
    }

    #[test]
    fn single_step_does_not_kill_healthy_cell() {
        let params = Params::default();
        let mut rng = StdRng::seed_from_u64(0);
        let cell = gen_cell(Phenotype::Glycolytic, &mut rng);
        let mut sim_rng = StdRng::seed_from_u64(1);
        let mut state = CellState::from_cell(&cell, Treatment::Control, &params, &mut sim_rng);
        let dead = sim_cell_step(&mut state, &cell, &params, 0, 0.0, &mut sim_rng);
        assert!(!dead, "One step should not kill a healthy glycolytic cell");
        assert!(state.lp < 1.0, "LP should be near zero after one step");
    }

    /// #270: a per-cell MUFA cap yields a per-cell durable steady state.
    /// Two MUFA levels start at the SAME low value and relax under
    /// `Params::spheroid()`: one toward a low per-cell cap (the spheroid core),
    /// one toward the global cap (`None` fallback = the OLD uniform behavior).
    /// The capped one stays low; the uncapped one climbs to the uniform M_ss —
    /// so the position-dependent MUFA persists instead of converging.
    #[test]
    fn per_cell_mufa_cap_gives_durable_per_cell_steady_state() {
        let params = Params::spheroid();
        let core_cap = 0.05; // a core cell's per-cell cap
        let global_cap = params.scd_mufa_max; // 0.25 — what `None` falls back to
        let (mut core, mut uncapped) = (0.05_f64, 0.05_f64); // same low start
        for _ in 0..300 {
            core = update_mufa_protection(core, core_cap, params.scd_mufa_rate, &params);
            uncapped = update_mufa_protection(uncapped, global_cap, params.scd_mufa_rate, &params);
        }
        // Core saturates near M_ss(0.05) ≈ 0.048; uncapped climbs to M_ss(0.25) ≈ 0.20.
        assert!(
            core < 0.08,
            "core MUFA stays low at its per-cell cap: {core}"
        );
        assert!(
            uncapped > 0.15,
            "uncapped (global) MUFA rises to the uniform M_ss: {uncapped}"
        );
        assert!(
            uncapped > 2.0 * core,
            "per-cell cap keeps a durable rim-vs-core spread: uncapped={uncapped}, core={core}"
        );
        // Each stays within its own cap (clamp invariant).
        assert!(core <= core_cap + 1e-9 && uncapped <= global_cap + 1e-9);
    }

    /// #270 wiring: `sim_cell_step` itself must read `cell.mufa_cap` (not just
    /// the `update_mufa_protection` helper). Routes two otherwise-identical
    /// cells through the full step — one with a low per-cell cap, one uncapped
    /// (`None` ⇒ global) — from the same low MUFA start under `Params::spheroid()`
    /// and Control (the cell stays alive so MUFA relaxes cleanly). The capped
    /// cell's late MUFA must stay below the uncapped one; this fails if the
    /// `sim_cell_step` call site ever stops threading the per-cell cap.
    #[test]
    fn sim_cell_step_reads_per_cell_mufa_cap() {
        let params = Params::spheroid();
        let mut gen_rng = StdRng::seed_from_u64(11);
        let base = gen_cell(Phenotype::Glycolytic, &mut gen_rng);

        let run = |mufa_cap: Option<f64>| -> f64 {
            let mut cell = base.clone();
            cell.mufa_cap = mufa_cap;
            let mut init_rng = StdRng::seed_from_u64(5);
            let mut state = CellState::from_cell(
                &cell,
                crate::cell::Treatment::Control,
                &params,
                &mut init_rng,
            );
            state.mufa_protection = 0.05; // same low start for both
            let mut step_rng = StdRng::seed_from_u64(77);
            for step in 0..150 {
                sim_cell_step(&mut state, &cell, &params, step, 0.0, &mut step_rng);
            }
            state.mufa_protection
        };

        let capped = run(Some(0.05)); // a core cell's low cap
        let uncapped = run(None); // global cap (0.25) — the old uniform target
        assert!(
            capped < uncapped,
            "sim_cell_step must honor the per-cell cap: capped={capped}, uncapped={uncapped}"
        );
        assert!(capped < 0.08, "capped cell stays MUFA-poor late: {capped}");
        assert!(
            uncapped > 0.15,
            "uncapped cell climbs to the global M_ss: {uncapped}"
        );
    }

    /// #363 wiring: `sim_cell_step` must read the per-cell `cell.mufa_rate`, so
    /// phenotype-specific SCD1/MUFA accumulation RATES actually diverge. Three
    /// otherwise-identical cells from the same naive MUFA start under
    /// `Params::spheroid()` + Control (alive, so MUFA accumulates cleanly), with a
    /// deliberately NON-binding per-cell cap so the cap does not mask the rate
    /// effect: (a) `None` ⇒ the global `scd_mufa_rate`, (b) `Some(global)` which
    /// must be byte-identical to (a), and (c) a FAST rate which must build more
    /// MUFA protection by the same step. This is the #363 acceptance ("two
    /// phenotypes diverge in MUFA build-up under the same dosing") at the engine
    /// level, and it guards the `None`-is-byte-identical default.
    #[test]
    fn sim_cell_step_reads_per_cell_mufa_rate() {
        let params = Params::spheroid();
        let mut gen_rng = StdRng::seed_from_u64(11);
        let base = gen_cell(Phenotype::Glycolytic, &mut gen_rng);

        let run = |mufa_rate: Option<f64>| -> f64 {
            let mut cell = base.clone();
            cell.mufa_rate = mufa_rate;
            cell.mufa_cap = Some(1.0); // non-binding cap: isolate the RATE effect
            let mut init_rng = StdRng::seed_from_u64(5);
            let mut state = CellState::from_cell(
                &cell,
                crate::cell::Treatment::Control,
                &params,
                &mut init_rng,
            );
            state.mufa_protection = 0.0; // same naive start for all runs
            let mut step_rng = StdRng::seed_from_u64(77);
            for step in 0..40 {
                sim_cell_step(&mut state, &cell, &params, step, 0.0, &mut step_rng);
            }
            state.mufa_protection
        };

        let global = run(None); // global scd_mufa_rate (the default path)
        let global_explicit = run(Some(params.scd_mufa_rate)); // == global, exactly
        let fast = run(Some(params.scd_mufa_rate * 3.0)); // a faster phenotype

        // `None` falls back to the global rate EXACTLY (the byte-identical default).
        assert_eq!(
            global, global_explicit,
            "mufa_rate=None must reproduce the global scd_mufa_rate bit-for-bit"
        );
        // A faster per-cell rate builds more MUFA protection by the same step —
        // the two phenotypes diverge in MUFA build-up under identical dosing.
        assert!(
            fast > global,
            "a faster per-cell mufa_rate must accumulate more MUFA: fast={fast}, global={global}"
        );
    }

    /// #446: `sim_cell_step` must respect both ALOX boosts. Isolates the
    /// propagation/substrate multipliers from the bistable switch by taking a
    /// SINGLE step from a fixed mid-LP state with the defenses zeroed (so repair
    /// and antioxidant quench are ~0 and the autocatalytic propagation term
    /// dominates). A positive `alox_propagation_boost` (faster enzymatic
    /// propagation) and a positive `mcfa_pufa_boost` (more oxidizable PUFA) each
    /// raise the post-step LP; `0.0` (the defaults) is the byte-identical
    /// baseline. The RNG seed is identical across runs, so the only difference is
    /// the boost.
    #[test]
    fn sim_cell_step_respects_alox_and_mcfa_boosts() {
        let mut gen_rng = StdRng::seed_from_u64(13);
        let mut cell = gen_cell(Phenotype::Glycolytic, &mut gen_rng);
        cell.lipid_unsat = 2.0;

        let step_once = |alox_boost: f64, mcfa_boost: f64| -> f64 {
            let mut params = Params::default();
            params.alox_propagation_boost = alox_boost;
            params.mcfa_pufa_boost = mcfa_boost;
            let mut init_rng = StdRng::seed_from_u64(7);
            let mut state = CellState::from_cell(
                &cell,
                crate::cell::Treatment::Control,
                &params,
                &mut init_rng,
            );
            // Zero defenses + seed a mid LP so the propagation term dominates.
            state.gsh = 0.0;
            state.gpx4 = 0.0;
            state.fsp1 = 0.0;
            state.lp = 3.0;
            let mut step_rng = StdRng::seed_from_u64(0); // identical noise draw
            sim_cell_step(&mut state, &cell, &params, 0, 0.0, &mut step_rng);
            state.lp
        };

        let baseline = step_once(0.0, 0.0);
        let alox_high = step_once(0.5, 0.0);
        let mcfa_high = step_once(0.0, 0.5);
        assert!(
            alox_high > baseline,
            "ALOX propagation boost must raise post-step LP: alox_high={alox_high}, baseline={baseline}"
        );
        assert!(
            mcfa_high > baseline,
            "MCFA PUFA boost must raise post-step LP: mcfa_high={mcfa_high}, baseline={baseline}"
        );
    }

    #[test]
    fn escrt_repair_brakes_death_execution() {
        // #465: an over-threshold cell dies by default, but ESCRT membrane repair
        // (rate > 0, budget remaining) reseals it and it survives the step; with no
        // budget the rescue cannot fire. Deterministic single-step A/B.
        let mut gen_rng = StdRng::seed_from_u64(21);
        let cell = gen_cell(Phenotype::Glycolytic, &mut gen_rng);

        let step_once = |escrt_rate: f64, escrt_budget: f64| -> (bool, f64) {
            let mut params = Params::default();
            params.escrt_repair_rate = escrt_rate;
            params.escrt_repair_budget = escrt_budget;
            let mut init_rng = StdRng::seed_from_u64(3);
            let mut state = CellState::from_cell(
                &cell,
                crate::cell::Treatment::Control,
                &params,
                &mut init_rng,
            );
            // Defenses off + LP already over the death threshold ⇒ this step's death
            // check fires.
            state.gsh = 0.0;
            state.gpx4 = 0.0;
            state.fsp1 = 0.0;
            state.lp = params.death_threshold + 0.5;
            let mut step_rng = StdRng::seed_from_u64(0);
            let died = sim_cell_step(&mut state, &cell, &params, 0, 0.0, &mut step_rng);
            (died, state.escrt_budget_used)
        };

        // Default (off): the over-threshold cell dies, no budget touched.
        let (died_off, used_off) = step_once(0.0, 0.0);
        assert!(died_off, "without ESCRT the over-threshold cell must die");
        assert_eq!(
            used_off, 0.0,
            "default path must not touch the repair budget"
        );

        // On (rate 1.0, budget available): resealed this step, one repair consumed.
        let (died_on, used_on) = step_once(1.0, 5.0);
        assert!(
            !died_on,
            "with ESCRT (rate 1, budget) the cell must be rescued"
        );
        assert_eq!(used_on, 1.0, "exactly one repair event consumed");

        // Enabled but zero budget ⇒ no rescue ⇒ death proceeds.
        let (died_spent, _) = step_once(1.0, 0.0);
        assert!(died_spent, "ESCRT with zero budget cannot rescue");
    }

    #[test]
    fn sim_cell_step_respects_por_h2o2_rate() {
        // #466: the POR/CYB5R1 enzymatic H2O2 source raises the oxidant pool, so a
        // positive `por_h2o2_rate` raises post-step LP. `0.0` (default) is the
        // byte-identical baseline. Identical seeds, so only the rate differs.
        let mut gen_rng = StdRng::seed_from_u64(31);
        let cell = gen_cell(Phenotype::Glycolytic, &mut gen_rng);

        let step_once = |por: f64| -> f64 {
            let mut params = Params::default();
            params.por_h2o2_rate = por;
            let mut init_rng = StdRng::seed_from_u64(9);
            let mut state = CellState::from_cell(
                &cell,
                crate::cell::Treatment::Control,
                &params,
                &mut init_rng,
            );
            state.gsh = 0.0; // remove scavenging so the extra oxidant reaches LP
            state.gpx4 = 0.0;
            state.fsp1 = 0.0;
            state.lp = 1.0;
            let mut step_rng = StdRng::seed_from_u64(0);
            sim_cell_step(&mut state, &cell, &params, 0, 0.0, &mut step_rng);
            state.lp
        };

        let baseline = step_once(0.0);
        let por_high = step_once(2.0);
        assert!(
            por_high > baseline,
            "POR H2O2 source must raise post-step LP: por_high={por_high}, baseline={baseline}"
        );
    }
}
