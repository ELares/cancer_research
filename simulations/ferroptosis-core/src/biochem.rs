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
            mufa_protection: params.initial_mufa_protection,
            lp: 0.0,
            dead: false,
            death_step: None,
            exo_ros_peak,
            persister_fraction: 0.0,
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
            mufa_protection: params.initial_mufa_protection,
            lp: 0.0,
            dead: false,
            death_step: None,
            exo_ros_peak,
            persister_fraction: 0.0,
        }
    }
}

/// Accumulate MUFA-style lipid-remodeling protection against peroxidation.
///
/// SCD1 (the enzyme converting SFA→MUFA) is regulated by SREBP1/mTORC1,
/// NOT by NRF2. In 3D culture and in vivo, SCD1 is constitutively active
/// and enriches membranes with MUFAs that displace PUFAs, reducing
/// ferroptosis susceptibility (Dixon/Park, Cancer Research 2025;
/// Magtanong et al., Mol Cell 2019; Tesfay et al., Cancer Res 2019).
///
/// The rate is context-dependent: zero in 2D culture (Params::default),
/// non-zero in in-vivo-like conditions (Params::invivo).
/// `mufa_max` is the MUFA carrying capacity: the per-cell `cell.mufa_cap` when
/// set (the spheroid's radial cap, #270), else the global `params.scd_mufa_max`.
/// The logistic growth saturates at — and the value is clamped to — `mufa_max`,
/// so a per-cell cap yields a per-cell steady state (durable position-dependent
/// MUFA). `cell.mufa_cap = None` ⇒ `params.scd_mufa_max` ⇒ byte-identical.
#[inline]
fn update_mufa_protection(current: f64, mufa_max: f64, params: &Params) -> f64 {
    let growth = params.scd_mufa_rate * (1.0 - current / (mufa_max + 1e-9));
    let decay = params.scd_mufa_decay * current;
    (current + growth - decay).clamp(0.0, mufa_max.max(0.0))
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
                let effective_iron = cell.iron + extra_iron;
                let fenton = effective_iron * params.fenton_rate * norm(rng, 1.0, 0.08).max(0.0);
                let exo = if step < 30 {
                    state.exo_ros_peak * norm(rng, 1.0, 0.1).max(0.0)
                } else {
                    state.exo_ros_peak * exo_decay_factor(step)
                };
                let total_ros = cell.basal_ros + exo + fenton;
                let effective_unsat = cell.lipid_unsat; // no MUFA protection
                let lp_direct = total_ros * effective_unsat * params.lp_rate;
                let lp_propagation = state.lp * effective_unsat * params.lp_propagation;
                state.lp += lp_direct + lp_propagation;
            }
        }
        return false;
    }

    // === ROS SOURCES ===
    let effective_iron = cell.iron + extra_iron;
    let fenton = effective_iron * params.fenton_rate * norm(rng, 1.0, 0.08).max(0.0);
    let exo = if step < 30 {
        state.exo_ros_peak * norm(rng, 1.0, 0.1).max(0.0)
    } else {
        state.exo_ros_peak * exo_decay_factor(step)
    };
    let total_ros = cell.basal_ros + exo + fenton;

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
        params,
    );

    let effective_unsat = (cell.lipid_unsat * (1.0 - state.mufa_protection)).max(0.05);
    let lp_direct = unscav * effective_unsat * params.lp_rate;
    // AUTOCATALYTIC PROPAGATION — GSH-gated bistable switch.
    // GCH1/BH4 (#338) adds GPX4-independent radical-trapping quench capacity
    // (`gch1_rate`, 0.0 by default ⇒ byte-identical).
    let antioxidant_quench =
        state.gpx4 * (state.gsh / (state.gsh + 0.5)) + state.fsp1 + params.gch1_rate;
    let propagation_rate = params.lp_propagation / (1.0 + antioxidant_quench * 5.0);
    let lp_propagation = state.lp * effective_unsat * propagation_rate;
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

    // Death check
    if state.lp > params.death_threshold {
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
    let mut mufa_protection = params.initial_mufa_protection;
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
        let fenton = cell.iron * params.fenton_rate * norm(rng, 1.0, 0.08).max(0.0);
        let exo = if step < 30 {
            exo_ros_peak * norm(rng, 1.0, 0.1).max(0.0)
        } else {
            exo_ros_peak * exo_decay_factor(step)
        };
        let total_ros = cell.basal_ros + exo + fenton;

        if death_step.is_some() {
            // Post-death: LP-only accumulation (no GSH, no repair, no GPX4).
            // Break check BEFORE accumulation to match sim_cell_step, which
            // only accumulates when step < death_step + post_death_steps.
            if step >= death_step.unwrap() + params.post_death_steps {
                break;
            }
            let effective_unsat = cell.lipid_unsat;
            let lp_direct = total_ros * effective_unsat * params.lp_rate;
            let lp_prop = lp * effective_unsat * params.lp_propagation;
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
            params,
        );
        let effective_unsat = (cell.lipid_unsat * (1.0 - mufa_protection)).max(0.05);
        let lp_direct = unscav * effective_unsat * params.lp_rate;
        // GCH1/BH4 (#338) adds GPX4-independent quench (`gch1_rate`, 0.0 default).
        let antioxidant_quench = gpx4 * (gsh / (gsh + 0.5)) + fsp1 + params.gch1_rate;
        let propagation_rate = params.lp_propagation / (1.0 + antioxidant_quench * 5.0);
        let lp_propagation = lp * effective_unsat * propagation_rate;
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

        // Death check
        if lp > params.death_threshold {
            death_step = Some(step);
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
            core = update_mufa_protection(core, core_cap, &params);
            uncapped = update_mufa_protection(uncapped, global_cap, &params);
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
}
