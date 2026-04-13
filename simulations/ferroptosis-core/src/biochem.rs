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
        }
    }

    /// Initialize with a custom exogenous ROS peak (for spatial model where
    /// ROS dose depends on depth/position).
    pub fn from_cell_with_ros(cell: &Cell, tx: Treatment, params: &Params, exo_ros_peak: f64) -> Self {
        let mut gpx4 = cell.gpx4;
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
#[inline]
fn update_mufa_protection(current: f64, params: &Params) -> f64 {
    let growth = params.scd_mufa_rate * (1.0 - current / (params.scd_mufa_max + 1e-9));
    let decay = params.scd_mufa_decay * current;
    (current + growth - decay).clamp(0.0, params.scd_mufa_max.max(0.0))
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
        return false;
    }

    // === ROS SOURCES ===
    let effective_iron = cell.iron + extra_iron;
    let fenton = effective_iron * params.fenton_rate * norm(rng, 1.0, 0.08).max(0.0);
    let exo = if step < 30 {
        state.exo_ros_peak * norm(rng, 1.0, 0.1).max(0.0)
    } else {
        state.exo_ros_peak * 0.5_f64.powf((step - 30) as f64 / 15.0)
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
    state.mufa_protection = update_mufa_protection(state.mufa_protection, params);

    let effective_unsat = (cell.lipid_unsat * (1.0 - state.mufa_protection)).max(0.05);
    let lp_direct = unscav * effective_unsat * params.lp_rate;
    // AUTOCATALYTIC PROPAGATION — GSH-gated bistable switch
    let antioxidant_quench = state.gpx4 * (state.gsh / (state.gsh + 0.5)) + state.fsp1;
    let propagation_rate = params.lp_propagation / (1.0 + antioxidant_quench * 5.0);
    let lp_propagation = state.lp * effective_unsat * propagation_rate;
    let lp_generation = lp_direct + lp_propagation;

    // === REPAIR ===
    let gpx4_repair = state.gpx4 * (state.gsh / (state.gsh + 1.0)) * params.gpx4_rate
        * (state.lp / (state.lp + 0.5));
    let fsp1_repair = state.fsp1 * params.fsp1_rate * (state.lp / (state.lp + 0.5));
    let total_repair = gpx4_repair + fsp1_repair;

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
        // === ROS SOURCES ===
        let fenton = cell.iron * params.fenton_rate * norm(rng, 1.0, 0.08).max(0.0);
        let exo = if step < 30 {
            exo_ros_peak * norm(rng, 1.0, 0.1).max(0.0)
        } else {
            exo_ros_peak * 0.5_f64.powf((step - 30) as f64 / 15.0)
        };
        let total_ros = cell.basal_ros + exo + fenton;

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
        mufa_protection = update_mufa_protection(mufa_protection, params);
        let effective_unsat = (cell.lipid_unsat * (1.0 - mufa_protection)).max(0.05);
        let lp_direct = unscav * effective_unsat * params.lp_rate;
        let antioxidant_quench = gpx4 * (gsh / (gsh + 0.5)) + fsp1;
        let propagation_rate = params.lp_propagation / (1.0 + antioxidant_quench * 5.0);
        let lp_propagation = lp * effective_unsat * propagation_rate;
        let lp_generation = lp_direct + lp_propagation;

        // === REPAIR ===
        let gpx4_repair = gpx4 * (gsh / (gsh + 1.0)) * params.gpx4_rate * (lp / (lp + 0.5));
        let fsp1_repair = fsp1 * params.fsp1_rate * (lp / (lp + 0.5));
        let total_repair = gpx4_repair + fsp1_repair;

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

        if lp > params.death_threshold {
            break;
        }
    }

    (lp > params.death_threshold, lp, gsh, gpx4)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cell::{gen_cell, Phenotype};
    use rand::SeedableRng;

    #[test]
    fn control_does_not_kill_glycolytic() {
        let params = Params::default();
        let mut rng = StdRng::seed_from_u64(0);
        let cell = gen_cell(Phenotype::Glycolytic, &mut rng);
        let mut sim_rng = StdRng::seed_from_u64(1);
        let (dead, lp, _, _) = sim_cell(&cell, Treatment::Control, &params, &mut sim_rng);
        assert!(!dead, "Glycolytic cell should survive Control");
        assert!(lp < params.death_threshold, "LP should stay below threshold");
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
        assert!((state.gpx4 - expected).abs() < 1e-10, "RSL3 should reduce GPX4 by {}%", params.rsl3_gpx4_inhib * 100.0);
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
            if sim_cell(&cell, Treatment::RSL3, &params_2d, &mut sr).0 { deaths_2d += 1; }
            let mut sr = StdRng::seed_from_u64(i * 2 + 1);
            if sim_cell(&cell, Treatment::RSL3, &params_vivo, &mut sr).0 { deaths_vivo += 1; }
        }
        assert!(deaths_vivo < deaths_2d, "In-vivo MUFA should reduce RSL3 deaths: 2D={deaths_2d}, vivo={deaths_vivo}");
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
}
