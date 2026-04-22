//! Tumor pharmacokinetic compartment model.
//!
//! Models time-dependent drug delivery from plasma to tumor interstitium
//! via a two-compartment ODE (vascular ↔ interstitial). Complements the
//! spatial `drug_transport` module (static Krogh cylinder C(r)) by adding
//! the temporal dimension C(t).
//!
//! The model captures five key tumor-specific barriers:
//! 1. Blood flow (Q): how much drug arrives per minute
//! 2. Vascular permeability (PS): how much crosses the vessel wall
//! 3. Partition coefficient (K_p): drug affinity for tumor tissue
//! 4. Interstitial fluid pressure (IFP): opposes drug inflow in tumors
//! 5. Cellular uptake and metabolism: drug consumed by cells
//!
//! References: Jain RK, Cancer Res 1988; Baxter & Jain, Microvascular Res 1989;
//! Boucher et al., Cancer Res 1990.
//!
//! All tumor-specific parameters are ESTIMATED. No textbook coverage for
//! tumor PK exists. Chemistry2e covers first-order kinetics (Ch.12);
//! Biology2e covers Michaelis-Menten kinetics conceptually.

use serde::{Deserialize, Serialize};

use crate::biochem::{sim_cell_step, CellState};
use crate::cell::{Cell, Treatment};
use crate::params::Params;
use rand::SeedableRng;

/// Two-compartment tumor pharmacokinetic parameters.
///
/// Models the vascular space (drug arrives from plasma) and interstitial
/// space (drug reaches cancer cells) as coupled compartments.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct TumorPKParams {
    /// Blood flow rate Q (mL/min/g tissue). Determines drug delivery rate.
    pub blood_flow_q: f64,
    /// Vascular volume fraction V_v (dimensionless, 0-1).
    pub vascular_fraction: f64,
    /// Interstitial volume fraction V_i (dimensionless, 0-1).
    pub interstitial_fraction: f64,
    /// Permeability-surface area product PS (mL/min/g). Governs vessel
    /// wall transport. Low for BBB tumors (GBM).
    pub ps_product: f64,
    /// Tissue:plasma partition coefficient K_p (dimensionless).
    pub partition_coeff: f64,
    /// Interstitial fluid pressure P_i (mmHg). High in desmoplastic tumors
    /// (pancreatic ~100 mmHg). Opposes drug inflow when P_i > P_v.
    pub ifp_mmhg: f64,
    /// Microvascular pressure P_v (mmHg). Typically ~25 mmHg.
    pub mvp_mmhg: f64,
    /// Hydraulic conductivity L_p (mL/min/mmHg/g).
    pub hydraulic_conductivity: f64,
    /// Osmotic reflection coefficient sigma (0-1).
    pub reflection_coeff: f64,
    /// Bulk cellular uptake rate (1/min). Fraction of interstitial drug
    /// consumed by cells per minute. Uses bulk rate (not per-cell × density)
    /// because the per-cell parameterization from Jain/Baxter produces
    /// unreasonably high clearance (k_up × ρ = 0.01 × 1e8 = 1e6/min).
    pub k_uptake_bulk: f64,
    /// Michaelis-Menten K_m for cellular uptake (concentration units).
    pub km_uptake: f64,
    /// Vascular metabolism rate (1/min).
    pub k_met_v: f64,
    /// Interstitial metabolism rate (1/min).
    pub k_met_i: f64,
    /// Human-readable tumor type name.
    pub name: &'static str,
}

/// Breast tumor: well-vascularized, moderate IFP.
/// Ref: Jain 1988 (blood flow), Boucher 1990 (IFP 10-30 mmHg).
pub fn breast_tumor() -> TumorPKParams {
    TumorPKParams {
        blood_flow_q: 0.25,
        vascular_fraction: 0.07,
        interstitial_fraction: 0.30,
        ps_product: 0.10,
        partition_coeff: 0.50,
        ifp_mmhg: 20.0,
        mvp_mmhg: 25.0,
        hydraulic_conductivity: 1e-7,
        reflection_coeff: 0.9,
        k_uptake_bulk: 0.02,
        km_uptake: 0.5,
        k_met_v: 0.001,
        k_met_i: 0.001,
        name: "Breast (well-vascularized)",
    }
}

/// Pancreatic tumor: desmoplastic stroma, very high IFP, poor perfusion.
/// Ref: Provenzano et al., Cancer Cell 2012 (IFP 75-130 mmHg).
pub fn pancreatic_tumor() -> TumorPKParams {
    TumorPKParams {
        blood_flow_q: 0.05,
        vascular_fraction: 0.02,
        interstitial_fraction: 0.15,
        ps_product: 0.03,
        partition_coeff: 0.30,
        ifp_mmhg: 100.0,
        mvp_mmhg: 25.0,
        hydraulic_conductivity: 1e-7,
        reflection_coeff: 0.9,
        k_uptake_bulk: 0.02,
        km_uptake: 0.5,
        k_met_v: 0.001,
        k_met_i: 0.001,
        name: "Pancreatic (desmoplastic)",
    }
}

/// Glioblastoma: behind blood-brain barrier, very low PS and K_p.
/// Ref: Sarkaria et al., Neuro-Oncology 2018 (BBB drug restriction).
pub fn glioblastoma_tumor() -> TumorPKParams {
    TumorPKParams {
        blood_flow_q: 0.10,
        vascular_fraction: 0.03,
        interstitial_fraction: 0.20,
        ps_product: 0.02,
        partition_coeff: 0.15,
        ifp_mmhg: 10.0,
        mvp_mmhg: 25.0,
        hydraulic_conductivity: 1e-7,
        reflection_coeff: 0.9,
        k_uptake_bulk: 0.02,
        km_uptake: 0.5,
        k_met_v: 0.001,
        k_met_i: 0.001,
        name: "GBM (blood-brain barrier)",
    }
}

// ============================================================
// Plasma concentration models
// ============================================================

/// Analytical plasma concentration models.
/// Phase 1 uses analytical solutions; Phase 2 will add PK-Sim CSV import.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum PlasmaModel {
    /// IV bolus: C(t) = C0 × exp(-k_el × t)
    IvBolus {
        /// Initial plasma concentration (normalized, typically 1.0).
        c0: f64,
        /// Elimination rate constant (1/min). k_el = ln(2) / t_half.
        k_el: f64,
    },
    /// Constant concentration (for 2D culture reference validation).
    Constant {
        concentration: f64,
    },
}

impl PlasmaModel {
    /// Plasma concentration at time t (minutes).
    pub fn concentration_at(&self, t_min: f64) -> f64 {
        match self {
            PlasmaModel::IvBolus { c0, k_el } => c0 * (-k_el * t_min).exp(),
            PlasmaModel::Constant { concentration } => *concentration,
        }
    }
}

/// RSL3-like small molecule IV bolus. t_half ≈ 30 min (estimated for
/// a chloroacetamide GPX4 inhibitor; actual RSL3 PK not well-characterized).
pub fn rsl3_iv_bolus() -> PlasmaModel {
    PlasmaModel::IvBolus {
        c0: 1.0,
        k_el: (2.0_f64).ln() / 30.0, // t_half = 30 min
    }
}

/// Constant max-exposure reference: concentration = 1.0 for all time.
/// NOT equivalent to the standard sim_cell RSL3 baseline (~42%): this model
/// applies continuous inhibition via GPX4 clamp (preventing NRF2 recovery),
/// while sim_cell applies one-time GPX4 reduction at init. Use as a
/// theoretical maximum for computing relative protection factors.
pub fn constant_reference() -> PlasmaModel {
    PlasmaModel::Constant { concentration: 1.0 }
}

// ============================================================
// ODE solver
// ============================================================

/// Result of tumor PK ODE integration.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct TumorPKResult {
    /// Time points in minutes (one per simulation step).
    pub time_min: Vec<f64>,
    /// Plasma concentration at each minute.
    pub c_plasma: Vec<f64>,
    /// Vascular compartment concentration.
    pub c_vascular: Vec<f64>,
    /// Interstitial compartment concentration (drug at cancer cells).
    pub c_interstitial: Vec<f64>,
}

/// Solve the two-compartment tumor PK ODE using forward Euler.
///
/// Returns concentration time-courses at 1-minute resolution (one value per
/// simulation step). Internal sub-stepping at `substeps_per_min` ensures
/// numerical stability for fast vascular equilibration dynamics.
///
/// The ODE:
/// ```text
/// dC_v/dt = Q/V_v × (C_p - C_v) - PS/V_v × (C_v - C_i/K_p) - k_met_v × C_v
/// dC_i/dt = PS/V_i × (C_v - C_i/K_p) - k_uptake × C_i/(C_i + K_m) - k_met_i × C_i
///           - σ × L_p × (P_v - P_i)/V_i × C_i
/// ```
pub fn solve_tumor_pk(
    plasma: &PlasmaModel,
    tumor: &TumorPKParams,
    n_steps: usize,
    substeps_per_min: usize,
) -> TumorPKResult {
    let dt = 1.0 / substeps_per_min as f64;
    let mut c_v = 0.0_f64;
    let mut c_i = 0.0_f64;

    let mut time_min = Vec::with_capacity(n_steps);
    let mut c_plasma_out = Vec::with_capacity(n_steps);
    let mut c_vascular_out = Vec::with_capacity(n_steps);
    let mut c_interstitial_out = Vec::with_capacity(n_steps);

    for minute in 0..n_steps {
        for sub in 0..substeps_per_min {
            let t = minute as f64 + sub as f64 * dt;
            let c_p = plasma.concentration_at(t);

            // Vascular ODE
            let dc_v = tumor.blood_flow_q / tumor.vascular_fraction * (c_p - c_v)
                - tumor.ps_product / tumor.vascular_fraction * (c_v - c_i / tumor.partition_coeff)
                - tumor.k_met_v * c_v;

            // Interstitial ODE
            let convection = tumor.reflection_coeff
                * tumor.hydraulic_conductivity
                * (tumor.mvp_mmhg - tumor.ifp_mmhg)
                / tumor.interstitial_fraction
                * c_i;
            let uptake = if c_i > 0.0 {
                tumor.k_uptake_bulk * c_i / (c_i + tumor.km_uptake)
            } else {
                0.0
            };
            let dc_i = tumor.ps_product / tumor.interstitial_fraction
                * (c_v - c_i / tumor.partition_coeff)
                - uptake
                - tumor.k_met_i * c_i
                - convection;

            c_v = (c_v + dc_v * dt).max(0.0);
            c_i = (c_i + dc_i * dt).max(0.0);
        }

        let c_p = plasma.concentration_at(minute as f64);
        time_min.push(minute as f64);
        c_plasma_out.push(c_p);
        c_vascular_out.push(c_v);
        c_interstitial_out.push(c_i);
    }

    TumorPKResult {
        time_min,
        c_plasma: c_plasma_out,
        c_vascular: c_vascular_out,
        c_interstitial: c_interstitial_out,
    }
}

// ============================================================
// Integration with ferroptosis-core
// ============================================================

/// Result of a single cell simulation with time-varying PK.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PKCellResult {
    pub dead: bool,
    pub death_step: Option<u32>,
    pub final_lp: f64,
    pub final_gsh: f64,
    pub final_gpx4: f64,
}

/// Simulate a single cell with time-varying drug concentration.
///
/// At each timestep, GPX4 inhibition tracks the interstitial drug
/// concentration: `effective_gpx4 = intrinsic_gpx4 × (1 - base_inhib × conc[step])`.
/// This models competitive inhibition where the drug and NRF2 upregulation
/// fight for GPX4 control at each moment.
///
/// The cell is initialized as Treatment::Control (no initial GPX4 reduction).
/// Drug effect is applied dynamically through the concentration schedule.
pub fn sim_cell_with_pk(
    cell: &Cell,
    params: &Params,
    conc_schedule: &[f64],
    base_gpx4_inhib: f64,
    seed: u64,
) -> PKCellResult {
    let n_steps = conc_schedule.len().min(180);
    let mut state = CellState::from_cell_with_ros(cell, Treatment::Control, params, 0.0);
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);

    for step in 0..n_steps as u32 {
        // Time-varying drug inhibition: CLAMP GPX4 BEFORE the step so the
        // repair calculation uses the drug-limited level. The drug is a
        // reversible competitive inhibitor that prevents GPX4 from exceeding
        // a concentration-dependent ceiling. NRF2 can upregulate GPX4 between
        // steps, but the drug caps it before the next repair calculation.
        // When drug washes out (conc→0), the cap lifts and GPX4 recovers.
        if !state.dead {
            let conc = conc_schedule[step as usize];
            let max_gpx4 = cell.gpx4 * (1.0 - base_gpx4_inhib * conc.clamp(0.0, 1.0));
            state.gpx4 = state.gpx4.min(max_gpx4);
        }

        let _died = sim_cell_step(&mut state, cell, params, step, 0.0, &mut rng);

        if state.dead {
            if let Some(ds) = state.death_step {
                if step >= ds + params.post_death_steps {
                    break;
                }
            }
        }
    }

    PKCellResult {
        dead: state.dead,
        death_step: state.death_step,
        final_lp: state.lp,
        final_gsh: state.gsh,
        final_gpx4: state.gpx4,
    }
}

// ============================================================
// Tests
// ============================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ode_reaches_steady_state_with_constant_plasma() {
        let plasma = PlasmaModel::Constant { concentration: 1.0 };
        let tumor = breast_tumor();
        let result = solve_tumor_pk(&plasma, &tumor, 500, 100);
        // After 500 minutes with constant input, C_i should be near steady state
        let last = *result.c_interstitial.last().unwrap();
        let second_last = result.c_interstitial[result.c_interstitial.len() - 2];
        assert!(
            (last - second_last).abs() < 1e-6,
            "C_i not at steady state: {last} vs {second_last}"
        );
        assert!(last > 0.0, "Steady-state C_i should be positive");
    }

    #[test]
    fn concentrations_never_negative() {
        let plasma = rsl3_iv_bolus();
        for tumor_fn in [breast_tumor, pancreatic_tumor, glioblastoma_tumor] {
            let tumor = tumor_fn();
            let result = solve_tumor_pk(&plasma, &tumor, 180, 100);
            for (i, &c) in result.c_interstitial.iter().enumerate() {
                assert!(c >= 0.0, "{}: C_i negative at step {i}: {c}", tumor.name);
            }
            for (i, &c) in result.c_vascular.iter().enumerate() {
                assert!(c >= 0.0, "{}: C_v negative at step {i}: {c}", tumor.name);
            }
        }
    }

    #[test]
    fn convergence_with_substep_doubling() {
        let plasma = rsl3_iv_bolus();
        let tumor = breast_tumor();
        let r100 = solve_tumor_pk(&plasma, &tumor, 180, 100);
        let r200 = solve_tumor_pk(&plasma, &tumor, 180, 200);
        // Peak C_i should agree within 1%
        let peak100: f64 = r100.c_interstitial.iter().cloned().fold(0.0, f64::max);
        let peak200: f64 = r200.c_interstitial.iter().cloned().fold(0.0, f64::max);
        let diff = (peak100 - peak200).abs() / peak200;
        assert!(diff < 0.01, "Convergence failed: {peak100} vs {peak200} ({diff:.4}%)");
    }

    #[test]
    fn breast_higher_exposure_than_gbm() {
        let plasma = rsl3_iv_bolus();
        let breast = solve_tumor_pk(&plasma, &breast_tumor(), 180, 100);
        let gbm = solve_tumor_pk(&plasma, &glioblastoma_tumor(), 180, 100);
        let auc_breast: f64 = breast.c_interstitial.iter().sum();
        let auc_gbm: f64 = gbm.c_interstitial.iter().sum();
        assert!(
            auc_breast > auc_gbm,
            "Breast AUC ({auc_breast}) should exceed GBM ({auc_gbm})"
        );
    }

    #[test]
    fn plasma_iv_bolus_decays() {
        let plasma = rsl3_iv_bolus();
        let c0 = plasma.concentration_at(0.0);
        let c30 = plasma.concentration_at(30.0);
        let c60 = plasma.concentration_at(60.0);
        assert!((c0 - 1.0).abs() < 1e-10, "C(0) should be 1.0");
        assert!((c30 - 0.5).abs() < 0.01, "C(t_half) should be ~0.5");
        assert!((c60 - 0.25).abs() < 0.01, "C(2×t_half) should be ~0.25");
    }
}
