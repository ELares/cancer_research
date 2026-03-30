//! Cell types, phenotypes, treatments, and cell generation.

use rand::prelude::*;
use rand_distr::Normal;
use serde::{Deserialize, Serialize};

/// Sample from a normal distribution. Does NOT clamp — callers must apply
/// `.max(threshold)` if a positive value is required.
pub fn norm(rng: &mut StdRng, mean: f64, sd: f64) -> f64 {
    Normal::new(mean, sd).unwrap().sample(rng)
}

/// A single cell's biochemical state.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Cell {
    /// Labile iron pool (µM). Ref: 0.2-1.5 normal, 2-6 overloaded
    /// (Kakhlon & Cabantchik, Free Radic Biol Med 2002)
    pub iron: f64,
    /// Glutathione (mM). Ref: 1-10 mM intracellular
    /// (Forman et al., Free Radic Biol Med 2009)
    pub gsh: f64,
    /// GPX4 activity (relative, 1.0 = normal). Ref: kcat ~40/s
    /// (Ursini et al., Free Radic Biol Med 1995)
    pub gpx4: f64,
    /// FSP1/DHODH activity (relative). GPX4-independent CoQ10 pathway.
    /// (Bersuker et al., Nature 2019; Mao et al., Nature 2021)
    pub fsp1: f64,
    /// Basal mitochondrial ROS production (relative).
    /// OXPHOS cells: ~2-3× higher due to active ETC (Murphy, Biochem J 2009)
    pub basal_ros: f64,
    /// Lipid unsaturation: PUFA content determines peroxidation susceptibility.
    /// OXPHOS cells have more mitochondrial membranes = more target.
    /// (Yang et al., Cell 2016 — PUFA requirement for ferroptosis)
    pub lipid_unsat: f64,
    /// NRF2 transcriptional activity. Master regulator of antioxidant response.
    /// Drives GSH synthesis (via GCL/GSS), GPX4 expression.
    /// (Dodson et al., Free Radic Biol Med 2019)
    pub nrf2: f64,
}

/// Treatment modalities.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Treatment {
    Control,
    RSL3,
    SDT,
    PDT,
}

/// Cell phenotypes.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Phenotype {
    Glycolytic,
    OXPHOS,
    Persister,
    PersisterNrf2,
    /// Non-tumor stromal cells (for spatial model). Low iron, high GSH, inert.
    Stromal,
}

/// Generate a cell with stochastic parameters for the given phenotype.
/// Identical to the original v3 simulation for the four base phenotypes.
pub fn gen_cell(pheno: Phenotype, rng: &mut StdRng) -> Cell {
    match pheno {
        Phenotype::Glycolytic => Cell {
            iron: norm(rng, 1.0, 0.25).max(0.3),
            gsh: norm(rng, 5.0, 1.0).max(1.5),
            gpx4: norm(rng, 1.0, 0.12).max(0.4),
            fsp1: norm(rng, 1.0, 0.12).max(0.4),
            basal_ros: norm(rng, 0.2, 0.05).max(0.05),
            lipid_unsat: norm(rng, 1.0, 0.12).max(0.5),
            nrf2: norm(rng, 1.0, 0.12).max(0.4),
        },
        Phenotype::OXPHOS => Cell {
            iron: norm(rng, 2.8, 0.6).max(0.8),
            gsh: norm(rng, 4.0, 0.8).max(1.0),
            gpx4: norm(rng, 1.0, 0.12).max(0.4),
            fsp1: norm(rng, 1.0, 0.12).max(0.4),
            basal_ros: norm(rng, 0.5, 0.12).max(0.1),
            lipid_unsat: norm(rng, 1.6, 0.2).max(0.7),
            nrf2: norm(rng, 1.2, 0.15).max(0.5),
        },
        Phenotype::Persister => Cell {
            iron: norm(rng, 1.5, 0.3).max(0.5),
            gsh: norm(rng, 4.8, 0.8).max(1.8),
            gpx4: norm(rng, 0.7, 0.15).max(0.15),
            fsp1: norm(rng, 0.15, 0.06).max(0.01),
            basal_ros: norm(rng, 0.25, 0.06).max(0.05),
            lipid_unsat: norm(rng, 1.4, 0.15).max(0.6),
            nrf2: norm(rng, 0.7, 0.15).max(0.2),
        },
        Phenotype::PersisterNrf2 => Cell {
            iron: norm(rng, 2.8, 0.6).max(0.8),
            gsh: norm(rng, 7.0, 1.2).max(3.0),
            gpx4: norm(rng, 1.3, 0.15).max(0.5),
            fsp1: norm(rng, 0.2, 0.08).max(0.02),
            basal_ros: norm(rng, 0.5, 0.12).max(0.1),
            lipid_unsat: norm(rng, 1.6, 0.2).max(0.7),
            nrf2: norm(rng, 3.0, 0.4).max(1.5),
        },
        Phenotype::Stromal => Cell {
            iron: norm(rng, 0.3, 0.08).max(0.1),
            gsh: norm(rng, 8.0, 1.0).max(4.0),
            gpx4: norm(rng, 1.5, 0.15).max(0.8),
            fsp1: norm(rng, 1.0, 0.12).max(0.4),
            basal_ros: norm(rng, 0.1, 0.03).max(0.02),
            lipid_unsat: norm(rng, 0.6, 0.1).max(0.3),
            nrf2: norm(rng, 1.0, 0.12).max(0.4),
        },
    }
}

/// Recovery rates for persister phenotype transitions post-chemotherapy.
/// Used by the vulnerability window simulation.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RecoveryRates {
    /// FSP1 re-expression half-time (days). Slowest — epigenetic.
    pub fsp1_half_recovery_days: f64,
    /// GPX4 re-expression half-time (days). Transcriptional.
    pub gpx4_half_recovery_days: f64,
    /// NRF2 re-activation half-time (days).
    pub nrf2_half_recovery_days: f64,
    /// GSH resynthesis half-time (days). Fastest — metabolic.
    pub gsh_half_recovery_days: f64,
}

impl Default for RecoveryRates {
    fn default() -> Self {
        RecoveryRates {
            fsp1_half_recovery_days: 7.0,
            gpx4_half_recovery_days: 3.0,
            nrf2_half_recovery_days: 5.0,
            gsh_half_recovery_days: 1.0,
        }
    }
}

/// Apply time-dependent recovery to a persister cell.
/// Returns new parameter means (not a full Cell — caller generates stochastic cell from these).
pub fn recovered_persister_means(days: f64, rates: &RecoveryRates) -> (f64, f64, f64, f64) {
    // Exponential recovery: fraction recovered = 1 - exp(-ln(2) * t / t_half)
    let frac = |t_half: f64| -> f64 {
        1.0 - (-(2.0_f64.ln()) * days / t_half).exp()
    };

    // Persister baseline → Glycolytic normal targets
    let fsp1 = 0.15 + (1.0 - 0.15) * frac(rates.fsp1_half_recovery_days);
    let gpx4 = 0.7 + (1.0 - 0.7) * frac(rates.gpx4_half_recovery_days);
    let nrf2 = 0.7 + (1.0 - 0.7) * frac(rates.nrf2_half_recovery_days);
    let gsh = 4.8 + (5.0 - 4.8) * frac(rates.gsh_half_recovery_days);

    (fsp1, gpx4, nrf2, gsh)
}

/// Generate a persister cell at a given number of days post-chemo withdrawal.
pub fn gen_recovered_persister(days: f64, rates: &RecoveryRates, rng: &mut StdRng) -> Cell {
    let (fsp1_mean, gpx4_mean, nrf2_mean, gsh_mean) = recovered_persister_means(days, rates);
    Cell {
        iron: norm(rng, 1.5, 0.3).max(0.5),
        gsh: norm(rng, gsh_mean, 0.8).max(1.8),
        gpx4: norm(rng, gpx4_mean, 0.15).max(0.15),
        fsp1: norm(rng, fsp1_mean, 0.06).max(0.01),
        basal_ros: norm(rng, 0.25, 0.06).max(0.05),
        lipid_unsat: norm(rng, 1.4, 0.15).max(0.6),
        nrf2: norm(rng, nrf2_mean, 0.15).max(0.2),
    }
}
