//! ICD-immune cascade model.
//!
//! Models: ferroptotic death → DAMP release → DC maturation → T cell priming → tumor killing.
//! Deliberately simple (3-stage linear cascade) to show qualitative ICD differences
//! between physical modalities and pharmacologic inducers.

use serde::Serialize;

use crate::params::ImmuneParams;

/// Result of the immune cascade for one treatment condition.
#[derive(Clone, Debug, Serialize)]
pub struct ImmuneResult {
    /// Total DAMP signal from all dead cells.
    pub total_damps: f64,
    /// Average DAMP per dead cell (higher = more immunogenic death).
    pub damp_per_dead_cell: f64,
    /// Number of dead cells that contributed DAMPs.
    pub n_dead_cells: usize,
    /// Fraction of DCs activated by DAMP signal.
    pub dc_activation_fraction: f64,
    /// Number of mature DCs.
    pub mature_dcs: f64,
    /// Number of primed T cells.
    pub primed_tcells: f64,
    /// Tumor cells killed by immune response.
    pub immune_kills: f64,
    /// Whether anti-PD-1 was applied.
    pub with_anti_pd1: bool,
}

/// Calculate total DAMP release from a set of dead cells.
///
/// DAMP is proportional to lipid peroxidation at death. Key biological insight:
/// - SDT/PDT kill via runaway LP cascade → LP at death is MUCH higher than threshold
/// - RSL3 kills via slow GPX4 inhibition → LP at death is NEAR threshold
/// Therefore SDT/PDT-killed cells release MORE DAMPs per cell.
///
/// Ref: Krysko et al., Nat Rev Cancer 2012 (ICD markers)
///      Berezhnoy et al., PLoS Comput Biol 2020 (Boolean ICD model)
pub fn calculate_damp_release(
    dead_cell_lps: &[f64],
    params: &ImmuneParams,
) -> (f64, f64) {
    if dead_cell_lps.is_empty() {
        return (0.0, 0.0);
    }
    let total: f64 = dead_cell_lps.iter().map(|lp| lp * params.damp_per_lp).sum();
    let per_cell = total / dead_cell_lps.len() as f64;
    (total, per_cell)
}

/// Run the DC maturation → T cell priming → killing cascade.
///
/// Model stages:
/// 1. DAMPs activate DCs: saturating Michaelis-Menten: DAMP / (DAMP + Kd)
/// 2. Activated DCs mature with probability dc_maturation_rate
/// 3. Each mature DC primes tcell_priming_rate T cells
/// 4. Each T cell kills tcell_kill_rate tumor cells
/// 5. PD-1 brake suppresses a fraction of T-cell killing
/// 6. Anti-PD-1 removes a fraction of the brake
///
/// This is a coarse-grained model showing QUALITATIVE differences, not absolute numbers.
pub fn immune_cascade(
    dead_cell_lps: &[f64],
    total_tumor_cells: usize,
    params: &ImmuneParams,
    with_anti_pd1: bool,
) -> ImmuneResult {
    let (total_damps, damp_per_dead_cell) = calculate_damp_release(dead_cell_lps, params);
    let n_dead = dead_cell_lps.len();

    // DC activation: saturating response to total DAMP
    let dc_activation_fraction = total_damps / (total_damps + params.dc_activation_kd);

    // Mature DCs (proportional to activation and number of dead cells presenting antigens)
    let mature_dcs = dc_activation_fraction * params.dc_maturation_rate * n_dead as f64;

    // T cell priming
    let primed_tcells = mature_dcs * params.tcell_priming_rate;

    // T cell killing (with PD-1 brake)
    let effective_brake = if with_anti_pd1 {
        params.pd1_brake * (1.0 - params.anti_pd1_efficacy)
    } else {
        params.pd1_brake
    };
    let kill_efficiency = 1.0 - effective_brake;
    let raw_kills = primed_tcells * params.tcell_kill_rate * kill_efficiency;

    // Cap at remaining alive tumor cells
    let remaining = total_tumor_cells.saturating_sub(n_dead) as f64;
    let immune_kills = raw_kills.min(remaining);

    ImmuneResult {
        total_damps,
        damp_per_dead_cell,
        n_dead_cells: n_dead,
        dc_activation_fraction,
        mature_dcs,
        primed_tcells,
        immune_kills,
        with_anti_pd1,
    }
}
