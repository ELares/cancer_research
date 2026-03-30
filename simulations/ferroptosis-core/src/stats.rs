//! Statistics, results aggregation, and parallel execution.

use rayon::prelude::*;
use rand::prelude::*;
use serde::Serialize;

use crate::cell::{gen_cell, Phenotype, Treatment};
use crate::biochem::sim_cell;
use crate::params::Params;

/// Aggregated result for one condition (phenotype × treatment).
#[derive(Serialize, Clone, Debug)]
pub struct SimResult {
    pub phenotype: String,
    pub treatment: String,
    pub n_cells: usize,
    pub n_dead: usize,
    pub death_rate: f64,
    pub ci_low: f64,
    pub ci_high: f64,
    pub mean_lipid_perox: f64,
    pub mean_gsh_final: f64,
    pub mean_gpx4_final: f64,
}

/// Wilson score confidence interval for binomial proportion.
pub fn wilson_ci(n: usize, k: usize) -> (f64, f64) {
    let (nf, p, z) = (n as f64, k as f64 / n as f64, 1.96);
    let d = 1.0 + z * z / nf;
    let c = (p + z * z / (2.0 * nf)) / d;
    let s = z * ((p * (1.0 - p) / nf + z * z / (4.0 * nf * nf)).sqrt()) / d;
    ((c - s).max(0.0), (c + s).min(1.0))
}

/// Run a single condition (phenotype × treatment) with n cells in parallel.
/// Uses the same RNG seeding as the original v3 for bitwise reproducibility.
pub fn run_condition(
    pheno: Phenotype,
    tx: Treatment,
    params: &Params,
    n: usize,
    pname: &str,
    tname: &str,
) -> SimResult {
    let outcomes: Vec<(bool, f64, f64, f64)> = (0..n)
        .into_par_iter()
        .map(|i| {
            let mut cell_rng = StdRng::seed_from_u64(i as u64 * 2);
            let mut sim_rng = StdRng::seed_from_u64(i as u64 * 2 + 1);
            let cell = gen_cell(pheno, &mut cell_rng);
            sim_cell(&cell, tx, params, &mut sim_rng)
        })
        .collect();

    let dead = outcomes.iter().filter(|(d, _, _, _)| *d).count();
    let dr = dead as f64 / n as f64;
    let (cl, ch) = wilson_ci(n, dead);

    SimResult {
        phenotype: pname.to_string(),
        treatment: tname.to_string(),
        n_cells: n,
        n_dead: dead,
        death_rate: dr,
        ci_low: cl,
        ci_high: ch,
        mean_lipid_perox: outcomes.iter().map(|(_, l, _, _)| l).sum::<f64>() / n as f64,
        mean_gsh_final: outcomes.iter().map(|(_, _, g, _)| g).sum::<f64>() / n as f64,
        mean_gpx4_final: outcomes.iter().map(|(_, _, _, p)| p).sum::<f64>() / n as f64,
    }
}
