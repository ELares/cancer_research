//! Combination Therapy Optimizer
//!
//! Combines temporal (vulnerability window) + immune (ICD cascade) models.
//! Sweeps SDT timing × anti-PD1 timing to find optimal treatment schedule.
//!
//! 3-phase model:
//! 1. Chemotherapy kills proliferating cells, persisters survive
//! 2. SDT at variable delay (0-28 days) targets persisters
//! 3. Immune cascade from ICD ± anti-PD1

use std::path::PathBuf;

use clap::Parser;
use rand::prelude::*;
use rayon::prelude::*;

use ferroptosis_core::biochem::sim_cell;
use ferroptosis_core::cell::{gen_recovered_persister, RecoveryRates, Treatment};
use ferroptosis_core::immune::immune_cascade;
use ferroptosis_core::io::write_json;
use ferroptosis_core::params::{ImmuneParams, Params};

#[derive(Parser)]
#[command(name = "sim-combo", about = "Combination therapy optimizer")]
struct Args {
    /// Cells per condition.
    #[arg(long, default_value_t = 50_000)]
    n_cells: usize,

    /// Random seed.
    #[arg(long, default_value_t = 42)]
    seed: u64,

    /// Output directory.
    #[arg(long, default_value = "output/combo")]
    output_dir: PathBuf,
}

// #585: SplitMix64 per-cell seed hash mix, ported verbatim from sim-tme-3d (#578),
// replacing the additive `seed + i*2(+1) + delay_days*1e6` scheme whose per-cell
// streams aliased across adjacent delay conditions once `2*(i-i') = 1e6` (n_cells >
// 500k). `cell_seed` fully decorrelates the (condition, cell, stream) triple, so
// there is no aliasing regime at ANY n_cells, and the two interleaved streams
// (cell-gen vs sim) are now separated by SALT instead of a `+1` offset. cond_seed
// still encodes ONLY delay_days (not tx / with_pd1) to preserve the deliberate
// paired design: at a fixed delay every treatment reuses the same cells + the same
// noise, a controlled comparison.
const COMBO_CELL_SALT: u64 = 0xC0B0_CE11_0000_0585;
const COMBO_SIM_SALT: u64 = 0xC0B0_5111_0000_0585;
const _: () = assert!(COMBO_CELL_SALT != COMBO_SIM_SALT);

#[inline]
fn splitmix64(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

#[inline]
fn cell_seed(cond_seed: u64, idx: usize, step: usize, salt: u64) -> u64 {
    splitmix64(
        cond_seed
            ^ salt
            ^ splitmix64(idx as u64).rotate_left(17)
            ^ splitmix64(step as u64).rotate_left(43),
    )
}

fn main() {
    let args = Args::parse();

    eprintln!("=== Combination Therapy Optimizer ===");
    eprintln!("Cells per condition: {}", args.n_cells);
    eprintln!("Seed: {}\n", args.seed);

    let params = Params::default();
    let recovery = RecoveryRates::default();
    let immune_params = ImmuneParams::default();

    // SDT delay: days after chemo withdrawal
    let sdt_delays_days: Vec<f64> = vec![0.0, 1.0, 3.0, 7.0, 14.0, 21.0, 28.0];

    // Anti-PD1 options
    let anti_pd1_options = [false, true];

    // Second-line treatment options
    let treatments = [
        (Treatment::SDT, "SDT"),
        (Treatment::PDT, "PDT"),
        (Treatment::RSL3, "RSL3"),
    ];

    std::fs::create_dir_all(&args.output_dir).expect("Failed to create output dir");

    let mut all_results = Vec::new();

    // Starting population: 1000 tumor cells post-chemo
    // Assume chemo killed 90% of proliferating cells, all survivors are persisters
    let initial_tumor_cells = 1000_usize;

    for &delay_days in &sdt_delays_days {
        for (tx, tx_name) in &treatments {
            for &with_pd1 in &anti_pd1_options {
                let n = args.n_cells;
                // #585: per-condition base seed (delay only — see the seed-helper
                // note above; the hash mix removes the old cross-condition aliasing).
                let cond_seed = args
                    .seed
                    .wrapping_add((delay_days as u64).wrapping_mul(1_000_000));

                // Phase 2: simulate ferroptosis on recovered persisters
                let outcomes: Vec<(bool, f64, f64, f64)> = (0..n)
                    .into_par_iter()
                    .map(|i| {
                        let mut cell_rng =
                            StdRng::seed_from_u64(cell_seed(cond_seed, i, 0, COMBO_CELL_SALT));
                        let mut sim_rng =
                            StdRng::seed_from_u64(cell_seed(cond_seed, i, 0, COMBO_SIM_SALT));
                        let cell = gen_recovered_persister(delay_days, &recovery, &mut cell_rng);
                        sim_cell(&cell, *tx, &params, &mut sim_rng)
                    })
                    .collect();

                let dead_cell_lps: Vec<f64> = outcomes
                    .iter()
                    .filter(|(dead, _, _, _)| *dead)
                    .map(|(_, lp, _, _)| *lp)
                    .collect();

                let ferroptosis_kill_rate = dead_cell_lps.len() as f64 / n as f64;

                // Phase 3: immune cascade
                // Scale dead_cell_lps to the biological population size (not simulation sample size)
                // to avoid saturating the immune cascade with 50K entries against 1K tumor cells.
                let ferroptosis_killed =
                    (ferroptosis_kill_rate * initial_tumor_cells as f64).round() as usize;
                let scaled_dead_lps: Vec<f64> = dead_cell_lps
                    .iter()
                    .take(ferroptosis_killed)
                    .copied()
                    .collect();
                let immune = immune_cascade(
                    &scaled_dead_lps,
                    initial_tumor_cells,
                    &immune_params,
                    with_pd1,
                );

                // Total tumor reduction: ferroptosis kills + immune kills
                let immune_killed = immune.immune_kills.round() as usize;
                let total_killed = (ferroptosis_killed + immune_killed).min(initial_tumor_cells);
                let survivors = initial_tumor_cells - total_killed;
                let survival_fraction = survivors as f64 / initial_tumor_cells as f64;

                eprintln!(
                    "  Day {:2.0} + {:<4} {}: ferro={:.1}%, immune={}, survivors={}/{} ({:.1}%)",
                    delay_days,
                    tx_name,
                    if with_pd1 { "+PD1" } else { "    " },
                    ferroptosis_kill_rate * 100.0,
                    immune_killed,
                    survivors,
                    initial_tumor_cells,
                    survival_fraction * 100.0,
                );

                all_results.push(serde_json::json!({
                    "sdt_delay_days": delay_days,
                    "treatment": tx_name,
                    "with_anti_pd1": with_pd1,
                    "initial_tumor_cells": initial_tumor_cells,
                    "ferroptosis_kill_rate": ferroptosis_kill_rate,
                    "ferroptosis_killed": ferroptosis_killed,
                    "immune_kills": immune.immune_kills,
                    "total_killed": total_killed,
                    "survivors": survivors,
                    "survival_fraction": survival_fraction,
                    "total_damps": immune.total_damps,
                    "damp_per_dead_cell": immune.damp_per_dead_cell,
                    "primed_tcells": immune.primed_tcells,
                }));
            }
        }
        eprintln!();
    }

    // Find optimal schedule
    let best = all_results.iter().min_by(|a, b| {
        a["survival_fraction"]
            .as_f64()
            .unwrap()
            .partial_cmp(&b["survival_fraction"].as_f64().unwrap())
            .unwrap()
    });

    if let Some(best) = best {
        eprintln!("=== Optimal Schedule ===");
        eprintln!(
            "  {} at day {:.0} {}: {:.1}% survival ({} survivors / {})",
            best["treatment"],
            best["sdt_delay_days"].as_f64().unwrap(),
            if best["with_anti_pd1"].as_bool().unwrap() {
                "+ anti-PD1"
            } else {
                ""
            },
            best["survival_fraction"].as_f64().unwrap() * 100.0,
            best["survivors"],
            best["initial_tumor_cells"],
        );
    }

    // Save results
    let json_path = args.output_dir.join("combo_results.json");
    write_json(&json_path, &all_results).expect("Failed to write JSON");

    eprintln!("\n=== Output saved to {} ===", args.output_dir.display());
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn cell_seed_has_no_additive_aliasing_585() {
        // The retired additive scheme aliased adjacent delay conditions once
        // `2*(i-i')` hit the 1e6 stride (n_cells > 500k). The hash mix has no such
        // regime: probe the old danger zone — seeds at and across the 500k onset
        // are all distinct, and the two interleaved streams never coincide (salt
        // separation replaces the old `+1` offset).
        let cond = 0xDEAD_BEEF_u64;
        let probe = [0usize, 1, 499_999, 500_000, 500_001, 1_000_000];
        let mut seen = HashSet::new();
        for &i in &probe {
            assert!(seen.insert(cell_seed(cond, i, 0, COMBO_CELL_SALT)));
            assert!(seen.insert(cell_seed(cond, i, 0, COMBO_SIM_SALT)));
        }
        assert_eq!(seen.len(), probe.len() * 2);
    }
}
