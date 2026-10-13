//! ICD-Immune Cascade Comparison
//!
//! Compares the per-dead-cell DAMP release of physical modalities against
//! pharmacologic ferroptosis induction.
//!
//! THIS BINARY'S KEY FINDING IS 1.1x, NOT THE 1.5-2x ITS README PROMISED
//! (#928). The header used to say "SDT/PDT kill via runaway LP cascade (LP >>
//! threshold at death), so they release MORE DAMPs per dead cell than RSL3
//! (LP ≈ threshold at death)" and the README's reproduction section told the
//! reader to expect ~1.5-2x. Run it: Persister SDT 20.0846 against RSL3
//! 17.7240 is **1.1332x**, and OXPHOS is 1.0748x. The direction holds in every
//! readable row; the magnitude does not reach the documented band.
//!
//! The "at death" half is also wrong, and this binary cannot see it. LP at the
//! moment of death is ~10.25 for BOTH arms against a `death_threshold` of 10.0
//! -- death IS the threshold crossing -- measured per cell in the committed 2D
//! archive (`analysis/immune-2d-measurements/observations.json.gz`: SDT
//! 10.2657, RSL3 10.2474). What this binary reports is `sim_cell`'s
//! `final_lp`, the LP after all 180 steps, which is why the field is now named
//! `avg_final_lp` rather than `avg_lp_at_death`.
//!
//! Analysis: `analysis/calibration/icd-damp-per-death.md`

use std::path::PathBuf;

use clap::Parser;
use rand::prelude::*;
use rayon::prelude::*;

use ferroptosis_core::biochem::sim_cell;
use ferroptosis_core::cell::{gen_cell, Phenotype, Treatment};
use ferroptosis_core::immune::immune_cascade;
use ferroptosis_core::io::write_json;
use ferroptosis_core::params::{ImmuneParams, Params};

#[derive(Parser)]
#[command(
    name = "sim-icd",
    about = "ICD-immune cascade comparison across treatments"
)]
struct Args {
    /// Cells per condition.
    #[arg(long, default_value_t = 100_000)]
    n_cells: usize,

    /// Random seed.
    #[arg(long, default_value_t = 42)]
    seed: u64,

    /// Output directory.
    #[arg(long, default_value = "output/icd")]
    output_dir: PathBuf,
}

fn main() {
    let args = Args::parse();

    eprintln!("=== ICD-Immune Cascade Simulation ===");
    eprintln!("Cells per condition: {}", args.n_cells);
    eprintln!("Seed: {}\n", args.seed);

    let params = Params::default();
    let immune_params = ImmuneParams::default();

    let phenotypes = [
        (Phenotype::Persister, "Persister"),
        (Phenotype::OXPHOS, "OXPHOS"),
        (Phenotype::Glycolytic, "Glycolytic"),
    ];
    let treatments = [
        (Treatment::Control, "Control"),
        (Treatment::RSL3, "RSL3"),
        (Treatment::SDT, "SDT"),
        (Treatment::PDT, "PDT"),
    ];

    std::fs::create_dir_all(&args.output_dir).expect("Failed to create output dir");

    let mut all_results = Vec::new();

    for (pheno, pheno_name) in &phenotypes {
        eprintln!("--- Phenotype: {} ---", pheno_name);

        for (tx, tx_name) in &treatments {
            let n = args.n_cells;

            // Run biochemistry simulation, collecting LP at death for dead cells
            let outcomes: Vec<(bool, f64, f64, f64)> = (0..n)
                .into_par_iter()
                .map(|i| {
                    let mut cell_rng = StdRng::seed_from_u64(args.seed.wrapping_add(i as u64 * 2));
                    let mut sim_rng =
                        StdRng::seed_from_u64(args.seed.wrapping_add(i as u64 * 2 + 1));
                    let cell = gen_cell(*pheno, &mut cell_rng);
                    sim_cell(&cell, *tx, &params, &mut sim_rng)
                })
                .collect();

            let dead_cell_lps: Vec<f64> = outcomes
                .iter()
                .filter(|(dead, _, _, _)| *dead)
                .map(|(_, lp, _, _)| *lp)
                .collect();

            let n_dead = dead_cell_lps.len();
            let death_rate = n_dead as f64 / n as f64;

            // NONE, NOT ZERO, WHEN NOTHING DIED (#928). This was `0.0`, and a
            // zero published under a column headed LP-per-dead-cell reads as a
            // measured potency of zero rather than as an absent cohort: at the
            // default config Glycolytic Control and Glycolytic RSL3 both kill
            // exactly 0 cells and both used to report 0.0. `damp_per_dead_cell`
            // carries the same sentinel from `calculate_damp_release`, so it is
            // withheld on the same condition rather than republished.
            let avg_final_lp = if n_dead > 0 {
                Some(dead_cell_lps.iter().sum::<f64>() / n_dead as f64)
            } else {
                None
            };

            // Run immune cascade (with and without anti-PD1)
            let immune_no_pd1 = immune_cascade(&dead_cell_lps, n, &immune_params, false);
            let immune_with_pd1 = immune_cascade(&dead_cell_lps, n, &immune_params, true);
            let damp_no_pd1 = (n_dead > 0).then_some(immune_no_pd1.damp_per_dead_cell);
            let damp_with_pd1 = (n_dead > 0).then_some(immune_with_pd1.damp_per_dead_cell);

            let fmt = |v: Option<f64>| match v {
                Some(x) => format!("{x:5.1}"),
                None => "  n/a".to_string(),
            };
            eprintln!(
                "  {:<8} → Dead: {:6.2}%, final LP/cell: {}, DAMP/cell: {}, T-kills: {:.0} / {:.0} (±PD1)",
                tx_name,
                death_rate * 100.0,
                fmt(avg_final_lp),
                fmt(damp_no_pd1),
                immune_no_pd1.immune_kills,
                immune_with_pd1.immune_kills,
            );

            all_results.push(serde_json::json!({
                "phenotype": pheno_name,
                "treatment": tx_name,
                "n_cells": n,
                "n_dead": n_dead,
                "death_rate": death_rate,
                // `final_lp` after all 180 steps, NOT LP at the threshold
                // crossing (which is ~10.25 for every arm). Renamed from
                // `avg_lp_at_death` in #928 because the old name and the
                // README's data dictionary both named the other quantity.
                "avg_final_lp": avg_final_lp,
                "immune_no_pd1": {
                    "total_damps": immune_no_pd1.total_damps,
                    "damp_per_dead_cell": damp_no_pd1,
                    "dc_activation_fraction": immune_no_pd1.dc_activation_fraction,
                    "mature_dcs": immune_no_pd1.mature_dcs,
                    "primed_tcells": immune_no_pd1.primed_tcells,
                    "immune_kills": immune_no_pd1.immune_kills,
                },
                "immune_with_pd1": {
                    "total_damps": immune_with_pd1.total_damps,
                    "damp_per_dead_cell": damp_with_pd1,
                    "dc_activation_fraction": immune_with_pd1.dc_activation_fraction,
                    "mature_dcs": immune_with_pd1.mature_dcs,
                    "primed_tcells": immune_with_pd1.primed_tcells,
                    "immune_kills": immune_with_pd1.immune_kills,
                },
            }));
        }
        eprintln!();
    }

    // Save
    let json_path = args.output_dir.join("icd_comparison.json");
    write_json(&json_path, &all_results).expect("Failed to write JSON");

    eprintln!("=== Key Finding ===");
    // Compare SDT vs RSL3 DAMP per dead cell for persisters
    let sdt_result = all_results
        .iter()
        .find(|r| r["phenotype"] == "Persister" && r["treatment"] == "SDT");
    let rsl3_result = all_results
        .iter()
        .find(|r| r["phenotype"] == "Persister" && r["treatment"] == "RSL3");
    if let (Some(sdt), Some(rsl3)) = (sdt_result, rsl3_result) {
        // ABSENT, NOT ZERO, ON BOTH SIDES. `damp_per_dead_cell` became
        // nullable in #928, and `unwrap_or(0.0)` on a nullable field is the
        // very defect that change removed -- it would print `0.000x` as a
        // measured ratio for an arm that killed nothing. Caught in self-review
        // of #928 itself, one level above the fix.
        let sdt_damp = sdt["immune_no_pd1"]["damp_per_dead_cell"].as_f64();
        let rsl3_damp = rsl3["immune_no_pd1"]["damp_per_dead_cell"].as_f64();
        if let (Some(sdt_damp), Some(rsl3_damp)) = (sdt_damp, rsl3_damp) {
            // THREE DECIMALS, DELIBERATELY (#928). At one decimal this printed
            // "1.1×" where the README promised 1.5-2x, and the two were near
            // enough in print that nobody read them against each other. The
            // documented band is named here so the line cannot be read as a
            // confirmation of it.
            eprintln!(
                "SDT produces {:.3}× the DAMP per dead cell of RSL3 ({:.4} vs {:.4}).",
                sdt_damp / rsl3_damp,
                sdt_damp,
                rsl3_damp,
            );
            eprintln!(
                "  DIRECTION only. sim-icd/README.md once said 1.5-2x and the v2 plan \
                 2-5x; neither is reached, here or in any other regime (#928)."
            );
        }
        let sdt_kills = sdt["immune_with_pd1"]["immune_kills"]
            .as_f64()
            .unwrap_or(0.0);
        let rsl3_kills = rsl3["immune_with_pd1"]["immune_kills"]
            .as_f64()
            .unwrap_or(0.0);
        let sdt_rate = sdt["death_rate"].as_f64().unwrap_or(0.0);
        let rsl3_rate = rsl3["death_rate"].as_f64().unwrap_or(0.0);
        eprintln!(
            "SDT+anti-PD1 immune kills: {:.0} vs RSL3+anti-PD1: {:.0}",
            sdt_kills, rsl3_kills,
        );
        // SAY WHY THE LOSER WON (#928). `immune_cascade` caps its kills at the
        // SURVIVING population, so an arm that kills everything biochemically
        // has nothing left to earn immune kills on. At the default config SDT
        // kills 99.997% of persisters and scores 3 against RSL3's 57,591 at
        // 42.4% -- the opposite of the README's "SDT >> RSL3", and decided by
        // who left survivors rather than by immunogenicity. sim-modality-panel
        // documents the same trap; this block printed the comparison without
        // it. The warning is DERIVED from the death rates, not a fixed string,
        // so it stops appearing if an arm stops saturating.
        if sdt_rate > 0.99 || rsl3_rate > 0.99 {
            eprintln!(
                "  NOT COMPARABLE: {} kills {:.3}% biochemically, and immune_cascade caps \
                 kills at the surviving population, so this ranks survivors not \
                 immunogenicity (#928).",
                if sdt_rate > rsl3_rate { "SDT" } else { "RSL3" },
                sdt_rate.max(rsl3_rate) * 100.0,
            );
        }
    }

    eprintln!("\n=== Output saved to {} ===", args.output_dir.display());
}
