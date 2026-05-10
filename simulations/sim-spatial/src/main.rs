//! Spatial Tumor Simulation with Energy Physics
//!
//! Validates the paper's core claims:
//! 1. "Physical modalities restrict ROS to a defined tumor volume"
//! 2. "Ultrasound penetrates centimeters; light penetrates millimeters"
//!
//! Simulates a 2D heterogeneous tumor with depth-dependent energy deposition
//! from PDT (Beer-Lambert), SDT (acoustic attenuation), and RSL3 (uniform).

use std::path::PathBuf;

use clap::Parser;
use rand::prelude::*;

use ferroptosis_core::biochem::{sim_cell_step, CellState};
use ferroptosis_core::cell::{norm, Treatment};
use ferroptosis_core::grid::{depth_kill_curve, death_heatmap, TumorGrid};
use ferroptosis_core::io::{write_depth_curves_csv, write_heatmap_csv, write_json};
use ferroptosis_core::params::{Params, SpatialParams};
use ferroptosis_core::photosensitizer_pk::{validate_dli_h, Photosensitizer};
use ferroptosis_core::physics::local_ros_multiplier;

#[derive(Parser)]
#[command(name = "sim-spatial", about = "Spatial tumor ferroptosis simulation")]
struct Args {
    /// Grid size (rows = cols).
    #[arg(long, default_value_t = 500)]
    grid_size: usize,

    /// Cell size in micrometers.
    #[arg(long, default_value_t = 20.0)]
    cell_size: f64,

    /// Random seed.
    #[arg(long, default_value_t = 42)]
    seed: u64,

    /// Output directory.
    #[arg(long, default_value = "output/spatial")]
    output_dir: PathBuf,

    /// Number of biochemistry timesteps per cell.
    #[arg(long, default_value_t = 180)]
    n_steps: u32,

    /// Photosensitizer PK model for PDT light scaling. Spec format
    /// (case-insensitive): `uniform` (= `uniform=1.0`, the default),
    /// `uniform=N` (constant fraction; values >1.0 represent enrichment
    /// rather than the typical [0, 1] drug-presence range, intentional
    /// forward-compat hook), `porfimer` (= `porfimer=504`, Bellnier 2006
    /// t½ in hours), or `porfimer=N` (custom t½ in hours).
    #[arg(long, default_value = "uniform")]
    photosensitizer: String,

    /// Drug-light interval in hours: time from photosensitizer
    /// post-distribution peak to light delivery. NOT the clinical DLI
    /// from injection — see ferroptosis_core::photosensitizer_pk for the
    /// distinction. Default 0.0 means light at peak.
    #[arg(long, default_value_t = 0.0)]
    dli_h: f64,
}


fn run_spatial(
    grid: &mut TumorGrid,
    tx: Treatment,
    params: &Params,
    spatial_params: &SpatialParams,
    n_steps: u32,
    seed: u64,
) {
    let base_ros = match tx {
        Treatment::SDT => params.sdt_ros,
        Treatment::PDT => params.pdt_ros,
        Treatment::RSL3 | Treatment::Control => 0.0,
    };

    let rows = grid.rows;
    let cols = grid.cols;
    let cell_size = grid.cell_size_um;

    // Initialize cell states with depth-dependent ROS
    for r in 0..rows {
        let ros_multiplier = local_ros_multiplier(r, cell_size, tx, spatial_params);
        for c in 0..cols {
            let exo_ros_peak = if tx == Treatment::Control || tx == Treatment::RSL3 {
                0.0
            } else {
                let mut rng = StdRng::seed_from_u64(seed.wrapping_add((r * cols + c) as u64));
                let peak = base_ros * ros_multiplier;
                norm(&mut rng, peak, peak * 0.2).max(0.0)
            };

            let gc = grid.get_mut(r, c);
            gc.state = CellState::from_cell_with_ros(&gc.cell, tx, params, exo_ros_peak);
            gc.extra_iron = 0.0;
            gc.newly_dead = false;
            gc.lp_at_death = 0.0;
        }
    }

    // Run simulation: 180 timesteps with interleaved diffusion
    for step in 0..n_steps {
        // Biochemistry step for each cell
        for r in 0..rows {
            for c in 0..cols {
                let idx = r * cols + c;
                if !grid.cells[idx].is_tumor {
                    continue;
                }
                if grid.cells[idx].state.dead {
                    if let Some(ds) = grid.cells[idx].state.death_step {
                        if step >= ds + params.post_death_steps { continue; }
                    } else {
                        continue;
                    }
                }

                // Seed incorporates a large offset to avoid collision with init seeds
                let mut rng = StdRng::seed_from_u64(
                    seed.wrapping_add(500_000)
                        .wrapping_add(idx as u64)
                        .wrapping_add(step as u64 * 1_000_000),
                );

                // Consume extra_iron this step, then reset (iron is consumed by Fenton reaction)
                let extra_iron = grid.cells[idx].extra_iron;
                grid.cells[idx].extra_iron = 0.0;

                let gc = &mut grid.cells[idx];
                let died = sim_cell_step(
                    &mut gc.state,
                    &gc.cell,
                    params,
                    step,
                    extra_iron,
                    &mut rng,
                );

                if died {
                    gc.newly_dead = true;
                    gc.lp_at_death = gc.state.lp;
                }
                // Update lp_at_death during grace period
                if gc.state.dead {
                    gc.lp_at_death = gc.state.lp;
                }
            }
        }

        // Diffusion step: distribute iron from newly dead cells
        grid.diffuse_iron(
            spatial_params.iron_release_per_death,
            spatial_params.neighbor_iron_fraction,
        );

        // Progress reporting every 30 steps
        if (step + 1) % 30 == 0 {
            let census = grid.census();
            eprintln!(
                "  Step {}/{}: {}/{} tumor cells dead ({:.1}%)",
                step + 1,
                n_steps,
                census.total_dead,
                census.total_tumor,
                if census.total_tumor > 0 {
                    census.total_dead as f64 / census.total_tumor as f64 * 100.0
                } else {
                    0.0
                }
            );
        }
    }
}

fn main() {
    let args = Args::parse();

    let photosensitizer = match Photosensitizer::from_cli_spec(&args.photosensitizer) {
        Ok(ps) => ps,
        Err(e) => {
            eprintln!("error: --photosensitizer {:?}: {}", args.photosensitizer, e);
            std::process::exit(2);
        }
    };
    if let Err(e) = validate_dli_h(args.dli_h) {
        eprintln!("error: {e}");
        std::process::exit(2);
    }

    eprintln!("=== Spatial Tumor Ferroptosis Simulation ===");
    eprintln!(
        "Grid: {}×{} cells ({:.1}mm × {:.1}mm tissue)",
        args.grid_size,
        args.grid_size,
        args.grid_size as f64 * args.cell_size / 1000.0,
        args.grid_size as f64 * args.cell_size / 1000.0,
    );
    eprintln!("Cell size: {} µm", args.cell_size);
    eprintln!("Seed: {}", args.seed);
    eprintln!("Photosensitizer: {photosensitizer}, DLI: {} h\n", args.dli_h);

    let params = Params::default();
    let spatial_params = SpatialParams {
        cell_size_um: args.cell_size,
        photosensitizer,
        t_drug_light_interval_h: args.dli_h,
        ..Default::default()
    };

    let treatments = [
        (Treatment::Control, "Control"),
        (Treatment::RSL3, "RSL3"),
        (Treatment::SDT, "SDT"),
        (Treatment::PDT, "PDT"),
    ];

    // Create output directory
    std::fs::create_dir_all(&args.output_dir).expect("Failed to create output directory");

    let mut all_depth_curves = Vec::new();
    let mut all_summaries = Vec::new();

    for (tx, tx_name) in &treatments {
        eprintln!("--- Treatment: {} ---", tx_name);

        // Generate fresh grid for each treatment (same seed = same tumor)
        let mut grid = TumorGrid::generate(args.grid_size, args.grid_size, args.cell_size, args.seed);

        let census_before = grid.census();
        eprintln!(
            "  Tumor composition: {} glycolytic, {} OXPHOS, {} persister, {} persister+NRF2, {} stromal",
            census_before.glycolytic,
            census_before.oxphos,
            census_before.persister,
            census_before.persister_nrf2,
            census_before.stromal,
        );

        // Run spatial simulation
        run_spatial(
            &mut grid,
            *tx,
            &params,
            &spatial_params,
            args.n_steps,
            args.seed.wrapping_add((*tx as u64) * 10_000_000),
        );

        let census_after = grid.census();
        eprintln!(
            "  Final: {}/{} dead ({:.1}%)",
            census_after.total_dead,
            census_after.total_tumor,
            if census_after.total_tumor > 0 {
                census_after.total_dead as f64 / census_after.total_tumor as f64 * 100.0
            } else {
                0.0
            }
        );
        eprintln!(
            "    Glycolytic: {}/{} dead, OXPHOS: {}/{} dead, Persister: {}/{} dead, NRF2: {}/{} dead",
            census_after.glycolytic_dead, census_after.glycolytic,
            census_after.oxphos_dead, census_after.oxphos,
            census_after.persister_dead, census_after.persister,
            census_after.persister_nrf2_dead, census_after.persister_nrf2,
        );
        eprintln!();

        // Save death heatmap
        let heatmap = death_heatmap(&grid);
        let heatmap_path = args.output_dir.join(format!("spatial_death_{}.csv", tx_name.to_lowercase()));
        write_heatmap_csv(&heatmap_path, &heatmap).expect("Failed to write heatmap");

        // Compute depth-kill curve
        let curve = depth_kill_curve(&grid);
        all_depth_curves.push((tx_name.to_string(), curve));

        // Summary
        all_summaries.push(serde_json::json!({
            "treatment": tx_name,
            "total_tumor": census_after.total_tumor,
            "total_dead": census_after.total_dead,
            "overall_death_rate": census_after.total_dead as f64 / census_after.total_tumor.max(1) as f64,
            "glycolytic": { "total": census_after.glycolytic, "dead": census_after.glycolytic_dead },
            "oxphos": { "total": census_after.oxphos, "dead": census_after.oxphos_dead },
            "persister": { "total": census_after.persister, "dead": census_after.persister_dead },
            "persister_nrf2": { "total": census_after.persister_nrf2, "dead": census_after.persister_nrf2_dead },
        }));
    }

    // Save depth-kill curves (all treatments in one CSV)
    let curves_path = args.output_dir.join("depth_kill_curves.csv");
    write_depth_curves_csv(&curves_path, &all_depth_curves).expect("Failed to write depth curves");

    // Save summary JSON
    let summary_path = args.output_dir.join("spatial_summary.json");
    write_json(&summary_path, &all_summaries).expect("Failed to write summary");

    eprintln!("=== Output saved to {} ===", args.output_dir.display());
    eprintln!("  depth_kill_curves.csv — death rate by depth for all treatments");
    eprintln!("  spatial_death_{{treatment}}.csv — 2D death heatmaps");
    eprintln!("  spatial_summary.json — aggregate statistics");
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Spec / DLI parsers and their edge-case rejection logic live in
    /// `ferroptosis_core::photosensitizer_pk`; their unit tests live with
    /// them. The remaining test here covers the only binary-specific
    /// concern: that clap's `default_value` attributes still map to
    /// `SpatialParams::default()`. Without this, the parser-level
    /// invariant could pass while the clap defaults silently drifted.
    #[test]
    fn default_args_match_default_spatial_params() {
        let args = Args::parse_from(["sim-spatial"]);
        let parsed_ps = Photosensitizer::from_cli_spec(&args.photosensitizer).unwrap();
        let default_sp = SpatialParams::default();
        assert_eq!(parsed_ps, default_sp.photosensitizer);
        assert_eq!(args.dli_h, default_sp.t_drug_light_interval_h);
        // Also pin the canonical CLI string so a rename would be visible
        // in stderr but not just in the parsed enum value.
        assert_eq!(args.photosensitizer, "uniform");
    }
}
