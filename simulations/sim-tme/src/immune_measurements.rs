//! Canonical 2D wrapper for the shared passive event ledger.
//!
//! Geometry and treatment seeds intentionally reproduce the historical 2D
//! matrix. They do not use the 3D condition-name hashing scheme.
//! The separate replicate entry point changes only the initialization seed
//! and corresponding treatment seeds, retaining the canonical observation mode.

use std::ffi::OsString;

use crate::passive_immune_measurements::Measurements;
use ferroptosis_core::grid::TUMOR_RADIUS_FRACTION;

use crate::*;

const GRID_SEED: u64 = 42;
const TREATMENT_SEED_STRIDE: u64 = 10_000_000;
const REPLICATE_BLOCKS: u32 = 20;
const REPLICATE_SEED_STRIDE: u64 = 1_u64 << 32;
const ARMS: [(Treatment, &str); 3] = [
    (Treatment::Control, "Control"),
    (Treatment::RSL3, "RSL3"),
    (Treatment::SDT, "SDT"),
];

pub(super) fn validate_args(args: &[String], ferro_env: &[OsString]) {
    assert!(
        args.len() == 2 && args[1] == "--immune-measurements",
        "--immune-measurements is a standalone canonical comparison; extra arguments are unsupported"
    );
    assert!(
        ferro_env.is_empty(),
        "--immune-measurements requires every FERRO_* environment variable to be absent"
    );
}

pub(super) fn replicate_requested(args: &[String]) -> bool {
    args.iter().any(|arg| arg.starts_with("--immune-replicate"))
}

pub(super) fn validate_replicate_args(args: &[String], ferro_env: &[OsString]) -> u32 {
    assert!(
        args.len() == 3 && args[1] == "--immune-replicate",
        "use --immune-replicate <block_id>; extra or mixed mode arguments are unsupported"
    );
    assert!(
        ferro_env.is_empty(),
        "--immune-replicate requires every FERRO_* environment variable to be absent"
    );
    let block: u32 = args[2]
        .parse()
        .expect("replicate block must be an integer from 1 through 20");
    assert!(
        (1..=REPLICATE_BLOCKS).contains(&block) && args[2] == block.to_string(),
        "replicate block must be written as an integer from 1 through 20"
    );
    block
}

fn replicate_seed(block: u32) -> u64 {
    assert!((1..=REPLICATE_BLOCKS).contains(&block));
    GRID_SEED + u64::from(block) * REPLICATE_SEED_STRIDE
}

fn condition_seed(tx: Treatment, grid_seed: u64) -> u64 {
    assert!(ARMS.iter().any(|&(arm, _)| arm == tx));
    grid_seed.wrapping_add((tx as u64) * TREATMENT_SEED_STRIDE)
}

fn spatial_config() -> SpatialParams {
    SpatialParams {
        cell_size_um: CELL_SIZE_UM,
        ..Default::default()
    }
}

fn configuration(grid_seed: u64) -> serde_json::Value {
    let immune = SpatialImmuneConfig::for_2d();
    serde_json::json!({
        "grid_rows": GRID_SIZE,
        "grid_cols": GRID_SIZE,
        "cell_size_um": CELL_SIZE_UM,
        "tumor_radius_um": GRID_SIZE as f64 * TUMOR_RADIUS_FRACTION * CELL_SIZE_UM,
        "n_steps": N_STEPS,
        "seed": grid_seed,
        "immune_start_step": IMMUNE_START_STEP,
        "damp_kill_threshold": DAMP_KILL_THRESHOLD,
        "o2_lambda_um": ZONE_REF_LAMBDA,
        "sdt_o2_dependence": 0.0,
        "params": Params::default(),
        "spatial_params": spatial_config(),
        "immune_config": {
            "damp_per_lp": immune.damp_per_lp,
            "damp_diffusion_fraction": immune.damp_diffusion_fraction,
            "damp_clearance_rate": immune.damp_clearance_rate,
            "dc_activation_kd": immune.dc_activation_kd,
            "immune_kill_rate": immune.immune_kill_rate,
            "pd1_brake": immune.pd1_brake,
            "anti_pd1_efficacy": immune.anti_pd1_efficacy,
            "exhaustion_rate": immune.exhaustion_rate,
            "ferro_immunosuppression_strength": immune.ferro_immunosuppression_strength,
        },
        "rng": {
            "scheme": "additive_per_cell",
            "treatment_seed_stride": TREATMENT_SEED_STRIDE,
            "init_offset": 0,
            "biochem_offset": 500_000,
            "biochem_step_stride": 1_000_000,
            "immune_offset": 900_000_000,
            "immune_step_stride": 2_000_000,
        },
    })
}

/// The measured path has no configurable death channels, parameter overrides,
/// or treatment labels. Tests reduce only grid size to check passivity cheaply.
fn run_condition(
    tx: Treatment,
    grid_size: usize,
    grid_seed: u64,
    diagnostics: ImmuneDiagnostics<'_>,
) -> (ConditionResult, Vec<f64>) {
    let tx_name = ARMS
        .iter()
        .find(|&&(arm, _)| arm == tx)
        .expect("canonical arm")
        .1;
    let params = Params::default();
    let spatial = spatial_config();
    let immune = SpatialImmuneConfig::for_2d();
    let mut grid = TumorGrid::generate(grid_size, grid_size, CELL_SIZE_UM, grid_seed);
    let stromal_mask = stromal_adjacency_mask_2d(&grid);
    let stromal_adj_count = stromal_mask.iter().filter(|&&b| b).count();
    let supply: Vec<f64> = apply_o2_gradient(&mut grid, ZONE_REF_LAMBDA)
        .iter()
        .map(|&(_, _, factor)| factor)
        .collect();
    let (ferro, immune_kills, damp) = run_spatial_with_immune_impl(
        &mut grid,
        tx,
        &params,
        &spatial,
        &immune,
        None,
        None,
        condition_seed(tx, grid_seed),
        Some(&supply),
        0.0,
        diagnostics,
    );
    let census = grid.census();
    let (norm, transition, hypoxic) = zone_kill_rates(&grid, ZONE_REF_LAMBDA);
    let result = ConditionResult {
        treatment: tx_name.to_string(),
        o2_condition: "gradient_120um".to_string(),
        o2_lambda_um: Some(ZONE_REF_LAMBDA),
        immune_mode: "immune_on".to_string(),
        total_tumor: census.total_tumor,
        total_dead: census.total_dead,
        ferroptosis_kills: Some(ferro),
        immune_kills: Some(immune_kills),
        overall_kill_rate: census.total_dead as f64 / census.total_tumor.max(1) as f64,
        normoxic_kill_rate: norm,
        transition_kill_rate: transition,
        hypoxic_kill_rate: hypoxic,
        stromal_mode: Some("off".to_string()),
        stromal_adjacent_kill_rate: Some(stromal_adjacent_kill_rate_2d(&grid, &stromal_mask)),
        stromal_adjacent_count: Some(stromal_adj_count),
        stromal_gsh_boost: None,
        stromal_mufa_boost: None,
        ph_mode: None,
        ph_iron_sensitivity: None,
        ph_ion_trap_sensitivity: None,
        ph_edge: None,
        ph_core: None,
        ph_lambda_um: None,
    };
    (result, damp)
}

fn observer(n_cells: usize) -> Measurements {
    let params = Params::default();
    let immune = SpatialImmuneConfig::for_2d();
    Measurements::new(
        n_cells,
        params.post_death_steps,
        immune.damp_per_lp,
        IMMUNE_START_STEP,
        DAMP_KILL_THRESHOLD,
        Some(immune.dc_activation_kd),
    )
}

pub(super) fn run(output_dir: &Path) {
    run_seeded(output_dir, GRID_SEED, None);
}

pub(super) fn run_replicate(output_dir: &Path, block: u32) {
    run_seeded(output_dir, replicate_seed(block), Some(block));
}

fn run_seeded(output_dir: &Path, grid_seed: u64, block: Option<u32>) {
    let conditions: Vec<_> = ARMS
        .into_iter()
        .map(|(tx, name)| {
            let mut measurements = observer(GRID_SIZE * GRID_SIZE);
            let (result, damp) = run_condition(
                tx,
                GRID_SIZE,
                grid_seed,
                ImmuneDiagnostics {
                    measurements: Some(&mut measurements),
                    #[cfg(test)]
                    snapshots: None,
                },
            );
            let total: f64 = damp.iter().sum();
            let peak = damp.iter().copied().fold(0.0_f64, f64::max);
            measurements.validate_totals(
                result.ferroptosis_kills.expect("ferroptosis count"),
                result.immune_kills.expect("immune count"),
                result.total_dead,
                total,
            );
            serde_json::json!({
                "condition_name": format!("immune_{name}"),
                "seed": condition_seed(tx, grid_seed),
                "result": result,
                "final_damp": {"total": total, "peak": peak},
                "measurements": measurements,
            })
        })
        .collect();
    let mut output = serde_json::json!({
        "schema_version": 2,
        "simulator": "sim-tme",
        "dimension": 2,
        "config": configuration(grid_seed),
        "conditions": conditions,
    });
    if let Some(block) = block {
        output["replicate_block"] = block.into();
    }
    fs::create_dir_all(output_dir).expect("create immune measurement directory");
    let path = output_dir.join(if block.is_some() {
        "immune_replicate.json"
    } else {
        "immune_measurements.json"
    });
    fs::write(
        &path,
        serde_json::to_vec(&output).expect("serialize immune measurements"),
    )
    .expect("write immune measurements");
    eprintln!("Wrote {}", path.display());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn immune_measurements_preserve_all_arm_results_and_every_step() {
        let size = 20;
        for (tx, _) in ARMS {
            let mut baseline = Vec::new();
            let (original, original_damp) = run_condition(
                tx,
                size,
                GRID_SEED,
                ImmuneDiagnostics {
                    measurements: None,
                    snapshots: Some(&mut baseline),
                },
            );
            let mut snapshots = Vec::new();
            let mut ledger = observer(size * size);
            let (measured, measured_damp) = run_condition(
                tx,
                size,
                GRID_SEED,
                ImmuneDiagnostics {
                    measurements: Some(&mut ledger),
                    snapshots: Some(&mut snapshots),
                },
            );
            assert_eq!(baseline.len(), N_STEPS as usize);
            assert_eq!(
                baseline, snapshots,
                "all cell state and LP/DAMP bits at each step"
            );
            assert_eq!(
                serde_json::to_vec(&original).unwrap(),
                serde_json::to_vec(&measured).unwrap()
            );
            assert_eq!(
                original_damp, measured_damp,
                "terminal field is also passive"
            );
            ledger.validate_totals(
                measured.ferroptosis_kills.unwrap(),
                measured.immune_kills.unwrap(),
                measured.total_dead,
                measured_damp.iter().sum(),
            );
            if tx == Treatment::SDT {
                assert!(!ledger.ferroptotic_events.is_empty());
                assert!(!ledger.eligible_cells.is_empty());
            }
            let mut repeated = observer(size * size);
            run_condition(
                tx,
                size,
                GRID_SEED,
                ImmuneDiagnostics {
                    measurements: Some(&mut repeated),
                    ..Default::default()
                },
            );
            assert_eq!(
                serde_json::to_vec(&ledger).unwrap(),
                serde_json::to_vec(&repeated).unwrap()
            );
        }
    }

    #[test]
    fn immune_measurements_reject_arguments_and_every_ferro_override() {
        let valid = vec!["sim-tme".into(), "--immune-measurements".into()];
        validate_args(&valid, &[]);
        for variable in [
            "FERRO_SEED",
            "FERRO_PARAM_OVERRIDES",
            "FERRO_SENSITIZER_ION_TRAP",
            "FERRO_FUTURE_KNOB",
        ] {
            assert!(
                std::panic::catch_unwind(|| validate_args(&valid, &[variable.into()])).is_err()
            );
        }
        for extra in ["--sdt-o2-dependence", "0", "--immune-measurements"] {
            let mut args = valid.clone();
            args.push(extra.into());
            assert!(std::panic::catch_unwind(|| validate_args(&args, &[])).is_err());
        }
        assert!(std::panic::catch_unwind(|| condition_seed(Treatment::PDT, GRID_SEED)).is_err());
        assert_eq!(
            ARMS.map(|(tx, _)| condition_seed(tx, GRID_SEED)),
            [42, 10_000_042, 20_000_042]
        );
    }

    #[test]
    fn immune_measurements_configuration_matches_runtime_defaults() {
        let config = configuration(GRID_SEED);
        let frozen: serde_json::Value = serde_json::from_str(include_str!(
            "../../../scripts/immune_2d_measurement_config_v1.json"
        ))
        .unwrap();
        assert_eq!(
            config, frozen["config"],
            "complete frozen 2D schema contract"
        );
        assert_eq!(
            config["params"],
            serde_json::to_value(Params::default()).unwrap()
        );
        assert_eq!(
            config["spatial_params"],
            serde_json::to_value(spatial_config()).unwrap()
        );
        assert_eq!(config["grid_rows"], 500);
        assert_eq!(config["tumor_radius_um"], 4500.0);
        assert_eq!(config["immune_config"]["dc_activation_kd"], 50.0);
        assert_eq!(config["immune_config"]["damp_diffusion_fraction"], 0.08);
    }

    #[test]
    fn immune_replicate_roots_have_disjoint_additive_seed_namespaces() {
        // Bound every seed used for initialization, biochemistry and immune
        // killing across all three arms, cells and simulated steps. This does
        // not claim that the treatment arms are independent within a block.
        let maximum_treatment_offset = 20_000_000;
        let last_cell = (GRID_SIZE * GRID_SIZE - 1) as u64;
        let last_step = u64::from(N_STEPS - 1);
        let stream_span = maximum_treatment_offset
            + last_cell
            + (500_000 + last_step * 1_000_000).max(900_000_000 + last_step * 2_000_000);
        assert_eq!(stream_span, 1_278_249_999);
        let mut previous_end = 61 + stream_span; // All historical roots 42–61.
        for block in 1..=20 {
            let root = replicate_seed(block);
            assert_eq!(root, 42 + u64::from(block) * 4_294_967_296);
            assert!(previous_end < root, "replicate namespaces must not overlap");
            previous_end = root.checked_add(stream_span).expect("no u64 seed wrapping");
            assert_eq!(
                ARMS.map(|(tx, _)| condition_seed(tx, root)),
                [root, root + 10_000_000, root + 20_000_000]
            );
        }
        for invalid in [0, 21, u32::MAX] {
            assert!(std::panic::catch_unwind(|| replicate_seed(invalid)).is_err());
        }
    }

    #[test]
    fn immune_replicate_configuration_changes_only_the_declared_seed() {
        let frozen: serde_json::Value = serde_json::from_str(include_str!(
            "../../../scripts/immune_2d_measurement_config_v1.json"
        ))
        .unwrap();
        for block in 1..=20 {
            let mut expected = frozen["config"].clone();
            expected["seed"] = replicate_seed(block).into();
            assert_eq!(configuration(replicate_seed(block)), expected);
        }
        assert_eq!(configuration(GRID_SEED), frozen["config"]);
    }

    #[test]
    fn immune_replicate_arguments_reject_malformed_mixed_or_overridden_runs() {
        for block in 1..=20 {
            let args = vec![
                "sim-tme".into(),
                "--immune-replicate".into(),
                block.to_string(),
            ];
            assert!(replicate_requested(&args));
            assert_eq!(validate_replicate_args(&args, &[]), block);
        }
        for value in [
            "0",
            "21",
            "-1",
            "+1",
            "01",
            "1.0",
            "1e1",
            " 1",
            "1 ",
            "",
            "１",
            "4294967296",
        ] {
            let args = vec!["sim-tme".into(), "--immune-replicate".into(), value.into()];
            assert!(std::panic::catch_unwind(|| validate_replicate_args(&args, &[])).is_err());
        }
        for args in [
            vec!["sim-tme", "--immune-replicate"],
            vec!["sim-tme", "--immune-replicate=1"],
            vec!["sim-tme", "--immune-replicates", "1"],
            vec!["sim-tme", "--immune-replicate", "1", "2"],
            vec![
                "sim-tme",
                "--immune-replicate",
                "1",
                "--immune-measurements",
            ],
            vec![
                "sim-tme",
                "--immune-measurements",
                "--immune-replicate",
                "1",
            ],
            vec!["sim-tme", "--other", "--immune-replicate", "1"],
        ] {
            let args: Vec<String> = args.into_iter().map(str::to_owned).collect();
            assert!(
                replicate_requested(&args),
                "reject malformed requests before default simulation"
            );
            assert!(std::panic::catch_unwind(|| validate_replicate_args(&args, &[])).is_err());
            if args.iter().any(|arg| arg == "--immune-measurements") {
                assert!(std::panic::catch_unwind(|| validate_args(&args, &[])).is_err());
            }
        }
        let valid = vec!["sim-tme".into(), "--immune-replicate".into(), "1".into()];
        for variable in ["FERRO_SEED", "FERRO_PARAM_OVERRIDES", "FERRO_FUTURE_KNOB"] {
            assert!(std::panic::catch_unwind(|| validate_replicate_args(
                &valid,
                &[variable.into()]
            ))
            .is_err());
        }
        assert!(!replicate_requested(&["sim-tme".into()]));
        assert!(!replicate_requested(&[
            "sim-tme".into(),
            "--immune-measurements".into()
        ]));
    }

    #[test]
    fn immune_replicate_observer_preserves_every_arm_and_step_at_extreme_blocks() {
        let size = 20;
        for block in [1, 20] {
            for (tx, _) in ARMS {
                let root = replicate_seed(block);
                let mut baseline = Vec::new();
                let (original, original_damp) = run_condition(
                    tx,
                    size,
                    root,
                    ImmuneDiagnostics {
                        measurements: None,
                        snapshots: Some(&mut baseline),
                    },
                );
                let mut snapshots = Vec::new();
                let mut ledger = observer(size * size);
                let (measured, measured_damp) = run_condition(
                    tx,
                    size,
                    root,
                    ImmuneDiagnostics {
                        measurements: Some(&mut ledger),
                        snapshots: Some(&mut snapshots),
                    },
                );
                assert_eq!(baseline.len(), N_STEPS as usize);
                assert_eq!(
                    baseline, snapshots,
                    "block {block}, arm {tx:?}: every phase-state snapshot"
                );
                assert_eq!(
                    serde_json::to_vec(&original).unwrap(),
                    serde_json::to_vec(&measured).unwrap()
                );
                assert_eq!(original_damp, measured_damp);
                ledger.validate_totals(
                    measured.ferroptosis_kills.unwrap(),
                    measured.immune_kills.unwrap(),
                    measured.total_dead,
                    measured_damp.iter().sum(),
                );
                if tx == Treatment::SDT {
                    assert!(!ledger.ferroptotic_events.is_empty());
                    assert!(!ledger.eligible_cells.is_empty());
                }
            }
        }
    }
}
