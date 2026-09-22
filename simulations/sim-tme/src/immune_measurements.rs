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

/// Deterministic controlled transport/activation comparison, isolated from the
/// canonical biochemical and stochastic immune loops. Kept in this source file
/// so the historical archives' Rust source inventory remains unchanged.
pub(super) mod immune_controlled {
    use std::ffi::OsString;
    use std::fs;
    use std::path::Path;

    use ferroptosis_core::biochem::CellState;
    use ferroptosis_core::cell::{Cell, Phenotype, Treatment};
    use ferroptosis_core::grid::{GridCell, TumorGrid};
    use ferroptosis_core::immune_spatial::{dc_activation, diffuse_damp_2d_step};
    use ferroptosis_core::params::Params;
    use serde::Serialize;

    const FNV_OFFSET: u64 = 0xcbf29ce484222325;
    const FNV_PRIME: u64 = 0x100000001b3;

    #[derive(Clone, Debug, Serialize)]
    struct Config {
        grid_rows: usize,
        grid_cols: usize,
        cell_size_um: f64,
        n_steps: u32,
        immune_start_step: u32,
        damp_diffusion_fraction: f64,
        damp_clearance_rate: f64,
        diffusion_source_cutoff: f64,
        damp_kill_threshold: f64,
        dc_activation_kd: f64,
    }

    impl Config {
        fn frozen() -> Self {
            Self {
                grid_rows: 500,
                grid_cols: 500,
                cell_size_um: 20.0,
                n_steps: 180,
                immune_start_step: 60,
                damp_diffusion_fraction: 0.08,
                damp_clearance_rate: 0.03,
                diffusion_source_cutoff: 0.001,
                damp_kill_threshold: 0.01,
                dc_activation_kd: 50.0,
            }
        }

        fn validate(&self) {
            assert!(self.grid_rows > 0 && self.grid_cols > 0 && self.n_steps > 0);
            assert!(self.cell_size_um.is_finite() && self.cell_size_um > 0.0);
            assert!(self.immune_start_step < self.n_steps);
            assert!(self.damp_diffusion_fraction.is_finite());
            assert!(self.damp_diffusion_fraction >= 0.0);
            assert!(8.0 * self.damp_diffusion_fraction < 1.0);
            assert!(self.damp_clearance_rate.is_finite());
            assert!((0.0..=1.0).contains(&self.damp_clearance_rate));
            assert_eq!(self.diffusion_source_cutoff, 0.001);
            assert!(self.damp_kill_threshold.is_finite() && self.damp_kill_threshold >= 0.0);
            assert!(self.dc_activation_kd.is_finite() && self.dc_activation_kd > 0.0);
        }
    }

    #[derive(Clone, Debug, Serialize)]
    struct MaskConfig {
        region_start: usize,
        region_stop: usize,
        tile_size: usize,
        source_offsets: [usize; 2],
        reference_recipient_count: usize,
        reduced_recipient_modulus: usize,
        reduced_recipient_remainder: usize,
    }

    impl MaskConfig {
        fn frozen() -> Self {
            Self {
                region_start: 218,
                region_stop: 282,
                tile_size: 8,
                source_offsets: [2, 5],
                reference_recipient_count: 3840,
                reduced_recipient_modulus: 2,
                reduced_recipient_remainder: 0,
            }
        }
    }

    #[derive(Debug, Serialize)]
    struct Masks {
        reserved_source_indices: Vec<usize>,
        source_indices_64: Vec<usize>,
        source_indices_256: Vec<usize>,
        reference_recipient_indices: Vec<usize>,
        available_recipient_indices_960: Vec<usize>,
    }

    fn frozen_masks(cfg: &Config, mask: &MaskConfig) -> Masks {
        assert!(mask.region_start < mask.region_stop);
        assert!(mask.region_stop <= cfg.grid_rows && mask.region_stop <= cfg.grid_cols);
        assert_eq!((mask.region_stop - mask.region_start) % mask.tile_size, 0);
        let mut reserved = Vec::new();
        let mut sparse = Vec::new();
        let mut reference = Vec::new();
        let mut reduced = Vec::new();
        for r in mask.region_start..mask.region_stop {
            for c in mask.region_start..mask.region_stop {
                let lr = r - mask.region_start;
                let lc = c - mask.region_start;
                let idx = r * cfg.grid_cols + c;
                let source = mask.source_offsets.contains(&(lr % mask.tile_size))
                    && mask.source_offsets.contains(&(lc % mask.tile_size));
                if source {
                    reserved.push(idx);
                    if lr % mask.tile_size == mask.source_offsets[(lr / mask.tile_size) % 2]
                        && lc % mask.tile_size == mask.source_offsets[(lc / mask.tile_size) % 2]
                    {
                        sparse.push(idx);
                    }
                } else {
                    reference.push(idx);
                    if lr % mask.reduced_recipient_modulus == mask.reduced_recipient_remainder
                        && lc % mask.reduced_recipient_modulus == mask.reduced_recipient_remainder
                    {
                        reduced.push(idx);
                    }
                }
            }
        }
        assert_eq!(reserved.len(), 256);
        assert_eq!(sparse.len(), 64);
        assert_eq!(reference.len(), mask.reference_recipient_count);
        assert_eq!(reduced.len(), 960);
        Masks {
            source_indices_256: reserved.clone(),
            reserved_source_indices: reserved,
            source_indices_64: sparse,
            reference_recipient_indices: reference,
            available_recipient_indices_960: reduced,
        }
    }

    #[derive(Clone, Debug, Serialize)]
    struct Condition {
        condition_id: String,
        source_count: usize,
        total_damp: f64,
        per_source_damp: f64,
        release_step: u32,
        recipient_count: usize,
    }

    fn frozen_conditions() -> Vec<Condition> {
        let mut result = Vec::new();
        for source_count in [64, 256] {
            for total_damp in [1280.0, 5120.0] {
                for release_step in [30, 60] {
                    for recipient_count in [3840, 960] {
                        result.push(Condition {
                            condition_id: format!(
                                "n{source_count}_m{total_damp:.0}_t{release_step}_r{recipient_count}"
                            ),
                            source_count,
                            total_damp,
                            per_source_damp: total_damp / source_count as f64,
                            release_step,
                            recipient_count,
                        });
                    }
                }
            }
        }
        for recipient_count in [3840, 960] {
            result.push(Condition {
                condition_id: format!("zero_r{recipient_count}"),
                source_count: 64,
                total_damp: 0.0,
                per_source_damp: 0.0,
                release_step: 60,
                recipient_count,
            });
        }
        result
    }

    // Only the topology is read by the shared diffusion helper. Construct inert
    // placeholders explicitly rather than sampling a biological grid or RNG.
    fn topology(cfg: &Config) -> TumorGrid {
        let cell = Cell {
            iron: 0.0,
            gsh: 0.0,
            gpx4: 0.0,
            fsp1: 0.0,
            basal_ros: 0.0,
            lipid_unsat: 0.0,
            nrf2: 0.0,
            mufa_cap: None,
            mufa_rate: None,
        };
        let state =
            CellState::from_cell_with_ros(&cell, Treatment::Control, &Params::default(), 0.0);
        let placeholder = GridCell {
            cell,
            state,
            phenotype: Phenotype::Stromal,
            is_tumor: false,
            extra_iron: 0.0,
            lp_at_grace_end: 0.0,
            newly_dead: false,
        };
        TumorGrid {
            cells: vec![placeholder; cfg.grid_rows * cfg.grid_cols],
            rows: cfg.grid_rows,
            cols: cfg.grid_cols,
            cell_size_um: cfg.cell_size_um,
        }
    }

    #[derive(Clone, Debug, Default, Serialize)]
    struct Activation {
        activation_sum: f64,
        max_eligible_damp: Option<f64>,
        damp_ge_kd_opportunities: usize,
        damp_ge_9kd_opportunities: usize,
    }

    impl Activation {
        fn observe(&mut self, damp: f64, kd: f64) {
            self.activation_sum += dc_activation(damp, kd);
            self.max_eligible_damp = Some(self.max_eligible_damp.map_or(damp, |old| old.max(damp)));
            self.damp_ge_kd_opportunities += usize::from(damp >= kd);
            self.damp_ge_9kd_opportunities += usize::from(damp >= 9.0 * kd);
        }
    }

    #[derive(Debug, Serialize)]
    struct Step {
        step: u32,
        injected_damp: f64,
        field_mass: f64,
        available_cells: usize,
        eligible_cells: usize,
        local_damp_sum: f64,
        #[serde(flatten)]
        activation: Activation,
        field_fingerprint: String,
    }

    #[derive(Debug, Serialize)]
    struct Recipient {
        cell_index: usize,
        available: bool,
        opportunities: usize,
        local_damp_sum: f64,
        #[serde(flatten)]
        activation: Activation,
    }

    #[derive(Debug, Serialize)]
    struct Summary {
        activation_sum: f64,
        activation_sum_per_reference_recipient: f64,
        opportunities: usize,
        unique_eligible_cells: usize,
        max_eligible_damp: Option<f64>,
        damp_ge_kd_opportunities: usize,
        damp_ge_9kd_opportunities: usize,
        final_field_mass: f64,
    }

    #[derive(Debug, Serialize)]
    struct Observation {
        #[serde(flatten)]
        condition: Condition,
        steps: Vec<Step>,
        recipients: Vec<Recipient>,
        summary: Summary,
        field_fingerprint: String,
    }

    fn fnv_bytes(mut hash: u64, bytes: &[u8]) -> u64 {
        for &byte in bytes {
            hash ^= u64::from(byte);
            hash = hash.wrapping_mul(FNV_PRIME);
        }
        hash
    }

    fn field_fingerprint(field: &[f64]) -> u64 {
        field.iter().fold(FNV_OFFSET, |hash, value| {
            assert!(
                value.is_finite() && *value >= 0.0,
                "finite nonnegative DAMP field"
            );
            fnv_bytes(hash, &value.to_bits().to_le_bytes())
        })
    }

    fn valid_indices(indices: &[usize], n: usize) {
        assert!(indices.iter().all(|&i| i < n), "mask index out of bounds");
        assert!(
            indices.windows(2).all(|w| w[0] < w[1]),
            "mask must be sorted and unique"
        );
    }

    fn run_condition(
        cfg: &Config,
        grid: &TumorGrid,
        condition: Condition,
        sources: &[usize],
        reference: &[usize],
        available: &[usize],
    ) -> Observation {
        cfg.validate();
        let n = cfg.grid_rows * cfg.grid_cols;
        assert_eq!((grid.rows, grid.cols), (cfg.grid_rows, cfg.grid_cols));
        assert_eq!(sources.len(), condition.source_count);
        assert!(condition.source_count > 0 && !reference.is_empty());
        assert_eq!(available.len(), condition.recipient_count);
        assert!(condition.total_damp.is_finite() && condition.total_damp >= 0.0);
        assert_eq!(
            condition.per_source_damp,
            condition.total_damp / condition.source_count as f64
        );
        assert!(condition.release_step < cfg.n_steps);
        for indices in [sources, reference, available] {
            valid_indices(indices, n);
        }
        assert!(
            sources.iter().all(|i| reference.binary_search(i).is_err()),
            "source/recipient overlap"
        );
        assert!(
            available.iter().all(|i| reference.binary_search(i).is_ok()),
            "availability outside reference mask"
        );

        let mut recipients: Vec<_> = reference
            .iter()
            .map(|&cell_index| Recipient {
                cell_index,
                available: available.binary_search(&cell_index).is_ok(),
                opportunities: 0,
                local_damp_sum: 0.0,
                activation: Activation::default(),
            })
            .collect();
        let mut damp = vec![0.0; n];
        let mut scratch = vec![0.0; n];
        let mut steps = Vec::with_capacity(cfg.n_steps as usize);
        let mut trajectory_hash = FNV_OFFSET;
        for step in 0..cfg.n_steps {
            let mut injected_damp = 0.0;
            if step == condition.release_step {
                for &idx in sources {
                    damp[idx] += condition.per_source_damp;
                    injected_damp += condition.per_source_damp;
                }
            }
            diffuse_damp_2d_step(
                &mut damp,
                &mut scratch,
                grid,
                cfg.damp_diffusion_fraction,
                cfg.damp_clearance_rate,
            );
            let fingerprint = field_fingerprint(&damp);
            trajectory_hash = fnv_bytes(trajectory_hash, &fingerprint.to_le_bytes());
            let mut row = Step {
                step,
                injected_damp,
                field_mass: damp.iter().sum(),
                available_cells: available.len(),
                eligible_cells: 0,
                local_damp_sum: 0.0,
                activation: Activation::default(),
                field_fingerprint: format!("{fingerprint:016x}"),
            };
            if step >= cfg.immune_start_step {
                for recipient in &mut recipients {
                    let local = damp[recipient.cell_index];
                    if !recipient.available || local < cfg.damp_kill_threshold {
                        continue;
                    }
                    recipient.opportunities += 1;
                    recipient.local_damp_sum += local;
                    recipient.activation.observe(local, cfg.dc_activation_kd);
                    row.eligible_cells += 1;
                    row.local_damp_sum += local;
                    row.activation.observe(local, cfg.dc_activation_kd);
                }
            }
            steps.push(row);
        }
        let activation_sum = recipients
            .iter()
            .map(|r| r.activation.activation_sum)
            .sum::<f64>();
        let summary = Summary {
            activation_sum,
            activation_sum_per_reference_recipient: activation_sum / reference.len() as f64,
            opportunities: recipients.iter().map(|r| r.opportunities).sum(),
            unique_eligible_cells: recipients.iter().filter(|r| r.opportunities > 0).count(),
            max_eligible_damp: recipients
                .iter()
                .filter_map(|r| r.activation.max_eligible_damp)
                .reduce(f64::max),
            damp_ge_kd_opportunities: recipients
                .iter()
                .map(|r| r.activation.damp_ge_kd_opportunities)
                .sum(),
            damp_ge_9kd_opportunities: recipients
                .iter()
                .map(|r| r.activation.damp_ge_9kd_opportunities)
                .sum(),
            final_field_mass: steps.last().expect("positive horizon").field_mass,
        };
        Observation {
            condition,
            steps,
            recipients,
            summary,
            field_fingerprint: format!("{trajectory_hash:016x}"),
        }
    }

    pub(crate) fn requested(args: &[String]) -> bool {
        args.iter()
            .any(|arg| arg.starts_with("--immune-controlled"))
    }

    pub(crate) fn validate_args(args: &[String], ferro_env: &[OsString]) {
        assert!(args.len() == 2 && args[1] == "--immune-controlled", "--immune-controlled is a standalone frozen comparison; extra or mixed arguments are unsupported");
        assert!(
            ferro_env.is_empty(),
            "--immune-controlled requires every FERRO_* environment variable to be absent"
        );
    }

    pub(crate) fn run(output_dir: &Path) {
        let cfg = Config::frozen();
        cfg.validate();
        let mask_cfg = MaskConfig::frozen();
        let masks = frozen_masks(&cfg, &mask_cfg);
        let grid = topology(&cfg);
        let observations: Vec<_> = frozen_conditions()
            .into_iter()
            .map(|condition| {
                let sources = match condition.source_count {
                    64 => &masks.source_indices_64,
                    256 => &masks.source_indices_256,
                    _ => unreachable!("frozen source count"),
                };
                let available = match condition.recipient_count {
                    3840 => &masks.reference_recipient_indices,
                    960 => &masks.available_recipient_indices_960,
                    _ => unreachable!("frozen recipient count"),
                };
                eprintln!("Controlled condition: {}", condition.condition_id);
                run_condition(
                    &cfg,
                    &grid,
                    condition,
                    sources,
                    &masks.reference_recipient_indices,
                    available,
                )
            })
            .collect();
        let output = serde_json::json!({
            "schema_version": 1,
            "study": "immune-2d-controlled-source-recipient",
            "simulator": "sim-tme",
            "dimension": 2,
            "config": cfg,
            "mask_config": mask_cfg,
            "masks": masks,
            "units": {
                "damp": "model field units",
                "time": "model steps; no physical-time conversion",
                "cell_size": "micrometers",
                "activation": "dimensionless",
                "activation_sum_per_reference_recipient": "activation-steps per fixed reference recipient"
            },
            "field_fingerprint_definition": {
                "step": "FNV-1a-64 of row-major finite nonnegative f64 bit patterns, each encoded little-endian",
                "condition": "FNV-1a-64 of the ordered step fingerprint u64 values, each encoded little-endian",
                "cryptographic": false
            },
            "conditions": observations,
        });
        fs::create_dir_all(output_dir).expect("create controlled output directory");
        fs::write(
            output_dir.join("immune_controlled.json"),
            serde_json::to_vec_pretty(&output).expect("serialize controlled observations"),
        )
        .expect("write controlled observations");
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        fn toy() -> Config {
            Config {
                grid_rows: 5,
                grid_cols: 5,
                n_steps: 6,
                immune_start_step: 2,
                ..Config::frozen()
            }
        }

        fn condition(total: f64, release: u32, recipient_count: usize) -> Condition {
            Condition {
                condition_id: "toy".to_owned(),
                source_count: 1,
                total_damp: total,
                per_source_damp: total,
                release_step: release,
                recipient_count,
            }
        }

        fn close(a: f64, b: f64) {
            assert!((a - b).abs() <= 1e-10 * b.abs() + 1e-8, "{a} != {b}");
        }

        #[test]
        fn frozen_design_masks_and_conditions_are_structural_only() {
            let cfg = Config::frozen();
            cfg.validate();
            let masks = frozen_masks(&cfg, &MaskConfig::frozen());
            assert_eq!(masks.reserved_source_indices, masks.source_indices_256);
            for indices in [
                &masks.source_indices_64,
                &masks.source_indices_256,
                &masks.reference_recipient_indices,
                &masks.available_recipient_indices_960,
            ] {
                valid_indices(indices, 500 * 500);
            }
            assert!(masks
                .source_indices_64
                .iter()
                .all(|i| masks.source_indices_256.binary_search(i).is_ok()));
            assert!(masks
                .source_indices_256
                .iter()
                .all(|i| masks.reference_recipient_indices.binary_search(i).is_err()));
            assert!(masks
                .available_recipient_indices_960
                .iter()
                .all(|i| masks.reference_recipient_indices.binary_search(i).is_ok()));
            let conditions = frozen_conditions();
            assert_eq!(conditions.len(), 18);
            assert_eq!(conditions[0].condition_id, "n64_m1280_t30_r3840");
            assert_eq!(conditions[15].condition_id, "n256_m5120_t60_r960");
            assert_eq!(conditions[16].condition_id, "zero_r3840");
            assert_eq!(conditions[17].condition_id, "zero_r960");
            for c in conditions {
                assert_eq!(c.total_damp, c.source_count as f64 * c.per_source_damp);
                assert!(c.release_step < cfg.n_steps);
            }
        }

        #[test]
        fn controlled_arguments_reject_mixed_modes_and_overrides() {
            let good = vec!["sim-tme".into(), "--immune-controlled".into()];
            validate_args(&good, &[]);
            for arg in ["--immune-controlled=1", "--immune-controlled-extra"] {
                let args = vec!["sim-tme".into(), arg.into()];
                assert!(requested(&args));
                assert!(std::panic::catch_unwind(|| validate_args(&args, &[])).is_err());
            }
            for extra in [
                "--immune-measurements",
                "--immune-replicate",
                "--other",
                "1",
            ] {
                let mut args = good.clone();
                args.push(extra.into());
                assert!(std::panic::catch_unwind(|| validate_args(&args, &[])).is_err());
            }
            assert!(std::panic::catch_unwind(|| validate_args(
                &good,
                &[OsString::from("FERRO_ANYTHING")]
            ))
            .is_err());
        }

        #[test]
        fn toy_zero_release_has_zero_sums_and_undefined_maxima() {
            let cfg = toy();
            let grid = topology(&cfg);
            let out = run_condition(&cfg, &grid, condition(0.0, 2, 2), &[12], &[6, 7], &[6, 7]);
            assert_eq!(out.summary.activation_sum, 0.0);
            assert_eq!(out.summary.opportunities, 0);
            assert!(out.summary.max_eligible_damp.is_none());
            for step in out.steps {
                assert_eq!(step.field_mass, 0.0);
                assert_eq!(step.injected_damp, 0.0);
                assert_eq!(step.eligible_cells, 0);
                assert!(step.activation.max_eligible_damp.is_none());
            }
        }

        #[test]
        fn toy_mass_timing_reconciliation_and_repeatability() {
            let cfg = toy();
            let grid = topology(&cfg);
            let run = || {
                run_condition(
                    &cfg,
                    &grid,
                    condition(50.0, 2, 3),
                    &[12],
                    &[6, 7, 8],
                    &[6, 7, 8],
                )
            };
            let out = run();
            assert_eq!(
                serde_json::to_vec(&out).unwrap(),
                serde_json::to_vec(&run()).unwrap()
            );
            let mut expected_mass = 0.0;
            for step in &out.steps {
                assert_eq!(step.injected_damp, if step.step == 2 { 50.0 } else { 0.0 });
                expected_mass = 0.97 * (expected_mass + step.injected_damp);
                close(step.field_mass, expected_mass);
                if step.step < 2 {
                    assert_eq!(step.eligible_cells, 0);
                    assert_eq!(step.activation.activation_sum, 0.0);
                }
            }
            close(
                out.summary.activation_sum,
                out.steps.iter().map(|s| s.activation.activation_sum).sum(),
            );
            close(
                out.recipients.iter().map(|r| r.local_damp_sum).sum(),
                out.steps.iter().map(|s| s.local_damp_sum).sum(),
            );
            assert_eq!(
                out.summary.opportunities,
                out.steps.iter().map(|s| s.eligible_cells).sum::<usize>()
            );
            close(
                out.summary.activation_sum_per_reference_recipient,
                out.summary.activation_sum / 3.0,
            );
        }

        #[test]
        fn toy_availability_changes_only_observation_not_transport() {
            let cfg = toy();
            let grid = topology(&cfg);
            let reference = &[6, 7, 8];
            let full = run_condition(
                &cfg,
                &grid,
                condition(50.0, 0, 3),
                &[12],
                reference,
                reference,
            );
            let reduced = run_condition(&cfg, &grid, condition(50.0, 0, 1), &[12], reference, &[7]);
            let empty = run_condition(&cfg, &grid, condition(50.0, 0, 0), &[12], reference, &[]);
            assert_eq!(full.field_fingerprint, reduced.field_fingerprint);
            assert_eq!(full.field_fingerprint, empty.field_fingerprint);
            for (f, r) in full.steps.iter().zip(&reduced.steps) {
                assert_eq!(f.field_mass.to_bits(), r.field_mass.to_bits());
                assert_eq!(f.field_fingerprint, r.field_fingerprint);
                if f.step < 2 {
                    assert_eq!(f.eligible_cells, 0);
                }
            }
            assert_eq!(
                serde_json::to_vec(&full.recipients[1]).unwrap(),
                serde_json::to_vec(&reduced.recipients[1]).unwrap()
            );
            assert_eq!(reduced.recipients[0].opportunities, 0);
            assert!(reduced.recipients[0].activation.max_eligible_damp.is_none());
            assert_eq!(empty.summary.activation_sum, 0.0);
            assert_eq!(empty.summary.opportunities, 0);
        }

        #[test]
        fn toy_matched_mass_with_different_sources_has_same_mass_history() {
            let cfg = toy();
            let grid = topology(&cfg);
            let one = run_condition(&cfg, &grid, condition(40.0, 1, 1), &[12], &[6], &[6]);
            let mut four = condition(40.0, 1, 1);
            four.source_count = 4;
            four.per_source_damp = 10.0;
            let four = run_condition(&cfg, &grid, four, &[11, 12, 13, 17], &[6], &[6]);
            for (a, b) in one.steps.iter().zip(&four.steps) {
                close(a.field_mass, b.field_mass);
            }
        }

        #[test]
        fn toy_eligibility_includes_exact_threshold_after_transport() {
            let mut cfg = toy();
            cfg.damp_clearance_rate = 0.0;
            cfg.immune_start_step = 0;
            cfg.n_steps = 1;
            let grid = topology(&cfg);
            let exact = run_condition(&cfg, &grid, condition(0.125, 0, 1), &[12], &[6], &[6]);
            assert_eq!(exact.steps[0].eligible_cells, 1);
            assert_eq!(exact.recipients[0].local_damp_sum, cfg.damp_kill_threshold);
            let below = run_condition(&cfg, &grid, condition(0.0625, 0, 1), &[12], &[6], &[6]);
            assert_eq!(below.steps[0].eligible_cells, 0);
            assert_eq!(dc_activation(0.0, 50.0), 0.0);
            assert_eq!(dc_activation(50.0, 50.0), 0.5);
            assert_eq!(dc_activation(450.0, 50.0), 0.9);
        }

        #[test]
        fn toy_invalid_masks_and_parameters_are_rejected() {
            let cfg = toy();
            let grid = topology(&cfg);
            assert!(std::panic::catch_unwind(|| run_condition(
                &cfg,
                &grid,
                condition(1.0, 0, 1),
                &[12],
                &[12],
                &[12]
            ))
            .is_err());
            assert!(std::panic::catch_unwind(|| run_condition(
                &cfg,
                &grid,
                condition(f64::NAN, 0, 1),
                &[12],
                &[6],
                &[6]
            ))
            .is_err());
            assert!(std::panic::catch_unwind(|| run_condition(
                &cfg,
                &grid,
                condition(1.0, 6, 1),
                &[12],
                &[6],
                &[6]
            ))
            .is_err());
            assert!(std::panic::catch_unwind(|| valid_indices(&[6, 6], 25)).is_err());
            assert!(std::panic::catch_unwind(|| valid_indices(&[25], 25)).is_err());
        }
    }
}

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
