//! Per-step trajectory capture for visualizing sim-tme-3d runs (#193).
//!
//! Activated only when the binary is invoked with `--snapshot`; the
//! default 24-condition matrix path doesn't touch any of this code,
//! preserving bit-identical output. Captured data is written as three
//! `.npy` files alongside a small JSON metadata sidecar, then consumed
//! by `scripts/render_tme_3d_trajectory.py` to produce an animated
//! axial-slice GIF/MP4.
//!
//! Memory budget at 60³ × 180 steps (216 000 cells/step × 180 =
//! 38.88 M elements per field):
//!
//! - `dead`: 38.88 M × 1 byte  ≈ 37 MB
//! - `damp`: 38.88 M × 4 bytes ≈ 148 MB (f64 source cast to f32)
//! - `lp`:   38.88 M × 4 bytes ≈ 148 MB (same)
//!
//! Total ≈ 333 MB held in RAM during the snapshot run, dropped at write
//! time (matches the on-disk size — the `.npy` payloads are uncompressed).

use std::path::Path;

use ferroptosis_core::grid::TumorGrid3D;
use serde::Serialize;

use crate::npy;

/// Schema version for `trajectory_meta.json`. Bump when the trajectory
/// output shape changes (file count, axis order, dtype).
pub const TRAJECTORY_SCHEMA_VERSION: u32 = 1;

/// In-memory buffers for per-step trajectory state. Flat layout:
/// `value(step, idx) = buf[step * n_cells + idx]`. The renderer
/// reshapes to `(n_steps, layers, rows, cols)` per the meta sidecar.
pub struct SnapshotBuffers {
    pub dead: Vec<u8>,
    pub damp: Vec<f32>,
    pub lp: Vec<f32>,
    grid_dim: usize,
    steps_captured: u32,
}

impl SnapshotBuffers {
    pub fn new(grid_dim: usize, n_steps: u32) -> Self {
        let n_cells = grid_dim.pow(3);
        let total = n_cells * n_steps as usize;
        Self {
            dead: Vec::with_capacity(total),
            damp: Vec::with_capacity(total),
            lp: Vec::with_capacity(total),
            grid_dim,
            steps_captured: 0,
        }
    }

    /// Append one step of state. Must be called exactly once per
    /// simulation step, in step order. `damp_field` length must match
    /// `grid.cells.len()`.
    pub fn capture_step(&mut self, grid: &TumorGrid3D, damp_field: &[f64]) {
        debug_assert_eq!(
            grid.cells.len(),
            damp_field.len(),
            "snapshot capture: grid/damp length mismatch"
        );
        for (idx, gc) in grid.cells.iter().enumerate() {
            self.dead.push(u8::from(gc.state.dead));
            self.damp.push(damp_field[idx] as f32);
            self.lp.push(gc.state.lp as f32);
        }
        self.steps_captured += 1;
    }

    /// Write the three trajectory `.npy` files to `output_dir`.
    /// Shape on disk: `(steps_captured, grid_dim, grid_dim, grid_dim)`.
    pub fn write(&self, output_dir: &Path) -> std::io::Result<()> {
        let shape = [
            self.steps_captured as usize,
            self.grid_dim,
            self.grid_dim,
            self.grid_dim,
        ];
        npy::write_u8_array(output_dir.join("trajectory_dead.npy"), &shape, &self.dead)?;
        npy::write_f32_array(output_dir.join("trajectory_damp.npy"), &shape, &self.damp)?;
        npy::write_f32_array(output_dir.join("trajectory_lp.npy"), &shape, &self.lp)?;
        Ok(())
    }

    pub fn steps_captured(&self) -> u32 {
        self.steps_captured
    }
}

/// Metadata sidecar (`trajectory_meta.json`) describing the run that
/// produced the `trajectory_*.npy` files. Lets the Python renderer
/// label the animation without inferring from the binary.
#[derive(Serialize)]
pub struct TrajectoryMeta {
    pub schema_version: u32,
    pub grid_dim: usize,
    pub cell_size_um: f64,
    pub tumor_radius_um: f64,
    pub n_steps: u32,
    /// Steps at which a drug dose is administered (#239). Empty for the
    /// steady-state `Constant` presets. The Python renderer draws a marker
    /// on these frames so the viewer can see death waves sync to doses.
    pub dose_steps: Vec<u32>,
    pub condition: TrajectoryCondition,
}

/// Subset of the run's `Condition` that's relevant for visualization
/// labels. Mirrors the shape of `sim-tme-3d`'s `ConditionResult` for
/// the immune/stromal/ph fields.
#[derive(Serialize)]
pub struct TrajectoryCondition {
    pub treatment: String,
    pub o2_condition: String,
    pub o2_lambda_um: Option<f64>,
    pub immune_mode: String,
    pub stromal_mode: Option<String>,
    pub ph_mode: Option<String>,
}
