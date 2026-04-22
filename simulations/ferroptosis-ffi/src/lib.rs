//! C FFI bindings for ferroptosis-core.
//!
//! Provides a C-compatible API for embedding the ferroptosis biochemistry
//! engine in C/C++ frameworks such as PhysiCell, CompuCell3D, or Chaste.
//!
//! ## Quick start (C)
//! ```c
//! #include "ferroptosis.h"
//!
//! FerroRng* rng = ferro_rng_new(42);
//! FerroParams params = ferro_params_default();
//! FerroCell cell = ferro_gen_cell(2, rng);  // Persister
//! FerroResult res = ferro_sim_cell(&cell, 1, &params, rng);  // RSL3
//! printf("Dead: %d, LP: %.2f\n", res.dead, res.final_lp);
//! ferro_rng_free(rng);
//! ```
//!
//! ## Enum values
//! - Phenotype: 0=Glycolytic, 1=OXPHOS, 2=Persister, 3=PersisterNrf2, 4=Stromal
//! - Treatment: 0=Control, 1=RSL3, 2=SDT, 3=PDT
//!
//! ## Safety
//! - Each thread must own its own FerroRng (not shared across threads)
//! - Caller must free FerroRng with ferro_rng_free (no double-free)
//! - Pointers must be valid and non-null (undefined behavior otherwise)

use ferroptosis_core::biochem::sim_cell;
use ferroptosis_core::cell::{gen_cell, Cell, Phenotype, Treatment};
use ferroptosis_core::params::Params;
use rand::rngs::StdRng;
use rand::SeedableRng;

// ============================================================
// C-compatible types
// ============================================================

/// Opaque random number generator. Create with ferro_rng_new, free with ferro_rng_free.
pub struct FerroRng {
    _private: [u8; 0],
}

/// Cell biochemical parameters (7 fields, all f64).
#[repr(C)]
pub struct FerroCell {
    pub iron: f64,
    pub gsh: f64,
    pub gpx4: f64,
    pub fsp1: f64,
    pub basal_ros: f64,
    pub lipid_unsat: f64,
    pub nrf2: f64,
}

/// Simulation parameters (20 f64 + 1 u32 = 21 fields).
/// Use ferro_params_default() or ferro_params_invivo() to create.
#[repr(C)]
pub struct FerroParams {
    pub fenton_rate: f64,
    pub gsh_scav_efficiency: f64,
    pub gsh_km: f64,
    pub nrf2_gsh_rate: f64,
    pub lp_rate: f64,
    pub lp_propagation: f64,
    pub gpx4_rate: f64,
    pub fsp1_rate: f64,
    pub scd_mufa_rate: f64,
    pub scd_mufa_max: f64,
    pub initial_mufa_protection: f64,
    pub scd_mufa_decay: f64,
    pub gpx4_degradation_by_ros: f64,
    pub gpx4_nrf2_upregulation: f64,
    pub sdt_ros: f64,
    pub pdt_ros: f64,
    pub rsl3_gpx4_inhib: f64,
    pub gsh_max: f64,
    pub gpx4_nrf2_target_multiplier: f64,
    pub death_threshold: f64,
    pub post_death_steps: u32,
}

/// Result of a single-cell ferroptosis simulation.
#[repr(C)]
pub struct FerroResult {
    pub dead: bool,
    pub final_lp: f64,
    pub final_gsh: f64,
    pub final_gpx4: f64,
}

// ============================================================
// Type conversions (internal)
// ============================================================

fn ferro_cell_to_cell(fc: &FerroCell) -> Cell {
    Cell {
        iron: fc.iron,
        gsh: fc.gsh,
        gpx4: fc.gpx4,
        fsp1: fc.fsp1,
        basal_ros: fc.basal_ros,
        lipid_unsat: fc.lipid_unsat,
        nrf2: fc.nrf2,
    }
}

fn cell_to_ferro_cell(c: &Cell) -> FerroCell {
    FerroCell {
        iron: c.iron,
        gsh: c.gsh,
        gpx4: c.gpx4,
        fsp1: c.fsp1,
        basal_ros: c.basal_ros,
        lipid_unsat: c.lipid_unsat,
        nrf2: c.nrf2,
    }
}

fn ferro_params_to_params(fp: &FerroParams) -> Params {
    Params {
        fenton_rate: fp.fenton_rate,
        gsh_scav_efficiency: fp.gsh_scav_efficiency,
        gsh_km: fp.gsh_km,
        nrf2_gsh_rate: fp.nrf2_gsh_rate,
        lp_rate: fp.lp_rate,
        lp_propagation: fp.lp_propagation,
        gpx4_rate: fp.gpx4_rate,
        fsp1_rate: fp.fsp1_rate,
        scd_mufa_rate: fp.scd_mufa_rate,
        scd_mufa_max: fp.scd_mufa_max,
        initial_mufa_protection: fp.initial_mufa_protection,
        scd_mufa_decay: fp.scd_mufa_decay,
        gpx4_degradation_by_ros: fp.gpx4_degradation_by_ros,
        gpx4_nrf2_upregulation: fp.gpx4_nrf2_upregulation,
        sdt_ros: fp.sdt_ros,
        pdt_ros: fp.pdt_ros,
        rsl3_gpx4_inhib: fp.rsl3_gpx4_inhib,
        gsh_max: fp.gsh_max,
        gpx4_nrf2_target_multiplier: fp.gpx4_nrf2_target_multiplier,
        death_threshold: fp.death_threshold,
        post_death_steps: fp.post_death_steps,
    }
}

fn params_to_ferro_params(p: &Params) -> FerroParams {
    FerroParams {
        fenton_rate: p.fenton_rate,
        gsh_scav_efficiency: p.gsh_scav_efficiency,
        gsh_km: p.gsh_km,
        nrf2_gsh_rate: p.nrf2_gsh_rate,
        lp_rate: p.lp_rate,
        lp_propagation: p.lp_propagation,
        gpx4_rate: p.gpx4_rate,
        fsp1_rate: p.fsp1_rate,
        scd_mufa_rate: p.scd_mufa_rate,
        scd_mufa_max: p.scd_mufa_max,
        initial_mufa_protection: p.initial_mufa_protection,
        scd_mufa_decay: p.scd_mufa_decay,
        gpx4_degradation_by_ros: p.gpx4_degradation_by_ros,
        gpx4_nrf2_upregulation: p.gpx4_nrf2_upregulation,
        sdt_ros: p.sdt_ros,
        pdt_ros: p.pdt_ros,
        rsl3_gpx4_inhib: p.rsl3_gpx4_inhib,
        gsh_max: p.gsh_max,
        gpx4_nrf2_target_multiplier: p.gpx4_nrf2_target_multiplier,
        death_threshold: p.death_threshold,
        post_death_steps: p.post_death_steps,
    }
}

fn i32_to_phenotype(v: i32) -> Phenotype {
    match v {
        0 => Phenotype::Glycolytic,
        1 => Phenotype::OXPHOS,
        2 => Phenotype::Persister,
        3 => Phenotype::PersisterNrf2,
        4 => Phenotype::Stromal,
        _ => Phenotype::Glycolytic, // default for invalid values
    }
}

fn i32_to_treatment(v: i32) -> Treatment {
    match v {
        0 => Treatment::Control,
        1 => Treatment::RSL3,
        2 => Treatment::SDT,
        3 => Treatment::PDT,
        _ => Treatment::Control, // default for invalid values
    }
}

// ============================================================
// FFI functions
// ============================================================

/// Create a new random number generator seeded with the given value.
/// The caller owns the returned pointer and MUST free it with ferro_rng_free.
/// Each thread should own its own FerroRng instance.
#[no_mangle]
pub extern "C" fn ferro_rng_new(seed: u64) -> *mut FerroRng {
    let rng = Box::new(StdRng::seed_from_u64(seed));
    Box::into_raw(rng) as *mut FerroRng
}

/// Free a FerroRng created by ferro_rng_new. Passing NULL is safe (no-op).
/// Do NOT double-free or use after free.
#[no_mangle]
pub extern "C" fn ferro_rng_free(rng: *mut FerroRng) {
    if !rng.is_null() {
        let _ = unsafe { Box::from_raw(rng as *mut StdRng) };
    }
}

/// Create default simulation parameters (2D culture).
/// Returns a FerroParams struct by value (all fields populated).
#[no_mangle]
pub extern "C" fn ferro_params_default() -> FerroParams {
    params_to_ferro_params(&Params::default())
}

/// Create in-vivo simulation parameters (3D/in-vivo with SCD1/MUFA protection).
/// Returns a FerroParams struct by value.
#[no_mangle]
pub extern "C" fn ferro_params_invivo() -> FerroParams {
    params_to_ferro_params(&Params::invivo())
}

/// Generate a stochastic cell with phenotype-specific biochemical parameters.
///
/// phenotype: 0=Glycolytic, 1=OXPHOS, 2=Persister, 3=PersisterNrf2, 4=Stromal.
/// Invalid values default to Glycolytic.
///
/// rng: Must be a valid FerroRng pointer (from ferro_rng_new). Must not be NULL.
#[no_mangle]
pub extern "C" fn ferro_gen_cell(phenotype: i32, rng: *mut FerroRng) -> FerroCell {
    if rng.is_null() {
        // Return zeroed cell on null RNG (defensive)
        return FerroCell {
            iron: 0.0, gsh: 0.0, gpx4: 0.0, fsp1: 0.0,
            basal_ros: 0.0, lipid_unsat: 0.0, nrf2: 0.0,
        };
    }
    let rng_ref = unsafe { &mut *(rng as *mut StdRng) };
    let pheno = i32_to_phenotype(phenotype);
    let cell = gen_cell(pheno, rng_ref);
    cell_to_ferro_cell(&cell)
}

/// Run a full 180-step ferroptosis simulation for one cell.
///
/// cell: Pointer to a FerroCell (from ferro_gen_cell). Must not be NULL.
/// treatment: 0=Control, 1=RSL3, 2=SDT, 3=PDT. Invalid values default to Control.
/// params: Pointer to FerroParams (from ferro_params_default/invivo). Must not be NULL.
/// rng: Must be a valid FerroRng pointer. Must not be NULL.
///
/// Returns a FerroResult with dead status and final LP/GSH/GPX4 values.
#[no_mangle]
pub extern "C" fn ferro_sim_cell(
    cell: *const FerroCell,
    treatment: i32,
    params: *const FerroParams,
    rng: *mut FerroRng,
) -> FerroResult {
    if cell.is_null() || params.is_null() || rng.is_null() {
        return FerroResult {
            dead: false, final_lp: 0.0, final_gsh: 0.0, final_gpx4: 0.0,
        };
    }

    let cell_ref = unsafe { &*cell };
    let params_ref = unsafe { &*params };
    let rng_ref = unsafe { &mut *(rng as *mut StdRng) };

    let rust_cell = ferro_cell_to_cell(cell_ref);
    let rust_params = ferro_params_to_params(params_ref);
    let tx = i32_to_treatment(treatment);

    let (dead, lp, gsh, gpx4) = sim_cell(&rust_cell, tx, &rust_params, rng_ref);

    FerroResult {
        dead,
        final_lp: lp,
        final_gsh: gsh,
        final_gpx4: gpx4,
    }
}
