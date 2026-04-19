//! Python bindings for ferroptosis-core.
//!
//! Exposes a functional, dict-based API designed for Python researchers.
//! No classes to learn — just functions that take strings and keyword
//! arguments, returning dicts.
//!
//! ```python
//! import ferroptosis_core as fc
//!
//! stats = fc.sim_batch("Persister", "RSL3", n=1000, seed=42)
//! print(f"Death rate: {stats['death_rate']:.1%}")
//! ```

use std::collections::HashMap;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use rand::prelude::*;
use rayon::prelude::*;

use ::ferroptosis_core::biochem::sim_cell as rust_sim_cell;
use ::ferroptosis_core::cell::{gen_cell, Phenotype, Treatment};
use ::ferroptosis_core::params::Params;
use ::ferroptosis_core::stats::wilson_ci;

// ============================================================
// Enum parsing
// ============================================================

fn parse_phenotype(s: &str) -> PyResult<Phenotype> {
    match s {
        "Glycolytic" => Ok(Phenotype::Glycolytic),
        "OXPHOS" => Ok(Phenotype::OXPHOS),
        "Persister" => Ok(Phenotype::Persister),
        "PersisterNrf2" | "Persister+NRF2" => Ok(Phenotype::PersisterNrf2),
        "Stromal" => Ok(Phenotype::Stromal),
        _ => Err(PyValueError::new_err(format!(
            "Unknown phenotype '{}'. Valid: Glycolytic, OXPHOS, Persister, PersisterNrf2, Stromal",
            s
        ))),
    }
}

fn parse_treatment(s: &str) -> PyResult<Treatment> {
    match s {
        "Control" => Ok(Treatment::Control),
        "RSL3" => Ok(Treatment::RSL3),
        "SDT" => Ok(Treatment::SDT),
        "PDT" => Ok(Treatment::PDT),
        _ => Err(PyValueError::new_err(format!(
            "Unknown treatment '{}'. Valid: Control, RSL3, SDT, PDT",
            s
        ))),
    }
}

fn parse_context(s: &str) -> PyResult<Params> {
    match s {
        "2d" | "default" => Ok(Params::default()),
        "invivo" | "in-vivo" | "in_vivo" => Ok(Params::invivo()),
        _ => Err(PyValueError::new_err(format!(
            "Unknown context '{}'. Valid: 2d, invivo",
            s
        ))),
    }
}

// ============================================================
// Params helpers
// ============================================================

fn params_to_dict(py: Python<'_>, params: &Params) -> PyResult<Py<PyDict>> {
    let dict = PyDict::new(py);
    dict.set_item("fenton_rate", params.fenton_rate)?;
    dict.set_item("gsh_scav_efficiency", params.gsh_scav_efficiency)?;
    dict.set_item("gsh_km", params.gsh_km)?;
    dict.set_item("nrf2_gsh_rate", params.nrf2_gsh_rate)?;
    dict.set_item("lp_rate", params.lp_rate)?;
    dict.set_item("lp_propagation", params.lp_propagation)?;
    dict.set_item("gpx4_rate", params.gpx4_rate)?;
    dict.set_item("fsp1_rate", params.fsp1_rate)?;
    dict.set_item("scd_mufa_rate", params.scd_mufa_rate)?;
    dict.set_item("scd_mufa_max", params.scd_mufa_max)?;
    dict.set_item("initial_mufa_protection", params.initial_mufa_protection)?;
    dict.set_item("scd_mufa_decay", params.scd_mufa_decay)?;
    dict.set_item("gpx4_degradation_by_ros", params.gpx4_degradation_by_ros)?;
    dict.set_item("gpx4_nrf2_upregulation", params.gpx4_nrf2_upregulation)?;
    dict.set_item("sdt_ros", params.sdt_ros)?;
    dict.set_item("pdt_ros", params.pdt_ros)?;
    dict.set_item("rsl3_gpx4_inhib", params.rsl3_gpx4_inhib)?;
    dict.set_item("gsh_max", params.gsh_max)?;
    dict.set_item("gpx4_nrf2_target_multiplier", params.gpx4_nrf2_target_multiplier)?;
    dict.set_item("death_threshold", params.death_threshold)?;
    Ok(dict.into())
}

fn apply_overrides(params: &mut Params, overrides: &HashMap<String, f64>) -> PyResult<()> {
    for (key, val) in overrides {
        match key.as_str() {
            "fenton_rate" => params.fenton_rate = *val,
            "gsh_scav_efficiency" => params.gsh_scav_efficiency = *val,
            "gsh_km" => params.gsh_km = *val,
            "nrf2_gsh_rate" => params.nrf2_gsh_rate = *val,
            "lp_rate" => params.lp_rate = *val,
            "lp_propagation" => params.lp_propagation = *val,
            "gpx4_rate" => params.gpx4_rate = *val,
            "fsp1_rate" => params.fsp1_rate = *val,
            "scd_mufa_rate" => params.scd_mufa_rate = *val,
            "scd_mufa_max" => params.scd_mufa_max = *val,
            "initial_mufa_protection" => params.initial_mufa_protection = *val,
            "scd_mufa_decay" => params.scd_mufa_decay = *val,
            "gpx4_degradation_by_ros" => params.gpx4_degradation_by_ros = *val,
            "gpx4_nrf2_upregulation" => params.gpx4_nrf2_upregulation = *val,
            "sdt_ros" => params.sdt_ros = *val,
            "pdt_ros" => params.pdt_ros = *val,
            "rsl3_gpx4_inhib" => params.rsl3_gpx4_inhib = *val,
            "gsh_max" => params.gsh_max = *val,
            "gpx4_nrf2_target_multiplier" => params.gpx4_nrf2_target_multiplier = *val,
            "death_threshold" => params.death_threshold = *val,
            _ => {
                return Err(PyValueError::new_err(format!(
                    "Unknown parameter '{}'. Use default_params() to see valid names.",
                    key
                )));
            }
        }
    }
    Ok(())
}

// ============================================================
// Python API
// ============================================================

/// Return default parameter values as a dict.
#[pyfunction]
fn default_params(py: Python<'_>) -> PyResult<Py<PyDict>> {
    params_to_dict(py, &Params::default())
}

/// Return in-vivo parameter values (with SCD1/MUFA protection) as a dict.
#[pyfunction]
fn invivo_params(py: Python<'_>) -> PyResult<Py<PyDict>> {
    params_to_dict(py, &Params::invivo())
}

/// Simulate a single cell and return the outcome as a dict.
///
/// Args:
///     phenotype: "Glycolytic", "OXPHOS", "Persister", "PersisterNrf2", or "Stromal"
///     treatment: "Control", "RSL3", "SDT", or "PDT"
///     seed: RNG seed for reproducibility
///     context: "2d" (default) or "invivo" (enables SCD1/MUFA protection)
///     **kwargs: parameter overrides (e.g., rsl3_gpx4_inhib=0.5)
///
/// Returns:
///     dict with keys: dead (bool), lp, gsh, gpx4 (floats)
#[pyfunction]
#[pyo3(signature = (phenotype, treatment, seed, context="2d", **kwargs))]
fn sim_cell(
    py: Python<'_>,
    phenotype: &str,
    treatment: &str,
    seed: u64,
    context: &str,
    kwargs: Option<HashMap<String, f64>>,
) -> PyResult<Py<PyDict>> {
    let pheno = parse_phenotype(phenotype)?;
    let tx = parse_treatment(treatment)?;
    let mut params = parse_context(context)?;
    if let Some(overrides) = &kwargs {
        apply_overrides(&mut params, overrides)?;
    }

    let mut cell_rng = StdRng::seed_from_u64(seed);
    let cell = gen_cell(pheno, &mut cell_rng);
    let mut sim_rng = StdRng::seed_from_u64(seed.wrapping_add(1));
    let (dead, lp, gsh, gpx4) = rust_sim_cell(&cell, tx, &params, &mut sim_rng);

    let dict = PyDict::new(py);
    dict.set_item("dead", dead)?;
    dict.set_item("lp", lp)?;
    dict.set_item("gsh", gsh)?;
    dict.set_item("gpx4", gpx4)?;
    Ok(dict.into())
}

/// Simulate a population of cells and return aggregate statistics.
///
/// Runs n cells in parallel (via rayon) and returns death rate with
/// Wilson confidence intervals and mean final pathway states.
///
/// Args:
///     phenotype: "Glycolytic", "OXPHOS", "Persister", "PersisterNrf2", or "Stromal"
///     treatment: "Control", "RSL3", "SDT", or "PDT"
///     n: number of cells to simulate
///     seed: RNG seed for reproducibility
///     context: "2d" (default) or "invivo" (enables SCD1/MUFA protection)
///     **kwargs: parameter overrides (e.g., rsl3_gpx4_inhib=0.5)
///
/// Returns:
///     dict with keys: death_rate, ci_low, ci_high, n_dead, n_cells,
///                     mean_lp, mean_gsh, mean_gpx4
#[pyfunction]
#[pyo3(signature = (phenotype, treatment, n, seed, context="2d", **kwargs))]
fn sim_batch(
    py: Python<'_>,
    phenotype: &str,
    treatment: &str,
    n: usize,
    seed: u64,
    context: &str,
    kwargs: Option<HashMap<String, f64>>,
) -> PyResult<Py<PyDict>> {
    let pheno = parse_phenotype(phenotype)?;
    let tx = parse_treatment(treatment)?;
    if n == 0 {
        return Err(PyValueError::new_err("n must be > 0"));
    }
    let mut params = parse_context(context)?;
    if let Some(overrides) = &kwargs {
        apply_overrides(&mut params, overrides)?;
    }

    // Run in parallel, releasing the GIL so Python threads aren't blocked
    let results: Vec<(bool, f64, f64, f64)> = py.allow_threads(|| {
        (0..n)
            .into_par_iter()
            .map(|i| {
                let cell_seed = seed.wrapping_add((i as u64) * 2);
                let mut cell_rng = StdRng::seed_from_u64(cell_seed);
                let cell = gen_cell(pheno, &mut cell_rng);
                let mut sim_rng = StdRng::seed_from_u64(cell_seed.wrapping_add(1));
                rust_sim_cell(&cell, tx, &params, &mut sim_rng)
            })
            .collect()
    });

    let n_dead = results.iter().filter(|(dead, _, _, _)| *dead).count();
    let death_rate = n_dead as f64 / n as f64;
    let (ci_low, ci_high) = wilson_ci(n, n_dead);
    let mean_lp = results.iter().map(|(_, lp, _, _)| lp).sum::<f64>() / n as f64;
    let mean_gsh = results.iter().map(|(_, _, gsh, _)| gsh).sum::<f64>() / n as f64;
    let mean_gpx4 = results.iter().map(|(_, _, _, gpx4)| gpx4).sum::<f64>() / n as f64;

    let dict = PyDict::new(py);
    dict.set_item("death_rate", death_rate)?;
    dict.set_item("ci_low", ci_low)?;
    dict.set_item("ci_high", ci_high)?;
    dict.set_item("n_dead", n_dead)?;
    dict.set_item("n_cells", n)?;
    dict.set_item("mean_lp", mean_lp)?;
    dict.set_item("mean_gsh", mean_gsh)?;
    dict.set_item("mean_gpx4", mean_gpx4)?;
    Ok(dict.into())
}

// ============================================================
// Module registration
// ============================================================

#[pymodule]
fn ferroptosis_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(default_params, m)?)?;
    m.add_function(wrap_pyfunction!(invivo_params, m)?)?;
    m.add_function(wrap_pyfunction!(sim_cell, m)?)?;
    m.add_function(wrap_pyfunction!(sim_batch, m)?)?;
    Ok(())
}
