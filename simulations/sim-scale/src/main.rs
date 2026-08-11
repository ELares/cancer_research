//! Large-n rare-event sweeps over the ferroptosis kill switch.
//!
//! WHY THIS EXISTS
//! ---------------
//! Several conditions in the manuscript's Figure 7 report a death rate of
//! exactly 0% at n = 1,000,000. That is not a measurement of zero: the per-cell
//! parameters are drawn from NORMAL distributions with unbounded upper support
//! (`gen_cell`), so every threshold retains a positive-probability tail. A
//! reported 0% is an upper bound set by the sample size, roughly 3e-6 at a
//! million cells.
//!
//! This resolves those bounds by sweeping n upward. The interesting output is
//! not any single number but the SHAPE: a rate that keeps falling as n grows is
//! still resolution-limited, while one that settles is a real estimate.
//!
//! WHY THE SCALES ARE WHAT THEY ARE
//! ---------------------------------
//! They are tumor burdens, not round numbers.
//!   1e9   about one gram, the smallest clinically detectable lesion;
//!   1e10  an intermediate burden;
//!   1e11  about a hundred grams, advanced metastatic disease.
//! 1e12 upward is typically lethal and is deliberately not swept. At a survival
//! rate of 1e-9 a 1e11-cell burden still harbours about a hundred survivors,
//! which is the resistance-escape question with a clinical denominator behind
//! it.
//!
//! WHAT IT IS NOT
//! --------------
//! 1e11 INDEPENDENT cells is not a 100-gram tumor. These cells never see each
//! other: no diffusion, no crowding, no vasculature, no clonal competition. The
//! defensible phrase is "a population the size of an advanced burden, sampled as
//! independent cells". The spatial model is the one with interaction, and it
//! does not reach these counts.
//!
//! PARAMETERS
//! ----------
//! `--params` takes the same JSON override map as `FERRO_PARAM_OVERRIDES`
//! (#331), so a run can be driven from the committed CTRPv2 in-vitro posterior
//! rather than from `Params::default()`. That matters: the default set is the
//! in-vivo one and `CALIBRATION_STATUS.md` records zero calibration targets for
//! it, while the in-vitro posterior is the repo's only data-anchored
//! parameterisation. The two are provably disjoint (#332, #500), so a run must
//! pick one and say which.
//!
//! SELF-CHECK
//! ----------
//! `--verify` reruns the committed Figure 7 conditions at n = 1e6 with default
//! parameters and compares against the values the figure reports. A large run
//! refuses to start unless that passes, because an engine that cannot reproduce
//! a known result should not be trusted with an expensive one.

use clap::Parser;
use ferroptosis_core::cell::{Phenotype, Treatment};
use ferroptosis_core::params::{apply_param_overrides, parse_param_overrides_json, Params};
use ferroptosis_core::stats::run_condition;
use std::time::Instant;

/// Conditions the manuscript's Figure 7 reports, with the death rates it
/// reports for them. Used by `--verify` as a fixed point: these come from the
/// committed figure, not from this binary.
const FIGURE_7: &[(&str, &str, f64)] = &[
    ("Glycolytic", "Control", 0.0),
    ("Glycolytic", "SDT", 0.872),
    ("OXPHOS", "SDT", 0.999),
    ("Persister", "RSL3", 0.425),
    ("Persister", "Control", 0.012),
    ("PersisterNrf2", "SDT", 0.995),
];

#[derive(Parser)]
#[command(about = "Large-n rare-event sweeps over the ferroptosis kill switch")]
struct Args {
    /// Cells per condition. Accepts 1e9 style shorthand.
    #[arg(long, default_value = "1e6")]
    cells: String,
    #[arg(long, default_value = "Persister")]
    phenotype: String,
    #[arg(long, default_value = "RSL3")]
    treatment: String,
    /// JSON map of parameter overrides, e.g. '{"lp_propagation":0.78}'
    #[arg(long)]
    params: Option<String>,
    /// A label carried into the output so a sweep's rows are self-describing.
    #[arg(long, default_value = "unlabelled")]
    label: String,
    /// Reproduce the committed Figure 7 values and exit.
    #[arg(long)]
    verify: bool,
}

fn parse_count(s: &str) -> Result<usize, String> {
    let t = s.trim().replace('_', "");
    if let Ok(v) = t.parse::<usize>() {
        return Ok(v);
    }
    // 1e9 / 2.5e10 shorthand: a sweep is unreadable written out in full.
    let f: f64 = t.parse().map_err(|_| format!("cannot parse count {s:?}"))?;
    if !f.is_finite() || f < 1.0 {
        return Err(format!("count {s:?} is not a positive number"));
    }
    Ok(f as usize)
}

fn parse_phenotype(s: &str) -> Result<Phenotype, String> {
    Ok(match s.to_ascii_lowercase().replace(['-', '_'], "").as_str() {
        "glycolytic" => Phenotype::Glycolytic,
        "oxphos" => Phenotype::OXPHOS,
        "persister" => Phenotype::Persister,
        "persisternrf2" => Phenotype::PersisterNrf2,
        _ => return Err(format!("unknown phenotype {s:?}")),
    })
}

fn parse_treatment(s: &str) -> Result<Treatment, String> {
    Ok(match s.to_ascii_lowercase().as_str() {
        "control" => Treatment::Control,
        "rsl3" => Treatment::RSL3,
        "sdt" => Treatment::SDT,
        "pdt" => Treatment::PDT,
        _ => return Err(format!("unknown treatment {s:?}")),
    })
}

/// Upper bound on a rate when zero events are observed, at 95%.
/// The rule of three: with 0 successes in n trials the one-sided 95% bound is
/// about 3/n. Reported because "0%" alone is the number that misleads.
fn zero_event_bound(n: usize) -> f64 {
    3.0 / n as f64
}

fn verify() -> i32 {
    println!("self-check: reproducing the committed Figure 7 at n=1e6, default params");
    let params = Params::default();
    let mut worst: f64 = 0.0;
    let mut failed = 0;
    for (p, t, expect) in FIGURE_7 {
        let r = run_condition(
            parse_phenotype(p).unwrap(),
            parse_treatment(t).unwrap(),
            &params,
            1_000_000,
            p,
            t,
        );
        // The figure prints to one decimal in percent, so agreement is checked
        // at that resolution rather than to the bit.
        let d = (r.death_rate - expect).abs();
        worst = worst.max(d);
        let ok = d <= 0.0006;
        if !ok {
            failed += 1;
        }
        println!(
            "  {:<14} {:<8} got {:>7.4}  figure {:>6.3}  {}",
            p,
            t,
            r.death_rate,
            expect,
            if ok { "ok" } else { "MISMATCH" }
        );
    }
    println!("  worst absolute difference: {worst:.2e}");
    if failed > 0 {
        eprintln!("self-check FAILED on {failed} condition(s); refusing to vouch for a large run");
        return 1;
    }
    println!("self-check passed");
    0
}

fn main() -> std::process::ExitCode {
    let args = Args::parse();
    if args.verify {
        return std::process::ExitCode::from(verify() as u8);
    }

    let n = match parse_count(&args.cells) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("{e}");
            return std::process::ExitCode::from(2);
        }
    };
    let pheno = parse_phenotype(&args.phenotype).unwrap();
    let tx = parse_treatment(&args.treatment).unwrap();

    let mut params = Params::default();
    let mut applied = serde_json::Map::new();
    if let Some(js) = &args.params {
        let ov = match parse_param_overrides_json(js) {
            Ok(v) => v,
            Err(e) => {
                eprintln!("bad --params: {e}");
                return std::process::ExitCode::from(2);
            }
        };
        for (k, v) in &ov {
            applied.insert(k.clone(), serde_json::json!(v));
        }
        if let Err(e) = apply_param_overrides(&mut params, ov) {
            eprintln!("bad --params: {e}");
            return std::process::ExitCode::from(2);
        }
    }

    let t0 = Instant::now();
    let r = run_condition(pheno, tx, &params, n, &args.phenotype, &args.treatment);
    let secs = t0.elapsed().as_secs_f64();

    let out = serde_json::json!({
        "label": args.label,
        "phenotype": args.phenotype,
        "treatment": args.treatment,
        "n_cells": r.n_cells,
        "n_dead": r.n_dead,
        "death_rate": r.death_rate,
        "ci_low": r.ci_low,
        "ci_high": r.ci_high,
        // The number a "0%" row actually licenses. Reported always, so a reader
        // never has to work out whether zero means zero.
        "zero_event_upper_bound_95": zero_event_bound(n),
        "observed_zero": r.n_dead == 0,
        "mean_lipid_perox": r.mean_lipid_perox,
        "mean_gsh_final": r.mean_gsh_final,
        "mean_gpx4_final": r.mean_gpx4_final,
        "param_overrides": applied,
        "params_source": if applied.is_empty() { "Params::default() (IN-VIVO, uncalibrated)" }
                         else { "override map supplied" },
        "wall_seconds": secs,
        "cells_per_cpu_second": n as f64 / secs.max(1e-9),
    });
    println!("{}", serde_json::to_string(&out).unwrap());
    std::process::ExitCode::SUCCESS
}
