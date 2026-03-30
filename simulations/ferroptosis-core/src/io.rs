//! Output helpers: JSON, CSV.

use std::io;
use std::path::Path;

use ndarray::Array2;

/// Write a 2D u8 array as CSV (for heatmap plotting).
pub fn write_heatmap_csv(path: &Path, data: &Array2<u8>) -> io::Result<()> {
    let mut wtr = csv::WriterBuilder::new()
        .has_headers(false)
        .from_path(path)?;
    for row in data.rows() {
        let record: Vec<String> = row.iter().map(|v| v.to_string()).collect();
        wtr.write_record(&record)?;
    }
    wtr.flush()?;
    Ok(())
}

/// Write depth-kill curves as CSV.
/// Format: depth_um, death_rate, n_cells, treatment
pub fn write_depth_curves_csv(
    path: &Path,
    curves: &[(String, Vec<(f64, f64, usize)>)],
) -> io::Result<()> {
    let mut wtr = csv::Writer::from_path(path)?;
    wtr.write_record(["depth_um", "death_rate", "n_cells", "treatment"])?;
    for (treatment, points) in curves {
        for (depth, rate, n) in points {
            wtr.write_record(&[
                format!("{:.1}", depth),
                format!("{:.6}", rate),
                n.to_string(),
                treatment.clone(),
            ])?;
        }
    }
    wtr.flush()?;
    Ok(())
}

/// Write vulnerability window results as CSV.
/// Format: timepoint_hours, treatment, death_rate, ci_low, ci_high
pub fn write_window_csv(
    path: &Path,
    results: &[(f64, String, f64, f64, f64)],
) -> io::Result<()> {
    let mut wtr = csv::Writer::from_path(path)?;
    wtr.write_record(["timepoint_hours", "treatment", "death_rate", "ci_low", "ci_high"])?;
    for (hours, tx, rate, ci_lo, ci_hi) in results {
        wtr.write_record(&[
            format!("{:.1}", hours),
            tx.clone(),
            format!("{:.6}", rate),
            format!("{:.6}", ci_lo),
            format!("{:.6}", ci_hi),
        ])?;
    }
    wtr.flush()?;
    Ok(())
}

/// Write any Serialize value as pretty-printed JSON to a file.
pub fn write_json<T: serde::Serialize>(path: &Path, data: &T) -> io::Result<()> {
    let json = serde_json::to_string_pretty(data)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
    std::fs::write(path, json)
}
