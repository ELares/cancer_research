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
pub fn write_window_csv(path: &Path, results: &[(f64, String, f64, f64, f64)]) -> io::Result<()> {
    let mut wtr = csv::Writer::from_path(path)?;
    wtr.write_record([
        "timepoint_hours",
        "treatment",
        "death_rate",
        "ci_low",
        "ci_high",
    ])?;
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
    let json = serde_json::to_string_pretty(data).map_err(io::Error::other)?;
    std::fs::write(path, json)
}

#[cfg(test)]
mod tests {
    //! #302: these writers feed the FIGURES.yaml-traced figure generators, so a
    //! silent header or float-precision drift would corrupt a downstream figure
    //! with nothing catching it. These tests pin the header row and the
    //! `{:.1}` / `{:.6}` formats. (`str::lines()` strips a trailing `\r`, so the
    //! exact-string asserts are robust to the csv crate's line terminator.)
    use super::*;
    use ndarray::arr2;
    use std::fs;

    /// Unique temp path per test (process id + a per-test tag) so parallel tests
    /// and repeat runs never collide.
    fn tmp_path(tag: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!("ferro_io_{}_{}.tmp", std::process::id(), tag))
    }

    #[test]
    fn heatmap_csv_is_headerless_rows_of_values() {
        let path = tmp_path("heatmap");
        let data = arr2(&[[0u8, 1, 2], [3, 4, 255]]);
        write_heatmap_csv(&path, &data).unwrap();
        let body = fs::read_to_string(&path).unwrap();
        let _ = fs::remove_file(&path);
        let lines: Vec<&str> = body.lines().collect();
        assert_eq!(lines.len(), 2, "one line per row, no header");
        assert_eq!(lines[0], "0,1,2");
        assert_eq!(lines[1], "3,4,255");
    }

    #[test]
    fn depth_curves_csv_header_and_precision() {
        let path = tmp_path("depth");
        // rate rounds to 6 dp, depth formats to 1 dp, n is raw.
        let curves = vec![(
            "SDT".to_string(),
            vec![(100.0_f64, 0.123_456_789_f64, 7usize)],
        )];
        write_depth_curves_csv(&path, &curves).unwrap();
        let body = fs::read_to_string(&path).unwrap();
        let _ = fs::remove_file(&path);
        let lines: Vec<&str> = body.lines().collect();
        assert_eq!(lines[0], "depth_um,death_rate,n_cells,treatment");
        assert_eq!(lines[1], "100.0,0.123457,7,SDT");
    }

    #[test]
    fn window_csv_header_and_precision() {
        let path = tmp_path("window");
        let results = vec![(24.0_f64, "RSL3".to_string(), 0.5_f64, 0.4_f64, 0.6_f64)];
        write_window_csv(&path, &results).unwrap();
        let body = fs::read_to_string(&path).unwrap();
        let _ = fs::remove_file(&path);
        let lines: Vec<&str> = body.lines().collect();
        assert_eq!(
            lines[0],
            "timepoint_hours,treatment,death_rate,ci_low,ci_high"
        );
        assert_eq!(lines[1], "24.0,RSL3,0.500000,0.400000,0.600000");
    }

    #[test]
    fn json_is_pretty_and_round_trips() {
        let path = tmp_path("json");
        let data = vec![1u32, 2, 3];
        write_json(&path, &data).unwrap();
        let body = fs::read_to_string(&path).unwrap();
        let _ = fs::remove_file(&path);
        assert!(
            body.contains('\n'),
            "pretty-printed JSON should have newlines"
        );
        let back: Vec<u32> = serde_json::from_str(&body).unwrap();
        assert_eq!(back, data, "JSON round-trips");
    }
}
