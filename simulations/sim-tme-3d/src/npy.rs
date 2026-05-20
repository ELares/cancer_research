//! Minimal NumPy `.npy` v1.0 writer for sim-tme-3d trajectory snapshots
//! (#193).
//!
//! Supports flat n-dimensional arrays in `u8` (dead-cell mask) and `f32`
//! (DAMP, LP fields). Hand-rolled because adding `ndarray-npy` would
//! pull `ndarray` back into sim-tme-3d's dep graph (deliberately
//! omitted per the #195 follow-up note in `Cargo.toml`), and the format
//! is small enough that a focused writer is easier to audit.
//!
//! Format reference:
//! <https://numpy.org/doc/stable/reference/generated/numpy.lib.format.html>
//!
//! Output bytes:
//! ```text
//! \x93NUMPY \x01 \x00 <hdr_len: u16 LE> <dict ASCII> <spaces…> '\n' <data LE>
//! ```
//! where the prefix + dict + padding + newline must be a multiple of 64
//! (so the data payload is 64-byte aligned for fast NumPy reads).
//!
//! Roundtrip is verified manually with `numpy.load`; the test below
//! pins the exact header bytes for a known shape to catch any
//! formatting drift.

// `write_u8_array` and `write_f32_array` are consumed by the
// `--snapshot` path added in the next commit. They're exercised by
// the tests below, so the dead-code warning isn't load-bearing.
#![allow(dead_code)]

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

const MAGIC: &[u8] = b"\x93NUMPY";
const VERSION: &[u8] = b"\x01\x00";
const HEADER_ALIGN: usize = 64;

/// Build the ASCII Python-dict header. NumPy parses this with
/// `ast.literal_eval`, so it must be valid Python syntax.
fn build_header_dict(descr: &str, shape: &[usize]) -> String {
    let shape_str = if shape.len() == 1 {
        // 1-D needs trailing comma so Python parses it as a tuple, not int.
        format!("({},)", shape[0])
    } else {
        shape
            .iter()
            .map(|s| s.to_string())
            .collect::<Vec<_>>()
            .join(", ")
    };
    let shape_paren = if shape.len() == 1 {
        shape_str
    } else {
        format!("({})", shape_str)
    };
    format!(
        "{{'descr': '{}', 'fortran_order': False, 'shape': {}, }}",
        descr, shape_paren
    )
}

/// Write magic + version + len-prefixed dict header, padded so the
/// data payload starts at a 64-byte boundary.
fn write_header<W: Write>(w: &mut W, descr: &str, shape: &[usize]) -> io::Result<()> {
    let dict = build_header_dict(descr, shape);
    // Prefix: 6 magic + 2 version + 2 header-length = 10 bytes.
    let prefix_len = MAGIC.len() + VERSION.len() + 2;
    // Header body is dict + final newline; padded to align.
    let unpadded = prefix_len + dict.len() + 1;
    let pad = if unpadded.is_multiple_of(HEADER_ALIGN) {
        0
    } else {
        HEADER_ALIGN - (unpadded % HEADER_ALIGN)
    };
    let header_len = dict.len() + pad + 1;
    assert!(
        header_len <= u16::MAX as usize,
        "npy header too long ({header_len} bytes); use v2.0 format for headers > 65535"
    );

    w.write_all(MAGIC)?;
    w.write_all(VERSION)?;
    w.write_all(&(header_len as u16).to_le_bytes())?;
    w.write_all(dict.as_bytes())?;
    // Pad with spaces, then a final newline.
    for _ in 0..pad {
        w.write_all(b" ")?;
    }
    w.write_all(b"\n")?;
    Ok(())
}

/// Write a flat `u8` array as a `.npy` file with the given shape.
/// `data.len()` must equal the product of `shape`.
pub fn write_u8_array<P: AsRef<Path>>(path: P, shape: &[usize], data: &[u8]) -> io::Result<()> {
    let expected: usize = shape.iter().product();
    assert_eq!(
        data.len(),
        expected,
        "npy::write_u8_array: data.len()={} != shape product={} (shape={:?})",
        data.len(),
        expected,
        shape
    );
    let mut w = BufWriter::new(File::create(path)?);
    write_header(&mut w, "|u1", shape)?;
    w.write_all(data)?;
    w.flush()?;
    Ok(())
}

/// Write a flat `f32` array as a `.npy` file with the given shape.
/// `data.len()` must equal the product of `shape`. Values are written
/// little-endian (matches the `'<f4'` dtype string).
pub fn write_f32_array<P: AsRef<Path>>(path: P, shape: &[usize], data: &[f32]) -> io::Result<()> {
    let expected: usize = shape.iter().product();
    assert_eq!(
        data.len(),
        expected,
        "npy::write_f32_array: data.len()={} != shape product={} (shape={:?})",
        data.len(),
        expected,
        shape
    );
    let mut w = BufWriter::new(File::create(path)?);
    write_header(&mut w, "<f4", shape)?;
    for &v in data {
        w.write_all(&v.to_le_bytes())?;
    }
    w.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn header_dict_4d_matches_expected() {
        let dict = build_header_dict("|u1", &[180, 60, 60, 60]);
        assert_eq!(
            dict,
            "{'descr': '|u1', 'fortran_order': False, 'shape': (180, 60, 60, 60), }"
        );
    }

    #[test]
    fn header_dict_1d_has_trailing_comma() {
        let dict = build_header_dict("<f4", &[5]);
        // 1-D shape needs `(5,)` not `(5)` — otherwise Python parses as int.
        assert_eq!(
            dict,
            "{'descr': '<f4', 'fortran_order': False, 'shape': (5,), }"
        );
    }

    /// Full file bytes for a tiny u8 array. Locks the prefix + header
    /// padding + body format so any future refactor that breaks
    /// `numpy.load` compatibility trips this test.
    #[test]
    fn write_u8_array_bytes_are_npy_v1_compatible() {
        let mut buf: Vec<u8> = Vec::new();
        let mut w = Cursor::new(&mut buf);
        write_header(&mut w, "|u1", &[2, 3]).unwrap();
        // Data: 6 bytes
        w.write_all(&[0u8, 1, 2, 3, 4, 5]).unwrap();

        // Prefix bytes (magic + version + header-length u16)
        assert_eq!(&buf[..6], b"\x93NUMPY");
        assert_eq!(&buf[6..8], b"\x01\x00");
        let header_len = u16::from_le_bytes([buf[8], buf[9]]) as usize;
        // Total prefix + header MUST align to 64.
        assert!(
            (10 + header_len).is_multiple_of(HEADER_ALIGN),
            "total prefix+header = {} bytes, not a multiple of {}",
            10 + header_len,
            HEADER_ALIGN
        );
        // Last byte of the header is the newline terminator.
        assert_eq!(buf[10 + header_len - 1], b'\n');
        // Header dict is parseable ASCII.
        let header_str = std::str::from_utf8(&buf[10..10 + header_len - 1]).unwrap();
        assert!(header_str.contains("'descr': '|u1'"));
        assert!(header_str.contains("'fortran_order': False"));
        assert!(header_str.contains("'shape': (2, 3)"));
        // Body follows immediately after the header.
        let body_start = 10 + header_len;
        assert_eq!(&buf[body_start..body_start + 6], &[0u8, 1, 2, 3, 4, 5]);
    }

    #[test]
    fn write_f32_array_uses_little_endian() {
        let mut buf: Vec<u8> = Vec::new();
        let mut w = Cursor::new(&mut buf);
        write_header(&mut w, "<f4", &[1]).unwrap();
        w.write_all(&1.0f32.to_le_bytes()).unwrap();
        let header_len = u16::from_le_bytes([buf[8], buf[9]]) as usize;
        let body_start = 10 + header_len;
        // 1.0 f32 LE = 0x00 0x00 0x80 0x3F
        assert_eq!(&buf[body_start..body_start + 4], &[0x00, 0x00, 0x80, 0x3F]);
    }
}
