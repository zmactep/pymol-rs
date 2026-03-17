//! GROMACS TRR trajectory reader
//!
//! Native Rust implementation using XDR (big-endian) decoding.
//! Coordinates are converted from nanometers to Angstroms (×10).
//!
//! TRR frame layout:
//! - Magic number (i32 = 1993)
//! - XDR string: version string (e.g. "GMX_trn_file")
//! - Header sizes: ir_size, e_size, box_size, vir_size, pres_size,
//!   top_size, sym_size, x_size, v_size, f_size (all i32)
//! - natoms (i32)
//! - step (i32)
//! - nre (i32)
//! - time (float32 or float64, determined by a size marker)
//! - lambda (float32 or float64)
//! - Box matrix (3×3 floats/doubles, box_size bytes)
//! - Positions (natoms × 3 floats/doubles, x_size bytes)
//! - Velocities (natoms × 3 floats/doubles, v_size bytes, optional)
//! - Forces (natoms × 3 floats/doubles, f_size bytes, optional)

use std::io::{BufReader, Read};

use byteorder::{BigEndian, ReadBytesExt};
use pymol_mol::CoordSet;

use crate::error::{IoError, IoResult};
use crate::traits::{TrajectoryReadOptions, TrajectoryReader};
use crate::units::NM_TO_ANGSTROM;

/// TRR magic number
const TRR_MAGIC: i32 = 1993;

/// Size of a double-precision float
const SIZE_DOUBLE: i32 = 8;

/// Number of spatial dimensions
const DIM: i32 = 3;

/// TRR trajectory reader
pub struct TrrReader<R> {
    reader: BufReader<R>,
}

impl<R: Read> TrrReader<R> {
    /// Create a new TRR reader
    pub fn new(reader: R) -> Self {
        TrrReader {
            reader: BufReader::new(reader),
        }
    }

    /// Read an XDR string: i32 length, then padded to 4-byte boundary
    fn read_xdr_string(&mut self) -> IoResult<Option<String>> {
        let len = match self.reader.read_i32::<BigEndian>() {
            Ok(l) => l as usize,
            Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => return Ok(None),
            Err(e) => return Err(IoError::Io(e)),
        };

        let mut buf = vec![0u8; len];
        self.reader.read_exact(&mut buf)?;

        // XDR pads strings to 4-byte boundary
        let padding = (4 - (len % 4)) % 4;
        if padding > 0 {
            let mut pad = vec![0u8; padding];
            self.reader.read_exact(&mut pad)?;
        }

        Ok(Some(
            String::from_utf8_lossy(&buf).trim_end_matches('\0').to_string(),
        ))
    }

    /// Skip a number of bytes
    fn skip(&mut self, n: usize) -> IoResult<()> {
        let mut remaining = n;
        let mut buf = [0u8; 4096];
        while remaining > 0 {
            let to_read = remaining.min(buf.len());
            self.reader.read_exact(&mut buf[..to_read])?;
            remaining -= to_read;
        }
        Ok(())
    }

    /// Read coordinates as f32 (single precision), converting nm → Å
    fn read_coords_f32(&mut self, natoms: usize) -> IoResult<Vec<f32>> {
        let n = natoms * DIM as usize;
        let mut coords = Vec::with_capacity(n);
        for _ in 0..n {
            let val = self.reader.read_f32::<BigEndian>()?;
            coords.push(val * NM_TO_ANGSTROM);
        }
        Ok(coords)
    }

    /// Read coordinates as f64 (double precision), converting nm → Å and downcasting to f32
    fn read_coords_f64(&mut self, natoms: usize) -> IoResult<Vec<f32>> {
        let n = natoms * DIM as usize;
        let mut coords = Vec::with_capacity(n);
        for _ in 0..n {
            let val = self.reader.read_f64::<BigEndian>()?;
            coords.push((val * NM_TO_ANGSTROM as f64) as f32);
        }
        Ok(coords)
    }

    /// Read a single TRR frame, returning coords or None at EOF
    fn read_frame(&mut self) -> IoResult<Option<Vec<f32>>> {
        // Magic number
        let magic = match self.reader.read_i32::<BigEndian>() {
            Ok(m) => m,
            Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => return Ok(None),
            Err(e) => return Err(IoError::Io(e)),
        };

        if magic != TRR_MAGIC {
            return Err(IoError::parse_msg(format!(
                "Invalid TRR magic number: {} (expected {})",
                magic, TRR_MAGIC
            )));
        }

        // Version string
        let _version = self.read_xdr_string()?;

        // Header size fields
        let _ir_size = self.reader.read_i32::<BigEndian>()?;
        let _e_size = self.reader.read_i32::<BigEndian>()?;
        let box_size = self.reader.read_i32::<BigEndian>()?;
        let _vir_size = self.reader.read_i32::<BigEndian>()?;
        let _pres_size = self.reader.read_i32::<BigEndian>()?;
        let _top_size = self.reader.read_i32::<BigEndian>()?;
        let _sym_size = self.reader.read_i32::<BigEndian>()?;
        let x_size = self.reader.read_i32::<BigEndian>()?;
        let v_size = self.reader.read_i32::<BigEndian>()?;
        let f_size = self.reader.read_i32::<BigEndian>()?;

        // natoms, step, nre
        let natoms = self.reader.read_i32::<BigEndian>()? as usize;
        let _step = self.reader.read_i32::<BigEndian>()?;
        let _nre = self.reader.read_i32::<BigEndian>()?;

        // Determine float precision from box_size or x_size
        let is_double = if box_size > 0 {
            box_size / (DIM * DIM) == SIZE_DOUBLE
        } else if x_size > 0 {
            x_size / (natoms as i32 * DIM) == SIZE_DOUBLE
        } else {
            false
        };

        // Time and lambda
        if is_double {
            let _time = self.reader.read_f64::<BigEndian>()?;
            let _lambda = self.reader.read_f64::<BigEndian>()?;
        } else {
            let _time = self.reader.read_f32::<BigEndian>()?;
            let _lambda = self.reader.read_f32::<BigEndian>()?;
        }

        // Skip box matrix
        if box_size > 0 {
            self.skip(box_size as usize)?;
        }

        // Read coordinates
        let coords = if x_size > 0 {
            if is_double {
                self.read_coords_f64(natoms)?
            } else {
                self.read_coords_f32(natoms)?
            }
        } else {
            return Ok(Some(Vec::new()));
        };

        // Skip velocities and forces
        if v_size > 0 {
            self.skip(v_size as usize)?;
        }
        if f_size > 0 {
            self.skip(f_size as usize)?;
        }

        Ok(Some(coords))
    }
}

impl<R: Read> TrajectoryReader for TrrReader<R> {
    fn read_frames(&mut self, opts: &TrajectoryReadOptions) -> IoResult<Vec<CoordSet>> {
        let interval = opts.interval.max(1);
        let mut frames = Vec::new();
        let mut frame_idx: usize = 0;

        while let Some(coords) = self.read_frame()? {
            if coords.is_empty() {
                frame_idx += 1;
                continue;
            }

            if frame_idx >= opts.start
                && opts.stop.map_or(true, |s| frame_idx < s)
                && (frame_idx - opts.start) % interval == 0
            {
                frames.push(CoordSet::from_coords(coords));
            }

            frame_idx += 1;

            if opts.stop.map_or(false, |s| frame_idx >= s) {
                break;
            }
        }

        if frames.is_empty() && frame_idx == 0 {
            return Err(IoError::EmptyFile);
        }

        Ok(frames)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a minimal valid TRR frame in XDR (big-endian) format
    fn build_trr_frame(natoms: usize, coords_nm: &[f32]) -> Vec<u8> {
        use byteorder::WriteBytesExt;

        let mut buf = Vec::new();

        // Magic
        buf.write_i32::<BigEndian>(TRR_MAGIC).unwrap();

        // Version string: "GMX_trn_file" (12 bytes, padded to 16)
        let version = b"GMX_trn_file";
        buf.write_i32::<BigEndian>(version.len() as i32).unwrap();
        buf.extend_from_slice(version);
        // Pad to 4-byte boundary (12 bytes -> no padding needed since 12 % 4 == 0)

        let x_size = (natoms as i32) * DIM * 4;
        let box_size = DIM * DIM * 4; // 36 bytes

        // Header sizes: ir, e, box, vir, pres, top, sym, x, v, f
        buf.write_i32::<BigEndian>(0).unwrap(); // ir_size
        buf.write_i32::<BigEndian>(0).unwrap(); // e_size
        buf.write_i32::<BigEndian>(box_size).unwrap(); // box_size
        buf.write_i32::<BigEndian>(0).unwrap(); // vir_size
        buf.write_i32::<BigEndian>(0).unwrap(); // pres_size
        buf.write_i32::<BigEndian>(0).unwrap(); // top_size
        buf.write_i32::<BigEndian>(0).unwrap(); // sym_size
        buf.write_i32::<BigEndian>(x_size).unwrap(); // x_size
        buf.write_i32::<BigEndian>(0).unwrap(); // v_size
        buf.write_i32::<BigEndian>(0).unwrap(); // f_size

        // natoms, step, nre
        buf.write_i32::<BigEndian>(natoms as i32).unwrap();
        buf.write_i32::<BigEndian>(0).unwrap(); // step
        buf.write_i32::<BigEndian>(0).unwrap(); // nre

        // Time and lambda (single precision since box uses single)
        buf.write_f32::<BigEndian>(0.0).unwrap(); // time
        buf.write_f32::<BigEndian>(0.0).unwrap(); // lambda

        // Box matrix (3x3 identity in nm)
        for i in 0..3 {
            for j in 0..3 {
                buf.write_f32::<BigEndian>(if i == j { 1.0 } else { 0.0 })
                    .unwrap();
            }
        }

        // Coordinates (in nm)
        for &c in coords_nm {
            buf.write_f32::<BigEndian>(c).unwrap();
        }

        buf
    }

    #[test]
    fn test_trr_single_frame() {
        // 2 atoms at (0.1, 0.2, 0.3) and (0.4, 0.5, 0.6) in nm
        let coords_nm = [0.1f32, 0.2, 0.3, 0.4, 0.5, 0.6];
        let data = build_trr_frame(2, &coords_nm);

        let mut reader = TrrReader::new(data.as_slice());
        let opts = TrajectoryReadOptions::default();
        let frames = reader.read_frames(&opts).unwrap();

        assert_eq!(frames.len(), 1);
        assert_eq!(frames[0].len(), 2);

        // Check nm → Å conversion
        let coord0 = frames[0].get_coord(pymol_mol::CoordIndex(0)).unwrap();
        assert!((coord0.x - 1.0).abs() < 0.01);
        assert!((coord0.y - 2.0).abs() < 0.01);
        assert!((coord0.z - 3.0).abs() < 0.01);

        let coord1 = frames[0].get_coord(pymol_mol::CoordIndex(1)).unwrap();
        assert!((coord1.x - 4.0).abs() < 0.01);
        assert!((coord1.y - 5.0).abs() < 0.01);
        assert!((coord1.z - 6.0).abs() < 0.01);
    }

    #[test]
    fn test_trr_multi_frame() {
        let frame1 = build_trr_frame(2, &[0.1, 0.2, 0.3, 0.4, 0.5, 0.6]);
        let frame2 = build_trr_frame(2, &[0.7, 0.8, 0.9, 1.0, 1.1, 1.2]);

        let mut data = Vec::new();
        data.extend_from_slice(&frame1);
        data.extend_from_slice(&frame2);

        let mut reader = TrrReader::new(data.as_slice());
        let opts = TrajectoryReadOptions::default();
        let frames = reader.read_frames(&opts).unwrap();

        assert_eq!(frames.len(), 2);

        // Frame 2: first atom at (7.0, 8.0, 9.0) Å
        let coord = frames[1].get_coord(pymol_mol::CoordIndex(0)).unwrap();
        assert!((coord.x - 7.0).abs() < 0.01);
    }

    #[test]
    fn test_trr_frame_filtering() {
        // Build 5 frames
        let mut data = Vec::new();
        for i in 0..5 {
            let x = i as f32 * 0.1;
            data.extend_from_slice(&build_trr_frame(1, &[x, 0.0, 0.0]));
        }

        // Read with start=1, stop=4, interval=2 → frames 1 and 3
        let mut reader = TrrReader::new(data.as_slice());
        let opts = TrajectoryReadOptions {
            start: 1,
            stop: Some(4),
            interval: 2,
        };
        let frames = reader.read_frames(&opts).unwrap();

        assert_eq!(frames.len(), 2);

        // Frame at index 1: x = 0.1 nm = 1.0 Å
        let c0 = frames[0].get_coord(pymol_mol::CoordIndex(0)).unwrap();
        assert!((c0.x - 1.0).abs() < 0.01);

        // Frame at index 3: x = 0.3 nm = 3.0 Å
        let c1 = frames[1].get_coord(pymol_mol::CoordIndex(0)).unwrap();
        assert!((c1.x - 3.0).abs() < 0.01);
    }

    #[test]
    fn test_trr_empty_file() {
        let data: Vec<u8> = Vec::new();
        let mut reader = TrrReader::new(data.as_slice());
        let opts = TrajectoryReadOptions::default();
        let result = reader.read_frames(&opts);
        assert!(result.is_err());
    }
}
