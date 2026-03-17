//! GROMACS XTC trajectory reader
//!
//! Uses the `molly` crate for XTC decompression.
//! Coordinates are converted from nanometers to Angstroms (×10).

use std::io::Read;

use pymol_mol::CoordSet;

use crate::error::{IoError, IoResult};
use crate::traits::{TrajectoryReadOptions, TrajectoryReader};
use crate::units::NM_TO_ANGSTROM;

/// XTC trajectory reader
pub struct XtcReader<R> {
    reader: molly::XTCReader<R>,
}

impl<R: Read> XtcReader<R> {
    /// Create a new XTC reader
    pub fn new(reader: R) -> Self {
        XtcReader {
            reader: molly::XTCReader::new(reader),
        }
    }
}

impl<R: Read> TrajectoryReader for XtcReader<R> {
    fn read_frames(&mut self, opts: &TrajectoryReadOptions) -> IoResult<Vec<CoordSet>> {
        let interval = opts.interval.max(1);
        let mut frames = Vec::new();
        let mut frame_idx: usize = 0;
        let mut frame = molly::Frame::default();

        loop {
            match self.reader.read_frame(&mut frame) {
                Ok(()) => {}
                Err(_) => break,
            }

            if frame_idx >= opts.start
                && opts.stop.map_or(true, |s| frame_idx < s)
                && (frame_idx - opts.start) % interval == 0
            {
                // Convert nm to Angstroms
                let mut coords = frame.positions.clone();
                for c in &mut coords {
                    *c *= NM_TO_ANGSTROM;
                }
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
