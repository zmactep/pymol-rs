//! XYZ file format parser and writer
//!
//! Provides reading and writing of XYZ coordinate format files.

mod parser;
mod writer;

pub use parser::XyzReader;
pub use writer::XyzWriter;

use std::io::Read;
use std::path::Path;

use pymol_mol::ObjectMolecule;

use crate::error::IoResult;
use crate::traits::MoleculeReader;

/// Read an XYZ file from a path
pub fn read_xyz(path: &Path) -> IoResult<ObjectMolecule> {
    let file = crate::compress::open_file(path)?;
    let mut reader = XyzReader::new(file);
    reader.read()
}

/// Read an XYZ file from a string
pub fn read_xyz_str(content: &str) -> IoResult<ObjectMolecule> {
    let mut reader = XyzReader::new(content.as_bytes());
    reader.read()
}

/// Read an XYZ file from a reader
pub fn read_xyz_from<R: Read>(reader: R) -> IoResult<ObjectMolecule> {
    let mut reader = XyzReader::new(reader);
    reader.read()
}
