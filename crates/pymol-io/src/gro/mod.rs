//! GROMACS GRO file format parser
//!
//! Provides reading of GRO coordinate format files.

mod parser;
mod writer;

pub use parser::GroReader;
pub use writer::GroWriter;

use std::io::Read;
use std::path::Path;

use pymol_mol::ObjectMolecule;

use crate::error::IoResult;
use crate::traits::MoleculeReader;

/// Read a GRO file from a path
pub fn read_gro(path: &Path) -> IoResult<ObjectMolecule> {
    let file = crate::compress::open_file(path)?;
    let mut reader = GroReader::new(file);
    reader.read()
}

/// Read a GRO file from a string
pub fn read_gro_str(content: &str) -> IoResult<ObjectMolecule> {
    let mut reader = GroReader::new(content.as_bytes());
    reader.read()
}

/// Read a GRO file from a reader
pub fn read_gro_from<R: Read>(reader: R) -> IoResult<ObjectMolecule> {
    let mut reader = GroReader::new(reader);
    reader.read()
}
