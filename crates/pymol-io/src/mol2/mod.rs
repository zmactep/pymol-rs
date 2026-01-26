//! MOL2 file format parser and writer
//!
//! Provides reading and writing of TRIPOS MOL2 format files.

mod parser;
mod writer;

pub use parser::Mol2Reader;
pub use writer::Mol2Writer;

use std::io::Read;
use std::path::Path;

use pymol_mol::ObjectMolecule;

use crate::error::IoResult;
use crate::traits::MoleculeReader;

/// Read a MOL2 file from a path
pub fn read_mol2(path: &Path) -> IoResult<Vec<ObjectMolecule>> {
    let file = crate::compress::open_file(path)?;
    let mut reader = Mol2Reader::new(file);
    reader.read_all()
}

/// Read a MOL2 file from a string
pub fn read_mol2_str(content: &str) -> IoResult<Vec<ObjectMolecule>> {
    let mut reader = Mol2Reader::new(content.as_bytes());
    reader.read_all()
}

/// Read a MOL2 file from a reader
pub fn read_mol2_from<R: Read>(reader: R) -> IoResult<Vec<ObjectMolecule>> {
    let mut reader = Mol2Reader::new(reader);
    reader.read_all()
}
