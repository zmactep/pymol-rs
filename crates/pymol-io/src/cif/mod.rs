//! mmCIF file format parser and writer
//!
//! Provides reading and writing of Macromolecular Crystallographic Information File format.

pub(crate) mod common;
mod lexer;
mod parser;
mod writer;

pub use parser::CifReader;
pub use writer::CifWriter;

use std::io::Read;
use std::path::Path;

use pymol_mol::ObjectMolecule;

use crate::error::IoResult;
use crate::traits::MoleculeReader;

/// Read an mmCIF file from a path
pub fn read_cif(path: &Path) -> IoResult<ObjectMolecule> {
    let file = crate::compress::open_file(path)?;
    let mut reader = CifReader::new(file);
    reader.read()
}

/// Read an mmCIF file from a string
pub fn read_cif_str(content: &str) -> IoResult<ObjectMolecule> {
    let mut reader = CifReader::new(content.as_bytes());
    reader.read()
}

/// Read an mmCIF file from a reader
pub fn read_cif_from<R: Read>(reader: R) -> IoResult<ObjectMolecule> {
    let mut reader = CifReader::new(reader);
    reader.read()
}
