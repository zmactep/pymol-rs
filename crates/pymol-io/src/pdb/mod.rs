//! PDB file format parser and writer
//!
//! Provides reading and writing of Protein Data Bank (PDB) format files.

pub mod hybrid36;
mod parser;
mod records;
mod writer;

pub use parser::PdbReader;
pub use records::*;
pub use writer::PdbWriter;

use std::io::Read;
use std::path::Path;

use pymol_mol::ObjectMolecule;

use crate::error::IoResult;
use crate::traits::MoleculeReader;

/// Read a PDB file from a path
pub fn read_pdb(path: &Path) -> IoResult<ObjectMolecule> {
    let file = crate::compress::open_file(path)?;
    let mut reader = PdbReader::new(file);
    reader.read()
}

/// Read a PDB file from a string
pub fn read_pdb_str(content: &str) -> IoResult<ObjectMolecule> {
    let mut reader = PdbReader::new(content.as_bytes());
    reader.read()
}

/// Read a PDB file from a reader
pub fn read_pdb_from<R: Read>(reader: R) -> IoResult<ObjectMolecule> {
    let mut reader = PdbReader::new(reader);
    reader.read()
}
