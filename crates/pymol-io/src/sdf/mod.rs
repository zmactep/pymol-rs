//! SDF/MOL file format parser and writer
//!
//! Provides reading and writing of MDL SDF/MOL V2000 format files.

mod parser;
mod writer;

pub use parser::SdfReader;
pub use writer::SdfWriter;

use std::io::Read;
use std::path::Path;

use pymol_mol::ObjectMolecule;

use crate::error::IoResult;
use crate::traits::MoleculeReader;

/// Read an SDF/MOL file from a path
pub fn read_sdf(path: &Path) -> IoResult<Vec<ObjectMolecule>> {
    let file = crate::compress::open_file(path)?;
    let mut reader = SdfReader::new(file);
    reader.read_all()
}

/// Read an SDF/MOL file from a string
pub fn read_sdf_str(content: &str) -> IoResult<Vec<ObjectMolecule>> {
    let mut reader = SdfReader::new(content.as_bytes());
    reader.read_all()
}

/// Read an SDF/MOL file from a reader
pub fn read_sdf_from<R: Read>(reader: R) -> IoResult<Vec<ObjectMolecule>> {
    let mut reader = SdfReader::new(reader);
    reader.read_all()
}
