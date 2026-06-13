//! GROMACS GRO file format parser
//!
//! Provides reading of GRO coordinate format files.

mod parser;
mod writer;

pub use parser::GroReader;
pub use writer::GroWriter;

use std::io::Read;
use std::path::Path;

use patinae_mol::ObjectMolecule;

use crate::error::IoResult;
use crate::traits::MoleculeReader;

/// Read a GRO file from a path
pub fn read_gro(path: &Path) -> IoResult<ObjectMolecule> {
    read_gro_with_bond_tolerance(path, patinae_mol::DEFAULT_BOND_TOLERANCE)
}

pub fn read_gro_with_bond_tolerance(path: &Path, bond_tolerance: f32) -> IoResult<ObjectMolecule> {
    let file = crate::compress::open_file(path)?;
    let mut reader = GroReader::new(file).with_bond_tolerance(bond_tolerance);
    reader.read()
}

/// Read a GRO file from a string
pub fn read_gro_str(content: &str) -> IoResult<ObjectMolecule> {
    let mut reader =
        GroReader::new(content.as_bytes()).with_bond_tolerance(patinae_mol::DEFAULT_BOND_TOLERANCE);
    reader.read()
}

pub fn read_gro_str_with_bond_tolerance(
    content: &str,
    bond_tolerance: f32,
) -> IoResult<ObjectMolecule> {
    let mut reader = GroReader::new(content.as_bytes()).with_bond_tolerance(bond_tolerance);
    reader.read()
}

/// Read a GRO file from a reader
pub fn read_gro_from<R: Read>(reader: R) -> IoResult<ObjectMolecule> {
    let mut reader =
        GroReader::new(reader).with_bond_tolerance(patinae_mol::DEFAULT_BOND_TOLERANCE);
    reader.read()
}

pub fn read_gro_from_with_bond_tolerance<R: Read>(
    reader: R,
    bond_tolerance: f32,
) -> IoResult<ObjectMolecule> {
    let mut reader = GroReader::new(reader).with_bond_tolerance(bond_tolerance);
    reader.read()
}
