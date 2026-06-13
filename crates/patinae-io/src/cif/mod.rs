//! mmCIF file format parser and writer
//!
//! Provides reading and writing of Macromolecular Crystallographic Information File format.

pub(crate) mod common;
pub mod lexer;
mod parser;
mod writer;

pub use parser::CifReader;
pub use writer::CifWriter;

use std::io::Read;
use std::path::Path;

use patinae_mol::ObjectMolecule;

use crate::error::IoResult;
use crate::traits::MoleculeReader;

/// Read an mmCIF file from a path
pub fn read_cif(path: &Path) -> IoResult<ObjectMolecule> {
    read_cif_with_bond_tolerance(path, patinae_mol::DEFAULT_BOND_TOLERANCE)
}

pub fn read_cif_with_bond_tolerance(path: &Path, bond_tolerance: f32) -> IoResult<ObjectMolecule> {
    let file = crate::compress::open_file(path)?;
    let mut reader = CifReader::new(file).with_bond_tolerance(bond_tolerance);
    reader.read()
}

/// Read an mmCIF file from a string
pub fn read_cif_str(content: &str) -> IoResult<ObjectMolecule> {
    let mut reader =
        CifReader::new(content.as_bytes()).with_bond_tolerance(patinae_mol::DEFAULT_BOND_TOLERANCE);
    reader.read()
}

pub fn read_cif_str_with_bond_tolerance(
    content: &str,
    bond_tolerance: f32,
) -> IoResult<ObjectMolecule> {
    let mut reader = CifReader::new(content.as_bytes()).with_bond_tolerance(bond_tolerance);
    reader.read()
}

/// Read an mmCIF file from a reader
pub fn read_cif_from<R: Read>(reader: R) -> IoResult<ObjectMolecule> {
    let mut reader =
        CifReader::new(reader).with_bond_tolerance(patinae_mol::DEFAULT_BOND_TOLERANCE);
    reader.read()
}

pub fn read_cif_from_with_bond_tolerance<R: Read>(
    reader: R,
    bond_tolerance: f32,
) -> IoResult<ObjectMolecule> {
    let mut reader = CifReader::new(reader).with_bond_tolerance(bond_tolerance);
    reader.read()
}
