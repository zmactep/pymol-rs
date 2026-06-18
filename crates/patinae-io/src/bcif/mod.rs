//! BinaryCIF file format parser
//!
//! Provides reading of BinaryCIF (bCIF) format files.
//! bCIF is a binary encoding of the CIF/mmCIF format using MessagePack.

mod decode;
mod parser;
#[cfg(test)]
pub(crate) mod test_support;
pub(crate) mod types;

pub use parser::BcifReader;

use std::io::Read;
use std::path::Path;

use patinae_mol::ObjectMolecule;

use crate::error::IoResult;
use crate::traits::MoleculeReader;

/// Read a bCIF file from a path
pub fn read_bcif(path: &Path) -> IoResult<ObjectMolecule> {
    read_bcif_with_bond_tolerance(path, patinae_mol::DEFAULT_BOND_TOLERANCE)
}

pub fn read_bcif_with_bond_tolerance(path: &Path, bond_tolerance: f32) -> IoResult<ObjectMolecule> {
    let file = crate::compress::open_file(path)?;
    let mut reader = BcifReader::new(file).with_bond_tolerance(bond_tolerance);
    reader.read()
}

pub(crate) fn read_all_bcif_with_bond_tolerance(
    path: &Path,
    bond_tolerance: f32,
) -> IoResult<Vec<ObjectMolecule>> {
    let file = crate::compress::open_file(path)?;
    let mut reader = BcifReader::new(file).with_bond_tolerance(bond_tolerance);
    reader.read_all()
}

/// Read a bCIF file from raw bytes
pub fn read_bcif_bytes(bytes: &[u8]) -> IoResult<ObjectMolecule> {
    let mut reader =
        BcifReader::new(bytes).with_bond_tolerance(patinae_mol::DEFAULT_BOND_TOLERANCE);
    reader.read()
}

pub fn read_bcif_bytes_with_bond_tolerance(
    bytes: &[u8],
    bond_tolerance: f32,
) -> IoResult<ObjectMolecule> {
    let mut reader = BcifReader::new(bytes).with_bond_tolerance(bond_tolerance);
    reader.read()
}

/// Read a bCIF file from a reader
pub fn read_bcif_from<R: Read>(reader: R) -> IoResult<ObjectMolecule> {
    let mut reader =
        BcifReader::new(reader).with_bond_tolerance(patinae_mol::DEFAULT_BOND_TOLERANCE);
    reader.read()
}

pub fn read_bcif_from_with_bond_tolerance<R: Read>(
    reader: R,
    bond_tolerance: f32,
) -> IoResult<ObjectMolecule> {
    let mut reader = BcifReader::new(reader).with_bond_tolerance(bond_tolerance);
    reader.read()
}
