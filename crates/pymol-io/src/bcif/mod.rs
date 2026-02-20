//! BinaryCIF file format parser
//!
//! Provides reading of BinaryCIF (bCIF) format files.
//! bCIF is a binary encoding of the CIF/mmCIF format using MessagePack.

mod decode;
mod parser;
pub(crate) mod types;

pub use parser::BcifReader;

use std::io::Read;
use std::path::Path;

use pymol_mol::ObjectMolecule;

use crate::error::IoResult;
use crate::traits::MoleculeReader;

/// Read a bCIF file from a path
pub fn read_bcif(path: &Path) -> IoResult<ObjectMolecule> {
    let file = crate::compress::open_file(path)?;
    let mut reader = BcifReader::new(file);
    reader.read()
}

/// Read a bCIF file from raw bytes
pub fn read_bcif_bytes(bytes: &[u8]) -> IoResult<ObjectMolecule> {
    let mut reader = BcifReader::new(bytes);
    reader.read()
}

/// Read a bCIF file from a reader
pub fn read_bcif_from<R: Read>(reader: R) -> IoResult<ObjectMolecule> {
    let mut reader = BcifReader::new(reader);
    reader.read()
}
