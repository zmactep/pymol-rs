//! Molecular file format I/O
//!
//! This crate provides parsers and writers for various molecular file formats:
//!
//! - **PDB** - Protein Data Bank format
//! - **SDF/MOL** - MDL Structure Data File / MOL V2000 format
//! - **MOL2** - TRIPOS MOL2 format
//! - **XYZ** - Simple XYZ coordinate format
//! - **mmCIF** - Macromolecular Crystallographic Information File
//!
//! # Quick Start
//!
//! ```no_run
//! use pymol_io::{read_file, write_file, FileFormat};
//! use std::path::Path;
//!
//! // Read a PDB file
//! let mol = pymol_io::pdb::read_pdb(Path::new("protein.pdb")).unwrap();
//!
//! // Write to a different format
//! // pymol_io::sdf::write_sdf(Path::new("protein.sdf"), &mol).unwrap();
//! ```
//!
//! # Fetching from RCSB PDB
//!
//! With the `fetch` or `fetch-async` feature enabled, you can download structures
//! directly from the RCSB Protein Data Bank:
//!
//! ```no_run
//! # #[cfg(feature = "fetch")]
//! # fn main() -> pymol_io::IoResult<()> {
//! use pymol_io::fetch::{fetch, FetchFormat};
//!
//! // Fetch ubiquitin in mmCIF format
//! let mol = fetch("1ubq", FetchFormat::default())?;
//! # Ok(())
//! # }
//! # #[cfg(not(feature = "fetch"))]
//! # fn main() {}
//! ```
//!
//! # Format-specific modules
//!
//! Each format has its own module with specialized readers and writers:
//!
//! - [`pdb`] - PDB format with ATOM/HETATM/CONECT records
//! - [`sdf`] - SDF/MOL format with connection tables
//! - [`mol2`] - MOL2 format with SYBYL atom types and charges
//! - [`xyz`] - Simple XYZ format (coordinates only)
//! - [`cif`] - mmCIF format with STAR syntax
//!
//! # Features
//!
//! - `fetch` - Enable synchronous fetching from RCSB PDB (uses `ureq`)
//! - `fetch-async` - Enable asynchronous fetching from RCSB PDB (uses `reqwest`)

pub mod cif;
pub mod compress;
pub mod detect;
pub mod error;
pub mod mol2;
pub mod pdb;
pub mod sdf;
pub mod traits;
pub mod xyz;

// Fetch module (requires feature)
#[cfg(any(feature = "fetch", feature = "fetch-async"))]
pub mod fetch;

// Re-exports
pub use error::{IoError, IoResult};
pub use traits::{
    create_reader, create_writer, FileFormat, MoleculeReader, MoleculeWriter, ReadOptions,
    WriteOptions,
};

// Fetch re-exports
#[cfg(any(feature = "fetch", feature = "fetch-async"))]
pub use fetch::FetchFormat;
#[cfg(feature = "fetch")]
pub use fetch::fetch;
#[cfg(feature = "fetch-async")]
pub use fetch::fetch_async;

use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use pymol_mol::ObjectMolecule;

/// Read a molecule from a file, auto-detecting the format
pub fn read_file(path: &Path) -> IoResult<ObjectMolecule> {
    let format = detect::detect_from_path(path);
    read_file_format(path, format)
}

/// Read a molecule from a file with a specific format
pub fn read_file_format(path: &Path, format: FileFormat) -> IoResult<ObjectMolecule> {
    match format {
        FileFormat::Pdb => pdb::read_pdb(path),
        FileFormat::Sdf => sdf::read_sdf(path).map(|v| {
            v.into_iter()
                .next()
                .ok_or(IoError::EmptyFile)
                .unwrap_or_default()
        }),
        FileFormat::Mol2 => mol2::read_mol2(path).map(|v| {
            v.into_iter()
                .next()
                .ok_or(IoError::EmptyFile)
                .unwrap_or_default()
        }),
        FileFormat::Xyz => xyz::read_xyz(path),
        FileFormat::Cif => cif::read_cif(path),
        FileFormat::Unknown => Err(IoError::UnknownFormat(
            path.to_string_lossy().into_owned(),
        )),
    }
}

/// Read all molecules from a file (for multi-molecule formats like SDF)
pub fn read_all(path: &Path) -> IoResult<Vec<ObjectMolecule>> {
    let format = detect::detect_from_path(path);
    read_all_format(path, format)
}

/// Read all molecules from a file with a specific format
pub fn read_all_format(path: &Path, format: FileFormat) -> IoResult<Vec<ObjectMolecule>> {
    match format {
        FileFormat::Pdb => pdb::read_pdb(path).map(|m| vec![m]),
        FileFormat::Sdf => sdf::read_sdf(path),
        FileFormat::Mol2 => mol2::read_mol2(path),
        FileFormat::Xyz => xyz::read_xyz(path).map(|m| vec![m]),
        FileFormat::Cif => cif::read_cif(path).map(|m| vec![m]),
        FileFormat::Unknown => Err(IoError::UnknownFormat(
            path.to_string_lossy().into_owned(),
        )),
    }
}

/// Write a molecule to a file, auto-detecting the format from extension
pub fn write_file(path: &Path, mol: &ObjectMolecule) -> IoResult<()> {
    let format = detect::detect_from_path(path);
    write_file_format(path, mol, format)
}

/// Write a molecule to a file with a specific format
pub fn write_file_format(path: &Path, mol: &ObjectMolecule, format: FileFormat) -> IoResult<()> {
    let file = File::create(path)?;
    let writer = BufWriter::new(file);
    let mut mol_writer = create_writer(writer, format)?;
    mol_writer.write(mol)?;
    mol_writer.flush()
}

/// Write multiple molecules to a file (for multi-molecule formats like SDF)
pub fn write_all(path: &Path, molecules: &[ObjectMolecule]) -> IoResult<()> {
    let format = detect::detect_from_path(path);
    write_all_format(path, molecules, format)
}

/// Write multiple molecules to a file with a specific format
pub fn write_all_format(
    path: &Path,
    molecules: &[ObjectMolecule],
    format: FileFormat,
) -> IoResult<()> {
    let file = File::create(path)?;
    let writer = BufWriter::new(file);
    let mut mol_writer = create_writer(writer, format)?;
    mol_writer.write_all(molecules)?;
    mol_writer.flush()
}

/// Parse a molecule from a string with the given format
pub fn parse_str(content: &str, format: FileFormat) -> IoResult<ObjectMolecule> {
    match format {
        FileFormat::Pdb => pdb::read_pdb_str(content),
        FileFormat::Sdf => sdf::read_sdf_str(content).map(|v| {
            v.into_iter()
                .next()
                .ok_or(IoError::EmptyFile)
                .unwrap_or_default()
        }),
        FileFormat::Mol2 => mol2::read_mol2_str(content).map(|v| {
            v.into_iter()
                .next()
                .ok_or(IoError::EmptyFile)
                .unwrap_or_default()
        }),
        FileFormat::Xyz => xyz::read_xyz_str(content),
        FileFormat::Cif => cif::read_cif_str(content),
        FileFormat::Unknown => Err(IoError::UnknownFormat("string input".to_string())),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_detection() {
        assert_eq!(
            FileFormat::from_path(Path::new("test.pdb")),
            FileFormat::Pdb
        );
        assert_eq!(
            FileFormat::from_path(Path::new("test.sdf")),
            FileFormat::Sdf
        );
        assert_eq!(
            FileFormat::from_path(Path::new("test.mol2")),
            FileFormat::Mol2
        );
        assert_eq!(
            FileFormat::from_path(Path::new("test.xyz")),
            FileFormat::Xyz
        );
        assert_eq!(
            FileFormat::from_path(Path::new("test.cif")),
            FileFormat::Cif
        );
    }

    #[test]
    fn test_gzip_format_detection() {
        assert_eq!(
            FileFormat::from_path(Path::new("test.pdb.gz")),
            FileFormat::Pdb
        );
    }
}
