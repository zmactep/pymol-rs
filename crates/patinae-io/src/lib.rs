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
//! use patinae_io::{read_file, write_file, FileFormat};
//! use std::path::Path;
//!
//! // Read a PDB file
//! let mol = patinae_io::pdb::read_pdb(Path::new("protein.pdb")).unwrap();
//!
//! // Write to a different format
//! // patinae_io::sdf::write_sdf(Path::new("protein.sdf"), &mol).unwrap();
//! ```
//!
//! # Fetching from RCSB PDB
//!
//! With the `fetch` or `fetch-async` feature enabled, you can download structures
//! directly from the RCSB Protein Data Bank:
//!
//! ```no_run
//! # #[cfg(feature = "fetch")]
//! # fn main() -> patinae_io::IoResult<()> {
//! use patinae_io::fetch::{fetch, FetchFormat};
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

pub mod bcif;
pub mod ccp4;
pub mod cif;
pub mod compress;
pub mod detect;
pub mod error;
pub mod gro;
pub mod mol2;
pub mod pdb;
pub mod sdf;
pub mod traits;
#[cfg(feature = "traj")]
pub mod traj;
pub(crate) mod units;
pub mod xyz;

// Fetch module (requires feature)
#[cfg(any(feature = "fetch", feature = "fetch-async"))]
pub mod fetch;

// Re-exports
pub use error::{IoError, IoResult};
#[cfg(feature = "traj")]
pub use traits::create_trajectory_reader;
pub use traits::{
    create_reader, create_writer, FileFormat, MoleculeReader, MoleculeWriter, ReadOptions,
    TrajectoryReadOptions, TrajectoryReader, WriteOptions,
};

// Fetch re-exports
#[cfg(any(feature = "fetch", feature = "fetch-async"))]
pub use fetch::FetchFormat;
#[cfg(feature = "fetch")]
pub use fetch::{fetch, fetch_with_bond_tolerance};
#[cfg(feature = "fetch-async")]
pub use fetch::{
    fetch_async, fetch_async_with_bond_tolerance, fetch_pdb_metadata, PdbMetadataPreview,
};

use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use patinae_mol::ObjectMolecule;

fn first_molecule(molecules: Vec<ObjectMolecule>) -> IoResult<ObjectMolecule> {
    molecules.into_iter().next().ok_or(IoError::empty_file())
}

/// Read a molecule from a file, auto-detecting the format
pub fn read_file(path: &Path) -> IoResult<ObjectMolecule> {
    let format = detect::detect_from_path(path);
    read_file_format(path, format)
}

pub fn read_file_with_bond_tolerance(path: &Path, bond_tolerance: f32) -> IoResult<ObjectMolecule> {
    let format = detect::detect_from_path(path);
    read_file_format_with_bond_tolerance(path, format, bond_tolerance)
}

/// Read a molecule from a file with a specific format
pub fn read_file_format(path: &Path, format: FileFormat) -> IoResult<ObjectMolecule> {
    read_file_format_with_bond_tolerance(path, format, patinae_mol::DEFAULT_BOND_TOLERANCE)
}

pub fn read_file_format_with_bond_tolerance(
    path: &Path,
    format: FileFormat,
    bond_tolerance: f32,
) -> IoResult<ObjectMolecule> {
    match format {
        FileFormat::Pdb => pdb::read_pdb_with_bond_tolerance(path, bond_tolerance),
        FileFormat::Sdf => sdf::read_sdf(path).and_then(first_molecule),
        FileFormat::Mol2 => mol2::read_mol2(path).and_then(first_molecule),
        FileFormat::Xyz => xyz::read_xyz(path),
        FileFormat::Cif => cif::read_cif_with_bond_tolerance(path, bond_tolerance),
        FileFormat::Bcif => bcif::read_bcif_with_bond_tolerance(path, bond_tolerance),
        FileFormat::Gro => gro::read_gro_with_bond_tolerance(path, bond_tolerance),
        FileFormat::Xtc | FileFormat::Trr => Err(IoError::unsupported(
            "Trajectory-only format; use load_traj instead".to_string(),
        )),
        FileFormat::Ccp4 => Err(IoError::unsupported(
            "CCP4 is a map format; use ccp4::read_ccp4 instead".to_string(),
        )),
        FileFormat::Unknown => Err(IoError::unknown_format(path.to_string_lossy().into_owned())),
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
        FileFormat::Bcif => bcif::read_bcif(path).map(|m| vec![m]),
        FileFormat::Gro => gro::read_gro(path).map(|m| vec![m]),
        FileFormat::Xtc | FileFormat::Trr => Err(IoError::unsupported(
            "Trajectory-only format; use load_traj instead".to_string(),
        )),
        FileFormat::Ccp4 => Err(IoError::unsupported(
            "CCP4 is a map format; use ccp4::read_ccp4 instead".to_string(),
        )),
        FileFormat::Unknown => Err(IoError::unknown_format(path.to_string_lossy().into_owned())),
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

#[cfg(feature = "traj")]
/// Read trajectory frames from a file, auto-detecting the format
pub fn read_trajectory(
    path: &Path,
    opts: &TrajectoryReadOptions,
) -> IoResult<Vec<patinae_mol::CoordSet>> {
    let format = detect::detect_from_path(path);
    read_trajectory_format(path, format, opts)
}

#[cfg(feature = "traj")]
/// Read trajectory frames from a file with a specific format
pub fn read_trajectory_format(
    path: &Path,
    format: FileFormat,
    opts: &TrajectoryReadOptions,
) -> IoResult<Vec<patinae_mol::CoordSet>> {
    let file = compress::open_file(path)?;
    let mut reader = create_trajectory_reader(file, format)?;
    reader.read_frames(opts)
}

/// Parse a molecule from a string with the given format
pub fn parse_str(content: &str, format: FileFormat) -> IoResult<ObjectMolecule> {
    match format {
        FileFormat::Pdb => pdb::read_pdb_str(content),
        FileFormat::Sdf => sdf::read_sdf_str(content).and_then(first_molecule),
        FileFormat::Mol2 => mol2::read_mol2_str(content).and_then(first_molecule),
        FileFormat::Xyz => xyz::read_xyz_str(content),
        FileFormat::Cif => cif::read_cif_str(content),
        FileFormat::Bcif => Err(IoError::unsupported(
            "bCIF is a binary format; use read_file instead".to_string(),
        )),
        FileFormat::Gro => gro::read_gro_str(content),
        FileFormat::Xtc | FileFormat::Trr => Err(IoError::unsupported(
            "Trajectory formats are binary; use load_traj with a file path".to_string(),
        )),
        FileFormat::Ccp4 => Err(IoError::unsupported(
            "CCP4 is a binary map format; use ccp4::read_ccp4 with a file path".to_string(),
        )),
        FileFormat::Unknown => Err(IoError::unknown_format("string input".to_string())),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

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
        assert_eq!(
            FileFormat::from_path(Path::new("test.cif.gz")),
            FileFormat::Cif
        );
        assert_eq!(
            FileFormat::from_path(Path::new("test.sdf.gz")),
            FileFormat::Sdf
        );
        assert_eq!(
            FileFormat::from_path(Path::new("test.mol2.gz")),
            FileFormat::Mol2
        );
        assert_eq!(
            FileFormat::from_path(Path::new("test.xyz.gz")),
            FileFormat::Xyz
        );
        assert_eq!(
            FileFormat::from_path(Path::new("test.gro.gz")),
            FileFormat::Gro
        );
        assert_eq!(
            FileFormat::from_path(Path::new("test.ent.gz")),
            FileFormat::Pdb
        );
    }

    #[test]
    fn parse_str_empty_sdf_returns_empty_file() {
        let result = parse_str("", FileFormat::Sdf);

        assert!(matches!(result, Err(err) if err.is_empty_file()));
    }

    #[test]
    fn parse_str_empty_mol2_returns_empty_file() {
        let result = parse_str("", FileFormat::Mol2);

        assert!(matches!(result, Err(err) if err.is_empty_file()));
    }

    #[test]
    fn read_file_format_empty_sdf_returns_empty_file() {
        let dir = tempdir().expect("temporary directory should be created");
        let path = dir.path().join("empty.sdf");
        std::fs::write(&path, "").expect("empty SDF fixture should be writable");

        let result = read_file_format(&path, FileFormat::Sdf);

        assert!(matches!(result, Err(err) if err.is_empty_file()));
    }

    #[test]
    fn read_file_format_empty_mol2_returns_empty_file() {
        let dir = tempdir().expect("temporary directory should be created");
        let path = dir.path().join("empty.mol2");
        std::fs::write(&path, "").expect("empty MOL2 fixture should be writable");

        let result = read_file_format(&path, FileFormat::Mol2);

        assert!(matches!(result, Err(err) if err.is_empty_file()));
    }

    #[test]
    fn read_file_detects_empty_sdf_as_empty_file() {
        let dir = tempdir().expect("temporary directory should be created");
        let path = dir.path().join("empty.sdf");
        std::fs::write(&path, "").expect("empty SDF fixture should be writable");

        let result = read_file(&path);

        assert!(matches!(result, Err(err) if err.is_empty_file()));
    }

    #[test]
    fn read_file_detects_empty_mol2_as_empty_file() {
        let dir = tempdir().expect("temporary directory should be created");
        let path = dir.path().join("empty.mol2");
        std::fs::write(&path, "").expect("empty MOL2 fixture should be writable");

        let result = read_file(&path);

        assert!(matches!(result, Err(err) if err.is_empty_file()));
    }

    #[test]
    fn read_all_format_empty_sdf_returns_empty_file() {
        let dir = tempdir().expect("temporary directory should be created");
        let path = dir.path().join("empty.sdf");
        std::fs::write(&path, "").expect("empty SDF fixture should be writable");

        let result = read_all_format(&path, FileFormat::Sdf);

        assert!(matches!(result, Err(err) if err.is_empty_file()));
    }

    #[test]
    fn read_all_format_empty_mol2_returns_empty_file() {
        let dir = tempdir().expect("temporary directory should be created");
        let path = dir.path().join("empty.mol2");
        std::fs::write(&path, "").expect("empty MOL2 fixture should be writable");

        let result = read_all_format(&path, FileFormat::Mol2);

        assert!(matches!(result, Err(err) if err.is_empty_file()));
    }
}
