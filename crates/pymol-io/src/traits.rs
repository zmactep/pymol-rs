//! Reader and writer traits for molecular file formats
//!
//! Defines the common interface for reading and writing molecular structures.

use std::io::{Read, Write};
use std::path::Path;

use pymol_mol::ObjectMolecule;

use crate::error::{IoError, IoResult};

/// Supported file formats
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FileFormat {
    /// Protein Data Bank format
    Pdb,
    /// MDL SDF/MOL format (V2000)
    Sdf,
    /// TRIPOS MOL2 format
    Mol2,
    /// Macromolecular Crystallographic Information File
    Cif,
    /// XYZ coordinate format
    Xyz,
    /// Unknown format
    Unknown,
}

impl FileFormat {
    /// Get the file format from a file extension
    pub fn from_extension(ext: &str) -> Self {
        match ext.to_lowercase().as_str() {
            "pdb" | "ent" => FileFormat::Pdb,
            "pdb.gz" | "ent.gz" => FileFormat::Pdb,
            "sdf" | "mol" | "sd" => FileFormat::Sdf,
            "mol2" | "ml2" => FileFormat::Mol2,
            "cif" | "mmcif" => FileFormat::Cif,
            "xyz" => FileFormat::Xyz,
            _ => FileFormat::Unknown,
        }
    }

    /// Get the file format from a path
    pub fn from_path(path: &Path) -> Self {
        // Check for double extensions like .pdb.gz
        let filename = path.file_name().and_then(|s| s.to_str()).unwrap_or("");
        if filename.ends_with(".pdb.gz") || filename.ends_with(".ent.gz") {
            return FileFormat::Pdb;
        }

        path.extension()
            .and_then(|s| s.to_str())
            .map(FileFormat::from_extension)
            .unwrap_or(FileFormat::Unknown)
    }

    /// Get the default file extension for this format
    pub fn extension(&self) -> &'static str {
        match self {
            FileFormat::Pdb => "pdb",
            FileFormat::Sdf => "sdf",
            FileFormat::Mol2 => "mol2",
            FileFormat::Cif => "cif",
            FileFormat::Xyz => "xyz",
            FileFormat::Unknown => "",
        }
    }

    /// Get a human-readable name for the format
    pub fn name(&self) -> &'static str {
        match self {
            FileFormat::Pdb => "PDB",
            FileFormat::Sdf => "SDF/MOL",
            FileFormat::Mol2 => "MOL2",
            FileFormat::Cif => "mmCIF",
            FileFormat::Xyz => "XYZ",
            FileFormat::Unknown => "Unknown",
        }
    }
}

/// Trait for reading molecular structures from a source
pub trait MoleculeReader {
    /// Read a single molecule from the source
    ///
    /// For formats that support multiple molecules (like SDF), this returns
    /// the first molecule. Use `read_all` to get all molecules.
    fn read(&mut self) -> IoResult<ObjectMolecule>;

    /// Read all molecules from the source
    ///
    /// For single-molecule formats like PDB, this returns a vector with one element.
    /// For multi-molecule formats like SDF, this returns all molecules.
    fn read_all(&mut self) -> IoResult<Vec<ObjectMolecule>> {
        let mut molecules = Vec::new();
        loop {
            match self.read() {
                Ok(mol) => molecules.push(mol),
                Err(IoError::EmptyFile) if !molecules.is_empty() => break,
                Err(e) if molecules.is_empty() => return Err(e),
                Err(_) => break,
            }
        }
        Ok(molecules)
    }
}

/// Trait for writing molecular structures to a destination
pub trait MoleculeWriter {
    /// Write a single molecule to the destination
    fn write(&mut self, mol: &ObjectMolecule) -> IoResult<()>;

    /// Write multiple molecules to the destination
    ///
    /// For formats that don't support multiple molecules, this writes
    /// only the first molecule.
    fn write_all(&mut self, molecules: &[ObjectMolecule]) -> IoResult<()> {
        for mol in molecules {
            self.write(mol)?;
        }
        Ok(())
    }

    /// Flush any buffered data
    fn flush(&mut self) -> IoResult<()>;
}

/// Options for reading molecular files
#[derive(Debug, Clone, Default)]
pub struct ReadOptions {
    /// Whether to read bond information (if available in format)
    pub read_bonds: bool,
    /// Whether to infer bonds from distances (for formats without bonds)
    pub infer_bonds: bool,
    /// Whether to read secondary structure annotations
    pub read_secondary_structure: bool,
    /// Which model/state to read (None = all models)
    pub model: Option<usize>,
}

impl ReadOptions {
    /// Create default read options
    pub fn new() -> Self {
        ReadOptions {
            read_bonds: true,
            read_secondary_structure: true,
            ..Default::default()
        }
    }

    /// Set whether to read bonds
    pub fn with_bonds(mut self, read_bonds: bool) -> Self {
        self.read_bonds = read_bonds;
        self
    }

    /// Set whether to infer bonds
    pub fn with_infer_bonds(mut self, infer_bonds: bool) -> Self {
        self.infer_bonds = infer_bonds;
        self
    }

    /// Set which model to read
    pub fn with_model(mut self, model: usize) -> Self {
        self.model = Some(model);
        self
    }
}

/// Options for writing molecular files
#[derive(Debug, Clone, Default)]
pub struct WriteOptions {
    /// Whether to write bond information (if supported by format)
    pub write_bonds: bool,
    /// Whether to write secondary structure annotations
    pub write_secondary_structure: bool,
    /// Which model/state to write (None = all states)
    pub state: Option<usize>,
    /// Title/header for the file
    pub title: Option<String>,
}

impl WriteOptions {
    /// Create default write options
    pub fn new() -> Self {
        WriteOptions {
            write_bonds: true,
            write_secondary_structure: true,
            ..Default::default()
        }
    }

    /// Set whether to write bonds
    pub fn with_bonds(mut self, write_bonds: bool) -> Self {
        self.write_bonds = write_bonds;
        self
    }

    /// Set which state to write
    pub fn with_state(mut self, state: usize) -> Self {
        self.state = Some(state);
        self
    }

    /// Set the title
    pub fn with_title(mut self, title: impl Into<String>) -> Self {
        self.title = Some(title.into());
        self
    }
}

/// Create a reader for the given format from a Read source
pub fn create_reader<R: Read + 'static>(
    reader: R,
    format: FileFormat,
) -> IoResult<Box<dyn MoleculeReader>> {
    match format {
        FileFormat::Pdb => Ok(Box::new(crate::pdb::PdbReader::new(reader))),
        FileFormat::Sdf => Ok(Box::new(crate::sdf::SdfReader::new(reader))),
        FileFormat::Mol2 => Ok(Box::new(crate::mol2::Mol2Reader::new(reader))),
        FileFormat::Xyz => Ok(Box::new(crate::xyz::XyzReader::new(reader))),
        FileFormat::Cif => Ok(Box::new(crate::cif::CifReader::new(reader))),
        FileFormat::Unknown => Err(IoError::UnknownFormat("Unknown format".to_string())),
    }
}

/// Create a writer for the given format to a Write destination
pub fn create_writer<W: Write + 'static>(
    writer: W,
    format: FileFormat,
) -> IoResult<Box<dyn MoleculeWriter>> {
    match format {
        FileFormat::Pdb => Ok(Box::new(crate::pdb::PdbWriter::new(writer))),
        FileFormat::Sdf => Ok(Box::new(crate::sdf::SdfWriter::new(writer))),
        FileFormat::Mol2 => Ok(Box::new(crate::mol2::Mol2Writer::new(writer))),
        FileFormat::Xyz => Ok(Box::new(crate::xyz::XyzWriter::new(writer))),
        FileFormat::Cif => Ok(Box::new(crate::cif::CifWriter::new(writer))),
        FileFormat::Unknown => Err(IoError::UnknownFormat("Unknown format".to_string())),
    }
}
