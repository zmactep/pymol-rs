//! File I/O module for Python
//!
//! Provides functions for reading and writing molecular files.

use pyo3::prelude::*;
use std::path::Path;

use crate::error::ResultExt;
use crate::mol::PyObjectMolecule;

/// Supported file formats
#[pyclass(name = "FileFormat", eq, eq_int)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PyFileFormat {
    /// PDB format
    Pdb,
    /// SDF/MOL format
    Sdf,
    /// MOL2 (Tripos) format
    Mol2,
    /// mmCIF format
    Cif,
    /// XYZ format
    Xyz,
    /// BinaryCIF format
    Bcif,
    /// GROMACS GRO format
    Gro,
    /// Unknown format
    Unknown,
}

impl From<pymol_io::FileFormat> for PyFileFormat {
    fn from(fmt: pymol_io::FileFormat) -> Self {
        match fmt {
            pymol_io::FileFormat::Pdb => PyFileFormat::Pdb,
            pymol_io::FileFormat::Sdf => PyFileFormat::Sdf,
            pymol_io::FileFormat::Mol2 => PyFileFormat::Mol2,
            pymol_io::FileFormat::Cif => PyFileFormat::Cif,
            pymol_io::FileFormat::Bcif => PyFileFormat::Bcif,
            pymol_io::FileFormat::Xyz => PyFileFormat::Xyz,
            pymol_io::FileFormat::Gro => PyFileFormat::Gro,
            pymol_io::FileFormat::Unknown => PyFileFormat::Unknown,
        }
    }
}

impl From<PyFileFormat> for pymol_io::FileFormat {
    fn from(fmt: PyFileFormat) -> Self {
        match fmt {
            PyFileFormat::Pdb => pymol_io::FileFormat::Pdb,
            PyFileFormat::Sdf => pymol_io::FileFormat::Sdf,
            PyFileFormat::Mol2 => pymol_io::FileFormat::Mol2,
            PyFileFormat::Cif => pymol_io::FileFormat::Cif,
            PyFileFormat::Bcif => pymol_io::FileFormat::Bcif,
            PyFileFormat::Xyz => pymol_io::FileFormat::Xyz,
            PyFileFormat::Gro => pymol_io::FileFormat::Gro,
            PyFileFormat::Unknown => pymol_io::FileFormat::Unknown,
        }
    }
}

#[pymethods]
impl PyFileFormat {
    fn __repr__(&self) -> &'static str {
        match self {
            PyFileFormat::Pdb => "FileFormat.Pdb",
            PyFileFormat::Sdf => "FileFormat.Sdf",
            PyFileFormat::Mol2 => "FileFormat.Mol2",
            PyFileFormat::Cif => "FileFormat.Cif",
            PyFileFormat::Bcif => "FileFormat.Bcif",
            PyFileFormat::Xyz => "FileFormat.Xyz",
            PyFileFormat::Gro => "FileFormat.Gro",
            PyFileFormat::Unknown => "FileFormat.Unknown",
        }
    }
}

/// Read a molecular file
///
/// Automatically detects the format from the file extension.
///
/// Args:
///     path: Path to the file
///     format: Optional format override (default: auto-detect)
///
/// Returns:
///     ObjectMolecule containing the molecular data
#[pyfunction]
#[pyo3(signature = (path, format=None))]
pub fn read_file(path: &str, format: Option<PyFileFormat>) -> PyResult<PyObjectMolecule> {
    let path = Path::new(path);
    
    let mol = if let Some(fmt) = format {
        pymol_io::read_file_format(path, fmt.into()).map_py_err()?
    } else {
        pymol_io::read_file(path).map_py_err()?
    };
    
    Ok(PyObjectMolecule::from(mol))
}

/// Read all molecules from a file
///
/// For multi-molecule formats (SDF, MOL2), returns all molecules.
/// For single-molecule formats, returns a list with one molecule.
///
/// Args:
///     path: Path to the file
///     format: Optional format override (default: auto-detect)
///
/// Returns:
///     List of ObjectMolecule objects
#[pyfunction]
#[pyo3(signature = (path, format=None))]
pub fn read_all(path: &str, format: Option<PyFileFormat>) -> PyResult<Vec<PyObjectMolecule>> {
    let path = Path::new(path);
    
    let mols = if let Some(fmt) = format {
        pymol_io::read_all_format(path, fmt.into()).map_py_err()?
    } else {
        pymol_io::read_all(path).map_py_err()?
    };
    
    Ok(mols.into_iter().map(PyObjectMolecule::from).collect())
}

/// Write a molecule to a file
///
/// Args:
///     path: Path to the output file
///     mol: The molecule to write
///     format: Optional format override (default: auto-detect from extension)
#[pyfunction]
#[pyo3(signature = (path, mol, format=None))]
pub fn write_file(path: &str, mol: &PyObjectMolecule, format: Option<PyFileFormat>) -> PyResult<()> {
    let path = Path::new(path);
    
    if let Some(fmt) = format {
        pymol_io::write_file_format(path, mol.inner(), fmt.into()).map_py_err()?;
    } else {
        pymol_io::write_file(path, mol.inner()).map_py_err()?;
    }
    
    Ok(())
}

/// Parse a molecule from a string
///
/// Args:
///     content: String content in the specified format
///     format: File format (required, cannot auto-detect from string)
///
/// Returns:
///     ObjectMolecule containing the molecular data
#[pyfunction]
pub fn parse_string(content: &str, format: PyFileFormat) -> PyResult<PyObjectMolecule> {
    let mol = pymol_io::parse_str(content, format.into()).map_py_err()?;
    Ok(PyObjectMolecule::from(mol))
}

/// Detect file format from path
///
/// Args:
///     path: Path to the file
///
/// Returns:
///     Detected FileFormat
#[pyfunction]
pub fn detect_format(path: &str) -> PyFileFormat {
    let path = Path::new(path);
    pymol_io::detect::detect_from_path(path).into()
}

/// Fetch a structure from RCSB PDB
///
/// Args:
///     pdb_id: 4-letter PDB ID (e.g., "1crn", "4hhb")
///     format: Download format ("pdb", "cif", or "bcif", default: "bcif")
///
/// Returns:
///     ObjectMolecule containing the fetched structure
#[pyfunction]
#[pyo3(signature = (pdb_id, format="bcif"))]
pub fn fetch(pdb_id: &str, format: &str) -> PyResult<PyObjectMolecule> {
    let fetch_format = match format.to_lowercase().as_str() {
        "pdb" => pymol_io::FetchFormat::Pdb,
        "cif" | "mmcif" => pymol_io::FetchFormat::Cif,
        "bcif" | "binarycif" => pymol_io::FetchFormat::Bcif,
        _ => {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "Unknown fetch format: '{}'. Use 'pdb', 'cif', or 'bcif'.",
                format
            )))
        }
    };
    
    let mol = pymol_io::fetch(pdb_id, fetch_format).map_py_err()?;
    Ok(PyObjectMolecule::from(mol))
}

/// Register the io submodule
pub fn register_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyFileFormat>()?;
    m.add_function(wrap_pyfunction!(read_file, m)?)?;
    m.add_function(wrap_pyfunction!(read_all, m)?)?;
    m.add_function(wrap_pyfunction!(write_file, m)?)?;
    m.add_function(wrap_pyfunction!(parse_string, m)?)?;
    m.add_function(wrap_pyfunction!(detect_format, m)?)?;
    m.add_function(wrap_pyfunction!(fetch, m)?)?;
    Ok(())
}
