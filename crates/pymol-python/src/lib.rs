//! PyMOL-RS Python Bindings
//!
//! This crate provides Python bindings for pymol-rs using PyO3.

use pyo3::prelude::*;
use pyo3::types::PyModule;

mod cmd;
mod convert;
mod error;
mod viewer;

pub mod color;
pub mod io;
pub mod mol;
pub mod selecting;
pub mod settings;

pub use cmd::PyCmd;
pub use error::{PymolError, SelectionError};
pub use mol::{PyAtom, PyBond, PyCoordSet, PyElement, PyObjectMolecule};
pub use color::PyColor;
pub use selecting::PySelectionResult;

/// Python module initialization
#[pymodule]
fn _pymol(py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Initialize logging
    let _ = env_logger::try_init();

    // Register exception types
    m.add("PymolError", py.get_type_bound::<error::PymolError>())?;
    m.add("SelectionError", py.get_type_bound::<error::SelectionError>())?;

    // Register main types
    m.add_class::<PyCmd>()?;
    m.add_class::<PyObjectMolecule>()?;
    m.add_class::<mol::PyAtom>()?;
    m.add_class::<mol::PyBond>()?;
    m.add_class::<mol::PyElement>()?;
    m.add_class::<mol::PyCoordSet>()?;
    m.add_class::<color::PyColor>()?;
    m.add_class::<selecting::PySelectionResult>()?;

    // Create and register the global cmd instance
    let cmd = PyCmd::new(true)?;
    m.add("cmd", cmd)?;

    // Register submodules
    register_submodule(py, m, "mol", mol::register_module)?;
    register_submodule(py, m, "io", io::register_module)?;
    register_submodule(py, m, "selecting", selecting::register_module)?;
    register_submodule(py, m, "color", color::register_module)?;
    register_submodule(py, m, "settings", settings::register_module)?;

    // Add version info
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    Ok(())
}

/// Helper function to register a submodule
fn register_submodule<F>(
    py: Python<'_>,
    parent: &Bound<'_, PyModule>,
    name: &str,
    register_fn: F,
) -> PyResult<()>
where
    F: FnOnce(&Bound<'_, PyModule>) -> PyResult<()>,
{
    let submodule = PyModule::new_bound(py, name)?;
    register_fn(&submodule)?;
    parent.add_submodule(&submodule)?;
    Ok(())
}
