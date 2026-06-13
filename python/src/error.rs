//! Error mapping from Rust errors to Python exceptions
//!
//! This module provides conversion from various patinae error types
//! to Python exceptions using helper traits.

use pyo3::create_exception;
use pyo3::exceptions::{PyIOError, PyKeyError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;

// Define custom exceptions
create_exception!(patinae, PatinaeError, pyo3::exceptions::PyException);
create_exception!(patinae, SelectionError, PatinaeError);

/// Trait for converting errors to PyErr
pub trait IntoPyErr {
    fn into_py_err(self) -> PyErr;
}

impl IntoPyErr for patinae_io::IoError {
    fn into_py_err(self) -> PyErr {
        PyIOError::new_err(self.to_string())
    }
}

impl IntoPyErr for patinae_mol::MolError {
    fn into_py_err(self) -> PyErr {
        PyValueError::new_err(self.to_string())
    }
}

impl IntoPyErr for patinae_select::SelectError {
    fn into_py_err(self) -> PyErr {
        SelectionError::new_err(self.to_string())
    }
}

impl IntoPyErr for patinae_cmd::CmdError {
    fn into_py_err(self) -> PyErr {
        PyRuntimeError::new_err(self.to_string())
    }
}

impl IntoPyErr for patinae_scene::ViewerError {
    fn into_py_err(self) -> PyErr {
        PyRuntimeError::new_err(self.to_string())
    }
}

impl IntoPyErr for patinae_scene::SceneError {
    fn into_py_err(self) -> PyErr {
        PatinaeError::new_err(self.to_string())
    }
}

impl IntoPyErr for patinae_settings::SettingError {
    fn into_py_err(self) -> PyErr {
        PyKeyError::new_err(self.to_string())
    }
}

/// Extension trait to convert Result<T, E> where E: IntoPyErr
pub trait ResultExt<T> {
    fn map_py_err(self) -> PyResult<T>;
}

impl<T, E: IntoPyErr> ResultExt<T> for Result<T, E> {
    fn map_py_err(self) -> PyResult<T> {
        self.map_err(|e| e.into_py_err())
    }
}
