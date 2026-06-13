use std::collections::HashMap;
use std::ffi::CString;
use std::sync::{Arc, Mutex};

use pyo3::prelude::*;
use pyo3::types::PyDict;

/// A single atom property value that can be changed by `alter`.
#[derive(Debug, Clone)]
pub enum PropertyValue {
    Str(String),
    F32(f32),
    I32(i32),
    I8(i8),
    Bool(bool),
}

/// A batch of property changes for a single atom.
#[derive(Debug, Clone)]
pub struct AtomChange {
    pub obj: String,
    pub idx: u32,
    pub changes: HashMap<String, PropertyValue>,
}

/// Thread-safe buffer for atom changes produced by `alter()` on the Python
/// worker thread and consumed via `queue_viewer_mutation` on the main thread.
pub type AlterBuffer = Arc<Mutex<Vec<Vec<AtomChange>>>>;

/// Build the globals dict for `iterate`/`alter` expression evaluation.
pub fn build_globals<'py>(
    py: Python<'py>,
    space: Option<&Bound<'py, PyDict>>,
) -> PyResult<Bound<'py, PyDict>> {
    let globals = py.import("__main__")?.dict();
    if let Some(s) = space {
        globals.update(s.as_mapping())?;
    }
    Ok(globals)
}

pub fn expression_to_cstring(expression: &str) -> PyResult<CString> {
    CString::new(expression)
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("Expression contains null byte"))
}
