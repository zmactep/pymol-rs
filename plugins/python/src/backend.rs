//! Plugin Backend for Embedded Mode
//!
//! Provides `PluginBackend` — a Python-visible object that the `pymol_rs`
//! package discovers via `sys._pymolrs_backend` when running inside the GUI.
//!
//! Unlike `StandaloneBackend` (which owns a Session), this backend
//! communicates with the host through shared snapshots and a command queue.

use std::sync::{Arc, Mutex};

use numpy::ndarray::Array3;
use numpy::{IntoPyArray, PyArray3, PyArrayMethods, PyUntypedArrayMethods};
use pyo3::prelude::*;
use pymol_mol::ObjectMolecule;

/// Shared state between the host (poll) and the Python backend.
pub struct SharedState {
    /// Object names currently loaded in the viewer.
    pub names: Vec<String>,
    /// Snapshot molecules (cloned from the viewer's registry).
    pub molecules: Vec<(String, ObjectMolecule)>,
    /// Queued commands to execute on the host side.
    pub cmd_queue: Vec<(String, bool)>,
    /// Viewport image snapshot for reading (RGBA data, width, height).
    pub viewport_image: Option<(Vec<u8>, u32, u32)>,
    /// Pending viewport image to set: `Some(Some(...))` = set, `Some(None)` = clear.
    pub set_image_queue: Option<Option<(Vec<u8>, u32, u32)>>,
}

impl SharedState {
    pub fn new() -> Self {
        Self {
            names: Vec::new(),
            molecules: Vec::new(),
            cmd_queue: Vec::new(),
            viewport_image: None,
            set_image_queue: None,
        }
    }
}

/// Thread-safe handle to shared state.
pub type SharedStateHandle = Arc<Mutex<SharedState>>;

/// Python-visible backend for embedded (plugin) mode.
///
/// Injected into `sys._pymolrs_backend` so that `pymol_rs.__init__`
/// auto-detects embedded mode and wraps this in a `Cmd` object.
#[pyclass(name = "PluginBackend")]
pub struct PluginBackend {
    shared: SharedStateHandle,
}

impl PluginBackend {
    pub fn new(shared: SharedStateHandle) -> Self {
        Self { shared }
    }
}

#[pymethods]
impl PluginBackend {
    /// Queue a command for execution by the host.
    ///
    /// The command is not executed immediately — it is drained during
    /// the next `poll()` cycle by the message handler.
    #[pyo3(signature = (command, silent=false))]
    fn execute(&self, command: &str, silent: bool) -> PyResult<()> {
        let mut state = self.shared.lock().unwrap();
        state.cmd_queue.push((command.to_string(), silent));
        Ok(())
    }

    /// Get a molecular object by name from the snapshot.
    fn get_model(&self, py: Python<'_>, name: &str) -> PyResult<Py<PyAny>> {
        let state = self.shared.lock().unwrap();
        for (n, mol) in &state.molecules {
            if n == name {
                // Return a basic dict with molecular data.
                // The pymol_rs package can wrap it if needed.
                let dict = pyo3::types::PyDict::new(py);
                dict.set_item("name", &mol.name)?;
                dict.set_item("atom_count", mol.atom_count())?;
                dict.set_item("bond_count", mol.bond_count())?;
                dict.set_item("state_count", mol.state_count())?;
                return Ok(dict.into_any().unbind());
            }
        }
        Err(pyo3::exceptions::PyKeyError::new_err(format!(
            "object '{}' not found",
            name
        )))
    }

    /// Get list of loaded object names.
    fn get_names(&self) -> Vec<String> {
        let state = self.shared.lock().unwrap();
        state.names.clone()
    }

    /// Count atoms matching a selection across all molecules.
    fn count_atoms(&self, selection: &str) -> PyResult<usize> {
        let state = self.shared.lock().unwrap();
        let mut count = 0;
        for (_, mol) in &state.molecules {
            if let Ok(result) = pymol_select::select(mol, selection) {
                count += result.count();
            }
        }
        Ok(count)
    }

    /// Get the current viewport image as a numpy array (H, W, 4) uint8.
    ///
    /// Returns `None` if no viewport image is set.
    /// The image is a snapshot from the last poll cycle.
    fn get_viewport_image<'py>(&self, py: Python<'py>) -> PyResult<Option<Py<PyAny>>> {
        let state = self.shared.lock().unwrap();
        match &state.viewport_image {
            Some((data, width, height)) => {
                let h = *height as usize;
                let w = *width as usize;
                let arr = Array3::from_shape_vec((h, w, 4), data.clone()).map_err(|e| {
                    pyo3::exceptions::PyRuntimeError::new_err(format!(
                        "Failed to create array: {}",
                        e
                    ))
                })?;
                Ok(Some(arr.into_pyarray(py).into_any().unbind()))
            }
            None => Ok(None),
        }
    }

    /// Set the viewport image from a numpy array (H, W, 4) uint8.
    ///
    /// The image is queued and applied on the next poll cycle.
    fn set_viewport_image(&self, array: &Bound<'_, PyAny>) -> PyResult<()> {
        let arr = array.cast::<PyArray3<u8>>().map_err(|_| {
            pyo3::exceptions::PyTypeError::new_err(
                "Expected numpy array with shape (H, W, 4) and dtype uint8",
            )
        })?;
        let shape = arr.shape();
        if shape.len() != 3 || shape[2] != 4 {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "Expected shape (H, W, 4), got ({}, {}, {})",
                shape[0], shape[1], shape[2]
            )));
        }
        let height = shape[0] as u32;
        let width = shape[1] as u32;
        let data = arr.to_vec()?;
        let mut state = self.shared.lock().unwrap();
        state.set_image_queue = Some(Some((data, width, height)));
        Ok(())
    }

    /// Clear the viewport image overlay.
    ///
    /// The clear is queued and applied on the next poll cycle.
    fn clear_viewport_image(&self) {
        let mut state = self.shared.lock().unwrap();
        state.set_image_queue = Some(None);
    }
}
