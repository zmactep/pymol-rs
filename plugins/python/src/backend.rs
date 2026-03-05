//! Plugin Backend for Embedded Mode
//!
//! Provides `PluginBackend` — a Python-visible object that the `pymol_rs`
//! package discovers via `sys._pymolrs_backend` when running inside the GUI.
//!
//! Unlike `StandaloneBackend` (which owns a Session), this backend
//! communicates with the host through shared snapshots and a command queue.

use std::sync::{Arc, Mutex};

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
}

impl SharedState {
    pub fn new() -> Self {
        Self {
            names: Vec::new(),
            molecules: Vec::new(),
            cmd_queue: Vec::new(),
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
}
