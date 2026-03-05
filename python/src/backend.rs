//! Standalone backend for the pymol_rs Python package.
//!
//! Provides direct access to a Session + CommandExecutor without IPC.
//! Used when running Python scripts standalone (not embedded in the GUI).

use pyo3::exceptions::{PyKeyError, PyRuntimeError};
use pyo3::prelude::*;

use pymol_cmd::CommandExecutor;
use pymol_scene::{Session, SessionAdapter};
use pymol_select::select;

use crate::mol::PyObjectMolecule;

/// Standalone backend — owns a Session and CommandExecutor.
///
/// All commands execute directly against the local session.
/// No IPC, no GUI — pure headless operation.
#[pyclass]
pub struct StandaloneBackend {
    session: Session,
    executor: CommandExecutor,
    needs_redraw: bool,
}

impl StandaloneBackend {
    pub fn create() -> Self {
        let mut session = Session::new();
        session.apply_default_settings();
        Self {
            session,
            executor: CommandExecutor::new(),
            needs_redraw: false,
        }
    }
}

#[pymethods]
impl StandaloneBackend {
    #[new]
    fn new() -> Self {
        let mut session = Session::new();
        session.apply_default_settings();
        Self {
            session,
            executor: CommandExecutor::new(),
            needs_redraw: false,
        }
    }

    /// Execute a PyMOL command string.
    #[pyo3(signature = (command, silent=false))]
    fn execute(&mut self, command: &str, silent: bool) -> PyResult<()> {
        let mut adapter = SessionAdapter {
            session: &mut self.session,
            render_context: None,
            default_size: (1024, 768),
            needs_redraw: &mut self.needs_redraw,
            async_fetch_fn: None,
        };
        self.executor
            .do_with_options(&mut adapter, command, true, silent)
            .map(|_| ())
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))
    }

    /// Get a molecular object by name.
    fn get_model(&self, name: &str) -> PyResult<PyObjectMolecule> {
        let mol_obj = self
            .session
            .registry
            .get_molecule(name)
            .ok_or_else(|| PyKeyError::new_err(format!("Object '{}' not found", name)))?;
        Ok(mol_obj.molecule().clone().into())
    }

    /// Get the names of all loaded objects.
    fn get_names(&self) -> Vec<String> {
        self.session
            .registry
            .names()
            .map(|s| s.to_string())
            .collect()
    }

    /// Count atoms matching a selection expression.
    #[pyo3(signature = (selection="all"))]
    fn count_atoms(&self, selection: &str) -> PyResult<usize> {
        let names: Vec<String> = self
            .session
            .registry
            .names()
            .map(|s| s.to_string())
            .collect();
        let mut total = 0;
        for name in &names {
            if let Some(mol_obj) = self.session.registry.get_molecule(name) {
                let mol = mol_obj.molecule();
                let result = select(mol, selection);
                match result {
                    Ok(mask) => total += mask.count(),
                    Err(e) => {
                        return Err(PyRuntimeError::new_err(format!(
                            "Selection error: {}",
                            e
                        )))
                    }
                }
            }
        }
        Ok(total)
    }
}
