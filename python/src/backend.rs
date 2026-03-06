//! Standalone backend for the pymol_rs Python package.
//!
//! Provides direct access to a Session + CommandExecutor without IPC.
//! Used when running Python scripts standalone (not embedded in the GUI).

use pyo3::exceptions::{PyKeyError, PyRuntimeError};
use pyo3::prelude::*;

use pymol_cmd::CommandExecutor;
use pymol_scene::{Session, SessionAdapter, ViewportImage};
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

    /// Get the current viewport image as a numpy array (H, W, 4) uint8.
    ///
    /// Returns `None` if no viewport image is set.
    #[cfg(feature = "numpy")]
    fn get_viewport_image<'py>(&self, py: Python<'py>) -> PyResult<Option<Py<PyAny>>> {
        use numpy::ndarray::Array3;
        use numpy::IntoPyArray;

        match &self.session.viewport_image {
            Some(img) => {
                let h = img.height as usize;
                let w = img.width as usize;
                let arr =
                    Array3::from_shape_vec((h, w, 4), img.data.clone()).map_err(|e| {
                        PyRuntimeError::new_err(format!("Failed to create array: {}", e))
                    })?;
                Ok(Some(arr.into_pyarray(py).into_any().unbind()))
            }
            None => Ok(None),
        }
    }

    /// Set the viewport image from a numpy array (H, W, 4) uint8.
    #[cfg(feature = "numpy")]
    fn set_viewport_image(&mut self, array: &Bound<'_, PyAny>) -> PyResult<()> {
        use numpy::{PyArray3, PyArrayMethods, PyUntypedArrayMethods};

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
        self.session.viewport_image = Some(ViewportImage {
            data,
            width,
            height,
        });
        Ok(())
    }

    /// Clear the viewport image.
    fn clear_viewport_image(&mut self) {
        self.session.viewport_image = None;
    }
}
