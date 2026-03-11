//! Plugin Backend for Embedded Mode
//!
//! Provides `PluginBackend` — a Python-visible object that the `pymol_rs`
//! package discovers via `sys._pymolrs_backend` when running inside the GUI.
//!
//! Unlike `StandaloneBackend` (which owns a Session), this backend
//! communicates with the host through shared snapshots and a command queue.

use std::collections::HashMap;
use std::ffi::CString;
use std::sync::{Arc, Mutex};

use numpy::ndarray::Array3;
use numpy::{IntoPyArray, PyArray3, PyArrayMethods, PyUntypedArrayMethods};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use pymol_mol::{Atom, AtomIndex, ObjectMolecule, SecondaryStructure};

use crate::commands::{AlterBuffer, AtomChange, PropertyValue};

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
    /// Shared buffer for atom mutations from `alter()` — drained by the handler's mutation queue.
    pub alter_buffer: AlterBuffer,
}

impl SharedState {
    pub fn new(alter_buffer: AlterBuffer) -> Self {
        Self {
            names: Vec::new(),
            molecules: Vec::new(),
            cmd_queue: Vec::new(),
            viewport_image: None,
            set_image_queue: None,
            alter_buffer,
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

    /// Execute a Python expression for each atom matching a selection (read-only).
    ///
    /// Operates on molecule snapshots from the last poll cycle.
    #[pyo3(signature = (selection, expression, space=None))]
    fn iterate(
        &self,
        py: Python<'_>,
        selection: &str,
        expression: &str,
        space: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<()> {
        let globals = match space {
            Some(s) => s.clone(),
            None => PyDict::new(py),
        };
        ensure_builtins(py, &globals)?;

        let code = CString::new(expression).map_err(|_| {
            pyo3::exceptions::PyRuntimeError::new_err("Expression contains null byte")
        })?;

        // Clone molecules to avoid holding lock during Python execution
        let molecules: Vec<(String, ObjectMolecule)> = {
            let state = self.shared.lock().unwrap();
            state.molecules.clone()
        };

        for (name, mol) in &molecules {
            let mask = pymol_select::select(mol, selection).map_err(|e| {
                pyo3::exceptions::PyRuntimeError::new_err(format!("Selection error: {}", e))
            })?;

            let cs = mol.current_coord_set();
            for idx in mask.raw_indices() {
                let atom = mol.get_atom(AtomIndex(idx as u32)).unwrap();
                let coord = cs
                    .and_then(|c| c.get_atom_coord(AtomIndex(idx as u32)))
                    .map(|v| (v.x, v.y, v.z));
                let locals = PyDict::new(py);
                set_atom_locals(&locals, atom, coord, idx, name)?;
                py.run(code.as_c_str(), Some(&globals), Some(&locals))?;
            }
        }
        Ok(())
    }

    /// Execute a Python expression for each atom matching a selection,
    /// allowing modification of atom properties.
    ///
    /// Changes are collected and queued for the host to apply on the next
    /// poll cycle via the handler's mutation queue.
    #[pyo3(signature = (selection, expression, space=None))]
    fn alter(
        &self,
        py: Python<'_>,
        selection: &str,
        expression: &str,
        space: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<()> {
        let globals = match space {
            Some(s) => s.clone(),
            None => PyDict::new(py),
        };
        ensure_builtins(py, &globals)?;

        let code = CString::new(expression).map_err(|_| {
            pyo3::exceptions::PyRuntimeError::new_err("Expression contains null byte")
        })?;

        // Clone molecules to avoid holding lock during Python execution
        let molecules: Vec<(String, ObjectMolecule)> = {
            let state = self.shared.lock().unwrap();
            state.molecules.clone()
        };

        let mut changes: Vec<AtomChange> = Vec::new();

        for (name, mol) in &molecules {
            let mask = pymol_select::select(mol, selection).map_err(|e| {
                pyo3::exceptions::PyRuntimeError::new_err(format!("Selection error: {}", e))
            })?;

            let cs = mol.current_coord_set();
            for idx in mask.raw_indices() {
                let atom = mol.get_atom(AtomIndex(idx as u32)).unwrap();
                let coord = cs
                    .and_then(|c| c.get_atom_coord(AtomIndex(idx as u32)))
                    .map(|v| (v.x, v.y, v.z));
                let locals = PyDict::new(py);
                set_atom_locals(&locals, atom, coord, idx, name)?;
                py.run(code.as_c_str(), Some(&globals), Some(&locals))?;

                // Diff locals vs original atom to collect changes
                let diff = diff_atom_locals(&locals, atom)?;
                if !diff.is_empty() {
                    changes.push(AtomChange {
                        obj: name.clone(),
                        idx: idx as u32,
                        changes: diff,
                    });
                }
            }
        }

        // Push changes to the shared alter buffer
        if !changes.is_empty() {
            let state = self.shared.lock().unwrap();
            let mut buf = state.alter_buffer.lock().unwrap();
            buf.push(changes);
        }

        Ok(())
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

// =============================================================================
// Helpers for iterate/alter
// =============================================================================

/// Populate a Python dict with atom properties.
fn set_atom_locals(
    locals: &Bound<'_, PyDict>,
    atom: &Atom,
    coord: Option<(f32, f32, f32)>,
    index: usize,
    model_name: &str,
) -> PyResult<()> {
    locals.set_item("name", &*atom.name)?;
    locals.set_item("resn", &atom.residue.resn)?;
    locals.set_item("resv", atom.residue.resv)?;
    locals.set_item("resi", atom.residue.resv)?;
    locals.set_item("chain", &atom.residue.chain)?;
    locals.set_item("segi", &atom.residue.segi)?;
    locals.set_item("alt", atom.alt.to_string())?;
    locals.set_item("elem", atom.element.symbol())?;
    locals.set_item("b", atom.b_factor)?;
    locals.set_item("q", atom.occupancy)?;
    locals.set_item("vdw", atom.effective_vdw())?;
    locals.set_item("partial_charge", atom.partial_charge)?;
    locals.set_item("formal_charge", atom.formal_charge)?;
    locals.set_item("ss", ss_to_str(atom.ss_type))?;
    locals.set_item("color", atom.repr.colors.base)?;
    locals.set_item("hetatm", atom.state.hetatm)?;
    locals.set_item(
        "type",
        if atom.state.hetatm { "HETATM" } else { "ATOM" },
    )?;
    locals.set_item("index", index + 1)?;
    locals.set_item("ID", atom.id)?;
    locals.set_item("rank", atom.rank)?;
    locals.set_item("model", model_name)?;
    if let Some((x, y, z)) = coord {
        locals.set_item("x", x)?;
        locals.set_item("y", y)?;
        locals.set_item("z", z)?;
    }
    Ok(())
}

/// Compare the Python locals dict against the original atom and return
/// a map of changed mutable properties.
fn diff_atom_locals(
    locals: &Bound<'_, PyDict>,
    atom: &Atom,
) -> PyResult<HashMap<String, PropertyValue>> {
    let mut changes = HashMap::new();

    if let Some(val) = locals.get_item("name")? {
        let v: String = val.extract()?;
        if &*atom.name != v.as_str() {
            changes.insert("name".into(), PropertyValue::Str(v));
        }
    }
    if let Some(val) = locals.get_item("b")? {
        let v: f32 = val.extract()?;
        if v != atom.b_factor {
            changes.insert("b".into(), PropertyValue::F32(v));
        }
    }
    if let Some(val) = locals.get_item("q")? {
        let v: f32 = val.extract()?;
        if v != atom.occupancy {
            changes.insert("q".into(), PropertyValue::F32(v));
        }
    }
    if let Some(val) = locals.get_item("vdw")? {
        let v: f32 = val.extract()?;
        if v != atom.effective_vdw() {
            changes.insert("vdw".into(), PropertyValue::F32(v));
        }
    }
    if let Some(val) = locals.get_item("partial_charge")? {
        let v: f32 = val.extract()?;
        if v != atom.partial_charge {
            changes.insert("partial_charge".into(), PropertyValue::F32(v));
        }
    }
    if let Some(val) = locals.get_item("formal_charge")? {
        let v: i8 = val.extract()?;
        if v != atom.formal_charge {
            changes.insert("formal_charge".into(), PropertyValue::I8(v));
        }
    }
    if let Some(val) = locals.get_item("color")? {
        let v: i32 = val.extract()?;
        if v != atom.repr.colors.base {
            changes.insert("color".into(), PropertyValue::I32(v));
        }
    }
    if let Some(val) = locals.get_item("elem")? {
        let v: String = val.extract()?;
        if v != atom.element.symbol() {
            changes.insert("elem".into(), PropertyValue::Str(v));
        }
    }
    if let Some(val) = locals.get_item("ss")? {
        let v: String = val.extract()?;
        if v != ss_to_str(atom.ss_type) {
            changes.insert("ss".into(), PropertyValue::Str(v));
        }
    }
    if let Some(val) = locals.get_item("type")? {
        let v: String = val.extract()?;
        let expected = if atom.state.hetatm { "HETATM" } else { "ATOM" };
        if v != expected {
            changes.insert("type".into(), PropertyValue::Str(v));
        }
    }
    if let Some(val) = locals.get_item("alt")? {
        let v: String = val.extract()?;
        let c = v.chars().next().unwrap_or(' ');
        if c != atom.alt {
            changes.insert("alt".into(), PropertyValue::Str(v));
        }
    }
    if let Some(val) = locals.get_item("chain")? {
        let v: String = val.extract()?;
        if v != atom.residue.chain {
            changes.insert("chain".into(), PropertyValue::Str(v));
        }
    }
    if let Some(val) = locals.get_item("resn")? {
        let v: String = val.extract()?;
        if v != atom.residue.resn {
            changes.insert("resn".into(), PropertyValue::Str(v));
        }
    }
    if let Some(val) = locals.get_item("resv")? {
        let v: i32 = val.extract()?;
        if v != atom.residue.resv {
            changes.insert("resv".into(), PropertyValue::I32(v));
        }
    }
    if let Some(val) = locals.get_item("segi")? {
        let v: String = val.extract()?;
        if v != atom.residue.segi {
            changes.insert("segi".into(), PropertyValue::Str(v));
        }
    }

    Ok(changes)
}

fn ss_to_str(ss: SecondaryStructure) -> &'static str {
    match ss {
        SecondaryStructure::Helix | SecondaryStructure::Helix310 | SecondaryStructure::HelixPi => {
            "H"
        }
        SecondaryStructure::Sheet => "S",
        _ => "L",
    }
}

fn ensure_builtins(py: Python<'_>, globals: &Bound<'_, PyDict>) -> PyResult<()> {
    if !globals.contains("__builtins__")? {
        let builtins = py.import("builtins")?;
        globals.set_item("__builtins__", builtins)?;
    }
    Ok(())
}
