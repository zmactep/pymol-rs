//! PyCmd - Main command interface for Python
//!
//! This module provides the `cmd` object that users interact with,
//! following PyMOL's familiar API pattern.

use std::path::Path;

use pyo3::prelude::*;
use pyo3::types::PyAny;

use pymol_cmd::CommandExecutor;
use pymol_mol::RepMask;

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;

use parking_lot::RwLock;

use crate::error::ResultExt;
use crate::settings::{get_setting_py, set_setting_py};
use crate::viewer::{HeadlessViewer, SharedViewerState, ViewerHandle, ViewerMode};

use pymol_scene::Viewer as GuiViewer;

/// Main command interface
///
/// Provides the familiar PyMOL `cmd.function()` API for molecular visualization.
///
/// Example:
///     from pymol import cmd
///     cmd.load("protein.pdb")
///     cmd.show("cartoon")
///     cmd.color("green", "chain A")
///     cmd.zoom()
#[pyclass(name = "Cmd")]
pub struct PyCmd {
    /// Viewer state (headless or GUI)
    viewer: ViewerMode,
    /// Command executor for text commands (reserved for future use)
    #[allow(dead_code)]
    executor: CommandExecutor,
    /// Handle to a running GUI viewer (if started with show_viewer)
    viewer_handle: Option<ViewerHandle>,
    /// Shared state with the GUI viewer
    shared_state: Option<Arc<RwLock<SharedViewerState>>>,
}

impl PyCmd {
    /// Create a new PyCmd instance
    pub fn new(headless: bool) -> PyResult<Self> {
        let viewer = if headless {
            ViewerMode::Headless(HeadlessViewer::new())
        } else {
            // For now, always use headless mode
            // GUI mode requires event loop integration
            ViewerMode::Headless(HeadlessViewer::new())
        };

        Ok(PyCmd {
            viewer,
            executor: CommandExecutor::new(),
            viewer_handle: None,
            shared_state: None,
        })
    }

    /// Get the shared state if a GUI viewer is running, otherwise return None
    ///
    /// Note: This is intended for future use when live sync between Python
    /// commands and the viewer is implemented.
    #[allow(dead_code)]
    fn get_shared_state(&self) -> Option<&Arc<RwLock<SharedViewerState>>> {
        if self.viewer_handle.as_ref().map(|h| h.is_running()).unwrap_or(false) {
            self.shared_state.as_ref()
        } else {
            None
        }
    }
}

#[allow(unused_variables)]
#[pymethods]
impl PyCmd {
    // =========================================================================
    // File I/O
    // =========================================================================

    /// Load a molecular structure file
    ///
    /// Args:
    ///     filename: Path to the file
    ///     object: Optional object name (default: derived from filename)
    ///     state: State to load into (default: append new state)
    ///     format: File format (default: auto-detect)
    ///     quiet: Suppress output (default: True)
    ///
    /// Returns:
    ///     Name of the loaded object
    #[pyo3(signature = (filename, object=None, state=0, format=None, quiet=true))]
    fn load(
        &mut self,
        filename: &str,
        object: Option<&str>,
        state: i32,
        format: Option<&str>,
        quiet: bool,
    ) -> PyResult<String> {
        let path = Path::new(filename);
        
        // Determine object name
        let obj_name = object
            .map(|s| s.to_string())
            .unwrap_or_else(|| {
                path.file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("obj")
                    .to_string()
            });

        // Read the file
        let mol = if let Some(fmt) = format {
            let file_format = match fmt.to_lowercase().as_str() {
                "pdb" => pymol_io::FileFormat::Pdb,
                "sdf" | "mol" => pymol_io::FileFormat::Sdf,
                "mol2" => pymol_io::FileFormat::Mol2,
                "cif" | "mmcif" => pymol_io::FileFormat::Cif,
                "xyz" => pymol_io::FileFormat::Xyz,
                _ => pymol_io::FileFormat::Unknown,
            };
            pymol_io::read_file_format(path, file_format).map_py_err()?
        } else {
            pymol_io::read_file(path).map_py_err()?
        };

        // Set the object name
        let mut mol = mol;
        mol.name = obj_name.clone();

        // Add to viewer
        self.viewer.add_molecule(mol);

        if !quiet {
            log::info!("Loaded {} atoms from {}", 
                self.viewer.objects().get_molecule(&obj_name)
                    .map(|m| m.molecule().atom_count())
                    .unwrap_or(0),
                filename);
        }

        Ok(obj_name)
    }

    /// Save molecular structure to a file
    ///
    /// Args:
    ///     filename: Output file path
    ///     selection: Atoms to save (default: all)
    ///     state: State to save (default: current)
    ///     format: File format (default: auto-detect from extension)
    #[pyo3(signature = (filename, selection="all", state=0, format=None))]
    fn save(
        &self,
        filename: &str,
        selection: &str,
        state: i32,
        format: Option<&str>,
    ) -> PyResult<()> {
        let path = Path::new(filename);
        
        // For now, save all atoms of first matching object
        // TODO: proper selection handling for save
        for name in self.viewer.objects().names() {
            if let Some(mol_obj) = self.viewer.objects().get_molecule(name) {
                if let Some(fmt) = format {
                    let file_format = match fmt.to_lowercase().as_str() {
                        "pdb" => pymol_io::FileFormat::Pdb,
                        "sdf" | "mol" => pymol_io::FileFormat::Sdf,
                        "mol2" => pymol_io::FileFormat::Mol2,
                        "cif" | "mmcif" => pymol_io::FileFormat::Cif,
                        "xyz" => pymol_io::FileFormat::Xyz,
                        _ => pymol_io::FileFormat::Unknown,
                    };
                    pymol_io::write_file_format(path, mol_obj.molecule(), file_format).map_py_err()?;
                } else {
                    pymol_io::write_file(path, mol_obj.molecule()).map_py_err()?;
                }
                return Ok(());
            }
        }

        Err(pyo3::exceptions::PyValueError::new_err("No objects to save"))
    }

    /// Fetch a structure from RCSB PDB
    ///
    /// Args:
    ///     code: 4-letter PDB ID
    ///     name: Optional object name (default: use PDB code)
    ///     type_: Download format ("pdb" or "cif", default: "cif")
    ///
    /// Returns:
    ///     Name of the fetched object
    #[pyo3(signature = (code, name=None, type_="cif"))]
    fn fetch(&mut self, code: &str, name: Option<&str>, type_: &str) -> PyResult<String> {
        let fetch_format = match type_.to_lowercase().as_str() {
            "pdb" => pymol_io::FetchFormat::Pdb,
            "cif" | "mmcif" => pymol_io::FetchFormat::Cif,
            _ => {
                return Err(pyo3::exceptions::PyValueError::new_err(format!(
                    "Unknown fetch format: '{}'. Use 'pdb' or 'cif'.",
                    type_
                )))
            }
        };

        let mol = pymol_io::fetch(code, fetch_format).map_py_err()?;

        let obj_name = name.unwrap_or(code).to_string();
        let mut mol = mol;
        mol.name = obj_name.clone();

        self.viewer.add_molecule(mol);

        Ok(obj_name)
    }

    // =========================================================================
    // Display Commands
    // =========================================================================

    /// Show a representation
    ///
    /// Args:
    ///     representation: Type of representation ("lines", "sticks", "spheres", 
    ///                     "surface", "cartoon", "ribbon", "mesh", "dots")
    ///     selection: Atoms to show (default: all)
    #[pyo3(signature = (representation="lines", selection="all"))]
    fn show(&mut self, representation: &str, selection: &str) -> PyResult<()> {
        let rep = parse_representation(representation)?;
        
        // Collect names first
        let names: Vec<String> = self.viewer.objects().names().map(|s| s.to_string()).collect();
        
        // Apply to all molecules
        for name in names {
            if let Some(mol_obj) = self.viewer.objects_mut().get_molecule_mut(&name) {
                mol_obj.show(rep);
            }
        }
        
        self.viewer.request_redraw();
        Ok(())
    }

    /// Hide a representation
    ///
    /// Args:
    ///     representation: Type of representation (or "everything" for all)
    ///     selection: Atoms to hide (default: all)
    #[pyo3(signature = (representation="everything", selection="all"))]
    fn hide(&mut self, representation: &str, selection: &str) -> PyResult<()> {
        // Collect names first
        let names: Vec<String> = self.viewer.objects().names().map(|s| s.to_string()).collect();
        
        if representation == "everything" {
            for name in names {
                if let Some(mol_obj) = self.viewer.objects_mut().get_molecule_mut(&name) {
                    mol_obj.hide_all();
                }
            }
        } else {
            let rep = parse_representation(representation)?;
            for name in names {
                if let Some(mol_obj) = self.viewer.objects_mut().get_molecule_mut(&name) {
                    mol_obj.hide(rep);
                }
            }
        }
        
        self.viewer.request_redraw();
        Ok(())
    }

    /// Color atoms or representations
    ///
    /// Args:
    ///     color: Color name (e.g., "red", "green", "blue") or hex code
    ///     selection: Atoms to color (default: all)
    #[pyo3(signature = (color, selection="all"))]
    fn color(&mut self, color: &str, selection: &str) -> PyResult<()> {
        // Look up color by name
        let color_idx = self.viewer.named_colors()
            .get_by_name(color)
            .map(|(idx, _)| idx as i32)
            .unwrap_or(-1); // -1 = by element

        // Collect names first
        let names: Vec<String> = self.viewer.objects().names().map(|s| s.to_string()).collect();

        // Apply color to matching atoms
        for name in names {
            if let Some(mol_obj) = self.viewer.objects_mut().get_molecule_mut(&name) {
                let mol = mol_obj.molecule_mut();
                
                // Parse selection and apply color
                if let Ok(result) = pymol_select::select(mol, selection) {
                    for idx in result.indices() {
                        if let Some(atom) = mol.get_atom_mut(idx) {
                            atom.color = color_idx;
                        }
                    }
                }
                
                // Mark molecule as needing re-rendering
                mol_obj.invalidate(pymol_scene::DirtyFlags::COLOR);
            }
        }

        self.viewer.request_redraw();
        Ok(())
    }

    // =========================================================================
    // Selection Commands
    // =========================================================================

    /// Create a named selection
    ///
    /// Args:
    ///     name: Name for the selection
    ///     selection: Selection expression
    ///     enable: Whether to enable the selection (default: auto)
    ///     quiet: Suppress output (default: True)
    ///
    /// Returns:
    ///     Number of atoms selected
    #[pyo3(signature = (name, selection, enable=-1, quiet=true))]
    fn select(&mut self, name: &str, selection: &str, enable: i32, quiet: bool) -> PyResult<i32> {
        let mut total = 0;

        // Store the selection
        if let ViewerMode::Headless(ref mut h) = self.viewer {
            h.define_selection(name, selection);
        }

        // Count matching atoms
        for obj_name in self.viewer.objects().names() {
            if let Some(mol_obj) = self.viewer.objects().get_molecule(obj_name) {
                if let Ok(result) = pymol_select::select(mol_obj.molecule(), selection) {
                    total += result.count() as i32;
                }
            }
        }

        if !quiet {
            log::info!("Selection '{}' created with {} atoms", name, total);
        }

        Ok(total)
    }

    /// Remove a named selection
    #[pyo3(signature = (name="sele"))]
    fn deselect(&mut self, name: &str) -> PyResult<()> {
        if let ViewerMode::Headless(ref mut h) = self.viewer {
            h.selections.remove(name);
        }
        Ok(())
    }

    /// Count atoms in a selection
    ///
    /// Args:
    ///     selection: Selection expression (default: all)
    ///
    /// Returns:
    ///     Number of atoms
    #[pyo3(signature = (selection="all"))]
    fn count_atoms(&self, selection: &str) -> PyResult<i32> {
        let mut total = 0;

        for name in self.viewer.objects().names() {
            if let Some(mol_obj) = self.viewer.objects().get_molecule(name) {
                if let Ok(result) = pymol_select::select(mol_obj.molecule(), selection) {
                    total += result.count() as i32;
                }
            }
        }

        Ok(total)
    }

    // =========================================================================
    // Viewing Commands
    // =========================================================================

    /// Zoom the camera to fit objects
    ///
    /// Args:
    ///     selection: Objects/atoms to zoom on (default: all)
    ///     buffer: Extra space around the selection (default: 0.0)
    ///     state: State to consider (default: all states)
    ///     complete: Whether to complete the zoom instantly (default: False)
    #[pyo3(signature = (selection="all", buffer=0.0, state=0, complete=false))]
    fn zoom(&mut self, selection: &str, buffer: f32, state: i32, complete: bool) -> PyResult<()> {
        self.viewer.zoom_all();
        Ok(())
    }

    /// Center the camera on objects
    ///
    /// Args:
    ///     selection: Objects/atoms to center on (default: all)
    ///     state: State to consider (default: all states)
    #[pyo3(signature = (selection="all", state=0))]
    fn center(&mut self, selection: &str, state: i32) -> PyResult<()> {
        self.viewer.center_all();
        Ok(())
    }

    /// Orient the camera to show objects optimally
    #[pyo3(signature = (selection="all", state=0))]
    fn orient(&mut self, selection: &str, state: i32) -> PyResult<()> {
        self.viewer.zoom_all();
        Ok(())
    }

    /// Reset the view to default
    fn reset(&mut self) -> PyResult<()> {
        self.viewer.reset_view();
        Ok(())
    }

    // =========================================================================
    // Object Commands
    // =========================================================================

    /// Delete an object or selection
    ///
    /// Args:
    ///     name: Object or selection name to delete
    #[pyo3(signature = (name))]
    fn delete(&mut self, name: &str) -> PyResult<()> {
        if name == "all" {
            match &mut self.viewer {
                ViewerMode::Headless(h) => h.delete_all(),
                ViewerMode::Gui(v) => v.objects_mut().clear(),
            }
        } else {
            self.viewer.objects_mut().remove(name);
        }
        Ok(())
    }

    /// List all object names
    fn get_names(&self) -> Vec<String> {
        self.viewer.objects().names().map(|s| s.to_string()).collect()
    }

    /// Check if an object exists
    fn object_exists(&self, name: &str) -> bool {
        self.viewer.objects().get_molecule(name).is_some()
    }

    /// Get object atom count
    fn get_object_atom_count(&self, name: &str) -> Option<usize> {
        self.viewer.objects().get_molecule(name).map(|mol_obj| {
            mol_obj.molecule().atom_count()
        })
    }

    /// Get object bond count
    fn get_object_bond_count(&self, name: &str) -> Option<usize> {
        self.viewer.objects().get_molecule(name).map(|mol_obj| {
            mol_obj.molecule().bond_count()
        })
    }

    // =========================================================================
    // Settings Commands
    // =========================================================================

    /// Set a PyMOL setting
    ///
    /// Args:
    ///     name: Setting name
    ///     value: New value
    ///     selection: Optional selection (for per-atom settings)
    ///     quiet: Suppress output
    #[pyo3(signature = (name, value, selection=None, quiet=true))]
    fn set(&mut self, name: &str, value: &Bound<'_, PyAny>, selection: Option<&str>, quiet: bool) -> PyResult<()> {
        set_setting_py(self.viewer.settings_mut(), name, value)
    }

    /// Get a PyMOL setting value
    ///
    /// Args:
    ///     name: Setting name
    ///     selection: Optional selection (for per-atom settings)
    ///
    /// Returns:
    ///     The setting value
    #[pyo3(signature = (name, selection=None))]
    fn get<'py>(&self, py: Python<'py>, name: &str, selection: Option<&str>) -> PyResult<PyObject> {
        get_setting_py(py, self.viewer.settings(), name)
    }

    // =========================================================================
    // State Commands
    // =========================================================================

    /// Get the current state number (1-based)
    fn get_state(&self) -> i32 {
        match &self.viewer {
            ViewerMode::Headless(h) => (h.current_state + 1) as i32,
            ViewerMode::Gui(_) => 1, // TODO: implement for GUI
        }
    }

    /// Set the current state (1-based)
    fn set_state(&mut self, state: i32) -> PyResult<()> {
        if state < 1 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "State must be >= 1",
            ));
        }
        match &mut self.viewer {
            ViewerMode::Headless(h) => h.current_state = (state - 1) as usize,
            ViewerMode::Gui(_) => {} // TODO: implement for GUI
        }
        Ok(())
    }

    // =========================================================================
    // Command Execution
    // =========================================================================

    /// Execute a PyMOL command string
    ///
    /// Args:
    ///     command: Command string (e.g., "load file.pdb; show cartoon")
    #[pyo3(name = "do")]
    fn do_(&mut self, command: &str) -> PyResult<()> {
        // For now, parse and execute simple commands
        // TODO: full command execution via CommandExecutor
        for cmd in command.split(';') {
            let cmd = cmd.trim();
            if cmd.is_empty() {
                continue;
            }

            let parts: Vec<&str> = cmd.splitn(2, ' ').collect();
            let name = parts[0];
            let args = parts.get(1).unwrap_or(&"");

            match name {
                "load" => {
                    self.load(args, None, 0, None, true)?;
                }
                "fetch" => {
                    self.fetch(args, None, "cif")?;
                }
                "show" => {
                    let parts: Vec<&str> = args.splitn(2, ',').collect();
                    let rep = parts.get(0).unwrap_or(&"lines").trim();
                    let sel = parts.get(1).unwrap_or(&"all").trim();
                    self.show(rep, sel)?;
                }
                "hide" => {
                    let parts: Vec<&str> = args.splitn(2, ',').collect();
                    let rep = parts.get(0).unwrap_or(&"everything").trim();
                    let sel = parts.get(1).unwrap_or(&"all").trim();
                    self.hide(rep, sel)?;
                }
                "color" => {
                    let parts: Vec<&str> = args.splitn(2, ',').collect();
                    let color = parts.get(0).unwrap_or(&"white").trim();
                    let sel = parts.get(1).unwrap_or(&"all").trim();
                    self.color(color, sel)?;
                }
                "zoom" => {
                    self.zoom(args, 0.0, 0, false)?;
                }
                "center" => {
                    self.center(args, 0)?;
                }
                "reset" => {
                    self.reset()?;
                }
                "delete" => {
                    self.delete(args)?;
                }
                "png" => {
                    // Parse: png filename [, width [, height]]
                    let parts: Vec<&str> = args.split(',').map(|s| s.trim()).collect();
                    let filename = parts.get(0).unwrap_or(&"output.png");
                    let width = parts.get(1).and_then(|s| s.parse().ok());
                    let height = parts.get(2).and_then(|s| s.parse().ok());
                    self.png(filename, width, height, 0, true)?;
                }
                _ => {
                    return Err(pyo3::exceptions::PyValueError::new_err(format!(
                        "Unknown command: {}",
                        name
                    )));
                }
            }
        }

        Ok(())
    }

    // =========================================================================
    // Utility
    // =========================================================================

    // =========================================================================
    // Image Output
    // =========================================================================

    /// Save current view as PNG image
    ///
    /// Args:
    ///     filename: Output file path (will add .png extension if missing)
    ///     width: Image width in pixels (default: current window size or 1024)
    ///     height: Image height in pixels (default: current window size or 768)
    ///     ray: Whether to ray-trace (not yet implemented, ignored)
    ///     quiet: Suppress output (default: True)
    ///
    /// Example:
    ///     cmd.png("output.png")
    ///     cmd.png("hires.png", 1920, 1080)
    ///     cmd.png("square.png", width=800, height=800)
    #[pyo3(signature = (filename, width=None, height=None, ray=0, quiet=true))]
    fn png(
        &mut self,
        filename: &str,
        width: Option<u32>,
        height: Option<u32>,
        ray: i32,
        quiet: bool,
    ) -> PyResult<()> {
        let path = Path::new(filename);
        
        // Ensure the path ends with .png
        let path = if path.extension().map(|e| e.to_ascii_lowercase()) != Some("png".into()) {
            path.with_extension("png")
        } else {
            path.to_path_buf()
        };

        match &mut self.viewer {
            ViewerMode::Gui(viewer) => {
                viewer.capture_png(&path, width, height).map_err(|e| {
                    pyo3::exceptions::PyRuntimeError::new_err(format!(
                        "Failed to capture PNG: {}",
                        e
                    ))
                })?;

                if !quiet {
                    let (w, h) = (width.unwrap_or(1024), height.unwrap_or(768));
                    log::info!("Saved PNG image to {} ({}x{})", path.display(), w, h);
                }

                Ok(())
            }
            ViewerMode::Headless(_) => {
                Err(pyo3::exceptions::PyRuntimeError::new_err(
                    "PNG capture requires GUI mode. Initialize PyMOL with headless=False or use ray() for offline rendering (not yet implemented)."
                ))
            }
        }
    }

    /// Refresh the display
    fn refresh(&mut self) {
        self.viewer.request_redraw();
    }

    /// Launch the interactive viewer window
    ///
    /// This opens a GUI window displaying the current molecular scene.
    /// By default, the viewer runs in the background (non-blocking) so you
    /// can continue to run commands that will update the display.
    ///
    /// Args:
    ///     width: Window width in pixels (default: 1024)
    ///     height: Window height in pixels (default: 768)
    ///     blocking: If True, block until window is closed (default: False)
    ///
    /// Example:
    ///     cmd.load("protein.pdb")
    ///     cmd.show("cartoon")
    ///     cmd.show_viewer()  # Opens window, returns immediately
    ///     cmd.color("red", "chain A")  # Updates the display
    ///     cmd.close_viewer()  # Close when done
    #[pyo3(signature = (width=1024, height=768, blocking=false))]
    fn show_viewer(&mut self, width: u32, height: u32, blocking: bool) -> PyResult<()> {
        // Check if viewer is already running
        if self.viewer_handle.as_ref().map(|h| h.is_running()).unwrap_or(false) {
            return Err(pyo3::exceptions::PyRuntimeError::new_err(
                "Viewer is already running. Use close_viewer() first."
            ));
        }

        // Create shared state
        let shared_state = Arc::new(RwLock::new(SharedViewerState::new()));
        
        // Transfer all molecules from the current viewer to the shared state
        {
            let mut state = shared_state.write();
            for name in self.viewer.objects().names() {
                if let Some(mol_obj) = self.viewer.objects().get_molecule(name) {
                    let mol = mol_obj.molecule().clone();
                    state.add_molecule(mol);
                }
            }
            // Copy settings
            state.settings = self.viewer.settings().clone();
        }

        // Create running flag
        let running = Arc::new(AtomicBool::new(true));

        if blocking {
            // Run in blocking mode on this thread
            let mut gui_viewer = GuiViewer::new();
            
            // Transfer from shared state to gui_viewer
            {
                let state = shared_state.read();
                for name in state.objects.names() {
                    if let Some(mol_obj) = state.objects.get_molecule(name) {
                        let mol = mol_obj.molecule().clone();
                        gui_viewer.add_molecule(mol);
                    }
                }
                *gui_viewer.settings_mut() = state.settings.clone();
            }
            
            gui_viewer.zoom_all();
            
            pymol_scene::run(gui_viewer).map_err(|e| {
                pyo3::exceptions::PyRuntimeError::new_err(format!(
                    "Failed to run viewer: {}",
                    e
                ))
            })?;
        } else {
            // Run in non-blocking mode in a background thread
            let state_clone = Arc::clone(&shared_state);
            let running_clone = Arc::clone(&running);
            let _width = width;
            let _height = height;
            
            let thread = std::thread::spawn(move || {
                let mut gui_viewer = GuiViewer::new();
                
                // Transfer from shared state to gui_viewer
                {
                    let state = state_clone.read();
                    for name in state.objects.names() {
                        if let Some(mol_obj) = state.objects.get_molecule(name) {
                            let mol = mol_obj.molecule().clone();
                            gui_viewer.add_molecule(mol);
                        }
                    }
                    *gui_viewer.settings_mut() = state.settings.clone();
                }
                
                gui_viewer.zoom_all();
                
                // Run the viewer - this blocks until the window is closed
                if let Err(e) = pymol_scene::run(gui_viewer) {
                    log::error!("Viewer error: {}", e);
                }
                
                // Mark as no longer running
                running_clone.store(false, Ordering::SeqCst);
            });
            
            // Store the handle
            self.shared_state = Some(shared_state);
            self.viewer_handle = Some(ViewerHandle::new(running, self.shared_state.clone().unwrap(), Some(thread)));
        }

        Ok(())
    }

    /// Check if the viewer window is open
    ///
    /// Returns:
    ///     True if the viewer is running, False otherwise
    fn viewer_is_open(&self) -> bool {
        self.viewer_handle.as_ref().map(|h| h.is_running()).unwrap_or(false)
    }

    /// Close the viewer window
    ///
    /// This requests the viewer to close. The window may not close immediately.
    fn close_viewer(&mut self) {
        if let Some(handle) = self.viewer_handle.take() {
            handle.request_close();
            // Note: We don't wait for the thread to finish to avoid blocking
        }
        self.shared_state = None;
    }

    /// Wait for the viewer window to close (blocking)
    ///
    /// This blocks until the user closes the viewer window.
    fn wait_viewer(&mut self) {
        if let Some(handle) = self.viewer_handle.take() {
            handle.wait();
        }
        self.shared_state = None;
    }

    /// Print version information
    fn version(&self) -> String {
        format!("PyMOL-RS {}", env!("CARGO_PKG_VERSION"))
    }

    fn __repr__(&self) -> String {
        format!(
            "Cmd(objects={}, mode={})",
            self.viewer.objects().names().count(),
            if self.viewer.is_headless() { "headless" } else { "gui" }
        )
    }
}

/// Parse a representation name to RepMask value
fn parse_representation(name: &str) -> PyResult<u32> {
    match name.to_lowercase().as_str() {
        "lines" | "line" => Ok(RepMask::LINES),
        "sticks" | "stick" => Ok(RepMask::STICKS),
        "spheres" | "sphere" => Ok(RepMask::SPHERES),
        "surface" | "surf" => Ok(RepMask::SURFACE),
        "cartoon" | "cart" => Ok(RepMask::CARTOON),
        "ribbon" | "rib" => Ok(RepMask::RIBBON),
        "mesh" => Ok(RepMask::MESH),
        "dots" | "dot" => Ok(RepMask::DOTS),
        "nonbonded" | "nb" => Ok(RepMask::NONBONDED),
        "labels" | "label" => Ok(RepMask::LABELS),
        _ => Err(pyo3::exceptions::PyValueError::new_err(format!(
            "Unknown representation: '{}'",
            name
        ))),
    }
}
