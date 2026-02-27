//! PyCmd - Main command interface for Python
//!
//! This module provides the `cmd` object that users interact with,
//! following PyMOL's familiar API pattern.
//!
//! All operations are performed via IPC to the pymol-rs server.

use std::path::Path;
use std::sync::{Arc, Mutex};

use pyo3::prelude::*;
use pyo3::types::PyAny;

use crate::connection::{establish_connection, Connection};
use crate::ipc::{CallbackListener, ExtendedCommands, IpcClient, IpcResponse};

/// Main command interface
///
/// Provides the familiar PyMOL `cmd.function()` API for molecular visualization.
/// All operations are performed via IPC to the pymol-rs server.
///
/// On module import, this class automatically connects to an existing IPC server
/// or spawns a new headless server.
///
/// Example:
///     from pymol_rs import cmd
///     cmd.load("protein.pdb")
///     cmd.show("cartoon")
///     cmd.color("green", "chain A")
///     cmd.show_gui()  # Show the window
///     cmd.zoom()
#[pyclass(name = "Cmd")]
pub struct PyCmd {
    /// Connection to the pymol-rs server (handles lifecycle management)
    connection: Option<Connection>,
    /// Background listener for callback requests (lazily started)
    listener: Arc<Mutex<Option<CallbackListener>>>,
    /// Extended command callbacks (shared with listener)
    extended_commands: Arc<Mutex<ExtendedCommands>>,
    /// Default silent mode for command execution
    silent: bool,
}

impl PyCmd {
    /// Create a new PyCmd instance by establishing an IPC connection
    ///
    /// This will connect to an existing server or spawn a new headless one.
    pub fn new() -> PyResult<Self> {
        log::info!("PyCmd::new() - establishing connection");
        
        let connection = establish_connection().map_err(|e| {
            pyo3::exceptions::PyRuntimeError::new_err(format!(
                "Failed to connect to pymol-rs server: {}",
                e
            ))
        })?;

        log::info!("PyCmd::new() - connection established");

        let extended_commands = Arc::new(Mutex::new(ExtendedCommands::new()));
        let listener = Arc::new(Mutex::new(None));

        // Note: Callback listener is started lazily when first needed,
        // to avoid issues with the listener thread competing with initial setup
        log::info!("PyCmd::new() - complete");

        Ok(PyCmd {
            connection: Some(connection),
            listener,
            extended_commands,
            silent: false,
        })
    }

    /// Get the client Arc for sharing with other components
    fn client_arc(&self) -> Option<Arc<Mutex<Option<IpcClient>>>> {
        self.connection.as_ref().map(|c| c.client_arc())
    }

    /// Execute a function with the locked IPC client
    ///
    /// This is the primary way to interact with the IPC client, ensuring
    /// proper error handling and lock management.
    fn with_client<F, R>(&self, f: F) -> PyResult<R>
    where
        F: FnOnce(&mut IpcClient) -> std::io::Result<R>,
    {
        let client_arc = self.client_arc().ok_or_else(|| {
            pyo3::exceptions::PyRuntimeError::new_err("Not connected to server")
        })?;

        let mut guard = client_arc.lock().map_err(|e| {
            pyo3::exceptions::PyRuntimeError::new_err(format!("Failed to lock client: {}", e))
        })?;

        let client = guard.as_mut().ok_or_else(|| {
            pyo3::exceptions::PyRuntimeError::new_err("IPC client not available")
        })?;

        f(client).map_err(|e| {
            pyo3::exceptions::PyRuntimeError::new_err(format!("IPC error: {}", e))
        })
    }
    
    /// Ensure the callback listener is running
    ///
    /// The listener runs in a background thread and handles callback requests
    /// from the GUI (e.g., for runpy and extended commands).
    fn ensure_listener_running(&self) {
        let Some(client_arc) = self.client_arc() else {
            return;
        };
        
        let mut listener_guard = self.listener.lock().unwrap();
        if listener_guard.is_none() {
            log::info!("Starting callback listener");
            *listener_guard = Some(CallbackListener::start(
                client_arc,
                self.extended_commands.clone(),
            ));
        }
    }

    /// Create from an existing connection
    pub fn from_connection(connection: Connection) -> PyResult<Self> {
        let extended_commands = Arc::new(Mutex::new(ExtendedCommands::new()));
        let listener = Arc::new(Mutex::new(None));

        Ok(PyCmd {
            connection: Some(connection),
            listener,
            extended_commands,
            silent: false,
        })
    }

    /// Execute a command via IPC and return the result
    ///
    /// Uses the default silent mode set by `set_silent()`.
    fn execute_cmd(&self, command: &str) -> PyResult<()> {
        self.execute_cmd_opts(command, self.silent)
    }

    /// Execute a command via IPC with options
    fn execute_cmd_opts(&self, command: &str, silent: bool) -> PyResult<()> {
        let response = self.with_client(|client| {
            if silent {
                client.execute_silent(command)
            } else {
                client.execute(command)
            }
        })?;

        match response {
            IpcResponse::Ok { .. } => Ok(()),
            IpcResponse::Error { message, .. } => {
                Err(pyo3::exceptions::PyRuntimeError::new_err(message))
            }
            _ => Ok(()), // Other responses are OK
        }
    }

    /// Register the runpy command with the server
    pub fn register_runpy(&self) -> PyResult<()> {
        // Start the listener since runpy is a callback command
        self.ensure_listener_running();
        
        self.with_client(|client| {
            client.register_command(
                "runpy",
                Some(r#"
DESCRIPTION

    "runpy" executes a Python script file (.py).

USAGE

    runpy filename [, namespace ]

ARGUMENTS

    filename = string: path to Python script file
    namespace = string: execution namespace (default: global)

EXAMPLES

    runpy setup.py
    runpy analysis.py, local
"#),
            )
        })
    }
}

impl Drop for PyCmd {
    fn drop(&mut self) {
        // Stop the listener first (before Connection cleanup sends quit)
        if let Ok(mut listener_guard) = self.listener.lock() {
            if let Some(ref mut listener) = *listener_guard {
                listener.stop();
            }
        }
        // Connection's Drop handles server cleanup automatically
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
        &self,
        filename: &str,
        object: Option<&str>,
        state: i32,
        format: Option<&str>,
        quiet: bool,
    ) -> PyResult<String> {
        let path = Path::new(filename);

        // Get absolute path
        let abs_path = if path.is_absolute() {
            path.to_path_buf()
        } else {
            std::env::current_dir()
                .map(|cwd| cwd.join(path))
                .unwrap_or_else(|_| path.to_path_buf())
        };

        // Determine object name
        let obj_name = object.map(|s| s.to_string()).unwrap_or_else(|| {
            path.file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("obj")
                .to_string()
        });

        let cmd = format!("load {}, {}", abs_path.display(), obj_name);
        self.execute_cmd_opts(&cmd, quiet)?;

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
        let cmd = format!("save {}, {}", filename, selection);
        self.execute_cmd(&cmd)
    }

    /// Fetch a structure from RCSB PDB
    ///
    /// Args:
    ///     code: 4-letter PDB ID
    ///     name: Optional object name (default: use PDB code)
    ///     type_: Download format ("pdb" or "cif", default: "cif")
    ///     wait: If True, block until the structure is loaded (default: False)
    ///     timeout: Maximum time to wait in seconds when wait=True (default: 30.0)
    ///
    /// Returns:
    ///     Name of the fetched object
    ///
    /// Example:
    ///     # Non-blocking (returns immediately, fetch runs in background)
    ///     cmd.fetch("1ubq")
    ///
    ///     # Blocking (waits until structure is loaded)
    ///     cmd.fetch("1ubq", sync=True)
    ///     cmd.show("cartoon")  # Safe to use immediately after
    #[pyo3(signature = (code, name=None, type_="cif", sync=false, timeout=30.0))]
    fn fetch(&self, code: &str, name: Option<&str>, type_: &str, sync: bool, timeout: f64) -> PyResult<String> {
        let obj_name = name.unwrap_or(code).to_string();
        let cmd = format!("fetch {}", code);
        self.execute_cmd(&cmd)?;

        if sync {
            // Poll get_names() until the object appears or timeout
            let start = std::time::Instant::now();
            let timeout_duration = std::time::Duration::from_secs_f64(timeout);
            let poll_interval = std::time::Duration::from_millis(100);

            loop {
                // Check if object exists
                let names = self.with_client(|client| client.get_names())?;
                if names.iter().any(|n| n == &obj_name) {
                    break;
                }

                // Check timeout
                if start.elapsed() >= timeout_duration {
                    return Err(pyo3::exceptions::PyTimeoutError::new_err(format!(
                        "Timeout waiting for fetch of '{}' after {:.1}s",
                        code, timeout
                    )));
                }

                // Sleep briefly before next poll
                std::thread::sleep(poll_interval);
            }
        }

        Ok(obj_name)
    }

    /// Execute a Python script file
    ///
    /// Args:
    ///     filename: Path to .py file
    ///     namespace: Execution namespace:
    ///         - "global": Execute in pymol_rs module namespace (default)
    ///         - "local": Execute in a fresh namespace with `cmd` pre-imported
    ///         - "module": Import as a Python module
    #[pyo3(signature = (filename, namespace="global"))]
    fn run(&self, py: Python<'_>, filename: &str, namespace: &str) -> PyResult<()> {
        crate::scripting::run_python_script(py, filename, namespace)
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
    fn show(&self, representation: &str, selection: &str) -> PyResult<()> {
        let cmd = format!("show {}, {}", representation, selection);
        self.execute_cmd(&cmd)
    }

    /// Hide a representation
    ///
    /// Args:
    ///     representation: Type of representation (or "everything" for all)
    ///     selection: Atoms to hide (default: all)
    #[pyo3(signature = (representation="everything", selection="all"))]
    fn hide(&self, representation: &str, selection: &str) -> PyResult<()> {
        let cmd = format!("hide {}, {}", representation, selection);
        self.execute_cmd(&cmd)
    }

    /// Show a representation as the only one (hide all others)
    ///
    /// This is equivalent to hiding everything and then showing the
    /// specified representation for the selection.
    ///
    /// Args:
    ///     representation: Type of representation ("lines", "sticks", "spheres",
    ///                     "surface", "cartoon", "ribbon", "mesh", "dots")
    ///     selection: Atoms to show (default: all)
    ///
    /// Example:
    ///     cmd.show_as("cartoon", "chain A")  # Show only cartoon for chain A
    #[pyo3(signature = (representation, selection="all"))]
    fn show_as(&self, representation: &str, selection: &str) -> PyResult<()> {
        let cmd = format!("as {}, {}", representation, selection);
        self.execute_cmd(&cmd)
    }

    /// Color atoms or representations
    ///
    /// Args:
    ///     color: Color name (e.g., "red", "green", "blue") or hex code
    ///     selection: Atoms to color (default: all)
    #[pyo3(signature = (color, selection="all"))]
    fn color(&self, color: &str, selection: &str) -> PyResult<()> {
        let cmd = format!("color {}, {}", color, selection);
        self.execute_cmd(&cmd)
    }

    /// Set background color
    ///
    /// Args:
    ///     color: Color name (e.g., "white", "black", "gray")
    #[pyo3(signature = (color))]
    fn bg_color(&self, color: &str) -> PyResult<()> {
        let cmd = format!("bg_color {}", color);
        self.execute_cmd(&cmd)
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
    ///     Number of atoms selected (0 when using IPC)
    #[pyo3(signature = (name, selection, enable=-1, quiet=true))]
    fn select(&self, name: &str, selection: &str, enable: i32, quiet: bool) -> PyResult<i32> {
        let cmd = format!("select {}, {}", name, selection);
        self.execute_cmd_opts(&cmd, quiet)?;
        // Can't return accurate count via IPC currently
        Ok(0)
    }

    /// Remove a named selection
    #[pyo3(signature = (name="sele"))]
    fn deselect(&self, name: &str) -> PyResult<()> {
        let cmd = format!("deselect {}", name);
        self.execute_cmd(&cmd)
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
        let count = self.with_client(|client| client.count_atoms(selection))?;
        Ok(count as i32)
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
    fn zoom(&self, selection: &str, buffer: f32, state: i32, complete: bool) -> PyResult<()> {
        let cmd = format!("zoom {}", selection);
        self.execute_cmd(&cmd)
    }

    /// Center the camera on objects
    ///
    /// Args:
    ///     selection: Objects/atoms to center on (default: all)
    ///     state: State to consider (default: all states)
    #[pyo3(signature = (selection="all", state=0))]
    fn center(&self, selection: &str, state: i32) -> PyResult<()> {
        let cmd = format!("center {}", selection);
        self.execute_cmd(&cmd)
    }

    /// Orient the camera to show objects optimally
    #[pyo3(signature = (selection="all", state=0))]
    fn orient(&self, selection: &str, state: i32) -> PyResult<()> {
        let cmd = format!("orient {}", selection);
        self.execute_cmd(&cmd)
    }

    /// Reset the view to default
    fn reset(&self) -> PyResult<()> {
        self.execute_cmd("reset")
    }

    /// Translate the camera along an axis
    ///
    /// Args:
    ///     axis: Axis to move along ("x", "y", or "z")
    ///     distance: Distance to move in Angstroms
    ///
    /// Notes:
    ///     Positive x moves right, positive y moves up, positive z moves
    ///     toward the viewer.
    ///
    /// See Also:
    ///     turn, zoom, center, clip
    #[pyo3(name = "move", signature = (axis, distance))]
    fn move_(&self, axis: &str, distance: f64) -> PyResult<()> {
        let cmd = format!("move {}, {}", axis, distance);
        self.execute_cmd(&cmd)
    }

    /// Rotate the camera about an axis
    ///
    /// Args:
    ///     axis: Axis to rotate about ("x", "y", or "z")
    ///     angle: Degrees of rotation
    ///
    /// Notes:
    ///     Rotations follow the right-hand rule. For example, a positive
    ///     rotation about the y-axis will rotate the view to the right.
    ///
    /// See Also:
    ///     move, rotate, zoom, center
    #[pyo3(signature = (axis, angle))]
    fn turn(&self, axis: &str, angle: f64) -> PyResult<()> {
        let cmd = format!("turn {}, {}", axis, angle);
        self.execute_cmd(&cmd)
    }

    /// Adjust near and far clipping planes
    ///
    /// Args:
    ///     mode: Clipping mode ("near", "far", "move", or "slab")
    ///     distance: Clipping distance adjustment
    ///     selection: Selection for reference (default: all)
    ///     state: State to consider (default: 0)
    ///
    /// Modes:
    ///     near: adjust near clipping plane
    ///     far: adjust far clipping plane
    ///     move: move both planes
    ///     slab: set slab thickness
    ///
    /// Examples:
    ///     cmd.clip("near", -5)
    ///     cmd.clip("far", 10)
    ///     cmd.clip("slab", 20)
    ///
    /// See Also:
    ///     zoom, move, turn
    #[pyo3(signature = (mode, distance, selection="all", state=0))]
    fn clip(&self, mode: &str, distance: f64, selection: &str, state: i32) -> PyResult<()> {
        let cmd = format!("clip {}, {}", mode, distance);
        self.execute_cmd(&cmd)
    }

    /// Change viewport size
    ///
    /// Args:
    ///     width: Width in pixels (optional)
    ///     height: Height in pixels (optional)
    ///
    /// Notes:
    ///     If only width is specified, height will be scaled to maintain
    ///     the current aspect ratio.
    ///     If no arguments are given, the current viewport size is displayed.
    ///
    /// Examples:
    ///     cmd.viewport()            # Show current size
    ///     cmd.viewport(800, 600)    # Set to 800x600
    ///     cmd.viewport(1920, 1080)  # Set to 1920x1080
    ///
    /// See Also:
    ///     full_screen, png
    #[pyo3(signature = (width=None, height=None))]
    fn viewport(&self, width: Option<u32>, height: Option<u32>) -> PyResult<()> {
        let cmd = match (width, height) {
            (Some(w), Some(h)) => format!("viewport {}, {}", w, h),
            (Some(w), None) => format!("viewport {}", w),
            _ => "viewport".to_string(),
        };
        self.execute_cmd(&cmd)
    }

    /// Enable or disable fullscreen mode
    ///
    /// Args:
    ///     toggle: Optional mode ("on", "off", "toggle", "1", "0", "-1")
    ///             If omitted, toggles the current state.
    ///
    /// Examples:
    ///     cmd.full_screen()       # Toggle fullscreen
    ///     cmd.full_screen("on")   # Enable fullscreen
    ///     cmd.full_screen("off")  # Disable fullscreen
    ///
    /// See Also:
    ///     viewport
    #[pyo3(signature = (toggle=None))]
    fn full_screen(&self, toggle: Option<&str>) -> PyResult<()> {
        let cmd = match toggle {
            Some(t) => format!("full_screen {}", t),
            None => "full_screen".to_string(),
        };
        self.execute_cmd(&cmd)
    }

    /// Set the center of rotation
    ///
    /// Args:
    ///     selection: Selection expression or name (default: all)
    ///     object: Object name (optional)
    ///     position: Explicit position [x, y, z] (optional)
    ///     state: State to use for coordinates (default: 0)
    ///
    /// Examples:
    ///     cmd.origin("chain A")
    ///     cmd.origin("resi 100")
    ///
    /// See Also:
    ///     zoom, orient, center, reset
    #[pyo3(signature = (selection="all", object=None, position=None, state=0))]
    fn origin(
        &self,
        selection: &str,
        object: Option<&str>,
        position: Option<Vec<f64>>,
        state: i32,
    ) -> PyResult<()> {
        let cmd = if let Some(pos) = position {
            if pos.len() >= 3 {
                format!("origin position=[{}, {}, {}]", pos[0], pos[1], pos[2])
            } else {
                format!("origin {}", selection)
            }
        } else {
            format!("origin {}", selection)
        };
        self.execute_cmd(&cmd)
    }

    /// Save or recall named camera views
    ///
    /// Unlike scenes, views only store camera state (rotation, position,
    /// origin, clipping planes, FOV) without colors, representations, or
    /// frame state.
    ///
    /// Args:
    ///     key: View name, or "*" for all views
    ///     action: "store", "recall", or "clear" (default: recall)
    ///     animate: Animation duration in seconds (default: 0)
    ///
    /// Examples:
    ///     cmd.view("F1", "store")      # Store current view as "F1"
    ///     cmd.view("F1")               # Recall view "F1"
    ///     cmd.view("F1", "recall", 1)  # Recall with 1 second animation
    ///     cmd.view("F1", "clear")      # Delete view "F1"
    ///     cmd.view("*", "clear")       # Delete all views
    ///
    /// See Also:
    ///     scene, get_view, set_view
    #[pyo3(signature = (key, action="recall", animate=0.0))]
    fn view(&self, key: &str, action: &str, animate: f64) -> PyResult<()> {
        let cmd = if animate > 0.0 {
            format!("view {}, {}, {}", key, action, animate)
        } else {
            format!("view {}, {}", key, action)
        };
        self.execute_cmd(&cmd)
    }

    /// Get the current view matrix
    ///
    /// Returns the current view as a list of 18 floats representing:
    /// - [0-8]: 3x3 rotation matrix (row-major)
    /// - [9-11]: Camera position (x, y, z)
    /// - [12-14]: Origin (center of rotation)
    /// - [15]: Front clipping plane
    /// - [16]: Back clipping plane
    /// - [17]: Field of view
    ///
    /// Args:
    ///     output: Output mode (default: 0, currently ignored)
    ///
    /// Returns:
    ///     List of 18 floats, or None if not available
    ///
    /// Examples:
    ///     view = cmd.get_view()
    ///     # Later restore with:
    ///     cmd.set_view(view)
    ///
    /// See Also:
    ///     set_view, view, scene
    #[pyo3(signature = (output=0))]
    fn get_view(&self, output: i32) -> PyResult<Option<Vec<f64>>> {
        // Try to get view via IPC
        self.with_client(|client| client.get_view())
    }

    /// Set the view matrix
    ///
    /// Sets the camera view from a tuple of 18 floats (as returned by get_view).
    ///
    /// Args:
    ///     view: Tuple or list of 18 floats:
    ///         - [0-8]: 3x3 rotation matrix (row-major)
    ///         - [9-11]: Camera position (x, y, z)
    ///         - [12-14]: Origin (center of rotation)
    ///         - [15]: Front clipping plane
    ///         - [16]: Back clipping plane
    ///         - [17]: Field of view
    ///
    /// Examples:
    ///     # Save and restore view
    ///     view = cmd.get_view()
    ///     # ... make changes ...
    ///     cmd.set_view(view)
    ///
    ///     # Set specific view
    ///     cmd.set_view([
    ///         1.0, 0.0, 0.0,
    ///         0.0, 1.0, 0.0,
    ///         0.0, 0.0, 1.0,
    ///         0.0, 0.0, -50.0,
    ///         0.0, 0.0, 0.0,
    ///         10.0, 100.0, 20.0
    ///     ])
    ///
    /// See Also:
    ///     get_view, view, scene
    #[pyo3(signature = (view))]
    fn set_view(&self, view: Vec<f64>) -> PyResult<()> {
        if view.len() < 18 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "view must have at least 18 elements",
            ));
        }

        // Format as PyMOL-style set_view command
        let cmd = format!(
            "set_view ({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {})",
            view[0], view[1], view[2],
            view[3], view[4], view[5],
            view[6], view[7], view[8],
            view[9], view[10], view[11],
            view[12], view[13], view[14],
            view[15], view[16], view[17],
        );
        self.execute_cmd(&cmd)
    }

    // =========================================================================
    // Object Commands
    // =========================================================================

    /// Delete an object or selection
    ///
    /// Args:
    ///     name: Object or selection name to delete
    #[pyo3(signature = (name))]
    fn delete(&self, name: &str) -> PyResult<()> {
        let cmd = format!("delete {}", name);
        self.execute_cmd(&cmd)
    }

    /// List all object names
    fn get_names(&self) -> PyResult<Vec<String>> {
        self.with_client(|client| client.get_names())
    }

    /// Enable an object (make visible)
    #[pyo3(signature = (name))]
    fn enable(&self, name: &str) -> PyResult<()> {
        let cmd = format!("enable {}", name);
        self.execute_cmd(&cmd)
    }

    /// Disable an object (make invisible)
    #[pyo3(signature = (name))]
    fn disable(&self, name: &str) -> PyResult<()> {
        let cmd = format!("disable {}", name);
        self.execute_cmd(&cmd)
    }

    /// Toggle object visibility
    ///
    /// Args:
    ///     name: Object name to toggle
    ///
    /// See Also:
    ///     enable, disable
    #[pyo3(signature = (name))]
    fn toggle(&self, name: &str) -> PyResult<()> {
        let cmd = format!("toggle {}", name);
        self.execute_cmd(&cmd)
    }

    /// Rename an object
    ///
    /// Args:
    ///     old_name: Current object name
    ///     new_name: New object name
    ///
    /// See Also:
    ///     delete, create, copy
    #[pyo3(signature = (old_name, new_name))]
    fn set_name(&self, old_name: &str, new_name: &str) -> PyResult<()> {
        let cmd = format!("set_name {}, {}", old_name, new_name);
        self.execute_cmd(&cmd)
    }

    /// Rename an object (alias for set_name)
    ///
    /// Args:
    ///     old_name: Current object name
    ///     new_name: New object name
    ///
    /// See Also:
    ///     set_name, delete
    #[pyo3(signature = (old_name, new_name))]
    fn rename(&self, old_name: &str, new_name: &str) -> PyResult<()> {
        self.set_name(old_name, new_name)
    }

    /// Show visual indicator for a selection
    ///
    /// Displays pink dots on atoms matching the selection to help
    /// visualize which atoms are selected.
    ///
    /// Args:
    ///     selection: Selection expression (default: all)
    ///
    /// See Also:
    ///     select, deselect
    #[pyo3(signature = (selection="all"))]
    fn indicate(&self, selection: &str) -> PyResult<()> {
        let cmd = format!("indicate {}", selection);
        self.execute_cmd(&cmd)
    }

    // =========================================================================
    // Transform Commands
    // =========================================================================

    /// Translate atomic coordinates
    ///
    /// Modifies the actual coordinates of atoms (unlike camera movement).
    ///
    /// Args:
    ///     vector: Translation vector [x, y, z]
    ///     selection: Atoms to translate (default: all)
    ///     state: State to modify (-1=current, 0=all, >0=specific)
    ///     camera: Is vector in camera coordinates? (1=yes, 0=no)
    ///
    /// Examples:
    ///     cmd.translate([1, 0, 0], "name CA")
    ///     cmd.translate([0, 5, 0], "chain A")
    ///     cmd.translate([0, 0, 10], "organic")
    ///
    /// See Also:
    ///     rotate, origin, move
    #[pyo3(signature = (vector, selection="all", state=-1, camera=1))]
    fn translate(
        &self,
        vector: Vec<f64>,
        selection: &str,
        state: i32,
        camera: i32,
    ) -> PyResult<()> {
        if vector.len() < 3 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "vector must have at least 3 elements [x, y, z]",
            ));
        }
        let cmd = format!(
            "translate [{}, {}, {}], {}, {}, {}",
            vector[0], vector[1], vector[2], selection, state, camera
        );
        self.execute_cmd(&cmd)
    }

    /// Rotate atomic coordinates
    ///
    /// Modifies the actual coordinates of atoms (unlike camera rotation).
    ///
    /// Args:
    ///     axis: Axis to rotate about - either:
    ///         - String: "x", "y", or "z" for principal axes
    ///         - List: [ax, ay, az] for arbitrary axis vector
    ///     angle: Degrees of rotation
    ///     selection: Atoms to rotate (default: all)
    ///     state: State to modify (-1=current, 0=all, >0=specific)
    ///     camera: Is axis in camera coordinates? (1=yes, 0=no)
    ///     origin: Center of rotation [x, y, z] (default: view origin)
    ///
    /// Examples:
    ///     cmd.rotate("x", 45, "all")
    ///     cmd.rotate("y", 90, "chain A")
    ///     cmd.rotate([1, 1, 0], 45, "organic")  # Rotate about diagonal axis
    ///
    /// See Also:
    ///     translate, turn, origin
    #[pyo3(signature = (axis, angle, selection="all", state=-1, camera=1, origin=None))]
    fn rotate(
        &self,
        axis: &Bound<'_, PyAny>,
        angle: f64,
        selection: &str,
        state: i32,
        camera: i32,
        origin: Option<Vec<f64>>,
    ) -> PyResult<()> {
        // Parse axis - either string ("x", "y", "z") or vector [ax, ay, az]
        let axis_str = if let Ok(s) = axis.extract::<String>() {
            s
        } else if let Ok(v) = axis.extract::<Vec<f64>>() {
            if v.len() < 3 {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "axis vector must have at least 3 elements [ax, ay, az]",
                ));
            }
            format!("[{}, {}, {}]", v[0], v[1], v[2])
        } else {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "axis must be a string ('x', 'y', 'z') or a list [ax, ay, az]",
            ));
        };

        let mut cmd = format!("rotate {}, {}, {}, {}, {}", axis_str, angle, selection, state, camera);

        if let Some(orig) = origin {
            if orig.len() >= 3 {
                cmd.push_str(&format!(", origin=[{}, {}, {}]", orig[0], orig[1], orig[2]));
            }
        }

        self.execute_cmd(&cmd)
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
    fn set(
        &self,
        name: &str,
        value: &Bound<'_, PyAny>,
        selection: Option<&str>,
        quiet: bool,
    ) -> PyResult<()> {
        // Convert value to string for IPC
        let value_str = if let Ok(b) = value.extract::<bool>() {
            if b { "on" } else { "off" }.to_string()
        } else if let Ok(i) = value.extract::<i64>() {
            i.to_string()
        } else if let Ok(f) = value.extract::<f64>() {
            f.to_string()
        } else if let Ok(s) = value.extract::<String>() {
            s
        } else {
            format!("{:?}", value)
        };

        let cmd = if let Some(sel) = selection {
            format!("set {}, {}, {}", name, value_str, sel)
        } else {
            format!("set {}, {}", name, value_str)
        };

        self.execute_cmd_opts(&cmd, quiet)
    }

    /// Get a setting value
    ///
    /// Args:
    ///     name: Setting name
    ///     selection: Optional selection (for per-atom settings)
    ///     state: State to query (default: 0)
    ///
    /// Note:
    ///     Currently returns None as structured return values require
    ///     IPC protocol enhancement. Use cmd.do("get setting_name") to
    ///     see the value printed to output.
    ///
    /// See Also:
    ///     set, unset
    #[pyo3(signature = (name, selection=None, state=0))]
    fn get(&self, name: &str, selection: Option<&str>, state: i32) -> PyResult<Option<String>> {
        let cmd = if let Some(sel) = selection {
            format!("get {}, {}", name, sel)
        } else {
            format!("get {}", name)
        };
        self.execute_cmd(&cmd)?;
        // TODO: Return actual value once IPC protocol supports it
        Ok(None)
    }

    /// Restore a setting to its default value
    ///
    /// Args:
    ///     name: Setting name
    ///     selection: Optional selection (for per-atom settings)
    ///     state: State to modify (default: 0)
    ///
    /// See Also:
    ///     set, get
    #[pyo3(signature = (name, selection=None, state=0))]
    fn unset(&self, name: &str, selection: Option<&str>, state: i32) -> PyResult<()> {
        let cmd = if let Some(sel) = selection {
            format!("unset {}, {}", name, sel)
        } else {
            format!("unset {}", name)
        };
        self.execute_cmd(&cmd)
    }

    // =========================================================================
    // Command Execution
    // =========================================================================

    /// Execute a PyMOL command string
    ///
    /// Args:
    ///     command: Command string (e.g., "load file.pdb; show cartoon")
    ///     quiet: Suppress command echo and info output (default: False)
    #[pyo3(name = "do", signature = (command, quiet=false))]
    fn do_(&self, command: &str, quiet: bool) -> PyResult<()> {
        self.execute_cmd_opts(command, quiet)
    }

    // =========================================================================
    // Image Output
    // =========================================================================

    /// Save current view as PNG image
    ///
    /// Args:
    ///     filename: Output file path (will add .png extension if missing)
    ///     width: Image width in pixels (default: current window size)
    ///     height: Image height in pixels (default: current window size)
    ///     ray: Whether to ray-trace (not yet implemented, ignored)
    ///     quiet: Suppress output (default: True)
    #[pyo3(signature = (filename, width=None, height=None, ray=0, quiet=true))]
    fn png(
        &self,
        filename: &str,
        width: Option<u32>,
        height: Option<u32>,
        ray: i32,
        quiet: bool,
    ) -> PyResult<()> {
        let cmd = match (width, height) {
            (Some(w), Some(h)) => format!("png {}, {}, {}", filename, w, h),
            (Some(w), None) => format!("png {}, {}", filename, w),
            _ => format!("png {}", filename),
        };
        self.execute_cmd_opts(&cmd, quiet)
    }

    /// Refresh the display
    fn refresh(&self) -> PyResult<()> {
        self.execute_cmd("refresh")
    }

    /// Reset PyMOL to initial state
    ///
    /// Args:
    ///     what: What to reinitialize (default: "everything")
    ///         "everything" - complete reset
    ///         "settings" - reset settings only
    ///     object: Object to reinitialize (optional)
    ///
    /// See Also:
    ///     delete, quit
    #[pyo3(signature = (what="everything", object=None))]
    fn reinitialize(&self, what: &str, object: Option<&str>) -> PyResult<()> {
        let cmd = if let Some(obj) = object {
            format!("reinitialize {}, {}", what, obj)
        } else {
            format!("reinitialize {}", what)
        };
        self.execute_cmd(&cmd)
    }

    /// Reset PyMOL to initial state (alias for reinitialize)
    ///
    /// See Also:
    ///     reinitialize
    #[pyo3(signature = (what="everything", object=None))]
    fn reinit(&self, what: &str, object: Option<&str>) -> PyResult<()> {
        self.reinitialize(what, object)
    }

    /// Force rebuilding of representations
    ///
    /// Args:
    ///     selection: Selection to rebuild (default: all)
    ///
    /// See Also:
    ///     refresh
    #[pyo3(signature = (selection="all"))]
    fn rebuild(&self, selection: &str) -> PyResult<()> {
        let cmd = format!("rebuild {}", selection);
        self.execute_cmd(&cmd)
    }

    /// Perform ray-tracing
    ///
    /// Ray tracing produces high-quality images with proper shadows, lighting,
    /// and transparency effects.
    ///
    /// Args:
    ///     width: Width in pixels (default: current window width)
    ///     height: Height in pixels (default: current window height)
    ///     antialias: Antialiasing level 1-4 (default: 1 = no AA)
    ///         1 = no antialiasing
    ///         2 = 2x2 supersampling
    ///         3 = 3x3 supersampling
    ///         4 = 4x4 supersampling
    ///     filename: Output file path (optional, displays if not provided)
    ///     quiet: Suppress feedback (default: False)
    ///
    /// Examples:
    ///     cmd.ray()                          # Raytrace at current resolution
    ///     cmd.ray(1920, 1080)                # Raytrace at 1080p
    ///     cmd.ray(1920, 1080, 2)             # Raytrace at 1080p with 2x2 AA
    ///     cmd.ray(filename="output.png")     # Raytrace and save to file
    ///
    /// See Also:
    ///     png
    #[pyo3(signature = (width=None, height=None, antialias=1, filename=None, quiet=false))]
    fn ray(
        &self,
        width: Option<u32>,
        height: Option<u32>,
        antialias: i32,
        filename: Option<&str>,
        quiet: bool,
    ) -> PyResult<()> {
        // Build command with named arguments for clarity
        let mut cmd = "ray".to_string();

        if let Some(w) = width {
            cmd.push_str(&format!(" width={}", w));
        }
        if let Some(h) = height {
            cmd.push_str(&format!(", height={}", h));
        }
        if antialias != 1 {
            cmd.push_str(&format!(", antialias={}", antialias));
        }
        if let Some(f) = filename {
            cmd.push_str(&format!(", filename={}", f));
        }

        self.execute_cmd(&cmd)
    }

    // =========================================================================
    // Movie Commands
    // =========================================================================

    /// Start movie playback
    ///
    /// Begins playing through the frames of a movie or multi-state object.
    ///
    /// See Also:
    ///     mstop, mpause, mtoggle, frame
    fn mplay(&self) -> PyResult<()> {
        self.execute_cmd("mplay")
    }

    /// Stop movie playback
    ///
    /// Stops playback and resets to frame 1.
    ///
    /// See Also:
    ///     mplay, mpause, mtoggle
    fn mstop(&self) -> PyResult<()> {
        self.execute_cmd("mstop")
    }

    /// Pause movie playback
    ///
    /// Pauses playback at the current frame without resetting position.
    ///
    /// See Also:
    ///     mplay, mstop, mtoggle
    fn mpause(&self) -> PyResult<()> {
        self.execute_cmd("mpause")
    }

    /// Toggle movie playback
    ///
    /// Toggles between play and pause states.
    ///
    /// See Also:
    ///     mplay, mstop, mpause
    fn mtoggle(&self) -> PyResult<()> {
        self.execute_cmd("mtoggle")
    }

    /// Advance to next frame
    ///
    /// Moves forward one frame in the movie.
    ///
    /// See Also:
    ///     backward, frame, rewind, ending
    fn forward(&self) -> PyResult<()> {
        self.execute_cmd("forward")
    }

    /// Go back to previous frame
    ///
    /// Moves backward one frame in the movie.
    ///
    /// See Also:
    ///     forward, frame, rewind, ending
    fn backward(&self) -> PyResult<()> {
        self.execute_cmd("backward")
    }

    /// Go to first frame
    ///
    /// Jumps to the beginning of the movie.
    ///
    /// See Also:
    ///     ending, middle, frame
    fn rewind(&self) -> PyResult<()> {
        self.execute_cmd("rewind")
    }

    /// Go to middle frame
    ///
    /// Jumps to the middle of the movie.
    ///
    /// See Also:
    ///     rewind, ending, frame
    fn middle(&self) -> PyResult<()> {
        self.execute_cmd("middle")
    }

    /// Go to last frame
    ///
    /// Jumps to the end of the movie.
    ///
    /// See Also:
    ///     rewind, middle, frame
    fn ending(&self) -> PyResult<()> {
        self.execute_cmd("ending")
    }

    /// Set number of frames in movie
    ///
    /// Args:
    ///     specification: Number of frames or frame specification string
    ///
    /// Examples:
    ///     cmd.mset("100")      # Create 100 frames
    ///     cmd.mset("1 x30")    # 30 frames of state 1
    ///     cmd.mset("1 -60")    # Frames for states 1-60
    ///
    /// See Also:
    ///     mplay, frame
    #[pyo3(signature = (specification))]
    fn mset(&self, specification: &str) -> PyResult<()> {
        let cmd = format!("mset {}", specification);
        self.execute_cmd(&cmd)
    }

    /// Go to specific frame
    ///
    /// Args:
    ///     frame_number: Frame number (1-indexed)
    ///
    /// See Also:
    ///     forward, backward, rewind, ending
    #[pyo3(signature = (frame_number))]
    fn frame(&self, frame_number: i32) -> PyResult<()> {
        let cmd = format!("frame {}", frame_number);
        self.execute_cmd(&cmd)
    }

    /// Toggle Y-axis rocking animation
    ///
    /// When enabled, the view continuously rocks back and forth around
    /// the Y-axis. This is useful for presentations and visual inspection.
    ///
    /// Args:
    ///     mode: Optional mode ("on", "off", "1", "0"). If omitted, toggles.
    ///
    /// Examples:
    ///     cmd.rock()       # Toggle rock mode
    ///     cmd.rock("on")   # Enable rock mode
    ///     cmd.rock("off")  # Disable rock mode
    ///
    /// See Also:
    ///     mplay, turn
    #[pyo3(signature = (mode=None))]
    fn rock(&self, mode: Option<&str>) -> PyResult<()> {
        let cmd = match mode {
            Some(m) => format!("rock {}", m),
            None => "rock".to_string(),
        };
        self.execute_cmd(&cmd)
    }

    // =========================================================================
    // Scene Commands
    // =========================================================================

    /// Save or recall named scenes
    ///
    /// Scenes store the complete visual state including:
    /// - Camera view (rotation, position, clipping)
    /// - Object visibility (enabled/disabled)
    /// - Representations (show/hide states)
    /// - Colors
    /// - Frame number
    ///
    /// Args:
    ///     key: Scene name, or "*" for all scenes
    ///     action: Action to perform (default: recall)
    ///         "store" - Save current state as scene
    ///         "recall" - Restore saved scene
    ///         "delete" - Delete a scene
    ///         "clear" - Delete all scenes (with key="*")
    ///         "list" - List all scenes
    ///         "rename" - Rename a scene (requires new_key)
    ///     message: Optional message to store with scene
    ///     view: Store/recall view (default: True)
    ///     color: Store/recall colors (default: True)
    ///     rep: Store/recall representations (default: True)
    ///     frame: Store/recall frame number (default: True)
    ///     animate: Animation duration in seconds (default: 0)
    ///
    /// Examples:
    ///     cmd.scene("F1", "store")       # Store current state as "F1"
    ///     cmd.scene("F1")                # Recall scene "F1"
    ///     cmd.scene("F1", "recall", animate=2)  # Recall with 2s animation
    ///     cmd.scene("F1", "delete")      # Delete scene "F1"
    ///     cmd.scene("*", "clear")        # Delete all scenes
    ///     cmd.scene("*", "list")         # List all scenes
    ///
    /// See Also:
    ///     view, get_view, set_view
    #[pyo3(signature = (key, action="recall", message=None, view=true, color=true, rep=true, frame=true, animate=0.0))]
    fn scene(
        &self,
        key: &str,
        action: &str,
        message: Option<&str>,
        view: bool,
        color: bool,
        rep: bool,
        frame: bool,
        animate: f64,
    ) -> PyResult<()> {
        let mut cmd = format!("scene {}, {}", key, action);

        if let Some(msg) = message {
            cmd.push_str(&format!(", message={}", msg));
        }
        if !view {
            cmd.push_str(", view=0");
        }
        if !color {
            cmd.push_str(", color=0");
        }
        if !rep {
            cmd.push_str(", rep=0");
        }
        if !frame {
            cmd.push_str(", frame=0");
        }
        if animate > 0.0 {
            cmd.push_str(&format!(", animate={}", animate));
        }

        self.execute_cmd(&cmd)
    }

    // =========================================================================
    // GUI Control
    // =========================================================================

    /// Show the GUI window (make it visible)
    ///
    /// This makes the pymol-rs window visible. The server continues running
    /// in the background regardless of window visibility.
    ///
    /// Example:
    ///     cmd.load("protein.pdb")
    ///     cmd.show_gui()  # Window appears
    ///     cmd.hide_gui()  # Window hidden, server still running
    fn show_gui(&self) -> PyResult<()> {
        self.with_client(|client| client.show_window())
    }

    /// Hide the GUI window (make it invisible)
    ///
    /// The server continues running in headless mode. Commands are still
    /// executed and state is preserved.
    fn hide_gui(&self) -> PyResult<()> {
        self.with_client(|client| client.hide_window())
    }

    /// Check if connected to the server
    ///
    /// Returns:
    ///     True if connected, False otherwise
    fn is_connected(&self) -> bool {
        self.with_client(|client| client.ping())
            .ok()
            .unwrap_or(false)
    }

    /// Quit the pymol-rs server
    ///
    /// This completely shuts down the server. After calling this,
    /// use cmd.start() to restart the server.
    fn quit(&mut self) -> PyResult<()> {
        // Stop the listener first
        if let Ok(mut listener_guard) = self.listener.lock() {
            if let Some(ref mut listener) = *listener_guard {
                listener.stop();
            }
            *listener_guard = None;
        }

        // Drop the connection (this triggers cleanup via Connection::Drop)
        self.connection = None;

        Ok(())
    }

    /// Start or reconnect to a pymol-rs server
    ///
    /// If a server is already running and connected, this does nothing.
    /// Otherwise, it spawns a new headless server and connects to it.
    ///
    /// This is useful after calling cmd.quit() or if the server was
    /// closed via the GUI.
    ///
    /// Example:
    ///     cmd.quit()  # Close the server
    ///     # ... do other things ...
    ///     cmd.start()  # Start a new server
    ///     cmd.load("protein.pdb")  # Works again!
    fn start(&mut self) -> PyResult<()> {
        // Check if already connected
        if self.is_connected() {
            log::info!("Server is already running");
            return Ok(());
        }

        log::info!("Starting new pymol-rs server");

        // Stop the old listener if any
        if let Ok(mut listener_guard) = self.listener.lock() {
            if let Some(ref mut listener) = *listener_guard {
                listener.stop();
            }
            *listener_guard = None;
        }

        // Establish a new connection
        let connection = establish_connection().map_err(|e| {
            pyo3::exceptions::PyRuntimeError::new_err(format!(
                "Failed to start pymol-rs server: {}",
                e
            ))
        })?;

        self.connection = Some(connection);

        // Re-register runpy command
        self.register_runpy()?;

        log::info!("Server started successfully");
        Ok(())
    }

    // =========================================================================
    // Extended Commands
    // =========================================================================

    /// Register a Python function as a command
    ///
    /// The function will be callable from:
    /// - Python: cmd.do_("highlight chain A")
    /// - GUI command line: type "highlight chain A" (with autocomplete!)
    ///
    /// Example:
    ///     def highlight(selection="all"):
    ///         cmd.show("sticks", selection)
    ///         cmd.color("yellow", selection)
    ///
    ///     cmd.extend("highlight", highlight)
    #[pyo3(signature = (name, function))]
    fn extend(&self, py: Python<'_>, name: &str, function: Py<PyAny>) -> PyResult<()> {
        // Start the listener if not already running (needed for callbacks)
        self.ensure_listener_running();
        
        // Extract help text from docstring
        let help = function
            .getattr(py, "__doc__")
            .ok()
            .and_then(|d| d.extract::<Option<String>>(py).ok())
            .flatten();

        // Store callback
        {
            let mut commands = self.extended_commands.lock().map_err(|e| {
                pyo3::exceptions::PyRuntimeError::new_err(format!(
                    "Failed to lock commands: {}",
                    e
                ))
            })?;
            commands.register(name.to_string(), function);
        }

        // Register with server for autocomplete
        self.with_client(|client| client.register_command(name, help.as_deref()))?;

        Ok(())
    }

    /// Unregister an extended command
    #[pyo3(signature = (name))]
    fn unextend(&self, name: &str) -> PyResult<()> {
        // Remove callback
        {
            let mut commands = self.extended_commands.lock().map_err(|e| {
                pyo3::exceptions::PyRuntimeError::new_err(format!(
                    "Failed to lock commands: {}",
                    e
                ))
            })?;
            commands.unregister(name);
        }

        // Unregister from server
        self.with_client(|client| client.unregister_command(name))?;

        Ok(())
    }

    // =========================================================================
    // Misc
    // =========================================================================

    /// Set default silent mode for all commands
    ///
    /// When silent mode is enabled, commands will not echo to the GUI
    /// output panel and info/warning messages will be suppressed.
    /// Individual commands can still override this with their `quiet` parameter.
    ///
    /// Args:
    ///     silent: True to suppress output, False to show output (default)
    ///
    /// Examples:
    ///     cmd.set_silent(True)   # All subsequent commands run silently
    ///     cmd.load("protein.pdb")  # No echo in GUI output
    ///     cmd.set_silent(False)  # Restore normal output
    #[pyo3(signature = (silent))]
    fn set_silent(&mut self, silent: bool) {
        self.silent = silent;
    }

    /// Get the current default silent mode
    ///
    /// Returns:
    ///     True if silent mode is enabled, False otherwise
    fn get_silent(&self) -> bool {
        self.silent
    }

    /// Print version information
    fn version(&self) -> String {
        format!("PyMOL-RS {}", env!("CARGO_PKG_VERSION"))
    }

    fn __repr__(&self) -> String {
        let connected = self.is_connected();
        let server_info = self.connection
            .as_ref()
            .map(|c| if c.owns_server() { "spawned" } else { "external" })
            .unwrap_or("disconnected");
        format!(
            "Cmd(connected={}, server={})",
            connected, server_info
        )
    }
}
