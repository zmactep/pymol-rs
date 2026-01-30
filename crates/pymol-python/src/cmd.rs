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
        })
    }

    /// Execute a command via IPC and return the result
    fn execute_cmd(&self, command: &str) -> PyResult<()> {
        let response = self.with_client(|client| client.execute(command))?;

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
                Some("Execute a Python script: runpy script.py [namespace]"),
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
        self.execute_cmd(&cmd)?;

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
        self.execute_cmd(&cmd)?;
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

        self.execute_cmd(&cmd)
    }

    // =========================================================================
    // Command Execution
    // =========================================================================

    /// Execute a PyMOL command string
    ///
    /// Args:
    ///     command: Command string (e.g., "load file.pdb; show cartoon")
    #[pyo3(name = "do")]
    fn do_(&self, command: &str) -> PyResult<()> {
        self.execute_cmd(command)
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
        self.execute_cmd(&cmd)
    }

    /// Refresh the display
    fn refresh(&self) -> PyResult<()> {
        self.execute_cmd("refresh")
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
