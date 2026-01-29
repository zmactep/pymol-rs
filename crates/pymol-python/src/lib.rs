//! PyMOL-RS Python Bindings
//!
//! This crate provides Python bindings for pymol-rs using PyO3.
//!
//! On module import, automatically connects to an existing pymol-rs IPC server
//! or spawns a new headless server if none is running.

use std::path::PathBuf;

use pyo3::prelude::*;
use pyo3::types::PyModule;

mod cmd;
pub mod connection;
mod convert;
mod error;
pub mod ipc;
mod scripting;

pub mod color;
pub mod io;
pub mod mol;
pub mod selecting;
pub mod settings;

pub use cmd::PyCmd;
pub use connection::{establish_connection, Connection};
pub use error::{PymolError, SelectionError};
pub use mol::{PyAtom, PyBond, PyCoordSet, PyElement, PyObjectMolecule};
pub use color::PyColor;
pub use selecting::PySelectionResult;

/// Initialize logging for the GUI, if not in quiet mode
fn init_gui_logging(quiet: bool, headless: bool, ipc_socket: Option<&PathBuf>) {
    if quiet {
        return;
    }
    
    // Set default log level to info if RUST_LOG is not set
    if std::env::var("RUST_LOG").is_err() {
        std::env::set_var("RUST_LOG", "info");
    }
    let _ = env_logger::try_init();
    
    log::info!("Starting PyMOL-RS GUI from Python");
    if headless {
        log::info!("Headless mode enabled");
    }
    if let Some(socket_path) = ipc_socket {
        log::info!("IPC mode enabled: {:?}", socket_path);
    }
}

/// Run the PyMOL-RS GUI application
///
/// This function blocks until the application is closed.
/// Note: This takes over the main thread and holds the Python GIL.
/// For interactive Python use, prefer using `cmd` with IPC to control
/// a running pymol-rs instance.
///
/// # Arguments
/// * `ipc_socket` - Optional path to the IPC socket for external control
/// * `headless` - If true, run without displaying a window
/// * `files` - List of files to load at startup
/// * `quiet` - If true, suppress log output (used when spawned by Python code)
#[pyfunction]
#[pyo3(signature = (ipc_socket=None, headless=false, files=None, quiet=false))]
fn run_gui(
    ipc_socket: Option<PathBuf>,
    headless: bool,
    files: Option<Vec<String>>,
    quiet: bool,
) -> PyResult<()> {
    use pymol_gui::App;
    use winit::event_loop::{ControlFlow, EventLoop};

    // Initialize logging (consolidated quiet flag handling)
    init_gui_logging(quiet, headless, ipc_socket.as_ref());

    // Create the application
    let mut app = if let Some(ref socket_path) = ipc_socket {
        App::with_ipc(socket_path, headless).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to create app with IPC: {}",
                e
            ))
        })?
    } else {
        App::new(headless)
    };

    // Queue files to load
    if let Some(file_list) = files {
        for file_path in file_list {
            app.queue_load_file(file_path);
        }
    }

    // Create event loop
    let event_loop = EventLoop::new().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to create event loop: {}",
            e
        ))
    })?;
    event_loop.set_control_flow(ControlFlow::Wait);

    // Run the application (this blocks until the window is closed)
    // Note: The GUI event loop must run on the main thread and cannot release
    // the GIL because EventLoop is not Send. This is fine for CLI usage where
    // Python is just the launcher.
    event_loop.run_app(&mut app).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Application error: {}",
            e
        ))
    })
}

/// Create and initialize the global cmd instance
///
/// This is called lazily when the user first accesses `cmd`.
/// It establishes the IPC connection (or spawns a headless server)
/// and registers built-in commands.
#[pyfunction]
fn _create_cmd() -> PyResult<PyCmd> {
    let cmd = PyCmd::new()?;
    cmd.register_runpy()?;
    Ok(cmd)
}

/// Python module initialization
///
/// Note: The `cmd` object is created lazily on first access to avoid
/// connection errors when the module is imported just for the `run_gui` function.
/// Note: Logging is NOT initialized here - it's done lazily in run_gui() or _create_cmd()
/// so that RUST_LOG can be set appropriately first.
#[pymodule]
fn _pymol_rs(py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Register exception types
    m.add("PymolError", py.get_type::<error::PymolError>())?;
    m.add("SelectionError", py.get_type::<error::SelectionError>())?;

    // Register main types
    m.add_class::<PyCmd>()?;
    m.add_class::<PyObjectMolecule>()?;
    m.add_class::<mol::PyAtom>()?;
    m.add_class::<mol::PyBond>()?;
    m.add_class::<mol::PyElement>()?;
    m.add_class::<mol::PyCoordSet>()?;
    m.add_class::<color::PyColor>()?;
    m.add_class::<selecting::PySelectionResult>()?;

    // Register GUI launcher function
    m.add_function(wrap_pyfunction!(run_gui, m)?)?;
    
    // Register the lazy cmd creator function
    m.add_function(wrap_pyfunction!(_create_cmd, m)?)?;

    // Register submodules
    register_submodule(py, m, "mol", mol::register_module)?;
    register_submodule(py, m, "io", io::register_module)?;
    register_submodule(py, m, "selecting", selecting::register_module)?;
    register_submodule(py, m, "color", color::register_module)?;
    register_submodule(py, m, "settings", settings::register_module)?;

    // Add version info
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    Ok(())
}

/// Helper function to register a submodule
fn register_submodule<F>(
    py: Python<'_>,
    parent: &Bound<'_, PyModule>,
    name: &str,
    register_fn: F,
) -> PyResult<()>
where
    F: FnOnce(&Bound<'_, PyModule>) -> PyResult<()>,
{
    let submodule = PyModule::new(py, name)?;
    register_fn(&submodule)?;
    parent.add_submodule(&submodule)?;
    Ok(())
}
