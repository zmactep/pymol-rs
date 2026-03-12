//! Python Worker Thread
//!
//! Runs the Python interpreter in a dedicated background thread so that
//! `eval()` / `exec_file()` never block the main (UI) thread.
//!
//! All callers submit [`WorkItem`]s via the [`WorkerHandle`] and results
//! are returned through an `mpsc` channel, drained in `PythonHandler::poll()`.

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{self, Receiver, Sender};
use std::sync::{Arc, Mutex};
use std::thread;

use pyo3::prelude::*;

use crate::backend::{PluginBackend, SharedStateHandle};
use crate::engine::PythonEngine;

// =============================================================================
// StreamingWriter — real-time stdout/stderr forwarding
// =============================================================================

/// A Python file-like object that sends each `write()` call through an mpsc
/// channel, enabling real-time output streaming from long-running scripts.
#[pyclass]
struct StreamingWriter {
    tx: Sender<WorkResult>,
    origin: WorkOrigin,
    /// Line buffer — accumulates text until a newline is seen, then sends
    /// complete lines. This avoids splitting `print("x")` into separate
    /// `"x"` and `"\n"` entries (which would render as extra blank lines).
    buffer: String,
}

#[pymethods]
impl StreamingWriter {
    fn write(&mut self, text: &str) -> PyResult<usize> {
        let len = text.len();
        self.buffer.push_str(text);

        // Flush complete lines
        while let Some(pos) = self.buffer.find('\n') {
            let line: String = self.buffer.drain(..=pos).collect();
            let _ = self.tx.send(WorkResult {
                origin: self.origin,
                result: Ok(line),
            });
        }

        Ok(len)
    }

    fn flush(&mut self) -> PyResult<()> {
        if !self.buffer.is_empty() {
            let text = std::mem::take(&mut self.buffer);
            let _ = self.tx.send(WorkResult {
                origin: self.origin,
                result: Ok(text),
            });
        }
        Ok(())
    }
}

// =============================================================================
// Types
// =============================================================================

/// Output entry from Python execution.
pub struct OutputEntry {
    pub text: String,
    pub is_error: bool,
}

/// Shared buffer for routing editor results from handler to editor component.
pub type EditorOutputHandle = Arc<Mutex<Vec<OutputEntry>>>;

/// Identifies the origin of a work request (for routing results).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WorkOrigin {
    /// From the `python` / `/` command line.
    Command,
    /// From the editor Run button.
    Editor,
    /// From `run script.py`.
    Script,
    /// Backend installation (one-time setup).
    Setup,
}

/// A unit of work to send to the Python worker thread.
pub enum WorkItem {
    /// Evaluate a code string.
    Eval { code: String, origin: WorkOrigin },
    /// Execute a file by path.
    ExecFile { path: String, origin: WorkOrigin },
    /// Install the PluginBackend into `sys._pymolrs_backend`.
    InstallBackend { shared: SharedStateHandle },
    /// Shut down the worker thread.
    #[allow(dead_code)]
    Shutdown,
}

/// Result sent back from the worker thread.
pub struct WorkResult {
    pub origin: WorkOrigin,
    pub result: Result<String, String>,
}

// =============================================================================
// WorkerHandle
// =============================================================================

/// Handle to the Python worker thread.
///
/// Cheaply cloneable — holds a channel sender and atomic busy flag.
#[derive(Clone)]
pub struct WorkerHandle {
    tx: Sender<WorkItem>,
    busy: Arc<AtomicBool>,
}

impl WorkerHandle {
    /// Submit a work item to the Python worker thread.
    ///
    /// Returns immediately; the result will arrive via the result channel.
    pub fn submit(&self, item: WorkItem) {
        let _ = self.tx.send(item);
    }

    /// Check if the worker is currently executing something.
    pub fn is_busy(&self) -> bool {
        self.busy.load(Ordering::Relaxed)
    }

    /// Request interruption of the currently running Python code.
    ///
    /// Calls `PyErr_SetInterrupt` which raises `KeyboardInterrupt` at the
    /// next bytecode instruction check. Safe to call from any thread.
    /// No-op if the worker is not busy.
    pub fn request_interrupt(&self) {
        if self.is_busy() {
            // Safety: PyErr_SetInterrupt is explicitly documented as safe
            // to call from any thread (it is async-signal-safe).
            unsafe {
                pyo3::ffi::PyErr_SetInterrupt();
            }
        }
    }
}

// =============================================================================
// spawn_worker
// =============================================================================

/// Spawn the Python worker thread.
///
/// Returns the handle (for submitting work) and the receiver (for collecting results).
/// The worker owns the `PythonEngine` exclusively — no `Arc<Mutex<>>` needed.
pub fn spawn_worker() -> (WorkerHandle, Receiver<WorkResult>) {
    let (work_tx, work_rx) = mpsc::channel::<WorkItem>();
    let (result_tx, result_rx) = mpsc::channel::<WorkResult>();
    let busy = Arc::new(AtomicBool::new(false));
    let busy_flag = busy.clone();

    thread::Builder::new()
        .name("python-worker".into())
        .spawn(move || {
            worker_loop(work_rx, result_tx, busy_flag);
        })
        .expect("Failed to spawn Python worker thread");

    let handle = WorkerHandle { tx: work_tx, busy };
    (handle, result_rx)
}

/// Main loop of the Python worker thread.
fn worker_loop(
    work_rx: Receiver<WorkItem>,
    result_tx: Sender<WorkResult>,
    busy: Arc<AtomicBool>,
) {
    let mut engine = PythonEngine::new();

    for item in work_rx {
        match item {
            WorkItem::Shutdown => break,

            WorkItem::InstallBackend { shared } => {
                busy.store(true, Ordering::Relaxed);
                let result = install_backend(&mut engine, &shared);
                let _ = result_tx.send(WorkResult {
                    origin: WorkOrigin::Setup,
                    result,
                });
                busy.store(false, Ordering::Relaxed);
            }

            WorkItem::Eval { code, origin } => {
                busy.store(true, Ordering::Relaxed);
                eval_streaming(&mut engine, &code, origin, &result_tx);
                busy.store(false, Ordering::Relaxed);
            }

            WorkItem::ExecFile { path, origin } => {
                busy.store(true, Ordering::Relaxed);
                match std::fs::read_to_string(&path) {
                    Ok(code) => eval_streaming(&mut engine, &code, origin, &result_tx),
                    Err(e) => {
                        let _ = result_tx.send(WorkResult {
                            origin,
                            result: Err(format!("failed to read {}: {}", path, e)),
                        });
                    }
                }
                busy.store(false, Ordering::Relaxed);
            }
        }
    }

    log::info!("Python worker thread exiting");
}

/// Run Python code with real-time output streaming via [`StreamingWriter`].
///
/// Output is sent through `result_tx` as it is produced. If the code raises
/// an exception, the error is sent as a final `Err` result.
fn eval_streaming(
    engine: &mut PythonEngine,
    code: &str,
    origin: WorkOrigin,
    result_tx: &Sender<WorkResult>,
) {
    let result = Python::attach(|py| -> Result<(), String> {
        let stdout_w = Py::new(
            py,
            StreamingWriter {
                tx: result_tx.clone(),
                origin,
                buffer: String::new(),
            },
        )
        .map_err(|e| e.to_string())?;

        let stderr_w = Py::new(
            py,
            StreamingWriter {
                tx: result_tx.clone(),
                origin,
                buffer: String::new(),
            },
        )
        .map_err(|e| e.to_string())?;

        let result = engine.eval_with_writers(code, stdout_w.bind(py), stderr_w.bind(py));

        // Flush any remaining buffered text (e.g. output without trailing newline)
        let _ = stdout_w.borrow_mut(py).flush();
        let _ = stderr_w.borrow_mut(py).flush();

        result
    });

    if let Err(e) = result {
        let _ = result_tx.send(WorkResult {
            origin,
            result: Err(e),
        });
    }
}

/// Install the PluginBackend into `sys._pymolrs_backend` and auto-import `cmd`.
fn install_backend(engine: &mut PythonEngine, shared: &SharedStateHandle) -> Result<String, String> {
    engine.ensure_init()?;

    let py_backend = Python::attach(|py| -> Result<Py<PyAny>, String> {
        let backend = PluginBackend::new(shared.clone());
        let py_obj = Py::new(py, backend).map_err(|e| e.to_string())?;
        Ok(py_obj.into_any())
    })?;

    engine.set_backend(py_backend)?;

    // Auto-import cmd and stored into __main__ so it's available in the REPL
    match engine.eval("from pymol_rs import cmd; from pymol_rs import stored") {
        Ok(_) => Ok("backend installed, cmd auto-imported".to_string()),
        Err(e) => {
            log::warn!("Python plugin: backend installed but auto-import cmd failed: {}", e);
            Ok("backend installed (cmd auto-import failed)".to_string())
        }
    }
}
