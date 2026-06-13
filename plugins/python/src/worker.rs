//! Python Worker Thread
//!
//! Runs the Python interpreter in a dedicated background thread so that
//! `eval()` / `exec_file()` never block the main (UI) thread.
//!
//! All callers submit [`WorkItem`]s via the [`WorkerHandle`] and results
//! are returned through an `mpsc` channel, drained in `PythonHandler::poll()`.

use std::ffi::CString;
use std::io;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{self, Receiver, Sender};
use std::sync::Arc;
use std::thread;

use pyo3::prelude::*;

use crate::backend::PluginBackend;
use crate::engine::PythonEngine;
use crate::shared::SharedStateHandle;

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
                payload: WorkResultPayload::Output(line),
            });
        }

        Ok(len)
    }

    fn flush(&mut self) -> PyResult<()> {
        if !self.buffer.is_empty() {
            let text = std::mem::take(&mut self.buffer);
            let _ = self.tx.send(WorkResult {
                origin: self.origin,
                payload: WorkResultPayload::Output(text),
            });
        }
        Ok(())
    }
}

#[pyclass]
struct CancellationToken {
    interrupt_requested: Arc<AtomicBool>,
}

#[pymethods]
impl CancellationToken {
    fn check(&self) -> PyResult<()> {
        check_interrupt_requested(&self.interrupt_requested)
    }
}

fn check_interrupt_requested(interrupt_requested: &AtomicBool) -> PyResult<()> {
    if interrupt_requested.load(Ordering::Acquire) {
        Err(pyo3::exceptions::PyRuntimeError::new_err(
            "Python script interrupted",
        ))
    } else {
        Ok(())
    }
}

// =============================================================================
// Types
// =============================================================================

/// Identifies the origin of a work request (for routing results).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WorkOrigin {
    /// From the `python` / `/` command line.
    Command,
    /// From `run script.py`.
    Script,
    /// From the scripting panel.
    Panel,
    /// Backend installation (one-time setup).
    Setup,
}

/// A unit of work to send to the Python worker thread.
pub enum WorkItem {
    /// Evaluate a code string.
    Eval { code: String, origin: WorkOrigin },
    /// Execute a file by path.
    ExecFile { path: String, origin: WorkOrigin },
    /// Install the PluginBackend into `sys._patinae_backend`.
    InstallBackend { shared: SharedStateHandle },
    /// Invoke a Python callable bound to a key (no arguments).
    InvokeKeybindCallback {
        callback: Py<PyAny>,
        origin: WorkOrigin,
    },
    /// Shut down the worker thread.
    #[allow(dead_code)]
    Shutdown,
}

/// Result sent back from the worker thread.
pub struct WorkResult {
    pub origin: WorkOrigin,
    pub payload: WorkResultPayload,
}

/// Typed worker event payload.
pub enum WorkResultPayload {
    /// A stdout/stderr chunk emitted by running Python code.
    Output(String),
    /// A work item completed, with an optional error.
    Finished(Result<(), String>),
    /// One-time backend setup result.
    Setup(Result<String, String>),
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
    interrupt_requested: Arc<AtomicBool>,
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

    /// Request cooperative cancellation of the running script.
    pub fn request_interrupt(&self) {
        self.interrupt_requested.store(true, Ordering::Release);
    }
}

// =============================================================================
// spawn_worker
// =============================================================================

/// Spawn the Python worker thread.
///
/// Returns the handle (for submitting work) and the receiver (for collecting results).
/// The worker owns the `PythonEngine` exclusively — no `Arc<Mutex<>>` needed.
pub fn spawn_worker(interrupt_requested: Arc<AtomicBool>) -> (WorkerHandle, Receiver<WorkResult>) {
    let (work_tx, work_rx) = mpsc::channel::<WorkItem>();
    let (result_tx, result_rx) = mpsc::channel::<WorkResult>();
    let busy = Arc::new(AtomicBool::new(false));
    let busy_flag = busy.clone();
    let worker_interrupt_requested = interrupt_requested.clone();

    thread::Builder::new()
        .name("python-worker".into())
        .spawn(move || {
            worker_loop(work_rx, result_tx, busy_flag, worker_interrupt_requested);
        })
        .expect("Failed to spawn Python worker thread");

    let handle = WorkerHandle {
        tx: work_tx,
        busy,
        interrupt_requested,
    };
    (handle, result_rx)
}

/// Main loop of the Python worker thread.
fn worker_loop(
    work_rx: Receiver<WorkItem>,
    result_tx: Sender<WorkResult>,
    busy: Arc<AtomicBool>,
    interrupt_requested: Arc<AtomicBool>,
) {
    let mut engine = PythonEngine::new();
    let script_reader = ProcessScriptReader;

    for item in work_rx {
        match item {
            WorkItem::Shutdown => break,

            WorkItem::InstallBackend { shared } => {
                busy.store(true, Ordering::Relaxed);
                let result = install_backend(&mut engine, &shared);
                let _ = result_tx.send(WorkResult {
                    origin: WorkOrigin::Setup,
                    payload: WorkResultPayload::Setup(result),
                });
                busy.store(false, Ordering::Relaxed);
            }

            WorkItem::InvokeKeybindCallback { callback, origin } => {
                busy.store(true, Ordering::Relaxed);
                let result = Python::attach(|py| -> Result<String, String> {
                    install_interrupt_trace(py, &interrupt_requested)?;
                    let result = callback.call0(py).map_err(|e| e.to_string());
                    clear_interrupt_trace(py);
                    result?;
                    Ok(String::new())
                });
                let _ = result_tx.send(WorkResult {
                    origin,
                    payload: WorkResultPayload::Finished(result.map(|_| ())),
                });
                busy.store(false, Ordering::Relaxed);
            }

            WorkItem::Eval { code, origin } => {
                busy.store(true, Ordering::Relaxed);
                eval_streaming(&mut engine, &code, origin, &result_tx, &interrupt_requested);
                busy.store(false, Ordering::Relaxed);
            }

            WorkItem::ExecFile { path, origin } => {
                busy.store(true, Ordering::Relaxed);
                match read_script_file(&script_reader, &path) {
                    Ok(code) => {
                        eval_streaming(&mut engine, &code, origin, &result_tx, &interrupt_requested)
                    }
                    Err(e) => {
                        let _ = result_tx.send(WorkResult {
                            origin,
                            payload: WorkResultPayload::Finished(Err(e)),
                        });
                    }
                }
                busy.store(false, Ordering::Relaxed);
            }
        }
    }

    log::info!("Python worker thread exiting");
}

trait ScriptReader {
    fn read_to_string(&self, path: &str) -> io::Result<String>;
}

struct ProcessScriptReader;

impl ScriptReader for ProcessScriptReader {
    fn read_to_string(&self, path: &str) -> io::Result<String> {
        std::fs::read_to_string(path)
    }
}

fn read_script_file(reader: &impl ScriptReader, path: &str) -> Result<String, String> {
    reader
        .read_to_string(path)
        .map_err(|e| format!("failed to read {}: {}", path, e))
}

fn install_interrupt_trace(
    py: Python<'_>,
    interrupt_requested: &Arc<AtomicBool>,
) -> Result<(), String> {
    let token = Py::new(
        py,
        CancellationToken {
            interrupt_requested: interrupt_requested.clone(),
        },
    )
    .map_err(|e| e.to_string())?;

    let globals = py.import("__main__").map_err(|e| e.to_string())?.dict();
    globals
        .set_item("__patinae_interrupt_token", token.bind(py))
        .map_err(|e| e.to_string())?;

    let trace_code = CString::new(
        r#"def __patinae_interrupt_trace(frame, event, arg):
    if frame.f_code.co_name == "<module>":
        __patinae_interrupt_token.check()
    return __patinae_interrupt_trace
"#,
    )
    .map_err(|e| e.to_string())?;
    py.run(trace_code.as_c_str(), Some(&globals), None)
        .map_err(|e| e.to_string())?;

    let trace_fn = globals
        .get_item("__patinae_interrupt_trace")
        .map_err(|e| e.to_string())?
        .ok_or_else(|| "failed to install Python interrupt trace".to_string())?;
    py.import("sys")
        .map_err(|e| e.to_string())?
        .call_method1("settrace", (trace_fn,))
        .map_err(|e| e.to_string())?;

    Ok(())
}

fn clear_interrupt_trace(py: Python<'_>) {
    if let Ok(sys) = py.import("sys") {
        let _ = sys.call_method1("settrace", (py.None(),));
    }
    if let Ok(main) = py.import("__main__") {
        let globals = main.dict();
        let _ = globals.del_item("__patinae_interrupt_token");
        let _ = globals.del_item("__patinae_interrupt_trace");
    }
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
    interrupt_requested: &Arc<AtomicBool>,
) {
    interrupt_requested.store(false, Ordering::Release);

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

        install_interrupt_trace(py, interrupt_requested)?;
        let result = engine.eval_with_writers(code, stdout_w.bind(py), stderr_w.bind(py));

        // Flush any remaining buffered text (e.g. output without trailing newline)
        let _ = stdout_w.borrow_mut(py).flush();
        let _ = stderr_w.borrow_mut(py).flush();
        clear_interrupt_trace(py);

        result
    });
    interrupt_requested.store(false, Ordering::Release);

    let _ = result_tx.send(WorkResult {
        origin,
        payload: WorkResultPayload::Finished(result),
    });
}

/// Install the PluginBackend into `sys._patinae_backend` and auto-import `cmd`.
fn install_backend(
    engine: &mut PythonEngine,
    shared: &SharedStateHandle,
) -> Result<String, String> {
    engine.ensure_init()?;

    let py_backend = Python::attach(|py| -> Result<Py<PyAny>, String> {
        let backend = PluginBackend::new(shared.clone());
        let py_obj = Py::new(py, backend).map_err(|e| e.to_string())?;
        Ok(py_obj.into_any())
    })?;

    engine.set_backend(py_backend)?;

    // Auto-import cmd and stored into __main__ so it's available in the REPL
    match engine.eval("from patinae import cmd; from patinae import stored") {
        Ok(_) => Ok("backend installed, cmd auto-imported".to_string()),
        Err(e) => {
            log::warn!(
                "Python plugin: backend installed but auto-import cmd failed: {}",
                e
            );
            Ok("backend installed (cmd auto-import failed)".to_string())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{Duration, Instant};

    #[derive(Clone, Copy)]
    enum FakeRead {
        Ok(&'static str),
        Err(io::ErrorKind),
    }

    struct FakeScriptReader {
        result: FakeRead,
    }

    impl ScriptReader for FakeScriptReader {
        fn read_to_string(&self, _path: &str) -> io::Result<String> {
            match self.result {
                FakeRead::Ok(code) => Ok(code.to_string()),
                FakeRead::Err(kind) => Err(io::Error::new(kind, "fake read failure")),
            }
        }
    }

    #[test]
    fn script_reader_returns_code_without_process_fs() {
        let reader = FakeScriptReader {
            result: FakeRead::Ok("print('ok')"),
        };

        let code = read_script_file(&reader, "/fake/script.py").expect("script should read");

        assert_eq!(code, "print('ok')");
    }

    #[test]
    fn script_reader_formats_read_failures() {
        let reader = FakeScriptReader {
            result: FakeRead::Err(io::ErrorKind::NotFound),
        };

        let err = read_script_file(&reader, "/fake/missing.py").expect_err("read should fail");

        assert!(
            err.contains("failed to read /fake/missing.py"),
            "unexpected error: {err}"
        );
        assert!(err.contains("fake read failure"), "unexpected error: {err}");
    }

    #[test]
    fn request_interrupt_stops_running_python_code() {
        let interrupt_requested = Arc::new(AtomicBool::new(false));
        let (worker, rx) = spawn_worker(interrupt_requested);
        worker.submit(WorkItem::Eval {
            code: "print('started')\nwhile True:\n    pass\n".to_string(),
            origin: WorkOrigin::Panel,
        });

        let startup_deadline = Instant::now() + Duration::from_secs(5);
        let mut saw_started = false;
        while !saw_started {
            match rx.recv_timeout(Duration::from_millis(100)) {
                Ok(result) => match result.payload {
                    WorkResultPayload::Output(output) => {
                        saw_started = output.contains("started");
                    }
                    WorkResultPayload::Finished(result) => {
                        panic!("Python worker finished before interrupt: {result:?}");
                    }
                    WorkResultPayload::Setup(_) => {}
                },
                Err(mpsc::RecvTimeoutError::Timeout) => {
                    assert!(
                        Instant::now() < startup_deadline,
                        "Python worker did not start executing before timeout"
                    );
                }
                Err(e) => panic!("Python worker result channel closed: {e}"),
            }
        }

        worker.request_interrupt();

        let finished_deadline = Instant::now() + Duration::from_secs(5);
        let mut interrupt_error = None;
        while Instant::now() < finished_deadline {
            match rx.recv_timeout(Duration::from_millis(100)) {
                Ok(result) => {
                    if let WorkResultPayload::Finished(result) = result.payload {
                        interrupt_error = Some(result.expect_err(
                            "interrupted Python code should finish with an interrupt error",
                        ));
                        break;
                    }
                }
                Err(mpsc::RecvTimeoutError::Timeout) => {}
                Err(e) => panic!("Python worker result channel closed: {e}"),
            }
        }

        worker.submit(WorkItem::Shutdown);

        let err = interrupt_error.expect("Python code did not stop after interrupt");
        assert!(
            err.contains("Python script interrupted"),
            "expected interrupt error, got {err:?}"
        );
    }
}
