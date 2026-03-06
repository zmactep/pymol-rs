//! Python Interpreter Engine
//!
//! Manages the embedded Python interpreter, providing `eval` and `exec_file` operations.
//! Uses PyO3 with `auto-initialize` — the interpreter starts on first GIL acquisition.
//!
//! Before the first GIL acquisition, we configure `PYTHONHOME` so that the
//! embedded interpreter can find its standard library (encodings, etc.).

use std::ffi::CString;
use std::sync::atomic::{AtomicU8, Ordering};
use std::sync::Once;

use pyo3::prelude::*;

// ---------------------------------------------------------------------------
// One-time Python environment configuration
// ---------------------------------------------------------------------------

static PYTHON_CONFIGURED: Once = Once::new();
static PYTHON_CONFIG_STATUS: AtomicU8 = AtomicU8::new(ConfigStatus::Pending as u8);

/// The Python major.minor version that PyO3 was compiled against.
/// Set by `build.rs` from PyO3's build config.
const PYO3_PYTHON_VERSION: &str = env!("PYMOLRS_PYO3_PYTHON_VERSION");

/// Result of the one-time Python environment configuration.
#[repr(u8)]
enum ConfigStatus {
    /// Configuration has not run yet.
    Pending = 0,
    /// A matching Python was found and PYTHONHOME is set.
    Ready = 1,
    /// No matching Python was found; the plugin is disabled.
    NoPython = 2,
}

impl ConfigStatus {
    fn load() -> Self {
        match PYTHON_CONFIG_STATUS.load(Ordering::Relaxed) {
            1 => Self::Ready,
            2 => Self::NoPython,
            _ => Self::Pending,
        }
    }

    fn store(self) {
        PYTHON_CONFIG_STATUS.store(self as u8, Ordering::Relaxed);
    }
}

/// Build a list of Python executables to probe, preferring venv if available.
fn python_candidates() -> Vec<std::path::PathBuf> {
    let mut candidates = Vec::new();

    // 1. VIRTUAL_ENV env var (activated venv)
    if let Some(venv) = std::env::var_os("VIRTUAL_ENV") {
        let venv = std::path::PathBuf::from(venv);
        candidates.push(venv.join("bin/python3"));
        candidates.push(venv.join("bin/python"));
    }

    // 2. .venv/ in current directory (common convention, even if not activated)
    let cwd_venv = std::path::PathBuf::from(".venv");
    if cwd_venv.is_dir() {
        candidates.push(cwd_venv.join("bin/python3"));
        candidates.push(cwd_venv.join("bin/python"));
    }

    // 3. System Python from PATH
    candidates.push("python3".into());
    candidates.push("python".into());

    candidates
}

/// Prepend a path to `PYTHONPATH`.
fn prepend_pythonpath(path: &str) {
    let existing = std::env::var("PYTHONPATH").unwrap_or_default();
    let new_path = if existing.is_empty() {
        path.to_string()
    } else {
        format!("{}:{}", path, existing)
    };
    log::debug!("Python plugin: setting PYTHONPATH={}", new_path);
    std::env::set_var("PYTHONPATH", &new_path);
}

/// Configure the Python interpreter's environment before initialization.
///
/// When embedding Python in a non-Python host binary, `Py_Initialize()` needs
/// `PYTHONHOME` to locate the standard library. We detect it from the best
/// available Python executable whose version matches the PyO3 build-time version.
///
/// If the matched interpreter lives in a venv, its `site-packages` is added
/// to `PYTHONPATH` so that installed packages (like `pymol_rs`) are importable.
fn configure_python_home() {
    // Skip if PYTHONHOME is already set by the user
    if std::env::var_os("PYTHONHOME").is_some() {
        ConfigStatus::Ready.store();
        return;
    }

    let required_version = PYO3_PYTHON_VERSION;
    log::info!(
        "Python plugin: compiled against Python {}",
        required_version
    );

    let version_check = required_version != "unknown";

    // Query each candidate for the paths we need:
    //   line 0: sys.version_info major.minor
    //   line 1: sys.base_prefix      (real stdlib root, for PYTHONHOME)
    //   line 2: sys.base_exec_prefix
    //   line 3: sys.prefix            (may differ in a venv)
    //   line 4: first site-packages   (for PYTHONPATH when in a venv)
    let script = "\
import sys, site
print(f'{sys.version_info.major}.{sys.version_info.minor}')
print(sys.base_prefix)
print(sys.base_exec_prefix)
print(sys.prefix)
sp = site.getsitepackages()
print(sp[0] if sp else '')
";

    for cmd in python_candidates() {
        let output = match std::process::Command::new(&cmd)
            .args(["-c", script])
            .output()
        {
            Ok(o) if o.status.success() => o,
            _ => continue,
        };

        let stdout = String::from_utf8_lossy(&output.stdout);
        let mut lines = stdout.lines();
        let version = match lines.next() { Some(s) => s, None => continue };
        let base_prefix = match lines.next() { Some(s) => s, None => continue };
        let base_exec_prefix = match lines.next() { Some(s) => s, None => continue };
        let prefix = match lines.next() { Some(s) => s, None => continue };
        let site_packages = lines.next().unwrap_or("");

        if version_check && version != required_version {
            log::debug!(
                "Python plugin: skipping {:?} (version {} != {})",
                cmd, version, required_version
            );
            continue;
        }

        log::debug!("Python plugin: using interpreter {:?} ({})", cmd, version);

        // PYTHONHOME must point to the real stdlib (base_prefix)
        let home = if base_prefix == base_exec_prefix {
            base_prefix.to_string()
        } else {
            format!("{}:{}", base_prefix, base_exec_prefix)
        };
        log::debug!("Python plugin: setting PYTHONHOME={}", home);
        std::env::set_var("PYTHONHOME", &home);

        // If this candidate lives in a venv, add its site-packages
        if prefix != base_prefix && !site_packages.is_empty() {
            prepend_pythonpath(site_packages);
        }

        ConfigStatus::Ready.store();
        return;
    }

    log::warn!(
        "Python plugin: no Python {} found; plugin disabled",
        required_version
    );
    ConfigStatus::NoPython.store();
}

// ---------------------------------------------------------------------------
// Public engine API
// ---------------------------------------------------------------------------

/// Embedded Python engine.
///
/// Wraps PyO3 GIL interactions for evaluating expressions and running scripts.
pub struct PythonEngine {
    initialized: bool,
    failed: bool,
}

impl PythonEngine {
    pub fn new() -> Self {
        Self {
            initialized: false,
            failed: false,
        }
    }

    /// Ensure the Python interpreter is ready and `pymol_rs` is importable.
    pub(crate) fn ensure_init(&mut self) -> Result<(), String> {
        if self.initialized {
            return Ok(());
        }
        if self.failed {
            return Err("Python interpreter failed to initialize".to_string());
        }

        // Configure PYTHONHOME once (before any Python::attach call)
        PYTHON_CONFIGURED.call_once(configure_python_home);

        // If no matching Python was found, fail gracefully instead of
        // letting Py_Initialize call abort() which kills the whole process.
        if matches!(ConfigStatus::load(), ConfigStatus::NoPython) {
            self.failed = true;
            return Err(format!(
                "Python {} not found on this system; Python plugin disabled",
                PYO3_PYTHON_VERSION
            ));
        }

        // Catch panics from Python initialization (Py_Initialize may abort,
        // but if it raises a catchable error we handle it gracefully).
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            Python::attach(|py| {
                match py.import("pymol_rs") {
                    Ok(_) => log::info!("Python plugin: pymol_rs package found"),
                    Err(_) => log::warn!(
                        "Python plugin: pymol_rs package not importable; \
                         cmd.* won't work in embedded mode"
                    ),
                }
            });
        }));

        match result {
            Ok(()) => {
                self.initialized = true;
                Ok(())
            }
            Err(_) => {
                self.failed = true;
                Err("Python interpreter panicked during initialization".to_string())
            }
        }
    }

    /// Evaluate a Python expression or statement string.
    ///
    /// Returns captured stdout/stderr output and any error message.
    pub fn eval(&mut self, code: &str) -> Result<String, String> {
        self.ensure_init()?;

        Python::attach(|py| {
            let io = py.import("io").map_err(|e| e.to_string())?;
            let sys = py.import("sys").map_err(|e| e.to_string())?;

            let string_io = io
                .getattr("StringIO")
                .map_err(|e| e.to_string())?
                .call0()
                .map_err(|e| e.to_string())?;

            let old_stdout = sys.getattr("stdout").map_err(|e| e.to_string())?;
            let old_stderr = sys.getattr("stderr").map_err(|e| e.to_string())?;

            sys.setattr("stdout", &string_io).map_err(|e| e.to_string())?;
            sys.setattr("stderr", &string_io).map_err(|e| e.to_string())?;

            let globals = py
                .import("__main__")
                .map_err(|e| e.to_string())?
                .dict();

            let c_code = CString::new(code).map_err(|e| e.to_string())?;
            let result = py.run(&c_code, Some(&globals), None);

            // Restore stdout/stderr before inspecting the result
            let _ = sys.setattr("stdout", old_stdout);
            let _ = sys.setattr("stderr", old_stderr);

            let output: String = string_io
                .call_method0("getvalue")
                .map_err(|e: PyErr| e.to_string())?
                .extract()
                .map_err(|e: PyErr| e.to_string())?;

            match result {
                Ok(()) => Ok(output),
                Err(e) => {
                    let err_msg = e.to_string();
                    if output.is_empty() {
                        Err(err_msg)
                    } else {
                        Err(format!("{}\n{}", output, err_msg))
                    }
                }
            }
        })
    }

    /// Execute a Python script file.
    pub fn exec_file(&mut self, path: &str) -> Result<String, String> {
        let code = std::fs::read_to_string(path)
            .map_err(|e| format!("failed to read {}: {}", path, e))?;
        self.eval(&code)
    }

    /// Install the plugin backend into `sys._pymolrs_backend`.
    ///
    /// Called once during initialization so that `pymol_rs` auto-detects
    /// embedded mode.
    pub fn set_backend(&mut self, backend: Py<PyAny>) -> Result<(), String> {
        self.ensure_init()?;

        Python::attach(|py| {
            let sys = py.import("sys").map_err(|e| e.to_string())?;
            sys.setattr("_pymolrs_backend", backend.bind(py))
                .map_err(|e| e.to_string())?;
            Ok(())
        })
    }
}
