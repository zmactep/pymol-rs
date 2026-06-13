//! Python Interpreter Engine
//!
//! Manages the embedded Python interpreter, providing `eval` and `exec_file` operations.
//! Uses PyO3 with `auto-initialize` — the interpreter starts on first GIL acquisition.
//!
//! Before the first GIL acquisition, we configure `PYTHONHOME` so that the
//! embedded interpreter can find its standard library (encodings, etc.).
//! PyO3 interpreter initialization and environment updates are process-global;
//! this module keeps the planning logic deterministic and applies the resulting
//! changes only at the initialization boundary.

use std::ffi::{CString, OsString};
use std::path::{Path, PathBuf};
use std::process::Command;
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
const PYO3_PYTHON_VERSION: &str = env!("PATINAE_PYO3_PYTHON_VERSION");

const PYTHON_PROBE_SCRIPT: &str = "\
import sys, site
print(f'{sys.version_info.major}.{sys.version_info.minor}')
print(sys.base_prefix)
print(sys.base_exec_prefix)
print(sys.prefix)
sp = [p for p in site.getsitepackages() if p.endswith('site-packages')]
print(sp[0] if sp else '')
";

/// Result of the one-time Python environment configuration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::{BTreeMap, BTreeSet};

    #[derive(Default)]
    struct FakePythonEnvironment {
        current_exe: Option<PathBuf>,
        vars: BTreeMap<String, OsString>,
        files: BTreeSet<PathBuf>,
        dirs: BTreeSet<PathBuf>,
        read_dirs: BTreeMap<PathBuf, Vec<PathBuf>>,
        probes: BTreeMap<PathBuf, PythonProbe>,
    }

    impl PythonEnvironment for FakePythonEnvironment {
        fn current_exe(&self) -> Option<PathBuf> {
            self.current_exe.clone()
        }

        fn var_os(&self, key: &str) -> Option<OsString> {
            self.vars.get(key).cloned()
        }

        fn is_file(&self, path: &Path) -> bool {
            self.files.contains(path)
        }

        fn is_dir(&self, path: &Path) -> bool {
            self.dirs.contains(path)
        }

        fn read_dir_paths(&self, path: &Path) -> Vec<PathBuf> {
            self.read_dirs.get(path).cloned().unwrap_or_default()
        }

        fn probe_python(&self, cmd: &Path) -> Option<PythonProbe> {
            self.probes.get(cmd).cloned()
        }
    }

    fn probe(
        version: &str,
        base_prefix: &str,
        base_exec_prefix: &str,
        prefix: &str,
        site_packages: Option<&str>,
    ) -> PythonProbe {
        PythonProbe {
            version: version.to_string(),
            base_prefix: PathBuf::from(base_prefix),
            base_exec_prefix: PathBuf::from(base_exec_prefix),
            prefix: PathBuf::from(prefix),
            site_packages: site_packages.map(PathBuf::from),
        }
    }

    fn venv_python(venv: &str) -> PathBuf {
        let venv = PathBuf::from(venv);
        if cfg!(target_os = "windows") {
            venv.join("Scripts").join("python.exe")
        } else {
            venv.join("bin/python3")
        }
    }

    fn expected_windows_path_prepend(path: &str) -> Option<PathBuf> {
        if cfg!(target_os = "windows") {
            Some(PathBuf::from(path))
        } else {
            None
        }
    }

    #[test]
    fn existing_pythonhome_returns_ready_without_updates() {
        let mut env = FakePythonEnvironment::default();
        env.vars
            .insert("PYTHONHOME".to_string(), OsString::from("/custom/python"));

        let plan = plan_python_environment("3.11", &env);

        assert_eq!(plan, PythonEnvPlan::ready_without_updates());
    }

    #[test]
    fn missing_python_returns_no_python() {
        let env = FakePythonEnvironment::default();

        let plan = plan_python_environment("3.11", &env);

        assert_eq!(plan.status, ConfigStatus::NoPython);
        assert!(!plan.has_process_updates());
    }

    #[test]
    fn wrong_python_version_is_skipped() {
        let mut env = FakePythonEnvironment::default();
        env.probes.insert(
            PathBuf::from("python3"),
            probe("3.10", "/python310", "/python310", "/python310", None),
        );

        let plan = plan_python_environment("3.11", &env);

        assert_eq!(plan.status, ConfigStatus::NoPython);
    }

    #[test]
    fn venv_site_packages_are_added_to_pythonpath() {
        let mut env = FakePythonEnvironment::default();
        env.vars
            .insert("VIRTUAL_ENV".to_string(), OsString::from("/venv"));
        env.probes.insert(
            venv_python("/venv"),
            probe(
                "3.11",
                "/python311",
                "/python311",
                "/venv",
                Some("/venv/lib/python3.11/site-packages"),
            ),
        );

        let plan = plan_python_environment("3.11", &env);

        assert_eq!(plan.status, ConfigStatus::Ready);
        assert_eq!(plan.python_home.as_deref(), Some("/python311"));
        assert_eq!(
            plan.pythonpath_prepend,
            vec![PathBuf::from("/venv/lib/python3.11/site-packages")]
        );
        assert_eq!(
            plan.windows_path_prepend,
            expected_windows_path_prepend("/python311")
        );
    }

    #[test]
    fn pythonhome_uses_base_and_exec_prefix() {
        let mut env = FakePythonEnvironment::default();
        env.probes.insert(
            PathBuf::from("python3"),
            probe("3.11", "/python311", "/python311-exec", "/python311", None),
        );

        let plan = plan_python_environment("3.11", &env);
        let sep = if cfg!(target_os = "windows") {
            ';'
        } else {
            ':'
        };

        assert_eq!(
            plan.python_home,
            Some(format!("/python311{sep}/python311-exec"))
        );
    }

    #[test]
    fn bundled_app_paths_have_highest_priority() {
        let mut env = FakePythonEnvironment {
            current_exe: Some(PathBuf::from(
                "/Applications/Patinae.app/Contents/MacOS/patinae",
            )),
            ..FakePythonEnvironment::default()
        };
        let bundled_python =
            PathBuf::from("/Applications/Patinae.app/Contents/Resources/python/bin/python3");
        let bundled_lib =
            PathBuf::from("/Applications/Patinae.app/Contents/Resources/python-venv/lib");
        let bundled_python_lib = bundled_lib.join("python3.11");
        let bundled_site_packages = bundled_python_lib.join("site-packages");
        env.files.insert(bundled_python.clone());
        env.dirs.insert(bundled_lib.clone());
        env.dirs.insert(bundled_site_packages.clone());
        env.read_dirs.insert(bundled_lib, vec![bundled_python_lib]);
        env.probes.insert(
            bundled_python,
            probe(
                "3.11",
                "/bundle/python",
                "/bundle/python",
                "/bundle/python",
                None,
            ),
        );
        env.probes.insert(
            PathBuf::from("python3"),
            probe(
                "3.11",
                "/system/python",
                "/system/python",
                "/system/python",
                None,
            ),
        );

        let plan = plan_python_environment("3.11", &env);

        assert_eq!(plan.python_home.as_deref(), Some("/bundle/python"));
        assert_eq!(plan.pythonpath_prepend, vec![bundled_site_packages]);
    }

    #[test]
    fn probe_stdout_is_parsed() {
        let stdout = b"3.11\n/base\n/exec\n/venv\n/venv/site-packages\n";

        let parsed = parse_python_probe_stdout(stdout).expect("probe output should parse");

        assert_eq!(
            parsed,
            probe(
                "3.11",
                "/base",
                "/exec",
                "/venv",
                Some("/venv/site-packages")
            )
        );
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct PythonProbe {
    version: String,
    base_prefix: PathBuf,
    base_exec_prefix: PathBuf,
    prefix: PathBuf,
    site_packages: Option<PathBuf>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct PythonEnvPlan {
    status: ConfigStatus,
    python_home: Option<String>,
    pythonpath_prepend: Vec<PathBuf>,
    windows_path_prepend: Option<PathBuf>,
}

impl PythonEnvPlan {
    fn ready_without_updates() -> Self {
        Self {
            status: ConfigStatus::Ready,
            python_home: None,
            pythonpath_prepend: Vec::new(),
            windows_path_prepend: None,
        }
    }

    fn no_python() -> Self {
        Self {
            status: ConfigStatus::NoPython,
            python_home: None,
            pythonpath_prepend: Vec::new(),
            windows_path_prepend: None,
        }
    }

    fn has_process_updates(&self) -> bool {
        self.python_home.is_some()
            || !self.pythonpath_prepend.is_empty()
            || self.windows_path_prepend.is_some()
    }
}

trait PythonEnvironment {
    fn current_exe(&self) -> Option<PathBuf>;
    fn var_os(&self, key: &str) -> Option<OsString>;
    fn is_file(&self, path: &Path) -> bool;
    fn is_dir(&self, path: &Path) -> bool;
    fn read_dir_paths(&self, path: &Path) -> Vec<PathBuf>;
    fn probe_python(&self, cmd: &Path) -> Option<PythonProbe>;
}

struct ProcessPythonEnvironment;

impl PythonEnvironment for ProcessPythonEnvironment {
    fn current_exe(&self) -> Option<PathBuf> {
        std::env::current_exe().ok()
    }

    fn var_os(&self, key: &str) -> Option<OsString> {
        std::env::var_os(key)
    }

    fn is_file(&self, path: &Path) -> bool {
        path.is_file()
    }

    fn is_dir(&self, path: &Path) -> bool {
        path.is_dir()
    }

    fn read_dir_paths(&self, path: &Path) -> Vec<PathBuf> {
        let Ok(entries) = std::fs::read_dir(path) else {
            return Vec::new();
        };
        entries.flatten().map(|entry| entry.path()).collect()
    }

    fn probe_python(&self, cmd: &Path) -> Option<PythonProbe> {
        let output = Command::new(cmd)
            .args(["-c", PYTHON_PROBE_SCRIPT])
            .output()
            .ok()?;
        if !output.status.success() {
            return None;
        }
        parse_python_probe_stdout(&output.stdout)
    }
}

fn parse_python_probe_stdout(stdout: &[u8]) -> Option<PythonProbe> {
    let stdout = String::from_utf8_lossy(stdout);
    let mut lines = stdout.lines();
    let version = lines.next()?.to_string();
    let base_prefix = PathBuf::from(lines.next()?);
    let base_exec_prefix = PathBuf::from(lines.next()?);
    let prefix = PathBuf::from(lines.next()?);
    let site_packages = lines
        .next()
        .filter(|path| !path.is_empty())
        .map(PathBuf::from);

    Some(PythonProbe {
        version,
        base_prefix,
        base_exec_prefix,
        prefix,
        site_packages,
    })
}

/// Check if running inside a macOS `.app` bundle and return the bundled Python path.
fn bundled_python_path(env: &impl PythonEnvironment) -> Option<PathBuf> {
    let exe = env.current_exe()?;
    let macos_dir = exe.parent()?;
    if macos_dir.file_name()?.to_str()? != "MacOS" {
        return None;
    }
    let contents = macos_dir.parent()?;
    let python = contents.join("Resources/python/bin/python3");
    if env.is_file(&python) {
        Some(python)
    } else {
        None
    }
}

/// Return the bundled venv's `site-packages` path if running inside a `.app` bundle.
fn bundled_venv_site_packages(env: &impl PythonEnvironment) -> Option<PathBuf> {
    let exe = env.current_exe()?;
    let macos_dir = exe.parent()?;
    if macos_dir.file_name()?.to_str()? != "MacOS" {
        return None;
    }
    let contents = macos_dir.parent()?;
    let venv_lib = contents.join("Resources/python-venv/lib");
    if !env.is_dir(&venv_lib) {
        return None;
    }
    // Find python3.X/site-packages inside the venv lib/
    for entry in env.read_dir_paths(&venv_lib) {
        let sp = entry.join("site-packages");
        if env.is_dir(&sp) {
            return Some(sp);
        }
    }
    None
}

/// Build a list of Python executables to probe, preferring venv if available.
fn python_candidates(env: &impl PythonEnvironment) -> Vec<PathBuf> {
    let mut candidates = Vec::new();

    // 0. Bundled Python inside .app bundle (highest priority)
    if let Some(bundled) = bundled_python_path(env) {
        candidates.push(bundled);
    }

    // 1. VIRTUAL_ENV env var (activated venv)
    if let Some(venv) = env.var_os("VIRTUAL_ENV") {
        let venv = PathBuf::from(venv);
        if cfg!(target_os = "windows") {
            candidates.push(venv.join("Scripts").join("python.exe"));
        } else {
            candidates.push(venv.join("bin/python3"));
            candidates.push(venv.join("bin/python"));
        }
    }

    // 2. .venv/ in current directory (common convention, even if not activated)
    let cwd_venv = PathBuf::from(".venv");
    if env.is_dir(&cwd_venv) {
        if cfg!(target_os = "windows") {
            candidates.push(cwd_venv.join("Scripts").join("python.exe"));
        } else {
            candidates.push(cwd_venv.join("bin/python3"));
            candidates.push(cwd_venv.join("bin/python"));
        }
    }

    // 3. System Python from PATH
    if cfg!(target_os = "windows") {
        candidates.push("python.exe".into());
        candidates.push("python3.exe".into());
    } else {
        candidates.push("python3".into());
        candidates.push("python".into());
    }

    candidates
}

/// Prepend a path to `PYTHONPATH`.
fn prepend_pythonpath(path: &Path) {
    let existing = std::env::var_os("PYTHONPATH").unwrap_or_default();
    let mut paths = vec![path.to_path_buf()];
    paths.extend(std::env::split_paths(&existing));
    if let Ok(new_path) = std::env::join_paths(&paths) {
        log::debug!(
            "Python plugin: setting PYTHONPATH={}",
            new_path.to_string_lossy()
        );
        std::env::set_var("PYTHONPATH", &new_path);
    }
}

fn plan_python_environment(required_version: &str, env: &impl PythonEnvironment) -> PythonEnvPlan {
    // Skip if PYTHONHOME is already set by the user.
    if env.var_os("PYTHONHOME").is_some() {
        return PythonEnvPlan::ready_without_updates();
    }

    let version_check = required_version != "unknown";
    let bundled_site_packages = bundled_venv_site_packages(env);

    for cmd in python_candidates(env) {
        let Some(probe) = env.probe_python(&cmd) else {
            continue;
        };

        if version_check && probe.version != required_version {
            log::debug!(
                "Python plugin: skipping {:?} (version {} != {})",
                cmd,
                probe.version,
                required_version
            );
            continue;
        }

        log::debug!(
            "Python plugin: using interpreter {:?} ({})",
            cmd,
            probe.version
        );

        let mut pythonpath_prepend = Vec::new();
        if probe.prefix != probe.base_prefix {
            if let Some(site_packages) = &probe.site_packages {
                pythonpath_prepend.push(site_packages.clone());
            }
        }
        if let Some(site_packages) = bundled_site_packages {
            pythonpath_prepend.push(site_packages);
        }

        let windows_path_prepend = if cfg!(target_os = "windows") {
            Some(probe.base_prefix.clone())
        } else {
            None
        };

        return PythonEnvPlan {
            status: ConfigStatus::Ready,
            python_home: Some(python_home_value(&probe)),
            pythonpath_prepend,
            windows_path_prepend,
        };
    }

    PythonEnvPlan::no_python()
}

fn python_home_value(probe: &PythonProbe) -> String {
    if probe.base_prefix == probe.base_exec_prefix {
        probe.base_prefix.to_string_lossy().into_owned()
    } else {
        // Python uses ';' as PYTHONHOME separator on Windows, ':' on Unix.
        let sep = if cfg!(target_os = "windows") {
            ';'
        } else {
            ':'
        };
        format!(
            "{}{}{}",
            probe.base_prefix.to_string_lossy(),
            sep,
            probe.base_exec_prefix.to_string_lossy()
        )
    }
}

fn apply_python_environment_plan(plan: &PythonEnvPlan) {
    if let Some(home) = &plan.python_home {
        log::debug!("Python plugin: setting PYTHONHOME={}", home);
        std::env::set_var("PYTHONHOME", home);
    }

    #[cfg(target_os = "windows")]
    if let Some(base_prefix) = &plan.windows_path_prepend {
        let current_path = std::env::var("PATH").unwrap_or_default();
        std::env::set_var(
            "PATH",
            format!("{};{}", base_prefix.to_string_lossy(), current_path),
        );
        log::debug!("Python plugin: added {} to PATH", base_prefix.display());
    }

    for path in &plan.pythonpath_prepend {
        prepend_pythonpath(path);
    }
}

/// Configure the Python interpreter's environment before initialization.
///
/// When embedding Python in a non-Python host binary, `Py_Initialize()` needs
/// `PYTHONHOME` to locate the standard library. We detect it from the best
/// available Python executable whose version matches the PyO3 build-time version.
///
/// If the matched interpreter lives in a venv, its `site-packages` is added
/// to `PYTHONPATH` so that installed packages (like `patinae`) are importable.
fn configure_python_home() {
    let required_version = PYO3_PYTHON_VERSION;
    let env = ProcessPythonEnvironment;
    let plan = plan_python_environment(required_version, &env);
    if matches!(plan.status, ConfigStatus::Ready) && !plan.has_process_updates() {
        ConfigStatus::Ready.store();
        return;
    }

    log::info!(
        "Python plugin: compiled against Python {}",
        required_version
    );

    match plan.status {
        ConfigStatus::Ready => {
            apply_python_environment_plan(&plan);
            ConfigStatus::Ready.store();
        }
        ConfigStatus::NoPython => {
            log::warn!(
                "Python plugin: no Python {} found; plugin disabled",
                required_version
            );
            ConfigStatus::NoPython.store();
        }
        ConfigStatus::Pending => ConfigStatus::Pending.store(),
    }
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

    /// Ensure the Python interpreter is ready and `patinae` is importable.
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
            Python::attach(|py| match py.import("patinae") {
                Ok(_) => log::info!("Python plugin: patinae package found"),
                Err(_) => log::warn!(
                    "Python plugin: patinae package not importable; \
                         cmd.* won't work in embedded mode"
                ),
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

            sys.setattr("stdout", &string_io)
                .map_err(|e| e.to_string())?;
            sys.setattr("stderr", &string_io)
                .map_err(|e| e.to_string())?;

            let globals = py.import("__main__").map_err(|e| e.to_string())?.dict();

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

    /// Evaluate Python code with custom stdout/stderr writers.
    ///
    /// Unlike [`eval`], this does not capture output into a string — the
    /// provided writer objects receive `write()` calls in real time.
    /// Returns `Ok(())` on success or `Err(error_message)` on exception.
    pub fn eval_with_writers(
        &mut self,
        code: &str,
        stdout: &Bound<'_, PyAny>,
        stderr: &Bound<'_, PyAny>,
    ) -> Result<(), String> {
        self.ensure_init()?;

        Python::attach(|py| {
            let sys = py.import("sys").map_err(|e| e.to_string())?;

            let old_stdout = sys.getattr("stdout").map_err(|e| e.to_string())?;
            let old_stderr = sys.getattr("stderr").map_err(|e| e.to_string())?;

            sys.setattr("stdout", stdout).map_err(|e| e.to_string())?;
            sys.setattr("stderr", stderr).map_err(|e| e.to_string())?;

            let globals = py.import("__main__").map_err(|e| e.to_string())?.dict();

            let c_code = CString::new(code).map_err(|e| e.to_string())?;
            let result = py.run(&c_code, Some(&globals), None);

            // Restore stdout/stderr before inspecting the result
            let _ = sys.setattr("stdout", old_stdout);
            let _ = sys.setattr("stderr", old_stderr);

            result.map_err(|e| e.to_string())
        })
    }

    /// Install the plugin backend into `sys._patinae_backend`.
    ///
    /// Called once during initialization so that `patinae` auto-detects
    /// embedded mode.
    pub fn set_backend(&mut self, backend: Py<PyAny>) -> Result<(), String> {
        self.ensure_init()?;

        Python::attach(|py| {
            let sys = py.import("sys").map_err(|e| e.to_string())?;
            sys.setattr("_patinae_backend", backend.bind(py))
                .map_err(|e| e.to_string())?;
            Ok(())
        })
    }
}
