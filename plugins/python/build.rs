//! Capture the Python version that PyO3 links against at compile time.

fn main() {
    // Try multiple env vars that PyO3's build chain might set
    let version = extract_version_from_dep_vars()
        .or_else(probe_pyo3_python)
        .unwrap_or_else(|| "unknown".to_string());

    println!("cargo:rustc-env=PYMOLRS_PYO3_PYTHON_VERSION={version}");

    // On Windows (MSVC), delay-load the Python DLL so the plugin can be loaded
    // even when python3XX.dll is not on PATH. The DLL will be resolved on first
    // use, after the plugin's init code has added Python to the DLL search path.
    let target_os = std::env::var("CARGO_CFG_TARGET_OS").unwrap_or_default();
    let target_env = std::env::var("CARGO_CFG_TARGET_ENV").unwrap_or_default();
    if target_os == "windows" && target_env == "msvc" && version != "unknown" {
        let dll_name = format!("python{}", version.replace('.', ""));
        println!("cargo:rustc-cdylib-link-arg=/DELAYLOAD:{dll_name}.dll");
        println!("cargo:rustc-link-lib=delayimp");
    }
}

/// Check DEP_*_PYO3_CONFIG env vars set by pyo3-ffi's build script.
fn extract_version_from_dep_vars() -> Option<String> {
    // The env var name depends on the dependency chain.
    // Try common patterns.
    for key in ["DEP_PYTHON_PYO3_CONFIG", "DEP_PYTHONXY_PYO3_CONFIG"] {
        if let Ok(config) = std::env::var(key) {
            for line in config.lines() {
                if let Some(v) = line.strip_prefix("version=") {
                    return Some(v.to_string());
                }
            }
        }
    }

    // Scan all env vars for any DEP_*PYO3_CONFIG
    for (key, val) in std::env::vars() {
        if key.contains("PYO3_CONFIG") {
            for line in val.lines() {
                if let Some(v) = line.strip_prefix("version=") {
                    return Some(v.to_string());
                }
            }
        }
    }

    None
}

/// Fallback: run the same Python that PYO3_PYTHON points to (or python3/python).
fn probe_pyo3_python() -> Option<String> {
    // On Windows "python" is the standard name; on Unix "python3".
    let candidates: &[&str] = if cfg!(target_os = "windows") {
        &["python", "python3"]
    } else {
        &["python3", "python"]
    };

    let python = std::env::var("PYO3_PYTHON").ok().or_else(|| {
        candidates.iter().find_map(|name| {
            std::process::Command::new(name)
                .arg("--version")
                .output()
                .ok()
                .filter(|o| o.status.success())
                .map(|_| name.to_string())
        })
    })?;

    let output = std::process::Command::new(&python)
        .args(["-c", "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')"])
        .output()
        .ok()?;

    if output.status.success() {
        Some(String::from_utf8_lossy(&output.stdout).trim().to_string())
    } else {
        None
    }
}
