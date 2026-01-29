//! Python Script Execution
//!
//! Provides functionality to execute Python scripts from within pymol-rs.

use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::ffi::CString;
use std::path::Path;

/// Execute a Python script file
///
/// # Arguments
/// * `py` - Python interpreter handle
/// * `filename` - Path to the Python script
/// * `namespace` - Execution namespace:
///   - "global": Execute in pymol_rs module namespace (default)
///   - "local": Execute in a fresh namespace with `cmd` pre-imported
///   - "module": Import as a Python module
///
/// # Errors
/// Returns an error if the file cannot be read or execution fails.
pub fn run_python_script(
    py: Python<'_>,
    filename: &str,
    namespace: &str,
) -> PyResult<()> {
    let path = Path::new(filename);
    
    // Expand ~ in path
    let expanded_path = if filename.starts_with('~') {
        if let Some(home) = std::env::var_os("HOME") {
            let home_str = home.to_string_lossy();
            Path::new(&filename.replacen('~', &home_str, 1)).to_path_buf()
        } else {
            path.to_path_buf()
        }
    } else {
        path.to_path_buf()
    };
    
    // Read the script file
    let code = std::fs::read_to_string(&expanded_path).map_err(|e| {
        pyo3::exceptions::PyIOError::new_err(format!(
            "Failed to read script '{}': {}",
            filename, e
        ))
    })?;
    
    // Convert to CString for py.run
    let code_cstr = CString::new(code.as_str()).map_err(|e| {
        pyo3::exceptions::PyValueError::new_err(format!(
            "Script contains null bytes: {}",
            e
        ))
    })?;
    
    match namespace {
        "global" | "" => {
            // Execute in pymol_rs module namespace
            let pymol_rs = py.import("pymol_rs")?;
            let globals = pymol_rs.dict();
            py.run(&code_cstr, Some(&globals), None)?;
        }
        "local" => {
            // Execute in a fresh local namespace
            let locals = PyDict::new(py);
            
            // Pre-import cmd for convenience
            let pymol_rs = py.import("pymol_rs")?;
            let cmd = pymol_rs.getattr("cmd")?;
            locals.set_item("cmd", cmd)?;
            
            // Also import common modules
            locals.set_item("pymol_rs", pymol_rs)?;
            
            py.run(&code_cstr, Some(&locals), Some(&locals))?;
        }
        "module" => {
            // Import as a module
            let module_name = expanded_path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("script");
            
            let module_name_cstr = CString::new(module_name).map_err(|e| {
                pyo3::exceptions::PyValueError::new_err(format!(
                    "Module name contains null bytes: {}",
                    e
                ))
            })?;
            let filename_cstr = CString::new(filename).map_err(|e| {
                pyo3::exceptions::PyValueError::new_err(format!(
                    "Filename contains null bytes: {}",
                    e
                ))
            })?;
            
            pyo3::types::PyModule::from_code(py, &code_cstr, &filename_cstr, &module_name_cstr)?;
        }
        _ => {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "Invalid namespace: '{}'. Use 'global', 'local', or 'module'.",
                namespace
            )));
        }
    }
    
    Ok(())
}

#[cfg(test)]
mod tests {
    // Tests would require a Python interpreter
}
