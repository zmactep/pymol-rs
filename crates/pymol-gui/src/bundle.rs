//! macOS `.app` bundle path resolution.
//!
//! When the binary runs inside `PyMOL-RS.app/Contents/MacOS/pymol-rs`,
//! this module locates bundled resources (plugins, Python, virtualenv)
//! relative to the executable.

use std::path::PathBuf;

/// Detect whether we are running inside a macOS `.app` bundle.
///
/// Returns the `Contents/` directory if the executable lives at
/// `<name>.app/Contents/MacOS/<binary>` and `Info.plist` exists.
pub fn bundle_contents_dir() -> Option<PathBuf> {
    let exe = std::env::current_exe().ok()?;
    let macos_dir = exe.parent()?;
    if macos_dir.file_name()?.to_str()? != "MacOS" {
        return None;
    }
    let contents_dir = macos_dir.parent()?;
    if contents_dir.file_name()?.to_str()? != "Contents" {
        return None;
    }
    if contents_dir.join("Info.plist").is_file() {
        Some(contents_dir.to_path_buf())
    } else {
        None
    }
}
