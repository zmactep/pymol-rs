//! Standard paths for PyMOL-RS configuration, plugins, and startup script.
//!
//! Each path can be overridden via an environment variable.  When the variable
//! is unset (or empty), the default under `~/.pymol-rs` is used.
//!
//! | Variable              | Default                    |
//! |-----------------------|----------------------------|
//! | `PYMOL_RS_CONFIG_DIR` | `~/.pymol-rs`              |
//! | `PYMOL_RS_PLUGIN_DIR` | `<config_dir>/plugins`     |
//! | `PYMOL_RS_PYMOLRC`    | `<config_dir>/pymolrc`     |

use std::path::PathBuf;

/// Default directory name under `$HOME`.
#[cfg(not(target_arch = "wasm32"))]
const DEFAULT_CONFIG_DIR_NAME: &str = ".pymol-rs";

/// Returns the application configuration directory.
///
/// Checks `PYMOL_RS_CONFIG_DIR` first, then falls back to `~/.pymol-rs`.
pub fn config_dir() -> PathBuf {
    if let Ok(val) = std::env::var("PYMOL_RS_CONFIG_DIR") {
        if !val.is_empty() {
            return PathBuf::from(val);
        }
    }
    default_config_dir()
}

/// Returns the plugin directory.
///
/// Checks `PYMOL_RS_PLUGIN_DIR` first, then falls back to `<config_dir>/plugins`.
pub fn plugin_dir() -> PathBuf {
    if let Ok(val) = std::env::var("PYMOL_RS_PLUGIN_DIR") {
        if !val.is_empty() {
            return PathBuf::from(val);
        }
    }
    config_dir().join("plugins")
}

/// Returns the path to the startup script.
///
/// Checks `PYMOL_RS_PYMOLRC` first, then falls back to `<config_dir>/pymolrc`.
pub fn pymolrc_path() -> PathBuf {
    if let Ok(val) = std::env::var("PYMOL_RS_PYMOLRC") {
        if !val.is_empty() {
            return PathBuf::from(val);
        }
    }
    config_dir().join("pymolrc")
}

#[cfg(not(target_arch = "wasm32"))]
fn default_config_dir() -> PathBuf {
    dirs::home_dir()
        .unwrap_or_default()
        .join(DEFAULT_CONFIG_DIR_NAME)
}

#[cfg(target_arch = "wasm32")]
fn default_config_dir() -> PathBuf {
    PathBuf::new()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn config_dir_env_override() {
        std::env::set_var("PYMOL_RS_CONFIG_DIR", "/tmp/test-pymol");
        assert_eq!(config_dir(), PathBuf::from("/tmp/test-pymol"));
        std::env::remove_var("PYMOL_RS_CONFIG_DIR");
    }

    #[test]
    fn plugin_dir_derives_from_config_dir() {
        std::env::remove_var("PYMOL_RS_PLUGIN_DIR");
        std::env::set_var("PYMOL_RS_CONFIG_DIR", "/tmp/test-pymol");
        assert_eq!(plugin_dir(), PathBuf::from("/tmp/test-pymol/plugins"));
        std::env::remove_var("PYMOL_RS_CONFIG_DIR");
    }

    #[test]
    fn plugin_dir_env_override_independent() {
        std::env::set_var("PYMOL_RS_CONFIG_DIR", "/tmp/test-pymol");
        std::env::set_var("PYMOL_RS_PLUGIN_DIR", "/opt/plugins");
        assert_eq!(plugin_dir(), PathBuf::from("/opt/plugins"));
        std::env::remove_var("PYMOL_RS_CONFIG_DIR");
        std::env::remove_var("PYMOL_RS_PLUGIN_DIR");
    }

    #[test]
    fn pymolrc_path_derives_from_config_dir() {
        std::env::remove_var("PYMOL_RS_PYMOLRC");
        std::env::set_var("PYMOL_RS_CONFIG_DIR", "/tmp/test-pymol");
        assert_eq!(pymolrc_path(), PathBuf::from("/tmp/test-pymol/pymolrc"));
        std::env::remove_var("PYMOL_RS_CONFIG_DIR");
    }
}
