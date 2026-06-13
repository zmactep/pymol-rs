use std::path::{Path, PathBuf};

use patinae_settings::paths::PathResolver;

/// Discovers plugin directories from settings and executable layout.
#[derive(Debug, Clone)]
pub struct PluginDiscovery {
    paths: PathResolver,
    executable_path: Option<PathBuf>,
}

impl PluginDiscovery {
    /// Creates plugin discovery from deterministic path inputs.
    pub fn new(paths: PathResolver) -> Self {
        Self {
            paths,
            executable_path: None,
        }
    }

    /// Creates plugin discovery from the current process.
    pub fn from_process_env() -> Self {
        Self::new(PathResolver::from_process_env())
            .with_optional_executable_path(std::env::current_exe().ok())
    }

    /// Sets the executable path used for adjacent plugin directories.
    pub fn with_executable_path(mut self, executable_path: impl Into<PathBuf>) -> Self {
        self.executable_path = Some(executable_path.into());
        self
    }

    /// Sets or clears the executable path.
    pub fn with_optional_executable_path(mut self, executable_path: Option<PathBuf>) -> Self {
        self.executable_path = executable_path;
        self
    }

    /// Returns the standard plugin search directories.
    pub fn standard_plugin_dirs(&self) -> Vec<PathBuf> {
        let mut dirs = Vec::new();
        push_unique(&mut dirs, self.paths.plugin_dir());

        if let Some(exe) = &self.executable_path {
            push_executable_plugin_dirs(&mut dirs, exe, cfg!(target_os = "macos"));
        }

        dirs
    }
}

pub fn is_plugin_library_path(path: &Path) -> bool {
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("");
    match std::env::consts::OS {
        "macos" => ext == "dylib",
        "linux" => ext == "so",
        "windows" => ext == "dll",
        _ => false,
    }
}

pub fn standard_plugin_dirs() -> Vec<PathBuf> {
    PluginDiscovery::from_process_env().standard_plugin_dirs()
}

fn push_unique(dirs: &mut Vec<PathBuf>, dir: PathBuf) {
    if !dirs.iter().any(|p| p == &dir) {
        dirs.push(dir);
    }
}

fn push_executable_plugin_dirs(dirs: &mut Vec<PathBuf>, exe: &Path, include_bundle_plugins: bool) {
    let Some(parent) = exe.parent() else {
        return;
    };
    push_unique(dirs, parent.join("plugins"));

    if include_bundle_plugins && parent.file_name().and_then(|s| s.to_str()) == Some("MacOS") {
        if let Some(contents) = parent.parent() {
            push_unique(dirs, contents.join("PlugIns"));
        }
    }
}

#[cfg(target_os = "windows")]
pub(crate) fn apply_deps_search_paths(plugin_path: &Path) {
    use std::os::windows::process::CommandExt;

    const CREATE_NO_WINDOW: u32 = 0x08000000;

    let deps_path = plugin_path.with_extension("deps");
    let content = match std::fs::read_to_string(&deps_path) {
        Ok(c) => c,
        Err(_) => return,
    };

    let current_path = std::env::var("PATH").unwrap_or_default();
    let added = deps_search_path_updates(&content, &current_path, |line, args| {
        let output = match std::process::Command::new(&args[0])
            .args(&args[1..])
            .creation_flags(CREATE_NO_WINDOW)
            .output()
        {
            Ok(o) if o.status.success() => o,
            Ok(o) => {
                log::debug!("Plugin deps: command exited {}: {}", o.status, line);
                return None;
            }
            Err(e) => {
                log::debug!("Plugin deps: failed to run: {} ({})", line, e);
                return None;
            }
        };

        Some(String::from_utf8_lossy(&output.stdout).trim().to_string())
    });

    if !added.is_empty() {
        let prepend = added.join(";");
        std::env::set_var("PATH", format!("{};{}", prepend, current_path));
        log::debug!(
            "Plugin deps ({:?}): added to PATH: {}",
            deps_path.file_name().unwrap_or_default(),
            prepend
        );
    }
}

#[cfg(any(target_os = "windows", test))]
fn deps_search_path_updates(
    content: &str,
    current_path: &str,
    mut run_command: impl FnMut(&str, &[String]) -> Option<String>,
) -> Vec<String> {
    let mut added = Vec::new();

    for line in content.lines().map(str::trim) {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let args = parse_shell_words(line);
        if args.is_empty() {
            continue;
        }

        let Some(dir) = run_command(line, &args) else {
            continue;
        };
        let dir = dir.trim();
        if !dir.is_empty()
            && !current_path.split(';').any(|p| p.eq_ignore_ascii_case(dir))
            && !added.iter().any(|p: &String| p.eq_ignore_ascii_case(dir))
        {
            added.push(dir.to_string());
        }
    }

    added
}

#[cfg(any(target_os = "windows", test))]
fn parse_shell_words(line: &str) -> Vec<String> {
    let mut args = Vec::new();
    let mut current = String::new();
    let mut in_quotes = false;

    for ch in line.chars() {
        match ch {
            '"' => in_quotes = !in_quotes,
            c if c.is_ascii_whitespace() && !in_quotes => {
                if !current.is_empty() {
                    args.push(std::mem::take(&mut current));
                }
            }
            c => current.push(c),
        }
    }
    if !current.is_empty() {
        args.push(current);
    }
    args
}

#[cfg(test)]
mod tests {
    use super::*;
    use patinae_settings::paths::PathResolverInput;

    fn discovery(config_dir: &str, plugin_dir: Option<&str>, exe: Option<&str>) -> PluginDiscovery {
        let input = PathResolverInput {
            home_dir: Some(PathBuf::from("/home/patinae")),
            config_dir: Some(PathBuf::from(config_dir)),
            plugin_dir: plugin_dir.map(PathBuf::from),
            ..PathResolverInput::default()
        };
        PluginDiscovery::new(PathResolver::new(input))
            .with_optional_executable_path(exe.map(PathBuf::from))
    }

    #[test]
    fn discovery_uses_plugin_dir_override() {
        let dirs = discovery("/cfg", Some("/plugins"), None).standard_plugin_dirs();
        assert_eq!(dirs, vec![PathBuf::from("/plugins")]);
    }

    #[test]
    fn discovery_derives_plugin_dir_from_config_dir() {
        let dirs = discovery("/cfg", None, None).standard_plugin_dirs();
        assert_eq!(dirs, vec![PathBuf::from("/cfg/plugins")]);
    }

    #[test]
    fn discovery_includes_executable_adjacent_plugins() {
        let dirs = discovery("/cfg", Some("/plugins"), Some("/opt/patinae/bin/patinae"))
            .standard_plugin_dirs();
        assert_eq!(
            dirs,
            vec![
                PathBuf::from("/plugins"),
                PathBuf::from("/opt/patinae/bin/plugins"),
            ]
        );
    }

    #[test]
    fn discovery_deduplicates_dirs() {
        let dirs = discovery(
            "/cfg",
            Some("/opt/patinae/bin/plugins"),
            Some("/opt/patinae/bin/patinae"),
        )
        .standard_plugin_dirs();
        assert_eq!(dirs, vec![PathBuf::from("/opt/patinae/bin/plugins")]);
    }

    #[test]
    fn executable_dirs_can_include_macos_bundle_plugins() {
        let mut dirs = Vec::new();
        push_executable_plugin_dirs(
            &mut dirs,
            Path::new("/Applications/Patinae.app/Contents/MacOS/patinae"),
            true,
        );
        assert_eq!(
            dirs,
            vec![
                PathBuf::from("/Applications/Patinae.app/Contents/MacOS/plugins"),
                PathBuf::from("/Applications/Patinae.app/Contents/PlugIns"),
            ]
        );
    }

    #[test]
    fn deps_search_path_updates_parse_commands_and_deduplicate() {
        let content = "\
# comment
python -c \"print('C:\\\\Dep A')\"
tool --dir
tool --dir
existing
";
        let mut calls = Vec::new();
        let updates = deps_search_path_updates(content, "C:\\Existing;C:\\Other", |_, args| {
            calls.push(args.to_vec());
            match args.first().map(String::as_str) {
                Some("python") => Some("C:\\Dep A".to_string()),
                Some("tool") => Some("C:\\Dep B".to_string()),
                Some("existing") => Some("c:\\existing".to_string()),
                _ => None,
            }
        });

        assert_eq!(
            calls,
            vec![
                vec![
                    "python".to_string(),
                    "-c".to_string(),
                    "print('C:\\\\Dep A')".to_string(),
                ],
                vec!["tool".to_string(), "--dir".to_string()],
                vec!["tool".to_string(), "--dir".to_string()],
                vec!["existing".to_string()],
            ]
        );
        assert_eq!(
            updates,
            vec!["C:\\Dep A".to_string(), "C:\\Dep B".to_string()]
        );
    }

    #[test]
    fn deps_search_path_updates_skip_failed_or_empty_commands() {
        let content = "missing\nempty\n";
        let updates = deps_search_path_updates(content, "", |_, args| match args[0].as_str() {
            "empty" => Some(String::new()),
            _ => None,
        });
        assert!(updates.is_empty());
    }
}
