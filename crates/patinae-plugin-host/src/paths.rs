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
    let deps_dir = deps_path.parent().unwrap_or_else(|| Path::new("."));
    let updates = deps_search_path_updates(
        &content,
        deps_dir,
        &current_path,
        |path| path.is_dir(),
        |line, args| {
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
        },
    );

    for (key, value) in &updates.env_updates {
        std::env::set_var(key, value);
        log::debug!(
            "Plugin deps ({:?}): set {}={}",
            deps_path.file_name().unwrap_or_default(),
            key,
            value
        );
    }

    if !updates.path_prepend.is_empty() {
        let prepend = updates.path_prepend.join(";");
        std::env::set_var("PATH", format!("{prepend};{current_path}"));
        log::debug!(
            "Plugin deps ({:?}): added to PATH: {}",
            deps_path.file_name().unwrap_or_default(),
            prepend
        );
    }
}

#[cfg(any(target_os = "windows", test))]
#[derive(Debug, Default, PartialEq, Eq)]
struct DepsSearchPathUpdates {
    path_prepend: Vec<String>,
    env_updates: Vec<(String, String)>,
}

#[cfg(any(target_os = "windows", test))]
fn deps_search_path_updates(
    content: &str,
    deps_dir: &Path,
    current_path: &str,
    mut path_exists: impl FnMut(&Path) -> bool,
    mut run_command: impl FnMut(&str, &[String]) -> Option<String>,
) -> DepsSearchPathUpdates {
    let mut updates = DepsSearchPathUpdates::default();

    for line in content.lines().map(str::trim) {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let args = parse_shell_words(line);
        if args.is_empty() {
            continue;
        }

        match args.first().map(String::as_str) {
            Some("path") => {
                let Some(raw_path) = args.get(1).filter(|_| args.len() == 2) else {
                    log::debug!("Plugin deps: invalid path directive: {}", line);
                    continue;
                };
                let path = resolve_deps_path(deps_dir, raw_path);
                if path_exists(&path) {
                    add_unique_path(&mut updates.path_prepend, current_path, path_string(&path));
                } else {
                    log::debug!("Plugin deps: path does not exist: {}", path.display());
                }
            }
            Some("env") => {
                let (Some(key), Some(raw_path)) =
                    (args.get(1), args.get(2).filter(|_| args.len() == 3))
                else {
                    log::debug!("Plugin deps: invalid env directive: {}", line);
                    continue;
                };
                let path = resolve_deps_path(deps_dir, raw_path);
                if path_exists(&path) {
                    updates.env_updates.push((key.clone(), path_string(&path)));
                } else {
                    log::debug!("Plugin deps: env path does not exist: {}", path.display());
                }
            }
            _ => {
                let Some(dir) = run_command(line, &args) else {
                    continue;
                };
                add_unique_path(
                    &mut updates.path_prepend,
                    current_path,
                    dir.trim().to_string(),
                );
            }
        }
    }

    updates
}

#[cfg(any(target_os = "windows", test))]
fn add_unique_path(added: &mut Vec<String>, current_path: &str, dir: String) {
    if !dir.is_empty()
        && !current_path
            .split(';')
            .any(|p| p.eq_ignore_ascii_case(&dir))
        && !added.iter().any(|p| p.eq_ignore_ascii_case(&dir))
    {
        added.push(dir);
    }
}

#[cfg(any(target_os = "windows", test))]
fn resolve_deps_path(deps_dir: &Path, raw_path: &str) -> PathBuf {
    let path = PathBuf::from(raw_path);
    if path.is_absolute() {
        normalize_deps_path(&path)
    } else {
        normalize_deps_path(&deps_dir.join(path))
    }
}

#[cfg(any(target_os = "windows", test))]
fn normalize_deps_path(path: &Path) -> PathBuf {
    let mut normalized = PathBuf::new();
    for component in path.components() {
        match component {
            std::path::Component::CurDir => {}
            std::path::Component::ParentDir => {
                if !normalized.pop() {
                    normalized.push("..");
                }
            }
            _ => normalized.push(component.as_os_str()),
        }
    }
    normalized
}

#[cfg(any(target_os = "windows", test))]
fn path_string(path: &Path) -> String {
    path.to_string_lossy().to_string()
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
        let updates = deps_search_path_updates(
            content,
            Path::new("/bundle/plugins"),
            "C:\\Existing;C:\\Other",
            |_| false,
            |_, args| {
                calls.push(args.to_vec());
                match args.first().map(String::as_str) {
                    Some("python") => Some("C:\\Dep A".to_string()),
                    Some("tool") => Some("C:\\Dep B".to_string()),
                    Some("existing") => Some("c:\\existing".to_string()),
                    _ => None,
                }
            },
        );

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
            updates.path_prepend,
            vec!["C:\\Dep A".to_string(), "C:\\Dep B".to_string()]
        );
        assert!(updates.env_updates.is_empty());
    }

    #[test]
    fn deps_search_path_updates_resolve_relative_path_and_env_directives() {
        let content = "\
path ../python
path ../python-venv/Scripts
env VIRTUAL_ENV ../python-venv
";
        let existing = [
            PathBuf::from("/bundle/python"),
            PathBuf::from("/bundle/python-venv"),
            PathBuf::from("/bundle/python-venv/Scripts"),
        ];

        let updates = deps_search_path_updates(
            content,
            Path::new("/bundle/plugins"),
            "",
            |path| existing.iter().any(|existing_path| existing_path == path),
            |_, _| None,
        );

        assert_eq!(
            updates.path_prepend,
            vec![
                "/bundle/python".to_string(),
                "/bundle/python-venv/Scripts".to_string(),
            ]
        );
        assert_eq!(
            updates.env_updates,
            vec![("VIRTUAL_ENV".to_string(), "/bundle/python-venv".to_string())]
        );
    }

    #[test]
    fn deps_search_path_updates_skip_missing_declarative_paths() {
        let content = "\
path ../missing-python
env VIRTUAL_ENV ../missing-venv
";
        let updates = deps_search_path_updates(
            content,
            Path::new("/bundle/plugins"),
            "",
            |_| false,
            |_, _| None,
        );

        assert!(updates.path_prepend.is_empty());
        assert!(updates.env_updates.is_empty());
    }

    #[test]
    fn deps_search_path_updates_deduplicate_declarative_paths_against_current_path() {
        let content = "\
path ../python
path ../python
";
        let updates = deps_search_path_updates(
            content,
            Path::new("/bundle/plugins"),
            "/bundle/python;/other",
            |path| path == Path::new("/bundle/python"),
            |_, _| None,
        );

        assert!(updates.path_prepend.is_empty());
    }

    #[test]
    fn deps_search_path_updates_skip_failed_or_empty_commands() {
        let content = "missing\nempty\n";
        let updates = deps_search_path_updates(
            content,
            Path::new("/bundle/plugins"),
            "",
            |_| false,
            |_, args| match args[0].as_str() {
                "empty" => Some(String::new()),
                _ => None,
            },
        );
        assert!(updates.path_prepend.is_empty());
        assert!(updates.env_updates.is_empty());
    }
}
