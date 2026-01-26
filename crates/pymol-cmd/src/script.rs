//! Script engine for executing .pml files
//!
//! Provides higher-level script execution with better error handling
//! and support for script-specific features like `@` file inclusion.

use std::path::{Path, PathBuf};

use pymol_scene::Viewer;

use crate::error::{CmdError, CmdResult};
use crate::executor::CommandExecutor;

/// Script engine for executing .pml files and command batches
pub struct ScriptEngine {
    /// The command executor
    executor: CommandExecutor,
    /// Stack of currently executing scripts (for include tracking)
    script_stack: Vec<PathBuf>,
    /// Maximum include depth (to prevent infinite recursion)
    max_include_depth: usize,
}

impl Default for ScriptEngine {
    fn default() -> Self {
        Self::new()
    }
}

impl ScriptEngine {
    /// Create a new script engine
    pub fn new() -> Self {
        Self {
            executor: CommandExecutor::new(),
            script_stack: Vec::new(),
            max_include_depth: 100,
        }
    }

    /// Create a script engine with a custom executor
    pub fn with_executor(executor: CommandExecutor) -> Self {
        Self {
            executor,
            script_stack: Vec::new(),
            max_include_depth: 100,
        }
    }

    /// Get a reference to the executor
    pub fn executor(&self) -> &CommandExecutor {
        &self.executor
    }

    /// Get a mutable reference to the executor
    pub fn executor_mut(&mut self) -> &mut CommandExecutor {
        &mut self.executor
    }

    /// Set the maximum include depth
    pub fn set_max_include_depth(&mut self, depth: usize) {
        self.max_include_depth = depth;
    }

    /// Run a .pml script file
    pub fn run_pml(&mut self, viewer: &mut Viewer, path: &Path) -> CmdResult {
        // Check include depth
        if self.script_stack.len() >= self.max_include_depth {
            return Err(CmdError::Script {
                line: 0,
                message: format!(
                    "maximum include depth ({}) exceeded",
                    self.max_include_depth
                ),
            });
        }

        // Canonicalize path for comparison
        let canonical = path.canonicalize().unwrap_or_else(|_| path.to_path_buf());

        // Check for circular includes
        if self.script_stack.contains(&canonical) {
            return Err(CmdError::Script {
                line: 0,
                message: format!("circular include detected: {}", path.display()),
            });
        }

        // Read the file
        let content = std::fs::read_to_string(path).map_err(|e| CmdError::Script {
            line: 0,
            message: format!("failed to read script: {}", e),
        })?;

        // Push onto stack
        self.script_stack.push(canonical);

        // Execute the script
        let result = self.run_string_with_base(viewer, &content, path.parent());

        // Pop from stack
        self.script_stack.pop();

        result
    }

    /// Run a script string
    pub fn run_string(&mut self, viewer: &mut Viewer, script: &str) -> CmdResult {
        self.run_string_with_base(viewer, script, None)
    }

    /// Run a script string with a base directory for relative paths
    fn run_string_with_base(
        &mut self,
        viewer: &mut Viewer,
        script: &str,
        base_dir: Option<&Path>,
    ) -> CmdResult {
        let lines: Vec<&str> = script.lines().collect();

        for (line_num, line) in lines.iter().enumerate() {
            let line = line.trim();

            // Skip empty lines and comments
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            // Handle @ include syntax
            if line.starts_with('@') {
                let include_path = line[1..].trim();
                let full_path = if let Some(base) = base_dir {
                    base.join(include_path)
                } else {
                    PathBuf::from(include_path)
                };

                self.run_pml(viewer, &full_path).map_err(|e| CmdError::Script {
                    line: line_num + 1,
                    message: format!("in @{}: {}", include_path, e),
                })?;
                continue;
            }

            // Execute the command
            if let Err(e) = self.executor.do_(viewer, line) {
                return Err(CmdError::Script {
                    line: line_num + 1,
                    message: format!("{}: {}", line, e),
                });
            }
        }

        Ok(())
    }

    /// Execute a single command through the engine
    pub fn do_(&mut self, viewer: &mut Viewer, cmd: &str) -> CmdResult {
        self.executor.do_(viewer, cmd)
    }

    /// Execute multiple commands through the engine
    pub fn do_multi(&mut self, viewer: &mut Viewer, cmds: &str) -> CmdResult {
        self.executor.do_multi(viewer, cmds)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_script_engine_creation() {
        let engine = ScriptEngine::new();
        assert!(engine.script_stack.is_empty());
    }
}
