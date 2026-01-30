//! Command line autocomplete/completion logic
//!
//! Provides tab completion for:
//! - Command names (when typing the first word)
//! - File/directory paths (for I/O commands like load, save, png, cd, ls)

use std::path::Path;

/// Result of completion generation
#[derive(Debug)]
pub struct CompletionResult {
    /// Byte position in input where the completion starts
    pub start_pos: usize,
    /// List of completion suggestions
    pub suggestions: Vec<String>,
}

impl CompletionResult {
    /// Create an empty result (no completions)
    pub fn empty(pos: usize) -> Self {
        Self {
            start_pos: pos,
            suggestions: Vec::new(),
        }
    }
}

/// Generate completions for the current command input
///
/// # Arguments
/// * `input` - The current command line text
/// * `cursor_pos` - Byte position of the cursor in the input
/// * `command_names` - List of available command names for completion
/// * `path_commands` - List of command names that take file paths as first argument
///
/// # Returns
/// A `CompletionResult` with the start position and list of suggestions
pub fn generate_completions(
    input: &str,
    cursor_pos: usize,
    command_names: &[&str],
    path_commands: &[&str],
) -> CompletionResult {
    // Ensure cursor position is valid
    let cursor_pos = cursor_pos.min(input.len());
    let text_before_cursor = &input[..cursor_pos];

    // Parse words before cursor
    let words: Vec<&str> = text_before_cursor.split_whitespace().collect();
    let ends_with_space = text_before_cursor.ends_with(char::is_whitespace);

    if words.is_empty() {
        // Empty input - show all commands
        complete_command("", cursor_pos, command_names)
    } else if words.len() == 1 && !ends_with_space {
        // Typing first word - complete command name
        complete_command(words[0], cursor_pos, command_names)
    } else if is_path_command(&words[0].to_lowercase(), path_commands) {
        // First word is a path command - complete file path
        complete_path(text_before_cursor, cursor_pos)
    } else {
        // No completion for other cases
        CompletionResult::empty(cursor_pos)
    }
}

/// Check if a command name is in the path commands list (case-insensitive)
fn is_path_command(cmd: &str, path_commands: &[&str]) -> bool {
    path_commands.iter().any(|pc| pc.eq_ignore_ascii_case(cmd))
}

/// Complete command names
fn complete_command(prefix: &str, cursor_pos: usize, command_names: &[&str]) -> CompletionResult {
    let prefix_lower = prefix.to_lowercase();
    let mut matches: Vec<String> = command_names
        .iter()
        .filter(|cmd| cmd.to_lowercase().starts_with(&prefix_lower))
        .map(|s| s.to_string())
        .collect();

    // Sort alphabetically for consistent display
    matches.sort();

    CompletionResult {
        start_pos: cursor_pos - prefix.len(),
        suggestions: matches,
    }
}

/// Complete file/directory paths
fn complete_path(input: &str, _cursor_pos: usize) -> CompletionResult {
    // Find where the path argument starts
    // Skip the command and find the start of the path portion
    let path_start = find_path_start(input);
    let path_portion = &input[path_start..];

    // Expand tilde if present
    let expanded = expand_tilde(path_portion);
    let expanded_path = Path::new(&expanded);

    // Determine the directory to search and the prefix to match
    let (search_dir, prefix) = if expanded_path.is_dir() && path_portion.ends_with('/') {
        // Path ends with / - list contents of that directory
        (expanded_path.to_path_buf(), String::new())
    } else if let Some(parent) = expanded_path.parent() {
        // Get parent directory and filename prefix
        let prefix = expanded_path
            .file_name()
            .map(|s| s.to_string_lossy().to_string())
            .unwrap_or_default();
        (parent.to_path_buf(), prefix)
    } else {
        // No parent - use current directory
        (std::env::current_dir().unwrap_or_else(|_| Path::new(".").to_path_buf()), expanded.clone())
    };

    // Read directory contents and filter by prefix
    let mut suggestions = Vec::new();

    if let Ok(entries) = std::fs::read_dir(&search_dir) {
        for entry in entries.filter_map(|e| e.ok()) {
            let name = entry.file_name().to_string_lossy().to_string();

            // Skip hidden files unless prefix starts with '.'
            if name.starts_with('.') && !prefix.starts_with('.') {
                continue;
            }

            if name.to_lowercase().starts_with(&prefix.to_lowercase()) {
                // Build the full suggestion
                let mut suggestion = if path_portion.starts_with('~') {
                    // Keep tilde in suggestion
                    format!("~/{}", 
                        search_dir
                            .strip_prefix(dirs::home_dir().unwrap_or_default())
                            .unwrap_or(&search_dir)
                            .join(&name)
                            .display()
                    )
                } else if search_dir.as_os_str().is_empty() || search_dir == Path::new(".") {
                    name.clone()
                } else {
                    search_dir.join(&name).display().to_string()
                };

                // Add trailing slash for directories
                if entry.path().is_dir() {
                    suggestion.push('/');
                }

                suggestions.push(suggestion);
            }
        }
    }

    // Sort alphabetically, directories first
    suggestions.sort_by(|a, b| {
        let a_is_dir = a.ends_with('/');
        let b_is_dir = b.ends_with('/');
        match (a_is_dir, b_is_dir) {
            (true, false) => std::cmp::Ordering::Less,
            (false, true) => std::cmp::Ordering::Greater,
            _ => a.to_lowercase().cmp(&b.to_lowercase()),
        }
    });

    CompletionResult {
        start_pos: path_start,
        suggestions,
    }
}

/// Find the byte position where the path argument starts in the input
fn find_path_start(input: &str) -> usize {
    // Skip the command name and any following whitespace
    let mut chars = input.char_indices().peekable();

    // Skip command (first word)
    while let Some((_, c)) = chars.peek() {
        if c.is_whitespace() {
            break;
        }
        chars.next();
    }

    // Skip whitespace after command
    while let Some((_, c)) = chars.peek() {
        if !c.is_whitespace() {
            break;
        }
        chars.next();
    }

    // Return position of path start, or end of input
    chars.peek().map(|(i, _)| *i).unwrap_or(input.len())
}

/// Expand tilde (~) to home directory
fn expand_tilde(path: &str) -> String {
    if path.starts_with('~') {
        if let Some(home) = dirs::home_dir() {
            return path.replacen('~', &home.display().to_string(), 1);
        }
    }
    path.to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_command_completion() {
        let commands = vec!["load", "label", "ls", "log", "zoom"];
        let path_commands = vec!["load", "ls"];

        // Complete "l"
        let result = generate_completions("l", 1, &commands, &path_commands);
        assert_eq!(result.start_pos, 0);
        assert!(result.suggestions.contains(&"load".to_string()));
        assert!(result.suggestions.contains(&"label".to_string()));
        assert!(result.suggestions.contains(&"ls".to_string()));
        assert!(result.suggestions.contains(&"log".to_string()));
        assert!(!result.suggestions.contains(&"zoom".to_string()));

        // Complete "lo"
        let result = generate_completions("lo", 2, &commands, &path_commands);
        assert_eq!(result.start_pos, 0);
        assert!(result.suggestions.contains(&"load".to_string()));
        assert!(result.suggestions.contains(&"log".to_string()));
        assert!(!result.suggestions.contains(&"label".to_string()));
    }

    #[test]
    fn test_find_path_start() {
        assert_eq!(find_path_start("load /path/to/file"), 5);
        assert_eq!(find_path_start("load  /path"), 6);
        assert_eq!(find_path_start("cd ~/"), 3);
        assert_eq!(find_path_start("ls"), 2);
    }

    #[test]
    fn test_expand_tilde() {
        let expanded = expand_tilde("~/test");
        assert!(!expanded.starts_with('~'));
        assert!(expanded.ends_with("/test"));
    }
}
