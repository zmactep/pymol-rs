//! Command line autocomplete/completion logic
//!
//! Provides tab completion for:
//! - Command names (when typing the first word)
//! - File/directory paths (for I/O commands like load, save, png, cd, ls)
//! - Setting names (for set, get, unset)
//! - Color names (for color, bg_color)
//! - Representation names (for show, hide, as)
//! - Object/selection names (for zoom, center, delete, etc.)

use std::path::Path;

use pymol_cmd::{ArgHint, CommandRegistry};

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

/// Static list of representation names for completion
const REPRESENTATION_NAMES: &[&str] = &[
    "cartoon", "sticks", "lines", "spheres", "surface", "mesh",
    "ribbon", "dots", "nonbonded", "nb_spheres", "labels", "cell",
    "cgo", "everything",
];

/// Context for generating completions, holding all available name lists
pub struct CompletionContext<'a> {
    pub command_names: &'a [&'a str],
    pub registry: &'a CommandRegistry,
    pub setting_names: &'a [&'a str],
    pub color_names: &'a [String],
    pub object_names: &'a [String],
    pub selection_names: &'a [String],
}

/// Generate completions for the current command input
pub fn generate_completions(
    input: &str,
    cursor_pos: usize,
    ctx: &CompletionContext,
) -> CompletionResult {
    let cursor_pos = cursor_pos.min(input.len());
    let text_before_cursor = &input[..cursor_pos];

    // Parse into command + segments (comma-separated, bracket-aware)
    let parsed = parse_input(text_before_cursor);

    match parsed {
        ParsedInput::Empty => {
            // Show all commands
            complete_from_list("", cursor_pos, ctx.command_names)
        }
        ParsedInput::CommandPrefix(prefix) => {
            // Completing command name
            complete_from_list(prefix, cursor_pos - prefix.len(), ctx.command_names)
        }
        ParsedInput::Argument { command, arg_index, prefix, prefix_start } => {
            // Look up arg hints for this command
            if let Some(cmd) = ctx.registry.get(&command.to_lowercase()) {
                let hints: &[ArgHint] = cmd.arg_hints();

                // Special case: "set" arg_index=1 → setting value completion
                if command.eq_ignore_ascii_case("set") && arg_index == 1 {
                    return complete_setting_value(text_before_cursor, prefix, prefix_start, ctx);
                }

                if let Some(hint) = hints.get(arg_index) {
                    match hint {
                        ArgHint::Path => complete_path(text_before_cursor, cursor_pos),
                        ArgHint::Setting => complete_from_list(prefix, prefix_start, ctx.setting_names),
                        ArgHint::Color => {
                            let mut suggestions = filter_prefix(prefix, ctx.color_names.iter().map(|s| s.as_str()));
                            suggestions.extend(filter_prefix(prefix, pymol_color::SCHEME_NAMES.iter().copied()));
                            suggestions.sort();
                            CompletionResult { start_pos: prefix_start, suggestions }
                        }
                        ArgHint::Representation => complete_from_static(prefix, prefix_start, REPRESENTATION_NAMES),
                        ArgHint::Object => complete_from_list(prefix, prefix_start, &ctx.object_names.iter().map(|s| s.as_str()).collect::<Vec<_>>()),
                        ArgHint::Selection => {
                            // Object names + selection names
                            let mut suggestions = filter_prefix(prefix, ctx.object_names.iter().map(|s| s.as_str()));
                            suggestions.extend(filter_prefix(prefix, ctx.selection_names.iter().map(|s| s.as_str())));
                            suggestions.sort();
                            suggestions.dedup();
                            CompletionResult { start_pos: prefix_start, suggestions }
                        }
                        ArgHint::NamedSelection => {
                            complete_from_list(prefix, prefix_start, &ctx.selection_names.iter().map(|s| s.as_str()).collect::<Vec<_>>())
                        }
                        ArgHint::None => CompletionResult::empty(cursor_pos),
                    }
                } else {
                    // No hint for this position - no completion
                    CompletionResult::empty(cursor_pos)
                }
            } else {
                // Unknown command
                CompletionResult::empty(cursor_pos)
            }
        }
    }
}

/// Parsed representation of the input before cursor
enum ParsedInput<'a> {
    /// Empty input
    Empty,
    /// Typing the first word (command name)
    CommandPrefix(&'a str),
    /// Typing an argument
    Argument {
        command: &'a str,
        arg_index: usize,
        prefix: &'a str,
        prefix_start: usize,
    },
}

/// Parse input into command and current argument position
fn parse_input(text: &str) -> ParsedInput<'_> {
    if text.is_empty() {
        return ParsedInput::Empty;
    }

    // Split into segments by commas (bracket-aware)
    let segments = split_by_commas(text);

    if segments.is_empty() {
        return ParsedInput::Empty;
    }

    // First segment: "command arg0_prefix" or just "command_prefix"
    let first_seg = segments[0];
    let first_words: Vec<&str> = first_seg.split_whitespace().collect();

    if first_words.is_empty() {
        return ParsedInput::Empty;
    }

    let ends_with_space = first_seg.ends_with(char::is_whitespace);

    if segments.len() == 1 && first_words.len() == 1 && !ends_with_space {
        // Still typing the command name
        return ParsedInput::CommandPrefix(first_words[0]);
    }

    let command = first_words[0];

    if segments.len() == 1 {
        // Still in first segment, after command name
        // arg_index = 0, prefix is everything after command
        let cmd_end = first_seg.find(char::is_whitespace).unwrap_or(first_seg.len());
        let after_cmd = &first_seg[cmd_end..];
        let prefix = after_cmd.trim_start();
        let prefix_start = text.len() - prefix.len();
        return ParsedInput::Argument {
            command,
            arg_index: 0,
            prefix,
            prefix_start,
        };
    }

    // We're in segment N (N >= 1), which means arg_index = N - 1... wait:
    // segment 0 = "command arg0", segment 1 = " arg1", segment 2 = " arg2"
    // So arg_index = segments.len() - 1
    let arg_index = segments.len() - 1;
    let last_seg = segments[arg_index];
    let prefix = last_seg.trim_start();
    let prefix_start = text.len() - prefix.len();

    ParsedInput::Argument {
        command,
        arg_index,
        prefix,
        prefix_start,
    }
}

/// Split text by commas, respecting brackets and quotes
fn split_by_commas(text: &str) -> Vec<&str> {
    let mut segments = Vec::new();
    let mut depth = 0i32;
    let mut in_quote = false;
    let mut last_start = 0;

    for (i, ch) in text.char_indices() {
        match ch {
            '"' | '\'' => in_quote = !in_quote,
            '(' | '[' | '{' if !in_quote => depth += 1,
            ')' | ']' | '}' if !in_quote => depth -= 1,
            ',' if !in_quote && depth <= 0 => {
                segments.push(&text[last_start..i]);
                last_start = i + 1;
            }
            _ => {}
        }
    }
    segments.push(&text[last_start..]);
    segments
}

/// Complete setting value (special case for `set name, <value>`)
fn complete_setting_value(
    text: &str,
    prefix: &str,
    prefix_start: usize,
    ctx: &CompletionContext,
) -> CompletionResult {
    // Extract setting name from arg 0
    let segments = split_by_commas(text);
    if segments.is_empty() {
        return CompletionResult::empty(prefix_start);
    }

    let first_seg = segments[0];
    let first_words: Vec<&str> = first_seg.split_whitespace().collect();
    let setting_name = if first_words.len() >= 2 { first_words[1] } else { return CompletionResult::empty(prefix_start); };

    // Look up setting type
    if let Some(id) = pymol_settings::get_setting_id(setting_name) {
        if let Some(setting) = pymol_settings::get_setting(id) {
            // Settings with named value variants (e.g., shading_mode: classic/skripkin)
            if setting.has_value_hints() {
                let hints: Vec<&str> = setting.hint_names().collect();
                return complete_from_static(prefix, prefix_start, &hints);
            }
            match setting.setting_type {
                pymol_settings::SettingType::Bool => {
                    return complete_from_static(prefix, prefix_start, &["on", "off"]);
                }
                pymol_settings::SettingType::Color => {
                    let mut suggestions = filter_prefix(prefix, ctx.color_names.iter().map(|s| s.as_str()));
                    suggestions.extend(filter_prefix(prefix, pymol_color::SCHEME_NAMES.iter().copied()));
                    suggestions.sort();
                    return CompletionResult { start_pos: prefix_start, suggestions };
                }
                _ => return CompletionResult::empty(prefix_start),
            };
        }
    }

    CompletionResult::empty(prefix_start)
}

/// Filter items by case-insensitive prefix
fn filter_prefix<'a>(prefix: &str, items: impl Iterator<Item = &'a str>) -> Vec<String> {
    let prefix_lower = prefix.to_lowercase();
    items
        .filter(|item| item.to_lowercase().starts_with(&prefix_lower))
        .map(|s| s.to_string())
        .collect()
}

/// Complete from a dynamic list of names
fn complete_from_list(prefix: &str, prefix_start: usize, names: &[&str]) -> CompletionResult {
    let prefix_lower = prefix.to_lowercase();
    let mut matches: Vec<String> = names
        .iter()
        .filter(|name| name.to_lowercase().starts_with(&prefix_lower))
        .map(|s| s.to_string())
        .collect();
    matches.sort();
    CompletionResult { start_pos: prefix_start, suggestions: matches }
}

/// Complete from a static list
fn complete_from_static(prefix: &str, start_pos: usize, names: &[&str]) -> CompletionResult {
    let prefix_lower = prefix.to_lowercase();
    let mut matches: Vec<String> = names
        .iter()
        .filter(|name| name.to_lowercase().starts_with(&prefix_lower))
        .map(|s| s.to_string())
        .collect();
    matches.sort();
    CompletionResult { start_pos, suggestions: matches }
}

/// Complete file/directory paths
fn complete_path(input: &str, _cursor_pos: usize) -> CompletionResult {
    // Find where the path argument starts
    let path_start = find_path_start(input);
    let path_portion = &input[path_start..];

    // Expand tilde if present
    let expanded = expand_tilde(path_portion);
    let expanded_path = Path::new(&expanded);

    // Determine the directory to search and the prefix to match
    let (search_dir, prefix) = if expanded_path.is_dir() && path_portion.ends_with('/') {
        (expanded_path.to_path_buf(), String::new())
    } else if let Some(parent) = expanded_path.parent() {
        let prefix = expanded_path
            .file_name()
            .map(|s| s.to_string_lossy().to_string())
            .unwrap_or_default();
        // Empty parent means bare filename (e.g. "1" without "./" prefix)
        // Use cwd so we get relative paths
        let dir = if parent.as_os_str().is_empty() {
            Path::new(".").to_path_buf()
        } else {
            parent.to_path_buf()
        };
        (dir, prefix)
    } else {
        (Path::new(".").to_path_buf(), expanded.clone())
    };

    let mut suggestions = Vec::new();

    if let Ok(entries) = std::fs::read_dir(&search_dir) {
        for entry in entries.filter_map(|e| e.ok()) {
            let name = entry.file_name().to_string_lossy().to_string();

            if name.starts_with('.') && !prefix.starts_with('.') {
                continue;
            }

            if name.to_lowercase().starts_with(&prefix.to_lowercase()) {
                let mut suggestion = if path_portion.starts_with('~') {
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

                if entry.path().is_dir() {
                    suggestion.push('/');
                }

                suggestions.push(suggestion);
            }
        }
    }

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
    let mut chars = input.char_indices().peekable();

    // Skip command (first word)
    while let Some((_, c)) = chars.peek() {
        if c.is_whitespace() { break; }
        chars.next();
    }

    // Skip whitespace after command
    while let Some((_, c)) = chars.peek() {
        if !c.is_whitespace() { break; }
        chars.next();
    }

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
    fn test_split_by_commas() {
        assert_eq!(split_by_commas("set sphere_scale, 0.5"), vec!["set sphere_scale", " 0.5"]);
        assert_eq!(split_by_commas("color green, chain A"), vec!["color green", " chain A"]);
        assert_eq!(split_by_commas("show cartoon"), vec!["show cartoon"]);
        // Brackets
        assert_eq!(split_by_commas("set bg_rgb, [1.0, 0.5, 0.0]"), vec!["set bg_rgb", " [1.0, 0.5, 0.0]"]);
    }

    #[test]
    fn test_parse_input_empty() {
        assert!(matches!(parse_input(""), ParsedInput::Empty));
    }

    #[test]
    fn test_parse_input_command_prefix() {
        match parse_input("lo") {
            ParsedInput::CommandPrefix("lo") => {}
            _ => panic!("Expected CommandPrefix"),
        }
    }

    #[test]
    fn test_parse_input_argument() {
        match parse_input("set sphere_sc") {
            ParsedInput::Argument { command: "set", arg_index: 0, prefix: "sphere_sc", .. } => {}
            _ => panic!("Expected Argument for arg 0"),
        }

        match parse_input("set sphere_scale, ") {
            ParsedInput::Argument { command: "set", arg_index: 1, prefix: "", .. } => {}
            _ => panic!("Expected Argument for arg 1"),
        }

        match parse_input("color green, ") {
            ParsedInput::Argument { command: "color", arg_index: 1, prefix: "", .. } => {}
            _ => panic!("Expected Argument for arg 1"),
        }
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
