//! Command parser using nom
//!
//! Parses PyMOL command strings into structured `ParsedCommand` objects.
//!
//! # Supported Syntax
//!
//! - Simple commands: `load file.pdb`
//! - Named arguments: `load file.pdb, object=mol`
//! - Mixed: `color red, selection=all`
//! - Quoted strings: `load "path with spaces.pdb"`
//! - Selections with parens: `select (name CA and chain A)`
//! - Multiple commands: `load file.pdb; zoom; show cartoon`
//! - Legacy syntax: `select name=CA` -> `select name, CA`

use nom::{
    branch::alt,
    bytes::complete::{escaped, tag, take_till, take_while, take_while1},
    character::complete::{char, multispace0, none_of, one_of},
    combinator::{map, recognize},
    multi::{many0, separated_list0},
    number::complete::recognize_float,
    sequence::{delimited, pair, preceded, tuple},
    IResult,
};

use crate::args::{ArgValue, ParsedCommand};
use crate::error::ParseError;

/// Normalize line continuations in a script/command string.
///
/// Handles:
/// - `\` followed by newline (standard line continuation)
/// - `\` followed by whitespace (collapsed line continuation, e.g. when newlines are removed)
///
/// This allows PyMOL-style multi-line commands like:
/// ```text
/// set_view (\
///   1.0, 0.0, 0.0,\
///   0.0, 1.0, 0.0)
/// ```
fn normalize_continuations(input: &str) -> String {
    let mut result = String::with_capacity(input.len());
    let mut chars = input.chars().peekable();

    while let Some(c) = chars.next() {
        if c == '\\' {
            // Check what follows the backslash
            if let Some(&next) = chars.peek() {
                if next == '\n' {
                    // Skip \ and newline (line continuation)
                    chars.next();
                    // Also skip any leading whitespace on next line
                    while chars.peek().map_or(false, |c| *c == ' ' || *c == '\t') {
                        chars.next();
                    }
                    continue;
                } else if next == ' ' || next == '\t' {
                    // Skip \ followed by whitespace (collapsed continuation)
                    chars.next();
                    // Skip any additional whitespace
                    while chars.peek().map_or(false, |c| *c == ' ' || *c == '\t') {
                        chars.next();
                    }
                    continue;
                }
            }
        }
        result.push(c);
    }
    result
}

/// Join lines that end with backslash (line continuation).
///
/// This handles PyMOL's line continuation syntax where a line ending with `\`
/// continues on the next line. This is useful for preprocessing script files
/// before line-by-line execution.
///
/// # Example
/// ```
/// use pymol_cmd::join_continued_lines;
///
/// let input = "set_view (\\\n  1.0, 2.0)";
/// let result = join_continued_lines(input);
/// assert!(result.contains("set_view"));
/// assert!(result.contains("1.0"));
/// ```
pub fn join_continued_lines(script: &str) -> String {
    let mut result = String::with_capacity(script.len());
    let mut continuation = false;

    for line in script.lines() {
        if continuation {
            // Continuing from previous line - add a space instead of newline
            result.push(' ');
        }

        if line.trim_end().ends_with('\\') {
            // Line continues - add without trailing backslash
            let trimmed = line.trim_end().trim_end_matches('\\');
            result.push_str(trimmed);
            continuation = true;
        } else {
            result.push_str(line);
            result.push('\n');
            continuation = false;
        }
    }
    result
}

/// Parse a single command from a string
///
/// # Example
/// ```
/// use pymol_cmd::parse_command;
///
/// let cmd = parse_command("load protein.pdb, object=mol").unwrap();
/// assert_eq!(cmd.name, "load");
/// assert_eq!(cmd.get_str(0), Some("protein.pdb"));
/// assert_eq!(cmd.get_named_str("object"), Some("mol"));
/// ```
pub fn parse_command(input: &str) -> Result<ParsedCommand, ParseError> {
    let input = normalize_continuations(input);
    let input = input.trim();
    if input.is_empty() {
        return Err(ParseError::EmptyCommand);
    }

    match parse_single_command(input) {
        Ok(("", cmd)) => Ok(cmd),
        Ok((remaining, _)) => Err(ParseError::Generic(format!(
            "unexpected trailing input: '{}'",
            remaining
        ))),
        Err(e) => Err(e.into()),
    }
}

/// Parse multiple commands separated by semicolons or newlines
///
/// # Example
/// ```
/// use pymol_cmd::parse_commands;
///
/// let cmds = parse_commands("load file.pdb; zoom; show cartoon").unwrap();
/// assert_eq!(cmds.len(), 3);
/// assert_eq!(cmds[0].name, "load");
/// assert_eq!(cmds[1].name, "zoom");
/// assert_eq!(cmds[2].name, "show");
/// ```
pub fn parse_commands(input: &str) -> Result<Vec<ParsedCommand>, ParseError> {
    let input = normalize_continuations(input);
    let input = input.trim();
    if input.is_empty() {
        return Ok(Vec::new());
    }

    let mut commands = Vec::new();
    let mut current: &str = &input;

    while !current.is_empty() {
        // Skip leading whitespace and semicolons
        current = current.trim_start();
        if current.starts_with(';') {
            current = &current[1..];
            continue;
        }

        // Handle comments
        if current.starts_with('#') {
            // Skip to end of line
            if let Some(newline_pos) = current.find('\n') {
                current = &current[newline_pos + 1..];
            } else {
                break;
            }
            continue;
        }

        if current.is_empty() {
            break;
        }

        // Find the end of this command (semicolon or newline, but not inside quotes/parens)
        let end = find_command_end(current);
        let cmd_str = &current[..end];
        current = &current[end..];

        let cmd_str = cmd_str.trim();
        if cmd_str.is_empty() {
            continue;
        }

        match parse_single_command(cmd_str) {
            Ok(("", cmd)) => commands.push(cmd),
            Ok((remaining, cmd)) => {
                // If there's remaining input, it might be another command
                if !remaining.trim().is_empty() {
                    commands.push(cmd);
                    current = remaining;
                } else {
                    commands.push(cmd);
                }
            }
            Err(e) => return Err(e.into()),
        }
    }

    Ok(commands)
}

/// Find the end of a command (semicolon or newline), respecting quotes and parens
fn find_command_end(input: &str) -> usize {
    let mut depth: i32 = 0; // Paren/bracket depth
    let mut in_string = false;
    let mut string_char = '"';
    let mut chars = input.char_indices().peekable();

    while let Some((i, c)) = chars.next() {
        if in_string {
            if c == '\\' {
                // Skip escaped character
                chars.next();
            } else if c == string_char {
                in_string = false;
            }
        } else {
            match c {
                '"' | '\'' => {
                    in_string = true;
                    string_char = c;
                }
                '(' | '[' | '{' => depth += 1,
                ')' | ']' | '}' => depth = depth.saturating_sub(1),
                ';' | '\n' if depth == 0 => return i,
                _ => {}
            }
        }
    }

    input.len()
}

/// Parse a single command (name and arguments)
fn parse_single_command(input: &str) -> IResult<&str, ParsedCommand> {
    let (input, _) = multispace0(input)?;
    let (input, name) = parse_command_name(input)?;
    let (input, _) = multispace0(input)?;

    // Parse arguments if present
    let (input, args) = if input.is_empty() || input.starts_with(';') || input.starts_with('\n') {
        (input, Vec::new())
    } else {
        parse_arguments(input)?
    };

    Ok((
        input,
        ParsedCommand {
            name: name.to_string(),
            args,
        },
    ))
}

/// Parse a command name (alphanumeric + underscore + dot for util.xxx, or @ for script execution)
fn parse_command_name(input: &str) -> IResult<&str, &str> {
    alt((
        // Special case: @ command for script execution
        tag("@"),
        // Regular command names: alphanumeric + underscore + dot
        recognize(pair(
            take_while1(|c: char| c.is_alphanumeric() || c == '_'),
            take_while(|c: char| c.is_alphanumeric() || c == '_' || c == '.'),
        )),
    ))(input)
}

/// Parse command arguments (comma-separated or space-separated list)
///
/// PyMOL accepts both comma and space as argument separators in many commands.
/// For example, `mset 1 x60` uses space separation while `color red, chain A`
/// uses comma separation. This parser handles both.
fn parse_arguments(input: &str) -> IResult<&str, Vec<(Option<String>, ArgValue)>> {
    let (input, first) = parse_argument(input)?;

    // Try comma-separated args first, then fall back to space-separated
    let (input, rest) = many0(preceded(
        alt((
            // Comma separator (with optional whitespace)
            map(tuple((multispace0, char(','), multispace0)), |_| ()),
            // Space separator (only whitespace, no comma) — must have at least one space
            // Only match if the next char is not a comma or end-of-input
            map(take_while1(|c: char| c == ' ' || c == '\t'), |_| ()),
        )),
        parse_argument,
    ))(input)?;

    let mut args = vec![first];
    args.extend(rest);
    Ok((input, args))
}

/// Parse a single argument (possibly named)
fn parse_argument(input: &str) -> IResult<&str, (Option<String>, ArgValue)> {
    let (input, _) = multispace0(input)?;

    // Try to parse as named argument (name=value)
    if let Ok((remaining, (name, _, _, value))) = tuple((
        parse_arg_name,
        multispace0,
        char('='),
        preceded(multispace0, parse_arg_value),
    ))(input)
    {
        return Ok((remaining, (Some(name.to_string()), value)));
    }

    // Parse as positional argument
    let (input, value) = parse_arg_value(input)?;
    Ok((input, (None, value)))
}

/// Parse an argument name (for named arguments)
fn parse_arg_name(input: &str) -> IResult<&str, &str> {
    take_while1(|c: char| c.is_alphanumeric() || c == '_')(input)
}

/// Parse an argument value
fn parse_arg_value(input: &str) -> IResult<&str, ArgValue> {
    alt((
        // Try parenthesized expression (selection or value list like set_view)
        map(parse_paren_expr, |s| ArgValue::String(s.to_string())),
        // Try bracketed list
        parse_list,
        // Try quoted string
        map(parse_quoted_string, |s| ArgValue::String(s)),
        // Try number
        parse_number,
        // Try boolean keywords
        parse_bool,
        // Try unquoted string/identifier
        map(parse_unquoted_value, |s| ArgValue::String(s.to_string())),
    ))(input)
}

/// Parse a parenthesized expression (typically a selection)
///
/// Handles both simple parenthesized expressions like `(1.0, 2.0)` and
/// selection algebra like `(chain A and resi 29) or (chain A and resi 31)`.
/// After matching the initial `(...)`, continues to consume selection operators
/// (`and`, `or`, `not`, `&`, `|`, `+`, `-`, `!`) and any subsequent terms.
fn parse_paren_expr(input: &str) -> IResult<&str, &str> {
    let (rest, _) = recognize(delimited(
        char('('),
        take_balanced_parens,
        char(')'),
    ))(input)?;

    // After the first (...), check if there are selection operators followed by more content.
    // This handles expressions like: (chain A and resi 29) or (chain A and resi 31)
    let mut end = input.len() - rest.len();
    let mut remaining = rest;

    loop {
        // Skip whitespace
        let trimmed = remaining.trim_start();
        if trimmed.is_empty() {
            break;
        }

        // Check for selection operators: and, or, not, in, like, &, |, +, -, !
        let has_operator = trimmed.starts_with("and ")
            || trimmed.starts_with("or ")
            || trimmed.starts_with("not ")
            || trimmed.starts_with("in ")
            || trimmed.starts_with("like ")
            || trimmed.starts_with("& ")
            || trimmed.starts_with("| ")
            || trimmed.starts_with("+ ")
            || trimmed.starts_with("- ")
            || trimmed.starts_with("! ");

        if !has_operator {
            break;
        }

        // Consume the rest of the selection expression after the operator
        let consumed_ws = remaining.len() - trimmed.len();
        let taken: usize = trimmed
            .chars()
            .take_while(|c| !matches!(c, ',' | ';' | '\n' | '\r' | '[' | ']'))
            .map(|c| c.len_utf8())
            .sum();
        if taken == 0 {
            break;
        }
        let taken_trimmed = trimmed[..taken].trim_end().len();
        end += consumed_ws + taken_trimmed;
        remaining = &input[end..];
    }

    Ok((&input[end..], &input[..end]))
}

/// Take content with balanced parentheses
fn take_balanced_parens(input: &str) -> IResult<&str, &str> {
    let mut depth = 1;
    let mut end = 0;
    let mut in_string = false;
    let mut string_char = '"';
    let mut chars = input.chars().peekable();

    while let Some(c) = chars.next() {
        if in_string {
            if c == '\\' {
                chars.next();
                end += 2;
                continue;
            } else if c == string_char {
                in_string = false;
            }
            end += c.len_utf8();
        } else {
            match c {
                '"' | '\'' => {
                    in_string = true;
                    string_char = c;
                    end += 1;
                }
                '(' => {
                    depth += 1;
                    end += 1;
                }
                ')' => {
                    depth -= 1;
                    if depth == 0 {
                        return Ok((&input[end..], &input[..end]));
                    }
                    end += 1;
                }
                _ => end += c.len_utf8(),
            }
        }
    }

    // Unbalanced - return all
    Ok(("", input))
}

/// Parse a bracketed list
fn parse_list(input: &str) -> IResult<&str, ArgValue> {
    let (input, _) = char('[')(input)?;
    let (input, _) = multispace0(input)?;
    let (input, items) = separated_list0(
        tuple((multispace0, char(','), multispace0)),
        parse_arg_value,
    )(input)?;
    let (input, _) = multispace0(input)?;
    let (input, _) = char(']')(input)?;

    Ok((input, ArgValue::List(items)))
}

/// Parse a quoted string (single or double quotes)
fn parse_quoted_string(input: &str) -> IResult<&str, String> {
    alt((
        // Double-quoted string
        map(
            delimited(
                char('"'),
                escaped(none_of("\"\\"), '\\', one_of("\"\\nrt")),
                char('"'),
            ),
            |s: &str| unescape_string(s),
        ),
        // Single-quoted string
        map(
            delimited(
                char('\''),
                escaped(none_of("'\\"), '\\', one_of("'\\nrt")),
                char('\''),
            ),
            |s: &str| unescape_string(s),
        ),
        // Triple-quoted strings (for multi-line)
        map(
            delimited(tag("'''"), take_till(|c| c == '\''), tag("'''")),
            |s: &str| s.to_string(),
        ),
        map(
            delimited(tag("\"\"\""), take_till(|c| c == '"'), tag("\"\"\"")),
            |s: &str| s.to_string(),
        ),
    ))(input)
}

/// Parse a number (int or float)
fn parse_number(input: &str) -> IResult<&str, ArgValue> {
    let (remaining, num_str) = recognize_float(input)?;

    // Check if it's followed by a valid separator or end
    if !remaining.is_empty() {
        let next_char = remaining.chars().next().unwrap();
        if next_char.is_alphanumeric() || next_char == '_' || next_char == '.' {
            // Not a number, it's part of an identifier like "chain1"
            return Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Float,
            )));
        }
    }

    // Try to parse as integer first
    if let Ok(i) = num_str.parse::<i64>() {
        return Ok((remaining, ArgValue::Int(i)));
    }

    // Parse as float
    if let Ok(f) = num_str.parse::<f64>() {
        return Ok((remaining, ArgValue::Float(f)));
    }

    Err(nom::Err::Error(nom::error::Error::new(
        input,
        nom::error::ErrorKind::Float,
    )))
}

/// Parse boolean keywords
fn parse_bool(input: &str) -> IResult<&str, ArgValue> {
    let (remaining, keyword) = alt((
        tag("true"),
        tag("false"),
        tag("yes"),
        tag("no"),
        tag("on"),
        tag("off"),
    ))(input)?;

    // Ensure word boundary - next char must not be alphanumeric or underscore
    if !remaining.is_empty() {
        let next_char = remaining.chars().next().unwrap();
        if next_char.is_alphanumeric() || next_char == '_' {
            return Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Tag,
            )));
        }
    }

    let value = match keyword {
        "true" | "yes" | "on" => ArgValue::Bool(true),
        "false" | "no" | "off" => ArgValue::Bool(false),
        _ => unreachable!(),
    };

    Ok((remaining, value))
}

/// Parse an unquoted value (stops at comma, semicolon, or whitespace before named arg)
fn parse_unquoted_value(input: &str) -> IResult<&str, &str> {
    // Take characters until we hit a delimiter
    take_while1(|c: char| {
        !matches!(c, ',' | ';' | '\n' | '\r' | '[' | ']')
    })(input)
    .map(|(remaining, value)| (remaining, value.trim_end()))
}

/// Unescape a string (handle \n, \t, etc.)
fn unescape_string(s: &str) -> String {
    let mut result = String::with_capacity(s.len());
    let mut chars = s.chars();

    while let Some(c) = chars.next() {
        if c == '\\' {
            match chars.next() {
                Some('n') => result.push('\n'),
                Some('r') => result.push('\r'),
                Some('t') => result.push('\t'),
                Some('\\') => result.push('\\'),
                Some('"') => result.push('"'),
                Some('\'') => result.push('\''),
                Some(c) => {
                    result.push('\\');
                    result.push(c);
                }
                None => result.push('\\'),
            }
        } else {
            result.push(c);
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_command() {
        let cmd = parse_command("zoom").unwrap();
        assert_eq!(cmd.name, "zoom");
        assert_eq!(cmd.args.len(), 0);
    }

    #[test]
    fn test_command_with_arg() {
        let cmd = parse_command("load protein.pdb").unwrap();
        assert_eq!(cmd.name, "load");
        assert_eq!(cmd.get_str(0), Some("protein.pdb"));
    }

    #[test]
    fn test_command_with_named_arg() {
        let cmd = parse_command("load protein.pdb, object=mol").unwrap();
        assert_eq!(cmd.name, "load");
        assert_eq!(cmd.get_str(0), Some("protein.pdb"));
        assert_eq!(cmd.get_named_str("object"), Some("mol"));
    }

    #[test]
    fn test_quoted_string() {
        let cmd = parse_command(r#"load "path with spaces.pdb""#).unwrap();
        assert_eq!(cmd.get_str(0), Some("path with spaces.pdb"));
    }

    #[test]
    fn test_selection_in_parens() {
        let cmd = parse_command("color red, (name CA and chain A)").unwrap();
        assert_eq!(cmd.name, "color");
        assert_eq!(cmd.get_str(0), Some("red"));
        assert_eq!(cmd.get_str(1), Some("(name CA and chain A)"));
    }

    #[test]
    fn test_numbers() {
        let cmd = parse_command("set sphere_scale, 0.5").unwrap();
        assert_eq!(cmd.get_float(1), Some(0.5));

        let cmd = parse_command("frame 10").unwrap();
        assert_eq!(cmd.get_int(0), Some(10));
    }

    #[test]
    fn test_multiple_commands() {
        let cmds = parse_commands("load file.pdb; zoom; show cartoon").unwrap();
        assert_eq!(cmds.len(), 3);
        assert_eq!(cmds[0].name, "load");
        assert_eq!(cmds[1].name, "zoom");
        assert_eq!(cmds[2].name, "show");
    }

    #[test]
    fn test_command_with_comment() {
        let cmds = parse_commands("# This is a comment\nload file.pdb").unwrap();
        assert_eq!(cmds.len(), 1);
        assert_eq!(cmds[0].name, "load");
    }

    #[test]
    fn test_util_command() {
        let cmd = parse_command("util.cbag all").unwrap();
        assert_eq!(cmd.name, "util.cbag");
        assert_eq!(cmd.get_str(0), Some("all"));
    }

    #[test]
    fn test_list_argument() {
        let cmd = parse_command("set_view [1.0, 0.0, 0.0]").unwrap();
        if let Some(ArgValue::List(items)) = cmd.get_arg(0) {
            assert_eq!(items.len(), 3);
        } else {
            panic!("Expected list argument");
        }
    }

    #[test]
    fn test_selection_with_not_operator() {
        // "not" should not be partially matched as "no" (boolean)
        let cmd = parse_command("select nca, not (chain A)").unwrap();
        assert_eq!(cmd.name, "select");
        assert_eq!(cmd.get_str(0), Some("nca"));
        assert_eq!(cmd.get_str(1), Some("not (chain A)"));
    }

    #[test]
    fn test_selection_with_parens_and_or() {
        let cmd = parse_command("select sele, (chain A and resi 29) or (chain A and resi 31)").unwrap();
        assert_eq!(cmd.name, "select");
        assert_eq!(cmd.get_str(0), Some("sele"));
        assert_eq!(cmd.get_str(1), Some("(chain A and resi 29) or (chain A and resi 31)"));
    }

    #[test]
    fn test_bool_word_boundary() {
        // "no" as standalone boolean should still work
        let cmd = parse_command("set valence, no").unwrap();
        assert_eq!(cmd.name, "set");
        assert_eq!(cmd.get_bool(1), Some(false));

        // "on" as standalone boolean should still work
        let cmd = parse_command("set valence, on").unwrap();
        assert_eq!(cmd.get_bool(1), Some(true));

        // "notable" should not match "no" - it's an identifier
        let cmd = parse_command("select notable").unwrap();
        assert_eq!(cmd.get_str(0), Some("notable"));
    }

    #[test]
    fn test_line_continuation_with_newline() {
        // Backslash followed by newline should join lines
        let cmd = parse_command("set_view (\\\n  1.0, 2.0, 3.0)").unwrap();
        assert_eq!(cmd.name, "set_view");
        // The argument should be the parenthesized expression
        assert!(cmd.get_str(0).is_some());
    }

    #[test]
    fn test_line_continuation_collapsed() {
        // Backslash followed by space (collapsed continuation from copy-paste)
        let cmd = parse_command("set_view (\\ 1.0, 2.0,\\ 3.0)").unwrap();
        assert_eq!(cmd.name, "set_view");
        let arg = cmd.get_str(0).unwrap();
        // The backslashes should be stripped, leaving just the values
        assert!(arg.contains("1.0"));
        assert!(arg.contains("2.0"));
        assert!(arg.contains("3.0"));
        // Backslashes should not be present
        assert!(!arg.contains('\\'));
    }

    #[test]
    fn test_set_view_with_18_values() {
        // Full set_view command with collapsed line continuations (like from test.pml)
        let cmd = parse_command(
            "set_view (\\ 0.381124, -0.381788, -0.842011,\\ -0.911104, -0.309719, -0.271964,\\ -0.156954, 0.870811, -0.465889,\\ 0.0, 0.0, 639.042053,\\ -1.574999, -16.596998, 9.544001,\\ 561.162231, 794.801453, 14.0 )"
        ).unwrap();
        assert_eq!(cmd.name, "set_view");
        let arg = cmd.get_str(0).unwrap();
        // Should have all 18 values parseable (no backslashes blocking them)
        let s = arg.trim_matches(|c| c == '(' || c == ')');
        let values: Vec<f32> = s
            .split(',')
            .filter_map(|v| v.trim().parse::<f32>().ok())
            .collect();
        assert_eq!(values.len(), 18, "Expected 18 values, got {}: {:?}", values.len(), values);
    }

    #[test]
    fn test_multiline_commands() {
        // Multiple commands with line continuation
        let cmds = parse_commands("set_view (\\\n  1.0, 2.0, 3.0)\nzoom").unwrap();
        assert_eq!(cmds.len(), 2);
        assert_eq!(cmds[0].name, "set_view");
        assert_eq!(cmds[1].name, "zoom");
    }

    #[test]
    fn test_normalize_continuations() {
        // Test the normalization function directly
        assert_eq!(normalize_continuations("a\\\nb"), "ab");
        assert_eq!(normalize_continuations("a\\ b"), "ab");
        assert_eq!(normalize_continuations("a\\  b"), "ab");
        assert_eq!(normalize_continuations("a\\\n  b"), "ab");
        // Regular backslash (not followed by newline or space) should be preserved
        assert_eq!(normalize_continuations("a\\xb"), "a\\xb");
    }

    #[test]
    fn test_set_view_multiline_exact_format() {
        // Exact format from test.pml with actual newlines
        let input = r#"set_view (\ 0.381124,    -0.381788,    -0.842011,\
           -0.911104,    -0.309719,    -0.271964,\
           -0.156954,     0.870811,    -0.465889,\
            0.000000,     0.000000,   639.042053,\
           -1.574999,   -16.596998,     9.544001,\
            561.162231,   794.801453,   14.000000 )"#;

        let cmd = parse_command(input).unwrap();
        assert_eq!(cmd.name, "set_view");
        let arg = cmd.get_str(0).unwrap();

        // Should have all 18 values parseable
        let s = arg.trim_matches(|c| c == '(' || c == ')');
        let values: Vec<f32> = s
            .split(',')
            .filter_map(|v| v.trim().parse::<f32>().ok())
            .collect();
        assert_eq!(values.len(), 18, "Expected 18 values, got {}: {:?}", values.len(), values);
    }

    #[test]
    fn test_parse_commands_full_test_pml() {
        // Full test.pml content
        let input = r#"load ~/Downloads/1IGT.cif
as sticks
set_view (\ 0.381124,    -0.381788,    -0.842011,\
           -0.911104,    -0.309719,    -0.271964,\
           -0.156954,     0.870811,    -0.465889,\
            0.000000,     0.000000,   639.042053,\
           -1.574999,   -16.596998,     9.544001,\
            561.162231,   794.801453,   14.000000 )
ray 1920, 1080, filename=test.png"#;

        let cmds = parse_commands(input).unwrap();
        assert_eq!(cmds.len(), 4, "Expected 4 commands");
        assert_eq!(cmds[0].name, "load");
        assert_eq!(cmds[1].name, "as");
        assert_eq!(cmds[2].name, "set_view");
        assert_eq!(cmds[3].name, "ray");
    }

    // ========================================================================
    // join_continued_lines tests
    // ========================================================================

    #[test]
    fn test_join_continued_lines_simple() {
        let input = "line1\\\nline2";
        let result = join_continued_lines(input);
        assert_eq!(result.trim(), "line1 line2");
    }

    #[test]
    fn test_join_continued_lines_multiple() {
        let input = "set_view (\\\n  1.0, 2.0,\\\n  3.0)";
        let result = join_continued_lines(input);
        // Should join all three lines
        assert!(result.contains("set_view"));
        assert!(result.contains("1.0"));
        assert!(result.contains("3.0"));
        // Should only have one line (newlines removed)
        assert_eq!(result.lines().count(), 1);
    }

    #[test]
    fn test_join_continued_lines_no_continuation() {
        let input = "line1\nline2\nline3";
        let result = join_continued_lines(input);
        assert_eq!(result.lines().count(), 3);
    }

    #[test]
    fn test_join_continued_lines_preserves_content() {
        let input = "load file.pdb\nzoom\nshow cartoon";
        let result = join_continued_lines(input);
        let lines: Vec<&str> = result.lines().collect();
        assert_eq!(lines.len(), 3);
        assert_eq!(lines[0], "load file.pdb");
        assert_eq!(lines[1], "zoom");
        assert_eq!(lines[2], "show cartoon");
    }

    #[test]
    fn test_join_continued_lines_test_pml_format() {
        // Exact format from test.pml
        let input = r#"load ~/Downloads/1IGT.cif
as sticks
set_view (\ 0.381124,    -0.381788,    -0.842011,\
           -0.911104,    -0.309719,    -0.271964,\
           -0.156954,     0.870811,    -0.465889,\
            0.000000,     0.000000,   639.042053,\
           -1.574999,   -16.596998,     9.544001,\
            561.162231,   794.801453,   14.000000 )
ray 1920, 1080, filename=test.png"#;

        let result = join_continued_lines(input);
        let lines: Vec<&str> = result.lines().collect();

        // Should have 4 lines: load, as, set_view (joined), ray
        assert_eq!(lines.len(), 4, "Expected 4 lines, got {}", lines.len());
        assert!(lines[0].starts_with("load"));
        assert!(lines[1].starts_with("as"));
        assert!(lines[2].starts_with("set_view"));
        assert!(lines[3].starts_with("ray"));

        // The set_view line should contain all values
        assert!(lines[2].contains("0.381124"));
        assert!(lines[2].contains("14.000000"));
    }
}
