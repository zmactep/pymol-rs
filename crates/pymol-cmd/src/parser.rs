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
    let input = input.trim();
    if input.is_empty() {
        return Ok(Vec::new());
    }

    let mut commands = Vec::new();
    let mut current = input;

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

/// Parse command arguments (comma-separated list)
fn parse_arguments(input: &str) -> IResult<&str, Vec<(Option<String>, ArgValue)>> {
    let (input, first) = parse_argument(input)?;
    let (input, rest) = many0(preceded(
        tuple((multispace0, char(','), multispace0)),
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
        // Try parenthesized expression (selection)
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
fn parse_paren_expr(input: &str) -> IResult<&str, &str> {
    recognize(delimited(
        char('('),
        take_balanced_parens,
        char(')'),
    ))(input)
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
}
