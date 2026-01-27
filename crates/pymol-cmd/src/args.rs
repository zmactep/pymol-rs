//! Command argument types and utilities
//!
//! This module defines the types used to represent parsed command arguments.

use std::fmt;

/// A command argument value
///
/// Arguments can be strings, numbers, booleans, lists, or None (for omitted optional args).
#[derive(Debug, Clone, PartialEq)]
pub enum ArgValue {
    /// String value (may be a selection, filename, object name, etc.)
    String(String),
    /// Integer value
    Int(i64),
    /// Floating-point value
    Float(f64),
    /// Boolean value (on/off, true/false, 1/0)
    Bool(bool),
    /// List of values (e.g., coordinates [1.0, 2.0, 3.0])
    List(Vec<ArgValue>),
    /// No value (omitted optional argument)
    None,
}

impl Default for ArgValue {
    fn default() -> Self {
        ArgValue::None
    }
}

impl fmt::Display for ArgValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ArgValue::String(s) => write!(f, "{}", s),
            ArgValue::Int(i) => write!(f, "{}", i),
            ArgValue::Float(n) => write!(f, "{}", n),
            ArgValue::Bool(b) => write!(f, "{}", if *b { "on" } else { "off" }),
            ArgValue::List(items) => {
                write!(f, "[")?;
                for (i, item) in items.iter().enumerate() {
                    if i > 0 {
                        write!(f, ", ")?;
                    }
                    write!(f, "{}", item)?;
                }
                write!(f, "]")
            }
            ArgValue::None => write!(f, ""),
        }
    }
}

impl ArgValue {
    /// Check if this value is None
    pub fn is_none(&self) -> bool {
        matches!(self, ArgValue::None)
    }

    /// Check if this value is present (not None)
    pub fn is_some(&self) -> bool {
        !self.is_none()
    }

    /// Try to get as a string
    pub fn as_str(&self) -> Option<&str> {
        match self {
            ArgValue::String(s) => Some(s),
            _ => None,
        }
    }

    /// Try to get as an integer
    pub fn as_int(&self) -> Option<i64> {
        match self {
            ArgValue::Int(i) => Some(*i),
            ArgValue::Float(f) => Some(*f as i64),
            ArgValue::String(s) => s.parse().ok(),
            _ => None,
        }
    }

    /// Try to get as a float
    pub fn as_float(&self) -> Option<f64> {
        match self {
            ArgValue::Float(f) => Some(*f),
            ArgValue::Int(i) => Some(*i as f64),
            ArgValue::String(s) => s.parse().ok(),
            _ => None,
        }
    }

    /// Try to get as a boolean
    pub fn as_bool(&self) -> Option<bool> {
        match self {
            ArgValue::Bool(b) => Some(*b),
            ArgValue::Int(i) => Some(*i != 0),
            ArgValue::String(s) => match s.to_lowercase().as_str() {
                "true" | "on" | "yes" | "1" => Some(true),
                "false" | "off" | "no" | "0" => Some(false),
                _ => None,
            },
            _ => None,
        }
    }

    /// Convert to string, consuming the value
    pub fn into_string(self) -> Option<String> {
        match self {
            ArgValue::String(s) => Some(s),
            ArgValue::Int(i) => Some(i.to_string()),
            ArgValue::Float(f) => Some(f.to_string()),
            ArgValue::Bool(b) => Some(if b { "on".to_string() } else { "off".to_string() }),
            _ => None,
        }
    }
}

impl From<&str> for ArgValue {
    fn from(s: &str) -> Self {
        ArgValue::String(s.to_string())
    }
}

impl From<String> for ArgValue {
    fn from(s: String) -> Self {
        ArgValue::String(s)
    }
}

impl From<i64> for ArgValue {
    fn from(i: i64) -> Self {
        ArgValue::Int(i)
    }
}

impl From<i32> for ArgValue {
    fn from(i: i32) -> Self {
        ArgValue::Int(i as i64)
    }
}

impl From<f64> for ArgValue {
    fn from(f: f64) -> Self {
        ArgValue::Float(f)
    }
}

impl From<f32> for ArgValue {
    fn from(f: f32) -> Self {
        ArgValue::Float(f as f64)
    }
}

impl From<bool> for ArgValue {
    fn from(b: bool) -> Self {
        ArgValue::Bool(b)
    }
}

/// A parsed command with its name and arguments
#[derive(Debug, Clone, PartialEq)]
pub struct ParsedCommand {
    /// The command name (e.g., "load", "zoom", "color")
    pub name: String,
    /// Arguments as (optional_name, value) pairs
    ///
    /// Positional arguments have `None` as the name.
    /// Named arguments have `Some(name)` as the name.
    pub args: Vec<(Option<String>, ArgValue)>,
}

impl ParsedCommand {
    /// Create a new parsed command with no arguments
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            args: Vec::new(),
        }
    }

    /// Add a positional argument
    pub fn with_arg(mut self, value: impl Into<ArgValue>) -> Self {
        self.args.push((None, value.into()));
        self
    }

    /// Add a named argument
    pub fn with_named_arg(mut self, name: impl Into<String>, value: impl Into<ArgValue>) -> Self {
        self.args.push((Some(name.into()), value.into()));
        self
    }

    /// Get the number of arguments
    pub fn arg_count(&self) -> usize {
        self.args.len()
    }

    /// Get a positional argument by index (0-based)
    pub fn get_arg(&self, index: usize) -> Option<&ArgValue> {
        self.args.get(index).map(|(_, v)| v)
    }

    /// Get a named argument by name
    pub fn get_named(&self, name: &str) -> Option<&ArgValue> {
        self.args
            .iter()
            .find(|(n, _)| n.as_deref() == Some(name))
            .map(|(_, v)| v)
    }

    /// Get positional argument as string
    pub fn get_str(&self, index: usize) -> Option<&str> {
        self.get_arg(index).and_then(|v| v.as_str())
    }

    /// Get positional argument as string, with default
    pub fn get_str_or<'a>(&'a self, index: usize, default: &'a str) -> &'a str {
        self.get_str(index).unwrap_or(default)
    }

    /// Get positional argument as int
    pub fn get_int(&self, index: usize) -> Option<i64> {
        self.get_arg(index).and_then(|v| v.as_int())
    }

    /// Get positional argument as int, with default
    pub fn get_int_or(&self, index: usize, default: i64) -> i64 {
        self.get_int(index).unwrap_or(default)
    }

    /// Get positional argument as float
    pub fn get_float(&self, index: usize) -> Option<f64> {
        self.get_arg(index).and_then(|v| v.as_float())
    }

    /// Get positional argument as float, with default
    pub fn get_float_or(&self, index: usize, default: f64) -> f64 {
        self.get_float(index).unwrap_or(default)
    }

    /// Get positional argument as bool
    pub fn get_bool(&self, index: usize) -> Option<bool> {
        self.get_arg(index).and_then(|v| v.as_bool())
    }

    /// Get positional argument as bool, with default
    pub fn get_bool_or(&self, index: usize, default: bool) -> bool {
        self.get_bool(index).unwrap_or(default)
    }

    /// Get named argument as string
    pub fn get_named_str(&self, name: &str) -> Option<&str> {
        self.get_named(name).and_then(|v| v.as_str())
    }

    /// Get named argument as string, with default
    pub fn get_named_str_or<'a>(&'a self, name: &str, default: &'a str) -> &'a str {
        self.get_named_str(name).unwrap_or(default)
    }

    /// Get named argument as int
    pub fn get_named_int(&self, name: &str) -> Option<i64> {
        self.get_named(name).and_then(|v| v.as_int())
    }

    /// Get named argument as int, with default
    pub fn get_named_int_or(&self, name: &str, default: i64) -> i64 {
        self.get_named_int(name).unwrap_or(default)
    }

    /// Get named argument as float
    pub fn get_named_float(&self, name: &str) -> Option<f64> {
        self.get_named(name).and_then(|v| v.as_float())
    }

    /// Get named argument as float, with default
    pub fn get_named_float_or(&self, name: &str, default: f64) -> f64 {
        self.get_named_float(name).unwrap_or(default)
    }

    /// Get named argument as bool
    pub fn get_named_bool(&self, name: &str) -> Option<bool> {
        self.get_named(name).and_then(|v| v.as_bool())
    }

    /// Get named argument as bool, with default
    pub fn get_named_bool_or(&self, name: &str, default: bool) -> bool {
        self.get_named_bool(name).unwrap_or(default)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_arg_value_conversions() {
        let s = ArgValue::String("hello".to_string());
        assert_eq!(s.as_str(), Some("hello"));

        let i = ArgValue::Int(42);
        assert_eq!(i.as_int(), Some(42));
        assert_eq!(i.as_float(), Some(42.0));

        let f = ArgValue::Float(3.14);
        assert_eq!(f.as_float(), Some(3.14));
        assert_eq!(f.as_int(), Some(3));

        let b = ArgValue::Bool(true);
        assert_eq!(b.as_bool(), Some(true));
    }

    #[test]
    fn test_string_to_bool() {
        assert_eq!(ArgValue::String("true".to_string()).as_bool(), Some(true));
        assert_eq!(ArgValue::String("on".to_string()).as_bool(), Some(true));
        assert_eq!(ArgValue::String("false".to_string()).as_bool(), Some(false));
        assert_eq!(ArgValue::String("off".to_string()).as_bool(), Some(false));
        assert_eq!(ArgValue::String("invalid".to_string()).as_bool(), None);
    }

    #[test]
    fn test_parsed_command() {
        let cmd = ParsedCommand::new("load")
            .with_arg("protein.pdb")
            .with_named_arg("object", "mol");

        assert_eq!(cmd.name, "load");
        assert_eq!(cmd.get_str(0), Some("protein.pdb"));
        assert_eq!(cmd.get_named_str("object"), Some("mol"));
        assert_eq!(cmd.get_named_str("missing"), None);
    }
}
