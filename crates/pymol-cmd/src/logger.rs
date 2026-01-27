//! Command logging for replay and scripting
//!
//! Records executed commands to a file for later replay.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::error::CmdResult;

/// Log file format
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum LogFormat {
    /// PyMOL command format (.pml)
    #[default]
    Pml,
    /// Python script format (.py)
    Python,
}

impl LogFormat {
    /// Detect format from file extension
    pub fn from_extension(path: &Path) -> Self {
        match path.extension().and_then(|e| e.to_str()) {
            Some("py") => LogFormat::Python,
            _ => LogFormat::Pml,
        }
    }
}

/// Command logger for recording executed commands
#[derive(Debug)]
pub struct CommandLogger {
    /// Output file writer
    writer: Option<BufWriter<File>>,
    /// Log format
    format: LogFormat,
    /// Whether logging is active
    active: bool,
}

impl Default for CommandLogger {
    fn default() -> Self {
        Self::new()
    }
}

impl CommandLogger {
    /// Create a new inactive logger
    pub fn new() -> Self {
        Self {
            writer: None,
            format: LogFormat::Pml,
            active: false,
        }
    }

    /// Open a log file for writing
    ///
    /// The format is auto-detected from the file extension:
    /// - `.py` -> Python format
    /// - anything else -> PML format
    pub fn log_open(&mut self, path: &Path) -> CmdResult<()> {
        self.log_close()?;

        let format = LogFormat::from_extension(path);
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        // Write header
        match format {
            LogFormat::Pml => {
                writeln!(writer, "# PyMOL-RS command log")?;
            }
            LogFormat::Python => {
                writeln!(writer, "# PyMOL-RS Python command log")?;
                writeln!(writer, "from pymol import cmd")?;
                writeln!(writer)?;
            }
        }

        self.writer = Some(writer);
        self.format = format;
        self.active = true;

        log::info!("Logging to {:?} ({:?} format)", path, format);
        Ok(())
    }

    /// Log a command
    ///
    /// Does nothing if logging is not active.
    pub fn log(&mut self, command: &str) {
        if !self.active {
            return;
        }

        if let Some(ref mut writer) = self.writer {
            let result = match self.format {
                LogFormat::Pml => writeln!(writer, "{}", command),
                LogFormat::Python => writeln!(writer, "cmd.do(\"{}\")", escape_python_string(command)),
            };

            if let Err(e) = result {
                log::warn!("Failed to write to log: {}", e);
            }
        }
    }

    /// Close the log file
    pub fn log_close(&mut self) -> CmdResult<()> {
        if let Some(mut writer) = self.writer.take() {
            writer.flush()?;
        }
        self.active = false;
        Ok(())
    }

    /// Check if logging is active
    pub fn is_active(&self) -> bool {
        self.active
    }
}

impl Drop for CommandLogger {
    fn drop(&mut self) {
        let _ = self.log_close();
    }
}

/// Escape a string for inclusion in a Python string literal
fn escape_python_string(s: &str) -> String {
    let mut result = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '\\' => result.push_str("\\\\"),
            '"' => result.push_str("\\\""),
            '\n' => result.push_str("\\n"),
            '\r' => result.push_str("\\r"),
            '\t' => result.push_str("\\t"),
            _ => result.push(c),
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_detection() {
        assert_eq!(
            LogFormat::from_extension(Path::new("log.pml")),
            LogFormat::Pml
        );
        assert_eq!(
            LogFormat::from_extension(Path::new("log.py")),
            LogFormat::Python
        );
        assert_eq!(
            LogFormat::from_extension(Path::new("log.txt")),
            LogFormat::Pml
        );
    }

    #[test]
    fn test_escape_python_string() {
        assert_eq!(escape_python_string("hello"), "hello");
        assert_eq!(escape_python_string("he\"llo"), "he\\\"llo");
        assert_eq!(escape_python_string("he\\llo"), "he\\\\llo");
        assert_eq!(escape_python_string("he\nllo"), "he\\nllo");
    }
}
