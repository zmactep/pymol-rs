//! Python Command
//!
//! Implements the `python` command (with `/` alias) for executing
//! Python code from the PyMOL-RS command line.

use std::sync::{Arc, Mutex};

use pymol_plugin::prelude::*;

use crate::engine::PythonEngine;

/// The `python` command — executes Python code inline.
///
/// Usage:
///   python print("hello")
///   /import math; print(math.pi)
pub struct PythonCommand {
    pub(crate) engine: Arc<Mutex<PythonEngine>>,
}

impl PythonCommand {
    pub fn new(engine: Arc<Mutex<PythonEngine>>) -> Self {
        Self { engine }
    }
}

impl Command for PythonCommand {
    fn name(&self) -> &str {
        "python"
    }

    fn aliases(&self) -> &[&str] {
        &["/"]
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        // Reconstruct the code from all positional args.
        // The parser splits on commas, so we rejoin them to reconstruct
        // the original Python expression (e.g., "print(1, 2)").
        let code: String = args
            .args
            .iter()
            .filter_map(|(name, val)| {
                if name.is_none() {
                    val.to_string_repr()
                } else {
                    None
                }
            })
            .collect::<Vec<_>>()
            .join(", ");
        let code = code.trim();
        if code.is_empty() {
            ctx.print("Usage: python <code>  or  /<code>");
            return Ok(());
        }

        let mut engine = self.engine.lock().unwrap();
        match engine.eval(code) {
            Ok(output) => {
                if !output.is_empty() {
                    for line in output.lines() {
                        ctx.print(line);
                    }
                }
                Ok(())
            }
            Err(err) => {
                ctx.print_error(&err);
                Err(CmdError::Execution(err))
            }
        }
    }

    fn help(&self) -> &str {
        "python <code>\n/code\n\n\
         Execute Python code.\n\n\
         Examples:\n\
             python print(\"hello\")\n\
             /import math; print(math.pi)\n\
             python from pymol_rs import cmd; cmd.color(\"red\", \"all\")"
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::None]
    }
}
