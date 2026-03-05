//! PyMOL-RS Python Plugin
//!
//! Embeds a Python interpreter in the GUI application, providing:
//! - `python` command (alias `/`) for inline Python execution
//! - `.py` file handler for `run script.py`
//! - `sys._pymolrs_backend` for embedded `pymol_rs` package access

mod backend;
mod commands;
mod engine;
mod handler;

use std::sync::{Arc, Mutex};

use pymol_plugin::pymol_plugin;

use backend::SharedState;
use commands::PythonCommand;
use engine::PythonEngine;
use handler::PythonHandler;

pymol_plugin! {
    name: "python",
    version: "0.2.0",
    description: "Embedded Python interpreter for scripting and automation",
    commands: [],
    register: |reg| {
        // Shared engine and state
        let engine = Arc::new(Mutex::new(PythonEngine::new()));
        let shared_state = Arc::new(Mutex::new(SharedState::new()));

        // Register the `python` command (with `/` alias)
        reg.register_command(PythonCommand::new(engine.clone()));

        // Register .py file handler
        let file_engine = engine.clone();
        reg.register_file_handler("py", move |path: &str| {
            let mut eng = file_engine.lock().unwrap();
            match eng.exec_file(path) {
                Ok(output) => {
                    if !output.is_empty() {
                        log::info!("{}", output);
                    }
                    Ok(())
                }
                Err(e) => Err(e),
            }
        });

        // Message handler for poll-based state sync and backend installation
        let handler = PythonHandler::new(engine, shared_state);
        reg.set_message_handler(handler);
    },
}
