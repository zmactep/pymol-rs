//! PyMOL-RS Python Plugin
//!
//! Embeds a Python interpreter in the GUI application, providing:
//! - `python` command (alias `/`) for inline Python execution
//! - `.py` file handler for `run script.py`
//! - `sys._pymolrs_backend` for embedded `pymol_rs` package access

mod backend;
mod commands;
mod editor;
mod engine;
mod handler;
mod highlight;
mod widgets;
mod worker;

use std::sync::{Arc, Mutex};

use pymol_plugin::prelude::PanelConfig;
use pymol_plugin::pymol_plugin;

use backend::SharedState;
use commands::{AlterBuffer, AlterCommand, IterateCommand, PythonCommand};
use editor::PythonEditorComponent;
use handler::PythonHandler;
use worker::{EditorOutputHandle, WorkItem, WorkOrigin};

pymol_plugin! {
    name: "python",
    description: "Embedded Python interpreter for scripting and automation",
    commands: [],
    register: |reg| {
        // Alter buffer shared between PluginBackend (producer) and handler (consumer)
        let alter_buffer: AlterBuffer = Arc::new(Mutex::new(Vec::new()));

        // Shared state for host ↔ Python communication
        let shared_state = Arc::new(Mutex::new(SharedState::new(alter_buffer)));

        // Spawn the Python worker thread
        let (worker_handle, result_rx) = worker::spawn_worker();

        // Shared buffer for routing editor results
        let editor_output: EditorOutputHandle = Arc::new(Mutex::new(Vec::new()));

        // Register commands
        reg.register_command(PythonCommand::new(worker_handle.clone()));
        reg.register_command(IterateCommand::new(worker_handle.clone(), shared_state.clone()));
        reg.register_command(AlterCommand::new(worker_handle.clone(), shared_state.clone()));

        // Register .py script handler (fire-and-forget: output via poll)
        let script_worker = worker_handle.clone();
        reg.register_script_handler("py", move |path: &str| {
            script_worker.submit(WorkItem::ExecFile {
                path: path.to_string(),
                origin: WorkOrigin::Script,
            });
            Ok(())
        });

        // Script editor panel (Top slot, tabbed with REPL)
        let editor = PythonEditorComponent::new(worker_handle.clone(), editor_output.clone());
        reg.register_component(editor, PanelConfig::top(150.0));

        // Message handler for poll-based state sync, result draining, and backend installation
        let handler = PythonHandler::new(worker_handle, shared_state, result_rx, editor_output);
        reg.set_message_handler(handler);
    },
}
