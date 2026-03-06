//! Python Message Handler
//!
//! Implements `MessageHandler` with `needs_poll() = true` for:
//! 1. Updating shared molecule/name snapshots from `SharedContext`
//! 2. Setting `sys._pymolrs_backend` on first poll
//! 3. Draining the command queue → `ctx.execute_command()`
//! 4. Syncing viewport image between Python and the host

use std::sync::{Arc, Mutex};

use pyo3::prelude::*;
use pymol_plugin::prelude::*;

use crate::backend::{PluginBackend, SharedStateHandle};
use crate::engine::PythonEngine;

/// Message handler for the Python plugin.
///
/// Polls each frame to synchronize state and drain queued commands.
pub struct PythonHandler {
    engine: Arc<Mutex<PythonEngine>>,
    shared: SharedStateHandle,
    backend_installed: bool,
    next_cmd_id: u64,
}

impl PythonHandler {
    pub fn new(engine: Arc<Mutex<PythonEngine>>, shared: SharedStateHandle) -> Self {
        Self {
            engine,
            shared,
            backend_installed: false,
            next_cmd_id: 1,
        }
    }

    /// Install the PluginBackend into `sys._pymolrs_backend`.
    fn install_backend(&mut self) {
        // Ensure Python is initialized (configures PYTHONHOME etc.)
        // before any Python::attach call.
        {
            let mut engine = self.engine.lock().unwrap();
            if let Err(e) = engine.ensure_init() {
                log::warn!("Python plugin: failed to initialize: {}", e);
                self.backend_installed = true; // don't retry
                return;
            }
        }

        let shared = self.shared.clone();

        let result = Python::attach(|py| -> Result<(), String> {
            let backend = PluginBackend::new(shared);
            let py_backend = Py::new(py, backend).map_err(|e| e.to_string())?;

            let mut engine = self.engine.lock().unwrap();
            engine.set_backend(py_backend.into_any())
        });

        match result {
            Ok(()) => {
                self.backend_installed = true;
                log::info!("Python plugin: backend installed into sys._pymolrs_backend");

                // Auto-import cmd into __main__ so it's available in the REPL
                let mut engine = self.engine.lock().unwrap();
                if let Err(e) = engine.eval("from pymol_rs import cmd") {
                    log::warn!("Python plugin: failed to auto-import cmd: {}", e);
                } else {
                    log::info!("Python plugin: cmd auto-imported into global namespace");
                }
            }
            Err(e) => {
                log::warn!("Python plugin: failed to install backend: {}", e);
            }
        }
    }

    /// Sync molecule snapshots and viewport image from SharedContext.
    fn update_snapshots(&self, shared_ctx: &SharedContext<'_>) {
        let mut state = self.shared.lock().unwrap();

        // Update names
        state.names = shared_ctx.registry.names().map(|s| s.to_string()).collect();

        // Update molecule snapshots
        state.molecules.clear();
        for name in shared_ctx.registry.names() {
            if let Some(mol_obj) = shared_ctx.registry.get_molecule(name) {
                state.molecules.push((name.to_string(), mol_obj.molecule().clone()));
            }
        }

        // Update viewport image snapshot
        state.viewport_image = shared_ctx.viewport_image.map(|img| {
            (img.data.clone(), img.width, img.height)
        });
    }

    /// Drain the command queue and schedule execution.
    fn drain_commands(&mut self, ctx: &mut PollContext<'_>) {
        let commands: Vec<(String, bool)> = {
            let mut state = self.shared.lock().unwrap();
            std::mem::take(&mut state.cmd_queue)
        };

        for (command, silent) in commands {
            let id = self.next_cmd_id;
            self.next_cmd_id += 1;
            ctx.execute_command(id, &command, silent);
        }
    }

    /// Drain queued viewport image set/clear and send via message bus.
    fn drain_image_queue(&self, ctx: &mut PollContext<'_>) {
        let pending = {
            let mut state = self.shared.lock().unwrap();
            state.set_image_queue.take()
        };

        if let Some(action) = pending {
            match action {
                Some((data, width, height)) => {
                    ctx.bus.set_viewport_image(data, width, height);
                }
                None => {
                    ctx.bus.clear_viewport_image();
                }
            }
        }
    }
}

impl MessageHandler for PythonHandler {
    fn on_message(&mut self, _msg: &AppMessage, _bus: &mut MessageBus) {
        // Python handler doesn't react to broadcast messages
    }

    fn needs_poll(&self) -> bool {
        true
    }

    fn poll(&mut self, ctx: &mut PollContext<'_>) {
        // Install backend on first poll
        if !self.backend_installed {
            self.install_backend();
        }

        // Sync state snapshots
        self.update_snapshots(ctx.shared);

        // Drain queued commands
        self.drain_commands(ctx);

        // Drain queued viewport image changes
        self.drain_image_queue(ctx);
    }
}
