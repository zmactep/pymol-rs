//! Python Message Handler
//!
//! Implements `MessageHandler` with `needs_poll() = true` for:
//! 1. Draining results from the Python worker thread
//! 2. Setting notification overlay while Python is busy
//! 3. Updating shared molecule/name snapshots from `SharedContext`
//! 4. Installing `sys._pymolrs_backend` on first poll (via worker)
//! 5. Draining the command queue → `ctx.execute_command()`
//! 6. Syncing viewport image between Python and the host

use std::sync::mpsc::Receiver;

use pymol_plugin::prelude::*;

use crate::backend::SharedStateHandle;
use crate::worker::{
    EditorOutputHandle, OutputEntry, WorkItem, WorkOrigin, WorkResult, WorkerHandle,
};

/// Message handler for the Python plugin.
///
/// Polls each frame to drain worker results, synchronize state,
/// and process queued commands.
pub struct PythonHandler {
    worker: WorkerHandle,
    shared: SharedStateHandle,
    result_rx: Receiver<WorkResult>,
    editor_output: EditorOutputHandle,
    backend_requested: bool,
    next_cmd_id: u64,
}

impl PythonHandler {
    pub fn new(
        worker: WorkerHandle,
        shared: SharedStateHandle,
        result_rx: Receiver<WorkResult>,
        editor_output: EditorOutputHandle,
    ) -> Self {
        Self {
            worker,
            shared,
            result_rx,
            editor_output,
            backend_requested: false,
            next_cmd_id: 1,
        }
    }

    /// Drain results from the Python worker thread and route by origin.
    fn drain_python_results(&mut self, ctx: &mut PollContext<'_>) {
        while let Ok(result) = self.result_rx.try_recv() {
            match result.origin {
                WorkOrigin::Command | WorkOrigin::Script => match &result.result {
                    Ok(output) if !output.is_empty() => {
                        ctx.bus.print_info(output);
                    }
                    Err(err) => {
                        ctx.bus.print_error(err);
                    }
                    _ => {}
                },
                WorkOrigin::Editor => {
                    let entry = match result.result {
                        Ok(output) if !output.is_empty() => Some(OutputEntry {
                            text: output,
                            is_error: false,
                        }),
                        Err(err) => Some(OutputEntry {
                            text: err,
                            is_error: true,
                        }),
                        _ => None,
                    };
                    if let Some(entry) = entry {
                        if let Ok(mut buf) = self.editor_output.lock() {
                            buf.push(entry);
                        }
                    }
                }
                WorkOrigin::Setup => match &result.result {
                    Ok(msg) => log::info!("Python plugin: {}", msg),
                    Err(e) => log::warn!("Python plugin: setup failed: {}", e),
                },
            }
        }
    }

    /// Sync molecule snapshots and viewport image from SharedContext.
    fn update_snapshots(&self, shared_ctx: &SharedContext<'_>) {
        let mut state = self.shared.lock().unwrap();

        // Update names
        state.names = shared_ctx
            .registry
            .names()
            .map(|s| s.to_string())
            .collect();

        // Update molecule snapshots
        state.molecules.clear();
        for name in shared_ctx.registry.names() {
            if let Some(mol_obj) = shared_ctx.registry.get_molecule(name) {
                state
                    .molecules
                    .push((name.to_string(), mol_obj.molecule().clone()));
            }
        }

        // Update viewport image snapshot
        state.viewport_image = shared_ctx
            .viewport_image
            .map(|img| (img.data.clone(), img.width, img.height));
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
        // Request backend installation on first poll (async via worker)
        if !self.backend_requested {
            self.worker.submit(WorkItem::InstallBackend {
                shared: self.shared.clone(),
            });
            self.backend_requested = true;
        }

        // Drain results from the Python worker thread
        self.drain_python_results(ctx);

        // Show notification overlay while Python is busy
        if self.worker.is_busy() {
            ctx.set_notification("Running Python...");
        }

        // Sync state snapshots
        self.update_snapshots(ctx.shared);

        // Drain queued commands from Python
        self.drain_commands(ctx);

        // Drain queued viewport image changes
        self.drain_image_queue(ctx);
    }
}
