use std::sync::{atomic::AtomicBool, Arc, Mutex};

use patinae_plugin::prelude::PluginRegistrar;

use crate::atom_ops::AlterBuffer;
use crate::commands::{AlterCommand, IterateCommand, PythonCommand};
use crate::handler::PythonHandler;
use crate::panel::PythonScriptPanel;
use crate::shared::SharedState;
use crate::worker::{self, WorkItem, WorkOrigin};

pub(crate) fn register(reg: &mut PluginRegistrar) {
    let alter_buffer: AlterBuffer = Arc::new(Mutex::new(Vec::new()));
    let interrupt_requested = Arc::new(AtomicBool::new(false));

    let shared_state = Arc::new(Mutex::new(SharedState::new(
        alter_buffer,
        interrupt_requested.clone(),
    )));
    let panel_state = crate::panel::shared_panel_state();

    let (worker_handle, result_rx) = worker::spawn_worker(interrupt_requested);

    reg.register_command(PythonCommand::new(worker_handle.clone()));
    reg.register_command(IterateCommand::new(
        worker_handle.clone(),
        shared_state.clone(),
    ));
    reg.register_command(AlterCommand::new(
        worker_handle.clone(),
        shared_state.clone(),
    ));
    reg.register_panel(PythonScriptPanel::new(
        worker_handle.clone(),
        panel_state.clone(),
    ));

    let script_worker = worker_handle.clone();
    reg.register_script_handler("py", move |path: &str| {
        script_worker.submit(WorkItem::ExecFile {
            path: path.to_string(),
            origin: WorkOrigin::Script,
        });
        Ok(())
    });

    reg.set_message_handler(PythonHandler::new(
        worker_handle,
        shared_state,
        result_rx,
        panel_state,
    ));
}
