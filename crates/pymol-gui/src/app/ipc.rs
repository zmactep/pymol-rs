//! IPC Request Handling
//!
//! Processing incoming IPC requests from external clients (e.g., pymol-python).

use crate::ipc::{
    IpcRequest, IpcResponse, IpcServer,
    OutputKind as IpcOutputKind,
};

use super::App;

impl App {
    /// Process IPC requests
    ///
    /// This should be called each frame to handle incoming IPC requests
    /// from external clients (e.g., pymol-python).
    pub(crate) fn process_ipc(&mut self) {
        // Take the server temporarily to avoid borrow conflicts
        let mut server = match self.ipc_server.take() {
            Some(s) => s,
            None => return,
        };

        // Process all pending requests
        while let Some(request) = server.poll() {
            let response = self.handle_ipc_request(&request, &mut server);
            if let Some(resp) = response {
                if let Err(e) = server.send(resp) {
                    log::error!("Failed to send IPC response: {}", e);
                }
            }
        }

        // Put the server back
        self.ipc_server = Some(server);
    }

    /// Handle a single IPC request
    fn handle_ipc_request(&mut self, request: &IpcRequest, server: &mut IpcServer) -> Option<IpcResponse> {
        match request {
            IpcRequest::Execute { id, command, silent } => {
                log::debug!("IPC Execute: {}", command);
                // Execute the command and propagate errors to the client
                match self.execute_command(command, *silent) {
                    Ok(()) => Some(IpcResponse::Ok { id: *id }),
                    Err(message) => Some(IpcResponse::Error { id: *id, message }),
                }
            }

            IpcRequest::RegisterCommand { name, help } => {
                log::info!("IPC RegisterCommand: {}", name);
                self.external_commands.register(name.clone());
                self.executor.registry_mut().register_external_help(name.clone(), help.clone());
                None
            }

            IpcRequest::UnregisterCommand { name } => {
                log::info!("IPC UnregisterCommand: {}", name);
                self.external_commands.unregister(name);
                self.executor.registry_mut().remove_external_help(name);
                None
            }

            IpcRequest::CallbackResponse { id, success, error, output } => {
                log::debug!("IPC CallbackResponse: id={}, success={}", id, success);
                // This is handled by the async task system via the server's pending callbacks
                let _ = server.handle_callback_response(request);

                // Also display output immediately
                for msg in output {
                    match msg.kind {
                        IpcOutputKind::Info => self.output.print_info(msg.text.clone()),
                        IpcOutputKind::Warning => self.output.print_warning(msg.text.clone()),
                        IpcOutputKind::Error => self.output.print_error(msg.text.clone()),
                    }
                }

                if !success {
                    if let Some(err) = error {
                        self.output.print_error(err.clone());
                    }
                }

                self.needs_redraw = true;
                None // No response needed
            }

            IpcRequest::GetState { id } => {
                // TODO: Implement state serialization
                Some(IpcResponse::Value {
                    id: *id,
                    value: serde_json::json!({})
                })
            }

            IpcRequest::GetNames { id } => {
                let names: Vec<String> = self.state.registry.names()
                    .map(|s| s.to_string())
                    .collect();
                Some(IpcResponse::Value {
                    id: *id,
                    value: serde_json::json!(names)
                })
            }

            IpcRequest::CountAtoms { id, selection } => {
                let mut count = 0;
                for name in self.state.registry.names() {
                    if let Some(mol_obj) = self.state.registry.get_molecule(name) {
                        if let Ok(result) = pymol_select::select(mol_obj.molecule(), selection) {
                            count += result.count();
                        }
                    }
                }
                Some(IpcResponse::Value {
                    id: *id,
                    value: serde_json::json!(count)
                })
            }

            IpcRequest::Hello { client_id } => {
                log::info!("IPC client identified as: {}", client_id);
                server.set_client_id(client_id.clone());
                Some(IpcResponse::Ok { id: 0 })
            }

            IpcRequest::Quit => {
                log::info!("IPC Quit received");
                self.quit_requested = true;
                Some(IpcResponse::Closing)
            }

            IpcRequest::Ping { id } => {
                Some(IpcResponse::Pong { id: *id })
            }

            IpcRequest::ShowWindow { id } => {
                log::info!("IPC ShowWindow received");
                self.show_window();
                self.headless = false;
                self.needs_redraw = true;
                Some(IpcResponse::Ok { id: *id })
            }

            IpcRequest::HideWindow { id } => {
                log::info!("IPC HideWindow received");
                self.hide_window();
                self.headless = true;
                Some(IpcResponse::Ok { id: *id })
            }

            IpcRequest::GetView { id } => {
                // Get current view as 18 floats
                let view = self.state.camera.current_view();
                let r = &view.rotation;

                // Build the 18-value array:
                // [0-8]: 3x3 rotation matrix (row-major)
                // [9-11]: Camera position
                // [12-14]: Origin
                // [15]: Front clip, [16]: Back clip, [17]: FOV
                let values: Vec<f64> = vec![
                    r.data[0] as f64, r.data[1] as f64, r.data[2] as f64,
                    r.data[4] as f64, r.data[5] as f64, r.data[6] as f64,
                    r.data[8] as f64, r.data[9] as f64, r.data[10] as f64,
                    view.position.x as f64, view.position.y as f64, view.position.z as f64,
                    view.origin.x as f64, view.origin.y as f64, view.origin.z as f64,
                    view.clip_front as f64, view.clip_back as f64, view.fov as f64,
                ];

                Some(IpcResponse::Value {
                    id: *id,
                    value: serde_json::json!(values),
                })
            }
        }
    }
}
