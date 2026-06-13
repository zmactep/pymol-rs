//! IPC Message Handler
//!
//! Implements the plugin `MessageHandler` trait, bridging the IPC server
//! with the host application through the `PollContext` polling API.

use patinae_plugin::prelude::*;
use patinae_plugin::wire::{WireHostQuery, WireHostQueryValue};

use crate::protocol::{IpcRequest, IpcResponse, OutputKind};
use crate::server::IpcServer;

/// IPC message handler — bridges the IPC server with the plugin system.
///
/// Uses `PollContext` for deferred command execution, dynamic command
/// registration, and reading application state.
pub struct IpcMessageHandler {
    server: IpcServer,
}

impl IpcMessageHandler {
    pub fn new(server: IpcServer) -> Self {
        Self { server }
    }

    /// Handle a single IPC request, returning an optional immediate response.
    ///
    /// Some requests (Execute, RegisterCommand, UnregisterCommand) are deferred
    /// through `PollContext` and have no immediate response. Others (Ping,
    /// GetNames, etc.) can be answered synchronously.
    fn handle_request(
        &mut self,
        request: &IpcRequest,
        ctx: &mut PollContext<'_>,
    ) -> Option<IpcResponse> {
        match request {
            IpcRequest::Execute {
                id,
                command,
                silent,
            } => {
                log::debug!("IPC Execute: {}", command);
                ctx.execute_command(*id, command, *silent);
                // Response sent later when result arrives
                None
            }

            IpcRequest::RegisterCommand {
                name,
                description,
                usage,
                arguments,
            } => {
                log::info!("IPC RegisterCommand: {}", name);
                ctx.register_dynamic_command(
                    name.clone(),
                    description.clone().unwrap_or_default(),
                    usage.clone().unwrap_or_default(),
                    arguments.clone().unwrap_or_default(),
                );
                None
            }

            IpcRequest::UnregisterCommand { name } => {
                log::info!("IPC UnregisterCommand: {}", name);
                ctx.unregister_dynamic_command(name);
                None
            }

            IpcRequest::CallbackResponse {
                id: _,
                success,
                error,
                output,
            } => {
                log::debug!("IPC CallbackResponse: success={}", success);
                // Display output via message bus
                for msg in output {
                    match msg.kind {
                        OutputKind::Info => ctx.bus.print_info(msg.text.clone()),
                        OutputKind::Warning => ctx.bus.print_warning(msg.text.clone()),
                        OutputKind::Error => ctx.bus.print_error(msg.text.clone()),
                        OutputKind::Timing => ctx.bus.print_info(msg.text.clone()),
                    }
                }
                if !success {
                    if let Some(err) = error {
                        ctx.bus.print_error(err.clone());
                    }
                }
                ctx.bus.request_redraw();
                None
            }

            IpcRequest::GetState { id } => {
                // TODO: Implement state serialization
                Some(IpcResponse::Value {
                    id: *id,
                    value: serde_json::json!({}),
                })
            }

            IpcRequest::GetNames { id } => Some(IpcResponse::Value {
                id: *id,
                value: serde_json::json!(&ctx.poll_shared.object_names),
            }),

            IpcRequest::CountAtoms { id, selection } => {
                ctx.query_host(WireHostQuery::CountAtoms {
                    id: *id,
                    selection: selection.clone(),
                });
                None
            }

            IpcRequest::Hello { client_id } => {
                log::info!("IPC client identified as: {}", client_id);
                self.server.set_client_id(client_id.clone());
                Some(IpcResponse::Ok { id: 0 })
            }

            IpcRequest::Quit => {
                log::info!("IPC Quit received");
                ctx.bus.send(AppMessage::Quit);
                Some(IpcResponse::Closing)
            }

            IpcRequest::Ping { id } => Some(IpcResponse::Pong { id: *id }),

            IpcRequest::ShowWindow { id } => {
                log::info!("IPC ShowWindow received");
                ctx.bus.send(AppMessage::ShowWindow);
                Some(IpcResponse::Ok { id: *id })
            }

            IpcRequest::HideWindow { id } => {
                log::info!("IPC HideWindow received");
                ctx.bus.send(AppMessage::HideWindow);
                Some(IpcResponse::Ok { id: *id })
            }

            IpcRequest::GetView { id } => {
                let view = ctx.poll_shared.camera.current_view();
                let r = &view.rotation;

                // Build the 18-value array:
                // [0-8]: 3x3 rotation matrix (row-major)
                // [9-11]: Camera position
                // [12-14]: Origin
                // [15]: Front clip, [16]: Back clip, [17]: FOV
                let values: Vec<f64> = vec![
                    r.data[0] as f64,
                    r.data[1] as f64,
                    r.data[2] as f64,
                    r.data[4] as f64,
                    r.data[5] as f64,
                    r.data[6] as f64,
                    r.data[8] as f64,
                    r.data[9] as f64,
                    r.data[10] as f64,
                    view.position.x as f64,
                    view.position.y as f64,
                    view.position.z as f64,
                    view.origin.x as f64,
                    view.origin.y as f64,
                    view.origin.z as f64,
                    view.clip_front as f64,
                    view.clip_back as f64,
                    view.fov as f64,
                ];

                Some(IpcResponse::Value {
                    id: *id,
                    value: serde_json::json!(values),
                })
            }
        }
    }
}

impl MessageHandler for IpcMessageHandler {
    fn on_message(&mut self, _msg: &AppMessage, _bus: &mut MessageBus) {
        // IPC handler doesn't react to broadcast messages
    }

    fn needs_poll(&self) -> bool {
        true
    }

    fn poll(&mut self, ctx: &mut PollContext<'_>) {
        for result in ctx.host_query_results {
            if let Ok(WireHostQueryValue::CountAtoms(count)) = &result.result {
                if let Err(e) = self.server.send(IpcResponse::Value {
                    id: result.id,
                    value: serde_json::json!(count),
                }) {
                    log::error!("Failed to send query result: {}", e);
                }
            }
        }

        // Deliver results from previously queued command executions
        for result in ctx.command_results {
            let response = match &result.result {
                Ok(()) => IpcResponse::Ok { id: result.id },
                Err(message) => IpcResponse::Error {
                    id: result.id,
                    message: message.clone(),
                },
            };
            if let Err(e) = self.server.send(response) {
                log::error!("Failed to send execution result: {}", e);
            }
        }

        // Forward dynamic command invocations to client as CallbackRequests
        for invocation in ctx.dynamic_invocations {
            let id = self.server.next_callback_id();
            let response = IpcResponse::CallbackRequest {
                id,
                name: invocation.name.clone(),
                args: invocation.args.clone(),
            };
            if let Err(e) = self.server.send(response) {
                log::error!("Failed to send callback request: {}", e);
            }
        }

        // Process incoming IPC requests
        while let Some(request) = self.server.poll() {
            if let Some(response) = self.handle_request(&request, ctx) {
                if let Err(e) = self.server.send(response) {
                    log::error!("Failed to send IPC response: {}", e);
                }
            }
        }
    }
}
