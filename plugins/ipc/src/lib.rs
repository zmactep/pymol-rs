//! PyMOL-RS IPC Plugin
//!
//! Provides a Unix domain socket IPC server for external control of
//! the application (e.g., from pymol-python).
//!
//! Activated by setting the `PYMOL_RS_IPC_SOCKET` environment variable
//! to the desired socket path.

mod handler;
mod protocol;
mod server;

use pymol_plugin::pymol_plugin;

use handler::IpcMessageHandler;
use server::IpcServer;

pymol_plugin! {
    name: "ipc",
    description: "IPC server for external control via Unix domain socket",
    commands: [],
    register: |reg| {
        if let Ok(socket_path) = std::env::var("PYMOL_RS_IPC_SOCKET") {
            let path = std::path::PathBuf::from(&socket_path);
            match IpcServer::bind(&path) {
                Ok(server) => {
                    log::info!("IPC plugin: server bound to {:?}", socket_path);
                    reg.set_message_handler(IpcMessageHandler::new(server));
                }
                Err(e) => {
                    log::warn!("IPC plugin: failed to bind server: {}", e);
                }
            }
        } else {
            log::info!("IPC plugin: PYMOL_RS_IPC_SOCKET not set, server disabled");
        }
    },
}
