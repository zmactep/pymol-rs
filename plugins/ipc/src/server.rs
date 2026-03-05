//! IPC Server implementation
//!
//! Provides a non-blocking IPC server using Unix domain sockets.
//! Simplified from the original pymol-gui implementation: no tokio deps,
//! no callback channels — all async handling goes through the plugin's
//! `PollContext` instead.

use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use super::protocol::{IpcRequest, IpcResponse};

/// IPC Server for external control
///
/// Listens on a Unix domain socket and processes incoming requests.
/// The server is non-blocking and should be polled each frame.
pub struct IpcServer {
    /// Path to the socket file
    socket_path: PathBuf,
    /// Listener for incoming connections
    #[cfg(unix)]
    listener: std::os::unix::net::UnixListener,
    /// Currently connected client stream
    #[cfg(unix)]
    client: Option<ClientConnection>,
    /// Next callback ID counter
    next_callback_id: AtomicU64,
}

#[cfg(unix)]
struct ClientConnection {
    stream: std::os::unix::net::UnixStream,
    reader: BufReader<std::os::unix::net::UnixStream>,
    client_id: String,
}

impl IpcServer {
    /// Create and bind a new IPC server
    #[cfg(unix)]
    pub fn bind(socket_path: &Path) -> std::io::Result<Self> {
        use std::os::unix::net::UnixListener;

        // Remove existing socket file if present
        let _ = std::fs::remove_file(socket_path);

        let listener = UnixListener::bind(socket_path)?;
        listener.set_nonblocking(true)?;

        log::info!("IPC server listening on {:?}", socket_path);

        Ok(Self {
            socket_path: socket_path.to_path_buf(),
            listener,
            client: None,
            next_callback_id: AtomicU64::new(1),
        })
    }

    #[cfg(not(unix))]
    pub fn bind(socket_path: &Path) -> std::io::Result<Self> {
        Err(std::io::Error::new(
            std::io::ErrorKind::Unsupported,
            "IPC server not yet supported on this platform",
        ))
    }

    /// Generate the next callback ID
    pub fn next_callback_id(&self) -> u64 {
        self.next_callback_id.fetch_add(1, Ordering::SeqCst)
    }

    /// Poll for incoming requests (non-blocking)
    #[cfg(unix)]
    pub fn poll(&mut self) -> Option<IpcRequest> {
        // Try to accept new connections
        match self.listener.accept() {
            Ok((stream, _addr)) => {
                if self.client.is_some() {
                    log::warn!("IPC: rejecting second client, only one connection allowed");
                    Self::reject_connection(stream);
                } else {
                    log::info!("IPC client connected");
                    if let Err(e) = stream.set_nonblocking(true) {
                        log::error!("Failed to set non-blocking: {}", e);
                        return None;
                    }
                    let reader_stream = match stream.try_clone() {
                        Ok(s) => s,
                        Err(e) => {
                            log::error!("Failed to clone stream: {}", e);
                            return None;
                        }
                    };
                    self.client = Some(ClientConnection {
                        stream,
                        reader: BufReader::new(reader_stream),
                        client_id: String::from("unknown"),
                    });
                }
            }
            Err(ref e) if e.kind() == std::io::ErrorKind::WouldBlock => {}
            Err(e) => {
                log::error!("Failed to accept connection: {}", e);
            }
        }

        // Try to read from connected client
        if let Some(ref mut client) = self.client {
            let mut line = String::new();
            match client.reader.read_line(&mut line) {
                Ok(0) => {
                    log::info!("IPC client disconnected");
                    self.client = None;
                    return None;
                }
                Ok(_) => {
                    match serde_json::from_str::<IpcRequest>(&line) {
                        Ok(request) => {
                            log::debug!("IPC request: {:?}", request);
                            return Some(request);
                        }
                        Err(e) => {
                            log::error!("Failed to parse IPC request: {}", e);
                            log::error!("Raw line: {}", line.trim());
                        }
                    }
                }
                Err(ref e) if e.kind() == std::io::ErrorKind::WouldBlock => {}
                Err(ref e) if e.kind() == std::io::ErrorKind::Interrupted => {
                    log::debug!("IPC read interrupted, will retry");
                }
                Err(ref e) => {
                    if e.raw_os_error() == Some(11) || e.raw_os_error() == Some(35) {
                        log::debug!("IPC read EAGAIN, no data available");
                    } else {
                        log::error!(
                            "Failed to read from client: {} (kind: {:?}, os_error: {:?})",
                            e,
                            e.kind(),
                            e.raw_os_error()
                        );
                        self.client = None;
                    }
                }
            }
        }

        None
    }

    #[cfg(not(unix))]
    pub fn poll(&mut self) -> Option<IpcRequest> {
        None
    }

    /// Send a response to the connected client
    #[cfg(unix)]
    pub fn send(&mut self, response: IpcResponse) -> std::io::Result<()> {
        if let Some(ref mut client) = self.client {
            let json = serde_json::to_string(&response)?;
            writeln!(client.stream, "{}", json)?;
            client.stream.flush()?;
            log::debug!("IPC response: {:?}", response);
            Ok(())
        } else {
            Err(std::io::Error::new(
                std::io::ErrorKind::NotConnected,
                "No client connected",
            ))
        }
    }

    #[cfg(not(unix))]
    pub fn send(&mut self, _response: IpcResponse) -> std::io::Result<()> {
        Err(std::io::Error::new(
            std::io::ErrorKind::Unsupported,
            "IPC server not yet supported on this platform",
        ))
    }

    /// Set the connected client's identifier (called when Hello is received)
    pub fn set_client_id(&mut self, id: String) {
        #[cfg(unix)]
        if let Some(ref mut client) = self.client {
            client.client_id = id;
        }
    }

    /// Reject an incoming connection by sending an error and closing
    #[cfg(unix)]
    fn reject_connection(mut stream: std::os::unix::net::UnixStream) {
        let response = IpcResponse::Error {
            id: 0,
            message: "Another client is already connected".to_string(),
        };
        if let Ok(json) = serde_json::to_string(&response) {
            let _ = writeln!(stream, "{}", json);
            let _ = stream.flush();
        }
    }
}

impl Drop for IpcServer {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.socket_path);
        log::info!("IPC server stopped");
    }
}
