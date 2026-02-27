//! IPC Server implementation
//!
//! Provides a non-blocking IPC server using Unix domain sockets (on Unix)
//! or named pipes (on Windows).

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use tokio::sync::oneshot;

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
    /// Next request ID for callbacks
    next_callback_id: AtomicU64,
    /// Pending callback responses (id -> oneshot sender)
    pending_callbacks: HashMap<u64, oneshot::Sender<IpcRequest>>,
}

#[cfg(unix)]
struct ClientConnection {
    stream: std::os::unix::net::UnixStream,
    reader: BufReader<std::os::unix::net::UnixStream>,
    client_id: String,
}

impl IpcServer {
    /// Create and bind a new IPC server
    ///
    /// # Arguments
    /// * `socket_path` - Path to the Unix domain socket
    ///
    /// # Errors
    /// Returns an error if the socket cannot be created
    #[cfg(unix)]
    pub fn bind(socket_path: &Path) -> std::io::Result<Self> {
        use std::os::unix::net::UnixListener;

        // Remove existing socket file if present
        let _ = std::fs::remove_file(socket_path);

        // Create the listener
        let listener = UnixListener::bind(socket_path)?;
        listener.set_nonblocking(true)?;

        log::info!("IPC server listening on {:?}", socket_path);

        Ok(Self {
            socket_path: socket_path.to_path_buf(),
            listener,
            client: None,
            next_callback_id: AtomicU64::new(1),
            pending_callbacks: HashMap::new(),
        })
    }

    /// Placeholder for Windows - not implemented yet
    #[cfg(not(unix))]
    pub fn bind(socket_path: &Path) -> std::io::Result<Self> {
        Err(std::io::Error::new(
            std::io::ErrorKind::Unsupported,
            "IPC server not yet supported on this platform",
        ))
    }

    /// Get the socket path
    pub fn socket_path(&self) -> &Path {
        &self.socket_path
    }

    /// Generate the next callback ID
    pub fn next_callback_id(&self) -> u64 {
        self.next_callback_id.fetch_add(1, Ordering::SeqCst)
    }

    /// Poll for incoming requests (non-blocking)
    ///
    /// Returns the next request if one is available.
    /// Call this each frame from the main loop.
    #[cfg(unix)]
    pub fn poll(&mut self) -> Option<IpcRequest> {
        // Try to accept new connections
        match self.listener.accept() {
            Ok((stream, _addr)) => {
                if self.client.is_some() {
                    // Already have a client — reject the new connection
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
            Err(ref e) if e.kind() == std::io::ErrorKind::WouldBlock => {
                // No pending connections
            }
            Err(e) => {
                log::error!("Failed to accept connection: {}", e);
            }
        }

        // Try to read from connected client
        if let Some(ref mut client) = self.client {
            let mut line = String::new();
            match client.reader.read_line(&mut line) {
                Ok(0) => {
                    // EOF - client disconnected
                    log::info!("IPC client disconnected");
                    self.client = None;
                    self.pending_callbacks.clear();
                    return None;
                }
                Ok(_) => {
                    // Parse the JSON request
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
                Err(ref e) if e.kind() == std::io::ErrorKind::WouldBlock => {
                    // No data available - this is expected for non-blocking sockets
                }
                Err(ref e) if e.kind() == std::io::ErrorKind::Interrupted => {
                    // Interrupted by signal, just retry next time
                    log::debug!("IPC read interrupted, will retry");
                }
                Err(ref e) => {
                    // Check for EAGAIN (11) which is the same as WouldBlock on Unix
                    // Some systems report this instead of WouldBlock
                    if e.raw_os_error() == Some(11) || e.raw_os_error() == Some(35) {
                        // EAGAIN (11 on Linux, 35 on macOS)
                        log::debug!("IPC read EAGAIN, no data available");
                    } else {
                        // Real error - close connection
                        log::error!("Failed to read from client: {} (kind: {:?}, os_error: {:?})", 
                            e, e.kind(), e.raw_os_error());
                        self.client = None;
                        self.pending_callbacks.clear();
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

    /// Get the connected client's identifier
    pub fn client_id(&self) -> Option<&str> {
        #[cfg(unix)]
        {
            self.client.as_ref().map(|c| c.client_id.as_str())
        }
        #[cfg(not(unix))]
        {
            None
        }
    }

    /// Set the connected client's identifier (called when Hello is received)
    pub fn set_client_id(&mut self, id: String) {
        #[cfg(unix)]
        if let Some(ref mut client) = self.client {
            client.client_id = id;
        }
    }

    /// Check if a client is connected
    pub fn is_connected(&self) -> bool {
        #[cfg(unix)]
        {
            self.client.is_some()
        }
        #[cfg(not(unix))]
        {
            false
        }
    }

    /// Check if there might be pending requests
    ///
    /// Returns true if a client is connected and might have data to read.
    /// This is used to trigger IPC processing in headless mode.
    pub fn has_pending(&self) -> bool {
        // If a client is connected, there might be pending requests
        // The actual check happens in poll()
        self.is_connected()
    }

    /// Register a pending callback
    ///
    /// When the GUI sends a CallbackRequest, it registers a sender here
    /// to receive the CallbackResponse when it arrives.
    pub fn register_callback(&mut self, id: u64, sender: oneshot::Sender<IpcRequest>) {
        self.pending_callbacks.insert(id, sender);
    }

    /// Handle a CallbackResponse by forwarding it to the waiting task
    pub fn handle_callback_response(&mut self, response: &IpcRequest) -> bool {
        if let IpcRequest::CallbackResponse { id, .. } = response {
            if let Some(sender) = self.pending_callbacks.remove(id) {
                // Send via oneshot (consumes sender)
                return sender.send(response.clone()).is_ok();
            }
        }
        false
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
        // stream is dropped here, closing the connection
    }
}

impl Drop for IpcServer {
    fn drop(&mut self) {
        // Clean up the socket file
        let _ = std::fs::remove_file(&self.socket_path);
        log::info!("IPC server stopped");
    }
}
