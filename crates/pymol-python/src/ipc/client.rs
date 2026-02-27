//! IPC Client implementation
//!
//! Provides a client that connects to the pymol-rs GUI via Unix domain socket.

use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use super::protocol::{IpcRequest, IpcResponse};

/// IPC Client for communicating with the GUI
pub struct IpcClient {
    /// Path to the socket file
    socket_path: PathBuf,
    /// Connected stream
    #[cfg(unix)]
    stream: std::os::unix::net::UnixStream,
    /// Reader for the stream
    #[cfg(unix)]
    reader: BufReader<std::os::unix::net::UnixStream>,
    /// Next request ID
    next_id: AtomicU64,
}

impl IpcClient {
    /// Connect to the GUI IPC server with the default client identifier
    #[cfg(unix)]
    pub fn connect(socket_path: &Path) -> std::io::Result<Self> {
        Self::connect_with_id(socket_path, "pymol-python")
    }

    /// Connect to the GUI IPC server with a specific client identifier
    ///
    /// # Arguments
    /// * `socket_path` - Path to the Unix domain socket
    /// * `client_id` - Identifier string sent to the server during handshake
    ///
    /// # Errors
    /// Returns an error if the connection fails or the server rejects the connection
    #[cfg(unix)]
    pub fn connect_with_id(socket_path: &Path, client_id: &str) -> std::io::Result<Self> {
        use std::os::unix::net::UnixStream;

        log::info!("Connecting to IPC server at {:?}", socket_path);

        let stream = UnixStream::connect(socket_path)?;
        let reader_stream = stream.try_clone()?;
        let reader = BufReader::new(reader_stream);

        let mut client = Self {
            socket_path: socket_path.to_path_buf(),
            stream,
            reader,
            next_id: AtomicU64::new(1),
        };

        // Send handshake
        let hello = IpcRequest::Hello {
            client_id: client_id.to_string(),
        };
        client.send(&hello)?;

        // Read acknowledgment
        let ack = client.recv()?;
        if let IpcResponse::Error { message, .. } = ack {
            return Err(std::io::Error::new(
                std::io::ErrorKind::ConnectionRefused,
                format!("Server rejected connection: {}", message),
            ));
        }

        Ok(client)
    }

    #[cfg(not(unix))]
    pub fn connect(_socket_path: &Path) -> std::io::Result<Self> {
        Err(std::io::Error::new(
            std::io::ErrorKind::Unsupported,
            "IPC client not yet supported on this platform",
        ))
    }

    #[cfg(not(unix))]
    pub fn connect_with_id(_socket_path: &Path, _client_id: &str) -> std::io::Result<Self> {
        Err(std::io::Error::new(
            std::io::ErrorKind::Unsupported,
            "IPC client not yet supported on this platform",
        ))
    }

    /// Get the socket path
    pub fn socket_path(&self) -> &Path {
        &self.socket_path
    }

    /// Generate the next request ID
    pub fn next_id(&self) -> u64 {
        self.next_id.fetch_add(1, Ordering::SeqCst)
    }

    /// Send a request to the GUI
    #[cfg(unix)]
    pub fn send(&mut self, request: &IpcRequest) -> std::io::Result<()> {
        let json = serde_json::to_string(request)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        writeln!(self.stream, "{}", json)?;
        self.stream.flush()?;
        log::debug!("IPC sent: {:?}", request);
        Ok(())
    }

    #[cfg(not(unix))]
    pub fn send(&mut self, _request: &IpcRequest) -> std::io::Result<()> {
        Err(std::io::Error::new(
            std::io::ErrorKind::Unsupported,
            "IPC client not yet supported on this platform",
        ))
    }

    /// Receive a response from the GUI (blocking)
    #[cfg(unix)]
    pub fn recv(&mut self) -> std::io::Result<IpcResponse> {
        let mut line = String::new();
        self.reader.read_line(&mut line)?;
        
        if line.is_empty() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::UnexpectedEof,
                "Connection closed",
            ));
        }

        let response: IpcResponse = serde_json::from_str(&line)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        
        log::debug!("IPC received: {:?}", response);
        Ok(response)
    }

    #[cfg(not(unix))]
    pub fn recv(&mut self) -> std::io::Result<IpcResponse> {
        Err(std::io::Error::new(
            std::io::ErrorKind::Unsupported,
            "IPC client not yet supported on this platform",
        ))
    }

    /// Try to receive a response without blocking
    #[cfg(unix)]
    pub fn try_recv(&mut self) -> std::io::Result<Option<IpcResponse>> {
        use std::io::ErrorKind;

        // Set non-blocking mode temporarily
        self.stream.set_nonblocking(true)?;
        
        let mut line = String::new();
        let result = self.reader.read_line(&mut line);
        
        // Restore blocking mode
        self.stream.set_nonblocking(false)?;

        match result {
            Ok(0) => Err(std::io::Error::new(
                ErrorKind::UnexpectedEof,
                "Connection closed",
            )),
            Ok(_) => {
                let response: IpcResponse = serde_json::from_str(&line)
                    .map_err(|e| std::io::Error::new(ErrorKind::InvalidData, e))?;
                Ok(Some(response))
            }
            Err(ref e) if e.kind() == ErrorKind::WouldBlock => Ok(None),
            Err(e) => Err(e),
        }
    }

    #[cfg(not(unix))]
    pub fn try_recv(&mut self) -> std::io::Result<Option<IpcResponse>> {
        Ok(None)
    }

    /// Send a request and wait for response
    pub fn request(&mut self, request: &IpcRequest) -> std::io::Result<IpcResponse> {
        self.send(request)?;
        self.recv()
    }

    /// Execute a command on the GUI
    pub fn execute(&mut self, command: &str) -> std::io::Result<IpcResponse> {
        let id = self.next_id();
        let request = IpcRequest::Execute {
            id,
            command: command.to_string(),
            silent: false,
        };
        self.request(&request)
    }

    /// Execute a command on the GUI in silent mode (no echo or info output)
    pub fn execute_silent(&mut self, command: &str) -> std::io::Result<IpcResponse> {
        let id = self.next_id();
        let request = IpcRequest::Execute {
            id,
            command: command.to_string(),
            silent: true,
        };
        self.request(&request)
    }

    /// Register an external command with the GUI
    pub fn register_command(&mut self, name: &str, help: Option<&str>) -> std::io::Result<()> {
        let request = IpcRequest::RegisterCommand {
            name: name.to_string(),
            help: help.map(|s| s.to_string()),
        };
        self.send(&request)
    }

    /// Unregister an external command from the GUI
    pub fn unregister_command(&mut self, name: &str) -> std::io::Result<()> {
        let request = IpcRequest::UnregisterCommand {
            name: name.to_string(),
        };
        self.send(&request)
    }

    /// Send callback response to the GUI
    pub fn send_callback_response(
        &mut self,
        id: u64,
        success: bool,
        error: Option<String>,
        output: Vec<super::protocol::OutputMessage>,
    ) -> std::io::Result<()> {
        let request = IpcRequest::CallbackResponse {
            id,
            success,
            error,
            output,
        };
        self.send(&request)
    }

    /// Get object names from the GUI
    pub fn get_names(&mut self) -> std::io::Result<Vec<String>> {
        let id = self.next_id();
        let request = IpcRequest::GetNames { id };
        let response = self.request(&request)?;
        
        match response {
            IpcResponse::Value { value, .. } => {
                serde_json::from_value(value)
                    .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))
            }
            IpcResponse::Error { message, .. } => {
                Err(std::io::Error::new(std::io::ErrorKind::Other, message))
            }
            _ => Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Unexpected response",
            )),
        }
    }

    /// Ping the GUI (health check)
    pub fn ping(&mut self) -> std::io::Result<bool> {
        let id = self.next_id();
        let request = IpcRequest::Ping { id };
        let response = self.request(&request)?;
        
        match response {
            IpcResponse::Pong { .. } => Ok(true),
            _ => Ok(false),
        }
    }

    /// Request the GUI to quit
    pub fn quit(&mut self) -> std::io::Result<()> {
        self.send(&IpcRequest::Quit)
    }

    /// Show the GUI window (make it visible)
    pub fn show_window(&mut self) -> std::io::Result<()> {
        let id = self.next_id();
        let request = IpcRequest::ShowWindow { id };
        let response = self.request(&request)?;

        match response {
            IpcResponse::Ok { .. } => Ok(()),
            IpcResponse::Error { message, .. } => {
                Err(std::io::Error::new(std::io::ErrorKind::Other, message))
            }
            _ => Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Unexpected response to ShowWindow",
            )),
        }
    }

    /// Hide the GUI window (make it invisible)
    pub fn hide_window(&mut self) -> std::io::Result<()> {
        let id = self.next_id();
        let request = IpcRequest::HideWindow { id };
        let response = self.request(&request)?;

        match response {
            IpcResponse::Ok { .. } => Ok(()),
            IpcResponse::Error { message, .. } => {
                Err(std::io::Error::new(std::io::ErrorKind::Other, message))
            }
            _ => Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Unexpected response to HideWindow",
            )),
        }
    }

    /// Count atoms matching a selection
    pub fn count_atoms(&mut self, selection: &str) -> std::io::Result<usize> {
        let id = self.next_id();
        let request = IpcRequest::CountAtoms {
            id,
            selection: selection.to_string(),
        };
        let response = self.request(&request)?;

        match response {
            IpcResponse::Value { value, .. } => {
                serde_json::from_value(value)
                    .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))
            }
            IpcResponse::Error { message, .. } => {
                Err(std::io::Error::new(std::io::ErrorKind::Other, message))
            }
            _ => Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Unexpected response to CountAtoms",
            )),
        }
    }

    /// Execute a command and check if it succeeded
    pub fn execute_ok(&mut self, command: &str) -> std::io::Result<()> {
        let response = self.execute(command)?;
        match response {
            IpcResponse::Ok { .. } => Ok(()),
            IpcResponse::Error { message, .. } => {
                Err(std::io::Error::new(std::io::ErrorKind::Other, message))
            }
            _ => Ok(()), // Consider other responses as success
        }
    }

    /// Get the current view matrix (18 floats)
    pub fn get_view(&mut self) -> std::io::Result<Option<Vec<f64>>> {
        let id = self.next_id();
        let request = IpcRequest::GetView { id };
        let response = self.request(&request)?;

        match response {
            IpcResponse::Value { value, .. } => {
                // Parse the JSON array of floats
                if let Some(arr) = value.as_array() {
                    let values: Vec<f64> = arr
                        .iter()
                        .filter_map(|v| v.as_f64())
                        .collect();
                    if values.len() >= 18 {
                        Ok(Some(values))
                    } else {
                        Ok(None)
                    }
                } else {
                    Ok(None)
                }
            }
            IpcResponse::Error { message, .. } => {
                Err(std::io::Error::new(std::io::ErrorKind::Other, message))
            }
            _ => Ok(None),
        }
    }
}
