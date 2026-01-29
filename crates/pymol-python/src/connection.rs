//! Connection management for pymol-rs IPC
//!
//! Handles connecting to an existing IPC server or spawning a new headless server.

use std::path::{Path, PathBuf};
use std::process::{Child, Command};
use std::sync::{Arc, Mutex};
use std::time::Duration;

use crate::ipc::IpcClient;

/// Default socket path when PYMOL_RS_IPC_SOCKET is not set
pub const DEFAULT_SOCKET_PATH: &str = "/tmp/pymol-rs.sock";

/// Environment variable for custom socket path
pub const SOCKET_ENV_VAR: &str = "PYMOL_RS_IPC_SOCKET";

/// Connection error types
#[derive(Debug)]
pub enum ConnectionError {
    /// Failed to find the pymol-rs binary
    BinaryNotFound(String),
    /// Failed to spawn the server process
    SpawnFailed(String),
    /// Failed to connect to the IPC server
    ConnectionFailed(String),
    /// Server health check failed
    HealthCheckFailed(String),
}

impl std::fmt::Display for ConnectionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::BinaryNotFound(msg) => write!(f, "pymol-rs binary not found: {}", msg),
            Self::SpawnFailed(msg) => write!(f, "Failed to spawn server: {}", msg),
            Self::ConnectionFailed(msg) => write!(f, "Failed to connect: {}", msg),
            Self::HealthCheckFailed(msg) => write!(f, "Health check failed: {}", msg),
        }
    }
}

impl std::error::Error for ConnectionError {}

/// Manages the connection to the pymol-rs IPC server
///
/// The client is stored in an `Arc<Mutex<...>>` to allow sharing with
/// background threads (e.g., the callback listener).
pub struct Connection {
    /// IPC client for communication (wrapped for thread-safe sharing)
    client: Arc<Mutex<Option<IpcClient>>>,
    /// Child process handle (if we spawned the server)
    child_process: Option<Child>,
    /// Path to the IPC socket
    socket_path: PathBuf,
    /// Whether we own (spawned) the server process
    owns_server: bool,
}

impl Connection {
    /// Get the shared client Arc for use with background threads
    ///
    /// This returns a clone of the Arc, allowing the client to be shared
    /// with callback listeners or other threads.
    pub fn client_arc(&self) -> Arc<Mutex<Option<IpcClient>>> {
        self.client.clone()
    }

    /// Execute a function with the locked client
    ///
    /// Returns None if the client has been taken or the lock fails.
    pub fn with_client<F, R>(&self, f: F) -> Option<R>
    where
        F: FnOnce(&mut IpcClient) -> R,
    {
        let mut guard = self.client.lock().ok()?;
        guard.as_mut().map(f)
    }

    /// Get the socket path
    pub fn socket_path(&self) -> &Path {
        &self.socket_path
    }

    /// Check if we own (spawned) the server process
    pub fn owns_server(&self) -> bool {
        self.owns_server
    }

    /// Take the IPC client out of the connection
    ///
    /// This transfers ownership of the client. The connection will no longer
    /// be able to communicate with the server after this call.
    pub fn take_client(&mut self) -> Option<IpcClient> {
        self.client.lock().ok().and_then(|mut guard| guard.take())
    }

    /// Check if the server is still running (only meaningful if we own it)
    pub fn is_server_running(&mut self) -> bool {
        if let Some(ref mut child) = self.child_process {
            match child.try_wait() {
                Ok(Some(_)) => false, // Process has exited
                Ok(None) => true,     // Process is still running
                Err(_) => false,      // Error checking - assume not running
            }
        } else {
            // We didn't spawn it, assume it's running if we're connected
            true
        }
    }

    /// Request the server to quit
    pub fn quit(&self) -> Result<(), std::io::Error> {
        if let Ok(mut guard) = self.client.lock() {
            if let Some(ref mut client) = *guard {
                return client.quit();
            }
        }
        Ok(())
    }
}

impl Drop for Connection {
    fn drop(&mut self) {
        // Only clean up if we own the server
        if self.owns_server {
            // Try to close gracefully
            let _ = self.quit();

            // Give it a moment to close
            std::thread::sleep(Duration::from_millis(100));

            // Kill if still running
            if let Some(ref mut child) = self.child_process {
                if let Ok(None) = child.try_wait() {
                    let _ = child.kill();
                }
            }

            // Clean up socket file
            let _ = std::fs::remove_file(&self.socket_path);
        }
    }
}

/// Get the socket path from environment or use default
pub fn get_socket_path() -> PathBuf {
    std::env::var(SOCKET_ENV_VAR)
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from(DEFAULT_SOCKET_PATH))
}

/// Establish a connection to the pymol-rs IPC server
///
/// This function tries to:
/// 1. Connect to an existing server at the socket path
/// 2. If no server is running, spawn a new headless server
///
/// # Returns
/// A `Connection` that manages the IPC client and optionally the server process
pub fn establish_connection() -> Result<Connection, ConnectionError> {
    let socket_path = get_socket_path();
    
    log::info!("Attempting to connect to IPC server at {:?}", socket_path);

    // Try connecting to existing server
    if let Ok(mut client) = IpcClient::connect(&socket_path) {
        // Verify the server is responsive
        match client.ping() {
            Ok(true) => {
                log::info!("Connected to existing IPC server at {:?}", socket_path);
                return Ok(Connection {
                    client: Arc::new(Mutex::new(Some(client))),
                    child_process: None,
                    socket_path,
                    owns_server: false,
                });
            }
            Ok(false) => {
                log::warn!("Server at {:?} did not respond to ping", socket_path);
            }
            Err(e) => {
                log::warn!("Failed to ping server at {:?}: {}", socket_path, e);
            }
        }
    }

    // No existing server, spawn a new headless one
    log::info!("No existing server found, spawning headless server");
    spawn_headless_server(&socket_path)
}

/// Spawn a new headless pymol-rs server
fn spawn_headless_server(socket_path: &Path) -> Result<Connection, ConnectionError> {
    // Remove existing socket file if present
    let _ = std::fs::remove_file(socket_path);

    // Spawn using the Python entry point (pymol-rs command from pip install)
    let mut child = spawn_via_python_entry_point(socket_path)
        .or_else(|e| {
            log::info!("Python entry point failed ({}), trying native binary", e);
            spawn_via_native_binary(socket_path)
        })?;

    // Wait for the server to start, checking if it's still running
    for i in 0..20 {
        std::thread::sleep(Duration::from_millis(100));
        
        // Check if the child process has exited
        match child.try_wait() {
            Ok(Some(status)) => {
                // Process has exited - this is an error
                return Err(ConnectionError::SpawnFailed(format!(
                    "Server process exited immediately with status: {}",
                    status
                )));
            }
            Ok(None) => {
                // Process is still running, good
                log::debug!("Server process still running after {}ms", (i + 1) * 100);
            }
            Err(e) => {
                return Err(ConnectionError::SpawnFailed(format!(
                    "Failed to check server process status: {}",
                    e
                )));
            }
        }
        
        // Check if socket file exists (server has started listening)
        if socket_path.exists() {
            log::info!("Socket file created after {}ms, attempting connection", (i + 1) * 100);
            break;
        }
    }

    // Connect with retry
    let mut client = connect_with_retry(socket_path, 10, Duration::from_millis(200))?;

    // Verify connection with ping before wrapping
    match client.ping() {
        Ok(true) => {
            log::info!("Connected to spawned headless server");
            Ok(Connection {
                client: Arc::new(Mutex::new(Some(client))),
                child_process: Some(child),
                socket_path: socket_path.to_path_buf(),
                owns_server: true,
            })
        }
        Ok(false) => Err(ConnectionError::HealthCheckFailed(
            "Server did not respond to ping".to_string(),
        )),
        Err(e) => Err(ConnectionError::HealthCheckFailed(format!(
            "Failed to ping server: {}",
            e
        ))),
    }
}

/// Spawn the server using the pymol-rs command (Python entry point)
///
/// This is the primary method when pymol-rs is installed from a wheel.
fn spawn_via_python_entry_point(socket_path: &Path) -> Result<Child, ConnectionError> {
    // Find the pymol-rs command in PATH
    let pymol_rs = which::which("pymol-rs").map_err(|_| {
        ConnectionError::SpawnFailed("pymol-rs command not found in PATH".to_string())
    })?;

    log::info!("Spawning headless server: {:?} --ipc {:?} --headless", pymol_rs, socket_path);

    // Spawn the process with --quiet to suppress logs in the Python console
    let child = Command::new(&pymol_rs)
        .arg("--ipc")
        .arg(socket_path)
        .arg("--headless")
        .arg("--quiet")
        .stdin(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .spawn()
        .map_err(|e| ConnectionError::SpawnFailed(format!(
            "Failed to spawn pymol-rs: {}", e
        )))?;

    Ok(child)
}

/// Spawn the server using a native binary
///
/// This is used when pymol-rs is installed via cargo.
fn spawn_via_native_binary(socket_path: &Path) -> Result<Child, ConnectionError> {
    let binary_path = find_pymol_rs_binary()?;

    log::info!("Spawning headless server: {:?}", binary_path);

    // Note: Native binary doesn't support --quiet flag, so we suppress I/O
    Command::new(&binary_path)
        .arg("--ipc")
        .arg(socket_path)
        .arg("--headless")
        .stdin(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .spawn()
        .map_err(|e| ConnectionError::SpawnFailed(format!("Failed to spawn pymol-rs: {}", e)))
}

/// Connect to an IPC server with retry logic
pub fn connect_with_retry(
    socket_path: &Path,
    max_attempts: u32,
    delay: Duration,
) -> Result<IpcClient, ConnectionError> {
    for attempt in 0..max_attempts {
        match IpcClient::connect(socket_path) {
            Ok(client) => {
                log::debug!("Connected on attempt {}", attempt + 1);
                return Ok(client);
            }
            Err(e) => {
                if attempt < max_attempts - 1 {
                    log::debug!(
                        "Connection attempt {} failed: {}, retrying...",
                        attempt + 1,
                        e
                    );
                    std::thread::sleep(delay);
                } else {
                    return Err(ConnectionError::ConnectionFailed(format!(
                        "Failed to connect after {} attempts: {}",
                        max_attempts, e
                    )));
                }
            }
        }
    }

    Err(ConnectionError::ConnectionFailed(
        "Max attempts reached".to_string(),
    ))
}

/// Find the pymol-rs binary
///
/// Searches in order:
/// 1. PYMOL_RS_BINARY environment variable
/// 2. In PATH using `which`
/// 3. Common install locations
pub fn find_pymol_rs_binary() -> Result<PathBuf, ConnectionError> {
    // 1. Check environment variable
    if let Ok(path) = std::env::var("PYMOL_RS_BINARY") {
        let path = PathBuf::from(path);
        if path.exists() {
            log::info!("Found pymol-rs from PYMOL_RS_BINARY: {:?}", path);
            return Ok(path);
        }
    }

    // 2. Check PATH using `which`
    if let Ok(path) = which::which("pymol-rs") {
        log::info!("Found pymol-rs in PATH: {:?}", path);
        return Ok(path);
    }

    // 3. Check common install locations
    let common_paths = [
        "/usr/local/bin/pymol-rs",
        "/usr/bin/pymol-rs",
        "~/.cargo/bin/pymol-rs",
    ];

    for path_str in &common_paths {
        let expanded = shellexpand::tilde(path_str);
        let path = PathBuf::from(expanded.as_ref());
        if path.exists() {
            log::info!("Found pymol-rs at: {:?}", path);
            return Ok(path);
        }
    }

    Err(ConnectionError::BinaryNotFound(
        "pymol-rs binary not found. Install it with 'cargo install pymol-gui' or set PYMOL_RS_BINARY environment variable.".to_string(),
    ))
}

