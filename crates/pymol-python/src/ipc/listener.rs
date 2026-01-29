//! Background Callback Listener
//!
//! Provides a background thread that listens for CallbackRequest messages from the GUI
//! and dispatches them to registered Python callbacks.

use std::collections::HashMap;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::thread::{self, JoinHandle};
use std::time::Duration;

use pyo3::prelude::*;
use pyo3::types::{PyString, PyTuple};

use super::protocol::{IpcResponse, OutputMessage};
use super::IpcClient;

/// Extended command callback storage
pub struct ExtendedCommands {
    /// Map of command name -> Python function
    callbacks: HashMap<String, Py<PyAny>>,
}

impl ExtendedCommands {
    /// Create a new empty registry
    pub fn new() -> Self {
        Self {
            callbacks: HashMap::new(),
        }
    }

    /// Register a callback for a command
    pub fn register(&mut self, name: String, func: Py<PyAny>) {
        self.callbacks.insert(name, func);
    }

    /// Unregister a callback
    pub fn unregister(&mut self, name: &str) -> Option<Py<PyAny>> {
        self.callbacks.remove(name)
    }

    /// Get a callback by name
    pub fn get(&self, name: &str) -> Option<&Py<PyAny>> {
        self.callbacks.get(name)
    }

    /// Check if a command is registered
    pub fn contains(&self, name: &str) -> bool {
        self.callbacks.contains_key(name)
    }

    /// Get all registered command names
    pub fn names(&self) -> impl Iterator<Item = &str> {
        self.callbacks.keys().map(|s| s.as_str())
    }

    /// Get the number of registered commands
    pub fn len(&self) -> usize {
        self.callbacks.len()
    }

    /// Check if the registry is empty
    pub fn is_empty(&self) -> bool {
        self.callbacks.is_empty()
    }
}

impl Default for ExtendedCommands {
    fn default() -> Self {
        Self::new()
    }
}

/// Background listener for callback requests from the GUI
pub struct CallbackListener {
    /// Background thread handle
    thread: Option<JoinHandle<()>>,
    /// Shutdown signal
    shutdown: Arc<AtomicBool>,
}

impl CallbackListener {
    /// Start the callback listener in a background thread
    ///
    /// # Arguments
    /// * `client` - Shared IPC client for communication (Option to handle disconnected state)
    /// * `callbacks` - Shared registry of extended command callbacks
    ///
    /// # Returns
    /// A CallbackListener that manages the background thread
    pub fn start(
        client: Arc<Mutex<Option<IpcClient>>>,
        callbacks: Arc<Mutex<ExtendedCommands>>,
    ) -> Self {
        let shutdown = Arc::new(AtomicBool::new(false));
        let shutdown_clone = shutdown.clone();

        let thread = thread::spawn(move || {
            listener_thread(client, callbacks, shutdown_clone);
        });

        Self {
            thread: Some(thread),
            shutdown,
        }
    }

    /// Signal the listener to stop and wait for it to finish
    pub fn stop(&mut self) {
        self.shutdown.store(true, Ordering::SeqCst);
        if let Some(thread) = self.thread.take() {
            let _ = thread.join();
        }
    }

    /// Check if the listener is still running
    pub fn is_running(&self) -> bool {
        self.thread.as_ref().map(|t| !t.is_finished()).unwrap_or(false)
    }
}

impl Drop for CallbackListener {
    fn drop(&mut self) {
        self.stop();
    }
}

/// The main listener thread function
fn listener_thread(
    client: Arc<Mutex<Option<IpcClient>>>,
    callbacks: Arc<Mutex<ExtendedCommands>>,
    shutdown: Arc<AtomicBool>,
) {
    log::info!("Callback listener thread started");

    while !shutdown.load(Ordering::SeqCst) {
        // Try to receive a message (non-blocking)
        let response = {
            let mut client_guard = match client.lock() {
                Ok(guard) => guard,
                Err(e) => {
                    log::error!("Failed to lock client: {}", e);
                    thread::sleep(Duration::from_millis(100));
                    continue;
                }
            };

            // Check if client is available
            let Some(ref mut ipc_client) = *client_guard else {
                thread::sleep(Duration::from_millis(100));
                continue;
            };

            match ipc_client.try_recv() {
                Ok(Some(response)) => Some(response),
                Ok(None) => None,
                Err(e) if e.kind() == std::io::ErrorKind::WouldBlock => None,
                Err(e) => {
                    log::warn!("Error receiving from IPC: {}", e);
                    None
                }
            }
        };

        // Handle the response if we got one
        if let Some(response) = response {
            match response {
                IpcResponse::CallbackRequest { id, name, args } => {
                    log::debug!("Received callback request: {} {:?}", name, args);
                    handle_callback_request(&client, &callbacks, id, &name, &args);
                }
                IpcResponse::Closing => {
                    log::info!("GUI is closing, stopping listener");
                    break;
                }
                _ => {
                    // Other responses are handled by the main thread
                    log::debug!("Ignoring response in listener: {:?}", response);
                }
            }
        } else {
            // No message, sleep briefly to avoid busy-waiting
            thread::sleep(Duration::from_millis(10));
        }
    }

    log::info!("Callback listener thread stopped");
}

/// Handle a callback request from the GUI
fn handle_callback_request(
    client: &Arc<Mutex<Option<IpcClient>>>,
    callbacks: &Arc<Mutex<ExtendedCommands>>,
    id: u64,
    name: &str,
    args: &[String],
) {
    let mut output = Vec::new();
    let (success, error) = if name == "runpy" {
        // Special case: runpy command for Python script execution
        handle_runpy(args, &mut output)
    } else {
        // User-defined extended command
        handle_extended_command(callbacks, name, args, &mut output)
    };

    // Send response back to GUI
    send_callback_response(client, id, success, error, output);
}

/// Handle the runpy command
fn handle_runpy(args: &[String], output: &mut Vec<OutputMessage>) -> (bool, Option<String>) {
    if let Some(script_path) = args.first() {
        let namespace = args.get(1).map(|s| s.as_str()).unwrap_or("global");

        // Attach to the Python interpreter and run the script
        Python::attach(|py| {
            match crate::scripting::run_python_script(py, script_path, namespace) {
                Ok(()) => {
                    output.push(OutputMessage::info(format!("Executed: {}", script_path)));
                    (true, None)
                }
                Err(e) => {
                    let err_msg = e.to_string();
                    output.push(OutputMessage::error(&err_msg));
                    (false, Some(err_msg))
                }
            }
        })
    } else {
        let err_msg = "runpy requires a script path".to_string();
        output.push(OutputMessage::error(&err_msg));
        (false, Some(err_msg))
    }
}

/// Handle a user-defined extended command
fn handle_extended_command(
    callbacks: &Arc<Mutex<ExtendedCommands>>,
    name: &str,
    args: &[String],
    output: &mut Vec<OutputMessage>,
) -> (bool, Option<String>) {
    // Get the callback function
    let func = {
        let callbacks_guard = match callbacks.lock() {
            Ok(guard) => guard,
            Err(e) => {
                let err_msg = format!("Failed to lock callbacks: {}", e);
                output.push(OutputMessage::error(&err_msg));
                return (false, Some(err_msg));
            }
        };
        // Clone the Py<PyAny> - requires attachment to Python interpreter
        Python::attach(|py| {
            callbacks_guard.get(name).map(|f| f.clone_ref(py))
        })
    };

    if let Some(func) = func {
        // Attach to the Python interpreter and call the Python function
        Python::attach(|py| {
            // Convert args to Python strings
            let py_args: Vec<Py<PyString>> = args
                .iter()
                .map(|s| PyString::new(py, s).into())
                .collect();
            
            // Create a tuple of arguments
            let py_tuple = match PyTuple::new(py, &py_args) {
                Ok(t) => t,
                Err(e) => {
                    let err_msg = format!("Failed to create argument tuple: {}", e);
                    output.push(OutputMessage::error(&err_msg));
                    return (false, Some(err_msg));
                }
            };

            // Call the Python function
            match func.call1(py, py_tuple) {
                Ok(_) => (true, None),
                Err(e) => {
                    let err_msg = e.to_string();
                    output.push(OutputMessage::error(&err_msg));
                    (false, Some(err_msg))
                }
            }
        })
    } else {
        let err_msg = format!("Unknown callback: {}", name);
        output.push(OutputMessage::error(&err_msg));
        (false, Some(err_msg))
    }
}

/// Send a callback response back to the GUI
fn send_callback_response(
    client: &Arc<Mutex<Option<IpcClient>>>,
    id: u64,
    success: bool,
    error: Option<String>,
    output: Vec<OutputMessage>,
) {
    let mut client_guard = match client.lock() {
        Ok(guard) => guard,
        Err(e) => {
            log::error!("Failed to lock client to send response: {}", e);
            return;
        }
    };

    let Some(ref mut ipc_client) = *client_guard else {
        log::error!("Client not available to send callback response");
        return;
    };

    if let Err(e) = ipc_client.send_callback_response(id, success, error, output) {
        log::error!("Failed to send callback response: {}", e);
    }
}
