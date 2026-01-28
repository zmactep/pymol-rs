//! Generic async task infrastructure for GUI
//!
//! Provides a reusable system for running async operations without blocking
//! the main thread. New task types are added by implementing the `AsyncTask` trait.
//!
//! # Architecture
//!
//! The task runner maintains a tokio runtime and an mpsc channel pair:
//! - Tasks are spawned on the tokio runtime
//! - Results are sent back via the channel
//! - The main thread polls for completed tasks each frame
//!
//! # Adding New Task Types
//!
//! 1. Create a new struct for your task with necessary parameters
//! 2. Implement `AsyncTask` trait for the struct
//! 3. Add a new variant to `AsyncTaskResult` for the result
//! 4. Handle the result in `App::process_async_tasks()`

use std::future::Future;
use std::pin::Pin;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::{self, Receiver, Sender};
use std::sync::{Arc, Mutex};
use std::time::Duration;
use tokio::runtime::Runtime;

use pymol_io::FetchFormat;
use pymol_mol::ObjectMolecule;

/// Default timeout for fetch operations (10 seconds)
const FETCH_TIMEOUT_SECS: u64 = 10;

// ============================================================================
// AsyncTask Trait
// ============================================================================

/// Trait for async tasks that can be spawned by the task runner
///
/// Implement this trait for each type of background operation. The trait
/// provides the notification message to display while the task is running
/// and the async execution logic.
///
/// # Example
///
/// ```ignore
/// pub struct MyTask {
///     pub param: String,
/// }
///
/// impl AsyncTask for MyTask {
///     fn notification_message(&self) -> String {
///         format!("Processing {}...", self.param)
///     }
///
///     fn execute(self) -> Pin<Box<dyn Future<Output = AsyncTaskResult> + Send>> {
///         Box::pin(async move {
///             // Do async work here
///             AsyncTaskResult::MyResult { ... }
///         })
///     }
/// }
/// ```
pub trait AsyncTask: Send + 'static {
    /// Get the notification message to display while this task is running
    ///
    /// This message will be shown in the notification overlay in the UI.
    fn notification_message(&self) -> String;

    /// Execute the async operation and return the result
    ///
    /// This method consumes the task and returns a boxed future that
    /// produces the task result.
    fn execute(self) -> Pin<Box<dyn Future<Output = AsyncTaskResult> + Send>>;
}

// ============================================================================
// Task Result Enum
// ============================================================================

/// Result of a completed async task
///
/// Each variant represents a different type of background operation.
/// Add new variants here for future async operations.
pub enum AsyncTaskResult {
    /// Fetch from RCSB PDB completed
    Fetch {
        /// Object name to use when adding to registry
        name: String,
        /// PDB code that was fetched
        code: String,
        /// Result: molecule data or error message
        result: Result<ObjectMolecule, String>,
    },
    // Future task results go here:
    // Save { path: String, result: Result<(), String> },
    // Calculate { name: String, result: Result<Data, String> },
}

// ============================================================================
// FetchTask - Fetches structures from RCSB PDB
// ============================================================================

/// Task to fetch a molecular structure from RCSB PDB
///
/// This task downloads a structure in the specified format (CIF or PDB)
/// and returns it as an `ObjectMolecule`.
pub struct FetchTask {
    /// PDB ID to fetch (e.g., "1ubq")
    pub code: String,
    /// Object name to use when adding to registry
    pub name: String,
    /// Format to fetch (CIF or PDB)
    pub format: FetchFormat,
}

impl FetchTask {
    /// Create a new fetch task
    pub fn new(code: String, name: String, format: FetchFormat) -> Self {
        Self { code, name, format }
    }
}

impl AsyncTask for FetchTask {
    fn notification_message(&self) -> String {
        format!("Fetching {}...", self.code)
    }

    fn execute(self) -> Pin<Box<dyn Future<Output = AsyncTaskResult> + Send>> {
        Box::pin(async move {
            let timeout_duration = Duration::from_secs(FETCH_TIMEOUT_SECS);

            let result = match tokio::time::timeout(
                timeout_duration,
                pymol_io::fetch_async(&self.code, self.format),
            )
            .await
            {
                Ok(Ok(mol)) => Ok(mol),
                Ok(Err(e)) => Err(format!("Network error: {}", e)),
                Err(_) => Err(format!(
                    "Timeout: fetch did not complete within {} seconds",
                    FETCH_TIMEOUT_SECS
                )),
            };

            AsyncTaskResult::Fetch {
                name: self.name,
                code: self.code,
                result,
            }
        })
    }
}

// ============================================================================
// AsyncTaskRunner
// ============================================================================

/// Async task runner - manages tokio runtime and result channel
///
/// This struct owns the tokio runtime and provides a generic `spawn` method
/// for running any task implementing the `AsyncTask` trait. Results are
/// collected via an mpsc channel that can be polled from the main thread.
///
/// The runner also tracks notification messages for all pending tasks,
/// which can be queried to display in the UI.
pub struct AsyncTaskRunner {
    /// Tokio runtime for executing async tasks
    runtime: Runtime,
    /// Sender for tasks to return results
    sender: Sender<AsyncTaskResult>,
    /// Receiver for main thread to poll for completed tasks
    receiver: Receiver<AsyncTaskResult>,
    /// Counter for pending (in-flight) tasks
    pending_count: Arc<AtomicUsize>,
    /// Messages for currently running tasks (for notification display)
    pending_messages: Arc<Mutex<Vec<String>>>,
}

impl Default for AsyncTaskRunner {
    fn default() -> Self {
        Self::new()
    }
}

impl AsyncTaskRunner {
    /// Create a new task runner with tokio runtime
    ///
    /// # Panics
    ///
    /// Panics if the tokio runtime cannot be created.
    pub fn new() -> Self {
        let runtime = Runtime::new().expect("Failed to create tokio runtime");
        let (sender, receiver) = mpsc::channel();
        Self {
            runtime,
            sender,
            receiver,
            pending_count: Arc::new(AtomicUsize::new(0)),
            pending_messages: Arc::new(Mutex::new(Vec::new())),
        }
    }

    /// Check if there are any pending (in-flight) tasks
    ///
    /// Returns `true` if at least one async task is currently running.
    pub fn has_pending_tasks(&self) -> bool {
        self.pending_count.load(Ordering::Relaxed) > 0
    }

    /// Get the number of pending (in-flight) tasks
    pub fn pending_task_count(&self) -> usize {
        self.pending_count.load(Ordering::Relaxed)
    }

    /// Get the notification messages for all pending tasks
    ///
    /// Returns a copy of all messages for currently running tasks.
    /// These can be displayed in the notification overlay.
    pub fn pending_messages(&self) -> Vec<String> {
        self.pending_messages
            .lock()
            .map(|guard| guard.clone())
            .unwrap_or_default()
    }

    /// Poll for completed tasks (non-blocking)
    ///
    /// Returns `Some(result)` if a task completed, `None` otherwise.
    /// Call this each frame to process completed tasks.
    pub fn poll(&self) -> Option<AsyncTaskResult> {
        self.receiver.try_recv().ok()
    }

    /// Spawn an async task
    ///
    /// The task runs asynchronously and sends the result back via the channel.
    /// Use `poll()` to check for completion.
    ///
    /// The notification message from the task is tracked while the task runs
    /// and can be retrieved via `pending_messages()`.
    ///
    /// # Arguments
    ///
    /// * `task` - Any type implementing the `AsyncTask` trait
    pub fn spawn<T: AsyncTask>(&self, task: T) {
        // Get the notification message before consuming the task
        let message = task.notification_message();

        // Add message to pending list
        if let Ok(mut messages) = self.pending_messages.lock() {
            messages.push(message.clone());
        }

        // Increment pending count before spawning
        self.pending_count.fetch_add(1, Ordering::Relaxed);

        let sender = self.sender.clone();
        let pending_count = self.pending_count.clone();
        let pending_messages = self.pending_messages.clone();

        // Get the future from the task
        let future = task.execute();

        self.runtime.spawn(async move {
            // Execute the task
            let result = future.await;

            // Decrement pending count when task completes
            pending_count.fetch_sub(1, Ordering::Relaxed);

            // Remove the message from pending list
            if let Ok(mut messages) = pending_messages.lock() {
                if let Some(pos) = messages.iter().position(|m| m == &message) {
                    messages.remove(pos);
                }
            }

            // Send result back to main thread (ignore send errors if receiver dropped)
            let _ = sender.send(result);
        });
    }

    /// Spawn a fetch task (convenience method)
    ///
    /// This is a convenience wrapper around `spawn` for the common case
    /// of fetching a structure from RCSB PDB.
    ///
    /// # Arguments
    ///
    /// * `code` - PDB ID to fetch (e.g., "1ubq")
    /// * `name` - Object name to use when adding to registry
    /// * `format` - Format to fetch (CIF or PDB)
    pub fn spawn_fetch(&self, code: String, name: String, format: FetchFormat) {
        self.spawn(FetchTask::new(code, name, format));
    }
}
