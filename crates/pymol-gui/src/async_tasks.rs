//! Generic async task infrastructure
//!
//! Provides traits for async tasks and their results, plus a TaskRunner
//! for executing tasks in the background.
//!
//! # Architecture
//!
//! - `TaskContext`: Abstract interface to application state (implemented by App)
//! - `TaskResult`: Trait for task results that know how to apply themselves
//! - `AsyncTask`: Trait for async tasks that produce TaskResult
//! - `TaskRunner`: Manages tokio runtime and channels for background tasks
//!
//! # Adding New Task Types
//!
//! 1. Create a new module with your task (e.g., `my_task.rs`)
//! 2. Define a result struct implementing `TaskResult`
//! 3. Define a task struct implementing `AsyncTask`
//! 4. No changes needed to this module or app.rs!

use std::future::Future;
use std::pin::Pin;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::{self, Receiver, Sender};
use std::sync::{Arc, Mutex};
use tokio::runtime::Runtime;

use pymol_mol::ObjectMolecule;

// ============================================================================
// TaskContext - Abstract interface to application state
// ============================================================================

/// Context trait providing access to application state for task results
///
/// This trait abstracts the application (App) so that task results don't
/// depend on App directly. App implements this trait.
pub trait TaskContext {
    /// Add a molecule to the registry
    fn add_molecule(&mut self, name: &str, mol: ObjectMolecule);

    /// Zoom the camera to focus on an object
    fn zoom_on(&mut self, name: &str);

    /// Print an informational message to the output
    fn print_info(&mut self, msg: String);

    /// Print a warning message to the output
    fn print_warning(&mut self, msg: String);

    /// Print an error message to the output
    fn print_error(&mut self, msg: String);

    /// Request a redraw of the viewport
    fn request_redraw(&mut self);
}

// ============================================================================
// TaskResult - Trait for async task results
// ============================================================================

/// Trait for async task results
///
/// Each result type implements this trait with its own apply() logic.
/// This allows adding new task types without modifying any central code.
pub trait TaskResult: Send + 'static {
    /// Apply this result to the application state
    ///
    /// Called when the async task completes. The result has full access
    /// to the application state via the TaskContext trait.
    fn apply(self: Box<Self>, ctx: &mut dyn TaskContext);
}

// ============================================================================
// AsyncTask - Trait for async tasks
// ============================================================================

/// Trait for async tasks that can be spawned by the TaskRunner
///
/// Tasks produce a boxed TaskResult that knows how to apply itself.
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
///     fn execute(self) -> Pin<Box<dyn Future<Output = Box<dyn TaskResult>> + Send>> {
///         Box::pin(async move {
///             // Do async work here
///             Box::new(MyResult { ... }) as Box<dyn TaskResult>
///         })
///     }
/// }
/// ```
pub trait AsyncTask: Send + 'static {
    /// Get the notification message to display while this task is running
    ///
    /// This message will be shown in the notification overlay in the UI.
    fn notification_message(&self) -> String;

    /// Execute the async operation and return a boxed result
    ///
    /// This method consumes the task and returns a boxed future that
    /// produces a boxed TaskResult.
    fn execute(self) -> Pin<Box<dyn Future<Output = Box<dyn TaskResult>> + Send>>;
}

// ============================================================================
// TaskRunner - Executes async tasks
// ============================================================================

/// Async task runner
///
/// Manages a tokio runtime and channels for running background tasks.
/// Results are polled from the main thread and applied via TaskContext.
///
/// The runner tracks notification messages for all pending tasks,
/// which can be queried to display in the UI.
pub struct TaskRunner {
    /// Tokio runtime for executing async tasks
    runtime: Runtime,
    /// Sender for tasks to return results
    sender: Sender<Box<dyn TaskResult>>,
    /// Receiver for main thread to poll for completed tasks
    receiver: Receiver<Box<dyn TaskResult>>,
    /// Counter for pending (in-flight) tasks
    pending_count: Arc<AtomicUsize>,
    /// Messages for currently running tasks (for notification display)
    pending_messages: Arc<Mutex<Vec<String>>>,
}

impl Default for TaskRunner {
    fn default() -> Self {
        Self::new()
    }
}

impl TaskRunner {
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
    pub fn poll(&self) -> Option<Box<dyn TaskResult>> {
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
}
