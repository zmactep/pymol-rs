//! Generic async task infrastructure for GUI hosts.
//!
//! Tasks run on a private Tokio runtime and return results to the main thread.
//! The main thread owns scene mutation, output, and UI state, so completed
//! task results are applied by [`AppKernel`].

use std::future::Future;
use std::pin::Pin;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::mpsc::{self, Receiver, Sender};
use std::sync::{Arc, Mutex};

use tokio::runtime::Runtime;

use crate::kernel::AppKernel;

/// Identifier assigned to a spawned async task.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct TaskId(pub u64);

/// User-visible status for a pending task.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PendingTaskInfo {
    /// Stable task identifier.
    pub id: TaskId,
    /// Message shown while the task is in flight.
    pub message: String,
}

/// Result returned by an async task and applied on the main thread.
pub trait TaskResult: Send + 'static {
    /// Apply this result to the application kernel.
    fn apply(self: Box<Self>, kernel: &mut AppKernel);
}

/// Builds a result for runtime-level task failures.
pub type TaskFailureHandler = Box<dyn FnOnce(String) -> Box<dyn TaskResult> + Send>;

/// Async task that can be spawned by [`TaskRunner`].
pub trait AsyncTask: Send + 'static {
    /// Message shown while this task is running.
    fn notification_message(&self) -> String;

    /// Result to apply if the task future panics or is cancelled by the runtime.
    fn failure_handler(&self) -> TaskFailureHandler {
        Box::new(|error| Box::new(TaskFailureResult { error }))
    }

    /// Execute the task in the background.
    fn execute(self: Box<Self>) -> Pin<Box<dyn Future<Output = Box<dyn TaskResult>> + Send>>;
}

struct TaskFailureResult {
    error: String,
}

impl TaskResult for TaskFailureResult {
    fn apply(self: Box<Self>, kernel: &mut AppKernel) {
        kernel
            .output
            .print_error(format!("Background task failed: {}", self.error));
    }
}

struct CompletedTask {
    result: Box<dyn TaskResult>,
}

/// Runs async tasks and exposes pending notifications plus completed results.
pub struct TaskRunner {
    runtime: Runtime,
    sender: Sender<CompletedTask>,
    receiver: Receiver<CompletedTask>,
    next_id: AtomicU64,
    pending_count: Arc<AtomicUsize>,
    pending: Arc<Mutex<Vec<PendingTaskInfo>>>,
}

impl TaskRunner {
    /// Create a task runner backed by a Tokio multi-thread runtime.
    pub fn new() -> Self {
        let runtime = Runtime::new().expect("failed to create async task runtime");
        let (sender, receiver) = mpsc::channel();
        Self {
            runtime,
            sender,
            receiver,
            next_id: AtomicU64::new(1),
            pending_count: Arc::new(AtomicUsize::new(0)),
            pending: Arc::new(Mutex::new(Vec::new())),
        }
    }

    /// Spawn a typed task.
    pub fn spawn<T: AsyncTask>(&self, task: T) -> TaskId {
        self.spawn_boxed(Box::new(task))
    }

    /// Spawn a boxed task.
    pub fn spawn_boxed(&self, task: Box<dyn AsyncTask>) -> TaskId {
        let id = TaskId(self.next_id.fetch_add(1, Ordering::Relaxed));
        let message = task.notification_message();

        self.pending_count.fetch_add(1, Ordering::Relaxed);
        if let Ok(mut pending) = self.pending.lock() {
            pending.push(PendingTaskInfo { id, message });
        }

        let failure_handler = task.failure_handler();
        let sender = self.sender.clone();
        let pending_count = self.pending_count.clone();
        let pending = self.pending.clone();

        self.runtime.spawn(async move {
            let task_join = tokio::spawn(async move { task.execute().await });
            let result = match task_join.await {
                Ok(result) => result,
                Err(err) => {
                    let detail = if err.is_panic() {
                        "task panicked".to_string()
                    } else if err.is_cancelled() {
                        "task was cancelled".to_string()
                    } else {
                        err.to_string()
                    };
                    failure_handler(detail)
                }
            };

            pending_count.fetch_sub(1, Ordering::Relaxed);
            if let Ok(mut pending) = pending.lock() {
                pending.retain(|info| info.id != id);
            }

            let _ = sender.send(CompletedTask { result });
        });

        id
    }

    /// Number of currently pending tasks.
    pub fn pending_task_count(&self) -> usize {
        self.pending_count.load(Ordering::Relaxed)
    }

    /// Whether any tasks are currently pending.
    pub fn has_pending_tasks(&self) -> bool {
        self.pending_task_count() > 0
    }

    /// Snapshot of pending task info.
    pub fn pending_tasks(&self) -> Vec<PendingTaskInfo> {
        self.pending
            .lock()
            .map(|pending| pending.clone())
            .unwrap_or_default()
    }

    /// User-visible messages for all pending tasks.
    pub fn pending_messages(&self) -> Vec<String> {
        self.pending_tasks()
            .into_iter()
            .map(|info| info.message)
            .collect()
    }

    /// Poll one completed task result without blocking.
    pub fn poll(&self) -> Option<Box<dyn TaskResult>> {
        self.receiver.try_recv().ok().map(|task| task.result)
    }
}

impl Default for TaskRunner {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::time::{Duration, Instant};

    use super::*;

    struct GatedTask {
        release: tokio::sync::oneshot::Receiver<()>,
        applied: Arc<AtomicUsize>,
    }

    struct GatedResult {
        applied: Arc<AtomicUsize>,
    }

    struct PanicTask {
        applied: Arc<AtomicUsize>,
    }

    struct PanicResult {
        applied: Arc<AtomicUsize>,
    }

    impl AsyncTask for GatedTask {
        fn notification_message(&self) -> String {
            "Testing task...".to_string()
        }

        fn execute(self: Box<Self>) -> Pin<Box<dyn Future<Output = Box<dyn TaskResult>> + Send>> {
            Box::pin(async move {
                let _ = self.release.await;
                Box::new(GatedResult {
                    applied: self.applied,
                }) as Box<dyn TaskResult>
            })
        }
    }

    impl TaskResult for GatedResult {
        fn apply(self: Box<Self>, _kernel: &mut AppKernel) {
            self.applied.fetch_add(1, Ordering::Relaxed);
        }
    }

    impl AsyncTask for PanicTask {
        fn notification_message(&self) -> String {
            "Panicking task...".to_string()
        }

        fn failure_handler(&self) -> TaskFailureHandler {
            let applied = self.applied.clone();
            Box::new(move |_error| Box::new(PanicResult { applied }))
        }

        fn execute(self: Box<Self>) -> Pin<Box<dyn Future<Output = Box<dyn TaskResult>> + Send>> {
            Box::pin(async move {
                panic!("intentional task panic");
            })
        }
    }

    impl TaskResult for PanicResult {
        fn apply(self: Box<Self>, _kernel: &mut AppKernel) {
            self.applied.fetch_add(1, Ordering::Relaxed);
        }
    }

    #[test]
    fn task_runner_tracks_pending_and_applies_completed_result() {
        let runner = TaskRunner::new();
        let applied = Arc::new(AtomicUsize::new(0));
        let (tx, rx) = tokio::sync::oneshot::channel();

        runner.spawn(GatedTask {
            release: rx,
            applied: applied.clone(),
        });

        assert_eq!(runner.pending_task_count(), 1);
        assert_eq!(runner.pending_messages(), vec!["Testing task..."]);

        tx.send(()).expect("release task");

        let deadline = Instant::now() + Duration::from_secs(1);
        let result = loop {
            if let Some(result) = runner.poll() {
                break result;
            }
            assert!(Instant::now() < deadline, "task did not complete");
            std::thread::sleep(Duration::from_millis(5));
        };

        let mut kernel = AppKernel::new();
        result.apply(&mut kernel);

        assert_eq!(applied.load(Ordering::Relaxed), 1);
        assert_eq!(runner.pending_task_count(), 0);
        assert!(runner.pending_messages().is_empty());
    }

    #[test]
    fn task_runner_clears_pending_and_returns_failure_result_after_panic() {
        let runner = TaskRunner::new();
        let applied = Arc::new(AtomicUsize::new(0));

        runner.spawn(PanicTask {
            applied: applied.clone(),
        });

        assert_eq!(runner.pending_task_count(), 1);
        assert_eq!(runner.pending_messages(), vec!["Panicking task..."]);

        let deadline = Instant::now() + Duration::from_secs(1);
        let result = loop {
            if let Some(result) = runner.poll() {
                break result;
            }
            assert!(Instant::now() < deadline, "panic result did not arrive");
            std::thread::sleep(Duration::from_millis(5));
        };

        let mut kernel = AppKernel::new();
        result.apply(&mut kernel);

        assert_eq!(applied.load(Ordering::Relaxed), 1);
        assert_eq!(runner.pending_task_count(), 0);
        assert!(runner.pending_messages().is_empty());
    }
}
