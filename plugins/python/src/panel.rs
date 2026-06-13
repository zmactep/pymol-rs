//! Declarative scripting panel for the embedded Python plugin.

use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use patinae_plugin::prelude::*;

use crate::highlight::PythonHighlightCache;
use crate::worker::{WorkItem, WorkOrigin, WorkerHandle};

const MAX_OUTPUT_CHARS: usize = 24_000;

pub(crate) type ScriptPanelStateHandle = Arc<Mutex<ScriptPanelState>>;

#[derive(Debug, Default)]
pub(crate) struct ScriptPanelState {
    input: String,
    output: String,
    started_at: Option<Instant>,
    stop_requested: bool,
    highlight_cache: PythonHighlightCache,
}

impl ScriptPanelState {
    fn edit_input(&mut self, input: String) {
        self.input = input;
    }

    pub(crate) fn append_output(&mut self, text: &str) {
        self.output.push_str(text);
        if !self.output.ends_with('\n') {
            self.output.push('\n');
        }
        self.trim_output();
    }

    pub(crate) fn append_error(&mut self, text: &str) {
        self.output.push_str("Error: ");
        self.output.push_str(text);
        if !self.output.ends_with('\n') {
            self.output.push('\n');
        }
        self.trim_output();
    }

    fn clear_output(&mut self) {
        self.output.clear();
    }

    fn start_run(&mut self) {
        self.stop_requested = false;
        self.started_at = Some(Instant::now());
    }

    fn request_stop(&mut self) -> bool {
        if self.stop_requested {
            false
        } else {
            self.stop_requested = true;
            true
        }
    }

    pub(crate) fn finish_run(&mut self) -> Option<String> {
        self.stop_requested = false;
        let started_at = self.started_at.take()?;

        Some(format_run_duration(started_at.elapsed()))
    }

    fn input_highlights(&mut self) -> Vec<PanelTextHighlight> {
        self.highlight_cache.highlights_for(&self.input)
    }

    fn trim_output(&mut self) {
        if self.output.len() <= MAX_OUTPUT_CHARS {
            return;
        }

        let mut start = self.output.len() - MAX_OUTPUT_CHARS;
        while !self.output.is_char_boundary(start) {
            start += 1;
        }
        self.output.drain(..start);
        self.output.insert_str(0, "[output truncated]\n");
    }
}

fn format_run_duration(duration: Duration) -> String {
    let seconds = duration.as_secs_f64();
    if seconds < 0.001 {
        format!("{:.0} us", seconds * 1_000_000.0)
    } else if seconds < 1.0 {
        format!("{:.1} ms", seconds * 1_000.0)
    } else {
        format!("{seconds:.3} s")
    }
}

pub(crate) fn shared_panel_state() -> ScriptPanelStateHandle {
    Arc::new(Mutex::new(ScriptPanelState::default()))
}

pub(crate) struct PythonScriptPanel {
    worker: WorkerHandle,
    state: ScriptPanelStateHandle,
}

impl PythonScriptPanel {
    pub(crate) fn new(worker: WorkerHandle, state: ScriptPanelStateHandle) -> Self {
        Self { worker, state }
    }

    fn run_script(&self) {
        let code = {
            let state = self.state.lock().unwrap();
            state.input.clone()
        };

        if code.trim().is_empty() {
            self.state
                .lock()
                .unwrap()
                .append_output("No script to run.");
            return;
        }

        if self.worker.is_busy() {
            self.state
                .lock()
                .unwrap()
                .append_output("Python worker is already running.");
            return;
        }

        self.state.lock().unwrap().start_run();
        self.worker.submit(WorkItem::Eval {
            code,
            origin: WorkOrigin::Panel,
        });
    }

    fn stop_script(&self) {
        let mut state = self.state.lock().unwrap();
        if !self.worker.is_busy() {
            state.append_output("Python worker is not running.");
            return;
        }

        if state.request_stop() {
            self.worker.request_interrupt();
            state.append_output("Interrupt requested.");
        }
    }
}

impl PluginPanel for PythonScriptPanel {
    fn descriptor(&self) -> PanelDescriptor {
        PanelDescriptor::bottom("python_scripting", "Python")
            .icon("Py")
            .default_visible(false)
    }

    fn runtime_requirements(&self) -> PanelRuntimeRequirements {
        PanelRuntimeRequirements::NONE
    }

    fn snapshot(&mut self, _ctx: &SharedContext<'_>) -> PanelSnapshot {
        let mut state = self.state.lock().unwrap();
        let busy = self.worker.is_busy();
        let has_output = !state.output.is_empty();
        let stop_requested = state.stop_requested;
        let input = state.input.clone();
        let input_highlights = state.input_highlights();
        let mut controls = vec![
            PanelControl::ButtonRow {
                id: "toolbar".into(),
                buttons: vec![
                    PanelButton::new("run", "Run script", "run", true).enabled(!busy),
                    PanelButton::new("stop", "Stop", "stop", false)
                        .enabled(busy && !stop_requested),
                    PanelButton::new("clear_output", "Clear output", "clear", false)
                        .enabled(has_output),
                ],
            },
            PanelControl::Row(
                PanelRow::new(
                    "script_output",
                    vec![
                        PanelControlNode::new(PanelControl::TextArea(
                            PanelTextArea::new(
                                "script",
                                "",
                                input,
                                "print('hello from Patinae')",
                                5,
                                false,
                            )
                            .with_highlights(input_highlights),
                        ))
                        .grow(1.0),
                        PanelControlNode::new(PanelControl::TextArea(PanelTextArea::new(
                            "output",
                            "",
                            state.output.clone(),
                            "Output will appear here.",
                            5,
                            true,
                        )))
                        .grow(1.0),
                    ],
                )
                .gap(8.0),
            ),
        ];

        if busy {
            controls.push(PanelControl::Text {
                id: "status".into(),
                text: "Running Python...".into(),
            });
        }

        PanelSnapshot::new(controls)
    }

    fn handle_event(
        &mut self,
        event: PanelEvent,
        _ctx: &SharedContext<'_>,
        _bus: &mut MessageBus,
    ) -> Vec<PanelAction> {
        match event.control_id.as_str() {
            "script" => {
                if let PanelValue::Text(text) = event.value {
                    self.state.lock().unwrap().edit_input(text);
                }
            }
            "run" => self.run_script(),
            "stop" => self.stop_script(),
            "clear_output" => self.state.lock().unwrap().clear_output(),
            _ => {}
        }
        Vec::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::AtomicBool;
    use std::sync::Arc;

    #[test]
    fn python_panel_uses_lightweight_runtime_input() {
        let (worker, _rx) = crate::worker::spawn_worker(Arc::new(AtomicBool::new(false)));
        let panel = PythonScriptPanel::new(worker, shared_panel_state());

        assert!(panel.runtime_requirements().is_empty());
    }

    #[test]
    fn edit_input_updates_cached_highlights_source() {
        let mut state = ScriptPanelState::default();
        state.edit_input("def f():\n    return 'hi'\n".into());

        let first = state.input_highlights();
        let second = state.input_highlights();

        assert!(!first.is_empty());
        assert_eq!(first, second);
    }

    #[test]
    fn output_can_be_cleared() {
        let mut state = ScriptPanelState::default();

        state.append_output("hello");
        assert!(!state.output.is_empty());

        state.clear_output();
        assert!(state.output.is_empty());
    }

    #[test]
    fn run_finish_returns_duration_and_resets_stop_request() {
        let mut state = ScriptPanelState::default();

        state.start_run();
        assert!(state.request_stop());

        let duration = state.finish_run();

        assert!(duration.is_some());
        assert!(!state.stop_requested);
    }

    #[test]
    fn stop_request_is_idempotent_until_finish() {
        let mut state = ScriptPanelState::default();

        assert!(state.request_stop());
        assert!(!state.request_stop());

        state.finish_run();
        assert!(state.request_stop());
    }
}
