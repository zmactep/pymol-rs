//! Python Script Editor Component
//!
//! A GUI panel with a multiline code editor (syntax-highlighted),
//! a Run button, and an output area for execution results.
//!
//! Execution is non-blocking: code is submitted to the Python worker
//! thread, and results are received via a shared output buffer.
//!
//! NOTE: This component runs inside a cdylib plugin, which has its own copy
//! of egui statics. Widgets that use TypeId-based context lookups (Label,
//! Button) will panic because their TypeIds differ from the host's.
//! We use cdylib-safe helpers from [`crate::widgets`] instead.

use egui::{Color32, FontId, Pos2, ScrollArea, Stroke, TextBuffer, TextEdit, Ui};
use pymol_plugin::prelude::*;

use crate::highlight::PythonHighlighter;
use crate::widgets::{painted_button, painted_text};
use crate::worker::{EditorOutputHandle, OutputEntry, WorkItem, WorkOrigin, WorkerHandle};

/// Python script editor with syntax highlighting and execution.
pub struct PythonEditorComponent {
    code: String,
    output: Vec<OutputEntry>,
    worker: WorkerHandle,
    editor_output: EditorOutputHandle,
    highlighter: PythonHighlighter,
}

impl PythonEditorComponent {
    pub fn new(worker: WorkerHandle, editor_output: EditorOutputHandle) -> Self {
        Self {
            code: String::new(),
            output: Vec::new(),
            worker,
            editor_output,
            highlighter: PythonHighlighter::new(),
        }
    }

    /// Submit the current editor contents for execution.
    fn run_code(&mut self) {
        if self.code.trim().is_empty() || self.worker.is_busy() {
            return;
        }

        self.worker.submit(WorkItem::Eval {
            code: self.code.clone(),
            origin: WorkOrigin::Editor,
        });
    }

    /// Drain any results that arrived from the worker thread.
    fn drain_results(&mut self) {
        if let Ok(mut buf) = self.editor_output.try_lock() {
            self.output.append(&mut buf);
        }
    }

    /// Calculate the height available for the output area.
    fn output_height(&self, available: f32) -> f32 {
        if self.output.is_empty() {
            0.0
        } else {
            (available * 0.3).clamp(40.0, 150.0)
        }
    }

    /// Render the toolbar row (Run button + shortcut hint + Clear).
    fn show_toolbar(&mut self, ui: &mut Ui) {
        ui.horizontal(|ui| {
            let busy = self.worker.is_busy();

            if busy {
                if painted_button(ui, "\u{25a0} Stop", Color32::from_rgb(255, 100, 100)) {
                    self.worker.request_interrupt();
                }
            } else if painted_button(ui, "\u{25b6} Run", Color32::from_rgb(150, 255, 150)) {
                self.run_code();
            }

            if !self.output.is_empty() {
                ui.add_space(4.0);
                if painted_button(ui, "Clear", Color32::GRAY) {
                    self.output.clear();
                }
            }

            ui.add_space(8.0);
            painted_text(
                ui,
                "Cmd+Enter to run",
                FontId::proportional(11.0),
                Color32::GRAY,
            );
        });
    }

    /// Render the syntax-highlighted code editor.
    fn show_editor(&mut self, ui: &mut Ui, height: f32) {
        let highlighter = &self.highlighter;
        let mut layouter = |ui: &Ui, text: &dyn TextBuffer, wrap_width: f32| {
            let font_id = FontId::monospace(13.0);
            let job = highlighter.highlight(text.as_str(), font_id, wrap_width);
            ui.painter().layout_job(job)
        };

        ScrollArea::vertical()
            .id_salt("python_editor_scroll")
            .max_height(height)
            .show(ui, |ui| {
                ui.add_sized(
                    [ui.available_width(), height],
                    TextEdit::multiline(&mut self.code)
                        .code_editor()
                        .desired_width(f32::INFINITY)
                        .layouter(&mut layouter),
                )
            });
    }

    /// Render the output area with a separator line.
    fn show_output(&self, ui: &mut Ui, height: f32) {
        if self.output.is_empty() {
            return;
        }

        ui.add_space(2.0);

        // Separator line
        let rect = ui.available_rect_before_wrap();
        let y = rect.min.y;
        ui.painter().line_segment(
            [Pos2::new(rect.min.x, y), Pos2::new(rect.max.x, y)],
            Stroke::new(1.0, Color32::from_gray(60)),
        );
        ui.add_space(3.0);

        ScrollArea::vertical()
            .id_salt("python_editor_output")
            .max_height(height)
            .stick_to_bottom(true)
            .show(ui, |ui| {
                ui.set_min_width(ui.available_width());
                let mono = FontId::monospace(12.0);
                for entry in &self.output {
                    let color = if entry.is_error {
                        Color32::from_rgb(255, 100, 100)
                    } else {
                        Color32::LIGHT_GRAY
                    };
                    for line in entry.text.lines() {
                        painted_text(ui, line, mono.clone(), color);
                    }
                }
            });
    }
}

impl Component for PythonEditorComponent {
    fn id(&self) -> &'static str {
        "python_editor"
    }

    fn title(&self) -> &str {
        "Python"
    }

    fn show(&mut self, ui: &mut Ui, _ctx: &SharedContext, _bus: &mut MessageBus) {
        // Drain any results from the worker thread
        self.drain_results();

        let available = ui.available_height();
        // Prevent the resizable panel from shrinking to fit content
        ui.set_min_height(available);

        let toolbar_height = 28.0;
        let output_height = self.output_height(available);
        let editor_height = (available - toolbar_height - output_height - 4.0).max(40.0);

        self.show_toolbar(ui);
        ui.add_space(2.0);

        // Consume Cmd+Enter / Ctrl+Enter BEFORE TextEdit sees it
        let run_shortcut = ui
            .ctx()
            .input_mut(|i| i.consume_key(egui::Modifiers::COMMAND, egui::Key::Enter));
        if run_shortcut {
            self.run_code();
        }

        self.show_editor(ui, editor_height);
        self.show_output(ui, output_height);
    }
}
