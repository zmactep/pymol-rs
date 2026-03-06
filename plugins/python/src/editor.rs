//! Python Script Editor Component
//!
//! A GUI panel with a multiline code editor (syntax-highlighted),
//! a Run button, and an output area for execution results.
//!
//! NOTE: This component runs inside a cdylib plugin, which has its own copy
//! of egui statics. Widgets that use TypeId-based context lookups (Label,
//! Button) will panic because their TypeIds differ from the host's.
//! We use cdylib-safe helpers from [`crate::widgets`] instead.

use std::sync::{Arc, Mutex};

use egui::{Color32, FontId, Pos2, ScrollArea, Stroke, TextBuffer, TextEdit, Ui};
use pymol_plugin::prelude::*;

use crate::engine::PythonEngine;
use crate::highlight::PythonHighlighter;
use crate::widgets::{painted_button, painted_text};

/// Output entry from script execution.
struct OutputEntry {
    text: String,
    is_error: bool,
}

/// Python script editor with syntax highlighting and execution.
pub struct PythonEditorComponent {
    code: String,
    output: Vec<OutputEntry>,
    engine: Arc<Mutex<PythonEngine>>,
    highlighter: PythonHighlighter,
}

impl PythonEditorComponent {
    pub fn new(engine: Arc<Mutex<PythonEngine>>) -> Self {
        Self {
            code: String::new(),
            output: Vec::new(),
            engine,
            highlighter: PythonHighlighter::new(),
        }
    }

    /// Execute the current editor contents and capture output.
    fn run_code(&mut self, bus: &mut MessageBus) {
        if self.code.trim().is_empty() {
            return;
        }

        let result = {
            let mut engine = self.engine.lock().unwrap();
            engine.eval(&self.code)
        };

        match result {
            Ok(output) => {
                if !output.is_empty() {
                    bus.print_info(&output);
                    self.output.push(OutputEntry {
                        text: output,
                        is_error: false,
                    });
                }
            }
            Err(err) => {
                bus.print_error(&err);
                self.output.push(OutputEntry {
                    text: err,
                    is_error: true,
                });
            }
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
    fn show_toolbar(&mut self, ui: &mut Ui, bus: &mut MessageBus) {
        ui.horizontal(|ui| {
            if painted_button(ui, "▶ Run", Color32::from_rgb(150, 255, 150)) {
                self.run_code(bus);
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

    fn show(&mut self, ui: &mut Ui, _ctx: &SharedContext, bus: &mut MessageBus) {
        let available = ui.available_height();
        // Prevent the resizable panel from shrinking to fit content
        ui.set_min_height(available);

        let toolbar_height = 28.0;
        let output_height = self.output_height(available);
        let editor_height = (available - toolbar_height - output_height - 4.0).max(40.0);

        self.show_toolbar(ui, bus);
        ui.add_space(2.0);

        // Consume Cmd+Enter / Ctrl+Enter BEFORE TextEdit sees it
        let run_shortcut = ui
            .ctx()
            .input_mut(|i| i.consume_key(egui::Modifiers::COMMAND, egui::Key::Enter));
        if run_shortcut {
            self.run_code(bus);
        }

        self.show_editor(ui, editor_height);
        self.show_output(ui, output_height);
    }
}
