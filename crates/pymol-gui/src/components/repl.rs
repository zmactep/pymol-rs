//! REPL Component (Read-Eval-Print Loop)
//!
//! Combines the output log and command line into a single panel.
//! Output is displayed above, command input below — always shown together.

use pymol_framework::component::{Component, SharedContext};
use pymol_framework::message::{AppMessage, MessageBus};
use crate::model::{CommandLineModel, OutputModel};
use crate::ui::command::CommandLineUiState;
use crate::ui::completion::CompletionContext;
use crate::ui::{CommandLinePanel, OutputPanel};

/// Self-contained REPL component: output log + command line input.
pub struct ReplComponent {
    pub output: OutputModel,
    pub command_line: CommandLineModel,
    pub command_line_ui: CommandLineUiState,
}

impl ReplComponent {
    pub fn new() -> Self {
        Self {
            output: OutputModel::new(),
            command_line: CommandLineModel::new(),
            command_line_ui: CommandLineUiState::new(),
        }
    }
}

impl Default for ReplComponent {
    fn default() -> Self {
        Self::new()
    }
}

impl Component for ReplComponent {
    fn id(&self) -> &'static str {
        "repl"
    }

    fn title(&self) -> &str {
        "REPL"
    }

    fn show(&mut self, ui: &mut egui::Ui, ctx: &SharedContext, bus: &mut MessageBus) {
        // Output log fills available space minus the command line row
        let available = ui.available_height();
        let cmd_line_height = 28.0;
        let output_height = (available - cmd_line_height - 2.0).max(0.0);

        // Fixed-height region for output so it doesn't resize with content
        ui.allocate_ui_with_layout(
            egui::vec2(ui.available_width(), output_height),
            *ui.layout(),
            |ui| {
                OutputPanel::show(ui, &self.output);
            },
        );

        ui.add_space(2.0);

        // Command line input at the bottom
        let command_name_refs: Vec<&str> =
            ctx.command_names.iter().map(|s| s.as_str()).collect();
        let color_names: Vec<String> = ctx
            .named_colors
            .names()
            .iter()
            .map(|s| s.to_string())
            .collect();
        let object_names: Vec<String> =
            ctx.registry.names().map(|s| s.to_string()).collect();
        let selection_names: Vec<String> = ctx.selections.names();

        let completion_ctx = CompletionContext {
            command_names: &command_name_refs,
            registry: ctx.command_registry,
            setting_names: ctx.setting_names,
            color_names: &color_names,
            object_names: &object_names,
            selection_names: &selection_names,
        };

        CommandLinePanel::show(
            ui,
            &mut self.command_line,
            &mut self.command_line_ui,
            &completion_ctx,
            bus,
        );
    }

    fn on_message(&mut self, msg: &AppMessage) {
        match msg {
            AppMessage::PrintInfo(text) => self.output.print_info(text),
            AppMessage::PrintWarning(text) => self.output.print_warning(text),
            AppMessage::PrintError(text) => self.output.print_error(text),
            AppMessage::PrintCommand(text) => self.output.print_command(text),
            _ => {}
        }
    }
}
