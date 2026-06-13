use std::panic::{catch_unwind, AssertUnwindSafe};

use patinae_framework::component::SharedContext;
use patinae_framework::message::{AppMessage, MessageBus};
use patinae_framework::plugin_ui::PanelAction;

use crate::host::PluginHost;
use crate::panic::panic_payload_to_string;

impl PluginHost {
    pub(crate) fn process_panel_events(
        &mut self,
        shared: &SharedContext<'_>,
        bus: &mut MessageBus,
    ) {
        let events = std::mem::take(&mut self.pending_panel_events);
        for event in events {
            let refresh = event.kind.refreshes_snapshot();
            let Some((plugin_idx, panel_idx)) = self.find_panel_indices(&event.panel_id) else {
                log::warn!("Unknown plugin panel event target '{}'", event.panel_id);
                continue;
            };
            if self.plugins[plugin_idx].faulted {
                continue;
            }
            let plugin_name = self.plugins[plugin_idx].metadata.name.clone();
            let result = {
                let panel = &mut self.plugins[plugin_idx].panels[panel_idx].panel;
                catch_unwind(AssertUnwindSafe(|| panel.handle_event(event, shared, bus)))
            };
            match result {
                Ok(actions) => {
                    apply_panel_actions(actions, bus);
                    if refresh {
                        self.bump_panel_ui_generation();
                    }
                }
                Err(panic_info) => {
                    log::error!(
                        "Plugin '{}' panicked during panel event: {}. Plugin disabled.",
                        plugin_name,
                        panic_payload_to_string(&panic_info),
                    );
                    self.plugins[plugin_idx].faulted = true;
                }
            }
        }
    }
}

fn apply_panel_actions(actions: Vec<PanelAction>, bus: &mut MessageBus) {
    for action in actions {
        match action {
            PanelAction::ExecuteCommand { command, silent } => {
                bus.send(AppMessage::ExecuteCommand { command, silent });
            }
            PanelAction::SetSetting { name, value } => {
                bus.execute_command_silent(format!("set {}, {}", name, value.as_command_value()));
            }
            PanelAction::Custom { topic, payload } => {
                bus.send(AppMessage::Custom { topic, payload });
            }
        }
    }
}
