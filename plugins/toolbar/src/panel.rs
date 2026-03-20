//! Toolbar Panel
//!
//! Rendering logic for the toolbar: groups of icon buttons with labels.

use std::collections::HashMap;

use pymol_plugin::prelude::{AppMessage, MessageBus};
use crate::model::{ToolbarAction, ToolbarButton, ToolbarGroup};

// ---------------------------------------------------------------------------
// Placeholder icon generation
// ---------------------------------------------------------------------------

/// Generate a simple placeholder icon: a colored square with a bright center.
///
/// Returns RGBA bytes for a `size x size` image.
fn generate_placeholder_icon(letter: char, color: [u8; 3], size: u32) -> Vec<u8> {
    let mut pixels = vec![0u8; (size * size * 4) as usize];
    let [r, g, b] = color;

    for y in 0..size {
        for x in 0..size {
            let idx = ((y * size + x) * 4) as usize;
            pixels[idx] = r;
            pixels[idx + 1] = g;
            pixels[idx + 2] = b;
            pixels[idx + 3] = 200;

            let cx = size / 2;
            let cy = size / 2;
            let dx = (x as i32 - cx as i32).unsigned_abs();
            let dy = (y as i32 - cy as i32).unsigned_abs();
            if dx < size / 5 && dy < size / 4 {
                pixels[idx] = 255;
                pixels[idx + 1] = 255;
                pixels[idx + 2] = 255;
                pixels[idx + 3] = 220;
            }
        }
    }

    let _ = letter;
    pixels
}

// ---------------------------------------------------------------------------
// Rendering
// ---------------------------------------------------------------------------

const ICON_SIZE: f32 = 32.0;
const BUTTON_SPACING: f32 = 4.0;
const GROUP_SPACING: f32 = 8.0;
/// Fixed height for button label area (fits up to 2 lines at 9pt).
const LABEL_HEIGHT: f32 = 26.0;
/// Fixed height for the group name row at 10pt.
const GROUP_LABEL_HEIGHT: f32 = 14.0;

/// Toolbar panel rendering.
pub struct ToolbarPanel;

impl ToolbarPanel {
    /// Render the toolbar into the given UI region.
    pub fn show(
        ui: &mut egui::Ui,
        groups: &[ToolbarGroup],
        textures: &mut HashMap<String, egui::TextureHandle>,
        bus: &mut MessageBus,
    ) {
        // Ensure all icons are loaded as textures
        for group in groups {
            for button in &group.buttons {
                if !textures.contains_key(button.icon_id) {
                    let texture = Self::load_icon_texture(ui.ctx(), button);
                    textures.insert(button.icon_id.to_string(), texture);
                }
            }
        }

        egui::ScrollArea::horizontal()
            .id_salt("toolbar_scroll")
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.spacing_mut().item_spacing.x = GROUP_SPACING;

                    for (i, group) in groups.iter().enumerate() {
                        if i > 0 {
                            ui.separator();
                        }
                        Self::show_group(ui, group, textures, bus);
                    }
                });
            });
    }

    fn show_group(
        ui: &mut egui::Ui,
        group: &ToolbarGroup,
        textures: &HashMap<String, egui::TextureHandle>,
        bus: &mut MessageBus,
    ) {
        let button_col_width = ICON_SIZE + 8.0;
        let n = group.buttons.len() as f32;
        let buttons_width = n * button_col_width + (n - 1.0).max(0.0) * BUTTON_SPACING;
        // Ensure group is wide enough for its label text.
        let label_width = {
            let font = egui::FontId::proportional(10.0);
            let galley = ui.painter().layout_no_wrap(group.name.to_string(), font, egui::Color32::WHITE);
            galley.size().x
        };
        let group_width = buttons_width.max(label_width + 4.0);
        let group_height = ICON_SIZE + LABEL_HEIGHT + GROUP_LABEL_HEIGHT;

        // Reserve a fixed rect for the entire group.
        let (group_rect, _) = ui.allocate_exact_size(
            egui::vec2(group_width, group_height),
            egui::Sense::hover(),
        );

        // Buttons row: placed at the top of the group rect.
        let buttons_rect = egui::Rect::from_min_size(
            group_rect.min,
            egui::vec2(group_width, ICON_SIZE + LABEL_HEIGHT),
        );
        let mut buttons_ui = ui.new_child(
            egui::UiBuilder::new()
                .max_rect(buttons_rect)
                .layout(egui::Layout::left_to_right(egui::Align::TOP)),
        );
        buttons_ui.spacing_mut().item_spacing.x = BUTTON_SPACING;
        for button in &group.buttons {
            Self::show_button(&mut buttons_ui, button, textures, bus);
        }

        // Group label: pinned to bottom of group rect, centered.
        let label_rect = egui::Rect::from_min_size(
            egui::pos2(group_rect.min.x, group_rect.max.y - GROUP_LABEL_HEIGHT),
            egui::vec2(group_width, GROUP_LABEL_HEIGHT),
        );
        ui.put(
            label_rect,
            egui::Label::new(
                egui::RichText::new(group.name)
                    .size(10.0)
                    .color(ui.visuals().weak_text_color()),
            )
            .selectable(false),
        );
    }

    fn show_button(
        ui: &mut egui::Ui,
        button: &ToolbarButton,
        textures: &HashMap<String, egui::TextureHandle>,
        bus: &mut MessageBus,
    ) {
        let button_col_width = ICON_SIZE + 8.0;
        let col_height = ICON_SIZE + LABEL_HEIGHT;

        let response = ui.allocate_ui_with_layout(
            egui::vec2(button_col_width, col_height),
            egui::Layout::top_down(egui::Align::Center),
            |ui| {
                let clicked = if let Some(texture) = textures.get(button.icon_id) {
                    let image = egui::Image::new(texture)
                        .fit_to_exact_size(egui::vec2(ICON_SIZE, ICON_SIZE));
                    ui.add(egui::Button::image(image))
                        .on_hover_text(button.tooltip)
                        .on_hover_cursor(egui::CursorIcon::PointingHand)
                        .clicked()
                } else {
                    ui.add_sized(
                        egui::vec2(ICON_SIZE, ICON_SIZE),
                        egui::Button::new(button.label),
                    )
                    .on_hover_text(button.tooltip)
                    .clicked()
                };

                ui.add(egui::Label::new(egui::RichText::new(button.label).size(9.0)).selectable(false));

                clicked
            },
        );

        if response.inner {
            Self::execute_action(&button.action, bus);
        }
    }

    fn execute_action(action: &ToolbarAction, bus: &mut MessageBus) {
        match action {
            ToolbarAction::Commands(cmds) => {
                for cmd in cmds {
                    bus.send(AppMessage::ExecuteCommand {
                        command: cmd.clone(),
                        silent: false,
                    });
                }
            }
        }
    }

    /// Load icon bytes into an egui texture.
    fn load_icon_texture(
        ctx: &egui::Context,
        button: &ToolbarButton,
    ) -> egui::TextureHandle {
        if let Ok(img) = image::load_from_memory(button.icon_bytes) {
            let rgba = img.to_rgba8();
            let size = [rgba.width() as usize, rgba.height() as usize];
            let pixels = rgba.into_raw();
            let color_image = egui::ColorImage::from_rgba_unmultiplied(size, &pixels);
            ctx.load_texture(button.icon_id, color_image, egui::TextureOptions::LINEAR)
        } else {
            let color = Self::label_to_color(button.label);
            let pixels = generate_placeholder_icon(
                button.label.chars().next().unwrap_or('?'),
                color,
                32,
            );
            let color_image = egui::ColorImage::from_rgba_unmultiplied([32, 32], &pixels);
            ctx.load_texture(button.icon_id, color_image, egui::TextureOptions::LINEAR)
        }
    }

    /// Derive a distinct color from a label string (simple hash).
    fn label_to_color(label: &str) -> [u8; 3] {
        let hash: u32 = label.bytes().fold(0u32, |acc, b| acc.wrapping_mul(31).wrapping_add(b as u32));
        [
            (60 + (hash % 150)) as u8,
            (60 + ((hash / 150) % 150)) as u8,
            (60 + ((hash / 22500) % 150)) as u8,
        ]
    }
}
