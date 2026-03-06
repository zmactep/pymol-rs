//! cdylib-Safe egui Widgets
//!
//! This plugin is compiled as a `cdylib`, which gets its own copy of egui
//! statics. Standard widgets like `Label` and `Button` use `TypeId`-based
//! context lookups that panic when the `TypeId` differs from the host's.
//!
//! These helpers use low-level painter APIs (`allocate_exact_size` +
//! `Painter`) to render text and buttons without any `TypeId` dependency.

use egui::{Color32, CornerRadius, FontId, Pos2, Sense, Ui, Vec2};

/// Paint a clickable button. Returns `true` if clicked.
pub fn painted_button(ui: &mut Ui, text: &str, text_color: Color32) -> bool {
    let font = FontId::proportional(13.0);
    let galley = ui.painter().layout_no_wrap(text.to_string(), font, text_color);
    let padding = Vec2::new(8.0, 4.0);
    let desired = galley.size() + padding * 2.0;

    let (rect, response) = ui.allocate_exact_size(desired, Sense::click());

    if ui.is_rect_visible(rect) {
        let bg = if response.hovered() {
            Color32::from_gray(60)
        } else {
            Color32::from_gray(45)
        };
        ui.painter()
            .rect_filled(rect, CornerRadius::same(4), bg);
        ui.painter().galley(
            Pos2::new(rect.min.x + padding.x, rect.min.y + padding.y),
            galley,
            text_color,
        );
    }

    response.clicked()
}

/// Paint a non-interactive text label.
pub fn painted_text(ui: &mut Ui, text: &str, font: FontId, color: Color32) {
    let galley = ui.painter().layout_no_wrap(text.to_string(), font, color);
    let size = galley.size();
    let (rect, _) = ui.allocate_exact_size(size, Sense::hover());
    if ui.is_rect_visible(rect) {
        ui.painter().galley(rect.min, galley, color);
    }
}
