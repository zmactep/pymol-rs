//! UI section helpers for the ray tracing toolbar.

use pymol_plugin::prelude::MessageBus;

use super::{background_label, trace_mode_label, ResolutionPreset, RtToolbarComponent};

/// Display the preview texture.
pub(crate) fn show_preview(comp: &RtToolbarComponent, ui: &mut egui::Ui) {
    if let Some(tex) = &comp.preview_texture {
        egui::Frame::new()
            .stroke(ui.visuals().widgets.noninteractive.bg_stroke)
            .corner_radius(2.0)
            .show(ui, |ui| {
                let w = ui.available_width();
                let h = w / (640.0 / 480.0);
                ui.image(egui::load::SizedTexture::new(tex.id(), egui::vec2(w, h)));
            });
        ui.add_space(4.0);
    }
}

/// "Use custom lighting" checkbox. Returns true if changed.
pub(crate) fn show_custom_lighting_toggle(
    comp: &mut RtToolbarComponent,
    ui: &mut egui::Ui,
    bus: &mut MessageBus,
) -> bool {
    if ui
        .checkbox(&mut comp.use_custom, "Use custom lighting")
        .changed()
    {
        bus.execute_command_silent(format!("set rt_use_custom, {}", comp.use_custom as i32));
        return true;
    }
    false
}

/// Resolution preset selector + width/height drag values.
pub(crate) fn show_resolution(comp: &mut RtToolbarComponent, ui: &mut egui::Ui) {
    ui.label("Resolution");
    ui.horizontal(|ui| {
        egui::ComboBox::from_id_salt("rt_preset")
            .selected_text(comp.preset.label())
            .width(70.0)
            .show_ui(ui, |ui| {
                for p in ResolutionPreset::ALL {
                    if ui
                        .selectable_value(&mut comp.preset, p, p.label())
                        .changed()
                    {
                        if let Some((w, h)) = p.dimensions() {
                            comp.width = w;
                            comp.height = h;
                        }
                    }
                }
            });

        let mut w = comp.width as f64;
        let mut h = comp.height as f64;
        ui.add(
            egui::DragValue::new(&mut w)
                .range(64.0..=7680.0)
                .speed(1.0),
        );
        ui.label("\u{00d7}"); // ×
        ui.add(
            egui::DragValue::new(&mut h)
                .range(64.0..=4320.0)
                .speed(1.0),
        );
        let new_w = w as u32;
        let new_h = h as u32;
        if new_w != comp.width || new_h != comp.height {
            comp.width = new_w;
            comp.height = new_h;
            comp.preset = ResolutionPreset::from_dimensions(new_w, new_h);
        }
    });
}

/// Antialias slider.
pub(crate) fn show_antialias(comp: &mut RtToolbarComponent, ui: &mut egui::Ui) {
    ui.horizontal(|ui| {
        ui.label("Antialias");
        ui.add(egui::Slider::new(&mut comp.antialias, 1..=4));
    });
}

/// Trace mode dropdown. Returns true if changed.
pub(crate) fn show_mode(
    comp: &mut RtToolbarComponent,
    ui: &mut egui::Ui,
    bus: &mut MessageBus,
) -> bool {
    let mut changed = false;
    ui.horizontal(|ui| {
        ui.label("Mode");
        let prev_mode = comp.mode;
        egui::ComboBox::from_id_salt("rt_mode")
            .selected_text(trace_mode_label(comp.mode))
            .show_ui(ui, |ui| {
                for mode in 0..=3 {
                    ui.selectable_value(&mut comp.mode, mode, trace_mode_label(mode));
                }
            });
        if comp.mode != prev_mode {
            bus.execute_command_silent(format!("set ray_trace_mode, {}", comp.mode));
            changed = true;
        }
    });
    changed
}

/// Lighting sliders (enabled when use_custom is on). Returns true if changed.
pub(crate) fn show_lighting(
    comp: &mut RtToolbarComponent,
    ui: &mut egui::Ui,
    bus: &mut MessageBus,
) -> bool {
    let mut changed = false;
    ui.add_enabled_ui(comp.use_custom, |ui| {
        egui::CollapsingHeader::new("Lighting")
            .default_open(false)
            .show(ui, |ui| {
                if ui
                    .add(egui::Slider::new(&mut comp.ambient, 0.0..=1.0).text("Ambient"))
                    .changed()
                {
                    bus.execute_command_silent(format!("set rt_ambient, {}", comp.ambient));
                    changed = true;
                }
                if ui
                    .add(egui::Slider::new(&mut comp.direct, 0.0..=1.0).text("Direct"))
                    .changed()
                {
                    bus.execute_command_silent(format!("set rt_direct, {}", comp.direct));
                    changed = true;
                }
                if ui
                    .add(egui::Slider::new(&mut comp.reflect, 0.0..=1.0).text("Reflect"))
                    .changed()
                {
                    bus.execute_command_silent(format!("set rt_reflect, {}", comp.reflect));
                    changed = true;
                }
                if ui
                    .add(egui::Slider::new(&mut comp.specular, 0.0..=1.0).text("Specular"))
                    .changed()
                {
                    bus.execute_command_silent(format!("set rt_specular, {}", comp.specular));
                    changed = true;
                }
                if ui
                    .add(
                        egui::Slider::new(&mut comp.shininess, 1.0..=128.0).text("Shininess"),
                    )
                    .changed()
                {
                    bus.execute_command_silent(format!("set rt_shininess, {}", comp.shininess));
                    changed = true;
                }
            });
    });
    changed
}

/// Shadow settings. Returns true if changed.
pub(crate) fn show_shadows(
    comp: &mut RtToolbarComponent,
    ui: &mut egui::Ui,
    bus: &mut MessageBus,
) -> bool {
    let mut changed = false;
    egui::CollapsingHeader::new("Shadows")
        .default_open(false)
        .show(ui, |ui| {
            if ui.checkbox(&mut comp.shadow, "Shadows").changed() {
                bus.execute_command_silent(format!("set ray_shadow, {}", comp.shadow as i32));
                changed = true;
            }
            if ui
                .checkbox(&mut comp.transparency_shadows, "Transparency shadows")
                .changed()
            {
                bus.execute_command_silent(format!(
                    "set ray_transparency_shadows, {}",
                    comp.transparency_shadows as i32
                ));
                changed = true;
            }
            if ui
                .add(egui::Slider::new(&mut comp.max_passes, 1..=100).text("Max passes"))
                .changed()
            {
                bus.execute_command_silent(format!("set ray_max_passes, {}", comp.max_passes));
                changed = true;
            }
        });
    changed
}

/// Edge detection settings. Returns true if changed.
pub(crate) fn show_edge_detection(
    comp: &mut RtToolbarComponent,
    ui: &mut egui::Ui,
    bus: &mut MessageBus,
) -> bool {
    let mut changed = false;
    egui::CollapsingHeader::new("Edge Detection")
        .default_open(false)
        .show(ui, |ui| {
            if ui
                .add(egui::Slider::new(&mut comp.slope_factor, 0.0..=2.0).text("Slope"))
                .changed()
            {
                bus.execute_command_silent(format!(
                    "set ray_trace_slope_factor, {}",
                    comp.slope_factor
                ));
                changed = true;
            }
            if ui
                .add(egui::Slider::new(&mut comp.depth_factor, 0.0..=1.0).text("Depth"))
                .changed()
            {
                bus.execute_command_silent(format!(
                    "set ray_trace_depth_factor, {}",
                    comp.depth_factor
                ));
                changed = true;
            }
            if ui
                .add(egui::Slider::new(&mut comp.disco_factor, 0.0..=1.0).text("Disco"))
                .changed()
            {
                bus.execute_command_silent(format!(
                    "set ray_trace_disco_factor, {}",
                    comp.disco_factor
                ));
                changed = true;
            }
            if ui
                .add(egui::Slider::new(&mut comp.gain, 0.0..=1.0).text("Gain"))
                .changed()
            {
                bus.execute_command_silent(format!("set ray_trace_gain, {}", comp.gain));
                changed = true;
            }
        });
    changed
}

/// Other settings (background, fog). Returns true if changed.
pub(crate) fn show_other(
    comp: &mut RtToolbarComponent,
    ui: &mut egui::Ui,
    bus: &mut MessageBus,
) -> bool {
    let mut changed = false;
    egui::CollapsingHeader::new("Other")
        .default_open(false)
        .show(ui, |ui| {
            ui.horizontal(|ui| {
                ui.label("Background");
                let prev_bg = comp.opaque_background;
                egui::ComboBox::from_id_salt("rt_bg")
                    .selected_text(background_label(comp.opaque_background))
                    .show_ui(ui, |ui| {
                        ui.selectable_value(&mut comp.opaque_background, -1, "Auto");
                        ui.selectable_value(&mut comp.opaque_background, 0, "Transparent");
                        ui.selectable_value(&mut comp.opaque_background, 1, "Opaque");
                    });
                if comp.opaque_background != prev_bg {
                    bus.execute_command_silent(format!(
                        "set ray_opaque_background, {}",
                        comp.opaque_background
                    ));
                    changed = true;
                }
            });
            if ui
                .add(egui::Slider::new(&mut comp.fog, -1.0..=1.0).text("Fog"))
                .changed()
            {
                bus.execute_command_silent(format!("set ray_trace_fog, {}", comp.fog));
                changed = true;
            }
        });
    changed
}

/// Render and Save As buttons.
pub(crate) fn show_action_buttons(
    comp: &RtToolbarComponent,
    ui: &mut egui::Ui,
    bus: &mut MessageBus,
) {
    ui.horizontal(|ui| {
        if ui.button("Render").clicked() {
            comp.send_render_commands(bus, None);
        }
        if ui.button("Save As\u{2026}").clicked() {
            let dialog = rfd::FileDialog::new()
                .set_title("Save Ray Trace Image")
                .add_filter("PNG Image", &["png"])
                .set_file_name("raytrace.png");
            if let Some(path) = dialog.save_file() {
                comp.send_render_commands(bus, Some(&path.display().to_string()));
            }
        }
    });
}
