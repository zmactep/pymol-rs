//! Preview rendering and scene change debouncing.

use pymol_plugin::prelude::SharedContext;

use crate::scene::raytrace_preview;
use crate::settings::read_ray_settings;

use super::RtToolbarComponent;

/// Number of frames to wait after the last scene change before re-rendering.
const DEBOUNCE_FRAMES: u32 = 10; // ~170ms at 60fps

/// Track scene generation changes and mark preview dirty after debounce.
pub(crate) fn update_scene_debounce(
    comp: &mut RtToolbarComponent,
    ui: &egui::Ui,
    ctx: &SharedContext,
) {
    if ctx.scene_generation != comp.last_scene_gen {
        comp.last_scene_gen = ctx.scene_generation;
        comp.stable_frames = 0;
    } else if comp.stable_frames < DEBOUNCE_FRAMES {
        comp.stable_frames += 1;
        if comp.stable_frames == DEBOUNCE_FRAMES {
            comp.preview_dirty = true;
        }
    }

    if comp.stable_frames < DEBOUNCE_FRAMES {
        ui.ctx().request_repaint();
    }
}

/// Render a 640x480 preview and update the texture.
pub(crate) fn render_preview(
    comp: &mut RtToolbarComponent,
    egui_ctx: &egui::Context,
    ctx: &SharedContext,
) {
    let (device, queue) = match (ctx.gpu_device, ctx.gpu_queue) {
        (Some(d), Some(q)) => (d, q),
        _ => return,
    };

    let ray_settings = read_ray_settings(|name| {
        ctx.dynamic_settings
            .and_then(|r| r.lookup(name))
            .and_then(|e| e.store.read().ok())
            .and_then(|s| s.get(name).cloned())
    });

    let result = raytrace_preview(
        ctx.registry,
        ctx.camera,
        ctx.settings,
        ctx.named_colors,
        ctx.clear_color,
        &ray_settings,
        device,
        queue,
        640,
        480,
    );

    match result {
        Ok((data, w, h)) => {
            let color_image =
                egui::ColorImage::from_rgba_unmultiplied([w as usize, h as usize], &data);
            match &mut comp.preview_texture {
                Some(tex) => tex.set(color_image, egui::TextureOptions::LINEAR),
                None => {
                    comp.preview_texture = Some(egui_ctx.load_texture(
                        "rt_preview",
                        color_image,
                        egui::TextureOptions::LINEAR,
                    ));
                }
            }
        }
        Err(e) => {
            log::debug!("preview: {e}");
        }
    }

    comp.preview_dirty = false;
}
