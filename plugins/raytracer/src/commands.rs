//! Command implementations for the raytracer plugin.
//!
//! Provides the `ray` command (perform raytracing).

use pymol_plugin::prelude::*;
use pymol_scene::ViewportImage;

use crate::scene::raytrace_scene;
use crate::settings::read_ray_settings;

// ---------------------------------------------------------------------------
// RayCommand
// ---------------------------------------------------------------------------

pub(crate) struct RayCommand;

impl Command for RayCommand {
    fn name(&self) -> &str {
        "ray"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "ray" performs ray-tracing and saves the resulting image to a file.
    Ray tracing produces high-quality images with proper shadows, lighting,
    and transparency effects.

USAGE

    ray [ width [, height [, antialias [, filename [, quiet ]]]]]

ARGUMENTS

    width = integer: width in pixels (default: current window width)
    height = integer: height in pixels (default: current window height)
    antialias = integer: antialiasing level 1-4 (default: antialias setting)
        1 = no antialiasing
        2 = 2x2 supersampling
        3 = 3x3 supersampling
        4 = 4x4 supersampling
    filename = string: output file path (default: no file saved, returns to display)
    quiet = 0/1: suppress feedback (default: 0)

EXAMPLES

    ray                          # Raytrace at current resolution
    ray 1920, 1080               # Raytrace at 1080p
    ray 1920, 1080, 2            # Raytrace at 1080p with 2x2 AA
    ray width=1920, height=1080, filename=output.png

SEE ALSO

    png
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        use std::time::Instant;

        let width = args
            .get_int(0)
            .or_else(|| args.get_named_int("width"))
            .map(|v| v as u32);

        let height = args
            .get_int(1)
            .or_else(|| args.get_named_int("height"))
            .map(|v| v as u32);

        let antialias = args
            .get_int(2)
            .or_else(|| args.get_named_int("antialias"))
            .unwrap_or_else(|| ctx.viewer.settings().ui.antialias as i64) as u32;

        let filename = args
            .get_str(3)
            .or_else(|| args.get_named_str("filename"));

        let quiet = args
            .get_bool(4)
            .or_else(|| args.get_named_bool("quiet"))
            .unwrap_or(false);

        // Ensure representations are built (mutable borrow, released immediately)
        ctx.viewer.prepare_render();

        // Read ray settings from plugin's dynamic settings
        let ray_settings = read_ray_settings(|name| {
            ctx.dynamic_setting(name)
                .and_then(|e| e.store.read().ok())
                .and_then(|s| s.get(name).cloned())
        });

        // Determine output dimensions
        let viewport = ctx.viewer.viewport_size();
        let final_width = width.unwrap_or(viewport.0.max(1024));
        let final_height = height.unwrap_or(viewport.1.max(768));

        let start = Instant::now();

        // Perform raytracing (handles borrow-checker dance internally)
        let (image_data, w, h) = raytrace_scene(
            ctx.viewer,
            &ray_settings,
            Some(final_width),
            Some(final_height),
            antialias,
        )
        .map_err(|e| CmdError::execution(format!("Ray tracing failed: {e}")))?;

        let elapsed = start.elapsed();

        if let Some(path) = filename {
            let path = expand_path(path);
            let path = if path.extension().map(|e| e.to_ascii_lowercase()) != Some("png".into()) {
                path.with_extension("png")
            } else {
                path.to_path_buf()
            };

            image::save_buffer(&path, &image_data, w, h, image::ColorType::Rgba8)
                .map_err(|e| CmdError::execution(format!("Failed to save PNG: {e}")))?;

            if !quiet {
                ctx.print(&format!(
                    " Ray: render time {:02}:{:02}:{:02}.{:02}  ({}x{})",
                    elapsed.as_secs() / 3600,
                    (elapsed.as_secs() % 3600) / 60,
                    elapsed.as_secs() % 60,
                    elapsed.subsec_millis() / 10,
                    w,
                    h
                ));
                ctx.print(&format!(" Saved \"{}\"", path.display()));
            }
        } else {
            ctx.viewer.set_viewport_image(Some(ViewportImage {
                data: image_data,
                width: w,
                height: h,
            }));

            if !quiet {
                ctx.print(&format!(
                    " Ray: render time {:02}:{:02}:{:02}.{:02}  ({}x{})",
                    elapsed.as_secs() / 3600,
                    (elapsed.as_secs() % 3600) / 60,
                    elapsed.as_secs() % 60,
                    elapsed.subsec_millis() / 10,
                    w,
                    h
                ));
                ctx.print(
                    " Ray trace complete. Use 'png filename' to save, or interact to dismiss.",
                );
            }
        }

        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Helper: expand ~ in paths
// ---------------------------------------------------------------------------

fn expand_path(path: &str) -> std::path::PathBuf {
    if path.starts_with('~') {
        if let Some(home) = dirs_hint() {
            return std::path::PathBuf::from(path.replacen('~', &home, 1));
        }
    }
    std::path::PathBuf::from(path)
}

fn dirs_hint() -> Option<String> {
    std::env::var("HOME").ok()
}
