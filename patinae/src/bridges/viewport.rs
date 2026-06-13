//! Viewport bridge: FPS tracking, texture push, frame dimensions.

use std::time::Instant;

use crate::ViewportState;

/// Tracks frame dimensions and FPS, pushes rendered images to Slint.
pub struct ViewportBridge {
    pub width: u32,
    pub height: u32,
    frame_count: u64,
    last_fps_time: Instant,
}

impl ViewportBridge {
    pub fn new(width: u32, height: u32) -> Self {
        Self {
            width,
            height,
            frame_count: 0,
            last_fps_time: Instant::now(),
        }
    }

    /// Update viewport dimensions.
    pub fn resize(&mut self, width: u32, height: u32) {
        self.width = width.max(1);
        self.height = height.max(1);
    }

    /// Push a rendered frame image and update FPS display in Slint.
    pub fn push_frame(&mut self, image: Option<slint::Image>, vp: &ViewportState) {
        if let Some(image) = image {
            vp.set_viewport_texture(image);
        }

        self.frame_count += 1;
        let elapsed = self.last_fps_time.elapsed();
        if elapsed.as_secs() >= 1 {
            let fps = self.frame_count as f64 / elapsed.as_secs_f64();
            vp.set_fps_text(format!("{:.0} fps", fps).into());
            vp.set_fps_value(fps as f32);
            self.frame_count = 0;
            self.last_fps_time = Instant::now();
        }
    }
}
