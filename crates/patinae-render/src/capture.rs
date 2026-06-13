//! Headless PNG capture via the live `RenderState`.
//!
//! Reuses the same frame targets, picking, WBOIT, marking, and silhouette
//! passes as the live viewport, so captured images match interactive frames
//! modulo the requested resolution.
//!
//! The function is synchronous: it submits the encoder, blocks on
//! `device.poll(Wait)` until the readback buffer is mapped, then encodes
//! the PNG via the `image` crate.

use std::path::Path;

use crate::frame::FrameTargets;
use crate::render_input::RenderInput;
use crate::render_state::RenderState;
use crate::uniforms::FrameUniforms;

/// Render `input` headlessly into an offscreen Rgba8Unorm texture sized
/// `(width, height)` and write it as PNG to `path`. Resizes the
/// `RenderState`'s frame targets to match — the caller is responsible for
/// resizing back to the live viewport size before the next live frame.
///
/// `frame_uniforms` should already carry view / proj / fog / clip for the
/// requested capture aspect (the caller builds these via the host bridges).
/// `silhouette_params` mirrors `RenderState::set_silhouette` — pass `None`
/// to skip the pass.
pub fn capture_png(
    state: &mut RenderState,
    path: &Path,
    width: u32,
    height: u32,
    frame_uniforms: &FrameUniforms,
    input: &RenderInput,
) -> Result<(), CaptureError> {
    let width = width.max(1);
    let height = height.max(1);

    // Resize internal targets to the requested aspect. We do not restore
    // the previous size — the caller owns the live viewport size and will
    // resize on its next frame anyway.
    state.targets = FrameTargets::new(&state.ctx.device, width, height, state.ctx.color_format);
    state.uniforms = *frame_uniforms;
    state.uniforms.set_viewport(width, height);
    state.ctx.upload_frame(&state.uniforms);
    state.sync(input);

    // Offscreen colour target — Rgba8Unorm so we can directly read it
    // back as RGBA bytes for PNG encoding.
    let format = wgpu::TextureFormat::Rgba8Unorm;
    let color_texture = state.ctx.device.create_texture(&wgpu::TextureDescriptor {
        label: Some("patinae.capture.color"),
        size: wgpu::Extent3d {
            width,
            height,
            depth_or_array_layers: 1,
        },
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format,
        usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC,
        view_formats: &[],
    });
    let color_view = color_texture.create_view(&wgpu::TextureViewDescriptor::default());

    // Read-back staging buffer. wgpu requires `bytes_per_row` aligned to
    // `COPY_BYTES_PER_ROW_ALIGNMENT` (256). Track the padded stride so we
    // can strip the padding while encoding.
    let unpadded_bpr = width * 4;
    let aligned_bpr = unpadded_bpr.div_ceil(256) * 256;
    let staging = state.ctx.device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("patinae.capture.readback"),
        size: (aligned_bpr * height) as u64,
        usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
        mapped_at_creation: false,
    });

    // Render + copy into staging.
    let mut encoder = state
        .ctx
        .device
        .create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("patinae.capture.encoder"),
        });
    state.render(&color_view, &mut encoder);
    encoder.copy_texture_to_buffer(
        wgpu::TexelCopyTextureInfo {
            texture: &color_texture,
            mip_level: 0,
            origin: wgpu::Origin3d::ZERO,
            aspect: wgpu::TextureAspect::All,
        },
        wgpu::TexelCopyBufferInfo {
            buffer: &staging,
            layout: wgpu::TexelCopyBufferLayout {
                offset: 0,
                bytes_per_row: Some(aligned_bpr),
                rows_per_image: Some(height),
            },
        },
        wgpu::Extent3d {
            width,
            height,
            depth_or_array_layers: 1,
        },
    );
    state.ctx.queue.submit(std::iter::once(encoder.finish()));

    // Block until the staging buffer is mapped.
    let slice = staging.slice(..);
    let (tx, rx) = std::sync::mpsc::channel();
    slice.map_async(wgpu::MapMode::Read, move |r| {
        let _ = tx.send(r);
    });
    state
        .ctx
        .device
        .poll(wgpu::PollType::Wait {
            submission_index: None,
            timeout: None,
        })
        .map_err(|e| CaptureError::Gpu(format!("device poll failed: {e:?}")))?;
    rx.recv()
        .map_err(|_| CaptureError::Gpu("staging buffer map channel closed".into()))?
        .map_err(|e| CaptureError::Gpu(format!("buffer map failed: {e:?}")))?;

    // Strip row padding while copying into a tightly packed RGBA8 buffer.
    let mapped = slice.get_mapped_range();
    let mut tightly_packed = Vec::with_capacity((unpadded_bpr * height) as usize);
    for row in 0..height as usize {
        let start = row * aligned_bpr as usize;
        let end = start + unpadded_bpr as usize;
        tightly_packed.extend_from_slice(&mapped[start..end]);
    }
    drop(mapped);
    staging.unmap();

    let buffer = image::RgbaImage::from_raw(width, height, tightly_packed).ok_or(
        CaptureError::Image("RgbaImage::from_raw size mismatch".into()),
    )?;
    buffer
        .save_with_format(path, image::ImageFormat::Png)
        .map_err(|e| CaptureError::Image(format!("PNG encode failed: {e}")))?;
    Ok(())
}

#[derive(Debug, thiserror::Error)]
pub enum CaptureError {
    #[error("GPU error during capture: {0}")]
    Gpu(String),
    #[error("Image encode error: {0}")]
    Image(String),
}
