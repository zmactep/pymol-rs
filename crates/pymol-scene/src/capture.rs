//! Screenshot capture functionality
//!
//! This module provides PNG capture capabilities for rendering the scene
//! to an image file.

use std::path::Path;

use pymol_color::{ChainColors, ElementColors, NamedColors};
use pymol_render::{ColorResolver, GlobalUniforms, RenderContext};
use pymol_settings::GlobalSettings;

use crate::camera::Camera;
use crate::error::ViewerError;
use crate::object::{Object, ObjectRegistry};
use crate::uniform::compute_reflect_scale;

/// Capture the current scene to a PNG file
///
/// This is the single source of truth for PNG capture logic, used by both
/// the `Viewer` and GUI's `ViewerAdapter`.
///
/// # Arguments
///
/// * `path` - Output file path
/// * `width` - Optional width in pixels
/// * `height` - Optional height in pixels  
/// * `context` - Render context with GPU device/queue
/// * `camera` - Camera for view/projection matrices (aspect ratio temporarily modified)
/// * `registry` - Object registry containing molecules to render
/// * `settings` - Global settings for lighting, fog, etc.
/// * `named_colors` - Named color definitions
/// * `element_colors` - Element color definitions
/// * `chain_colors` - Chain color definitions
/// * `clear_color` - Background color RGB
/// * `default_size` - Fallback size when width/height not specified
///
/// # Example
///
/// ```ignore
/// capture_png_to_file(
///     Path::new("screenshot.png"),
///     Some(1920), Some(1080),
///     &context, &mut camera, &mut registry,
///     &settings, &named_colors, &element_colors, &chain_colors,
///     [0.0, 0.0, 0.0], (1024, 768),
/// )?;
/// ```
pub fn capture_png_to_file(
    path: &Path,
    width: Option<u32>,
    height: Option<u32>,
    context: &RenderContext,
    camera: &mut Camera,
    registry: &mut ObjectRegistry,
    settings: &GlobalSettings,
    named_colors: &NamedColors,
    element_colors: &ElementColors,
    chain_colors: &ChainColors,
    clear_color: [f32; 3],
    default_size: (u32, u32),
) -> Result<(), ViewerError> {
    // Determine output size
    let (output_width, output_height) = match (width, height) {
        (Some(w), Some(h)) => (w, h),
        (Some(w), None) => {
            // Maintain aspect ratio based on width
            let aspect = camera.aspect();
            (w, (w as f32 / aspect) as u32)
        }
        (None, Some(h)) => {
            // Maintain aspect ratio based on height
            let aspect = camera.aspect();
            ((h as f32 * aspect) as u32, h)
        }
        (None, None) => default_size,
    };

    let device = context.device();
    let queue = context.queue();

    // Create offscreen color texture using the same format as render pipelines
    let texture_format = context.surface_format();
    let color_texture = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("Screenshot Color Texture"),
        size: wgpu::Extent3d {
            width: output_width,
            height: output_height,
            depth_or_array_layers: 1,
        },
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format: texture_format,
        usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC,
        view_formats: &[],
    });
    let color_view = color_texture.create_view(&wgpu::TextureViewDescriptor::default());

    // Create offscreen depth texture
    let depth_texture = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("Screenshot Depth Texture"),
        size: wgpu::Extent3d {
            width: output_width,
            height: output_height,
            depth_or_array_layers: 1,
        },
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format: wgpu::TextureFormat::Depth32Float,
        usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::TEXTURE_BINDING,
        view_formats: &[],
    });
    let depth_view = depth_texture.create_view(&wgpu::TextureViewDescriptor::default());

    // Update uniforms for the capture resolution
    let mut uniforms = GlobalUniforms::new();

    // Temporarily update camera aspect ratio for correct projection
    let original_aspect = camera.aspect();
    camera.set_aspect(output_width as f32 / output_height as f32);

    uniforms.set_camera(camera.view_matrix(), camera.projection_matrix());
    uniforms.set_background(clear_color);
    uniforms.set_viewport(output_width as f32, output_height as f32);

    let camera_pos = camera.world_position();
    uniforms.camera_pos = [camera_pos.x, camera_pos.y, camera_pos.z, 1.0];

    // Multi-light support (PyMOL light_count setting)
    let light_count = settings.get_int(pymol_settings::id::light_count);

    // spec_count controls how many positional lights contribute specular
    // -1 (default) means all positional lights contribute specular
    let spec_count = settings.get_int(pymol_settings::id::spec_count);

    // Gather light directions from settings (light, light2, ..., light9)
    let light_setting_ids = [
        pymol_settings::id::light,
        pymol_settings::id::light2,
        pymol_settings::id::light3,
        pymol_settings::id::light4,
        pymol_settings::id::light5,
        pymol_settings::id::light6,
        pymol_settings::id::light7,
        pymol_settings::id::light8,
        pymol_settings::id::light9,
    ];

    let light_dirs: Vec<[f32; 3]> = light_setting_ids
        .iter()
        .map(|&id| settings.get_float3(id))
        .collect();

    uniforms.set_lights(light_count, spec_count, &light_dirs);

    // Lighting settings
    let ambient = settings.get_float(pymol_settings::id::ambient);
    let direct = settings.get_float(pymol_settings::id::direct);
    let reflect = settings.get_float(pymol_settings::id::reflect);
    let specular = settings.get_float(pymol_settings::id::specular);
    let shininess = settings.get_float(pymol_settings::id::shininess);
    let spec_direct = settings.get_float(pymol_settings::id::spec_direct);
    let spec_direct_power = settings.get_float(pymol_settings::id::spec_direct_power);

    // PyMOL brightness consistency adjustments:
    // 1. Scale reflect based on light directions to maintain consistent brightness
    // 2. When light_count < 2, redirect reflect energy to direct (headlight)
    let reflect_scale = compute_reflect_scale(light_count, &light_dirs);
    let reflect_scaled = reflect * reflect_scale;

    let (direct_adjusted, reflect_final) = if light_count < 2 {
        // No positional lights - add reflect to direct (PyMOL behavior)
        ((direct + reflect_scaled).min(1.0), 0.0)
    } else {
        (direct, reflect_scaled)
    };

    uniforms.set_lighting(
        ambient,
        direct_adjusted,
        reflect_final,
        specular,
        shininess,
        spec_direct,
        spec_direct_power,
    );

    // Clip planes
    let current_view = camera.current_view();
    uniforms.set_clip_planes(current_view.clip_front, current_view.clip_back);

    // Fog parameters
    let depth_cue_enabled = settings.get_bool(pymol_settings::id::depth_cue);
    let fog_density = settings.get_float(pymol_settings::id::fog);
    if depth_cue_enabled && fog_density > 0.0 {
        let fog_start_ratio = settings.get_float(pymol_settings::id::fog_start);
        let fog_start_actual = (current_view.clip_back - current_view.clip_front) * fog_start_ratio + current_view.clip_front;
        let fog_end_actual = if (fog_density - 1.0).abs() < 0.001 {
            current_view.clip_back
        } else {
            fog_start_actual + (current_view.clip_back - fog_start_actual) / fog_density
        };
        uniforms.set_fog(fog_start_actual, fog_end_actual, fog_density, clear_color);
        uniforms.set_depth_cue(1.0);
    }

    context.update_uniforms(&uniforms);

    // Restore original aspect ratio
    camera.set_aspect(original_aspect);

    // Prepare molecules for rendering
    let names: Vec<_> = registry.names().map(|s| s.to_string()).collect();
    for name in &names {
        let color_resolver = ColorResolver::new(named_colors, element_colors, chain_colors);
        if let Some(mol_obj) = registry.get_molecule_mut(name) {
            mol_obj.prepare_render(context, &color_resolver, settings);
        }
    }

    // Create command encoder
    let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
        label: Some("Screenshot Encoder"),
    });

    // Render pass
    {
        let mut render_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
            label: Some("Screenshot Render Pass"),
            color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                view: &color_view,
                resolve_target: None,
                ops: wgpu::Operations {
                    load: wgpu::LoadOp::Clear(wgpu::Color {
                        r: clear_color[0] as f64,
                        g: clear_color[1] as f64,
                        b: clear_color[2] as f64,
                        a: 1.0,
                    }),
                    store: wgpu::StoreOp::Store,
                },
                depth_slice: None,
            })],
            depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                view: &depth_view,
                depth_ops: Some(wgpu::Operations {
                    load: wgpu::LoadOp::Clear(1.0),
                    store: wgpu::StoreOp::Store,
                }),
                stencil_ops: None,
            }),
            timestamp_writes: None,
            occlusion_query_set: None,
        });

        // Render all enabled objects
        for name in &names {
            if let Some(mol_obj) = registry.get_molecule(name) {
                if mol_obj.is_enabled() {
                    mol_obj.render(&mut render_pass, context);
                }
            }
        }
    }

    // Calculate buffer size with proper alignment
    // wgpu requires rows to be aligned to 256 bytes
    let bytes_per_pixel = 4u32; // RGBA8
    let unpadded_bytes_per_row = output_width * bytes_per_pixel;
    let align = wgpu::COPY_BYTES_PER_ROW_ALIGNMENT;
    let padded_bytes_per_row = (unpadded_bytes_per_row + align - 1) / align * align;
    let buffer_size = (padded_bytes_per_row * output_height) as u64;

    // Create output buffer
    let output_buffer = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("Screenshot Buffer"),
        size: buffer_size,
        usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
        mapped_at_creation: false,
    });

    // Copy texture to buffer
    encoder.copy_texture_to_buffer(
        wgpu::TexelCopyTextureInfo {
            texture: &color_texture,
            mip_level: 0,
            origin: wgpu::Origin3d::ZERO,
            aspect: wgpu::TextureAspect::All,
        },
        wgpu::TexelCopyBufferInfo {
            buffer: &output_buffer,
            layout: wgpu::TexelCopyBufferLayout {
                offset: 0,
                bytes_per_row: Some(padded_bytes_per_row),
                rows_per_image: Some(output_height),
            },
        },
        wgpu::Extent3d {
            width: output_width,
            height: output_height,
            depth_or_array_layers: 1,
        },
    );

    // Submit and wait
    queue.submit(std::iter::once(encoder.finish()));

    // Map buffer and read data
    let buffer_slice = output_buffer.slice(..);
    let (tx, rx) = std::sync::mpsc::channel();
    buffer_slice.map_async(wgpu::MapMode::Read, move |result| {
        tx.send(result).unwrap();
    });
    device.poll(wgpu::PollType::Wait { submission_index: None, timeout: None }).ok();
    rx.recv()
        .map_err(|e| ViewerError::GpuInitFailed(format!("Failed to receive map result: {}", e)))?
        .map_err(|e| ViewerError::GpuInitFailed(format!("Failed to map buffer: {:?}", e)))?;

    // Read pixel data (removing padding)
    let data = buffer_slice.get_mapped_range();
    let mut pixels = Vec::with_capacity((output_width * output_height * bytes_per_pixel) as usize);
    for row in 0..output_height {
        let start = (row * padded_bytes_per_row) as usize;
        let end = start + (output_width * bytes_per_pixel) as usize;
        pixels.extend_from_slice(&data[start..end]);
    }
    drop(data);
    output_buffer.unmap();

    // Convert BGRA to RGBA if needed (common on many platforms)
    let is_bgra = matches!(
        texture_format,
        wgpu::TextureFormat::Bgra8Unorm | wgpu::TextureFormat::Bgra8UnormSrgb
    );
    if is_bgra {
        // Swap B and R channels: BGRA -> RGBA
        for chunk in pixels.chunks_exact_mut(4) {
            chunk.swap(0, 2);
        }
    }

    // Save as PNG using image crate
    let img: image::RgbaImage = image::ImageBuffer::from_raw(output_width, output_height, pixels)
        .ok_or_else(|| ViewerError::GpuInitFailed("Failed to create image buffer".to_string()))?;

    img.save(path)
        .map_err(|e| ViewerError::IoError(format!("Failed to save PNG: {}", e)))?;

    Ok(())
}
