//! GPU raytracing orchestration.
//!
//! Coordinates buffer creation, compute shader dispatch, and pixel readback
//! for the multi-pass raytracing pipeline.

pub mod buffers;
pub mod dispatch;
pub mod textures;
pub mod uniforms;

use crate::bvh::Bvh;
use crate::error::{RaytraceError, RaytraceResult};
use crate::pipeline::RaytracePipeline;
use crate::primitive::Primitives;
use crate::settings::RaytraceSettings;

use wgpu::util::DeviceExt;

// ---------------------------------------------------------------------------
// RaytraceParams
// ---------------------------------------------------------------------------

/// Raytracing parameters: output size, camera matrices, and settings.
#[derive(Clone, Debug)]
pub struct RaytraceParams {
    /// Output width in pixels
    pub width: u32,
    /// Output height in pixels
    pub height: u32,
    /// Antialiasing level (1 = no AA, 2 = 2x2 supersampling, etc.)
    pub antialias: u32,
    /// View matrix (world to camera)
    pub view_matrix: [[f32; 4]; 4],
    /// Projection matrix
    pub proj_matrix: [[f32; 4]; 4],
    /// Inverse view matrix
    pub view_inv_matrix: [[f32; 4]; 4],
    /// Inverse projection matrix
    pub proj_inv_matrix: [[f32; 4]; 4],
    /// Camera position in world space
    pub camera_pos: [f32; 4],
    /// Raytracing settings
    pub settings: RaytraceSettings,
}

impl RaytraceParams {
    /// Create new raytracing parameters with default settings.
    pub fn new(width: u32, height: u32) -> Self {
        Self {
            width,
            height,
            antialias: 1,
            view_matrix: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
            proj_matrix: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
            view_inv_matrix: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
            proj_inv_matrix: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
            camera_pos: [0.0, 0.0, 10.0, 1.0],
            settings: RaytraceSettings::default(),
        }
    }

    /// Set antialiasing level.
    pub fn with_antialias(mut self, level: u32) -> Self {
        self.antialias = level.clamp(1, 4);
        self
    }

    /// Set camera matrices.
    pub fn with_camera(
        mut self,
        view: [[f32; 4]; 4],
        proj: [[f32; 4]; 4],
        view_inv: [[f32; 4]; 4],
        proj_inv: [[f32; 4]; 4],
        camera_pos: [f32; 4],
    ) -> Self {
        self.view_matrix = view;
        self.proj_matrix = proj;
        self.view_inv_matrix = view_inv;
        self.proj_inv_matrix = proj_inv;
        self.camera_pos = camera_pos;
        self
    }

    /// Set raytracing settings.
    pub fn with_settings(mut self, settings: RaytraceSettings) -> Self {
        self.settings = settings;
        self
    }
}

// ---------------------------------------------------------------------------
// Main raytrace entry point
// ---------------------------------------------------------------------------

/// Perform raytracing and return the image data as RGBA bytes.
///
/// This is a standalone function that takes device/queue references,
/// avoiding lifetime and ownership issues.
pub fn raytrace(
    device: &wgpu::Device,
    queue: &wgpu::Queue,
    primitives: &Primitives,
    bvh: &Bvh,
    params: &RaytraceParams,
) -> RaytraceResult<Vec<u8>> {
    if primitives.is_empty() {
        return Err(RaytraceError::NoPrimitives);
    }

    log::debug!(
        "Raytrace primitives: {} spheres, {} cylinders, {} triangles",
        primitives.spheres.len(),
        primitives.cylinders.len(),
        primitives.triangles.len()
    );
    if let Some((min, max)) = primitives.aabb() {
        log::debug!(
            "Raytrace scene bounds: ({:.2}, {:.2}, {:.2}) - ({:.2}, {:.2}, {:.2})",
            min[0], min[1], min[2], max[0], max[1], max[2]
        );
    }
    log::debug!(
        "Raytrace camera position: ({:.2}, {:.2}, {:.2})",
        params.camera_pos[0], params.camera_pos[1], params.camera_pos[2]
    );

    // Validate buffer sizes against GPU limits
    buffers::validate_buffer_sizes(device, primitives)?;

    // Pipeline and render dimensions
    let pipeline = RaytracePipeline::new(device)?;
    let supersample = params.antialias.max(1);
    let render_width = params.width * supersample;
    let render_height = params.height * supersample;

    // Create GPU buffers for primitives and BVH
    let sphere_buffer = buffers::create_storage_buffer(device, "Spheres", &primitives.spheres);
    let cylinder_buffer =
        buffers::create_storage_buffer(device, "Cylinders", &primitives.cylinders);
    let triangle_buffer =
        buffers::create_storage_buffer(device, "Triangles", &primitives.triangles);
    let bvh_node_buffer = buffers::create_storage_buffer(device, "BVH Nodes", &bvh.nodes);
    let bvh_index_buffer =
        buffers::create_storage_buffer(device, "BVH Indices", &bvh.primitive_indices);

    // Create uniform buffer
    let gpu_uniforms = uniforms::RaytraceUniforms::from_params(params, primitives, bvh);
    let uniform_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("Raytrace Uniforms"),
        contents: bytemuck::bytes_of(&gpu_uniforms),
        usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
    });

    // Create render textures
    let render_textures = textures::RenderTextures::new(device, render_width, render_height);

    // Bind group for pass 1
    let raytrace_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("Raytrace Bind Group"),
        layout: pipeline.bind_group_layout(),
        entries: &[
            wgpu::BindGroupEntry {
                binding: 0,
                resource: uniform_buffer.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 1,
                resource: sphere_buffer.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 2,
                resource: cylinder_buffer.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 3,
                resource: triangle_buffer.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 4,
                resource: bvh_node_buffer.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 5,
                resource: bvh_index_buffer.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 6,
                resource: wgpu::BindingResource::TextureView(&render_textures.color_view),
            },
            wgpu::BindGroupEntry {
                binding: 7,
                resource: wgpu::BindingResource::TextureView(&render_textures.depth_view),
            },
            wgpu::BindGroupEntry {
                binding: 8,
                resource: wgpu::BindingResource::TextureView(&render_textures.normal_view),
            },
        ],
    });

    // Record all compute passes on a single encoder
    let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
        label: Some("Raytrace Encoder"),
    });

    let workgroups_x = render_width.div_ceil(8);
    let workgroups_y = render_height.div_ceil(8);

    // Pass 1: Raytrace
    dispatch::dispatch_raytrace(
        &mut encoder,
        &pipeline,
        &raytrace_bind_group,
        workgroups_x,
        workgroups_y,
    );

    // Passes 2+3: Edge detection + composite (modes 1, 2, 3)
    let ray_trace_mode = params.settings.ray_trace_mode as u32;
    let composite_texture = if ray_trace_mode > 0 {
        dispatch::dispatch_edge_and_composite(
            device,
            &mut encoder,
            &render_textures,
            params,
            render_width,
            render_height,
            workgroups_x,
            workgroups_y,
        )?
    } else {
        None
    };

    let output_texture = composite_texture
        .as_ref()
        .unwrap_or(&render_textures.color_texture);

    // Record texture-to-buffer copy, then submit and read back
    let (readback_buffer, padded_bytes_per_row) =
        buffers::record_texture_copy(device, &mut encoder, output_texture, render_width, render_height);

    queue.submit(std::iter::once(encoder.finish()));

    let render_pixels = buffers::map_readback_buffer(
        device,
        &readback_buffer,
        padded_bytes_per_row,
        render_width,
        render_height,
    )?;

    // Downsample if supersampling was used
    if supersample > 1 {
        Ok(downsample(
            &render_pixels,
            render_width,
            render_height,
            params.width,
            params.height,
            supersample,
        ))
    } else {
        Ok(render_pixels)
    }
}

// ---------------------------------------------------------------------------
// Downsample
// ---------------------------------------------------------------------------

/// Downsample supersampled image using box filter.
fn downsample(
    src: &[u8],
    src_width: u32,
    _src_height: u32,
    dst_width: u32,
    dst_height: u32,
    factor: u32,
) -> Vec<u8> {
    let mut dst = vec![0u8; (dst_width * dst_height * 4) as usize];
    let factor_sq = (factor * factor) as f32;

    for y in 0..dst_height {
        for x in 0..dst_width {
            let mut r = 0.0f32;
            let mut g = 0.0f32;
            let mut b = 0.0f32;
            let mut a = 0.0f32;

            for sy in 0..factor {
                for sx in 0..factor {
                    let src_x = x * factor + sx;
                    let src_y = y * factor + sy;
                    let src_idx = ((src_y * src_width + src_x) * 4) as usize;
                    r += src[src_idx] as f32;
                    g += src[src_idx + 1] as f32;
                    b += src[src_idx + 2] as f32;
                    a += src[src_idx + 3] as f32;
                }
            }

            let dst_idx = ((y * dst_width + x) * 4) as usize;
            dst[dst_idx] = (r / factor_sq) as u8;
            dst[dst_idx + 1] = (g / factor_sq) as u8;
            dst[dst_idx + 2] = (b / factor_sq) as u8;
            dst[dst_idx + 3] = (a / factor_sq) as u8;
        }
    }

    dst
}
