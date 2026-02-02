//! Raytracing context and GPU resource management

use wgpu::util::DeviceExt;

use crate::bvh::Bvh;
use crate::edge_pipeline::{CompositePipeline, CompositeParams, EdgeDetectPipeline, EdgeParams};
use crate::error::{RaytraceError, RaytraceResult};
use crate::pipeline::RaytracePipeline;
use crate::primitive::Primitives;
use crate::render::RaytraceParams;

/// Get the maximum storage buffer size from device limits
fn get_max_storage_buffer_size(device: &wgpu::Device) -> usize {
    device.limits().max_storage_buffer_binding_size as usize
}

/// Perform raytracing and return the image data as RGBA bytes
///
/// This is a standalone function that takes device/queue references,
/// avoiding lifetime and ownership issues.
///
/// # Arguments
/// * `device` - wgpu device reference
/// * `queue` - wgpu queue reference
/// * `primitives` - Scene primitives (spheres, cylinders, triangles)
/// * `bvh` - Bounding volume hierarchy for acceleration
/// * `params` - Rendering parameters (size, camera, settings)
///
/// # Returns
/// * `Ok(Vec<u8>)` - RGBA image data, row-major, top-to-bottom
/// * `Err(RaytraceError)` - If rendering fails
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

    // Log primitive counts for debugging
    log::debug!(
        "Raytrace primitives: {} spheres, {} cylinders, {} triangles",
        primitives.spheres.len(),
        primitives.cylinders.len(),
        primitives.triangles.len()
    );

    // Log scene bounds
    if let Some((min, max)) = primitives.aabb() {
        log::debug!(
            "Raytrace scene bounds: ({:.2}, {:.2}, {:.2}) - ({:.2}, {:.2}, {:.2})",
            min[0], min[1], min[2], max[0], max[1], max[2]
        );
    }

    // Log camera position
    log::debug!(
        "Raytrace camera position: ({:.2}, {:.2}, {:.2})",
        params.camera_pos[0], params.camera_pos[1], params.camera_pos[2]
    );

    // Check buffer sizes don't exceed GPU limits
    let max_buffer_size = get_max_storage_buffer_size(device);
    let sphere_size = primitives.spheres.len() * std::mem::size_of::<crate::primitive::GpuSphere>();
    let cylinder_size = primitives.cylinders.len() * std::mem::size_of::<crate::primitive::GpuCylinder>();
    let triangle_size = primitives.triangles.len() * std::mem::size_of::<crate::primitive::GpuTriangle>();

    log::debug!(
        "Buffer sizes: spheres={:.1}MB, cylinders={:.1}MB, triangles={:.1}MB (limit={:.1}MB)",
        sphere_size as f64 / (1024.0 * 1024.0),
        cylinder_size as f64 / (1024.0 * 1024.0),
        triangle_size as f64 / (1024.0 * 1024.0),
        max_buffer_size as f64 / (1024.0 * 1024.0)
    );
    
    if triangle_size > max_buffer_size {
        return Err(RaytraceError::RenderFailed(format!(
            "Triangle buffer too large ({:.1} MB, limit {:.1} MB). Scene has {} triangles. \
             Try reducing geometry detail with 'set cartoon_sampling' or use a simpler representation.",
            triangle_size as f64 / (1024.0 * 1024.0),
            max_buffer_size as f64 / (1024.0 * 1024.0),
            primitives.triangles.len()
        )));
    }
    if sphere_size > max_buffer_size {
        return Err(RaytraceError::RenderFailed(format!(
            "Sphere buffer too large ({:.1} MB, limit {:.1} MB). Scene has {} spheres.",
            sphere_size as f64 / (1024.0 * 1024.0),
            max_buffer_size as f64 / (1024.0 * 1024.0),
            primitives.spheres.len()
        )));
    }
    if cylinder_size > max_buffer_size {
        return Err(RaytraceError::RenderFailed(format!(
            "Cylinder buffer too large ({:.1} MB, limit {:.1} MB). Scene has {} cylinders.",
            cylinder_size as f64 / (1024.0 * 1024.0),
            max_buffer_size as f64 / (1024.0 * 1024.0),
            primitives.cylinders.len()
        )));
    }

    // Create pipeline
    let pipeline = RaytracePipeline::new(device)?;

    // Calculate actual render size (with supersampling)
    let supersample = params.antialias.max(1);
    let render_width = params.width * supersample;
    let render_height = params.height * supersample;

    // Create GPU buffers for primitives
    let sphere_buffer = create_storage_buffer(device, "Spheres", &primitives.spheres);
    let cylinder_buffer = create_storage_buffer(device, "Cylinders", &primitives.cylinders);
    let triangle_buffer = create_storage_buffer(device, "Triangles", &primitives.triangles);
    let bvh_node_buffer = create_storage_buffer(device, "BVH Nodes", &bvh.nodes);
    let bvh_index_buffer = create_storage_buffer(device, "BVH Indices", &bvh.primitive_indices);

    // Create uniform buffer with raytracing parameters
    let uniforms = RaytraceUniforms::from_params(params, primitives, bvh);
    let uniform_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("Raytrace Uniforms"),
        contents: bytemuck::bytes_of(&uniforms),
        usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
    });

    // Get ray_trace_mode from settings
    let ray_trace_mode = params.settings.ray_trace_mode as u32;
    let use_multipass = ray_trace_mode > 0;

    // Create output texture (color)
    let color_texture = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("Raytrace Color Output"),
        size: wgpu::Extent3d {
            width: render_width,
            height: render_height,
            depth_or_array_layers: 1,
        },
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format: wgpu::TextureFormat::Rgba8Unorm,
        usage: wgpu::TextureUsages::STORAGE_BINDING
            | wgpu::TextureUsages::TEXTURE_BINDING
            | wgpu::TextureUsages::COPY_SRC,
        view_formats: &[],
    });
    let color_view = color_texture.create_view(&wgpu::TextureViewDescriptor::default());

    // Create depth texture (for edge detection)
    let depth_texture = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("Raytrace Depth Output"),
        size: wgpu::Extent3d {
            width: render_width,
            height: render_height,
            depth_or_array_layers: 1,
        },
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format: wgpu::TextureFormat::R32Float,
        usage: wgpu::TextureUsages::STORAGE_BINDING | wgpu::TextureUsages::TEXTURE_BINDING,
        view_formats: &[],
    });
    let depth_view = depth_texture.create_view(&wgpu::TextureViewDescriptor::default());

    // Create normal texture (for edge detection)
    let normal_texture = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("Raytrace Normal Output"),
        size: wgpu::Extent3d {
            width: render_width,
            height: render_height,
            depth_or_array_layers: 1,
        },
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format: wgpu::TextureFormat::Rgba16Float,
        usage: wgpu::TextureUsages::STORAGE_BINDING | wgpu::TextureUsages::TEXTURE_BINDING,
        view_formats: &[],
    });
    let normal_view = normal_texture.create_view(&wgpu::TextureViewDescriptor::default());

    // Create bind group for pass 1 (raytrace)
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
                resource: wgpu::BindingResource::TextureView(&color_view),
            },
            wgpu::BindGroupEntry {
                binding: 7,
                resource: wgpu::BindingResource::TextureView(&depth_view),
            },
            wgpu::BindGroupEntry {
                binding: 8,
                resource: wgpu::BindingResource::TextureView(&normal_view),
            },
        ],
    });

    // Create encoder
    let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
        label: Some("Raytrace Encoder"),
    });

    // Workgroup size is 8x8
    let workgroups_x = (render_width + 7) / 8;
    let workgroups_y = (render_height + 7) / 8;

    // Pass 1: Raytrace
    {
        let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
            label: Some("Raytrace Pass 1"),
            timestamp_writes: None,
        });
        pass.set_pipeline(pipeline.compute_pipeline());
        pass.set_bind_group(0, &raytrace_bind_group, &[]);
        pass.dispatch_workgroups(workgroups_x, workgroups_y, 1);
    }

    // For modes 1, 2, 3: run edge detection and composite passes
    // The final output texture depends on whether we're doing multipass
    let output_texture = if use_multipass {
        // Create edge detection pipeline
        let edge_pipeline = EdgeDetectPipeline::new(device)?;

        // Create edge texture (use R32Float for storage binding support)
        let edge_texture = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("Edge Detection Output"),
            size: wgpu::Extent3d {
                width: render_width,
                height: render_height,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: wgpu::TextureFormat::R32Float,
            usage: wgpu::TextureUsages::STORAGE_BINDING | wgpu::TextureUsages::TEXTURE_BINDING,
            view_formats: &[],
        });
        let edge_view = edge_texture.create_view(&wgpu::TextureViewDescriptor::default());

        // Create edge params uniform using PyMOL's gradient-of-gradient algorithm
        // Scale PyMOL parameters for our normalized depth buffer (0-1 range)
        // PyMOL uses screen-space depth, so their thresholds need scaling
        //
        // Balanced thresholds for both surface (smooth) and cartoon (sharp) geometry
        // slope_factor: PyMOL 0.6 -> 0.004 for our depth range
        // depth_factor: PyMOL 0.1 -> ~0.000005 for our depth range
        // gain: Controls gradient amplification
        let edge_params = EdgeParams {
            viewport: [render_width as f32, render_height as f32],
            slope_factor: params.settings.ray_trace_slope_factor / 150.0, // 0.6 -> 0.004 (balanced)
            depth_factor: (params.settings.ray_trace_depth_factor / 45.0).powi(2), // 0.1 -> ~0.000005
            disco_factor: params.settings.ray_trace_disco_factor * 5.0, // 0.05 -> 0.25
            gain: 1.0 / (params.settings.ray_trace_gain + 0.001) * 12.0, // 0.12 -> ~100
            _pad: [0.0; 2],
        };
        let edge_params_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Edge Params"),
            contents: bytemuck::bytes_of(&edge_params),
            usage: wgpu::BufferUsages::UNIFORM,
        });

        // Create edge detection bind group
        let edge_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("Edge Detect Bind Group"),
            layout: edge_pipeline.bind_group_layout(),
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: wgpu::BindingResource::TextureView(&depth_view),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: wgpu::BindingResource::TextureView(&normal_view),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: wgpu::BindingResource::TextureView(&edge_view),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: edge_params_buffer.as_entire_binding(),
                },
            ],
        });

        // Pass 2: Edge detection
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("Edge Detect Pass 2"),
                timestamp_writes: None,
            });
            pass.set_pipeline(edge_pipeline.compute_pipeline());
            pass.set_bind_group(0, &edge_bind_group, &[]);
            pass.dispatch_workgroups(workgroups_x, workgroups_y, 1);
        }

        // Create composite pipeline
        let composite_pipeline = CompositePipeline::new(device)?;

        // Create final output texture
        let final_texture = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("Composite Output"),
            size: wgpu::Extent3d {
                width: render_width,
                height: render_height,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: wgpu::TextureFormat::Rgba8Unorm,
            usage: wgpu::TextureUsages::STORAGE_BINDING | wgpu::TextureUsages::COPY_SRC,
            view_formats: &[],
        });
        let final_view = final_texture.create_view(&wgpu::TextureViewDescriptor::default());

        // Create composite params uniform
        let use_transparent_bg = params.settings.ray_opaque_background == 0;
        let composite_params = CompositeParams {
            viewport: [render_width as f32, render_height as f32],
            mode: ray_trace_mode,
            use_transparent_bg: if use_transparent_bg { 1 } else { 0 },
            edge_color: [
                params.settings.ray_trace_color[0],
                params.settings.ray_trace_color[1],
                params.settings.ray_trace_color[2],
            ],
            quantize_levels: 4.0,
            bg_color: [
                params.settings.bg_color[0],
                params.settings.bg_color[1],
                params.settings.bg_color[2],
            ],
            _pad: 0.0,
        };
        let composite_params_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Composite Params"),
            contents: bytemuck::bytes_of(&composite_params),
            usage: wgpu::BufferUsages::UNIFORM,
        });

        // Create composite bind group
        let composite_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("Composite Bind Group"),
            layout: composite_pipeline.bind_group_layout(),
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: wgpu::BindingResource::TextureView(&color_view),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: wgpu::BindingResource::TextureView(&edge_view),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: wgpu::BindingResource::TextureView(&depth_view),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: wgpu::BindingResource::TextureView(&final_view),
                },
                wgpu::BindGroupEntry {
                    binding: 4,
                    resource: composite_params_buffer.as_entire_binding(),
                },
            ],
        });

        // Pass 3: Composite
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("Composite Pass 3"),
                timestamp_writes: None,
            });
            pass.set_pipeline(composite_pipeline.compute_pipeline());
            pass.set_bind_group(0, &composite_bind_group, &[]);
            pass.dispatch_workgroups(workgroups_x, workgroups_y, 1);
        }

        final_texture
    } else {
        // Mode 0: use color texture directly as output
        color_texture
    };

    // Copy texture to buffer for readback
    let bytes_per_pixel = 4u32;
    let unpadded_bytes_per_row = render_width * bytes_per_pixel;
    let align = wgpu::COPY_BYTES_PER_ROW_ALIGNMENT;
    let padded_bytes_per_row = (unpadded_bytes_per_row + align - 1) / align * align;
    let buffer_size = (padded_bytes_per_row * render_height) as u64;

    let readback_buffer = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("Raytrace Readback"),
        size: buffer_size,
        usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
        mapped_at_creation: false,
    });

    encoder.copy_texture_to_buffer(
        wgpu::TexelCopyTextureInfo {
            texture: &output_texture,
            mip_level: 0,
            origin: wgpu::Origin3d::ZERO,
            aspect: wgpu::TextureAspect::All,
        },
        wgpu::TexelCopyBufferInfo {
            buffer: &readback_buffer,
            layout: wgpu::TexelCopyBufferLayout {
                offset: 0,
                bytes_per_row: Some(padded_bytes_per_row),
                rows_per_image: Some(render_height),
            },
        },
        wgpu::Extent3d {
            width: render_width,
            height: render_height,
            depth_or_array_layers: 1,
        },
    );

    queue.submit(std::iter::once(encoder.finish()));

    // Read back results
    let buffer_slice = readback_buffer.slice(..);
    let (tx, rx) = std::sync::mpsc::channel();
    buffer_slice.map_async(wgpu::MapMode::Read, move |result| {
        // Ignore send error - receiver may have been dropped
        let _ = tx.send(result);
    });
    let _ = device.poll(wgpu::PollType::Wait {
        submission_index: None,
        timeout: None,
    });
    rx.recv()
        .map_err(|e| RaytraceError::RenderFailed(format!("Failed to receive map result: {}", e)))?
        .map_err(|e| RaytraceError::RenderFailed(format!("Failed to map buffer: {:?}", e)))?;

    // Extract pixel data (removing row padding)
    let data = buffer_slice.get_mapped_range();
    let mut render_pixels =
        Vec::with_capacity((render_width * render_height * bytes_per_pixel) as usize);
    for row in 0..render_height {
        let start = (row * padded_bytes_per_row) as usize;
        let end = start + (render_width * bytes_per_pixel) as usize;
        render_pixels.extend_from_slice(&data[start..end]);
    }
    drop(data);
    readback_buffer.unmap();

    // Downsample if supersampling was used
    if supersample > 1 {
        let final_pixels = downsample(
            &render_pixels,
            render_width,
            render_height,
            params.width,
            params.height,
            supersample,
        );
        Ok(final_pixels)
    } else {
        Ok(render_pixels)
    }
}

/// Create a storage buffer from data
///
/// When data is empty, creates a buffer with one zeroed element to satisfy
/// WGSL's minimum binding size requirements.
fn create_storage_buffer<T: bytemuck::Pod + bytemuck::Zeroable>(
    device: &wgpu::Device,
    label: &str,
    data: &[T],
) -> wgpu::Buffer {
    if data.is_empty() {
        // Create a buffer with one zeroed element - WGSL requires the buffer
        // to be at least the size of one struct element for array types
        let dummy = T::zeroed();
        device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some(label),
            contents: bytemuck::bytes_of(&dummy),
            usage: wgpu::BufferUsages::STORAGE,
        })
    } else {
        device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some(label),
            contents: bytemuck::cast_slice(data),
            usage: wgpu::BufferUsages::STORAGE,
        })
    }
}

/// Downsample supersampled image using box filter
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

/// GPU uniforms for raytracing
#[repr(C)]
#[derive(Copy, Clone, Debug, bytemuck::Pod, bytemuck::Zeroable)]
pub struct RaytraceUniforms {
    // Camera
    pub view_matrix: [[f32; 4]; 4],
    pub proj_matrix: [[f32; 4]; 4],
    pub view_inv_matrix: [[f32; 4]; 4],
    pub proj_inv_matrix: [[f32; 4]; 4],
    pub camera_pos: [f32; 4],

    // Viewport
    pub viewport: [f32; 4], // width, height, 1/width, 1/height

    // Light
    pub light_dir: [f32; 4],
    pub ambient: f32,
    pub direct: f32,
    pub reflect: f32,
    pub specular: f32,
    pub shininess: f32,
    pub _pad_light1: f32,
    pub _pad_light2: f32,
    pub _pad_light3: f32,

    // Background
    pub bg_color: [f32; 4],

    // Fog
    pub fog_start: f32,
    pub fog_end: f32,
    pub fog_density: f32,
    pub _pad0: f32,
    pub fog_color: [f32; 4],

    // Primitive counts
    pub sphere_count: u32,
    pub cylinder_count: u32,
    pub triangle_count: u32,
    pub bvh_node_count: u32,

    // Ray settings
    pub ray_shadow: u32,
    pub ray_max_passes: u32,
    pub ray_trace_fog: u32,
    pub ray_transparency_shadows: u32,

    // Ray trace mode settings
    pub ray_trace_mode: u32,
    pub ray_opaque_background: i32,
    pub _pad1: u32,
    pub _pad2: u32,
    pub ray_trace_color: [f32; 4],
}

impl RaytraceUniforms {
    /// Create uniforms from raytracing parameters
    pub fn from_params(params: &RaytraceParams, primitives: &Primitives, bvh: &Bvh) -> Self {
        let settings = &params.settings;
        let supersample = params.antialias.max(1);
        let width = (params.width * supersample) as f32;
        let height = (params.height * supersample) as f32;

        Self {
            view_matrix: params.view_matrix,
            proj_matrix: params.proj_matrix,
            view_inv_matrix: params.view_inv_matrix,
            proj_inv_matrix: params.proj_inv_matrix,
            camera_pos: params.camera_pos,
            viewport: [width, height, 1.0 / width, 1.0 / height],
            light_dir: settings.light_dir,
            ambient: settings.ambient,
            direct: settings.direct,
            reflect: settings.reflect,
            specular: settings.specular,
            shininess: settings.shininess,
            _pad_light1: 0.0,
            _pad_light2: 0.0,
            _pad_light3: 0.0,
            bg_color: settings.bg_color,
            fog_start: settings.fog_start,
            fog_end: settings.fog_end,
            fog_density: settings.fog_density,
            _pad0: 0.0,
            fog_color: settings.fog_color,
            sphere_count: primitives.spheres.len() as u32,
            cylinder_count: primitives.cylinders.len() as u32,
            triangle_count: primitives.triangles.len() as u32,
            bvh_node_count: bvh.nodes.len() as u32,
            ray_shadow: if settings.ray_shadow { 1 } else { 0 },
            ray_max_passes: settings.ray_max_passes,
            ray_trace_fog: if settings.ray_trace_fog { 1 } else { 0 },
            ray_transparency_shadows: if settings.ray_transparency_shadows { 1 } else { 0 },
            ray_trace_mode: settings.ray_trace_mode as u32,
            ray_opaque_background: settings.ray_opaque_background,
            _pad1: 0,
            _pad2: 0,
            ray_trace_color: settings.ray_trace_color,
        }
    }
}
