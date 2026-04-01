//! Compute shader dispatch for the raytracing pipeline passes.

use wgpu::util::DeviceExt;

use crate::edge_pipeline::{CompositePipeline, CompositeParams, EdgeDetectPipeline, EdgeParams};
use crate::error::RaytraceResult;

use super::textures::RenderTextures;
use super::RaytraceParams;

/// Dispatch pass 1: main raytracing compute shader.
pub(crate) fn dispatch_raytrace(
    encoder: &mut wgpu::CommandEncoder,
    pipeline: &crate::pipeline::RaytracePipeline,
    bind_group: &wgpu::BindGroup,
    workgroups_x: u32,
    workgroups_y: u32,
) {
    let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
        label: Some("Raytrace Pass 1"),
        timestamp_writes: None,
    });
    pass.set_pipeline(pipeline.compute_pipeline());
    pass.set_bind_group(0, bind_group, &[]);
    pass.dispatch_workgroups(workgroups_x, workgroups_y, 1);
}

/// Dispatch passes 2+3: edge detection and composite.
///
/// Returns the final composited texture, or `None` if multipass is not needed
/// (this should only be called when `ray_trace_mode > 0`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn dispatch_edge_and_composite(
    device: &wgpu::Device,
    encoder: &mut wgpu::CommandEncoder,
    render_textures: &RenderTextures,
    params: &RaytraceParams,
    render_width: u32,
    render_height: u32,
    workgroups_x: u32,
    workgroups_y: u32,
) -> RaytraceResult<Option<wgpu::Texture>> {
    let ray_trace_mode = params.settings.ray_trace_mode as u32;

    // --- Pass 2: Edge detection ---

    let edge_pipeline = EdgeDetectPipeline::new(device)?;

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

    let edge_params = EdgeParams {
        viewport: [render_width as f32, render_height as f32],
        thickness: params.settings.silhouette_thickness,
        depth_jump: params.settings.silhouette_depth_jump,
    };
    let edge_params_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("Edge Params"),
        contents: bytemuck::bytes_of(&edge_params),
        usage: wgpu::BufferUsages::UNIFORM,
    });

    let edge_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("Edge Detect Bind Group"),
        layout: edge_pipeline.bind_group_layout(),
        entries: &[
            wgpu::BindGroupEntry {
                binding: 0,
                resource: wgpu::BindingResource::TextureView(&render_textures.depth_view),
            },
            wgpu::BindGroupEntry {
                binding: 1,
                resource: wgpu::BindingResource::TextureView(&edge_view),
            },
            wgpu::BindGroupEntry {
                binding: 2,
                resource: edge_params_buffer.as_entire_binding(),
            },
        ],
    });

    {
        let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
            label: Some("Edge Detect Pass 2"),
            timestamp_writes: None,
        });
        pass.set_pipeline(edge_pipeline.compute_pipeline());
        pass.set_bind_group(0, &edge_bind_group, &[]);
        pass.dispatch_workgroups(workgroups_x, workgroups_y, 1);
    }

    // --- Pass 3: Composite ---

    let composite_pipeline = CompositePipeline::new(device)?;

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

    let use_transparent_bg = !params.settings.ray_opaque_background;
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

    let composite_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("Composite Bind Group"),
        layout: composite_pipeline.bind_group_layout(),
        entries: &[
            wgpu::BindGroupEntry {
                binding: 0,
                resource: wgpu::BindingResource::TextureView(&render_textures.color_view),
            },
            wgpu::BindGroupEntry {
                binding: 1,
                resource: wgpu::BindingResource::TextureView(&edge_view),
            },
            wgpu::BindGroupEntry {
                binding: 2,
                resource: wgpu::BindingResource::TextureView(&render_textures.depth_view),
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

    {
        let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
            label: Some("Composite Pass 3"),
            timestamp_writes: None,
        });
        pass.set_pipeline(composite_pipeline.compute_pipeline());
        pass.set_bind_group(0, &composite_bind_group, &[]);
        pass.dispatch_workgroups(workgroups_x, workgroups_y, 1);
    }

    Ok(Some(final_texture))
}
