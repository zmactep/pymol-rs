//! Compute shader dispatch for the raytracing pipeline passes.

use wgpu::util::DeviceExt;

use crate::edge_pipeline::{
    CompositeParams, CompositePipeline, DownsampleParams, DownsamplePipeline, EdgeDetectPipeline,
    EdgeParams,
};

use super::textures::{RenderTextures, TextureWithView};
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
#[allow(clippy::too_many_arguments)]
pub(crate) fn dispatch_edge_and_composite(
    device: &wgpu::Device,
    encoder: &mut wgpu::CommandEncoder,
    render_textures: &RenderTextures,
    edge_texture: &TextureWithView,
    composite_texture: &TextureWithView,
    edge_pipeline: &EdgeDetectPipeline,
    composite_pipeline: &CompositePipeline,
    params: &RaytraceParams,
    render_width: u32,
    render_height: u32,
    workgroups_x: u32,
    workgroups_y: u32,
) {
    let ray_trace_mode = params.settings.ray_trace_mode as u32;

    // --- Pass 2: Edge detection ---

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
                resource: wgpu::BindingResource::TextureView(&edge_texture.view),
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
                resource: wgpu::BindingResource::TextureView(&edge_texture.view),
            },
            wgpu::BindGroupEntry {
                binding: 2,
                resource: wgpu::BindingResource::TextureView(&render_textures.depth_view),
            },
            wgpu::BindGroupEntry {
                binding: 3,
                resource: wgpu::BindingResource::TextureView(&composite_texture.view),
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
}

/// Dispatch the standalone GPU antialiasing downsample pass.
#[allow(clippy::too_many_arguments)]
pub(crate) fn dispatch_downsample(
    device: &wgpu::Device,
    encoder: &mut wgpu::CommandEncoder,
    pipeline: &DownsamplePipeline,
    source_view: &wgpu::TextureView,
    output_view: &wgpu::TextureView,
    src_width: u32,
    dst_width: u32,
    dst_height: u32,
    factor: u32,
) {
    let params = DownsampleParams {
        src_width,
        dst_width,
        dst_height,
        factor,
    };
    let params_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("Standalone Downsample Params"),
        contents: bytemuck::bytes_of(&params),
        usage: wgpu::BufferUsages::UNIFORM,
    });
    let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("Standalone Downsample Bind Group"),
        layout: pipeline.bind_group_layout(),
        entries: &[
            wgpu::BindGroupEntry {
                binding: 0,
                resource: wgpu::BindingResource::TextureView(source_view),
            },
            wgpu::BindGroupEntry {
                binding: 1,
                resource: wgpu::BindingResource::TextureView(output_view),
            },
            wgpu::BindGroupEntry {
                binding: 2,
                resource: params_buffer.as_entire_binding(),
            },
        ],
    });

    let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
        label: Some("Standalone Downsample Pass"),
        timestamp_writes: None,
    });
    pass.set_pipeline(pipeline.compute_pipeline());
    pass.set_bind_group(0, &bind_group, &[]);
    pass.dispatch_workgroups(dst_width.div_ceil(8), dst_height.div_ceil(8), 1);
}
