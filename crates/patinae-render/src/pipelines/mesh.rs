//! Mesh color pipeline — atom-driven wireframe of the molecular surface.
//! `StdVertex` LineList topology emitted by the shared `surface_mc.wgsl`
//! kernel with `emit_lines = 1`. See `representations/surface/mod.rs`'s
//! `SurfaceMode::Mesh` branch.
//!
//! Bind groups:
//! - 0: `FrameUniforms` (`RenderContext`)
//! - 1: shadow placeholder (`RenderContext`)
//! - 2: [`SceneStoreLayout`] — scene-wide LUTs + `obj` dynamic-uniform
//! - 3: [`MeshParamsLayout`] — per-rep alpha multiplier
//!
//! Reuses `StdVertex` from the mesh module; only difference vs
//! `SurfacePipeline` is `PrimitiveTopology::LineList`.

use std::num::NonZeroU64;

use bytemuck::{Pod, Zeroable};

use crate::context::RenderContext;
use crate::pipelines::{build_draw_pair, build_fast_overlay_pipeline};
use crate::representations::mesh::StdVertex;
use crate::scene_store::SceneStoreLayout;
use crate::shader_source::MESH_WGSL;

/// Mirror of `MeshParams` in `shaders/representations/mesh.wgsl`. Kept at
/// the same size and alpha offset as `SurfaceParams`.
#[repr(C, align(16))]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct MeshParams {
    pub alpha_mul: f32,
    pub _pad0: f32,
    pub _pad1: f32,
    pub _pad2: f32,
}

impl MeshParams {
    pub const SIZE: u64 = std::mem::size_of::<Self>() as u64;
}

impl Default for MeshParams {
    fn default() -> Self {
        Self {
            alpha_mul: 1.0,
            _pad0: 0.0,
            _pad1: 0.0,
            _pad2: 0.0,
        }
    }
}

pub struct MeshParamsLayout {
    pub bind_group_layout: wgpu::BindGroupLayout,
}

impl MeshParamsLayout {
    pub fn new(device: &wgpu::Device) -> Self {
        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("patinae.mesh.params_layout"),
            entries: &[wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::VERTEX | wgpu::ShaderStages::FRAGMENT,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: NonZeroU64::new(MeshParams::SIZE),
                },
                count: None,
            }],
        });
        Self { bind_group_layout }
    }
}

pub struct MeshPipeline {
    pub pipeline: wgpu::RenderPipeline,
    pub pipeline_opaque: wgpu::RenderPipeline,
    pub pipeline_fast_overlay: wgpu::RenderPipeline,
}

impl MeshPipeline {
    pub fn new(
        ctx: &RenderContext,
        scene_layout: &SceneStoreLayout,
        mesh_layout: &MeshParamsLayout,
    ) -> Self {
        let layouts = [
            &ctx.frame.bind_group_layout,
            &ctx.lighting.bind_group_layout,
            &scene_layout.bind_group_layout,
            &mesh_layout.bind_group_layout,
        ];
        let pair = build_draw_pair(
            ctx,
            "patinae.mesh",
            MESH_WGSL,
            &layouts,
            StdVertex::vertex_layout(),
            wgpu::PrimitiveTopology::LineList,
        );
        let pipeline_fast_overlay = build_fast_overlay_pipeline(
            ctx,
            "patinae.mesh.fast_overlay",
            MESH_WGSL,
            &layouts,
            StdVertex::vertex_layout(),
            wgpu::PrimitiveTopology::LineList,
        );
        Self {
            pipeline: pair.translucent,
            pipeline_opaque: pair.opaque,
            pipeline_fast_overlay,
        }
    }
}
