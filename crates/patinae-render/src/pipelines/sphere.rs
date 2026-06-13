//! Sphere color pipeline — atom-impostor billboard + ray-sphere FS.
//!
//! Bind groups:
//! - 0: `FrameUniforms` (`RenderContext`)
//! - 1: shadow placeholder (`RenderContext`)
//! - 2: [`SceneStoreLayout`] — scene-wide LUTs + `obj` dynamic-uniform
//! - 3: [`SphereParamsLayout`] — per-rep alpha multiplier
//!
//! Opaque spheres write directly to the frame colour/depth targets; spheres
//! with alpha overrides route through WBOIT.

use std::num::NonZeroU64;

use bytemuck::{Pod, Zeroable};

use crate::context::RenderContext;
use crate::pipelines::build_draw_pair;
use crate::representations::sphere::SphereInstance;
use crate::scene_store::SceneStoreLayout;
use crate::shader_source::SPHERE_WGSL;

/// Mirror of `SphereParams` in `shaders/representations/sphere.wgsl`. 16 B.
#[repr(C, align(16))]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct SphereParams {
    pub alpha_mul: f32,
    pub _pad0: f32,
    pub _pad1: f32,
    pub _pad2: f32,
}

impl SphereParams {
    pub const SIZE: u64 = std::mem::size_of::<Self>() as u64;
}

impl Default for SphereParams {
    fn default() -> Self {
        Self {
            alpha_mul: 1.0,
            _pad0: 0.0,
            _pad1: 0.0,
            _pad2: 0.0,
        }
    }
}

/// Group-3 BGL for the sphere render pipeline. Just a uniform block.
pub struct SphereParamsLayout {
    pub bind_group_layout: wgpu::BindGroupLayout,
}

impl SphereParamsLayout {
    pub fn new(device: &wgpu::Device) -> Self {
        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("patinae.sphere.params_layout"),
            entries: &[wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::VERTEX | wgpu::ShaderStages::FRAGMENT,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: NonZeroU64::new(SphereParams::SIZE),
                },
                count: None,
            }],
        });
        Self { bind_group_layout }
    }
}

pub struct SpherePipeline {
    /// WBOIT pipeline — used when the rep contains any per-atom α < 1.
    /// Writes (accum, reveal); depth is read-only against the prepass.
    pub pipeline: wgpu::RenderPipeline,
    /// Opaque fast-path — single host-format target, writes depth.
    /// Selected when every per-atom alpha is fully opaque.
    pub pipeline_opaque: wgpu::RenderPipeline,
}

impl SpherePipeline {
    pub fn new(
        ctx: &RenderContext,
        scene_layout: &SceneStoreLayout,
        sphere_layout: &SphereParamsLayout,
    ) -> Self {
        let pair = build_draw_pair(
            ctx,
            "patinae.sphere",
            SPHERE_WGSL,
            &[
                &ctx.frame.bind_group_layout,
                &ctx.lighting.bind_group_layout,
                &scene_layout.bind_group_layout,
                &sphere_layout.bind_group_layout,
            ],
            SphereInstance::vertex_layout(),
            wgpu::PrimitiveTopology::TriangleList,
        );
        Self {
            pipeline: pair.translucent,
            pipeline_opaque: pair.opaque,
        }
    }
}
