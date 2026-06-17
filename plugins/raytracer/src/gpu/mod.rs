//! GPU raytracing orchestration.
//!
//! Coordinates buffer creation, compute shader dispatch, and pixel readback
//! for the multi-pass raytracing pipeline.

pub mod buffers;
pub mod dispatch;
pub mod textures;
pub mod uniforms;

use std::time::Instant;

use crate::bvh::Bvh;
use crate::edge_pipeline::{CompositePipeline, DownsamplePipeline, EdgeDetectPipeline};
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
// RaytraceProfile
// ---------------------------------------------------------------------------

/// Timing breakdown for one raytrace call.
#[derive(Clone, Debug, Default)]
pub struct RaytraceProfile {
    /// Total wall time measured by the caller-side raytrace function.
    pub total_ms: f64,
    /// Buffer-size validation time.
    pub validate_ms: f64,
    /// Main raytrace pipeline creation time.
    pub pipeline_ms: f64,
    /// Primitive and BVH storage-buffer upload time.
    pub upload_buffers_ms: f64,
    /// Uniform-buffer creation time.
    pub uniform_ms: f64,
    /// Output texture allocation time.
    pub textures_ms: f64,
    /// Main bind-group creation time.
    pub bind_group_ms: f64,
    /// Command encoder creation and compute/copy recording time.
    pub encode_ms: f64,
    /// Queue submission call time, excluding GPU completion.
    pub submit_ms: f64,
    /// Readback map and device-poll wait time.
    pub readback_ms: f64,
    /// CPU downsample time when antialiasing is greater than one.
    pub downsample_ms: f64,
    /// Output width requested by the caller.
    pub output_width: u32,
    /// Output height requested by the caller.
    pub output_height: u32,
    /// Supersampled render width used internally.
    pub render_width: u32,
    /// Supersampled render height used internally.
    pub render_height: u32,
    /// Antialiasing factor used by the render.
    pub antialias: u32,
    /// Ray trace mode used by the render.
    pub ray_trace_mode: i32,
    /// Number of sphere primitives.
    pub sphere_count: usize,
    /// Number of cylinder primitives.
    pub cylinder_count: usize,
    /// Number of capsule primitives.
    pub capsule_count: usize,
    /// Number of triangle primitives.
    pub triangle_count: usize,
    /// Number of flattened BVH nodes.
    pub bvh_node_count: usize,
    /// Final readback byte count.
    pub readback_bytes: u64,
}

impl RaytraceProfile {
    /// Return the total primitive count.
    pub fn total_primitives(&self) -> usize {
        self.sphere_count + self.cylinder_count + self.capsule_count + self.triangle_count
    }

    /// Return final output pixel count.
    pub fn output_pixels(&self) -> u64 {
        u64::from(self.output_width) * u64::from(self.output_height)
    }

    /// Return supersampled render pixel count.
    pub fn render_pixels(&self) -> u64 {
        u64::from(self.render_width) * u64::from(self.render_height)
    }
}

/// Image bytes plus timing details from a profiled raytrace.
#[derive(Clone, Debug)]
pub struct RaytraceProfiledImage {
    /// Final RGBA image bytes.
    pub image: Vec<u8>,
    /// Timing breakdown for the render.
    pub profile: RaytraceProfile,
}

// ---------------------------------------------------------------------------
// Main raytrace entry point
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct FrameKey {
    width: u32,
    height: u32,
    antialias: u32,
    ray_trace_mode: i32,
}

impl FrameKey {
    fn from_params(params: &RaytraceParams) -> Self {
        Self {
            width: params.width,
            height: params.height,
            antialias: params.antialias.max(1),
            ray_trace_mode: params.settings.ray_trace_mode,
        }
    }

    fn render_width(self) -> u32 {
        self.width * self.antialias
    }

    fn render_height(self) -> u32 {
        self.height * self.antialias
    }

    fn needs_edge(self) -> bool {
        (self.ray_trace_mode as u32) > 0
    }

    fn needs_downsample(self) -> bool {
        self.antialias > 1
    }
}

struct FrameResources {
    key: FrameKey,
    render_textures: textures::RenderTextures,
    edge_texture: Option<textures::TextureWithView>,
    composite_texture: Option<textures::TextureWithView>,
    downsample_texture: Option<textures::TextureWithView>,
    readback_buffer: buffers::ReadbackBuffer,
}

impl FrameResources {
    fn new(device: &wgpu::Device, key: FrameKey) -> Self {
        let render_width = key.render_width();
        let render_height = key.render_height();
        Self {
            key,
            render_textures: textures::RenderTextures::new(device, render_width, render_height),
            edge_texture: key
                .needs_edge()
                .then(|| textures::edge_texture(device, render_width, render_height)),
            composite_texture: key
                .needs_edge()
                .then(|| textures::composite_texture(device, render_width, render_height)),
            downsample_texture: key
                .needs_downsample()
                .then(|| textures::downsample_texture(device, key.width, key.height)),
            readback_buffer: buffers::ReadbackBuffer::new(device, key.width, key.height),
        }
    }
}

/// Reusable standalone raytrace context.
///
/// A context owns GPU pipelines and same-shape frame resources for repeated
/// standalone renders. Use it only with the same [`wgpu::Device`] that created
/// it; `wgpu` resources are device-owned and cannot be moved between devices.
#[must_use]
pub struct RaytraceContext {
    pipeline: RaytracePipeline,
    edge_pipeline: Option<EdgeDetectPipeline>,
    composite_pipeline: Option<CompositePipeline>,
    downsample_pipeline: Option<DownsamplePipeline>,
    frame_resources: Option<FrameResources>,
}

impl std::fmt::Debug for RaytraceContext {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("RaytraceContext")
            .field("has_edge_pipeline", &self.edge_pipeline.is_some())
            .field("has_composite_pipeline", &self.composite_pipeline.is_some())
            .field(
                "has_downsample_pipeline",
                &self.downsample_pipeline.is_some(),
            )
            .field(
                "frame_key",
                &self.frame_resources.as_ref().map(|resources| resources.key),
            )
            .finish()
    }
}

impl RaytraceContext {
    /// Create a reusable standalone raytrace context.
    ///
    /// # Errors
    /// Returns an error if the main raytrace compute pipeline cannot be
    /// created for the provided device.
    pub fn new(device: &wgpu::Device) -> RaytraceResult<Self> {
        Ok(Self {
            pipeline: RaytracePipeline::new(device)?,
            edge_pipeline: None,
            composite_pipeline: None,
            downsample_pipeline: None,
            frame_resources: None,
        })
    }

    /// Perform raytracing and return RGBA bytes.
    ///
    /// Reuses cached pipelines and frame resources while the output shape,
    /// antialiasing factor, and ray trace mode remain unchanged.
    ///
    /// # Errors
    /// Returns an error for empty scenes, oversized GPU buffers, pipeline
    /// creation failures, or readback mapping failures.
    pub fn raytrace(
        &mut self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        primitives: &Primitives,
        bvh: &Bvh,
        params: &RaytraceParams,
    ) -> RaytraceResult<Vec<u8>> {
        self.raytrace_inner(device, queue, primitives, bvh, params, None)
    }

    /// Perform raytracing and return profiling details.
    ///
    /// Warm calls report `pipeline_ms` only for optional pipeline cache misses;
    /// the main pipeline is paid once in [`RaytraceContext::new`].
    ///
    /// # Errors
    /// Returns an error for empty scenes, oversized GPU buffers, pipeline
    /// creation failures, or readback mapping failures.
    pub fn raytrace_profiled(
        &mut self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        primitives: &Primitives,
        bvh: &Bvh,
        params: &RaytraceParams,
    ) -> RaytraceResult<RaytraceProfiledImage> {
        let mut profile = RaytraceProfile::default();
        let image =
            self.raytrace_inner(device, queue, primitives, bvh, params, Some(&mut profile))?;
        Ok(RaytraceProfiledImage { image, profile })
    }

    fn raytrace_inner(
        &mut self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        primitives: &Primitives,
        bvh: &Bvh,
        params: &RaytraceParams,
        mut profile: Option<&mut RaytraceProfile>,
    ) -> RaytraceResult<Vec<u8>> {
        if primitives.is_empty() {
            return Err(RaytraceError::NoPrimitives);
        }
        let total_start = Instant::now();

        log::debug!(
            "Raytrace primitives: {} spheres, {} cylinders, {} capsules, {} triangles",
            primitives.spheres.len(),
            primitives.cylinders.len(),
            primitives.capsules.len(),
            primitives.triangles.len()
        );
        if let Some((min, max)) = primitives.aabb() {
            log::debug!(
                "Raytrace scene bounds: ({:.2}, {:.2}, {:.2}) - ({:.2}, {:.2}, {:.2})",
                min[0],
                min[1],
                min[2],
                max[0],
                max[1],
                max[2]
            );
        }
        log::debug!(
            "Raytrace camera position: ({:.2}, {:.2}, {:.2})",
            params.camera_pos[0],
            params.camera_pos[1],
            params.camera_pos[2]
        );

        let stage_start = Instant::now();
        buffers::validate_buffer_sizes(device, primitives)?;
        set_profile(&mut profile, |profile| {
            profile.validate_ms = elapsed_ms(stage_start);
        });

        let frame_key = FrameKey::from_params(params);
        let render_width = frame_key.render_width();
        let render_height = frame_key.render_height();
        set_profile(&mut profile, |profile| {
            profile.output_width = params.width;
            profile.output_height = params.height;
            profile.render_width = render_width;
            profile.render_height = render_height;
            profile.antialias = frame_key.antialias;
            profile.ray_trace_mode = params.settings.ray_trace_mode;
            profile.sphere_count = primitives.spheres.len();
            profile.cylinder_count = primitives.cylinders.len();
            profile.capsule_count = primitives.capsules.len();
            profile.triangle_count = primitives.triangles.len();
            profile.bvh_node_count = bvh.nodes.len();
            profile.readback_bytes = u64::from(params.width) * u64::from(params.height) * 4;
        });

        let pipeline_ms = self.ensure_optional_pipelines(
            device,
            frame_key.needs_edge(),
            frame_key.needs_downsample(),
        )?;
        set_profile(&mut profile, |profile| {
            profile.pipeline_ms += pipeline_ms;
        });

        let stage_start = Instant::now();
        let sphere_buffer = buffers::create_storage_buffer(device, "Spheres", &primitives.spheres);
        let cylinder_buffer =
            buffers::create_storage_buffer(device, "Cylinders", &primitives.cylinders);
        let capsule_buffer =
            buffers::create_storage_buffer(device, "Capsules", &primitives.capsules);
        let triangle_buffer =
            buffers::create_storage_buffer(device, "Triangles", &primitives.triangles);
        let bvh_node_buffer = buffers::create_storage_buffer(device, "BVH Nodes", &bvh.nodes);
        let bvh_index_buffer =
            buffers::create_storage_buffer(device, "BVH Indices", &bvh.primitive_indices);
        set_profile(&mut profile, |profile| {
            profile.upload_buffers_ms = elapsed_ms(stage_start);
        });

        let stage_start = Instant::now();
        let gpu_uniforms = uniforms::RaytraceUniforms::from_params(params, primitives, bvh);
        let uniform_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Raytrace Uniforms"),
            contents: bytemuck::bytes_of(&gpu_uniforms),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });
        set_profile(&mut profile, |profile| {
            profile.uniform_ms = elapsed_ms(stage_start);
        });

        let stage_start = Instant::now();
        let reused_frame_resources = self.ensure_frame_resources(device, frame_key);
        set_profile(&mut profile, |profile| {
            profile.textures_ms = if reused_frame_resources {
                0.0
            } else {
                elapsed_ms(stage_start)
            };
        });
        let frame_resources = self
            .frame_resources
            .as_ref()
            .expect("frame resources were initialized");

        let stage_start = Instant::now();
        let raytrace_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("Raytrace Bind Group"),
            layout: self.pipeline.bind_group_layout(),
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
                    resource: capsule_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 4,
                    resource: triangle_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 5,
                    resource: bvh_node_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 6,
                    resource: bvh_index_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 7,
                    resource: wgpu::BindingResource::TextureView(
                        &frame_resources.render_textures.color_view,
                    ),
                },
                wgpu::BindGroupEntry {
                    binding: 8,
                    resource: wgpu::BindingResource::TextureView(
                        &frame_resources.render_textures.depth_view,
                    ),
                },
                wgpu::BindGroupEntry {
                    binding: 9,
                    resource: wgpu::BindingResource::TextureView(
                        &frame_resources.render_textures.normal_view,
                    ),
                },
            ],
        });
        set_profile(&mut profile, |profile| {
            profile.bind_group_ms = elapsed_ms(stage_start);
        });

        let stage_start = Instant::now();
        let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("Raytrace Encoder"),
        });

        let workgroups_x = render_width.div_ceil(8);
        let workgroups_y = render_height.div_ceil(8);

        dispatch::dispatch_raytrace(
            &mut encoder,
            &self.pipeline,
            &raytrace_bind_group,
            workgroups_x,
            workgroups_y,
        );

        if frame_key.needs_edge() {
            dispatch::dispatch_edge_and_composite(
                device,
                &mut encoder,
                &frame_resources.render_textures,
                frame_resources
                    .edge_texture
                    .as_ref()
                    .expect("edge texture exists for edge mode"),
                frame_resources
                    .composite_texture
                    .as_ref()
                    .expect("composite texture exists for edge mode"),
                self.edge_pipeline
                    .as_ref()
                    .expect("edge pipeline exists for edge mode"),
                self.composite_pipeline
                    .as_ref()
                    .expect("composite pipeline exists for edge mode"),
                params,
                render_width,
                render_height,
                workgroups_x,
                workgroups_y,
            );
        }

        let final_texture = if frame_key.needs_downsample() {
            let source_view = if frame_key.needs_edge() {
                &frame_resources
                    .composite_texture
                    .as_ref()
                    .expect("composite texture exists for downsample")
                    .view
            } else {
                &frame_resources.render_textures.color_view
            };
            let downsample_texture = frame_resources
                .downsample_texture
                .as_ref()
                .expect("downsample texture exists for antialiasing");
            dispatch::dispatch_downsample(
                device,
                &mut encoder,
                self.downsample_pipeline
                    .as_ref()
                    .expect("downsample pipeline exists for antialiasing"),
                source_view,
                &downsample_texture.view,
                render_width,
                params.width,
                params.height,
                frame_key.antialias,
            );
            &downsample_texture.texture
        } else if frame_key.needs_edge() {
            &frame_resources
                .composite_texture
                .as_ref()
                .expect("composite texture exists for edge mode")
                .texture
        } else {
            &frame_resources.render_textures.color_texture
        };

        buffers::record_texture_copy_to_buffer(
            &mut encoder,
            final_texture,
            &frame_resources.readback_buffer,
            params.width,
            params.height,
        );
        set_profile(&mut profile, |profile| {
            profile.encode_ms = elapsed_ms(stage_start);
        });

        let stage_start = Instant::now();
        queue.submit(std::iter::once(encoder.finish()));
        set_profile(&mut profile, |profile| {
            profile.submit_ms = elapsed_ms(stage_start);
        });

        let stage_start = Instant::now();
        let image = buffers::map_readback_buffer(
            device,
            &frame_resources.readback_buffer.buffer,
            frame_resources.readback_buffer.layout.padded_bytes_per_row,
            params.width,
            params.height,
        )?;
        set_profile(&mut profile, |profile| {
            profile.readback_ms = elapsed_ms(stage_start);
            profile.total_ms = elapsed_ms(total_start);
        });

        Ok(image)
    }

    fn ensure_optional_pipelines(
        &mut self,
        device: &wgpu::Device,
        needs_edge: bool,
        needs_downsample: bool,
    ) -> RaytraceResult<f64> {
        let mut pipeline_ms = 0.0;
        if needs_edge && self.edge_pipeline.is_none() {
            let stage_start = Instant::now();
            self.edge_pipeline = Some(EdgeDetectPipeline::new(device)?);
            pipeline_ms += elapsed_ms(stage_start);
        }
        if needs_edge && self.composite_pipeline.is_none() {
            let stage_start = Instant::now();
            self.composite_pipeline = Some(CompositePipeline::new(device)?);
            pipeline_ms += elapsed_ms(stage_start);
        }
        if needs_downsample && self.downsample_pipeline.is_none() {
            let stage_start = Instant::now();
            self.downsample_pipeline = Some(DownsamplePipeline::new(device)?);
            pipeline_ms += elapsed_ms(stage_start);
        }
        Ok(pipeline_ms)
    }

    fn ensure_frame_resources(&mut self, device: &wgpu::Device, key: FrameKey) -> bool {
        if self
            .frame_resources
            .as_ref()
            .is_some_and(|resources| resources.key == key)
        {
            return true;
        }
        self.frame_resources = Some(FrameResources::new(device, key));
        false
    }
}

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
    let mut context = RaytraceContext::new(device)?;
    context.raytrace(device, queue, primitives, bvh, params)
}

/// Perform raytracing and return image data with timing details.
pub fn raytrace_profiled(
    device: &wgpu::Device,
    queue: &wgpu::Queue,
    primitives: &Primitives,
    bvh: &Bvh,
    params: &RaytraceParams,
) -> RaytraceResult<RaytraceProfiledImage> {
    let total_start = Instant::now();
    let stage_start = Instant::now();
    let mut context = RaytraceContext::new(device)?;
    let setup_ms = elapsed_ms(stage_start);
    let mut result = context.raytrace_profiled(device, queue, primitives, bvh, params)?;
    result.profile.pipeline_ms += setup_ms;
    result.profile.total_ms = elapsed_ms(total_start);
    Ok(result)
}

fn elapsed_ms(start: Instant) -> f64 {
    start.elapsed().as_secs_f64() * 1000.0
}

fn set_profile(
    profile: &mut Option<&mut RaytraceProfile>,
    update: impl FnOnce(&mut RaytraceProfile),
) {
    if let Some(profile) = profile.as_deref_mut() {
        update(profile);
    }
}
