//! Multi-directional shadow map ambient occlusion (Skripkin shading mode)
//!
//! Renders the scene from N uniformly distributed directions into a shadow map
//! atlas, then samples the atlas during the main lighting pass to compute
//! world-space ambient occlusion. Since shadows are computed in world space,
//! the atlas only needs re-rendering when geometry changes — not on camera moves.

use bytemuck::{Pod, Zeroable};

use crate::pipeline::depth_stencil_state;
use crate::vertex::{CylinderVertex, MeshVertex, SphereVertex};

/// Maximum number of shadow directions supported
pub const MAX_SHADOW_DIRECTIONS: usize = 128;

/// Shadow map uniform passed to depth-only shaders
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct ShadowUniform {
    /// View-projection matrix for the current shadow direction
    pub view_proj: [[f32; 4]; 4],
}

/// Parameters passed to the main fragment shaders for shadow sampling
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct ShadowParams {
    /// Number of active shadow directions
    pub shadow_count: u32,
    /// Grid size (tiles per row in atlas)
    pub grid_size: u32,
    /// Depth bias (fraction of scene diameter)
    pub bias: f32,
    /// AO intensity (0 = no AO, 1 = full)
    pub intensity: f32,
}

/// Generate N uniformly distributed directions on a unit sphere
/// using Fibonacci spiral (golden angle method)
pub fn fibonacci_sphere_directions(n: usize) -> Vec<[f32; 3]> {
    let golden_ratio = (1.0 + 5.0_f32.sqrt()) / 2.0;
    let angle_increment = std::f32::consts::TAU / golden_ratio;

    (0..n)
        .map(|i| {
            let t = if n > 1 {
                i as f32 / (n as f32 - 1.0)
            } else {
                0.5
            };
            let inclination = (1.0 - 2.0 * t).acos();
            let azimuth = angle_increment * i as f32;

            [
                inclination.sin() * azimuth.cos(),
                inclination.sin() * azimuth.sin(),
                inclination.cos(),
            ]
        })
        .collect()
}

/// Compute an orthographic view-projection matrix looking from `direction`
/// at the scene center, covering the entire bounding sphere.
///
/// Returns a column-major 4x4 matrix as `[[f32; 4]; 4]`.
pub fn compute_shadow_matrix(
    direction: [f32; 3],
    scene_center: [f32; 3],
    scene_radius: f32,
) -> [[f32; 4]; 4] {
    let dx = direction[0];
    let dy = direction[1];
    let dz = direction[2];

    // Eye position: center + direction * radius * 2
    let dist = scene_radius * 2.0;
    let eye = [
        scene_center[0] + dx * dist,
        scene_center[1] + dy * dist,
        scene_center[2] + dz * dist,
    ];

    // Up vector: pick any vector not parallel to direction
    let up = if dy.abs() < 0.99 {
        [0.0, 1.0, 0.0]
    } else {
        [1.0, 0.0, 0.0]
    };

    // Look-at matrix
    let view = look_at(eye, scene_center, up);

    // Orthographic projection covering the bounding sphere
    let r = scene_radius;
    let ortho = ortho_matrix(-r, r, -r, r, 0.0, scene_radius * 4.0);

    // view_proj = ortho * view (column-major multiply)
    mat4_mul(&ortho, &view)
}

/// Look-at view matrix (column-major)
fn look_at(eye: [f32; 3], target: [f32; 3], up: [f32; 3]) -> [[f32; 4]; 4] {
    // Forward = normalize(eye - target) — camera looks toward -Z
    let f = normalize([
        target[0] - eye[0],
        target[1] - eye[1],
        target[2] - eye[2],
    ]);
    let s = normalize(cross(f, up));
    let u = cross(s, f);

    // Column-major
    [
        [s[0], u[0], -f[0], 0.0],
        [s[1], u[1], -f[1], 0.0],
        [s[2], u[2], -f[2], 0.0],
        [
            -(s[0] * eye[0] + s[1] * eye[1] + s[2] * eye[2]),
            -(u[0] * eye[0] + u[1] * eye[1] + u[2] * eye[2]),
            f[0] * eye[0] + f[1] * eye[1] + f[2] * eye[2],
            1.0,
        ],
    ]
}

/// Orthographic projection matrix (column-major, depth range [0, 1])
fn ortho_matrix(
    left: f32,
    right: f32,
    bottom: f32,
    top: f32,
    near: f32,
    far: f32,
) -> [[f32; 4]; 4] {
    [
        [2.0 / (right - left), 0.0, 0.0, 0.0],
        [0.0, 2.0 / (top - bottom), 0.0, 0.0],
        [0.0, 0.0, -1.0 / (far - near), 0.0],
        [
            -(right + left) / (right - left),
            -(top + bottom) / (top - bottom),
            -near / (far - near),
            1.0,
        ],
    ]
}

fn normalize(v: [f32; 3]) -> [f32; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len > 1e-10 {
        [v[0] / len, v[1] / len, v[2] / len]
    } else {
        [0.0, 0.0, 1.0]
    }
}

fn cross(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn mat4_mul(a: &[[f32; 4]; 4], b: &[[f32; 4]; 4]) -> [[f32; 4]; 4] {
    let mut result = [[0.0f32; 4]; 4];
    for col in 0..4 {
        for row in 0..4 {
            result[col][row] = a[0][row] * b[col][0]
                + a[1][row] * b[col][1]
                + a[2][row] * b[col][2]
                + a[3][row] * b[col][3];
        }
    }
    result
}

/// Shadow map atlas texture and associated resources
pub struct MultishadowAtlas {
    pub texture: wgpu::Texture,
    pub depth_view: wgpu::TextureView,
    pub sample_view: wgpu::TextureView,
    pub tile_size: u32,
    pub grid_size: u32,
    pub direction_count: u32,
}

impl MultishadowAtlas {
    /// Create a new shadow map atlas
    pub fn new(device: &wgpu::Device, direction_count: u32, tile_size: u32) -> Self {
        let grid_size = (direction_count as f32).sqrt().ceil() as u32;
        let total_size = grid_size * tile_size;

        let texture = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("Multishadow Atlas"),
            size: wgpu::Extent3d {
                width: total_size,
                height: total_size,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: wgpu::TextureFormat::Depth32Float,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::TEXTURE_BINDING,
            view_formats: &[],
        });

        let depth_view = texture.create_view(&wgpu::TextureViewDescriptor {
            label: Some("Multishadow Atlas Depth View"),
            ..Default::default()
        });

        let sample_view = texture.create_view(&wgpu::TextureViewDescriptor {
            label: Some("Multishadow Atlas Sample View"),
            ..Default::default()
        });

        Self {
            texture,
            depth_view,
            sample_view,
            tile_size,
            grid_size,
            direction_count,
        }
    }

    /// Get viewport (x, y, w, h) for shadow map tile at index i
    pub fn tile_viewport(&self, index: u32) -> (u32, u32, u32, u32) {
        let col = index % self.grid_size;
        let row = index / self.grid_size;
        (
            col * self.tile_size,
            row * self.tile_size,
            self.tile_size,
            self.tile_size,
        )
    }
}

/// Depth-only render pipelines for each geometry type
pub struct ShadowPipelines {
    pub mesh_pipeline: wgpu::RenderPipeline,
    pub sphere_pipeline: wgpu::RenderPipeline,
    pub cylinder_pipeline: wgpu::RenderPipeline,
    pub shadow_bind_group_layout: wgpu::BindGroupLayout,
    pub shadow_uniform_buffer: wgpu::Buffer,
}

impl ShadowPipelines {
    pub fn new(device: &wgpu::Device) -> Self {
        let shadow_bind_group_layout =
            device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("Shadow Uniform Layout"),
                entries: &[wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::VERTEX | wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                }],
            });

        let shadow_uniform_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Shadow Uniform Buffer"),
            size: std::mem::size_of::<ShadowUniform>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("Shadow Pipeline Layout"),
            bind_group_layouts: &[&shadow_bind_group_layout],
            push_constant_ranges: &[],
        });

        // Compile depth-only shaders
        let mesh_shadow_shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("Mesh Shadow Shader"),
            source: wgpu::ShaderSource::Wgsl(
                include_str!("shaders/shadow_mesh.wgsl").into(),
            ),
        });

        let sphere_shadow_shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("Sphere Shadow Shader"),
            source: wgpu::ShaderSource::Wgsl(
                include_str!("shaders/shadow_sphere.wgsl").into(),
            ),
        });

        let cylinder_shadow_shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("Cylinder Shadow Shader"),
            source: wgpu::ShaderSource::Wgsl(
                include_str!("shaders/shadow_cylinder.wgsl").into(),
            ),
        });

        let depth_stencil = Some(depth_stencil_state(true, wgpu::CompareFunction::Less));

        // Mesh depth-only pipeline
        let mesh_pipeline =
            device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
                label: Some("Shadow Mesh Pipeline"),
                layout: Some(&pipeline_layout),
                vertex: wgpu::VertexState {
                    module: &mesh_shadow_shader,
                    entry_point: Some("vs_main"),
                    buffers: &[MeshVertex::layout()],
                    compilation_options: Default::default(),
                },
                fragment: Some(wgpu::FragmentState {
                    module: &mesh_shadow_shader,
                    entry_point: Some("fs_main"),
                    targets: &[],
                    compilation_options: Default::default(),
                }),
                primitive: wgpu::PrimitiveState {
                    topology: wgpu::PrimitiveTopology::TriangleList,
                    cull_mode: None,
                    ..Default::default()
                },
                depth_stencil: depth_stencil.clone(),
                multisample: wgpu::MultisampleState::default(),
                multiview: None,
                cache: None,
            });

        // Billboard vertex layout (shared with sphere and cylinder)
        let billboard_layout = wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<[f32; 2]>() as wgpu::BufferAddress,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: &[wgpu::VertexAttribute {
                format: wgpu::VertexFormat::Float32x2,
                offset: 0,
                shader_location: 10,
            }],
        };

        // Sphere depth-only pipeline
        let sphere_pipeline =
            device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
                label: Some("Shadow Sphere Pipeline"),
                layout: Some(&pipeline_layout),
                vertex: wgpu::VertexState {
                    module: &sphere_shadow_shader,
                    entry_point: Some("vs_main"),
                    buffers: &[billboard_layout.clone(), SphereVertex::layout()],
                    compilation_options: Default::default(),
                },
                fragment: Some(wgpu::FragmentState {
                    module: &sphere_shadow_shader,
                    entry_point: Some("fs_main"),
                    targets: &[],
                    compilation_options: Default::default(),
                }),
                primitive: wgpu::PrimitiveState {
                    topology: wgpu::PrimitiveTopology::TriangleList,
                    cull_mode: None,
                    ..Default::default()
                },
                depth_stencil: depth_stencil.clone(),
                multisample: wgpu::MultisampleState::default(),
                multiview: None,
                cache: None,
            });

        // Cylinder depth-only pipeline
        let cylinder_pipeline =
            device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
                label: Some("Shadow Cylinder Pipeline"),
                layout: Some(&pipeline_layout),
                vertex: wgpu::VertexState {
                    module: &cylinder_shadow_shader,
                    entry_point: Some("vs_main"),
                    buffers: &[billboard_layout, CylinderVertex::layout()],
                    compilation_options: Default::default(),
                },
                fragment: Some(wgpu::FragmentState {
                    module: &cylinder_shadow_shader,
                    entry_point: Some("fs_main"),
                    targets: &[],
                    compilation_options: Default::default(),
                }),
                primitive: wgpu::PrimitiveState {
                    topology: wgpu::PrimitiveTopology::TriangleList,
                    cull_mode: None,
                    ..Default::default()
                },
                depth_stencil,
                multisample: wgpu::MultisampleState::default(),
                multiview: None,
                cache: None,
            });

        Self {
            mesh_pipeline,
            sphere_pipeline,
            cylinder_pipeline,
            shadow_bind_group_layout,
            shadow_uniform_buffer,
        }
    }

    /// Create a bind group for the shadow uniform buffer
    pub fn create_bind_group(&self, device: &wgpu::Device) -> wgpu::BindGroup {
        device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("Shadow Uniform Bind Group"),
            layout: &self.shadow_bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: self.shadow_uniform_buffer.as_entire_binding(),
            }],
        })
    }

    /// Update shadow uniform with a new view-projection matrix
    pub fn update_uniform(&self, queue: &wgpu::Queue, view_proj: &[[f32; 4]; 4]) {
        let uniform = ShadowUniform {
            view_proj: *view_proj,
        };
        queue.write_buffer(
            &self.shadow_uniform_buffer,
            0,
            bytemuck::bytes_of(&uniform),
        );
    }
}

/// Bind group layout and resources for sampling the shadow atlas in main shaders
pub struct ShadowSampling {
    pub bind_group_layout: wgpu::BindGroupLayout,
    pub params_buffer: wgpu::Buffer,
    pub matrices_buffer: wgpu::Buffer,
    pub sampler: wgpu::Sampler,
}

impl ShadowSampling {
    pub fn new(device: &wgpu::Device) -> Self {
        let bind_group_layout =
            device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("Shadow Sampling Layout"),
                entries: &[
                    // binding 0: shadow atlas depth texture
                    wgpu::BindGroupLayoutEntry {
                        binding: 0,
                        visibility: wgpu::ShaderStages::FRAGMENT,
                        ty: wgpu::BindingType::Texture {
                            sample_type: wgpu::TextureSampleType::Float { filterable: false },
                            view_dimension: wgpu::TextureViewDimension::D2,
                            multisampled: false,
                        },
                        count: None,
                    },
                    // binding 1: sampler
                    wgpu::BindGroupLayoutEntry {
                        binding: 1,
                        visibility: wgpu::ShaderStages::FRAGMENT,
                        ty: wgpu::BindingType::Sampler(wgpu::SamplerBindingType::NonFiltering),
                        count: None,
                    },
                    // binding 2: shadow params
                    wgpu::BindGroupLayoutEntry {
                        binding: 2,
                        visibility: wgpu::ShaderStages::FRAGMENT,
                        ty: wgpu::BindingType::Buffer {
                            ty: wgpu::BufferBindingType::Uniform,
                            has_dynamic_offset: false,
                            min_binding_size: None,
                        },
                        count: None,
                    },
                    // binding 3: shadow matrices (storage buffer)
                    wgpu::BindGroupLayoutEntry {
                        binding: 3,
                        visibility: wgpu::ShaderStages::FRAGMENT,
                        ty: wgpu::BindingType::Buffer {
                            ty: wgpu::BufferBindingType::Storage { read_only: true },
                            has_dynamic_offset: false,
                            min_binding_size: None,
                        },
                        count: None,
                    },
                ],
            });

        let params_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Shadow Params Buffer"),
            size: std::mem::size_of::<ShadowParams>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        // Storage buffer for shadow matrices (MAX_SHADOW_DIRECTIONS * mat4x4)
        let matrices_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Shadow Matrices Buffer"),
            size: (MAX_SHADOW_DIRECTIONS * 16 * 4) as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let sampler = device.create_sampler(&wgpu::SamplerDescriptor {
            label: Some("Shadow Atlas Sampler"),
            mag_filter: wgpu::FilterMode::Nearest,
            min_filter: wgpu::FilterMode::Nearest,
            ..Default::default()
        });

        Self {
            bind_group_layout,
            params_buffer,
            matrices_buffer,
            sampler,
        }
    }

    /// Create a bind group for the given shadow atlas
    pub fn create_bind_group(
        &self,
        device: &wgpu::Device,
        atlas_view: &wgpu::TextureView,
    ) -> wgpu::BindGroup {
        device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("Shadow Sampling Bind Group"),
            layout: &self.bind_group_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: wgpu::BindingResource::TextureView(atlas_view),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: wgpu::BindingResource::Sampler(&self.sampler),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: self.params_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: self.matrices_buffer.as_entire_binding(),
                },
            ],
        })
    }

    /// Upload shadow parameters and matrices
    pub fn update(
        &self,
        queue: &wgpu::Queue,
        params: &ShadowParams,
        matrices: &[[[f32; 4]; 4]],
    ) {
        queue.write_buffer(&self.params_buffer, 0, bytemuck::bytes_of(params));
        queue.write_buffer(&self.matrices_buffer, 0, bytemuck::cast_slice(matrices));
    }
}

/// Dummy (disabled) bind group for when shading mode is not Skripkin.
/// All zeros / identity — the shader checks shadow_count == 0 and skips.
pub fn create_disabled_shadow_bind_group(
    device: &wgpu::Device,
    sampling: &ShadowSampling,
) -> wgpu::BindGroup {
    // Create a tiny 1x1 depth texture as placeholder
    let placeholder = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("Shadow Placeholder Texture"),
        size: wgpu::Extent3d {
            width: 1,
            height: 1,
            depth_or_array_layers: 1,
        },
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format: wgpu::TextureFormat::Depth32Float,
        usage: wgpu::TextureUsages::TEXTURE_BINDING | wgpu::TextureUsages::RENDER_ATTACHMENT,
        view_formats: &[],
    });
    let placeholder_view = placeholder.create_view(&wgpu::TextureViewDescriptor::default());

    sampling.create_bind_group(device, &placeholder_view)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fibonacci_sphere_directions() {
        let dirs = fibonacci_sphere_directions(64);
        assert_eq!(dirs.len(), 64);
        // All directions should be unit vectors
        for d in &dirs {
            let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            assert!((len - 1.0).abs() < 0.01, "Direction not unit length: {}", len);
        }
    }

    #[test]
    fn test_shadow_matrix() {
        let mat = compute_shadow_matrix([0.0, 0.0, 1.0], [0.0, 0.0, 0.0], 10.0);
        // Should be a valid matrix (not all zeros)
        let sum: f32 = mat.iter().flat_map(|col| col.iter()).map(|v| v.abs()).sum();
        assert!(sum > 0.0);
    }

    #[test]
    fn test_shadow_params_size() {
        assert_eq!(std::mem::size_of::<ShadowParams>() % 4, 0);
    }

    #[test]
    fn test_shadow_uniform_size() {
        assert_eq!(std::mem::size_of::<ShadowUniform>(), 64);
    }
}
