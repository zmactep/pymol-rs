//! Render context and GPU resource management
//!
//! The `RenderContext` manages the wgpu device, queue, and shared resources
//! needed for rendering molecular representations.

use crate::buffer::{create_uniform_buffer, update_uniform_buffer};
use crate::pipeline::{depth_stencil_state, get_blend_state, BlendMode, PipelineKey, PipelineType};
use crate::uniforms::GlobalUniforms;
use crate::vertex::{CylinderVertex, DotVertex, LineVertex, MeshVertex, SphereVertex, BILLBOARD_VERTICES, QUAD_INDICES};
use wgpu::util::DeviceExt;

/// Main render context for molecular visualization
///
/// This struct holds all shared GPU resources and provides methods for
/// creating pipelines and managing the render state.
pub struct RenderContext {
    /// The wgpu device for creating GPU resources
    device: wgpu::Device,
    /// The wgpu queue for submitting commands
    queue: wgpu::Queue,
    /// Surface texture format for color targets
    surface_format: wgpu::TextureFormat,
    /// Depth texture format
    depth_format: wgpu::TextureFormat,
    /// Global uniforms buffer
    uniform_buffer: wgpu::Buffer,
    /// Bind group for global uniforms
    uniform_bind_group: wgpu::BindGroup,
    /// Bind group layout for global uniforms
    uniform_bind_group_layout: wgpu::BindGroupLayout,
    /// Billboard vertex buffer (shared by sphere and cylinder impostors)
    billboard_vertex_buffer: wgpu::Buffer,
    /// Quad index buffer (shared by impostors)
    quad_index_buffer: wgpu::Buffer,
    /// Compiled shaders
    shaders: Shaders,
}

/// Compiled shader modules
struct Shaders {
    sphere: wgpu::ShaderModule,
    cylinder: wgpu::ShaderModule,
    line: wgpu::ShaderModule,
    mesh: wgpu::ShaderModule,
    dot: wgpu::ShaderModule,
}

impl RenderContext {
    /// Create a new render context
    ///
    /// Takes ownership of the device and queue. The surface_format should match
    /// the target surface or texture format.
    pub fn new(
        device: wgpu::Device,
        queue: wgpu::Queue,
        surface_format: wgpu::TextureFormat,
    ) -> Self {
        let depth_format = wgpu::TextureFormat::Depth32Float;

        // Create uniform buffer and bind group layout
        let uniforms = GlobalUniforms::default();
        let uniform_buffer = create_uniform_buffer(&device, "Global Uniforms", &uniforms);

        let uniform_bind_group_layout =
            device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("Global Uniforms Layout"),
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

        let uniform_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("Global Uniforms Bind Group"),
            layout: &uniform_bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: uniform_buffer.as_entire_binding(),
            }],
        });

        // Create shared billboard vertices
        let billboard_vertex_buffer =
            device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("Billboard Vertices"),
                contents: bytemuck::cast_slice(&BILLBOARD_VERTICES),
                usage: wgpu::BufferUsages::VERTEX,
            });

        // Create shared quad indices
        let quad_index_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Quad Indices"),
            contents: bytemuck::cast_slice(&QUAD_INDICES),
            usage: wgpu::BufferUsages::INDEX,
        });

        // Compile shaders
        let shaders = Shaders {
            sphere: device.create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some("Sphere Shader"),
                source: wgpu::ShaderSource::Wgsl(include_str!("shaders/sphere.wgsl").into()),
            }),
            cylinder: device.create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some("Cylinder Shader"),
                source: wgpu::ShaderSource::Wgsl(include_str!("shaders/cylinder.wgsl").into()),
            }),
            line: device.create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some("Line Shader"),
                source: wgpu::ShaderSource::Wgsl(include_str!("shaders/line.wgsl").into()),
            }),
            mesh: device.create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some("Mesh Shader"),
                source: wgpu::ShaderSource::Wgsl(include_str!("shaders/mesh.wgsl").into()),
            }),
            dot: device.create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some("Dot Shader"),
                source: wgpu::ShaderSource::Wgsl(include_str!("shaders/dot.wgsl").into()),
            }),
        };

        Self {
            device,
            queue,
            surface_format,
            depth_format,
            uniform_buffer,
            uniform_bind_group,
            uniform_bind_group_layout,
            billboard_vertex_buffer,
            quad_index_buffer,
            shaders,
        }
    }

    /// Get a reference to the wgpu device
    pub fn device(&self) -> &wgpu::Device {
        &self.device
    }

    /// Get a reference to the wgpu queue
    pub fn queue(&self) -> &wgpu::Queue {
        &self.queue
    }

    /// Get the surface format
    pub fn surface_format(&self) -> wgpu::TextureFormat {
        self.surface_format
    }

    /// Get the depth format
    pub fn depth_format(&self) -> wgpu::TextureFormat {
        self.depth_format
    }

    /// Update global uniforms
    pub fn update_uniforms(&self, uniforms: &GlobalUniforms) {
        update_uniform_buffer(&self.queue, &self.uniform_buffer, uniforms);
    }

    /// Get the uniform bind group for shaders
    pub fn uniform_bind_group(&self) -> &wgpu::BindGroup {
        &self.uniform_bind_group
    }

    /// Get the billboard vertex buffer
    pub fn billboard_vertex_buffer(&self) -> &wgpu::Buffer {
        &self.billboard_vertex_buffer
    }

    /// Get the quad index buffer
    pub fn quad_index_buffer(&self) -> &wgpu::Buffer {
        &self.quad_index_buffer
    }

    /// Get or create the sphere render pipeline
    pub fn sphere_pipeline(&self, blend_mode: BlendMode) -> wgpu::RenderPipeline {
        let key = PipelineKey {
            pipeline_type: PipelineType::Sphere,
            blend_mode,
            depth_write: blend_mode == BlendMode::Opaque,
        };

        self.create_pipeline(key, &self.shaders.sphere, SphereVertex::layout())
    }

    /// Get or create the cylinder render pipeline
    pub fn cylinder_pipeline(&self, blend_mode: BlendMode) -> wgpu::RenderPipeline {
        let key = PipelineKey {
            pipeline_type: PipelineType::Cylinder,
            blend_mode,
            depth_write: blend_mode == BlendMode::Opaque,
        };

        self.create_pipeline(key, &self.shaders.cylinder, CylinderVertex::layout())
    }

    /// Get or create the line render pipeline
    pub fn line_pipeline(&self, blend_mode: BlendMode) -> wgpu::RenderPipeline {
        let key = PipelineKey {
            pipeline_type: PipelineType::Line,
            blend_mode,
            depth_write: blend_mode == BlendMode::Opaque,
        };

        self.create_line_pipeline(key)
    }

    /// Get or create the mesh render pipeline
    pub fn mesh_pipeline(&self, blend_mode: BlendMode) -> wgpu::RenderPipeline {
        let key = PipelineKey {
            pipeline_type: PipelineType::Mesh,
            blend_mode,
            depth_write: blend_mode == BlendMode::Opaque,
        };

        self.create_mesh_pipeline(key)
    }

    /// Get or create the dot render pipeline
    pub fn dot_pipeline(&self, blend_mode: BlendMode) -> wgpu::RenderPipeline {
        let key = PipelineKey {
            pipeline_type: PipelineType::Dot,
            blend_mode,
            depth_write: blend_mode == BlendMode::Opaque,
        };

        self.create_dot_pipeline(key)
    }

    /// Create a render pipeline for impostor rendering (sphere, cylinder)
    fn create_pipeline(
        &self,
        key: PipelineKey,
        shader: &wgpu::ShaderModule,
        instance_layout: wgpu::VertexBufferLayout<'static>,
    ) -> wgpu::RenderPipeline {
        let pipeline_layout = self.device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("Impostor Pipeline Layout"),
            bind_group_layouts: &[&self.uniform_bind_group_layout],
            push_constant_ranges: &[],
        });

        // Billboard vertex layout
        let billboard_layout = wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<[f32; 2]>() as wgpu::BufferAddress,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: &[wgpu::VertexAttribute {
                format: wgpu::VertexFormat::Float32x2,
                offset: 0,
                shader_location: 10, // High location to avoid conflicts
            }],
        };

        self.device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("Impostor Pipeline"),
            layout: Some(&pipeline_layout),
            vertex: wgpu::VertexState {
                module: shader,
                entry_point: Some("vs_main"),
                buffers: &[billboard_layout, instance_layout],
                compilation_options: Default::default(),
            },
            fragment: Some(wgpu::FragmentState {
                module: shader,
                entry_point: Some("fs_main"),
                targets: &[Some(wgpu::ColorTargetState {
                    format: self.surface_format,
                    blend: get_blend_state(key.blend_mode),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
                compilation_options: Default::default(),
            }),
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::TriangleList,
                strip_index_format: None,
                front_face: wgpu::FrontFace::Ccw,
                cull_mode: None, // Impostors need both sides
                unclipped_depth: false,
                polygon_mode: wgpu::PolygonMode::Fill,
                conservative: false,
            },
            depth_stencil: Some(depth_stencil_state(key.depth_write)),
            multisample: wgpu::MultisampleState {
                count: 1,
                mask: !0,
                alpha_to_coverage_enabled: false,
            },
            multiview: None,
            cache: None,
        })
    }

    /// Create the line render pipeline
    fn create_line_pipeline(&self, key: PipelineKey) -> wgpu::RenderPipeline {
        let pipeline_layout = self.device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("Line Pipeline Layout"),
            bind_group_layouts: &[&self.uniform_bind_group_layout],
            push_constant_ranges: &[],
        });

        self.device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("Line Pipeline"),
            layout: Some(&pipeline_layout),
            vertex: wgpu::VertexState {
                module: &self.shaders.line,
                entry_point: Some("vs_main"),
                buffers: &[LineVertex::layout()],
                compilation_options: Default::default(),
            },
            fragment: Some(wgpu::FragmentState {
                module: &self.shaders.line,
                entry_point: Some("fs_main"),
                targets: &[Some(wgpu::ColorTargetState {
                    format: self.surface_format,
                    blend: get_blend_state(key.blend_mode),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
                compilation_options: Default::default(),
            }),
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::LineList,
                strip_index_format: None,
                front_face: wgpu::FrontFace::Ccw,
                cull_mode: None,
                unclipped_depth: false,
                polygon_mode: wgpu::PolygonMode::Fill,
                conservative: false,
            },
            depth_stencil: Some(depth_stencil_state(key.depth_write)),
            multisample: wgpu::MultisampleState {
                count: 1,
                mask: !0,
                alpha_to_coverage_enabled: false,
            },
            multiview: None,
            cache: None,
        })
    }

    /// Create the mesh render pipeline
    fn create_mesh_pipeline(&self, key: PipelineKey) -> wgpu::RenderPipeline {
        let pipeline_layout = self.device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("Mesh Pipeline Layout"),
            bind_group_layouts: &[&self.uniform_bind_group_layout],
            push_constant_ranges: &[],
        });

        self.device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("Mesh Pipeline"),
            layout: Some(&pipeline_layout),
            vertex: wgpu::VertexState {
                module: &self.shaders.mesh,
                entry_point: Some("vs_main"),
                buffers: &[MeshVertex::layout()],
                compilation_options: Default::default(),
            },
            fragment: Some(wgpu::FragmentState {
                module: &self.shaders.mesh,
                entry_point: Some("fs_main"),
                targets: &[Some(wgpu::ColorTargetState {
                    format: self.surface_format,
                    blend: get_blend_state(key.blend_mode),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
                compilation_options: Default::default(),
            }),
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::TriangleList,
                strip_index_format: None,
                front_face: wgpu::FrontFace::Ccw,
                cull_mode: None, // Two-sided rendering for surfaces (shader handles normal flipping)
                unclipped_depth: false,
                polygon_mode: wgpu::PolygonMode::Fill,
                conservative: false,
            },
            depth_stencil: Some(depth_stencil_state(key.depth_write)),
            multisample: wgpu::MultisampleState {
                count: 1,
                mask: !0,
                alpha_to_coverage_enabled: false,
            },
            multiview: None,
            cache: None,
        })
    }

    /// Create the dot render pipeline
    fn create_dot_pipeline(&self, key: PipelineKey) -> wgpu::RenderPipeline {
        let pipeline_layout = self.device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("Dot Pipeline Layout"),
            bind_group_layouts: &[&self.uniform_bind_group_layout],
            push_constant_ranges: &[],
        });

        // Billboard vertex layout for dots
        let billboard_layout = wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<[f32; 2]>() as wgpu::BufferAddress,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: &[wgpu::VertexAttribute {
                format: wgpu::VertexFormat::Float32x2,
                offset: 0,
                shader_location: 10,
            }],
        };

        self.device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("Dot Pipeline"),
            layout: Some(&pipeline_layout),
            vertex: wgpu::VertexState {
                module: &self.shaders.dot,
                entry_point: Some("vs_main"),
                buffers: &[billboard_layout, DotVertex::layout()],
                compilation_options: Default::default(),
            },
            fragment: Some(wgpu::FragmentState {
                module: &self.shaders.dot,
                entry_point: Some("fs_main"),
                targets: &[Some(wgpu::ColorTargetState {
                    format: self.surface_format,
                    blend: get_blend_state(key.blend_mode),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
                compilation_options: Default::default(),
            }),
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::TriangleList,
                strip_index_format: None,
                front_face: wgpu::FrontFace::Ccw,
                cull_mode: None,
                unclipped_depth: false,
                polygon_mode: wgpu::PolygonMode::Fill,
                conservative: false,
            },
            depth_stencil: Some(depth_stencil_state(key.depth_write)),
            multisample: wgpu::MultisampleState {
                count: 1,
                mask: !0,
                alpha_to_coverage_enabled: false,
            },
            multiview: None,
            cache: None,
        })
    }

    /// Create a depth texture for the given dimensions
    pub fn create_depth_texture(&self, width: u32, height: u32) -> wgpu::TextureView {
        let texture = self.device.create_texture(&wgpu::TextureDescriptor {
            label: Some("Depth Texture"),
            size: wgpu::Extent3d {
                width,
                height,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: self.depth_format,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::TEXTURE_BINDING,
            view_formats: &[],
        });

        texture.create_view(&wgpu::TextureViewDescriptor::default())
    }
}
