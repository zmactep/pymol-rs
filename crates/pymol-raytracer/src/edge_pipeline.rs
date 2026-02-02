//! Edge detection and compositing pipelines for ray trace modes 1, 2, 3

use crate::error::RaytraceResult;
use bytemuck::{Pod, Zeroable};

/// Uniform parameters for edge detection shader
/// 
/// Uses PyMOL's gradient-of-gradient algorithm parameters:
/// - slope_factor: Max gradient magnitude difference threshold (default 0.6)
/// - depth_factor: Max gradient direction difference threshold (default 0.1)  
/// - disco_factor: Gradient discontinuity threshold (default 0.05)
/// - gain: Pixel radius adjustment factor
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct EdgeParams {
    pub viewport: [f32; 2],
    pub slope_factor: f32,
    pub depth_factor: f32,
    pub disco_factor: f32,
    pub gain: f32,
    pub _pad: [f32; 2],
}

impl Default for EdgeParams {
    fn default() -> Self {
        Self {
            viewport: [0.0, 0.0],
            // Thresholds balanced for both surface (smooth) and cartoon (sharp) geometry
            slope_factor: 0.004,     // Gradient magnitude difference threshold
            depth_factor: 0.000005,  // Gradient direction difference threshold
            disco_factor: 0.25,      // Gradient discontinuity dot product threshold
            gain: 100.0,             // Gradient amplification
            _pad: [0.0; 2],
        }
    }
}

/// Uniform parameters for composite shader
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct CompositeParams {
    pub viewport: [f32; 2],
    pub mode: u32,
    pub use_transparent_bg: u32,
    pub edge_color: [f32; 3],
    pub quantize_levels: f32,
    pub bg_color: [f32; 3],
    pub _pad: f32,
}

impl Default for CompositeParams {
    fn default() -> Self {
        Self {
            viewport: [0.0, 0.0],
            mode: 0,
            use_transparent_bg: 0,
            edge_color: [0.0, 0.0, 0.0], // Black edges
            quantize_levels: 4.0,        // 4 levels for posterization
            bg_color: [0.0, 0.0, 0.0],   // Black background
            _pad: 0.0,
        }
    }
}

/// Edge detection pipeline (Pass 2)
pub struct EdgeDetectPipeline {
    bind_group_layout: wgpu::BindGroupLayout,
    compute_pipeline: wgpu::ComputePipeline,
}

impl EdgeDetectPipeline {
    /// Create the edge detection pipeline
    pub fn new(device: &wgpu::Device) -> RaytraceResult<Self> {
        // Create bind group layout for edge detection
        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("Edge Detect Bind Group Layout"),
            entries: &[
                // Depth texture (input)
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Texture {
                        sample_type: wgpu::TextureSampleType::Float { filterable: false },
                        view_dimension: wgpu::TextureViewDimension::D2,
                        multisampled: false,
                    },
                    count: None,
                },
                // Normal texture (input)
                wgpu::BindGroupLayoutEntry {
                    binding: 1,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Texture {
                        sample_type: wgpu::TextureSampleType::Float { filterable: false },
                        view_dimension: wgpu::TextureViewDimension::D2,
                        multisampled: false,
                    },
                    count: None,
                },
                // Edge texture (output) - use R32Float for storage binding support
                wgpu::BindGroupLayoutEntry {
                    binding: 2,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::StorageTexture {
                        access: wgpu::StorageTextureAccess::WriteOnly,
                        format: wgpu::TextureFormat::R32Float,
                        view_dimension: wgpu::TextureViewDimension::D2,
                    },
                    count: None,
                },
                // Edge params uniform
                wgpu::BindGroupLayoutEntry {
                    binding: 3,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
            ],
        });

        // Create shader module
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("Edge Detect Shader"),
            source: wgpu::ShaderSource::Wgsl(include_str!("shaders/edge_detect.wgsl").into()),
        });

        // Create pipeline layout
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("Edge Detect Pipeline Layout"),
            bind_group_layouts: &[&bind_group_layout],
            push_constant_ranges: &[],
        });

        // Create compute pipeline
        let compute_pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("Edge Detect Compute Pipeline"),
            layout: Some(&pipeline_layout),
            module: &shader,
            entry_point: Some("main"),
            compilation_options: Default::default(),
            cache: None,
        });

        Ok(Self {
            bind_group_layout,
            compute_pipeline,
        })
    }

    /// Get the bind group layout
    pub fn bind_group_layout(&self) -> &wgpu::BindGroupLayout {
        &self.bind_group_layout
    }

    /// Get the compute pipeline
    pub fn compute_pipeline(&self) -> &wgpu::ComputePipeline {
        &self.compute_pipeline
    }
}

/// Composite pipeline (Pass 3)
pub struct CompositePipeline {
    bind_group_layout: wgpu::BindGroupLayout,
    compute_pipeline: wgpu::ComputePipeline,
}

impl CompositePipeline {
    /// Create the composite pipeline
    pub fn new(device: &wgpu::Device) -> RaytraceResult<Self> {
        // Create bind group layout for compositing
        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("Composite Bind Group Layout"),
            entries: &[
                // Color texture (input)
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Texture {
                        sample_type: wgpu::TextureSampleType::Float { filterable: false },
                        view_dimension: wgpu::TextureViewDimension::D2,
                        multisampled: false,
                    },
                    count: None,
                },
                // Edge texture (input)
                wgpu::BindGroupLayoutEntry {
                    binding: 1,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Texture {
                        sample_type: wgpu::TextureSampleType::Float { filterable: false },
                        view_dimension: wgpu::TextureViewDimension::D2,
                        multisampled: false,
                    },
                    count: None,
                },
                // Depth texture (input - for background detection)
                wgpu::BindGroupLayoutEntry {
                    binding: 2,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Texture {
                        sample_type: wgpu::TextureSampleType::Float { filterable: false },
                        view_dimension: wgpu::TextureViewDimension::D2,
                        multisampled: false,
                    },
                    count: None,
                },
                // Output texture
                wgpu::BindGroupLayoutEntry {
                    binding: 3,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::StorageTexture {
                        access: wgpu::StorageTextureAccess::WriteOnly,
                        format: wgpu::TextureFormat::Rgba8Unorm,
                        view_dimension: wgpu::TextureViewDimension::D2,
                    },
                    count: None,
                },
                // Composite params uniform
                wgpu::BindGroupLayoutEntry {
                    binding: 4,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
            ],
        });

        // Create shader module
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("Composite Shader"),
            source: wgpu::ShaderSource::Wgsl(include_str!("shaders/composite.wgsl").into()),
        });

        // Create pipeline layout
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("Composite Pipeline Layout"),
            bind_group_layouts: &[&bind_group_layout],
            push_constant_ranges: &[],
        });

        // Create compute pipeline
        let compute_pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("Composite Compute Pipeline"),
            layout: Some(&pipeline_layout),
            module: &shader,
            entry_point: Some("main"),
            compilation_options: Default::default(),
            cache: None,
        });

        Ok(Self {
            bind_group_layout,
            compute_pipeline,
        })
    }

    /// Get the bind group layout
    pub fn bind_group_layout(&self) -> &wgpu::BindGroupLayout {
        &self.bind_group_layout
    }

    /// Get the compute pipeline
    pub fn compute_pipeline(&self) -> &wgpu::ComputePipeline {
        &self.compute_pipeline
    }
}
