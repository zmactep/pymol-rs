//! Pipeline management

/// Key for identifying a pipeline configuration
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct PipelineKey {
    pub pipeline_type: PipelineType,
    pub blend_mode: BlendMode,
    pub depth_write: bool,
}

/// Types of render pipelines
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PipelineType {
    Sphere,
    Cylinder,
    Line,
    Mesh,
    Dot,
}

/// Blend modes for rendering
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BlendMode {
    Opaque,
    Transparent,
}

impl Default for BlendMode {
    fn default() -> Self {
        BlendMode::Opaque
    }
}

/// Get the blend state for a given blend mode
pub fn get_blend_state(mode: BlendMode) -> Option<wgpu::BlendState> {
    match mode {
        BlendMode::Opaque => None,
        BlendMode::Transparent => Some(wgpu::BlendState {
            color: wgpu::BlendComponent {
                src_factor: wgpu::BlendFactor::SrcAlpha,
                dst_factor: wgpu::BlendFactor::OneMinusSrcAlpha,
                operation: wgpu::BlendOperation::Add,
            },
            alpha: wgpu::BlendComponent {
                src_factor: wgpu::BlendFactor::One,
                dst_factor: wgpu::BlendFactor::OneMinusSrcAlpha,
                operation: wgpu::BlendOperation::Add,
            },
        }),
    }
}

/// Common depth stencil state for 3D rendering
pub fn depth_stencil_state(depth_write: bool) -> wgpu::DepthStencilState {
    wgpu::DepthStencilState {
        format: wgpu::TextureFormat::Depth32Float,
        depth_write_enabled: depth_write,
        depth_compare: wgpu::CompareFunction::Less,
        stencil: wgpu::StencilState::default(),
        bias: wgpu::DepthBiasState::default(),
    }
}
