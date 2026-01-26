//! Pipeline management and caching

use parking_lot::RwLock;
use std::collections::HashMap;

/// Cache for render pipelines
pub struct PipelineCache {
    #[allow(dead_code)]
    pipelines: RwLock<HashMap<PipelineKey, wgpu::RenderPipeline>>,
}

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

impl PipelineCache {
    /// Create a new empty pipeline cache
    pub fn new() -> Self {
        Self {
            pipelines: RwLock::new(HashMap::new()),
        }
    }

    /// Get a pipeline from the cache, or create it if it doesn't exist
    #[allow(dead_code)]
    pub fn get_or_create<F>(&self, key: PipelineKey, create_fn: F) -> wgpu::RenderPipeline
    where
        F: FnOnce() -> wgpu::RenderPipeline,
    {
        // Check if pipeline exists
        {
            let pipelines = self.pipelines.read();
            if let Some(pipeline) = pipelines.get(&key) {
                // Clone the pipeline handle (it's an Arc internally)
                return unsafe { std::ptr::read(pipeline) };
            }
        }

        // Create new pipeline
        let pipeline = create_fn();

        // Store in cache
        {
            let mut pipelines = self.pipelines.write();
            pipelines.insert(key, unsafe { std::ptr::read(&pipeline) });
        }

        pipeline
    }

    /// Clear the pipeline cache
    #[allow(dead_code)]
    pub fn clear(&self) {
        self.pipelines.write().clear();
    }
}

impl Default for PipelineCache {
    fn default() -> Self {
        Self::new()
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
