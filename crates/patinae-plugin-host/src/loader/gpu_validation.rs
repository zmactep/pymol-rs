use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use super::*;

pub(super) fn semantic_hash<T: Hash>(value: &T) -> u64 {
    let mut hasher = DefaultHasher::new();
    value.hash(&mut hasher);
    hasher.finish()
}

pub(super) fn fingerprint_for(kind: GpuHandleKind, value: impl Hash) -> GpuResourceFingerprint {
    let mut hasher = DefaultHasher::new();
    GPU_CACHE_LAYOUT_VERSION.hash(&mut hasher);
    kind.hash(&mut hasher);
    value.hash(&mut hasher);
    GpuResourceFingerprint {
        kind,
        hash: hasher.finish(),
    }
}

pub(super) fn shader_module_fingerprint(
    descriptor: &GpuShaderModuleDescriptor,
) -> GpuResourceFingerprint {
    fingerprint_for(GpuHandleKind::ShaderModule, &descriptor.wgsl)
}

pub(super) fn bind_group_layout_fingerprint(
    descriptor: &GpuBindGroupLayoutDescriptor,
) -> GpuResourceFingerprint {
    fingerprint_for(GpuHandleKind::BindGroupLayout, &descriptor.entries)
}

pub(super) fn pipeline_layout_fingerprint(
    bind_group_layouts: &[GpuResourceFingerprint],
) -> GpuResourceFingerprint {
    fingerprint_for(GpuHandleKind::PipelineLayout, bind_group_layouts)
}

pub(super) fn compute_pipeline_fingerprint(
    layout: GpuResourceFingerprint,
    module: GpuResourceFingerprint,
    entry_point: &str,
) -> GpuResourceFingerprint {
    fingerprint_for(
        GpuHandleKind::ComputePipeline,
        (layout, module, entry_point),
    )
}

pub(super) struct RenderPipelineFingerprintInput<'a> {
    pub(super) layout: GpuResourceFingerprint,
    pub(super) vertex_module: GpuResourceFingerprint,
    pub(super) vertex_entry_point: &'a str,
    pub(super) vertex_buffers: &'a [GpuVertexBufferLayout],
    pub(super) fragment: Option<(
        GpuResourceFingerprint,
        &'a str,
        &'a Vec<Option<GpuColorTargetState>>,
    )>,
    pub(super) primitive: GpuPrimitiveState,
    pub(super) depth_stencil: Option<GpuDepthStencilState>,
    pub(super) multisample: GpuMultisampleState,
}

pub(super) fn render_pipeline_fingerprint(
    input: RenderPipelineFingerprintInput<'_>,
) -> GpuResourceFingerprint {
    let RenderPipelineFingerprintInput {
        layout,
        vertex_module,
        vertex_entry_point,
        vertex_buffers,
        fragment,
        primitive,
        depth_stencil,
        multisample,
    } = input;

    fingerprint_for(
        GpuHandleKind::RenderPipeline,
        (
            layout,
            vertex_module,
            vertex_entry_point,
            vertex_buffers,
            fragment,
            primitive,
            depth_stencil,
            multisample,
        ),
    )
}

pub(super) fn render_pipeline_metadata(
    descriptor: &GpuRenderPipelineDescriptor,
) -> RenderPipelineMetadata {
    RenderPipelineMetadata {
        color_targets: descriptor
            .fragment
            .as_ref()
            .map(|fragment| {
                fragment
                    .targets
                    .iter()
                    .map(|target| target.map(|target| target.format))
                    .collect()
            })
            .unwrap_or_default(),
        depth_stencil: descriptor.depth_stencil.as_ref().map(|depth| depth.format),
    }
}

pub(super) fn validate_render_pipeline_descriptor(
    descriptor: &GpuRenderPipelineDescriptor,
) -> Result<(), String> {
    if descriptor.vertex.entry_point.is_empty() {
        return Err("GPU render pipeline vertex entry point must not be empty".to_string());
    }
    validate_vertex_state(&descriptor.vertex)?;
    for layout in &descriptor.vertex.buffers {
        validate_vertex_buffer_layout(layout)?;
    }
    if let Some(fragment) = &descriptor.fragment {
        if fragment.entry_point.is_empty() {
            return Err("GPU render pipeline fragment entry point must not be empty".to_string());
        }
        for target in fragment.targets.iter().flatten() {
            validate_color_target_state(*target)?;
        }
    }
    if let Some(depth) = descriptor.depth_stencil {
        validate_depth_format(depth.format)?;
    }
    if descriptor.multisample.count != 1 {
        return Err("GPU render pipeline v1 requires multisample count 1".to_string());
    }
    if descriptor.multisample.alpha_to_coverage_enabled {
        return Err("GPU render pipeline v1 does not support alpha-to-coverage".to_string());
    }
    Ok(())
}

pub(super) fn validate_vertex_state(state: &GpuVertexState) -> Result<(), String> {
    if state.module.kind != GpuHandleKind::ShaderModule {
        return Err("GPU render pipeline vertex module must be a shader module handle".to_string());
    }
    Ok(())
}

pub(super) fn validate_vertex_buffer_layout(layout: &GpuVertexBufferLayout) -> Result<(), String> {
    if layout.array_stride == 0 {
        return Err("GPU render pipeline vertex array_stride must be non-zero".to_string());
    }
    if !layout.array_stride.is_multiple_of(wgpu::VERTEX_ALIGNMENT) {
        return Err(format!(
            "GPU render pipeline vertex array_stride must be aligned to {}",
            wgpu::VERTEX_ALIGNMENT
        ));
    }
    for attribute in &layout.attributes {
        let size = gpu_vertex_format_size(attribute.format);
        let end = attribute
            .offset
            .checked_add(size)
            .ok_or_else(|| "GPU render pipeline vertex attribute range overflow".to_string())?;
        if end > layout.array_stride {
            return Err(format!(
                "GPU render pipeline vertex attribute at location {} exceeds array_stride",
                attribute.shader_location
            ));
        }
    }
    Ok(())
}

pub(super) fn validate_color_target_state(target: GpuColorTargetState) -> Result<(), String> {
    validate_color_format(target.format)?;
    validate_color_write_mask(target.write_mask)
}

pub(super) fn validate_color_format(format: GpuTextureFormat) -> Result<(), String> {
    if matches!(format, GpuTextureFormat::Depth32Float) {
        return Err("GPU render color targets cannot use depth formats".to_string());
    }
    Ok(())
}

pub(super) fn validate_depth_format(format: GpuTextureFormat) -> Result<(), String> {
    if format != GpuTextureFormat::Depth32Float {
        return Err("GPU render depth targets require Depth32Float format".to_string());
    }
    Ok(())
}

pub(super) fn validate_color_write_mask(mask: GpuColorWriteMask) -> Result<(), String> {
    if mask.bits & !GpuColorWriteMask::ALL.bits != 0 {
        return Err(format!("unknown GPU color write mask bits: {}", mask.bits));
    }
    Ok(())
}

pub(super) fn bind_group_layout_has_dynamic_offsets(entries: &[GpuBindGroupLayoutEntry]) -> bool {
    entries.iter().any(|entry| {
        matches!(
            entry.ty,
            GpuBindingType::Buffer {
                has_dynamic_offset: true,
                ..
            }
        )
    })
}
