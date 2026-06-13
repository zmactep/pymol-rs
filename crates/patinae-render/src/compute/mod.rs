//! GPU compute pipelines used by representations whose builds satisfy the
//! GPU-first criteria. Each submodule exposes a
//! `dispatch(...)` helper called from a representation's `build()` via the
//! encoder it already holds — no separate buffer marshalling.

pub mod cartoon_extrude;
pub mod cull;
pub mod dot_build;
pub mod ellipsoid_build;
pub mod line_build;
pub mod sphere_build;
pub mod sphere_lod_count;
pub mod ssao;
pub mod stick_build;
pub mod stick_lod_count;
pub mod surface_density;
pub mod surface_mc;
pub mod surface_ses_morph;
pub mod surface_vdw_sdf;

/// WebGPU caps each `dispatch_workgroups` dimension at 65535. For 1D
/// dispatches with > 65535 workgroups (assemblies with millions of atoms /
/// bonds / vertices), split into a 2D grid.
///
/// Resulting dispatch dimensions cover at least `total_workgroups`
/// workgroups; shaders must guard with a linear-index bounds check since
/// `wg_x * wg_y` may exceed `total_workgroups` at the tail. Shader-side
/// pattern (replacing `let i = gid.x;`):
///
/// ```wgsl
/// @compute @workgroup_size(64)
/// fn cs_main(
///     @builtin(global_invocation_id) gid: vec3<u32>,
///     @builtin(num_workgroups) nwg: vec3<u32>,
/// ) {
///     let i = gid.x + gid.y * nwg.x * 64u;
///     if i >= count { return; }
///     ...
/// }
/// ```
pub const MAX_WORKGROUPS_PER_DIM: u32 = 65_535;

/// Returns `(wg_x, wg_y)` covering `total_workgroups` workgroups within
/// the WebGPU `max_compute_workgroups_per_dimension = 65535` limit.
#[inline]
pub fn split_1d_dispatch(total_workgroups: u32) -> (u32, u32) {
    if total_workgroups <= MAX_WORKGROUPS_PER_DIM {
        return (total_workgroups.max(1), 1);
    }
    let wg_x = MAX_WORKGROUPS_PER_DIM;
    let wg_y = total_workgroups.div_ceil(MAX_WORKGROUPS_PER_DIM);
    (wg_x, wg_y)
}

/// Bytes laid out as [`wgpu::util::DrawIndirectArgs`]. The host primes
/// `vertex_count` plus three zero words before each rep's build /
/// cull dispatch; the GPU writes `instance_count` via atomic.
#[inline]
pub fn indirect_seed(vertex_count: u32) -> [u32; 4] {
    [vertex_count, 0, 0, 0]
}

/// Bundle of compute pipeline + its build-side BGL. Used by every per-rep
/// build kernel (sphere / stick / line / dot / ellipsoid).
pub(crate) struct ComputeBuildPipeline {
    pub pipeline: wgpu::ComputePipeline,
    pub build_layout: wgpu::BindGroupLayout,
}

/// The standard 3-binding BGL shared by sphere/stick/line/dot build
/// pipelines: `params` uniform, `raw_instance` storage (RW), `raw_count`
/// atomic counter (RW).
///
/// `instance_min_binding` is the WGSL-aligned size of one instance entry
/// (sphere: 32, stick/line: 48, dot: 32 — driven by struct layout).
pub(crate) fn build_3binding_build_layout(
    device: &wgpu::Device,
    label: &str,
    params_size: u64,
    instance_min_binding: u64,
) -> wgpu::BindGroupLayout {
    use std::num::NonZeroU64;
    device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        label: Some(label),
        entries: &[
            wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::COMPUTE,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: NonZeroU64::new(params_size),
                },
                count: None,
            },
            wgpu::BindGroupLayoutEntry {
                binding: 1,
                visibility: wgpu::ShaderStages::COMPUTE,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Storage { read_only: false },
                    has_dynamic_offset: false,
                    min_binding_size: NonZeroU64::new(instance_min_binding),
                },
                count: None,
            },
            wgpu::BindGroupLayoutEntry {
                binding: 2,
                visibility: wgpu::ShaderStages::COMPUTE,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Storage { read_only: false },
                    has_dynamic_offset: false,
                    // raw_count: array<atomic<u32>, 1> — 4 B
                    min_binding_size: NonZeroU64::new(4),
                },
                count: None,
            },
        ],
    })
}

/// Compile a per-rep compute build pipeline. `label_prefix` follows the
/// `patinae.<rep>_build` convention; downstream labels are
/// `<prefix>.pipeline_layout`, `<prefix>.wgsl`, `<prefix>.pipeline`.
/// The pipeline is bound over `[scene_layout, build_layout]`.
pub(crate) fn make_compute_build_pipeline(
    device: &wgpu::Device,
    scene_layout: &crate::scene_store::SceneStoreLayout,
    label_prefix: &str,
    wgsl: &str,
    build_layout: wgpu::BindGroupLayout,
) -> ComputeBuildPipeline {
    let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
        label: Some(&format!("{label_prefix}.pipeline_layout")),
        bind_group_layouts: &[&scene_layout.bind_group_layout, &build_layout],
        immediate_size: 0,
    });

    let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some(&format!("{label_prefix}.wgsl")),
        source: wgpu::ShaderSource::Wgsl(crate::shader_source::expand(wgsl).into()),
    });

    let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
        label: Some(&format!("{label_prefix}.pipeline")),
        layout: Some(&pipeline_layout),
        module: &module,
        entry_point: Some("cs_main"),
        compilation_options: wgpu::PipelineCompilationOptions::default(),
        cache: None,
    });

    ComputeBuildPipeline {
        pipeline,
        build_layout,
    }
}
