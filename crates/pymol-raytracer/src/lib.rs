//! GPU compute shader-based raytracer for PyMOL-RS
//!
//! This crate provides high-quality offline rendering of molecular structures
//! using GPU compute shaders for raytracing.
//!
//! # Architecture
//!
//! The raytracer uses a multi-pass compute shader pipeline:
//! 1. Primary ray generation and BVH traversal
//! 2. Shading with shadows and lighting
//! 3. Transparency handling with depth sorting
//! 4. Optional antialiasing via supersampling
//!
//! # Example
//!
//! ```ignore
//! use pymol_raytracer::{RayContext, RaytraceParams, PrimitiveCollector};
//!
//! // Create raytracing context from existing wgpu device/queue
//! let ray_ctx = RayContext::new(&device, &queue);
//!
//! // Collect primitives from scene
//! let primitives = PrimitiveCollector::new()
//!     .collect_spheres(&sphere_data)
//!     .collect_cylinders(&cylinder_data)
//!     .build();
//!
//! // Raytrace
//! let params = RaytraceParams::new(1920, 1080);
//! let image_data = ray_ctx.raytrace(&primitives, &params)?;
//! ```

pub mod bvh;
pub mod collect;
pub mod context;
pub mod edge_pipeline;
pub mod error;
pub mod pipeline;
pub mod primitive;
pub mod render;

// Re-exports
pub use bvh::{Bvh, BvhNode};
pub use collect::{collect_from_molecule, CollectOptions, RayColorResolver};
pub use context::raytrace;
pub use error::{RaytraceError, RaytraceResult};
pub use primitive::{GpuCylinder, GpuSphere, GpuTriangle, PrimitiveCollector, Primitives};
pub use render::{RaytraceParams, RaytraceSettings};

/// Convert a lin_alg Mat4 to a [[f32; 4]; 4] array for raytracing uniforms
///
/// Uses the same conversion as pymol-render's GlobalUniforms for consistency.
#[inline]
pub fn matrix_to_array(m: lin_alg::f32::Mat4) -> [[f32; 4]; 4] {
    <[[f32; 4]; 4]>::from(m)
}
