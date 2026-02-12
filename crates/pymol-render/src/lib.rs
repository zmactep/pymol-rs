//! pymol-render: wgpu-based rendering engine for PyMOL-RS
//!
//! This crate provides GPU-accelerated molecular visualization using modern
//! graphics APIs via wgpu. It supports various molecular representations:
//!
//! - **Spheres**: Van der Waals spheres using impostor shaders
//! - **Sticks**: Cylinder-based bond representation using impostor shaders
//! - **Lines**: Simple line-based bond representation
//! - **Dots**: Dot surface representation
//! - **Mesh**: Triangle mesh representation with lighting
//! - **Cartoon**: Secondary structure visualization (helices, sheets, loops)
//! - **Surface**: Molecular surface (VdW, SAS, SES) via marching cubes
//!
//! # Architecture
//!
//! The crate is organized around these key concepts:
//!
//! - [`RenderContext`]: Manages wgpu device, queue, and shared GPU resources
//! - [`Representation`]: Trait for generating and rendering molecular representations
//! - [`ColorResolver`]: Resolves color indices to final RGBA values
//! - [`GlobalUniforms`]: Camera, lighting, and fog parameters for shaders
//!
//! # Example
//!
//! ```ignore
//! use pymol_render::{RenderContext, SphereRep, CartoonRep, ColorResolver};
//! use pymol_mol::ObjectMolecule;
//!
//! // Create render context
//! let context = RenderContext::new(&device, &queue, surface_format);
//!
//! // Build sphere representation
//! let mut spheres = SphereRep::new();
//! spheres.build(&molecule, &coord_set, &color_resolver, &settings);
//! spheres.upload(&device, &queue);
//!
//! // Build cartoon representation
//! let mut cartoon = CartoonRep::new();
//! cartoon.build(&molecule, &coord_set, &color_resolver, &settings);
//! cartoon.upload(&device, &queue);
//!
//! // Render in a render pass
//! spheres.render(&mut render_pass);
//! cartoon.render(&mut render_pass);
//! ```

mod buffer;
mod color_resolver;
mod context;
mod error;
pub mod picking;
pub mod pipeline;
mod representation;
mod uniforms;
mod vertex;

// Re-export main types
pub use color_resolver::ColorResolver;
pub use context::RenderContext;
pub use error::RenderError;
pub use representation::{RepType, Representation};
pub use uniforms::GlobalUniforms;
pub use vertex::{CylinderVertex, DotVertex, LineVertex, MeshVertex, SphereVertex};

// Re-export specific representations
pub use representation::cartoon::CartoonRep;
pub use representation::dot::DotRep;
pub use representation::line::LineRep;
pub use representation::mesh::MeshRep;
pub use representation::ribbon::RibbonRep;
pub use representation::selection_indicator::{
    SelectionIndicatorRep, DEFAULT_INDICATOR_SIZE, SELECTION_INDICATOR_COLOR,
};
pub use representation::sphere::SphereRep;
pub use representation::stick::StickRep;
pub use representation::surface::{Grid3D, SurfaceRep, SurfaceType};
pub use representation::wire_surface::WireSurfaceRep;
