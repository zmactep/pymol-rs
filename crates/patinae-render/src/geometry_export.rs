//! Renderer-neutral export of currently displayed molecular geometry.
//!
//! This module deliberately exports display geometry, not raytracer
//! primitives. Downstream consumers can adapt it to ray tracing, Blender,
//! VR/AR, debug visualisation, or file formats without coupling
//! `patinae-render` to any one consumer.

use thiserror::Error;

use crate::picking::{ObjectId, RepKind};

pub(crate) mod analytic;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GeometryExportOptions {
    /// Include renderer mesh outputs such as cartoon/ribbon/surface.
    pub include_meshes: bool,
    /// Include analytic atoms/bonds where the renderer semantics are
    /// naturally analytic.
    pub include_analytic: bool,
    /// Include semantic samples for screen-space representations.
    pub include_semantic_samples: bool,
}

impl Default for GeometryExportOptions {
    fn default() -> Self {
        Self {
            include_meshes: true,
            include_analytic: true,
            include_semantic_samples: true,
        }
    }
}

/// Whole-scene displayed geometry grouped by renderer object id.
#[derive(Debug, Clone, Default)]
pub struct DisplayedGeometry {
    pub objects: Vec<DisplayedObjectGeometry>,
}

impl DisplayedGeometry {
    pub fn is_empty(&self) -> bool {
        self.objects.iter().all(|obj| obj.primitives.is_empty())
    }

    pub fn primitive_count(&self) -> usize {
        self.objects.iter().map(|obj| obj.primitives.len()).sum()
    }
}

/// Display geometry for one renderer object.
#[derive(Debug, Clone)]
pub struct DisplayedObjectGeometry {
    pub object_id: ObjectId,
    pub primitives: Vec<DisplayedPrimitive>,
}

/// A renderer-neutral displayed primitive.
#[derive(Debug, Clone)]
pub enum DisplayedPrimitive {
    /// Non-indexed triangle-list mesh. Used for Cartoon, Ribbon, Surface.
    Mesh { rep: RepKind, mesh: DisplayedMesh },
    /// Analytic sphere. Used for Spheres and also stick end caps.
    Sphere {
        rep: RepKind,
        owner_atom_id: u32,
        center: [f32; 3],
        radius: f32,
        material: DisplayedMaterial,
    },
    /// Analytic cylinder body. Consumers that need capsules can pair stick
    /// cylinders with the exported stick cap spheres.
    Cylinder {
        rep: RepKind,
        owner_atom_ids: [u32; 2],
        start: [f32; 3],
        end: [f32; 3],
        radius: f32,
        material_start: DisplayedMaterial,
        material_end: DisplayedMaterial,
    },
    /// Semantic line segment for screen-space line or mesh-wire displays.
    LineSegment {
        rep: RepKind,
        owner_atom_ids: [u32; 2],
        start: [f32; 3],
        end: [f32; 3],
        width_px: f32,
        material_start: DisplayedMaterial,
        material_end: DisplayedMaterial,
    },
    /// Semantic dot sample for screen-space dot displays.
    PointSample {
        rep: RepKind,
        owner_atom_id: u32,
        position: [f32; 3],
        radius_px: f32,
        material: DisplayedMaterial,
    },
}

impl DisplayedPrimitive {
    pub fn rep_kind(&self) -> RepKind {
        match self {
            DisplayedPrimitive::Mesh { rep, .. }
            | DisplayedPrimitive::Sphere { rep, .. }
            | DisplayedPrimitive::Cylinder { rep, .. }
            | DisplayedPrimitive::LineSegment { rep, .. }
            | DisplayedPrimitive::PointSample { rep, .. } => *rep,
        }
    }
}

/// Non-indexed triangle-list mesh data.
#[derive(Debug, Clone, Default)]
pub struct DisplayedMesh {
    pub vertices: Vec<DisplayedMeshVertex>,
}

/// Vertex emitted by a displayed mesh.
#[derive(Debug, Clone)]
pub struct DisplayedMeshVertex {
    pub position: [f32; 3],
    pub normal: [f32; 3],
    pub owner_atom_id: u32,
    pub material: DisplayedMaterial,
    pub flags: u32,
}

/// Renderer-resolved material for an exported primitive or vertex.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DisplayedMaterial {
    /// Host-resolved atom base colour before representation override/alpha.
    pub base_rgba: [f32; 4],
    /// Representation colour after rep-specific override, before alpha.
    pub rep_rgba: [f32; 4],
    /// Final display colour after per-rep/per-atom alpha resolution.
    pub rgba: [f32; 4],
    /// Transparency in PyMOL terms: 0 = opaque, 1 = fully transparent.
    pub transparency: f32,
}

impl DisplayedMaterial {
    pub fn from_rgba(rgba: [f32; 4]) -> Self {
        Self {
            base_rgba: rgba,
            rep_rgba: rgba,
            rgba,
            transparency: 1.0 - rgba[3].clamp(0.0, 1.0),
        }
    }
}

/// Errors from blocking GPU readback or export conversion.
#[derive(Debug, Error)]
pub enum GeometryExportError {
    #[error("GPU readback failed: {0}")]
    Gpu(String),
}
