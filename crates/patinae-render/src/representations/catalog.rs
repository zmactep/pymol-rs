//! Central registry for built-in render representations.
//!
//! The renderer does not expose a custom-representation extension point. This
//! table is the crate-local source of truth for built-in representation
//! construction and per-pass pipeline routing.

use patinae_mol::RepMask;

use crate::compute::cull::CullPipeline;
use crate::picking::pass::PickingPass;
use crate::picking::RepKind;
use crate::render_state::state::GeometryRuntime;
use crate::representations::cartoon::CartoonRep;
use crate::representations::dot::DotRep;
use crate::representations::ellipsoid::EllipsoidRep;
use crate::representations::line::LineRep;
use crate::representations::sphere::SphereRep;
use crate::representations::stick::StickRep;
use crate::representations::surface::SurfaceRep;
use crate::representations::{DrawPhase, Representation};

type RepConstructor = fn(&wgpu::Device) -> Box<dyn Representation>;

pub(crate) struct RepCatalogEntry {
    pub(crate) kind: RepKind,
    pub(crate) mask: RepMask,
    constructor: RepConstructor,
    pub(crate) cullable: bool,
    picking_scene_group: bool,
}

impl RepCatalogEntry {
    pub(crate) fn construct(&self, device: &wgpu::Device) -> Box<dyn Representation> {
        (self.constructor)(device)
    }

    pub(crate) fn color_pipeline<'a>(
        &self,
        geometry: &'a GeometryRuntime,
        phase: DrawPhase,
    ) -> Option<&'a wgpu::RenderPipeline> {
        match phase {
            DrawPhase::Opaque => self.opaque_pipeline(geometry),
            DrawPhase::Wboit => self.wboit_pipeline(geometry),
            DrawPhase::FastOverlay => self.fast_overlay_pipeline(geometry),
        }
    }

    pub(crate) fn picking_pipeline<'a>(
        &self,
        id_pass: &'a PickingPass,
    ) -> Option<&'a wgpu::RenderPipeline> {
        match self.kind {
            RepKind::Sphere => Some(&id_pass.sphere_pipeline),
            RepKind::Stick => Some(&id_pass.stick_pipeline),
            RepKind::Line => Some(&id_pass.line_pipeline),
            RepKind::Dot => Some(&id_pass.dot_pipeline),
            RepKind::Mesh => Some(&id_pass.std_vertex_line_pipeline),
            RepKind::Surface | RepKind::Cartoon | RepKind::Ribbon => {
                Some(&id_pass.std_vertex_pipeline)
            }
            RepKind::Ellipsoid => Some(&id_pass.ellipsoid_pipeline),
            _ => None,
        }
    }

    pub(crate) fn shadow_pipeline<'a>(
        &self,
        geometry: &'a GeometryRuntime,
    ) -> Option<&'a wgpu::RenderPipeline> {
        match self.kind {
            RepKind::Sphere => Some(&geometry.depth_prepass.sphere),
            RepKind::Stick => Some(&geometry.depth_prepass.stick),
            RepKind::Line => Some(&geometry.depth_prepass.line),
            RepKind::Dot => Some(&geometry.depth_prepass.dot),
            RepKind::Mesh => Some(&geometry.depth_prepass.mesh),
            RepKind::Surface => Some(&geometry.depth_prepass.surface),
            RepKind::Cartoon | RepKind::Ribbon => Some(&geometry.depth_prepass.cartoon),
            RepKind::Ellipsoid => Some(&geometry.depth_prepass.ellipsoid),
            _ => None,
        }
    }

    pub(crate) fn cull_pipeline<'a>(
        &self,
        cull: &'a CullPipeline,
    ) -> Option<&'a wgpu::ComputePipeline> {
        match self.kind {
            RepKind::Sphere => Some(&cull.pipeline_sphere),
            RepKind::Stick => Some(&cull.pipeline_stick),
            RepKind::Line => Some(&cull.pipeline_line),
            RepKind::Dot => Some(&cull.pipeline_dot),
            RepKind::Ellipsoid => Some(&cull.pipeline_ellipsoid),
            _ => None,
        }
    }

    pub(crate) fn cull_label(&self) -> &'static str {
        match self.kind {
            RepKind::Sphere => "patinae.cull.sphere",
            RepKind::Stick => "patinae.cull.stick",
            RepKind::Line => "patinae.cull.line",
            RepKind::Dot => "patinae.cull.dot",
            RepKind::Ellipsoid => "patinae.cull.ellipsoid",
            _ => "patinae.cull.unknown",
        }
    }

    pub(crate) fn picking_needs_scene_group(&self) -> bool {
        self.picking_scene_group
    }

    fn opaque_pipeline<'a>(
        &self,
        geometry: &'a GeometryRuntime,
    ) -> Option<&'a wgpu::RenderPipeline> {
        match self.kind {
            RepKind::Sphere => Some(&geometry.sphere_pipeline.pipeline_opaque),
            RepKind::Stick => Some(&geometry.stick_pipeline.pipeline_opaque),
            RepKind::Line => Some(&geometry.line_pipeline.pipeline_opaque),
            RepKind::Dot => Some(&geometry.dot_pipeline.pipeline_opaque),
            RepKind::Mesh => Some(&geometry.mesh_pipeline.pipeline_opaque),
            RepKind::Surface => Some(&geometry.surface_pipeline.pipeline_opaque),
            RepKind::Cartoon | RepKind::Ribbon => Some(&geometry.cartoon_pipeline.pipeline_opaque),
            RepKind::Ellipsoid => Some(&geometry.ellipsoid_pipeline.pipeline_opaque),
            _ => None,
        }
    }

    fn wboit_pipeline<'a>(
        &self,
        geometry: &'a GeometryRuntime,
    ) -> Option<&'a wgpu::RenderPipeline> {
        match self.kind {
            RepKind::Sphere => Some(&geometry.sphere_pipeline.pipeline),
            RepKind::Stick => Some(&geometry.stick_pipeline.pipeline),
            RepKind::Line => Some(&geometry.line_pipeline.pipeline),
            RepKind::Dot => Some(&geometry.dot_pipeline.pipeline),
            RepKind::Mesh => Some(&geometry.mesh_pipeline.pipeline),
            RepKind::Surface => Some(&geometry.surface_pipeline.pipeline),
            RepKind::Cartoon | RepKind::Ribbon => Some(&geometry.cartoon_pipeline.pipeline),
            RepKind::Ellipsoid => Some(&geometry.ellipsoid_pipeline.pipeline),
            _ => None,
        }
    }

    fn fast_overlay_pipeline<'a>(
        &self,
        geometry: &'a GeometryRuntime,
    ) -> Option<&'a wgpu::RenderPipeline> {
        match self.kind {
            RepKind::Dot => Some(&geometry.dot_pipeline.pipeline_fast_overlay),
            RepKind::Mesh => Some(&geometry.mesh_pipeline.pipeline_fast_overlay),
            _ => None,
        }
    }
}

pub(crate) fn active_entries(
    visible_reps: RepMask,
) -> impl Iterator<Item = &'static RepCatalogEntry> {
    REPS.iter()
        .filter(move |entry| visible_reps.is_visible(entry.mask))
}

pub(crate) fn entry(kind: RepKind) -> Option<&'static RepCatalogEntry> {
    REPS.iter().find(|entry| entry.kind == kind)
}

fn sphere(device: &wgpu::Device) -> Box<dyn Representation> {
    Box::new(SphereRep::new(device))
}

fn stick(device: &wgpu::Device) -> Box<dyn Representation> {
    Box::new(StickRep::new(device))
}

fn line(device: &wgpu::Device) -> Box<dyn Representation> {
    Box::new(LineRep::new(device))
}

fn dot(device: &wgpu::Device) -> Box<dyn Representation> {
    Box::new(DotRep::new(device))
}

fn surface(device: &wgpu::Device) -> Box<dyn Representation> {
    Box::new(SurfaceRep::new(device))
}

fn mesh(device: &wgpu::Device) -> Box<dyn Representation> {
    Box::new(SurfaceRep::new_mesh(device))
}

fn cartoon(device: &wgpu::Device) -> Box<dyn Representation> {
    Box::new(CartoonRep::new(device))
}

fn ribbon(device: &wgpu::Device) -> Box<dyn Representation> {
    Box::new(CartoonRep::new_ribbon(device))
}

fn ellipsoid(device: &wgpu::Device) -> Box<dyn Representation> {
    Box::new(EllipsoidRep::new(device))
}

pub(crate) const REPS: &[RepCatalogEntry] = &[
    RepCatalogEntry {
        kind: RepKind::Sphere,
        mask: RepMask::SPHERES,
        constructor: sphere,
        cullable: true,
        picking_scene_group: true,
    },
    RepCatalogEntry {
        kind: RepKind::Stick,
        mask: RepMask::STICKS,
        constructor: stick,
        cullable: true,
        picking_scene_group: false,
    },
    RepCatalogEntry {
        kind: RepKind::Line,
        mask: RepMask::LINES,
        constructor: line,
        cullable: true,
        picking_scene_group: false,
    },
    RepCatalogEntry {
        kind: RepKind::Dot,
        mask: RepMask::DOTS,
        constructor: dot,
        cullable: true,
        picking_scene_group: false,
    },
    RepCatalogEntry {
        kind: RepKind::Surface,
        mask: RepMask::SURFACE,
        constructor: surface,
        cullable: false,
        picking_scene_group: false,
    },
    RepCatalogEntry {
        kind: RepKind::Mesh,
        mask: RepMask::MESH,
        constructor: mesh,
        cullable: false,
        picking_scene_group: false,
    },
    RepCatalogEntry {
        kind: RepKind::Cartoon,
        mask: RepMask::CARTOON,
        constructor: cartoon,
        cullable: false,
        picking_scene_group: false,
    },
    RepCatalogEntry {
        kind: RepKind::Ribbon,
        mask: RepMask::RIBBON,
        constructor: ribbon,
        cullable: false,
        picking_scene_group: false,
    },
    RepCatalogEntry {
        kind: RepKind::Ellipsoid,
        mask: RepMask::ELLIPSOIDS,
        constructor: ellipsoid,
        cullable: true,
        picking_scene_group: false,
    },
];
