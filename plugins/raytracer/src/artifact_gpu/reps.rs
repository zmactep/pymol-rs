use patinae_plugin::prelude::ViewerLike;
use patinae_scene::{GpuHandle, RenderArtifactRepDescriptor, RenderArtifactRepKind};

use super::layout::PrimitiveCounts;
use super::resources::direct_draw_args_buffer;

pub(super) struct SphereArtifactRep<'a> {
    pub(super) rep: &'a RenderArtifactRepDescriptor,
    pub(super) sphere_offset: u32,
    pub(super) instance_capacity: u32,
    pub(super) geometry_binding_size: u64,
    pub(super) rep_slot: u32,
}

impl SphereArtifactRep<'_> {
    pub(super) fn indirect_or_direct_args(
        &self,
        viewer: &mut dyn ViewerLike,
        label: &str,
    ) -> Result<GpuHandle, String> {
        if let Some(indirect) = self.rep.indirect {
            return Ok(indirect);
        }
        direct_draw_args_buffer(viewer, label, [0, self.instance_capacity, 0, 0])
    }
}

pub(super) struct CylinderArtifactRep<'a> {
    pub(super) rep: &'a RenderArtifactRepDescriptor,
    pub(super) cylinder_offset: u32,
    pub(super) instance_capacity: u32,
    pub(super) geometry_binding_size: u64,
    pub(super) rep_slot: u32,
    pub(super) radius: f32,
}

impl CylinderArtifactRep<'_> {
    pub(super) fn indirect_or_direct_args(
        &self,
        viewer: &mut dyn ViewerLike,
        label: &str,
    ) -> Result<GpuHandle, String> {
        if let Some(indirect) = self.rep.indirect {
            return Ok(indirect);
        }
        direct_draw_args_buffer(viewer, label, [0, self.instance_capacity, 0, 0])
    }
}

pub(super) struct CapsuleArtifactRep<'a> {
    pub(super) rep: &'a RenderArtifactRepDescriptor,
    pub(super) capsule_offset: u32,
    pub(super) instance_capacity: u32,
    pub(super) geometry_binding_size: u64,
    pub(super) rep_slot: u32,
}

impl CapsuleArtifactRep<'_> {
    pub(super) fn indirect_or_direct_args(
        &self,
        viewer: &mut dyn ViewerLike,
        label: &str,
    ) -> Result<GpuHandle, String> {
        if let Some(indirect) = self.rep.indirect {
            return Ok(indirect);
        }
        direct_draw_args_buffer(viewer, label, [0, self.instance_capacity, 0, 0])
    }
}

pub(super) struct TriangleArtifactRep<'a> {
    pub(super) rep: &'a RenderArtifactRepDescriptor,
    pub(super) triangle_offset: u32,
    pub(super) triangle_count: u32,
    pub(super) vertex_count: u32,
    pub(super) geometry_binding_size: u64,
    pub(super) rep_slot: u32,
    pub(super) visibility_counter_index: Option<u32>,
}

impl TriangleArtifactRep<'_> {
    pub(super) fn uses_surface_visibility_culling(&self) -> bool {
        self.rep.rep_kind == RenderArtifactRepKind::Surface
    }

    pub(super) fn source_triangle_count(&self) -> u32 {
        self.vertex_count / 3
    }
}

pub(super) struct ArtifactPlan<'a> {
    pub(super) color_lut: GpuHandle,
    pub(super) scene_atoms: Option<GpuHandle>,
    pub(super) sphere_reps: Vec<SphereArtifactRep<'a>>,
    pub(super) cylinder_reps: Vec<CylinderArtifactRep<'a>>,
    pub(super) capsule_reps: Vec<CapsuleArtifactRep<'a>>,
    pub(super) triangle_reps: Vec<TriangleArtifactRep<'a>>,
    pub(super) sphere_count: u32,
    pub(super) cylinder_count: u32,
    pub(super) capsule_count: u32,
    pub(super) triangle_count: u32,
}

impl ArtifactPlan<'_> {
    pub(super) fn primitive_counts(&self) -> PrimitiveCounts {
        PrimitiveCounts {
            spheres: self.sphere_count,
            cylinders: self.cylinder_count,
            capsules: self.capsule_count,
            triangles: self.triangle_count,
        }
    }

    pub(super) fn primitive_count(&self) -> Result<u32, String> {
        self.sphere_count
            .checked_add(self.cylinder_count)
            .and_then(|count| count.checked_add(self.capsule_count))
            .and_then(|count| count.checked_add(self.triangle_count))
            .ok_or_else(|| "artifact primitive count overflow".to_string())
    }
}
