use patinae_plugin::prelude::ViewerLike;
use patinae_scene::{
    GpuBatchCommand, GpuHandle, RenderArtifactRepDescriptor, RenderArtifactRepKind,
};

use super::layout::PrimitiveCounts;
use super::resources::{create_buffer, direct_draw_args_buffer, storage_usage};

const DRAW_ARGS_BYTES: u64 = std::mem::size_of::<[u32; 4]>() as u64;

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
        batch_commands: &mut Vec<GpuBatchCommand>,
    ) -> Result<GpuHandle, String> {
        if let Some(indirect) = self.rep.indirect {
            return shader_draw_args_buffer(viewer, label, indirect, batch_commands);
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
        batch_commands: &mut Vec<GpuBatchCommand>,
    ) -> Result<GpuHandle, String> {
        if let Some(indirect) = self.rep.indirect {
            return shader_draw_args_buffer(viewer, label, indirect, batch_commands);
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
        batch_commands: &mut Vec<GpuBatchCommand>,
    ) -> Result<GpuHandle, String> {
        if let Some(indirect) = self.rep.indirect {
            return shader_draw_args_buffer(viewer, label, indirect, batch_commands);
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
    pub(super) fn indirect_or_direct_args(
        &self,
        viewer: &mut dyn ViewerLike,
        label: &str,
        batch_commands: &mut Vec<GpuBatchCommand>,
    ) -> Result<GpuHandle, String> {
        if let Some(indirect) = self.rep.indirect {
            return shader_draw_args_buffer(viewer, label, indirect, batch_commands);
        }
        direct_draw_args_buffer(viewer, label, [self.vertex_count, 1, 0, 0])
    }

    pub(super) fn uses_surface_visibility_culling(&self) -> bool {
        self.rep.rep_kind == RenderArtifactRepKind::Surface
    }

    pub(super) fn source_triangle_count(&self) -> u32 {
        self.vertex_count / 3
    }
}

fn shader_draw_args_buffer(
    viewer: &mut dyn ViewerLike,
    label: &str,
    indirect: GpuHandle,
    batch_commands: &mut Vec<GpuBatchCommand>,
) -> Result<GpuHandle, String> {
    let shader_args = create_buffer(viewer, label, DRAW_ARGS_BYTES, storage_usage(), None)?;
    append_shader_draw_args_copy(indirect, shader_args, batch_commands);
    Ok(shader_args)
}

fn append_shader_draw_args_copy(
    indirect: GpuHandle,
    shader_args: GpuHandle,
    batch_commands: &mut Vec<GpuBatchCommand>,
) {
    batch_commands.push(GpuBatchCommand::CopyBufferToBuffer {
        source: indirect,
        source_offset: 0,
        destination: shader_args,
        destination_offset: 0,
        size: DRAW_ARGS_BYTES,
    });
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

#[cfg(test)]
mod tests {
    use super::*;
    use patinae_scene::GpuHandleKind;

    #[test]
    fn indirect_draw_args_copy_is_recorded_for_shader_buffer() {
        let indirect = GpuHandle {
            id: 99,
            kind: GpuHandleKind::Buffer,
            generation: 1,
        };
        let shader_args = GpuHandle {
            id: 7,
            kind: GpuHandleKind::Buffer,
            generation: 1,
        };
        let mut commands = Vec::new();

        append_shader_draw_args_copy(indirect, shader_args, &mut commands);

        assert!(matches!(
            commands.as_slice(),
            [GpuBatchCommand::CopyBufferToBuffer {
                source,
                source_offset: 0,
                destination,
                destination_offset: 0,
                size: DRAW_ARGS_BYTES,
            }] if *source == indirect && *destination == shader_args
        ));
    }
}
