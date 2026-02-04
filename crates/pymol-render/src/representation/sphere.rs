//! Sphere representation for VdW sphere visualization
//!
//! Renders atoms as spheres using impostor shaders.

use crate::buffer::GrowableBuffer;
use crate::color_resolver::ColorResolver;
use crate::representation::Representation;
use crate::vertex::SphereVertex;

use pymol_mol::{CoordSet, ObjectMolecule, RepMask};
use pymol_settings::SettingResolver;

/// Sphere representation for atoms
///
/// Renders atoms as Van der Waals spheres using impostor shaders.
/// Each sphere is rendered as a billboard quad with ray-sphere intersection
/// computed in the fragment shader.
pub struct SphereRep {
    /// Instance data (CPU side)
    instances: Vec<SphereVertex>,
    /// GPU instance buffer
    instance_buffer: GrowableBuffer,
    /// Pipeline for rendering
    pipeline: Option<wgpu::RenderPipeline>,
    /// Billboard vertex buffer (shared)
    billboard_buffer: Option<wgpu::Buffer>,
    /// Quad index buffer (shared)
    index_buffer: Option<wgpu::Buffer>,
    /// Uniform bind group
    bind_group: Option<wgpu::BindGroup>,
    /// Whether the representation needs to be rebuilt
    dirty: bool,
    /// Number of sphere instances to render
    instance_count: u32,
    /// Sphere scale factor
    sphere_scale: f32,
}

impl SphereRep {
    /// Create a new sphere representation
    pub fn new() -> Self {
        Self {
            instances: Vec::new(),
            instance_buffer: GrowableBuffer::new("Sphere Instances", wgpu::BufferUsages::VERTEX),
            pipeline: None,
            billboard_buffer: None,
            index_buffer: None,
            bind_group: None,
            dirty: true,
            instance_count: 0,
            sphere_scale: 1.0,
        }
    }

    /// Set the render pipeline and shared buffers
    pub fn set_pipeline(
        &mut self,
        pipeline: wgpu::RenderPipeline,
        billboard_buffer: wgpu::Buffer,
        index_buffer: wgpu::Buffer,
        bind_group: wgpu::BindGroup,
    ) {
        self.pipeline = Some(pipeline);
        self.billboard_buffer = Some(billboard_buffer);
        self.index_buffer = Some(index_buffer);
        self.bind_group = Some(bind_group);
    }

    /// Set the sphere scale factor
    pub fn set_scale(&mut self, scale: f32) {
        self.sphere_scale = scale;
    }
}

impl Default for SphereRep {
    fn default() -> Self {
        Self::new()
    }
}

impl Representation for SphereRep {
    fn build(
        &mut self,
        molecule: &ObjectMolecule,
        coord_set: &CoordSet,
        colors: &ColorResolver,
        settings: &SettingResolver,
    ) {
        self.instances.clear();

        // Get sphere scale from settings if available
        // Setting ID 155 is sphere_scale, returns f32 directly
        let sphere_scale = settings
            .get_float_if_defined(155)
            .unwrap_or(self.sphere_scale);

        // Iterate over atoms
        for (atom_idx, coord) in coord_set.iter_with_atoms() {
            let atom = match molecule.get_atom(atom_idx) {
                Some(a) => a,
                None => continue,
            };

            // Check if spheres representation is visible for this atom
            if !atom.repr.visible_reps.is_visible(RepMask::SPHERES) {
                continue;
            }

            // Get VdW radius - use per-atom sphere_scale if set, otherwise global
            let scale = atom.repr.sphere_scale.unwrap_or(sphere_scale);
            let radius = atom.effective_vdw() * scale;

            // Get color
            let color = colors.resolve_sphere(atom, molecule);

            // Add sphere instance
            self.instances.push(SphereVertex {
                center: [coord.x, coord.y, coord.z],
                radius,
                color,
            });
        }

        self.instance_count = self.instances.len() as u32;
        self.dirty = false;
    }

    fn upload(&mut self, device: &wgpu::Device, queue: &wgpu::Queue) {
        if !self.instances.is_empty() {
            self.instance_buffer.write(device, queue, &self.instances);
        }
    }

    fn render<'a>(&'a self, render_pass: &mut wgpu::RenderPass<'a>) {
        if self.instance_count == 0 {
            return;
        }

        // The caller is responsible for setting pipeline, bind_group, billboard vertex buffer,
        // and index buffer. We just need to set the instance buffer (slot 1) and draw.
        if let Some(instances) = self.instance_buffer.buffer() {
            render_pass.set_vertex_buffer(1, instances.slice(..));
            render_pass.draw_indexed(0..6, 0, 0..self.instance_count);
        }
    }

    fn is_dirty(&self) -> bool {
        self.dirty
    }

    fn set_dirty(&mut self) {
        self.dirty = true;
    }

    fn primitive_count(&self) -> usize {
        self.instance_count as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sphere_rep_new() {
        let rep = SphereRep::new();
        assert!(rep.is_dirty());
        assert!(rep.is_empty());
    }
}
