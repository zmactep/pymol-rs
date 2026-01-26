//! Stick representation for bond visualization
//!
//! Renders bonds as cylinders using impostor shaders, with sphere caps at atoms.

use crate::buffer::GrowableBuffer;
use crate::color_resolver::ColorResolver;
use crate::representation::Representation;
use crate::vertex::{CylinderVertex, SphereVertex};

use std::collections::HashMap;
use pymol_mol::{AtomIndex, CoordSet, ObjectMolecule, RepMask};
use pymol_settings::SettingResolver;

/// Stick representation for bonds
///
/// Renders bonds as cylinders using impostor shaders. Each cylinder is
/// rendered as an oriented billboard quad with ray-cylinder intersection
/// computed in the fragment shader. Sphere caps are rendered at atom positions
/// to give a smooth, continuous appearance.
pub struct StickRep {
    /// Cylinder instance data (CPU side)
    cylinder_instances: Vec<CylinderVertex>,
    /// Sphere cap instance data (CPU side)
    sphere_instances: Vec<SphereVertex>,
    /// GPU cylinder instance buffer
    cylinder_buffer: GrowableBuffer,
    /// GPU sphere instance buffer
    sphere_buffer: GrowableBuffer,
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
    /// Number of cylinder instances to render
    cylinder_count: u32,
    /// Number of sphere cap instances to render
    sphere_count: u32,
    /// Default stick radius
    stick_radius: f32,
}

impl StickRep {
    /// Create a new stick representation
    pub fn new() -> Self {
        Self {
            cylinder_instances: Vec::new(),
            sphere_instances: Vec::new(),
            cylinder_buffer: GrowableBuffer::new("Stick Cylinder Instances", wgpu::BufferUsages::VERTEX),
            sphere_buffer: GrowableBuffer::new("Stick Sphere Caps", wgpu::BufferUsages::VERTEX),
            pipeline: None,
            billboard_buffer: None,
            index_buffer: None,
            bind_group: None,
            dirty: true,
            cylinder_count: 0,
            sphere_count: 0,
            stick_radius: 0.25,
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

    /// Set the default stick radius
    pub fn set_radius(&mut self, radius: f32) {
        self.stick_radius = radius;
    }

    /// Get the sphere instance buffer for rendering caps
    pub fn sphere_buffer(&self) -> Option<&wgpu::Buffer> {
        self.sphere_buffer.buffer()
    }

    /// Get the number of sphere cap instances
    pub fn sphere_count(&self) -> u32 {
        self.sphere_count
    }
}

impl Default for StickRep {
    fn default() -> Self {
        Self::new()
    }
}

impl Representation for StickRep {
    fn build(
        &mut self,
        molecule: &ObjectMolecule,
        coord_set: &CoordSet,
        colors: &ColorResolver,
        settings: &SettingResolver,
    ) {
        self.cylinder_instances.clear();
        self.sphere_instances.clear();

        // Get stick radius from settings if available
        // Setting ID 21 is stick_radius, returns f32 directly
        let stick_radius = settings
            .get_float_if_defined(21)
            .unwrap_or(self.stick_radius);

        // Track atoms with their maximum connecting stick radius
        // The sphere cap needs to be large enough to cover all connecting cylinders
        let mut atom_max_radius: HashMap<AtomIndex, f32> = HashMap::new();

        // Iterate over bonds to create cylinders
        for bond in molecule.bonds() {
            let atom1_idx = bond.atom1;
            let atom2_idx = bond.atom2;

            // Get atom data
            let atom1 = match molecule.get_atom(atom1_idx) {
                Some(a) => a,
                None => continue,
            };
            let atom2 = match molecule.get_atom(atom2_idx) {
                Some(a) => a,
                None => continue,
            };

            // Check if sticks representation is visible for either atom
            if !atom1.visible_reps.is_visible(RepMask::STICKS)
                && !atom2.visible_reps.is_visible(RepMask::STICKS)
            {
                continue;
            }

            // Get coordinates
            let pos1 = match coord_set.get_atom_coord(atom1_idx) {
                Some(p) => p,
                None => continue,
            };
            let pos2 = match coord_set.get_atom_coord(atom2_idx) {
                Some(p) => p,
                None => continue,
            };

            // Get colors
            let color1 = colors.resolve_atom(atom1, molecule);
            let color2 = colors.resolve_atom(atom2, molecule);

            // Adjust radius based on bond order
            let radius = stick_radius * bond.order.as_float().sqrt();

            // Add cylinder instance
            self.cylinder_instances.push(CylinderVertex {
                start: [pos1.x, pos1.y, pos1.z],
                radius,
                end: [pos2.x, pos2.y, pos2.z],
                flags: CylinderVertex::CAP_BOTH,
                color1,
                color2,
            });

            // Track maximum radius for each atom's sphere cap
            atom_max_radius
                .entry(atom1_idx)
                .and_modify(|r| *r = r.max(radius))
                .or_insert(radius);
            atom_max_radius
                .entry(atom2_idx)
                .and_modify(|r| *r = r.max(radius))
                .or_insert(radius);
        }

        // Create sphere caps at atom positions with the maximum radius
        for (atom_idx, max_radius) in atom_max_radius {
            let atom = match molecule.get_atom(atom_idx) {
                Some(a) => a,
                None => continue,
            };

            let coord = match coord_set.get_atom_coord(atom_idx) {
                Some(c) => c,
                None => continue,
            };

            let color = colors.resolve_atom(atom, molecule);

            // Sphere radius is the maximum of all connecting stick radii
            self.sphere_instances.push(SphereVertex {
                center: [coord.x, coord.y, coord.z],
                radius: max_radius,
                color,
            });
        }

        self.cylinder_count = self.cylinder_instances.len() as u32;
        self.sphere_count = self.sphere_instances.len() as u32;
        self.dirty = false;
    }

    fn upload(&mut self, device: &wgpu::Device, queue: &wgpu::Queue) {
        if !self.cylinder_instances.is_empty() {
            self.cylinder_buffer.write(device, queue, &self.cylinder_instances);
        }
        if !self.sphere_instances.is_empty() {
            self.sphere_buffer.write(device, queue, &self.sphere_instances);
        }
    }

    fn render<'a>(&'a self, render_pass: &mut wgpu::RenderPass<'a>) {
        if self.cylinder_count == 0 {
            return;
        }

        // The caller is responsible for setting pipeline, bind_group, billboard vertex buffer,
        // and index buffer. We just need to set the instance buffer (slot 1) and draw.
        // Note: Sphere caps are rendered separately by the caller using sphere_buffer()
        if let Some(instances) = self.cylinder_buffer.buffer() {
            render_pass.set_vertex_buffer(1, instances.slice(..));
            render_pass.draw_indexed(0..6, 0, 0..self.cylinder_count);
        }
    }

    fn is_dirty(&self) -> bool {
        self.dirty
    }

    fn set_dirty(&mut self) {
        self.dirty = true;
    }

    fn primitive_count(&self) -> usize {
        (self.cylinder_count + self.sphere_count) as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stick_rep_new() {
        let rep = StickRep::new();
        assert!(rep.is_dirty());
        assert!(rep.is_empty());
    }
}
