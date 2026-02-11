//! Stick representation for bond visualization
//!
//! Renders bonds as cylinders using impostor shaders, with sphere caps at atoms.
//! Supports multiple bond visualization (double, triple, aromatic) when valence display is enabled.

use crate::buffer::GrowableBuffer;
use crate::color_resolver::ColorResolver;
use crate::representation::Representation;
use crate::vertex::{CylinderVertex, SphereVertex};

use pymol_mol::{CoordSet, ObjectMolecule, RepMask};
use pymol_settings::SettingResolver;

use super::bond_utils::{
    calculate_perpendicular_with_neighbor, find_neighbor_position, get_bond_offsets, normalize,
};

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

        // Get valence display settings
        // Setting ID 64 = valence (bool), ID 512 = stick_valence_scale (float)
        let valence_enabled = settings.get_bool_if_defined(64).unwrap_or(true);
        let stick_valence_scale = settings.get_float_if_defined(512).unwrap_or(1.0);

        // For sticks, the offset needs to be based on stick_radius to prevent overlap.
        // Using stick_radius * 1.5 as base gives a compact appearance for double/triple bonds.
        // The offset factors (-0.5, 0.5 for double bonds) will be multiplied by this.
        let stick_valence_offset = stick_radius * 1.5 * stick_valence_scale;

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
            if !atom1.repr.visible_reps.is_visible(RepMask::STICKS)
                && !atom2.repr.visible_reps.is_visible(RepMask::STICKS)
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
            let color1 = colors.resolve_stick(atom1, molecule);
            let color2 = colors.resolve_stick(atom2, molecule);

            // Use smaller radius for multiple bonds (like PyMOL)
            let radius = if bond.order.is_multiple() && valence_enabled {
                stick_radius * 0.5
            } else {
                stick_radius
            };

            // Get offsets for multiple bonds (or single offset of 0.0 for single bonds)
            let offsets = if valence_enabled {
                get_bond_offsets(bond.order)
            } else {
                &[0.0_f32] as &[f32]
            };

            // Calculate perpendicular direction for offsetting multiple bonds
            let bond_dir = normalize([
                pos2.x - pos1.x,
                pos2.y - pos1.y,
                pos2.z - pos1.z,
            ]);

            // Find a neighbor atom to determine the molecular plane
            // This ensures double bonds stay in-plane (e.g., in aromatic rings)
            let neighbor_pos = find_neighbor_position(
                molecule,
                coord_set,
                atom1_idx,
                atom2_idx,
            );
            let perp = calculate_perpendicular_with_neighbor(
                bond_dir,
                [pos1.x, pos1.y, pos1.z],
                [pos2.x, pos2.y, pos2.z],
                neighbor_pos,
            );

            // Generate cylinder(s) for this bond, with sphere caps at each endpoint
            for &offset_factor in offsets {
                let offset = offset_factor * stick_valence_offset;
                let off = [perp[0] * offset, perp[1] * offset, perp[2] * offset];

                // Offset positions
                let p1 = [pos1.x + off[0], pos1.y + off[1], pos1.z + off[2]];
                let p2 = [pos2.x + off[0], pos2.y + off[1], pos2.z + off[2]];

                // Add cylinder instance (no built-in caps, we use spheres)
                self.cylinder_instances.push(CylinderVertex {
                    start: p1,
                    radius,
                    end: p2,
                    flags: 0, // No caps - we use sphere caps instead
                    color1,
                    color2,
                });

                // Add sphere caps at each cylinder endpoint
                self.sphere_instances.push(SphereVertex {
                    center: p1,
                    radius,
                    color: color1,
                });
                self.sphere_instances.push(SphereVertex {
                    center: p2,
                    radius,
                    color: color2,
                });
            }
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

    fn clear_cpu_data(&mut self) {
        self.cylinder_instances = Vec::new();
        self.sphere_instances = Vec::new();
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
