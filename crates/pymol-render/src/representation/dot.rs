//! Dot representation for dot surface visualization
//!
//! Renders atoms as small dots on their VdW surface.

use crate::buffer::GrowableBuffer;
use crate::color_resolver::ColorResolver;
use crate::representation::Representation;
use crate::vertex::DotVertex;

use pymol_mol::{CoordSet, ObjectMolecule, RepMask};
use pymol_settings::SettingResolver;

/// Dot representation for atoms
///
/// Renders atoms as small dots, typically distributed on the VdW surface.
/// This is a simplified implementation that renders a single dot per atom.
pub struct DotRep {
    /// Instance data (CPU side)
    instances: Vec<DotVertex>,
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
    /// Number of dot instances to render
    instance_count: u32,
    /// Dot size
    dot_width: f32,
    /// Dot density (points per Angstrom^2)
    dot_density: i32,
}

impl DotRep {
    /// Create a new dot representation
    pub fn new() -> Self {
        Self {
            instances: Vec::new(),
            instance_buffer: GrowableBuffer::new("Dot Instances", wgpu::BufferUsages::VERTEX),
            pipeline: None,
            billboard_buffer: None,
            index_buffer: None,
            bind_group: None,
            dirty: true,
            instance_count: 0,
            dot_width: 2.0,
            dot_density: 2,
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

    /// Set dot width (size)
    pub fn set_width(&mut self, width: f32) {
        self.dot_width = width;
    }

    /// Set dot density
    pub fn set_density(&mut self, density: i32) {
        self.dot_density = density;
    }

    /// Generate points on a sphere surface using golden spiral
    fn generate_sphere_points(n: usize, radius: f32, center: [f32; 3]) -> Vec<[f32; 3]> {
        let mut points = Vec::with_capacity(n);
        let phi = std::f32::consts::PI * (3.0 - 5.0_f32.sqrt()); // Golden angle

        for i in 0..n {
            let y = 1.0 - (i as f32 / (n - 1) as f32) * 2.0; // y goes from 1 to -1
            let r_at_y = (1.0 - y * y).sqrt();
            let theta = phi * i as f32;

            let x = theta.cos() * r_at_y;
            let z = theta.sin() * r_at_y;

            points.push([
                center[0] + x * radius,
                center[1] + y * radius,
                center[2] + z * radius,
            ]);
        }

        points
    }
}

impl Default for DotRep {
    fn default() -> Self {
        Self::new()
    }
}

impl Representation for DotRep {
    fn build(
        &mut self,
        molecule: &ObjectMolecule,
        coord_set: &CoordSet,
        colors: &ColorResolver,
        settings: &SettingResolver,
    ) {
        self.instances.clear();

        // Get settings
        // Setting ID 77 is dot_width, Setting ID 2 is dot_density
        let dot_width = settings
            .get_float_if_defined(77)
            .unwrap_or(self.dot_width);
        let dot_density = settings
            .get_int_if_defined(2)
            .unwrap_or(self.dot_density);

        // Calculate number of dots based on density
        // density 1 = ~10 dots, 2 = ~30 dots, 3 = ~100 dots, 4 = ~300 dots
        let base_dots = match dot_density {
            0 => 5,
            1 => 10,
            2 => 30,
            3 => 100,
            _ => 300,
        };

        // Iterate over atoms
        for (atom_idx, coord) in coord_set.iter_with_atoms() {
            let atom = match molecule.get_atom(atom_idx) {
                Some(a) => a,
                None => continue,
            };

            // Check if dots representation is visible for this atom
            if !atom.repr.visible_reps.is_visible(RepMask::DOTS) {
                continue;
            }

            // Get VdW radius
            let radius = atom.effective_vdw();

            // Get color
            let color = colors.resolve_atom(atom, molecule);

            // Generate dots on surface
            let center = [coord.x, coord.y, coord.z];
            let points = Self::generate_sphere_points(base_dots, radius, center);

            for point in points {
                self.instances.push(DotVertex {
                    position: point,
                    size: dot_width,
                    color,
                });
            }
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
    fn test_dot_rep_new() {
        let rep = DotRep::new();
        assert!(rep.is_dirty());
        assert!(rep.is_empty());
    }

    #[test]
    fn test_sphere_points() {
        let points = DotRep::generate_sphere_points(10, 1.0, [0.0, 0.0, 0.0]);
        assert_eq!(points.len(), 10);

        // All points should be approximately on unit sphere
        for p in &points {
            let dist = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            assert!((dist - 1.0).abs() < 0.01);
        }
    }
}
