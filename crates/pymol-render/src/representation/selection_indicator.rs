//! Selection indicator representation
//!
//! Renders pink/magenta indicators at selected atom positions,
//! similar to PyMOL's selection visualization.

use crate::buffer::GrowableBuffer;
use crate::vertex::DotVertex;

use pymol_mol::{CoordSet, ObjectMolecule};
use pymol_select::SelectionResult;

/// Default selection indicator color (pink/magenta like PyMOL)
pub const SELECTION_INDICATOR_COLOR: [f32; 4] = [1.0, 0.467, 1.0, 1.0]; // RGB(255, 119, 255)

/// Default selection indicator size (larger for visibility)
pub const DEFAULT_INDICATOR_SIZE: f32 = 8.0;

/// Selection indicator representation
///
/// Renders small pink/magenta dots at selected atom positions.
/// This provides visual feedback for the current selection in the 3D viewport.
pub struct SelectionIndicatorRep {
    /// Instance data (CPU side)
    instances: Vec<DotVertex>,
    /// GPU instance buffer
    instance_buffer: GrowableBuffer,
    /// Number of indicator instances to render
    instance_count: u32,
    /// Indicator dot size
    indicator_size: f32,
    /// Indicator color
    indicator_color: [f32; 4],
    /// Whether the representation is dirty and needs rebuild
    dirty: bool,
}

impl SelectionIndicatorRep {
    /// Create a new selection indicator representation
    pub fn new() -> Self {
        Self {
            instances: Vec::new(),
            instance_buffer: GrowableBuffer::new(
                "Selection Indicator Instances",
                wgpu::BufferUsages::VERTEX,
            ),
            instance_count: 0,
            indicator_size: DEFAULT_INDICATOR_SIZE,
            indicator_color: SELECTION_INDICATOR_COLOR,
            dirty: true,
        }
    }

    /// Set the indicator dot size
    pub fn set_size(&mut self, size: f32) {
        if (self.indicator_size - size).abs() > f32::EPSILON {
            self.indicator_size = size;
            self.dirty = true;
        }
    }

    /// Get the current indicator size
    pub fn size(&self) -> f32 {
        self.indicator_size
    }

    /// Set the indicator color
    pub fn set_color(&mut self, color: [f32; 4]) {
        self.indicator_color = color;
        self.dirty = true;
    }

    /// Get the current indicator color
    pub fn color(&self) -> [f32; 4] {
        self.indicator_color
    }

    /// Build the indicator geometry for a selection
    ///
    /// This generates dot instances at each selected atom's position.
    ///
    /// # Arguments
    /// * `molecule` - The molecule containing the atoms
    /// * `coord_set` - The coordinate set to use for positions
    /// * `selection` - The selection result indicating which atoms to show
    pub fn build_for_selection(
        &mut self,
        molecule: &ObjectMolecule,
        coord_set: &CoordSet,
        selection: &SelectionResult,
    ) {
        self.instances.clear();

        // Iterate over selected atom indices
        for atom_idx in selection.indices() {
            // Get the coordinate for this atom
            if let Some(coord) = coord_set.get_atom_coord(atom_idx) {
                // Verify atom exists
                if molecule.get_atom(atom_idx).is_some() {
                    self.instances.push(DotVertex {
                        position: [coord.x, coord.y, coord.z],
                        size: self.indicator_size,
                        color: self.indicator_color,
                    });
                }
            }
        }

        self.instance_count = self.instances.len() as u32;
        self.dirty = false;
    }

    /// Clear all indicator instances
    pub fn clear(&mut self) {
        self.instances.clear();
        self.instance_count = 0;
        self.dirty = false;
    }

    /// Upload vertex data to GPU buffers
    pub fn upload(&mut self, device: &wgpu::Device, queue: &wgpu::Queue) {
        if !self.instances.is_empty() {
            self.instance_buffer.write(device, queue, &self.instances);
        }
    }

    /// Render the selection indicators
    ///
    /// Note: The caller must set up the pipeline, bind groups, billboard vertex buffer,
    /// and index buffer before calling this. This method only sets the instance buffer
    /// and issues the draw call.
    pub fn render<'a>(&'a self, render_pass: &mut wgpu::RenderPass<'a>) {
        if self.instance_count == 0 {
            return;
        }

        if let Some(instances) = self.instance_buffer.buffer() {
            render_pass.set_vertex_buffer(1, instances.slice(..));
            render_pass.draw_indexed(0..6, 0, 0..self.instance_count);
        }
    }

    /// Check if the representation is dirty
    pub fn is_dirty(&self) -> bool {
        self.dirty
    }

    /// Mark the representation as dirty
    pub fn set_dirty(&mut self) {
        self.dirty = true;
    }

    /// Get the number of indicator instances
    pub fn instance_count(&self) -> u32 {
        self.instance_count
    }

    /// Check if there are no indicators to render
    pub fn is_empty(&self) -> bool {
        self.instance_count == 0
    }
}

impl Default for SelectionIndicatorRep {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_selection_indicator_new() {
        let rep = SelectionIndicatorRep::new();
        assert!(rep.is_dirty());
        assert!(rep.is_empty());
        assert_eq!(rep.size(), DEFAULT_INDICATOR_SIZE);
        assert_eq!(rep.color(), SELECTION_INDICATOR_COLOR);
    }

    #[test]
    fn test_selection_indicator_set_size() {
        let mut rep = SelectionIndicatorRep::new();
        rep.set_size(8.0);
        assert_eq!(rep.size(), 8.0);
    }

    #[test]
    fn test_selection_indicator_set_color() {
        let mut rep = SelectionIndicatorRep::new();
        let new_color = [1.0, 0.0, 0.0, 1.0];
        rep.set_color(new_color);
        assert_eq!(rep.color(), new_color);
    }

    #[test]
    fn test_selection_indicator_clear() {
        let mut rep = SelectionIndicatorRep::new();
        rep.clear();
        assert!(rep.is_empty());
        assert!(!rep.is_dirty());
    }
}
