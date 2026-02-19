//! Spatial hash grid for efficient neighbor queries
//!
//! Used by bond generation, H-bond detection, and other distance-based algorithms.

use lin_alg::f32::Vec3;
use std::collections::HashMap;

/// Simple spatial hash grid for efficient neighbor queries
///
/// Divides 3D space into uniform cubic cells. Each cell stores indices of
/// points that fall within it. Neighbor queries check the 3×3×3 neighborhood
/// of cells around the query point, guaranteeing all points within `cell_size`
/// distance are found.
pub(crate) struct SpatialGrid {
    cells: HashMap<(i32, i32, i32), Vec<usize>>,
    cell_size: f32,
}

impl SpatialGrid {
    pub fn new(cell_size: f32) -> Self {
        Self {
            cells: HashMap::new(),
            cell_size,
        }
    }

    fn cell_key(&self, pos: Vec3) -> (i32, i32, i32) {
        (
            (pos.x / self.cell_size).floor() as i32,
            (pos.y / self.cell_size).floor() as i32,
            (pos.z / self.cell_size).floor() as i32,
        )
    }

    pub fn insert(&mut self, pos: Vec3, idx: usize) {
        let key = self.cell_key(pos);
        self.cells.entry(key).or_default().push(idx);
    }

    /// Iterate over all indices in the 3×3×3 neighborhood of the given position
    pub fn query_neighbors(&self, pos: Vec3, out: &mut Vec<usize>) {
        out.clear();
        let (cx, cy, cz) = self.cell_key(pos);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    if let Some(indices) = self.cells.get(&(cx + dx, cy + dy, cz + dz)) {
                        out.extend_from_slice(indices);
                    }
                }
            }
        }
    }
}
