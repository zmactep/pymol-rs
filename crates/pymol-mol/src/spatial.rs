//! Spatial hash grid for efficient neighbor queries
//!
//! Used by bond generation, H-bond detection, and other distance-based algorithms.

use ahash::AHashMap;
use lin_alg::f32::Vec3;

/// Simple spatial hash grid for efficient neighbor queries
///
/// Divides 3D space into uniform cubic cells. Each cell stores indices of
/// points that fall within it. Neighbor queries check the 3×3×3 neighborhood
/// of cells around the query point, guaranteeing all points within `cell_size`
/// distance are found.
pub(crate) struct SpatialGrid {
    cells: AHashMap<(i32, i32, i32), Vec<usize>>,
    cell_size: f32,
}

impl SpatialGrid {
    pub fn with_capacity(cell_size: f32, expected_atoms: usize) -> Self {
        Self {
            cells: AHashMap::with_capacity(expected_atoms),
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_insert_and_query() {
        let mut grid = SpatialGrid::with_capacity(2.0, 10);
        grid.insert(Vec3::new(0.0, 0.0, 0.0), 0);
        grid.insert(Vec3::new(1.0, 0.0, 0.0), 1);
        grid.insert(Vec3::new(0.0, 1.0, 0.0), 2);

        let mut neighbors = Vec::new();
        grid.query_neighbors(Vec3::new(0.5, 0.5, 0.0), &mut neighbors);
        neighbors.sort();

        assert_eq!(neighbors, vec![0, 1, 2]);
    }

    #[test]
    fn test_query_no_neighbors() {
        let mut grid = SpatialGrid::with_capacity(1.0, 10);
        grid.insert(Vec3::new(0.0, 0.0, 0.0), 0);

        let mut neighbors = Vec::new();
        // Query far away — outside 3x3x3 neighborhood
        grid.query_neighbors(Vec3::new(100.0, 100.0, 100.0), &mut neighbors);

        assert!(neighbors.is_empty());
    }

    #[test]
    fn test_all_atoms_same_cell() {
        let mut grid = SpatialGrid::with_capacity(10.0, 5);
        for i in 0..5 {
            grid.insert(Vec3::new(0.1 * i as f32, 0.0, 0.0), i);
        }

        let mut neighbors = Vec::new();
        grid.query_neighbors(Vec3::new(0.0, 0.0, 0.0), &mut neighbors);
        neighbors.sort();

        assert_eq!(neighbors, vec![0, 1, 2, 3, 4]);
    }

    #[test]
    fn test_negative_coordinates() {
        let mut grid = SpatialGrid::with_capacity(2.0, 4);
        grid.insert(Vec3::new(-1.0, -1.0, -1.0), 0);
        grid.insert(Vec3::new(1.0, 1.0, 1.0), 1);

        let mut neighbors = Vec::new();
        grid.query_neighbors(Vec3::new(0.0, 0.0, 0.0), &mut neighbors);
        neighbors.sort();

        assert_eq!(neighbors, vec![0, 1]);
    }
}
