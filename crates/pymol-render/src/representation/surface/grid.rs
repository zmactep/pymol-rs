//! 3D grid structure for surface generation
//!
//! Provides spatial data structures for distance field computation and
//! marching cubes isosurface extraction.

/// 3D uniform grid for scalar field storage
#[derive(Debug, Clone)]
pub struct Grid3D {
    /// Grid origin (minimum corner in world coordinates)
    pub origin: [f32; 3],
    /// Grid cell spacing (uniform in all dimensions)
    pub spacing: f32,
    /// Grid dimensions (number of cells in x, y, z)
    pub dims: [usize; 3],
    /// Scalar field values at grid vertices
    /// Layout: values[x + y * (dims[0]+1) + z * (dims[0]+1) * (dims[1]+1)]
    pub values: Vec<f32>,
}

impl Grid3D {
    /// Create a new grid from bounding box and spacing
    ///
    /// The grid will encompass the bounding box with the given padding.
    pub fn from_bounds(min: [f32; 3], max: [f32; 3], spacing: f32, padding: f32) -> Self {
        // Add padding to bounds
        let padded_min = [
            min[0] - padding,
            min[1] - padding,
            min[2] - padding,
        ];
        let padded_max = [
            max[0] + padding,
            max[1] + padding,
            max[2] + padding,
        ];

        // Calculate grid dimensions
        let dims = [
            ((padded_max[0] - padded_min[0]) / spacing).ceil() as usize + 1,
            ((padded_max[1] - padded_min[1]) / spacing).ceil() as usize + 1,
            ((padded_max[2] - padded_min[2]) / spacing).ceil() as usize + 1,
        ];

        // Number of vertices = dims + 1 in each dimension
        let vertex_count = (dims[0] + 1) * (dims[1] + 1) * (dims[2] + 1);

        Self {
            origin: padded_min,
            spacing,
            dims,
            values: vec![f32::MAX; vertex_count],
        }
    }

    /// Get the number of vertices in each dimension
    #[inline]
    pub fn vertex_dims(&self) -> [usize; 3] {
        [self.dims[0] + 1, self.dims[1] + 1, self.dims[2] + 1]
    }

    /// Get the total number of vertices
    #[inline]
    #[allow(dead_code)]
    pub fn vertex_count(&self) -> usize {
        let vd = self.vertex_dims();
        vd[0] * vd[1] * vd[2]
    }

    /// Convert 3D vertex indices to linear index
    #[inline]
    pub fn vertex_index(&self, x: usize, y: usize, z: usize) -> usize {
        let vd = self.vertex_dims();
        x + y * vd[0] + z * vd[0] * vd[1]
    }

    /// Get world position of a vertex
    #[inline]
    pub fn vertex_position(&self, x: usize, y: usize, z: usize) -> [f32; 3] {
        [
            self.origin[0] + x as f32 * self.spacing,
            self.origin[1] + y as f32 * self.spacing,
            self.origin[2] + z as f32 * self.spacing,
        ]
    }

    /// Get the value at a vertex
    #[inline]
    pub fn get(&self, x: usize, y: usize, z: usize) -> f32 {
        self.values[self.vertex_index(x, y, z)]
    }

    /// Set the value at a vertex
    #[inline]
    #[allow(dead_code)]
    pub fn set(&mut self, x: usize, y: usize, z: usize, value: f32) {
        let idx = self.vertex_index(x, y, z);
        self.values[idx] = value;
    }

    /// Update the value at a vertex to the minimum of current and new value
    #[inline]
    #[allow(dead_code)]
    pub fn update_min(&mut self, x: usize, y: usize, z: usize, value: f32) {
        let idx = self.vertex_index(x, y, z);
        if value < self.values[idx] {
            self.values[idx] = value;
        }
    }

    /// Get the 8 corner values of a cube at the given cell position
    #[inline]
    pub fn cube_values(&self, x: usize, y: usize, z: usize) -> [f32; 8] {
        [
            self.get(x, y, z),         // 0: (0,0,0)
            self.get(x + 1, y, z),     // 1: (1,0,0)
            self.get(x + 1, y + 1, z), // 2: (1,1,0)
            self.get(x, y + 1, z),     // 3: (0,1,0)
            self.get(x, y, z + 1),     // 4: (0,0,1)
            self.get(x + 1, y, z + 1), // 5: (1,0,1)
            self.get(x + 1, y + 1, z + 1), // 6: (1,1,1)
            self.get(x, y + 1, z + 1), // 7: (0,1,1)
        ]
    }

    /// Get the 8 corner positions of a cube at the given cell position
    #[inline]
    pub fn cube_positions(&self, x: usize, y: usize, z: usize) -> [[f32; 3]; 8] {
        [
            self.vertex_position(x, y, z),         // 0
            self.vertex_position(x + 1, y, z),     // 1
            self.vertex_position(x + 1, y + 1, z), // 2
            self.vertex_position(x, y + 1, z),     // 3
            self.vertex_position(x, y, z + 1),     // 4
            self.vertex_position(x + 1, y, z + 1), // 5
            self.vertex_position(x + 1, y + 1, z + 1), // 6
            self.vertex_position(x, y + 1, z + 1), // 7
        ]
    }
}

/// Spatial hash grid for accelerated atom lookups
pub struct SpatialHash {
    /// Cell size for hashing
    cell_size: f32,
    /// Origin of the spatial hash grid
    origin: [f32; 3],
    /// Grid dimensions
    dims: [usize; 3],
    /// Cells containing atom indices
    cells: Vec<Vec<usize>>,
}

impl SpatialHash {
    /// Create a new spatial hash from atom positions
    ///
    /// The caller provides the bounds that the spatial hash should cover.
    /// The spatial hash will use these bounds directly (no additional padding).
    pub fn new(
        positions: &[[f32; 3]],
        radii: &[f32],
        min: [f32; 3],
        max: [f32; 3],
    ) -> Self {
        // Use a cell size that balances memory and query efficiency
        // Typical VdW radii are 1-2 Å, so 3 Å cells work well
        let cell_size = 3.0_f32;

        // Use the provided bounds directly - the caller is responsible for
        // ensuring the bounds cover the region of interest
        let origin = min;

        let dims = [
            ((max[0] - min[0]) / cell_size).ceil() as usize + 1,
            ((max[1] - min[1]) / cell_size).ceil() as usize + 1,
            ((max[2] - min[2]) / cell_size).ceil() as usize + 1,
        ];

        let cell_count = dims[0] * dims[1] * dims[2];
        let mut cells: Vec<Vec<usize>> = vec![Vec::new(); cell_count];

        // Insert atoms into cells, accounting for their radii
        for (atom_idx, (pos, radius)) in positions.iter().zip(radii.iter()).enumerate() {
            // Calculate the cell range this atom could influence
            let r = *radius;
            let min_cell = [
                ((pos[0] - r - origin[0]) / cell_size).floor() as i32,
                ((pos[1] - r - origin[1]) / cell_size).floor() as i32,
                ((pos[2] - r - origin[2]) / cell_size).floor() as i32,
            ];
            let max_cell = [
                ((pos[0] + r - origin[0]) / cell_size).ceil() as i32,
                ((pos[1] + r - origin[1]) / cell_size).ceil() as i32,
                ((pos[2] + r - origin[2]) / cell_size).ceil() as i32,
            ];

            // Insert into all cells the atom could influence
            for z in min_cell[2].max(0)..=max_cell[2].min(dims[2] as i32 - 1) {
                for y in min_cell[1].max(0)..=max_cell[1].min(dims[1] as i32 - 1) {
                    for x in min_cell[0].max(0)..=max_cell[0].min(dims[0] as i32 - 1) {
                        let cell_idx = x as usize + y as usize * dims[0] + z as usize * dims[0] * dims[1];
                        cells[cell_idx].push(atom_idx);
                    }
                }
            }
        }

        Self {
            cell_size,
            origin,
            dims,
            cells,
        }
    }

    /// Get atoms that could be near the given point (single cell - fast but may miss edge cases)
    #[allow(dead_code)]
    pub fn query(&self, point: [f32; 3]) -> &[usize] {
        let x = ((point[0] - self.origin[0]) / self.cell_size).floor() as i32;
        let y = ((point[1] - self.origin[1]) / self.cell_size).floor() as i32;
        let z = ((point[2] - self.origin[2]) / self.cell_size).floor() as i32;

        if x < 0 || y < 0 || z < 0 
            || x >= self.dims[0] as i32 
            || y >= self.dims[1] as i32 
            || z >= self.dims[2] as i32 
        {
            return &[];
        }

        let cell_idx = x as usize + y as usize * self.dims[0] + z as usize * self.dims[0] * self.dims[1];
        &self.cells[cell_idx]
    }

    /// Get atoms from a 3x3x3 neighborhood of cells around the given point
    /// 
    /// This is more robust than single-cell query, as it handles cases where
    /// atoms near cell boundaries might influence points in adjacent cells.
    /// 
    /// Optimized for performance with pre-computed strides and minimized bounds checks.
    pub fn query_neighborhood(&self, point: [f32; 3], result: &mut Vec<usize>) {
        result.clear();
        
        let inv_cell_size = 1.0 / self.cell_size;
        let cx = ((point[0] - self.origin[0]) * inv_cell_size).floor() as i32;
        let cy = ((point[1] - self.origin[1]) * inv_cell_size).floor() as i32;
        let cz = ((point[2] - self.origin[2]) * inv_cell_size).floor() as i32;

        // Pre-compute bounds for the neighborhood
        let x_min = (cx - 1).max(0) as usize;
        let x_max = ((cx + 1) as usize).min(self.dims[0] - 1);
        let y_min = (cy - 1).max(0) as usize;
        let y_max = ((cy + 1) as usize).min(self.dims[1] - 1);
        let z_min = (cz - 1).max(0) as usize;
        let z_max = ((cz + 1) as usize).min(self.dims[2] - 1);

        // Pre-compute strides
        let stride_y = self.dims[0];
        let stride_z = self.dims[0] * self.dims[1];

        // Iterate with pre-computed bounds (no inner-loop bounds checks)
        for z in z_min..=z_max {
            let z_offset = z * stride_z;
            for y in y_min..=y_max {
                let zy_offset = z_offset + y * stride_y;
                for x in x_min..=x_max {
                    let cell_idx = zy_offset + x;
                    let cell = &self.cells[cell_idx];
                    if !cell.is_empty() {
                        result.extend_from_slice(cell);
                    }
                }
            }
        }
    }

    /// Get cell index for a point (used for cache-friendly iteration patterns)
    #[inline]
    #[allow(dead_code)]
    pub fn cell_index(&self, point: [f32; 3]) -> Option<usize> {
        let inv_cell_size = 1.0 / self.cell_size;
        let x = ((point[0] - self.origin[0]) * inv_cell_size).floor() as i32;
        let y = ((point[1] - self.origin[1]) * inv_cell_size).floor() as i32;
        let z = ((point[2] - self.origin[2]) * inv_cell_size).floor() as i32;

        if x < 0 || y < 0 || z < 0 
            || x >= self.dims[0] as i32 
            || y >= self.dims[1] as i32 
            || z >= self.dims[2] as i32 
        {
            return None;
        }

        Some(x as usize + y as usize * self.dims[0] + z as usize * self.dims[0] * self.dims[1])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_grid_creation() {
        let grid = Grid3D::from_bounds(
            [0.0, 0.0, 0.0],
            [10.0, 10.0, 10.0],
            1.0,
            2.0,
        );

        // With padding of 2.0, bounds become [-2, -2, -2] to [12, 12, 12]
        // That's 14 units, at 1.0 spacing = 15 cells, 16 vertices per dimension
        assert_eq!(grid.dims[0], 15);
        assert_eq!(grid.dims[1], 15);
        assert_eq!(grid.dims[2], 15);
        assert_eq!(grid.vertex_dims(), [16, 16, 16]);
    }

    #[test]
    fn test_vertex_position() {
        let grid = Grid3D::from_bounds(
            [0.0, 0.0, 0.0],
            [10.0, 10.0, 10.0],
            0.5,
            1.0,
        );

        let pos = grid.vertex_position(0, 0, 0);
        assert!((pos[0] - (-1.0)).abs() < 0.001);
        assert!((pos[1] - (-1.0)).abs() < 0.001);
        assert!((pos[2] - (-1.0)).abs() < 0.001);

        let pos2 = grid.vertex_position(2, 2, 2);
        assert!((pos2[0] - 0.0).abs() < 0.001);
        assert!((pos2[1] - 0.0).abs() < 0.001);
        assert!((pos2[2] - 0.0).abs() < 0.001);
    }

    #[test]
    fn test_spatial_hash() {
        let positions = vec![
            [0.0, 0.0, 0.0],
            [5.0, 5.0, 5.0],
            [10.0, 10.0, 10.0],
        ];
        let radii = vec![1.5, 1.5, 1.5];

        let hash = SpatialHash::new(&positions, &radii, [0.0, 0.0, 0.0], [10.0, 10.0, 10.0]);

        // Query near first atom should return it
        let near_first = hash.query([0.5, 0.5, 0.5]);
        assert!(near_first.contains(&0));

        // Query near second atom
        let near_second = hash.query([5.0, 5.0, 5.0]);
        assert!(near_second.contains(&1));
    }
}
