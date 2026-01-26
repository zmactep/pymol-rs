//! Triangulation of surface dots.
//!
//! Provides two triangulation approaches:
//! 1. Fast per-atom triangulation using sphere connectivity (default)
//! 2. Global edge-based triangulation for merging atom surfaces (slower)

use super::dot_surface::SurfaceDot;
use std::collections::{HashMap, HashSet};

/// Result of triangulation
pub struct TriangulationResult {
    /// Triangle indices (groups of 3)
    pub indices: Vec<u32>,
    /// Number of triangles created (for debugging/statistics)
    #[allow(dead_code)]
    pub num_triangles: usize,
}

/// Spatial hash for dot lookups during triangulation
struct DotSpatialHash {
    cell_size: f32,
    origin: [f32; 3],
    dims: [usize; 3],
    cells: Vec<Vec<u32>>,
}

impl DotSpatialHash {
    fn new(dots: &[SurfaceDot], cutoff: f32) -> Self {
        if dots.is_empty() {
            return Self {
                cell_size: cutoff,
                origin: [0.0; 3],
                dims: [1, 1, 1],
                cells: vec![Vec::new()],
            };
        }

        // Find bounds
        let mut min = [f32::MAX; 3];
        let mut max = [f32::MIN; 3];
        for dot in dots {
            for i in 0..3 {
                min[i] = min[i].min(dot.position[i]);
                max[i] = max[i].max(dot.position[i]);
            }
        }

        // Add padding
        for i in 0..3 {
            min[i] -= cutoff;
            max[i] += cutoff;
        }

        let cell_size = cutoff;
        let dims = [
            ((max[0] - min[0]) / cell_size).ceil() as usize + 1,
            ((max[1] - min[1]) / cell_size).ceil() as usize + 1,
            ((max[2] - min[2]) / cell_size).ceil() as usize + 1,
        ];

        let cell_count = dims[0] * dims[1] * dims[2];
        let mut cells: Vec<Vec<u32>> = vec![Vec::new(); cell_count];

        // Insert dots into cells
        for (idx, dot) in dots.iter().enumerate() {
            let cx = ((dot.position[0] - min[0]) / cell_size).floor() as usize;
            let cy = ((dot.position[1] - min[1]) / cell_size).floor() as usize;
            let cz = ((dot.position[2] - min[2]) / cell_size).floor() as usize;

            let cell_idx = cx.min(dims[0] - 1)
                + cy.min(dims[1] - 1) * dims[0]
                + cz.min(dims[2] - 1) * dims[0] * dims[1];
            cells[cell_idx].push(idx as u32);
        }

        Self {
            cell_size,
            origin: min,
            dims,
            cells,
        }
    }

    /// Query neighbors within the cutoff distance
    fn query_neighbors(&self, pos: &[f32; 3], result: &mut Vec<u32>) {
        result.clear();

        let cx = ((pos[0] - self.origin[0]) / self.cell_size).floor() as i32;
        let cy = ((pos[1] - self.origin[1]) / self.cell_size).floor() as i32;
        let cz = ((pos[2] - self.origin[2]) / self.cell_size).floor() as i32;

        // Check 3x3x3 neighborhood
        for dz in -1..=1 {
            let z = cz + dz;
            if z < 0 || z >= self.dims[2] as i32 {
                continue;
            }
            for dy in -1..=1 {
                let y = cy + dy;
                if y < 0 || y >= self.dims[1] as i32 {
                    continue;
                }
                for dx in -1..=1 {
                    let x = cx + dx;
                    if x < 0 || x >= self.dims[0] as i32 {
                        continue;
                    }

                    let cell_idx = x as usize
                        + y as usize * self.dims[0]
                        + z as usize * self.dims[0] * self.dims[1];
                    result.extend_from_slice(&self.cells[cell_idx]);
                }
            }
        }
    }
}

/// Triangulation state tracking edges and triangles
struct TriangulationState {
    /// Edge -> (triangle count, vertex on one side, vertex on other side)
    /// Key is ordered pair (min, max) of vertex indices
    edge_info: HashMap<(u32, u32), EdgeInfo>,
    /// Active edges (edges with exactly 1 triangle)
    active_edges: Vec<(u32, u32)>,
    /// Output triangles
    triangles: Vec<[u32; 3]>,
    /// Vertices that are part of active edges
    active_vertices: HashSet<u32>,
}

#[derive(Clone, Copy, Default)]
struct EdgeInfo {
    /// Number of triangles sharing this edge (1 = active, 2 = closed)
    tri_count: u8,
    /// Third vertex of first triangle (for orientation check)
    third_vertex: u32,
}

impl TriangulationState {
    fn new(estimated_triangles: usize) -> Self {
        Self {
            edge_info: HashMap::with_capacity(estimated_triangles * 3),
            active_edges: Vec::with_capacity(estimated_triangles),
            triangles: Vec::with_capacity(estimated_triangles),
            active_vertices: HashSet::with_capacity(estimated_triangles),
        }
    }

    /// Create ordered edge key
    #[inline]
    fn edge_key(v1: u32, v2: u32) -> (u32, u32) {
        if v1 < v2 {
            (v1, v2)
        } else {
            (v2, v1)
        }
    }

    /// Add a triangle, updating edge status
    fn add_triangle(&mut self, v0: u32, v1: u32, v2: u32) {
        self.triangles.push([v0, v1, v2]);

        // Update edge info for each edge
        self.update_edge(v0, v1, v2);
        self.update_edge(v1, v2, v0);
        self.update_edge(v2, v0, v1);
    }

    /// Update edge info when adding a triangle
    fn update_edge(&mut self, v1: u32, v2: u32, third: u32) {
        let key = Self::edge_key(v1, v2);
        let info = self.edge_info.entry(key).or_default();
        info.tri_count += 1;

        if info.tri_count == 1 {
            // New active edge
            info.third_vertex = third;
            self.active_edges.push(key);
            self.active_vertices.insert(v1);
            self.active_vertices.insert(v2);
        } else if info.tri_count == 2 {
            // Edge is now closed - remove from active edges
            if let Some(pos) = self.active_edges.iter().position(|&e| e == key) {
                self.active_edges.swap_remove(pos);
            }
        }
    }

    /// Check if an edge can accept another triangle
    #[inline]
    fn edge_can_add_triangle(&self, v1: u32, v2: u32) -> bool {
        let key = Self::edge_key(v1, v2);
        match self.edge_info.get(&key) {
            None => true,
            Some(info) => info.tri_count < 2,
        }
    }

    /// Get the third vertex of the existing triangle on an edge (if any)
    fn get_edge_third_vertex(&self, v1: u32, v2: u32) -> Option<u32> {
        let key = Self::edge_key(v1, v2);
        self.edge_info.get(&key).map(|info| info.third_vertex)
    }
}

/// Squared distance between two points
#[inline]
fn dist_sq(a: &[f32; 3], b: &[f32; 3]) -> f32 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}

/// Check if a triangle has consistent normal orientation with vertex normals
fn triangle_orientation_ok(
    dots: &[SurfaceDot],
    v0: u32,
    v1: u32,
    v2: u32,
) -> bool {
    let p0 = &dots[v0 as usize].position;
    let p1 = &dots[v1 as usize].position;
    let p2 = &dots[v2 as usize].position;

    let n0 = &dots[v0 as usize].normal;
    let n1 = &dots[v1 as usize].normal;
    let n2 = &dots[v2 as usize].normal;

    // Compute triangle normal
    let edge1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
    let edge2 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];

    let tri_normal = [
        edge1[1] * edge2[2] - edge1[2] * edge2[1],
        edge1[2] * edge2[0] - edge1[0] * edge2[2],
        edge1[0] * edge2[1] - edge1[1] * edge2[0],
    ];

    // Average vertex normal
    let avg_normal = [
        n0[0] + n1[0] + n2[0],
        n0[1] + n1[1] + n2[1],
        n0[2] + n1[2] + n2[2],
    ];

    // Dot product should be positive (normals pointing same direction)
    let dot = tri_normal[0] * avg_normal[0]
        + tri_normal[1] * avg_normal[1]
        + tri_normal[2] * avg_normal[2];

    dot > 0.0
}

/// Check if all vertex normals point to the same side of the triangle
fn normals_consistent(dots: &[SurfaceDot], v0: u32, v1: u32, v2: u32) -> bool {
    let p0 = &dots[v0 as usize].position;
    let p1 = &dots[v1 as usize].position;
    let p2 = &dots[v2 as usize].position;

    let n0 = &dots[v0 as usize].normal;
    let n1 = &dots[v1 as usize].normal;
    let n2 = &dots[v2 as usize].normal;

    // Compute triangle normal
    let edge1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
    let edge2 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];

    let tri_normal = [
        edge1[1] * edge2[2] - edge1[2] * edge2[1],
        edge1[2] * edge2[0] - edge1[0] * edge2[2],
        edge1[0] * edge2[1] - edge1[1] * edge2[0],
    ];

    // Check each vertex normal against triangle normal
    let dot0 = tri_normal[0] * n0[0] + tri_normal[1] * n0[1] + tri_normal[2] * n0[2];
    let dot1 = tri_normal[0] * n1[0] + tri_normal[1] * n1[1] + tri_normal[2] * n1[2];
    let dot2 = tri_normal[0] * n2[0] + tri_normal[1] * n2[1] + tri_normal[2] * n2[2];

    // All should have the same sign
    (dot0 > 0.0 && dot1 > 0.0 && dot2 > 0.0) || (dot0 < 0.0 && dot1 < 0.0 && dot2 < 0.0)
}

/// Find seed triangles to start the triangulation
fn find_seed_triangles(
    dots: &[SurfaceDot],
    spatial_hash: &DotSpatialHash,
    cutoff_sq: f32,
    state: &mut TriangulationState,
) {
    let mut neighbors = Vec::with_capacity(64);
    let mut used = vec![false; dots.len()];

    for i in 0..dots.len() {
        if used[i] {
            continue;
        }

        let pi = &dots[i].position;
        spatial_hash.query_neighbors(pi, &mut neighbors);

        // Find two other vertices to form a seed triangle
        for (idx_j, &j) in neighbors.iter().enumerate() {
            if j as usize <= i || used[j as usize] {
                continue;
            }

            let pj = &dots[j as usize].position;
            let d_ij = dist_sq(pi, pj);
            if d_ij > cutoff_sq {
                continue;
            }

            for &k in neighbors.iter().skip(idx_j + 1) {
                if k as usize <= j as usize || used[k as usize] {
                    continue;
                }

                let pk = &dots[k as usize].position;
                let d_ik = dist_sq(pi, pk);
                let d_jk = dist_sq(pj, pk);

                if d_ik > cutoff_sq || d_jk > cutoff_sq {
                    continue;
                }

                // Check orientation
                if !normals_consistent(dots, i as u32, j, k) {
                    continue;
                }

                // Determine correct winding order
                if triangle_orientation_ok(dots, i as u32, j, k) {
                    state.add_triangle(i as u32, j, k);
                } else {
                    state.add_triangle(i as u32, k, j);
                }

                used[i] = true;
                used[j as usize] = true;
                used[k as usize] = true;
                break;
            }

            if used[i] {
                break;
            }
        }
    }
}

/// Expand the mesh by adding triangles to active edges
fn expand_mesh(
    dots: &[SurfaceDot],
    spatial_hash: &DotSpatialHash,
    cutoff_sq: f32,
    state: &mut TriangulationState,
    max_iterations: usize,
) {
    let mut neighbors = Vec::with_capacity(64);
    let mut iteration = 0;

    while !state.active_edges.is_empty() && iteration < max_iterations {
        iteration += 1;

        // Take the last active edge
        let (v1, v2) = match state.active_edges.pop() {
            Some(edge) => edge,
            None => break,
        };

        // Check if edge is still active (might have been closed by another triangle)
        if !state.edge_can_add_triangle(v1, v2) {
            continue;
        }

        let p1 = &dots[v1 as usize].position;
        let p2 = &dots[v2 as usize].position;

        // Get the third vertex of the existing triangle (to avoid duplicating it)
        let existing_third = state.get_edge_third_vertex(v1, v2);

        // Find midpoint of edge for neighbor query
        let mid = [
            (p1[0] + p2[0]) * 0.5,
            (p1[1] + p2[1]) * 0.5,
            (p1[2] + p2[2]) * 0.5,
        ];
        spatial_hash.query_neighbors(&mid, &mut neighbors);

        // Also query from both endpoints to ensure coverage
        let mut more_neighbors = Vec::new();
        spatial_hash.query_neighbors(p1, &mut more_neighbors);
        neighbors.extend(more_neighbors.iter());
        spatial_hash.query_neighbors(p2, &mut more_neighbors);
        neighbors.extend(more_neighbors.iter());

        // Remove duplicates
        neighbors.sort_unstable();
        neighbors.dedup();

        // Find the best third vertex
        let mut best_vertex: Option<u32> = None;
        let mut best_score = f32::MAX;

        for &candidate in &neighbors {
            // Skip endpoints of the edge
            if candidate == v1 || candidate == v2 {
                continue;
            }

            // Skip the existing third vertex
            if Some(candidate) == existing_third {
                continue;
            }

            let pc = &dots[candidate as usize].position;

            // Check edge lengths
            let d1 = dist_sq(p1, pc);
            let d2 = dist_sq(p2, pc);
            if d1 > cutoff_sq || d2 > cutoff_sq {
                continue;
            }

            // Check that the other two edges can accept triangles
            if !state.edge_can_add_triangle(v1, candidate)
                || !state.edge_can_add_triangle(v2, candidate)
            {
                continue;
            }

            // Check normal consistency
            if !normals_consistent(dots, v1, v2, candidate) {
                continue;
            }

            // Score: prefer smaller triangles (sum of edge lengths)
            let score = d1 + d2;
            if score < best_score {
                best_score = score;
                best_vertex = Some(candidate);
            }
        }

        // Add the best triangle
        if let Some(v3) = best_vertex {
            // Determine correct winding order based on existing triangle
            // The new triangle should be on the opposite side of the edge
            if let Some(existing) = existing_third {
                // Check which side the new vertex is on compared to existing
                let p_existing = &dots[existing as usize].position;
                
                // Cross product of edge with vector to existing third
                let edge = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]];
                let to_existing = [
                    p_existing[0] - p1[0],
                    p_existing[1] - p1[1],
                    p_existing[2] - p1[2],
                ];
                let cross_existing = [
                    edge[1] * to_existing[2] - edge[2] * to_existing[1],
                    edge[2] * to_existing[0] - edge[0] * to_existing[2],
                    edge[0] * to_existing[1] - edge[1] * to_existing[0],
                ];

                let p3 = &dots[v3 as usize].position;
                let to_new = [p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2]];
                let cross_new = [
                    edge[1] * to_new[2] - edge[2] * to_new[1],
                    edge[2] * to_new[0] - edge[0] * to_new[2],
                    edge[0] * to_new[1] - edge[1] * to_new[0],
                ];

                // Dot product of cross products
                let dot = cross_existing[0] * cross_new[0]
                    + cross_existing[1] * cross_new[1]
                    + cross_existing[2] * cross_new[2];

                if dot < 0.0 {
                    // Opposite side - use reversed winding
                    state.add_triangle(v1, v3, v2);
                } else {
                    // Same side (shouldn't happen in normal manifold)
                    state.add_triangle(v1, v2, v3);
                }
            } else {
                // No existing triangle - use orientation check
                if triangle_orientation_ok(dots, v1, v2, v3) {
                    state.add_triangle(v1, v2, v3);
                } else {
                    state.add_triangle(v1, v3, v2);
                }
            }
        } else {
            // Couldn't find a valid third vertex - re-add edge if still active
            if state.edge_can_add_triangle(v1, v2) {
                // Don't re-add, just leave the edge open
            }
        }
    }
}

/// Triangulate surface dots into a triangle mesh using fast local triangulation.
///
/// This uses a greedy local approach that's much faster than global edge-based
/// triangulation for large dot counts.
///
/// # Arguments
/// * `dots` - Surface dots with positions and normals
/// * `cutoff` - Maximum edge length for triangles
///
/// # Returns
/// Triangle indices as a flat vector (groups of 3)
pub fn triangulate_dots(dots: &[SurfaceDot], cutoff: f32) -> TriangulationResult {
    if dots.len() < 3 {
        return TriangulationResult {
            indices: Vec::new(),
            num_triangles: 0,
        };
    }

    // Use fast local triangulation for large dot counts
    if dots.len() > 1000 {
        return triangulate_dots_fast(dots, cutoff);
    }

    // Use detailed triangulation for small dot counts
    triangulate_dots_detailed(dots, cutoff)
}

/// Fast local triangulation - O(n) complexity
/// 
/// For each dot, finds nearby dots and creates triangles locally.
/// Much faster than global edge-based approach but may have some gaps.
fn triangulate_dots_fast(dots: &[SurfaceDot], cutoff: f32) -> TriangulationResult {
    let cutoff_sq = cutoff * cutoff;
    let spatial_hash = DotSpatialHash::new(dots, cutoff);
    
    let mut triangles: Vec<[u32; 3]> = Vec::with_capacity(dots.len() * 2);
    let mut edge_used: HashSet<(u32, u32)> = HashSet::with_capacity(dots.len() * 3);
    let mut neighbors = Vec::with_capacity(64);

    for i in 0..dots.len() {
        let pi = &dots[i].position;
        let ni = &dots[i].normal;
        spatial_hash.query_neighbors(pi, &mut neighbors);

        // Sort neighbors by distance for better triangle quality
        neighbors.sort_by(|&a, &b| {
            let da = dist_sq(pi, &dots[a as usize].position);
            let db = dist_sq(pi, &dots[b as usize].position);
            da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
        });

        // Try to form triangles with nearby pairs
        for (idx_j, &j) in neighbors.iter().enumerate() {
            if j as usize == i {
                continue;
            }

            let pj = &dots[j as usize].position;
            let d_ij = dist_sq(pi, pj);
            if d_ij > cutoff_sq {
                continue;
            }

            // Look for a third vertex to complete a triangle
            for &k in neighbors.iter().skip(idx_j + 1) {
                if k as usize == i || k == j {
                    continue;
                }

                let pk = &dots[k as usize].position;
                let d_ik = dist_sq(pi, pk);
                let d_jk = dist_sq(pj, pk);

                if d_ik > cutoff_sq || d_jk > cutoff_sq {
                    continue;
                }

                // Create ordered triangle indices
                let mut tri = [i as u32, j, k];
                tri.sort();

                // Check if we've already created this triangle
                let edge1 = (tri[0], tri[1]);
                let edge2 = (tri[1], tri[2]);
                let edge3 = (tri[0], tri[2]);

                // Only add if at least one edge is new (prevents duplicate triangles)
                if edge_used.contains(&edge1) && edge_used.contains(&edge2) && edge_used.contains(&edge3) {
                    continue;
                }

                // Check normal consistency
                let nj = &dots[j as usize].normal;
                let nk = &dots[k as usize].normal;
                
                // Compute triangle normal
                let e1 = [pj[0] - pi[0], pj[1] - pi[1], pj[2] - pi[2]];
                let e2 = [pk[0] - pi[0], pk[1] - pi[1], pk[2] - pi[2]];
                let tri_normal = [
                    e1[1] * e2[2] - e1[2] * e2[1],
                    e1[2] * e2[0] - e1[0] * e2[2],
                    e1[0] * e2[1] - e1[1] * e2[0],
                ];

                // Check all vertex normals point same side as triangle normal
                let avg_normal = [ni[0] + nj[0] + nk[0], ni[1] + nj[1] + nk[1], ni[2] + nj[2] + nk[2]];
                let dot = tri_normal[0] * avg_normal[0] + tri_normal[1] * avg_normal[1] + tri_normal[2] * avg_normal[2];

                if dot.abs() < 0.001 {
                    continue; // Degenerate triangle
                }

                // Add triangle with correct winding
                if dot > 0.0 {
                    triangles.push([i as u32, j, k]);
                } else {
                    triangles.push([i as u32, k, j]);
                }

                edge_used.insert(edge1);
                edge_used.insert(edge2);
                edge_used.insert(edge3);

                // Limit triangles per vertex for performance
                if triangles.len() > dots.len() * 3 {
                    break;
                }
            }
        }
    }

    let indices: Vec<u32> = triangles.iter().flat_map(|t| [t[0], t[1], t[2]]).collect();
    TriangulationResult {
        num_triangles: triangles.len(),
        indices,
    }
}

/// Detailed edge-based triangulation - higher quality but slower
fn triangulate_dots_detailed(dots: &[SurfaceDot], cutoff: f32) -> TriangulationResult {
    let cutoff_sq = cutoff * cutoff;

    // Build spatial hash for dot queries
    let spatial_hash = DotSpatialHash::new(dots, cutoff);

    // Estimate number of triangles (roughly 2x number of dots for a closed surface)
    let estimated_triangles = dots.len() * 2;

    let mut state = TriangulationState::new(estimated_triangles);

    // Phase 1: Find seed triangles
    find_seed_triangles(dots, &spatial_hash, cutoff_sq, &mut state);

    // Phase 2: Expand mesh from active edges (with conservative limit)
    let max_iterations = dots.len().min(10000) * 5;
    expand_mesh(dots, &spatial_hash, cutoff_sq, &mut state, max_iterations);

    // Convert triangles to flat index buffer
    let indices: Vec<u32> = state
        .triangles
        .iter()
        .flat_map(|tri| [tri[0], tri[1], tri[2]])
        .collect();

    TriangulationResult {
        num_triangles: state.triangles.len(),
        indices,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_dot(x: f32, y: f32, z: f32, nx: f32, ny: f32, nz: f32) -> SurfaceDot {
        SurfaceDot {
            position: [x, y, z],
            normal: [nx, ny, nz],
            atom_index: 0,
        }
    }

    #[test]
    fn test_triangulate_tetrahedron() {
        // Four dots forming a tetrahedron
        let dots = vec![
            make_dot(0.0, 0.0, 0.0, 0.0, 0.0, -1.0),
            make_dot(1.0, 0.0, 0.0, 1.0, 0.0, 0.0),
            make_dot(0.5, 0.866, 0.0, 0.0, 1.0, 0.0),
            make_dot(0.5, 0.289, 0.816, 0.0, 0.0, 1.0),
        ];

        let result = triangulate_dots(&dots, 2.0);

        // Should have created some triangles
        assert!(!result.indices.is_empty(), "Should create triangles");
        assert_eq!(result.indices.len() % 3, 0, "Indices should be multiple of 3");
    }

    #[test]
    fn test_triangulate_plane() {
        // Dots on a plane - should triangulate well
        let mut dots = Vec::new();
        for i in 0..5 {
            for j in 0..5 {
                dots.push(make_dot(
                    i as f32 * 0.5,
                    j as f32 * 0.5,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                ));
            }
        }

        let result = triangulate_dots(&dots, 1.0);

        // Should have created triangles
        assert!(!result.indices.is_empty(), "Should create triangles for plane");
    }

    #[test]
    fn test_triangulate_empty() {
        let dots: Vec<SurfaceDot> = vec![];
        let result = triangulate_dots(&dots, 1.0);

        assert!(result.indices.is_empty());
        assert_eq!(result.num_triangles, 0);
    }

    #[test]
    fn test_triangulate_too_few_dots() {
        let dots = vec![
            make_dot(0.0, 0.0, 0.0, 0.0, 0.0, 1.0),
            make_dot(1.0, 0.0, 0.0, 0.0, 0.0, 1.0),
        ];

        let result = triangulate_dots(&dots, 2.0);

        assert!(result.indices.is_empty());
        assert_eq!(result.num_triangles, 0);
    }
}
