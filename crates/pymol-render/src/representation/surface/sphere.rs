//! Pre-tessellated sphere data for dot-based surface generation.
//!
//! This module provides icosahedron-based sphere tessellations at various
//! subdivision levels. Each level quadruples the number of triangles.
//!
//! | Level | Dots  | Triangles |
//! |-------|-------|-----------|
//! | 0     | 12    | 20        |
//! | 1     | 42    | 80        |
//! | 2     | 162   | 320       |
//! | 3     | 642   | 1280      |
//! | 4     | 2562  | 5120      |

use std::collections::HashMap;
use std::sync::OnceLock;

/// Pre-tessellated sphere with dots on unit sphere surface.
#[derive(Debug, Clone)]
pub struct SphereTessellation {
    /// Unit sphere dot positions (also serve as normals)
    pub dots: Vec<[f32; 3]>,
    /// Triangle indices for the tessellation (for direct rendering if needed)
    #[allow(dead_code)]
    pub triangles: Vec<[u32; 3]>,
}

impl SphereTessellation {
    /// Number of dots in this tessellation
    #[inline]
    pub fn num_dots(&self) -> usize {
        self.dots.len()
    }

    /// Number of triangles in this tessellation
    #[inline]
    #[allow(dead_code)]
    pub fn num_triangles(&self) -> usize {
        self.triangles.len()
    }
}

/// Golden ratio constant for icosahedron construction
const PHI: f32 = 1.618033988749895;

/// Generate the base icosahedron (12 vertices, 20 faces)
fn create_icosahedron() -> (Vec<[f32; 3]>, Vec<[u32; 3]>) {
    // Normalize factor for vertices at distance 1 from origin
    let norm = (1.0 + PHI * PHI).sqrt();
    let a = 1.0 / norm;
    let b = PHI / norm;

    // 12 vertices of icosahedron
    let vertices = vec![
        [-a, b, 0.0],
        [a, b, 0.0],
        [-a, -b, 0.0],
        [a, -b, 0.0],
        [0.0, -a, b],
        [0.0, a, b],
        [0.0, -a, -b],
        [0.0, a, -b],
        [b, 0.0, -a],
        [b, 0.0, a],
        [-b, 0.0, -a],
        [-b, 0.0, a],
    ];

    // 20 triangular faces
    let triangles = vec![
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ];

    (vertices, triangles)
}

/// Normalize a 3D vector to unit length
#[inline]
fn normalize(v: [f32; 3]) -> [f32; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len > 0.0 {
        [v[0] / len, v[1] / len, v[2] / len]
    } else {
        v
    }
}

/// Get or create a midpoint vertex, returning its index
fn get_midpoint(
    v1: u32,
    v2: u32,
    vertices: &mut Vec<[f32; 3]>,
    midpoint_cache: &mut HashMap<(u32, u32), u32>,
) -> u32 {
    // Use ordered pair as key to avoid duplicates
    let key = if v1 < v2 { (v1, v2) } else { (v2, v1) };

    if let Some(&idx) = midpoint_cache.get(&key) {
        return idx;
    }

    // Calculate midpoint and normalize to unit sphere
    let p1 = vertices[v1 as usize];
    let p2 = vertices[v2 as usize];
    let mid = normalize([
        (p1[0] + p2[0]) * 0.5,
        (p1[1] + p2[1]) * 0.5,
        (p1[2] + p2[2]) * 0.5,
    ]);

    let idx = vertices.len() as u32;
    vertices.push(mid);
    midpoint_cache.insert(key, idx);
    idx
}

/// Subdivide a sphere tessellation by one level.
/// Each triangle is split into 4 triangles.
fn subdivide(vertices: &mut Vec<[f32; 3]>, triangles: Vec<[u32; 3]>) -> Vec<[u32; 3]> {
    let mut midpoint_cache: HashMap<(u32, u32), u32> = HashMap::new();
    let mut new_triangles = Vec::with_capacity(triangles.len() * 4);

    for tri in triangles {
        let v0 = tri[0];
        let v1 = tri[1];
        let v2 = tri[2];

        // Get midpoints of each edge
        let m01 = get_midpoint(v0, v1, vertices, &mut midpoint_cache);
        let m12 = get_midpoint(v1, v2, vertices, &mut midpoint_cache);
        let m20 = get_midpoint(v2, v0, vertices, &mut midpoint_cache);

        // Create 4 new triangles
        new_triangles.push([v0, m01, m20]);
        new_triangles.push([v1, m12, m01]);
        new_triangles.push([v2, m20, m12]);
        new_triangles.push([m01, m12, m20]);
    }

    new_triangles
}

/// Generate a tessellated sphere at the given subdivision level.
pub fn generate_sphere(level: usize) -> SphereTessellation {
    let (mut vertices, mut triangles) = create_icosahedron();

    for _ in 0..level {
        triangles = subdivide(&mut vertices, triangles);
    }

    SphereTessellation {
        dots: vertices,
        triangles,
    }
}

// Static storage for pre-computed spheres (lazily initialized)
static SPHERE_LEVEL_0: OnceLock<SphereTessellation> = OnceLock::new();
static SPHERE_LEVEL_1: OnceLock<SphereTessellation> = OnceLock::new();
static SPHERE_LEVEL_2: OnceLock<SphereTessellation> = OnceLock::new();
static SPHERE_LEVEL_3: OnceLock<SphereTessellation> = OnceLock::new();
static SPHERE_LEVEL_4: OnceLock<SphereTessellation> = OnceLock::new();

/// Get a pre-computed sphere tessellation at the given level (0-4).
///
/// Higher levels have more dots and better surface quality but slower performance.
/// - Level 0: 12 dots (very coarse)
/// - Level 1: 42 dots (coarse)
/// - Level 2: 162 dots (medium - recommended for most uses)
/// - Level 3: 642 dots (fine)
/// - Level 4: 2562 dots (very fine)
///
/// Levels > 4 are clamped to 4.
pub fn get_sphere(level: usize) -> &'static SphereTessellation {
    match level.min(4) {
        0 => SPHERE_LEVEL_0.get_or_init(|| generate_sphere(0)),
        1 => SPHERE_LEVEL_1.get_or_init(|| generate_sphere(1)),
        2 => SPHERE_LEVEL_2.get_or_init(|| generate_sphere(2)),
        3 => SPHERE_LEVEL_3.get_or_init(|| generate_sphere(3)),
        _ => SPHERE_LEVEL_4.get_or_init(|| generate_sphere(4)),
    }
}

/// Convert quality setting (-4 to 4) to sphere tessellation level (0 to 4).
///
/// Lower quality uses fewer dots for faster computation.
pub fn quality_to_level(quality: i32) -> usize {
    match quality {
        q if q <= -3 => 0,
        -2 => 1,
        -1 | 0 => 2,
        1 | 2 => 3,
        _ => 4,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_icosahedron() {
        let sphere = get_sphere(0);
        assert_eq!(sphere.num_dots(), 12);
        assert_eq!(sphere.num_triangles(), 20);
    }

    #[test]
    fn test_sphere_levels() {
        // Verify dot counts match PyMOL's SphereData.h
        assert_eq!(get_sphere(0).num_dots(), 12);
        assert_eq!(get_sphere(1).num_dots(), 42);
        assert_eq!(get_sphere(2).num_dots(), 162);
        assert_eq!(get_sphere(3).num_dots(), 642);
        assert_eq!(get_sphere(4).num_dots(), 2562);
    }

    #[test]
    fn test_unit_sphere() {
        let sphere = get_sphere(2);
        for dot in &sphere.dots {
            let len = (dot[0] * dot[0] + dot[1] * dot[1] + dot[2] * dot[2]).sqrt();
            assert!((len - 1.0).abs() < 1e-6, "Dot should be on unit sphere");
        }
    }

    #[test]
    fn test_triangle_indices_valid() {
        let sphere = get_sphere(2);
        let n = sphere.num_dots() as u32;
        for tri in &sphere.triangles {
            assert!(tri[0] < n);
            assert!(tri[1] < n);
            assert!(tri[2] < n);
        }
    }
}
