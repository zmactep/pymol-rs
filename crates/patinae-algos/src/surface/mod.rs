//! Surface support types shared by CPU helpers and GPU render code.

pub mod contour;
pub mod grid;

pub use contour::{
    extract_isomesh, extract_isosurface, ContourGeometry, ContourTopology, ContourVertex,
};
pub use grid::{Grid3D, SpatialHash};

/// Sentinel for vertices with no owning atom (e.g. outside any narrowband).
pub const NO_OWNER: u32 = u32::MAX;

/// Atom data for surface generation.
#[derive(Debug, Clone)]
pub struct SurfaceAtom {
    pub position: [f32; 3],
    pub radius: f32,
    pub atom_index: usize,
}

/// Output of surface generation — rendering-independent triangle mesh.
#[derive(Debug, Default, Clone)]
pub struct SurfaceMesh {
    pub positions: Vec<[f32; 3]>,
    pub normals: Vec<[f32; 3]>,
    pub atom_ids: Vec<u32>,
    pub indices: Vec<u32>,
}

/// Build a molecular surface mesh.
///
/// Stub: always returns an empty mesh until a CPU implementation is needed by
/// callers outside the renderer.
pub fn build_surface(_atoms: &[SurfaceAtom], _probe_radius: f32, _quality: i32) -> SurfaceMesh {
    SurfaceMesh::default()
}

/// Build a `SpatialHash` over `atoms` using their VDW radii. Surface ownership
/// helpers use it to find nearby atoms during color and owner assignment.
pub fn create_spatial_hash(atoms: &[SurfaceAtom]) -> SpatialHash {
    if atoms.is_empty() {
        return SpatialHash::new(&[], &[], [0.0; 3], [0.0; 3]);
    }
    let positions: Vec<[f32; 3]> = atoms.iter().map(|a| a.position).collect();
    let radii: Vec<f32> = atoms.iter().map(|a| a.radius).collect();
    let mut min = positions[0];
    let mut max = positions[0];
    for (p, r) in positions.iter().zip(&radii) {
        for i in 0..3 {
            if p[i] - r < min[i] {
                min[i] = p[i] - r;
            }
            if p[i] + r > max[i] {
                max[i] = p[i] + r;
            }
        }
    }
    SpatialHash::new(&positions, &radii, min, max)
}
