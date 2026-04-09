//! Lightweight BVH for CPU ray-picking.
//!
//! Flat-array AABB tree built with median splits. Leaf nodes store ranges
//! into a sorted atom buffer. Ray traversal is iterative (explicit stack).

use pymol_mol::{AtomIndex, ObjectMolecule, RepMask};
use pymol_render::picking::Ray;

/// Maximum atoms in a leaf node before we split.
const LEAF_THRESHOLD: usize = 8;

/// A sphere (atom) stored in the BVH leaf buffer.
#[derive(Clone, Copy)]
struct Sphere {
    index: AtomIndex,
    center: [f32; 3],
    radius: f32,
}

/// Axis-aligned bounding box.
#[derive(Clone, Copy)]
struct Aabb {
    min: [f32; 3],
    max: [f32; 3],
}

impl Aabb {
    fn empty() -> Self {
        Self {
            min: [f32::INFINITY; 3],
            max: [f32::NEG_INFINITY; 3],
        }
    }

    fn expand_sphere(&mut self, center: [f32; 3], radius: f32) {
        for (i, &c) in center.iter().enumerate() {
            self.min[i] = self.min[i].min(c - radius);
            self.max[i] = self.max[i].max(c + radius);
        }
    }

    fn longest_axis(&self) -> usize {
        let dx = self.max[0] - self.min[0];
        let dy = self.max[1] - self.min[1];
        let dz = self.max[2] - self.min[2];
        if dx >= dy && dx >= dz {
            0
        } else if dy >= dz {
            1
        } else {
            2
        }
    }
}

/// A node in the flat BVH array.
///
/// Internal nodes: `count == 0`, `offset` = index of right child (left child is `self + 1`).
/// Leaf nodes: `count > 0`, `offset` = start index into `spheres`.
#[derive(Clone, Copy)]
struct BvhNode {
    aabb: Aabb,
    /// For leaves: number of spheres. For internal nodes: 0.
    count: u32,
    /// For leaves: start index in `spheres`. For internal: right-child index in `nodes`.
    offset: u32,
}

/// Cached BVH for fast ray-picking on a molecule.
pub struct PickBvh {
    nodes: Vec<BvhNode>,
    spheres: Vec<Sphere>,
}

impl PickBvh {
    /// Build a BVH from a molecule's current display state.
    ///
    /// Only atoms whose per-atom `visible_reps` intersect `obj_reps` are included.
    pub fn build(molecule: &ObjectMolecule, state: usize, obj_reps: RepMask) -> Self {
        let coord_set = match molecule.get_coord_set(state) {
            Some(cs) => cs,
            None => {
                return Self {
                    nodes: Vec::new(),
                    spheres: Vec::new(),
                };
            }
        };

        // Collect pickable atoms
        let mut spheres: Vec<Sphere> = Vec::with_capacity(molecule.atom_count());
        for (atom_idx, coord) in coord_set.iter_with_atoms() {
            if let Some(atom) = molecule.get_atom(atom_idx) {
                if atom.repr.visible_reps.intersection(obj_reps) != RepMask::NONE {
                    spheres.push(Sphere {
                        index: atom_idx,
                        center: [coord.x, coord.y, coord.z],
                        radius: atom.effective_vdw(),
                    });
                }
            }
        }

        if spheres.is_empty() {
            return Self {
                nodes: Vec::new(),
                spheres,
            };
        }

        // Pre-allocate nodes (at most 2n-1 for n leaves)
        let max_nodes = 2 * spheres.len();
        let mut nodes: Vec<BvhNode> = Vec::with_capacity(max_nodes.min(4096));

        // Build recursively
        let len = spheres.len();
        Self::build_recursive(&mut nodes, &mut spheres, 0, len);

        Self { nodes, spheres }
    }

    fn build_recursive(
        nodes: &mut Vec<BvhNode>,
        spheres: &mut [Sphere],
        start: usize,
        end: usize,
    ) -> usize {
        let node_idx = nodes.len();

        // Compute AABB for this range
        let mut aabb = Aabb::empty();
        for s in &spheres[start..end] {
            aabb.expand_sphere(s.center, s.radius);
        }

        let count = end - start;

        if count <= LEAF_THRESHOLD {
            // Leaf
            nodes.push(BvhNode {
                aabb,
                count: count as u32,
                offset: start as u32,
            });
            return node_idx;
        }

        // Split on longest axis at median
        let axis = aabb.longest_axis();
        let mid = start + count / 2;

        // Partial sort to find median — O(n) average
        spheres[start..end].select_nth_unstable_by(mid - start, |a, b| {
            a.center[axis]
                .partial_cmp(&b.center[axis])
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        // Reserve space for this internal node (will fill in right_child later)
        nodes.push(BvhNode {
            aabb,
            count: 0,
            offset: 0, // placeholder
        });

        // Build left child (immediately follows this node)
        Self::build_recursive(nodes, spheres, start, mid);

        // Build right child
        let right_idx = nodes.len();
        nodes[node_idx].offset = right_idx as u32;
        Self::build_recursive(nodes, spheres, mid, end);

        node_idx
    }

    /// Find the closest atom hit by a ray.
    ///
    /// Returns `(atom_index, distance)` for the closest intersection, or `None`.
    pub fn intersect_ray(&self, ray: &Ray) -> Option<(AtomIndex, f32)> {
        if self.nodes.is_empty() {
            return None;
        }

        let mut closest_t = f32::MAX;
        let mut closest_atom: Option<AtomIndex> = None;

        // Iterative traversal with explicit stack
        let mut stack = [0u32; 64]; // 64 levels is more than enough
        let mut stack_ptr = 0usize;
        stack[0] = 0;
        stack_ptr += 1;

        while stack_ptr > 0 {
            stack_ptr -= 1;
            let node_idx = stack[stack_ptr] as usize;
            let node = &self.nodes[node_idx];

            // Test ray against node AABB
            if !ray_aabb_hit(ray, &node.aabb, closest_t) {
                continue;
            }

            if node.count > 0 {
                // Leaf: test all spheres
                let start = node.offset as usize;
                let end = start + node.count as usize;
                for sphere in &self.spheres[start..end] {
                    if let Some(t) = ray.intersect_sphere(sphere.center, sphere.radius) {
                        if t < closest_t {
                            closest_t = t;
                            closest_atom = Some(sphere.index);
                        }
                    }
                }
            } else {
                // Internal: push children (right first so left is processed first)
                let right = node.offset;
                let left = (node_idx + 1) as u32;
                stack[stack_ptr] = right;
                stack_ptr += 1;
                stack[stack_ptr] = left;
                stack_ptr += 1;
            }
        }

        closest_atom.map(|idx| (idx, closest_t))
    }
}

/// Fast ray-AABB slab test. Returns true if the ray hits the box
/// at a distance less than `max_t`.
fn ray_aabb_hit(ray: &Ray, aabb: &Aabb, max_t: f32) -> bool {
    let mut tmin = 0.0_f32;
    let mut tmax = max_t;

    for i in 0..3 {
        let inv_d = 1.0 / ray.direction[i];
        let mut t0 = (aabb.min[i] - ray.origin[i]) * inv_d;
        let mut t1 = (aabb.max[i] - ray.origin[i]) * inv_d;
        if inv_d < 0.0 {
            std::mem::swap(&mut t0, &mut t1);
        }
        tmin = tmin.max(t0);
        tmax = tmax.min(t1);
        if tmax < tmin {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use lin_alg::f32::Vec3;
    use pymol_mol::{AtomBuilder, Element, MoleculeBuilder, RepMask as MolRepMask};

    fn test_molecule() -> ObjectMolecule {
        let mut mol = MoleculeBuilder::new("test")
            .add_atom(
                AtomBuilder::new()
                    .name("CA")
                    .element(Element::Carbon)
                    .chain("A")
                    .resn("ALA")
                    .resv(1)
                    .build(),
                Vec3::new(0.0, 0.0, 0.0),
            )
            .add_atom(
                AtomBuilder::new()
                    .name("CB")
                    .element(Element::Carbon)
                    .chain("A")
                    .resn("ALA")
                    .resv(1)
                    .build(),
                Vec3::new(10.0, 0.0, 0.0),
            )
            .add_atom(
                AtomBuilder::new()
                    .name("N")
                    .element(Element::Nitrogen)
                    .chain("A")
                    .resn("ALA")
                    .resv(1)
                    .build(),
                Vec3::new(0.0, 10.0, 0.0),
            )
            .build();

        // Make all atoms visible
        for atom in mol.atoms_mut() {
            atom.repr.visible_reps = MolRepMask::ALL;
        }
        mol
    }

    #[test]
    fn test_bvh_build_empty() {
        let mol = MoleculeBuilder::new("empty").build();
        let bvh = PickBvh::build(&mol, 0, RepMask::ALL);
        assert!(bvh.nodes.is_empty());
        assert!(bvh.spheres.is_empty());
    }

    #[test]
    fn test_bvh_ray_hit() {
        let mol = test_molecule();
        let bvh = PickBvh::build(&mol, 0, RepMask::ALL);

        // Ray pointing at first atom (0, 0, 0)
        let ray = Ray::new([0.0, 0.0, 10.0], [0.0, 0.0, -1.0]);
        let result = bvh.intersect_ray(&ray);
        assert!(result.is_some());
        let (idx, t) = result.unwrap();
        assert_eq!(idx, AtomIndex::from(0usize));
        assert!(t > 0.0);
    }

    #[test]
    fn test_bvh_ray_miss() {
        let mol = test_molecule();
        let bvh = PickBvh::build(&mol, 0, RepMask::ALL);

        // Ray far from all atoms
        let ray = Ray::new([100.0, 100.0, 10.0], [0.0, 0.0, -1.0]);
        assert!(bvh.intersect_ray(&ray).is_none());
    }

    #[test]
    fn test_bvh_finds_closest() {
        let mol = test_molecule();
        let bvh = PickBvh::build(&mol, 0, RepMask::ALL);

        // Ray aimed at atom at (10, 0, 0) — should hit that one, not (0,0,0)
        let ray = Ray::new([10.0, 0.0, 10.0], [0.0, 0.0, -1.0]);
        let result = bvh.intersect_ray(&ray);
        assert!(result.is_some());
        let (idx, _t) = result.unwrap();
        assert_eq!(idx, AtomIndex::from(1usize));
    }
}
