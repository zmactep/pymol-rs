//! Primitive collection from molecular data
//!
//! This module provides functions to extract raytracing primitives from
//! molecular structures, similar to how the render representations work.

use std::collections::HashMap;

use pymol_color::{ChainColors, ElementColors, NamedColors};
use pymol_mol::{AtomIndex, CoordSet, ObjectMolecule, RepMask};

use crate::primitive::{GpuCylinder, GpuSphere, PrimitiveCollector, Primitives};

/// Color resolver for raytracing (simplified version)
pub struct RayColorResolver<'a> {
    named_colors: &'a NamedColors,
    element_colors: &'a ElementColors,
}

impl<'a> RayColorResolver<'a> {
    /// Create a new color resolver
    pub fn new(
        named_colors: &'a NamedColors,
        element_colors: &'a ElementColors,
        _chain_colors: &ChainColors, // Kept for API compatibility
    ) -> Self {
        Self {
            named_colors,
            element_colors,
        }
    }

    /// Resolve color for an atom
    pub fn resolve_atom(
        &self,
        atom: &pymol_mol::Atom,
        _molecule: &ObjectMolecule,
    ) -> [f32; 4] {
        // Check for explicit atom color first (positive index means custom color)
        if atom.colors.base >= 0 {
            if let Some(color) = self.named_colors.get_by_index(atom.colors.base as u32) {
                return [color.r, color.g, color.b, 1.0];
            }
        }

        // Try element color
        let color = self.element_colors.get(atom.element.atomic_number());
        if color.r != 1.0 || color.g != 0.0 || color.b != 1.0 {
            // Not the default magenta, so it's a valid element color
            return [color.r, color.g, color.b, 1.0];
        }

        // Try chain color
        if !atom.chain.is_empty() {
            let color = ChainColors::get(&atom.chain);
            return [color.r, color.g, color.b, 1.0];
        }

        // Default to gray (carbon-like)
        [0.56, 0.56, 0.56, 1.0]
    }
}

/// Collect sphere primitives from a molecule
///
/// Extracts atoms visible as spheres representation.
pub fn collect_spheres(
    molecule: &ObjectMolecule,
    coord_set: &CoordSet,
    colors: &RayColorResolver,
    sphere_scale: f32,
    sphere_transparency: f32,
) -> Vec<GpuSphere> {
    let mut spheres = Vec::new();

    for (atom_idx, coord) in coord_set.iter_with_atoms() {
        let atom = match molecule.get_atom(atom_idx) {
            Some(a) => a,
            None => continue,
        };

        // Check if spheres representation is visible
        if !atom.visible_reps.is_visible(RepMask::SPHERES) {
            continue;
        }

        // Use per-atom sphere_scale if set, otherwise global
        let scale = atom.sphere_scale.unwrap_or(sphere_scale);
        let radius = atom.effective_vdw() * scale;
        let color = colors.resolve_atom(atom, molecule);

        spheres.push(GpuSphere::new(
            [coord.x, coord.y, coord.z],
            radius,
            color,
            sphere_transparency,
        ));
    }

    spheres
}

/// Collect cylinder primitives from a molecule (sticks representation)
///
/// Extracts bonds visible as sticks representation.
pub fn collect_cylinders(
    molecule: &ObjectMolecule,
    coord_set: &CoordSet,
    colors: &RayColorResolver,
    stick_radius: f32,
    stick_transparency: f32,
) -> (Vec<GpuCylinder>, Vec<GpuSphere>) {
    let mut cylinders = Vec::new();
    let mut atom_max_radius: HashMap<AtomIndex, (f32, [f32; 4])> = HashMap::new();

    for bond in molecule.bonds() {
        let atom1_idx = bond.atom1;
        let atom2_idx = bond.atom2;

        let atom1 = match molecule.get_atom(atom1_idx) {
            Some(a) => a,
            None => continue,
        };
        let atom2 = match molecule.get_atom(atom2_idx) {
            Some(a) => a,
            None => continue,
        };

        // Check if sticks representation is visible
        if !atom1.visible_reps.is_visible(RepMask::STICKS)
            && !atom2.visible_reps.is_visible(RepMask::STICKS)
        {
            continue;
        }

        let pos1 = match coord_set.get_atom_coord(atom1_idx) {
            Some(p) => p,
            None => continue,
        };
        let pos2 = match coord_set.get_atom_coord(atom2_idx) {
            Some(p) => p,
            None => continue,
        };

        let color1 = colors.resolve_atom(atom1, molecule);
        let color2 = colors.resolve_atom(atom2, molecule);
        let radius = stick_radius * bond.order.as_float().sqrt();

        cylinders.push(GpuCylinder::new(
            [pos1.x, pos1.y, pos1.z],
            [pos2.x, pos2.y, pos2.z],
            radius,
            color1,
            color2,
            stick_transparency,
        ));

        // Track atoms for sphere caps
        atom_max_radius
            .entry(atom1_idx)
            .and_modify(|(r, _)| *r = r.max(radius))
            .or_insert((radius, color1));
        atom_max_radius
            .entry(atom2_idx)
            .and_modify(|(r, _)| *r = r.max(radius))
            .or_insert((radius, color2));
    }

    // Create sphere caps
    let mut sphere_caps = Vec::new();
    for (atom_idx, (max_radius, color)) in atom_max_radius {
        if let Some(coord) = coord_set.get_atom_coord(atom_idx) {
            sphere_caps.push(GpuSphere::new(
                [coord.x, coord.y, coord.z],
                max_radius,
                color,
                stick_transparency,
            ));
        }
    }

    (cylinders, sphere_caps)
}

/// Options for primitive collection
#[derive(Clone, Debug)]
pub struct CollectOptions {
    /// Sphere scale factor
    pub sphere_scale: f32,
    /// Stick radius
    pub stick_radius: f32,
    /// Default sphere transparency
    pub sphere_transparency: f32,
    /// Default stick transparency
    pub stick_transparency: f32,
    /// Collect spheres representation
    pub collect_spheres: bool,
    /// Collect sticks representation  
    pub collect_sticks: bool,
    /// Collect cartoon representation
    pub collect_cartoon: bool,
    /// Collect surface representation
    pub collect_surface: bool,
}

impl Default for CollectOptions {
    fn default() -> Self {
        Self {
            sphere_scale: 1.0,
            stick_radius: 0.25,
            sphere_transparency: 0.0,
            stick_transparency: 0.0,
            collect_spheres: true,
            collect_sticks: true,
            collect_cartoon: true,
            collect_surface: true,
        }
    }
}

/// Collect all primitives from a molecule
pub fn collect_from_molecule(
    molecule: &ObjectMolecule,
    coord_set: &CoordSet,
    colors: &RayColorResolver,
    options: &CollectOptions,
) -> Primitives {
    let mut collector = PrimitiveCollector::new();

    // Collect spheres
    if options.collect_spheres {
        let spheres = collect_spheres(
            molecule,
            coord_set,
            colors,
            options.sphere_scale,
            options.sphere_transparency,
        );
        collector.add_spheres(spheres);
    }

    // Collect sticks (cylinders + sphere caps)
    if options.collect_sticks {
        let (cylinders, caps) = collect_cylinders(
            molecule,
            coord_set,
            colors,
            options.stick_radius,
            options.stick_transparency,
        );
        collector.add_cylinders(cylinders);
        collector.add_spheres(caps);
    }

    // TODO: Collect cartoon representation (triangles)
    // TODO: Collect surface representation (triangles)

    collector.build()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_collect_options_default() {
        let opts = CollectOptions::default();
        assert_eq!(opts.sphere_scale, 1.0);
        assert_eq!(opts.stick_radius, 0.25);
    }
}
