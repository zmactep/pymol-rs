//! Color resolution for molecular visualization.
//!
//! Resolves `ColorIndex` values from atoms to final RGBA values. Reads only
//! `NamedPalette`, `ThemedPalette`, the optional per-chain
//! `patinae_mol::ResidueRankTable` (for spectrum sampling) and `patinae_mol::Atom`
//! itself — no rendering crate is touched. Lives in in this crate so the
//! viewport renderer and headless code paths can both use it without pulling
//! in a concrete rendering backend.

use crate::{Color, ColorIndex, NamedPalette, ThemedPalette};
use patinae_mol::{Atom, ResidueRankTable};

/// Resolves color indices to final RGBA values.
pub struct ColorResolver<'a> {
    named_palette: &'a NamedPalette,
    palette: &'a ThemedPalette,
    default_alpha: f32,
    b_factor_range: (f32, f32),
    residue_ranks: Option<&'a ResidueRankTable>,
}

impl<'a> ColorResolver<'a> {
    pub fn new(named_palette: &'a NamedPalette, palette: &'a ThemedPalette) -> Self {
        Self {
            named_palette,
            palette,
            default_alpha: 1.0,
            b_factor_range: (0.0, 100.0),
            residue_ranks: None,
        }
    }

    pub fn palette(&self) -> &ThemedPalette {
        self.palette
    }

    pub fn with_alpha(mut self, alpha: f32) -> Self {
        self.default_alpha = alpha;
        self
    }

    pub fn with_b_factor_range(mut self, min: f32, max: f32) -> Self {
        self.b_factor_range = (min, max);
        self
    }

    /// Provide a per-chain polymer-residue rank table; `ColorIndex::ByResidueIndex`
    /// will sample the spectrum at `rank / (chain_total - 1)` so each chain
    /// reads blue → red end-to-end regardless of resv numbering or gaps.
    /// Pass `None` to leave the table unset (non-polymer atoms always fall
    /// back to element color in that case).
    pub fn with_residue_ranks(mut self, ranks: Option<&'a ResidueRankTable>) -> Self {
        self.residue_ranks = ranks;
        self
    }

    pub fn resolve_atom(&self, atom: &Atom) -> [f32; 4] {
        let color = self.resolve_color_index(atom);
        color.to_rgba(self.default_alpha)
    }

    pub fn resolve_atom_with_alpha(&self, atom: &Atom, alpha: f32) -> [f32; 4] {
        let color = self.resolve_color_index(atom);
        color.to_rgba(alpha)
    }

    fn resolve_color_index(&self, atom: &Atom) -> Color {
        self.resolve_color_index_value(atom.repr.colors.base, atom)
    }

    fn resolve_color_index_value(&self, color_idx: i32, atom: &Atom) -> Color {
        match ColorIndex::from(color_idx) {
            ColorIndex::ByElement | ColorIndex::Atomic => {
                self.palette.element.get(atom.element as u8)
            }
            ColorIndex::ByChain => self.resolve_by_chain(atom),
            ColorIndex::BySS => self.palette.ss.get(atom.ss_type as u8),
            ColorIndex::ByBFactor => self.resolve_by_b_factor(atom),
            ColorIndex::ByResidueType => self.resolve_by_residue_type(atom),
            ColorIndex::ByResidueIndex => self.resolve_by_residue_index(atom),
            ColorIndex::Named(idx) => self.named_palette.get_by_index(idx).unwrap_or(Color::WHITE),
        }
    }

    fn resolve_by_chain(&self, atom: &Atom) -> Color {
        if atom.state.flags.is_biomolecule() {
            self.palette.chains.get(&atom.residue.chain)
        } else {
            self.palette.element.get(atom.element as u8)
        }
    }

    pub fn resolve_cartoon(&self, atom: &Atom) -> [f32; 4] {
        let color_idx = atom.repr.colors.cartoon_or_base();
        let color = self.resolve_color_index_value(color_idx, atom);
        color.to_rgba(self.default_alpha)
    }

    pub fn resolve_rep_color(
        &self,
        atom: &Atom,
        per_atom_color: i32,
        default_color: i32,
    ) -> [f32; 4] {
        let color_idx = if per_atom_color != patinae_mol::COLOR_UNSET {
            per_atom_color
        } else if default_color >= 0 {
            default_color
        } else {
            atom.repr.colors.base
        };
        let color = self.resolve_color_index_value(color_idx, atom);
        color.to_rgba(self.default_alpha)
    }

    pub fn resolve_cartoon_with_default(&self, atom: &Atom, default_color: i32) -> [f32; 4] {
        self.resolve_rep_color(atom, atom.repr.colors.cartoon, default_color)
    }

    fn resolve_by_residue_type(&self, atom: &Atom) -> Color {
        self.palette
            .residue
            .get(&atom.residue.resn)
            .unwrap_or_else(|| self.palette.element.get(atom.element as u8))
    }

    fn resolve_by_residue_index(&self, atom: &Atom) -> Color {
        // Polymer residues are coloured by per-chain rank (blue to red
        // end-to-end). Non-polymer atoms fall back to element colour so
        // ligands, waters, and ions do not inherit the polymer gradient.
        match self.residue_ranks.and_then(|t| t.t_for(atom)) {
            Some(t) => self.palette.spectrum.sample(t),
            None => self.palette.element.get(atom.element as u8),
        }
    }

    fn resolve_by_b_factor(&self, atom: &Atom) -> Color {
        let (min, max) = self.b_factor_range;
        let range = max - min;
        if range <= 0.0 {
            return Color::WHITE;
        }
        let t = ((atom.b_factor - min) / range).clamp(0.0, 1.0);
        self.palette.b_factor.sample(t)
    }

    pub fn get_named_color(&self, name: &str) -> Option<Color> {
        self.named_palette.get_by_name(name).map(|(_, c)| c)
    }

    pub fn get_color_by_index(&self, index: u32) -> Option<Color> {
        self.named_palette.get_by_index(index)
    }

    pub fn get_element_color(&self, atomic_number: u8) -> Color {
        self.palette.element.get(atomic_number)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use patinae_mol::Element;

    fn test_resolver() -> (NamedPalette, ThemedPalette) {
        (NamedPalette::default(), ThemedPalette::dark())
    }

    #[test]
    fn test_element_coloring() {
        let (named, palette) = test_resolver();
        let resolver = ColorResolver::new(&named, &palette);

        let mut atom = Atom {
            element: Element::Carbon,
            ..Default::default()
        };
        atom.repr.colors.base = -1;

        let color = resolver.resolve_atom(&atom);
        assert!(color[0] > 0.4 && color[0] < 0.7);
        assert!(color[1] > 0.4 && color[1] < 0.7);
        assert!(color[2] > 0.4 && color[2] < 0.7);
        assert_eq!(color[3], 1.0);
    }

    #[test]
    fn test_chain_coloring() {
        use patinae_mol::AtomFlags;

        let (named, palette) = test_resolver();
        let resolver = ColorResolver::new(&named, &palette);

        let mut atom = Atom {
            element: Element::Carbon,
            ..Default::default()
        };
        atom.set_residue("ALA", 1, "A");
        atom.state.flags = AtomFlags::PROTEIN | AtomFlags::POLYMER;
        atom.repr.colors.base = -2;
        let color_a = resolver.resolve_atom(&atom);
        let expected_a = palette.chains.get("A").to_rgba(1.0);
        assert_eq!(color_a, expected_a);

        atom.set_residue("ALA", 1, "B");
        atom.state.flags = AtomFlags::PROTEIN | AtomFlags::POLYMER;
        let color_b = resolver.resolve_atom(&atom);
        assert_ne!(color_a, color_b);

        let mut lig = Atom {
            element: Element::Carbon,
            ..Default::default()
        };
        lig.set_residue("LIG", 1, "Z");
        lig.state.flags = AtomFlags::ORGANIC;
        lig.repr.colors.base = -2;
        let lig_color = resolver.resolve_atom(&lig);
        let cpk_c = palette.element.get(Element::Carbon as u8).to_rgba(1.0);
        assert_eq!(lig_color, cpk_c);
    }

    #[test]
    fn test_residue_type_coloring() {
        let (named, palette) = test_resolver();
        let resolver = ColorResolver::new(&named, &palette);

        let mut atom = Atom {
            element: Element::Carbon,
            ..Default::default()
        };
        atom.set_residue("ALA", 1, "A");
        atom.repr.colors.base = -5;
        let color = resolver.resolve_atom(&atom);
        let expected = crate::RES_HYDROPHOBIC;
        assert_eq!(
            [color[0], color[1], color[2]],
            [expected.r, expected.g, expected.b]
        );

        let mut water = Atom {
            element: Element::Carbon,
            ..Default::default()
        };
        water.set_residue("HOH", 1, "A");
        water.repr.colors.base = -5;
        let w_color = resolver.resolve_atom(&water);
        let cpk_c = palette.element.get(Element::Carbon as u8);
        assert_eq!(w_color[0], cpk_c.r);
    }

    #[test]
    fn test_residue_index_coloring() {
        use patinae_mol::{polymer_residue_ranks, AtomFlags, ObjectMolecule};
        use std::sync::Arc;

        let mut mol = ObjectMolecule::new("t");
        for (resn, resv) in &[("ALA", 1), ("ALA", 50), ("ALA", 100)] {
            let mut atom = Atom {
                residue: Arc::new(patinae_mol::AtomResidue::from_parts(
                    "A", *resn, *resv, ' ', "",
                )),
                ..Default::default()
            };
            atom.state.flags = AtomFlags::PROTEIN | AtomFlags::POLYMER;
            atom.repr.colors.base = -6;
            mol.add_atom(atom);
        }
        let ranks = polymer_residue_ranks(&mol);

        let (named, palette) = test_resolver();
        let resolver = ColorResolver::new(&named, &palette).with_residue_ranks(Some(&ranks));

        let blue = resolver.resolve_atom(mol.get_atom(patinae_mol::AtomIndex(0)).unwrap());
        assert!(blue[2] > 0.9 && blue[0] < 0.1, "rank 0 should be blue");

        let green = resolver.resolve_atom(mol.get_atom(patinae_mol::AtomIndex(1)).unwrap());
        assert!(green[1] > 0.9, "rank 1 should be green");

        let red = resolver.resolve_atom(mol.get_atom(patinae_mol::AtomIndex(2)).unwrap());
        assert!(red[0] > 0.9 && red[2] < 0.1, "rank 2 should be red");
    }

    #[test]
    fn test_cartoon_with_default_priority() {
        use patinae_mol::COLOR_UNSET;

        let (named, palette) = test_resolver();
        let resolver = ColorResolver::new(&named, &palette);

        let mut atom = Atom {
            element: Element::Carbon,
            ..Default::default()
        };
        atom.repr.colors.base = -1;
        atom.repr.colors.cartoon = COLOR_UNSET;

        let color = resolver.resolve_cartoon_with_default(&atom, -1);
        assert!(color[0] > 0.4 && color[0] < 0.7);

        let color = resolver.resolve_cartoon_with_default(&atom, 2);
        assert!(color[0] > 0.9);
        assert!(color[1] < 0.1);

        atom.repr.colors.cartoon = 3;
        let color = resolver.resolve_cartoon_with_default(&atom, 2);
        assert!(color[1] > 0.9);
        assert!(color[0] < 0.1);

        atom.repr.colors.cartoon = COLOR_UNSET;
        let color = resolver.resolve_cartoon_with_default(&atom, 2);
        assert!(color[0] > 0.9);
    }
}
