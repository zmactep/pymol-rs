use lin_alg::f32::Vec3;
use patinae_color::{NamedPalette, ThemedPalette};
use patinae_mol::{AtomBuilder, MoleculeBuilder};
use patinae_render::{pack_rep_rgb8, REP_COLOR_INHERIT};
use patinae_scene::{bridge::ResolvedSceneColors, MoleculeObject, ObjectRegistry};
use patinae_settings::{Color as SettingColor, Settings};

fn single_atom_registry(atom: patinae_mol::Atom) -> ObjectRegistry {
    let mol = MoleculeBuilder::new("obj")
        .add_atom(atom, Vec3::new(0.0, 0.0, 0.0))
        .build();
    let mut registry = ObjectRegistry::new();
    registry.add(MoleculeObject::with_name(mol, "obj"));
    registry
}

#[test]
fn bridge_resolves_rep_specific_atom_colors() {
    let named = NamedPalette::new();
    let themed = ThemedPalette::dark();
    let red = named.get_by_name("red").unwrap().0 as i32;
    let green = named.get_by_name("green").unwrap().0 as i32;
    let blue = named.get_by_name("blue").unwrap().0 as i32;
    let yellow = named.get_by_name("yellow").unwrap().0 as i32;
    let cyan = named.get_by_name("cyan").unwrap().0 as i32;
    let magenta = named.get_by_name("magenta").unwrap().0 as i32;
    let white = named.get_by_name("white").unwrap().0 as i32;
    let gray = named.get_by_name("gray").unwrap().0 as i32;
    let black = named.get_by_name("black").unwrap().0 as i32;

    let mut atom = AtomBuilder::new().name("CA").element_symbol("C").build();
    atom.repr.colors.sphere = red;
    atom.repr.colors.stick = green;
    atom.repr.colors.line = blue;
    atom.repr.colors.dot = yellow;
    atom.repr.colors.cartoon = cyan;
    atom.repr.colors.ribbon = magenta;
    atom.repr.colors.surface = white;
    atom.repr.colors.mesh = gray;
    atom.repr.colors.ellipsoid = black;

    let registry = single_atom_registry(atom);
    let colors = ResolvedSceneColors::build(&registry, &Settings::default(), &named, &themed);
    let rep = colors.get_rep("obj").unwrap()[0];
    let expected = |idx: i32| pack_rep_rgb8(named.get_by_index(idx as u32).unwrap().to_rgba(1.0));

    assert_eq!(rep.sphere, expected(red));
    assert_eq!(rep.stick, expected(green));
    assert_eq!(rep.line, expected(blue));
    assert_eq!(rep.dot, expected(yellow));
    assert_eq!(rep.cartoon, expected(cyan));
    assert_eq!(rep.ribbon, expected(magenta));
    assert_eq!(rep.surface, expected(white));
    assert_eq!(rep.mesh, expected(gray));
    assert_eq!(rep.ellipsoid, expected(black));
}

#[test]
fn bridge_resolves_global_rep_color_defaults() {
    let named = NamedPalette::new();
    let themed = ThemedPalette::dark();
    let red = named.get_by_name("red").unwrap().0 as i32;

    let atom = AtomBuilder::new().name("CA").element_symbol("C").build();
    let registry = single_atom_registry(atom);
    let mut settings = Settings::default();
    settings.sphere.color = SettingColor(red);

    let colors = ResolvedSceneColors::build(&registry, &settings, &named, &themed);
    let rep = colors.get_rep("obj").unwrap()[0];

    assert_eq!(
        rep.sphere,
        pack_rep_rgb8(named.get_by_index(red as u32).unwrap().to_rgba(1.0))
    );
    assert_eq!(rep.stick, REP_COLOR_INHERIT);
}
