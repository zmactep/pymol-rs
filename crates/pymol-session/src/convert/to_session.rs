//! Convert [`PseSession`] → [`Session`]

use std::sync::Arc;

use ahash::AHashMap;
use lin_alg::f32::{Mat4, Vec3};

use pymol_color::{Color, ColorIndex, NamedColors};
use pymol_mol::{
    Atom, AtomColors, AtomRepresentation, AtomResidue, AtomIndex, BondOrder, BondStereo,
    CoordSet, Element, ObjectMolecule, RepMask, ResidueKey, SecondaryStructure, Symmetry,
};
use pymol_scene::{MoleculeObject, Object, SceneView, Session};
use pymol_settings::{SerializedSetting, SettingValue};

use super::pymol_colors::pymol_color_rgb;
use crate::pse::{
    PseAtom, PseCoordSet, PseMolecule, PseNameEntry, PseObject, PseObjectData,
    PseSelection, PseSession, PseSettingValue, PseSymmetry,
};
use crate::SessionError;

// =============================================================================
// PyMOL → pymol-rs color index remapping
// =============================================================================

/// Caches the mapping from PyMOL color indices to our NamedColors indices.
struct PseColorMapper {
    cache: AHashMap<i64, i32>,
}

impl PseColorMapper {
    fn new() -> Self {
        Self { cache: AHashMap::new() }
    }

    /// Map a PyMOL color index to our color index.
    ///
    /// Negative values have special meaning (-1=ByElement, etc.) and pass through.
    /// Positive values are looked up in the PyMOL color table, registered in
    /// `NamedColors` if needed, and the resulting local index is returned.
    fn map(&mut self, pymol_index: i64, named_colors: &mut NamedColors) -> i32 {
        if pymol_index < 0 {
            return pymol_index as i32;
        }
        if let Some(&local) = self.cache.get(&pymol_index) {
            return local;
        }
        let local = if let Some((r, g, b)) = pymol_color_rgb(pymol_index) {
            let color_name = format!("_pse_{}", pymol_index);
            named_colors.register(&color_name, Color::new(r, g, b)) as i32
        } else {
            // Unknown PyMOL index — fall back to element coloring
            -1
        };
        self.cache.insert(pymol_index, local);
        local
    }
}

/// Convert a parsed PSE session into a live [`Session`].
pub fn pse_to_session(pse: &PseSession) -> Result<Session, SessionError> {
    let mut session = Session::new();

    // --- Settings ---
    // Skip PSE global settings for now — they can cause rendering issues
    // convert_settings(pse, &mut session);

    // --- Camera / View ---
    let view = convert_view(&pse.view);
    session.camera.set_view(view);

    // --- Named views ---
    for (name, view_arr) in &pse.view_dict {
        session.views.store_view(name, convert_view(view_arr));
    }

    // --- Custom colors ---
    for color in &pse.colors {
        session.named_colors.register(
            &color.name,
            Color::new(color.rgb[0], color.rgb[1], color.rgb[2]),
        );
    }

    // --- Color mapper ---
    let mut color_mapper = PseColorMapper::new();

    // --- Objects and Selections ---
    let mut object_names: Vec<Option<String>> = Vec::new();

    for entry in &pse.names {
        match entry {
            None => object_names.push(None),
            Some(PseNameEntry::Object(obj)) => {
                object_names.push(Some(obj.name.clone()));
                convert_object(obj, pse, &mut session, &mut color_mapper)?;
            }
            Some(PseNameEntry::Selection(sel)) => {
                object_names.push(Some(sel.name.clone()));
            }
        }
    }

    // --- Selections (after all objects are registered) ---
    for entry in &pse.names {
        if let Some(PseNameEntry::Selection(sel)) = entry {
            convert_selection(sel, &object_names, &mut session);
        }
    }

    Ok(session)
}

// =============================================================================
// Settings
// =============================================================================

#[allow(dead_code)]
fn convert_settings(pse: &PseSession, session: &mut Session) {
    let serialized: Vec<SerializedSetting> = pse
        .settings
        .iter()
        .filter_map(|s| {
            let value = convert_setting_value(&s.value)?;
            Some(SerializedSetting { id: s.id, value })
        })
        .collect();
    session.settings.from_session_list(&serialized);
}

fn convert_setting_value(pse_val: &PseSettingValue) -> Option<SettingValue> {
    Some(match pse_val {
        PseSettingValue::Bool(b) => SettingValue::Bool(*b),
        PseSettingValue::Int(i) => SettingValue::Int(*i as i32),
        PseSettingValue::Float(f) => SettingValue::Float(*f as f32),
        PseSettingValue::Float3(a) => SettingValue::Float3([a[0] as f32, a[1] as f32, a[2] as f32]),
        PseSettingValue::Color(c) => SettingValue::Color(*c as i32),
        PseSettingValue::String(s) => SettingValue::String(s.clone()),
    })
}

// =============================================================================
// View / Camera
// =============================================================================

/// Convert PSE 18-float view to [`SceneView`].
///
/// PSE layout:
/// - `[0..9]`   3×3 rotation (row-major)
/// - `[9..12]`  camera position
/// - `[12..15]` rotation origin
/// - `[15]`     front clip
/// - `[16]`     back clip
/// - `[17]`     orthographic flag
fn convert_view(view: &[f64; 18]) -> SceneView {
    let mut rotation = Mat4::new_identity();
    rotation.data[0] = view[0] as f32;
    rotation.data[1] = view[1] as f32;
    rotation.data[2] = view[2] as f32;
    rotation.data[4] = view[3] as f32;
    rotation.data[5] = view[4] as f32;
    rotation.data[6] = view[5] as f32;
    rotation.data[8] = view[6] as f32;
    rotation.data[9] = view[7] as f32;
    rotation.data[10] = view[8] as f32;

    SceneView {
        rotation,
        position: Vec3::new(view[9] as f32, view[10] as f32, view[11] as f32),
        origin: Vec3::new(view[12] as f32, view[13] as f32, view[14] as f32),
        clip_front: view[15] as f32,
        clip_back: view[16] as f32,
        fov: 14.0, // PSE doesn't store FOV in the view tuple; default
    }
}

// =============================================================================
// Objects
// =============================================================================

fn convert_object(
    obj: &PseObject,
    pse: &PseSession,
    session: &mut Session,
    color_mapper: &mut PseColorMapper,
) -> Result<(), SessionError> {
    match &obj.data {
        PseObjectData::Molecule(mol) => {
            let mol_obj = convert_molecule(&obj.name, obj, mol, pse, color_mapper, &mut session.named_colors)?;
            session.registry.add(mol_obj);
        }
        PseObjectData::Unsupported => {
            log::debug!("Skipping unsupported object '{}' (type {})", obj.name, obj.type_code);
        }
    }
    Ok(())
}

fn convert_molecule(
    name: &str,
    obj: &PseObject,
    pse_mol: &PseMolecule,
    pse: &PseSession,
    color_mapper: &mut PseColorMapper,
    named_colors: &mut NamedColors,
) -> Result<MoleculeObject, SessionError> {
    let mut mol = ObjectMolecule::with_capacity(name, pse_mol.atoms.len(), pse_mol.bonds.len());
    mol.discrete = pse_mol.discrete;

    // Symmetry
    if let Some(sym) = &pse_mol.symmetry {
        mol.symmetry = Some(convert_symmetry(sym));
    }

    // Atoms
    for pse_atom in &pse_mol.atoms {
        mol.add_atom(convert_atom(pse_atom, color_mapper, named_colors));
    }

    // Bonds
    for pse_bond in &pse_mol.bonds {
        let a1 = AtomIndex(pse_bond.index0 as u32);
        let a2 = AtomIndex(pse_bond.index1 as u32);
        let order = BondOrder::from_raw(pse_bond.order as i8);
        if let Ok(bi) = mol.add_bond(a1, a2, order) {
            if let Some(bond) = mol.get_bond_mut(bi) {
                bond.stereo = match pse_bond.stereo {
                    1 => BondStereo::Up,
                    2 => BondStereo::Down,
                    3 => BondStereo::Either,
                    _ => BondStereo::None,
                };
                if pse_bond.unique_id != 0 {
                    bond.unique_id = Some(pse_bond.unique_id as i32);
                    bond.has_setting = true;
                }
            }
        }
    }

    // Coordinate sets
    for pse_cs in &pse_mol.coord_sets {
        if let Some(cs) = pse_cs {
            mol.add_coord_set(convert_coord_set(cs, pse_mol.atoms.len()));
        }
    }

    // Unique settings
    let unique_list: Vec<_> = pse
        .unique_settings
        .iter()
        .filter_map(|us| {
            let value = convert_setting_value(&us.value)?;
            Some((us.unique_id as i32, us.setting_id, value))
        })
        .collect();
    mol.unique_settings.from_session_list(&unique_list);

    // Classify atoms
    mol.classify_atoms();

    // Compute object-level rep_mask as union of all per-atom reps
    let obj_rep_mask = pse_mol.atoms.iter().fold(0u32, |acc, a| acc | a.visible_reps);

    // Use normal constructor (which sets up rendering state),
    // then override per-atom visible_reps from PSE data.
    let mut mol_obj = MoleculeObject::new(mol);
    // Restore per-atom reps and remapped colors from PSE
    for (i, pse_atom) in pse_mol.atoms.iter().enumerate() {
        if let Some(atom) = mol_obj.molecule_mut().get_atom_mut(AtomIndex(i as u32)) {
            atom.repr.visible_reps = RepMask(pse_atom.visible_reps);
            atom.repr.colors.base = color_mapper.map(pse_atom.color, named_colors);
        }
    }
    let state = mol_obj.state_mut();
    state.enabled = obj.visible;
    state.visible_reps = RepMask(obj_rep_mask);
    state.color = color_index_from_pse(obj.color);

    Ok(mol_obj)
}

fn convert_atom(pse: &PseAtom, color_mapper: &mut PseColorMapper, named_colors: &mut NamedColors) -> Atom {
    let element = Element::from_symbol(&pse.elem).unwrap_or(Element::Unknown);
    let inscode = pse.resi.chars().last().filter(|c| c.is_ascii_alphabetic()).unwrap_or(' ');
    let residue = Arc::new(AtomResidue::new(
        ResidueKey::new(&pse.chain, &pse.resn, pse.resv as i32, inscode),
        pse.segi.clone(),
    ));
    let ss_type = match pse.ss.as_str() {
        "H" => SecondaryStructure::Helix,
        "S" => SecondaryStructure::Sheet,
        _ => SecondaryStructure::Loop,
    };

    let mut atom = Atom::default();
    atom.name = Arc::from(pse.name.as_str());
    atom.element = element;
    atom.residue = residue;
    atom.alt = pse.alt.chars().next().unwrap_or(' ');
    atom.b_factor = pse.b as f32;
    atom.occupancy = pse.q as f32;
    atom.vdw = pse.vdw as f32;
    atom.partial_charge = pse.partial_charge as f32;
    atom.formal_charge = pse.formal_charge as i8;
    atom.ss_type = ss_type;
    atom.id = pse.id as i32;
    atom.state.hetatm = pse.hetatm;
    atom.repr = AtomRepresentation {
        colors: AtomColors::with_base(color_mapper.map(pse.color, named_colors)),
        visible_reps: RepMask(pse.visible_reps),
        text_type: pse.text_type.clone(),
        label: pse.label.clone(),
        unique_id: if pse.unique_id != 0 { Some(pse.unique_id as i32) } else { None },
        has_setting: pse.unique_id != 0,
        ..Default::default()
    };
    atom
}

fn convert_coord_set(pse: &PseCoordSet, n_atoms: usize) -> CoordSet {
    if pse.idx_to_atm.is_empty() || pse.idx_to_atm.len() == n_atoms {
        // Check if idx_to_atm is identity (or empty) — no reordering needed
        let is_identity = pse.idx_to_atm.is_empty()
            || pse.idx_to_atm.iter().enumerate().all(|(i, &v)| v == i);

        if is_identity {
            let coords: Vec<f32> = pse.coords.iter().map(|&c| c as f32).collect();
            return CoordSet::from_coords(coords);
        }
    }

    // Reorder coordinates: PSE stores coords in index order,
    // but we need them in atom order. idx_to_atm[idx] = atom_index.
    let mut atom_coords = vec![0.0f32; n_atoms * 3];
    for (idx, &atm) in pse.idx_to_atm.iter().enumerate() {
        if atm < n_atoms && idx * 3 + 2 < pse.coords.len() {
            atom_coords[atm * 3] = pse.coords[idx * 3] as f32;
            atom_coords[atm * 3 + 1] = pse.coords[idx * 3 + 1] as f32;
            atom_coords[atm * 3 + 2] = pse.coords[idx * 3 + 2] as f32;
        }
    }
    CoordSet::from_coords(atom_coords)
}

fn convert_symmetry(pse: &PseSymmetry) -> Symmetry {
    Symmetry::new(
        &pse.space_group,
        [pse.cell[0] as f32, pse.cell[1] as f32, pse.cell[2] as f32],
        [pse.cell[3] as f32, pse.cell[4] as f32, pse.cell[5] as f32],
    )
}

// =============================================================================
// Selections
// =============================================================================

fn convert_selection(
    sel: &PseSelection,
    object_names: &[Option<String>],
    session: &mut Session,
) {
    // Build a union expression from referenced object names.
    // This is approximate — PSE selections store per-atom membership,
    // but SelectionManager uses text expressions.
    let obj_refs: Vec<&str> = sel
        .members
        .iter()
        .filter_map(|(idx, _)| object_names.get(*idx)?.as_deref())
        .collect();

    let expression = if obj_refs.is_empty() {
        "none".to_string()
    } else {
        obj_refs.join(" or ")
    };

    session.selections.define(&sel.name, &expression);
    session.selections.set_visible(&sel.name, sel.visible);
}

// =============================================================================
// Color helpers
// =============================================================================

fn color_index_from_pse(color: i64) -> ColorIndex {
    match color {
        -1 => ColorIndex::ByElement,
        -2 => ColorIndex::ByChain,
        -3 => ColorIndex::BySS,
        -4 => ColorIndex::ByBFactor,
        c if c >= 0 => ColorIndex::Named(c as u32),
        _ => ColorIndex::default(),
    }
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pse::*;
    use std::collections::HashMap;

    fn minimal_pse() -> PseSession {
        PseSession {
            version: 1830000,
            settings: vec![],
            view: [
                1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
                0.0, 0.0, -50.0,
                0.0, 0.0, 0.0,
                40.0, 100.0, 0.0,
            ],
            names: vec![],
            colors: vec![],
            view_dict: HashMap::new(),
            scene_order: vec![],
            unique_settings: vec![],
        }
    }

    #[test]
    fn test_empty_session() {
        let session = pse_to_session(&minimal_pse()).unwrap();
        assert_eq!(session.registry.len(), 0);
    }

    #[test]
    fn test_view_conversion() {
        let session = pse_to_session(&minimal_pse()).unwrap();
        let view = session.camera.current_view();
        assert!((view.position.z - (-50.0)).abs() < 0.001);
    }

    #[test]
    fn test_custom_colors() {
        let mut pse = minimal_pse();
        pse.colors.push(PseColor { name: "mycolor".into(), rgb: [0.5, 0.3, 0.1] });
        let session = pse_to_session(&pse).unwrap();
        let (_, c) = session.named_colors.get_by_name("mycolor").unwrap();
        assert!((c.r - 0.5).abs() < 0.001);
    }

    #[test]
    fn test_molecule_conversion() {
        let mut pse = minimal_pse();
        pse.names.push(Some(PseNameEntry::Object(PseObject {
            name: "mol1".into(),
            type_code: 1,
            visible: true,
            rep_mask: 0x01,
            color: -1,
            data: PseObjectData::Molecule(PseMolecule {
                atoms: vec![PseAtom {
                    resv: 1, chain: "A".into(), alt: String::new(), resi: "1".into(),
                    segi: String::new(), resn: "ALA".into(), name: "CA".into(),
                    elem: "C".into(), text_type: String::new(), label: String::new(),
                    ss: "H".into(), b: 20.0, q: 1.0, vdw: 1.7,
                    partial_charge: 0.0, formal_charge: 0, color: -1, id: 1, unique_id: 0,
                    visible_reps: 0x20, hetatm: false,
                }],
                bonds: vec![],
                coord_sets: vec![Some(PseCoordSet { n_atom: 1, coords: vec![1.0, 2.0, 3.0], idx_to_atm: vec![] })],
                discrete: false, n_discrete: 0, symmetry: None,
            }),
        })));

        let session = pse_to_session(&pse).unwrap();
        assert_eq!(session.registry.len(), 1);
        let mol_obj = session.registry.get_molecule("mol1").unwrap();
        assert_eq!(mol_obj.molecule().atom_count(), 1);
        assert_eq!(&*mol_obj.molecule().get_atom(AtomIndex(0)).unwrap().name, "CA");
    }

    #[test]
    fn test_selection_conversion() {
        let mut pse = minimal_pse();
        pse.names.push(Some(PseNameEntry::Object(PseObject {
            name: "mol1".into(), type_code: 1, visible: true, rep_mask: 1, color: -1,
            data: PseObjectData::Molecule(PseMolecule {
                atoms: vec![], bonds: vec![], coord_sets: vec![],
                discrete: false, n_discrete: 0, symmetry: None,
            }),
        })));
        pse.names.push(Some(PseNameEntry::Selection(PseSelection {
            name: "sele".into(), visible: true, members: vec![(0, vec![0])],
        })));

        let session = pse_to_session(&pse).unwrap();
        assert!(session.selections.contains("sele"));
    }

    #[test]
    fn test_named_views() {
        let mut pse = minimal_pse();
        pse.view_dict.insert("front".into(), [
            1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
            0.0, 0.0, -30.0, 0.0, 0.0, 0.0, 20.0, 80.0, 0.0,
        ]);
        let session = pse_to_session(&pse).unwrap();
        assert!(session.views.get("front").is_some());
    }

    #[test]
    fn test_color_index_mapping() {
        assert!(matches!(color_index_from_pse(-1), ColorIndex::ByElement));
        assert!(matches!(color_index_from_pse(-2), ColorIndex::ByChain));
        assert!(matches!(color_index_from_pse(5), ColorIndex::Named(5)));
    }
}
