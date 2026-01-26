//! Simple molecular viewer example
//!
//! Demonstrates basic usage of pymol-scene:
//! - Creating a viewer
//! - Loading a molecule (programmatically or from file)
//! - Setting up representations
//! - Running the interactive viewer
//!
//! ## Controls
//!
//! ### Camera
//! - Left mouse drag: Rotate
//! - Middle mouse drag: Pan
//! - Right mouse drag: Zoom
//! - Scroll wheel: Zoom
//! - R key: Reset view
//! - O key: Toggle orthographic/perspective
//! - Escape: Exit
//!
//! ### Representations (toggle with number keys)
//! - 1: Lines (wireframe bonds)
//! - 2: Sticks (cylinder bonds)
//! - 3: Spheres (VdW spheres)
//! - 4: Cartoon (secondary structure with arrows)
//! - 5: Surface (molecular surface)
//! - 6: Mesh (mesh surface)
//! - 7: Dots (dot surface)
//! - 8: Ribbon (secondary structure without arrows)
//! - H: Hide all representations
//! - A: Show default (lines + sticks)
//!
//! ### Surface Quality
//! - +/=: Increase quality (finer mesh, slower)
//! - -: Decrease quality (coarser mesh, faster)

use std::path::Path;

use lin_alg::f32::Vec3;
use pymol_mol::{Atom, BondOrder, CoordSet, Element, ObjectMolecule, RepMask};
use pymol_scene::{run, KeyBinding, KeyCode, MoleculeObject, Viewer};

fn main() {
    // Initialize logging
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    // Create the viewer
    let mut viewer = Viewer::new();

    // Register default key bindings
    // Camera controls
    viewer.bind_key(KeyCode::KeyR, |v| v.reset_view());
    viewer.bind_key(KeyCode::KeyO, |v| {
        v.camera_mut().toggle_projection();
        v.request_redraw();
    });

    // Representation hotkeys (number keys)
    viewer.bind_key(KeyCode::Digit1, |v| v.toggle_representation(RepMask::LINES));
    viewer.bind_key(KeyCode::Digit2, |v| v.toggle_representation(RepMask::STICKS));
    viewer.bind_key(KeyCode::Digit3, |v| v.toggle_representation(RepMask::SPHERES));
    viewer.bind_key(KeyCode::Digit4, |v| v.toggle_representation(RepMask::CARTOON));
    viewer.bind_key(KeyCode::Digit5, |v| v.toggle_representation(RepMask::SURFACE));
    viewer.bind_key(KeyCode::Digit6, |v| v.toggle_representation(RepMask::MESH));
    viewer.bind_key(KeyCode::Digit7, |v| v.toggle_representation(RepMask::DOTS));
    viewer.bind_key(KeyCode::Digit8, |v| v.toggle_representation(RepMask::RIBBON));

    // Representation management
    viewer.bind_key(KeyCode::KeyH, |v| v.hide_all_representations());
    viewer.bind_key(KeyCode::KeyA, |v| v.show_default_representations());

    // Surface quality controls
    viewer.bind_key(KeyCode::Equal, |v| v.increase_surface_quality()); // + key
    viewer.bind_key(KeyCode::Minus, |v| v.decrease_surface_quality()); // - key

    // Example: Key combinations with modifiers
    // Ctrl+R: Reset view (alternative binding)
    viewer.bind_key(KeyBinding::new(KeyCode::KeyR).ctrl(), |v| {
        log::info!("Ctrl+R: Reset view");
        v.reset_view();
    });

    // Check for file argument
    let args: Vec<String> = std::env::args().collect();
    
    // Find the first argument that looks like a file path (has an extension)
    let file_arg = args.iter().skip(1).find(|arg| {
        let path = Path::new(arg);
        // Check if it has a recognized extension or exists as a file
        path.extension().is_some() || path.exists()
    });

    if let Some(file_path) = file_arg {
        // Load molecule from file
        let path = Path::new(file_path);
        match pymol_io::read_file(path) {
            Ok(mut mol) => {
                log::info!("Loaded {} atoms from {}", mol.atom_count(), path.display());

                // Generate bonds if none were present (e.g., for CIF files)
                if mol.bond_count() == 0 {
                    log::info!("No bonds found, generating from distances...");
                    mol.generate_bonds(0.6);
                    log::info!("Generated {} bonds", mol.bond_count());
                }

                // Create molecule object and show as sticks
                let mut obj = MoleculeObject::new(mol);
                obj.show(RepMask::LINES);
                obj.show(RepMask::STICKS);
                viewer.objects_mut().add(obj);
            }
            Err(e) => {
                log::error!("Failed to load {}: {}", path.display(), e);
                log::info!("Using sample peptide instead");

                let peptide = create_simple_peptide();
                let mut obj = MoleculeObject::new(peptide);
                obj.show(RepMask::STICKS);
                obj.show(RepMask::LINES);
                viewer.objects_mut().add(obj);
            }
        }
    } else {
        // Use sample molecule
        log::info!("No file provided, using sample peptide");
        log::info!("Usage: viewer [molecule_file.pdb]");

        // Create a simple peptide backbone
        let peptide = create_simple_peptide();
        log::info!("Created peptide with {} atoms, {} bonds", peptide.atom_count(), peptide.bond_count());
        
        let mut obj = MoleculeObject::new(peptide);
        obj.show(RepMask::STICKS);
        obj.show(RepMask::LINES);
        viewer.objects_mut().add(obj);
    }

    // Set a nice dark blue background
    viewer.set_background_color(0.0, 0.0, 0.1);

    // Center the view on all objects
    viewer.center_all();

    // Run the viewer
    log::info!("Starting viewer...");
    log::info!("Controls:");
    log::info!("  Camera:");
    log::info!("    Left drag: Rotate");
    log::info!("    Middle drag: Pan");
    log::info!("    Right drag / Scroll: Zoom");
    log::info!("    R: Reset view");
    log::info!("    O: Toggle ortho/perspective");
    log::info!("  Representations:");
    log::info!("    1: Lines    2: Sticks   3: Spheres  4: Cartoon");
    log::info!("    5: Surface  6: Mesh     7: Dots     8: Ribbon");
    log::info!("    H: Hide all  A: Default (lines+sticks)");
    log::info!("  Surface Quality:");
    log::info!("    +/-: Increase/Decrease quality");

    if let Err(e) = run(viewer) {
        log::error!("Viewer error: {}", e);
    }
}

/// Create a simple water molecule (example of programmatic molecule creation)
#[allow(dead_code)]
fn create_water_molecule() -> ObjectMolecule {
    let mut mol = ObjectMolecule::new("water");

    // Add atoms
    let o = mol.add_atom(Atom::new("O", Element::Oxygen));
    let h1 = mol.add_atom(Atom::new("H1", Element::Hydrogen));
    let h2 = mol.add_atom(Atom::new("H2", Element::Hydrogen));

    // Add bonds
    mol.add_bond(o, h1, BondOrder::Single).unwrap();
    mol.add_bond(o, h2, BondOrder::Single).unwrap();

    // Add coordinates (water geometry)
    let coords = vec![
        Vec3::new(0.0, 0.0, 0.0),      // O at origin
        Vec3::new(0.96, 0.0, 0.0),     // H1 along x-axis
        Vec3::new(-0.24, 0.93, 0.0),   // H2 at ~104.5 degree angle
    ];
    mol.add_coord_set(CoordSet::from_vec3(&coords));

    mol
}

/// Create a simple peptide backbone (Ala-Gly-Ala)
fn create_simple_peptide() -> ObjectMolecule {
    let mut mol = ObjectMolecule::new("peptide");

    // Define atoms with coordinates
    struct AtomDef {
        name: &'static str,
        element: Element,
        resn: &'static str,
        resv: i32,
        coord: Vec3,
    }

    let atoms = vec![
        // Residue 1: ALA
        AtomDef { name: "N", element: Element::Nitrogen, resn: "ALA", resv: 1, coord: Vec3::new(0.0, 0.0, 0.0) },
        AtomDef { name: "CA", element: Element::Carbon, resn: "ALA", resv: 1, coord: Vec3::new(1.45, 0.0, 0.0) },
        AtomDef { name: "C", element: Element::Carbon, resn: "ALA", resv: 1, coord: Vec3::new(2.0, 1.4, 0.0) },
        AtomDef { name: "O", element: Element::Oxygen, resn: "ALA", resv: 1, coord: Vec3::new(1.4, 2.4, 0.0) },
        AtomDef { name: "CB", element: Element::Carbon, resn: "ALA", resv: 1, coord: Vec3::new(1.9, -0.7, 1.2) },
        // Residue 2: GLY
        AtomDef { name: "N", element: Element::Nitrogen, resn: "GLY", resv: 2, coord: Vec3::new(3.3, 1.5, 0.0) },
        AtomDef { name: "CA", element: Element::Carbon, resn: "GLY", resv: 2, coord: Vec3::new(4.0, 2.8, 0.0) },
        AtomDef { name: "C", element: Element::Carbon, resn: "GLY", resv: 2, coord: Vec3::new(5.5, 2.7, 0.0) },
        AtomDef { name: "O", element: Element::Oxygen, resn: "GLY", resv: 2, coord: Vec3::new(6.1, 1.6, 0.0) },
        // Residue 3: ALA
        AtomDef { name: "N", element: Element::Nitrogen, resn: "ALA", resv: 3, coord: Vec3::new(6.1, 3.9, 0.0) },
        AtomDef { name: "CA", element: Element::Carbon, resn: "ALA", resv: 3, coord: Vec3::new(7.5, 4.1, 0.0) },
        AtomDef { name: "C", element: Element::Carbon, resn: "ALA", resv: 3, coord: Vec3::new(8.0, 5.5, 0.0) },
        AtomDef { name: "O", element: Element::Oxygen, resn: "ALA", resv: 3, coord: Vec3::new(7.2, 6.4, 0.0) },
        AtomDef { name: "CB", element: Element::Carbon, resn: "ALA", resv: 3, coord: Vec3::new(8.3, 3.2, 1.0) },
    ];

    // Collect coordinates
    let coords: Vec<Vec3> = atoms.iter().map(|a| a.coord.clone()).collect();

    // Add atoms
    let mut indices = Vec::new();
    for def in &atoms {
        let mut atom = Atom::new(def.name, def.element);
        atom.resn = def.resn.to_string();
        atom.resv = def.resv;
        atom.chain = "A".to_string();
        atom.color = -1; // Color by element
        indices.push(mol.add_atom(atom));
    }

    // Add coordinate set
    mol.add_coord_set(CoordSet::from_vec3(&coords));

    // Add backbone bonds (using indices)
    // Residue 1: N(0)-CA(1), CA(1)-C(2), C(2)=O(3), CA(1)-CB(4), C(2)-N(5)
    mol.add_bond(indices[0], indices[1], BondOrder::Single).unwrap();
    mol.add_bond(indices[1], indices[2], BondOrder::Single).unwrap();
    mol.add_bond(indices[2], indices[3], BondOrder::Double).unwrap();
    mol.add_bond(indices[1], indices[4], BondOrder::Single).unwrap();
    mol.add_bond(indices[2], indices[5], BondOrder::Single).unwrap();  // Peptide bond

    // Residue 2: N(5)-CA(6), CA(6)-C(7), C(7)=O(8), C(7)-N(9)
    mol.add_bond(indices[5], indices[6], BondOrder::Single).unwrap();
    mol.add_bond(indices[6], indices[7], BondOrder::Single).unwrap();
    mol.add_bond(indices[7], indices[8], BondOrder::Double).unwrap();
    mol.add_bond(indices[7], indices[9], BondOrder::Single).unwrap();  // Peptide bond

    // Residue 3: N(9)-CA(10), CA(10)-C(11), C(11)=O(12), CA(10)-CB(13)
    mol.add_bond(indices[9], indices[10], BondOrder::Single).unwrap();
    mol.add_bond(indices[10], indices[11], BondOrder::Single).unwrap();
    mol.add_bond(indices[11], indices[12], BondOrder::Double).unwrap();
    mol.add_bond(indices[10], indices[13], BondOrder::Single).unwrap();

    // Classify atoms as protein so cartoon representation works
    mol.classify_atoms();

    mol
}
