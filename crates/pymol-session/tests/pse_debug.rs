use pymol_session::{load_pse, convert::pse_to_session};

#[test]
fn debug_1fsd_pse_atoms() {
    let pse = load_pse(std::path::Path::new("../../_tests/1fsd.pse")).unwrap();
    let session = pse_to_session(&pse).unwrap();

    println!("Named colors count: {}", session.named_colors.len());
    for i in 50..session.named_colors.len().min(70) {
        if let Some(c) = session.named_colors.get_by_index(i as u32) {
            println!("  color[{i}]: ({:.2}, {:.2}, {:.2})", c.r, c.g, c.b);
        }
    }

    for name in session.registry.names() {
        println!("Object: {name}");
        if let Some(mol) = session.registry.get_molecule(name) {
            let mol_data = mol.molecule();
            let atom_vec: Vec<_> = mol_data.atoms().collect();
            println!("  Atoms: {}", atom_vec.len());
            for (i, atom) in atom_vec.iter().take(10).enumerate() {
                println!("  Atom[{i}]: name='{}' elem='{}' colors.base={} hetatm={} flags={:?}",
                    atom.name, atom.element.symbol(), atom.repr.colors.base,
                    atom.state.hetatm, atom.state.flags);
            }
        }
    }
}
