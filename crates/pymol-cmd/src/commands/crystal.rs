//! Crystal symmetry commands: symexp

use std::sync::Arc;

use lin_alg::f32::{Mat4, Vec3};

use pymol_algos::crystal::CrystalCell;
use pymol_algos::linalg::{is_identity_mat4, left_multiply_mat4, transform_mat4};
use pymol_mol::spatial::SpatialGrid;
use pymol_mol::RepMask;
use pymol_scene::MoleculeObject;

use crate::args::ParsedCommand;
use crate::command::{ArgHint, Command, CommandContext, CommandRegistry, ViewerLike};
use crate::commands::selecting::evaluate_selection;
use crate::error::{CmdError, CmdResult};

pub fn register(registry: &mut CommandRegistry) {
    registry.register(SymExpCommand);
}

// ============================================================================
// symexp command
// ============================================================================

struct SymExpCommand;

impl Command for SymExpCommand {
    fn name(&self) -> &str {
        "symexp"
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[
            ArgHint::None,      // prefix
            ArgHint::Object,    // object
            ArgHint::Selection, // selection
            ArgHint::None,      // cutoff
        ]
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "symexp" generates symmetry-related copies of an object within a
    specified distance cutoff. The new objects are named using the
    given prefix followed by symmetry operation and translation codes.

USAGE

    symexp prefix, object, selection, cutoff [, segi [, quiet]]

ARGUMENTS

    prefix = string: prefix for new object names
    object = string: source object with crystallographic symmetry
    selection = string: atom selection defining the region of interest
    cutoff = float: distance cutoff in Angstroms (default: 10.0)
    segi = 0/1: store symmetry info in segment identifier (default: 0)
    quiet = 0/1: suppress output messages (default: 1)

EXAMPLES

    symexp sym, 1oky, chain A, 20.0
    symexp s, myprotein, all, 15.0, segi=1
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        // Parse arguments
        let prefix = args
            .get_str(0)
            .or_else(|| args.get_named_str("prefix"))
            .ok_or_else(|| CmdError::MissingArgument("prefix".to_string()))?;
        let object = args
            .get_str(1)
            .or_else(|| args.get_named_str("object"))
            .ok_or_else(|| CmdError::MissingArgument("object".to_string()))?;
        let selection = args
            .get_str(2)
            .or_else(|| args.get_named_str("selection"))
            .ok_or_else(|| CmdError::MissingArgument("selection".to_string()))?;
        let cutoff = args
            .get_float(3)
            .or_else(|| args.get_named_float("cutoff"))
            .unwrap_or(10.0) as f32;
        let segi = args
            .get_int(4)
            .or_else(|| args.get_named_int("segi"))
            .unwrap_or(0)
            != 0;
        let quiet = args
            .get_int(5)
            .or_else(|| args.get_named_int("quiet"))
            .map(|v| v != 0)
            .unwrap_or(ctx.quiet);

        // Get source molecule data (clone to release borrow on ctx.viewer)
        let (src_mol, source_reps, symop_matrices, cell, state_centers) = {
            let mol_obj = ctx
                .viewer
                .objects()
                .get_molecule(object)
                .ok_or_else(|| CmdError::ObjectNotFound(object.to_string()))?;

            let mol = mol_obj.molecule();

            // Find symmetry: try molecule-level first, then first coord set with symmetry
            let symmetry = mol
                .symmetry
                .as_ref()
                .or_else(|| {
                    (0..mol.state_count()).find_map(|s| {
                        mol.get_coord_set(s)
                            .and_then(|cs| cs.symmetry.as_ref())
                    })
                })
                .ok_or_else(|| CmdError::execution("No symmetry loaded!"))?;

            // Get space group operations
            let symop_matrices = pymol_algos::space_groups::get_symops(&symmetry.space_group)
                .ok_or_else(|| {
                    CmdError::execution(format!(
                        "Unknown space group: '{}'",
                        symmetry.space_group
                    ))
                })?;

            let cell = CrystalCell::new(symmetry.cell_lengths, symmetry.cell_angles);

            let state_centers: Vec<Vec3> = (0..mol.state_count())
                .map(|s| {
                    mol.get_coord_set(s)
                        .map(|cs| cs.center())
                        .unwrap_or(Vec3::new(0.0, 0.0, 0.0))
                })
                .collect();

            let source_reps = mol_obj.visible_reps();

            (mol.clone(), source_reps, symop_matrices, cell, state_centers)
        };

        // Evaluate selection and collect all selected atom coordinates
        let results = evaluate_selection(ctx.viewer, selection)?;
        let (sel_center, sel_coords) = collect_selection_coords(ctx.viewer, &results);

        if sel_coords.is_empty() {
            return Err(CmdError::execution("No atoms in selection!"));
        }

        // Convert selection center to fractional coordinates
        let tc = cell.to_fractional(sel_center);

        // Build spatial hash of selection atoms for distance checks
        let mut grid = SpatialGrid::with_capacity(cutoff, sel_coords.len());
        for (i, coord) in sel_coords.iter().enumerate() {
            grid.insert(*coord, i);
        }

        let cutoff_sq = cutoff * cutoff;

        // Generate symmetry mates
        if !quiet {
            ctx.print(" SymExp: Generating symmetry mates...");
        }

        let mut objects_to_add: Vec<(String, pymol_mol::ObjectMolecule, RepMask)> = Vec::new();

        for x in -1..=1_i32 {
            for y in -1..=1_i32 {
                for z in -1..=1_i32 {
                    for (a, symop) in symop_matrices.iter().enumerate() {
                        // Compute per-state matrices and check if any state has atoms in range
                        let matrices =
                            compute_per_state_matrices(&cell, symop, &state_centers, tc, x, y, z);

                        // Skip if all matrices are identity (this is the original molecule)
                        if matrices.iter().all(|m| is_identity_mat4(m)) {
                            continue;
                        }

                        // Check distance: does any transformed atom fall within cutoff?
                        let keep = check_distance(
                            &src_mol,
                            &matrices,
                            &grid,
                            &sel_coords,
                            cutoff_sq,
                        );

                        if !keep {
                            continue;
                        }

                        // Clone and transform
                        let mut new_mol = src_mol.clone();
                        for (state_idx, mat) in matrices.iter().enumerate() {
                            if !is_identity_mat4(mat) {
                                new_mol.transform(state_idx, mat);
                            }
                        }

                        // Set segment identifiers if requested
                        if segi {
                            let label = make_segi_label(a, x, y, z);
                            for atom in new_mol.atoms_mut() {
                                let res = Arc::make_mut(&mut atom.residue);
                                res.segi = label.clone();
                            }
                        }

                        // Generate object name matching PyMOL convention
                        let name = format!(
                            "{}{:02}{:02}{:02}{:02}",
                            prefix,
                            a,
                            x + 1,
                            y + 1,
                            z + 1
                        );

                        objects_to_add.push((name, new_mol, source_reps));
                    }
                }
            }
        }

        // Add all new objects to the scene
        let count = objects_to_add.len();
        for (name, mol, reps) in objects_to_add {
            let mut obj = MoleculeObject::from_raw_with_name(mol, &name);
            obj.set_visible_reps(reps);
            ctx.viewer.objects_mut().add(obj);
        }

        ctx.viewer.request_redraw();

        if !quiet {
            ctx.print(&format!(
                " SymExp: Created {} symmetry mate objects.",
                count
            ));
        }

        Ok(())
    }
}

/// Collect all selected atom coordinates across all objects and states,
/// and compute their center of mass.
fn collect_selection_coords(
    viewer: &dyn ViewerLike,
    results: &[(String, pymol_select::SelectionResult)],
) -> (Vec3, Vec<Vec3>) {
    let mut coords = Vec::new();
    let mut sum = Vec3::new(0.0, 0.0, 0.0);

    for (obj_name, selection) in results {
        if selection.count() == 0 {
            continue;
        }
        let mol_obj = match viewer.objects().get_molecule(obj_name) {
            Some(m) => m,
            None => continue,
        };
        let mol = mol_obj.molecule();

        // Collect from all states (matching PyMOL's behavior: OMOP_SUMC + OMOP_VERT)
        for state in 0..mol.state_count() {
            if let Some(cs) = mol.get_coord_set(state) {
                for idx in selection.indices() {
                    if let Some(coord) = cs.get_atom_coord(idx) {
                        sum += coord;
                        coords.push(coord);
                    }
                }
            }
        }
    }

    let center = if coords.is_empty() {
        Vec3::new(0.0, 0.0, 0.0)
    } else {
        sum * (1.0 / coords.len() as f32)
    };

    (center, coords)
}

/// Compute the transformation matrix for each state of the molecule.
///
/// Following PyMOL's algorithm in `ExecutiveSymExp`:
/// 1. mat = realToFrac (3x3→4x4)
/// 2. left_multiply by symOp
/// 3. Transform state center through mat → compute rounding shift
/// 4. Add unit cell translation (x, y, z)
/// 5. Build shift matrix and left_multiply into mat
/// 6. left_multiply fracToReal into mat
fn compute_per_state_matrices(
    cell: &CrystalCell,
    symop: &[f32; 16],
    state_centers: &[Vec3],
    tc: Vec3,
    x: i32,
    y: i32,
    z: i32,
) -> Vec<Mat4> {
    let symop_mat = Mat4 { data: *symop };
    let r2f_4x4 = cell.real_to_frac_4x4();
    let f2r_4x4 = cell.frac_to_real_4x4();

    state_centers
        .iter()
        .map(|center| {
            // Step 1-2: mat = symOp * realToFrac
            let mut mat = left_multiply_mat4(&symop_mat, &r2f_4x4);

            // Step 3: Transform state center to get fractional position
            let ts = transform_mat4(&mat, *center);

            // Compute rounding shift to bring into same unit cell as selection center
            let shift = [
                (tc.x - ts.x).round() + x as f32,
                (tc.y - ts.y).round() + y as f32,
                (tc.z - ts.z).round() + z as f32,
            ];

            // Step 4-5: Build shift matrix and left-multiply
            let shift_mat = Mat4 {
                data: [
                    1.0, 0.0, 0.0, shift[0],
                    0.0, 1.0, 0.0, shift[1],
                    0.0, 0.0, 1.0, shift[2],
                    0.0, 0.0, 0.0, 1.0,
                ],
            };
            mat = left_multiply_mat4(&shift_mat, &mat);

            // Step 6: Left-multiply fracToReal
            left_multiply_mat4(&f2r_4x4, &mat)
        })
        .collect()
}

/// Check if any atom from the source molecule, when transformed, falls within
/// cutoff distance of any selection atom.
fn check_distance(
    mol: &pymol_mol::ObjectMolecule,
    matrices: &[Mat4],
    grid: &SpatialGrid,
    sel_coords: &[Vec3],
    cutoff_sq: f32,
) -> bool {
    let mut neighbors = Vec::new();

    for (state_idx, mat) in matrices.iter().enumerate() {
        let cs = match mol.get_coord_set(state_idx) {
            Some(cs) => cs,
            None => continue,
        };

        for coord in cs.iter() {
            // Transform the coordinate by the symmetry matrix
            let transformed = transform_mat4(mat, coord);

            // Check spatial grid for nearby selection atoms
            grid.query_neighbors(transformed, &mut neighbors);
            for &idx in &neighbors {
                let sel = sel_coords[idx];
                let dx = transformed.x - sel.x;
                let dy = transformed.y - sel.y;
                let dz = transformed.z - sel.z;
                if dx * dx + dy * dy + dz * dz <= cutoff_sq {
                    return true;
                }
            }
        }
    }

    false
}

/// Generate segment identifier label for symmetry operation encoding.
///
/// Follows PyMOL's `make_symexp_segi_label`:
/// - Char 0: symop index (A=0, B=1, ..., Z=25)
/// - Chars 1-3: x,y,z translations encoded as digit (0=−1, 1=0, 2=+1)
fn make_segi_label(symop_idx: usize, x: i32, y: i32, z: i32) -> String {
    let op_char = if symop_idx < 26 {
        (b'A' + symop_idx as u8) as char
    } else if symop_idx < 36 {
        (b'0' + (symop_idx - 26) as u8) as char
    } else {
        (b'a' + (symop_idx - 36) as u8) as char
    };
    format!("{}{}{}{}", op_char, x + 1, y + 1, z + 1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_make_segi_label() {
        assert_eq!(make_segi_label(0, 0, 0, 0), "A111");
        assert_eq!(make_segi_label(3, -1, 0, 1), "D012");
        assert_eq!(make_segi_label(25, 1, 1, 1), "Z222");
    }
}
