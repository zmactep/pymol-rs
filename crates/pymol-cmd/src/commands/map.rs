//! Map visualization commands: isomesh, isodot

use pymol_mol::AtomIndex;
use pymol_scene::MapDisplayMode;

use crate::args::ParsedCommand;
use crate::command::{ArgHint, Command, CommandContext, CommandRegistry, ViewerLike};
use crate::error::{CmdError, CmdResult};

use super::selecting::evaluate_selection;

/// Register map commands
pub fn register(registry: &mut CommandRegistry) {
    registry.register(IsomeshCommand);
    registry.register(IsodotCommand);
    registry.register(IsosurfaceCommand);
}

/// Shared implementation for isomesh/isodot/isosurface commands
fn execute_isocontour(
    ctx: &mut CommandContext<'_, '_, dyn ViewerLike + '_>,
    args: &ParsedCommand,
    mode: MapDisplayMode,
) -> CmdResult {
    // isomesh name, map, level [, selection [, buffer [, carve ]]]
    let name = args
        .str_arg(0, "name")
        .ok_or_else(|| CmdError::MissingArgument("name".to_string()))?;

    let map_name = args
        .str_arg(1, "map")
        .ok_or_else(|| CmdError::MissingArgument("map".to_string()))?;

    let level = args.float_arg_or(2, "level", 1.0) as f32;

    let selection = args.str_arg(3, "selection");

    let _buffer = args.float_arg_or(4, "buffer", 0.0);

    let carve = args.float_arg_or(5, "carve", 0.0) as f32;

    // Get the source map's grid data
    let grid = {
        let map_obj = ctx
            .viewer
            .objects()
            .get_map(map_name)
            .ok_or_else(|| CmdError::execution(format!("'{}' is not a map object", map_name)))?;

        map_obj
            .grid()
            .ok_or_else(|| CmdError::execution(format!("Map '{}' has no data", map_name)))?
            .clone()
    };

    // Create a new MapObject for the visualization
    let mut vis_obj = pymol_scene::MapObject::new(name, grid);
    vis_obj.set_level(level);
    vis_obj.set_display_mode(mode);

    // If carve is specified with a selection, get atom positions
    if carve > 0.0 {
        if let Some(sel_str) = selection {
            let results = evaluate_selection(ctx.viewer, sel_str)?;
            let mut positions = Vec::new();
            for (obj_name, sel) in &results {
                if let Some(mol_obj) = ctx.viewer.objects().get_molecule(obj_name) {
                    let mol = mol_obj.molecule();
                    if let Some(cs) = mol.current_coord_set() {
                        for (i, _) in mol.atoms().enumerate() {
                            let atom_idx = AtomIndex::new(i as u32);
                            if sel.contains(atom_idx) {
                                if let Some(pos) = cs.get_atom_coord(atom_idx) {
                                    positions.push([pos.x, pos.y, pos.z]);
                                }
                            }
                        }
                    }
                }
            }
            if !positions.is_empty() {
                vis_obj.set_carve_radius(carve);
                vis_obj.set_carve_positions(positions);
            }
        }
    }

    ctx.viewer.objects_mut().add(vis_obj);

    let mode_name = match mode {
        MapDisplayMode::Isomesh => "isomesh",
        MapDisplayMode::Isodot => "isodot",
        MapDisplayMode::Isosurface => "isosurface",
        MapDisplayMode::Volume => "volume",
    };

    ctx.print(&format!(
        " Created {} \"{}\" from map \"{}\" at level {:.2}",
        mode_name, name, map_name, level
    ));

    Ok(())
}

// ============================================================================
// isomesh command
// ============================================================================

struct IsomeshCommand;

impl Command for IsomeshCommand {
    fn name(&self) -> &str {
        "isomesh"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "isomesh" creates a mesh representation of an electron density map
    at a specified contour level.

USAGE

    isomesh name, map, level [, selection [, buffer [, carve ]]]

ARGUMENTS

    name = string: name for the new mesh object
    map = string: name of the map object
    level = float: contour level in sigma (default: 1.0)
    selection = string: atoms around which to show density
    buffer = float: buffer distance around selection (default: 0)
    carve = float: carve radius - hide density beyond this distance (default: 0)

EXAMPLES

    load density.ccp4
    isomesh mesh1, density, 1.0
    isomesh mesh2, density, 1.5, organic, carve=1.6
"#
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::Object, ArgHint::Object]
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        execute_isocontour(ctx, args, MapDisplayMode::Isomesh)
    }
}

// ============================================================================
// isodot command
// ============================================================================

struct IsodotCommand;

impl Command for IsodotCommand {
    fn name(&self) -> &str {
        "isodot"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "isodot" creates a dot representation of an electron density map
    at a specified contour level.

USAGE

    isodot name, map, level [, selection [, buffer [, carve ]]]

ARGUMENTS

    name = string: name for the new dot object
    map = string: name of the map object
    level = float: contour level in sigma (default: 1.0)
    selection = string: atoms around which to show density
    buffer = float: buffer distance around selection (default: 0)
    carve = float: carve radius (default: 0)

EXAMPLES

    load density.ccp4
    isodot dots1, density, 1.0
    isodot dots2, density, 2.0, chain A, carve=2.0
"#
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::Object, ArgHint::Object]
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        execute_isocontour(ctx, args, MapDisplayMode::Isodot)
    }
}

// ============================================================================
// isosurface command
// ============================================================================

struct IsosurfaceCommand;

impl Command for IsosurfaceCommand {
    fn name(&self) -> &str {
        "isosurface"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "isosurface" creates a solid surface representation of an electron
    density map at a specified contour level.

USAGE

    isosurface name, map, level [, selection [, buffer [, carve ]]]

ARGUMENTS

    name = string: name for the new surface object
    map = string: name of the map object
    level = float: contour level in sigma (default: 1.0)
    selection = string: atoms around which to show density
    buffer = float: buffer distance around selection (default: 0)
    carve = float: carve radius (default: 0)

EXAMPLES

    load density.ccp4
    isosurface surf1, density, 1.0
"#
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::Object, ArgHint::Object]
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        execute_isocontour(ctx, args, MapDisplayMode::Isosurface)
    }
}
