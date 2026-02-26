//! Measurement commands: distance, angle, dihedral
//!
//! Create measurement objects that display dashed lines between atoms
//! with labels showing the measured value.

use lin_alg::f32::Vec3;

use crate::args::ParsedCommand;
use crate::command::{ArgHint, Command, CommandContext, CommandRegistry, ViewerLike};
use crate::commands::selecting::evaluate_selection;
use crate::error::{CmdError, CmdResult};

use pymol_scene::{Measurement, MeasurementObject};

/// Default measurement line color (yellow)
const DEFAULT_COLOR: [f32; 4] = [1.0, 1.0, 0.0, 1.0];

pub fn register(registry: &mut CommandRegistry) {
    registry.register(DistanceCommand);
    registry.register(AngleCommand);
    registry.register(DihedralCommand);
}

/// Resolve the first atom position from a selection expression.
fn resolve_atom_position(viewer: &dyn ViewerLike, selection: &str) -> CmdResult<Vec3> {
    let results = evaluate_selection(viewer, selection)?;
    for (obj_name, selected) in &results {
        if let Some(mol_obj) = viewer.objects().get_molecule(obj_name) {
            if let Some(idx) = selected.indices().next() {
                if let Some(coord) =
                    mol_obj
                        .molecule()
                        .get_coord_set(mol_obj.display_state())
                        .and_then(|cs| cs.get_atom_coord(idx))
                {
                    return Ok(coord);
                }
            }
        }
    }
    Err(CmdError::selection(format!(
        "no atoms found in selection '{}'",
        selection
    )))
}

/// Add a measurement to an existing or new measurement object.
fn add_measurement_to_scene(
    viewer: &mut dyn ViewerLike,
    name: &str,
    measurement: Measurement,
) {
    if let Some(meas_obj) = viewer.objects_mut().get_measurement_mut(name) {
        meas_obj.add_measurement(measurement);
    } else {
        let mut meas_obj = MeasurementObject::new(name);
        meas_obj.add_measurement(measurement);
        viewer.objects_mut().add(meas_obj);
    }
    viewer.request_redraw();
}

// ============================================================================
// distance command
// ============================================================================

struct DistanceCommand;

impl Command for DistanceCommand {
    fn name(&self) -> &str {
        "distance"
    }

    fn aliases(&self) -> &[&str] {
        &["dist"]
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "distance" creates a distance measurement between two atom selections.

USAGE

    distance name, selection1, selection2

ARGUMENTS

    name = string: name for the measurement object
    selection1 = string: first atom selection
    selection2 = string: second atom selection

EXAMPLES

    distance dist1, /1hpx///A/1/CA, /1hpx///A/10/CA
    distance d1, chain A and name CA and resi 1, chain A and name CA and resi 10
"#
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::None, ArgHint::Selection, ArgHint::Selection]
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        let name = args
            .get_str(0)
            .ok_or_else(|| CmdError::MissingArgument("name".into()))?;
        let sel1 = args
            .get_str(1)
            .ok_or_else(|| CmdError::MissingArgument("selection1".into()))?;
        let sel2 = args
            .get_str(2)
            .ok_or_else(|| CmdError::MissingArgument("selection2".into()))?;

        let p1 = resolve_atom_position(ctx.viewer, sel1)?;
        let p2 = resolve_atom_position(ctx.viewer, sel2)?;

        let measurement = Measurement::distance(p1, p2, DEFAULT_COLOR);
        let value = measurement.value;

        add_measurement_to_scene(ctx.viewer, name, measurement);
        ctx.print(&format!(" distance: {:.3} Angstroms", value));

        Ok(())
    }
}

// ============================================================================
// angle command
// ============================================================================

struct AngleCommand;

impl Command for AngleCommand {
    fn name(&self) -> &str {
        "angle"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "angle" creates an angle measurement between three atom selections.
    The angle is measured at the second atom (vertex).

USAGE

    angle name, selection1, selection2, selection3

ARGUMENTS

    name = string: name for the measurement object
    selection1 = string: first atom selection
    selection2 = string: second atom (vertex) selection
    selection3 = string: third atom selection

EXAMPLES

    angle ang1, /1hpx///A/1/CA, /1hpx///A/5/CA, /1hpx///A/10/CA
"#
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[
            ArgHint::None,
            ArgHint::Selection,
            ArgHint::Selection,
            ArgHint::Selection,
        ]
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        let name = args
            .get_str(0)
            .ok_or_else(|| CmdError::MissingArgument("name".into()))?;
        let sel1 = args
            .get_str(1)
            .ok_or_else(|| CmdError::MissingArgument("selection1".into()))?;
        let sel2 = args
            .get_str(2)
            .ok_or_else(|| CmdError::MissingArgument("selection2".into()))?;
        let sel3 = args
            .get_str(3)
            .ok_or_else(|| CmdError::MissingArgument("selection3".into()))?;

        let p1 = resolve_atom_position(ctx.viewer, sel1)?;
        let p2 = resolve_atom_position(ctx.viewer, sel2)?;
        let p3 = resolve_atom_position(ctx.viewer, sel3)?;

        let measurement = Measurement::angle(p1, p2, p3, DEFAULT_COLOR);
        let value = measurement.value;

        add_measurement_to_scene(ctx.viewer, name, measurement);
        ctx.print(&format!(" angle: {:.1} degrees", value));

        Ok(())
    }
}

// ============================================================================
// dihedral command
// ============================================================================

struct DihedralCommand;

impl Command for DihedralCommand {
    fn name(&self) -> &str {
        "dihedral"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "dihedral" creates a dihedral angle measurement between four atom selections.
    The dihedral is measured around the bond between atoms 2 and 3.

USAGE

    dihedral name, selection1, selection2, selection3, selection4

ARGUMENTS

    name = string: name for the measurement object
    selection1 = string: first atom selection
    selection2 = string: second atom selection
    selection3 = string: third atom selection
    selection4 = string: fourth atom selection

EXAMPLES

    dihedral dih1, /1hpx///A/1/N, /1hpx///A/1/CA, /1hpx///A/1/C, /1hpx///A/2/N
"#
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[
            ArgHint::None,
            ArgHint::Selection,
            ArgHint::Selection,
            ArgHint::Selection,
            ArgHint::Selection,
        ]
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        let name = args
            .get_str(0)
            .ok_or_else(|| CmdError::MissingArgument("name".into()))?;
        let sel1 = args
            .get_str(1)
            .ok_or_else(|| CmdError::MissingArgument("selection1".into()))?;
        let sel2 = args
            .get_str(2)
            .ok_or_else(|| CmdError::MissingArgument("selection2".into()))?;
        let sel3 = args
            .get_str(3)
            .ok_or_else(|| CmdError::MissingArgument("selection3".into()))?;
        let sel4 = args
            .get_str(4)
            .ok_or_else(|| CmdError::MissingArgument("selection4".into()))?;

        let p1 = resolve_atom_position(ctx.viewer, sel1)?;
        let p2 = resolve_atom_position(ctx.viewer, sel2)?;
        let p3 = resolve_atom_position(ctx.viewer, sel3)?;
        let p4 = resolve_atom_position(ctx.viewer, sel4)?;

        let measurement = Measurement::dihedral(p1, p2, p3, p4, DEFAULT_COLOR);
        let value = measurement.value;

        add_measurement_to_scene(ctx.viewer, name, measurement);
        ctx.print(&format!(" dihedral: {:.1} degrees", value));

        Ok(())
    }
}
