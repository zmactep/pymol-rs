//! Python Plugin Commands
//!
//! Implements commands for the Python plugin:
//! - `python` (alias `/`) — execute Python code inline
//! - `iterate` — execute a Python expression for each atom in a selection (read-only)
//! - `alter` — like iterate, but allows modifying atom properties
//!
//! All commands are non-blocking: code is submitted to the Python worker
//! thread, and output appears asynchronously via the poll cycle.

use std::collections::HashMap;
use std::sync::{Arc, Mutex};

use pymol_mol::{Element, SecondaryStructure};
use pymol_plugin::prelude::*;

use crate::backend::SharedStateHandle;
use crate::worker::{WorkItem, WorkOrigin, WorkerHandle};

// =============================================================================
// Alter types
// =============================================================================

/// A single atom property value that can be changed by `alter`.
#[derive(Debug, Clone)]
pub enum PropertyValue {
    Str(String),
    F32(f32),
    I32(i32),
    I8(i8),
    Bool(bool),
}

/// A batch of property changes for a single atom.
#[derive(Debug, Clone)]
pub struct AtomChange {
    pub obj: String,
    pub idx: u32,
    pub changes: HashMap<String, PropertyValue>,
}

/// Thread-safe buffer for atom changes produced by `alter()` on the Python
/// worker thread and consumed via `queue_viewer_mutation` on the main thread.
pub type AlterBuffer = Arc<Mutex<Vec<Vec<AtomChange>>>>;

/// Apply property changes from a HashMap to a mutable Atom.
pub fn apply_property_changes(atom: &mut Atom, changes: &HashMap<String, PropertyValue>) {
    for (key, value) in changes {
        match key.as_str() {
            "name" => {
                if let PropertyValue::Str(v) = value {
                    atom.name = Arc::from(v.as_str());
                }
            }
            "b" => {
                if let PropertyValue::F32(v) = value {
                    atom.b_factor = *v;
                }
            }
            "q" => {
                if let PropertyValue::F32(v) = value {
                    atom.occupancy = *v;
                }
            }
            "vdw" => {
                if let PropertyValue::F32(v) = value {
                    atom.vdw = *v;
                }
            }
            "partial_charge" => {
                if let PropertyValue::F32(v) = value {
                    atom.partial_charge = *v;
                }
            }
            "formal_charge" => {
                if let PropertyValue::I8(v) = value {
                    atom.formal_charge = *v;
                }
            }
            "color" => {
                if let PropertyValue::I32(v) = value {
                    atom.repr.colors.base = *v;
                }
            }
            "elem" => {
                if let PropertyValue::Str(v) = value {
                    if let Some(el) = Element::from_symbol(v) {
                        atom.element = el;
                    }
                }
            }
            "ss" => {
                if let PropertyValue::Str(v) = value {
                    atom.ss_type = match v.as_str() {
                        "H" => SecondaryStructure::Helix,
                        "S" => SecondaryStructure::Sheet,
                        _ => SecondaryStructure::Loop,
                    };
                }
            }
            "type" => {
                if let PropertyValue::Str(v) = value {
                    atom.state.hetatm = v == "HETATM";
                }
            }
            "alt" => {
                if let PropertyValue::Str(v) = value {
                    atom.alt = v.chars().next().unwrap_or(' ');
                }
            }
            // Residue fields require Arc clone-on-write
            "chain" | "resn" | "resv" | "segi" => {
                let mut res = (*atom.residue).clone();
                match key.as_str() {
                    "chain" => {
                        if let PropertyValue::Str(v) = value {
                            res.key.chain = v.clone();
                        }
                    }
                    "resn" => {
                        if let PropertyValue::Str(v) = value {
                            res.key.resn = v.clone();
                        }
                    }
                    "resv" => {
                        if let PropertyValue::I32(v) = value {
                            res.key.resv = *v;
                        }
                    }
                    "segi" => {
                        if let PropertyValue::Str(v) = value {
                            res.segi = v.clone();
                        }
                    }
                    _ => {}
                }
                atom.residue = Arc::new(res);
            }
            _ => {} // Unknown property — ignore
        }
    }
}

// =============================================================================
// Shared state sync
// =============================================================================

/// Sync molecule snapshots from the viewer into shared state.
///
/// Called before submitting `iterate`/`alter` work to the Python worker so
/// the worker sees up-to-date molecule data, even when no GUI poll cycle
/// has run (e.g., during .pml script execution).
fn sync_shared_molecules(shared: &SharedStateHandle, viewer: &dyn ViewerLike) {
    let mut state = shared.lock().unwrap();
    state.names = viewer.objects().names().map(|s| s.to_string()).collect();
    state.molecules.clear();
    for name in viewer.objects().names() {
        if let Some(mol_obj) = viewer.objects().get_molecule(name) {
            state
                .molecules
                .push((name.to_string(), mol_obj.molecule().clone()));
        }
    }
}

// =============================================================================
// Helpers
// =============================================================================

/// Collect all arguments as strings, preserving order.
///
/// Named arguments are reconstructed as `key=value`. String values in
/// named args are re-quoted (the parser strips quotes, but we need them
/// for Python expressions like `chain="C"`).
fn collect_all_args(args: &ParsedCommand) -> Vec<String> {
    args.args
        .iter()
        .filter_map(|(name, val)| {
            let repr = val.to_string_repr()?;
            match name {
                Some(key) => {
                    // Re-quote string values: the parser strips quotes from
                    // `chain="C"` → String("C"), but Python needs `chain="C"`.
                    let rhs = if val.as_str().is_some() {
                        format!("\"{}\"", repr)
                    } else {
                        repr
                    };
                    Some(format!("{}={}", key, rhs))
                }
                None => Some(repr),
            }
        })
        .collect()
}

/// Escape a string for inclusion in Python source code (single-quoted).
fn python_escape(s: &str) -> String {
    s.replace('\\', "\\\\").replace('\'', "\\'")
}

/// Shared logic for `iterate` and `alter` commands.
///
/// Parses `<selection>, <expression>` from args, wraps into
/// `cmd.<method>('<selection>', '<expression>')`, and submits to the worker.
fn submit_atom_command(
    worker: &WorkerHandle,
    ctx: &mut CommandContext<'_, '_, dyn ViewerLike + '_>,
    args: &ParsedCommand,
    method: &str,
) -> CmdResult {
    let all_args = collect_all_args(args);
    if all_args.len() < 2 {
        ctx.print(&format!("Usage: {} <selection>, <expression>", method));
        return Ok(());
    }

    let selection = &all_args[0];
    let expression = all_args[1..].join(", ");

    let code = format!(
        "cmd.{}('{}', '{}')",
        method,
        python_escape(selection),
        python_escape(&expression),
    );

    worker.submit(WorkItem::Eval {
        code,
        origin: WorkOrigin::Command,
    });

    Ok(())
}

// =============================================================================
// PythonCommand
// =============================================================================

/// The `python` command — executes Python code inline.
///
/// Usage:
///   python print("hello")
///   /import math; print(math.pi)
pub struct PythonCommand {
    pub(crate) worker: WorkerHandle,
}

impl PythonCommand {
    pub fn new(worker: WorkerHandle) -> Self {
        Self { worker }
    }
}

impl Command for PythonCommand {
    fn name(&self) -> &str {
        "python"
    }

    fn aliases(&self) -> &[&str] {
        &["/"]
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        // Reconstruct the code from all args.
        // The parser splits on commas, so we rejoin them to reconstruct
        // the original Python expression (e.g., "print(1, 2)").
        let code = collect_all_args(args).join(", ");
        let code = code.trim();
        if code.is_empty() {
            ctx.print("Usage: python <code>  or  /<code>");
            return Ok(());
        }

        self.worker.submit(WorkItem::Eval {
            code: code.to_string(),
            origin: WorkOrigin::Command,
        });

        Ok(())
    }

    fn help(&self) -> &str {
        "python <code>\n/code\n\n\
         Execute Python code.\n\n\
         Examples:\n\
             python print(\"hello\")\n\
             /import math; print(math.pi)\n\
             python from pymol_rs import cmd; cmd.color(\"red\", \"all\")"
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::None]
    }
}

// =============================================================================
// IterateCommand
// =============================================================================

/// The `iterate` command — executes a Python expression for each atom
/// in a selection (read-only).
///
/// Usage:
///   iterate all, print(name, resn, b)
///   iterate chain A, mylist.append(b)
pub struct IterateCommand {
    worker: WorkerHandle,
    shared: SharedStateHandle,
}

impl IterateCommand {
    pub fn new(worker: WorkerHandle, shared: SharedStateHandle) -> Self {
        Self { worker, shared }
    }
}

impl Command for IterateCommand {
    fn name(&self) -> &str {
        "iterate"
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        sync_shared_molecules(&self.shared, ctx.viewer);
        submit_atom_command(&self.worker, ctx, args, "iterate")
    }

    fn help(&self) -> &str {
        "iterate <selection>, <expression>\n\n\
         Execute a Python expression for each atom in a selection (read-only).\n\n\
         Available atom properties:\n\
             name, resn, resv, resi, chain, segi, alt, elem,\n\
             b, q, vdw, partial_charge, formal_charge,\n\
             ss, color, type, hetatm, index, ID, rank, model,\n\
             x, y, z\n\n\
         Examples:\n\
             iterate all, print(name, resn, b)\n\
             iterate chain A, mylist.append(b)"
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::Selection, ArgHint::None]
    }
}

// =============================================================================
// AlterCommand
// =============================================================================

/// The `alter` command — executes a Python expression for each atom
/// in a selection, allowing modification of atom properties.
///
/// Usage:
///   alter all, b=0.0
///   alter chain A, chain='B'
pub struct AlterCommand {
    worker: WorkerHandle,
    shared: SharedStateHandle,
}

impl AlterCommand {
    pub fn new(worker: WorkerHandle, shared: SharedStateHandle) -> Self {
        Self { worker, shared }
    }
}

impl Command for AlterCommand {
    fn name(&self) -> &str {
        "alter"
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        sync_shared_molecules(&self.shared, ctx.viewer);
        submit_atom_command(&self.worker, ctx, args, "alter")
    }

    fn help(&self) -> &str {
        "alter <selection>, <expression>\n\n\
         Execute a Python expression for each atom in a selection,\n\
         allowing modification of atom properties.\n\n\
         Mutable properties:\n\
             name, resn, resv, chain, segi, alt, elem,\n\
             b, q, vdw, partial_charge, formal_charge,\n\
             ss, color, type\n\n\
         Read-only: x, y, z, index, ID, rank, model, hetatm\n\n\
         Examples:\n\
             alter all, b=0.0\n\
             alter chain A, chain='B'\n\
             alter name CA, vdw=2.0"
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::Selection, ArgHint::None]
    }
}
