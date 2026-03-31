//! File I/O commands: load, save, png, fetch, cd, pwd, ls

use std::env;
use std::path::{Path, PathBuf};

use pymol_io::FileFormat;
use pymol_mol::dss::{assign_secondary_structure, DssSettings};
use pymol_mol::ObjectMolecule;
use pymol_render::Grid3D;
use pymol_scene::{MapData, MapObject, MoleculeObject};

use crate::args::ParsedCommand;
use crate::command::{ArgHint, Command, CommandContext, CommandRegistry, ViewerLike};
use crate::error::{CmdError, CmdResult};

use super::objects::extract_molecule;
use super::selecting::evaluate_selection;

/// Expand shell-style paths: ~ to home directory, $VAR to environment variables
pub fn expand_path(path: &str) -> PathBuf {
    #[cfg(feature = "native-io")]
    {
        match shellexpand::full(path) {
            Ok(expanded) => PathBuf::from(expanded.as_ref()),
            Err(_) => PathBuf::from(path),
        }
    }
    #[cfg(not(feature = "native-io"))]
    {
        PathBuf::from(path)
    }
}

/// Apply DSS (Define Secondary Structure) algorithm to a molecule
///
/// This recalculates secondary structure based on backbone phi/psi angles,
/// similar to PyMOL's auto_dss functionality. It overwrites the secondary
/// structure assignments from the PDB/CIF file with computed values.
fn apply_dss(mol: &mut ObjectMolecule) {
    let settings = DssSettings::default();
    let state = 0; // First state
    assign_secondary_structure(mol, state, &settings);
}

/// Register I/O commands
pub fn register(registry: &mut CommandRegistry) {
    registry.register(LoadCommand);
    #[cfg(feature = "traj")]
    registry.register(LoadTrajCommand);
    registry.register(SaveCommand);
    registry.register(PngCommand);
    registry.register(CdCommand);
    registry.register(PwdCommand);
    registry.register(LsCommand);
    #[cfg(any(feature = "fetch", feature = "fetch-async"))]
    registry.register(FetchCommand);
}

// ============================================================================
// load command
// ============================================================================

/// Load a molecular structure file
struct LoadCommand;

impl Command for LoadCommand {
    fn name(&self) -> &str {
        "load"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "load" reads several molecular file formats and creates a new
    object or adds states to an existing object.

USAGE

    load filename [, object [, state [, format [, quiet ]]]]

ARGUMENTS

    filename = string: file path or URL
    object = string: object name (default: filename stem)
    state = integer: state to load into (0 = append, default)
    format = string: file format (auto-detected if omitted)
    quiet = 0/1: suppress feedback (default: 0)

EXAMPLES

    load protein.pdb
    load ligand.sdf, object=lig
    load structure.cif, format=cif
    load session.pse         → import PyMOL session
    load session.prs         → load native session
"#
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::Path, ArgHint::Object]
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        // Get filename (required)
        let filename = args
            .get_str(0)
            .or_else(|| args.get_named_str("filename"))
            .ok_or_else(|| CmdError::MissingArgument("filename".to_string()))?;

        // Get object name (optional, defaults to filename stem)
        let object_name = args
            .get_str(1)
            .or_else(|| args.get_named_str("object"))
            .map(|s| s.to_string())
            .unwrap_or_else(|| {
                let p = Path::new(filename);
                let stem = p.file_stem().and_then(|s| s.to_str()).unwrap_or("obj");
                // Strip double extension for .gz files (e.g. "1AOI.cif.gz" → "1AOI")
                if p.extension().and_then(|s| s.to_str()) == Some("gz") {
                    Path::new(stem)
                        .file_stem()
                        .and_then(|s| s.to_str())
                        .unwrap_or(stem)
                        .to_string()
                } else {
                    stem.to_string()
                }
            });

        // Get state (optional, 0 = append)
        let state = args
            .get_int(2)
            .or_else(|| args.get_named_int("state"))
            .unwrap_or(0);

        if state != 0 {
            ctx.print(" Warning: state= parameter not yet fully supported, loading all states.");
        }

        // Get format (optional, auto-detect)
        let format = args
            .get_str(3)
            .or_else(|| args.get_named_str("format"))
            .map(|s| match s.to_lowercase().as_str() {
                "pdb" => FileFormat::Pdb,
                "sdf" | "mol" => FileFormat::Sdf,
                "mol2" => FileFormat::Mol2,
                "xyz" => FileFormat::Xyz,
                "cif" | "mmcif" => FileFormat::Cif,
                "gro" => FileFormat::Gro,
                "xtc" => FileFormat::Xtc,
                "trr" => FileFormat::Trr,
                "ccp4" | "map" | "mrc" => FileFormat::Ccp4,
                _ => FileFormat::Unknown,
            });

        // Get quiet flag
        let quiet = args
            .get_bool(4)
            .or_else(|| args.get_named_bool("quiet"))
            .unwrap_or(false);

        // Expand path and check for session file formats
        let path = expand_path(filename);
        let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("").to_lowercase();
        match ext.as_str() {
            "pse" | "pze" => {
                use pymol_session::{load_pse, pse_to_session};
                let pse = load_pse(&path).map_err(|e| CmdError::FileFormat(e.to_string()))?;
                let mut session = pse_to_session(&pse).map_err(|e| CmdError::FileFormat(e.to_string()))?;
                // Merge PSE named colors into the viewer's color table,
                // remapping atom color indices to match the viewer's registry.
                let index_map = ctx.viewer.named_colors_mut().merge_from(&session.named_colors);
                let mol_names: Vec<String> = session.registry.names().map(|s| s.to_string()).collect();
                let n_objs = mol_names.len();
                // Move objects into viewer, then remap atom color indices
                for name in &mol_names {
                    if let Some(obj) = session.registry.remove(name) {
                        ctx.viewer.objects_mut().insert_boxed(name, obj);
                    }
                }
                // Remap color indices on molecule atoms to match viewer's color table
                for name in &mol_names {
                    if let Some(mol_obj) = ctx.viewer.objects_mut().get_molecule_mut(name) {
                        for atom in mol_obj.molecule_mut().atoms_mut() {
                            let base = atom.repr.colors.base;
                            if base >= 0 {
                                if let Some(&new_idx) = index_map.get(&(base as u32)) {
                                    atom.repr.colors.base = new_idx as i32;
                                }
                            }
                        }
                    }
                }
                if !quiet {
                    ctx.print(&format!(" Loaded PSE session \"{}\" ({} objects)", filename, n_objs));
                }
                // Zoom to show everything
                ctx.viewer.zoom_all(0.0);
                ctx.viewer.update_movie_state_count();
                return Ok(());
            }
            "prs" => {
                use pymol_session::load_prs;
                let session = load_prs(&path).map_err(|e| CmdError::FileFormat(e.to_string()))?;
                ctx.viewer.replace_session(session);
                ctx.viewer.update_movie_state_count();
                if !quiet {
                    ctx.print(&format!(" Loaded PRS session \"{}\"", filename));
                }
                return Ok(());
            }
            _ => {}
        }

        // Detect format
        let builtin_format = format.unwrap_or_else(|| FileFormat::from_path(&path));

        // Reject trajectory-only formats (no topology)
        if builtin_format.is_trajectory_only() {
            return Err(CmdError::execution(format!(
                "{} is a trajectory-only format with no topology. Use 'load_traj' instead.",
                builtin_format.name()
            )));
        }

        // Handle CCP4/MRC map files as MapObject (not MoleculeObject)
        if builtin_format.is_map_format() {
            let ccp4 = pymol_io::ccp4::read_ccp4(&path)
                .map_err(|e| CmdError::FileFormat(e.to_string()))?;

            if !quiet {
                ctx.print(&format!(
                    " CCP4 map: {}x{}x{} grid, spacing ({:.2}, {:.2}, {:.2}) Å",
                    ccp4.dims[0] + 1, ccp4.dims[1] + 1, ccp4.dims[2] + 1,
                    ccp4.spacing[0], ccp4.spacing[1], ccp4.spacing[2]
                ));
                ctx.print(&format!(
                    " Density: mean {:.4}, rms {:.4} (normalized to sigma)",
                    ccp4.stats.dmean, ccp4.stats.rms
                ));
                ctx.print(&format!(
                    " Origin: ({:.2}, {:.2}, {:.2})",
                    ccp4.origin[0], ccp4.origin[1], ccp4.origin[2]
                ));
            }

            let grid = Grid3D::from_dims(ccp4.origin, ccp4.spacing, ccp4.dims, ccp4.values);
            let map_data = MapData::new(grid);
            let map_obj = MapObject::from_map_data(&object_name, map_data);
            ctx.viewer.objects_mut().add(map_obj);

            if !quiet {
                ctx.print(&format!(" Loaded map \"{}\" as \"{}\"", filename, object_name));
            }

            ctx.viewer.zoom_on(&object_name, 0.0);
            return Ok(());
        }

        // Load molecular file
        let mut mol = if builtin_format != FileFormat::Unknown {
            // Built-in format
            pymol_io::read_file_format(&path, builtin_format)
                .map_err(|e| CmdError::FileFormat(e.to_string()))?
        } else if let Some(handler) = ctx.format_handler(&ext) {
            // Plugin-provided format
            let reader_fn = handler.reader.as_ref().ok_or_else(|| {
                CmdError::execution(format!(
                    "Format '{}' does not support reading",
                    handler.name
                ))
            })?;
            let file = std::fs::File::open(&path)
                .map_err(|e| CmdError::FileFormat(e.to_string()))?;
            let buf_reader: Box<dyn std::io::Read> = Box::new(std::io::BufReader::new(file));
            let molecules = (reader_fn)(buf_reader)
                .map_err(CmdError::FileFormat)?;
            molecules.into_iter().next().ok_or_else(|| {
                CmdError::FileFormat("No molecules in file".to_string())
            })?
        } else {
            return Err(CmdError::execution(format!(
                "Unknown file format: .{} (no built-in or plugin handler found)",
                ext
            )));
        };

        // Apply DSS (Define Secondary Structure) if auto_dss is enabled
        // This recalculates secondary structure based on backbone geometry
        let auto_dss = ctx.viewer.settings().behavior.auto_dss;
        if auto_dss {
            apply_dss(&mut mol);
        }

        // Create molecule object and add to viewer
        let mol_obj = MoleculeObject::with_name(mol, &object_name);
        ctx.viewer.objects_mut().add(mol_obj);

        if !quiet {
            ctx.print(&format!(" Loaded \"{}\" as \"{}\"", filename, object_name));
        }

        // Zoom to loaded molecule (preserves rotation)
        ctx.viewer.zoom_on(&object_name, 0.0);
        ctx.viewer.update_movie_state_count();

        Ok(())
    }
}

// ============================================================================
// load_traj command
// ============================================================================

#[cfg(feature = "traj")]
/// Load trajectory coordinates onto an existing molecular object
struct LoadTrajCommand;

#[cfg(feature = "traj")]
impl Command for LoadTrajCommand {
    fn name(&self) -> &str {
        "load_traj"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "load_traj" loads trajectory coordinates from a file onto an
    existing molecular object. The trajectory must have the same
    number of atoms as the target object.

USAGE

    load_traj filename [, object [, start [, stop [, interval [, format [, append]]]]]]

ARGUMENTS

    filename = string: trajectory file path (XTC, TRR)
    object = string: target object name (default: first molecule)
    start = integer: first frame to load, 1-based (default: 1)
    stop = integer: last frame to load, 1-based (default: all)
    interval = integer: load every Nth frame (default: 1)
    format = string: file format (auto-detected if omitted)
    append = bool: append loaded states to existed (default: false)

EXAMPLES

    load protein.gro
    load_traj trajectory.xtc, protein
    load_traj long_sim.xtc, protein, start=1, stop=1000, interval=10
    load_traj trajectory.trr, protein
"#
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::Path, ArgHint::Object]
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        // Get filename (required)
        let filename = args
            .get_str(0)
            .or_else(|| args.get_named_str("filename"))
            .ok_or_else(|| CmdError::MissingArgument("filename".to_string()))?;

        // Get object name (optional, defaults to first loaded molecule)
        let object_name = args
            .get_str(1)
            .or_else(|| args.get_named_str("object"));

        // PyMOL uses 1-based frame numbers; convert to 0-based
        let start = args
            .get_int(2)
            .or_else(|| args.get_named_int("start"))
            .map(|s| (s.max(1) - 1) as usize)
            .unwrap_or(0);

        let stop = args
            .get_int(3)
            .or_else(|| args.get_named_int("stop"))
            .map(|s| s.max(0) as usize);

        let interval = args
            .get_int(4)
            .or_else(|| args.get_named_int("interval"))
            .map(|i| i.max(1) as usize)
            .unwrap_or(1);

        let append = args
            .get_bool(5)
            .unwrap_or(args.get_named_bool_or("append", false));

        let quiet = args
            .get_named_bool("quiet")
            .unwrap_or(false);

        // Resolve target object name
        let target_name = if let Some(name) = object_name {
            name.to_string()
        } else {
            // Find first molecule object
            let mut first_mol = None;
            for name in ctx.viewer.objects().names() {
                if ctx.viewer.objects().get_molecule(name).is_some() {
                    first_mol = Some(name.to_string());
                    break;
                }
            }
            first_mol.ok_or_else(|| {
                CmdError::execution("No molecule objects loaded. Load a structure first.")
            })?
        };

        // Validate target exists and get atom count
        let expected_atoms = {
            let mol_obj = ctx
                .viewer
                .objects()
                .get_molecule(&target_name)
                .ok_or_else(|| CmdError::ObjectNotFound(target_name.clone()))?;
            mol_obj.molecule().atom_count()
        };

        // Detect format
        let path = expand_path(filename);
        let format = args
            .get_named_str("format")
            .map(|s| match s.to_lowercase().as_str() {
                "xtc" => FileFormat::Xtc,
                "trr" => FileFormat::Trr,
                _ => FileFormat::from_path(&path),
            })
            .unwrap_or_else(|| FileFormat::from_path(&path));

        // Read trajectory frames
        let opts = pymol_io::TrajectoryReadOptions {
            start,
            stop,
            interval,
        };

        let frames = pymol_io::read_trajectory_format(&path, format, &opts)
            .map_err(|e| CmdError::FileFormat(e.to_string()))?;

        if frames.is_empty() {
            return Err(CmdError::execution("No frames found in trajectory"));
        }

        // Validate atom count
        let traj_atoms = frames[0].len();
        if traj_atoms != expected_atoms {
            return Err(CmdError::execution(format!(
                "Atom count mismatch: object '{}' has {} atoms, trajectory has {}",
                target_name, expected_atoms, traj_atoms
            )));
        }

        // Replace existing coordinates with trajectory frames
        let n_frames = frames.len();
        let mol_obj = ctx
            .viewer
            .objects_mut()
            .get_molecule_mut(&target_name)
            .ok_or_else(|| CmdError::ObjectNotFound(target_name.clone()))?;

        if !append {
            mol_obj.molecule_mut().clear_coord_sets();
        }
        for cs in frames {
            mol_obj.molecule_mut().add_coord_set(cs);
        }

        let total_states = mol_obj.molecule().state_count();

        if !quiet {
            ctx.print(&format!(
                " Loaded {} frames from \"{}\" onto \"{}\" ({} total states)",
                n_frames, filename, target_name, total_states
            ));
        }

        ctx.viewer.update_movie_state_count();

        Ok(())
    }
}

// ============================================================================
// save command
// ============================================================================

/// Save a molecular structure to file
struct SaveCommand;

impl Command for SaveCommand {
    fn name(&self) -> &str {
        "save"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "save" writes molecular data to a file.

USAGE

    save filename [, selection [, state [, format ]]]

ARGUMENTS

    filename = string: output file path
    selection = string: atoms to save (default: all)
    state = integer: state to save (default: -1 = current)
    format = string: file format (auto-detected from extension)

EXAMPLES

    save output.pdb
    save ligand.sdf, organic
    save protein.mol2, polymer
    save session.prs         → save native session
"#
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::Path, ArgHint::Selection]
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let filename = args
            .get_str(0)
            .or_else(|| args.get_named_str("filename"))
            .ok_or_else(|| CmdError::MissingArgument("filename".to_string()))?;

        let selection = args
            .get_str(1)
            .or_else(|| args.get_named_str("selection"))
            .unwrap_or("all");

        let state_arg = args
            .get_int(2)
            .or_else(|| args.get_named_int("state"))
            .unwrap_or(-1);

        let path = expand_path(filename);
        let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("").to_lowercase();

        // Session formats — handled before anything else
        match ext.as_str() {
            "prs" => {
                use pymol_session::save_prs;
                let session = ctx.viewer.session();
                save_prs(session, &path).map_err(|e| CmdError::FileFormat(e.to_string()))?;
                ctx.print(&format!(" Saved session to \"{}\"", filename));
                return Ok(());
            }
            "pse" => {
                return Err(CmdError::execution("PSE format is read-only. Use .prs for saving sessions."));
            }
            "bcif" => {
                return Err(CmdError::execution(
                    "BinaryCIF writing is not supported. Save as .cif instead.",
                ));
            }
            "pml" => {
                let mut file = std::fs::File::create(&path)
                    .map_err(|e| CmdError::FileFormat(e.to_string()))?;
                use std::io::Write;
                writeln!(file, "# PyMOL-RS command log")
                    .map_err(|e| CmdError::FileFormat(e.to_string()))?;
                let count = if let Some(history) = ctx.history() {
                    let mut n = 0;
                    for cmd in history.iter() {
                        writeln!(file, "{}", cmd)
                            .map_err(|e| CmdError::FileFormat(e.to_string()))?;
                        n += 1;
                    }
                    n
                } else {
                    0
                };
                ctx.print(&format!(
                    " Saved {} commands to \"{}\"",
                    count,
                    filename
                ));
                return Ok(());
            }
            _ => {}
        }

        let format = FileFormat::from_extension(&ext);
        if format == FileFormat::Unknown && ctx.format_handler(&ext).is_none() {
            return Err(CmdError::execution(format!(
                "Unknown file format: .{} (supported: pdb, cif, mol2, sdf, xyz, gro, prs + plugins)",
                ext
            )));
        }

        // Evaluate selection across all objects
        let results = evaluate_selection(ctx.viewer, selection)?;
        let matches: Vec<(String, _)> = results
            .into_iter()
            .filter(|(_, sel)| sel.count() > 0)
            .collect();

        if matches.is_empty() {
            return Err(CmdError::Selection(format!(
                "No atoms match '{}'",
                selection
            )));
        }

        // Resolve state index (PyMOL: -1/0 = current, 1-based otherwise)
        let state_idx: Option<usize> = match state_arg {
            n if n <= 0 => None,
            n => Some((n - 1) as usize),
        };

        // Build filtered molecules
        let filtered: Vec<ObjectMolecule> = matches
            .iter()
            .filter_map(|(name, sel)| {
                let mol_obj = ctx.viewer.objects().get_molecule(name)?;
                Some(extract_molecule(mol_obj.molecule(), sel, name, state_idx))
            })
            .filter(|m| m.atom_count() > 0)
            .collect();

        if filtered.is_empty() {
            return Err(CmdError::execution(
                "No atoms to save after applying selection",
            ));
        }

        let total_atoms: usize = filtered.iter().map(|m| m.atom_count()).sum();

        // Write
        let format_name = if format != FileFormat::Unknown {
            // Built-in format
            pymol_io::write_all(&path, &filtered)
                .map_err(|e| CmdError::FileFormat(e.to_string()))?;
            format.name().to_string()
        } else {
            // Plugin-provided format
            let handler = ctx.format_handler(&ext).unwrap();
            let writer_fn = handler.writer.as_ref().ok_or_else(|| {
                CmdError::execution(format!(
                    "Format '{}' does not support writing",
                    handler.name
                ))
            })?;
            let file = std::fs::File::create(&path)
                .map_err(|e| CmdError::FileFormat(e.to_string()))?;
            let buf_writer: Box<dyn std::io::Write> = Box::new(std::io::BufWriter::new(file));
            (writer_fn)(buf_writer, &filtered)
                .map_err(CmdError::FileFormat)?;
            handler.name.clone()
        };

        ctx.print(&format!(
            " Saved {} atoms to \"{}\" ({} format)",
            total_atoms,
            filename,
            format_name
        ));

        Ok(())
    }
}

// ============================================================================
// png command
// ============================================================================

/// Save current view as PNG image
struct PngCommand;

impl Command for PngCommand {
    fn name(&self) -> &str {
        "png"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "png" saves a PNG format image file of the current display.

USAGE

    png filename [, width [, height [, dpi [, quiet]]]]

ARGUMENTS

    filename = string: file path to be written
    width = integer: width in pixels (default: current window width)
    height = integer: height in pixels (default: current window height)
    dpi = float: dots-per-inch (not yet implemented)
    quiet = 0/1: suppress feedback (default: 0)

NOTES

    If only width or height is specified, the aspect ratio is preserved.
    Use the "ray" command (from the raytracer plugin) for ray-traced images.

EXAMPLES

    png ~/Desktop/screenshot.png
    png output.png, 1920, 1080
    png output.png, width=800
    png hires.png, 4096, 4096, quiet=1
"#
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::Path]
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        // Get filename (required)
        let filename = args
            .get_str(0)
            .or_else(|| args.get_named_str("filename"))
            .ok_or_else(|| CmdError::MissingArgument("filename".to_string()))?;

        // Get width (optional)
        let width = args
            .get_int(1)
            .or_else(|| args.get_named_int("width"))
            .map(|v| v as u32);

        // Get height (optional)
        let height = args
            .get_int(2)
            .or_else(|| args.get_named_int("height"))
            .map(|v| v as u32);

        // Get dpi (optional, not yet used)
        let _dpi = args
            .get_float(3)
            .or_else(|| args.get_named_float("dpi"));

        // Get quiet flag
        let quiet = args
            .get_bool(4)
            .or_else(|| args.get_named_bool("quiet"))
            .unwrap_or(false);

        // Ensure the path ends with .png
        let path = expand_path(filename);
        let path = if path.extension().map(|e| e.to_ascii_lowercase()) != Some("png".into()) {
            path.with_extension("png")
        } else {
            path.to_path_buf()
        };

        if let Some(vp_img) = ctx.viewer.get_viewport_image() {
            // Export stored viewport image overlay
            #[cfg(feature = "native-io")]
            {
                use image::RgbaImage;

                let img =
                    RgbaImage::from_raw(vp_img.width, vp_img.height, vp_img.data.clone())
                        .ok_or_else(|| {
                            CmdError::execution("Invalid viewport image data".to_string())
                        })?;

                img.save(&path)
                    .map_err(|e| CmdError::execution(format!("Failed to save PNG: {}", e)))?;

                if !quiet {
                    ctx.print(&format!(" ({}x{})", vp_img.width, vp_img.height));
                    ctx.print(&format!(
                        " Saved viewport image to \"{}\"",
                        path.display()
                    ));
                }
            }
            #[cfg(not(feature = "native-io"))]
            {
                let _ = vp_img;
                return Err(CmdError::execution(
                    "PNG export not available in this build".to_string(),
                ));
            }
        } else {
            // Capture rasterized screenshot
            ctx.viewer
                .capture_png(&path, width, height)
                .map_err(|e| CmdError::execution(format!("Failed to capture PNG: {}", e)))?;

            if !quiet {
                let (w, h) = match (width, height) {
                    (Some(w), Some(h)) => (w, h),
                    _ => {
                        // Report actual captured dimensions
                        // For now just report the requested or window size
                        (width.unwrap_or(0), height.unwrap_or(0))
                    }
                };
                if w > 0 && h > 0 {
                    ctx.print(&format!(" ({}x{})", w, h));
                }
                ctx.print(&format!(" Saved \"{}\"", path.display()));
            }
        }

        Ok(())
    }
}

// ============================================================================
// cd command
// ============================================================================

/// Change current working directory
struct CdCommand;

impl Command for CdCommand {
    fn name(&self) -> &str {
        "cd"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "cd" changes the current working directory.

USAGE

    cd path

ARGUMENTS

    path = string: directory path (default: home directory)

EXAMPLES

    cd /home/user/molecules
    cd ..
    cd
"#
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::Path]
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let path = args.get_str(0).or_else(|| args.get_named_str("path"));

        let target = if let Some(p) = path {
            expand_path(p)
        } else {
            // Default to home directory
            dirs_path_or_current()
        };

        env::set_current_dir(&target).map_err(|e| {
            CmdError::execution(format!("Cannot change to '{}': {}", target.display(), e))
        })?;

        ctx.print(&format!(" Changed to {}", target.display()));

        Ok(())
    }
}

/// Get home directory or current directory
fn dirs_path_or_current() -> std::path::PathBuf {
    env::var("HOME")
        .map(|h| Path::new(&h).to_path_buf())
        .unwrap_or_else(|_| env::current_dir().unwrap_or_else(|_| Path::new(".").to_path_buf()))
}

// ============================================================================
// pwd command
// ============================================================================

/// Print current working directory
struct PwdCommand;

impl Command for PwdCommand {
    fn name(&self) -> &str {
        "pwd"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "pwd" prints the current working directory.

USAGE

    pwd
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, _args: &ParsedCommand) -> CmdResult {
        let cwd = env::current_dir().map_err(|e| CmdError::execution(e.to_string()))?;

        ctx.print(&format!(" {}", cwd.display()));

        Ok(())
    }
}

// ============================================================================
// ls command
// ============================================================================

/// List directory contents
struct LsCommand;

impl Command for LsCommand {
    fn name(&self) -> &str {
        "ls"
    }

    fn aliases(&self) -> &[&str] {
        &["dir"]
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "ls" lists the contents of a directory.

USAGE

    ls [path]

ARGUMENTS

    path = string: directory path (default: current directory)

EXAMPLES

    ls
    ls /home/user/molecules
"#
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::Path]
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let path_str = args
            .get_str(0)
            .or_else(|| args.get_named_str("path"))
            .unwrap_or(".");

        let dir = expand_path(path_str);

        if !dir.is_dir() {
            return Err(CmdError::execution(format!(
                "'{}' is not a directory",
                path_str
            )));
        }

        let mut entries: Vec<String> = std::fs::read_dir(dir)
            .map_err(|e| CmdError::execution(e.to_string()))?
            .filter_map(|e| e.ok())
            .map(|e| {
                let name = e.file_name().to_string_lossy().to_string();
                if e.file_type().map(|t| t.is_dir()).unwrap_or(false) {
                    format!("{}/", name)
                } else {
                    name
                }
            })
            .collect();

        entries.sort();

        // Format in columns
        for entry in entries {
            ctx.print(&format!("   {}", entry));
        }

        Ok(())
    }
}

// ============================================================================
// fetch command (optional, requires fetch feature)
// ============================================================================

#[cfg(any(feature = "fetch", feature = "fetch-async"))]
struct FetchCommand;

#[cfg(any(feature = "fetch", feature = "fetch-async"))]
impl Command for FetchCommand {
    fn name(&self) -> &str {
        "fetch"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "fetch" downloads a structure from the RCSB PDB.

USAGE

    fetch code [, name [, type [, async ]]]

ARGUMENTS

    code = string: PDB ID (e.g., "1ubq")
    name = string: object name (default: PDB ID)
    type = string: file type to fetch (default: cif)
           Supported: "pdb", "cif" (or "mmcif")
    async = 0/1: asynchronous fetch (default: 1)
            When 0, forces synchronous (blocking) fetch

EXAMPLES

    fetch 1ubq
    fetch 4hhb, name=hemoglobin
    fetch 1crn, type=pdb
    fetch 1igt, async=0
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let code = args
            .get_str(0)
            .or_else(|| args.get_named_str("code"))
            .ok_or_else(|| CmdError::MissingArgument("code".to_string()))?;

        let name = args
            .get_str(1)
            .or_else(|| args.get_named_str("name"))
            .unwrap_or(code);

        // Parse the type argument (default: bcif)
        let format = args
            .get_str(2)
            .or_else(|| args.get_named_str("type"))
            .map(|s| match s.to_lowercase().as_str() {
                "pdb" => pymol_io::FetchFormat::Pdb,
                "cif" | "mmcif" => pymol_io::FetchFormat::Cif,
                "bcif" | "binarycif" => pymol_io::FetchFormat::Bcif,
                _ => pymol_io::FetchFormat::default(),
            })
            .unwrap_or_default();

        // Parse async flag (default: true = non-blocking)
        let use_async = args.get_named_bool_or("async", true);

        // Convert format to code for async API
        let format_code = match format {
            pymol_io::FetchFormat::Pdb => 1u8,
            pymol_io::FetchFormat::Cif => 0u8,
            pymol_io::FetchFormat::Bcif => 2u8,
        };

        // Try async path first (GUI supports this)
        if use_async && ctx.viewer.request_async_fetch(code, name, format_code) {
            ctx.print(&format!(" Fetching {}...", code));
            return Ok(());
        }

        // Sync fallback for non-GUI viewers (headless, scripts, etc.)
        #[cfg(feature = "fetch")]
        {
            let mut mol = pymol_io::fetch(code, format)
                .map_err(|e| CmdError::FileFormat(e.to_string()))?;

            // Apply DSS if auto_dss is enabled
            let auto_dss = ctx.viewer.settings().behavior.auto_dss;
            if auto_dss {
                apply_dss(&mut mol);
            }

            // Add to viewer
            ctx.viewer
                .objects_mut()
                .add(MoleculeObject::with_name(mol, name));

            ctx.print(&format!(" Fetched {} as \"{}\"", code, name));

            // Zoom to fetched molecule (preserves rotation)
            ctx.viewer.zoom_on(name, 0.0);
            ctx.viewer.update_movie_state_count();

            Ok(())
        }

        // No fetch implementation available
        #[cfg(not(feature = "fetch"))]
        {
            return Err(CmdError::execution(
                "Fetch not available: neither async nor sync fetch is enabled",
            ));
        }
    }
}
