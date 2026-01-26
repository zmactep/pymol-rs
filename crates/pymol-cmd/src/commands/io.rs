//! File I/O commands: load, save, png, fetch, cd, pwd, ls

use std::env;
use std::path::Path;

use pymol_io::{read_file, write_file, FileFormat};
use pymol_scene::MoleculeObject;

use crate::args::{ArgDef, ParsedCommand};
use crate::command::{Command, CommandContext, CommandRegistry};
use crate::error::{CmdError, CmdResult};

/// Register I/O commands
pub fn register(registry: &mut CommandRegistry) {
    registry.register(LoadCommand);
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
"#
    }

    fn args(&self) -> &[ArgDef] {
        static ARGS: &[ArgDef] = &[];
        ARGS
    }

    fn execute(&self, ctx: &mut CommandContext, args: &ParsedCommand) -> CmdResult {
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
                Path::new(filename)
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("obj")
                    .to_string()
            });

        // Get state (optional, 0 = append)
        let _state = args
            .get_int(2)
            .or_else(|| args.get_named_int("state"))
            .unwrap_or(0);

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
                _ => FileFormat::Unknown,
            });

        // Get quiet flag
        let quiet = args
            .get_bool(4)
            .or_else(|| args.get_named_bool("quiet"))
            .unwrap_or(false);

        // Load the file
        let path = Path::new(filename);
        let mol = if let Some(fmt) = format {
            pymol_io::read_file_format(path, fmt)
        } else {
            read_file(path)
        }
        .map_err(|e| CmdError::FileFormat(e.to_string()))?;

        // Create molecule object and add to viewer
        let mol_obj = MoleculeObject::with_name(mol, &object_name);
        ctx.viewer.objects_mut().add(mol_obj);

        if !quiet {
            ctx.print(&format!(" Loaded \"{}\" as \"{}\"", filename, object_name));
        }

        // Center on loaded molecule
        ctx.viewer.center_on(&object_name);

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
"#
    }

    fn execute(&self, ctx: &mut CommandContext, args: &ParsedCommand) -> CmdResult {
        let filename = args
            .get_str(0)
            .or_else(|| args.get_named_str("filename"))
            .ok_or_else(|| CmdError::MissingArgument("filename".to_string()))?;

        let _selection = args
            .get_str(1)
            .or_else(|| args.get_named_str("selection"))
            .unwrap_or("all");

        let _state = args
            .get_int(2)
            .or_else(|| args.get_named_int("state"))
            .unwrap_or(-1);

        // For now, save the first molecule object
        // TODO: Implement proper selection support
        let path = Path::new(filename);

        // Find first molecule to save
        let mol = ctx
            .viewer
            .objects()
            .names()
            .find_map(|name| ctx.viewer.objects().get_molecule(name))
            .ok_or_else(|| CmdError::execution("No molecule objects to save"))?;

        write_file(path, mol.molecule()).map_err(|e| CmdError::FileFormat(e.to_string()))?;

        ctx.print(&format!(" Saved \"{}\"", filename));

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

    png filename [, width [, height [, dpi [, ray [, quiet]]]]]

ARGUMENTS

    filename = string: file path to be written
    width = integer: width in pixels (default: current window width)
    height = integer: height in pixels (default: current window height)
    dpi = float: dots-per-inch (not yet implemented)
    ray = 0/1: run ray tracing first (not yet implemented)
    quiet = 0/1: suppress feedback (default: 0)

NOTES

    If only width or height is specified, the aspect ratio is preserved.
    Ray tracing is not yet implemented in this version.

EXAMPLES

    png ~/Desktop/screenshot.png
    png output.png, 1920, 1080
    png output.png, width=800
    png hires.png, 4096, 4096, quiet=1
"#
    }

    fn execute(&self, ctx: &mut CommandContext, args: &ParsedCommand) -> CmdResult {
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

        // Get ray flag (optional, not yet implemented)
        let _ray = args
            .get_bool(4)
            .or_else(|| args.get_named_bool("ray"))
            .unwrap_or(false);

        // Get quiet flag
        let quiet = args
            .get_bool(5)
            .or_else(|| args.get_named_bool("quiet"))
            .unwrap_or(false);

        // Ensure the path ends with .png
        let path = Path::new(filename);
        let path = if path.extension().map(|e| e.to_ascii_lowercase()) != Some("png".into()) {
            path.with_extension("png")
        } else {
            path.to_path_buf()
        };

        // Capture the screenshot
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
                ctx.print(&format!(" Ray: render time --:--:--  ({}x{})", w, h));
            }
            ctx.print(&format!(" Saved \"{}\"", path.display()));
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

    fn execute(&self, ctx: &mut CommandContext, args: &ParsedCommand) -> CmdResult {
        let path = args.get_str(0).or_else(|| args.get_named_str("path"));

        let target = if let Some(p) = path {
            Path::new(p).to_path_buf()
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

    fn execute(&self, ctx: &mut CommandContext, _args: &ParsedCommand) -> CmdResult {
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

    fn execute(&self, ctx: &mut CommandContext, args: &ParsedCommand) -> CmdResult {
        let path = args
            .get_str(0)
            .or_else(|| args.get_named_str("path"))
            .unwrap_or(".");

        let dir = Path::new(path);

        if !dir.is_dir() {
            return Err(CmdError::execution(format!(
                "'{}' is not a directory",
                path
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
    type = string: file type to fetch (default: pdb)

EXAMPLES

    fetch 1ubq
    fetch 4hhb, name=hemoglobin
"#
    }

    fn execute(&self, ctx: &mut CommandContext, args: &ParsedCommand) -> CmdResult {
        let code = args
            .get_str(0)
            .or_else(|| args.get_named_str("code"))
            .ok_or_else(|| CmdError::MissingArgument("code".to_string()))?;

        let name = args
            .get_str(1)
            .or_else(|| args.get_named_str("name"))
            .unwrap_or(code);

        // Fetch the structure
        #[cfg(feature = "fetch")]
        let mol = pymol_io::fetch(code, pymol_io::FetchFormat::default())
            .map_err(|e| CmdError::FileFormat(e.to_string()))?;

        #[cfg(all(feature = "fetch-async", not(feature = "fetch")))]
        let mol = {
            // For async, we'd need a runtime - for now just error
            return Err(CmdError::execution(
                "Async fetch not supported in synchronous context",
            ));
        };

        // Add to viewer
        ctx.viewer
            .objects_mut()
            .add(MoleculeObject::with_name(mol, name));

        ctx.print(&format!(" Fetched {} as \"{}\"", code, name));

        // Center on fetched molecule
        ctx.viewer.center_on(name);

        Ok(())
    }
}
