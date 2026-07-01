use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::{Path, PathBuf};

use thiserror::Error;

use crate::topology::{
    AngleType, AtomType, BondType, CombinationRule, ForceField, ForceFieldDefaults, ResidueAtom,
    ResidueTemplate, TerminalPatchSummary,
};

#[derive(Debug, Error)]
pub enum GromacsError {
    #[error("force-field path does not exist: {0}")]
    MissingPath(String),
    #[error("expected a .ff directory or forcefield.itp path, got: {0}")]
    InvalidRoot(String),
    #[error("I/O error reading {path}: {source}")]
    Io {
        path: String,
        #[source]
        source: std::io::Error,
    },
    #[error("unsupported GROMACS preprocessor directive in {path}:{line}: {directive}")]
    UnsupportedDirective {
        path: String,
        line: usize,
        directive: String,
    },
    #[error("unresolved include in {path}:{line}: {include}")]
    UnresolvedInclude {
        path: String,
        line: usize,
        include: String,
    },
    #[error("parse error in {path}:{line}: {message}")]
    Parse {
        path: String,
        line: usize,
        message: String,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TopSection {
    Defaults,
    AtomTypes,
    BondTypes,
    AngleTypes,
    NonbondParams,
    PairTypes,
    Other,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RtpSubsection {
    Atoms,
    Bonds,
    Exclusions,
    Impropers,
    Other,
}

#[derive(Debug, Default)]
struct ForceFieldBuilder {
    defaults: ForceFieldDefaults,
    atom_types: HashMap<String, AtomType>,
    bond_types: Vec<BondType>,
    angle_types: Vec<AngleType>,
}

#[derive(Debug, Clone, Copy)]
struct ConditionalFrame {
    parent_active: bool,
    condition_active: bool,
    else_seen: bool,
}

struct TopologyParseContext<'a> {
    root: &'a Path,
    seen: &'a mut HashSet<PathBuf>,
    defines: &'a mut HashSet<String>,
    builder: &'a mut ForceFieldBuilder,
}

pub fn load_force_field(path: impl AsRef<Path>) -> Result<ForceField, GromacsError> {
    let input = path.as_ref();
    if !input.exists() {
        return Err(GromacsError::MissingPath(input.display().to_string()));
    }

    let (root, forcefield_itp) = forcefield_root(input)?;
    let name = root
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("gromacs-force-field")
        .to_string();

    let mut builder = ForceFieldBuilder::default();
    let mut seen = HashSet::new();
    let mut defines = HashSet::new();
    let mut context = TopologyParseContext {
        root: &root,
        seen: &mut seen,
        defines: &mut defines,
        builder: &mut builder,
    };
    parse_topology_file(&forcefield_itp, &mut context)?;
    let residues = load_rtp_files(&root)?;
    let hydrogen_rules = load_hdb_files(&root)?;
    let terminal_patches = load_terminal_patches(&root)?;

    if residues.is_empty() {
        return Err(GromacsError::Parse {
            path: root.display().to_string(),
            line: 0,
            message: "no .rtp residue topology files found".into(),
        });
    }

    Ok(ForceField {
        name,
        defaults: builder.defaults,
        atom_types: builder.atom_types,
        bond_types: builder.bond_types,
        angle_types: builder.angle_types,
        residues,
        hydrogen_rules,
        terminal_patches,
    })
}

fn forcefield_root(path: &Path) -> Result<(PathBuf, PathBuf), GromacsError> {
    if path.is_dir() {
        let itp = path.join("forcefield.itp");
        if itp.exists() {
            return Ok((path.to_path_buf(), itp));
        }
        return Err(GromacsError::InvalidRoot(path.display().to_string()));
    }

    if path.file_name().and_then(|name| name.to_str()) == Some("forcefield.itp") {
        let root = path.parent().unwrap_or_else(|| Path::new("."));
        return Ok((root.to_path_buf(), path.to_path_buf()));
    }

    Err(GromacsError::InvalidRoot(path.display().to_string()))
}

fn parse_topology_file(
    path: &Path,
    context: &mut TopologyParseContext<'_>,
) -> Result<(), GromacsError> {
    let canonical = path.canonicalize().unwrap_or_else(|_| path.to_path_buf());
    if !context.seen.insert(canonical) {
        return Ok(());
    }

    let text = fs::read_to_string(path).map_err(|source| GromacsError::Io {
        path: path.display().to_string(),
        source,
    })?;

    let mut section = TopSection::Other;
    let mut conditionals = Vec::new();
    for (idx, raw_line) in text.lines().enumerate() {
        let line_no = idx + 1;
        let line = strip_comment(raw_line);
        if line.is_empty() {
            continue;
        }

        if line.starts_with('#') {
            handle_directive(path, line_no, line, context, &mut conditionals)?;
            continue;
        }

        if !conditionals_active(&conditionals) {
            continue;
        }

        if let Some(name) = bracket_name(line) {
            section = match name {
                "defaults" => TopSection::Defaults,
                "atomtypes" => TopSection::AtomTypes,
                "bondtypes" => TopSection::BondTypes,
                "angletypes" => TopSection::AngleTypes,
                "nonbond_params" => TopSection::NonbondParams,
                "pairtypes" => TopSection::PairTypes,
                _ => TopSection::Other,
            };
            continue;
        }

        let tokens: Vec<&str> = line.split_whitespace().collect();
        match section {
            TopSection::Defaults => parse_defaults(path, line_no, &tokens, context.builder)?,
            TopSection::AtomTypes => parse_atom_type(path, line_no, &tokens, context.builder)?,
            TopSection::BondTypes => parse_bond_type(&tokens, context.builder),
            TopSection::AngleTypes => parse_angle_type(&tokens, context.builder),
            TopSection::NonbondParams | TopSection::PairTypes | TopSection::Other => {}
        }
    }

    if !conditionals.is_empty() {
        return Err(parse_error(
            path,
            text.lines().count(),
            "unterminated #if block",
        ));
    }
    Ok(())
}

fn handle_directive(
    path: &Path,
    line_no: usize,
    line: &str,
    context: &mut TopologyParseContext<'_>,
    conditionals: &mut Vec<ConditionalFrame>,
) -> Result<(), GromacsError> {
    if let Some(name) = directive_arg(line, "#ifdef") {
        push_conditional(conditionals, context.defines.contains(name));
        return Ok(());
    }
    if let Some(name) = directive_arg(line, "#ifndef") {
        push_conditional(conditionals, !context.defines.contains(name));
        return Ok(());
    }
    if directive_exact(line, "#else") {
        flip_conditional(path, line_no, conditionals)?;
        return Ok(());
    }
    if directive_exact(line, "#endif") {
        pop_conditional(path, line_no, conditionals)?;
        return Ok(());
    }
    if !conditionals_active(conditionals) {
        return Ok(());
    }

    if let Some(include) = parse_include(line) {
        let include_path = resolve_include(path, context.root, include).ok_or_else(|| {
            GromacsError::UnresolvedInclude {
                path: path.display().to_string(),
                line: line_no,
                include: include.to_string(),
            }
        })?;
        parse_topology_file(&include_path, context)?;
        return Ok(());
    }

    if let Some(name) = directive_arg(line, "#define") {
        context.defines.insert(name.to_string());
        return Ok(());
    }

    Err(GromacsError::UnsupportedDirective {
        path: path.display().to_string(),
        line: line_no,
        directive: line.to_string(),
    })
}

fn conditionals_active(conditionals: &[ConditionalFrame]) -> bool {
    conditionals
        .last()
        .is_none_or(|frame| frame.parent_active && frame.condition_active)
}

fn push_conditional(conditionals: &mut Vec<ConditionalFrame>, condition_active: bool) {
    conditionals.push(ConditionalFrame {
        parent_active: conditionals_active(conditionals),
        condition_active,
        else_seen: false,
    });
}

fn flip_conditional(
    path: &Path,
    line_no: usize,
    conditionals: &mut [ConditionalFrame],
) -> Result<(), GromacsError> {
    let Some(frame) = conditionals.last_mut() else {
        return Err(parse_error(
            path,
            line_no,
            "#else without matching #ifdef/#ifndef",
        ));
    };
    if frame.else_seen {
        return Err(parse_error(
            path,
            line_no,
            "duplicate #else in conditional block",
        ));
    }
    frame.condition_active = !frame.condition_active;
    frame.else_seen = true;
    Ok(())
}

fn pop_conditional(
    path: &Path,
    line_no: usize,
    conditionals: &mut Vec<ConditionalFrame>,
) -> Result<(), GromacsError> {
    if conditionals.pop().is_none() {
        return Err(parse_error(
            path,
            line_no,
            "#endif without matching #ifdef/#ifndef",
        ));
    }
    Ok(())
}

fn directive_arg<'a>(line: &'a str, directive: &str) -> Option<&'a str> {
    let rest = line.strip_prefix(directive)?.trim();
    rest.split_whitespace().next()
}

fn directive_exact(line: &str, directive: &str) -> bool {
    line.strip_prefix(directive)
        .is_some_and(|rest| rest.trim().is_empty())
}

fn parse_include(line: &str) -> Option<&str> {
    let rest = line.strip_prefix("#include")?.trim();
    rest.strip_prefix('"')?.split('"').next()
}

fn resolve_include(current: &Path, root: &Path, include: &str) -> Option<PathBuf> {
    let direct = current.parent().unwrap_or(root).join(include);
    if direct.exists() {
        return Some(direct);
    }
    let rooted = root.join(include);
    rooted.exists().then_some(rooted)
}

fn parse_defaults(
    path: &Path,
    line: usize,
    tokens: &[&str],
    builder: &mut ForceFieldBuilder,
) -> Result<(), GromacsError> {
    if tokens.len() < 2 {
        return Err(parse_error(
            path,
            line,
            "expected at least two [ defaults ] fields",
        ));
    }
    let comb_rule = tokens
        .get(1)
        .and_then(|value| value.parse::<i32>().ok())
        .unwrap_or(2);
    let fudge_lj = tokens
        .get(3)
        .and_then(|value| value.parse::<f64>().ok())
        .unwrap_or(1.0);
    let fudge_qq = tokens
        .get(4)
        .and_then(|value| value.parse::<f64>().ok())
        .unwrap_or(1.0);

    builder.defaults = ForceFieldDefaults {
        combination_rule: CombinationRule::from_i32(comb_rule),
        fudge_lj,
        fudge_qq,
    };
    Ok(())
}

fn parse_atom_type(
    path: &Path,
    line: usize,
    tokens: &[&str],
    builder: &mut ForceFieldBuilder,
) -> Result<(), GromacsError> {
    if tokens.len() < 7 {
        return Err(parse_error(
            path,
            line,
            "expected at least seven [ atomtypes ] fields",
        ));
    }
    let name = tokens[0].to_string();
    let mass = parse_f64_at(path, line, tokens[tokens.len() - 5], "atom type mass")?;
    let _charge = parse_f64_at(path, line, tokens[tokens.len() - 4], "atom type charge")?;
    let v1 = parse_f64_at(
        path,
        line,
        tokens[tokens.len() - 2],
        "atom type V parameter",
    )?;
    let v2 = parse_f64_at(
        path,
        line,
        tokens[tokens.len() - 1],
        "atom type W parameter",
    )?;
    let (sigma_a, epsilon_kj) = match builder.defaults.combination_rule {
        CombinationRule::C6C12 => c6_c12_to_sigma_epsilon(v1, v2),
        CombinationRule::LorentzBerthelot | CombinationRule::Geometric => (v1 * 10.0, v2),
    };

    builder.atom_types.insert(
        name.clone(),
        AtomType {
            mass,
            sigma_a,
            epsilon_kj,
        },
    );
    Ok(())
}

fn parse_bond_type(tokens: &[&str], builder: &mut ForceFieldBuilder) {
    if tokens.len() < 5 {
        return;
    }
    let Some(length_nm) = tokens[tokens.len() - 2].parse::<f64>().ok() else {
        return;
    };
    let Some(force) = tokens[tokens.len() - 1].parse::<f64>().ok() else {
        return;
    };
    builder.bond_types.push(BondType {
        atoms: [tokens[0].to_string(), tokens[1].to_string()],
        length_nm,
        force_kj_mol_nm2: force,
    });
}

fn parse_angle_type(tokens: &[&str], builder: &mut ForceFieldBuilder) {
    if tokens.len() < 6 {
        return;
    }
    let Some(angle_deg) = tokens[tokens.len() - 2].parse::<f64>().ok() else {
        return;
    };
    let Some(force) = tokens[tokens.len() - 1].parse::<f64>().ok() else {
        return;
    };
    builder.angle_types.push(AngleType {
        atoms: [
            tokens[0].to_string(),
            tokens[1].to_string(),
            tokens[2].to_string(),
        ],
        angle_deg,
        force_kj_mol_rad2: force,
    });
}

fn c6_c12_to_sigma_epsilon(c6: f64, c12: f64) -> (f64, f64) {
    if c6 <= 0.0 || c12 <= 0.0 {
        return (0.0, 0.0);
    }
    let sigma_nm = (c12 / c6).powf(1.0 / 6.0);
    let epsilon = c6 * c6 / (4.0 * c12);
    (sigma_nm * 10.0, epsilon)
}

fn load_rtp_files(root: &Path) -> Result<HashMap<String, ResidueTemplate>, GromacsError> {
    let mut residues = HashMap::new();
    for entry in fs::read_dir(root).map_err(|source| GromacsError::Io {
        path: root.display().to_string(),
        source,
    })? {
        let entry = entry.map_err(|source| GromacsError::Io {
            path: root.display().to_string(),
            source,
        })?;
        let path = entry.path();
        if path.extension().and_then(|ext| ext.to_str()) == Some("rtp") {
            parse_rtp_file(&path, &mut residues)?;
        }
    }
    Ok(residues)
}

fn parse_rtp_file(
    path: &Path,
    residues: &mut HashMap<String, ResidueTemplate>,
) -> Result<(), GromacsError> {
    let text = fs::read_to_string(path).map_err(|source| GromacsError::Io {
        path: path.display().to_string(),
        source,
    })?;

    let mut current: Option<ResidueTemplate> = None;
    let mut subsection = RtpSubsection::Other;

    for (idx, raw_line) in text.lines().enumerate() {
        let line_no = idx + 1;
        let line = strip_comment(raw_line);
        if line.is_empty() {
            continue;
        }
        if line.starts_with('#') {
            return Err(GromacsError::UnsupportedDirective {
                path: path.display().to_string(),
                line: line_no,
                directive: line.to_string(),
            });
        }

        if let Some(name) = bracket_name(line) {
            match name {
                "bondedtypes" => {
                    commit_residue(&mut current, residues);
                    subsection = RtpSubsection::Other;
                }
                "atoms" if current.is_some() => subsection = RtpSubsection::Atoms,
                "bonds" if current.is_some() => subsection = RtpSubsection::Bonds,
                "exclusions" if current.is_some() => subsection = RtpSubsection::Exclusions,
                "impropers" if current.is_some() => subsection = RtpSubsection::Impropers,
                "angles" | "dihedrals" if current.is_some() => subsection = RtpSubsection::Other,
                _ => {
                    commit_residue(&mut current, residues);
                    current = Some(ResidueTemplate {
                        name: name.to_string(),
                        atoms: Vec::new(),
                        bonds: Vec::new(),
                        exclusions: Vec::new(),
                        impropers: Vec::new(),
                    });
                    subsection = RtpSubsection::Other;
                }
            }
            continue;
        }

        let tokens: Vec<&str> = line.split_whitespace().collect();
        let Some(residue) = current.as_mut() else {
            continue;
        };
        match subsection {
            RtpSubsection::Atoms => parse_rtp_atom(path, line_no, &tokens, residue)?,
            RtpSubsection::Bonds => push_pair(&tokens, &mut residue.bonds),
            RtpSubsection::Exclusions => push_pair(&tokens, &mut residue.exclusions),
            RtpSubsection::Impropers => push_quad(&tokens, &mut residue.impropers),
            RtpSubsection::Other => {}
        }
    }

    commit_residue(&mut current, residues);
    Ok(())
}

fn parse_rtp_atom(
    path: &Path,
    line: usize,
    tokens: &[&str],
    residue: &mut ResidueTemplate,
) -> Result<(), GromacsError> {
    if tokens.len() < 3 {
        return Err(parse_error(
            path,
            line,
            "expected RTP atom name, type, and charge",
        ));
    }
    residue.atoms.push(ResidueAtom {
        name: tokens[0].to_string(),
        atom_type: tokens[1].to_string(),
        charge: parse_f64_at(path, line, tokens[2], "RTP atom charge")?,
    });
    Ok(())
}

fn push_pair(tokens: &[&str], pairs: &mut Vec<[String; 2]>) {
    if tokens.len() >= 2 && is_local_atom(tokens[0]) && is_local_atom(tokens[1]) {
        pairs.push([tokens[0].to_string(), tokens[1].to_string()]);
    }
}

fn push_quad(tokens: &[&str], quads: &mut Vec<[String; 4]>) {
    if tokens.len() >= 4 && tokens.iter().take(4).all(|atom| is_local_atom(atom)) {
        quads.push([
            tokens[0].to_string(),
            tokens[1].to_string(),
            tokens[2].to_string(),
            tokens[3].to_string(),
        ]);
    }
}

fn is_local_atom(name: &str) -> bool {
    !name.starts_with('+') && !name.starts_with('-')
}

fn commit_residue(
    current: &mut Option<ResidueTemplate>,
    residues: &mut HashMap<String, ResidueTemplate>,
) {
    if let Some(residue) = current.take() {
        if !residue.atoms.is_empty() {
            residues.insert(residue.name.clone(), residue);
        }
    }
}

fn load_terminal_patches(root: &Path) -> Result<TerminalPatchSummary, GromacsError> {
    let mut summary = TerminalPatchSummary::default();
    for entry in fs::read_dir(root).map_err(|source| GromacsError::Io {
        path: root.display().to_string(),
        source,
    })? {
        let entry = entry.map_err(|source| GromacsError::Io {
            path: root.display().to_string(),
            source,
        })?;
        let path = entry.path();
        let Some(name) = path.file_name().and_then(|name| name.to_str()) else {
            continue;
        };
        if name.ends_with(".n.tdb") {
            summary.n_terminal_files += 1;
            collect_tdb_entries(&path, &mut summary.entries)?;
        } else if name.ends_with(".c.tdb") {
            summary.c_terminal_files += 1;
            collect_tdb_entries(&path, &mut summary.entries)?;
        }
    }
    summary.entries.sort();
    summary.entries.dedup();
    Ok(summary)
}

fn load_hdb_files(
    root: &Path,
) -> Result<HashMap<String, Vec<crate::topology::HydrogenRule>>, GromacsError> {
    let mut rules = HashMap::new();
    for entry in fs::read_dir(root).map_err(|source| GromacsError::Io {
        path: root.display().to_string(),
        source,
    })? {
        let entry = entry.map_err(|source| GromacsError::Io {
            path: root.display().to_string(),
            source,
        })?;
        let path = entry.path();
        if path.extension().and_then(|ext| ext.to_str()) == Some("hdb") {
            parse_hdb_file(&path, &mut rules)?;
        }
    }
    Ok(rules)
}

fn parse_hdb_file(
    path: &Path,
    rules: &mut HashMap<String, Vec<crate::topology::HydrogenRule>>,
) -> Result<(), GromacsError> {
    let text = fs::read_to_string(path).map_err(|source| GromacsError::Io {
        path: path.display().to_string(),
        source,
    })?;

    let mut current_residue = String::new();
    let mut remaining = 0usize;
    for (idx, raw_line) in text.lines().enumerate() {
        let line_no = idx + 1;
        let line = strip_comment(raw_line);
        if line.is_empty() {
            continue;
        }
        if line.starts_with('#') {
            return Err(GromacsError::UnsupportedDirective {
                path: path.display().to_string(),
                line: line_no,
                directive: line.to_string(),
            });
        }
        let tokens: Vec<&str> = line.split_whitespace().collect();
        if tokens.len() == 2
            && tokens[1].parse::<usize>().is_ok()
            && tokens[0].parse::<usize>().is_err()
        {
            current_residue = tokens[0].to_string();
            remaining = tokens[1]
                .parse::<usize>()
                .map_err(|_| parse_error(path, line_no, "invalid HDB residue entry count"))?;
            rules.entry(current_residue.clone()).or_default();
            continue;
        }
        if current_residue.is_empty() || remaining == 0 {
            return Err(parse_error(
                path,
                line_no,
                "HDB hydrogen entry appears outside a residue block",
            ));
        }
        if tokens.len() < 4 {
            return Err(parse_error(
                path,
                line_no,
                "expected HDB count, method, name, and control atoms",
            ));
        }
        let count = tokens[0]
            .parse::<usize>()
            .map_err(|_| parse_error(path, line_no, "invalid HDB hydrogen count"))?;
        let method = tokens[1]
            .parse::<i32>()
            .map_err(|_| parse_error(path, line_no, "invalid HDB hydrogen method"))?;
        rules
            .entry(current_residue.clone())
            .or_default()
            .push(crate::topology::HydrogenRule {
                count,
                method,
                name: tokens[2].to_string(),
                controls: tokens[3..]
                    .iter()
                    .map(|token| (*token).to_string())
                    .collect(),
            });
        remaining -= 1;
    }
    Ok(())
}

fn collect_tdb_entries(path: &Path, entries: &mut Vec<String>) -> Result<(), GromacsError> {
    let text = fs::read_to_string(path).map_err(|source| GromacsError::Io {
        path: path.display().to_string(),
        source,
    })?;
    for raw_line in text.lines() {
        let line = strip_comment(raw_line);
        if let Some(name) = bracket_name(line) {
            entries.push(name.to_string());
        }
    }
    Ok(())
}

fn strip_comment(line: &str) -> &str {
    line.split(';').next().unwrap_or("").trim()
}

fn bracket_name(line: &str) -> Option<&str> {
    line.strip_prefix('[')?.split(']').next().map(str::trim)
}

fn parse_f64_at(path: &Path, line: usize, value: &str, label: &str) -> Result<f64, GromacsError> {
    value
        .parse::<f64>()
        .map_err(|_| parse_error(path, line, format!("invalid {label}: {value}")))
}

fn parse_error(path: &Path, line: usize, message: impl Into<String>) -> GromacsError {
    GromacsError::Parse {
        path: path.display().to_string(),
        line,
        message: message.into(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_file(path: &Path, text: &str) {
        let mut file = fs::File::create(path).unwrap();
        file.write_all(text.as_bytes()).unwrap();
    }

    #[test]
    fn imports_forcefield_includes_atomtypes_and_rtp() {
        let dir = tempfile::tempdir().unwrap();
        let ff = dir.path().join("tiny.ff");
        fs::create_dir(&ff).unwrap();
        write_file(
            &ff.join("forcefield.itp"),
            r#"
#define _FF_TINY
[ defaults ]
1 2 yes 0.5 0.8333
#include "ffnonbonded.itp"
#include "ffbonded.itp"
"#,
        );
        write_file(
            &ff.join("ffnonbonded.itp"),
            r#"
[ atomtypes ]
C 6 12.011 0.0 A 0.34 0.45
H 1 1.008 0.0 A 0.25 0.10
"#,
        );
        write_file(
            &ff.join("ffbonded.itp"),
            r#"
[ bondtypes ]
C H 1 0.109 284512
[ angletypes ]
H C H 1 109.5 276.144
"#,
        );
        write_file(
            &ff.join("aminoacids.rtp"),
            r#"
[ bondedtypes ]
1 1 1 2
[ LIG ]
 [ atoms ]
 C1 C 0.1 0
 H1 H -0.1 0
 [ bonds ]
 C1 H1
"#,
        );
        write_file(
            &ff.join("aminoacids.hdb"),
            r#"
LIG 1
1 1 H1 C1
"#,
        );

        let loaded = load_force_field(&ff).unwrap();
        assert_eq!(loaded.defaults.fudge_qq, 0.8333);
        assert!(loaded.atom_type("C").is_some());
        assert_eq!(loaded.bond_types.len(), 1);
        assert_eq!(loaded.angle_types.len(), 1);
        assert_eq!(loaded.residue("LIG").unwrap().atoms.len(), 2);
        assert_eq!(loaded.hydrogen_rules["LIG"][0].name, "H1");
    }

    #[test]
    fn rejects_unresolved_include() {
        let dir = tempfile::tempdir().unwrap();
        let ff = dir.path().join("bad.ff");
        fs::create_dir(&ff).unwrap();
        write_file(&ff.join("forcefield.itp"), "#include \"missing.itp\"\n");

        let err = load_force_field(&ff).unwrap_err();
        assert!(err.to_string().contains("unresolved include"));
    }

    #[test]
    fn honors_supported_conditionals() {
        let dir = tempfile::tempdir().unwrap();
        let ff = dir.path().join("conditional.ff");
        fs::create_dir(&ff).unwrap();
        write_file(
            &ff.join("forcefield.itp"),
            r#"
#define USE_REAL
[ defaults ]
1 2 yes 1.0 1.0
#ifdef HEAVY_H
#include "missing-heavy.itp"
#else
#include "ffnonbonded.itp"
#endif
#ifndef USE_REAL
#include "missing-disabled.itp"
#endif
"#,
        );
        write_file(
            &ff.join("ffnonbonded.itp"),
            r#"
[ atomtypes ]
C 6 12.011 0.0 A 0.34 0.45
"#,
        );
        write_file(
            &ff.join("aminoacids.rtp"),
            r#"
[ LIG ]
 [ atoms ]
 C1 C 0.0
"#,
        );

        let loaded = load_force_field(&ff).unwrap();
        assert!(loaded.atom_type("C").is_some());
    }

    #[test]
    fn rejects_unsupported_active_preprocessor_directive() {
        let dir = tempfile::tempdir().unwrap();
        let ff = dir.path().join("bad.ff");
        fs::create_dir(&ff).unwrap();
        write_file(&ff.join("forcefield.itp"), "#if 1\n");

        let err = load_force_field(&ff).unwrap_err();
        assert!(err.to_string().contains("unsupported"));
    }

    #[test]
    fn loads_bundled_forcefields() {
        for alias in ["AMBER", "CHARMM", "OPLS-AA"] {
            let path = crate::forcefields::resolve_force_field_path(alias, &[]).unwrap();
            let loaded = load_force_field(path).unwrap();
            assert!(
                !loaded.residues.is_empty(),
                "{alias} did not load residue templates"
            );
            assert!(
                !loaded.atom_types.is_empty(),
                "{alias} did not load atom types"
            );
            assert!(
                !loaded.hydrogen_rules.is_empty(),
                "{alias} did not load hydrogen rules"
            );
        }
    }
}
