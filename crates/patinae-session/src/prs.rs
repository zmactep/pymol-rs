//! Native PRS session format: MessagePack + gzip.

use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;

use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use patinae_scene::Session;
use serde::{Deserialize, Serialize};

/// Current native PRS document format version.
pub const PRS_FORMAT_VERSION: u32 = 2;

/// Format version assigned to legacy raw [`Session`] files.
pub const PRS_LEGACY_FORMAT_VERSION: u32 = 1;

/// Native PRS producer name.
pub const PRS_PRODUCER: &str = "patinae";

/// Native PRS producer crate version.
pub const PRS_PRODUCER_VERSION: &str = env!("CARGO_PKG_VERSION");

/// A loaded PRS document plus its format metadata.
#[derive(Serialize, Deserialize)]
pub struct PrsDocument {
    /// PRS envelope format version.
    pub prs_format_version: u32,
    /// Producer name. Missing for legacy raw `Session` files.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub producer: Option<String>,
    /// Producer version. Missing for legacy raw `Session` files.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub producer_version: Option<String>,
    /// Scene session payload.
    pub session: Session,
}

impl std::fmt::Debug for PrsDocument {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PrsDocument")
            .field("prs_format_version", &self.prs_format_version)
            .field("producer", &self.producer)
            .field("producer_version", &self.producer_version)
            .field("session", &"<Session>")
            .finish()
    }
}

impl PrsDocument {
    /// Return user-facing PRS compatibility warnings.
    pub fn warning_messages(&self) -> Vec<String> {
        prs_document_warning_messages(self)
    }

    fn legacy(session: Session) -> Self {
        Self {
            prs_format_version: PRS_LEGACY_FORMAT_VERSION,
            producer: None,
            producer_version: None,
            session,
        }
    }
}

#[derive(Serialize)]
struct PrsDocumentRef<'a> {
    prs_format_version: u32,
    producer: &'static str,
    producer_version: &'static str,
    session: &'a Session,
}

/// Errors for PRS save/load operations.
#[derive(Debug, thiserror::Error)]
pub enum PrsError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("serialization error: {0}")]
    Serialize(#[from] rmp_serde::encode::Error),

    #[error("deserialization error: {0}")]
    Deserialize(#[from] rmp_serde::decode::Error),
}

/// Save a session to a `.prs` file (MessagePack + gzip).
pub fn save_prs(session: &Session, path: &Path) -> Result<(), PrsError> {
    let document = PrsDocumentRef {
        prs_format_version: PRS_FORMAT_VERSION,
        producer: PRS_PRODUCER,
        producer_version: PRS_PRODUCER_VERSION,
        session,
    };
    let data = rmp_serde::to_vec_named(&document)?;
    write_prs_bytes(path, &data)?;
    Ok(())
}

/// Load a session from a `.prs` file (MessagePack + gzip).
pub fn load_prs(path: &Path) -> Result<Session, PrsError> {
    let document = load_prs_document(path)?;
    log_prs_document_warnings(&document);
    Ok(document.session)
}

/// Load a `.prs` file together with PRS format and producer metadata.
pub fn load_prs_document(path: &Path) -> Result<PrsDocument, PrsError> {
    let bytes = read_prs_bytes(path)?;
    decode_prs_document(&bytes)
}

fn decode_prs_document(bytes: &[u8]) -> Result<PrsDocument, PrsError> {
    match rmp_serde::from_slice::<PrsDocument>(bytes) {
        Ok(document) => Ok(document),
        Err(envelope_error) => match rmp_serde::from_slice::<Session>(bytes) {
            Ok(session) => Ok(PrsDocument::legacy(session)),
            Err(_) => Err(PrsError::Deserialize(envelope_error)),
        },
    }
}

fn log_prs_document_warnings(document: &PrsDocument) {
    for warning in document.warning_messages() {
        log::warn!("{}", warning);
    }
}

fn prs_document_warning_messages(document: &PrsDocument) -> Vec<String> {
    let mut warnings = Vec::new();

    if document.prs_format_version == PRS_LEGACY_FORMAT_VERSION
        && document.producer.is_none()
        && document.producer_version.is_none()
    {
        warnings.push(format!(
            "Loaded legacy PRS session (format v{}). Re-save it with this Patinae version to upgrade PRS metadata.",
            PRS_LEGACY_FORMAT_VERSION,
        ));
    }

    let Some(producer_version) = document.producer_version.as_deref() else {
        return warnings;
    };
    if document.producer.as_deref() == Some(PRS_PRODUCER)
        && version_is_newer(producer_version, PRS_PRODUCER_VERSION)
    {
        warnings.push(format!(
            "Loaded PRS session produced by newer Patinae version {}. Current Patinae version is {}; some session data may not be interpreted exactly as saved.",
            producer_version,
            PRS_PRODUCER_VERSION,
        ));
    }

    warnings
}

fn version_is_newer(candidate: &str, current: &str) -> bool {
    match (semver_core(candidate), semver_core(current)) {
        (Some(candidate), Some(current)) => candidate > current,
        _ => false,
    }
}

fn semver_core(version: &str) -> Option<[u64; 3]> {
    let core = version
        .split_once(['-', '+'])
        .map_or(version, |(core, _)| core);
    let mut parsed = [0; 3];
    let mut parts = core.split('.');

    for slot in &mut parsed {
        let Some(part) = parts.next() else {
            break;
        };
        if part.is_empty() {
            return None;
        }
        *slot = part.parse().ok()?;
    }

    if parts.next().is_some() {
        return None;
    }

    Some(parsed)
}

fn write_prs_bytes(path: &Path, data: &[u8]) -> Result<(), PrsError> {
    let file = File::create(path)?;
    let mut encoder = GzEncoder::new(file, Compression::default());
    encoder.write_all(data)?;
    encoder.finish()?;
    Ok(())
}

fn read_prs_bytes(path: &Path) -> Result<Vec<u8>, PrsError> {
    let file = File::open(path)?;
    let mut decoder = GzDecoder::new(file);
    let mut bytes = Vec::new();
    decoder.read_to_end(&mut bytes)?;
    Ok(bytes)
}

#[cfg(test)]
mod tests {
    use super::*;
    use lin_alg::f32::{Mat4, Vec3};
    use patinae_color::ColorIndex;
    use patinae_mol::{Atom, CoordSet, Element, ObjectMolecule, RepMask};
    use patinae_scene::{
        Camera, GroupObject, MoleculeObject, Object, SceneManager, SelectionManager, ViewManager,
    };
    use patinae_settings::{ObjectOverrides, Settings};
    use serde::Serialize;
    use std::fs;
    use std::path::PathBuf;

    fn temp_prs_path(name: &str) -> (PathBuf, PathBuf) {
        let dir =
            std::env::temp_dir().join(format!("patinae_prs_test_{}_{}", name, std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let path = dir.join("test.prs");
        (dir, path)
    }

    fn cartoon_session_with_restore() -> Session {
        let mut session = Session::new();
        let mut mol = ObjectMolecule::new("mol");
        mol.add_atom(Atom::new("CA", Element::Carbon));
        mol.add_coord_set(CoordSet::from_vec3(&[Vec3::new(0.0, 0.0, 0.0)]));
        let mut obj = MoleculeObject::with_name(mol, "mol");
        obj.state_mut().hide_draw_rep(RepMask::CARTOON);
        session.registry.add(obj);
        session
    }

    #[test]
    fn test_prs_round_trip() {
        let session = cartoon_session_with_restore();
        let (dir, path) = temp_prs_path("roundtrip");

        save_prs(&session, &path).unwrap();
        assert!(path.exists());

        let document = load_prs_document(&path).unwrap();
        assert_eq!(document.prs_format_version, PRS_FORMAT_VERSION);
        assert_eq!(document.producer.as_deref(), Some(PRS_PRODUCER));
        assert_eq!(
            document.producer_version.as_deref(),
            Some(PRS_PRODUCER_VERSION)
        );
        assert!(document
            .session
            .registry
            .get_molecule("mol")
            .unwrap()
            .draw_mask_restorable_reps()
            .is_visible(RepMask::CARTOON));

        let loaded = load_prs(&path).unwrap();
        assert_eq!(loaded.clear_color, session.clear_color);
        assert!(loaded
            .registry
            .get_molecule("mol")
            .unwrap()
            .draw_mask_restorable_reps()
            .is_visible(RepMask::CARTOON));

        // Cleanup
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn legacy_raw_session_loads_as_format_v1_without_producer_metadata() {
        let session = Session::new();
        let (dir, path) = temp_prs_path("legacy_empty");
        let data = rmp_serde::to_vec(&session).unwrap();
        write_prs_bytes(&path, &data).unwrap();

        let document = load_prs_document(&path).unwrap();

        assert_eq!(document.prs_format_version, PRS_LEGACY_FORMAT_VERSION);
        assert_eq!(document.producer, None);
        assert_eq!(document.producer_version, None);
        assert_eq!(document.session.clear_color, session.clear_color);
        assert!(document
            .warning_messages()
            .iter()
            .any(|warning| warning.contains("legacy PRS session")));

        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn producer_version_warning_detects_newer_semver_core() {
        assert!(version_is_newer("0.10.0", "0.9.9"));
        assert!(version_is_newer("1.0.0-alpha.1", "0.9.9"));
        assert!(!version_is_newer("0.4.1", "0.4.1"));
        assert!(!version_is_newer("0.4.1+local", "0.4.1"));
        assert!(!version_is_newer("0.4.0", "0.4.1"));
        assert!(!version_is_newer("not-a-version", "0.4.1"));
    }

    #[test]
    fn warning_messages_detect_newer_producer_version() {
        let document = PrsDocument {
            prs_format_version: PRS_FORMAT_VERSION,
            producer: Some(PRS_PRODUCER.to_string()),
            producer_version: Some("999.0.0".to_string()),
            session: Session::new(),
        };

        assert!(document.warning_messages().iter().any(|warning| {
            warning.contains("newer Patinae version 999.0.0")
                && warning.contains(PRS_PRODUCER_VERSION)
        }));
    }

    #[derive(Serialize)]
    struct LegacyObjectState {
        enabled: bool,
        color: ColorIndex,
        visible_reps: RepMask,
        draw_reps: Option<RepMask>,
        #[serde(with = "patinae_scene::serde_helpers::mat4_serde")]
        transform: Mat4,
    }

    #[derive(Serialize)]
    struct LegacyMoleculeObjectSnapshot {
        molecule: ObjectMolecule,
        state: LegacyObjectState,
        display_state: usize,
        overrides: Option<ObjectOverrides>,
        surface_quality: i32,
    }

    #[derive(Serialize)]
    struct LegacyObjectRegistrySnapshot {
        molecules: Vec<(String, LegacyMoleculeObjectSnapshot)>,
        groups: Vec<(String, GroupObject)>,
        render_order: Vec<String>,
        object_states: Vec<(String, LegacyObjectState)>,
        render_ids: Vec<(String, u32)>,
        next_render_id: u32,
        next_id: u32,
        generation: u64,
    }

    #[derive(Serialize)]
    struct LegacySessionRef<'a> {
        registry: LegacyObjectRegistrySnapshot,
        camera: &'a Camera,
        selections: &'a SelectionManager,
        scenes: &'a SceneManager,
        views: &'a ViewManager,
        movie: &'a patinae_scene::Movie,
        settings: &'a Settings,
        named_palette: &'a patinae_scene::NamedPalette,
        palette: &'a patinae_scene::ThemedPalette,
        clear_color: [f32; 3],
        clear_color_set: bool,
    }

    fn legacy_cartoon_state() -> LegacyObjectState {
        LegacyObjectState {
            enabled: true,
            color: ColorIndex::default(),
            visible_reps: RepMask::CARTOON,
            draw_reps: Some(RepMask::NONE),
            transform: Mat4::new_identity(),
        }
    }

    #[test]
    fn legacy_object_state_sequence_defaults_restore_mask_to_none() {
        let data = rmp_serde::to_vec(&legacy_cartoon_state()).unwrap();
        let state: patinae_scene::ObjectState = rmp_serde::from_slice(&data).unwrap();

        assert!(state.visible_reps.is_visible(RepMask::CARTOON));
        assert!(!state.draw_reps.is_visible(RepMask::CARTOON));
        assert_eq!(state.draw_mask_restorable_reps, RepMask::NONE);
    }

    #[test]
    fn legacy_raw_session_object_state_defaults_restore_mask_to_none() {
        let mut mol = ObjectMolecule::new("legacy");
        mol.add_atom(Atom::new("CA", Element::Carbon));
        mol.add_coord_set(CoordSet::from_vec3(&[Vec3::new(0.0, 0.0, 0.0)]));
        let legacy_registry = LegacyObjectRegistrySnapshot {
            molecules: vec![(
                "legacy".to_string(),
                LegacyMoleculeObjectSnapshot {
                    molecule: mol,
                    state: legacy_cartoon_state(),
                    display_state: 0,
                    overrides: None,
                    surface_quality: 0,
                },
            )],
            groups: Vec::new(),
            render_order: vec!["legacy".to_string()],
            object_states: vec![("legacy".to_string(), legacy_cartoon_state())],
            render_ids: Vec::new(),
            next_render_id: 1,
            next_id: 1,
            generation: 0,
        };
        let base = Session::new();
        let legacy_session = LegacySessionRef {
            registry: legacy_registry,
            camera: &base.camera,
            selections: &base.selections,
            scenes: &base.scenes,
            views: &base.views,
            movie: &base.movie,
            settings: &base.settings,
            named_palette: &base.named_palette,
            palette: &base.palette,
            clear_color: base.clear_color,
            clear_color_set: base.clear_color_set,
        };
        let (dir, path) = temp_prs_path("legacy_restore_flag");
        let data = rmp_serde::to_vec_named(&legacy_session).unwrap();
        write_prs_bytes(&path, &data).unwrap();

        let document = load_prs_document(&path).unwrap();
        let obj = document.session.registry.get_molecule("legacy").unwrap();

        assert_eq!(document.prs_format_version, PRS_LEGACY_FORMAT_VERSION);
        assert!(obj.visible_reps().is_visible(RepMask::CARTOON));
        assert!(!obj.draw_reps().is_visible(RepMask::CARTOON));
        assert_eq!(obj.draw_mask_restorable_reps(), RepMask::NONE);

        let _ = fs::remove_dir_all(&dir);
    }
}
