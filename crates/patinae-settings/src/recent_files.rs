//! Persists the recent-files list.
//!
//! The store keeps only local filesystem paths. Loading is intentionally
//! forgiving: missing, malformed, or unsupported files return an empty list so
//! startup is not blocked by user-editable configuration.

use std::ffi::{OsStr, OsString};
use std::io;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

/// Maximum number of recent files kept on disk.
pub const RECENT_FILES_LIMIT: usize = 10;

const RECENT_FILES_VERSION: u32 = 1;

/// Ordered recent-files store.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RecentFiles {
    paths: Vec<PathBuf>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct RecentFilesDocument {
    version: u32,
    files: Vec<PathBuf>,
}

impl RecentFiles {
    /// Creates an empty recent-files store.
    pub fn new() -> Self {
        Self::default()
    }

    /// Loads a recent-files store from disk.
    ///
    /// Missing, malformed, or unsupported files return an empty store.
    pub fn load(path: &Path) -> Self {
        let Ok(bytes) = std::fs::read(path) else {
            return Self::new();
        };
        Self::from_json_slice(&bytes).unwrap_or_default()
    }

    /// Saves the recent-files store atomically.
    ///
    /// The file is written next to the destination and then renamed into place.
    ///
    /// # Errors
    /// Returns filesystem or serialization errors from directory creation,
    /// temporary-file writing, or final rename.
    pub fn save(&self, path: &Path) -> io::Result<()> {
        if let Some(parent) = path
            .parent()
            .filter(|parent| !parent.as_os_str().is_empty())
        {
            std::fs::create_dir_all(parent)?;
        }

        let temp_path = temporary_path(path);
        let bytes = serde_json::to_vec_pretty(&self.to_document()).map_err(io::Error::other)?;
        std::fs::write(&temp_path, bytes)?;
        if let Err(err) = std::fs::rename(&temp_path, path) {
            let _ = std::fs::remove_file(&temp_path);
            return Err(err);
        }
        Ok(())
    }

    /// Returns the recent paths newest first.
    pub fn paths(&self) -> &[PathBuf] {
        &self.paths
    }

    /// Records a path as the newest recent file.
    pub fn record(&mut self, path: PathBuf) {
        if path.as_os_str().is_empty() {
            return;
        }
        self.paths.retain(|existing| existing != &path);
        self.paths.insert(0, path);
        self.paths.truncate(RECENT_FILES_LIMIT);
    }

    /// Removes a path from the recent-files store.
    ///
    /// Returns `true` when an entry was removed.
    pub fn remove(&mut self, path: &Path) -> bool {
        let old_len = self.paths.len();
        self.paths.retain(|existing| existing != path);
        self.paths.len() != old_len
    }

    fn from_json_slice(bytes: &[u8]) -> Option<Self> {
        let document: RecentFilesDocument = serde_json::from_slice(bytes).ok()?;
        if document.version != RECENT_FILES_VERSION {
            return None;
        }

        let mut paths = Vec::new();
        for path in document.files {
            if path.as_os_str().is_empty() || paths.contains(&path) {
                continue;
            }
            paths.push(path);
            if paths.len() == RECENT_FILES_LIMIT {
                break;
            }
        }
        Some(Self { paths })
    }

    fn to_document(&self) -> RecentFilesDocument {
        RecentFilesDocument {
            version: RECENT_FILES_VERSION,
            files: self.paths.clone(),
        }
    }
}

fn temporary_path(path: &Path) -> PathBuf {
    let mut name = OsString::from(".");
    name.push(
        path.file_name()
            .unwrap_or_else(|| OsStr::new("recent-files.json")),
    );
    name.push(".tmp");
    path.with_file_name(name)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn temp_path(name: &str) -> PathBuf {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("system clock should be after Unix epoch")
            .as_nanos();
        std::env::temp_dir().join(format!("patinae-recent-files-{name}-{nonce}.json"))
    }

    #[test]
    fn record_dedupes_and_moves_to_front() {
        let mut recent = RecentFiles::new();

        recent.record(PathBuf::from("/tmp/a.pdb"));
        recent.record(PathBuf::from("/tmp/b.pdb"));
        recent.record(PathBuf::from("/tmp/a.pdb"));

        assert_eq!(
            recent.paths(),
            &[PathBuf::from("/tmp/a.pdb"), PathBuf::from("/tmp/b.pdb")]
        );
    }

    #[test]
    fn record_trims_to_limit() {
        let mut recent = RecentFiles::new();

        for idx in 0..(RECENT_FILES_LIMIT + 2) {
            recent.record(PathBuf::from(format!("/tmp/{idx}.pdb")));
        }

        assert_eq!(recent.paths().len(), RECENT_FILES_LIMIT);
        assert_eq!(recent.paths()[0], PathBuf::from("/tmp/11.pdb"));
        assert_eq!(
            recent.paths()[RECENT_FILES_LIMIT - 1],
            PathBuf::from("/tmp/2.pdb")
        );
    }

    #[test]
    fn missing_and_corrupt_files_load_empty() {
        let missing = temp_path("missing");
        assert!(RecentFiles::load(&missing).paths().is_empty());

        let corrupt = temp_path("corrupt");
        std::fs::write(&corrupt, b"{not json").expect("corrupt fixture should be writable");
        assert!(RecentFiles::load(&corrupt).paths().is_empty());
        let _ = std::fs::remove_file(corrupt);
    }

    #[test]
    fn save_and_load_roundtrip() {
        let path = temp_path("roundtrip");
        let mut recent = RecentFiles::new();
        recent.record(PathBuf::from("/tmp/a.pdb"));
        recent.record(PathBuf::from("/tmp/b.pdb"));

        recent.save(&path).expect("recent store should save");
        let restored = RecentFiles::load(&path);

        assert_eq!(restored, recent);
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn unsupported_version_loads_empty() {
        let path = temp_path("version");
        std::fs::write(&path, br#"{"version":999,"files":["/tmp/a.pdb"]}"#)
            .expect("version fixture should be writable");

        assert!(RecentFiles::load(&path).paths().is_empty());
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn remove_reports_whether_path_was_present() {
        let mut recent = RecentFiles::new();
        recent.record(PathBuf::from("/tmp/a.pdb"));

        assert!(recent.remove(Path::new("/tmp/a.pdb")));
        assert!(!recent.remove(Path::new("/tmp/a.pdb")));
        assert!(recent.paths().is_empty());
    }
}
