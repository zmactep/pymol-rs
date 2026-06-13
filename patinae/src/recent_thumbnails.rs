use std::ffi::{OsStr, OsString};
use std::io;
use std::path::{Path, PathBuf};
use std::time::UNIX_EPOCH;

use slint::Image;

/// Width of a cached recent-card thumbnail in pixels.
pub(crate) const RECENT_THUMBNAIL_WIDTH: u32 = 88;
/// Height of a cached recent-card thumbnail in pixels.
///
/// Keep recent captures square: wide captures make molecular previews look
/// vertically compressed inside the compact recent-card artwork slot.
pub(crate) const RECENT_THUMBNAIL_HEIGHT: u32 = 88;

/// Filesystem-backed cache for recent-file thumbnails.
#[derive(Debug, Clone)]
pub(crate) struct RecentThumbnailCache {
    dir: PathBuf,
}

impl RecentThumbnailCache {
    pub(crate) fn new(dir: PathBuf) -> Self {
        Self { dir }
    }

    pub(crate) fn load(&self, source: &Path) -> Option<Image> {
        let path = self.thumbnail_path(source);
        if !path.is_file() {
            return None;
        }
        Image::load_from_path(&path).ok()
    }

    pub(crate) fn write_paths(&self, source: &Path) -> io::Result<ThumbnailWritePaths> {
        std::fs::create_dir_all(&self.dir)?;
        let path = self.thumbnail_path(source);
        let temp_path = temporary_path(&path);
        Ok(ThumbnailWritePaths { path, temp_path })
    }

    fn thumbnail_path(&self, source: &Path) -> PathBuf {
        self.dir.join(format!("{}.png", thumbnail_key(source)))
    }
}

/// Final and temporary paths for one atomic thumbnail write.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ThumbnailWritePaths {
    pub(crate) path: PathBuf,
    pub(crate) temp_path: PathBuf,
}

impl ThumbnailWritePaths {
    pub(crate) fn commit(&self) -> io::Result<()> {
        std::fs::rename(&self.temp_path, &self.path)
    }

    pub(crate) fn cleanup_temp(&self) {
        let _ = std::fs::remove_file(&self.temp_path);
    }
}

pub(crate) fn extension_badge(path: &Path) -> String {
    let Some(extension) = normalized_extension(path) else {
        return String::new();
    };
    match extension.as_str() {
        "pdb" | "cif" => String::new(),
        _ => format!(".{extension}"),
    }
}

fn normalized_extension(path: &Path) -> Option<String> {
    path.extension()
        .and_then(OsStr::to_str)
        .map(str::trim)
        .filter(|extension| !extension.is_empty())
        .map(|extension| extension.to_ascii_lowercase())
}

fn thumbnail_key(source: &Path) -> String {
    let identity = thumbnail_identity(source);
    format!("{:016x}", fnv1a64(identity.as_bytes()))
}

fn thumbnail_identity(source: &Path) -> String {
    let resolved = source
        .canonicalize()
        .unwrap_or_else(|_| source.to_path_buf());
    let mut identity = resolved.to_string_lossy().into_owned();
    match std::fs::metadata(&resolved).or_else(|_| std::fs::metadata(source)) {
        Ok(metadata) => {
            identity.push_str("|len=");
            identity.push_str(&metadata.len().to_string());
            if let Ok(modified) = metadata.modified() {
                identity.push_str("|mtime-ns=");
                let nanos = modified
                    .duration_since(UNIX_EPOCH)
                    .map(|duration| duration.as_nanos())
                    .unwrap_or_default();
                identity.push_str(&nanos.to_string());
            }
        }
        Err(_) => identity.push_str("|missing"),
    }
    identity
}

fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut hash = 0xcbf2_9ce4_8422_2325u64;
    for byte in bytes {
        hash ^= u64::from(*byte);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    hash
}

fn temporary_path(path: &Path) -> PathBuf {
    let mut name = OsString::from(".");
    name.push(
        path.file_name()
            .unwrap_or_else(|| OsStr::new("recent-thumbnail.png")),
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
        std::env::temp_dir().join(format!("patinae-thumbnail-{name}-{nonce}"))
    }

    #[test]
    fn extension_badge_hides_default_structure_formats() {
        assert_eq!(extension_badge(Path::new("/tmp/a.pdb")), "");
        assert_eq!(extension_badge(Path::new("/tmp/a.CIF")), "");
    }

    #[test]
    fn extension_badge_keeps_session_and_other_formats() {
        assert_eq!(extension_badge(Path::new("/tmp/a.pse")), ".pse");
        assert_eq!(extension_badge(Path::new("/tmp/a.mol2")), ".mol2");
    }

    #[test]
    fn write_paths_use_recent_cache_directory_and_temp_sibling() {
        let root = temp_path("write-paths-cache");
        let dir = root.join("recent");
        let cache = RecentThumbnailCache::new(dir.clone());
        let paths = cache
            .write_paths(Path::new("/tmp/model.pdb"))
            .expect("temp cache directory should be creatable");

        assert_eq!(paths.path.parent(), Some(dir.as_path()));
        assert_eq!(paths.temp_path.parent(), Some(dir.as_path()));
        assert_eq!(paths.path.extension().and_then(OsStr::to_str), Some("png"));

        let _ = std::fs::remove_dir_all(root);
    }

    #[test]
    fn thumbnail_key_changes_when_file_metadata_changes() {
        let path = temp_path("metadata.pdb");
        std::fs::write(&path, b"first").expect("fixture should be writable");
        let first = thumbnail_key(&path);

        std::fs::write(&path, b"second content").expect("fixture should update");
        let second = thumbnail_key(&path);

        assert_ne!(first, second);
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn missing_thumbnail_load_returns_none() {
        let cache = RecentThumbnailCache::new(temp_path("missing-cache"));
        assert!(cache.load(Path::new("/tmp/absent.pdb")).is_none());
    }

    #[test]
    fn recent_thumbnail_capture_is_square() {
        assert_eq!(RECENT_THUMBNAIL_WIDTH, RECENT_THUMBNAIL_HEIGHT);
    }
}
