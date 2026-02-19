//! Native PRS (PyMOL-RS Session) format: MessagePack + gzip.

use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;

use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use pymol_scene::Session;

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
    let data = rmp_serde::to_vec(session)?;
    let file = File::create(path)?;
    let mut encoder = GzEncoder::new(file, Compression::default());
    encoder.write_all(&data)?;
    encoder.finish()?;
    Ok(())
}

/// Load a session from a `.prs` file (MessagePack + gzip).
pub fn load_prs(path: &Path) -> Result<Session, PrsError> {
    let file = File::open(path)?;
    let mut decoder = GzDecoder::new(file);
    let mut bytes = Vec::new();
    decoder.read_to_end(&mut bytes)?;
    let session: Session = rmp_serde::from_slice(&bytes)?;
    Ok(session)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn test_prs_round_trip() {
        let session = Session::new();
        let dir = std::env::temp_dir().join("pymol_prs_test");
        fs::create_dir_all(&dir).unwrap();
        let path = dir.join("test.prs");

        save_prs(&session, &path).unwrap();
        assert!(path.exists());

        let loaded = load_prs(&path).unwrap();
        assert_eq!(loaded.clear_color, session.clear_color);

        // Cleanup
        let _ = fs::remove_dir_all(&dir);
    }
}
