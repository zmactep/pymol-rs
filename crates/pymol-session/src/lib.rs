pub mod convert;
pub mod pickle;
pub mod prs;
pub mod pse;

pub use convert::pse_to_session;
pub use prs::{load_prs, save_prs};

use std::io::Read;
use std::path::Path;

use flate2::read::GzDecoder;

use crate::pse::PseSession;

/// Errors that can occur when loading a PSE session.
#[derive(Debug, thiserror::Error)]
pub enum SessionError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("pickle error: {0}")]
    Pickle(String),

    #[error("PSE format error: {0}")]
    Pse(String),
}

/// Load a `.pse` file and return the parsed session data.
///
/// PSE files are typically gzip-compressed Python pickle (protocol 0–2),
/// but older PyMOL versions may save uncompressed pickle. This function
/// tries gzip first, then falls back to raw pickle.
pub fn load_pse(path: &Path) -> Result<PseSession, SessionError> {
    let raw = std::fs::read(path)?;

    // Try gzip decompression first
    let bytes = match decompress_gzip(&raw) {
        Ok(decompressed) => decompressed,
        Err(_) => raw, // Fall back to raw bytes (uncompressed pickle)
    };

    let value = pickle::read(&bytes)?;
    pse::reader::read_session(&value)
}

/// Attempt gzip decompression; returns Err if not valid gzip.
fn decompress_gzip(data: &[u8]) -> Result<Vec<u8>, std::io::Error> {
    let mut decoder = GzDecoder::new(data);
    let mut bytes = Vec::new();
    decoder.read_to_end(&mut bytes)?;
    Ok(bytes)
}

/// Load a PSE session from raw (already decompressed) pickle bytes.
pub fn load_pse_bytes(bytes: &[u8]) -> Result<PseSession, SessionError> {
    let value = pickle::read(bytes)?;
    pse::reader::read_session(&value)
}
