//! Fetch molecular structures from RCSB PDB
//!
//! This module provides functionality to download molecular structures from the
//! RCSB Protein Data Bank by their PDB ID.
//!
//! # Features
//!
//! - `fetch` - Enables synchronous fetching using `ureq`
//! - `fetch-async` - Enables asynchronous fetching using `reqwest`
//!
//! # Example (synchronous)
//!
//! ```no_run
//! # #[cfg(feature = "fetch")]
//! # fn main() -> pymol_io::IoResult<()> {
//! use pymol_io::fetch::{fetch, FetchFormat};
//!
//! // Fetch structure in BinaryCIF format (default)
//! let mol = fetch("1ubq", FetchFormat::default())?;
//!
//! // Fetch in PDB format
//! let mol = fetch("4hhb", FetchFormat::Pdb)?;
//! # Ok(())
//! # }
//! # #[cfg(not(feature = "fetch"))]
//! # fn main() {}
//! ```
//!
//! # Example (asynchronous)
//!
//! ```no_run
//! # #[cfg(feature = "fetch-async")]
//! # async fn example() -> pymol_io::IoResult<()> {
//! use pymol_io::fetch::{fetch_async, FetchFormat};
//!
//! // Async fetch
//! let mol = fetch_async("1ubq", FetchFormat::Cif).await?;
//! # Ok(())
//! # }
//! ```

use std::io::Cursor;

use crate::compress::gzip_reader;
use crate::error::{IoError, IoResult};
use pymol_mol::ObjectMolecule;

/// RCSB PDB base URL for text file downloads (PDB, mmCIF)
const RCSB_BASE_URL: &str = "https://files.rcsb.org/download";

/// RCSB models API base URL for BinaryCIF downloads
const RCSB_MODELS_URL: &str = "https://models.rcsb.org";

/// User-Agent header for HTTP requests
const USER_AGENT: &str = concat!("pymol-rs/", env!("CARGO_PKG_VERSION"));

/// Format to fetch from RCSB PDB
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum FetchFormat {
    /// PDB format (.pdb)
    ///
    /// Traditional PDB format. Note that PDB format has limitations for
    /// large structures (>99,999 atoms or >9,999 residues).
    Pdb,

    /// mmCIF format (.cif)
    ///
    /// Macromolecular Crystallographic Information File format.
    /// Handles large structures and contains more complete information
    /// than PDB format.
    Cif,

    /// BinaryCIF format (.bcif) - default
    ///
    /// Binary encoding of mmCIF using MessagePack. Smallest file size
    /// and fastest parsing. The primary format served by RCSB PDB.
    #[default]
    Bcif,
}

impl FetchFormat {
    /// Get the file extension for this format
    pub fn extension(&self) -> &'static str {
        match self {
            FetchFormat::Pdb => "pdb",
            FetchFormat::Cif => "cif",
            FetchFormat::Bcif => "bcif",
        }
    }

    /// Get a human-readable name for this format
    pub fn name(&self) -> &'static str {
        match self {
            FetchFormat::Pdb => "PDB",
            FetchFormat::Cif => "mmCIF",
            FetchFormat::Bcif => "BinaryCIF",
        }
    }
}

/// Validate a PDB ID
///
/// PDB IDs are 4-character alphanumeric codes (e.g., "1ubq", "4hhb").
/// This function validates the format but does not check if the ID exists.
///
/// # Arguments
///
/// * `pdb_id` - The PDB ID to validate
///
/// # Returns
///
/// Returns `Ok(())` if valid, or an error if the PDB ID is invalid.
pub fn validate_pdb_id(pdb_id: &str) -> IoResult<()> {
    if pdb_id.len() != 4 {
        return Err(IoError::invalid_pdb_id(format!(
            "'{}' - must be exactly 4 characters",
            pdb_id
        )));
    }

    if !pdb_id.chars().all(|c| c.is_ascii_alphanumeric()) {
        return Err(IoError::invalid_pdb_id(format!(
            "'{}' - must contain only alphanumeric characters",
            pdb_id
        )));
    }

    Ok(())
}

/// Build the RCSB download URL for a PDB ID
///
/// # Arguments
///
/// * `pdb_id` - The 4-character PDB ID (case-insensitive)
/// * `format` - The format to download
///
/// # Returns
///
/// The full URL to download the structure file (gzip-compressed).
pub fn build_rcsb_url(pdb_id: &str, format: FetchFormat) -> String {
    // RCSB uses lowercase PDB IDs in URLs
    let pdb_id_lower = pdb_id.to_lowercase();
    let url = match format {
        FetchFormat::Bcif => RCSB_MODELS_URL,
        _ => RCSB_BASE_URL
    };
    format!(
        "{}/{}.{}.gz",
        url,
        pdb_id_lower,
        format.extension()
    )
}

/// Parse the fetched content based on format
fn parse_content(content: &[u8], format: FetchFormat) -> IoResult<ObjectMolecule> {
    let decoder = gzip_reader(Cursor::new(content));
    match format {
        FetchFormat::Pdb => crate::pdb::read_pdb_from(decoder),
        FetchFormat::Cif => crate::cif::read_cif_from(decoder),
        FetchFormat::Bcif => crate::bcif::read_bcif_from(decoder)
    }
}

/// Fetch a molecular structure from RCSB PDB (synchronous)
///
/// Downloads and parses a molecular structure from the RCSB Protein Data Bank.
///
/// # Arguments
///
/// * `pdb_id` - The 4-character PDB ID (e.g., "1ubq", "4hhb"). Case-insensitive.
/// * `format` - The format to fetch ([`FetchFormat::Pdb`] or [`FetchFormat::Cif`])
///
/// # Returns
///
/// The parsed molecular structure, or an error if the fetch or parse failed.
///
/// # Example
///
/// ```no_run
/// # #[cfg(feature = "fetch")]
/// # fn main() -> pymol_io::IoResult<()> {
/// use pymol_io::fetch::{fetch, FetchFormat};
///
/// // Fetch ubiquitin in mmCIF format
/// let mol = fetch("1ubq", FetchFormat::Cif)?;
/// println!("Fetched {} atoms", mol.atom_count());
/// # Ok(())
/// # }
/// # #[cfg(not(feature = "fetch"))]
/// # fn main() {}
/// ```
///
/// # Errors
///
/// Returns an error if:
/// - The PDB ID is invalid (not 4 alphanumeric characters)
/// - The network request fails
/// - The PDB ID does not exist (404 error)
/// - The downloaded content cannot be parsed
#[cfg(feature = "fetch")]
pub fn fetch(pdb_id: &str, format: FetchFormat) -> IoResult<ObjectMolecule> {
    use std::io::Read;

    validate_pdb_id(pdb_id)?;

    let url = build_rcsb_url(pdb_id, format);

    let response = ureq::get(&url)
        .header("User-Agent", USER_AGENT)
        .call()
        .map_err(|e| IoError::fetch(format!("Network error: {}", e)))?;

    let status = response.status().as_u16();
    if status == 404 {
        return Err(IoError::fetch(format!("PDB ID '{}' not found", pdb_id)));
    }
    if status >= 400 {
        return Err(IoError::fetch(format!("HTTP error {}: {}", status, url)));
    }

    let mut content = Vec::new();
    response
        .into_body()
        .into_reader()
        .read_to_end(&mut content)
        .map_err(|e| IoError::fetch(format!("Failed to read response: {}", e)))?;

    parse_content(&content, format)
}

/// Fetch a molecular structure from RCSB PDB (asynchronous)
///
/// Downloads and parses a molecular structure from the RCSB Protein Data Bank.
///
/// # Arguments
///
/// * `pdb_id` - The 4-character PDB ID (e.g., "1ubq", "4hhb"). Case-insensitive.
/// * `format` - The format to fetch ([`FetchFormat::Pdb`] or [`FetchFormat::Cif`])
///
/// # Returns
///
/// The parsed molecular structure, or an error if the fetch or parse failed.
///
/// # Example
///
/// ```no_run
/// # #[cfg(feature = "fetch-async")]
/// # async fn example() -> pymol_io::IoResult<()> {
/// use pymol_io::fetch::{fetch_async, FetchFormat};
///
/// // Fetch hemoglobin in PDB format
/// let mol = fetch_async("4hhb", FetchFormat::Pdb).await?;
/// println!("Fetched {} atoms", mol.atom_count());
/// # Ok(())
/// # }
/// ```
///
/// # Errors
///
/// Returns an error if:
/// - The PDB ID is invalid (not 4 alphanumeric characters)
/// - The network request fails
/// - The PDB ID does not exist (404 error)
/// - The downloaded content cannot be parsed
#[cfg(feature = "fetch-async")]
pub async fn fetch_async(pdb_id: &str, format: FetchFormat) -> IoResult<ObjectMolecule> {
    validate_pdb_id(pdb_id)?;

    let url = build_rcsb_url(pdb_id, format);

    let client = reqwest::Client::new();
    let response = client
        .get(&url)
        .header("User-Agent", USER_AGENT)
        .send()
        .await
        .map_err(|e| IoError::fetch(format!("Network error: {}", e)))?;

    let status = response.status();
    if status == reqwest::StatusCode::NOT_FOUND {
        return Err(IoError::fetch(format!("PDB ID '{}' not found", pdb_id)));
    }
    if !status.is_success() {
        return Err(IoError::fetch(format!(
            "HTTP error {}: {}",
            status.as_u16(),
            url
        )));
    }

    let content = response
        .bytes()
        .await
        .map_err(|e| IoError::fetch(format!("Failed to read response: {}", e)))?;

    parse_content(&content, format)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_validate_pdb_id_valid() {
        assert!(validate_pdb_id("1ubq").is_ok());
        assert!(validate_pdb_id("4HHB").is_ok());
        assert!(validate_pdb_id("1A2B").is_ok());
        assert!(validate_pdb_id("9xyz").is_ok());
    }

    #[test]
    fn test_validate_pdb_id_invalid_length() {
        assert!(validate_pdb_id("").is_err());
        assert!(validate_pdb_id("1ub").is_err());
        assert!(validate_pdb_id("1ubqx").is_err());
        assert!(validate_pdb_id("12345").is_err());
    }

    #[test]
    fn test_validate_pdb_id_invalid_chars() {
        assert!(validate_pdb_id("1ub!").is_err());
        assert!(validate_pdb_id("1u-q").is_err());
        assert!(validate_pdb_id("1u q").is_err());
        assert!(validate_pdb_id("1u.q").is_err());
    }

    #[test]
    fn test_build_rcsb_url() {
        assert_eq!(
            build_rcsb_url("1UBQ", FetchFormat::Pdb),
            "https://files.rcsb.org/download/1ubq.pdb.gz"
        );
        assert_eq!(
            build_rcsb_url("4hhb", FetchFormat::Cif),
            "https://files.rcsb.org/download/4hhb.cif.gz"
        );
        assert_eq!(
            build_rcsb_url("1UBQ", FetchFormat::Bcif),
            "https://models.rcsb.org/1ubq.bcif.gz"
        );
    }

    #[test]
    fn test_fetch_format_extension() {
        assert_eq!(FetchFormat::Pdb.extension(), "pdb");
        assert_eq!(FetchFormat::Cif.extension(), "cif");
        assert_eq!(FetchFormat::Bcif.extension(), "bcif");
    }

    #[test]
    fn test_fetch_format_name() {
        assert_eq!(FetchFormat::Pdb.name(), "PDB");
        assert_eq!(FetchFormat::Cif.name(), "mmCIF");
        assert_eq!(FetchFormat::Bcif.name(), "BinaryCIF");
    }

    #[test]
    fn test_fetch_format_default() {
        assert_eq!(FetchFormat::default(), FetchFormat::Bcif);
    }

    // Integration tests - require network access
    // Run with: cargo test -p pymol-io --features fetch -- --ignored

    #[cfg(feature = "fetch")]
    #[test]
    #[ignore = "requires network access"]
    fn test_fetch_pdb_format() {
        let mol = fetch("1ubq", FetchFormat::Pdb).expect("Failed to fetch 1ubq in PDB format");
        // Ubiquitin has 76 residues, ~600 atoms
        assert!(mol.atom_count() > 500, "Expected >500 atoms, got {}", mol.atom_count());
    }

    #[cfg(feature = "fetch")]
    #[test]
    #[ignore = "requires network access"]
    fn test_fetch_cif_format() {
        let mol = fetch("1ubq", FetchFormat::Cif).expect("Failed to fetch 1ubq in mmCIF format");
        // Ubiquitin has 76 residues, ~600 atoms
        assert!(mol.atom_count() > 500, "Expected >500 atoms, got {}", mol.atom_count());
    }

    #[cfg(feature = "fetch")]
    #[test]
    #[ignore = "requires network access"]
    fn test_fetch_bcif_format() {
        let mol = fetch("1ubq", FetchFormat::Bcif).expect("Failed to fetch 1ubq in bCIF format");
        // Ubiquitin has 76 residues, ~600 atoms
        assert!(mol.atom_count() > 500, "Expected >500 atoms, got {}", mol.atom_count());
    }

    #[cfg(feature = "fetch")]
    #[test]
    #[ignore = "requires network access"]
    fn test_fetch_nonexistent_pdb() {
        let result = fetch("9999", FetchFormat::Pdb);
        assert!(result.is_err(), "Expected error for nonexistent PDB ID");
        let err = result.unwrap_err();
        assert!(
            matches!(err, IoError::Fetch(_)),
            "Expected Fetch error, got {:?}",
            err
        );
    }

    #[cfg(feature = "fetch")]
    #[test]
    fn test_fetch_invalid_pdb_id() {
        let result = fetch("invalid", FetchFormat::Pdb);
        assert!(result.is_err(), "Expected error for invalid PDB ID");
        let err = result.unwrap_err();
        assert!(
            matches!(err, IoError::InvalidPdbId(_)),
            "Expected InvalidPdbId error, got {:?}",
            err
        );
    }

    // Async integration tests
    #[cfg(feature = "fetch-async")]
    #[tokio::test]
    #[ignore = "requires network access"]
    async fn test_fetch_async_pdb_format() {
        let mol = fetch_async("1ubq", FetchFormat::Pdb)
            .await
            .expect("Failed to fetch 1ubq in PDB format");
        assert!(mol.atom_count() > 500, "Expected >500 atoms, got {}", mol.atom_count());
    }

    #[cfg(feature = "fetch-async")]
    #[tokio::test]
    #[ignore = "requires network access"]
    async fn test_fetch_async_cif_format() {
        let mol = fetch_async("1ubq", FetchFormat::Cif)
            .await
            .expect("Failed to fetch 1ubq in mmCIF format");
        assert!(mol.atom_count() > 500, "Expected >500 atoms, got {}", mol.atom_count());
    }

    #[cfg(feature = "fetch-async")]
    #[tokio::test]
    #[ignore = "requires network access"]
    async fn test_fetch_async_bcif_format() {
        let mol = fetch_async("1ubq", FetchFormat::Bcif)
            .await
            .expect("Failed to fetch 1ubq in bCIF format");
        assert!(mol.atom_count() > 500, "Expected >500 atoms, got {}", mol.atom_count());
    }

    #[cfg(feature = "fetch-async")]
    #[tokio::test]
    #[ignore = "requires network access"]
    async fn test_fetch_async_nonexistent_pdb() {
        let result = fetch_async("9999", FetchFormat::Pdb).await;
        assert!(result.is_err(), "Expected error for nonexistent PDB ID");
        let err = result.unwrap_err();
        assert!(
            matches!(err, IoError::Fetch(_)),
            "Expected Fetch error, got {:?}",
            err
        );
    }

    #[cfg(feature = "fetch-async")]
    #[tokio::test]
    async fn test_fetch_async_invalid_pdb_id() {
        let result = fetch_async("invalid", FetchFormat::Pdb).await;
        assert!(result.is_err(), "Expected error for invalid PDB ID");
        let err = result.unwrap_err();
        assert!(
            matches!(err, IoError::InvalidPdbId(_)),
            "Expected InvalidPdbId error, got {:?}",
            err
        );
    }
}
