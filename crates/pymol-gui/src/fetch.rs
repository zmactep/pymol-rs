//! Fetch task for downloading structures from RCSB PDB
//!
//! This module provides async task support for fetching molecular structures
//! from the RCSB Protein Data Bank.

use std::future::Future;
use std::pin::Pin;
use std::time::Duration;

use pymol_io::FetchFormat;
use pymol_mol::ObjectMolecule;

use crate::async_tasks::{AsyncTask, TaskContext, TaskResult};

/// Default timeout for fetch operations (10 seconds)
const FETCH_TIMEOUT_SECS: u64 = 10;

// ============================================================================
// FetchResult
// ============================================================================

/// Result of a fetch operation
pub struct FetchResult {
    /// Object name to use when adding to registry
    pub name: String,
    /// PDB code that was fetched
    pub code: String,
    /// Result: molecule data or error message
    pub result: Result<ObjectMolecule, String>,
}

impl TaskResult for FetchResult {
    fn apply(self: Box<Self>, ctx: &mut dyn TaskContext) {
        match self.result {
            Ok(mol) => {
                ctx.add_molecule(&self.name, mol);
                ctx.print_info(format!(" Fetched {} as \"{}\"", self.code, self.name));
                // Use command system for zoom (also triggers redraw)
                ctx.execute_command(&format!("zoom {}", self.name));
            }
            Err(e) => {
                ctx.print_error(format!("Fetch failed for {}: {}", self.code, e));
            }
        }
    }
}

// ============================================================================
// FetchTask
// ============================================================================

/// Task to fetch a molecular structure from RCSB PDB
///
/// This task downloads a structure in the specified format (CIF or PDB)
/// and returns it as an `ObjectMolecule`.
pub struct FetchTask {
    /// PDB ID to fetch (e.g., "1ubq")
    pub code: String,
    /// Object name to use when adding to registry
    pub name: String,
    /// Format to fetch (CIF or PDB)
    pub format: FetchFormat,
}

impl FetchTask {
    /// Create a new fetch task
    pub fn new(code: String, name: String, format: FetchFormat) -> Self {
        Self { code, name, format }
    }
}

impl AsyncTask for FetchTask {
    fn notification_message(&self) -> String {
        format!("Fetching {}...", self.code)
    }

    fn execute(self) -> Pin<Box<dyn Future<Output = Box<dyn TaskResult>> + Send>> {
        Box::pin(async move {
            let timeout_duration = Duration::from_secs(FETCH_TIMEOUT_SECS);

            let result = match tokio::time::timeout(
                timeout_duration,
                pymol_io::fetch_async(&self.code, self.format),
            )
            .await
            {
                Ok(Ok(mol)) => Ok(mol),
                Ok(Err(e)) => Err(format!("Network error: {}", e)),
                Err(_) => Err(format!(
                    "Timeout: fetch did not complete within {} seconds",
                    FETCH_TIMEOUT_SECS
                )),
            };

            Box::new(FetchResult {
                name: self.name,
                code: self.code,
                result,
            }) as Box<dyn TaskResult>
        })
    }
}