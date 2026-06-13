//! Async fetch task integration for Patinae.

use std::future::Future;
use std::pin::Pin;
use std::time::Duration;

use patinae_cmd::{AsyncCommandRequest, FetchFormatCode, FetchRequest};
use patinae_framework::kernel::AppKernel;
use patinae_framework::tasks::{AsyncTask, TaskFailureHandler, TaskResult};
use patinae_framework::topics;
use patinae_mol::ObjectMolecule;

const FETCH_TIMEOUT_SECS: u64 = 10;
const METADATA_TIMEOUT_SECS: u64 = 6;

pub const PDB_METADATA_TOPIC: &str = "patinae.ui.fetch_pdb.metadata";

#[derive(Clone, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct PdbMetadataMessage {
    pub pdb_id: String,
    pub title: String,
    pub details: String,
    pub method: String,
    pub deposit_date: String,
    pub release_date: String,
    pub doi: String,
    pub error: String,
}

/// Convert a command async request into a Patinae task.
pub fn task_from_request(request: AsyncCommandRequest) -> Option<Box<dyn AsyncTask>> {
    match request {
        AsyncCommandRequest::Fetch(request) => Some(Box::new(FetchTask { request })),
    }
}

struct FetchTask {
    request: FetchRequest,
}

pub struct PdbMetadataTask {
    pdb_id: String,
}

impl PdbMetadataTask {
    pub fn new(pdb_id: impl Into<String>) -> Self {
        Self {
            pdb_id: pdb_id.into(),
        }
    }
}

impl AsyncTask for FetchTask {
    fn notification_message(&self) -> String {
        format!("Fetching {}...", self.request.code)
    }

    fn failure_handler(&self) -> TaskFailureHandler {
        let request = self.request.clone();
        Box::new(move |error| {
            Box::new(FetchResult {
                request,
                result: Err(error),
            })
        })
    }

    fn execute(self: Box<Self>) -> Pin<Box<dyn Future<Output = Box<dyn TaskResult>> + Send>> {
        Box::pin(async move {
            let request = self.request;
            let format = to_io_format(request.format);
            let result = match tokio::time::timeout(
                Duration::from_secs(FETCH_TIMEOUT_SECS),
                patinae_io::fetch_async_with_bond_tolerance(
                    &request.code,
                    format,
                    request.bond_tolerance,
                ),
            )
            .await
            {
                Ok(Ok(mol)) => Ok(mol),
                Ok(Err(err)) => Err(err.to_string()),
                Err(_) => Err(format!("timeout after {} seconds", FETCH_TIMEOUT_SECS)),
            };

            Box::new(FetchResult { request, result }) as Box<dyn TaskResult>
        })
    }
}

struct FetchResult {
    request: FetchRequest,
    result: Result<ObjectMolecule, String>,
}

impl TaskResult for FetchResult {
    fn apply(self: Box<Self>, kernel: &mut AppKernel) {
        match self.result {
            Ok(mol) => kernel.apply_fetched_molecule(&self.request, mol),
            Err(err) => kernel.print_fetch_error(&self.request, err),
        }
    }
}

impl AsyncTask for PdbMetadataTask {
    fn notification_message(&self) -> String {
        format!("Looking up {} metadata...", self.pdb_id)
    }

    fn failure_handler(&self) -> TaskFailureHandler {
        let pdb_id = self.pdb_id.clone();
        Box::new(move |error| {
            Box::new(PdbMetadataResult {
                pdb_id,
                result: Err(error),
            })
        })
    }

    fn execute(self: Box<Self>) -> Pin<Box<dyn Future<Output = Box<dyn TaskResult>> + Send>> {
        Box::pin(async move {
            let pdb_id = self.pdb_id;
            let result = match tokio::time::timeout(
                Duration::from_secs(METADATA_TIMEOUT_SECS),
                patinae_io::fetch_pdb_metadata(&pdb_id),
            )
            .await
            {
                Ok(Ok(preview)) => Ok(preview),
                Ok(Err(err)) => Err(err.to_string()),
                Err(_) => Err(format!("timeout after {} seconds", METADATA_TIMEOUT_SECS)),
            };

            Box::new(PdbMetadataResult { pdb_id, result }) as Box<dyn TaskResult>
        })
    }
}

struct PdbMetadataResult {
    pdb_id: String,
    result: Result<patinae_io::PdbMetadataPreview, String>,
}

impl TaskResult for PdbMetadataResult {
    fn apply(self: Box<Self>, kernel: &mut AppKernel) {
        let message = match self.result {
            Ok(preview) => PdbMetadataMessage {
                pdb_id: self.pdb_id,
                title: format!("{} - {}", preview.rcsb_id, preview.title),
                details: String::new(),
                method: preview.methods.join(", "),
                deposit_date: preview.deposit_date.unwrap_or_default(),
                release_date: preview.release_date.unwrap_or_default(),
                doi: preview.doi.unwrap_or_default(),
                error: String::new(),
            },
            Err(error) => PdbMetadataMessage {
                pdb_id: self.pdb_id,
                title: String::new(),
                details: String::new(),
                method: String::new(),
                deposit_date: String::new(),
                release_date: String::new(),
                doi: String::new(),
                error,
            },
        };
        topics::publish(&mut kernel.bus, PDB_METADATA_TOPIC, &message);
    }
}

fn to_io_format(format: FetchFormatCode) -> patinae_io::FetchFormat {
    match format {
        FetchFormatCode::Pdb => patinae_io::FetchFormat::Pdb,
        FetchFormatCode::Cif => patinae_io::FetchFormat::Cif,
        FetchFormatCode::Bcif => patinae_io::FetchFormat::Bcif,
    }
}
