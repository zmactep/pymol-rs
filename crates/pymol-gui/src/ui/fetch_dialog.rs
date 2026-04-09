//! Fetch PDB Dialog
//!
//! Modal dialog for fetching structures from RCSB PDB with metadata preview.

use std::sync::mpsc;

use egui::{Align2, Color32, Id, RichText};

// ============================================================================
// GraphQL query & response types
// ============================================================================

const RCSB_GRAPHQL_URL: &str = "https://data.rcsb.org/graphql";

const RCSB_QUERY: &str = r#"
query structure($id: String!) {
  entry(entry_id: $id) {
    rcsb_id
    struct {
      title
    }
    pubmed {
      rcsb_pubmed_doi
      rcsb_pubmed_abstract_text
    }
    rcsb_accession_info {
      deposit_date
      initial_release_date
    }
    rcsb_entry_info {
      structure_determination_methodology
    }
  }
}
"#;

/// Parsed metadata from RCSB GraphQL response.
pub struct PdbPreview {
    pub rcsb_id: String,
    pub title: String,
    pub deposit_date: Option<String>,
    pub release_date: Option<String>,
    pub methodology: Option<String>,
    pub pubmed_doi: Option<String>,
    pub pubmed_abstract: Option<String>,
}

/// Parse the RCSB GraphQL JSON response into a `PdbPreview`.
fn parse_preview(json: &serde_json::Value) -> Result<PdbPreview, String> {
    let entry = json
        .pointer("/data/entry")
        .ok_or("Structure not found in RCSB")?;

    if entry.is_null() {
        return Err("PDB ID not found".into());
    }

    let rcsb_id = entry
        .pointer("/rcsb_id")
        .and_then(|v| v.as_str())
        .unwrap_or("")
        .to_string();

    let title = entry
        .pointer("/struct/title")
        .and_then(|v| v.as_str())
        .unwrap_or("(no title)")
        .to_string();

    let deposit_date = entry
        .pointer("/rcsb_accession_info/deposit_date")
        .and_then(|v| v.as_str())
        .map(|s| s.to_string());

    let release_date = entry
        .pointer("/rcsb_accession_info/initial_release_date")
        .and_then(|v| v.as_str())
        .map(|s| s.to_string());

    let methodology = entry
        .pointer("/rcsb_entry_info/structure_determination_methodology")
        .and_then(|v| v.as_str())
        .map(|s| s.to_string());

    let pubmed_doi = entry
        .pointer("/pubmed/rcsb_pubmed_doi")
        .and_then(|v| v.as_str())
        .map(|s| s.to_string());

    let pubmed_abstract = entry
        .pointer("/pubmed/rcsb_pubmed_abstract_text")
        .and_then(|v| v.as_str())
        .map(|s| s.to_string());

    Ok(PdbPreview {
        rcsb_id,
        title,
        deposit_date,
        release_date,
        methodology,
        pubmed_doi,
        pubmed_abstract,
    })
}

/// Fetch preview metadata from RCSB (async).
async fn fetch_preview(pdb_id: &str) -> Result<PdbPreview, String> {
    let body = serde_json::json!({
        "query": RCSB_QUERY,
        "variables": { "id": pdb_id }
    });

    let resp = reqwest::Client::new()
        .post(RCSB_GRAPHQL_URL)
        .json(&body)
        .send()
        .await
        .map_err(|e| format!("Network error: {}", e))?;

    let json: serde_json::Value = resp
        .json()
        .await
        .map_err(|e| format!("Invalid response: {}", e))?;

    if let Some(errors) = json.get("errors") {
        if let Some(msg) = errors.pointer("/0/message").and_then(|v| v.as_str()) {
            return Err(msg.to_string());
        }
    }

    parse_preview(&json)
}

// ============================================================================
// FetchDialog
// ============================================================================

/// Modal dialog for fetching PDB structures with RCSB metadata preview.
pub struct FetchDialog {
    visible: bool,
    pdb_id: String,
    preview: Option<PdbPreview>,
    loading: bool,
    error: Option<String>,
    /// Command to execute when Fetch is clicked (consumed by caller).
    pending_command: Option<String>,
    preview_tx: mpsc::Sender<Result<PdbPreview, String>>,
    preview_rx: mpsc::Receiver<Result<PdbPreview, String>>,
    /// Whether to focus the text input on next frame.
    focus_input: bool,
}

impl Default for FetchDialog {
    fn default() -> Self {
        Self::new()
    }
}

impl FetchDialog {
    pub fn new() -> Self {
        let (tx, rx) = mpsc::channel();
        Self {
            visible: false,
            pdb_id: String::new(),
            preview: None,
            loading: false,
            error: None,
            pending_command: None,
            preview_tx: tx,
            preview_rx: rx,
            focus_input: false,
        }
    }

    /// Open the dialog.
    pub fn open(&mut self) {
        self.visible = true;
        self.pdb_id.clear();
        self.preview = None;
        self.loading = false;
        self.error = None;
        self.focus_input = true;
    }

    /// Take a pending fetch command, if any (called by App after show).
    pub fn take_command(&mut self) -> Option<String> {
        self.pending_command.take()
    }

    /// Render the dialog. Returns true if visible.
    pub fn show(&mut self, ctx: &egui::Context, rt_handle: &tokio::runtime::Handle) -> bool {
        if !self.visible {
            return false;
        }

        // Poll for preview results
        if let Ok(result) = self.preview_rx.try_recv() {
            self.loading = false;
            match result {
                Ok(preview) => {
                    self.error = None;
                    self.preview = Some(preview);
                }
                Err(e) => {
                    self.preview = None;
                    self.error = Some(e);
                }
            }
        }

        let mut should_fetch = false;
        let mut should_close = false;

        egui::Window::new("Fetch PDB")
            .anchor(Align2::CENTER_CENTER, [0.0, 0.0])
            .resizable(false)
            .collapsible(false)
            .default_width(420.0)
            .show(ctx, |ui| {
                ui.spacing_mut().item_spacing.y = 8.0;

                // PDB ID input
                ui.horizontal(|ui| {
                    ui.label("PDB ID:");
                    let input = ui.add(
                        egui::TextEdit::singleline(&mut self.pdb_id)
                            .desired_width(80.0)
                            .hint_text("e.g. 1UBQ")
                            .char_limit(10),
                    );

                    if self.focus_input {
                        input.request_focus();
                        self.focus_input = false;
                    }

                    // Auto-uppercase
                    self.pdb_id = self.pdb_id.to_uppercase();

                    // Enter triggers preview
                    let enter_pressed = input.lost_focus()
                        && ui.input(|i| i.key_pressed(egui::Key::Enter));

                    let preview_enabled = self.pdb_id.len() >= 4 && !self.loading;

                    if ui.add_enabled(preview_enabled, egui::Button::new("Preview")).clicked()
                        || (enter_pressed && preview_enabled)
                    {
                        self.fire_preview(rt_handle);
                    }

                    if self.loading {
                        ui.spinner();
                    }
                });

                ui.separator();

                // Preview content
                if let Some(ref preview) = self.preview {
                    self.render_preview(ui, preview);
                } else if let Some(ref err) = self.error {
                    ui.colored_label(Color32::from_rgb(255, 100, 100), err);
                } else {
                    ui.colored_label(
                        Color32::from_gray(120),
                        "Enter a PDB ID and press Enter to preview.",
                    );
                }

                ui.separator();

                // Buttons
                ui.horizontal(|ui| {
                    let fetch_enabled = self.preview.is_some() && !self.loading;
                    if ui.add_enabled(fetch_enabled, egui::Button::new("Fetch")).clicked() {
                        should_fetch = true;
                    }
                    if ui.button("Cancel").clicked() {
                        should_close = true;
                    }
                });

                // Esc to close
                if ui.input(|i| i.key_pressed(egui::Key::Escape)) {
                    should_close = true;
                }
            });

        if should_fetch {
            self.pending_command = Some(format!("fetch {}", self.pdb_id));
            self.visible = false;
        }

        if should_close {
            self.visible = false;
        }

        self.visible
    }

    /// Fire an async RCSB preview query.
    fn fire_preview(&mut self, rt_handle: &tokio::runtime::Handle) {
        self.loading = true;
        self.preview = None;
        self.error = None;

        let pdb_id = self.pdb_id.clone();
        let tx = self.preview_tx.clone();

        rt_handle.spawn(async move {
            let result = fetch_preview(&pdb_id).await;
            let _ = tx.send(result);
        });
    }

    /// Render the preview card.
    fn render_preview(&self, ui: &mut egui::Ui, preview: &PdbPreview) {
        ui.label(RichText::new(&preview.rcsb_id).strong().size(16.0));
        ui.label(
            RichText::new(&preview.title)
                .italics()
                .color(Color32::from_gray(200)),
        );

        ui.add_space(4.0);

        egui::Grid::new("pdb_preview_grid")
            .num_columns(2)
            .spacing([12.0, 4.0])
            .show(ui, |ui| {
                if let Some(ref meth) = preview.methodology {
                    ui.label("Method:");
                    ui.label(meth);
                    ui.end_row();
                }
                if let Some(ref date) = preview.deposit_date {
                    ui.label("Deposited:");
                    ui.label(format_date(date));
                    ui.end_row();
                }
                if let Some(ref date) = preview.release_date {
                    ui.label("Released:");
                    ui.label(format_date(date));
                    ui.end_row();
                }
                if let Some(ref doi) = preview.pubmed_doi {
                    ui.label("DOI:");
                    ui.label(doi);
                    ui.end_row();
                }
            });

        if let Some(ref abs) = preview.pubmed_abstract {
            ui.add_space(4.0);
            egui::CollapsingHeader::new("Abstract")
                .id_salt(Id::new("fetch_abstract"))
                .show(ui, |ui| {
                    egui::ScrollArea::vertical()
                        .max_height(150.0)
                        .show(ui, |ui| {
                            ui.label(abs);
                        });
                });
        }
    }
}

/// Format an ISO date string (e.g. "2024-01-15T00:00:00Z") to "2024-01-15".
fn format_date(iso: &str) -> &str {
    iso.get(..10).unwrap_or(iso)
}
