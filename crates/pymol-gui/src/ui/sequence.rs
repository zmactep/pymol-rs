//! Sequence Viewer Panel
//!
//! Displays amino acid / nucleotide sequences of loaded structures with
//! interactive selection synchronized with the 3D view.

use std::collections::HashSet;
use std::ops::Range;

use egui::{Color32, RichText, Ui};
use pymol_color::NamedColors;
use pymol_mol::{
    is_capping_group, is_ion, is_standard_amino_acid, is_standard_nucleotide, is_water,
    nucleotide_to_char, residue_to_char, three_to_one, ObjectMolecule,
};
use pymol_scene::{Object, ObjectRegistry};
use pymol_select::build_sele_command;

/// Compress a sorted list of residue numbers into range notation.
/// E.g., `[74, 75, 76, 77, 80, 85, 86, 87]` → `"74-77+80+85-87"`
fn compress_resi_list(values: &[i32]) -> String {
    if values.is_empty() {
        return String::new();
    }
    let mut parts = Vec::new();
    let mut start = values[0];
    let mut end = values[0];
    for &v in &values[1..] {
        if v == end + 1 {
            end = v;
        } else {
            if start == end {
                parts.push(start.to_string());
            } else {
                parts.push(format!("{}-{}", start, end));
            }
            start = v;
            end = v;
        }
    }
    if start == end {
        parts.push(start.to_string());
    } else {
        parts.push(format!("{}-{}", start, end));
    }
    parts.join("+")
}

/// Pale magenta — standard DNA/RNA nucleotides
const NUCLEOTIDE_COLOR: Color32 = Color32::from_rgb(205, 150, 205);

/// Orange — non-canonical AAs, AA variants (HIP/CYX/MSE/…), modified bases, caps
const NONCANONICAL_COLOR: Color32 = Color32::from_rgb(230, 160, 60);

/// Dim gray — non-polymer HETATM ligands
const LIGAND_COLOR: Color32 = Color32::from_rgb(140, 140, 140);

/// Semantic classification of a residue for sequence viewer coloring.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ResidueKind {
    /// One of the 20 standard amino acids — PyMOL atom color
    AminoAcidCanonical,
    /// Known AA variant (HIP, CYX, MSE, SEP, ACE, NME, …) or unknown protein residue
    /// — single letter if parent known, brackets if truly unknown; pale yellow
    AminoAcidNonCanonical,
    /// Standard unmodified nucleotide — pale magenta
    NucleotideCanonical,
    /// Modified or unknown nucleotide — pale yellow
    NucleotideNonCanonical,
    /// Non-polymer HETATM: ligand, lipid, cofactor — gray
    Ligand,
}

impl ResidueKind {
    pub fn is_polymer(self) -> bool {
        !matches!(self, ResidueKind::Ligand)
    }
}

/// Single residue entry in the sequence viewer
#[derive(Debug, Clone)]
pub struct SeqResidue {
    pub chain: String,
    pub resn: String,
    pub resv: i32,
    pub display_char: char,
    /// Full display label: single char for polymer residues, "[ATP]" for ligands
    pub display_label: String,
    /// Number of character cells this residue occupies in the sequence row
    pub char_width: usize,
    pub atom_range: Range<usize>,
    /// Base color index of the CA atom (or first atom), for display coloring
    pub color_index: i32,
    /// Drives color and font size in the sequence viewer
    pub kind: ResidueKind,
}

/// Sequence data for one chain
#[derive(Debug, Clone)]
pub struct SeqChain {
    pub chain_id: String,
    pub residues: Vec<SeqResidue>,
}

/// Sequence data for one object
#[derive(Debug, Clone)]
pub struct SeqObject {
    pub object_name: String,
    pub chains: Vec<SeqChain>,
}

/// Reference to a specific residue in a loaded object
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ResidueRef {
    pub object_name: String,
    pub chain_id: String,
    pub resv: i32,
}

/// Action produced by the sequence panel
pub enum SequenceAction {
    /// Execute a command string (e.g., "select sele, ...")
    Execute(String),
    /// Display a notification message in the output panel
    Notify(String),
    /// Hover over a residue, or None when leaving
    Hover(Option<ResidueRef>),
}

/// Persistent state for the sequence viewer
pub struct SequencePanel {
    /// Cached sequence data per object
    sequences: Vec<SeqObject>,
    /// Number of objects when cache was last built
    cached_object_count: usize,
    /// Number of enabled objects when cache was last built
    cached_enabled_count: usize,
    /// Registry generation when cache was last built
    cached_generation: u64,
    /// Currently highlighted residues from 3D selection
    highlighted: HashSet<ResidueRef>,
    /// Active drag start: (object_index, chain_index, start_residue_index)
    drag_start: Option<(usize, usize, usize)>,
    /// Current drag end residue index (updated each frame while dragging)
    drag_end: Option<usize>,
    /// Current sequence hover state for change detection
    current_hover: Option<ResidueRef>,
}

impl Default for SequencePanel {
    fn default() -> Self {
        Self::new()
    }
}

impl SequencePanel {
    pub fn new() -> Self {
        Self {
            sequences: Vec::new(),
            cached_object_count: 0,
            cached_enabled_count: 0,
            cached_generation: 0,
            highlighted: HashSet::new(),
            drag_start: None,
            drag_end: None,
            current_hover: None,
        }
    }

    /// Clear any active drag state (called on dock/float transition).
    pub fn clear_drag(&mut self) {
        self.drag_start = None;
        self.drag_end = None;
    }

    /// Update the highlighted residues from the current 3D selection.
    pub fn update_highlights(&mut self, highlighted: HashSet<ResidueRef>) {
        self.highlighted = highlighted;
    }

    /// Rebuild the sequence cache from the registry
    fn rebuild_cache(&mut self, registry: &ObjectRegistry) {
        self.sequences.clear();

        for name in registry.names() {
            if let Some(mol_obj) = registry.get_molecule(name) {
                if !mol_obj.is_enabled() {
                    continue;
                }
                let mol = mol_obj.molecule();
                let seq_obj = build_seq_object(name, mol);
                if !seq_obj.chains.is_empty() {
                    self.sequences.push(seq_obj);
                }
            }
        }

        self.cached_object_count = registry.names().count();
        self.cached_enabled_count = registry.enabled_objects().count();
        self.cached_generation = registry.generation();
    }

    /// Check if cache needs rebuilding (object count, enabled count, or generation changed)
    fn needs_rebuild(&self, registry: &ObjectRegistry) -> bool {
        registry.names().count() != self.cached_object_count
            || registry.enabled_objects().count() != self.cached_enabled_count
            || registry.generation() != self.cached_generation
    }

    /// Compute the desired panel height based on the number of objects and chains.
    ///
    /// When the total chain count exceeds 4, the panel defaults to showing 4 chains
    /// and the rest are accessible via vertical scroll.
    pub fn desired_height(&mut self, registry: &ObjectRegistry) -> f32 {
        if self.needs_rebuild(registry) {
            self.rebuild_cache(registry);
        }
        if self.sequences.is_empty() {
            return 30.0;
        }
        // These must match the values used in show().
        // row_height: egui Body text style is typically ~18px (14pt font).
        // ruler_height: 10px for the residue-number ruler above each chain row.
        let row_height = 18.0;
        let ruler_height = 10.0;
        let spacing = 2.0;
        let chain_height = ruler_height + row_height;
        let total_chains: usize = self.sequences.iter().map(|s| s.chains.len()).sum();
        let visible_chains = total_chains.min(4) as f32;
        let obj_spacing = self.sequences.len() as f32 * spacing;
        // 16.0 accounts for panel padding + scroll area chrome
        visible_chains * chain_height + obj_spacing + 16.0
    }

    /// Show the sequence panel and return any actions
    pub fn show(
        &mut self,
        ui: &mut Ui,
        registry: &ObjectRegistry,
        named_colors: &NamedColors,
        has_sele: bool,
    ) -> Vec<SequenceAction> {
        if self.needs_rebuild(registry) {
            self.rebuild_cache(registry);
        }

        let mut actions = Vec::new();

        if self.sequences.is_empty() {
            ui.centered_and_justified(|ui| {
                ui.label(
                    RichText::new("No sequences loaded")
                        .color(Color32::GRAY)
                        .italics(),
                );
            });
            return actions;
        }

        let highlighted = &self.highlighted;
        let drag_start = &mut self.drag_start;
        let drag_end = &mut self.drag_end;

        let cw = char_width(ui);
        let ruler_height = 10.0;
        let row_height = ui.text_style_height(&egui::TextStyle::Body);
        let obj_name_font = egui::FontId::new(11.0, egui::FontFamily::Proportional);
        let obj_name_color = Color32::from_rgb(100, 255, 100);
        let label_font = egui::FontId::new(12.0, egui::FontFamily::Monospace);
        let label_color = Color32::from_rgb(180, 180, 180);
        let separator_color = Color32::from_rgb(60, 60, 60);
        let obj_name_width = 40.0;
        let label_width = cw * 2.5; // enough for "D:" etc.
        let spacing = 2.0;
        let n_objects = self.sequences.len();

        // Pre-compute per-object heights
        let obj_heights: Vec<f32> = self
            .sequences
            .iter()
            .map(|seq_obj| seq_obj.chains.len() as f32 * (ruler_height + row_height))
            .collect();

        // Track which residue the mouse hovers over (across all chain rows)
        let mut hover_out: Option<ResidueRef> = None;

        // Outer vertical scroll — keeps all three columns in sync when the
        // total chain height exceeds the panel (e.g. nucleosome with 12+ chains).
        egui::ScrollArea::vertical()
            .id_salt("sequence_vscroll")
            .auto_shrink([false, true])
            .show(ui, |ui| {
        // Three-column layout: vertical object names | chain labels | scrollable sequences
        ui.horizontal(|ui| {
            // Column 1: horizontal object names (fixed, truncated)
            ui.vertical(|ui| {
                ui.spacing_mut().item_spacing.y = 0.0;
                for (i, seq_obj) in self.sequences.iter().enumerate() {
                    let obj_height = obj_heights[i];
                    let (rect, response) = ui.allocate_exact_size(
                        egui::vec2(obj_name_width, obj_height),
                        egui::Sense::hover(),
                    );
                    let painter = ui.painter_at(rect);
                    let name = &seq_obj.object_name;
                    let display_name = if name.chars().count() > 4 {
                        let truncated: String = name.chars().take(3).collect();
                        format!("{truncated}...")
                    } else {
                        name.clone()
                    };
                    painter.text(
                        egui::pos2(rect.min.x, rect.center().y),
                        egui::Align2::LEFT_CENTER,
                        &display_name,
                        obj_name_font.clone(),
                        obj_name_color,
                    );
                    if name.chars().count() > 4 {
                        response.on_hover_text(name);
                    }
                    if i + 1 < n_objects {
                        let y = rect.max.y + spacing / 2.0;
                        painter.line_segment(
                            [egui::pos2(rect.min.x, y), egui::pos2(rect.max.x + 1000.0, y)],
                            egui::Stroke::new(1.0, separator_color),
                        );
                    }
                    ui.add_space(spacing);
                }
            });

            // Column 2: chain labels (fixed)
            ui.vertical(|ui| {
                ui.spacing_mut().item_spacing.y = 0.0;
                for seq_obj in self.sequences.iter() {
                    for chain in seq_obj.chains.iter() {
                        // Ruler spacer
                        ui.allocate_exact_size(
                            egui::vec2(label_width, ruler_height),
                            egui::Sense::hover(),
                        );
                        // Chain label aligned with sequence row
                        let (rect, _) = ui.allocate_exact_size(
                            egui::vec2(label_width, row_height),
                            egui::Sense::hover(),
                        );
                        let painter = ui.painter_at(rect);
                        painter.text(
                            rect.right_center(),
                            egui::Align2::RIGHT_CENTER,
                            format!("{}:", chain.chain_id),
                            label_font.clone(),
                            label_color,
                        );
                    }
                    ui.add_space(spacing);
                }
            });

            // Column 3: scrollable sequences (horizontal)
            egui::ScrollArea::horizontal()
                .id_salt("sequence_hscroll")
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.vertical(|ui| {
                        ui.spacing_mut().item_spacing.y = 0.0;
                        for (obj_idx, seq_obj) in self.sequences.iter().enumerate() {
                            for (chain_idx, chain) in seq_obj.chains.iter().enumerate() {
                                render_chain_sequence(
                                    ui,
                                    &seq_obj.object_name,
                                    chain,
                                    obj_idx,
                                    chain_idx,
                                    cw,
                                    named_colors,
                                    highlighted,
                                    &self.sequences,
                                    drag_start,
                                    drag_end,
                                    &mut actions,
                                    has_sele,
                                    &mut hover_out,
                                );
                            }
                            ui.add_space(spacing);
                        }
                    });
                });
        });
            });

        // Emit hover action only when the hovered residue changes
        if hover_out != self.current_hover {
            actions.push(SequenceAction::Hover(hover_out.clone()));
            self.current_hover = hover_out;
        }

        actions
    }
}

/// Measure the width of a single monospace character at the sequence font size
fn char_width(ui: &Ui) -> f32 {
    let font_id = egui::FontId::new(12.0, egui::FontFamily::Monospace);
    let galley = ui
        .painter()
        .layout_no_wrap("M".to_string(), font_id, Color32::WHITE);
    galley.size().x
}

/// Render the sequence of residues for one chain: ruler row + sequence row.
fn render_chain_sequence(
    ui: &mut Ui,
    object_name: &str,
    chain: &SeqChain,
    obj_idx: usize,
    chain_idx: usize,
    cw: f32,
    named_colors: &NamedColors,
    highlighted: &HashSet<ResidueRef>,
    all_sequences: &[SeqObject],
    drag_start: &mut Option<(usize, usize, usize)>,
    drag_end: &mut Option<usize>,
    actions: &mut Vec<SequenceAction>,
    has_sele: bool,
    hover_out: &mut Option<ResidueRef>,
) {
    let n = chain.residues.len();
    if n == 0 {
        return;
    }

    // Precompute cumulative x-offsets (in char cells) for each residue
    let mut x_offsets: Vec<usize> = Vec::with_capacity(n);
    let mut cumulative = 0usize;
    for residue in &chain.residues {
        x_offsets.push(cumulative);
        cumulative += residue.char_width;
    }
    let total_cells = cumulative;
    let total_width = cw * total_cells as f32;

    // --- Ruler row: residue numbers above the sequence ---
    let (ruler_rect, _) =
        ui.allocate_exact_size(egui::vec2(total_width, 10.0), egui::Sense::hover());
    let ruler_painter = ui.painter_at(ruler_rect);
    let ruler_font = egui::FontId::new(8.0, egui::FontFamily::Monospace);
    let ruler_color = Color32::from_rgb(120, 120, 120);

    for (res_idx, residue) in chain.residues.iter().enumerate() {
        if res_idx == 0 || res_idx % 10 == 0 {
            let x = ruler_rect.min.x + cw * x_offsets[res_idx] as f32;
            let label = format!("{}", residue.resv);
            // Skip if the label would overflow past the ruler width
            let label_width = label.len() as f32 * cw * 0.7; // ruler font is 8pt vs 12pt char width
            if x + label_width <= ruler_rect.max.x {
                ruler_painter.text(
                    egui::pos2(x, ruler_rect.min.y),
                    egui::Align2::LEFT_TOP,
                    label,
                    ruler_font.clone(),
                    ruler_color,
                );
            }
        }
    }

    // --- Sequence row: single interactive rect with painted characters ---
    let row_height = ui.text_style_height(&egui::TextStyle::Body);
    let (seq_rect, response) = ui.allocate_exact_size(
        egui::vec2(total_width, row_height),
        egui::Sense::click_and_drag(),
    );

    let painter = ui.painter_at(seq_rect);
    let seq_font = egui::FontId::new(12.0, egui::FontFamily::Monospace);
    let ligand_font = egui::FontId::new(10.0, egui::FontFamily::Monospace);

    // Ctrl/Cmd held = exclusion mode (remove from selection instead of adding)
    let ctrl_held = ui.input(|i| i.modifiers.ctrl || i.modifiers.command);

    // Helper: convert pointer x position to residue index using cumulative offsets
    let pointer_to_res_idx = |pos: egui::Pos2| -> usize {
        let rel_x = pos.x - seq_rect.min.x;
        let cell = (rel_x / cw).floor().max(0.0) as usize;
        // Binary search: find the last residue whose x_offset <= cell
        match x_offsets.binary_search(&cell) {
            Ok(idx) => idx.min(n - 1),
            Err(idx) => idx.saturating_sub(1).min(n - 1),
        }
    };

    // Determine the active drag range for this chain (if any)
    let drag_range = if let Some((d_obj, d_chain, d_start)) = *drag_start {
        if d_obj == obj_idx && d_chain == chain_idx {
            drag_end.map(|d_end| {
                let lo = d_start.min(d_end);
                let hi = d_start.max(d_end);
                lo..=hi
            })
        } else {
            None
        }
    } else {
        None
    };

    // Draw drag highlight background
    if let Some(ref range) = drag_range {
        let x_start = seq_rect.min.x + cw * x_offsets[*range.start()] as f32;
        let end_res = *range.end();
        let x_end = seq_rect.min.x
            + cw * (x_offsets[end_res] + chain.residues[end_res].char_width) as f32;
        let highlight_rect = egui::Rect::from_min_max(
            egui::pos2(x_start, seq_rect.min.y),
            egui::pos2(x_end, seq_rect.max.y),
        );
        let drag_color = if ctrl_held {
            Color32::from_rgba_premultiplied(255, 100, 50, 80)
        } else {
            Color32::from_rgba_premultiplied(255, 180, 200, 80)
        };
        painter.rect_filled(highlight_rect, 0.0, drag_color);
    }

    // Detect hovered residue (when not dragging) for Ctrl+hover visual feedback
    let hover_res_idx = if drag_start.is_none() && ctrl_held {
        response.hover_pos().map(|pos| pointer_to_res_idx(pos))
    } else {
        None
    };

    // Pre-filter highlighted residues for this chain to avoid per-residue String allocation
    let highlighted_resvs: HashSet<i32> = highlighted
        .iter()
        .filter(|r| r.object_name == object_name && r.chain_id == chain.chain_id)
        .map(|r| r.resv)
        .collect();

    // Paint residue labels
    for (res_idx, residue) in chain.residues.iter().enumerate() {
        let is_highlighted = highlighted_resvs.contains(&residue.resv);

        let is_in_drag = drag_range
            .as_ref()
            .is_some_and(|range| range.contains(&res_idx));

        let is_ctrl_hovered = hover_res_idx == Some(res_idx);

        let x = seq_rect.min.x + cw * x_offsets[res_idx] as f32;

        // Ctrl+hover background: orange-red rect behind hovered residue
        if is_ctrl_hovered {
            let w = cw * residue.char_width as f32;
            let hover_rect = egui::Rect::from_min_size(
                egui::pos2(x, seq_rect.min.y),
                egui::vec2(w, seq_rect.height()),
            );
            painter.rect_filled(
                hover_rect,
                0.0,
                Color32::from_rgba_premultiplied(255, 80, 30, 120),
            );
        }

        // Font: 12pt for single-character, 10pt for bracket notation
        let font = if residue.char_width == 1 {
            seq_font.clone()
        } else {
            ligand_font.clone()
        };

        let color = if is_ctrl_hovered {
            Color32::WHITE
        } else if is_highlighted || is_in_drag {
            Color32::from_rgb(255, 119, 255)
        } else {
            match residue.kind {
                ResidueKind::AminoAcidCanonical => {
                    resolve_residue_color(residue, named_colors)
                }
                ResidueKind::NucleotideCanonical => NUCLEOTIDE_COLOR,
                ResidueKind::AminoAcidNonCanonical
                | ResidueKind::NucleotideNonCanonical => NONCANONICAL_COLOR,
                ResidueKind::Ligand => LIGAND_COLOR,
            }
        };

        // Bracket labels get a 1px vertical nudge for visual alignment at smaller font size
        let y = if residue.char_width == 1 {
            seq_rect.min.y
        } else {
            seq_rect.min.y + 1.0
        };

        painter.text(
            egui::pos2(x, y),
            egui::Align2::LEFT_TOP,
            &residue.display_label,
            font,
            color,
        );
    }

    // --- Interaction logic ---

    // Drag started: record start position
    if response.drag_started() {
        if let Some(pos) = response.interact_pointer_pos() {
            let res_idx = pointer_to_res_idx(pos);
            *drag_start = Some((obj_idx, chain_idx, res_idx));
            *drag_end = Some(res_idx);
        }
    }

    // Dragging: update end position (only if drag originated in this chain)
    if response.dragged() {
        if let Some((d_obj, d_chain, _)) = *drag_start {
            if d_obj == obj_idx && d_chain == chain_idx {
                if let Some(pos) = response.interact_pointer_pos() {
                    *drag_end = Some(pointer_to_res_idx(pos));
                }
            }
        }
    }

    // Drag stopped: emit selection command and clear drag state
    if response.drag_stopped() {
        if let Some((d_obj, d_chain, d_start_idx)) = *drag_start {
            if d_obj == obj_idx && d_chain == chain_idx {
                if let Some(d_end_idx) = *drag_end {
                    let lo = d_start_idx.min(d_end_idx);
                    let hi = d_start_idx.max(d_end_idx);
                    let resv_values: Vec<i32> =
                        chain.residues[lo..=hi].iter().map(|r| r.resv).collect();
                    let expr = format!(
                        "model {} and chain {} and resi {}",
                        object_name,
                        chain.chain_id,
                        compress_resi_list(&resv_values)
                    );
                    if let Some(cmd) = build_sele_command(&expr, ctrl_held, has_sele) {
                        actions.push(SequenceAction::Execute(cmd));
                    }
                }
                *drag_start = None;
                *drag_end = None;
            }
        }
    }

    // Click (no drag): select or exclude single residue
    if response.clicked() {
        if let Some(pos) = response.interact_pointer_pos() {
            let res_idx = pointer_to_res_idx(pos);
            let residue = &chain.residues[res_idx];
            let expr = format!(
                "model {} and chain {} and resi {}",
                object_name, chain.chain_id, residue.resv
            );
            if let Some(cmd) = build_sele_command(&expr, ctrl_held, has_sele) {
                actions.push(SequenceAction::Execute(cmd));
            }
        }
    }

    // Tooltip on hover (when not dragging) — shown above the sequence row
    // Also signals hover to the 3D viewport for highlighting
    if drag_start.is_none() {
        if let Some(pos) = response.hover_pos() {
            let res_idx = pointer_to_res_idx(pos);
            let residue = &chain.residues[res_idx];
            let tooltip_text = format!("{}.{} {}", chain.chain_id, residue.resn, residue.resv);
            egui::Tooltip::always_open(
                ui.ctx().clone(),
                ui.layer_id(),
                response.id.with("seq_tip"),
                egui::pos2(pos.x, seq_rect.min.y),
            )
            .show(|ui| {
                ui.set_max_width(f32::INFINITY);
                ui.label(&tooltip_text);
            });

            // Signal hover to the 3D viewport
            *hover_out = Some(ResidueRef {
                object_name: object_name.to_string(),
                chain_id: chain.chain_id.clone(),
                resv: residue.resv,
            });
        }
    }

    // Right-click context menu
    response.context_menu(|ui| {
        let has_selection = !highlighted_resvs.is_empty();

        if has_selection {
            if ui.button("Copy selected sequence").clicked() {
                let seq = chain_highlighted_sequence(chain, &highlighted_resvs);
                if !seq.is_empty() {
                    ui.ctx().output_mut(|o| {
                        o.commands
                            .push(egui::OutputCommand::CopyText(seq.clone()));
                    });
                    actions.push(SequenceAction::Notify(format!(
                        "Copied {} residues from chain {}",
                        seq.len(),
                        chain.chain_id
                    )));
                }
                ui.close();
            }
        }

        if ui.button("Copy full chain sequence").clicked() {
            let seq = chain_to_sequence(chain);
            if !seq.is_empty() {
                ui.ctx().output_mut(|o| {
                    o.commands
                        .push(egui::OutputCommand::CopyText(seq.clone()));
                });
                actions.push(SequenceAction::Notify(format!(
                    "Copied {} residues from chain {}",
                    seq.len(),
                    chain.chain_id
                )));
            }
            ui.close();
        }

        ui.separator();

        if ui.button("Copy all chains").clicked() {
            let (text, count) = collect_all_sequences(all_sequences);
            if !text.is_empty() {
                ui.ctx().output_mut(|o| {
                    o.commands
                        .push(egui::OutputCommand::CopyText(text));
                });
                actions.push(SequenceAction::Notify(format!(
                    "Copied {} residues to clipboard",
                    count
                )));
            }
            ui.close();
        }
    });
}

/// Collect all polymer sequences from all objects/chains.
/// Returns `(text, residue_count)` — single chain = plain, multiple = FASTA.
fn collect_all_sequences(sequences: &[SeqObject]) -> (String, usize) {
    let mut entries: Vec<(String, String)> = Vec::new();
    let mut total = 0;
    for seq_obj in sequences {
        for chain in &seq_obj.chains {
            let seq = chain_to_sequence(chain);
            if !seq.is_empty() {
                total += seq.len();
                entries.push((
                    format!("{}:{}", seq_obj.object_name, chain.chain_id),
                    seq,
                ));
            }
        }
    }
    if entries.is_empty() {
        return (String::new(), 0);
    }
    let text = if entries.len() == 1 {
        entries.into_iter().next().unwrap().1
    } else {
        format_fasta(&entries)
    };
    (text, total)
}

/// Extract a one-letter-per-residue sequence string from a chain, skipping ligands.
/// Unknown polymer residues (`display_char == '?'`) become `'X'`.
fn chain_to_sequence(chain: &SeqChain) -> String {
    chain
        .residues
        .iter()
        .filter(|r| r.kind.is_polymer())
        .map(|r| if r.display_char == '?' { 'X' } else { r.display_char })
        .collect()
}

/// Extract sequence of only the highlighted residues from a chain.
fn chain_highlighted_sequence(
    chain: &SeqChain,
    highlighted_resvs: &HashSet<i32>,
) -> String {
    chain
        .residues
        .iter()
        .filter(|r| r.kind.is_polymer())
        .filter(|r| highlighted_resvs.contains(&r.resv))
        .map(|r| if r.display_char == '?' { 'X' } else { r.display_char })
        .collect()
}

/// Format `(header, sequence)` pairs as FASTA with 60-character line wrapping.
fn format_fasta(entries: &[(String, String)]) -> String {
    let mut result = String::new();
    for (header, seq) in entries {
        result.push('>');
        result.push_str(header);
        result.push('\n');
        for chunk in seq.as_bytes().chunks(60) {
            result.push_str(std::str::from_utf8(chunk).unwrap());
            result.push('\n');
        }
    }
    if result.ends_with('\n') {
        result.pop();
    }
    result
}

/// Build sequence data for one molecule object
///
/// ChainIterator groups *consecutive* atoms by chain ID. In mmCIF/PDB files,
/// HETATM records often follow all ATOM records, so the same chain ID (e.g. "A")
/// can appear multiple times — once for protein ATOMs and again for HETATM
/// ligands. We merge these into a single SeqChain so that ligands appear on the
/// same row as the polymer they belong to.
fn build_seq_object(object_name: &str, mol: &ObjectMolecule) -> SeqObject {
    // chain_id → index into `chains`
    let mut chain_index: std::collections::HashMap<String, usize> = Default::default();
    let mut chains: Vec<SeqChain> = Vec::new();

    for chain_view in mol.chains() {
        for residue_view in chain_view.residues() {
            let resn = residue_view.resn();

            let (display_label, kind) = if is_water(resn) || is_ion(resn) {
                continue;
            } else if is_capping_group(resn) {
                // ACE, NME — attached to polymer, no parent letter
                (format!("[{}]", resn), ResidueKind::AminoAcidNonCanonical)
            } else if residue_view.is_protein() {
                if is_standard_amino_acid(resn) {
                    let ch = three_to_one(resn).expect("standard AA must have 1-letter code");
                    (ch.to_string(), ResidueKind::AminoAcidCanonical)
                } else if let Some(ch) = three_to_one(resn) {
                    // Known variant: HID/HIE/HIP → 'H', CYX → 'C', MSE → 'M', etc.
                    (ch.to_string(), ResidueKind::AminoAcidNonCanonical)
                } else {
                    // Truly unknown modified residue with no mapping
                    (format!("[{}]", resn), ResidueKind::AminoAcidNonCanonical)
                }
            } else if residue_view.is_nucleic() {
                if is_standard_nucleotide(resn) {
                    let ch = nucleotide_to_char(resn).expect("standard nucleotide must have 1-letter code");
                    (ch.to_string(), ResidueKind::NucleotideCanonical)
                } else if let Some(ch) = nucleotide_to_char(resn) {
                    // Modified base with known mapping (lowercase)
                    (ch.to_string(), ResidueKind::NucleotideNonCanonical)
                } else {
                    // Unknown modified base
                    (format!("[{}]", resn), ResidueKind::NucleotideNonCanonical)
                }
            } else {
                // Non-polymer: ligand, lipid, unknown HETATM
                (format!("[{}]", resn), ResidueKind::Ligand)
            };

            let display_char = if kind.is_polymer() {
                residue_to_char(resn)
            } else {
                '?'
            };
            let char_width = display_label.len();

            let color_index = residue_view
                .ca()
                .map(|(_, atom)| atom.repr.colors.base)
                .unwrap_or_else(|| {
                    residue_view
                        .iter()
                        .next()
                        .map(|a| a.repr.colors.base)
                        .unwrap_or(-1)
                });

            let seq_residue = SeqResidue {
                chain: residue_view.chain().to_string(),
                resn: resn.to_string(),
                resv: residue_view.resv(),
                display_char,
                display_label,
                char_width,
                atom_range: residue_view.atom_range.clone(),
                color_index,
                kind,
            };

            let chain_id = chain_view.id().to_string();
            if let Some(&idx) = chain_index.get(&chain_id) {
                chains[idx].residues.push(seq_residue);
            } else {
                let idx = chains.len();
                chain_index.insert(chain_id.clone(), idx);
                chains.push(SeqChain {
                    chain_id,
                    residues: vec![seq_residue],
                });
            }
        }
    }

    SeqObject {
        object_name: object_name.to_string(),
        chains,
    }
}

/// Resolve a residue's display color from its color index
fn resolve_residue_color(residue: &SeqResidue, named_colors: &NamedColors) -> Color32 {
    if residue.color_index >= 0 {
        if let Some(color) = named_colors.get_by_index(residue.color_index as u32) {
            return Color32::from_rgb(
                (color.r * 255.0) as u8,
                (color.g * 255.0) as u8,
                (color.b * 255.0) as u8,
            );
        }
    }
    Color32::from_rgb(200, 200, 200)
}
