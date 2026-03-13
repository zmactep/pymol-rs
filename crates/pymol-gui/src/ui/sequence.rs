//! Sequence Viewer Panel
//!
//! Displays amino acid / nucleotide sequences of loaded structures with
//! interactive selection synchronized with the 3D view.
//! Sends `AppMessage` directly to the `MessageBus`.

use std::collections::HashSet;

use egui::{Color32, RichText, Ui};
use pymol_color::NamedColors;
use pymol_scene::ObjectRegistry;
use pymol_select::build_sele_command;

use pymol_framework::message::MessageBus;
use crate::model::sequence::{
    ResidueRef, ResidueKind, SeqChain, SeqObject, SequenceModel,
    chain_highlighted_sequence, chain_to_sequence, collect_all_sequences,
    compress_resi_list,
};

/// Pale magenta — standard DNA/RNA nucleotides
const NUCLEOTIDE_COLOR: Color32 = Color32::from_rgb(205, 150, 205);

/// Orange — non-canonical AAs, AA variants (HIP/CYX/MSE/…), modified bases, caps
const NONCANONICAL_COLOR: Color32 = Color32::from_rgb(230, 160, 60);

/// Dim gray — non-polymer HETATM ligands
const LIGAND_COLOR: Color32 = Color32::from_rgb(140, 140, 140);

// =============================================================================
// UI State
// =============================================================================

/// Sequence panel UI state (egui-specific: drag, hover)
pub struct SequenceUiState {
    /// Active drag start: (object_index, chain_index, start_residue_index)
    pub drag_start: Option<(usize, usize, usize)>,
    /// Current drag end residue index (updated each frame while dragging)
    pub drag_end: Option<usize>,
    /// Current sequence hover state for change detection
    pub current_hover: Option<ResidueRef>,
}

impl Default for SequenceUiState {
    fn default() -> Self {
        Self::new()
    }
}

impl SequenceUiState {
    pub fn new() -> Self {
        Self {
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
}

// =============================================================================
// View
// =============================================================================

/// Sequence viewer panel
pub struct SequencePanel;

impl SequencePanel {
    /// Show the sequence panel, sending messages directly to the bus.
    pub fn show(
        ui: &mut Ui,
        model: &mut SequenceModel,
        ui_state: &mut SequenceUiState,
        registry: &ObjectRegistry,
        named_colors: &NamedColors,
        has_sele: bool,
        bus: &mut MessageBus,
    ) {
        if model.needs_rebuild(registry) {
            model.rebuild_cache(registry);
        }

        if model.sequences.is_empty() {
            ui.centered_and_justified(|ui| {
                ui.label(
                    RichText::new("No sequences loaded")
                        .color(Color32::GRAY)
                        .italics(),
                );
            });
            return;
        }

        let highlighted = &model.highlighted;
        let drag_start = &mut ui_state.drag_start;
        let drag_end = &mut ui_state.drag_end;

        let cw = char_width(ui);
        let ruler_height = 10.0;
        let row_height = ui.text_style_height(&egui::TextStyle::Body);
        let obj_name_font = egui::FontId::new(11.0, egui::FontFamily::Proportional);
        let obj_name_color = Color32::from_rgb(100, 255, 100);
        let label_font = egui::FontId::new(12.0, egui::FontFamily::Monospace);
        let label_color = Color32::from_rgb(180, 180, 180);
        let separator_color = Color32::from_rgb(60, 60, 60);
        let obj_name_width = 40.0;
        let label_width = cw * 2.5;
        let spacing = 2.0;
        let n_objects = model.sequences.len();

        // Pre-compute per-object heights
        let obj_heights: Vec<f32> = model
            .sequences
            .iter()
            .map(|seq_obj| seq_obj.chains.len() as f32 * (ruler_height + row_height))
            .collect();

        // Track which residue the mouse hovers over
        let mut hover_out: Option<ResidueRef> = None;

        egui::ScrollArea::vertical()
            .id_salt("sequence_vscroll")
            .auto_shrink([false, true])
            .show(ui, |ui| {
        // Three-column layout
        ui.horizontal(|ui| {
            // Column 1: object names
            ui.vertical(|ui| {
                ui.spacing_mut().item_spacing.y = 0.0;
                for (i, seq_obj) in model.sequences.iter().enumerate() {
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

            // Column 2: chain labels
            ui.vertical(|ui| {
                ui.spacing_mut().item_spacing.y = 0.0;
                for seq_obj in model.sequences.iter() {
                    for chain in seq_obj.chains.iter() {
                        ui.allocate_exact_size(
                            egui::vec2(label_width, ruler_height),
                            egui::Sense::hover(),
                        );
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

            // Column 3: scrollable sequences
            egui::ScrollArea::horizontal()
                .id_salt("sequence_hscroll")
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.vertical(|ui| {
                        ui.spacing_mut().item_spacing.y = 0.0;
                        for (obj_idx, seq_obj) in model.sequences.iter().enumerate() {
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
                                    &model.sequences,
                                    drag_start,
                                    drag_end,
                                    bus,
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

        // Track hover state changes (synced to viewport after rendering)
        if hover_out != ui_state.current_hover {
            ui_state.current_hover = hover_out;
        }
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
#[allow(clippy::too_many_arguments)]
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
    bus: &mut MessageBus,
    has_sele: bool,
    hover_out: &mut Option<ResidueRef>,
) {
    let n = chain.residues.len();
    if n == 0 {
        return;
    }

    // Precompute cumulative x-offsets
    let mut x_offsets: Vec<usize> = Vec::with_capacity(n);
    let mut cumulative = 0usize;
    for residue in &chain.residues {
        x_offsets.push(cumulative);
        cumulative += residue.char_width;
    }
    let total_cells = cumulative;
    let total_width = cw * total_cells as f32;

    // --- Ruler row ---
    let (ruler_rect, _) =
        ui.allocate_exact_size(egui::vec2(total_width, 10.0), egui::Sense::hover());
    let ruler_painter = ui.painter_at(ruler_rect);
    let ruler_font = egui::FontId::new(8.0, egui::FontFamily::Monospace);
    let ruler_color = Color32::from_rgb(120, 120, 120);

    for (res_idx, residue) in chain.residues.iter().enumerate() {
        if res_idx == 0 || res_idx % 10 == 0 {
            let x = ruler_rect.min.x + cw * x_offsets[res_idx] as f32;
            let label = format!("{}", residue.resv);
            let label_width = label.len() as f32 * cw * 0.7;
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

    // --- Sequence row ---
    let row_height = ui.text_style_height(&egui::TextStyle::Body);
    let (seq_rect, response) = ui.allocate_exact_size(
        egui::vec2(total_width, row_height),
        egui::Sense::click_and_drag(),
    );

    let painter = ui.painter_at(seq_rect);
    let seq_font = egui::FontId::new(12.0, egui::FontFamily::Monospace);
    let ligand_font = egui::FontId::new(10.0, egui::FontFamily::Monospace);

    let ctrl_held = ui.input(|i| i.modifiers.ctrl || i.modifiers.command);

    let pointer_to_res_idx = |pos: egui::Pos2| -> usize {
        let rel_x = pos.x - seq_rect.min.x;
        let cell = (rel_x / cw).floor().max(0.0) as usize;
        match x_offsets.binary_search(&cell) {
            Ok(idx) => idx.min(n - 1),
            Err(idx) => idx.saturating_sub(1).min(n - 1),
        }
    };

    // Determine the active drag range for this chain
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

    // Detect hovered residue
    let hover_res_idx = if drag_start.is_none() && ctrl_held {
        response.hover_pos().map(pointer_to_res_idx)
    } else {
        None
    };

    // Pre-filter highlighted residues for this chain
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

    if response.drag_started() {
        if let Some(pos) = response.interact_pointer_pos() {
            let res_idx = pointer_to_res_idx(pos);
            *drag_start = Some((obj_idx, chain_idx, res_idx));
            *drag_end = Some(res_idx);
        }
    }

    if response.dragged() {
        if let Some((d_obj, d_chain, _)) = *drag_start {
            if d_obj == obj_idx && d_chain == chain_idx {
                if let Some(pos) = response.interact_pointer_pos() {
                    *drag_end = Some(pointer_to_res_idx(pos));
                }
            }
        }
    }

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
                        bus.execute_command(cmd);
                    }
                }
                *drag_start = None;
                *drag_end = None;
            }
        }
    }

    if response.clicked() {
        if let Some(pos) = response.interact_pointer_pos() {
            let res_idx = pointer_to_res_idx(pos);
            let residue = &chain.residues[res_idx];
            let expr = format!(
                "model {} and chain {} and resi {}",
                object_name, chain.chain_id, residue.resv
            );
            if let Some(cmd) = build_sele_command(&expr, ctrl_held, has_sele) {
                bus.execute_command(cmd);
            }
        }
    }

    // Tooltip on hover
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

        if has_selection
            && ui.button("Copy selected sequence").clicked()
        {
            let seq = chain_highlighted_sequence(chain, &highlighted_resvs);
            if !seq.is_empty() {
                ui.ctx().output_mut(|o| {
                    o.commands
                        .push(egui::OutputCommand::CopyText(seq.clone()));
                });
                bus.print_info(format!(
                    "Copied {} residues from chain {}",
                    seq.len(),
                    chain.chain_id
                ));
            }
            ui.close();
        }

        if ui.button("Copy full chain sequence").clicked() {
            let seq = chain_to_sequence(chain);
            if !seq.is_empty() {
                ui.ctx().output_mut(|o| {
                    o.commands
                        .push(egui::OutputCommand::CopyText(seq.clone()));
                });
                bus.print_info(format!(
                    "Copied {} residues from chain {}",
                    seq.len(),
                    chain.chain_id
                ));
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
                bus.print_info(format!(
                    "Copied {} residues to clipboard",
                    count
                ));
            }
            ui.close();
        }
    });
}

/// Resolve a residue's display color from its color index
fn resolve_residue_color(residue: &crate::model::SeqResidue, named_colors: &NamedColors) -> Color32 {
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
