//! Sequence Viewer Panel
//!
//! Displays amino acid / nucleotide sequences of loaded structures with
//! interactive selection synchronized with the 3D view.

use std::collections::HashSet;
use std::ops::Range;

use egui::{Color32, RichText, Ui};
use pymol_color::NamedColors;
use pymol_mol::{residue_to_char, ObjectMolecule};
use pymol_scene::ObjectRegistry;

/// Single residue entry in the sequence viewer
#[derive(Debug, Clone)]
pub struct SeqResidue {
    pub chain: String,
    pub resn: String,
    pub resv: i32,
    pub display_char: char,
    pub atom_range: Range<usize>,
    /// Base color index of the CA atom (or first atom), for display coloring
    pub color_index: i32,
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

/// Action produced by the sequence panel
pub enum SequenceAction {
    /// Execute a command string (e.g., "select sele, ...")
    Execute(String),
}

/// Persistent state for the sequence viewer
pub struct SequencePanel {
    /// Cached sequence data per object
    sequences: Vec<SeqObject>,
    /// Number of objects when cache was last built
    cached_object_count: usize,
    /// Currently highlighted residues from 3D selection: (object_name, chain_id, resv)
    highlighted: HashSet<(String, String, i32)>,
    /// Last clicked residue for shift-click range selection: (object_index, chain_index, residue_index)
    last_clicked: Option<(usize, usize, usize)>,
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
            highlighted: HashSet::new(),
            last_clicked: None,
        }
    }

    /// Update the highlighted residues from the current 3D selection.
    pub fn update_highlights(&mut self, highlighted: HashSet<(String, String, i32)>) {
        self.highlighted = highlighted;
    }

    /// Rebuild the sequence cache from the registry
    fn rebuild_cache(&mut self, registry: &ObjectRegistry) {
        self.sequences.clear();

        for name in registry.names() {
            if let Some(mol_obj) = registry.get_molecule(name) {
                let mol = mol_obj.molecule();
                let seq_obj = build_seq_object(name, mol);
                if !seq_obj.chains.is_empty() {
                    self.sequences.push(seq_obj);
                }
            }
        }

        self.cached_object_count = registry.names().count();
    }

    /// Check if cache needs rebuilding (simple heuristic: object count changed)
    fn needs_rebuild(&self, registry: &ObjectRegistry) -> bool {
        let count = registry.names().count();
        count != self.cached_object_count
    }

    /// Show the sequence panel and return any actions
    pub fn show(
        &mut self,
        ui: &mut Ui,
        registry: &ObjectRegistry,
        named_colors: &NamedColors,
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
        let last_clicked = &mut self.last_clicked;

        // Find the longest chain to size the scrollable content
        let max_residues = self
            .sequences
            .iter()
            .flat_map(|s| s.chains.iter())
            .map(|c| c.residues.len())
            .max()
            .unwrap_or(0);

        let cw = char_width(ui);

        // Two-column layout: fixed labels on the left, single shared horizontal
        // scroll area on the right containing all chains' sequences.
        ui.horizontal(|ui| {
            // Left column: object/chain labels (fixed, not scrolled)
            ui.vertical(|ui| {
                ui.spacing_mut().item_spacing.y = 0.0;
                for seq_obj in self.sequences.iter() {
                    // Object name row
                    ui.label(
                        RichText::new(&seq_obj.object_name)
                            .strong()
                            .color(Color32::from_rgb(100, 255, 100))
                            .size(11.0),
                    );
                    for chain in seq_obj.chains.iter() {
                        // Ruler placeholder row (same height as the ruler)
                        ui.allocate_exact_size(egui::vec2(0.0, 10.0), egui::Sense::hover());
                        // Sequence row label
                        ui.label(
                            RichText::new(format!("{}:", chain.chain_id))
                                .color(Color32::from_rgb(180, 180, 180))
                                .family(egui::FontFamily::Monospace)
                                .size(12.0),
                        );
                    }
                    ui.add_space(2.0);
                }
            });

            // Right column: single horizontal scroll for all sequences
            egui::ScrollArea::horizontal()
                .id_salt("sequence_hscroll")
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.vertical(|ui| {
                        ui.spacing_mut().item_spacing.y = 0.0;
                        for (obj_idx, seq_obj) in self.sequences.iter().enumerate() {
                            // Object header spacer (match left column)
                            ui.allocate_exact_size(
                                egui::vec2(cw * max_residues as f32, ui.text_style_height(&egui::TextStyle::Body)),
                                egui::Sense::hover(),
                            );

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
                                    last_clicked,
                                    &mut actions,
                                );
                            }
                            ui.add_space(2.0);
                        }
                    });
                });
        });

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
    highlighted: &HashSet<(String, String, i32)>,
    last_clicked: &mut Option<(usize, usize, usize)>,
    actions: &mut Vec<SequenceAction>,
) {
    let n = chain.residues.len();

    // --- Ruler row: residue numbers above the sequence ---
    let (ruler_rect, _) =
        ui.allocate_exact_size(egui::vec2(cw * n as f32, 10.0), egui::Sense::hover());
    let painter = ui.painter_at(ruler_rect);
    let ruler_font = egui::FontId::new(8.0, egui::FontFamily::Monospace);
    let ruler_color = Color32::from_rgb(120, 120, 120);

    for (res_idx, residue) in chain.residues.iter().enumerate() {
        if res_idx == 0 || res_idx % 10 == 0 {
            let x = ruler_rect.min.x + cw * res_idx as f32;
            painter.text(
                egui::pos2(x, ruler_rect.min.y),
                egui::Align2::LEFT_TOP,
                format!("{}", residue.resv),
                ruler_font.clone(),
                ruler_color,
            );
        }
    }

    // --- Sequence row: clickable residue letters ---
    ui.horizontal(|ui| {
        ui.spacing_mut().item_spacing.x = 0.0;
        ui.spacing_mut().button_padding = egui::Vec2::ZERO;

        for (res_idx, residue) in chain.residues.iter().enumerate() {
            let is_highlighted = highlighted.contains(&(
                object_name.to_string(),
                chain.chain_id.clone(),
                residue.resv,
            ));

            let text_color = resolve_residue_color(residue, named_colors);

            let mut text = RichText::new(residue.display_char.to_string())
                .family(egui::FontFamily::Monospace)
                .size(12.0);

            if is_highlighted {
                text = text.color(Color32::from_rgb(255, 80, 80)).strong();
            } else {
                text = text.color(text_color);
            }

            let fill = Color32::TRANSPARENT;

            let button = egui::Button::new(text)
                .fill(fill)
                .frame(false)
                .min_size(egui::vec2(cw, 0.0));

            let response = ui.add(button);

            response.clone().on_hover_text(format!(
                "{}.{} {}",
                chain.chain_id, residue.resn, residue.resv
            ));

            if response.clicked() {
                let modifiers = ui.input(|i| i.modifiers);

                if modifiers.shift {
                    if let Some((prev_obj, prev_chain, prev_res)) = *last_clicked {
                        if prev_obj == obj_idx && prev_chain == chain_idx {
                            let start = prev_res.min(res_idx);
                            let end = prev_res.max(res_idx);
                            let resi_list: Vec<String> = chain.residues[start..=end]
                                .iter()
                                .map(|r| r.resv.to_string())
                                .collect();
                            let expr = format!(
                                "model {} and chain {} and resi {}",
                                object_name,
                                chain.chain_id,
                                resi_list.join("+")
                            );
                            actions.push(SequenceAction::Execute(format!(
                                "select sele, {}",
                                expr
                            )));
                        }
                    }
                } else if modifiers.command || modifiers.ctrl {
                    let expr = format!(
                        "model {} and chain {} and resi {}",
                        object_name, chain.chain_id, residue.resv
                    );
                    if is_highlighted {
                        actions.push(SequenceAction::Execute(format!(
                            "select sele, sele and not ({})",
                            expr
                        )));
                    } else {
                        actions.push(SequenceAction::Execute(format!(
                            "select sele, sele or ({})",
                            expr
                        )));
                    }
                } else {
                    let expr = format!(
                        "model {} and chain {} and resi {}",
                        object_name, chain.chain_id, residue.resv
                    );
                    actions.push(SequenceAction::Execute(format!(
                        "select sele, {}",
                        expr
                    )));
                }

                *last_clicked = Some((obj_idx, chain_idx, res_idx));
            }
        }
    });
}

/// Build sequence data for one molecule object
fn build_seq_object(object_name: &str, mol: &ObjectMolecule) -> SeqObject {
    let mut chains: Vec<SeqChain> = Vec::new();

    for chain_view in mol.chains() {
        let mut residues = Vec::new();

        for residue_view in chain_view.residues() {
            if !residue_view.is_protein() && !residue_view.is_nucleic() {
                continue;
            }

            let display_char = residue_to_char(residue_view.resn());

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

            residues.push(SeqResidue {
                chain: residue_view.chain().to_string(),
                resn: residue_view.resn().to_string(),
                resv: residue_view.resv(),
                display_char,
                atom_range: residue_view.atom_range.clone(),
                color_index,
            });
        }

        if !residues.is_empty() {
            chains.push(SeqChain {
                chain_id: chain_view.id().to_string(),
                residues,
            });
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
