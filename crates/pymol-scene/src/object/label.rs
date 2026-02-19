//! Label object for 3D text rendering
//!
//! Labels display text in 3D space, positioned at specific coordinates.
//! They can be attached to atoms or free-floating.

use lin_alg::f32::Vec3;
use pymol_render::{LineRep, LineVertex, RenderContext, Representation};
use serde::{Deserialize, Serialize};

use super::{Object, ObjectState, ObjectType};

/// A single label entry
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Label {
    /// Label text
    pub text: String,
    /// Position in world space
    #[serde(with = "crate::serde_helpers::vec3_serde")]
    pub position: Vec3,
    /// Text color (RGBA)
    pub color: [f32; 4],
    /// Font size (in Angstroms)
    pub size: f32,
    /// Screen-space offset (x, y in pixels)
    pub offset: [f32; 2],
    /// Whether to draw a connector line to the position
    pub connector: bool,
    /// Background color (None = transparent)
    pub background: Option<[f32; 4]>,
    /// Associated atom index (if atom-attached)
    pub atom_index: Option<usize>,
    /// Associated object name (if object-attached)
    pub object_name: Option<String>,
}

impl Label {
    /// Create a new label at a position
    pub fn new(text: &str, position: Vec3) -> Self {
        Self {
            text: text.to_string(),
            position,
            color: [1.0, 1.0, 1.0, 1.0],
            size: 14.0,
            offset: [0.0, 0.0],
            connector: false,
            background: None,
            atom_index: None,
            object_name: None,
        }
    }

    /// Create a label attached to an atom
    pub fn atom_label(text: &str, position: Vec3, atom_index: usize, object_name: &str) -> Self {
        Self {
            text: text.to_string(),
            position,
            color: [1.0, 1.0, 1.0, 1.0],
            size: 14.0,
            offset: [0.0, 0.0],
            connector: false,
            background: None,
            atom_index: Some(atom_index),
            object_name: Some(object_name.to_string()),
        }
    }

    /// Set the label color
    pub fn with_color(mut self, color: [f32; 4]) -> Self {
        self.color = color;
        self
    }

    /// Set the font size
    pub fn with_size(mut self, size: f32) -> Self {
        self.size = size;
        self
    }

    /// Set the screen offset
    pub fn with_offset(mut self, offset: [f32; 2]) -> Self {
        self.offset = offset;
        self
    }

    /// Enable connector line
    pub fn with_connector(mut self) -> Self {
        self.connector = true;
        self
    }

    /// Set background color
    pub fn with_background(mut self, color: [f32; 4]) -> Self {
        self.background = Some(color);
        self
    }
}

/// Label anchor point
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum LabelAnchor {
    /// Anchor at bottom-left of text
    #[default]
    BottomLeft,
    /// Anchor at bottom-center
    BottomCenter,
    /// Anchor at bottom-right
    BottomRight,
    /// Anchor at center-left
    CenterLeft,
    /// Anchor at center
    Center,
    /// Anchor at center-right
    CenterRight,
    /// Anchor at top-left
    TopLeft,
    /// Anchor at top-center
    TopCenter,
    /// Anchor at top-right
    TopRight,
}

/// A collection of labels as a scene object
///
/// LabelObject stores multiple labels that can be rendered together.
/// Labels can be atom-attached (following atom positions) or free-floating.
#[derive(Serialize, Deserialize)]
pub struct LabelObject {
    /// Object name
    name: String,
    /// Visual state
    state: ObjectState,
    /// List of labels
    labels: Vec<Label>,
    /// Default anchor point
    anchor: LabelAnchor,
    /// Cached line representation for connectors
    #[serde(skip)]
    connector_lines: Option<LineRep>,
    /// Whether cache needs rebuilding
    #[serde(skip, default = "default_dirty_true")]
    dirty: bool,
}

fn default_dirty_true() -> bool {
    true
}

impl LabelObject {
    /// Create a new empty label object
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            state: ObjectState::default(),
            labels: Vec::new(),
            anchor: LabelAnchor::default(),
            connector_lines: None,
            dirty: true,
        }
    }

    /// Create a label object with initial labels
    pub fn with_labels(name: &str, labels: Vec<Label>) -> Self {
        Self {
            name: name.to_string(),
            state: ObjectState::default(),
            labels,
            anchor: LabelAnchor::default(),
            connector_lines: None,
            dirty: true,
        }
    }

    /// Get the labels
    pub fn labels(&self) -> &[Label] {
        &self.labels
    }

    /// Get mutable access to labels (marks dirty)
    pub fn labels_mut(&mut self) -> &mut Vec<Label> {
        self.dirty = true;
        &mut self.labels
    }

    /// Add a label
    pub fn add_label(&mut self, label: Label) {
        self.labels.push(label);
        self.dirty = true;
    }

    /// Add a simple text label at a position
    pub fn add_text(&mut self, text: &str, position: Vec3) {
        self.labels.push(Label::new(text, position));
        self.dirty = true;
    }

    /// Remove a label by index
    pub fn remove_label(&mut self, index: usize) -> Option<Label> {
        if index < self.labels.len() {
            self.dirty = true;
            Some(self.labels.remove(index))
        } else {
            None
        }
    }

    /// Clear all labels
    pub fn clear(&mut self) {
        self.labels.clear();
        self.dirty = true;
    }

    /// Get the number of labels
    pub fn len(&self) -> usize {
        self.labels.len()
    }

    /// Check if there are no labels
    pub fn is_empty(&self) -> bool {
        self.labels.is_empty()
    }

    /// Get the anchor point
    pub fn anchor(&self) -> LabelAnchor {
        self.anchor
    }

    /// Set the anchor point for all labels
    pub fn set_anchor(&mut self, anchor: LabelAnchor) {
        self.anchor = anchor;
        self.dirty = true;
    }

    /// Check if dirty
    pub fn is_dirty(&self) -> bool {
        self.dirty
    }

    /// Mark as dirty
    pub fn invalidate(&mut self) {
        self.dirty = true;
    }

    /// Build connector lines for labels that have them enabled
    fn build_connectors(&mut self) {
        let _lines = self.connector_lines.get_or_insert_with(LineRep::new);

        // Build line vertices for connectors
        let mut vertices: Vec<LineVertex> = Vec::new();

        for label in &self.labels {
            if label.connector {
                // Draw a small cross at the label position to indicate connector
                // In a full implementation, this would draw from the text to the position
                let p = label.position;
                let size = label.size * 0.1;
                let color = label.color;

                // Horizontal line
                vertices.push(LineVertex {
                    position: [p.x - size, p.y, p.z],
                    color,
                });
                vertices.push(LineVertex {
                    position: [p.x + size, p.y, p.z],
                    color,
                });

                // Vertical line
                vertices.push(LineVertex {
                    position: [p.x, p.y - size, p.z],
                    color,
                });
                vertices.push(LineVertex {
                    position: [p.x, p.y + size, p.z],
                    color,
                });
            }
        }

        // Set the line data directly using internal method
        // Note: LineRep doesn't have set_lines, so we'll need to handle this differently
        // For now, we'll skip the connector rendering until LineRep is extended
    }

    /// Prepare for rendering
    pub fn prepare_render(&mut self, context: &RenderContext) {
        if !self.dirty {
            return;
        }

        self.build_connectors();

        if let Some(ref mut lines) = self.connector_lines {
            lines.upload(context.device(), context.queue());
        }

        self.dirty = false;
    }

    /// Render the labels
    ///
    /// Note: Full text rendering requires a font atlas and text shader.
    /// This implementation renders connector lines; text rendering would
    /// be added with a font atlas system.
    pub fn render<'a>(&'a self, render_pass: &mut wgpu::RenderPass<'a>, context: &'a RenderContext) {
        if !self.state.enabled || self.is_empty() {
            return;
        }

        // Render connector lines
        if let Some(ref lines) = self.connector_lines {
            if !lines.is_empty() {
                let pipeline = context.line_pipeline(pymol_render::pipeline::BlendMode::Opaque);
                render_pass.set_pipeline(&pipeline);
                render_pass.set_bind_group(0, context.uniform_bind_group(), &[]);
                render_pass.set_bind_group(1, context.shadow_bind_group(), &[]);
                lines.render(render_pass);
            }
        }

        // Note: Full text rendering would happen here with a text shader
        // and font atlas. For now, labels are stored but not rendered as text.
    }

    /// Compute bounding box of all label positions
    fn compute_extent(&self) -> Option<(Vec3, Vec3)> {
        if self.labels.is_empty() {
            return None;
        }

        let mut min = Vec3::new(f32::MAX, f32::MAX, f32::MAX);
        let mut max = Vec3::new(f32::MIN, f32::MIN, f32::MIN);

        for label in &self.labels {
            let p = label.position;
            min.x = min.x.min(p.x);
            min.y = min.y.min(p.y);
            min.z = min.z.min(p.z);
            max.x = max.x.max(p.x);
            max.y = max.y.max(p.y);
            max.z = max.z.max(p.z);
        }

        Some((min, max))
    }
}

impl Object for LabelObject {
    fn name(&self) -> &str {
        &self.name
    }

    fn object_type(&self) -> ObjectType {
        // Labels don't have a dedicated ObjectType in the enum yet
        // Using Measurement as it's semantically closest
        ObjectType::Measurement
    }

    fn state(&self) -> &ObjectState {
        &self.state
    }

    fn state_mut(&mut self) -> &mut ObjectState {
        &mut self.state
    }

    fn extent(&self) -> Option<(Vec3, Vec3)> {
        self.compute_extent()
    }

    fn n_states(&self) -> usize {
        1
    }

    fn current_state(&self) -> usize {
        0
    }

    fn set_current_state(&mut self, _state: usize) -> bool {
        false
    }

    fn set_name(&mut self, name: String) {
        self.name = name;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_label_creation() {
        let label = Label::new("Test", Vec3::new(0.0, 0.0, 0.0));
        assert_eq!(label.text, "Test");
        assert_eq!(label.position, Vec3::new(0.0, 0.0, 0.0));
        assert!(!label.connector);
    }

    #[test]
    fn test_label_builder() {
        let label = Label::new("Test", Vec3::new(1.0, 2.0, 3.0))
            .with_color([1.0, 0.0, 0.0, 1.0])
            .with_size(20.0)
            .with_connector()
            .with_background([0.0, 0.0, 0.0, 0.5]);

        assert_eq!(label.color, [1.0, 0.0, 0.0, 1.0]);
        assert_eq!(label.size, 20.0);
        assert!(label.connector);
        assert!(label.background.is_some());
    }

    #[test]
    fn test_atom_label() {
        let label = Label::atom_label("CA", Vec3::new(0.0, 0.0, 0.0), 42, "protein");
        assert_eq!(label.atom_index, Some(42));
        assert_eq!(label.object_name, Some("protein".to_string()));
    }

    #[test]
    fn test_label_object_creation() {
        let obj = LabelObject::new("labels");
        assert_eq!(obj.name(), "labels");
        assert!(obj.is_empty());
        assert!(obj.is_enabled());
    }

    #[test]
    fn test_label_object_add_remove() {
        let mut obj = LabelObject::new("labels");

        obj.add_text("Label 1", Vec3::new(0.0, 0.0, 0.0));
        obj.add_text("Label 2", Vec3::new(1.0, 0.0, 0.0));
        assert_eq!(obj.len(), 2);

        obj.remove_label(0);
        assert_eq!(obj.len(), 1);
        assert_eq!(obj.labels()[0].text, "Label 2");
    }

    #[test]
    fn test_label_object_extent() {
        let mut obj = LabelObject::new("labels");
        obj.add_text("A", Vec3::new(0.0, 0.0, 0.0));
        obj.add_text("B", Vec3::new(10.0, 5.0, 2.0));

        let (min, max) = obj.extent().expect("Should have extent");
        assert_eq!(min, Vec3::new(0.0, 0.0, 0.0));
        assert_eq!(max, Vec3::new(10.0, 5.0, 2.0));
    }

    #[test]
    fn test_label_object_clear() {
        let mut obj = LabelObject::new("labels");
        obj.add_text("A", Vec3::new(0.0, 0.0, 0.0));
        obj.add_text("B", Vec3::new(1.0, 0.0, 0.0));

        obj.clear();
        assert!(obj.is_empty());
    }

    #[test]
    fn test_empty_label_object_extent() {
        let obj = LabelObject::new("labels");
        assert!(obj.extent().is_none());
    }
}
