//! Declarative plugin UI contracts.
//!
//! Plugins describe panel contents with these data types. Frontends own the
//! actual rendering, so plugin UI stays independent of egui, Slint, or any
//! other widget toolkit.

use crate::component::SharedContext;
use crate::message::MessageBus;
use serde::{Deserialize, Serialize};

/// Dock location supported by Patinae plugin panels.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PanelPlacement {
    Right,
    Bottom,
}

/// Static metadata for a plugin panel.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PanelDescriptor {
    pub id: String,
    pub title: String,
    /// Short text/icon shown in the left toolbar.
    pub icon: String,
    pub placement: PanelPlacement,
    pub default_visible: bool,
}

impl PanelDescriptor {
    pub fn right(id: impl Into<String>, title: impl Into<String>) -> Self {
        let title = title.into();
        Self {
            id: id.into(),
            icon: default_icon(&title),
            title,
            placement: PanelPlacement::Right,
            default_visible: false,
        }
    }

    pub fn bottom(id: impl Into<String>, title: impl Into<String>) -> Self {
        let title = title.into();
        Self {
            id: id.into(),
            icon: default_icon(&title),
            title,
            placement: PanelPlacement::Bottom,
            default_visible: false,
        }
    }

    pub fn icon(mut self, icon: impl Into<String>) -> Self {
        self.icon = icon.into();
        self
    }

    pub fn default_visible(mut self, visible: bool) -> Self {
        self.default_visible = visible;
        self
    }
}

fn default_icon(title: &str) -> String {
    title
        .chars()
        .find(|c| !c.is_whitespace())
        .map(|c| c.to_uppercase().to_string())
        .unwrap_or_else(|| "+".to_string())
}

/// One selectable option for segmented/dropdown-style controls.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PanelOption {
    pub label: String,
    pub value: String,
}

impl PanelOption {
    pub fn new(label: impl Into<String>, value: impl Into<String>) -> Self {
        Self {
            label: label.into(),
            value: value.into(),
        }
    }
}

/// A compact icon button used inside toolbar-style rows.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PanelButton {
    pub id: String,
    pub label: String,
    pub icon: String,
    pub primary: bool,
    pub enabled: bool,
}

impl PanelButton {
    pub fn new(
        id: impl Into<String>,
        label: impl Into<String>,
        icon: impl Into<String>,
        primary: bool,
    ) -> Self {
        Self {
            id: id.into(),
            label: label.into(),
            icon: icon.into(),
            primary,
            enabled: true,
        }
    }

    pub fn enabled(mut self, enabled: bool) -> Self {
        self.enabled = enabled;
        self
    }
}

/// Semantic style for highlighted text ranges in plugin text areas.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PanelTextStyle {
    Keyword,
    String,
    Comment,
    Number,
    Function,
    Type,
    Constant,
    Operator,
    Punctuation,
    Builtin,
}

impl PanelTextStyle {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Keyword => "keyword",
            Self::String => "string",
            Self::Comment => "comment",
            Self::Number => "number",
            Self::Function => "function",
            Self::Type => "type",
            Self::Constant => "constant",
            Self::Operator => "operator",
            Self::Punctuation => "punctuation",
            Self::Builtin => "builtin",
        }
    }
}

/// A highlighted byte range inside a plugin text area.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct PanelTextHighlight {
    pub start: usize,
    pub end: usize,
    pub style: PanelTextStyle,
}

impl PanelTextHighlight {
    pub fn new(start: usize, end: usize, style: PanelTextStyle) -> Self {
        Self { start, end, style }
    }
}

/// A multi-line text surface used by scripting/editor panels.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PanelTextArea {
    pub id: String,
    pub label: String,
    pub value: String,
    pub placeholder: String,
    pub rows: u32,
    pub read_only: bool,
    pub highlights: Vec<PanelTextHighlight>,
}

impl PanelTextArea {
    pub fn new(
        id: impl Into<String>,
        label: impl Into<String>,
        value: impl Into<String>,
        placeholder: impl Into<String>,
        rows: u32,
        read_only: bool,
    ) -> Self {
        Self {
            id: id.into(),
            label: label.into(),
            value: value.into(),
            placeholder: placeholder.into(),
            rows,
            read_only,
            highlights: Vec::new(),
        }
    }

    pub fn with_highlights(mut self, highlights: Vec<PanelTextHighlight>) -> Self {
        self.highlights = highlights;
        self
    }
}

/// A child control inside a generic panel layout container.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PanelControlNode {
    pub control: PanelControl,
    pub grow: f32,
}

impl PanelControlNode {
    pub fn new(control: PanelControl) -> Self {
        Self { control, grow: 0.0 }
    }

    pub fn grow(mut self, grow: f32) -> Self {
        self.grow = grow;
        self
    }
}

/// A horizontal group of plugin controls.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PanelRow {
    pub id: String,
    pub children: Vec<PanelControlNode>,
    pub gap: f32,
}

impl PanelRow {
    pub fn new(id: impl Into<String>, children: Vec<PanelControlNode>) -> Self {
        Self {
            id: id.into(),
            children,
            gap: 8.0,
        }
    }

    pub fn gap(mut self, gap: f32) -> Self {
        self.gap = gap;
        self
    }
}

/// A vertical group of plugin controls.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PanelColumn {
    pub id: String,
    pub children: Vec<PanelControlNode>,
    pub gap: f32,
}

impl PanelColumn {
    pub fn new(id: impl Into<String>, children: Vec<PanelControlNode>) -> Self {
        Self {
            id: id.into(),
            children,
            gap: 4.0,
        }
    }

    pub fn gap(mut self, gap: f32) -> Self {
        self.gap = gap;
        self
    }
}

/// A titled visual group of plugin controls.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PanelGroup {
    pub id: String,
    pub title: String,
    pub open: bool,
    pub children: Vec<PanelControlNode>,
    pub gap: f32,
}

impl PanelGroup {
    pub fn new(
        id: impl Into<String>,
        title: impl Into<String>,
        children: Vec<PanelControlNode>,
    ) -> Self {
        Self {
            id: id.into(),
            title: title.into(),
            open: true,
            children,
            gap: 8.0,
        }
    }

    pub fn open(mut self, open: bool) -> Self {
        self.open = open;
        self
    }

    pub fn gap(mut self, gap: f32) -> Self {
        self.gap = gap;
        self
    }
}

/// Controls supported by the v1 declarative renderer.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum PanelControl {
    Text {
        id: String,
        text: String,
    },
    /// A bold title with a muted description beneath — a compact header block
    /// that pairs well with an action button to its right.
    TitleDesc {
        id: String,
        title: String,
        desc: String,
    },
    Heading {
        id: String,
        text: String,
    },
    Section {
        id: String,
        title: String,
        open: bool,
    },
    Button {
        id: String,
        label: String,
        primary: bool,
    },
    ButtonRow {
        id: String,
        buttons: Vec<PanelButton>,
    },
    Toggle {
        id: String,
        label: String,
        value: bool,
    },
    Slider {
        id: String,
        label: String,
        value: f32,
        min: f32,
        max: f32,
        step: f32,
    },
    Number {
        id: String,
        label: String,
        value: f32,
        min: f32,
        max: f32,
        step: f32,
    },
    Select {
        id: String,
        label: String,
        value: String,
        options: Vec<PanelOption>,
    },
    TextInput {
        id: String,
        label: String,
        value: String,
        placeholder: String,
    },
    TextArea(PanelTextArea),
    Row(PanelRow),
    Column(PanelColumn),
    Group(PanelGroup),
    Image {
        id: String,
        width: u32,
        height: u32,
        rgba: Vec<u8>,
    },
    Spacer {
        id: String,
        height: f32,
    },
}

/// A complete render snapshot for one panel.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct PanelSnapshot {
    pub controls: Vec<PanelControl>,
}

impl PanelSnapshot {
    pub fn new(controls: Vec<PanelControl>) -> Self {
        Self { controls }
    }
}

/// Runtime value sent from the frontend to a plugin panel.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum PanelValue {
    None,
    Bool(bool),
    Number(f32),
    Text(String),
}

impl PanelValue {
    pub fn as_command_value(&self) -> String {
        match self {
            Self::None => String::new(),
            Self::Bool(v) => {
                if *v {
                    "1".to_string()
                } else {
                    "0".to_string()
                }
            }
            Self::Number(v) => format!("{v:.6}")
                .trim_end_matches('0')
                .trim_end_matches('.')
                .to_string(),
            Self::Text(v) => v.clone(),
        }
    }
}

/// Frontend event kind emitted by rendered plugin controls.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PanelEventKind {
    Click,
    Toggle,
    NumberChange,
    Select,
    TextEdit,
    TextAreaEdit,
    TextCommit,
}

impl PanelEventKind {
    pub fn refreshes_snapshot(self) -> bool {
        !matches!(self, Self::TextEdit)
    }
}

/// An event emitted by the rendered panel.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PanelEvent {
    pub panel_id: String,
    pub control_id: String,
    pub kind: PanelEventKind,
    pub value: PanelValue,
}

/// High-level actions returned by plugin panels.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum PanelAction {
    ExecuteCommand { command: String, silent: bool },
    SetSetting { name: String, value: PanelValue },
    Custom { topic: String, payload: Vec<u8> },
}

/// Runtime host inputs requested by a plugin panel.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PanelRuntimeRequirements {
    bits: u64,
}

impl PanelRuntimeRequirements {
    /// No extra host runtime inputs are required.
    pub const NONE: Self = Self { bits: 0 };

    /// The panel needs a full serialized session.
    pub const FULL_SESSION: Self = Self { bits: 1 << 0 };

    /// Builds requirements from raw ABI bits.
    pub const fn from_bits(bits: u64) -> Self {
        Self { bits }
    }

    /// Returns the raw ABI bitset.
    pub const fn bits(self) -> u64 {
        self.bits
    }

    /// Returns true when all requested bits are present.
    pub const fn contains(self, other: Self) -> bool {
        (self.bits & other.bits) == other.bits
    }

    /// Returns the union of two requirement sets.
    pub const fn union(self, other: Self) -> Self {
        Self {
            bits: self.bits | other.bits,
        }
    }

    /// Returns true when no extra runtime inputs are required.
    pub const fn is_empty(self) -> bool {
        self.bits == 0
    }
}

/// A UI panel implemented by a plugin.
pub trait PluginPanel: Send {
    fn descriptor(&self) -> PanelDescriptor;

    fn runtime_requirements(&self) -> PanelRuntimeRequirements {
        PanelRuntimeRequirements::FULL_SESSION
    }

    fn snapshot(&mut self, ctx: &SharedContext<'_>) -> PanelSnapshot;

    fn handle_event(
        &mut self,
        _event: PanelEvent,
        _ctx: &SharedContext<'_>,
        _bus: &mut MessageBus,
    ) -> Vec<PanelAction> {
        Vec::new()
    }
}

#[cfg(test)]
mod tests {
    use super::PanelEventKind;

    #[test]
    fn text_edit_is_the_only_non_refreshing_panel_event() {
        assert!(!PanelEventKind::TextEdit.refreshes_snapshot());
        assert!(PanelEventKind::Click.refreshes_snapshot());
        assert!(PanelEventKind::Toggle.refreshes_snapshot());
        assert!(PanelEventKind::NumberChange.refreshes_snapshot());
        assert!(PanelEventKind::Select.refreshes_snapshot());
        assert!(PanelEventKind::TextAreaEdit.refreshes_snapshot());
        assert!(PanelEventKind::TextCommit.refreshes_snapshot());
    }
}
