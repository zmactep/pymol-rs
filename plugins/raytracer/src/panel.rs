//! Declarative Patinae panel for ray tracing controls.

use patinae_plugin::prelude::*;
use patinae_settings::SettingValue;

const PANEL_ID: &str = "rt_toolbar";
const SAVE_REPLY_CONTROL_ID: &str = "save_file_selected";
const DEFAULT_SAVE_FILE: &str = "raytrace.png";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ResolutionPreset {
    Hd720,
    Hd1080,
    Uhd4K,
    Custom,
}

impl ResolutionPreset {
    fn label(self) -> &'static str {
        match self {
            Self::Hd720 => "720p",
            Self::Hd1080 => "1080p",
            Self::Uhd4K => "4K",
            Self::Custom => "Custom",
        }
    }

    fn value(self) -> &'static str {
        match self {
            Self::Hd720 => "720",
            Self::Hd1080 => "1080",
            Self::Uhd4K => "4k",
            Self::Custom => "custom",
        }
    }

    fn dimensions(self) -> Option<(u32, u32)> {
        match self {
            Self::Hd720 => Some((1280, 720)),
            Self::Hd1080 => Some((1920, 1080)),
            Self::Uhd4K => Some((3840, 2160)),
            Self::Custom => None,
        }
    }

    fn from_dimensions(w: u32, h: u32) -> Self {
        match (w, h) {
            (1280, 720) => Self::Hd720,
            (1920, 1080) => Self::Hd1080,
            (3840, 2160) => Self::Uhd4K,
            _ => Self::Custom,
        }
    }

    fn from_value(value: &str) -> Self {
        match value {
            "720" => Self::Hd720,
            "1080" => Self::Hd1080,
            "4k" => Self::Uhd4K,
            _ => Self::Custom,
        }
    }
}

pub(crate) struct RtPanel {
    width: u32,
    height: u32,
    preset: ResolutionPreset,
    antialias: u32,
    output_open: bool,
    trace_open: bool,
    quality_open: bool,
    lighting_open: bool,
    edge_open: bool,
    use_custom: bool,
    ambient: f32,
    direct: f32,
    reflect: f32,
    specular: f32,
    shininess: f32,
    shadow: bool,
    transparency_shadows: bool,
    max_passes: i32,
    mode: i32,
    fog: f32,
    slope_factor: f32,
    depth_factor: f32,
    disco_factor: f32,
    gain: f32,
    opaque_background: i32,
}

impl Default for RtPanel {
    fn default() -> Self {
        Self {
            width: 1920,
            height: 1080,
            preset: ResolutionPreset::Hd1080,
            antialias: 2,
            output_open: true,
            trace_open: true,
            quality_open: true,
            lighting_open: false,
            edge_open: false,
            use_custom: false,
            ambient: 0.14,
            direct: 0.45,
            reflect: 0.45,
            specular: 0.5,
            shininess: 40.0,
            shadow: true,
            transparency_shadows: true,
            max_passes: 25,
            mode: 0,
            fog: -1.0,
            slope_factor: 0.6,
            depth_factor: 0.1,
            disco_factor: 0.05,
            gain: 0.12,
            opaque_background: -1,
        }
    }
}

impl RtPanel {
    pub(crate) fn new() -> Self {
        Self::default()
    }

    fn sync_from_settings(&mut self, ctx: &SharedContext<'_>) {
        let read = |name: &str| -> Option<SettingValue> {
            let entry = ctx.dynamic_settings?.lookup(name)?;
            let store = entry.store.read().ok()?;
            store.get(name).cloned()
        };

        macro_rules! sync_bool {
            ($field:ident, $name:expr) => {
                if let Some(SettingValue::Bool(v)) = read($name) {
                    self.$field = v;
                }
            };
        }
        macro_rules! sync_int {
            ($field:ident, $name:expr) => {
                if let Some(SettingValue::Int(v)) = read($name) {
                    self.$field = v;
                }
            };
        }
        macro_rules! sync_float {
            ($field:ident, $name:expr) => {
                if let Some(SettingValue::Float(v)) = read($name) {
                    self.$field = v;
                }
            };
        }

        sync_bool!(use_custom, "rt_use_custom");
        sync_float!(ambient, "rt_ambient");
        sync_float!(direct, "rt_direct");
        sync_float!(reflect, "rt_reflect");
        sync_float!(specular, "rt_specular");
        sync_float!(shininess, "rt_shininess");
        sync_bool!(shadow, "ray_shadow");
        sync_bool!(transparency_shadows, "ray_transparency_shadows");
        sync_int!(max_passes, "ray_max_passes");
        sync_int!(mode, "ray_trace_mode");
        sync_float!(fog, "ray_trace_fog");
        sync_float!(slope_factor, "ray_trace_slope_factor");
        sync_float!(depth_factor, "ray_trace_depth_factor");
        sync_float!(disco_factor, "ray_trace_disco_factor");
        sync_float!(gain, "ray_trace_gain");
        sync_int!(opaque_background, "ray_opaque_background");
    }

    fn controls(&mut self) -> Vec<PanelControl> {
        let mut output_controls = vec![node(PanelControl::Select {
            id: "preset".into(),
            label: "Resolution".into(),
            value: self.preset.value().into(),
            options: ResolutionPreset::ALL
                .iter()
                .map(|p| PanelOption::new(p.label(), p.value()))
                .collect(),
        })];
        if self.preset == ResolutionPreset::Custom {
            output_controls.push(node(row_with_gap(
                "custom_dimensions",
                vec![
                    node(number(
                        "width",
                        "Width",
                        self.width as f32,
                        64.0,
                        7680.0,
                        1.0,
                    ))
                    .grow(1.0),
                    node(number(
                        "height",
                        "Height",
                        self.height as f32,
                        64.0,
                        4320.0,
                        1.0,
                    ))
                    .grow(1.0),
                ],
                8.0,
            )));
        }
        output_controls.push(node(slider(
            "antialias",
            "Antialias",
            self.antialias as f32,
            1.0,
            4.0,
            1.0,
        )));

        let mut controls = vec![
            PanelControl::Heading {
                id: "title".into(),
                text: "Ray Tracing".into(),
            },
            PanelControl::Text {
                id: "summary".into(),
                text: format!("{} | AA {}", self.resolution_summary(), self.antialias),
            },
            row_with_gap(
                "actions",
                vec![
                    node(button("render", "Render", true)).grow(1.0),
                    node(button("save", "Save", false)).grow(0.35),
                ],
                10.0,
            ),
            group("output", "Output", self.output_open, output_controls),
            group(
                "trace",
                "Trace",
                self.trace_open,
                vec![
                    node(PanelControl::Select {
                        id: "mode".into(),
                        label: "Mode".into(),
                        value: self.mode.to_string(),
                        options: vec![
                            PanelOption::new("Normal", "0"),
                            PanelOption::new("Outline", "1"),
                            PanelOption::new("Edges", "2"),
                            PanelOption::new("Quantized", "3"),
                        ],
                    }),
                    node(slider("fog", "Fog", self.fog, -1.0, 1.0, 0.01)),
                    node(PanelControl::Select {
                        id: "background".into(),
                        label: "Background".into(),
                        value: self.opaque_background.to_string(),
                        options: vec![
                            PanelOption::new("Auto", "-1"),
                            PanelOption::new("Transparent", "0"),
                            PanelOption::new("Opaque", "1"),
                        ],
                    }),
                ],
            ),
            group(
                "quality",
                "Quality",
                self.quality_open,
                vec![
                    node(PanelControl::Toggle {
                        id: "shadow".into(),
                        label: "Shadows".into(),
                        value: self.shadow,
                    }),
                    node(PanelControl::Toggle {
                        id: "transparency_shadows".into(),
                        label: "Transparency shadows".into(),
                        value: self.transparency_shadows,
                    }),
                    node(slider(
                        "max_passes",
                        "Max passes",
                        self.max_passes as f32,
                        1.0,
                        100.0,
                        1.0,
                    )),
                ],
            ),
        ];

        let mut lighting_controls = vec![node(PanelControl::Toggle {
            id: "use_custom".into(),
            label: "Use custom lighting".into(),
            value: self.use_custom,
        })];
        if self.use_custom {
            lighting_controls.extend([
                node(slider("ambient", "Ambient", self.ambient, 0.0, 1.0, 0.01)),
                node(slider("direct", "Direct", self.direct, 0.0, 1.0, 0.01)),
                node(slider("reflect", "Reflect", self.reflect, 0.0, 1.0, 0.01)),
                node(slider(
                    "specular",
                    "Specular",
                    self.specular,
                    0.0,
                    1.0,
                    0.01,
                )),
                node(slider(
                    "shininess",
                    "Shininess",
                    self.shininess,
                    1.0,
                    128.0,
                    1.0,
                )),
            ]);
        }
        controls.push(group(
            "lighting",
            "Lighting",
            self.lighting_open,
            lighting_controls,
        ));

        controls.push(group(
            "edge",
            "Edge Detection",
            self.edge_open,
            vec![
                node(slider("slope", "Slope", self.slope_factor, 0.0, 2.0, 0.01)),
                node(slider("depth", "Depth", self.depth_factor, 0.0, 1.0, 0.01)),
                node(slider("disco", "Disco", self.disco_factor, 0.0, 1.0, 0.01)),
                node(slider("gain", "Gain", self.gain, 0.0, 1.0, 0.01)),
            ],
        ));

        controls
    }

    fn resolution_summary(&self) -> String {
        match self.preset {
            ResolutionPreset::Custom => format!("Custom {}x{}", self.width, self.height),
            preset => format!("{} {}x{}", preset.label(), self.width, self.height),
        }
    }

    fn set_section(&mut self, id: &str, open: bool) -> bool {
        match id {
            "output" => self.output_open = open,
            "trace" => self.trace_open = open,
            "quality" => self.quality_open = open,
            "lighting" => self.lighting_open = open,
            "edge" => self.edge_open = open,
            _ => return false,
        }
        true
    }

    fn apply_event(&mut self, event: &PanelEvent) -> Vec<PanelAction> {
        if let PanelValue::Bool(open) = &event.value {
            if self.set_section(&event.control_id, *open) {
                return Vec::new();
            }
        }

        match event.control_id.as_str() {
            "preset" => {
                if let Some(value) = text_value(&event.value) {
                    self.preset = ResolutionPreset::from_value(&value);
                    if let Some((w, h)) = self.preset.dimensions() {
                        self.width = w;
                        self.height = h;
                    }
                }
                Vec::new()
            }
            "width" => {
                if let Some(v) = numeric_value(&event.value) {
                    self.width = v.round().clamp(64.0, 7680.0) as u32;
                    self.preset = ResolutionPreset::from_dimensions(self.width, self.height);
                }
                Vec::new()
            }
            "height" => {
                if let Some(v) = numeric_value(&event.value) {
                    self.height = v.round().clamp(64.0, 4320.0) as u32;
                    self.preset = ResolutionPreset::from_dimensions(self.width, self.height);
                }
                Vec::new()
            }
            "antialias" => {
                if let Some(v) = numeric_value(&event.value) {
                    self.antialias = v.round().clamp(1.0, 4.0) as u32;
                }
                Vec::new()
            }
            "render" => self.render_actions(),
            "save" => vec![self.save_file_request_action()],
            SAVE_REPLY_CONTROL_ID => {
                let Some(path) = text_value(&event.value) else {
                    return Vec::new();
                };
                let path = path.trim();
                if path.is_empty() {
                    Vec::new()
                } else {
                    vec![PanelAction::ExecuteCommand {
                        command: format!("png {}", quote_command_arg(path)),
                        silent: false,
                    }]
                }
            }
            "mode" => self.set_int_from_event(event, "ray_trace_mode", |s, v| s.mode = v),
            "background" => self.set_int_from_event(event, "ray_opaque_background", |s, v| {
                s.opaque_background = v;
            }),
            "fog" => self.set_float_from_event(event, "ray_trace_fog", |s, v| s.fog = v),
            "use_custom" => {
                self.set_bool_from_event(event, "rt_use_custom", |s, v| s.use_custom = v)
            }
            "ambient" => self.set_float_from_event(event, "rt_ambient", |s, v| s.ambient = v),
            "direct" => self.set_float_from_event(event, "rt_direct", |s, v| s.direct = v),
            "reflect" => self.set_float_from_event(event, "rt_reflect", |s, v| s.reflect = v),
            "specular" => self.set_float_from_event(event, "rt_specular", |s, v| s.specular = v),
            "shininess" => self.set_float_from_event(event, "rt_shininess", |s, v| s.shininess = v),
            "shadow" => self.set_bool_from_event(event, "ray_shadow", |s, v| s.shadow = v),
            "transparency_shadows" => {
                self.set_bool_from_event(event, "ray_transparency_shadows", |s, v| {
                    s.transparency_shadows = v
                })
            }
            "max_passes" => {
                self.set_int_from_event(event, "ray_max_passes", |s, v| s.max_passes = v)
            }
            "slope" => self.set_float_from_event(event, "ray_trace_slope_factor", |s, v| {
                s.slope_factor = v;
            }),
            "depth" => self.set_float_from_event(event, "ray_trace_depth_factor", |s, v| {
                s.depth_factor = v;
            }),
            "disco" => self.set_float_from_event(event, "ray_trace_disco_factor", |s, v| {
                s.disco_factor = v;
            }),
            "gain" => self.set_float_from_event(event, "ray_trace_gain", |s, v| s.gain = v),
            _ => Vec::new(),
        }
    }

    fn set_bool_from_event(
        &mut self,
        event: &PanelEvent,
        setting: &str,
        set: impl FnOnce(&mut Self, bool),
    ) -> Vec<PanelAction> {
        let Some(value) = bool_value(&event.value) else {
            return Vec::new();
        };
        set(self, value);
        vec![set_setting(setting, PanelValue::Bool(value))]
    }

    fn set_int_from_event(
        &mut self,
        event: &PanelEvent,
        setting: &str,
        set: impl FnOnce(&mut Self, i32),
    ) -> Vec<PanelAction> {
        let Some(value) = numeric_value(&event.value) else {
            return Vec::new();
        };
        let value = value.round() as i32;
        set(self, value);
        vec![set_setting(setting, PanelValue::Text(value.to_string()))]
    }

    fn set_float_from_event(
        &mut self,
        event: &PanelEvent,
        setting: &str,
        set: impl FnOnce(&mut Self, f32),
    ) -> Vec<PanelAction> {
        let Some(value) = numeric_value(&event.value) else {
            return Vec::new();
        };
        set(self, value);
        vec![set_setting(setting, PanelValue::Number(value))]
    }

    fn render_actions(&self) -> Vec<PanelAction> {
        let mut actions = vec![
            set_cmd("ray_shadow", self.shadow as i32),
            set_cmd("ray_transparency_shadows", self.transparency_shadows as i32),
            set_cmd("ray_max_passes", self.max_passes),
            set_cmd("ray_trace_mode", self.mode),
            set_cmd("ray_trace_fog", self.fog),
            set_cmd("ray_trace_slope_factor", self.slope_factor),
            set_cmd("ray_trace_depth_factor", self.depth_factor),
            set_cmd("ray_trace_disco_factor", self.disco_factor),
            set_cmd("ray_trace_gain", self.gain),
            set_cmd("ray_opaque_background", self.opaque_background),
            set_cmd("rt_use_custom", self.use_custom as i32),
        ];
        if self.use_custom {
            actions.extend([
                set_cmd("rt_ambient", self.ambient),
                set_cmd("rt_direct", self.direct),
                set_cmd("rt_reflect", self.reflect),
                set_cmd("rt_specular", self.specular),
                set_cmd("rt_shininess", self.shininess),
            ]);
        }

        actions.push(PanelAction::ExecuteCommand {
            command: format!("ray {}, {}, {}", self.width, self.height, self.antialias),
            silent: false,
        });
        actions
    }

    fn save_file_request_action(&self) -> PanelAction {
        custom_action(
            SAVE_FILE_REQUEST_TOPIC,
            &SaveFileRequest {
                panel_id: PANEL_ID.into(),
                reply_control_id: SAVE_REPLY_CONTROL_ID.into(),
                title: "Save ray-traced image".into(),
                default_file_name: DEFAULT_SAVE_FILE.into(),
                allowed_extensions: vec!["png".into()],
            },
        )
    }
}

impl PluginPanel for RtPanel {
    fn descriptor(&self) -> PanelDescriptor {
        PanelDescriptor::right(PANEL_ID, "Ray Tracing")
            .icon("RT")
            .default_visible(false)
    }

    fn snapshot(&mut self, ctx: &SharedContext<'_>) -> PanelSnapshot {
        self.sync_from_settings(ctx);
        PanelSnapshot::new(self.controls())
    }

    fn handle_event(
        &mut self,
        event: PanelEvent,
        _ctx: &SharedContext<'_>,
        _bus: &mut MessageBus,
    ) -> Vec<PanelAction> {
        self.apply_event(&event)
    }
}

impl ResolutionPreset {
    const ALL: [Self; 4] = [Self::Hd720, Self::Hd1080, Self::Uhd4K, Self::Custom];
}

fn button(id: &str, label: &str, primary: bool) -> PanelControl {
    PanelControl::Button {
        id: id.into(),
        label: label.into(),
        primary,
    }
}

fn row_with_gap(id: &str, children: Vec<PanelControlNode>, gap: f32) -> PanelControl {
    PanelControl::Row(PanelRow::new(id, children).gap(gap))
}

fn group(id: &str, title: &str, open: bool, children: Vec<PanelControlNode>) -> PanelControl {
    PanelControl::Group(PanelGroup::new(id, title, children).open(open))
}

fn node(control: PanelControl) -> PanelControlNode {
    PanelControlNode::new(control)
}

fn number(id: &str, label: &str, value: f32, min: f32, max: f32, step: f32) -> PanelControl {
    PanelControl::Number {
        id: id.into(),
        label: label.into(),
        value,
        min,
        max,
        step,
    }
}

fn slider(id: &str, label: &str, value: f32, min: f32, max: f32, step: f32) -> PanelControl {
    PanelControl::Slider {
        id: id.into(),
        label: label.into(),
        value,
        min,
        max,
        step,
    }
}

fn set_setting(name: &str, value: PanelValue) -> PanelAction {
    PanelAction::SetSetting {
        name: name.into(),
        value,
    }
}

fn set_cmd(name: &str, value: impl std::fmt::Display) -> PanelAction {
    PanelAction::ExecuteCommand {
        command: format!("set {name}, {value}"),
        silent: true,
    }
}

fn text_value(value: &PanelValue) -> Option<String> {
    match value {
        PanelValue::Text(v) => Some(v.clone()),
        PanelValue::Number(v) => Some(v.to_string()),
        PanelValue::Bool(v) => Some((*v as i32).to_string()),
        PanelValue::None => None,
    }
}

fn numeric_value(value: &PanelValue) -> Option<f32> {
    match value {
        PanelValue::Number(v) => Some(*v),
        PanelValue::Text(v) => v.trim().parse().ok(),
        PanelValue::Bool(v) => Some(if *v { 1.0 } else { 0.0 }),
        PanelValue::None => None,
    }
}

fn bool_value(value: &PanelValue) -> Option<bool> {
    match value {
        PanelValue::Bool(v) => Some(*v),
        PanelValue::Number(v) => Some(*v != 0.0),
        PanelValue::Text(v) => match v.trim() {
            "1" | "true" | "on" => Some(true),
            "0" | "false" | "off" => Some(false),
            _ => None,
        },
        PanelValue::None => None,
    }
}

fn quote_command_arg(s: &str) -> String {
    if s.chars()
        .any(|c| c.is_whitespace() || matches!(c, ',' | '"' | '\''))
    {
        format!("\"{}\"", s.replace('\\', "\\\\").replace('"', "\\\""))
    } else {
        s.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn event(control_id: &str, value: PanelValue) -> PanelEvent {
        PanelEvent {
            panel_id: PANEL_ID.into(),
            control_id: control_id.into(),
            kind: PanelEventKind::TextCommit,
            value,
        }
    }

    fn control_ids(controls: &[PanelControl]) -> Vec<String> {
        let mut ids = Vec::new();
        for control in controls {
            collect_control_ids(control, &mut ids);
        }
        ids
    }

    fn collect_control_ids(control: &PanelControl, ids: &mut Vec<String>) {
        match control {
            PanelControl::Text { id, .. }
            | PanelControl::Heading { id, .. }
            | PanelControl::Section { id, .. }
            | PanelControl::Button { id, .. }
            | PanelControl::ButtonRow { id, .. }
            | PanelControl::Toggle { id, .. }
            | PanelControl::Slider { id, .. }
            | PanelControl::Number { id, .. }
            | PanelControl::Select { id, .. }
            | PanelControl::TextInput { id, .. }
            | PanelControl::Image { id, .. }
            | PanelControl::Spacer { id, .. } => ids.push(id.clone()),
            PanelControl::TextArea(area) => ids.push(area.id.clone()),
            PanelControl::Row(row) => {
                ids.push(row.id.clone());
                for child in &row.children {
                    collect_control_ids(&child.control, ids);
                }
            }
            PanelControl::Column(column) => {
                ids.push(column.id.clone());
                for child in &column.children {
                    collect_control_ids(&child.control, ids);
                }
            }
            PanelControl::Group(group) => {
                ids.push(group.id.clone());
                for child in &group.children {
                    collect_control_ids(&child.control, ids);
                }
            }
        }
    }

    fn execute_commands(actions: &[PanelAction]) -> Vec<String> {
        actions
            .iter()
            .filter_map(|action| match action {
                PanelAction::ExecuteCommand { command, .. } => Some(command.clone()),
                PanelAction::SetSetting { .. } | PanelAction::Custom { .. } => None,
            })
            .collect()
    }

    fn group_open(controls: &[PanelControl], id: &str) -> Option<bool> {
        controls.iter().find_map(|control| match control {
            PanelControl::Group(group) if group.id == id => Some(group.open),
            _ => None,
        })
    }

    fn find_control<'a>(controls: &'a [PanelControl], id: &str) -> Option<&'a PanelControl> {
        controls
            .iter()
            .find_map(|control| find_control_in(control, id))
    }

    fn find_control_in<'a>(control: &'a PanelControl, id: &str) -> Option<&'a PanelControl> {
        if control_id(control).is_some_and(|control_id| control_id == id) {
            return Some(control);
        }

        match control {
            PanelControl::Row(row) => row
                .children
                .iter()
                .find_map(|child| find_control_in(&child.control, id)),
            PanelControl::Column(column) => column
                .children
                .iter()
                .find_map(|child| find_control_in(&child.control, id)),
            PanelControl::Group(group) => group
                .children
                .iter()
                .find_map(|child| find_control_in(&child.control, id)),
            _ => None,
        }
    }

    fn control_id(control: &PanelControl) -> Option<&str> {
        match control {
            PanelControl::Text { id, .. }
            | PanelControl::Heading { id, .. }
            | PanelControl::Section { id, .. }
            | PanelControl::Button { id, .. }
            | PanelControl::ButtonRow { id, .. }
            | PanelControl::Toggle { id, .. }
            | PanelControl::Slider { id, .. }
            | PanelControl::Number { id, .. }
            | PanelControl::Select { id, .. }
            | PanelControl::TextInput { id, .. }
            | PanelControl::Image { id, .. }
            | PanelControl::Spacer { id, .. }
            | PanelControl::Row(PanelRow { id, .. })
            | PanelControl::Column(PanelColumn { id, .. })
            | PanelControl::Group(PanelGroup { id, .. }) => Some(id),
            PanelControl::TextArea(area) => Some(area.id.as_str()),
        }
    }

    #[test]
    fn snapshot_omits_save_path() {
        let mut panel = RtPanel::default();

        let ids = control_ids(&panel.controls());

        assert!(!ids.iter().any(|id| id == "save_path"));
    }

    #[test]
    fn custom_dimensions_are_only_visible_for_custom_preset() {
        let mut panel = RtPanel::default();

        let ids = control_ids(&panel.controls());
        assert!(!ids.iter().any(|id| id == "width"));
        assert!(!ids.iter().any(|id| id == "height"));

        panel.apply_event(&event("preset", PanelValue::Text("custom".into())));
        let ids = control_ids(&panel.controls());

        assert!(ids.iter().any(|id| id == "width"));
        assert!(ids.iter().any(|id| id == "height"));
    }

    #[test]
    fn custom_dimensions_render_as_one_row() {
        let mut panel = RtPanel::default();
        panel.apply_event(&event("preset", PanelValue::Text("custom".into())));
        let controls = panel.controls();

        let Some(PanelControl::Row(row)) = find_control(&controls, "custom_dimensions") else {
            panic!("expected custom dimensions row");
        };

        assert_eq!(row.children.len(), 2);
        assert_eq!(control_id(&row.children[0].control), Some("width"));
        assert_eq!(control_id(&row.children[1].control), Some("height"));
    }

    #[test]
    fn integer_slider_controls_use_integral_steps() {
        let mut panel = RtPanel::default();
        let controls = panel.controls();

        for id in ["antialias", "max_passes"] {
            let Some(PanelControl::Slider { step, .. }) = find_control(&controls, id) else {
                panic!("expected {id} slider");
            };

            assert_eq!(*step, 1.0);
        }
    }

    #[test]
    fn fixed_preset_updates_dimensions() {
        let mut panel = RtPanel::default();
        panel.apply_event(&event("preset", PanelValue::Text("720".into())));

        assert_eq!(panel.width, 1280);
        assert_eq!(panel.height, 720);
    }

    #[test]
    fn panel_groups_are_collapsible() {
        let mut panel = RtPanel::default();

        for id in ["output", "trace", "quality", "lighting", "edge"] {
            panel.apply_event(&event(id, PanelValue::Bool(false)));
            assert_eq!(group_open(&panel.controls(), id), Some(false));

            panel.apply_event(&event(id, PanelValue::Bool(true)));
            assert_eq!(group_open(&panel.controls(), id), Some(true));
        }
    }

    #[test]
    fn save_emits_typed_save_file_request() {
        let mut panel = RtPanel::default();

        let actions = panel.apply_event(&event("save", PanelValue::None));
        let [PanelAction::Custom { topic, payload }] = actions.as_slice() else {
            panic!("expected one custom action");
        };
        let msg = AppMessage::Custom {
            topic: topic.clone(),
            payload: payload.clone(),
        };
        let request: Option<SaveFileRequest> = subscribe(&msg, SAVE_FILE_REQUEST_TOPIC);

        assert_eq!(
            request,
            Some(SaveFileRequest {
                panel_id: PANEL_ID.into(),
                reply_control_id: SAVE_REPLY_CONTROL_ID.into(),
                title: "Save ray-traced image".into(),
                default_file_name: DEFAULT_SAVE_FILE.into(),
                allowed_extensions: vec!["png".into()],
            })
        );
    }

    #[test]
    fn selected_save_file_saves_current_viewport_image() {
        let mut panel = RtPanel::default();

        let actions = panel.apply_event(&event(
            SAVE_REPLY_CONTROL_ID,
            PanelValue::Text("/tmp/ray trace.png".into()),
        ));
        let commands = execute_commands(&actions);

        assert!(commands
            .iter()
            .any(|command| command == "png \"/tmp/ray trace.png\""));
        assert!(!commands.iter().any(|command| command.starts_with("ray ")));
    }
}
