//! Ray Tracing Toolbar Component
//!
//! An egui panel that provides interactive controls for ray tracing parameters.
//! Supports both system settings and custom overrides for Phong lighting.

mod preview;
mod sections;

use pymol_plugin::prelude::{Component, EguiComponent, MessageBus, SharedContext};
use pymol_settings::SettingValue;

// ---------------------------------------------------------------------------
// Resolution presets
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ResolutionPreset {
    Hd720,
    Hd1080,
    Uhd4K,
    Custom,
}

impl ResolutionPreset {
    pub fn label(self) -> &'static str {
        match self {
            Self::Hd720 => "720p",
            Self::Hd1080 => "1080p",
            Self::Uhd4K => "4K",
            Self::Custom => "Custom",
        }
    }

    pub fn dimensions(self) -> Option<(u32, u32)> {
        match self {
            Self::Hd720 => Some((1280, 720)),
            Self::Hd1080 => Some((1920, 1080)),
            Self::Uhd4K => Some((3840, 2160)),
            Self::Custom => None,
        }
    }

    pub fn from_dimensions(w: u32, h: u32) -> Self {
        match (w, h) {
            (1280, 720) => Self::Hd720,
            (1920, 1080) => Self::Hd1080,
            (3840, 2160) => Self::Uhd4K,
            _ => Self::Custom,
        }
    }

    pub const ALL: [Self; 4] = [Self::Hd720, Self::Hd1080, Self::Uhd4K, Self::Custom];
}

// ---------------------------------------------------------------------------
// Label helpers
// ---------------------------------------------------------------------------

pub(crate) fn trace_mode_label(mode: i32) -> &'static str {
    match mode {
        0 => "Normal",
        1 => "Normal + Outline",
        2 => "Outline Only",
        3 => "Quantized + Outline",
        _ => "Unknown",
    }
}

pub(crate) fn background_label(val: i32) -> &'static str {
    match val {
        -1 => "Auto",
        0 => "Transparent",
        1 => "Opaque",
        _ => "Auto",
    }
}

// ---------------------------------------------------------------------------
// Component state
// ---------------------------------------------------------------------------

pub struct RtToolbarComponent {
    // Resolution
    pub(crate) width: u32,
    pub(crate) height: u32,
    pub(crate) preset: ResolutionPreset,
    pub(crate) antialias: u32,

    // Use custom Phong values
    pub(crate) use_custom: bool,

    // Phong lighting
    pub(crate) ambient: f32,
    pub(crate) direct: f32,
    pub(crate) reflect: f32,
    pub(crate) specular: f32,
    pub(crate) shininess: f32,

    // Shadows
    pub(crate) shadow: bool,
    pub(crate) transparency_shadows: bool,
    pub(crate) max_passes: i32,

    // Trace mode
    pub(crate) mode: i32,

    // Fog
    pub(crate) fog: f32,

    // Edge detection
    pub(crate) slope_factor: f32,
    pub(crate) depth_factor: f32,
    pub(crate) disco_factor: f32,
    pub(crate) gain: f32,

    // Background
    pub(crate) opaque_background: i32,

    // Preview
    pub(crate) preview_texture: Option<egui::TextureHandle>,
    pub(crate) preview_dirty: bool,

    // Scene change detection (debounced)
    pub(crate) last_scene_gen: u64,
    pub(crate) stable_frames: u32,
}

impl Default for RtToolbarComponent {
    fn default() -> Self {
        Self {
            width: 1920,
            height: 1080,
            preset: ResolutionPreset::Hd1080,
            antialias: 2,
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
            preview_texture: None,
            preview_dirty: true,
            last_scene_gen: 0,
            stable_frames: 0,
        }
    }
}

impl RtToolbarComponent {
    pub fn new() -> Self {
        Self::default()
    }

    /// Sync component fields from the dynamic settings store.
    fn sync_from_settings(&mut self, ctx: &SharedContext) {
        let dyn_reg = match ctx.dynamic_settings {
            Some(r) => r,
            None => return,
        };

        let read = |name: &str| -> Option<SettingValue> {
            let entry = dyn_reg.lookup(name)?;
            let store = entry.store.read().ok()?;
            store.get(name).cloned()
        };

        macro_rules! sync {
            ($target:expr, $name:expr, Bool) => {
                if let Some(SettingValue::Bool(v)) = read($name) { $target = v; }
            };
            ($target:expr, $name:expr, Float) => {
                if let Some(SettingValue::Float(v)) = read($name) { $target = v; }
            };
            ($target:expr, $name:expr, Int) => {
                if let Some(SettingValue::Int(v)) = read($name) { $target = v; }
            };
        }

        sync!(self.use_custom, "rt_use_custom", Bool);
        sync!(self.ambient, "rt_ambient", Float);
        sync!(self.direct, "rt_direct", Float);
        sync!(self.reflect, "rt_reflect", Float);
        sync!(self.specular, "rt_specular", Float);
        sync!(self.shininess, "rt_shininess", Float);
        sync!(self.shadow, "ray_shadow", Bool);
        sync!(self.transparency_shadows, "ray_transparency_shadows", Bool);
        sync!(self.max_passes, "ray_max_passes", Int);
        sync!(self.mode, "ray_trace_mode", Int);
        sync!(self.fog, "ray_trace_fog", Float);
        sync!(self.slope_factor, "ray_trace_slope_factor", Float);
        sync!(self.depth_factor, "ray_trace_depth_factor", Float);
        sync!(self.disco_factor, "ray_trace_disco_factor", Float);
        sync!(self.gain, "ray_trace_gain", Float);
        sync!(self.opaque_background, "ray_opaque_background", Int);
    }

    /// Send `set` commands for all ray settings, then execute the ray command.
    fn send_render_commands(&self, bus: &mut MessageBus, save_path: Option<&str>) {
        bus.execute_command_silent(format!("set ray_shadow, {}", self.shadow as i32));
        bus.execute_command_silent(format!(
            "set ray_transparency_shadows, {}",
            self.transparency_shadows as i32
        ));
        bus.execute_command_silent(format!("set ray_max_passes, {}", self.max_passes));
        bus.execute_command_silent(format!("set ray_trace_mode, {}", self.mode));
        bus.execute_command_silent(format!("set ray_trace_fog, {}", self.fog));
        bus.execute_command_silent(format!(
            "set ray_trace_slope_factor, {}",
            self.slope_factor
        ));
        bus.execute_command_silent(format!(
            "set ray_trace_depth_factor, {}",
            self.depth_factor
        ));
        bus.execute_command_silent(format!(
            "set ray_trace_disco_factor, {}",
            self.disco_factor
        ));
        bus.execute_command_silent(format!("set ray_trace_gain, {}", self.gain));
        bus.execute_command_silent(format!(
            "set ray_opaque_background, {}",
            self.opaque_background
        ));

        bus.execute_command_silent(format!(
            "set rt_use_custom, {}",
            self.use_custom as i32
        ));
        if self.use_custom {
            bus.execute_command_silent(format!("set rt_ambient, {}", self.ambient));
            bus.execute_command_silent(format!("set rt_direct, {}", self.direct));
            bus.execute_command_silent(format!("set rt_reflect, {}", self.reflect));
            bus.execute_command_silent(format!("set rt_specular, {}", self.specular));
            bus.execute_command_silent(format!("set rt_shininess, {}", self.shininess));
        }

        let cmd = if let Some(path) = save_path {
            format!(
                "ray {}, {}, {}, {}",
                self.width, self.height, self.antialias, path
            )
        } else {
            format!("ray {}, {}, {}", self.width, self.height, self.antialias)
        };
        bus.execute_command(cmd);
    }
}

impl Component for RtToolbarComponent {
    fn id(&self) -> &'static str {
        "rt_toolbar"
    }

    fn title(&self) -> &str {
        "Ray Tracing"
    }
}

impl EguiComponent for RtToolbarComponent {
    fn show(&mut self, ui: &mut egui::Ui, ctx: &SharedContext, bus: &mut MessageBus) {
        // Disable selectable labels for the entire plugin UI.
        ui.style_mut().interaction.selectable_labels = false;

        self.sync_from_settings(ctx);

        // Debounced scene change detection for preview updates
        preview::update_scene_debounce(self, ui, ctx);

        if self.preview_dirty {
            preview::render_preview(self, ui.ctx(), ctx);
        }

        let mut setting_changed = false;

        egui::ScrollArea::vertical().show(ui, |ui| {
            ui.spacing_mut().item_spacing.y = 4.0;

            sections::show_preview(self, ui);

            setting_changed |= sections::show_custom_lighting_toggle(self, ui, bus);
            ui.separator();

            sections::show_resolution(self, ui);
            sections::show_antialias(self, ui);
            setting_changed |= sections::show_mode(self, ui, bus);
            ui.separator();

            setting_changed |= sections::show_lighting(self, ui, bus);
            ui.separator();

            setting_changed |= sections::show_shadows(self, ui, bus);
            setting_changed |= sections::show_edge_detection(self, ui, bus);
            setting_changed |= sections::show_other(self, ui, bus);
            ui.separator();

            sections::show_action_buttons(self, ui, bus);
        });

        if setting_changed {
            self.preview_dirty = true;
        }
    }
}
