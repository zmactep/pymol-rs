//! Session — pure scene state
//!
//! [`Session`] holds all domain/scene data that can exist independently of
//! rendering: object registry, camera, selections, settings, colors, etc.
//!
//! This is the single source of truth that both the GUI (`pymol-gui`) and
//! any headless adapter can own. GPU resources live elsewhere (e.g., in
//! `pymol-render::RenderContext`).

use pymol_color::{ChainColors, ElementColors, NamedColors};
use pymol_render::{ColorResolver, RenderContext};
use pymol_settings::GlobalSettings;
use serde::{Deserialize, Serialize};

use crate::camera::Camera;
use crate::movie::Movie;
use crate::object::{ObjectRegistry, ObjectRegistrySnapshot};
use crate::raytrace::RaytraceInput;
use crate::scene::SceneManager;
use crate::selection::SelectionManager;
use crate::view::ViewManager;
use crate::viewer_trait::RaytracedImage;

/// Pure scene state — no GPU resources, no window, no event loop.
///
/// Owns all molecular objects, camera state, named selections, scenes,
/// views, animation, settings, and color tables.
///
/// Implements `Serialize` and `Deserialize` via a proxy that converts
/// the [`ObjectRegistry`] to/from an [`ObjectRegistrySnapshot`].
pub struct Session {
    // =========================================================================
    // Scene
    // =========================================================================
    /// Object registry (molecules, surfaces, maps, CGO, etc.)
    pub registry: ObjectRegistry,
    /// Camera for view control
    pub camera: Camera,
    /// Named selections manager
    pub selections: SelectionManager,
    /// Scene manager for named snapshots (camera + object state)
    pub scenes: SceneManager,
    /// Named views (camera state only — simpler than scenes)
    pub views: ViewManager,
    /// Movie player for frame-based animation
    pub movie: Movie,

    // =========================================================================
    // Settings and Colors
    // =========================================================================
    /// Global rendering settings
    pub settings: GlobalSettings,
    /// Named colors table (e.g., "red", "carbon")
    pub named_colors: NamedColors,
    /// Per-element color defaults
    pub element_colors: ElementColors,
    /// Chain-based coloring (unit struct)
    pub chain_colors: ChainColors,

    // =========================================================================
    // Visual Properties
    // =========================================================================
    /// Background (clear) color as linear RGB floats
    pub clear_color: [f32; 3],

    // =========================================================================
    // Raytraced Image Overlay
    // =========================================================================
    /// Stored raytraced image for display (from `ray` command without filename)
    pub raytraced_image: Option<RaytracedImage>,
}

/// Serializable proxy for [`Session`] (for deserialization).
#[derive(Deserialize)]
struct SessionProxy {
    registry: ObjectRegistrySnapshot,
    camera: Camera,
    selections: SelectionManager,
    scenes: SceneManager,
    views: ViewManager,
    movie: Movie,
    settings: GlobalSettings,
    named_colors: NamedColors,
    element_colors: ElementColors,
    chain_colors: ChainColors,
    clear_color: [f32; 3],
}

impl Serialize for Session {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        use serde::ser::SerializeStruct;
        let mut s = serializer.serialize_struct("Session", 11)?;
        s.serialize_field("registry", &self.registry.to_snapshot())?;
        s.serialize_field("camera", &self.camera)?;
        s.serialize_field("selections", &self.selections)?;
        s.serialize_field("scenes", &self.scenes)?;
        s.serialize_field("views", &self.views)?;
        s.serialize_field("movie", &self.movie)?;
        s.serialize_field("settings", &self.settings)?;
        s.serialize_field("named_colors", &self.named_colors)?;
        s.serialize_field("element_colors", &self.element_colors)?;
        s.serialize_field("chain_colors", &self.chain_colors)?;
        s.serialize_field("clear_color", &self.clear_color)?;
        s.end()
    }
}

impl<'de> Deserialize<'de> for Session {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let proxy = SessionProxy::deserialize(deserializer)?;
        Ok(Session {
            registry: ObjectRegistry::from_snapshot(proxy.registry),
            camera: proxy.camera,
            selections: proxy.selections,
            scenes: proxy.scenes,
            views: proxy.views,
            movie: proxy.movie,
            settings: proxy.settings,
            named_colors: proxy.named_colors,
            element_colors: proxy.element_colors,
            chain_colors: proxy.chain_colors,
            clear_color: proxy.clear_color,
            raytraced_image: None,
        })
    }
}

impl Default for Session {
    fn default() -> Self {
        Self::new()
    }
}

impl Session {
    /// Ensure all molecule representations are built (needed before raytrace
    /// when the render loop hasn't run yet, e.g. in headless/script mode).
    pub fn prepare_render_all(&mut self, context: &RenderContext) {
        let names: Vec<String> = self.registry.names().map(|s| s.to_string()).collect();
        for name in &names {
            let color_resolver = ColorResolver::new(
                &self.named_colors,
                &self.element_colors,
                &self.chain_colors,
            );
            if let Some(mol_obj) = self.registry.get_molecule_mut(name) {
                mol_obj.prepare_render(context, &color_resolver, &self.settings);
            }
        }
    }

    /// Build a [`RaytraceInput`] borrowing from this session.
    pub fn raytrace_input<'a>(
        &'a mut self,
        context: &'a RenderContext,
        default_size: (u32, u32),
    ) -> RaytraceInput<'a> {
        RaytraceInput {
            device: context.device(),
            queue: context.queue(),
            camera: &mut self.camera,
            registry: &self.registry,
            settings: &self.settings,
            named_colors: &self.named_colors,
            element_colors: &self.element_colors,
            chain_colors: &self.chain_colors,
            clear_color: self.clear_color,
            default_size,
        }
    }

    /// Create a new session with default values.
    pub fn new() -> Self {
        Self {
            registry: ObjectRegistry::new(),
            camera: Camera::new(),
            selections: SelectionManager::new(),
            scenes: SceneManager::new(),
            views: ViewManager::new(),
            movie: Movie::new(),
            settings: GlobalSettings::new(),
            named_colors: NamedColors::default(),
            element_colors: ElementColors::default(),
            chain_colors: ChainColors,
            clear_color: [0.0, 0.0, 0.0],
            raytraced_image: None,
        }
    }
}
