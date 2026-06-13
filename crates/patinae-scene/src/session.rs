//! Session — pure scene state
//!
//! [`Session`] holds all domain/scene data that can exist independently of
//! rendering: object registry, camera, selections, settings, colors, etc.
//!
//! This is the single source of truth that GUI and headless adapters can own.
//! GPU resources live behind a host-provided render target.

use patinae_color::{NamedPalette, ThemedPalette};
use patinae_select::SelectionResult;
use patinae_settings::Settings;
use serde::{Deserialize, Serialize};

use crate::camera::Camera;
use crate::highlight_state::HighlightState;
use crate::movie::Movie;
use crate::object::{DirtyFlags, Object, ObjectRegistry, ObjectRegistrySnapshot};
use crate::scene::SceneManager;
use crate::selection::SelectionManager;
use crate::view::ViewManager;
use crate::viewer_trait::ViewportImage;

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
    pub settings: Settings,
    /// Named colors table (e.g., "red", "carbon")
    pub named_palette: NamedPalette,
    /// Theme-aware palette (element, chain, SS, residue, gradients, etc.)
    pub palette: ThemedPalette,

    // =========================================================================
    // Visual Properties
    // =========================================================================
    /// Background (clear) color as linear RGB floats
    pub clear_color: [f32; 3],
    /// Whether clear_color has been explicitly set by the user
    pub clear_color_set: bool,

    // =========================================================================
    // Viewport Image Overlay
    // =========================================================================
    /// Image overlay for display in the viewport (e.g. from `ray` command or plugins)
    pub viewport_image: Option<ViewportImage>,

    // =========================================================================
    // Highlight state (transient — not serialized)
    // =========================================================================
    /// GPU selection / hover bitmask state for the screen-space highlight pass.
    /// Rebuilt every frame in `prepare_scene` from `selections.evaluate_visible`
    /// and `hover_target`.
    pub highlight_state: HighlightState,
    /// Active hover target. `None` when nothing is hovered.
    pub hover_target: Option<HoverTarget>,
}

/// Atoms currently under the cursor, fed into the screen-space highlight pass.
///
/// The bridge layer (patinae / web) sets this on cursor-move; the highlight
/// state reads `(object, selection)` to set hover bits in its bitmap.
#[derive(Debug, Clone)]
pub struct HoverTarget {
    /// Object whose coords resolve `selection`'s indices.
    pub object: String,
    /// Atom indices to mark.
    pub selection: SelectionResult,
}

/// Result of advancing session-owned animations for one host frame.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct AnimationUpdate {
    /// Movie playback advanced to a different frame.
    pub movie_frame_changed: bool,
    /// The current movie frame was applied to scene objects/camera.
    pub movie_synced: bool,
    /// Rock animation changed the camera.
    pub rock_changed: bool,
    /// Camera interpolation changed the camera.
    pub camera_changed: bool,
    /// Any visible state changed and the host should render.
    pub needs_redraw: bool,
}

/// Lightweight movie state for host UI/API layers.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
pub struct MovieStateSnapshot {
    /// Effective frame count: explicit movie frames or max object states.
    pub frame_count: usize,
    /// Current 0-based movie frame.
    pub current_frame: usize,
    /// Whether playback is active.
    pub is_playing: bool,
    /// Whether rock animation is active.
    pub rock_enabled: bool,
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
    settings: Settings,
    #[serde(alias = "named_colors")]
    named_palette: NamedPalette,
    #[serde(alias = "element_palette", alias = "element_colors")]
    palette: ThemedPalette,
    clear_color: [f32; 3],
    #[serde(default)]
    clear_color_set: bool,
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
        s.serialize_field("named_palette", &self.named_palette)?;
        s.serialize_field("palette", &self.palette)?;
        s.serialize_field("clear_color", &self.clear_color)?;
        s.serialize_field("clear_color_set", &self.clear_color_set)?;
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
            named_palette: proxy.named_palette,
            palette: proxy.palette,
            clear_color: proxy.clear_color,
            clear_color_set: proxy.clear_color_set,
            viewport_image: None,
            highlight_state: HighlightState::new(),
            hover_target: None,
        })
    }
}

impl Default for Session {
    fn default() -> Self {
        Self::new()
    }
}

impl Session {
    /// Create a new session with default values.
    pub fn new() -> Self {
        Self {
            registry: ObjectRegistry::new(),
            camera: Camera::new(),
            selections: SelectionManager::new(),
            scenes: SceneManager::new(),
            views: ViewManager::new(),
            movie: Movie::new(),
            settings: Settings::default(),
            named_palette: NamedPalette::default(),
            palette: ThemedPalette::dark(),
            clear_color: [0.0, 0.0, 0.0],
            clear_color_set: false,
            viewport_image: None,
            highlight_state: HighlightState::new(),
            hover_target: None,
        }
    }

    /// Set the active hover target. Renders next frame.
    pub fn set_hover(&mut self, target: HoverTarget) {
        self.hover_target = Some(target);
    }

    /// Clear the active hover target.
    pub fn clear_hover(&mut self) {
        self.hover_target = None;
    }

    /// Recompute the effective trajectory/movie frame count from loaded objects.
    pub fn refresh_movie_state_count(&mut self) {
        let max_states = self
            .registry
            .iter()
            .map(|obj| obj.n_states())
            .max()
            .unwrap_or(1);
        self.movie.set_n_object_states(max_states);
    }

    /// Apply the current movie frame to scenes, object states, transforms, and camera.
    ///
    /// Frame state is applied after scene recall, and movie camera view is applied
    /// last so explicit movie keyframes win over scene camera snapshots.
    pub fn sync_movie_frame(&mut self) -> bool {
        let current_frame = self.movie.current_frame();
        let state_index = self.movie.frame_to_state(current_frame);
        let scene_name = self.movie.current_scene_name().map(ToOwned::to_owned);
        let view = self.movie.interpolated_view();
        let (object_states, object_keyframe_states) = self
            .movie
            .frame()
            .map(|frame| {
                let object_states = frame
                    .object_states
                    .iter()
                    .map(|(name, state)| (name.clone(), *state))
                    .collect::<Vec<_>>();
                let object_keyframe_states = frame
                    .object_keyframes
                    .iter()
                    .filter_map(|(name, keyframe)| {
                        keyframe.state.map(|state| (name.clone(), state))
                    })
                    .collect::<Vec<_>>();
                (object_states, object_keyframe_states)
            })
            .unwrap_or_default();

        let mut changed = false;

        if let Some(scene_name) = scene_name {
            if let Some(scene) = self.scenes.get(&scene_name) {
                scene.apply(&mut self.camera, &mut self.registry, false, 0.0);
                changed = true;
            }
        }

        let names: Vec<String> = self.registry.names().map(ToOwned::to_owned).collect();
        for name in &names {
            if let Some(obj) = self.registry.get_molecule_mut(name) {
                let before = obj.display_state();
                if obj.set_display_state(state_index) && obj.display_state() != before {
                    changed = true;
                }
            }
        }

        for (name, state) in object_states.into_iter().chain(object_keyframe_states) {
            if let Some(obj) = self.registry.get_mut(&name) {
                let before = obj.current_state();
                let state = state.saturating_sub(1);
                if obj.set_current_state(state) && obj.current_state() != before {
                    changed = true;
                }
            }
        }

        changed |= self.apply_movie_object_transforms();

        if let Some(view) = view {
            self.camera.set_view(view);
            changed = true;
        }

        changed
    }

    /// Advance movie playback, rock animation, and camera interpolation.
    pub fn update_animations(&mut self, dt: f32) -> AnimationUpdate {
        let dt = if dt.is_finite() {
            dt.clamp(0.0, 0.25)
        } else {
            0.0
        };

        self.movie.set_fps(self.settings.movie.movie_fps);

        let was_playing = self.movie.is_playing();
        let movie_frame_changed = self.movie.update(dt);
        let movie_synced = if movie_frame_changed {
            self.sync_movie_frame()
        } else {
            false
        };
        let playback_state_changed = was_playing != self.movie.is_playing();

        let rock_delta = if self.movie.is_rock_enabled() {
            let amplitude = 45.0_f32.to_radians();
            let speed = 5.0;
            self.movie.update_rock(dt, amplitude, speed)
        } else {
            0.0
        };
        let rock_changed = rock_delta != 0.0;
        if rock_changed {
            self.camera.rotate_y(rock_delta);
        }

        let camera_changed = self.camera.update(dt);

        AnimationUpdate {
            movie_frame_changed,
            movie_synced,
            rock_changed,
            camera_changed,
            needs_redraw: movie_frame_changed
                || movie_synced
                || playback_state_changed
                || rock_changed
                || camera_changed,
        }
    }

    /// Get a host-facing snapshot of the current movie state.
    pub fn movie_state_snapshot(&self) -> MovieStateSnapshot {
        MovieStateSnapshot {
            frame_count: self.movie.effective_frame_count(),
            current_frame: self.movie.current_frame(),
            is_playing: self.movie.is_playing(),
            rock_enabled: self.movie.is_rock_enabled(),
        }
    }

    /// Apply interpolated object transforms from the movie's current frame to the registry.
    pub fn apply_movie_object_transforms(&mut self) -> bool {
        let mut changed = false;
        for (name, transform) in self.movie.objects_with_transforms() {
            if let Some(obj) = self.registry.get_molecule_mut(&name) {
                obj.state_mut().set_transform(transform);
                obj.invalidate(DirtyFlags::COORDS);
                changed = true;
            }
        }
        changed
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::movie::LoopMode;
    use crate::object::MoleculeObject;
    use crate::scene::SceneStoreMask;
    use lin_alg::f32::Vec3;
    use patinae_mol::{Atom, CoordSet, Element, ObjectMolecule};

    fn multi_state_molecule(name: &str, states: usize) -> ObjectMolecule {
        let mut mol = ObjectMolecule::new(name);
        mol.add_atom(Atom::new("C", Element::Carbon));
        for state in 0..states {
            mol.add_coord_set(CoordSet::from_vec3(&[Vec3::new(state as f32, 0.0, 0.0)]));
        }
        mol
    }

    #[test]
    fn sync_movie_frame_applies_mset_display_state() {
        let mut session = Session::new();
        session.registry.add(MoleculeObject::with_name(
            multi_state_molecule("mol", 3),
            "mol",
        ));
        session.movie.set_from_spec(vec![1, 2, 3]);
        session.movie.goto_frame(2);

        assert!(session.sync_movie_frame());

        let obj = session.registry.get_molecule("mol").unwrap();
        assert_eq!(obj.display_state(), 2);
    }

    #[test]
    fn update_animations_advances_and_syncs_state() {
        let mut session = Session::new();
        session.registry.add(MoleculeObject::with_name(
            multi_state_molecule("mol", 3),
            "mol",
        ));
        session.refresh_movie_state_count();
        session.movie.set_loop_mode(LoopMode::Loop);
        session.settings.movie.movie_fps = 10.0;
        session.movie.play();

        let update = session.update_animations(0.11);

        assert!(update.movie_frame_changed);
        assert!(update.movie_synced);
        assert!(update.needs_redraw);
        let obj = session.registry.get_molecule("mol").unwrap();
        assert_eq!(session.movie.current_frame(), 1);
        assert_eq!(obj.display_state(), 1);
    }

    #[test]
    fn sync_movie_frame_applies_movie_view_after_scene_view() {
        let mut session = Session::new();
        session.camera.set_fov(20.0);
        session.scenes.store(
            "scene_a",
            SceneStoreMask::VIEW,
            &session.camera,
            &session.registry,
        );

        session.movie.set_frame_count(1);
        let mut movie_view = session.camera.current_view();
        movie_view.fov = 35.0;
        let frame = session.movie.frame_mut(0).unwrap();
        frame.scene_name = Some("scene_a".to_string());
        frame.set_view(movie_view);

        session.camera.set_fov(10.0);
        assert!(session.sync_movie_frame());

        assert_eq!(session.camera.fov(), 35.0);
    }
}
