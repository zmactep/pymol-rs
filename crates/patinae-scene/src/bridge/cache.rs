//! Cached render-scene bridge state.
//!
//! Hosts keep this cache alive across frames so per-atom colors, marker bits,
//! and picking name lookups are rebuilt only when scene state changes.

use std::cell::RefCell;

use patinae_render::{RenderInput, RenderMapInput, RenderObjectInput, SceneLod};
use patinae_settings::ResolvedSettings;

use crate::session::Session;

use super::{
    picking::render_id_slot_index, visit_render_scene, ResolvedSceneColors, ResolvedSceneMarkers,
};

/// Persistent host-side cache for renderer input.
///
/// The renderer input itself borrows from the current [`Session`], so each
/// frame still owns short-lived input vectors. The expensive per-atom buffers
/// and the sparse render-id-to-name picking lookup persist here.
#[derive(Default)]
pub struct CachedRenderScene {
    colors: ResolvedSceneColors,
    markers: ResolvedSceneMarkers,
    object_names: Vec<Option<String>>,
}

impl CachedRenderScene {
    /// Builds a frame input using cached color and marker buffers.
    pub fn prepare<'a>(&'a mut self, session: &'a mut Session) -> CachedRenderFrame<'a> {
        if self.colors.needs_rebuild(&session.registry) {
            self.colors.rebuild(
                &session.registry,
                &session.settings,
                &session.named_palette,
                &session.palette,
            );
        }

        self.markers.rebuild(
            &mut session.selections,
            &session.registry,
            session.hover_target.as_ref(),
        );

        let mut objects = Vec::new();
        let mut maps = Vec::new();
        {
            let names = RefCell::new(&mut self.object_names);
            names.borrow_mut().clear();
            visit_render_scene(
                &session.registry,
                &session.settings,
                &self.colors,
                &self.markers,
                &mut |name, obj| {
                    record_object_name(&mut names.borrow_mut(), obj.object_id.0, name);
                    objects.push(obj);
                },
                &mut |name, map| {
                    record_object_name(&mut names.borrow_mut(), map.object_id.0, name);
                    maps.push(map);
                },
            );
        }

        let settings = ResolvedSettings::resolve(&session.settings, None);
        let lod = objects.first().map(|o| o.lod).unwrap_or(SceneLod::Auto);

        CachedRenderFrame {
            objects,
            maps,
            settings,
            lod,
        }
    }

    /// Returns sparse names indexed by `RenderObjectId::slot_index()`.
    pub fn object_names(&self) -> &[Option<String>] {
        &self.object_names
    }
}

fn record_object_name(names: &mut Vec<Option<String>>, object_id: u32, name: &str) {
    let Some(idx) = render_id_slot_index(object_id) else {
        return;
    };
    if names.len() <= idx {
        names.resize_with(idx + 1, || None);
    }
    names[idx] = Some(name.to_string());
}

/// Short-lived renderer input for a single frame.
pub struct CachedRenderFrame<'a> {
    objects: Vec<RenderObjectInput<'a>>,
    maps: Vec<RenderMapInput<'a>>,
    settings: ResolvedSettings,
    lod: SceneLod,
}

impl<'a> CachedRenderFrame<'a> {
    /// Returns borrowed render input for [`patinae_render::RenderState::sync`].
    pub fn render_input(&self) -> RenderInput<'_> {
        RenderInput {
            objects: &self.objects,
            maps: &self.maps,
            settings: &self.settings,
            lod: self.lod,
        }
    }
}
