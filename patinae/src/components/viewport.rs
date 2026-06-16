//! 3D viewport renderer for the Patinae (Slint) GUI.
//!
//! Pure GPU renderer — owns wgpu resources and shading pipelines.
//! Input handling and viewport state live in the application bridge layer.
//!
//! Drives `patinae_render::RenderState` for both live frames and headless
//! PNG capture (via `CaptureRenderer`).

use std::path::Path;
use std::sync::Arc;
use std::time::{Duration, Instant};

use patinae_color::{NamedPalette, ThemedPalette};
use patinae_render::picking::readback::PendingPick;
use patinae_render::{
    DisplayedGeometry, GeometryExportOptions, RenderArtifactSnapshot, RenderInput, RenderState,
    SceneLod, TraceGeometryChunk,
};
use patinae_scene::{Camera, CaptureRenderer, ObjectRegistry, PickHit, Session, ViewerError};
use patinae_settings::{ResolvedSettings, Settings, ShadingMode};

use patinae_scene::bridge::{
    frame_uniforms_from_camera, frame_uniforms_from_session, resolve_pick, resolve_setting_color,
    visit_render_scene, CachedRenderScene, ResolvedSceneColors, ResolvedSceneMarkers,
};

const VIEWPORT_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Rgba8Unorm;
const LARGE_LOD_HOVER_PICK_INTERVAL: Duration = Duration::from_millis(33);

pub struct AdapterInfo {
    pub device_name: String,
    pub backend: String,
}

/// GPU viewport that renders a [`Session`] to an offscreen wgpu texture.
pub struct ViewportRenderer {
    /// Live-render path — patinae-render.
    state: RenderState,
    /// Last viewport size — used to detect resize before each frame.
    last_viewport_size: Option<(u32, u32)>,
    /// Persistent host-side render cache shared with the web viewer.
    render_scene: CachedRenderScene,
    /// In-flight async hover pick. At most one outstanding handle.
    pending_hover: Option<PendingPick>,
    /// Last hover readback submission, used to pace 3J3Q-class sphere scenes.
    last_hover_submit: Option<Instant>,
    /// Ring of offscreen color textures handed to Slint. Rotating across
    /// 3 textures lets us write the next frame's image while Slint is
    /// still reading the previous frame — otherwise the GPU serialises
    /// our write on Slint's read fence, capping fps independent of CPU.
    color_textures: Vec<wgpu::Texture>,
    color_textures_next: usize,
}

impl ViewportRenderer {
    pub fn new(device: wgpu::Device, queue: wgpu::Queue) -> Self {
        // wgpu::Device / wgpu::Queue clone bumps an internal Arc refcount.
        let device_arc = Arc::new(device);
        let queue_arc = Arc::new(queue);
        let state = RenderState::new(device_arc, queue_arc, VIEWPORT_FORMAT, (1, 1));
        Self {
            state,
            last_viewport_size: None,
            render_scene: CachedRenderScene::default(),
            pending_hover: None,
            last_hover_submit: None,
            color_textures: Vec::new(),
            color_textures_next: 0,
        }
    }

    /// Create from a Slint graphics API handle. Returns the renderer and
    /// optional adapter info (for HUD display).
    pub fn setup(api: &slint::GraphicsAPI) -> Option<(Self, Option<AdapterInfo>)> {
        let slint::GraphicsAPI::WGPU28 {
            instance,
            device,
            queue,
            ..
        } = api
        else {
            return None;
        };

        let info = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            ..Default::default()
        }))
        .ok()
        .map(|adapter| {
            let info = adapter.get_info();
            AdapterInfo {
                device_name: info.name,
                backend: format!("{:?}", info.backend),
            }
        });

        let renderer = Self::new(device.clone(), queue.clone());
        Some((renderer, info))
    }

    /// Render the scene and return the result as a Slint image.
    /// Drain whatever FrameStats are ready right now (gated by the
    /// `patinae-render/stats` feature, which patinae enables by default).
    pub fn take_frame_stats_for_log(&self) -> Option<patinae_render::FrameStats> {
        self.state.take_frame_stats()
    }

    pub fn render(
        &mut self,
        session: &mut Session,
        width: u32,
        height: u32,
    ) -> Option<slint::Image> {
        let texture = self.render_frame(session, width, height);
        match slint::Image::try_from(texture) {
            Ok(image) => Some(image),
            Err(e) => {
                log::error!("Failed to import texture: {:?}", e);
                None
            }
        }
    }

    fn render_frame(&mut self, session: &mut Session, width: u32, height: u32) -> wgpu::Texture {
        // Resize state's internal targets if the viewport changed.
        let resized = self.last_viewport_size != Some((width, height));
        if resized {
            self.state.resize((width, height));
            self.last_viewport_size = Some((width, height));
            // Drop the cached ring so it's re-created at the new size.
            self.color_textures.clear();
            self.color_textures_next = 0;
        }

        // 1) Frame uniforms.
        let frame = frame_uniforms_from_session(session, (width, height), session.clear_color);
        self.state.uniforms = frame;
        self.state.ctx.upload_frame(&self.state.uniforms);
        // Live viewports should show the RGB background; PNG capture keeps
        // using `opaque_background` to decide exported alpha.
        self.state.set_clear_color(session.clear_color, true);
        self.state
            .set_marking_width(session.settings.ui.selection_width);

        // 2) Build cached host-side render input, then sync the renderer.
        // The cache owns expensive per-atom color/marker buffers and the
        // picking name lookup. The frame input borrows from `session` only
        // until `RenderState::sync` returns.
        let frame = self.render_scene.prepare(session);
        self.state.sync(&frame.render_input());
        drop(frame);
        // patinae-render reps track dirty internally; consume the host-side
        // flag set by commands so SceneModel::sync stops perpetually
        // rebuilding the objects panel (which would destroy panel
        // TouchAreas every frame and break hover/click detection).
        session.registry.clear_all_dirty_molecules();
        session.registry.clear_all_dirty_maps();

        // 5) Marker bits already flow through `render_scene` →
        // `RenderObjectInput.atom_markers` → `marker_lut`. The legacy
        // `apply_highlight` / `selection_mut` path is gone; the `max_atoms`
        // value is no longer consumed here.

        // 6) Silhouette pass on/off from settings.
        let common = &session.settings.shading.common;
        let silhouette_ink = resolve_setting_color(
            common.silhouette_color,
            &session.named_palette,
            [0.0, 0.0, 0.0, 1.0],
        );
        if common.silhouettes {
            self.state
                .set_silhouette(true, common.silhouette_width, silhouette_ink);
        } else {
            self.state.set_silhouette(false, 1.0, silhouette_ink);
        }

        let full = &session.settings.shading.full;
        self.state.set_shadows(
            session.settings.shading.mode == ShadingMode::Full,
            full.shadow_map_size.max(64) as u32,
            full.shadow_bias,
            full.shadow_intensity,
            full.shadow_pcf.max(1) as u32,
        );
        // 6b) Ambient-occlusion config. Skripkin uses a multi-directional
        // shadow atlas; Classic / Full use the explicit screen-space AO
        // settings.
        configure_ambient_occlusion(&mut self.state, &session.settings);
        // 6c) FXAA postprocess. Default-enabled; toggle via
        // `set fxaa_enabled, 0`.
        self.state.set_fxaa(effective_fxaa(&session.settings));

        // 7) Build (or reuse) the offscreen colour texture handed back to
        // Slint each frame and run the FrameGraph. We keep a ring of
        // textures (default 3) and rotate. Slint may still be sampling
        // the previous frame's texture when we kick off this frame —
        // writing into a fresh one lets the GPU pipeline both passes
        // instead of serialising on the read fence.
        const COLOR_TEXTURE_RING: usize = 3;
        if self.color_textures.is_empty() {
            self.color_textures = (0..COLOR_TEXTURE_RING)
                .map(|i| {
                    self.state
                        .ctx
                        .device
                        .create_texture(&wgpu::TextureDescriptor {
                            label: Some(match i {
                                0 => "Patinae Viewport 0",
                                1 => "Patinae Viewport 1",
                                _ => "Patinae Viewport 2",
                            }),
                            size: wgpu::Extent3d {
                                width,
                                height,
                                depth_or_array_layers: 1,
                            },
                            mip_level_count: 1,
                            sample_count: 1,
                            dimension: wgpu::TextureDimension::D2,
                            format: VIEWPORT_FORMAT,
                            usage: wgpu::TextureUsages::RENDER_ATTACHMENT
                                | wgpu::TextureUsages::TEXTURE_BINDING,
                            view_formats: &[],
                        })
                })
                .collect();
            self.color_textures_next = 0;
        }
        let idx = self.color_textures_next;
        self.color_textures_next = (idx + 1) % self.color_textures.len();
        let color_texture = &self.color_textures[idx];
        let color_view = color_texture.create_view(&wgpu::TextureViewDescriptor::default());

        let mut encoder =
            self.state
                .ctx
                .device
                .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                    label: Some("Patinae Render Encoder"),
                });
        self.state.render(&color_view, &mut encoder);

        // PATINAE_GPU_SYNC=1: force device.poll(Wait) after submit and
        // measure real GPU wall-clock per frame. Slint's wgpu device is
        // created without TIMESTAMP_QUERY, so per-pass timing is gone —
        // this is the only way to see whether GPU work alone busts the
        // vsync budget on big assemblies.
        let gpu_sync = std::env::var("PATINAE_GPU_SYNC")
            .ok()
            .map(|v| v == "1")
            .unwrap_or(false);
        let t_submit = std::time::Instant::now();
        self.state
            .ctx
            .queue
            .submit(std::iter::once(encoder.finish()));
        if gpu_sync {
            self.state
                .ctx
                .device
                .poll(wgpu::PollType::Wait {
                    submission_index: None,
                    timeout: None,
                })
                .ok();
            let dt_ms = t_submit.elapsed().as_secs_f64() * 1000.0;
            eprintln!("[patinae] gpu submit+wait: {:.2} ms", dt_ms);
        }

        self.color_textures[idx].clone()
    }

    /// Async hover pick: skip while a previous hover readback is in flight.
    /// - `Some(Some(hit))` — fresh hit landed this frame.
    /// - `Some(None)` — fresh miss landed this frame.
    /// - `None` — no fresh result (one is in flight or being submitted).
    pub fn poll_hover_pick(
        &mut self,
        session: &Session,
        x: u32,
        y: u32,
    ) -> Option<Option<PickHit>> {
        let (width, height) = self.last_viewport_size?;
        if x >= width || y >= height {
            return None;
        }
        self.refresh_pick_frame_uniforms(session, (width, height));

        if let Some(pending) = self.pending_hover.as_ref() {
            // Wait for the GPU to finish; never submit a second pick.
            let result = self.state.try_collect_pick(pending)?;
            self.pending_hover = None;
            return Some(
                result.and_then(|hit| resolve_pick(hit, self.render_scene.object_names(), session)),
            );
        }

        // No pick in flight — submit a fresh one. `submit_pick` returns
        // `None` if the renderer was built with `PickingMode::Disabled`;
        // in that case there's nothing to wait on.
        if self.large_lod_hover_pick_waiting() {
            return None;
        }
        self.pending_hover = self.state.submit_pick(x, y);
        if self.pending_hover.is_some() {
            self.last_hover_submit = Some(Instant::now());
        }
        None
    }

    fn large_lod_hover_pick_waiting(&self) -> bool {
        (self.state.sphere_lod_diagnostics().active || self.state.stick_lod_diagnostics().active)
            && self
                .last_hover_submit
                .is_some_and(|last| last.elapsed() < LARGE_LOD_HOVER_PICK_INTERVAL)
    }

    /// Synchronous click pick. Uses a separate renderer readback from hover,
    /// so it never waits for a pending hover query before resolving a click.
    pub fn pick_at_click(&mut self, session: &Session, x: u32, y: u32) -> Option<PickHit> {
        let (width, height) = self.last_viewport_size?;
        if x >= width || y >= height {
            return None;
        }
        self.refresh_pick_frame_uniforms(session, (width, height));

        let hit = self.state.pick(x, y)?;
        resolve_pick(hit, self.render_scene.object_names(), session)
    }

    fn refresh_pick_frame_uniforms(&mut self, session: &Session, viewport: (u32, u32)) {
        let frame = frame_uniforms_from_session(session, viewport, session.clear_color);
        self.state.uniforms = frame;
        self.state.ctx.upload_frame(&self.state.uniforms);
    }
}

impl CaptureRenderer for ViewportRenderer {
    fn gpu_device(&self) -> &Arc<wgpu::Device> {
        &self.state.ctx.device
    }

    fn gpu_queue(&self) -> &Arc<wgpu::Queue> {
        &self.state.ctx.queue
    }

    fn capture_png(
        &mut self,
        path: &Path,
        width: u32,
        height: u32,
        camera: &mut Camera,
        registry: &mut ObjectRegistry,
        settings: &Settings,
        named: &NamedPalette,
        themed: &ThemedPalette,
        clear_color: [f32; 3],
    ) -> Result<(), ViewerError> {
        // Build per-atom colours + render objects via the same bridges the
        // live frame uses, so the captured PNG matches the on-screen scene
        // (modulo resolution).
        let colors = ResolvedSceneColors::build(registry, settings, named, themed);
        // Capture path doesn't have access to a `Session`; render with no
        // markers (selection / hover overlays are interactive UI affordances
        // and don't belong in offscreen captures).
        let markers = ResolvedSceneMarkers::default();
        let mut object_inputs = Vec::new();
        let mut map_inputs = Vec::new();
        visit_render_scene(
            registry,
            settings,
            &colors,
            &markers,
            &mut |_name, obj| object_inputs.push(obj),
            &mut |_name, map| map_inputs.push(map),
        );

        let frame = frame_uniforms_from_camera(camera, settings, (width, height), clear_color);

        let resolved = ResolvedSettings::resolve(settings, None);
        let lod = object_inputs
            .first()
            .map(|o| o.lod)
            .unwrap_or(SceneLod::Auto);
        let input = RenderInput {
            objects: &object_inputs,
            maps: &map_inputs,
            settings: &resolved,
            lod,
        };

        // Apply the silhouette setting before render.
        let common = &settings.shading.common;
        let silhouette_ink =
            resolve_setting_color(common.silhouette_color, named, [0.0, 0.0, 0.0, 1.0]);
        if common.silhouettes {
            self.state
                .set_silhouette(true, common.silhouette_width, silhouette_ink);
        } else {
            self.state.set_silhouette(false, 1.0, silhouette_ink);
        }
        let full = &settings.shading.full;
        self.state.set_shadows(
            settings.shading.mode == ShadingMode::Full,
            full.shadow_map_size.max(64) as u32,
            full.shadow_bias,
            full.shadow_intensity,
            full.shadow_pcf.max(1) as u32,
        );
        configure_ambient_occlusion(&mut self.state, settings);
        self.state
            .set_clear_color(clear_color, settings.ui.opaque_background);
        self.state.set_fxaa(effective_fxaa(settings));

        patinae_render::capture::capture_png(&mut self.state, path, width, height, &frame, &input)
            .map_err(|e| ViewerError::capture_error(e.to_string()))?;
        registry.clear_all_dirty_molecules();
        registry.clear_all_dirty_maps();

        // Restore live viewport size on the next frame — `render_frame`
        // detects the dimension mismatch and resizes back.
        self.last_viewport_size = None;
        Ok(())
    }

    fn export_displayed_geometry(
        &mut self,
        _camera: &mut Camera,
        registry: &mut ObjectRegistry,
        settings: &Settings,
        named: &NamedPalette,
        themed: &ThemedPalette,
        _clear_color: [f32; 3],
        options: &GeometryExportOptions,
    ) -> Result<DisplayedGeometry, ViewerError> {
        let colors = ResolvedSceneColors::build(registry, settings, named, themed);
        let markers = ResolvedSceneMarkers::default();
        let mut object_inputs = Vec::new();
        let mut map_inputs = Vec::new();
        visit_render_scene(
            registry,
            settings,
            &colors,
            &markers,
            &mut |_name, obj| object_inputs.push(obj),
            &mut |_name, map| map_inputs.push(map),
        );

        let resolved = ResolvedSettings::resolve(settings, None);
        let lod = object_inputs
            .first()
            .map(|o| o.lod)
            .unwrap_or(SceneLod::Auto);
        let input = RenderInput {
            objects: &object_inputs,
            maps: &map_inputs,
            settings: &resolved,
            lod,
        };

        let geometry = self
            .state
            .export_displayed_geometry(&input, options)
            .map_err(|e| ViewerError::capture_error(e.to_string()))?;
        registry.clear_all_dirty_molecules();
        registry.clear_all_dirty_maps();

        // Restore live viewport size on the next frame, matching capture.
        self.last_viewport_size = None;
        Ok(geometry)
    }

    fn for_each_trace_geometry_chunk(
        &mut self,
        _camera: &mut Camera,
        registry: &mut ObjectRegistry,
        settings: &Settings,
        named: &NamedPalette,
        themed: &ThemedPalette,
        _clear_color: [f32; 3],
        options: &GeometryExportOptions,
        visitor: &mut dyn FnMut(TraceGeometryChunk) -> Result<(), String>,
    ) -> Result<(), ViewerError> {
        let colors = ResolvedSceneColors::build(registry, settings, named, themed);
        let markers = ResolvedSceneMarkers::default();
        let mut object_inputs = Vec::new();
        let mut map_inputs = Vec::new();
        visit_render_scene(
            registry,
            settings,
            &colors,
            &markers,
            &mut |_name, obj| object_inputs.push(obj),
            &mut |_name, map| map_inputs.push(map),
        );

        let resolved = ResolvedSettings::resolve(settings, None);
        let lod = object_inputs
            .first()
            .map(|o| o.lod)
            .unwrap_or(SceneLod::Auto);
        let input = RenderInput {
            objects: &object_inputs,
            maps: &map_inputs,
            settings: &resolved,
            lod,
        };

        self.state
            .for_each_trace_geometry_chunk(&input, options, visitor)
            .map_err(|e| ViewerError::capture_error(e.to_string()))?;
        registry.clear_all_dirty_molecules();
        registry.clear_all_dirty_maps();

        self.last_viewport_size = None;
        Ok(())
    }

    fn visit_render_artifacts(
        &mut self,
        _camera: &mut Camera,
        registry: &mut ObjectRegistry,
        settings: &Settings,
        named: &NamedPalette,
        themed: &ThemedPalette,
        _clear_color: [f32; 3],
        visitor: &mut dyn FnMut(RenderArtifactSnapshot<'_>) -> Result<(), String>,
    ) -> Result<(), ViewerError> {
        let colors = ResolvedSceneColors::build(registry, settings, named, themed);
        let markers = ResolvedSceneMarkers::default();
        let mut object_inputs = Vec::new();
        let mut map_inputs = Vec::new();
        visit_render_scene(
            registry,
            settings,
            &colors,
            &markers,
            &mut |_name, obj| object_inputs.push(obj),
            &mut |_name, map| map_inputs.push(map),
        );

        let resolved = ResolvedSettings::resolve(settings, None);
        let lod = object_inputs
            .first()
            .map(|o| o.lod)
            .unwrap_or(SceneLod::Auto);
        let input = RenderInput {
            objects: &object_inputs,
            maps: &map_inputs,
            settings: &resolved,
            lod,
        };

        let snapshot = self.state.render_artifact_snapshot(&input);
        visitor(snapshot).map_err(ViewerError::capture_error)?;
        registry.clear_all_dirty_molecules();
        registry.clear_all_dirty_maps();

        self.last_viewport_size = None;
        Ok(())
    }
}

fn effective_fxaa(settings: &Settings) -> bool {
    settings.ui.antialias > 0 && settings.fxaa.enabled
}

fn configure_ambient_occlusion(state: &mut RenderState, settings: &Settings) {
    let ssao = &settings.ssao;
    if settings.shading.mode == ShadingMode::Skripkin {
        let skripkin = &settings.shading.skripkin;
        state.set_ssao(false, ssao.radius, ssao.intensity, ssao.bias);
        state.set_skripkin_ao(
            true,
            skripkin.directions.max(1) as u32,
            skripkin.map_size.max(32) as u32,
            skripkin.bias,
            skripkin.intensity,
        );
    } else {
        state.set_skripkin_ao(false, 1, 32, 0.0, 0.0);
        state.set_ssao(ssao.enabled, ssao.radius, ssao.intensity, ssao.bias);
    }
}
