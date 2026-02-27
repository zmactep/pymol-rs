//! GPU Rendering Pipeline
//!
//! Frame orchestration, shadow passes, molecule preparation,
//! 3D scene rendering, silhouette post-processing, and egui overlay.

use pymol_render::ColorResolver;
use pymol_scene::{setup_uniforms, Object};
use pymol_select::SelectionResult;
use winit::dpi::PhysicalSize;

use super::App;

impl App {
    /// Initialize GPU resources
    pub(crate) async fn init_gpu(
        &mut self,
        window: std::sync::Arc<winit::window::Window>,
    ) -> Result<(String, String), String> {
        // We need to get adapter info before init_gpu consumes it
        // For now, we'll do a quick adapter request just for info
        let adapter = self.view.instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                compatible_surface: None,
                force_fallback_adapter: false,
            })
            .await
            .map_err(|e| format!("No suitable GPU adapter found: {}", e))?;

        let info = adapter.get_info();
        let device_name = info.name.clone();
        let backend = format!("{:?}", info.backend);

        // Drop adapter so init_gpu can create its own
        drop(adapter);

        // Initialize the view's GPU resources
        self.view.init_gpu(window.clone()).await?;

        // Set camera aspect ratio
        let size = window.inner_size();
        self.state.camera.set_aspect(size.width as f32 / size.height as f32);

        Ok((device_name, backend))
    }

    /// Handle window resize
    pub(crate) fn resize(&mut self, new_size: PhysicalSize<u32>) {
        if new_size.width == 0 || new_size.height == 0 {
            return;
        }

        self.view.resize(new_size);

        // Update camera aspect ratio
        self.state.camera.set_aspect(new_size.width as f32 / new_size.height as f32);
        self.needs_redraw = true;
    }

    /// Render a frame — thin orchestrator.
    pub(crate) fn render(&mut self) -> Result<(), wgpu::SurfaceError> {
        // Extract config values before any mutable borrows.
        let (width, height) = {
            let config = self.view.surface_config.as_ref().unwrap();
            (config.width, config.height)
        };
        let scale_factor = self.view.scale_factor();

        // Run egui UI first (mutably borrows self to process commands/actions).
        let egui_output = self.run_egui_ui();

        // Update camera aspect ratio based on viewport dimensions.
        let (viewport_width, viewport_height) = if let Some(vp) = &self.view.viewport_rect {
            (vp.width().max(1.0), vp.height().max(1.0))
        } else {
            (width as f32, height as f32)
        };
        self.state.camera.set_aspect(viewport_width / viewport_height);

        // Get surface texture.
        let surface = self.view.surface.as_ref().unwrap();
        let output = surface.get_current_texture()?;
        let output_view = output.texture.create_view(&wgpu::TextureViewDescriptor::default());

        // Upload scene uniforms and set the active shading mode.
        let shading_mode = pymol_settings::ShadingMode::from_settings(&self.state.settings);
        {
            let context = self.view.render_context.as_mut().unwrap();
            self.shading.set_mode(shading_mode, context);
            let uniforms = setup_uniforms(
                &self.state.camera,
                &self.state.settings,
                self.state.clear_color,
                (viewport_width, viewport_height),
            );
            context.update_uniforms(&uniforms);
        }

        // Prepare molecule GPU geometry; detect geometry changes.
        let names: Vec<String> = self.state.registry.names().map(|s: &str| s.to_string()).collect();
        self.prepare_molecules(&names);

        // Create command encoder.
        let mut encoder = {
            let context = self.view.render_context.as_ref().unwrap();
            context.device().create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("Render Encoder"),
            })
        };

        // Give the active pipeline the current scene bounding sphere so it can
        // compute pre-pass matrices (e.g. shadow projections). Modes that don't
        // need it ignore this call.
        if let Some((min, max)) = self.state.registry.extent() {
            let center = [
                (min.x + max.x) * 0.5,
                (min.y + max.y) * 0.5,
                (min.z + max.z) * 0.5,
            ];
            let dx = max.x - min.x;
            let dy = max.y - min.y;
            let dz = max.z - min.z;
            let radius = (dx * dx + dy * dy + dz * dz).sqrt() * 0.5;
            self.shading.set_scene_bounds(center, radius.max(1.0));
        }

        // Run the active pipeline's prepare step. Returns true when the pipeline
        // requires additional render passes this frame (e.g. shadow depth passes).
        let need_shadow_passes = {
            let context = self.view.render_context.as_mut().unwrap();
            self.shading.prepare(context, &self.state.settings)
        };

        // Execute pre-passes if required, then notify the pipeline they're done.
        if need_shadow_passes {
            self.render_shadow_passes(&mut encoder, &names);
            self.shading.finish_shadow_passes();
        }

        // 3D scene pass (opaque + transparent).
        self.render_molecules(&mut encoder, &output_view, &names);

        // Post-process: silhouette edge detection.
        self.render_silhouettes(&mut encoder, &output_view, width, height);

        // egui UI overlay.
        self.render_egui(&mut encoder, &output_view, egui_output, width, height, scale_factor);

        // Submit and present.
        let context = self.view.render_context.as_ref().unwrap();
        context.queue().submit(std::iter::once(encoder.finish()));
        output.present();

        Ok(())
    }

    /// Render shadow depth passes requested by the active shading pipeline.
    ///
    /// Uses `ShadowPassState` from the pipeline — fully mode-agnostic.
    fn render_shadow_passes(&mut self, encoder: &mut wgpu::CommandEncoder, names: &[String]) {
        let context = self.view.render_context.as_ref().unwrap();

        let state = match self.shading.shadow_pass_state() {
            Some(s) => s,
            None => return,
        };

        let bind_group = state.pipelines.create_bind_group(context.device());

        for (i, matrix) in state.matrices.iter().enumerate() {
            state.pipelines.update_uniform(context.queue(), matrix);
            let (vp_x, vp_y, vp_w, vp_h) = state.atlas.tile_viewport(i as u32);

            let mut shadow_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("Shadow Depth Pass"),
                color_attachments: &[],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &state.atlas.depth_view,
                    depth_ops: Some(wgpu::Operations {
                        load: if i == 0 { wgpu::LoadOp::Clear(1.0) } else { wgpu::LoadOp::Load },
                        store: wgpu::StoreOp::Store,
                    }),
                    stencil_ops: None,
                }),
                timestamp_writes: None,
                occlusion_query_set: None,
            });

            shadow_pass.set_viewport(vp_x as f32, vp_y as f32, vp_w as f32, vp_h as f32, 0.0, 1.0);
            shadow_pass.set_scissor_rect(vp_x, vp_y, vp_w, vp_h);

            for name in names {
                if let Some(mol_obj) = self.state.registry.get_molecule(name) {
                    mol_obj.render_shadow_depth(&mut shadow_pass, state.pipelines, &bind_group, context);
                }
            }
        }
    }

    /// Prepare molecule GPU geometry for the frame.
    ///
    /// Calls `prepare_render` on each molecule, updates selection indicators,
    /// and marks shadows dirty if any geometry changed.
    fn prepare_molecules(&mut self, names: &[String]) {
        let selection_results = self.evaluate_visible_selections();
        let selection_width = self.state.settings.get_float(pymol_settings::id::selection_width).max(6.0);
        let mouse_selection_mode = self.state.settings.get_int(pymol_settings::id::mouse_selection_mode);

        // Snapshot hover state to avoid borrow issues
        let hover_hit = self.hover_hit.clone();
        let sequence_hover = self.sequence_hover.clone();

        // Sequence hover uses Ctrl/Cmd for exclusion mode
        let ctrl_held = self.input.ctrl_or_cmd_held();

        const COLOR_YELLOW: [f32; 4] = [1.0, 1.0, 0.5, 0.7];
        const COLOR_RED: [f32; 4] = [1.0, 0.0, 0.0, 1.0];

        let mut geometry_changed = false;
        let context = self.view.render_context.as_ref().unwrap();
        for name in names {
            let color_resolver = ColorResolver::new(
                &self.state.named_colors,
                &self.state.element_colors,
            );
            if let Some(mol_obj) = self.state.registry.get_molecule_mut(name) {
                if mol_obj.is_dirty() {
                    geometry_changed = true;
                }
                mol_obj.prepare_render(context, color_resolver, &self.state.settings);

                // Selection indicator
                let sele_for_obj = selection_results.iter().find(|(n, _)| n == name).map(|(_, s)| s);
                if let Some(sel) = sele_for_obj {
                    log::debug!("Setting selection indicator for '{}' with {} atoms", name, sel.count());
                    mol_obj.set_selection_indicator_with_size(sel, context, Some(selection_width));
                } else {
                    mol_obj.clear_selection_indicator();
                }

                // Hover indicator: 3D viewport pick takes priority, then sequence hover.
                if let Some(ref hit) = hover_hit {
                    if hit.object_name == *name {
                        let sel = pymol_scene::expand_pick_to_selection(hit, mouse_selection_mode, mol_obj.molecule());
                        // 3D viewport: red if hovered atoms overlap sele, yellow otherwise
                        let overlaps_sele = sele_for_obj
                            .is_some_and(|sele| sel.intersection(sele).any());
                        let color = if overlaps_sele { COLOR_RED } else { COLOR_YELLOW };
                        mol_obj.set_hover_indicator(&sel, context, selection_width, color);
                    } else {
                        mol_obj.clear_hover_indicator();
                    }
                } else if let Some(ref hover) = sequence_hover {
                    if hover.object_name == *name {
                        let mol = mol_obj.molecule();
                        let sel = SelectionResult::from_indices(
                            mol.atom_count(),
                            mol.atoms_indexed().filter_map(|(idx, atom)| {
                                if atom.residue.key.chain == hover.chain_id
                                    && atom.residue.key.resv == hover.resv
                                {
                                    Some(idx)
                                } else {
                                    None
                                }
                            }),
                        );
                        // Sequence hover: red with Ctrl (exclusion), yellow without
                        let color = if ctrl_held { COLOR_RED } else { COLOR_YELLOW };
                        mol_obj.set_hover_indicator(&sel, context, selection_width, color);
                    } else {
                        mol_obj.clear_hover_indicator();
                    }
                } else {
                    mol_obj.clear_hover_indicator();
                }
            }
        }

        // Prepare measurement objects
        for name in names {
            if let Some(meas_obj) = self.state.registry.get_measurement_mut(name) {
                meas_obj.prepare_render(context);
            }
        }

        if geometry_changed {
            self.shading.invalidate_shadows();
        }
    }

    /// 3D scene render pass — opaque + transparent molecules.
    fn render_molecules(
        &mut self,
        encoder: &mut wgpu::CommandEncoder,
        output_view: &wgpu::TextureView,
        names: &[String],
    ) {
        let context = self.view.render_context.as_ref().unwrap();
        let depth_view = self.view.depth_view.as_ref().unwrap();

        let mut render_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
            label: Some("3D Render Pass"),
            color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                view: output_view,
                resolve_target: None,
                ops: wgpu::Operations {
                    load: wgpu::LoadOp::Clear(wgpu::Color {
                        r: self.state.clear_color[0] as f64,
                        g: self.state.clear_color[1] as f64,
                        b: self.state.clear_color[2] as f64,
                        a: 1.0,
                    }),
                    store: wgpu::StoreOp::Store,
                },
                depth_slice: None,
            })],
            depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                view: depth_view,
                depth_ops: Some(wgpu::Operations {
                    load: wgpu::LoadOp::Clear(1.0),
                    store: wgpu::StoreOp::Store,
                }),
                stencil_ops: None,
            }),
            timestamp_writes: None,
            occlusion_query_set: None,
        });

        // Restrict 3D rendering to the central viewport area, clamped to surface bounds.
        if let Some(vp) = &self.view.viewport_rect {
            let (surf_w, surf_h) = self
                .view
                .surface_config
                .as_ref()
                .map(|c| (c.width, c.height))
                .unwrap_or((1, 1));
            let vp_x = (vp.min.x.max(0.0) as u32).min(surf_w.saturating_sub(1));
            let vp_y = (vp.min.y.max(0.0) as u32).min(surf_h.saturating_sub(1));
            let vp_w = (vp.width().max(1.0) as u32).min(surf_w - vp_x);
            let vp_h = (vp.height().max(1.0) as u32).min(surf_h - vp_y);
            render_pass.set_viewport(vp_x as f32, vp_y as f32, vp_w as f32, vp_h as f32, 0.0, 1.0);
            render_pass.set_scissor_rect(vp_x, vp_y, vp_w, vp_h);
        }

        for name in names {
            if let Some(mol_obj) = self.state.registry.get_molecule(name) {
                if mol_obj.is_enabled() {
                    mol_obj.render(&mut render_pass, context);
                }
            }
        }

        // Render measurement objects (dashed lines)
        for name in names {
            if let Some(meas_obj) = self.state.registry.get_measurement(name) {
                if meas_obj.is_enabled() {
                    meas_obj.render(&mut render_pass, context);
                }
            }
        }
    }

    /// Post-process silhouette edge detection pass.
    fn render_silhouettes(
        &mut self,
        encoder: &mut wgpu::CommandEncoder,
        output_view: &wgpu::TextureView,
        width: u32,
        height: u32,
    ) {
        if !self.state.settings.get_bool(pymol_settings::id::silhouettes) {
            return;
        }
        if let (Some(silhouette), Some(depth_view)) =
            (&self.view.silhouette_pipeline, &self.view.depth_view)
        {
            let context = self.view.render_context.as_ref().unwrap();
            let thickness = self.state.settings.get_float(pymol_settings::id::silhouette_width);
            let depth_jump = self.state.settings.get_float(pymol_settings::id::silhouette_depth_jump);
            let color_int = self.state.settings.get_color(pymol_settings::id::silhouette_color);
            let color = pymol_color::Color::from_packed_rgb(color_int).to_rgba(1.0);
            silhouette.render(
                encoder,
                context.queue(),
                context.device(),
                output_view,
                depth_view,
                width,
                height,
                thickness,
                depth_jump,
                color,
                None,
            );
        }
    }

    /// egui UI overlay pass.
    fn render_egui(
        &mut self,
        encoder: &mut wgpu::CommandEncoder,
        output_view: &wgpu::TextureView,
        egui_output: Option<(Vec<egui::ClippedPrimitive>, egui::TexturesDelta)>,
        width: u32,
        height: u32,
        scale_factor: f32,
    ) {
        let (clipped_primitives, textures_delta) = match egui_output {
            Some(o) => o,
            None => return,
        };
        let egui_renderer = match &mut self.view.egui_renderer {
            Some(r) => r,
            None => return,
        };

        let context = self.view.render_context.as_ref().unwrap();
        let device = context.device();
        let queue = context.queue();

        for (id, image_delta) in &textures_delta.set {
            egui_renderer.update_texture(device, queue, *id, image_delta);
        }

        let screen_descriptor = egui_wgpu::ScreenDescriptor {
            size_in_pixels: [width, height],
            pixels_per_point: scale_factor,
        };
        egui_renderer.update_buffers(device, queue, encoder, &clipped_primitives, &screen_descriptor);

        {
            let render_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("egui Render Pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: output_view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Load,
                        store: wgpu::StoreOp::Store,
                    },
                    depth_slice: None,
                })],
                depth_stencil_attachment: None,
                timestamp_writes: None,
                occlusion_query_set: None,
            });
            let mut render_pass = render_pass.forget_lifetime();
            egui_renderer.render(&mut render_pass, &clipped_primitives, &screen_descriptor);
        }

        for id in &textures_delta.free {
            egui_renderer.free_texture(id);
        }
    }
}
