//! Frame rendering pipeline — port of pymol-gui render.rs without egui.

use pymol_render::{ColorResolver, ShadingManager};
use pymol_scene::{setup_uniforms, Object, Session};

use crate::gpu::GpuState;

/// Render one frame to the canvas surface.
///
/// This is the core rendering pipeline, ported from `pymol-gui/src/app/render.rs`
/// with egui overlay removed.
pub fn render_frame(
    gpu: &mut GpuState,
    session: &mut Session,
    shading: &mut ShadingManager,
) -> Result<(), wgpu::SurfaceError> {
    let width = gpu.surface_config.width;
    let height = gpu.surface_config.height;

    // Setup uniforms
    let shading_mode = pymol_settings::ShadingMode::from_settings(&session.settings);
    shading.set_mode(shading_mode, &mut gpu.render_context);
    let uniforms = setup_uniforms(
        &session.camera,
        &session.settings,
        session.clear_color,
        (width as f32, height as f32),
    );
    gpu.render_context.update_uniforms(&uniforms);

    // Update selection indicators
    let names: Vec<String> = session.registry.names().map(|s| s.to_string()).collect();
    let selection_results = session.selections.evaluate_visible(&session.registry, Default::default());
    let selection_width = session.settings.get_float(pymol_settings::id::selection_width).max(6.0);

    for name in &names {
        if let Some(mol_obj) = session.registry.get_molecule_mut(name) {
            let sele_for_obj = selection_results.iter().find(|(n, _)| n == name).map(|(_, s)| s);
            if let Some(sel) = sele_for_obj {
                mol_obj.set_selection_indicator_with_size(sel, &gpu.render_context, Some(selection_width));
            } else {
                mol_obj.clear_selection_indicator();
            }
        }
    }

    // Prepare GPU geometry for all objects
    let mut geometry_changed = false;

    for name in &names {
        let color_resolver = ColorResolver::new(
            &session.named_colors,
            &session.element_colors,
        );
        if let Some(mol_obj) = session.registry.get_molecule_mut(name) {
            if mol_obj.is_dirty() {
                geometry_changed = true;
            }
            mol_obj.prepare_render(&gpu.render_context, color_resolver, &session.settings);
        }
    }

    for name in &names {
        if let Some(map_obj) = session.registry.get_map_mut(name) {
            if map_obj.is_dirty() {
                geometry_changed = true;
            }
            map_obj.prepare_render(&gpu.render_context);
        }
    }

    for name in &names {
        if let Some(meas_obj) = session.registry.get_measurement_mut(name) {
            meas_obj.prepare_render(&gpu.render_context);
        }
    }

    if geometry_changed {
        shading.invalidate_shadows();
    }

    // Scene bounding sphere for shadow projection
    if let Some((min, max)) = session.registry.extent() {
        let center = [
            (min.x + max.x) * 0.5,
            (min.y + max.y) * 0.5,
            (min.z + max.z) * 0.5,
        ];
        let dx = max.x - min.x;
        let dy = max.y - min.y;
        let dz = max.z - min.z;
        let radius = (dx * dx + dy * dy + dz * dz).sqrt() * 0.5;
        shading.set_scene_bounds(center, radius.max(1.0));
    }

    // Shadow pre-passes
    let need_shadow_passes = shading.prepare(&mut gpu.render_context, &session.settings);

    // Get surface texture
    let output = gpu.surface.get_current_texture()?;
    let output_view = output
        .texture
        .create_view(&wgpu::TextureViewDescriptor::default());

    // Create command encoder
    let mut encoder = gpu
        .render_context
        .device()
        .create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("Render Encoder"),
        });

    // Shadow depth passes
    if need_shadow_passes {
        render_shadow_passes(&mut encoder, &gpu.render_context, shading, session, &names);
        shading.finish_shadow_passes();
    }

    // 3D scene pass
    {
        let mut render_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
            label: Some("3D Render Pass"),
            color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                view: &output_view,
                resolve_target: None,
                ops: wgpu::Operations {
                    load: wgpu::LoadOp::Clear(wgpu::Color {
                        r: session.clear_color[0] as f64,
                        g: session.clear_color[1] as f64,
                        b: session.clear_color[2] as f64,
                        a: 1.0,
                    }),
                    store: wgpu::StoreOp::Store,
                },
                depth_slice: None,
            })],
            depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                view: &gpu.depth_view,
                depth_ops: Some(wgpu::Operations {
                    load: wgpu::LoadOp::Clear(1.0),
                    store: wgpu::StoreOp::Store,
                }),
                stencil_ops: None,
            }),
            timestamp_writes: None,
            occlusion_query_set: None,
        });

        for name in &names {
            if let Some(mol_obj) = session.registry.get_molecule(name) {
                if mol_obj.is_enabled() {
                    mol_obj.render(&mut render_pass, &gpu.render_context);
                }
            }
        }
        for name in &names {
            if let Some(meas_obj) = session.registry.get_measurement(name) {
                if meas_obj.is_enabled() {
                    meas_obj.render(&mut render_pass, &gpu.render_context);
                }
            }
        }
        for name in &names {
            if let Some(map_obj) = session.registry.get_map(name) {
                if map_obj.is_enabled() {
                    map_obj.render(&mut render_pass, &gpu.render_context);
                }
            }
        }
    }

    // Post-process: silhouette edge detection
    if session.settings.get_bool(pymol_settings::id::silhouettes) {
        let thickness = session.settings.get_float(pymol_settings::id::silhouette_width);
        let depth_jump = session.settings.get_float(pymol_settings::id::silhouette_depth_jump);
        let color_int = session.settings.get_color(pymol_settings::id::silhouette_color);
        let color = pymol_color::Color::from_packed_rgb(color_int).to_rgba(1.0);
        gpu.silhouette.render(
            &mut encoder,
            gpu.render_context.queue(),
            gpu.render_context.device(),
            &output_view,
            &gpu.depth_view,
            width,
            height,
            thickness,
            depth_jump,
            color,
            None,
        );
    }

    // Submit and present
    gpu.render_context
        .queue()
        .submit(std::iter::once(encoder.finish()));
    output.present();

    Ok(())
}

/// Shadow depth passes (mirrors pymol-gui render_shadow_passes).
fn render_shadow_passes(
    encoder: &mut wgpu::CommandEncoder,
    context: &pymol_render::RenderContext,
    shading: &ShadingManager,
    session: &Session,
    names: &[String],
) {
    let state = match shading.shadow_pass_state() {
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
                    load: if i == 0 {
                        wgpu::LoadOp::Clear(1.0)
                    } else {
                        wgpu::LoadOp::Load
                    },
                    store: wgpu::StoreOp::Store,
                }),
                stencil_ops: None,
            }),
            timestamp_writes: None,
            occlusion_query_set: None,
        });

        shadow_pass.set_viewport(
            vp_x as f32,
            vp_y as f32,
            vp_w as f32,
            vp_h as f32,
            0.0,
            1.0,
        );
        shadow_pass.set_scissor_rect(vp_x, vp_y, vp_w, vp_h);

        for name in names {
            if let Some(mol_obj) = session.registry.get_molecule(name) {
                mol_obj.render_shadow_depth(
                    &mut shadow_pass,
                    state.pipelines,
                    &bind_group,
                    context,
                );
            }
        }
    }
}
