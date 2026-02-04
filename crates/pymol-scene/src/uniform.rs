//! Shader uniform setup
//!
//! This module provides utilities for setting up global uniforms used in rendering,
//! including camera matrices, lighting parameters, fog settings, and clip planes.

use crate::camera::Camera;
use pymol_render::GlobalUniforms;
use pymol_settings::GlobalSettings;

/// Compute reflect scale factor to maintain consistent brightness (PyMOL algorithm)
///
/// This scales the `reflect` value inversely based on the light directions to ensure
/// that the total light contribution stays constant regardless of light count.
///
/// From PyMOL's `SceneGetReflectScaleValue` in Scene.cpp:
/// - For each light, adds (1 - z_normalized) where z is the light direction's z component
/// - Multiplies sum by 0.5
/// - Returns 1.0 / sum
pub fn compute_reflect_scale(light_count: i32, light_dirs: &[[f32; 3]]) -> f32 {
    let num_pos_lights = (light_count - 1).max(0) as usize;
    if num_pos_lights > 0 {
        let mut sum = 0.0f32;
        for i in 0..num_pos_lights.min(light_dirs.len()) {
            let dir = light_dirs[i];
            let len = (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]).sqrt();
            if len > 0.0 {
                let z_normalized = dir[2] / len;
                sum += 1.0 - z_normalized;
            }
        }
        sum *= 0.5;
        if sum > 0.0 {
            1.0 / sum
        } else {
            1.0
        }
    } else {
        1.0
    }
}

/// Setup global uniforms for rendering
///
/// This is the single source of truth for uniform setup logic, used by both
/// the `Viewer` and GUI's `App::render()`.
///
/// # Arguments
///
/// * `camera` - Camera for view/projection matrices
/// * `settings` - Global settings for lighting, fog, etc.
/// * `clear_color` - Background color RGB
/// * `viewport_size` - Viewport dimensions (width, height)
///
/// # Example
///
/// ```ignore
/// let uniforms = setup_uniforms(&camera, &settings, clear_color, (width, height));
/// context.update_uniforms(&uniforms);
/// ```
pub fn setup_uniforms(
    camera: &Camera,
    settings: &GlobalSettings,
    clear_color: [f32; 3],
    viewport_size: (f32, f32),
) -> GlobalUniforms {
    let mut uniforms = GlobalUniforms::new();
    uniforms.set_camera(camera.view_matrix(), camera.projection_matrix());
    uniforms.set_background(clear_color);
    uniforms.set_viewport(viewport_size.0, viewport_size.1);

    let camera_pos = camera.world_position();
    uniforms.camera_pos = [camera_pos.x, camera_pos.y, camera_pos.z, 1.0];

    // Multi-light support (PyMOL light_count setting)
    // light_count = 1: ambient only
    // light_count = 2: ambient + 1 directional light
    // light_count = 3-10: ambient + (N-1) directional lights
    let light_count = settings.get_int(pymol_settings::id::light_count);

    // spec_count controls how many positional lights contribute specular
    // -1 (default) means all positional lights contribute specular
    let spec_count = settings.get_int(pymol_settings::id::spec_count);

    // Gather light directions from settings (light, light2, ..., light9)
    let light_setting_ids = [
        pymol_settings::id::light,
        pymol_settings::id::light2,
        pymol_settings::id::light3,
        pymol_settings::id::light4,
        pymol_settings::id::light5,
        pymol_settings::id::light6,
        pymol_settings::id::light7,
        pymol_settings::id::light8,
        pymol_settings::id::light9,
    ];

    let light_dirs: Vec<[f32; 3]> = light_setting_ids
        .iter()
        .map(|&id| settings.get_float3(id))
        .collect();

    uniforms.set_lights(light_count, spec_count, &light_dirs);

    // Lighting settings (PyMOL dual-light model)
    let ambient = settings.get_float(pymol_settings::id::ambient);
    let direct = settings.get_float(pymol_settings::id::direct);
    let reflect = settings.get_float(pymol_settings::id::reflect);
    let specular = settings.get_float(pymol_settings::id::specular);
    let shininess = settings.get_float(pymol_settings::id::shininess);
    let spec_direct = settings.get_float(pymol_settings::id::spec_direct);
    let spec_direct_power = settings.get_float(pymol_settings::id::spec_direct_power);

    // PyMOL brightness consistency adjustments:
    // 1. Scale reflect based on light directions to maintain consistent brightness
    // 2. When light_count < 2, redirect reflect energy to direct (headlight)
    let reflect_scale = compute_reflect_scale(light_count, &light_dirs);
    let reflect_scaled = reflect * reflect_scale;

    let (direct_adjusted, reflect_final) = if light_count < 2 {
        // No positional lights - add reflect to direct (PyMOL behavior)
        ((direct + reflect_scaled).min(1.0), 0.0)
    } else {
        (direct, reflect_scaled)
    };

    uniforms.set_lighting(
        ambient,
        direct_adjusted,
        reflect_final,
        specular,
        shininess,
        spec_direct,
        spec_direct_power,
    );

    // Clip planes from camera
    let current_view = camera.current_view();
    let clip_front = current_view.clip_front;
    let clip_back = current_view.clip_back;
    uniforms.set_clip_planes(clip_front, clip_back);

    // Fog parameters (PyMOL-compatible)
    // depth_cue controls whether fog is enabled
    let depth_cue_enabled = settings.get_bool(pymol_settings::id::depth_cue);
    let fog_density = settings.get_float(pymol_settings::id::fog);

    if depth_cue_enabled && fog_density > 0.0 {
        // fog_start setting is a ratio (0.0-1.0) between clip planes
        let fog_start_ratio = settings.get_float(pymol_settings::id::fog_start);

        // Compute actual fog distances based on clip planes
        // FogStart = (back - front) * fog_start_ratio + front
        let fog_start_actual = (clip_back - clip_front) * fog_start_ratio + clip_front;

        // FogEnd depends on fog density
        let fog_end_actual = if (fog_density - 1.0).abs() < 0.001 {
            // If density is ~1.0, fog ends at back clip plane
            clip_back
        } else {
            // Otherwise scale the fog range by density
            fog_start_actual + (clip_back - fog_start_actual) / fog_density
        };

        // Set fog with background color (so objects fade to background)
        uniforms.set_fog(fog_start_actual, fog_end_actual, fog_density, clear_color);

        // Enable depth cue (darkening based on depth)
        uniforms.set_depth_cue(1.0);
    }
    // else: fog_density and depth_cue remain at 0.0 (disabled) from defaults

    uniforms
}
