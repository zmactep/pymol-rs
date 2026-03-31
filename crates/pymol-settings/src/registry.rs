//! Setting descriptor registry — O(1) lookup by name.
//!
//! Provides the bridge between string-based `set`/`get` commands and the typed
//! `Settings` / `ObjectOverrides` structs. Each setting has a descriptor with
//! function pointers for reading/writing values.

use std::sync::LazyLock;

use ahash::AHashMap;

use crate::error::SettingError;
use crate::macros::SettingDescriptor;
use crate::setting::{SettingType, SettingValue};
use crate::side_effects::SideEffectCategory;

/// Build the full list of setting descriptors with real get/set function pointers.
fn build_descriptors() -> Vec<SettingDescriptor> {
    let mut descs = Vec::new();

    // Helper macro to reduce boilerplate for registering settings.
    // `global_path` is the path from Settings to the field (e.g., shading.classic.ambient).
    // `override_path` (optional) is the path from ObjectOverrides to the field.
    macro_rules! reg {
        // Global-only setting (no overrides)
        (
            $name:expr, $ty:ty,
            global = |s| s.$($gpath:ident).+
            $(, min = $min:expr, max = $max:expr)?
            $(, hints = $hints_ty:ty)?
            $(, side_effects = [$($se:ident),*])?
        ) => {
            descs.push(SettingDescriptor {
                name: $name,
                setting_type: <$ty as crate::macros::SettingKind>::setting_type(),
                get: |s| crate::macros::SettingKind::wrap(&s.$($gpath).+),
                set: |s, v| {
                    s.$($gpath).+ = <$ty as crate::macros::SettingKind>::unwrap_value(v)?;
                    Ok(())
                },
                get_override: None,
                set_override: None,
                unset_override: None,
                min: reg!(@opt_f32 $($min)?),
                max: reg!(@opt_f32 $($max)?),
                value_hints: reg!(@hints $($hints_ty)?),
                side_effects: &[$($($crate::SideEffectCategory::$se),*)?],
            });
        };
        // Object-overridable setting (with overrides)
        (
            $name:expr, $ty:ty,
            global = |s| s.$($gpath:ident).+,
            override = |o| o.$($opath:ident).+
            $(, min = $min:expr, max = $max:expr)?
            $(, hints = $hints_ty:ty)?
            $(, side_effects = [$($se:ident),*])?
        ) => {
            descs.push(SettingDescriptor {
                name: $name,
                setting_type: <$ty as crate::macros::SettingKind>::setting_type(),
                get: |s| crate::macros::SettingKind::wrap(&s.$($gpath).+),
                set: |s, v| {
                    s.$($gpath).+ = <$ty as crate::macros::SettingKind>::unwrap_value(v)?;
                    Ok(())
                },
                get_override: Some(|o| {
                    o.$($opath).+.as_ref().map(|v| crate::macros::SettingKind::wrap(v))
                }),
                set_override: Some(|o, v| {
                    o.$($opath).+ = Some(<$ty as crate::macros::SettingKind>::unwrap_value(v)?);
                    Ok(())
                }),
                unset_override: Some(|o| { o.$($opath).+ = None; }),
                min: reg!(@opt_f32 $($min)?),
                max: reg!(@opt_f32 $($max)?),
                value_hints: reg!(@hints $($hints_ty)?),
                side_effects: &[$($($crate::SideEffectCategory::$se),*)?],
            });
        };
        (@opt_f32) => { None };
        (@opt_f32 $val:expr) => { Some($val as f32) };
        (@hints) => { &[] };
        (@hints $ty:ty) => { <$ty as crate::SettingEnum>::value_hints() };
    }

    // =========================================================================
    // Shading — Common
    // =========================================================================
    reg!("silhouettes", bool, global = |s| s.shading.common.silhouettes, side_effects = [ShaderReload, SceneInvalidate]);
    reg!("silhouette_width", f32, global = |s| s.shading.common.silhouette_width, min = 0.1, max = 10.0, side_effects = [SceneInvalidate]);
    reg!("silhouette_color", i32, global = |s| s.shading.common.silhouette_color, side_effects = [SceneInvalidate]);
    reg!("silhouette_depth_jump", f32, global = |s| s.shading.common.silhouette_depth_jump, min = 0.0, max = 1.0, side_effects = [SceneInvalidate]);
    reg!("depth_cue", bool, global = |s| s.shading.common.depth_cue, side_effects = [ShaderReload, SceneInvalidate]);
    reg!("fog", f32, global = |s| s.shading.common.fog, min = 0.0, max = 1.0, side_effects = [SceneInvalidate]);
    reg!("fog_start", f32, global = |s| s.shading.common.fog_start, min = 0.0, max = 1.0, side_effects = [SceneInvalidate]);

    // =========================================================================
    // Shading — Mode
    // =========================================================================
    reg!("shading_mode", crate::ShadingMode, global = |s| s.shading.mode, hints = crate::ShadingMode, side_effects = [ShaderReload, ShaderComputeLighting, SceneInvalidate]);

    // =========================================================================
    // Shading — Classic
    // =========================================================================
    reg!("ambient", f32, global = |s| s.shading.classic.ambient, min = 0.0, max = 1.0, side_effects = [ShaderComputeLighting, SceneInvalidate]);
    reg!("direct", f32, global = |s| s.shading.classic.direct, min = 0.0, max = 1.0, side_effects = [ShaderComputeLighting, SceneInvalidate]);
    reg!("reflect", f32, global = |s| s.shading.classic.reflect, min = 0.0, max = 1.0, side_effects = [ShaderComputeLighting, SceneInvalidate]);
    reg!("specular", f32, global = |s| s.shading.classic.specular, min = 0.0, max = 1.0, side_effects = [ShaderComputeLighting, SceneInvalidate]);
    reg!("shininess", f32, global = |s| s.shading.classic.shininess, min = 0.0, max = 1000.0, side_effects = [ShaderComputeLighting, SceneInvalidate]);
    reg!("spec_direct", f32, global = |s| s.shading.classic.spec_direct, min = 0.0, max = 1.0, side_effects = [ShaderComputeLighting, SceneInvalidate]);
    reg!("spec_direct_power", f32, global = |s| s.shading.classic.spec_direct_power, min = 0.0, max = 1000.0, side_effects = [ShaderComputeLighting, SceneInvalidate]);
    reg!("spec_count", i32, global = |s| s.shading.classic.spec_count, side_effects = [ShaderComputeLighting, SceneInvalidate]);
    reg!("light_count", i32, global = |s| s.shading.classic.light_count, min = 1.0, max = 9.0, side_effects = [ShaderComputeLighting, SceneInvalidate]);

    // Lights — [f32; 3] doesn't implement SettingKind the same way; register manually
    macro_rules! reg_light {
        ($name:expr, $field:ident) => {
            descs.push(SettingDescriptor {
                name: $name,
                setting_type: SettingType::Float3,
                get: |s| SettingValue::Float3(s.shading.classic.$field),
                set: |s, v| {
                    s.shading.classic.$field = v.as_float3().ok_or(SettingError::type_mismatch("float3"))?;
                    Ok(())
                },
                get_override: None,
                set_override: None,
                unset_override: None,
                min: None,
                max: None,
                value_hints: &[],
                side_effects: &[SideEffectCategory::ShaderComputeLighting, SideEffectCategory::SceneInvalidate],
            });
        };
    }
    reg_light!("light", light);
    reg_light!("light2", light2);
    reg_light!("light3", light3);
    reg_light!("light4", light4);
    reg_light!("light5", light5);
    reg_light!("light6", light6);
    reg_light!("light7", light7);
    reg_light!("light8", light8);
    reg_light!("light9", light9);

    // =========================================================================
    // Shading — Skripkin
    // =========================================================================
    reg!("skripkin_directions", i32, global = |s| s.shading.skripkin.directions, min = 1.0, max = 256.0, side_effects = [ShaderReload, SceneInvalidate]);
    reg!("skripkin_map_size", i32, global = |s| s.shading.skripkin.map_size, min = 32.0, max = 4096.0, side_effects = [ShaderReload, SceneInvalidate]);
    reg!("skripkin_bias", f32, global = |s| s.shading.skripkin.bias, min = 0.0, max = 1.0, side_effects = [SceneInvalidate]);
    reg!("skripkin_intensity", f32, global = |s| s.shading.skripkin.intensity, min = 0.0, max = 2.0, side_effects = [SceneInvalidate]);
    reg!("skripkin_ambient", f32, global = |s| s.shading.skripkin.ambient, min = 0.0, max = 1.0, side_effects = [ShaderComputeLighting, SceneInvalidate]);
    reg!("skripkin_specular", f32, global = |s| s.shading.skripkin.specular, min = 0.0, max = 1.0, side_effects = [ShaderComputeLighting, SceneInvalidate]);
    reg!("skripkin_shininess", f32, global = |s| s.shading.skripkin.shininess, min = 0.0, max = 1000.0, side_effects = [ShaderComputeLighting, SceneInvalidate]);

    // =========================================================================
    // Shading — Full shadows
    // =========================================================================
    reg!("shadow_map_size", i32, global = |s| s.shading.full.shadow_map_size, min = 64.0, max = 4096.0, side_effects = [ShaderReload, SceneInvalidate]);
    reg!("shadow_bias", f32, global = |s| s.shading.full.shadow_bias, min = 0.0, max = 1.0, side_effects = [SceneInvalidate]);
    reg!("shadow_intensity", f32, global = |s| s.shading.full.shadow_intensity, min = 0.0, max = 1.0, side_effects = [SceneInvalidate]);
    reg!("shadow_pcf", i32, global = |s| s.shading.full.shadow_pcf, min = 1.0, max = 16.0, side_effects = [SceneInvalidate]);

    // =========================================================================
    // UI
    // =========================================================================
    reg!("mouse_selection_mode", crate::MouseSelectionMode, global = |s| s.ui.mouse_selection_mode, hints = crate::MouseSelectionMode, side_effects = [OrthoDirty]);
    reg!("mouse_wheel_scale", f32, global = |s| s.ui.mouse_wheel_scale, min = 0.01, max = 10.0);
    reg!("antialias", i32, global = |s| s.ui.antialias, side_effects = [SceneInvalidate]);
    reg!("selection_width", f32, global = |s| s.ui.selection_width, min = 0.5, max = 20.0, side_effects = [SceneInvalidate]);
    reg!("transparent_panels", bool, global = |s| s.ui.transparent_panels, side_effects = [OrthoDirty]);
    reg!("opaque_background", bool, global = |s| s.ui.opaque_background, side_effects = [SceneInvalidate]);
    reg!("orthoscopic", bool, global = |s| s.ui.orthoscopic, side_effects = [SceneInvalidate]);
    reg!("bg_rgb", i32, global = |s| s.ui.bg_rgb, side_effects = [ViewportUpdate]);
    reg!("bg_rgb_top", i32, global = |s| s.ui.bg_rgb_top, side_effects = [ViewportUpdate]);
    reg!("bg_rgb_bottom", i32, global = |s| s.ui.bg_rgb_bottom, side_effects = [ViewportUpdate]);

    // =========================================================================
    // Movie
    // =========================================================================
    reg!("movie_loop", bool, global = |s| s.movie.movie_loop);
    reg!("movie_fps", f32, global = |s| s.movie.movie_fps, min = 1.0, max = 120.0);
    reg!("movie_auto_interpolate", bool, global = |s| s.movie.movie_auto_interpolate);

    // =========================================================================
    // Behavior
    // =========================================================================
    reg!("ignore_case", bool, global = |s| s.behavior.ignore_case);
    reg!("ignore_case_chain", bool, global = |s| s.behavior.ignore_case_chain);
    reg!("auto_dss", bool, global = |s| s.behavior.auto_dss);
    reg!("bonding_vdw_cutoff", f32, global = |s| s.behavior.bonding_vdw_cutoff, min = 0.0, max = 1.0);

    // =========================================================================
    // Cartoon (object-overridable)
    // =========================================================================
    reg!("cartoon_sampling", i32, global = |s| s.cartoon.sampling, override = |o| o.cartoon.sampling, side_effects = [RepresentationRebuild]);
    reg!("cartoon_power", f32, global = |s| s.cartoon.power, override = |o| o.cartoon.power, min = 0.0, max = 10.0, side_effects = [RepresentationRebuild]);
    reg!("cartoon_power_b", f32, global = |s| s.cartoon.power_b, override = |o| o.cartoon.power_b, min = 0.0, max = 10.0, side_effects = [RepresentationRebuild]);
    reg!("cartoon_throw", f32, global = |s| s.cartoon.throw, override = |o| o.cartoon.throw, min = 0.0, max = 10.0, side_effects = [RepresentationRebuild]);
    reg!("cartoon_smooth_first", i32, global = |s| s.cartoon.smooth_first, override = |o| o.cartoon.smooth_first, side_effects = [RepresentationRebuild]);
    reg!("cartoon_smooth_last", i32, global = |s| s.cartoon.smooth_last, override = |o| o.cartoon.smooth_last, side_effects = [RepresentationRebuild]);
    reg!("cartoon_smooth_cycles", i32, global = |s| s.cartoon.smooth_cycles, override = |o| o.cartoon.smooth_cycles, side_effects = [RepresentationRebuild]);
    reg!("cartoon_smooth_loops", bool, global = |s| s.cartoon.smooth_loops, override = |o| o.cartoon.smooth_loops, side_effects = [RepresentationRebuild]);
    reg!("cartoon_flat_cycles", i32, global = |s| s.cartoon.flat_cycles, override = |o| o.cartoon.flat_cycles, side_effects = [RepresentationRebuild]);
    reg!("cartoon_gap_cutoff", i32, global = |s| s.cartoon.gap_cutoff, override = |o| o.cartoon.gap_cutoff, side_effects = [RepresentationRebuild]);
    reg!("cartoon_refine", i32, global = |s| s.cartoon.refine, override = |o| o.cartoon.refine, side_effects = [RepresentationRebuild]);
    reg!("cartoon_refine_normals", i32, global = |s| s.cartoon.refine_normals, override = |o| o.cartoon.refine_normals, side_effects = [RepresentationRebuild]);
    reg!("cartoon_round_helices", bool, global = |s| s.cartoon.round_helices, override = |o| o.cartoon.round_helices, side_effects = [RepresentationRebuild]);
    reg!("cartoon_color", i32, global = |s| s.cartoon.color, override = |o| o.cartoon.color, side_effects = [ColorRebuild]);
    reg!("cartoon_nucleic_acid_color", i32, global = |s| s.cartoon.nucleic_acid_color, override = |o| o.cartoon.nucleic_acid_color, side_effects = [ColorRebuild]);

    // =========================================================================
    // Stick (object-overridable)
    // =========================================================================
    reg!("stick_radius", f32, global = |s| s.stick.radius, override = |o| o.stick.radius, min = 0.01, max = 5.0, side_effects = [RepresentationRebuild]);
    reg!("stick_color", i32, global = |s| s.stick.color, override = |o| o.stick.color, side_effects = [ColorRebuild]);
    reg!("stick_valence_scale", f32, global = |s| s.stick.valence_scale, override = |o| o.stick.valence_scale, min = 0.0, max = 5.0, side_effects = [RepresentationRebuild]);
    reg!("valence", bool, global = |s| s.stick.valence, override = |o| o.stick.valence, side_effects = [RepresentationRebuild]);
    reg!("valence_size", f32, global = |s| s.stick.valence_size, override = |o| o.stick.valence_size, min = 0.0, max = 1.0, side_effects = [RepresentationRebuild]);

    // =========================================================================
    // Sphere (object-overridable)
    // =========================================================================
    reg!("sphere_scale", f32, global = |s| s.sphere.scale, override = |o| o.sphere.scale, min = 0.01, max = 10.0, side_effects = [RepresentationRebuild]);
    reg!("sphere_quality", i32, global = |s| s.sphere.quality, override = |o| o.sphere.quality, min = 0.0, max = 4.0, side_effects = [RepresentationRebuild]);
    reg!("sphere_color", i32, global = |s| s.sphere.color, override = |o| o.sphere.color, side_effects = [ColorRebuild]);

    // =========================================================================
    // Surface (object-overridable)
    // =========================================================================
    reg!("surface_type", i32, global = |s| s.surface.surface_type, override = |o| o.surface.surface_type, side_effects = [RepresentationRebuild]);
    reg!("surface_solvent", bool, global = |s| s.surface.solvent, override = |o| o.surface.solvent, side_effects = [RepresentationRebuild]);
    reg!("surface_quality", i32, global = |s| s.surface.quality, override = |o| o.surface.quality, side_effects = [RepresentationRebuild]);
    reg!("surface_color", i32, global = |s| s.surface.color, override = |o| o.surface.color, side_effects = [ColorRebuild]);
    reg!("surface_individual_chains", bool, global = |s| s.surface.individual_chains, override = |o| o.surface.individual_chains, side_effects = [RepresentationRebuild]);
    reg!("transparency", f32, global = |s| s.surface.transparency, override = |o| o.surface.transparency, min = 0.0, max = 1.0, side_effects = [SceneInvalidate]);
    reg!("transparency_mode", i32, global = |s| s.surface.transparency_mode, override = |o| o.surface.transparency_mode, side_effects = [SceneInvalidate]);

    // =========================================================================
    // Ribbon (object-overridable)
    // =========================================================================
    reg!("ribbon_sampling", i32, global = |s| s.ribbon.sampling, override = |o| o.ribbon.sampling, side_effects = [RepresentationRebuild]);
    reg!("ribbon_power", f32, global = |s| s.ribbon.power, override = |o| o.ribbon.power, min = 0.0, max = 10.0, side_effects = [RepresentationRebuild]);
    reg!("ribbon_power_b", f32, global = |s| s.ribbon.power_b, override = |o| o.ribbon.power_b, min = 0.0, max = 10.0, side_effects = [RepresentationRebuild]);
    reg!("ribbon_throw", f32, global = |s| s.ribbon.throw, override = |o| o.ribbon.throw, min = 0.0, max = 10.0, side_effects = [RepresentationRebuild]);
    reg!("ribbon_radius", f32, global = |s| s.ribbon.radius, override = |o| o.ribbon.radius, min = 0.0, max = 5.0, side_effects = [RepresentationRebuild]);
    reg!("ribbon_color", i32, global = |s| s.ribbon.color, override = |o| o.ribbon.color, side_effects = [ColorRebuild]);

    // =========================================================================
    // Line (object-overridable)
    // =========================================================================
    reg!("line_color", i32, global = |s| s.line.color, override = |o| o.line.color, side_effects = [ColorRebuild]);

    // =========================================================================
    // Dot (object-overridable)
    // =========================================================================
    reg!("dot_density", i32, global = |s| s.dot.density, override = |o| o.dot.density, min = 0.0, max = 4.0, side_effects = [RepresentationRebuild]);
    reg!("dot_width", f32, global = |s| s.dot.width, override = |o| o.dot.width, min = 0.5, max = 20.0, side_effects = [SceneInvalidate]);
    reg!("dot_lighting", bool, global = |s| s.dot.lighting, override = |o| o.dot.lighting, side_effects = [SceneInvalidate]);

    // =========================================================================
    // Mesh (object-overridable)
    // =========================================================================
    reg!("mesh_color", i32, global = |s| s.mesh.color, override = |o| o.mesh.color, side_effects = [ColorRebuild]);
    reg!("mesh_lighting", bool, global = |s| s.mesh.lighting, override = |o| o.mesh.lighting, side_effects = [SceneInvalidate]);

    descs
}

/// All setting descriptors, built once at startup.
static DESCRIPTORS: LazyLock<Vec<SettingDescriptor>> = LazyLock::new(build_descriptors);

/// Name → index into DESCRIPTORS.
static NAME_INDEX: LazyLock<AHashMap<&'static str, usize>> = LazyLock::new(|| {
    DESCRIPTORS
        .iter()
        .enumerate()
        .map(|(i, d)| (d.name, i))
        .collect()
});

/// Look up a descriptor by string name. O(1).
pub fn lookup_by_name(name: &str) -> Option<&'static SettingDescriptor> {
    NAME_INDEX.get(name).map(|&i| &DESCRIPTORS[i])
}

/// Iterator over all descriptors.
pub fn all_descriptors() -> &'static [SettingDescriptor] {
    &DESCRIPTORS
}

/// Total number of registered settings.
pub fn descriptor_count() -> usize {
    DESCRIPTORS.len()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::groups::Settings;
    use crate::overrides::ObjectOverrides;

    #[test]
    fn test_descriptor_count() {
        // We registered ~94 settings
        let count = descriptor_count();
        assert!(count > 80, "expected at least 80 descriptors, got {count}");
        assert!(count < 120, "expected fewer than 120 descriptors, got {count}");
    }

    #[test]
    fn test_lookup_by_name() {
        let d = lookup_by_name("ambient").expect("ambient not found");
        assert_eq!(d.setting_type, SettingType::Float);
    }

    #[test]
    fn test_get_set_roundtrip() {
        let mut settings = Settings::default();
        let d = lookup_by_name("ambient").unwrap();

        // Read default
        let val = (d.get)(&settings);
        assert_eq!(val, SettingValue::Float(0.14));

        // Write new value
        (d.set)(&mut settings, SettingValue::Float(0.3)).unwrap();
        let val = (d.get)(&settings);
        assert_eq!(val, SettingValue::Float(0.3));
    }

    #[test]
    fn test_override_get_set() {
        let mut overrides = ObjectOverrides::default();
        let d = lookup_by_name("stick_radius").unwrap();

        // Override starts as None
        let get_ov = d.get_override.expect("stick_radius should be overridable");
        assert!(get_ov(&overrides).is_none());

        // Set override
        let set_ov = d.set_override.unwrap();
        set_ov(&mut overrides, SettingValue::Float(0.5)).unwrap();
        assert_eq!(get_ov(&overrides), Some(SettingValue::Float(0.5)));

        // Unset override
        let unset_ov = d.unset_override.unwrap();
        unset_ov(&mut overrides);
        assert!(get_ov(&overrides).is_none());
    }

    #[test]
    fn test_global_only_has_no_override() {
        let d = lookup_by_name("ambient").unwrap();
        assert!(d.get_override.is_none());
        assert!(d.set_override.is_none());
        assert!(d.unset_override.is_none());
    }

    #[test]
    fn test_enum_hints() {
        let d = lookup_by_name("shading_mode").unwrap();
        assert_eq!(d.value_hints.len(), 3);
        assert_eq!(d.value_hints[0].0, "classic");

        let d = lookup_by_name("mouse_selection_mode").unwrap();
        assert_eq!(d.value_hints.len(), 7);
    }

    #[test]
    fn test_enum_get_set() {
        let mut settings = Settings::default();
        let d = lookup_by_name("shading_mode").unwrap();

        let val = (d.get)(&settings);
        assert_eq!(val, SettingValue::Int(0)); // Classic

        (d.set)(&mut settings, SettingValue::Int(1)).unwrap(); // Skripkin
        assert_eq!(settings.shading.mode, crate::ShadingMode::Skripkin);
    }

    #[test]
    fn test_no_duplicate_names() {
        let names: Vec<&str> = DESCRIPTORS.iter().map(|d| d.name).collect();
        let mut seen = std::collections::HashSet::new();
        for name in &names {
            assert!(seen.insert(name), "duplicate descriptor name: {name}");
        }
    }

    #[test]
    fn test_light_float3() {
        let mut settings = Settings::default();
        let d = lookup_by_name("light").unwrap();
        assert_eq!(d.setting_type, SettingType::Float3);

        let val = (d.get)(&settings);
        assert_eq!(val, SettingValue::Float3([-0.4, -0.4, -1.0]));

        (d.set)(&mut settings, SettingValue::Float3([1.0, 0.0, 0.0])).unwrap();
        assert_eq!(settings.shading.classic.light, [1.0, 0.0, 0.0]);
    }
}
