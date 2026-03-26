// PyMOL multi-light lighting model.
// Used by sphere, cylinder, and mesh shaders (geometry with normals).
// Requires common.wgsl to be prepended.

// PyMOL multi-light model
// - Headlight (from camera direction) controlled by 'direct' setting
// - Positional lights controlled by light_count and 'reflect'/'specular' settings
// - spec_count controls how many positional lights contribute specular (-1 = all)
fn pymol_lighting(normal: vec3<f32>, view_dir: vec3<f32>, base_color: vec3<f32>) -> vec3<f32> {
    let n = normalize(normal);
    let v = normalize(view_dir);

    var color = base_color * uniforms.ambient;

    // Headlight (camera light)
    let headlight_ndotl = max(dot(n, v), 0.0);
    color += base_color * headlight_ndotl * uniforms.direct;

    if uniforms.spec_direct > 0.0 && headlight_ndotl > 0.0 {
        let spec = pow(headlight_ndotl, uniforms.spec_direct_power);
        color += vec3<f32>(1.0) * spec * uniforms.spec_direct;
    }

    // Positional lights
    let num_pos_lights = max(uniforms.light_count - 1, 0);
    let effective_spec_count = select(uniforms.spec_count, num_pos_lights, uniforms.spec_count < 0);

    for (var i = 0; i < num_pos_lights; i++) {
        let l = normalize(-uniforms.light_dirs[i].xyz);
        let ndotl = max(dot(n, l), 0.0);
        color += base_color * ndotl * uniforms.reflect;

        if uniforms.specular > 0.0 && ndotl > 0.0 && i < effective_spec_count {
            let h = normalize(l + v);
            let ndoth = max(dot(n, h), 0.0);
            let spec = pow(ndoth, uniforms.shininess);
            color += vec3<f32>(1.0) * spec * uniforms.specular;
        }
    }

    return color;
}

// PyMOL multi-light model with per-light directional shadow sampling
fn pymol_lighting_shadowed(normal: vec3<f32>, view_dir: vec3<f32>, base_color: vec3<f32>, world_pos: vec3<f32>) -> vec3<f32> {
    let n = normalize(normal);
    let v = normalize(view_dir);

    var color = base_color * uniforms.ambient;

    // Headlight — no shadow (camera light)
    let headlight_ndotl = max(dot(n, v), 0.0);
    color += base_color * headlight_ndotl * uniforms.direct;
    if uniforms.spec_direct > 0.0 && headlight_ndotl > 0.0 {
        let spec = pow(headlight_ndotl, uniforms.spec_direct_power);
        color += vec3<f32>(1.0) * spec * uniforms.spec_direct;
    }

    // Positional lights with per-light shadow
    let num_pos_lights = max(uniforms.light_count - 1, 0);
    let effective_spec_count = select(uniforms.spec_count, num_pos_lights, uniforms.spec_count < 0);

    for (var i = 0; i < num_pos_lights; i++) {
        var shadow = 1.0;
        if u32(i) < shadow_params.shadow_count {
            shadow = sample_light_shadow(world_pos, u32(i));
        }

        let l = normalize(-uniforms.light_dirs[i].xyz);
        let ndotl = max(dot(n, l), 0.0);
        color += base_color * ndotl * uniforms.reflect * shadow;

        if uniforms.specular > 0.0 && ndotl > 0.0 && i < effective_spec_count {
            let h = normalize(l + v);
            let ndoth = max(dot(n, h), 0.0);
            let spec = pow(ndoth, uniforms.shininess);
            color += vec3<f32>(1.0) * spec * uniforms.specular * shadow;
        }
    }

    return color;
}
