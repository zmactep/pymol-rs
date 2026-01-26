// Sphere impostor shader with PyMOL dual-light model
// Renders spheres as billboard quads with ray-sphere intersection in fragment shader
//
// PyMOL Lighting Model:
// - Headlight (direct): Always from camera direction, ensures front-facing surfaces are lit
// - Positional light (reflect): World-space directional light for depth/shadow cues

struct GlobalUniforms {
    view_proj: mat4x4<f32>,
    view: mat4x4<f32>,
    view_inv: mat4x4<f32>,
    proj: mat4x4<f32>,
    camera_pos: vec4<f32>,
    light_dir: vec4<f32>,
    // Headlight (camera light) parameters
    ambient: f32,
    direct: f32,
    spec_direct: f32,
    spec_direct_power: f32,
    // Positional light parameters
    reflect: f32,
    specular: f32,
    shininess: f32,
    _pad0: f32,
    // Fog parameters
    fog_start: f32,
    fog_end: f32,
    fog_density: f32,
    depth_cue: f32,
    fog_color: vec4<f32>,
    bg_color: vec4<f32>,
    viewport: vec4<f32>,
    clip_planes: vec4<f32>,
}

@group(0) @binding(0)
var<uniform> uniforms: GlobalUniforms;

// Billboard vertex (screen-space offset)
struct BillboardVertex {
    @location(10) offset: vec2<f32>,
}

// Instance data (per sphere)
struct SphereInstance {
    @location(0) center: vec3<f32>,
    @location(1) radius: f32,
    @location(2) color: vec4<f32>,
}

struct VertexOutput {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) color: vec4<f32>,
    @location(1) center_view: vec3<f32>,
    @location(2) ray_dir: vec3<f32>,
    @location(3) radius: f32,
}

@vertex
fn vs_main(billboard: BillboardVertex, instance: SphereInstance) -> VertexOutput {
    var output: VertexOutput;
    
    // Transform sphere center to view space
    let center_world = vec4<f32>(instance.center, 1.0);
    let center_view = (uniforms.view * center_world).xyz;
    
    // Create billboard quad in view space
    // Scale by radius and a factor to ensure sphere fits in quad
    let scale = instance.radius * 1.5; // Extra margin for perspective
    let billboard_pos = center_view + vec3<f32>(billboard.offset * scale, 0.0);
    
    // Project to clip space
    output.clip_position = uniforms.proj * vec4<f32>(billboard_pos, 1.0);
    output.color = instance.color;
    output.center_view = center_view;
    output.radius = instance.radius;
    
    // Ray direction in view space (from camera through this vertex)
    // In view space, camera is at origin looking down -Z
    output.ray_dir = normalize(billboard_pos);
    
    return output;
}

// Simplified camera-centric lighting model
// Uses only headlight (from camera) for rotation-independent illumination
fn pymol_lighting(normal: vec3<f32>, view_dir: vec3<f32>, light_dir_view: vec3<f32>, base_color: vec3<f32>) -> vec3<f32> {
    let n = normalize(normal);
    let v = normalize(view_dir);
    
    // Higher ambient for better base illumination
    let ambient_color = base_color * (uniforms.ambient + 0.15);
    
    // === HEADLIGHT (camera light) ===
    // Light comes from camera direction - provides 3D depth while being rotation-independent
    // Combine direct and reflect into a single camera light for consistency
    let headlight_intensity = uniforms.direct + uniforms.reflect;
    let headlight_ndotl = max(dot(n, v), 0.0);
    let headlight_diffuse = base_color * headlight_ndotl * headlight_intensity;
    
    // Headlight specular (Blinn-Phong)
    // For headlight, half vector H = normalize(L + V) = normalize(V + V) = V
    var headlight_specular = vec3<f32>(0.0);
    let spec_intensity = uniforms.spec_direct + uniforms.specular;
    if spec_intensity > 0.0 && headlight_ndotl > 0.0 {
        let spec_power = max(uniforms.spec_direct_power, uniforms.shininess);
        let spec = pow(headlight_ndotl, spec_power);
        headlight_specular = vec3<f32>(1.0) * spec * spec_intensity;
    }
    
    // Combine contributions
    return ambient_color + headlight_diffuse + headlight_specular;
}

// Apply fog
fn apply_fog(color: vec3<f32>, depth: f32) -> vec3<f32> {
    if uniforms.fog_density <= 0.0 {
        return color;
    }
    let fog_factor = clamp((uniforms.fog_end - depth) / (uniforms.fog_end - uniforms.fog_start), 0.0, 1.0);
    return mix(uniforms.fog_color.rgb, color, fog_factor);
}

// Apply depth cue
fn apply_depth_cue(color: vec3<f32>, depth: f32) -> vec3<f32> {
    if uniforms.depth_cue <= 0.0 {
        return color;
    }
    let near = uniforms.clip_planes.x;
    let far = uniforms.clip_planes.y;
    let normalized_depth = clamp((depth - near) / (far - near), 0.0, 1.0);
    let cue = 1.0 - uniforms.depth_cue * normalized_depth * 0.5;
    return color * cue;
}

struct FragmentOutput {
    @location(0) color: vec4<f32>,
    @builtin(frag_depth) depth: f32,
}

@fragment
fn fs_main(input: VertexOutput) -> FragmentOutput {
    var output: FragmentOutput;
    
    // Ray-sphere intersection in view space
    // Ray: P = O + t * D where O is camera (origin in view space)
    // Sphere: |P - C|^2 = r^2
    
    let ray_origin = vec3<f32>(0.0, 0.0, 0.0); // Camera at origin in view space
    let ray_dir = normalize(input.ray_dir);
    let sphere_center = input.center_view;
    let radius = input.radius;
    
    // Solve quadratic: |O + t*D - C|^2 = r^2
    let oc = ray_origin - sphere_center;
    let a = dot(ray_dir, ray_dir);
    let b = 2.0 * dot(oc, ray_dir);
    let c = dot(oc, oc) - radius * radius;
    let discriminant = b * b - 4.0 * a * c;
    
    if discriminant < 0.0 {
        discard;
    }
    
    // Find nearest intersection
    let t = (-b - sqrt(discriminant)) / (2.0 * a);
    
    if t < 0.0 {
        discard;
    }
    
    // Calculate hit point and normal in view space
    let hit_point = ray_origin + t * ray_dir;
    var normal = normalize(hit_point - sphere_center);
    
    // View direction (from hit point to camera, which is at origin)
    let view_dir = normalize(-hit_point);
    
    // Ensure normal faces the camera for proper lighting
    // This handles edge cases where the normal might point away from the viewer
    if dot(normal, view_dir) < 0.0 {
        normal = -normal;
    }
    
    // Transform light direction to view space
    let light_dir_view = (uniforms.view * vec4<f32>(uniforms.light_dir.xyz, 0.0)).xyz;
    
    // Calculate lighting using PyMOL dual-light model
    var color = pymol_lighting(normal, view_dir, light_dir_view, input.color.rgb);
    
    // Calculate depth
    let depth = -hit_point.z;
    
    // Apply depth cue and fog
    color = apply_depth_cue(color, depth);
    color = apply_fog(color, depth);
    
    // Calculate clip-space depth
    let clip_pos = uniforms.proj * vec4<f32>(hit_point, 1.0);
    output.depth = clip_pos.z / clip_pos.w;
    
    output.color = vec4<f32>(color, input.color.a);
    
    return output;
}
