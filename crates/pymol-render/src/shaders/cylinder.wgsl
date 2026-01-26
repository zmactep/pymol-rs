// Cylinder impostor shader with PyMOL dual-light model for stick/bond rendering
// Renders cylinders as oriented quads with ray-cylinder intersection
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

// Billboard vertex
struct BillboardVertex {
    @location(10) offset: vec2<f32>,
}

// Instance data (per cylinder)
struct CylinderInstance {
    @location(0) start: vec3<f32>,
    @location(1) radius: f32,
    @location(2) end: vec3<f32>,
    @location(3) flags: u32,
    @location(4) color1: vec4<f32>,
    @location(5) color2: vec4<f32>,
}

struct VertexOutput {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) color1: vec4<f32>,
    @location(1) color2: vec4<f32>,
    @location(2) start_view: vec3<f32>,
    @location(3) end_view: vec3<f32>,
    @location(4) ray_origin: vec3<f32>,
    @location(5) radius: f32,
    @location(6) flags: u32,
}

@vertex
fn vs_main(billboard: BillboardVertex, instance: CylinderInstance) -> VertexOutput {
    var output: VertexOutput;
    
    // Transform cylinder endpoints to view space
    let start_view = (uniforms.view * vec4<f32>(instance.start, 1.0)).xyz;
    let end_view = (uniforms.view * vec4<f32>(instance.end, 1.0)).xyz;
    
    // Cylinder center in view space
    let center_view = (start_view + end_view) * 0.5;
    
    // Calculate the bounding sphere radius for the cylinder
    // This ensures the billboard covers the cylinder from any angle
    let axis = end_view - start_view;
    let axis_len = length(axis);
    let half_len = axis_len * 0.5;
    
    // Bounding sphere radius = sqrt(half_length^2 + radius^2) + margin
    let bound_radius = sqrt(half_len * half_len + instance.radius * instance.radius) + instance.radius;
    
    // Create camera-facing billboard quad
    // In view space, X is right, Y is up, Z is toward camera
    // Simply offset in X and Y from the center
    let scale = bound_radius * 1.5; // Extra margin for ray intersection
    let billboard_pos = center_view + vec3<f32>(billboard.offset.x * scale, billboard.offset.y * scale, 0.0);
    
    output.clip_position = uniforms.proj * vec4<f32>(billboard_pos, 1.0);
    output.color1 = instance.color1;
    output.color2 = instance.color2;
    output.start_view = start_view;
    output.end_view = end_view;
    output.ray_origin = billboard_pos;
    output.radius = instance.radius;
    output.flags = instance.flags;
    
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

// Ray-cylinder intersection
fn ray_cylinder_intersect(
    ray_origin: vec3<f32>,
    ray_dir: vec3<f32>,
    cyl_start: vec3<f32>,
    cyl_end: vec3<f32>,
    radius: f32
) -> vec4<f32> {
    // Returns vec4(t, u, nx, ny) where:
    // t = ray parameter, u = position along cylinder [0,1]
    // (nx, ny) can be used to reconstruct normal
    // Returns t < 0 if no hit
    
    let axis = cyl_end - cyl_start;
    let axis_len = length(axis);
    let axis_dir = axis / axis_len;
    
    let oc = ray_origin - cyl_start;
    
    // Project ray and oc onto plane perpendicular to cylinder axis
    let ray_proj = ray_dir - dot(ray_dir, axis_dir) * axis_dir;
    let oc_proj = oc - dot(oc, axis_dir) * axis_dir;
    
    let a = dot(ray_proj, ray_proj);
    let b = 2.0 * dot(ray_proj, oc_proj);
    let c = dot(oc_proj, oc_proj) - radius * radius;
    
    let discriminant = b * b - 4.0 * a * c;
    
    if discriminant < 0.0 {
        return vec4<f32>(-1.0, 0.0, 0.0, 0.0);
    }
    
    let t = (-b - sqrt(discriminant)) / (2.0 * a);
    
    if t < 0.0 {
        return vec4<f32>(-1.0, 0.0, 0.0, 0.0);
    }
    
    // Check if hit is within cylinder bounds
    let hit_point = ray_origin + t * ray_dir;
    let hit_proj = dot(hit_point - cyl_start, axis_dir);
    let u = hit_proj / axis_len;
    
    if u < 0.0 || u > 1.0 {
        return vec4<f32>(-1.0, 0.0, 0.0, 0.0);
    }
    
    return vec4<f32>(t, u, 0.0, 0.0);
}

struct FragmentOutput {
    @location(0) color: vec4<f32>,
    @builtin(frag_depth) depth: f32,
}

@fragment
fn fs_main(input: VertexOutput) -> FragmentOutput {
    var output: FragmentOutput;
    
    // Ray from camera through this fragment
    let ray_origin = vec3<f32>(0.0, 0.0, 0.0);
    let ray_dir = normalize(input.ray_origin);
    
    // Ray-cylinder intersection
    let hit = ray_cylinder_intersect(
        ray_origin,
        ray_dir,
        input.start_view,
        input.end_view,
        input.radius
    );
    
    if hit.x < 0.0 {
        discard;
    }
    
    let t = hit.x;
    let u = hit.y; // Position along cylinder [0, 1]
    
    let hit_point = ray_origin + t * ray_dir;
    
    // Calculate normal (perpendicular to cylinder axis, pointing outward)
    let axis = input.end_view - input.start_view;
    let axis_dir = normalize(axis);
    let center_on_axis = input.start_view + u * axis;
    var normal = normalize(hit_point - center_on_axis);
    
    // View direction (from hit point to camera, which is at origin)
    let view_dir = normalize(-hit_point);
    
    // Ensure normal faces the camera for proper lighting
    // This handles edge cases where the normal might point away from the viewer
    if dot(normal, view_dir) < 0.0 {
        normal = -normal;
    }
    
    // Interpolate color based on position along cylinder
    let base_color = mix(input.color1.rgb, input.color2.rgb, u);
    let alpha = mix(input.color1.a, input.color2.a, u);
    
    // Transform light direction to view space
    let light_dir_view = (uniforms.view * vec4<f32>(uniforms.light_dir.xyz, 0.0)).xyz;
    
    // Calculate lighting using PyMOL dual-light model
    var color = pymol_lighting(normal, view_dir, light_dir_view, base_color);
    
    // Calculate depth
    let depth = -hit_point.z;
    
    // Apply depth cue and fog
    color = apply_depth_cue(color, depth);
    color = apply_fog(color, depth);
    
    // Calculate clip-space depth
    let clip_pos = uniforms.proj * vec4<f32>(hit_point, 1.0);
    output.depth = clip_pos.z / clip_pos.w;
    
    output.color = vec4<f32>(color, alpha);
    
    return output;
}
