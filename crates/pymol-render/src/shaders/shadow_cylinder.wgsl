// Depth-only cylinder impostor shader for shadow map rendering
// Simplified ray-cylinder intersection from the light's viewpoint.
// Writes corrected depth — no color output.

struct ShadowUniform {
    view_proj: mat4x4<f32>,
}

@group(0) @binding(0)
var<uniform> shadow: ShadowUniform;

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
    @location(0) start_world: vec3<f32>,
    @location(1) end_world: vec3<f32>,
    @location(2) radius: f32,
}

@vertex
fn vs_main(billboard: BillboardVertex, instance: CylinderInstance) -> VertexOutput {
    var output: VertexOutput;

    // Cylinder center in world space
    let center = (instance.start + instance.end) * 0.5;
    let axis = instance.end - instance.start;
    let axis_len = length(axis);
    let half_len = axis_len * 0.5;

    // Bounding sphere radius
    let bound_radius = sqrt(half_len * half_len + instance.radius * instance.radius) + instance.radius;

    // Transform center to clip space
    let center_clip = shadow.view_proj * vec4<f32>(center, 1.0);

    // Scale billboard offset in clip space
    let sx = length(vec3<f32>(shadow.view_proj[0][0], shadow.view_proj[0][1], shadow.view_proj[0][2]));
    let sy = length(vec3<f32>(shadow.view_proj[1][0], shadow.view_proj[1][1], shadow.view_proj[1][2]));
    let scale_x = bound_radius * sx * 1.5;
    let scale_y = bound_radius * sy * 1.5;

    output.clip_position = center_clip + vec4<f32>(
        billboard.offset.x * scale_x,
        billboard.offset.y * scale_y,
        0.0, 0.0
    );

    output.start_world = instance.start;
    output.end_world = instance.end;
    output.radius = instance.radius;

    return output;
}

struct FragmentOutput {
    @builtin(frag_depth) depth: f32,
}

@fragment
fn fs_main(input: VertexOutput) -> FragmentOutput {
    var output: FragmentOutput;

    // For orthographic shadow maps, approximate the cylinder depth
    // by projecting the closest point on the cylinder axis to the light.

    // Light forward direction from the matrix
    let fwd = normalize(vec3<f32>(
        shadow.view_proj[0][2],
        shadow.view_proj[1][2],
        shadow.view_proj[2][2]
    ));

    // Unproject fragment position to get the ray in world space
    // For orthographic: ray direction = -fwd (constant for all fragments)
    // We need the ray origin for this fragment. Use the fragment's clip position
    // to find the world-space point on the near plane.

    // Simpler approach: find the closest point on the cylinder surface
    // along the light direction at the cylinder axis center line.

    let axis = input.end_world - input.start_world;
    let axis_len = length(axis);
    let axis_dir = axis / axis_len;

    // Project the midpoint along the cylinder
    let center = (input.start_world + input.end_world) * 0.5;

    // Fragment NDC position
    let frag_ndc = input.clip_position.xy / input.clip_position.w;

    // Find where on the axis this fragment corresponds to
    // Use the clip-space x,y to determine the closest axis point
    let center_clip = shadow.view_proj * vec4<f32>(center, 1.0);
    let start_clip = shadow.view_proj * vec4<f32>(input.start_world, 1.0);
    let end_clip = shadow.view_proj * vec4<f32>(input.end_world, 1.0);

    let center_ndc = center_clip.xy / center_clip.w;
    let start_ndc = start_clip.xy / start_clip.w;
    let end_ndc = end_clip.xy / end_clip.w;

    // Project fragment position onto cylinder axis in NDC
    let axis_ndc = end_ndc - start_ndc;
    let axis_ndc_len_sq = dot(axis_ndc, axis_ndc);
    var t = 0.5;
    if axis_ndc_len_sq > 1e-8 {
        t = clamp(dot(frag_ndc - start_ndc, axis_ndc) / axis_ndc_len_sq, 0.0, 1.0);
    }

    // Point on axis in world space
    let axis_point = input.start_world + axis_dir * (t * axis_len);

    // The closest point on cylinder surface facing the light
    let to_light = -fwd;
    let proj_on_axis = dot(to_light, axis_dir) * axis_dir;
    let perp = normalize(to_light - proj_on_axis);
    let surface_point = axis_point + perp * input.radius;

    // Check if fragment is within cylinder bounds in NDC
    let axis_point_clip = shadow.view_proj * vec4<f32>(axis_point, 1.0);
    let axis_ndc_point = axis_point_clip.xy / axis_point_clip.w;
    let dist_from_axis_ndc = length(frag_ndc - axis_ndc_point);

    let inv_sx = 1.0 / length(vec3<f32>(shadow.view_proj[0][0], shadow.view_proj[0][1], shadow.view_proj[0][2]));
    let inv_sy = 1.0 / length(vec3<f32>(shadow.view_proj[1][0], shadow.view_proj[1][1], shadow.view_proj[1][2]));

    // Approximate world-space distance from axis (using average scale)
    let avg_inv_s = (inv_sx + inv_sy) * 0.5;
    let world_dist_from_axis = dist_from_axis_ndc * avg_inv_s;

    if world_dist_from_axis > input.radius * 1.2 {
        discard;
    }

    // Write depth of the surface point
    let clip = shadow.view_proj * vec4<f32>(surface_point, 1.0);
    output.depth = clip.z / clip.w;

    return output;
}
