// Depth-only sphere impostor shader for shadow map rendering
// Ray-sphere intersection from the shadow light's viewpoint.
// Writes corrected depth — no color output.

struct ShadowUniform {
    view_proj: mat4x4<f32>,
}

@group(0) @binding(0)
var<uniform> shadow: ShadowUniform;

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
    @location(0) center_view: vec3<f32>,
    @location(1) ray_dir: vec3<f32>,
    @location(2) radius: f32,
}

// Extract the view matrix from view_proj is non-trivial for ortho,
// so we decompose: for orthographic shadow maps, the view part is the
// upper-left 3x3 (with translation in col 3). We compute the view-space
// positions using the full matrix's linear part.

@vertex
fn vs_main(billboard: BillboardVertex, instance: SphereInstance) -> VertexOutput {
    var output: VertexOutput;

    // For orthographic projection, we need to extract view transform.
    // shadow.view_proj = ortho * view. For ortho with symmetric bounds,
    // the inverse projection is simple scaling. We'll work in clip space
    // and use the ortho projection properties.

    // Transform sphere center to clip space
    let center_clip = shadow.view_proj * vec4<f32>(instance.center, 1.0);

    // For orthographic: clip.xy maps linearly to NDC.
    // Billboard offset in clip space needs to be scaled by radius in clip units.
    // We need to figure out the scale. The ortho matrix maps world-space radius
    // to clip-space offset. ortho[0][0] = 2/(right-left), which is the x scale.
    // So clip_radius = radius * ortho_scale. We use view_proj[0][0] and [1][1]
    // as the combined scale.
    let sx = length(vec3<f32>(shadow.view_proj[0][0], shadow.view_proj[0][1], shadow.view_proj[0][2]));
    let sy = length(vec3<f32>(shadow.view_proj[1][0], shadow.view_proj[1][1], shadow.view_proj[1][2]));
    let scale_x = instance.radius * sx * 1.5;
    let scale_y = instance.radius * sy * 1.5;

    output.clip_position = center_clip + vec4<f32>(
        billboard.offset.x * scale_x,
        billboard.offset.y * scale_y,
        0.0, 0.0
    );

    // For fragment shader, pass world-space center and radius.
    // The fragment shader will reconstruct the world position and do
    // ray-sphere intersection from the light direction.
    output.center_view = instance.center;
    output.ray_dir = vec3<f32>(billboard.offset, 0.0); // placeholder
    output.radius = instance.radius;

    return output;
}

struct FragmentOutput {
    @builtin(frag_depth) depth: f32,
}

@fragment
fn fs_main(input: VertexOutput) -> FragmentOutput {
    var output: FragmentOutput;

    // For orthographic shadow maps, the ray direction is constant:
    // it's the -Z axis of the view (the "forward" direction of the light).
    // We can extract it from the view_proj matrix.
    // The view part of view_proj: row 2 of the view matrix = forward direction.
    // In column-major: view_proj[0][2], view_proj[1][2], view_proj[2][2] is
    // the third row of the linear part (possibly scaled by ortho).
    // For orthographic, the ray is parallel — all fragments have the same direction.

    // Extract the light's forward direction (row 2 of the view matrix part)
    // For ortho * view, the z-row of the combined matrix gives us the depth mapping.
    // We need the view direction in world space. It's the 3rd row of the inverse view.
    // Simpler: the light direction is encoded in row 2 of the view_proj matrix.

    // Get forward direction from the matrix's 3rd row (columns 0-2)
    let fwd = normalize(vec3<f32>(
        shadow.view_proj[0][2],
        shadow.view_proj[1][2],
        shadow.view_proj[2][2]
    ));

    // Ray origin: for orthographic, each fragment has its own ray origin
    // We need to unproject the fragment position back to world space.
    // Instead of full unprojection, we use a simpler approach:
    // Project the sphere center and adjust depth based on the intersection.

    // The closest point on the sphere along the light direction is:
    // sphere_center + radius * (-fwd)  (the point facing the light)
    let front_point = input.center_view - fwd * input.radius;

    // Project to get the correct depth
    let clip = shadow.view_proj * vec4<f32>(front_point, 1.0);
    output.depth = clip.z / clip.w;

    // Check if this fragment is actually within the sphere's projected circle
    // by comparing the clip-space distance from the center
    let center_clip = shadow.view_proj * vec4<f32>(input.center_view, 1.0);
    let frag_ndc = input.clip_position.xy / input.clip_position.w;
    let center_ndc = center_clip.xy / center_clip.w;
    let dist_ndc = frag_ndc - center_ndc;

    // Scale factors to convert NDC distance to world distance
    let inv_sx = 1.0 / length(vec3<f32>(shadow.view_proj[0][0], shadow.view_proj[0][1], shadow.view_proj[0][2]));
    let inv_sy = 1.0 / length(vec3<f32>(shadow.view_proj[1][0], shadow.view_proj[1][1], shadow.view_proj[1][2]));
    let world_dist_sq = (dist_ndc.x * inv_sx) * (dist_ndc.x * inv_sx) +
                        (dist_ndc.y * inv_sy) * (dist_ndc.y * inv_sy);

    if world_dist_sq > input.radius * input.radius {
        discard;
    }

    // More accurate depth: use actual ray-sphere intersection depth
    let penetration = sqrt(max(input.radius * input.radius - world_dist_sq, 0.0));
    let hit_point = input.center_view - fwd * penetration;
    let hit_clip = shadow.view_proj * vec4<f32>(hit_point, 1.0);
    output.depth = hit_clip.z / hit_clip.w;

    return output;
}
