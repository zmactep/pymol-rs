// Cylinder impostor shader — ray-cylinder intersection for stick/bond rendering.
// Requires: common.wgsl + lighting.wgsl prepended via concat!(include_str!(...))

struct BillboardVertex {
    @location(10) offset: vec2<f32>,
}

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
    @location(6) @interpolate(flat) flags: u32,
}

@vertex
fn vs_main(billboard: BillboardVertex, instance: CylinderInstance) -> VertexOutput {
    var output: VertexOutput;

    let start_view = (uniforms.view * vec4<f32>(instance.start, 1.0)).xyz;
    let end_view = (uniforms.view * vec4<f32>(instance.end, 1.0)).xyz;
    let center_view = (start_view + end_view) * 0.5;

    let axis = end_view - start_view;
    let axis_len = length(axis);
    let half_len = axis_len * 0.5;
    let bound_radius = sqrt(half_len * half_len + instance.radius * instance.radius) + instance.radius;

    let scale = bound_radius * 1.5;
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

fn ray_cylinder_intersect(
    ray_origin: vec3<f32>,
    ray_dir: vec3<f32>,
    cyl_start: vec3<f32>,
    cyl_end: vec3<f32>,
    radius: f32
) -> vec4<f32> {
    let axis = cyl_end - cyl_start;
    let axis_len = length(axis);
    let axis_dir = axis / axis_len;

    let oc = ray_origin - cyl_start;

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

    let ray_origin = vec3<f32>(0.0, 0.0, 0.0);
    let ray_dir = normalize(input.ray_origin);

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
    let u = hit.y;

    let hit_point = ray_origin + t * ray_dir;

    let axis = input.end_view - input.start_view;
    let axis_dir = normalize(axis);
    let center_on_axis = input.start_view + u * axis;
    var normal = normalize(hit_point - center_on_axis);

    let view_dir = normalize(-hit_point);

    if dot(normal, view_dir) < 0.0 {
        normal = -normal;
    }

    let base_color = mix(input.color1.rgb, input.color2.rgb, u);
    let alpha = mix(input.color1.a, input.color2.a, u);

    // Calculate lighting — branch on shadow mode
    let world_hit = (uniforms.view_inv * vec4<f32>(hit_point, 1.0)).xyz;
    var color: vec3<f32>;
    if shadow_params.mode == 2u {
        color = pymol_lighting_shadowed(normal, view_dir, base_color, world_hit);
    } else {
        color = pymol_lighting(normal, view_dir, base_color);
        let ao = compute_shadow_ao(world_hit);
        color *= ao;
    }

    let depth = -hit_point.z;
    color = apply_depth_cue(color, depth);
    color = apply_fog(color, depth);

    let clip_pos = uniforms.proj * vec4<f32>(hit_point, 1.0);
    output.depth = clip_pos.z / clip_pos.w;
    output.color = vec4<f32>(color, alpha);

    return output;
}
