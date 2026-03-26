// Sphere impostor shader — ray-sphere intersection in fragment shader.
// Requires: common.wgsl + lighting.wgsl prepended via concat!(include_str!(...))

struct BillboardVertex {
    @location(10) offset: vec2<f32>,
}

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

    let center_world = vec4<f32>(instance.center, 1.0);
    let center_view = (uniforms.view * center_world).xyz;

    let scale = instance.radius * 1.5;
    let billboard_pos = center_view + vec3<f32>(billboard.offset * scale, 0.0);

    output.clip_position = uniforms.proj * vec4<f32>(billboard_pos, 1.0);
    output.color = instance.color;
    output.center_view = center_view;
    output.radius = instance.radius;
    output.ray_dir = normalize(billboard_pos);

    return output;
}

struct FragmentOutput {
    @location(0) color: vec4<f32>,
    @builtin(frag_depth) depth: f32,
}

@fragment
fn fs_main(input: VertexOutput) -> FragmentOutput {
    var output: FragmentOutput;

    let ray_origin = vec3<f32>(0.0, 0.0, 0.0);
    let ray_dir = normalize(input.ray_dir);
    let sphere_center = input.center_view;
    let radius = input.radius;

    let oc = ray_origin - sphere_center;
    let a = dot(ray_dir, ray_dir);
    let b = 2.0 * dot(oc, ray_dir);
    let c = dot(oc, oc) - radius * radius;
    let discriminant = b * b - 4.0 * a * c;

    if discriminant < 0.0 {
        discard;
    }

    let t = (-b - sqrt(discriminant)) / (2.0 * a);

    if t < 0.0 {
        discard;
    }

    let hit_point = ray_origin + t * ray_dir;
    var normal = normalize(hit_point - sphere_center);
    let view_dir = normalize(-hit_point);

    if dot(normal, view_dir) < 0.0 {
        normal = -normal;
    }

    // Calculate lighting — branch on shadow mode
    let world_hit = (uniforms.view_inv * vec4<f32>(hit_point, 1.0)).xyz;
    var color: vec3<f32>;
    if shadow_params.mode == 2u {
        color = pymol_lighting_shadowed(normal, view_dir, input.color.rgb, world_hit);
    } else {
        color = pymol_lighting(normal, view_dir, input.color.rgb);
        let ao = compute_shadow_ao(world_hit);
        color *= ao;
    }

    let depth = -hit_point.z;
    color = apply_depth_cue(color, depth);
    color = apply_fog(color, depth);

    let clip_pos = uniforms.proj * vec4<f32>(hit_point, 1.0);
    output.depth = clip_pos.z / clip_pos.w;
    output.color = vec4<f32>(color, input.color.a);

    return output;
}
