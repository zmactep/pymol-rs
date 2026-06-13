// Sphere picking — same impostor + ray-sphere as the colour shader, writes
// `PackedId(rep_kind=Sphere, object_id, atom_id=group_id)` into Rg32Uint.

// {{INCLUDE_FRAME}}
// {{INCLUDE_PICKING}}

@group(2) @binding(0) var<uniform> picking: PickingParams;

struct ObjectEntry {
    atom_offset: u32,
    atom_count:  u32,
    bond_offset: u32,
    bond_count:  u32,
    object_id:   u32,
    flags:       u32,
    _pad0_a:     u32,
    _pad0_b:     u32,
    model_matrix: mat4x4<f32>,
};

@group(3) @binding(0) var<uniform> obj: ObjectEntry;
@group(3) @binding(5) var<storage, read> mask_lut: array<u32>;

struct SphereInstance {
    @location(0) center:   vec3<f32>,
    @location(1) radius:   f32,
    @location(2) group_id: u32,
};

struct VsOut {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) center_view: vec3<f32>,
    @location(1) ray_dir:     vec3<f32>,
    @location(2) radius:      f32,
    @location(3) @interpolate(flat) group_id: u32,
};

fn billboard_offset(vid: u32) -> vec2<f32> {
    let id = vid % 6u;
    var x: f32 = -1.0;
    var y: f32 = -1.0;
    if id == 1u || id == 4u || id == 5u { x = 1.0; }
    if id == 2u || id == 3u || id == 5u { y = 1.0; }
    return vec2<f32>(x, y);
}

fn scene_visible(gid: u32) -> bool {
    let word = gid >> 5u;
    let bit  = gid & 31u;
    return (mask_lut[word] & (1u << bit)) != 0u;
}

@vertex
fn vs_main(
    @builtin(vertex_index) vid: u32,
    instance: SphereInstance,
) -> VsOut {
    var out: VsOut;
    let off = billboard_offset(vid);
    let center_view = (frame.view * vec4<f32>(instance.center, 1.0)).xyz;
    let scale = instance.radius * 1.5;
    let billboard_pos = center_view + vec3<f32>(off * scale, 0.0);
    out.clip_position = frame.proj * vec4<f32>(billboard_pos, 1.0);
    out.center_view = center_view;
    out.radius      = instance.radius;
    out.ray_dir     = billboard_pos;
    out.group_id    = instance.group_id;
    return out;
}

struct PickOut {
    @location(0) id:    vec2<u32>,
    @builtin(frag_depth) depth: f32,
};

@fragment
fn fs_main(input: VsOut) -> PickOut {
    let global_id = obj.atom_offset + input.group_id;
    if !scene_visible(global_id) {
        discard;
    }

    let ray_dir = normalize(input.ray_dir);
    let oc = -input.center_view;
    let a = dot(ray_dir, ray_dir);
    let b = 2.0 * dot(oc, ray_dir);
    let c = dot(oc, oc) - input.radius * input.radius;
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        discard;
    }
    let t = (-b - sqrt(disc)) / (2.0 * a);
    if t < 0.0 {
        discard;
    }
    // Override the rasterized depth (which came from the bounding billboard
    // at z = sphere centre) with the actual ray-sphere hit. This lets the
    // picking depth-test resolve correctly between impostors and meshes.
    let hit_view = ray_dir * t;
    let clip = frame.proj * vec4<f32>(hit_view, 1.0);
    var out: PickOut;
    out.id    = pack_id(picking.rep_object, input.group_id);
    out.depth = clip.z / clip.w;
    return out;
}
