struct ConvertParams {
    instance_capacity: u32,
    sphere_offset: u32,
    atom_offset: u32,
    rep_slot: u32,
    transparency: f32,
    dispatch_width: u32,
    _pad0: u32,
    _pad1: u32,
}

struct SphereInstance {
    center: vec3<f32>,
    radius: f32,
    group_id: u32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
}

struct Sphere {
    center: vec3<f32>,
    radius: f32,
    color: vec4<f32>,
    transparency: f32,
    _pad_a: f32, _pad_b: f32, _pad_c: f32,
}

// {{INCLUDE_ARTIFACT_COLOR}}

@group(0) @binding(0) var<storage, read> instances: array<SphereInstance>;
@group(0) @binding(1) var<storage, read> color_lut: array<ColorLutEntry>;
@group(0) @binding(2) var<storage, read> draw_args: array<u32>;
@group(0) @binding(3) var<storage, read_write> spheres: array<Sphere>;
@group(0) @binding(4) var<uniform> params: ConvertParams;

fn material_for_group(group_id: u32) -> vec4<f32> {
    let entry = color_lut[params.atom_offset + group_id];
    let packed = rep_packed_color(entry.reps, params.rep_slot);
    var rgba = select(unpack_rep_rgb8(packed), entry.base, packed == REP_COLOR_INHERIT);
    rgba.a = clamp(rgba.a * (1.0 - params.transparency), 0.0, 1.0);
    return rgba;
}

fn invalid_sphere() -> Sphere {
    return Sphere(vec3<f32>(0.0), 0.0, vec4<f32>(0.0), 1.0, 0.0, 0.0, 0.0);
}

@compute @workgroup_size(128)
fn build_spheres(@builtin(global_invocation_id) gid: vec3<u32>) {
    let instance_index = gid.y * params.dispatch_width + gid.x;
    if instance_index >= params.instance_capacity {
        return;
    }
    let out_index = params.sphere_offset + instance_index;
    if instance_index >= draw_args[1] {
        spheres[out_index] = invalid_sphere();
        return;
    }

    let inst = instances[instance_index];
    let material = material_for_group(inst.group_id);
    spheres[out_index] = Sphere(
        inst.center,
        inst.radius,
        material,
        1.0 - material.a,
        0.0, 0.0, 0.0,
    );
}
