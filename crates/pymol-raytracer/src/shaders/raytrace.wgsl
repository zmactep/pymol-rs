// GPU compute shader for molecular raytracing
//
// This shader performs raytracing of molecular structures including:
// - Spheres (atoms)
// - Cylinders (bonds)
// - Triangles (surfaces, cartoons)
//
// Features:
// - BVH acceleration for efficient ray traversal
// - Shadows
// - Transparency with depth sorting
// - Fog/depth cue effects

// Uniform buffer
struct Uniforms {
    view_matrix: mat4x4<f32>,
    proj_matrix: mat4x4<f32>,
    view_inv_matrix: mat4x4<f32>,
    proj_inv_matrix: mat4x4<f32>,
    camera_pos: vec4<f32>,
    viewport: vec4<f32>,  // width, height, 1/width, 1/height
    light_dir: vec4<f32>,
    ambient: f32,
    direct: f32,
    reflect: f32,
    specular: f32,
    shininess: f32,
    _pad_light1: f32,
    _pad_light2: f32,
    _pad_light3: f32,
    bg_color: vec4<f32>,
    fog_start: f32,
    fog_end: f32,
    fog_density: f32,
    _pad0: f32,
    fog_color: vec4<f32>,
    sphere_count: u32,
    cylinder_count: u32,
    triangle_count: u32,
    bvh_node_count: u32,
    ray_shadow: u32,
    ray_max_passes: u32,
    ray_trace_fog: u32,
    ray_transparency_shadows: u32,
}

// Sphere primitive
// Note: Using individual f32 for padding to match Rust struct layout
// (vec3<f32> has 16-byte alignment in WGSL which would cause layout mismatch)
struct Sphere {
    center: vec3<f32>,
    radius: f32,
    color: vec4<f32>,
    transparency: f32,
    _pad_a: f32, _pad_b: f32, _pad_c: f32,
}

// Cylinder primitive
struct Cylinder {
    start: vec3<f32>,
    radius: f32,
    end: vec3<f32>,
    _pad0: f32,
    color1: vec4<f32>,
    color2: vec4<f32>,
    transparency: f32,
    _pad_a: f32, _pad_b: f32, _pad_c: f32,
}

// Triangle primitive
struct Triangle {
    v0: vec3<f32>, _pad0: f32,
    v1: vec3<f32>, _pad1: f32,
    v2: vec3<f32>, _pad2: f32,
    n0: vec3<f32>, _pad3: f32,
    n1: vec3<f32>, _pad4: f32,
    n2: vec3<f32>, _pad5: f32,
    color: vec4<f32>,
    transparency: f32,
    _pad_a: f32, _pad_b: f32, _pad_c: f32,
}

// BVH node
struct BvhNode {
    min: vec3<f32>,
    left_or_first: u32,
    max: vec3<f32>,
    count: u32,
}

// Ray structure
struct Ray {
    origin: vec3<f32>,
    direction: vec3<f32>,
}

// Hit information
struct HitInfo {
    t: f32,
    normal: vec3<f32>,
    color: vec4<f32>,
    transparency: f32,
    hit: bool,
}

// Bindings
@group(0) @binding(0) var<uniform> uniforms: Uniforms;
@group(0) @binding(1) var<storage, read> spheres: array<Sphere>;
@group(0) @binding(2) var<storage, read> cylinders: array<Cylinder>;
@group(0) @binding(3) var<storage, read> triangles: array<Triangle>;
@group(0) @binding(4) var<storage, read> bvh_nodes: array<BvhNode>;
@group(0) @binding(5) var<storage, read> bvh_indices: array<u32>;
@group(0) @binding(6) var output: texture_storage_2d<rgba8unorm, write>;

// Constants
const EPSILON: f32 = 0.0001;
const MAX_T: f32 = 1e30;
const MAX_STACK_SIZE: u32 = 64u;

// Ray-sphere intersection
fn intersect_sphere(ray: Ray, sphere: Sphere) -> HitInfo {
    var hit: HitInfo;
    hit.hit = false;
    hit.t = MAX_T;
    
    let oc = ray.origin - sphere.center;
    let a = dot(ray.direction, ray.direction);
    let b = 2.0 * dot(oc, ray.direction);
    let c = dot(oc, oc) - sphere.radius * sphere.radius;
    let discriminant = b * b - 4.0 * a * c;
    
    if discriminant < 0.0 {
        return hit;
    }
    
    let sqrt_d = sqrt(discriminant);
    var t = (-b - sqrt_d) / (2.0 * a);
    
    if t < EPSILON {
        t = (-b + sqrt_d) / (2.0 * a);
    }
    
    if t < EPSILON {
        return hit;
    }
    
    hit.hit = true;
    hit.t = t;
    let hit_point = ray.origin + t * ray.direction;
    hit.normal = normalize(hit_point - sphere.center);
    hit.color = sphere.color;
    hit.transparency = sphere.transparency;
    
    return hit;
}

// Ray-cylinder intersection
fn intersect_cylinder(ray: Ray, cyl: Cylinder) -> HitInfo {
    var hit: HitInfo;
    hit.hit = false;
    hit.t = MAX_T;
    
    let axis = cyl.end - cyl.start;
    let axis_len = length(axis);
    if axis_len < EPSILON {
        return hit;
    }
    let axis_norm = axis / axis_len;
    
    // Transform ray to cylinder-local space
    let oc = ray.origin - cyl.start;
    let d_perp = ray.direction - dot(ray.direction, axis_norm) * axis_norm;
    let o_perp = oc - dot(oc, axis_norm) * axis_norm;
    
    let a = dot(d_perp, d_perp);
    let b = 2.0 * dot(o_perp, d_perp);
    let c = dot(o_perp, o_perp) - cyl.radius * cyl.radius;
    
    let discriminant = b * b - 4.0 * a * c;
    
    if discriminant < 0.0 {
        return hit;
    }
    
    let sqrt_d = sqrt(discriminant);
    var t = (-b - sqrt_d) / (2.0 * a);
    
    if t < EPSILON {
        t = (-b + sqrt_d) / (2.0 * a);
    }
    
    if t < EPSILON {
        return hit;
    }
    
    // Check if hit is within cylinder bounds
    let hit_point = ray.origin + t * ray.direction;
    let height = dot(hit_point - cyl.start, axis_norm);
    
    if height < 0.0 || height > axis_len {
        return hit;
    }
    
    hit.hit = true;
    hit.t = t;
    
    // Calculate normal (perpendicular to axis)
    let closest_on_axis = cyl.start + height * axis_norm;
    hit.normal = normalize(hit_point - closest_on_axis);
    
    // Interpolate color based on position along cylinder
    let t_along = height / axis_len;
    hit.color = mix(cyl.color1, cyl.color2, t_along);
    hit.transparency = cyl.transparency;
    
    return hit;
}

// Ray-triangle intersection (Möller-Trumbore)
fn intersect_triangle(ray: Ray, tri: Triangle) -> HitInfo {
    var hit: HitInfo;
    hit.hit = false;
    hit.t = MAX_T;
    
    let e1 = tri.v1 - tri.v0;
    let e2 = tri.v2 - tri.v0;
    let h = cross(ray.direction, e2);
    let a = dot(e1, h);
    
    if abs(a) < EPSILON {
        return hit;
    }
    
    let f = 1.0 / a;
    let s = ray.origin - tri.v0;
    let u = f * dot(s, h);
    
    if u < 0.0 || u > 1.0 {
        return hit;
    }
    
    let q = cross(s, e1);
    let v = f * dot(ray.direction, q);
    
    if v < 0.0 || u + v > 1.0 {
        return hit;
    }
    
    let t = f * dot(e2, q);
    
    if t < EPSILON {
        return hit;
    }
    
    hit.hit = true;
    hit.t = t;
    
    // Interpolate normal
    let w = 1.0 - u - v;
    hit.normal = normalize(w * tri.n0 + u * tri.n1 + v * tri.n2);
    hit.color = tri.color;
    hit.transparency = tri.transparency;
    
    return hit;
}

// Ray-AABB intersection
fn intersect_aabb(ray: Ray, box_min: vec3<f32>, box_max: vec3<f32>) -> bool {
    let inv_dir = 1.0 / ray.direction;
    
    let t1 = (box_min - ray.origin) * inv_dir;
    let t2 = (box_max - ray.origin) * inv_dir;
    
    let tmin_v = min(t1, t2);
    let tmax_v = max(t1, t2);
    
    let tmin = max(max(tmin_v.x, tmin_v.y), tmin_v.z);
    let tmax = min(min(tmax_v.x, tmax_v.y), tmax_v.z);
    
    return tmax >= max(tmin, 0.0);
}

// Decode primitive index
fn decode_index(encoded: u32) -> vec2<u32> {
    let prim_type = encoded >> 30u;
    let index = encoded & 0x3FFFFFFFu;
    return vec2<u32>(prim_type, index);
}

// Trace ray through BVH
fn trace_ray(ray: Ray) -> HitInfo {
    var closest: HitInfo;
    closest.hit = false;
    closest.t = MAX_T;
    
    // If no BVH nodes, do brute force
    if arrayLength(&bvh_nodes) == 0u {
        // Brute force spheres
        for (var i = 0u; i < uniforms.sphere_count; i++) {
            let hit = intersect_sphere(ray, spheres[i]);
            if hit.hit && hit.t < closest.t {
                closest = hit;
            }
        }
        // Brute force cylinders
        for (var i = 0u; i < uniforms.cylinder_count; i++) {
            let hit = intersect_cylinder(ray, cylinders[i]);
            if hit.hit && hit.t < closest.t {
                closest = hit;
            }
        }
        // Brute force triangles
        for (var i = 0u; i < uniforms.triangle_count; i++) {
            let hit = intersect_triangle(ray, triangles[i]);
            if hit.hit && hit.t < closest.t {
                closest = hit;
            }
        }
        return closest;
    }
    
    // BVH traversal with stack
    var stack: array<u32, 64>;
    var stack_ptr = 0u;
    stack[stack_ptr] = 0u;
    stack_ptr++;
    
    while stack_ptr > 0u {
        stack_ptr--;
        let node_idx = stack[stack_ptr];
        let node = bvh_nodes[node_idx];
        
        // Test AABB
        if !intersect_aabb(ray, node.min, node.max) {
            continue;
        }
        
        if node.count > 0u {
            // Leaf node - test primitives
            let first = node.left_or_first;
            let count = node.count;
            
            for (var i = 0u; i < count; i++) {
                let encoded = bvh_indices[first + i];
                let decoded = decode_index(encoded);
                let prim_type = decoded.x;
                let prim_idx = decoded.y;
                
                var hit: HitInfo;
                hit.hit = false;
                
                if prim_type == 0u && prim_idx < uniforms.sphere_count {
                    hit = intersect_sphere(ray, spheres[prim_idx]);
                } else if prim_type == 1u && prim_idx < uniforms.cylinder_count {
                    hit = intersect_cylinder(ray, cylinders[prim_idx]);
                } else if prim_type == 2u && prim_idx < uniforms.triangle_count {
                    hit = intersect_triangle(ray, triangles[prim_idx]);
                }
                
                if hit.hit && hit.t < closest.t {
                    closest = hit;
                }
            }
        } else {
            // Internal node - push children
            let left = node.left_or_first;
            let right = left + 1u;
            
            if stack_ptr < MAX_STACK_SIZE - 1u {
                stack[stack_ptr] = left;
                stack_ptr++;
                stack[stack_ptr] = right;
                stack_ptr++;
            }
        }
    }
    
    return closest;
}

// Check shadow ray using BVH acceleration
fn trace_shadow(ray: Ray, max_t: f32) -> bool {
    // If no BVH nodes, do brute force (fallback for small scenes)
    if arrayLength(&bvh_nodes) == 0u {
        for (var i = 0u; i < uniforms.sphere_count; i++) {
            let hit = intersect_sphere(ray, spheres[i]);
            if hit.hit && hit.t > EPSILON && hit.t < max_t && hit.transparency < 0.5 {
                return true;
            }
        }
        for (var i = 0u; i < uniforms.cylinder_count; i++) {
            let hit = intersect_cylinder(ray, cylinders[i]);
            if hit.hit && hit.t > EPSILON && hit.t < max_t && hit.transparency < 0.5 {
                return true;
            }
        }
        for (var i = 0u; i < uniforms.triangle_count; i++) {
            let hit = intersect_triangle(ray, triangles[i]);
            if hit.hit && hit.t > EPSILON && hit.t < max_t && hit.transparency < 0.5 {
                return true;
            }
        }
        return false;
    }
    
    // BVH traversal - early exit on any hit
    var stack: array<u32, 64>;
    var stack_ptr = 0u;
    stack[stack_ptr] = 0u;
    stack_ptr++;
    
    while stack_ptr > 0u {
        stack_ptr--;
        let node_idx = stack[stack_ptr];
        let node = bvh_nodes[node_idx];
        
        // Test AABB
        if !intersect_aabb(ray, node.min, node.max) {
            continue;
        }
        
        if node.count > 0u {
            // Leaf node - test primitives
            let first = node.left_or_first;
            let count = node.count;
            
            for (var i = 0u; i < count; i++) {
                let encoded = bvh_indices[first + i];
                let decoded = decode_index(encoded);
                let prim_type = decoded.x;
                let prim_idx = decoded.y;
                
                var hit: HitInfo;
                hit.hit = false;
                
                if prim_type == 0u && prim_idx < uniforms.sphere_count {
                    hit = intersect_sphere(ray, spheres[prim_idx]);
                } else if prim_type == 1u && prim_idx < uniforms.cylinder_count {
                    hit = intersect_cylinder(ray, cylinders[prim_idx]);
                } else if prim_type == 2u && prim_idx < uniforms.triangle_count {
                    hit = intersect_triangle(ray, triangles[prim_idx]);
                }
                
                // Early exit on any valid shadow hit
                if hit.hit && hit.t > EPSILON && hit.t < max_t && hit.transparency < 0.5 {
                    return true;
                }
            }
        } else {
            // Internal node - push children
            let left = node.left_or_first;
            let right = left + 1u;
            
            if stack_ptr < MAX_STACK_SIZE - 1u {
                stack[stack_ptr] = left;
                stack_ptr++;
                stack[stack_ptr] = right;
                stack_ptr++;
            }
        }
    }
    
    return false;
}

// Shade a hit point - PyMOL-compatible lighting
fn shade(ray: Ray, hit: HitInfo) -> vec3<f32> {
    let hit_point = ray.origin + hit.t * ray.direction;
    let normal = hit.normal;
    let view_dir = normalize(-ray.direction);
    
    // Ensure normal faces the viewer
    var n = normal;
    if dot(n, view_dir) < 0.0 {
        n = -n;
    }
    
    // PyMOL light directions are in VIEW SPACE (relative to camera)
    // Transform from view space to world space so light rotates with camera
    // light setting is direction FROM light TO scene, negate to get FROM surface TO light
    let light_view = -uniforms.light_dir.xyz;
    let light_dir = normalize((uniforms.view_inv_matrix * vec4<f32>(light_view, 0.0)).xyz);
    
    // Secondary light (light2) for "reflect" component
    // light2 default in view space: (-0.55, -0.7, 0.15), negate for surface-to-light
    let light2_view = vec3<f32>(0.55, 0.7, -0.15);
    let light2_dir = normalize((uniforms.view_inv_matrix * vec4<f32>(light2_view, 0.0)).xyz);
    
    // === AMBIENT ===
    // PyMOL: ambient = 0.14
    var color = hit.color.rgb * uniforms.ambient;
    
    // === DIRECT (primary light) ===
    let ndotl = max(dot(n, light_dir), 0.0);
    
    // Shadow check for primary light (reduced bias for tighter shadows)
    var shadow_factor = 1.0;
    if uniforms.ray_shadow != 0u && ndotl > 0.001 {
        let shadow_bias = max(hit.t * 0.0005, 0.05);
        let shadow_origin = hit_point + n * shadow_bias;
        let shadow_ray = Ray(shadow_origin, light_dir);
        if trace_shadow(shadow_ray, MAX_T) {
            shadow_factor = 0.0;  // PyMOL shadows are hard
        }
    }
    
    // Direct diffuse: PyMOL direct = 0.45
    color += hit.color.rgb * ndotl * uniforms.direct * shadow_factor;
    
    // === REFLECT (secondary fill light, no shadows) ===
    // PyMOL reflect default = 0.45
    let ndotl2 = max(dot(n, light2_dir), 0.0);
    color += hit.color.rgb * ndotl2 * uniforms.reflect;
    
    // === SPECULAR ===
    // PyMOL: specular = 1.0, specular_intensity = 0.5, shininess = 55
    if uniforms.specular > 0.0 && shadow_factor > 0.5 {
        // Primary light specular
        let h1 = normalize(light_dir + view_dir);
        let spec1 = pow(max(dot(n, h1), 0.0), uniforms.shininess);
        
        // Secondary light specular (weaker)
        let h2 = normalize(light2_dir + view_dir);
        let spec2 = pow(max(dot(n, h2), 0.0), uniforms.shininess);
        
        let spec_intensity = 0.5;  // specular_intensity default
        color += vec3<f32>(1.0) * (spec1 + spec2 * 0.3) * uniforms.specular * spec_intensity;
    }
    
    return color;
}

// Apply fog
fn apply_fog(color: vec3<f32>, depth: f32) -> vec3<f32> {
    if uniforms.fog_density <= 0.0 || uniforms.ray_trace_fog == 0u {
        return color;
    }
    let fog_factor = clamp((uniforms.fog_end - depth) / (uniforms.fog_end - uniforms.fog_start), 0.0, 1.0);
    return mix(uniforms.fog_color.rgb, color, fog_factor);
}

// Generate primary ray from pixel coordinates
fn generate_ray(pixel: vec2<f32>) -> Ray {
    // Convert pixel to NDC (x,y: -1 to 1)
    let ndc = vec2<f32>(
        (pixel.x + 0.5) / uniforms.viewport.x * 2.0 - 1.0,
        1.0 - (pixel.y + 0.5) / uniforms.viewport.y * 2.0
    );
    
    // Extract FOV parameters from projection matrix
    // For perspective: proj[0][0] = 1/(aspect*tan(fov/2)), proj[1][1] = 1/tan(fov/2)
    // So tan(fov/2) = 1/proj[1][1], and aspect-scaled x = ndc.x * aspect * tan(fov/2)
    let tan_half_fov = 1.0 / uniforms.proj_matrix[1][1];
    let aspect = uniforms.proj_matrix[1][1] / uniforms.proj_matrix[0][0];
    
    // Ray direction in view space: looking down -Z, with X right and Y up
    let ray_dir_view = normalize(vec3<f32>(
        ndc.x * aspect * tan_half_fov,
        ndc.y * tan_half_fov,
        -1.0  // Looking down negative Z in view space
    ));
    
    // Camera is at origin in view space, transform to world space
    let origin = (uniforms.view_inv_matrix * vec4<f32>(0.0, 0.0, 0.0, 1.0)).xyz;
    
    // Transform direction to world space (w=0 for direction)
    let direction = normalize((uniforms.view_inv_matrix * vec4<f32>(ray_dir_view, 0.0)).xyz);
    
    return Ray(origin, direction);
}

@compute @workgroup_size(8, 8, 1)
fn main(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let pixel = vec2<f32>(f32(global_id.x), f32(global_id.y));
    
    // Bounds check
    if pixel.x >= uniforms.viewport.x || pixel.y >= uniforms.viewport.y {
        return;
    }
    
    // Generate ray
    let ray = generate_ray(pixel);
    
    // Trace primary ray
    let hit = trace_ray(ray);
    
    var color: vec3<f32>;
    var alpha: f32;
    
    if hit.hit {
        // Shade hit point
        color = shade(ray, hit);
        alpha = 1.0;
        
        // Apply fog based on depth
        let depth = hit.t;
        color = apply_fog(color, depth);
        
        // Handle transparency (simplified)
        if hit.transparency > 0.001 {
            alpha = 1.0 - hit.transparency;
            color = mix(uniforms.bg_color.rgb, color, alpha);
        }
    } else {
        // Transparent background (alpha=0) like PyMOL
        color = vec3<f32>(0.0, 0.0, 0.0);
        alpha = 0.0;
    }
    
    // Clamp and output with alpha
    color = clamp(color, vec3<f32>(0.0), vec3<f32>(1.0));
    textureStore(output, vec2<i32>(global_id.xy), vec4<f32>(color, alpha));
}
