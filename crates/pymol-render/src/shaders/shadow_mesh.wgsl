// Depth-only mesh shader for shadow map rendering
// Transforms mesh vertices using the shadow view-projection matrix.
// No color output — writes only depth.

struct ShadowUniform {
    view_proj: mat4x4<f32>,
}

@group(0) @binding(0)
var<uniform> shadow: ShadowUniform;

struct VertexInput {
    @location(0) position: vec3<f32>,
    @location(1) normal: vec3<f32>,
    @location(2) color: vec4<f32>,
}

@vertex
fn vs_main(input: VertexInput) -> @builtin(position) vec4<f32> {
    return shadow.view_proj * vec4<f32>(input.position, 1.0);
}

@fragment
fn fs_main() {
    // Depth-only: no color output
}
