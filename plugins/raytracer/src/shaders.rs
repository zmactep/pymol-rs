//! Embedded WGSL sources for the raytracer plugin.

pub(crate) const RAYTRACE: &str = include_str!("shaders/raytrace.wgsl");
pub(crate) const EDGE_DETECT: &str = include_str!("shaders/edge_detect.wgsl");
pub(crate) const COMPOSITE: &str = include_str!("shaders/composite.wgsl");

const ARTIFACT_COLOR: &str = include_str!("shaders/artifacts/color.wgsl");
const ARTIFACT_TRIANGLE_STORAGE: &str = include_str!("shaders/artifacts/triangle_storage.wgsl");
const ARTIFACT_STD_VERTEX: &str = include_str!("shaders/artifacts/std_vertex.wgsl");
const ARTIFACT_ATOM_ALPHA: &str = include_str!("shaders/artifacts/atom_alpha.wgsl");
const ARTIFACT_SPHERES: &str = include_str!("shaders/artifacts/build_spheres.wgsl");
const ARTIFACT_STICK_CAPSULES: &str = include_str!("shaders/artifacts/build_stick_capsules.wgsl");
const ARTIFACT_LINE_CYLINDERS: &str = include_str!("shaders/artifacts/build_line_cylinders.wgsl");
const ARTIFACT_TRIANGLES: &str = include_str!("shaders/artifacts/build_triangles.wgsl");
const ARTIFACT_COUNT_VISIBLE_TRIANGLES: &str =
    include_str!("shaders/artifacts/count_visible_triangles.wgsl");
const ARTIFACT_VISIBLE_TRIANGLES: &str =
    include_str!("shaders/artifacts/build_visible_triangles.wgsl");

const ARTIFACT_BVH: &str = include_str!("shaders/artifacts/bvh.wgsl");
const ARTIFACT_RAYTRACE_OUTPUT_BUFFER: &str =
    include_str!("shaders/artifacts/raytrace_output_buffer.wgsl");
pub(crate) const ARTIFACT_DOWNSAMPLE: &str = include_str!("shaders/artifacts/downsample.wgsl");

const INCLUDE_ARTIFACT_COLOR: &str = "// {{INCLUDE_ARTIFACT_COLOR}}";
const INCLUDE_ARTIFACT_TRIANGLE: &str = "// {{INCLUDE_ARTIFACT_TRIANGLE}}";
const INCLUDE_ARTIFACT_STD_VERTEX: &str = "// {{INCLUDE_ARTIFACT_STD_VERTEX}}";
const INCLUDE_ARTIFACT_ATOM_ALPHA: &str = "// {{INCLUDE_ARTIFACT_ATOM_ALPHA}}";

pub(crate) fn artifact_spheres() -> String {
    expand_artifact(ARTIFACT_SPHERES)
}

pub(crate) fn artifact_stick_capsules() -> String {
    expand_artifact(ARTIFACT_STICK_CAPSULES)
}

pub(crate) fn artifact_line_cylinders() -> String {
    expand_artifact(ARTIFACT_LINE_CYLINDERS)
}

pub(crate) fn artifact_triangles() -> String {
    expand_artifact(ARTIFACT_TRIANGLES)
}

pub(crate) fn artifact_count_visible_triangles() -> String {
    ARTIFACT_COUNT_VISIBLE_TRIANGLES.to_string()
}

pub(crate) fn artifact_visible_triangles() -> String {
    expand_artifact(ARTIFACT_VISIBLE_TRIANGLES)
}

pub(crate) fn artifact_bvh() -> String {
    expand_artifact(ARTIFACT_BVH)
}

pub(crate) fn artifact_raytrace_output_buffer() -> String {
    expand_artifact(ARTIFACT_RAYTRACE_OUTPUT_BUFFER)
}

fn expand_artifact(source: &str) -> String {
    source
        .replace(INCLUDE_ARTIFACT_COLOR, ARTIFACT_COLOR)
        .replace(INCLUDE_ARTIFACT_TRIANGLE, ARTIFACT_TRIANGLE_STORAGE)
        .replace(INCLUDE_ARTIFACT_STD_VERTEX, ARTIFACT_STD_VERTEX)
        .replace(INCLUDE_ARTIFACT_ATOM_ALPHA, ARTIFACT_ATOM_ALPHA)
}
