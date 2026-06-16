//! Embedded WGSL sources for the raytracer plugin.

pub(crate) const RAYTRACE: &str = include_str!("shaders/raytrace.wgsl");
pub(crate) const EDGE_DETECT: &str = include_str!("shaders/edge_detect.wgsl");
pub(crate) const COMPOSITE: &str = include_str!("shaders/composite.wgsl");

const ARTIFACT_COLOR: &str = include_str!("shaders/artifacts/color.wgsl");
const ARTIFACT_SPHERES: &str = include_str!("shaders/artifacts/build_spheres.wgsl");
const ARTIFACT_STICK_CYLINDERS: &str = include_str!("shaders/artifacts/build_stick_cylinders.wgsl");
const ARTIFACT_LINE_CYLINDERS: &str = include_str!("shaders/artifacts/build_line_cylinders.wgsl");
const ARTIFACT_TRIANGLES: &str = include_str!("shaders/artifacts/build_triangles.wgsl");

pub(crate) const ARTIFACT_BVH: &str = include_str!("shaders/artifacts/bvh.wgsl");
pub(crate) const ARTIFACT_RAYTRACE_OUTPUT_BUFFER: &str =
    include_str!("shaders/artifacts/raytrace_output_buffer.wgsl");
pub(crate) const ARTIFACT_DOWNSAMPLE: &str = include_str!("shaders/artifacts/downsample.wgsl");

const INCLUDE_ARTIFACT_COLOR: &str = "// {{INCLUDE_ARTIFACT_COLOR}}";

pub(crate) fn artifact_spheres() -> String {
    expand_artifact(ARTIFACT_SPHERES)
}

pub(crate) fn artifact_stick_cylinders() -> String {
    expand_artifact(ARTIFACT_STICK_CYLINDERS)
}

pub(crate) fn artifact_line_cylinders() -> String {
    expand_artifact(ARTIFACT_LINE_CYLINDERS)
}

pub(crate) fn artifact_triangles() -> String {
    expand_artifact(ARTIFACT_TRIANGLES)
}

fn expand_artifact(source: &str) -> String {
    source.replace(INCLUDE_ARTIFACT_COLOR, ARTIFACT_COLOR)
}
