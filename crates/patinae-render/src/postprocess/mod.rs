//! Full-screen postprocess passes.
//!
//! `silhouette` runs after the composite for edge detect. All passes use the
//! standard fullscreen-triangle trick (`@builtin(vertex_index)` with no vertex
//! buffer).

pub mod fxaa;
pub mod marking;
pub mod silhouette;
pub mod ssao_compose;
