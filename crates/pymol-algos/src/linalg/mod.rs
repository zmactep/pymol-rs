//! Linear algebra utilities
//!
//! General-purpose matrix and vector operations used across the crate:
//!
//! - [`svd3`] — Analytical 3×3 SVD decomposition (Jacobi eigenvalue method)
//! - [`mat4`] — 4×4 / 3×3 row-major matrix operations (multiply, transform, invert)

pub mod svd3;
pub mod mat4;

pub use svd3::{svd3, Svd3};
pub use mat4::{
    mat3x3_to_mat4, left_multiply_mat4, transform_mat4, is_identity_mat4,
    transform_3x3, invert_3x3,
};
