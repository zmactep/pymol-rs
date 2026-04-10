//! Linear algebra utilities
//!
//! General-purpose matrix and vector operations used across the crate:
//!
//! - [`svd3`] — Analytical 3×3 SVD decomposition (Jacobi eigenvalue method)

pub mod svd3;

pub use svd3::{svd3, Svd3};
