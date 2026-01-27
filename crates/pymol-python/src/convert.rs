//! Type conversions between Rust and Python types
//!
//! This module provides utilities for converting between:
//! - `lin_alg::f32::Vec3` <-> Python tuple `(f32, f32, f32)`
//! - `pymol_color::Color` <-> Python tuple `(f32, f32, f32)`
//! - Various array/list conversions

// Allow unused functions - these are utility functions for the complete API
#![allow(dead_code)]

use lin_alg::f32::Vec3;
use pymol_color::Color;

/// Convert a Vec3 to a Python tuple
pub fn vec3_to_tuple(v: Vec3) -> (f32, f32, f32) {
    (v.x, v.y, v.z)
}

/// Convert a Python tuple to Vec3
pub fn tuple_to_vec3(t: (f32, f32, f32)) -> Vec3 {
    Vec3::new(t.0, t.1, t.2)
}

/// Convert a Color to a Python tuple (r, g, b) with values 0.0-1.0
pub fn color_to_tuple(c: Color) -> (f32, f32, f32) {
    (c.r, c.g, c.b)
}

/// Convert a Python tuple (r, g, b) to Color
pub fn tuple_to_color(t: (f32, f32, f32)) -> Color {
    Color::new(t.0, t.1, t.2)
}

/// Convert a Color to a Python tuple (r, g, b) with values 0-255
pub fn color_to_rgb8(c: Color) -> (u8, u8, u8) {
    let rgb = c.to_rgb8();
    (rgb[0], rgb[1], rgb[2])
}

/// Convert a Python tuple (r, g, b) with values 0-255 to Color
pub fn rgb8_to_color(t: (u8, u8, u8)) -> Color {
    Color::from_rgb8(t.0, t.1, t.2)
}

/// Trait for Python-compatible coordinate extraction
pub trait ToPyCoord {
    fn to_py_coord(&self) -> (f32, f32, f32);
}

impl ToPyCoord for Vec3 {
    fn to_py_coord(&self) -> (f32, f32, f32) {
        (self.x, self.y, self.z)
    }
}

impl ToPyCoord for [f32; 3] {
    fn to_py_coord(&self) -> (f32, f32, f32) {
        (self[0], self[1], self[2])
    }
}

/// Trait for converting from Python coordinates
pub trait FromPyCoord: Sized {
    fn from_py_coord(coord: (f32, f32, f32)) -> Self;
}

impl FromPyCoord for Vec3 {
    fn from_py_coord(coord: (f32, f32, f32)) -> Self {
        Vec3::new(coord.0, coord.1, coord.2)
    }
}

impl FromPyCoord for [f32; 3] {
    fn from_py_coord(coord: (f32, f32, f32)) -> Self {
        [coord.0, coord.1, coord.2]
    }
}

/// Convert a slice of Vec3 to a list of tuples
pub fn coords_to_py_list(coords: &[Vec3]) -> Vec<(f32, f32, f32)> {
    coords.iter().map(|v| v.to_py_coord()).collect()
}

/// Convert a list of Python tuples to Vec<Vec3>
pub fn py_list_to_coords(coords: Vec<(f32, f32, f32)>) -> Vec<Vec3> {
    coords.into_iter().map(Vec3::from_py_coord).collect()
}

/// Extract coordinates from CoordSet as Python list
pub fn coordset_to_py_list(coordset: &pymol_mol::CoordSet) -> Vec<(f32, f32, f32)> {
    coordset.iter()
        .map(|v| (v.x, v.y, v.z))
        .collect()
}

/// Convert bounding box to Python tuple of tuples
pub fn bbox_to_py(min: Vec3, max: Vec3) -> ((f32, f32, f32), (f32, f32, f32)) {
    (vec3_to_tuple(min), vec3_to_tuple(max))
}

/// Helper to extract optional selection string
pub fn extract_selection(selection: Option<&str>) -> &str {
    selection.unwrap_or("all")
}

/// Helper to extract optional object name
pub fn extract_object_name(name: Option<&str>, default: &str) -> String {
    name.map(|s| s.to_string()).unwrap_or_else(|| default.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vec3_conversion() {
        let v = Vec3::new(1.0, 2.0, 3.0);
        let t = vec3_to_tuple(v);
        assert_eq!(t, (1.0, 2.0, 3.0));

        let v2 = tuple_to_vec3(t);
        assert_eq!(v2.x, 1.0);
        assert_eq!(v2.y, 2.0);
        assert_eq!(v2.z, 3.0);
    }

    #[test]
    fn test_color_conversion() {
        let c = Color::new(0.5, 0.25, 0.75);
        let t = color_to_tuple(c);
        assert!((t.0 - 0.5).abs() < 0.001);
        assert!((t.1 - 0.25).abs() < 0.001);
        assert!((t.2 - 0.75).abs() < 0.001);
    }

    #[test]
    fn test_rgb8_conversion() {
        let c = Color::from_rgb8(255, 128, 0);
        let rgb = color_to_rgb8(c);
        assert_eq!(rgb, (255, 128, 0));
    }
}
