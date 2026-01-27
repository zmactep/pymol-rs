//! Color module for Python
//!
//! Provides color types, named colors, and color ramps.

use pyo3::prelude::*;
use pymol_color::{Color, ColorRamp, NamedColors};

/// RGB Color with floating-point components (0.0-1.0)
#[pyclass(name = "Color")]
#[derive(Debug, Clone, Copy)]
pub struct PyColor {
    pub(crate) inner: Color,
}

impl From<Color> for PyColor {
    fn from(c: Color) -> Self {
        PyColor { inner: c }
    }
}

impl From<PyColor> for Color {
    fn from(c: PyColor) -> Self {
        c.inner
    }
}

#[pymethods]
impl PyColor {
    /// Create a new color from RGB components (0.0-1.0)
    #[new]
    #[pyo3(signature = (r=0.0, g=0.0, b=0.0))]
    fn new(r: f32, g: f32, b: f32) -> Self {
        PyColor {
            inner: Color::new(r, g, b),
        }
    }

    /// Create a color from RGB components (0-255)
    #[staticmethod]
    fn from_rgb8(r: u8, g: u8, b: u8) -> Self {
        PyColor {
            inner: Color::from_rgb8(r, g, b),
        }
    }

    /// Create a color from a hex string (e.g., "#ff0000" or "ff0000")
    #[staticmethod]
    fn from_hex(hex: &str) -> PyResult<Self> {
        Color::from_hex(hex)
            .map(|c| PyColor { inner: c })
            .ok_or_else(|| {
                pyo3::exceptions::PyValueError::new_err(format!("Invalid hex color: {}", hex))
            })
    }

    /// Red component (0.0-1.0)
    #[getter]
    fn r(&self) -> f32 {
        self.inner.r
    }

    /// Green component (0.0-1.0)
    #[getter]
    fn g(&self) -> f32 {
        self.inner.g
    }

    /// Blue component (0.0-1.0)
    #[getter]
    fn b(&self) -> f32 {
        self.inner.b
    }

    /// Get as tuple (r, g, b) with values 0.0-1.0
    fn to_tuple(&self) -> (f32, f32, f32) {
        (self.inner.r, self.inner.g, self.inner.b)
    }

    /// Get as tuple (r, g, b) with values 0-255
    fn to_rgb8(&self) -> (u8, u8, u8) {
        let rgb = self.inner.to_rgb8();
        (rgb[0], rgb[1], rgb[2])
    }

    /// Get as RGBA tuple with alpha=1.0
    fn to_rgba(&self) -> (f32, f32, f32, f32) {
        let rgba = self.inner.to_rgba(1.0);
        (rgba[0], rgba[1], rgba[2], rgba[3])
    }

    /// Linearly interpolate between this color and another
    fn lerp(&self, other: &PyColor, t: f32) -> PyColor {
        PyColor {
            inner: self.inner.lerp(&other.inner, t),
        }
    }

    fn __repr__(&self) -> String {
        format!("Color({:.3}, {:.3}, {:.3})", self.inner.r, self.inner.g, self.inner.b)
    }

    fn __str__(&self) -> String {
        let rgb = self.inner.to_rgb8();
        format!("#{:02x}{:02x}{:02x}", rgb[0], rgb[1], rgb[2])
    }

    fn __eq__(&self, other: &PyColor) -> bool {
        (self.inner.r - other.inner.r).abs() < 0.001
            && (self.inner.g - other.inner.g).abs() < 0.001
            && (self.inner.b - other.inner.b).abs() < 0.001
    }
}

// Common color constants
impl PyColor {
    /// Black color
    pub const BLACK: PyColor = PyColor {
        inner: Color { r: 0.0, g: 0.0, b: 0.0 },
    };
    /// White color
    pub const WHITE: PyColor = PyColor {
        inner: Color { r: 1.0, g: 1.0, b: 1.0 },
    };
    /// Red color
    pub const RED: PyColor = PyColor {
        inner: Color { r: 1.0, g: 0.0, b: 0.0 },
    };
    /// Green color
    pub const GREEN: PyColor = PyColor {
        inner: Color { r: 0.0, g: 1.0, b: 0.0 },
    };
    /// Blue color
    pub const BLUE: PyColor = PyColor {
        inner: Color { r: 0.0, g: 0.0, b: 1.0 },
    };
}

/// Color ramp for continuous color mapping
#[pyclass(name = "ColorRamp")]
#[derive(Clone)]
pub struct PyColorRamp {
    inner: ColorRamp,
}

#[pymethods]
impl PyColorRamp {
    /// Create a new empty color ramp
    #[new]
    fn new() -> Self {
        PyColorRamp {
            inner: ColorRamp::new("custom"),
        }
    }

    /// Create a blue-white-red ramp (for negative to positive values)
    #[staticmethod]
    fn blue_white_red() -> Self {
        PyColorRamp {
            inner: ColorRamp::blue_white_red(),
        }
    }

    /// Create a rainbow ramp
    #[staticmethod]
    fn rainbow() -> Self {
        PyColorRamp {
            inner: ColorRamp::rainbow(),
        }
    }

    /// Create a grayscale ramp
    #[staticmethod]
    fn grayscale() -> Self {
        PyColorRamp {
            inner: ColorRamp::grayscale(),
        }
    }

    /// Create a hot ramp (black -> red -> yellow -> white)
    #[staticmethod]
    fn hot() -> Self {
        PyColorRamp {
            inner: ColorRamp::hot(),
        }
    }

    /// Add a color point to the ramp
    fn add_point(&mut self, position: f32, color: &PyColor) {
        self.inner.add_point(position, color.inner);
    }

    /// Get the color at a position (0.0-1.0)
    fn get_color(&self, position: f32) -> PyColor {
        PyColor {
            inner: self.inner.get_color(position),
        }
    }

    /// Map a value to a color within a range
    #[pyo3(signature = (value, min_val=0.0, max_val=1.0))]
    fn map_value(&self, value: f32, min_val: f32, max_val: f32) -> PyColor {
        let normalized = if (max_val - min_val).abs() < 0.0001 {
            0.5
        } else {
            (value - min_val) / (max_val - min_val)
        };
        self.get_color(normalized.clamp(0.0, 1.0))
    }

    fn __repr__(&self) -> String {
        "ColorRamp(...)".to_string()
    }
}

/// Named colors registry
#[pyclass(name = "NamedColors")]
pub struct PyNamedColors {
    inner: NamedColors,
}

#[pymethods]
impl PyNamedColors {
    /// Create a new named colors registry with default colors
    #[new]
    fn new() -> Self {
        PyNamedColors {
            inner: NamedColors::default(),
        }
    }

    /// Get a color by name
    fn get(&self, name: &str) -> Option<PyColor> {
        self.inner
            .get_by_name(name)
            .map(|(_, c)| PyColor { inner: c })
    }

    /// Get a color by index
    fn get_by_index(&self, index: u32) -> Option<PyColor> {
        self.inner
            .get_by_index(index)
            .map(|c| PyColor { inner: c })
    }

    /// Register a new named color
    fn register(&mut self, name: &str, color: &PyColor) -> u32 {
        self.inner.register(name, color.inner)
    }

    fn __getitem__(&self, name: &str) -> PyResult<PyColor> {
        self.get(name).ok_or_else(|| {
            pyo3::exceptions::PyKeyError::new_err(format!("Unknown color: {}", name))
        })
    }

    fn __repr__(&self) -> String {
        "NamedColors(...)".to_string()
    }
}

// Module-level helper functions

/// Get a named color
#[pyfunction]
pub fn get_color(name: &str) -> Option<PyColor> {
    let colors = NamedColors::default();
    colors.get_by_name(name).map(|(_, c)| PyColor { inner: c })
}

/// Get element color (CPK coloring)
#[pyfunction]
pub fn element_color(symbol: &str) -> PyColor {
    use pymol_color::ElementColors;
    let colors = ElementColors::default();
    PyColor {
        inner: colors.get_by_symbol(symbol),
    }
}

/// Get chain color
#[pyfunction]
pub fn chain_color(chain_id: &str) -> PyColor {
    use pymol_color::ChainColors;
    PyColor {
        inner: ChainColors::get(chain_id),
    }
}

/// Register the color submodule
pub fn register_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyColor>()?;
    m.add_class::<PyColorRamp>()?;
    m.add_class::<PyNamedColors>()?;
    m.add_function(wrap_pyfunction!(get_color, m)?)?;
    m.add_function(wrap_pyfunction!(element_color, m)?)?;
    m.add_function(wrap_pyfunction!(chain_color, m)?)?;

    // Add color constants
    m.add("BLACK", PyColor::BLACK)?;
    m.add("WHITE", PyColor::WHITE)?;
    m.add("RED", PyColor::RED)?;
    m.add("GREEN", PyColor::GREEN)?;
    m.add("BLUE", PyColor::BLUE)?;

    Ok(())
}
