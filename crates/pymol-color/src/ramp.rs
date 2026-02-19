//! Color ramps for continuous coloring

use serde::{Deserialize, Serialize};

use crate::Color;

/// Type of color ramp interpolation
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum RampType {
    /// Linear interpolation in RGB space
    Linear,
    /// Smooth step interpolation
    Smooth,
}

/// A color ramp for mapping values to colors
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ColorRamp {
    /// Name of the ramp
    pub name: String,
    /// Control points: (value, color) pairs, sorted by value
    points: Vec<(f32, Color)>,
    /// Interpolation type
    pub ramp_type: RampType,
}

impl ColorRamp {
    /// Create a new color ramp
    pub fn new(name: impl Into<String>) -> Self {
        ColorRamp {
            name: name.into(),
            points: Vec::new(),
            ramp_type: RampType::Linear,
        }
    }

    /// Add a control point
    pub fn add_point(&mut self, value: f32, color: Color) -> &mut Self {
        self.points.push((value, color));
        self.points.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        self
    }

    /// Set the interpolation type
    pub fn with_type(mut self, ramp_type: RampType) -> Self {
        self.ramp_type = ramp_type;
        self
    }

    /// Get the color for a value
    pub fn get_color(&self, value: f32) -> Color {
        if self.points.is_empty() {
            return Color::WHITE;
        }

        if self.points.len() == 1 {
            return self.points[0].1;
        }

        // Find the two control points that bracket the value
        let mut lower_idx = 0;
        let mut upper_idx = self.points.len() - 1;

        for (i, &(v, _)) in self.points.iter().enumerate() {
            if v <= value {
                lower_idx = i;
            }
            if v >= value {
                upper_idx = i;
                break;
            }
        }

        if lower_idx == upper_idx {
            return self.points[lower_idx].1;
        }

        let (v0, c0) = self.points[lower_idx];
        let (v1, c1) = self.points[upper_idx];

        let t = if (v1 - v0).abs() < f32::EPSILON {
            0.5
        } else {
            (value - v0) / (v1 - v0)
        };

        let t = match self.ramp_type {
            RampType::Linear => t,
            RampType::Smooth => t * t * (3.0 - 2.0 * t), // smoothstep
        };

        c0.lerp(&c1, t)
    }

    /// Create a blue-white-red ramp (for B-factors, etc.)
    pub fn blue_white_red() -> Self {
        let mut ramp = ColorRamp::new("bwr");
        ramp.add_point(0.0, Color::BLUE)
            .add_point(0.5, Color::WHITE)
            .add_point(1.0, Color::RED);
        ramp
    }

    /// Create a rainbow ramp
    pub fn rainbow() -> Self {
        let mut ramp = ColorRamp::new("rainbow");
        ramp.add_point(0.0, Color::BLUE)
            .add_point(0.25, Color::CYAN)
            .add_point(0.5, Color::GREEN)
            .add_point(0.75, Color::YELLOW)
            .add_point(1.0, Color::RED);
        ramp
    }

    /// Create a grayscale ramp
    pub fn grayscale() -> Self {
        let mut ramp = ColorRamp::new("grayscale");
        ramp.add_point(0.0, Color::BLACK)
            .add_point(1.0, Color::WHITE);
        ramp
    }

    /// Create a hot ramp (black -> red -> yellow -> white)
    pub fn hot() -> Self {
        let mut ramp = ColorRamp::new("hot");
        ramp.add_point(0.0, Color::BLACK)
            .add_point(0.33, Color::RED)
            .add_point(0.66, Color::YELLOW)
            .add_point(1.0, Color::WHITE);
        ramp
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_color_ramp() {
        let ramp = ColorRamp::blue_white_red();

        let c0 = ramp.get_color(0.0);
        assert_eq!(c0, Color::BLUE);

        let c1 = ramp.get_color(1.0);
        assert_eq!(c1, Color::RED);

        let c_mid = ramp.get_color(0.5);
        assert_eq!(c_mid, Color::WHITE);
    }
}
