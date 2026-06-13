//! Color gradients for continuous coloring (spectrum, b-factor, etc.)

use serde::{Deserialize, Serialize};

use crate::Color;

/// A multi-stop color gradient mapping t ∈ [0, 1] to a Color.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Gradient {
    stops: Vec<(f32, Color)>,
}

impl Gradient {
    /// Create a gradient from a list of (position, color) stops.
    /// Positions should be in [0.0, 1.0] and will be sorted.
    pub fn new(mut stops: Vec<(f32, Color)>) -> Self {
        stops.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        Gradient { stops }
    }

    /// Sample the gradient at position t ∈ [0, 1].
    pub fn sample(&self, t: f32) -> Color {
        if self.stops.is_empty() {
            return Color::WHITE;
        }
        if self.stops.len() == 1 {
            return self.stops[0].1;
        }

        let t = t.clamp(0.0, 1.0);

        // Find bracketing stops
        if t <= self.stops[0].0 {
            return self.stops[0].1;
        }
        if t >= self.stops[self.stops.len() - 1].0 {
            return self.stops[self.stops.len() - 1].1;
        }

        for i in 0..self.stops.len() - 1 {
            let (t0, c0) = self.stops[i];
            let (t1, c1) = self.stops[i + 1];
            if t >= t0 && t <= t1 {
                let span = t1 - t0;
                if span < f32::EPSILON {
                    return c0;
                }
                let s = (t - t0) / span;
                return c0.lerp(&c1, s);
            }
        }

        self.stops[self.stops.len() - 1].1
    }

    // --- Predefined gradients -----------------------------------------------

    /// Rainbow: blue → cyan → green → yellow → red
    pub fn rainbow() -> Self {
        Gradient::new(vec![
            (0.0, Color::BLUE),
            (0.25, Color::CYAN),
            (0.5, Color::GREEN),
            (0.75, Color::YELLOW),
            (1.0, Color::RED),
        ])
    }

    /// Blue → white → red (for B-factors, diverging data)
    pub fn blue_white_red() -> Self {
        Gradient::new(vec![
            (0.0, Color::BLUE),
            (0.5, Color::WHITE),
            (1.0, Color::RED),
        ])
    }

    /// Viridis perceptual colormap (16-stop approximation)
    pub fn viridis() -> Self {
        Gradient::new(vec![
            (0.000, Color::new(0.267, 0.004, 0.329)),
            (0.067, Color::new(0.283, 0.141, 0.458)),
            (0.133, Color::new(0.254, 0.265, 0.530)),
            (0.200, Color::new(0.207, 0.372, 0.553)),
            (0.267, Color::new(0.164, 0.471, 0.558)),
            (0.333, Color::new(0.128, 0.567, 0.551)),
            (0.400, Color::new(0.135, 0.659, 0.518)),
            (0.467, Color::new(0.209, 0.747, 0.452)),
            (0.533, Color::new(0.328, 0.816, 0.374)),
            (0.600, Color::new(0.477, 0.870, 0.290)),
            (0.667, Color::new(0.648, 0.905, 0.210)),
            (0.733, Color::new(0.815, 0.922, 0.165)),
            (0.800, Color::new(0.937, 0.923, 0.176)),
            (0.867, Color::new(0.993, 0.906, 0.144)),
            (0.933, Color::new(0.998, 0.853, 0.109)),
            (1.000, Color::new(0.993, 0.906, 0.144)),
        ])
    }

    /// Plasma perceptual colormap (10-stop approximation)
    pub fn plasma() -> Self {
        Gradient::new(vec![
            (0.000, Color::new(0.050, 0.030, 0.528)),
            (0.111, Color::new(0.294, 0.012, 0.631)),
            (0.222, Color::new(0.492, 0.012, 0.658)),
            (0.333, Color::new(0.658, 0.134, 0.588)),
            (0.444, Color::new(0.797, 0.280, 0.470)),
            (0.556, Color::new(0.899, 0.424, 0.360)),
            (0.667, Color::new(0.964, 0.575, 0.258)),
            (0.778, Color::new(0.993, 0.738, 0.154)),
            (0.889, Color::new(0.973, 0.907, 0.131)),
            (1.000, Color::new(0.940, 0.975, 0.131)),
        ])
    }

    /// Inferno perceptual colormap (10-stop approximation)
    pub fn inferno() -> Self {
        Gradient::new(vec![
            (0.000, Color::new(0.001, 0.000, 0.014)),
            (0.111, Color::new(0.121, 0.047, 0.282)),
            (0.222, Color::new(0.316, 0.072, 0.485)),
            (0.333, Color::new(0.526, 0.120, 0.507)),
            (0.444, Color::new(0.735, 0.215, 0.330)),
            (0.556, Color::new(0.891, 0.350, 0.140)),
            (0.667, Color::new(0.975, 0.549, 0.040)),
            (0.778, Color::new(0.993, 0.765, 0.156)),
            (0.889, Color::new(0.955, 0.935, 0.389)),
            (1.000, Color::new(0.988, 0.998, 0.645)),
        ])
    }

    /// Magma perceptual colormap (10-stop approximation)
    pub fn magma() -> Self {
        Gradient::new(vec![
            (0.000, Color::new(0.001, 0.000, 0.014)),
            (0.111, Color::new(0.113, 0.065, 0.276)),
            (0.222, Color::new(0.302, 0.084, 0.498)),
            (0.333, Color::new(0.506, 0.145, 0.506)),
            (0.444, Color::new(0.710, 0.212, 0.476)),
            (0.556, Color::new(0.899, 0.323, 0.390)),
            (0.667, Color::new(0.987, 0.510, 0.374)),
            (0.778, Color::new(0.997, 0.707, 0.500)),
            (0.889, Color::new(0.994, 0.877, 0.718)),
            (1.000, Color::new(0.987, 0.991, 0.750)),
        ])
    }

    /// Cool-warm diverging colormap (blue → white → red, perceptual)
    pub fn coolwarm() -> Self {
        Gradient::new(vec![
            (0.000, Color::new(0.230, 0.299, 0.754)),
            (0.250, Color::new(0.553, 0.627, 0.934)),
            (0.500, Color::new(0.866, 0.866, 0.866)),
            (0.750, Color::new(0.906, 0.533, 0.451)),
            (1.000, Color::new(0.706, 0.016, 0.150)),
        ])
    }

    /// Grayscale: black → white
    pub fn grayscale() -> Self {
        Gradient::new(vec![(0.0, Color::BLACK), (1.0, Color::WHITE)])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sample_endpoints() {
        let g = Gradient::rainbow();
        let blue = g.sample(0.0);
        assert!(blue.b > 0.9 && blue.r < 0.1);
        let red = g.sample(1.0);
        assert!(red.r > 0.9 && red.b < 0.1);
    }

    #[test]
    fn test_sample_midpoint() {
        let g = Gradient::rainbow();
        let green = g.sample(0.5);
        assert!(green.g > 0.9 && green.r < 0.1 && green.b < 0.1);
    }

    #[test]
    fn test_sample_clamping() {
        let g = Gradient::rainbow();
        assert_eq!(g.sample(-1.0), g.sample(0.0));
        assert_eq!(g.sample(2.0), g.sample(1.0));
    }

    #[test]
    fn test_blue_white_red() {
        let g = Gradient::blue_white_red();
        let mid = g.sample(0.5);
        assert!((mid.r - 1.0).abs() < 0.01);
        assert!((mid.g - 1.0).abs() < 0.01);
        assert!((mid.b - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_empty_gradient() {
        let g = Gradient::new(vec![]);
        assert_eq!(g.sample(0.5), Color::WHITE);
    }

    #[test]
    fn test_single_stop() {
        let g = Gradient::new(vec![(0.5, Color::RED)]);
        assert_eq!(g.sample(0.0), Color::RED);
        assert_eq!(g.sample(1.0), Color::RED);
    }
}
