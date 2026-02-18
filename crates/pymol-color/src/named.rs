//! Named color registry

use ahash::AHashMap;
use crate::Color;

/// Registry of named colors
#[derive(Debug)]
pub struct NamedColors {
    colors: Vec<Color>,
    by_name: AHashMap<String, u32>,
}

impl NamedColors {
    /// Create a new named color registry with PyMOL's default colors
    pub fn new() -> Self {
        let mut registry = NamedColors {
            colors: Vec::with_capacity(256),
            by_name: AHashMap::new(),
        };

        // Register PyMOL's standard colors
        registry.register_defaults();
        registry
    }

    fn register_defaults(&mut self) {
        // Basic colors
        self.register("white", Color::new(1.0, 1.0, 1.0));
        self.register("black", Color::new(0.0, 0.0, 0.0));
        self.register("red", Color::new(1.0, 0.0, 0.0));
        self.register("green", Color::new(0.0, 1.0, 0.0));
        self.register("blue", Color::new(0.0, 0.0, 1.0));
        self.register("yellow", Color::new(1.0, 1.0, 0.0));
        self.register("cyan", Color::new(0.0, 1.0, 1.0));
        self.register("magenta", Color::new(1.0, 0.0, 1.0));

        // Gray scale
        self.register("gray", Color::new(0.5, 0.5, 0.5));
        self.register("grey", Color::new(0.5, 0.5, 0.5));
        self.register("gray10", Color::new(0.1, 0.1, 0.1));
        self.register("gray20", Color::new(0.2, 0.2, 0.2));
        self.register("gray30", Color::new(0.3, 0.3, 0.3));
        self.register("gray40", Color::new(0.4, 0.4, 0.4));
        self.register("gray50", Color::new(0.5, 0.5, 0.5));
        self.register("gray60", Color::new(0.6, 0.6, 0.6));
        self.register("gray70", Color::new(0.7, 0.7, 0.7));
        self.register("gray80", Color::new(0.8, 0.8, 0.8));
        self.register("gray90", Color::new(0.9, 0.9, 0.9));

        // PyMOL specific colors
        self.register("carbon", Color::new(0.2, 1.0, 0.2));
        self.register("nitrogen", Color::new(0.2, 0.2, 1.0));
        self.register("oxygen", Color::new(1.0, 0.3, 0.3));
        self.register("hydrogen", Color::new(0.9, 0.9, 0.9));
        self.register("sulfur", Color::new(0.9, 0.775, 0.25));

        // Extended palette
        self.register("orange", Color::new(1.0, 0.5, 0.0));
        self.register("pink", Color::new(1.0, 0.65, 0.85));
        self.register("purple", Color::new(0.75, 0.0, 0.75));
        self.register("brown", Color::new(0.65, 0.32, 0.17));
        self.register("salmon", Color::new(1.0, 0.6, 0.6));
        self.register("lime", Color::new(0.5, 1.0, 0.5));
        self.register("slate", Color::new(0.5, 0.5, 1.0));
        self.register("hotpink", Color::new(1.0, 0.0, 0.5));
        self.register("teal", Color::new(0.0, 0.75, 0.75));
        self.register("olive", Color::new(0.77, 0.7, 0.0));
        self.register("marine", Color::new(0.0, 0.5, 1.0));
        self.register("forest", Color::new(0.2, 0.6, 0.2));
        self.register("firebrick", Color::new(0.7, 0.13, 0.13));
        self.register("chocolate", Color::new(0.55, 0.27, 0.07));
        self.register("wheat", Color::new(0.99, 0.82, 0.65));
        self.register("violet", Color::new(1.0, 0.5, 1.0));
        self.register("lightblue", Color::new(0.75, 0.75, 1.0));
        self.register("lightgreen", Color::new(0.75, 1.0, 0.75));
        self.register("lightorange", Color::new(1.0, 0.8, 0.5));
        self.register("lightpink", Color::new(1.0, 0.75, 0.87));
        self.register("palecyan", Color::new(0.8, 1.0, 1.0));
        self.register("paleyellow", Color::new(1.0, 1.0, 0.8));
        self.register("palegreen", Color::new(0.65, 0.9, 0.65));
        self.register("aquamarine", Color::new(0.5, 1.0, 0.83));
        self.register("deepteal", Color::new(0.1, 0.6, 0.6));
        self.register("deeppurple", Color::new(0.6, 0.1, 0.6));
        self.register("deepolive", Color::new(0.6, 0.6, 0.1));
        self.register("deepsalmon", Color::new(1.0, 0.42, 0.42));
        self.register("tv_red", Color::new(1.0, 0.2, 0.2));
        self.register("tv_green", Color::new(0.2, 1.0, 0.2));
        self.register("tv_blue", Color::new(0.3, 0.3, 1.0));
        self.register("tv_yellow", Color::new(1.0, 1.0, 0.2));
        self.register("tv_orange", Color::new(1.0, 0.55, 0.15));
    }

    /// Register a new named color
    pub fn register(&mut self, name: &str, color: Color) -> u32 {
        let index = self.colors.len() as u32;
        self.colors.push(color);
        self.by_name.insert(name.to_lowercase(), index);
        index
    }

    /// Get a color by name
    pub fn get_by_name(&self, name: &str) -> Option<(u32, Color)> {
        self.by_name
            .get(&name.to_lowercase())
            .map(|&idx| (idx, self.colors[idx as usize]))
    }

    /// Get a color by index
    pub fn get_by_index(&self, index: u32) -> Option<Color> {
        self.colors.get(index as usize).copied()
    }

    /// Get the number of registered colors
    pub fn len(&self) -> usize {
        self.colors.len()
    }

    /// Check if the registry is empty
    pub fn is_empty(&self) -> bool {
        self.colors.is_empty()
    }

    /// Iterate over all named colors
    /// Get all registered color names
    pub fn names(&self) -> Vec<&str> {
        let mut names: Vec<&str> = self.by_name.keys().map(|s| s.as_str()).collect();
        names.sort();
        names
    }

    pub fn iter(&self) -> impl Iterator<Item = (&str, Color)> + '_ {
        self.by_name
            .iter()
            .map(|(name, &idx)| (name.as_str(), self.colors[idx as usize]))
    }
}

impl Default for NamedColors {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_named_colors() {
        let colors = NamedColors::new();

        let (_, red) = colors.get_by_name("red").unwrap();
        assert_eq!(red.r, 1.0);
        assert_eq!(red.g, 0.0);
        assert_eq!(red.b, 0.0);

        // Test case insensitivity
        let (_, red2) = colors.get_by_name("RED").unwrap();
        assert_eq!(red, red2);
    }
}
