//! PyMOL standard color index → RGB mapping table.
//!
//! PyMOL assigns fixed integer indices to its built-in colors.  This table
//! lets us convert those indices to RGB values during PSE import so that we
//! can register them in our own [`NamedColors`] registry.

/// (pymol_index, r, g, b)  — values taken from PyMOL's `Color.cpp`.
pub static PYMOL_COLORS: &[(i64, f32, f32, f32)] = &[
    // Basic colors (0-7)
    (0, 0.2, 1.0, 0.2),      // default (carbon green)
    (1, 1.0, 0.0, 0.0),      // red
    (2, 0.0, 1.0, 0.0),      // green
    (3, 0.0, 0.0, 1.0),      // blue
    (4, 1.0, 1.0, 0.0),      // yellow
    (5, 1.0, 0.0, 1.0),      // magenta
    (6, 0.0, 1.0, 1.0),      // cyan
    (7, 1.0, 0.5, 0.5),      // salmon
    (8, 0.6, 0.6, 0.1),      // lime
    (9, 1.0, 0.65, 0.85),    // slate → pink
    (10, 1.0, 1.0, 1.0),     // white
    (11, 0.0, 0.0, 0.0),     // black (actually dash)
    (12, 0.5, 0.5, 0.5),     // grey/gray

    // PyMOL extended named colors
    (13, 1.0, 0.5, 0.0),     // orange
    (14, 1.0, 1.0, 0.0),     // brightyellow (same as yellow)
    (15, 0.65, 0.32, 0.17),  // brown
    (16, 0.2, 0.6, 0.2),     // forest
    (17, 0.75, 0.0, 0.75),   // purple
    (18, 0.0, 0.75, 0.75),   // teal
    (19, 0.0, 0.0, 0.5),     // navy
    (20, 0.5, 0.0, 0.0),     // firebrick / maroon
    (21, 0.5, 0.5, 0.0),     // olive
    (22, 0.85, 0.85, 1.0),   // lightblue
    (23, 0.8, 1.0, 0.8),     // lightgreen
    (24, 0.6, 0.6, 1.0),     // slate
    (25, 1.0, 0.0, 0.5),     // hotpink
    (26, 0.2, 1.0, 0.2),     // carbon (limegreen)
    (27, 0.2, 0.2, 1.0),     // nitrogen
    (28, 1.0, 0.3, 0.3),     // oxygen
    (29, 0.9, 0.775, 0.25),  // sulfur
    (30, 0.9, 0.9, 0.9),     // hydrogen

    // PyMOL color spectrum: greens, cyans, blues (50-99)
    (50, 0.0, 0.0, 1.0),     // tv_blue
    (51, 0.3, 0.3, 1.0),     // tv_green → actually tv_blue variant
    (52, 0.2, 0.6, 0.2),     // smudge
    (53, 1.0, 0.2, 0.2),     // tv_red
    (54, 0.0, 0.9, 0.04),    // tv_green
    (55, 1.0, 1.0, 0.2),     // tv_yellow
    (56, 1.0, 0.5, 0.0),     // tv_orange

    // Deep colors
    (100, 0.1, 0.1, 0.6),    // deepblue
    (101, 0.6, 0.1, 0.6),    // deeppurple
    (102, 0.1, 0.6, 0.1),    // deepgreen
    (103, 0.1, 0.6, 0.6),    // deepteal
    (104, 0.6, 0.1, 0.1),    // deepred

    // Light colors
    (150, 0.75, 0.75, 1.0),  // lightblue
    (151, 1.0, 0.75, 0.87),  // lightpink
    (152, 0.75, 1.0, 0.75),  // lightgreen
    (153, 0.75, 1.0, 1.0),   // lightcyan
    (154, 1.0, 0.87, 0.37),  // lightyellow
    (155, 1.0, 0.75, 0.5),   // lightorange
    (156, 0.87, 0.75, 1.0),  // lightpurple
    (157, 1.0, 1.0, 0.75),   // paleyellow
    (158, 0.65, 0.9, 0.65),  // limegreen
    (159, 0.6, 0.6, 1.0),    // slate

    // Warm colors
    (200, 1.0, 0.6, 0.6),    // lightsalmon / warmpink
    (201, 0.85, 0.2, 0.5),   // raspberry
    (202, 0.6, 0.0, 0.0),    // darkred
    (203, 1.0, 0.5, 0.5),    // salmon
    (204, 0.85, 0.65, 0.125),// dirtyviolet → wheat
    (205, 0.55, 0.25, 0.6),  // violet / deeppurple
    (206, 0.5, 0.5, 0.0),    // olive
    (207, 1.0, 0.85, 0.7),   // wheat

    // Chain colors (auto coloring)
    (5000, 0.2, 1.0, 0.2),   // auto color 0 (green)
    (5001, 0.0, 1.0, 1.0),   // auto color 1 (cyan)
    (5002, 1.0, 0.5, 1.0),   // auto color 2 (lightmagenta)
    (5003, 1.0, 1.0, 0.0),   // auto color 3 (yellow)
    (5004, 0.65, 0.32, 0.17),// auto color 4 (salmon→brown)
    (5005, 0.5, 0.5, 0.0),   // auto color 5 (olive)
    (5006, 0.6, 0.6, 0.1),   // auto color 6 (smudge)
    (5007, 0.0, 0.46, 0.74), // auto color 7 (steel blue)
    (5008, 1.0, 0.65, 0.85), // auto color 8 (pink)
    (5009, 0.55, 0.25, 0.6), // auto color 9 (violet)

    // Element colors (negative-offset in PyMOL, but stored as positive in PSE)
    // These are typically encoded as negative indices in the atom color field
    // and handled by the ByElement color scheme, not this table.

    // Gray scale
    (5050, 0.0, 0.0, 0.0),   // gray0
    (5051, 0.1, 0.1, 0.1),   // gray10
    (5052, 0.2, 0.2, 0.2),   // gray20
    (5053, 0.3, 0.3, 0.3),   // gray30
    (5054, 0.4, 0.4, 0.4),   // gray40
    (5055, 0.5, 0.5, 0.5),   // gray50
    (5056, 0.6, 0.6, 0.6),   // gray60
    (5057, 0.7, 0.7, 0.7),   // gray70
    (5058, 0.8, 0.8, 0.8),   // gray80
    (5059, 0.9, 0.9, 0.9),   // gray90
];

use ahash::AHashMap;
use std::sync::LazyLock;

/// Fast lookup: PyMOL color index → (R, G, B)
static PYMOL_COLOR_MAP: LazyLock<AHashMap<i64, (f32, f32, f32)>> = LazyLock::new(|| {
    PYMOL_COLORS.iter().map(|&(idx, r, g, b)| (idx, (r, g, b))).collect()
});

/// Look up a PyMOL color index and return its RGB values.
pub fn pymol_color_rgb(index: i64) -> Option<(f32, f32, f32)> {
    PYMOL_COLOR_MAP.get(&index).copied()
}
