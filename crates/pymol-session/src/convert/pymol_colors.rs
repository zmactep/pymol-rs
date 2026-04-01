//! PyMOL standard color index → RGB mapping table.
//!
//! PyMOL assigns sequential integer indices to its built-in colors during
//! initialization in `Color.cpp`.  The indices are determined by the
//! registration order, starting from 0.  This table lets us convert those
//! indices to RGB values during PSE import.

/// (pymol_index, r, g, b) — values and order taken from PyMOL's `Color.cpp` `ColorInit`.
///
/// Only the first 57 named colors are listed here (indices 0–56).
/// After that come grey00–grey99 (indices 57–156), lightmagenta (157),
/// then large generated spectrum ranges (s000–s999, r000–r999, etc.).
pub static PYMOL_COLORS: &[(i64, f32, f32, f32)] = &[
    // Named colors (registered sequentially in ColorInit)
    (0,  1.0,   1.0,   1.0),    // white
    (1,  0.0,   0.0,   0.0),    // black
    (2,  0.0,   0.0,   1.0),    // blue
    (3,  0.0,   1.0,   0.0),    // green
    (4,  1.0,   0.0,   0.0),    // red
    (5,  0.0,   1.0,   1.0),    // cyan
    (6,  1.0,   1.0,   0.0),    // yellow
    (7,  1.0,   1.0,   0.0),    // dash
    (8,  1.0,   0.0,   1.0),    // magenta
    (9,  1.0,   0.6,   0.6),    // salmon
    (10, 0.5,   1.0,   0.5),    // lime
    (11, 0.5,   0.5,   1.0),    // slate
    (12, 1.0,   0.0,   0.5),    // hotpink
    (13, 1.0,   0.5,   0.0),    // orange
    (14, 0.5,   1.0,   0.0),    // chartreuse
    (15, 0.0,   1.0,   0.5),    // limegreen
    (16, 0.5,   0.0,   1.0),    // purpleblue
    (17, 0.0,   0.5,   1.0),    // marine
    (18, 0.77,  0.7,   0.0),    // olive
    (19, 0.75,  0.0,   0.75),   // purple
    (20, 0.0,   0.75,  0.75),   // teal
    (21, 0.6,   0.2,   0.2),    // ruby
    (22, 0.2,   0.6,   0.2),    // forest
    (23, 0.25,  0.25,  0.65),   // deepblue
    (24, 0.5,   0.5,   0.5),    // grey
    (25, 0.5,   0.5,   0.5),    // gray
    (26, 0.2,   1.0,   0.2),    // carbon
    (27, 0.2,   0.2,   1.0),    // nitrogen
    (28, 1.0,   0.3,   0.3),    // oxygen
    (29, 0.9,   0.9,   0.9),    // hydrogen
    (30, 1.0,   0.7,   0.2),    // brightorange
    (31, 0.9,   0.775, 0.25),   // sulfur
    (32, 1.0,   0.2,   0.2),    // tv_red
    (33, 0.2,   1.0,   0.2),    // tv_green
    (34, 0.3,   0.3,   1.0),    // tv_blue
    (35, 1.0,   1.0,   0.2),    // tv_yellow
    (36, 1.0,   0.87,  0.37),   // yelloworange
    (37, 1.0,   0.55,  0.15),   // tv_orange
    (38, 0.1,   0.1,   1.0),    // br0
    (39, 0.2,   0.1,   0.9),    // br1
    (40, 0.3,   0.1,   0.8),    // br2
    (41, 0.4,   0.1,   0.7),    // br3
    (42, 0.5,   0.1,   0.6),    // br4
    (43, 0.6,   0.1,   0.5),    // br5
    (44, 0.7,   0.1,   0.4),    // br6
    (45, 0.8,   0.1,   0.3),    // br7
    (46, 0.9,   0.1,   0.2),    // br8
    (47, 1.0,   0.1,   0.1),    // br9
    (48, 1.0,   0.65,  0.85),   // pink
    (49, 0.698, 0.13,  0.13),   // firebrick
    (50, 0.555, 0.222, 0.111),  // chocolate
    (51, 0.65,  0.32,  0.17),   // brown
    (52, 0.99,  0.82,  0.65),   // wheat
    (53, 1.0,   0.5,   1.0),    // violet

    // grey00–grey99 (indices 54–153)
    // grey00 = index 54, grey01 = 55, ..., grey99 = 153
    // formula: grey_N → (N/99.0, N/99.0, N/99.0)
    // (only a few key entries; the full range is generated at runtime via the formula)

    // lightmagenta (index 154 after 100 greys = 54+100)
    (154, 1.0, 0.2, 0.8),       // lightmagenta
];

use ahash::AHashMap;
use std::sync::LazyLock;

/// Fast lookup: PyMOL color index → (R, G, B)
static PYMOL_COLOR_MAP: LazyLock<AHashMap<i64, (f32, f32, f32)>> = LazyLock::new(|| {
    let mut map: AHashMap<i64, (f32, f32, f32)> = PYMOL_COLORS
        .iter()
        .map(|&(idx, r, g, b)| (idx, (r, g, b)))
        .collect();

    // Generate grey00–grey99 (indices 54–153)
    for a in 0..100i64 {
        let v = a as f32 / 99.0;
        map.insert(54 + a, (v, v, v));
    }

    map
});

/// Look up a PyMOL color index and return its RGB values.
pub fn pymol_color_rgb(index: i64) -> Option<(f32, f32, f32)> {
    PYMOL_COLOR_MAP.get(&index).copied()
}
