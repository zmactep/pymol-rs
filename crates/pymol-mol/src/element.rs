//! Chemical element definitions
//!
//! Provides the `Element` enum representing chemical elements with their
//! atomic numbers, symbols, names, VdW radii, and atomic masses.

use ahash::AHashMap;
use std::fmt;
use std::sync::OnceLock;

/// Chemical element with atomic properties
///
/// Elements are represented by their atomic number, with special handling
/// for the "lonepair" pseudo-element (atomic number 0) used in some
/// molecular modeling contexts.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[repr(u8)]
pub enum Element {
    /// Lone pair / pseudo-atom (atomic number 0)
    #[default]
    Unknown = 0,
    /// Hydrogen (atomic number 1)
    Hydrogen = 1,
    /// Helium (atomic number 2)
    Helium = 2,
    /// Lithium (atomic number 3)
    Lithium = 3,
    /// Beryllium (atomic number 4)
    Beryllium = 4,
    /// Boron (atomic number 5)
    Boron = 5,
    /// Carbon (atomic number 6)
    Carbon = 6,
    /// Nitrogen (atomic number 7)
    Nitrogen = 7,
    /// Oxygen (atomic number 8)
    Oxygen = 8,
    /// Fluorine (atomic number 9)
    Fluorine = 9,
    /// Neon (atomic number 10)
    Neon = 10,
    /// Sodium (atomic number 11)
    Sodium = 11,
    /// Magnesium (atomic number 12)
    Magnesium = 12,
    /// Aluminum (atomic number 13)
    Aluminum = 13,
    /// Silicon (atomic number 14)
    Silicon = 14,
    /// Phosphorus (atomic number 15)
    Phosphorus = 15,
    /// Sulfur (atomic number 16)
    Sulfur = 16,
    /// Chlorine (atomic number 17)
    Chlorine = 17,
    /// Argon (atomic number 18)
    Argon = 18,
    /// Potassium (atomic number 19)
    Potassium = 19,
    /// Calcium (atomic number 20)
    Calcium = 20,
    /// Scandium (atomic number 21)
    Scandium = 21,
    /// Titanium (atomic number 22)
    Titanium = 22,
    /// Vanadium (atomic number 23)
    Vanadium = 23,
    /// Chromium (atomic number 24)
    Chromium = 24,
    /// Manganese (atomic number 25)
    Manganese = 25,
    /// Iron (atomic number 26)
    Iron = 26,
    /// Cobalt (atomic number 27)
    Cobalt = 27,
    /// Nickel (atomic number 28)
    Nickel = 28,
    /// Copper (atomic number 29)
    Copper = 29,
    /// Zinc (atomic number 30)
    Zinc = 30,
    /// Gallium (atomic number 31)
    Gallium = 31,
    /// Germanium (atomic number 32)
    Germanium = 32,
    /// Arsenic (atomic number 33)
    Arsenic = 33,
    /// Selenium (atomic number 34)
    Selenium = 34,
    /// Bromine (atomic number 35)
    Bromine = 35,
    /// Krypton (atomic number 36)
    Krypton = 36,
    /// Rubidium (atomic number 37)
    Rubidium = 37,
    /// Strontium (atomic number 38)
    Strontium = 38,
    /// Yttrium (atomic number 39)
    Yttrium = 39,
    /// Zirconium (atomic number 40)
    Zirconium = 40,
    /// Niobium (atomic number 41)
    Niobium = 41,
    /// Molybdenum (atomic number 42)
    Molybdenum = 42,
    /// Technetium (atomic number 43)
    Technetium = 43,
    /// Ruthenium (atomic number 44)
    Ruthenium = 44,
    /// Rhodium (atomic number 45)
    Rhodium = 45,
    /// Palladium (atomic number 46)
    Palladium = 46,
    /// Silver (atomic number 47)
    Silver = 47,
    /// Cadmium (atomic number 48)
    Cadmium = 48,
    /// Indium (atomic number 49)
    Indium = 49,
    /// Tin (atomic number 50)
    Tin = 50,
    /// Antimony (atomic number 51)
    Antimony = 51,
    /// Tellurium (atomic number 52)
    Tellurium = 52,
    /// Iodine (atomic number 53)
    Iodine = 53,
    /// Xenon (atomic number 54)
    Xenon = 54,
    /// Cesium (atomic number 55)
    Cesium = 55,
    /// Barium (atomic number 56)
    Barium = 56,
    /// Lanthanum (atomic number 57)
    Lanthanum = 57,
    /// Cerium (atomic number 58)
    Cerium = 58,
    /// Praseodymium (atomic number 59)
    Praseodymium = 59,
    /// Neodymium (atomic number 60)
    Neodymium = 60,
    /// Promethium (atomic number 61)
    Promethium = 61,
    /// Samarium (atomic number 62)
    Samarium = 62,
    /// Europium (atomic number 63)
    Europium = 63,
    /// Gadolinium (atomic number 64)
    Gadolinium = 64,
    /// Terbium (atomic number 65)
    Terbium = 65,
    /// Dysprosium (atomic number 66)
    Dysprosium = 66,
    /// Holmium (atomic number 67)
    Holmium = 67,
    /// Erbium (atomic number 68)
    Erbium = 68,
    /// Thulium (atomic number 69)
    Thulium = 69,
    /// Ytterbium (atomic number 70)
    Ytterbium = 70,
    /// Lutetium (atomic number 71)
    Lutetium = 71,
    /// Hafnium (atomic number 72)
    Hafnium = 72,
    /// Tantalum (atomic number 73)
    Tantalum = 73,
    /// Tungsten (atomic number 74)
    Tungsten = 74,
    /// Rhenium (atomic number 75)
    Rhenium = 75,
    /// Osmium (atomic number 76)
    Osmium = 76,
    /// Iridium (atomic number 77)
    Iridium = 77,
    /// Platinum (atomic number 78)
    Platinum = 78,
    /// Gold (atomic number 79)
    Gold = 79,
    /// Mercury (atomic number 80)
    Mercury = 80,
    /// Thallium (atomic number 81)
    Thallium = 81,
    /// Lead (atomic number 82)
    Lead = 82,
    /// Bismuth (atomic number 83)
    Bismuth = 83,
    /// Polonium (atomic number 84)
    Polonium = 84,
    /// Astatine (atomic number 85)
    Astatine = 85,
    /// Radon (atomic number 86)
    Radon = 86,
    /// Francium (atomic number 87)
    Francium = 87,
    /// Radium (atomic number 88)
    Radium = 88,
    /// Actinium (atomic number 89)
    Actinium = 89,
    /// Thorium (atomic number 90)
    Thorium = 90,
    /// Protactinium (atomic number 91)
    Protactinium = 91,
    /// Uranium (atomic number 92)
    Uranium = 92,
    /// Neptunium (atomic number 93)
    Neptunium = 93,
    /// Plutonium (atomic number 94)
    Plutonium = 94,
    /// Americium (atomic number 95)
    Americium = 95,
    /// Curium (atomic number 96)
    Curium = 96,
    /// Berkelium (atomic number 97)
    Berkelium = 97,
    /// Californium (atomic number 98)
    Californium = 98,
    /// Einsteinium (atomic number 99)
    Einsteinium = 99,
    /// Fermium (atomic number 100)
    Fermium = 100,
    /// Mendelevium (atomic number 101)
    Mendelevium = 101,
    /// Nobelium (atomic number 102)
    Nobelium = 102,
    /// Lawrencium (atomic number 103)
    Lawrencium = 103,
    /// Rutherfordium (atomic number 104)
    Rutherfordium = 104,
    /// Dubnium (atomic number 105)
    Dubnium = 105,
    /// Seaborgium (atomic number 106)
    Seaborgium = 106,
    /// Bohrium (atomic number 107)
    Bohrium = 107,
    /// Hassium (atomic number 108)
    Hassium = 108,
    /// Meitnerium (atomic number 109)
    Meitnerium = 109,
    /// Darmstadtium (atomic number 110)
    Darmstadtium = 110,
    /// Roentgenium (atomic number 111)
    Roentgenium = 111,
    /// Copernicium (atomic number 112)
    Copernicium = 112,
    /// Nihonium (atomic number 113)
    Nihonium = 113,
    /// Flerovium (atomic number 114)
    Flerovium = 114,
    /// Moscovium (atomic number 115)
    Moscovium = 115,
    /// Livermorium (atomic number 116)
    Livermorium = 116,
    /// Tennessine (atomic number 117)
    Tennessine = 117,
    /// Oganesson (atomic number 118)
    Oganesson = 118,
}

/// Element data table entry
struct ElementData {
    name: &'static str,
    symbol: &'static str,
    vdw: f32,
    mass: f32,
}

/// Element data table indexed by atomic number
/// Data from PyMOL's ElementTable in AtomInfo.cpp
static ELEMENT_DATA: &[ElementData] = &[
    ElementData { name: "lonepair", symbol: "LP", vdw: 0.50, mass: 0.0 },           // 0
    ElementData { name: "hydrogen", symbol: "H", vdw: 1.20, mass: 1.00794 },        // 1
    ElementData { name: "helium", symbol: "He", vdw: 1.40, mass: 4.002602 },        // 2
    ElementData { name: "lithium", symbol: "Li", vdw: 1.82, mass: 6.941 },          // 3
    ElementData { name: "beryllium", symbol: "Be", vdw: 1.80, mass: 9.012182 },     // 4
    ElementData { name: "boron", symbol: "B", vdw: 1.85, mass: 10.811 },            // 5
    ElementData { name: "carbon", symbol: "C", vdw: 1.70, mass: 12.0107 },          // 6
    ElementData { name: "nitrogen", symbol: "N", vdw: 1.55, mass: 14.0067 },        // 7
    ElementData { name: "oxygen", symbol: "O", vdw: 1.52, mass: 15.9994 },          // 8
    ElementData { name: "fluorine", symbol: "F", vdw: 1.47, mass: 18.998403 },      // 9
    ElementData { name: "neon", symbol: "Ne", vdw: 1.54, mass: 20.1797 },           // 10
    ElementData { name: "sodium", symbol: "Na", vdw: 2.27, mass: 22.98977 },        // 11
    ElementData { name: "magnesium", symbol: "Mg", vdw: 1.73, mass: 24.305 },       // 12
    ElementData { name: "aluminum", symbol: "Al", vdw: 2.00, mass: 26.981538 },     // 13
    ElementData { name: "silicon", symbol: "Si", vdw: 2.10, mass: 28.0855 },        // 14
    ElementData { name: "phosphorus", symbol: "P", vdw: 1.80, mass: 30.973761 },    // 15
    ElementData { name: "sulfur", symbol: "S", vdw: 1.80, mass: 32.065 },           // 16
    ElementData { name: "chlorine", symbol: "Cl", vdw: 1.75, mass: 35.453 },        // 17
    ElementData { name: "argon", symbol: "Ar", vdw: 1.88, mass: 39.948 },           // 18
    ElementData { name: "potassium", symbol: "K", vdw: 2.75, mass: 39.0983 },       // 19
    ElementData { name: "calcium", symbol: "Ca", vdw: 1.80, mass: 40.078 },         // 20
    ElementData { name: "scandium", symbol: "Sc", vdw: 1.80, mass: 44.95591 },      // 21
    ElementData { name: "titanium", symbol: "Ti", vdw: 1.80, mass: 47.867 },        // 22
    ElementData { name: "vanadium", symbol: "V", vdw: 1.80, mass: 50.9415 },        // 23
    ElementData { name: "chromium", symbol: "Cr", vdw: 1.80, mass: 51.9961 },       // 24
    ElementData { name: "manganese", symbol: "Mn", vdw: 1.73, mass: 54.938049 },    // 25
    ElementData { name: "iron", symbol: "Fe", vdw: 1.80, mass: 55.845 },            // 26
    ElementData { name: "cobalt", symbol: "Co", vdw: 1.80, mass: 58.9332 },         // 27
    ElementData { name: "nickel", symbol: "Ni", vdw: 1.63, mass: 58.6934 },         // 28
    ElementData { name: "copper", symbol: "Cu", vdw: 1.40, mass: 63.546 },          // 29
    ElementData { name: "zinc", symbol: "Zn", vdw: 1.39, mass: 65.39 },             // 30
    ElementData { name: "gallium", symbol: "Ga", vdw: 1.87, mass: 69.723 },         // 31
    ElementData { name: "germanium", symbol: "Ge", vdw: 1.80, mass: 72.64 },        // 32
    ElementData { name: "arsenic", symbol: "As", vdw: 1.85, mass: 74.9216 },        // 33
    ElementData { name: "selenium", symbol: "Se", vdw: 1.90, mass: 78.96 },         // 34
    ElementData { name: "bromine", symbol: "Br", vdw: 1.85, mass: 79.904 },         // 35
    ElementData { name: "krypton", symbol: "Kr", vdw: 2.02, mass: 83.8 },           // 36
    ElementData { name: "rubidium", symbol: "Rb", vdw: 1.80, mass: 85.4678 },       // 37
    ElementData { name: "strontium", symbol: "Sr", vdw: 1.80, mass: 87.62 },        // 38
    ElementData { name: "yttrium", symbol: "Y", vdw: 1.80, mass: 88.90585 },        // 39
    ElementData { name: "zirconium", symbol: "Zr", vdw: 1.80, mass: 91.224 },       // 40
    ElementData { name: "niobium", symbol: "Nb", vdw: 1.80, mass: 92.90638 },       // 41
    ElementData { name: "molybdenum", symbol: "Mo", vdw: 1.80, mass: 95.94 },       // 42
    ElementData { name: "technetium", symbol: "Tc", vdw: 1.80, mass: 98.0 },        // 43
    ElementData { name: "ruthenium", symbol: "Ru", vdw: 1.80, mass: 101.07 },       // 44
    ElementData { name: "rhodium", symbol: "Rh", vdw: 1.80, mass: 102.9055 },       // 45
    ElementData { name: "palladium", symbol: "Pd", vdw: 1.63, mass: 106.42 },       // 46
    ElementData { name: "silver", symbol: "Ag", vdw: 1.72, mass: 107.8682 },        // 47
    ElementData { name: "cadmium", symbol: "Cd", vdw: 1.58, mass: 112.411 },        // 48
    ElementData { name: "indium", symbol: "In", vdw: 1.93, mass: 114.818 },         // 49
    ElementData { name: "tin", symbol: "Sn", vdw: 2.17, mass: 118.71 },             // 50
    ElementData { name: "antimony", symbol: "Sb", vdw: 1.80, mass: 121.76 },        // 51
    ElementData { name: "tellurium", symbol: "Te", vdw: 2.06, mass: 127.6 },        // 52
    ElementData { name: "iodine", symbol: "I", vdw: 1.98, mass: 126.90447 },        // 53
    ElementData { name: "xenon", symbol: "Xe", vdw: 2.16, mass: 131.293 },          // 54
    ElementData { name: "cesium", symbol: "Cs", vdw: 1.80, mass: 132.90545 },       // 55
    ElementData { name: "barium", symbol: "Ba", vdw: 1.80, mass: 137.327 },         // 56
    ElementData { name: "lanthanum", symbol: "La", vdw: 1.80, mass: 138.9055 },     // 57
    ElementData { name: "cerium", symbol: "Ce", vdw: 1.80, mass: 140.116 },         // 58
    ElementData { name: "praseodymium", symbol: "Pr", vdw: 1.80, mass: 140.90765 }, // 59
    ElementData { name: "neodymium", symbol: "Nd", vdw: 1.80, mass: 144.24 },       // 60
    ElementData { name: "promethium", symbol: "Pm", vdw: 1.80, mass: 145.0 },       // 61
    ElementData { name: "samarium", symbol: "Sm", vdw: 1.80, mass: 150.36 },        // 62
    ElementData { name: "europium", symbol: "Eu", vdw: 1.80, mass: 151.964 },       // 63
    ElementData { name: "gadolinium", symbol: "Gd", vdw: 1.80, mass: 157.25 },      // 64
    ElementData { name: "terbium", symbol: "Tb", vdw: 1.80, mass: 158.92534 },      // 65
    ElementData { name: "dysprosium", symbol: "Dy", vdw: 1.80, mass: 162.5 },       // 66
    ElementData { name: "holmium", symbol: "Ho", vdw: 1.80, mass: 164.93032 },      // 67
    ElementData { name: "erbium", symbol: "Er", vdw: 1.80, mass: 167.259 },         // 68
    ElementData { name: "thulium", symbol: "Tm", vdw: 1.80, mass: 168.93421 },      // 69
    ElementData { name: "ytterbium", symbol: "Yb", vdw: 1.80, mass: 173.04 },       // 70
    ElementData { name: "lutetium", symbol: "Lu", vdw: 1.80, mass: 174.967 },       // 71
    ElementData { name: "hafnium", symbol: "Hf", vdw: 1.80, mass: 178.49 },         // 72
    ElementData { name: "tantalum", symbol: "Ta", vdw: 1.80, mass: 180.9479 },      // 73
    ElementData { name: "tungsten", symbol: "W", vdw: 1.80, mass: 183.84 },         // 74
    ElementData { name: "rhenium", symbol: "Re", vdw: 1.80, mass: 186.207 },        // 75
    ElementData { name: "osmium", symbol: "Os", vdw: 1.80, mass: 190.23 },          // 76
    ElementData { name: "iridium", symbol: "Ir", vdw: 1.80, mass: 192.217 },        // 77
    ElementData { name: "platinum", symbol: "Pt", vdw: 1.75, mass: 195.078 },       // 78
    ElementData { name: "gold", symbol: "Au", vdw: 1.66, mass: 196.96655 },         // 79
    ElementData { name: "mercury", symbol: "Hg", vdw: 1.55, mass: 200.59 },         // 80
    ElementData { name: "thallium", symbol: "Tl", vdw: 1.96, mass: 204.3833 },      // 81
    ElementData { name: "lead", symbol: "Pb", vdw: 2.02, mass: 207.2 },             // 82
    ElementData { name: "bismuth", symbol: "Bi", vdw: 1.80, mass: 208.98038 },      // 83
    ElementData { name: "polonium", symbol: "Po", vdw: 1.80, mass: 208.98 },        // 84
    ElementData { name: "astatine", symbol: "At", vdw: 1.80, mass: 209.99 },        // 85
    ElementData { name: "radon", symbol: "Rn", vdw: 1.80, mass: 222.02 },           // 86
    ElementData { name: "francium", symbol: "Fr", vdw: 1.80, mass: 223.02 },        // 87
    ElementData { name: "radium", symbol: "Ra", vdw: 1.80, mass: 226.03 },          // 88
    ElementData { name: "actinium", symbol: "Ac", vdw: 1.80, mass: 227.03 },        // 89
    ElementData { name: "thorium", symbol: "Th", vdw: 1.80, mass: 232.0381 },       // 90
    ElementData { name: "protactinium", symbol: "Pa", vdw: 1.80, mass: 231.03588 }, // 91
    ElementData { name: "uranium", symbol: "U", vdw: 1.86, mass: 238.02891 },       // 92
    ElementData { name: "neptunium", symbol: "Np", vdw: 1.80, mass: 237.05 },       // 93
    ElementData { name: "plutonium", symbol: "Pu", vdw: 1.80, mass: 244.06 },       // 94
    ElementData { name: "americium", symbol: "Am", vdw: 1.80, mass: 243.06 },       // 95
    ElementData { name: "curium", symbol: "Cm", vdw: 1.80, mass: 247.07 },          // 96
    ElementData { name: "berkelium", symbol: "Bk", vdw: 1.80, mass: 247.07 },       // 97
    ElementData { name: "californium", symbol: "Cf", vdw: 1.80, mass: 251.08 },     // 98
    ElementData { name: "einsteinium", symbol: "Es", vdw: 1.80, mass: 252.08 },     // 99
    ElementData { name: "fermium", symbol: "Fm", vdw: 1.80, mass: 257.1 },          // 100
    ElementData { name: "mendelevium", symbol: "Md", vdw: 1.80, mass: 258.1 },      // 101
    ElementData { name: "nobelium", symbol: "No", vdw: 1.80, mass: 259.1 },         // 102
    ElementData { name: "lawrencium", symbol: "Lr", vdw: 1.80, mass: 262.11 },      // 103
    ElementData { name: "rutherfordium", symbol: "Rf", vdw: 1.80, mass: 261.11 },   // 104
    ElementData { name: "dubnium", symbol: "Db", vdw: 1.80, mass: 262.11 },         // 105
    ElementData { name: "seaborgium", symbol: "Sg", vdw: 1.80, mass: 266.12 },      // 106
    ElementData { name: "bohrium", symbol: "Bh", vdw: 1.80, mass: 264.12 },         // 107
    ElementData { name: "hassium", symbol: "Hs", vdw: 1.80, mass: 269.13 },         // 108
    ElementData { name: "meitnerium", symbol: "Mt", vdw: 1.80, mass: 268.14 },      // 109
    ElementData { name: "darmstadtium", symbol: "Ds", vdw: 1.80, mass: 281.0 },     // 110
    ElementData { name: "roentgenium", symbol: "Rg", vdw: 1.80, mass: 281.0 },      // 111
    ElementData { name: "copernicium", symbol: "Cn", vdw: 1.80, mass: 285.0 },      // 112
    ElementData { name: "nihonium", symbol: "Nh", vdw: 1.80, mass: 286.0 },         // 113
    ElementData { name: "flerovium", symbol: "Fl", vdw: 1.80, mass: 289.0 },        // 114
    ElementData { name: "moscovium", symbol: "Mc", vdw: 1.80, mass: 290.0 },        // 115
    ElementData { name: "livermorium", symbol: "Lv", vdw: 1.80, mass: 293.0 },      // 116
    ElementData { name: "tennessine", symbol: "Ts", vdw: 1.80, mass: 294.0 },       // 117
    ElementData { name: "oganesson", symbol: "Og", vdw: 1.80, mass: 294.0 },        // 118
];

/// Default VdW radius for unknown elements
pub const DEFAULT_VDW_RADIUS: f32 = 1.80;

/// Number of elements in the table (including Unknown/LP at index 0)
pub const ELEMENT_COUNT: usize = 119;

/// Static symbol lookup map for O(1) symbol-to-element conversion.
/// Initialized lazily on first access.
static SYMBOL_MAP: OnceLock<AHashMap<&'static str, Element>> = OnceLock::new();

/// Initialize the symbol lookup map
fn init_symbol_map() -> AHashMap<&'static str, Element> {
    let mut map = AHashMap::with_capacity(ELEMENT_COUNT + 10); // Extra for aliases

    // Add all elements from the data table
    for (i, data) in ELEMENT_DATA.iter().enumerate() {
        if let Some(elem) = Element::from_atomic_number(i as u8) {
            map.insert(data.symbol, elem);
        }
    }

    // Add special aliases
    map.insert("D", Element::Hydrogen); // Deuterium
    map.insert("Q", Element::Hydrogen); // Unspecified hydrogen
    map.insert("Lp", Element::Unknown); // Lone pair
    map.insert("LP", Element::Unknown); // Lone pair (uppercase)

    map
}

/// Get the symbol lookup map (initializes on first access)
fn symbol_map() -> &'static AHashMap<&'static str, Element> {
    SYMBOL_MAP.get_or_init(init_symbol_map)
}

impl Element {
    /// Create an element from its atomic number
    ///
    /// Returns `None` if the atomic number is out of range (> 118)
    #[inline]
    pub fn from_atomic_number(n: u8) -> Option<Self> {
        if (n as usize) < ELEMENT_COUNT {
            // SAFETY: All values 0-118 are valid Element discriminants
            Some(unsafe { std::mem::transmute(n) })
        } else {
            None
        }
    }

    /// Create an element from its symbol (case-insensitive)
    ///
    /// Supports common variations like "D" for deuterium (returns Hydrogen)
    pub fn from_symbol(symbol: &str) -> Option<Self> {
        let symbol = symbol.trim();
        if symbol.is_empty() {
            return None;
        }

        // Fast path for common single-letter elements
        if symbol.len() == 1 {
            return match symbol.chars().next()?.to_ascii_uppercase() {
                'H' | 'D' | 'Q' => Some(Element::Hydrogen), // D = deuterium, Q = unspecified
                'C' => Some(Element::Carbon),
                'N' => Some(Element::Nitrogen),
                'O' => Some(Element::Oxygen),
                'F' => Some(Element::Fluorine),
                'P' => Some(Element::Phosphorus),
                'S' => Some(Element::Sulfur),
                'K' => Some(Element::Potassium),
                'V' => Some(Element::Vanadium),
                'Y' => Some(Element::Yttrium),
                'I' => Some(Element::Iodine),
                'W' => Some(Element::Tungsten),
                'U' => Some(Element::Uranium),
                'B' => Some(Element::Boron),
                _ => None,
            };
        }

        // Titlecase the symbol for lookup (e.g., "HE" -> "He", "ca" -> "Ca")
        let mut titlecase = String::with_capacity(symbol.len());
        for (i, c) in symbol.chars().enumerate() {
            if i == 0 {
                titlecase.push(c.to_ascii_uppercase());
            } else {
                titlecase.push(c.to_ascii_lowercase());
            }
        }

        // O(1) lookup using the hash map
        symbol_map().get(titlecase.as_str()).copied()
    }

    /// Get the atomic number
    #[inline]
    pub const fn atomic_number(&self) -> u8 {
        *self as u8
    }

    /// Get the element symbol (e.g., "H", "He", "Li")
    #[inline]
    pub fn symbol(&self) -> &'static str {
        ELEMENT_DATA[*self as usize].symbol
    }

    /// Get the full element name (e.g., "hydrogen", "helium", "lithium")
    #[inline]
    pub fn name(&self) -> &'static str {
        ELEMENT_DATA[*self as usize].name
    }

    /// Get the default Van der Waals radius in Angstroms
    #[inline]
    pub fn vdw_radius(&self) -> f32 {
        ELEMENT_DATA[*self as usize].vdw
    }

    /// Get the atomic mass in Daltons (g/mol)
    #[inline]
    pub fn mass(&self) -> f32 {
        ELEMENT_DATA[*self as usize].mass
    }

    /// Check if this is a common organic element (H, C, N, O, P, S)
    #[inline]
    pub fn is_organic(&self) -> bool {
        matches!(
            self,
            Element::Hydrogen
                | Element::Carbon
                | Element::Nitrogen
                | Element::Oxygen
                | Element::Phosphorus
                | Element::Sulfur
        )
    }

    /// Check if this is a halogen (F, Cl, Br, I, At)
    #[inline]
    pub fn is_halogen(&self) -> bool {
        matches!(
            self,
            Element::Fluorine
                | Element::Chlorine
                | Element::Bromine
                | Element::Iodine
                | Element::Astatine
        )
    }

    /// Check if this is a metal
    #[inline]
    pub fn is_metal(&self) -> bool {
        let n = *self as u8;
        // Alkali metals
        matches!(n, 3 | 11 | 19 | 37 | 55 | 87)
            // Alkaline earth metals
            || matches!(n, 4 | 12 | 20 | 38 | 56 | 88)
            // Transition metals
            || (21..=30).contains(&n)
            || (39..=48).contains(&n)
            || (72..=80).contains(&n)
            || (104..=112).contains(&n)
            // Lanthanides
            || (57..=71).contains(&n)
            // Actinides
            || (89..=103).contains(&n)
            // Post-transition metals
            || matches!(n, 13 | 31 | 49 | 50 | 81 | 82 | 83 | 113 | 114 | 115 | 116)
    }

    /// Check if this is a noble gas
    #[inline]
    pub fn is_noble_gas(&self) -> bool {
        matches!(
            self,
            Element::Helium
                | Element::Neon
                | Element::Argon
                | Element::Krypton
                | Element::Xenon
                | Element::Radon
                | Element::Oganesson
        )
    }

    /// Check if this is hydrogen
    #[inline]
    pub fn is_hydrogen(&self) -> bool {
        *self == Element::Hydrogen
    }

    /// Check if this is carbon
    #[inline]
    pub fn is_carbon(&self) -> bool {
        *self == Element::Carbon
    }

    /// Check if this element is unknown/lonepair
    #[inline]
    pub fn is_unknown(&self) -> bool {
        *self == Element::Unknown
    }
}

impl fmt::Display for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.symbol())
    }
}

impl From<u8> for Element {
    fn from(n: u8) -> Self {
        Element::from_atomic_number(n).unwrap_or(Element::Unknown)
    }
}

impl From<Element> for u8 {
    fn from(e: Element) -> Self {
        e.atomic_number()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_element_from_atomic_number() {
        assert_eq!(Element::from_atomic_number(1), Some(Element::Hydrogen));
        assert_eq!(Element::from_atomic_number(6), Some(Element::Carbon));
        assert_eq!(Element::from_atomic_number(118), Some(Element::Oganesson));
        assert_eq!(Element::from_atomic_number(119), None);
    }

    #[test]
    fn test_element_from_symbol() {
        assert_eq!(Element::from_symbol("H"), Some(Element::Hydrogen));
        assert_eq!(Element::from_symbol("h"), Some(Element::Hydrogen));
        assert_eq!(Element::from_symbol("He"), Some(Element::Helium));
        assert_eq!(Element::from_symbol("HE"), Some(Element::Helium));
        assert_eq!(Element::from_symbol("C"), Some(Element::Carbon));
        assert_eq!(Element::from_symbol("Ca"), Some(Element::Calcium));
        assert_eq!(Element::from_symbol("CA"), Some(Element::Calcium));
        assert_eq!(Element::from_symbol("D"), Some(Element::Hydrogen)); // Deuterium
        assert_eq!(Element::from_symbol("Xx"), None);
    }

    #[test]
    fn test_element_properties() {
        let carbon = Element::Carbon;
        assert_eq!(carbon.atomic_number(), 6);
        assert_eq!(carbon.symbol(), "C");
        assert_eq!(carbon.name(), "carbon");
        assert!((carbon.vdw_radius() - 1.70).abs() < 0.01);
        assert!((carbon.mass() - 12.0107).abs() < 0.01);
    }

    #[test]
    fn test_element_classification() {
        assert!(Element::Carbon.is_organic());
        assert!(Element::Chlorine.is_halogen());
        assert!(Element::Iron.is_metal());
        assert!(Element::Argon.is_noble_gas());
        assert!(Element::Hydrogen.is_hydrogen());
    }

    #[test]
    fn test_element_display() {
        assert_eq!(format!("{}", Element::Carbon), "C");
        assert_eq!(format!("{}", Element::Helium), "He");
    }

    #[test]
    fn test_element_default() {
        assert_eq!(Element::default(), Element::Unknown);
    }
}
