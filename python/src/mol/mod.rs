//! Molecular data types for Python
//!
//! This module provides Python bindings for molecular data structures:
//! - `ObjectMolecule` - Main molecular container
//! - `Atom` - Atom properties
//! - `Bond` - Bond connectivity
//! - `Element` - Chemical elements
//! - `CoordSet` - Coordinate sets

mod atom;
mod bond;
mod coordset;
mod element;
mod molecule;

pub use atom::{PyAtom, PyAtomIter};
pub use bond::PyBond;
pub use coordset::PyCoordSet;
pub use element::PyElement;
pub use molecule::PyObjectMolecule;

use pyo3::prelude::*;

/// Register the mol submodule
pub fn register_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyObjectMolecule>()?;
    m.add_class::<PyAtom>()?;
    m.add_class::<PyBond>()?;
    m.add_class::<PyElement>()?;
    m.add_class::<PyCoordSet>()?;
    Ok(())
}
