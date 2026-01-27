"""
PyMOL-RS: A modern Rust implementation of PyMOL molecular visualization.

This package provides Python bindings to the pymol-rs Rust library,
offering a familiar PyMOL-like API for molecular visualization.

Example usage:
    >>> from pymol import cmd
    >>> cmd.load("protein.pdb")
    >>> cmd.show("cartoon")
    >>> cmd.color("green", "chain A")
    >>> cmd.zoom()
"""

from pymol._pymol import (
    # Main command interface
    cmd,
    # Molecular data types
    ObjectMolecule,
    Atom,
    Bond,
    Element,
    CoordSet,
    # Color types
    Color,
    # Selection
    SelectionResult,
    # Exceptions
    PymolError,
    SelectionError,
)

# Re-export submodules
from pymol._pymol import mol
from pymol._pymol import io
from pymol._pymol import selecting
from pymol._pymol import color
from pymol._pymol import settings

__version__ = "0.1.0"
__all__ = [
    # Main interface
    "cmd",
    # Types
    "ObjectMolecule",
    "Atom",
    "Bond",
    "Element",
    "CoordSet",
    "Color",
    "SelectionResult",
    # Exceptions
    "PymolError",
    "SelectionError",
    # Submodules
    "mol",
    "io",
    "selecting",
    "color",
    "settings",
]
