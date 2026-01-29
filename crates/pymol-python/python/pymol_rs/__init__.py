"""
PyMOL-RS: A modern Rust implementation of PyMOL molecular visualization.

This package provides Python bindings to the pymol-rs Rust library,
offering a familiar PyMOL-like API for molecular visualization.

Example usage:
    >>> from pymol_rs import cmd
    >>> cmd.load("protein.pdb")
    >>> cmd.show("cartoon")
    >>> cmd.color("green", "chain A")
    >>> cmd.show_gui()  # Show the visualization window
    >>> cmd.zoom()

To run the GUI directly from the command line:
    $ pymol-rs protein.pdb
"""

from pymol_rs._pymol_rs import (
    # GUI launcher function
    run_gui,
    # Lazy cmd creator (internal)
    _create_cmd,
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
from pymol_rs._pymol_rs import mol
from pymol_rs._pymol_rs import io
from pymol_rs._pymol_rs import selecting
from pymol_rs._pymol_rs import color
from pymol_rs._pymol_rs import settings


# Lazy initialization of cmd
# The cmd object is only created when first accessed, to avoid connection
# errors when the module is imported just for run_gui (e.g., CLI usage)
_cmd_instance = None


def __getattr__(name: str):
    """Lazy attribute access for the cmd object."""
    global _cmd_instance
    if name == "cmd":
        if _cmd_instance is None:
            _cmd_instance = _create_cmd()
        return _cmd_instance
    raise AttributeError(f"module 'pymol_rs' has no attribute {name!r}")


__version__ = "0.1.0"
__all__ = [
    # Main interface
    "cmd",
    # GUI launcher
    "run_gui",
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
