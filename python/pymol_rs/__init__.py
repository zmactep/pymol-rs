"""
PyMOL-RS: Python bindings for the Rust reimplementation of PyMOL.

Unified package that works in two modes:
- Standalone: ``python -c "from pymol_rs import cmd; cmd.load('x.pdb')"``
- Embedded: ``run script.py`` inside the PyMOL-RS GUI (plugin sets backend)

The ``cmd`` object is created lazily on first access.
"""

import sys

from ._cmd import Cmd
from ._pymol_rs import (
    Atom,
    Bond,
    Color,
    CoordSet,
    Element,
    # Types
    ObjectMolecule,
    # Exceptions
    PymolError,
    SelectionError,
    SelectionResult,
    color,
    io,
    # Submodules
    mol,
    selecting,
    settings,
)

__version__ = "0.2.0"

__all__ = [
    "cmd",
    "Cmd",
    "ObjectMolecule",
    "Atom",
    "Bond",
    "Element",
    "CoordSet",
    "Color",
    "SelectionResult",
    "PymolError",
    "SelectionError",
    "mol",
    "io",
    "selecting",
    "color",
    "settings",
]

# Lazy cmd singleton
_cmd = None


class Store:
    def __getattr__(self, name):
        return self.__dict__.get(name, None)

    def __setattr__(self, name, value):
        object.__setattr__(self, name, value)


stored = Store()


def _get_cmd():
    """Get or create the global cmd instance."""
    global _cmd
    if _cmd is None:
        if hasattr(sys, "_pymolrs_backend"):
            # Embedded mode: running inside PyMOL-RS GUI (plugin set this up)
            _cmd = Cmd(sys._pymolrs_backend)
        else:
            # Standalone mode: own Session + CommandExecutor
            from ._pymol_rs import _create_standalone_backend

            _cmd = Cmd(_create_standalone_backend())
    return _cmd


def __getattr__(name):
    if name == "cmd":
        return _get_cmd()
    raise AttributeError(f"module 'pymol_rs' has no attribute {name!r}")
