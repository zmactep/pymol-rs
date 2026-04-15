# PyMOL-RS Python Bindings

Python bindings for [PyMOL-RS](https://github.com/zmactep/pymol-rs), a modern molecular visualization tool built in Rust with WebGPU rendering.

## Installation

```bash
pip install pymol-rs
```

## Usage

```python
import pymol_rs

# Load a structure
viewer = pymol_rs.Viewer()
viewer.load("structure.pdb")
```

## Jupyter Widget

PyMOL-RS includes a Jupyter widget for interactive visualization in notebooks:

```bash
pip install pymol-rs[widget]
```

```python
from pymol_rs.widget import MoleculeViewer

viewer = MoleculeViewer()
viewer.load("structure.pdb")
viewer
```

## License

BSD-3-Clause
