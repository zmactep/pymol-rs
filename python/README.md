# PyMOL-RS Python Bindings

Python bindings for [PyMOL-RS](https://github.com/zmactep/pymol-rs), a modern molecular visualization tool built in Rust with WebGPU rendering.

## Installation

```bash
pip install pymol-rs
```

## Jupyter Widget

PyMOL-RS includes an [anywidget](https://anywidget.dev/)-based widget for interactive visualization in notebooks. Works in JupyterLab, Jupyter Notebook, VS Code, and Google Colab.

```bash
pip install pymol-rs[widget]
```

```python
from pymol_rs.widget import Viewer

view = Viewer()
view.show()
cmd = view.get_cmd()
cmd.fetch("1CRN")
cmd.show("cartoon")
cmd.color("green", "chain A")
```

Features:
- **Fire-and-forget commands** — `cmd.fetch()`, `cmd.show()`, `cmd.color()`, etc.
- **Synchronous queries** — request/response channel for commands that return data
- **Local file loading** — load structures from the local filesystem into the browser viewer
- **Picking support** — optional click-to-select atoms (`Viewer(picking=True)`)
- **Configurable layout** — `width` and `height` parameters with sensible defaults

## License

BSD-3-Clause
