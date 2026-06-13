# Patinae Python Bindings

Python bindings for the Rust/WebGPU molecular visualization workspace.

## Installation

```bash
pip install patinae
```

## Jupyter Widget

The package includes an [anywidget](https://anywidget.dev/)-based widget for interactive visualization in notebooks. Works in JupyterLab, Jupyter Notebook, VS Code, and Google Colab.

```bash
pip install patinae[widget]
```

```python
from patinae.widget import Viewer

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
