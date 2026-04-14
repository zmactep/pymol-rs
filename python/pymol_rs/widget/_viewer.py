"""Viewer widget — anywidget subclass that embeds the WASM molecular viewer."""

import base64
import pathlib

import anywidget
import traitlets

_STATIC = pathlib.Path(__file__).parent / "static"


class Viewer(anywidget.AnyWidget):
    """Interactive 3D molecular viewer for Jupyter notebooks.

    Uses the PyMOL-RS WASM + WebGPU viewer running in the browser.
    Requires a WebGPU-capable browser (Chrome, Edge, or Firefox Nightly).

    Usage::

        from pymol_rs.widget import Viewer

        view = Viewer()
        view.show()
        cmd = view.get_cmd()
        cmd.fetch("1CRN")
        cmd.show("cartoon")
        cmd.color("green", "chain A")
    """

    _esm = pathlib.Path(__file__).parent / "_frontend.js"

    # --- Synced traitlets ---

    # WASM glue JS source (~52KB string)
    _glue_js = traitlets.Unicode("").tag(sync=True)

    # WASM binary as base64 (~4MB string, sent once on init)
    _wasm_b64 = traitlets.Unicode("").tag(sync=True)

    # Fire-and-forget command channel
    _command = traitlets.Unicode("").tag(sync=True)
    _command_id = traitlets.Int(0).tag(sync=True)

    # Synchronous query request/response
    _query_request = traitlets.Dict({}).tag(sync=True)
    _query_response = traitlets.Dict({}).tag(sync=True)
    _query_id = traitlets.Int(0).tag(sync=True)

    # Layout
    _width = traitlets.Unicode("100%").tag(sync=True)
    _height = traitlets.Unicode("500px").tag(sync=True)

    # Picking (click-to-select atoms)
    _picking = traitlets.Bool(False).tag(sync=True)

    def __init__(self, width="100%", height="500px", picking=False, **kwargs):
        glue_js = (_STATIC / "pymol_web_glue.js").read_text()
        wasm_b64 = base64.b64encode(
            (_STATIC / "pymol_web_bg.wasm").read_bytes()
        ).decode("ascii")
        super().__init__(_glue_js=glue_js, _wasm_b64=wasm_b64, **kwargs)
        self._width = width
        self._height = height
        self._picking = picking

    def show(self):
        """Display the widget in the notebook."""
        from IPython.display import display

        display(self)

    def get_cmd(self):
        """Return a Cmd object that proxies commands to the browser WASM viewer.

        The returned object has the same API as ``pymol_rs.cmd``.
        """
        from .._cmd import Cmd
        from ._backend import WidgetBackend

        return Cmd(WidgetBackend(self))
