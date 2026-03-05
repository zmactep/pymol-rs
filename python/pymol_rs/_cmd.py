"""
Unified Cmd class for PyMOL-RS.

Works with any backend (StandaloneBackend or PluginBackend).
All commands are thin wrappers that build command strings
and delegate to ``backend.execute()``.
"""


class Cmd:
    """PyMOL-RS command interface.

    Wraps a backend (standalone or embedded) and provides
    a PyMOL-compatible API.
    """

    def __init__(self, backend):
        self._backend = backend

    # =====================================================================
    # File I/O
    # =====================================================================

    def load(self, filename, object=None, state=0, format=None, quiet=True):
        """Load a molecular file."""
        cmd_str = f"load {filename}"
        if object:
            cmd_str += f", {object}"
        if state:
            cmd_str += f", state={state}"
        if format:
            cmd_str += f", format={format}"
        self._backend.execute(cmd_str, quiet)

    def save(self, filename, selection="all", state=-1, format=None, quiet=True):
        """Save molecular data to a file."""
        cmd_str = f"save {filename}, {selection}"
        if state != -1:
            cmd_str += f", state={state}"
        if format:
            cmd_str += f", format={format}"
        self._backend.execute(cmd_str, quiet)

    def fetch(self, code, name=None, state=0, type_="cif", quiet=True):
        """Fetch a structure from the PDB."""
        cmd_str = f"fetch {code}"
        if name:
            cmd_str += f", {name}"
        if type_ != "cif":
            cmd_str += f", type={type_}"
        self._backend.execute(cmd_str, quiet)

    # =====================================================================
    # Display
    # =====================================================================

    def show(self, representation, selection="all"):
        """Show a representation."""
        self._backend.execute(f"show {representation}, {selection}")

    def hide(self, representation, selection="all"):
        """Hide a representation."""
        self._backend.execute(f"hide {representation}, {selection}")

    def show_as(self, representation, selection="all"):
        """Show only the specified representation (hide others)."""
        self._backend.execute(f"as {representation}, {selection}")

    def color(self, color, selection="all"):
        """Color a selection."""
        self._backend.execute(f"color {color}, {selection}")

    def bg_color(self, color):
        """Set background color."""
        self._backend.execute(f"bg_color {color}")

    # =====================================================================
    # Selections
    # =====================================================================

    def select(self, name, selection):
        """Create a named selection."""
        self._backend.execute(f"select {name}, {selection}")

    def deselect(self):
        """Deselect all."""
        self._backend.execute("deselect")

    def count_atoms(self, selection="all"):
        """Count atoms in a selection."""
        return self._backend.count_atoms(selection)

    # =====================================================================
    # Viewing
    # =====================================================================

    def zoom(self, selection="all", buffer=0.0, complete=0):
        """Zoom to a selection."""
        self._backend.execute(f"zoom {selection}, {buffer}, {complete}")

    def center(self, selection="all"):
        """Center on a selection."""
        self._backend.execute(f"center {selection}")

    def orient(self, selection="all"):
        """Orient on a selection."""
        self._backend.execute(f"orient {selection}")

    def reset(self):
        """Reset the view."""
        self._backend.execute("reset")

    # =====================================================================
    # Objects
    # =====================================================================

    def delete(self, name):
        """Delete an object or selection."""
        self._backend.execute(f"delete {name}")

    def get_names(self):
        """Get the names of all loaded objects."""
        return self._backend.get_names()

    def enable(self, name):
        """Enable (show) an object."""
        self._backend.execute(f"enable {name}")

    def disable(self, name):
        """Disable (hide) an object."""
        self._backend.execute(f"disable {name}")

    # =====================================================================
    # Data access
    # =====================================================================

    def get_model(self, name):
        """Get a molecular object by name."""
        return self._backend.get_model(name)

    # =====================================================================
    # Settings
    # =====================================================================

    def set(self, name, value, selection=None, quiet=True):
        """Set a setting."""
        cmd_str = f"set {name}, {value}"
        if selection:
            cmd_str += f", {selection}"
        self._backend.execute(cmd_str, quiet)

    def get(self, name, selection=None, quiet=True):
        """Get a setting value."""
        cmd_str = f"get {name}"
        if selection:
            cmd_str += f", {selection}"
        self._backend.execute(cmd_str, quiet)

    # =====================================================================
    # Image output
    # =====================================================================

    def png(self, filename, width=0, height=0, ray=0, quiet=True):
        """Save a PNG image."""
        cmd_str = f"png {filename}"
        if width:
            cmd_str += f", {width}"
        if height:
            cmd_str += f", {height}"
        if ray:
            cmd_str += f", ray={ray}"
        self._backend.execute(cmd_str, quiet)

    def ray(self, width=0, height=0, quiet=True):
        """Ray trace the scene."""
        cmd_str = "ray"
        if width:
            cmd_str += f" {width}"
        if height:
            cmd_str += f", {height}"
        self._backend.execute(cmd_str, quiet)

    # =====================================================================
    # Control
    # =====================================================================

    def refresh(self):
        """Refresh the scene."""
        self._backend.execute("refresh")

    def rebuild(self, selection="all"):
        """Rebuild representations."""
        self._backend.execute(f"rebuild {selection}")

    def reinitialize(self, what="everything"):
        """Reinitialize the scene."""
        self._backend.execute(f"reinitialize {what}")

    reinit = reinitialize

    # =====================================================================
    # General execution
    # =====================================================================

    def do(self, command, quiet=False):
        """Execute an arbitrary command string."""
        self._backend.execute(command, quiet)

    def __repr__(self):
        return "<pymol_rs.Cmd>"
