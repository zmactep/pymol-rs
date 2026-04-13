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

    def fetch(self, code, name=None, state=0, type_="cif", sync=False, quiet=True):
        """Fetch a structure from the PDB."""
        cmd_str = f"fetch {code}"
        if name:
            cmd_str += f", {name}"
        if type_ != "cif":
            cmd_str += f", type={type_}"
        if sync:
            cmd_str += ", async=0"
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

    def get_viewport_image(self):
        """Get viewport image as numpy array (H, W, 4) uint8, or None."""
        return self._backend.get_viewport_image()

    def set_viewport_image(self, array):
        """Set viewport image from numpy array (H, W, 4) uint8."""
        self._backend.set_viewport_image(array)

    def clear_viewport_image(self):
        """Clear viewport image overlay."""
        self._backend.clear_viewport_image()

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
    # Iteration
    # =====================================================================

    def iterate(self, selection, expression, space=None):
        """Execute expression for each atom in selection (read-only).

        Atom properties available as local variables:
            name, resn, resv, resi, chain, segi, alt, elem,
            b, q, vdw, partial_charge, formal_charge,
            ss, color, type, hetatm, index, ID, rank, model,
            x, y, z

        Use ``space`` to accumulate data across atoms::

            mylist = []
            cmd.iterate("name CA", "mylist.append(b)", space=locals())
        """
        return self._backend.iterate(selection, expression, space)

    def alter(self, selection, expression, space=None):
        """Execute expression for each atom in selection (read-write).

        Modifiable properties:
            name, resn, resv, chain, segi, alt, elem,
            b, q, vdw, partial_charge, formal_charge,
            ss, color, type

        Coordinates (x, y, z) are read-only — use alter_state for those.

        Example::

            cmd.alter("all", "b=0.0")
            cmd.alter("chain A", "chain='B'")
        """
        return self._backend.alter(selection, expression, space)

    # =====================================================================
    # Keybindings
    # =====================================================================

    def set_key(self, key, callback):
        """Bind a key or key combination to a Python callback.

        The callback is called with no arguments when the key is pressed.
        Rebinding the same key replaces the previous callback.

        Examples::

            cmd.set_key("F1", my_help_fn)
            cmd.set_key("ctrl+s", lambda: cmd.save("output.pdb"))
            cmd.set_key("ctrl+shift+r", reload_fn)
        """
        if not callable(callback):
            raise TypeError("callback must be callable")
        self._backend.set_key(key, callback)

    def unset_key(self, key):
        """Unbind a key or key combination.

        No error if the key was not bound.

        Examples::

            cmd.unset_key("ctrl+s")
        """
        self._backend.unset_key(key)

    # =====================================================================
    # General execution
    # =====================================================================

    def do(self, command, quiet=False):
        """Execute an arbitrary command string."""
        self._backend.execute(command, quiet)

    def __repr__(self):
        return "<pymol_rs.Cmd>"
