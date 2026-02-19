## Backlog

### Bugs
- [x] **Multistate Transformation Bug**: The `translate` and `rotate` commands currently only affect state 0 (the first state) of an object. They need to respect the `state` argument or the current global state to allow transforming specific conformations in a trajectory/ensemble.
- [ ] **Reinitialize Does Not Reset Settings**: `reinitialize` currently preserves some settings (e.g., background color) instead of resetting everything to defaults. It should fully reset all settings, representations, objects, and state to a clean initial configuration.

### Features
- [ ] **GROMACS .gro Support**: Implement a parser for `.gro` coordinate files (GROMACS format). This should handle atom names, residues, and box vectors if present.
- [ ] **MMTF Format Support (Read + Write)**: Implement full support for MMTF (Macromolecular Transmission Format) in `pymol-io`. Both reading and writing, since MMTF is a compact binary format useful for fast loading and efficient storage/transfer of large structures.
- [ ] **Mouse/Touchpad UX Improvements**: Align mouse behavior with PyMOL standard:
    - **Scroll + Click (Drag)**: Zoom (Z-translation of camera).
    - **Scroll (Wheel/Touchpad) only**: Adjust Slab/Clipping planes (move the front/back clipping planes).
- [ ] **Selection UX Logic**: Revisit selection behavior, especially in Sequence Viewer.
    - **Additive by Default**: Clicking residues should *add* to the current selection (like Ctrl/Cmd+Click currently does), rather than replace it.
    - **Review needed**: This deviates from PyMOL standard (click = replace, shift+click = range, ctrl+click = add). We need to decide if this should be a global setting (`mouse_selection_mode`) or specific to the sequence viewer.
- [ ] **Interactive 3D Picking**: Implement full mouse interaction for 3D objects.
    - **Picking Modes**: Support `picking_mode` setting (atom, CA, residue, chain, molecule, object).
    - **Hover Highlight**: Temporarily highlight the entity under the mouse cursor before clicking.
    - **Click Behavior**: Add the hovered entity to the current selection `sele` (respecting `mouse_selection_mode`).
