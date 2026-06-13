//! `DirtyFlags` — bitmask of which aspects of a molecule have changed since
//! the last render sync.
//!
//! Produced by mutation logic (commands, animations, picking) and consumed by
//! the renderer to decide which work to redo. Lives in `patinae-mol` because
//! both the data layer (mutating commands in `patinae-cmd`) and the renderer
//! (`patinae-render`) need the same vocabulary, and `patinae-mol` is already a
//! workspace dependency of both.
//!
//! Flags are short-lived: cleared after each `RenderState::sync` walk. They
//! are intentionally not part of any snapshot — restoring a session always
//! re-renders from scratch with `DirtyFlags::ALL`.

use bitflags::bitflags;

bitflags! {
    /// What changed for an object since the last render sync.
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub struct DirtyFlags: u32 {
        /// Atomic positions changed (animation, dynamics, state switch).
        const COORDS       = 1 << 0;
        /// Per-atom base colour changed (`color red, polymer`).
        const COLOR        = 1 << 1;
        /// Visible-rep mask or other rep config changed (`show / hide`).
        const REPS         = 1 << 2;
        /// Selection bitset changed (`select foo, …`).
        const SELECTION    = 1 << 3;
        /// Per-rep alpha or per-atom transparency override changed.
        const TRANSPARENCY = 1 << 4;
        /// Per-atom visibility changed: both coarse `mask_lut` bits and
        /// per-representation visible bits in `AtomGpu::repr_flags`.
        const VISIBILITY   = 1 << 5;
        /// Hover atom changed (single-atom highlight). Set by picking on
        /// mouse-move; consumed by the marking pass.
        const HOVER        = 1 << 6;
        /// Bonds added or removed; rep instance topology must rebuild.
        const TOPOLOGY     = 1 << 7;
        /// Scene-wide LOD bucket changed (atoms crossed a threshold).
        const LOD          = 1 << 8;
        /// Object-level draw mask changed without changing per-atom rep bits.
        const DRAW_MASK    = 1 << 9;

        /// Convenience: every defined bit.
        const ALL = Self::COORDS.bits()
                  | Self::COLOR.bits()
                  | Self::REPS.bits()
                  | Self::SELECTION.bits()
                  | Self::TRANSPARENCY.bits()
                  | Self::VISIBILITY.bits()
                  | Self::HOVER.bits()
                  | Self::TOPOLOGY.bits()
                  | Self::LOD.bits()
                  | Self::DRAW_MASK.bits();
    }
}

impl DirtyFlags {
    /// Bits that affect *only* per-atom LUTs (color / mask / params), not
    /// instance geometry. When `dirty & !LUT_ONLY_MASK == 0` and the rep
    /// already exists, a representation can take the LUT-only fast path:
    /// flush colours / masks / alpha and skip atom iteration entirely.
    pub const LUT_ONLY_MASK: Self = Self::COLOR
        .union(Self::TRANSPARENCY)
        .union(Self::VISIBILITY)
        .union(Self::SELECTION)
        .union(Self::HOVER);

    /// True iff every set bit lies inside `LUT_ONLY_MASK`.
    pub fn is_lut_only(self) -> bool {
        !self.is_empty() && (self - Self::LUT_ONLY_MASK).is_empty()
    }
}
