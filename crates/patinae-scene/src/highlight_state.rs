//! Per-frame GPU state for the selection / hover highlight pass.
//!
//! Holds two storage buffers — `ranges[(offset, nwords); n_slots]` and
//! `bits[u32]` — that the highlight shader uses to test atom membership in
//! the union of visible selections. Slot indices match the picking pass:
//! `object_id == slot + 1`, with 0 reserved for the cleared-pixel sentinel.
//!
//! Rebuilt every frame in `prepare_scene` from
//! `SelectionManager::evaluate_visible` and the active hover target.

use patinae_select::SelectionResult;
use wgpu::util::DeviceExt;

use crate::object::ObjectRegistry;
use crate::session::HoverTarget;

/// At least 1 word in each storage buffer; wgpu rejects zero-sized bindings.
const MIN_WORDS: usize = 1;
const MIN_RANGES: usize = 1;

/// Per-slot bit range: `(offset_in_words, nwords)`. `nwords == 0` means the
/// slot has no selection bits.
#[repr(C)]
#[derive(Debug, Copy, Clone, Default, bytemuck::Pod, bytemuck::Zeroable)]
struct Range {
    offset_words: u32,
    n_words: u32,
}

/// GPU-side selection / hover state consumed by `HighlightPipeline`.
///
/// Buffers are allocated lazily in `rebuild`, so `Session` can construct an
/// empty state without a device. After the first non-empty rebuild the
/// buffers are valid for binding.
pub struct HighlightState {
    ranges_buffer: Option<wgpu::Buffer>,
    selection_bits_buffer: Option<wgpu::Buffer>,
    hover_bits_buffer: Option<wgpu::Buffer>,
    /// Allocated capacity, in elements.
    ranges_capacity: usize,
    bits_capacity: usize,

    /// True when at least one selection or hover bit is set. The render loop
    /// skips the highlight pass entirely when this is false.
    active: bool,
}

impl Default for HighlightState {
    fn default() -> Self {
        Self::new()
    }
}

impl HighlightState {
    pub fn new() -> Self {
        Self {
            ranges_buffer: None,
            selection_bits_buffer: None,
            hover_bits_buffer: None,
            ranges_capacity: 0,
            bits_capacity: 0,
            active: false,
        }
    }

    /// GPU buffers — `None` until the first `rebuild` call. Callers should
    /// only consult these when `is_active()` is true.
    pub fn ranges_buffer(&self) -> Option<&wgpu::Buffer> {
        self.ranges_buffer.as_ref()
    }

    pub fn selection_bits_buffer(&self) -> Option<&wgpu::Buffer> {
        self.selection_bits_buffer.as_ref()
    }

    pub fn hover_bits_buffer(&self) -> Option<&wgpu::Buffer> {
        self.hover_bits_buffer.as_ref()
    }

    /// True when a highlight should be drawn this frame.
    pub fn is_active(&self) -> bool {
        self.active
    }

    /// Rebuild CPU bit blocks and upload to GPU. Allocates new buffers when
    /// the previous capacity is too small; otherwise reuses the existing ones.
    ///
    /// `slot_names` must be ordered identically to the picking pass — i.e.
    /// `pickable_object_names(session)`. `selections` is the union returned
    /// by `SelectionManager::evaluate_visible`. Selections for objects not in
    /// `slot_names` (disabled molecules) are ignored.
    pub fn rebuild(
        &mut self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        slot_names: &[String],
        selections: &[(String, SelectionResult)],
        registry: &ObjectRegistry,
        hover: Option<&HoverTarget>,
    ) {
        let n_slots = slot_names.len();

        // Per-slot lookup of selection result and (optional) hover selection.
        // Both share the same `ranges` layout so the shader can do two
        // bitmap probes with a single offset/n_words read.
        let mut sel_lookup: Vec<Option<&SelectionResult>> = vec![None; n_slots];
        for (name, sel) in selections {
            if let Some(slot) = slot_names.iter().position(|n| n == name) {
                if !sel.is_empty() {
                    sel_lookup[slot] = Some(sel);
                }
            }
        }
        let hover_slot: Option<usize> =
            hover.and_then(|h| slot_names.iter().position(|n| n == &h.object));
        let hover_sel: Option<&SelectionResult> = match (hover, hover_slot) {
            (Some(h), Some(_)) if !h.selection.is_empty() => Some(&h.selection),
            _ => None,
        };

        // Pass 1: per-slot word count = atom_count.div_ceil(32) when either
        // selection or hover touches this slot. Caps at the molecule's
        // atom_count so trailing bits never address out-of-range atoms.
        let mut ranges: Vec<Range> = Vec::with_capacity(n_slots.max(MIN_RANGES));
        let mut total_words: usize = 0;
        for (slot, name) in slot_names.iter().enumerate() {
            let needs_slot =
                sel_lookup[slot].is_some() || (hover_slot == Some(slot) && hover_sel.is_some());
            let n_words = if needs_slot {
                let atom_count = registry
                    .get_molecule(name)
                    .map(|m| m.molecule().atom_count())
                    .unwrap_or(0);
                atom_count.div_ceil(32)
            } else {
                0
            };
            ranges.push(Range {
                offset_words: total_words as u32,
                n_words: n_words as u32,
            });
            total_words += n_words;
        }

        // Always allocate at least MIN_WORDS to satisfy wgpu's nonzero-size rule.
        let bits_len = total_words.max(MIN_WORDS);
        let mut sel_bits: Vec<u32> = vec![0u32; bits_len];
        let mut hov_bits: Vec<u32> = vec![0u32; bits_len];

        // Pass 2: fill selection bits.
        let mut any_selected = false;
        for slot in 0..n_slots {
            let r = ranges[slot];
            if r.n_words == 0 {
                continue;
            }
            let sel = match sel_lookup[slot] {
                Some(s) => s,
                None => continue,
            };
            let offset = r.offset_words as usize;
            for atom in sel.indices() {
                let i = atom.as_usize();
                let w = i >> 5;
                if w >= r.n_words as usize {
                    continue;
                }
                sel_bits[offset + w] |= 1u32 << (i & 31);
                any_selected = true;
            }
        }

        // Pass 3: fill hover bits at the hover slot.
        let mut any_hover = false;
        if let (Some(slot), Some(hsel)) = (hover_slot, hover_sel) {
            let r = ranges[slot];
            if r.n_words > 0 {
                let offset = r.offset_words as usize;
                for atom in hsel.indices() {
                    let i = atom.as_usize();
                    let w = i >> 5;
                    if w >= r.n_words as usize {
                        continue;
                    }
                    hov_bits[offset + w] |= 1u32 << (i & 31);
                    any_hover = true;
                }
            }
        }

        self.active = any_selected || any_hover;

        // Always need ranges entries; pad to MIN_RANGES if no slots.
        let ranges_to_upload: Vec<Range> = if ranges.is_empty() {
            vec![Range::default(); MIN_RANGES]
        } else {
            ranges
        };

        self.upload(device, queue, &ranges_to_upload, &sel_bits, &hov_bits);
    }

    fn upload(
        &mut self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        ranges: &[Range],
        sel_bits: &[u32],
        hov_bits: &[u32],
    ) {
        debug_assert_eq!(sel_bits.len(), hov_bits.len());

        let need_new_ranges = self.ranges_buffer.is_none() || ranges.len() > self.ranges_capacity;
        if need_new_ranges {
            self.ranges_buffer = Some(device.create_buffer_init(
                &wgpu::util::BufferInitDescriptor {
                    label: Some("Highlight Ranges Buffer"),
                    contents: bytemuck::cast_slice(ranges),
                    usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
                },
            ));
            self.ranges_capacity = ranges.len();
        } else if let Some(buf) = &self.ranges_buffer {
            queue.write_buffer(buf, 0, bytemuck::cast_slice(ranges));
        }

        let need_new_bits =
            self.selection_bits_buffer.is_none() || sel_bits.len() > self.bits_capacity;
        if need_new_bits {
            self.selection_bits_buffer = Some(device.create_buffer_init(
                &wgpu::util::BufferInitDescriptor {
                    label: Some("Highlight Selection Bits Buffer"),
                    contents: bytemuck::cast_slice(sel_bits),
                    usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
                },
            ));
            self.hover_bits_buffer = Some(device.create_buffer_init(
                &wgpu::util::BufferInitDescriptor {
                    label: Some("Highlight Hover Bits Buffer"),
                    contents: bytemuck::cast_slice(hov_bits),
                    usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
                },
            ));
            self.bits_capacity = sel_bits.len();
        } else {
            if let Some(buf) = &self.selection_bits_buffer {
                queue.write_buffer(buf, 0, bytemuck::cast_slice(sel_bits));
            }
            if let Some(buf) = &self.hover_bits_buffer {
                queue.write_buffer(buf, 0, bytemuck::cast_slice(hov_bits));
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn range_size_matches_wgsl_vec2_u32() {
        // WGSL `vec2<u32>` is 8 bytes with std430 alignment.
        assert_eq!(std::mem::size_of::<Range>(), 8);
    }
}
