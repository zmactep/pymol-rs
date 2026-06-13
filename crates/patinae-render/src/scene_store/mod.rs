//! `SceneStore` — group 2, scene-wide. The single source of truth for
//! per-atom data shared across every representation.
//!
//! One `SceneStore` per `RenderState`. Holds **every** object's atom soup,
//! bond graph, coordinates, and per-atom LUTs in linearly-packed storage
//! buffers; one [`ObjectEntry`] per object indexes those slabs via
//! `atom_offset / bond_offset`. Drawing rebinds via `set_bind_group(2,
//! &scene.bg, &[obj.dynamic_offset])` once per object — N×M
//! bind-group churn collapses to N.
//!
//! Layout invariants:
//!
//! - `atoms`, `coords`, `color_lut`, `mask_lut`, `marker_lut` are all
//!   indexed by **global** atom id
//!   `gid = obj.atom_offset + atom_index`.
//! - `bonds` and `atom_bond_csr` are indexed by global ids the same way
//!   (`obj.bond_offset + bond_index`).
//! - `obj_table` is a dynamic-uniform array; the host advances
//!   `dynamic_offset` per object inside a draw loop. Each entry is padded
//!   to `min_uniform_buffer_offset_alignment` (256 B on every WebGPU
//!   adapter we target).
//! - All per-rep alpha (`sphere_transparency`, etc.) flows through the
//!   per-rep `material_params.alpha_mul` uniform on group 3 — never baked
//!   into `color_lut`. Per-atom per-rep alpha overrides live on
//!   [`AtomGpu`] (`sphere_alpha`, `stick_alpha`) with `f32::NAN` as the
//!   "no override" sentinel.

pub mod layout;
pub mod marker;
pub mod sync;

use std::collections::HashMap;

use bytemuck::{Pod, Zeroable};

use crate::lut_buffer::{BufferFlushStats, GrowableStorageBuffer};
use crate::picking::ObjectId;
use crate::render_input::ColorLutEntry;

pub use layout::SceneStoreLayout;

/// Required dynamic-uniform-buffer offset alignment. WebGPU spec guarantees
/// adapters expose at least 256; we hard-code that and assert at runtime if
/// the adapter reports a larger value.
pub const OBJECT_ENTRY_STRIDE: u64 = 256;

/// Per-atom GPU soup. **32 B**. Mirror of `AtomGpu` in
/// `shaders/common/scene.wgsl`.
///
/// Host-side baking semantics (populated by `sync::write_atoms_for_object`):
///
/// - `vdw` = `atom.effective_vdw() * atom.repr.sphere_scale.unwrap_or(1.0)`.
///   Per-atom `sphere_scale` override is folded in here so changing
///   `settings.sphere.scale` slider only rewrites the per-rep params
///   uniform — not the entire atom buffer. The compute shader applies
///   `params.scale` on top.
/// - `repr_flags` = `atom.repr.visible_reps.0`. These bits are rewritten on
///   both `REPS` and `VISIBILITY`, because surface/mesh/cartoon shaders use
///   them for per-representation visibility without rebuilding geometry.
///   Reserved high bits will carry override-present flags as more reps gain
///   per-atom overrides.
/// - `alpha_pack_a` / `alpha_pack_b` — packed per-atom per-rep alpha
///   overrides, **one byte per rep kind**. Sentinel `0xFF` means
///   "no override → use the per-rep `<rep>_params.alpha_mul` from
///   group 3". Real overrides span `0..=254` mapping to alpha
///   `0.0..=1.0` (precision ≈ 0.4%, perceptually indistinguishable).
///   Layout (LE byte order, byte 0 = LSB):
///   `alpha_pack_a` bytes [sphere | stick | ellipsoid | cartoon]
///   `alpha_pack_b` bytes [surface |   _   |   _    |   _   ]
///   3 spare slots in `alpha_pack_b` for future per-rep alphas without
///   widening the AtomGpu struct.
/// - `element_id` / `chain_id` / `residue_id` are placeholder slots for
///   future culling / grouping. Currently zero.
#[repr(C, align(16))]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct AtomGpu {
    pub vdw: f32,
    pub repr_flags: u32,
    pub alpha_pack_a: u32,
    pub alpha_pack_b: u32,
    pub element_id: u32,
    pub chain_id: u32,
    pub residue_id: u32,
    pub _pad: u32,
}

/// Sentinel byte in `AtomGpu::alpha_pack_*`. Means "no override for this
/// rep kind — use the per-rep `alpha_mul` uniform".
pub const ALPHA_NO_OVERRIDE: u8 = 0xFF;

/// Pack a per-atom transparency override (`atom.repr.<rep>_transparency`)
/// into a byte slot. `None` ⇒ `ALPHA_NO_OVERRIDE`. `Some(t)` is converted
/// to alpha `1.0 - t`, scaled to `0..=254` (rounded). Value `255` is
/// reserved for the sentinel; actual alpha `1.0` clamps to `254`
/// (≈0.996, perceptually identical).
#[inline]
pub fn pack_alpha_override(override_transparency: Option<f32>) -> u8 {
    match override_transparency {
        None => ALPHA_NO_OVERRIDE,
        Some(t) => {
            let alpha = (1.0 - t.clamp(0.0, 1.0)).clamp(0.0, 1.0);
            (alpha * 254.0).round().clamp(0.0, 254.0) as u8
        }
    }
}

const _: () = assert!(std::mem::size_of::<AtomGpu>() == 32);

/// Per-bond GPU entry. **32 B**. Mirror of `BondGpu` in
/// `shaders/common/scene.wgsl`.
///
/// `atom1` / `atom2` are **object-local** atom indices. The compute build
/// shader resolves them to global ids via `obj.atom_offset` so per-object
/// bond data stays relocatable. `valence_perp_pad.xyz` is the precomputed
/// local-plane offset direction used by sticks for multiple bonds.
#[repr(C, align(16))]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct BondGpu {
    pub atom1: u32,
    pub atom2: u32,
    pub order: u32,
    pub flags: u32,
    pub valence_perp_pad: [f32; 4],
}

const _: () = assert!(std::mem::size_of::<BondGpu>() == 32);

/// Per-object table entry. Bound as a `<uniform>` with `dynamic_offset`,
/// so each entry is padded to `OBJECT_ENTRY_STRIDE` (256 B). The Pod
/// payload is 96 B; the trailing 160 B are pad.
#[repr(C, align(16))]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct ObjectEntry {
    pub atom_offset: u32,
    pub atom_count: u32,
    pub bond_offset: u32,
    pub bond_count: u32,
    pub object_id: u32,
    pub flags: u32,
    pub _pad0: [u32; 2],
    pub model_matrix: [[f32; 4]; 4],
}

const _: () = assert!(std::mem::size_of::<ObjectEntry>() == 96);

impl ObjectEntry {
    /// Padded byte size used as `dynamic_offset` step. 256 B per WebGPU.
    pub const STRIDE: u64 = OBJECT_ENTRY_STRIDE;
}

/// Per-object slot bookkeeping. The store keeps one slot per `ObjectId`,
/// allocated linearly the first time the object is seen. Slot reuse on
/// object removal is not implemented yet; reps churning objects
/// would leak buffer space until a future compaction pass.
#[derive(Debug, Clone, Copy)]
pub struct ObjectSlot {
    pub atom_offset: u32,
    pub atom_count: u32,
    pub bond_offset: u32,
    pub bond_count: u32,
    /// Index into the dynamic-uniform table. `dynamic_offset = idx *
    /// ObjectEntry::STRIDE`.
    pub table_index: u32,
}

impl ObjectSlot {
    pub fn dynamic_offset(&self) -> u32 {
        self.table_index * ObjectEntry::STRIDE as u32
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SceneStoreFragmentationStats {
    pub live_atoms: u64,
    pub allocated_atoms: u64,
    pub orphaned_atoms: u64,
    pub live_bonds: u64,
    pub allocated_bonds: u64,
    pub orphaned_bonds: u64,
    pub live_table_slots: u64,
    pub allocated_table_slots: u64,
    pub orphaned_table_slots: u64,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SceneStoreFlushStats {
    pub marker_lut_upload_ranges: u32,
    pub marker_lut_upload_bytes: u64,
    pub marker_lut_reallocated: bool,
    pub fragmentation: SceneStoreFragmentationStats,
}

/// Group-2 backing store. One per `RenderState`.
pub struct SceneStore {
    pub atoms: GrowableStorageBuffer<AtomGpu>,
    pub coords: GrowableStorageBuffer<[f32; 4]>,
    pub bonds: GrowableStorageBuffer<BondGpu>,
    pub color_lut: GrowableStorageBuffer<ColorLutEntry>,
    pub mask_lut: GrowableStorageBuffer<u32>,
    pub marker_lut: GrowableStorageBuffer<u32>,
    /// CSR-style atom→bond adjacency. Two parallel buffers:
    /// - `csr_offsets[gid] .. csr_offsets[gid+1]` — slice into
    ///   `csr_indices` of bond global ids touching the atom.
    /// - `csr_indices[k]` = global bond id.
    pub csr_offsets: GrowableStorageBuffer<u32>,
    pub csr_indices: GrowableStorageBuffer<u32>,
    /// Dynamic-uniform `ObjectEntry` table. Each entry occupies
    /// [`ObjectEntry::STRIDE`] bytes; the payload sits at the start of
    /// each slot, the rest is pad.
    obj_table_cpu: Vec<u8>,
    obj_table_gpu: Option<wgpu::Buffer>,
    obj_table_capacity_bytes: u64,
    obj_table_dirty: Option<(usize, usize)>,

    slots: HashMap<u32, ObjectSlot>,
    /// Linear allocation watermarks. Slot reclamation is left for a future
    /// compaction pass; object adds/removes are infrequent on the current UI path.
    next_atom_offset: u32,
    next_bond_offset: u32,
    next_table_index: u32,
    /// Ranges made unreachable by replacing an existing `ObjectId` with a
    /// differently-sized object. A future compaction pass may reclaim these
    /// ranges, but draw/pick paths already iterate only entries in `slots`.
    orphaned_atom_count: u64,
    orphaned_bond_count: u64,
    orphaned_table_slots: u64,

    bind_group: Option<wgpu::BindGroup>,
    object_coords_bind_group: Option<wgpu::BindGroup>,
    /// `(atoms_cap, coords_cap, bonds_cap, color_cap, mask_cap, marker_cap,
    ///   csr_off_cap, csr_idx_cap, obj_table_cap)` snapshot from the last
    /// bind-group rebuild. Any mismatch triggers a rebuild.
    bind_capacities: [usize; 9],
    object_coords_bind_capacities: [usize; 2],
}

impl SceneStore {
    pub fn new() -> Self {
        Self {
            atoms: GrowableStorageBuffer::new("patinae.scene_store.atoms"),
            coords: GrowableStorageBuffer::new("patinae.scene_store.coords"),
            bonds: GrowableStorageBuffer::new("patinae.scene_store.bonds"),
            color_lut: GrowableStorageBuffer::new("patinae.scene_store.color_lut"),
            mask_lut: GrowableStorageBuffer::new("patinae.scene_store.mask_lut"),
            marker_lut: GrowableStorageBuffer::new("patinae.scene_store.marker_lut"),
            csr_offsets: GrowableStorageBuffer::new("patinae.scene_store.csr_offsets"),
            csr_indices: GrowableStorageBuffer::new("patinae.scene_store.csr_indices"),
            obj_table_cpu: Vec::new(),
            obj_table_gpu: None,
            obj_table_capacity_bytes: 0,
            obj_table_dirty: None,
            slots: HashMap::new(),
            next_atom_offset: 0,
            next_bond_offset: 0,
            next_table_index: 0,
            orphaned_atom_count: 0,
            orphaned_bond_count: 0,
            orphaned_table_slots: 0,
            bind_group: None,
            object_coords_bind_group: None,
            bind_capacities: [0; 9],
            object_coords_bind_capacities: [0; 2],
        }
    }

    pub fn slot(&self, object_id: ObjectId) -> Option<&ObjectSlot> {
        self.slots.get(&object_id.0)
    }

    pub fn has_slot(&self, object_id: ObjectId) -> bool {
        self.slots.contains_key(&object_id.0)
    }

    pub fn bind_group(&self) -> Option<&wgpu::BindGroup> {
        self.bind_group.as_ref()
    }

    pub fn object_coords_bind_group(&self) -> Option<&wgpu::BindGroup> {
        self.object_coords_bind_group.as_ref()
    }

    pub fn obj_table_buffer(&self) -> Option<&wgpu::Buffer> {
        self.obj_table_gpu.as_ref()
    }

    pub fn marker_lut_buffer(&self) -> Option<&wgpu::Buffer> {
        self.marker_lut.buffer()
    }

    pub fn iter_atom_offsets(&self) -> impl Iterator<Item = (u32, u32)> + '_ {
        self.slots
            .iter()
            .map(|(&object_id, slot)| (object_id, slot.atom_offset))
    }

    pub(crate) fn fragmentation_stats(&self) -> SceneStoreFragmentationStats {
        let live_atoms = self
            .slots
            .values()
            .map(|slot| u64::from(slot.atom_count))
            .sum();
        let live_bonds = self
            .slots
            .values()
            .map(|slot| u64::from(slot.bond_count))
            .sum();
        SceneStoreFragmentationStats {
            live_atoms,
            allocated_atoms: u64::from(self.next_atom_offset),
            orphaned_atoms: self.orphaned_atom_count,
            live_bonds,
            allocated_bonds: u64::from(self.next_bond_offset),
            orphaned_bonds: self.orphaned_bond_count,
            live_table_slots: self.slots.len() as u64,
            allocated_table_slots: u64::from(self.next_table_index),
            orphaned_table_slots: self.orphaned_table_slots,
        }
    }

    /// Allocate (or look up) a slot for `object_id` of the given size.
    ///
    /// If a slot already exists with the same size, return it as-is. If the
    /// size changed (host reloaded the molecule, swapped objects under the
    /// same id, etc.) the old slot is **orphaned** — its `atom_offset` /
    /// `bond_offset` range stays as dead space in the linear buffers — and
    /// a fresh slot is appended at the current tail. The buffer-level
    /// orphan space is acceptable for now: a future compaction pass will
    /// reclaim it (`slots.iter` then sees only live ranges so picking and
    /// draw loops are unaffected).
    pub(crate) fn ensure_slot(
        &mut self,
        object_id: ObjectId,
        atom_count: u32,
        bond_count: u32,
    ) -> ObjectSlot {
        if let Some(existing) = self.slots.get(&object_id.0).copied() {
            if existing.atom_count == atom_count && existing.bond_count == bond_count {
                return existing;
            }
            // Size mismatch: drop the stale slot, fall through to a fresh
            // allocation at the buffer tail. The previous range
            // (`existing.atom_offset .. +existing.atom_count`) is now
            // orphaned. `table_index` is also re-assigned; the old
            // dynamic-uniform slot becomes unused (kept at its last value,
            // harmless since no rep references it anymore).
            self.orphaned_atom_count = self
                .orphaned_atom_count
                .checked_add(u64::from(existing.atom_count))
                .expect("scene_store: orphaned atom counter overflow");
            self.orphaned_bond_count = self
                .orphaned_bond_count
                .checked_add(u64::from(existing.bond_count))
                .expect("scene_store: orphaned bond counter overflow");
            self.orphaned_table_slots = self
                .orphaned_table_slots
                .checked_add(1)
                .expect("scene_store: orphaned table counter overflow");
            self.slots.remove(&object_id.0);
        }
        let slot = ObjectSlot {
            atom_offset: self.next_atom_offset,
            atom_count,
            bond_offset: self.next_bond_offset,
            bond_count,
            table_index: self.next_table_index,
        };
        self.next_atom_offset = self
            .next_atom_offset
            .checked_add(atom_count)
            .expect("scene_store: atom offset overflow");
        self.next_bond_offset = self
            .next_bond_offset
            .checked_add(bond_count)
            .expect("scene_store: bond offset overflow");
        self.next_table_index += 1;

        // Grow the per-buffer CPU staging vectors so subsequent writes can
        // index into them without bounds checks.
        let atom_total = self.next_atom_offset as usize;
        let bond_total = self.next_bond_offset as usize;
        self.atoms.resize(atom_total, AtomGpu::zeroed());
        self.coords.resize(atom_total, [0.0; 4]);
        self.color_lut.resize(atom_total, ColorLutEntry::default());
        let n_mask_words = atom_total.div_ceil(32);
        self.mask_lut.resize(n_mask_words, u32::MAX);
        self.marker_lut.resize(atom_total, 0);
        self.bonds.resize(bond_total, BondGpu::zeroed());
        // CSR offsets has one extra entry (sentinel) so `csr[gid+1] -
        // csr[gid]` always works.
        self.csr_offsets.resize(atom_total + 1, 0);
        // CSR indices length is twice the bond count (each bond touches
        // two atoms). Capacity grows; actual size set at write time.
        self.csr_indices.resize(2 * bond_total, 0);

        self.slots.insert(object_id.0, slot);
        self.grow_obj_table(self.next_table_index as usize);
        slot
    }

    fn grow_obj_table(&mut self, count: usize) {
        let needed = count * ObjectEntry::STRIDE as usize;
        if needed > self.obj_table_cpu.len() {
            let prev = self.obj_table_cpu.len();
            self.obj_table_cpu.resize(needed, 0);
            self.mark_obj_table_dirty(prev, needed);
        }
    }

    pub(crate) fn write_obj_entry(&mut self, table_index: u32, entry: ObjectEntry) {
        let offset = table_index as usize * ObjectEntry::STRIDE as usize;
        let payload = bytemuck::bytes_of(&entry);
        self.obj_table_cpu[offset..offset + payload.len()].copy_from_slice(payload);
        self.mark_obj_table_dirty(offset, offset + ObjectEntry::STRIDE as usize);
    }

    fn mark_obj_table_dirty(&mut self, lo: usize, hi: usize) {
        self.obj_table_dirty = Some(match self.obj_table_dirty {
            Some((a, b)) => (a.min(lo), b.max(hi)),
            None => (lo, hi),
        });
    }

    /// Upload all pending writes and (re)build the bind group. Idempotent
    /// when nothing changed since the last call.
    pub fn flush(
        &mut self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        layout: &SceneStoreLayout,
    ) -> SceneStoreFlushStats {
        let atoms_grew = self.atoms.ensure_capacity(device);
        let coords_grew = self.coords.ensure_capacity(device);
        let bonds_grew = self.bonds.ensure_capacity(device);
        let color_grew = self.color_lut.ensure_capacity(device);
        let mask_grew = self.mask_lut.ensure_capacity(device);
        let marker_grew = self.marker_lut.ensure_capacity(device);
        let csr_off_grew = self.csr_offsets.ensure_capacity(device);
        let csr_idx_grew = self.csr_indices.ensure_capacity(device);
        let obj_table_grew = self.ensure_obj_table_capacity(device);

        self.atoms.flush(queue);
        self.coords.flush(queue);
        self.bonds.flush(queue);
        self.color_lut.flush(queue);
        self.mask_lut.flush(queue);
        let marker_flush = self.marker_lut.flush(queue);
        self.csr_offsets.flush(queue);
        self.csr_indices.flush(queue);
        if let Some((lo, hi)) = self.obj_table_dirty.take() {
            if let Some(buf) = self.obj_table_gpu.as_ref() {
                queue.write_buffer(buf, lo as u64, &self.obj_table_cpu[lo..hi]);
            }
        }

        let caps = [
            self.atoms.capacity_entries(),
            self.coords.capacity_entries(),
            self.bonds.capacity_entries(),
            self.color_lut.capacity_entries(),
            self.mask_lut.capacity_entries(),
            self.marker_lut.capacity_entries(),
            self.csr_offsets.capacity_entries(),
            self.csr_indices.capacity_entries(),
            self.obj_table_capacity_bytes as usize,
        ];
        let object_coords_caps = [
            self.coords.capacity_entries(),
            self.obj_table_capacity_bytes as usize,
        ];
        let any_grew = atoms_grew
            || coords_grew
            || bonds_grew
            || color_grew
            || mask_grew
            || marker_grew
            || csr_off_grew
            || csr_idx_grew
            || obj_table_grew;
        if self.bind_group.is_none() || any_grew || caps != self.bind_capacities {
            self.rebuild_bind_group(device, layout);
            self.bind_capacities = caps;
        }
        if self.object_coords_bind_group.is_none()
            || coords_grew
            || obj_table_grew
            || object_coords_caps != self.object_coords_bind_capacities
        {
            self.rebuild_object_coords_bind_group(device, layout);
            self.object_coords_bind_capacities = object_coords_caps;
        }

        SceneStoreFlushStats::from_marker_flush(
            marker_flush,
            marker_grew,
            self.fragmentation_stats(),
        )
    }

    fn ensure_obj_table_capacity(&mut self, device: &wgpu::Device) -> bool {
        let needed = self.obj_table_cpu.len().max(ObjectEntry::STRIDE as usize) as u64;
        if needed <= self.obj_table_capacity_bytes && self.obj_table_gpu.is_some() {
            return false;
        }
        let new_cap = needed.next_power_of_two().max(ObjectEntry::STRIDE * 8);
        self.obj_table_gpu = Some(device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("patinae.scene_store.obj_table"),
            size: new_cap,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        }));
        self.obj_table_capacity_bytes = new_cap;
        // Whole CPU range becomes dirty after realloc.
        if !self.obj_table_cpu.is_empty() {
            self.obj_table_dirty = Some((0, self.obj_table_cpu.len()));
        }
        true
    }

    fn rebuild_bind_group(&mut self, device: &wgpu::Device, layout: &SceneStoreLayout) {
        let (
            Some(atoms),
            Some(coords),
            Some(bonds),
            Some(color),
            Some(mask),
            Some(marker),
            Some(csr_off),
            Some(csr_idx),
            Some(obj_table),
        ) = (
            self.atoms.buffer(),
            self.coords.buffer(),
            self.bonds.buffer(),
            self.color_lut.buffer(),
            self.mask_lut.buffer(),
            self.marker_lut.buffer(),
            self.csr_offsets.buffer(),
            self.csr_indices.buffer(),
            self.obj_table_gpu.as_ref(),
        )
        else {
            return;
        };
        self.bind_group = Some(device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("patinae.scene_store.bind_group"),
            layout: &layout.bind_group_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: wgpu::BindingResource::Buffer(wgpu::BufferBinding {
                        buffer: obj_table,
                        offset: 0,
                        size: std::num::NonZeroU64::new(std::mem::size_of::<ObjectEntry>() as u64),
                    }),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: atoms.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: coords.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: bonds.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 4,
                    resource: color.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 5,
                    resource: mask.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 6,
                    resource: marker.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 7,
                    resource: csr_off.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 8,
                    resource: csr_idx.as_entire_binding(),
                },
            ],
        }));
    }

    fn rebuild_object_coords_bind_group(
        &mut self,
        device: &wgpu::Device,
        layout: &SceneStoreLayout,
    ) {
        let (Some(coords), Some(obj_table)) = (self.coords.buffer(), self.obj_table_gpu.as_ref())
        else {
            return;
        };
        self.object_coords_bind_group =
            Some(device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("patinae.scene_store.object_coords_bind_group"),
                layout: &layout.object_coords_bind_group_layout,
                entries: &[
                    wgpu::BindGroupEntry {
                        binding: 0,
                        resource: wgpu::BindingResource::Buffer(wgpu::BufferBinding {
                            buffer: obj_table,
                            offset: 0,
                            size: std::num::NonZeroU64::new(
                                std::mem::size_of::<ObjectEntry>() as u64
                            ),
                        }),
                    },
                    wgpu::BindGroupEntry {
                        binding: 2,
                        resource: coords.as_entire_binding(),
                    },
                ],
            }));
    }
}

impl SceneStoreFlushStats {
    fn from_marker_flush(
        marker: BufferFlushStats,
        marker_lut_reallocated: bool,
        fragmentation: SceneStoreFragmentationStats,
    ) -> Self {
        Self {
            marker_lut_upload_ranges: marker.ranges,
            marker_lut_upload_bytes: marker.bytes,
            marker_lut_reallocated,
            fragmentation,
        }
    }
}

impl Default for SceneStore {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resizing_same_object_id_records_orphaned_slots() {
        let mut store = SceneStore::new();
        let object = ObjectId(7);

        let first = store.ensure_slot(object, 3, 2);
        assert_eq!(first.atom_offset, 0);
        assert_eq!(first.bond_offset, 0);

        let second = store.ensure_slot(object, 5, 4);
        assert_eq!(second.atom_offset, 3);
        assert_eq!(second.bond_offset, 2);

        let stats = store.fragmentation_stats();
        assert_eq!(stats.live_atoms, 5);
        assert_eq!(stats.allocated_atoms, 8);
        assert_eq!(stats.orphaned_atoms, 3);
        assert_eq!(stats.live_bonds, 4);
        assert_eq!(stats.allocated_bonds, 6);
        assert_eq!(stats.orphaned_bonds, 2);
        assert_eq!(stats.live_table_slots, 1);
        assert_eq!(stats.allocated_table_slots, 2);
        assert_eq!(stats.orphaned_table_slots, 1);
    }
}
