//! Lite selected atom dots.

use std::collections::{BTreeMap, HashSet};
use std::num::NonZeroU64;

use bytemuck::{Pod, Zeroable};
use patinae_mol::DirtyFlags;

use crate::context::RenderContext;
use crate::frame::DEPTH_FORMAT;
use crate::memory::{buffer_usage, GpuMemoryUsage};
use crate::memory_policy::{RenderMemoryPolicy, RenderMemoryProfile};
use crate::picking::ObjectId;
use crate::render_input::MarkerUpdate;
use crate::scene_store::marker::MARKER_SELECTED;
use crate::scene_store::{ObjectSlot, SceneStore, SceneStoreLayout};
use crate::shader_source;

const SELECTED_INDEX_MIN_CAPACITY: usize = 16;
const DOT_RADIUS_BASE_PX: f32 = 7.0;
const DOT_RADIUS_PER_MARKING_PX: f32 = 3.0;
const DOT_RADIUS_MIN_PX: f32 = 8.0;
const DOT_RADIUS_MAX_PX: f32 = 20.0;
// Push the dot center slightly toward the camera so selected atoms remain
// visible on top of their own molecular surface without becoming X-ray marks.
const DOT_VIEW_BIAS_ANGSTROM: f32 = 0.25;

#[repr(C, align(16))]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
struct SelectionDotsParams {
    radius_px: f32,
    view_bias: f32,
    _pad0: f32,
    _pad1: f32,
}

impl SelectionDotsParams {
    const SIZE: u64 = std::mem::size_of::<Self>() as u64;
}

struct SelectionDotsObject {
    selected_indices: Vec<u32>,
    buffer: wgpu::Buffer,
    bind_group: wgpu::BindGroup,
    capacity: usize,
}

impl SelectionDotsObject {
    fn selected_count(&self) -> u32 {
        self.selected_indices.len() as u32
    }
}

pub(crate) struct SelectionDotsPass {
    pipeline: wgpu::RenderPipeline,
    bind_group_layout: wgpu::BindGroupLayout,
    params_buffer: wgpu::Buffer,
    objects: BTreeMap<u32, SelectionDotsObject>,
}

impl SelectionDotsPass {
    pub(crate) fn new(ctx: &RenderContext, scene_layout: &SceneStoreLayout) -> Self {
        let device = &ctx.device;
        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("patinae.selection_dots.bgl"),
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::VERTEX,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: NonZeroU64::new(SelectionDotsParams::SIZE),
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 1,
                    visibility: wgpu::ShaderStages::VERTEX,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: true },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
            ],
        });
        let params_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("patinae.selection_dots.params"),
            size: SelectionDotsParams::SIZE,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("patinae.selection_dots.shader"),
            source: wgpu::ShaderSource::Wgsl(
                shader_source::expand(shader_source::SELECTION_DOTS_WGSL).into(),
            ),
        });
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("patinae.selection_dots.pipeline_layout"),
            bind_group_layouts: &[
                Some(&ctx.frame.bind_group_layout),
                Some(&ctx.lighting.bind_group_layout),
                Some(&scene_layout.bind_group_layout),
                Some(&bind_group_layout),
            ],
            immediate_size: 0,
        });
        let color_targets = [Some(wgpu::ColorTargetState {
            format: ctx.color_format,
            blend: Some(wgpu::BlendState::ALPHA_BLENDING),
            write_mask: wgpu::ColorWrites::ALL,
        })];
        let pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("patinae.selection_dots.pipeline"),
            layout: Some(&pipeline_layout),
            vertex: wgpu::VertexState {
                module: &module,
                entry_point: Some("vs_main"),
                compilation_options: Default::default(),
                buffers: &[],
            },
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::TriangleList,
                cull_mode: None,
                ..Default::default()
            },
            depth_stencil: Some(wgpu::DepthStencilState {
                format: DEPTH_FORMAT,
                depth_write_enabled: Some(false),
                depth_compare: Some(wgpu::CompareFunction::LessEqual),
                stencil: wgpu::StencilState::default(),
                bias: wgpu::DepthBiasState::default(),
            }),
            multisample: wgpu::MultisampleState::default(),
            fragment: Some(wgpu::FragmentState {
                module: &module,
                entry_point: Some("fs_main"),
                compilation_options: Default::default(),
                targets: &color_targets,
            }),
            multiview_mask: None,
            cache: None,
        });

        Self {
            pipeline,
            bind_group_layout,
            params_buffer,
            objects: BTreeMap::new(),
        }
    }

    pub(crate) fn upload_params(&self, queue: &wgpu::Queue, marking_width: f32) {
        let radius_px = selection_dot_radius_px(marking_width);
        queue.write_buffer(
            &self.params_buffer,
            0,
            bytemuck::bytes_of(&SelectionDotsParams {
                radius_px,
                view_bias: DOT_VIEW_BIAS_ANGSTROM,
                _pad0: 0.0,
                _pad1: 0.0,
            }),
        );
    }

    pub(crate) fn selected_indices(&self, object_id: u32) -> &[u32] {
        self.objects
            .get(&object_id)
            .map(|object| object.selected_indices.as_slice())
            .unwrap_or(&[])
    }

    pub(crate) fn sync_object(
        &mut self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        object_id: u32,
        marker_bits: &[u32],
    ) -> bool {
        let selected_indices = collect_selected_indices(marker_bits);
        let Some(existing) = self.objects.get_mut(&object_id) else {
            if selected_indices.is_empty() {
                return false;
            }
            let object = make_selection_dots_object(
                device,
                queue,
                &self.bind_group_layout,
                &self.params_buffer,
                &selected_indices,
            );
            self.objects.insert(object_id, object);
            return true;
        };

        if existing.selected_indices == selected_indices {
            return false;
        }
        if selected_indices.is_empty() {
            self.objects.remove(&object_id);
            return true;
        }

        if selected_indices.len() > existing.capacity {
            *existing = make_selection_dots_object(
                device,
                queue,
                &self.bind_group_layout,
                &self.params_buffer,
                &selected_indices,
            );
        } else {
            queue.write_buffer(&existing.buffer, 0, bytemuck::cast_slice(&selected_indices));
            existing.selected_indices = selected_indices;
        }
        true
    }

    pub(crate) fn retain_objects<I>(&mut self, live_object_ids: I) -> bool
    where
        I: IntoIterator<Item = u32>,
    {
        let live: HashSet<u32> = live_object_ids.into_iter().collect();
        let before = self.objects.len();
        self.objects.retain(|object_id, _| live.contains(object_id));
        self.objects.len() != before
    }

    pub(crate) fn clear(&mut self) {
        self.objects.clear();
    }

    pub(crate) fn has_selected_atoms(&self) -> bool {
        !self.objects.is_empty()
    }

    pub(crate) fn memory_usage(&self) -> GpuMemoryUsage {
        let mut usage = buffer_usage(&self.params_buffer);
        for object in self.objects.values() {
            usage.add(buffer_usage(&object.buffer));
        }
        usage
    }

    pub(crate) fn record(
        &self,
        encoder: &mut wgpu::CommandEncoder,
        target: &wgpu::TextureView,
        depth: &wgpu::TextureView,
        frame_bind_group: &wgpu::BindGroup,
        lighting_bind_group: &wgpu::BindGroup,
        scene_store: &SceneStore,
    ) {
        if !self.has_selected_atoms() {
            return;
        }
        let Some(scene_bind_group) = scene_store.bind_group() else {
            return;
        };

        let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
            label: Some("patinae.selection_dots_pass"),
            color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                view: target,
                depth_slice: None,
                resolve_target: None,
                ops: wgpu::Operations {
                    load: wgpu::LoadOp::Load,
                    store: wgpu::StoreOp::Store,
                },
            })],
            depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                view: depth,
                depth_ops: Some(wgpu::Operations {
                    load: wgpu::LoadOp::Load,
                    store: wgpu::StoreOp::Store,
                }),
                stencil_ops: None,
            }),
            timestamp_writes: None,
            occlusion_query_set: None,
            multiview_mask: None,
        });
        pass.set_pipeline(&self.pipeline);
        pass.set_bind_group(0, frame_bind_group, &[]);
        pass.set_bind_group(1, lighting_bind_group, &[]);
        for (&object_id, object) in &self.objects {
            let Some(slot) = scene_store.slot(ObjectId(object_id)) else {
                continue;
            };
            pass.set_bind_group(2, scene_bind_group, &[slot.dynamic_offset()]);
            pass.set_bind_group(3, &object.bind_group, &[]);
            pass.draw(0..6, 0..object.selected_count());
        }
    }
}

fn make_selection_dots_object(
    device: &wgpu::Device,
    queue: &wgpu::Queue,
    bind_group_layout: &wgpu::BindGroupLayout,
    params_buffer: &wgpu::Buffer,
    selected_indices: &[u32],
) -> SelectionDotsObject {
    let capacity = selected_indices
        .len()
        .next_power_of_two()
        .max(SELECTED_INDEX_MIN_CAPACITY);
    let buffer = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("patinae.selection_dots.indices"),
        size: (capacity * std::mem::size_of::<u32>()) as u64,
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    queue.write_buffer(&buffer, 0, bytemuck::cast_slice(selected_indices));
    let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("patinae.selection_dots.bg"),
        layout: bind_group_layout,
        entries: &[
            wgpu::BindGroupEntry {
                binding: 0,
                resource: params_buffer.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 1,
                resource: buffer.as_entire_binding(),
            },
        ],
    });
    SelectionDotsObject {
        selected_indices: selected_indices.to_vec(),
        buffer,
        bind_group,
        capacity,
    }
}

pub(crate) fn uses_selection_dots_fallback(policy: RenderMemoryPolicy) -> bool {
    matches!(
        policy.profile,
        RenderMemoryProfile::Lite | RenderMemoryProfile::Manual { .. }
    ) && !policy.overlays.selection_enabled
}

pub(crate) fn should_rebuild_selected_indices(
    dirty: DirtyFlags,
    old_selected_indices: &[u32],
    marker_updates: &[MarkerUpdate],
) -> bool {
    if dirty.intersects(DirtyFlags::TOPOLOGY | DirtyFlags::SELECTION) {
        return true;
    }
    dirty.intersects(DirtyFlags::HOVER)
        && marker_updates_change_selected(old_selected_indices, marker_updates)
}

pub(crate) fn object_marker_bits(marker_lut: &[u32], slot: ObjectSlot) -> &[u32] {
    let start = slot.atom_offset as usize;
    if start >= marker_lut.len() {
        return &[];
    }
    let end = start
        .saturating_add(slot.atom_count as usize)
        .min(marker_lut.len());
    &marker_lut[start..end]
}

fn collect_selected_indices(marker_bits: &[u32]) -> Vec<u32> {
    marker_bits
        .iter()
        .enumerate()
        .filter_map(|(index, bits)| {
            if bits & MARKER_SELECTED != 0 {
                Some(index as u32)
            } else {
                None
            }
        })
        .collect()
}

fn marker_updates_change_selected(
    old_selected_indices: &[u32],
    marker_updates: &[MarkerUpdate],
) -> bool {
    marker_updates.iter().any(|update| {
        let was_selected = old_selected_indices
            .binary_search(&update.atom_index)
            .is_ok();
        let is_selected = update.bits & MARKER_SELECTED != 0;
        was_selected != is_selected
    })
}

fn selection_dot_radius_px(marking_width: f32) -> f32 {
    (marking_width.clamp(0.5, 20.0) * DOT_RADIUS_PER_MARKING_PX + DOT_RADIUS_BASE_PX)
        .clamp(DOT_RADIUS_MIN_PX, DOT_RADIUS_MAX_PX)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::byte_units::mib_to_bytes;
    use crate::memory_policy::RenderMemoryPolicy;
    use crate::scene_store::marker::MARKER_HOVER;

    #[test]
    fn selected_indices_ignore_hover_bits() {
        assert_eq!(
            collect_selected_indices(&[
                0,
                MARKER_SELECTED,
                MARKER_HOVER,
                MARKER_SELECTED | MARKER_HOVER,
            ]),
            vec![1, 3]
        );
    }

    #[test]
    fn hover_only_dirty_reuses_selected_indices_when_selection_bits_stay_put() {
        let old_selected = [1, 3];
        let updates = [
            MarkerUpdate {
                atom_index: 0,
                bits: MARKER_HOVER,
            },
            MarkerUpdate {
                atom_index: 1,
                bits: MARKER_SELECTED | MARKER_HOVER,
            },
        ];

        assert!(!should_rebuild_selected_indices(
            DirtyFlags::HOVER,
            &old_selected,
            &updates,
        ));
    }

    #[test]
    fn hover_dirty_rebuilds_if_sparse_update_changes_selection_bit() {
        let add_selected = [MarkerUpdate {
            atom_index: 2,
            bits: MARKER_SELECTED | MARKER_HOVER,
        }];
        let clear_selected = [MarkerUpdate {
            atom_index: 3,
            bits: MARKER_HOVER,
        }];

        assert!(should_rebuild_selected_indices(
            DirtyFlags::HOVER,
            &[],
            &add_selected,
        ));
        assert!(should_rebuild_selected_indices(
            DirtyFlags::HOVER,
            &[3],
            &clear_selected,
        ));
    }

    #[test]
    fn selection_and_topology_dirty_rebuild_selected_indices() {
        assert!(should_rebuild_selected_indices(
            DirtyFlags::SELECTION,
            &[],
            &[],
        ));
        assert!(should_rebuild_selected_indices(
            DirtyFlags::TOPOLOGY,
            &[],
            &[],
        ));
        assert!(!should_rebuild_selected_indices(
            DirtyFlags::COLOR,
            &[],
            &[]
        ));
        assert!(!should_rebuild_selected_indices(
            DirtyFlags::empty(),
            &[],
            &[]
        ));
    }

    #[test]
    fn object_marker_slice_clamps_to_marker_lut_capacity() {
        let slot = ObjectSlot {
            atom_offset: 2,
            atom_count: 4,
            bond_offset: 0,
            bond_count: 0,
            table_index: 0,
        };

        assert_eq!(object_marker_bits(&[0, 0, 1, 2], slot), &[1, 2]);
    }

    #[test]
    fn selection_dots_fallback_is_only_lite_and_manual() {
        assert!(!uses_selection_dots_fallback(
            RenderMemoryPolicy::performance()
        ));
        assert!(!uses_selection_dots_fallback(RenderMemoryPolicy::balanced()));
        assert!(uses_selection_dots_fallback(RenderMemoryPolicy::lite()));
        assert!(uses_selection_dots_fallback(RenderMemoryPolicy::manual(
            mib_to_bytes(2048)
        )));
    }

    #[test]
    fn default_selection_dot_radius_is_readable() {
        assert_eq!(selection_dot_radius_px(1.0), 10.0);
        assert!(selection_dot_radius_px(0.0) >= DOT_RADIUS_MIN_PX);
        assert_eq!(selection_dot_radius_px(20.0), DOT_RADIUS_MAX_PX);
    }
}
