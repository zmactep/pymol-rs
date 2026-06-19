//! Shared GPU buffer owner for compute-built, camera-culled reps.

use crate::compute::cull::{CullParams, CullPipeline};
use crate::memory::{buffer_usage, GpuMemoryUsage};
use crate::representations::{prepare_raw_shadow_indirect, CullPlan, CullPlanCtx};

pub(crate) struct CullableBuffers {
    label: &'static str,
    instance_size: u64,
    raw_instance_buffer: Option<wgpu::Buffer>,
    compacted_instance_buffer: Option<wgpu::Buffer>,
    instance_capacity: u64,
    instance_capacity_count: u32,
    raw_count_buffer: Option<wgpu::Buffer>,
    indirect_buffer: Option<wgpu::Buffer>,
    shadow_indirect_buffer: Option<wgpu::Buffer>,
    cull_params_buffer: wgpu::Buffer,
    cull_bind_group: Option<wgpu::BindGroup>,
    has_shadow_indirect: bool,
}

impl CullableBuffers {
    pub(crate) fn new(
        device: &wgpu::Device,
        label: &'static str,
        instance_size: u64,
        has_shadow_indirect: bool,
    ) -> Self {
        let cull_params_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some(&buffer_label(label, "cull_params")),
            size: CullParams::SIZE,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        Self {
            label,
            instance_size,
            raw_instance_buffer: None,
            compacted_instance_buffer: None,
            instance_capacity: 0,
            instance_capacity_count: 0,
            raw_count_buffer: None,
            indirect_buffer: None,
            shadow_indirect_buffer: None,
            cull_params_buffer,
            cull_bind_group: None,
            has_shadow_indirect,
        }
    }

    pub(crate) fn ensure(
        &mut self,
        device: &wgpu::Device,
        instance_capacity_bytes: u64,
        cull_pipeline: Option<&CullPipeline>,
    ) -> bool {
        let mut storage_changed = false;
        if self.instance_capacity < instance_capacity_bytes || self.raw_instance_buffer.is_none() {
            self.raw_instance_buffer = Some(device.create_buffer(&wgpu::BufferDescriptor {
                label: Some(&buffer_label(self.label, "raw_instances")),
                size: instance_capacity_bytes,
                usage: wgpu::BufferUsages::VERTEX
                    | wgpu::BufferUsages::STORAGE
                    | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            }));
            self.compacted_instance_buffer = Some(device.create_buffer(&wgpu::BufferDescriptor {
                label: Some(&buffer_label(self.label, "compacted_instances")),
                size: instance_capacity_bytes,
                usage: wgpu::BufferUsages::VERTEX
                    | wgpu::BufferUsages::STORAGE
                    | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            }));
            self.instance_capacity = instance_capacity_bytes;
            self.instance_capacity_count = (instance_capacity_bytes / self.instance_size) as u32;
            storage_changed = true;
        }

        if self.raw_count_buffer.is_none() {
            self.raw_count_buffer = Some(device.create_buffer(&wgpu::BufferDescriptor {
                label: Some(&buffer_label(self.label, "raw_count")),
                size: 4,
                usage: wgpu::BufferUsages::STORAGE
                    | wgpu::BufferUsages::COPY_DST
                    | wgpu::BufferUsages::COPY_SRC,
                mapped_at_creation: false,
            }));
            storage_changed = true;
        }

        if self.indirect_buffer.is_none() {
            self.indirect_buffer = Some(device.create_buffer(&wgpu::BufferDescriptor {
                label: Some(&buffer_label(self.label, "indirect")),
                size: 16,
                usage: wgpu::BufferUsages::INDIRECT
                    | wgpu::BufferUsages::STORAGE
                    | wgpu::BufferUsages::COPY_DST
                    | wgpu::BufferUsages::COPY_SRC,
                mapped_at_creation: false,
            }));
            storage_changed = true;
        }

        if self.has_shadow_indirect && self.shadow_indirect_buffer.is_none() {
            self.shadow_indirect_buffer = Some(device.create_buffer(&wgpu::BufferDescriptor {
                label: Some(&buffer_label(self.label, "shadow_indirect")),
                size: 16,
                usage: wgpu::BufferUsages::INDIRECT | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            }));
        }

        if let Some(cull) = cull_pipeline {
            if storage_changed || self.cull_bind_group.is_none() {
                self.cull_bind_group = Some(
                    cull.make_rep_bind_group(
                        device,
                        &self.cull_params_buffer,
                        self.raw_instance_buffer
                            .as_ref()
                            .expect("raw instance buffer"),
                        self.raw_count_buffer.as_ref().expect("raw count buffer"),
                        self.compacted_instance_buffer
                            .as_ref()
                            .expect("compacted instance buffer"),
                        self.indirect_buffer.as_ref().expect("indirect buffer"),
                    ),
                );
            }
        }

        storage_changed
    }

    pub(crate) fn raw_instance_buffer(&self) -> Option<&wgpu::Buffer> {
        self.raw_instance_buffer.as_ref()
    }

    pub(crate) fn compacted_instance_buffer(&self) -> Option<&wgpu::Buffer> {
        self.compacted_instance_buffer.as_ref()
    }

    pub(crate) fn instance_stride(&self) -> u64 {
        self.instance_size
    }

    pub(crate) fn instance_capacity(&self) -> u64 {
        self.instance_capacity
    }

    pub(crate) fn raw_count_buffer(&self) -> Option<&wgpu::Buffer> {
        self.raw_count_buffer.as_ref()
    }

    pub(crate) fn indirect_buffer(&self) -> Option<&wgpu::Buffer> {
        self.indirect_buffer.as_ref()
    }

    pub(crate) fn shadow_indirect_buffer(&self) -> Option<&wgpu::Buffer> {
        self.shadow_indirect_buffer.as_ref()
    }

    pub(crate) fn cull_bind_group(&self) -> Option<&wgpu::BindGroup> {
        self.cull_bind_group.as_ref()
    }

    pub(crate) fn cull_params_buffer(&self) -> &wgpu::Buffer {
        &self.cull_params_buffer
    }

    pub(crate) fn instance_capacity_count(&self) -> u32 {
        self.instance_capacity_count
    }

    pub(crate) fn has_raw_instances(&self) -> bool {
        self.raw_instance_buffer.is_some()
    }

    pub(crate) fn memory_usage(&self) -> GpuMemoryUsage {
        let mut usage = GpuMemoryUsage::default();
        usage.add(buffer_usage(&self.cull_params_buffer));
        if let Some(buffer) = self.raw_instance_buffer.as_ref() {
            usage.add(buffer_usage(buffer));
        }
        if let Some(buffer) = self.compacted_instance_buffer.as_ref() {
            usage.add(buffer_usage(buffer));
        }
        if let Some(buffer) = self.raw_count_buffer.as_ref() {
            usage.add(buffer_usage(buffer));
        }
        if let Some(buffer) = self.indirect_buffer.as_ref() {
            usage.add(buffer_usage(buffer));
        }
        if let Some(buffer) = self.shadow_indirect_buffer.as_ref() {
            usage.add(buffer_usage(buffer));
        }
        usage
    }

    pub(crate) fn reset_raw_count(&self, queue: &wgpu::Queue) {
        if let Some(buf) = self.raw_count_buffer.as_ref() {
            queue.write_buffer(buf, 0, bytemuck::bytes_of(&0u32));
        }
    }

    pub(crate) fn prepare_raw_shadow_indirect(
        &self,
        encoder: &mut wgpu::CommandEncoder,
        queue: &wgpu::Queue,
        seed: &[u32; 4],
    ) {
        prepare_raw_shadow_indirect(
            encoder,
            queue,
            self.raw_count_buffer.as_ref(),
            self.shadow_indirect_buffer.as_ref(),
            seed,
        );
    }

    pub(crate) fn plan_cull(
        &self,
        ctx: &CullPlanCtx<'_>,
        upper: u32,
        kind_radius: f32,
        seed: &[u32; 4],
    ) -> Option<CullPlan<'_>> {
        if upper == 0 || self.cull_bind_group.is_none() {
            return None;
        }
        let params = CullParams {
            view_proj: ctx.view_proj,
            frustum_planes: ctx.frustum_planes,
            raw_capacity: self.instance_capacity_count,
            kind_radius,
            _pad0: 0,
            _pad1: 0,
        };
        ctx.queue
            .write_buffer(&self.cull_params_buffer, 0, bytemuck::bytes_of(&params));
        if let Some(buf) = self.indirect_buffer.as_ref() {
            ctx.queue.write_buffer(buf, 0, bytemuck::cast_slice(seed));
        }
        Some(CullPlan {
            bind_group: self.cull_bind_group.as_ref().expect("cull bind group"),
            upper,
        })
    }
}

fn buffer_label(rep: &str, name: &str) -> String {
    format!("patinae.{rep}.{name}")
}
