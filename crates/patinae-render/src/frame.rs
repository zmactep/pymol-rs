//! Per-frame GPU targets used by the WBOIT pipeline.
//!
//! - `accum`   — Rgba16Float, additive blend, holds Σ (color·α·w, α·w)
//! - `reveal`  — R16Float, multiplicative blend, holds Π (1-α)
//! - `depth`   — Depth32Float, read-only during the translucent pass; written
//!   by the opaque pass
//! - `picking` — Rg32Uint @ low-res, packed `(rep_kind, object_id, atom_id)`
//!
//! The WBOIT targets are always available, while picking uses a low-resolution
//! target to reduce per-frame bandwidth.

pub const ACCUM_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Rgba16Float;
pub const REVEAL_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::R16Float;
pub const DEPTH_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Depth32Float;
pub const PICKING_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Rg32Uint;
pub const MARKING_MASK_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Rgba8Unorm;
/// Single-channel occlusion factor written by the SSAO compute pass and
/// consumed by the bilateral blur + final compose. R32Float is in the
/// WebGPU 1.0 mandatory storage-write set; R8Unorm / R16Float are not.
/// 4 B/px is more than AO needs but stays portable without an
/// adapter-specific feature gate.
pub const SSAO_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::R32Float;

/// Fraction of the main viewport used for the picking texture.
///
/// This uses 0.5× to cut bandwidth after profiling showed the picking pass
/// costs as much as the translucent pass on assemblies (60+ chains, 86M
/// cartoon vertices on 7KP3). At 0.5× the picking pass runs at
/// quarter the fragment cost. Thin-feature misses (cartoon ribbons,
/// sticks) are tolerated; click handlers can dilate the lookup to nearby
/// taps if a single-pixel sample misses.
pub const PICKING_SCALE: f32 = 0.5;

/// Per-frame GPU targets. The picking-related targets are `None` when the
/// host constructs `RenderState` with `PickingMode::Disabled` — saving the
/// half-res Rg32Uint + matching depth surface. Visual overlays use their own
/// full-resolution id target so silhouettes and selection highlighting stay
/// independent from hit-test picking.
pub struct FrameTargets {
    pub width: u32,
    pub height: u32,
    pub accum: wgpu::TextureView,
    pub reveal: wgpu::TextureView,
    pub depth: wgpu::TextureView,
    pub picking: Option<wgpu::TextureView>,
    /// Low-res depth that matches `picking` 1:1 (same w×h × `PICKING_SCALE`).
    /// Lets the picking pass resolve overlapping representations to the
    /// nearest hit instead of "last rasterized wins".
    pub picking_depth: Option<wgpu::TextureView>,
    /// Snapshot of the previous full re-record's picking texture. Read by
    /// the reprojection compute pass to warp into the current frame. `None`
    /// when picking is disabled.
    pub picking_prev: Option<wgpu::TextureView>,
    /// Snapshot of the previous full re-record's picking depth. Read by the
    /// reprojection compute pass to unproject source pixels back into world
    /// space. `None` when picking is disabled.
    pub picking_depth_prev: Option<wgpu::TextureView>,
    /// SSAO occlusion-factor texture. Written by `compute/ssao.rs`,
    /// blurred in-place by `compute/ssao_blur.rs` (ping-pong with
    /// `ssao_blurred_texture`), read by `postprocess/ssao_compose`.
    /// `None` until SSAO is first dispatched.
    pub ssao_texture: wgpu::Texture,
    pub ssao_view: wgpu::TextureView,
    pub ssao_blurred_texture: wgpu::Texture,
    pub ssao_blurred_view: wgpu::TextureView,
    /// Off-screen colour target used when FXAA is enabled. Every
    /// render pass writes here; the final FXAA fragment shader reads
    /// it + writes to the host target. Always allocated (cost = 1
    /// host-sized RGBA8 ≈ 8 MB at 1080p), making `set_fxaa(true)` a
    /// zero-rebuild toggle.
    pub color_scratch_texture: wgpu::Texture,
    pub color_scratch_view: wgpu::TextureView,
    /// Full-resolution id target used by visual overlays (silhouettes and
    /// selection highlighting). Allocated lazily only when an overlay needs it.
    pub overlay_id: Option<wgpu::TextureView>,
    pub overlay_id_depth: Option<wgpu::TextureView>,
    pub marking_mask: Option<wgpu::TextureView>,
    pub overlay_color_scratch: Option<wgpu::TextureView>,
    /// Underlying textures kept alive — exposed for advanced callers (e.g.
    /// readback). Most code should bind through the views above.
    pub accum_texture: wgpu::Texture,
    pub reveal_texture: wgpu::Texture,
    pub depth_texture: wgpu::Texture,
    pub picking_texture: Option<wgpu::Texture>,
    pub picking_depth_texture: Option<wgpu::Texture>,
    pub picking_prev_texture: Option<wgpu::Texture>,
    pub picking_depth_prev_texture: Option<wgpu::Texture>,
    pub overlay_id_texture: Option<wgpu::Texture>,
    pub overlay_id_depth_texture: Option<wgpu::Texture>,
    pub marking_mask_texture: Option<wgpu::Texture>,
    pub overlay_color_scratch_texture: Option<wgpu::Texture>,
}

impl FrameTargets {
    /// Default constructor — allocates the picking textures. Most hosts go
    /// through `RenderState::new` which calls `new_with_picking(true)`.
    pub fn new(
        device: &wgpu::Device,
        width: u32,
        height: u32,
        color_format: wgpu::TextureFormat,
    ) -> Self {
        Self::new_with_picking(device, width, height, true, color_format)
    }

    /// Allocate frame targets with picking textures conditionally present.
    /// `with_picking = false` skips the Rg32Uint + matching depth surface.
    /// `color_format` matches the host swap-chain target so the FXAA
    /// scratch can be written to via the same blend states.
    pub fn new_with_picking(
        device: &wgpu::Device,
        width: u32,
        height: u32,
        with_picking: bool,
        color_format: wgpu::TextureFormat,
    ) -> Self {
        let width = width.max(1);
        let height = height.max(1);

        let accum_texture = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("patinae.frame.accum"),
            size: wgpu::Extent3d {
                width,
                height,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: ACCUM_FORMAT,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::TEXTURE_BINDING,
            view_formats: &[],
        });

        let reveal_texture = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("patinae.frame.reveal"),
            size: wgpu::Extent3d {
                width,
                height,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: REVEAL_FORMAT,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::TEXTURE_BINDING,
            view_formats: &[],
        });

        let depth_texture = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("patinae.frame.depth"),
            size: wgpu::Extent3d {
                width,
                height,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: DEPTH_FORMAT,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::TEXTURE_BINDING,
            view_formats: &[],
        });

        let (
            picking_texture,
            picking_depth_texture,
            picking_prev_texture,
            picking_depth_prev_texture,
            picking,
            picking_depth,
            picking_prev,
            picking_depth_prev,
        ) = if with_picking {
            let pick_w = ((width as f32 * PICKING_SCALE) as u32).max(1);
            let pick_h = ((height as f32 * PICKING_SCALE) as u32).max(1);
            // `picking` and `picking_prev` share the same descriptor so we
            // can `copy_texture_to_texture` from current → prev after a
            // full re-record. The current texture also needs
            // `STORAGE_BINDING` because the reprojection write pass binds
            // it as `texture_storage_2d<…, write>`.
            let picking_desc = wgpu::TextureDescriptor {
                label: Some("patinae.frame.picking"),
                size: wgpu::Extent3d {
                    width: pick_w,
                    height: pick_h,
                    depth_or_array_layers: 1,
                },
                mip_level_count: 1,
                sample_count: 1,
                dimension: wgpu::TextureDimension::D2,
                format: PICKING_FORMAT,
                usage: wgpu::TextureUsages::RENDER_ATTACHMENT
                    | wgpu::TextureUsages::TEXTURE_BINDING
                    | wgpu::TextureUsages::STORAGE_BINDING
                    | wgpu::TextureUsages::COPY_SRC
                    | wgpu::TextureUsages::COPY_DST,
                view_formats: &[],
            };
            let picking_prev_desc = wgpu::TextureDescriptor {
                label: Some("patinae.frame.picking_prev"),
                ..picking_desc.clone()
            };
            let picking_texture = device.create_texture(&picking_desc);
            let picking_prev_texture = device.create_texture(&picking_prev_desc);

            // Picking depth needs TEXTURE_BINDING so reproject can sample
            // it (depth_in). And StoreOp::Store on the render pass — set
            // in render_state. Both depth textures must be COPY_SRC/DST
            // so we can snapshot current → prev after full re-record.
            let depth_desc = wgpu::TextureDescriptor {
                label: Some("patinae.frame.picking_depth"),
                size: wgpu::Extent3d {
                    width: pick_w,
                    height: pick_h,
                    depth_or_array_layers: 1,
                },
                mip_level_count: 1,
                sample_count: 1,
                dimension: wgpu::TextureDimension::D2,
                format: DEPTH_FORMAT,
                usage: wgpu::TextureUsages::RENDER_ATTACHMENT
                    | wgpu::TextureUsages::TEXTURE_BINDING
                    | wgpu::TextureUsages::COPY_SRC
                    | wgpu::TextureUsages::COPY_DST,
                view_formats: &[],
            };
            let depth_prev_desc = wgpu::TextureDescriptor {
                label: Some("patinae.frame.picking_depth_prev"),
                ..depth_desc.clone()
            };
            let picking_depth_texture = device.create_texture(&depth_desc);
            let picking_depth_prev_texture = device.create_texture(&depth_prev_desc);

            let picking = picking_texture.create_view(&wgpu::TextureViewDescriptor::default());
            let picking_prev =
                picking_prev_texture.create_view(&wgpu::TextureViewDescriptor::default());
            let picking_depth =
                picking_depth_texture.create_view(&wgpu::TextureViewDescriptor::default());
            let picking_depth_prev =
                picking_depth_prev_texture.create_view(&wgpu::TextureViewDescriptor::default());
            (
                Some(picking_texture),
                Some(picking_depth_texture),
                Some(picking_prev_texture),
                Some(picking_depth_prev_texture),
                Some(picking),
                Some(picking_depth),
                Some(picking_prev),
                Some(picking_depth_prev),
            )
        } else {
            (None, None, None, None, None, None, None, None)
        };

        let accum = accum_texture.create_view(&wgpu::TextureViewDescriptor::default());
        let reveal = reveal_texture.create_view(&wgpu::TextureViewDescriptor::default());
        let depth = depth_texture.create_view(&wgpu::TextureViewDescriptor::default());

        let ssao_desc = wgpu::TextureDescriptor {
            label: Some("patinae.frame.ssao"),
            size: wgpu::Extent3d {
                width,
                height,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: SSAO_FORMAT,
            usage: wgpu::TextureUsages::STORAGE_BINDING | wgpu::TextureUsages::TEXTURE_BINDING,
            view_formats: &[],
        };
        let ssao_texture = device.create_texture(&ssao_desc);
        let ssao_view = ssao_texture.create_view(&wgpu::TextureViewDescriptor::default());
        let ssao_blurred_texture = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("patinae.frame.ssao_blurred"),
            ..ssao_desc
        });
        let ssao_blurred_view =
            ssao_blurred_texture.create_view(&wgpu::TextureViewDescriptor::default());

        // FXAA needs a separate texture to sample from (a render pass
        // cannot bind its color attachment as a sampled texture). When
        // FXAA is enabled, every render pass writes into
        // `color_scratch_view` and the final FXAA pass reads it +
        // outputs to the host target.
        let color_scratch_texture = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("patinae.frame.color_scratch"),
            size: wgpu::Extent3d {
                width,
                height,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: color_format,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::TEXTURE_BINDING,
            view_formats: &[],
        });
        let color_scratch_view =
            color_scratch_texture.create_view(&wgpu::TextureViewDescriptor::default());

        Self {
            width,
            height,
            accum,
            reveal,
            depth,
            picking,
            picking_depth,
            picking_prev,
            picking_depth_prev,
            ssao_texture,
            ssao_view,
            ssao_blurred_texture,
            ssao_blurred_view,
            color_scratch_texture,
            color_scratch_view,
            overlay_id: None,
            overlay_id_depth: None,
            marking_mask: None,
            overlay_color_scratch: None,
            accum_texture,
            reveal_texture,
            depth_texture,
            picking_texture,
            picking_depth_texture,
            picking_prev_texture,
            picking_depth_prev_texture,
            overlay_id_texture: None,
            overlay_id_depth_texture: None,
            marking_mask_texture: None,
            overlay_color_scratch_texture: None,
        }
    }

    /// Half-resolution picking texture size used for both `picking` and
    /// `picking_prev`. Returns `(width, height)`.
    pub fn picking_dims(&self) -> (u32, u32) {
        let pick_w = ((self.width as f32 * PICKING_SCALE) as u32).max(1);
        let pick_h = ((self.height as f32 * PICKING_SCALE) as u32).max(1);
        (pick_w, pick_h)
    }

    pub fn overlay_dims(&self) -> (u32, u32) {
        (self.width.max(1), self.height.max(1))
    }

    /// Lazily allocate the full-resolution id resources used by visual
    /// overlays. Returns `true` when new views were created and bind groups
    /// must be rebuilt.
    pub fn ensure_overlay_id_targets(&mut self, device: &wgpu::Device) -> bool {
        if self.overlay_id.is_some() && self.overlay_id_depth.is_some() {
            return false;
        }

        let (width, height) = self.overlay_dims();
        let size = wgpu::Extent3d {
            width,
            height,
            depth_or_array_layers: 1,
        };

        let overlay_id_texture = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("patinae.frame.overlay_id"),
            size,
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: PICKING_FORMAT,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::TEXTURE_BINDING,
            view_formats: &[],
        });
        let overlay_id_depth_texture = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("patinae.frame.overlay_id_depth"),
            size,
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: DEPTH_FORMAT,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            view_formats: &[],
        });

        self.overlay_id = Some(overlay_id_texture.create_view(&Default::default()));
        self.overlay_id_depth = Some(overlay_id_depth_texture.create_view(&Default::default()));
        self.overlay_id_texture = Some(overlay_id_texture);
        self.overlay_id_depth_texture = Some(overlay_id_depth_texture);
        true
    }

    /// Lazily allocate resources used only by the selection / hover overlay.
    /// Returns `true` when new views were created and bind groups must be
    /// rebuilt.
    pub fn ensure_marking_targets(
        &mut self,
        device: &wgpu::Device,
        color_format: wgpu::TextureFormat,
        needs_overlay_color_scratch: bool,
    ) -> bool {
        if self.marking_mask.is_some()
            && (!needs_overlay_color_scratch || self.overlay_color_scratch.is_some())
        {
            return false;
        }

        let (width, height) = self.overlay_dims();
        let size = wgpu::Extent3d {
            width,
            height,
            depth_or_array_layers: 1,
        };

        let mut allocated = false;
        let mask_desc = wgpu::TextureDescriptor {
            label: Some("patinae.frame.marking_mask"),
            size,
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: MARKING_MASK_FORMAT,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::TEXTURE_BINDING,
            view_formats: &[],
        };

        if self.marking_mask.is_none() {
            let marking_mask_texture = device.create_texture(&mask_desc);
            self.marking_mask = Some(marking_mask_texture.create_view(&Default::default()));
            self.marking_mask_texture = Some(marking_mask_texture);
            allocated = true;
        }
        if needs_overlay_color_scratch && self.overlay_color_scratch.is_none() {
            let overlay_color_scratch_texture = device.create_texture(&wgpu::TextureDescriptor {
                label: Some("patinae.frame.overlay_color_scratch"),
                size,
                mip_level_count: 1,
                sample_count: 1,
                dimension: wgpu::TextureDimension::D2,
                format: color_format,
                usage: wgpu::TextureUsages::RENDER_ATTACHMENT
                    | wgpu::TextureUsages::TEXTURE_BINDING,
                view_formats: &[],
            });
            self.overlay_color_scratch =
                Some(overlay_color_scratch_texture.create_view(&Default::default()));
            self.overlay_color_scratch_texture = Some(overlay_color_scratch_texture);
            allocated = true;
        }
        allocated
    }

    pub fn clear_marking_targets(&mut self) {
        self.marking_mask = None;
        self.overlay_color_scratch = None;
        self.marking_mask_texture = None;
        self.overlay_color_scratch_texture = None;
    }

    pub fn clear_overlay_targets(&mut self) {
        self.overlay_id = None;
        self.overlay_id_depth = None;
        self.overlay_id_texture = None;
        self.overlay_id_depth_texture = None;
        self.clear_marking_targets();
    }
}
