//! Render texture creation for the raytracing pipeline.

/// Holds a texture and its default view.
pub(crate) struct TextureWithView {
    pub texture: wgpu::Texture,
    pub view: wgpu::TextureView,
}

impl TextureWithView {
    /// Create a texture and its default view.
    fn new(
        device: &wgpu::Device,
        label: &'static str,
        width: u32,
        height: u32,
        format: wgpu::TextureFormat,
        usage: wgpu::TextureUsages,
    ) -> Self {
        let texture = device.create_texture(&wgpu::TextureDescriptor {
            label: Some(label),
            size: texture_size(width, height),
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format,
            usage,
            view_formats: &[],
        });
        let view = texture.create_view(&wgpu::TextureViewDescriptor::default());
        Self { texture, view }
    }
}

/// Holds the color, depth, and normal textures used by the raytracing passes.
pub(crate) struct RenderTextures {
    pub color_texture: wgpu::Texture,
    pub color_view: wgpu::TextureView,
    pub depth_view: wgpu::TextureView,
    pub normal_view: wgpu::TextureView,
}

impl RenderTextures {
    /// Create the three render textures (color, depth, normal) at the given resolution.
    pub fn new(device: &wgpu::Device, width: u32, height: u32) -> Self {
        let color_texture = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("Raytrace Color Output"),
            size: texture_size(width, height),
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: wgpu::TextureFormat::Rgba8Unorm,
            usage: wgpu::TextureUsages::STORAGE_BINDING
                | wgpu::TextureUsages::TEXTURE_BINDING
                | wgpu::TextureUsages::COPY_SRC,
            view_formats: &[],
        });
        let color_view = color_texture.create_view(&wgpu::TextureViewDescriptor::default());

        let depth_texture = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("Raytrace Depth Output"),
            size: texture_size(width, height),
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: wgpu::TextureFormat::R32Float,
            usage: wgpu::TextureUsages::STORAGE_BINDING | wgpu::TextureUsages::TEXTURE_BINDING,
            view_formats: &[],
        });
        let depth_view = depth_texture.create_view(&wgpu::TextureViewDescriptor::default());

        let normal_texture = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("Raytrace Normal Output"),
            size: texture_size(width, height),
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: wgpu::TextureFormat::Rgba16Float,
            usage: wgpu::TextureUsages::STORAGE_BINDING | wgpu::TextureUsages::TEXTURE_BINDING,
            view_formats: &[],
        });
        let normal_view = normal_texture.create_view(&wgpu::TextureViewDescriptor::default());

        Self {
            color_texture,
            color_view,
            depth_view,
            normal_view,
        }
    }
}

/// Create the intermediate edge-detection texture.
pub(crate) fn edge_texture(device: &wgpu::Device, width: u32, height: u32) -> TextureWithView {
    TextureWithView::new(
        device,
        "Edge Detection Output",
        width,
        height,
        wgpu::TextureFormat::R32Float,
        wgpu::TextureUsages::STORAGE_BINDING | wgpu::TextureUsages::TEXTURE_BINDING,
    )
}

/// Create the ray trace mode composite texture.
pub(crate) fn composite_texture(device: &wgpu::Device, width: u32, height: u32) -> TextureWithView {
    TextureWithView::new(
        device,
        "Composite Output",
        width,
        height,
        wgpu::TextureFormat::Rgba8Unorm,
        wgpu::TextureUsages::STORAGE_BINDING
            | wgpu::TextureUsages::TEXTURE_BINDING
            | wgpu::TextureUsages::COPY_SRC,
    )
}

/// Create the final antialiasing downsample texture.
pub(crate) fn downsample_texture(
    device: &wgpu::Device,
    width: u32,
    height: u32,
) -> TextureWithView {
    TextureWithView::new(
        device,
        "Standalone Downsample Output",
        width,
        height,
        wgpu::TextureFormat::Rgba8Unorm,
        wgpu::TextureUsages::STORAGE_BINDING | wgpu::TextureUsages::COPY_SRC,
    )
}

fn texture_size(width: u32, height: u32) -> wgpu::Extent3d {
    wgpu::Extent3d {
        width,
        height,
        depth_or_array_layers: 1,
    }
}
