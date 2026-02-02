//! High-level raytracing interface

/// Raytracing parameters
#[derive(Clone, Debug)]
pub struct RaytraceParams {
    /// Output width in pixels
    pub width: u32,
    /// Output height in pixels
    pub height: u32,
    /// Antialiasing level (1 = no AA, 2 = 2x2 supersampling, etc.)
    pub antialias: u32,
    /// View matrix (world to camera)
    pub view_matrix: [[f32; 4]; 4],
    /// Projection matrix
    pub proj_matrix: [[f32; 4]; 4],
    /// Inverse view matrix
    pub view_inv_matrix: [[f32; 4]; 4],
    /// Inverse projection matrix
    pub proj_inv_matrix: [[f32; 4]; 4],
    /// Camera position in world space
    pub camera_pos: [f32; 4],
    /// Raytracing settings
    pub settings: RaytraceSettings,
}

impl RaytraceParams {
    /// Create new raytracing parameters with default settings
    pub fn new(width: u32, height: u32) -> Self {
        Self {
            width,
            height,
            antialias: 1,
            view_matrix: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
            proj_matrix: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
            view_inv_matrix: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
            proj_inv_matrix: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
            camera_pos: [0.0, 0.0, 10.0, 1.0],
            settings: RaytraceSettings::default(),
        }
    }

    /// Set antialiasing level
    pub fn with_antialias(mut self, level: u32) -> Self {
        self.antialias = level.clamp(1, 4);
        self
    }

    /// Set camera matrices
    pub fn with_camera(
        mut self,
        view: [[f32; 4]; 4],
        proj: [[f32; 4]; 4],
        view_inv: [[f32; 4]; 4],
        proj_inv: [[f32; 4]; 4],
        camera_pos: [f32; 4],
    ) -> Self {
        self.view_matrix = view;
        self.proj_matrix = proj;
        self.view_inv_matrix = view_inv;
        self.proj_inv_matrix = proj_inv;
        self.camera_pos = camera_pos;
        self
    }

    /// Set raytracing settings
    pub fn with_settings(mut self, settings: RaytraceSettings) -> Self {
        self.settings = settings;
        self
    }
}

/// Raytracing-specific settings
#[derive(Clone, Debug)]
pub struct RaytraceSettings {
    /// Light direction (normalized)
    pub light_dir: [f32; 4],
    /// Ambient light intensity
    pub ambient: f32,
    /// Direct light intensity
    pub direct: f32,
    /// Reflect light intensity (secondary fill light)
    pub reflect: f32,
    /// Specular intensity
    pub specular: f32,
    /// Specular shininess
    pub shininess: f32,
    /// Background color (RGBA)
    pub bg_color: [f32; 4],
    /// Fog start distance
    pub fog_start: f32,
    /// Fog end distance
    pub fog_end: f32,
    /// Fog density (0 = no fog)
    pub fog_density: f32,
    /// Fog color (RGBA)
    pub fog_color: [f32; 4],
    /// Enable shadows
    pub ray_shadow: bool,
    /// Maximum passes for transparency
    pub ray_max_passes: u32,
    /// Apply fog in raytracing
    pub ray_trace_fog: bool,
    /// Shadows through transparent objects
    pub ray_transparency_shadows: bool,
    /// Interior color index (-1 = use surface color)
    pub ray_interior_color: i32,
    /// Enable interior shadows
    pub ray_interior_shadows: bool,
    /// Direct shading factor
    pub ray_direct_shade: f32,
    /// Orthographic mode (-1 = auto)
    pub ray_orthoscopic: i32,
    /// Ray trace mode (0=normal, 1=normal+outline, 2=outline only, 3=quantized+outline)
    pub ray_trace_mode: i32,
    /// Color for outlines in modes 1, 2, and 3 (RGBA)
    pub ray_trace_color: [f32; 4],
    /// Opaque background (-1=auto, 0=transparent, 1=opaque)
    pub ray_opaque_background: i32,
    /// Edge detection: gradient magnitude difference threshold (PyMOL default: 0.6)
    pub ray_trace_slope_factor: f32,
    /// Edge detection: gradient direction difference threshold (PyMOL default: 0.1)
    pub ray_trace_depth_factor: f32,
    /// Edge detection: gradient discontinuity threshold (PyMOL default: 0.05)
    pub ray_trace_disco_factor: f32,
    /// Edge detection: pixel radius adjustment factor (PyMOL default: 0.12)
    pub ray_trace_gain: f32,
}

impl Default for RaytraceSettings {
    fn default() -> Self {
        Self {
            light_dir: [-0.4, -0.4, -1.0, 0.0],
            ambient: 0.14,
            direct: 0.45,
            reflect: 0.45,
            specular: 0.5,
            shininess: 40.0,
            bg_color: [0.0, 0.0, 0.0, 1.0],
            fog_start: 0.45,
            fog_end: 1.0,
            fog_density: 0.0,
            fog_color: [0.0, 0.0, 0.0, 1.0],
            ray_shadow: true,
            ray_max_passes: 25,
            ray_trace_fog: false,
            ray_transparency_shadows: true,
            ray_interior_color: -1,
            ray_interior_shadows: false,
            ray_direct_shade: 0.0,
            ray_orthoscopic: -1,
            ray_trace_mode: 0,
            ray_trace_color: [0.0, 0.0, 0.0, 1.0], // Black outline by default
            ray_opaque_background: -1, // Auto
            ray_trace_slope_factor: 0.6, // PyMOL default
            ray_trace_depth_factor: 0.1, // PyMOL default
            ray_trace_disco_factor: 0.05, // PyMOL default
            ray_trace_gain: 0.12, // PyMOL default
        }
    }
}
