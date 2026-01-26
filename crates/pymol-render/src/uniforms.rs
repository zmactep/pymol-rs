//! Global uniforms for shaders
//!
//! Contains camera matrices, lighting parameters, and fog settings.

use bytemuck::{Pod, Zeroable};
use lin_alg::f32::{Mat4, Vec3};

/// Global uniforms passed to all shaders
///
/// This structure contains all the global rendering parameters including
/// camera matrices, lighting, and fog settings. It's uploaded to a uniform
/// buffer and bound to all render pipelines.
///
/// ## PyMOL Lighting Model
///
/// PyMOL uses a dual-light system:
/// - **Headlight/Direct**: A light that always comes from the camera direction,
///   ensuring front-facing surfaces are always illuminated regardless of view angle.
/// - **Positional/Reflect**: Directional lights in world space that provide depth
///   and shadow cues.
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct GlobalUniforms {
    /// Combined view-projection matrix
    pub view_proj: [[f32; 4]; 4],
    /// View matrix (camera transform)
    pub view: [[f32; 4]; 4],
    /// Inverse view matrix
    pub view_inv: [[f32; 4]; 4],
    /// Projection matrix
    pub proj: [[f32; 4]; 4],
    /// Camera position in world space (w unused)
    pub camera_pos: [f32; 4],
    /// Primary light direction (normalized, w unused)
    pub light_dir: [f32; 4],

    // === Headlight (camera light) parameters ===
    /// Ambient light intensity (PyMOL default: 0.14)
    pub ambient: f32,
    /// Headlight/camera diffuse intensity (PyMOL default: 0.45)
    pub direct: f32,
    /// Headlight specular intensity (PyMOL default: 0.0)
    pub spec_direct: f32,
    /// Headlight specular exponent (PyMOL default: 55.0)
    pub spec_direct_power: f32,

    // === Positional light parameters ===
    /// Positional light diffuse intensity (PyMOL default: 0.45)
    pub reflect: f32,
    /// Positional light specular intensity (PyMOL default: 1.0)
    pub specular: f32,
    /// Positional light specular exponent (PyMOL default: 55.0)
    pub shininess: f32,
    /// Padding for 16-byte alignment
    pub _pad0: f32,

    // === Fog parameters ===
    /// Fog start distance
    pub fog_start: f32,
    /// Fog end distance (full fog)
    pub fog_end: f32,
    /// Fog density (0 = disabled)
    pub fog_density: f32,
    /// Depth cue factor (0 = disabled)
    pub depth_cue: f32,

    /// Fog color (w unused)
    pub fog_color: [f32; 4],
    /// Background color (w unused)
    pub bg_color: [f32; 4],
    /// Viewport size (width, height) and inverse (1/width, 1/height)
    pub viewport: [f32; 4],
    /// Near and far clip planes (near, far, unused, unused)
    pub clip_planes: [f32; 4],
}

/// Helper to convert Mat4 to uniform-compatible array
fn mat4_to_array(m: Mat4) -> [[f32; 4]; 4] {
    <[[f32; 4]; 4]>::from(m)
}

impl Default for GlobalUniforms {
    fn default() -> Self {
        let identity = mat4_to_array(Mat4::new_identity());
        Self {
            view_proj: identity,
            view: identity,
            view_inv: identity,
            proj: identity,
            camera_pos: [0.0, 0.0, 10.0, 1.0],
            light_dir: [-0.4, -0.4, -1.0, 0.0],
            // Headlight (camera light) - PyMOL defaults
            ambient: 0.14,
            direct: 0.45,
            spec_direct: 0.0,
            spec_direct_power: 55.0,
            // Positional light - PyMOL defaults
            reflect: 0.45,
            specular: 1.0,
            shininess: 55.0,
            _pad0: 0.0,
            // Fog
            fog_start: 0.0,
            fog_end: 1.0,
            fog_density: 0.0,  // Disabled by default - requires proper computation based on clip planes
            depth_cue: 0.0,    // Disabled by default
            fog_color: [0.0, 0.0, 0.0, 1.0],
            bg_color: [0.0, 0.0, 0.0, 1.0],
            viewport: [800.0, 600.0, 1.0 / 800.0, 1.0 / 600.0],
            clip_planes: [0.1, 1000.0, 0.0, 0.0],
        }
    }
}

impl GlobalUniforms {
    /// Create new uniforms with default lighting
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the camera matrices
    pub fn set_camera(&mut self, view: Mat4, proj: Mat4) {
        // Calculate view_proj first (need to clone since Mul takes ownership)
        self.view_proj = mat4_to_array(proj.clone() * view.clone());

        // Compute inverse view matrix for world-space calculations
        if let Some(inv) = view.inverse() {
            // Extract camera position from inverse view matrix
            self.camera_pos = [inv.data[12], inv.data[13], inv.data[14], 1.0];
            self.view_inv = mat4_to_array(inv);
        }

        self.view = mat4_to_array(view);
        self.proj = mat4_to_array(proj);
    }

    /// Set the primary light direction (will be normalized)
    pub fn set_light_direction(&mut self, dir: Vec3) {
        let normalized = dir.to_normalized();
        self.light_dir = [normalized.x, normalized.y, normalized.z, 0.0];
    }

    /// Set lighting parameters (PyMOL dual-light model)
    ///
    /// # Arguments
    /// * `ambient` - Ambient light intensity
    /// * `direct` - Headlight (camera) diffuse intensity
    /// * `reflect` - Positional light diffuse intensity
    /// * `specular` - Positional light specular intensity
    /// * `shininess` - Positional light specular exponent
    /// * `spec_direct` - Headlight specular intensity
    /// * `spec_direct_power` - Headlight specular exponent
    pub fn set_lighting(
        &mut self,
        ambient: f32,
        direct: f32,
        reflect: f32,
        specular: f32,
        shininess: f32,
        spec_direct: f32,
        spec_direct_power: f32,
    ) {
        self.ambient = ambient;
        self.direct = direct;
        self.reflect = reflect;
        self.specular = specular;
        self.shininess = shininess;
        self.spec_direct = spec_direct;
        self.spec_direct_power = spec_direct_power;
    }

    /// Set fog parameters
    pub fn set_fog(&mut self, start: f32, end: f32, density: f32, color: [f32; 3]) {
        self.fog_start = start;
        self.fog_end = end;
        self.fog_density = density;
        self.fog_color = [color[0], color[1], color[2], 1.0];
    }

    /// Set viewport size
    pub fn set_viewport(&mut self, width: f32, height: f32) {
        self.viewport = [width, height, 1.0 / width, 1.0 / height];
    }

    /// Set clip planes
    pub fn set_clip_planes(&mut self, near: f32, far: f32) {
        self.clip_planes = [near, far, 0.0, 0.0];
    }

    /// Set background color
    pub fn set_background(&mut self, color: [f32; 3]) {
        self.bg_color = [color[0], color[1], color[2], 1.0];
    }

    /// Set depth cue factor
    pub fn set_depth_cue(&mut self, factor: f32) {
        self.depth_cue = factor;
    }
}

/// Builder for creating global uniforms from settings
#[allow(dead_code)]
pub struct UniformsBuilder {
    uniforms: GlobalUniforms,
}

#[allow(dead_code)]
impl UniformsBuilder {
    /// Create a new builder with default values
    pub fn new() -> Self {
        Self {
            uniforms: GlobalUniforms::default(),
        }
    }

    /// Set camera matrices
    pub fn camera(mut self, view: Mat4, proj: Mat4) -> Self {
        self.uniforms.set_camera(view, proj);
        self
    }

    /// Set light direction
    pub fn light_direction(mut self, dir: Vec3) -> Self {
        self.uniforms.set_light_direction(dir);
        self
    }

    /// Set lighting parameters (PyMOL dual-light model)
    pub fn lighting(
        mut self,
        ambient: f32,
        direct: f32,
        reflect: f32,
        specular: f32,
        shininess: f32,
        spec_direct: f32,
        spec_direct_power: f32,
    ) -> Self {
        self.uniforms.set_lighting(ambient, direct, reflect, specular, shininess, spec_direct, spec_direct_power);
        self
    }

    /// Set fog parameters
    pub fn fog(mut self, start: f32, end: f32, density: f32, color: [f32; 3]) -> Self {
        self.uniforms.set_fog(start, end, density, color);
        self
    }

    /// Set viewport size
    pub fn viewport(mut self, width: f32, height: f32) -> Self {
        self.uniforms.set_viewport(width, height);
        self
    }

    /// Build the uniforms
    pub fn build(self) -> GlobalUniforms {
        self.uniforms
    }
}

impl Default for UniformsBuilder {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_uniforms_size() {
        // Ensure the struct is properly aligned for GPU
        assert_eq!(std::mem::size_of::<GlobalUniforms>() % 16, 0);
    }

    #[test]
    fn test_uniforms_default() {
        let uniforms = GlobalUniforms::default();
        assert!(uniforms.ambient > 0.0);
        assert!(uniforms.shininess > 0.0);
    }

    #[test]
    fn test_builder() {
        let uniforms = UniformsBuilder::new()
            .lighting(0.2, 0.5, 0.4, 0.8, 32.0, 0.1, 40.0)
            .fog(0.5, 1.0, 0.8, [0.1, 0.1, 0.1])
            .build();

        assert_eq!(uniforms.ambient, 0.2);
        assert_eq!(uniforms.direct, 0.5);
        assert_eq!(uniforms.reflect, 0.4);
        assert_eq!(uniforms.specular, 0.8);
        assert_eq!(uniforms.shininess, 32.0);
        assert_eq!(uniforms.spec_direct, 0.1);
        assert_eq!(uniforms.spec_direct_power, 40.0);
    }

    #[test]
    fn test_pymol_defaults() {
        let uniforms = GlobalUniforms::default();
        // PyMOL default lighting values
        assert_eq!(uniforms.ambient, 0.14);
        assert_eq!(uniforms.direct, 0.45);
        assert_eq!(uniforms.reflect, 0.45);
        assert_eq!(uniforms.specular, 1.0);
        assert_eq!(uniforms.shininess, 55.0);
        assert_eq!(uniforms.spec_direct, 0.0);
        assert_eq!(uniforms.spec_direct_power, 55.0);
    }
}
