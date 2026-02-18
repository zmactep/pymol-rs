//! Classic shading pipeline — no shadow passes.

use crate::RenderContext;
use pymol_settings::GlobalSettings;

use super::ShadingPipeline;

pub struct ClassicPipeline;

impl ShadingPipeline for ClassicPipeline {
    fn prepare(&mut self, _context: &mut RenderContext, _settings: &GlobalSettings) -> bool {
        false // No shadow passes needed.
    }

    fn deactivate(&mut self, _context: &mut RenderContext) {}

    fn needs_shadow_update(&self) -> bool { false }

    fn invalidate_shadows(&mut self) {}
}
