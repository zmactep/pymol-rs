//! Renderer memory settings.

use crate::define_settings_group;
use crate::enums::RenderMemoryProfileSetting;

const MAX_RENDER_MEMORY_BUDGET_MIB: f32 = 1_048_576.0;

define_settings_group! {
    /// Renderer memory controls.
    group_global RendererSettings {
        memory_profile: RenderMemoryProfileSetting = RenderMemoryProfileSetting::Auto,
            name = "render_memory_profile",
            hints = RenderMemoryProfileSetting,
            side_effects = [ViewportUpdate];
        memory_budget_mib: i32 = 0,
            name = "render_memory_budget",
            min = 0.0, max = MAX_RENDER_MEMORY_BUDGET_MIB,
            side_effects = [ViewportUpdate];
    }
}
