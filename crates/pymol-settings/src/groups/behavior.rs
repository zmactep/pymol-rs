//! Behavioral / algorithmic settings (global-only).

use crate::define_settings_group;

define_settings_group! {
    /// Behavioral and algorithmic defaults.
    group_global BehaviorSettings {
        ignore_case: bool = true,
            name = "ignore_case";
        ignore_case_chain: bool = false,
            name = "ignore_case_chain";
        auto_dss: bool = true,
            name = "auto_dss";
        dss_algorithm: crate::DssAlgorithm = crate::DssAlgorithm::PyMol,
            name = "dss_algorithm",
            hints = crate::DssAlgorithm;
        bonding_vdw_cutoff: f32 = 0.2,
            name = "bonding_vdw_cutoff",
            min = 0.0, max = 1.0;
    }
}
