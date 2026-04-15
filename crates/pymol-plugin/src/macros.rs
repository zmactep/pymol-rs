//! Plugin Declaration Macro
//!
//! The `pymol_plugin!` macro generates the FFI entry point, metadata,
//! and registration logic in a single declarative invocation.

/// Define a typed plugin settings struct.
///
/// Generates a struct with typed fields, a `Default` impl, and helper methods
/// for creating descriptors and syncing with a shared store.
///
/// # Syntax
///
/// ```rust,ignore
/// define_plugin_settings! {
///     MySettings {
///         threshold: f32 = 0.5, name = "my_threshold",
///             min = 0.0, max = 1.0,
///             side_effects = [SceneInvalidate];
///         mode: i32 = 0, name = "my_mode";
///         enabled: bool = true, name = "my_enabled";
///     }
/// }
/// ```
///
/// # Generated API
///
/// - `MySettings` struct with `pub` fields + `Default` impl
/// - `MySettings::descriptors()` → `Vec<DynamicSettingDescriptor>`
/// - `MySettings::init_store(&self)` → `SharedSettingStore` (populated with defaults)
/// - `MySettings::sync_from_store(&mut self, store: &DynamicSettingStore)` — reads store into typed fields
///
/// Settings are global-only by default. To make a setting object-overridable,
/// modify the descriptor after calling `descriptors()`:
///
/// ```rust,ignore
/// let mut descs = MySettings::descriptors();
/// descs[0].object_overridable = true;
/// ```
#[macro_export]
macro_rules! define_plugin_settings {
    (
        $(#[$meta:meta])*
        $Name:ident {
            $(
                $(#[$field_meta:meta])*
                $field:ident : $ty:tt = $default:expr,
                    name = $sname:expr
                    $(, min = $min:expr, max = $max:expr)?
                    $(, side_effects = [$($se:ident),* $(,)?])?
                    ;
            )*
        }
    ) => {
        $(#[$meta])*
        #[derive(Debug, Clone)]
        pub struct $Name {
            $(
                $(#[$field_meta])*
                pub $field: $ty,
            )*
        }

        impl Default for $Name {
            fn default() -> Self {
                Self {
                    $( $field: $default, )*
                }
            }
        }

        impl $Name {
            /// Build descriptors for all settings in this group.
            pub fn descriptors() -> Vec<$crate::registrar::DynamicSettingDescriptor> {
                vec![
                    $(
                        $crate::registrar::DynamicSettingDescriptor {
                            name: $sname.to_string(),
                            setting_type: $crate::define_plugin_settings!(@setting_type $ty),
                            default: $crate::define_plugin_settings!(@wrap $ty, $default),
                            min: $crate::define_plugin_settings!(@opt_f32 $($min)?),
                            max: $crate::define_plugin_settings!(@opt_f32 $($max)?),
                            value_hints: vec![],
                            side_effects: vec![$($($crate::define_plugin_settings!(@side_effect $se)),*)?],
                            object_overridable: false,
                        },
                    )*
                ]
            }

            /// Create a shared store populated with default values.
            pub fn init_store(&self) -> $crate::registrar::SharedSettingStore {
                use std::sync::{Arc, RwLock};
                let mut store = $crate::registrar::DynamicSettingStore::new();
                $(
                    store.set($sname, $crate::define_plugin_settings!(@wrap $ty, self.$field));
                )*
                Arc::new(RwLock::new(store))
            }

            /// Read current values from the store into typed fields.
            ///
            /// Call this at the start of `poll()` to pick up changes
            /// made by the `set` command.
            pub fn sync_from_store(&mut self, store: &$crate::registrar::DynamicSettingStore) {
                $(
                    if let Some(val) = store.get($sname) {
                        $crate::define_plugin_settings!(@unwrap self.$field, $ty, val);
                    }
                )*
            }
        }
    };

    // =========================================================================
    // Internal helpers
    // =========================================================================

    // --- Setting type mapping ---
    (@setting_type bool) => { $crate::__private::SettingType::Bool };
    (@setting_type i32) => { $crate::__private::SettingType::Int };
    (@setting_type f32) => { $crate::__private::SettingType::Float };

    // --- Wrap value into SettingValue ---
    (@wrap bool, $val:expr) => { $crate::__private::SettingValue::Bool($val) };
    (@wrap i32, $val:expr) => { $crate::__private::SettingValue::Int($val) };
    (@wrap f32, $val:expr) => { $crate::__private::SettingValue::Float($val) };

    // --- Unwrap SettingValue into typed field ---
    (@unwrap $field:expr, bool, $val:expr) => {
        if let Some(v) = $val.as_bool() { $field = v; }
    };
    (@unwrap $field:expr, i32, $val:expr) => {
        if let Some(v) = $val.as_int() { $field = v; }
    };
    (@unwrap $field:expr, f32, $val:expr) => {
        if let Some(v) = $val.as_float() { $field = v; }
    };

    // --- Optional f32 ---
    (@opt_f32) => { None };
    (@opt_f32 $val:expr) => { Some($val as f32) };

    // --- Side effect mapping ---
    (@side_effect RepresentationRebuild) => { $crate::__private::SideEffectCategory::RepresentationRebuild };
    (@side_effect ColorRebuild) => { $crate::__private::SideEffectCategory::ColorRebuild };
    (@side_effect SceneInvalidate) => { $crate::__private::SideEffectCategory::SceneInvalidate };
    (@side_effect SceneChanged) => { $crate::__private::SideEffectCategory::SceneChanged };
    (@side_effect ShaderReload) => { $crate::__private::SideEffectCategory::ShaderReload };
    (@side_effect ShaderComputeLighting) => { $crate::__private::SideEffectCategory::ShaderComputeLighting };
    (@side_effect ViewportUpdate) => { $crate::__private::SideEffectCategory::ViewportUpdate };
    (@side_effect OrthoDirty) => { $crate::__private::SideEffectCategory::OrthoDirty };
    (@side_effect FullRebuild) => { $crate::__private::SideEffectCategory::FullRebuild };
}

/// Declare a PyMOL-RS plugin.
///
/// Generates a `#[no_mangle] static PYMOL_PLUGIN_DECLARATION` with the correct
/// ABI version, SDK version, and a registration function that populates a
/// [`PluginRegistrar`](crate::registrar::PluginRegistrar) with commands,
/// components, and optional custom logic.
///
/// # Syntax
///
/// ```rust,ignore
/// pymol_plugin! {
///     name: "my-plugin",
///     // version is optional — defaults to env!("CARGO_PKG_VERSION")
///     description: "What this plugin does",
///     commands: [Cmd1, Cmd2],
///     // Optional: register GUI components with panel configs
///     components: [(MyPanel::new(), PanelConfig::right(250.0))],
///     // Optional: register settings structs (from define_plugin_settings!)
///     settings: [MySettings],
///     // Optional: custom registration logic
///     register: |reg| {
///         reg.set_message_handler(MyHandler);
///     },
/// }
/// ```
///
/// # Panic Safety
///
/// The generated registration function wraps all user code in
/// `std::panic::catch_unwind` to prevent panics from crossing the FFI boundary.
#[macro_export]
macro_rules! pymol_plugin {
    (
        name: $name:expr,
        version: $version:expr,
        description: $desc:expr,
        commands: [$($cmd:expr),* $(,)?]
        $(, components: [$( ($comp:expr, $config:expr) ),* $(,)?] )?
        $(, settings: [$( $settings_ty:ty ),* $(,)?] )?
        $(, hotkeys: [$( ($hk_key:expr, $hk_action:expr) ),* $(,)?] )?
        $(, register: |$reg:ident| $body:block )?
        $(,)?
    ) => {
        #[no_mangle]
        pub static PYMOL_PLUGIN_DECLARATION: $crate::ffi::PluginDeclaration =
            $crate::ffi::PluginDeclaration {
                abi_version: $crate::ffi::ABI_VERSION,
                sdk_version_ptr: $crate::ffi::SDK_VERSION.as_ptr(),
                sdk_version_len: $crate::ffi::SDK_VERSION.len(),
                init: {
                    #[allow(improper_ctypes_definitions)]
                    unsafe extern "C" fn __pymol_init(
                        logger: &'static dyn $crate::log::Log,
                        level: $crate::log::LevelFilter,
                    ) {
                        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                            let _ = $crate::log::set_logger(logger);
                            $crate::log::set_max_level(level);
                        }));
                    }
                    __pymol_init
                },
                register: {
                    unsafe extern "C" fn __pymol_register(
                        __registrar: *mut $crate::registrar::PluginRegistrar,
                    ) {
                        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                            let __registrar = unsafe { &mut *__registrar };
                            __registrar.set_metadata($crate::registrar::PluginMetadata {
                                name: $name,
                                version: $version,
                                description: $desc,
                            });
                            $( __registrar.register_command($cmd); )*
                            $( $( __registrar.register_component($comp, $config); )* )?
                            $( $( __registrar.register_hotkey($hk_key, $hk_action); )* )?
                            $( $(
                                {
                                    let __s = <$settings_ty>::default();
                                    let __store = __s.init_store();
                                    __registrar.register_settings(<$settings_ty>::descriptors(), __store);
                                }
                            )* )?
                            $( let $reg = __registrar; $body )?
                        }));
                    }
                    __pymol_register
                },
            };
    };

    // Arm without `version:` — defaults to the crate version from Cargo.toml.
    (
        name: $name:expr,
        description: $desc:expr,
        commands: [$($cmd:expr),* $(,)?]
        $(, components: [$( ($comp:expr, $config:expr) ),* $(,)?] )?
        $(, settings: [$( $settings_ty:ty ),* $(,)?] )?
        $(, hotkeys: [$( ($hk_key:expr, $hk_action:expr) ),* $(,)?] )?
        $(, register: |$reg:ident| $body:block )?
        $(,)?
    ) => {
        pymol_plugin! {
            name: $name,
            version: env!("CARGO_PKG_VERSION"),
            description: $desc,
            commands: [$($cmd),*]
            $(, components: [$( ($comp, $config) ),*] )?
            $(, settings: [$( $settings_ty ),*] )?
            $(, hotkeys: [$( ($hk_key, $hk_action) ),*] )?
            $(, register: |$reg| $body )?
        }
    };
}
