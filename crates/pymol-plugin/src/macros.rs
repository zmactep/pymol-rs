//! Plugin Declaration Macro
//!
//! The `pymol_plugin!` macro generates the FFI entry point, metadata,
//! and registration logic in a single declarative invocation.

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
///     version: "0.1.0",
///     description: "What this plugin does",
///     commands: [Cmd1, Cmd2],
///     // Optional: register GUI components with panel configs
///     components: [(MyPanel::new(), PanelConfig::right(250.0))],
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
                        let _ = $crate::log::set_logger(logger);
                        $crate::log::set_max_level(level);
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
                            $( let $reg = __registrar; $body )?
                        }));
                    }
                    __pymol_register
                },
            };
    };
}
