//! FFI Contract
//!
//! Defines the C-compatible entry point that plugins export and the host
//! discovers via `libloading`.

/// ABI version — bumped when the `PluginDeclaration` layout changes.
pub const ABI_VERSION: u32 = 1;

/// SDK version — must match between host and plugin at load time.
pub const SDK_VERSION: &str = env!("CARGO_PKG_VERSION");

/// C-compatible plugin declaration exported by every plugin.
///
/// The host loads the shared library, looks up the `PYMOL_PLUGIN_DECLARATION`
/// symbol, checks ABI and SDK versions, then calls `register` with a
/// `PluginRegistrar` to collect commands, components, and handlers.
#[repr(C)]
pub struct PluginDeclaration {
    /// Must equal [`ABI_VERSION`].
    pub abi_version: u32,
    /// Pointer to the SDK version string.
    pub sdk_version_ptr: *const u8,
    /// Length of the SDK version string (no null terminator).
    pub sdk_version_len: usize,
    /// Registration function. Called once at load time.
    pub register: unsafe extern "C" fn(registrar: *mut crate::registrar::PluginRegistrar),
}

// SAFETY: PluginDeclaration is a static read-only descriptor.
// The register function pointer is called on the main thread.
unsafe impl Send for PluginDeclaration {}
unsafe impl Sync for PluginDeclaration {}
