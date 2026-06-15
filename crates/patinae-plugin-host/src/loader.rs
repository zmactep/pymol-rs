use std::collections::{hash_map::DefaultHasher, HashMap, VecDeque};
use std::ffi::c_void;
use std::fs::OpenOptions;
use std::hash::{Hash, Hasher};
use std::io::{Read, Write};
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex, RwLock};

use libloading::Library;
use wgpu::util::DeviceExt;

use patinae_cmd::{
    ArgHint, CmdError, CmdResult, Command, CommandAction, CommandContext, CommandExecutor,
    CommandRuntimeRequirements, DynamicSettingRegistry, FormatHandler, MessageKind, ParsedCommand,
    PluginReaderFn, PluginWriterFn, ScriptHandler,
};
use patinae_framework::component::SharedContext;
use patinae_framework::message::{AppMessage, MessageBus};
use patinae_framework::plugin_ui::{
    PanelControl, PanelDescriptor, PanelEvent, PanelPlacement, PanelRuntimeRequirements,
    PanelSnapshot, PluginPanel,
};
use patinae_mol::{Atom, AtomIndex, Element, ObjectMolecule, SecondaryStructure};
use patinae_plugin::ffi::{
    AbiBytesSinkFn, AbiCommandDescriptor, AbiCommandVTable, AbiFormatHandlerDescriptor,
    AbiFormatHandlerVTable, AbiHotkeyDescriptor, AbiHotkeyVTable, AbiMessageHandlerDescriptor,
    AbiMessageHandlerVTable, AbiPanelDescriptor, AbiPanelVTable, AbiReadError,
    AbiScriptHandlerDescriptor, AbiScriptHandlerVTable, AbiSettingDescriptor, AbiSettingValue,
    AbiSettingValueHint, AbiStatus, AbiStr, AbiStrSlice, AbiU8Slice, HostCallbacks,
    HostCommandRuntimeCallbacks, HostCommandRuntimeHandle, HostRegistrarHandle,
    PluginCommandHandle, PluginDeclaration, PluginFormatHandlerHandle, PluginHotkeyHandle,
    PluginMessageHandlerHandle, PluginPanelHandle, PluginScriptHandlerHandle,
    ABI_STATUS_HOST_ERROR, ABI_STATUS_INVALID, ABI_STATUS_PANIC, ABI_STATUS_UNSUPPORTED,
    ABI_VERSION, ARG_HINT_COLOR, ARG_HINT_COMMAND, ARG_HINT_LABEL_PROPERTY,
    ARG_HINT_NAMED_SELECTION, ARG_HINT_OBJECT, ARG_HINT_PATH, ARG_HINT_REPRESENTATION,
    ARG_HINT_SELECTION, ARG_HINT_SETTING, ARG_HINT_SETTING_VALUE, CAPABILITY_REGISTRATION,
    HOST_CALLBACKS_VERSION, HOST_COMMAND_RUNTIME_CALLBACKS_VERSION, KNOWN_CAPABILITIES,
    KNOWN_COMMAND_RUNTIME_REQUIREMENTS, KNOWN_PANEL_RUNTIME_REQUIREMENTS, MAX_ABI_SLICE_LEN,
    MAX_ABI_STRING_LEN, PANEL_PLACEMENT_BOTTOM, PANEL_PLACEMENT_RIGHT, SDK_VERSION,
    SETTING_TYPE_BLANK, SETTING_TYPE_BOOL, SETTING_TYPE_COLOR, SETTING_TYPE_FLOAT,
    SETTING_TYPE_FLOAT3, SETTING_TYPE_INT, SETTING_TYPE_STRING, SETTING_VALUE_BOOL,
    SETTING_VALUE_COLOR, SETTING_VALUE_FLOAT, SETTING_VALUE_FLOAT3, SETTING_VALUE_INT,
    SETTING_VALUE_NONE, SETTING_VALUE_STRING, SIDE_EFFECT_COLOR_REBUILD, SIDE_EFFECT_FULL_REBUILD,
    SIDE_EFFECT_ORTHO_DIRTY, SIDE_EFFECT_REPRESENTATION_REBUILD, SIDE_EFFECT_SCENE_CHANGED,
    SIDE_EFFECT_SCENE_INVALIDATE, SIDE_EFFECT_SEQ_CHANGED, SIDE_EFFECT_SHADER_COMPUTE_LIGHTING,
    SIDE_EFFECT_SHADER_RELOAD, SIDE_EFFECT_STEREO_UPDATE, SIDE_EFFECT_SURFACE_TRANSPARENCY,
    SIDE_EFFECT_VIEWPORT_UPDATE,
};
use patinae_plugin::registrar::{MessageHandler, PluginKeyAction, PluginMetadata, PollContext};
use patinae_plugin::wire::{
    self, WireCommandInput, WireCommandOutput, WireCommandRuntimeRequest,
    WireCommandRuntimeResponse, WireCommandRuntimeValue, WireFormatReadInput, WireFormatReadOutput,
    WireFormatWriteInput, WireFormatWriteOutput, WireHotkeyAction, WireMessageInput,
    WireMessageOutput, WirePanelEventInput, WirePanelEventOutput, WirePanelSnapshotOutput,
    WirePollInput, WirePollOutput, WireScriptInput, WireScriptOutput, WireSharedInput,
    WireTraceGeometryOpened, WireViewerAction, RUNTIME_WIRE_VERSION,
};
use patinae_scene::{
    parse_key_string, GpuBatchCommand, GpuBatchResult, GpuBindGroupDescriptor, GpuBindGroupEntry,
    GpuBindGroupLayoutDescriptor, GpuBindingResource, GpuBindingType, GpuBufferBindingType,
    GpuBufferDescriptor, GpuBufferUsage, GpuCacheStats, GpuCacheStatus, GpuCachedHandle,
    GpuComputePipelineDescriptor, GpuDeviceLimits, GpuHandle, GpuHandleKind,
    GpuPipelineLayoutDescriptor, GpuShaderModuleDescriptor, GpuShaderStages, GpuSubmitBatch,
    KeyBindings, RenderArtifactBufferDescriptor, RenderArtifactBufferRole,
    RenderArtifactPrimitiveTopology, RenderArtifactRepDescriptor, RenderArtifactRepKind,
    RenderArtifactSnapshotDescriptor, Session, ViewerLike,
};
use patinae_settings::{
    DynamicSettingDescriptor, DynamicSettingStore, SettingType, SettingValue, SharedSettingStore,
    SideEffectCategory,
};
use serde::{de::DeserializeOwned, Serialize};

use crate::host::PluginHost;
use crate::panic::panic_payload_to_string;
use crate::paths::{is_plugin_library_path, PluginDiscovery};
use crate::plugin::{LibraryHandle, LoadedPanel, LoadedPlugin};

/// Vertices per spooled mesh chunk.
///
/// `DisplayedMesh` is a triangle list. The value is divisible by three so
/// chunks never split a triangle, and it keeps encoded chunks far below the
/// 64 MiB ABI payload ceiling.
const DISPLAYED_GEOMETRY_SPOOL_MESH_VERTICES: usize = 24_576;

/// Non-mesh primitives per spooled chunk.
///
/// Analytic primitives encode much smaller than mesh vertices; this still
/// leaves wide headroom below the ABI payload ceiling.
const DISPLAYED_GEOMETRY_SPOOL_PRIMITIVES: usize = 32_768;

static DISPLAYED_GEOMETRY_SPOOL_COUNTER: AtomicU64 = AtomicU64::new(0);

type SharedGpuPluginCache = Arc<Mutex<GpuPluginCache>>;

const GPU_CACHE_LAYOUT_VERSION: u32 = 1;

fn rt_profile_enabled() -> bool {
    std::env::var_os("PATINAE_RT_PROFILE").is_some()
}

static HOST_CALLBACKS: HostCallbacks = HostCallbacks {
    table_version: HOST_CALLBACKS_VERSION,
    register_metadata: Some(host_register_metadata),
    register_command: Some(host_register_command),
    register_panel: Some(host_register_panel),
    register_setting: Some(host_register_setting),
    register_message_handler: Some(host_register_message_handler),
    register_script_handler: Some(host_register_script_handler),
    register_format_handler: Some(host_register_format_handler),
    register_hotkey: Some(host_register_hotkey),
    report_unsupported: Some(host_report_unsupported),
};

static HOST_COMMAND_RUNTIME_CALLBACKS: HostCommandRuntimeCallbacks = HostCommandRuntimeCallbacks {
    table_version: HOST_COMMAND_RUNTIME_CALLBACKS_VERSION,
    request: Some(host_command_runtime_request),
};

impl PluginHost {
    pub fn load_standard_dirs(&mut self, executor: &mut CommandExecutor) {
        self.load_discovered_dirs(&PluginDiscovery::from_process_env(), executor);
    }

    pub fn load_discovered_dirs(
        &mut self,
        discovery: &PluginDiscovery,
        executor: &mut CommandExecutor,
    ) {
        for dir in discovery.standard_plugin_dirs() {
            self.load_dir(&dir, executor);
        }
    }

    pub fn load_dir(&mut self, dir: &Path, executor: &mut CommandExecutor) {
        self.add_plugin_dir(dir);

        let entries = match std::fs::read_dir(dir) {
            Ok(entries) => entries,
            Err(e) => {
                log::debug!("Plugin directory {:?}: {}", dir, e);
                return;
            }
        };

        for entry in entries.flatten() {
            let path = entry.path();
            if !is_plugin_library_path(&path) {
                continue;
            }
            match self.load_library(&path, executor) {
                Ok(name) => log::info!("Loaded plugin: {}", name),
                Err(e) => log::warn!("Failed to load plugin {:?}: {}", path, e),
            }
        }
    }

    pub fn load_library(
        &mut self,
        path: &Path,
        executor: &mut CommandExecutor,
    ) -> Result<String, String> {
        #[cfg(target_os = "windows")]
        crate::paths::apply_deps_search_paths(path);

        // SAFETY: Loading a plugin is the explicit native-code extension point.
        // The resulting handle is kept alive in LoadedPlugin while needed.
        let library =
            unsafe { Library::new(path).map_err(|e| format!("Failed to load library: {e}"))? };

        let declaration = load_declaration(&library)?;
        validate_declaration(&declaration)?;
        initialize_plugin(&declaration)?;

        let mut registration = RegistrationSink::new();
        register_plugin(&declaration, &mut registration)?;

        finish_registration(
            self,
            executor,
            LibraryHandle::Dynamic(library),
            registration,
        )
    }
}

#[cfg(test)]
pub(crate) fn load_declaration_for_test(
    host: &mut PluginHost,
    executor: &mut CommandExecutor,
    declaration: PluginDeclaration,
) -> Result<String, String> {
    validate_declaration(&declaration)?;
    initialize_plugin(&declaration)?;

    let mut registration = RegistrationSink::new();
    register_plugin(&declaration, &mut registration)?;

    finish_registration(host, executor, LibraryHandle::Static, registration)
}

fn finish_registration(
    host: &mut PluginHost,
    executor: &mut CommandExecutor,
    library_handle: LibraryHandle,
    mut registration: RegistrationSink,
) -> Result<String, String> {
    let metadata = registration
        .metadata
        .take()
        .ok_or("Plugin did not set metadata")?;
    let plugin_name = format!("{} v{}", metadata.name, metadata.version);
    let plugin_id = format!("{}@{}", metadata.name, metadata.version);
    let gpu_cache = Arc::new(Mutex::new(GpuPluginCache::new(plugin_id)));
    let library = Arc::new(library_handle);

    install_registration_assets(executor, &mut registration, library.clone(), gpu_cache);

    let mut panels = Vec::new();
    for panel in registration.panels {
        if host.panel_exists(&panel.descriptor.id)
            || panels
                .iter()
                .any(|p: &LoadedPanel| p.descriptor.id == panel.descriptor.id)
        {
            log::warn!(
                "Skipping duplicate plugin panel id '{}'",
                panel.descriptor.id
            );
            continue;
        }
        let descriptor = panel.descriptor.clone();
        panels.push(LoadedPanel::new(
            descriptor,
            Box::new(AbiPanelProxy::new(panel, library.clone())),
        ));
    }

    let mut hotkeys = KeyBindings::new();
    for hotkey in registration.hotkeys {
        let key_str = hotkey.key.clone();
        let action = hotkey.into_action(library.clone());
        match parse_key_string(&key_str) {
            Ok(key) => hotkeys.bind(key, action),
            Err(error) => log::warn!("Invalid plugin hotkey string '{}': {}", key_str, error),
        }
    }

    let message_handler = registration.message_handler.map(|handler| {
        Box::new(AbiMessageHandlerProxy::new(handler, library.clone())) as Box<dyn MessageHandler>
    });

    for diagnostic in &registration.diagnostics {
        log::warn!("Plugin '{}': {}", metadata.name, diagnostic);
    }

    host.plugins.push(LoadedPlugin {
        _library: library,
        metadata,
        message_handler,
        panels,
        hotkeys,
        atom_streams: Default::default(),
        faulted: false,
    });
    host.ensure_single_active_per_placement();
    host.bump_panel_ui_generation();

    Ok(plugin_name)
}

fn load_declaration(library: &Library) -> Result<PluginDeclaration, String> {
    // SAFETY: The plugin is expected to export PATINAE_PLUGIN_DECLARATION as a
    // data symbol containing a pointer to a static PluginDeclaration.
    unsafe {
        let symbol = library
            .get::<*const PluginDeclaration>(b"PATINAE_PLUGIN_DECLARATION")
            .map_err(|e| format!("Missing PATINAE_PLUGIN_DECLARATION symbol: {e}"))?;
        let ptr = *symbol;
        if ptr.is_null() {
            return Err("PATINAE_PLUGIN_DECLARATION symbol was null".to_string());
        }
        Ok(*ptr)
    }
}

fn initialize_plugin(declaration: &PluginDeclaration) -> Result<(), String> {
    let Some(init) = declaration.init else {
        return Ok(());
    };
    let init_result = catch_unwind(AssertUnwindSafe(|| {
        // SAFETY: The declaration passed validation, and HOST_CALLBACKS is a
        // process-static callback table valid for the duration of this call.
        unsafe { init(&HOST_CALLBACKS) }
    }));
    match init_result {
        Ok(status) => status_to_result("init", status),
        Err(panic_info) => Err(format!(
            "Plugin panicked during init: {}",
            panic_payload_to_string(&panic_info)
        )),
    }
}

fn register_plugin(
    declaration: &PluginDeclaration,
    registration: &mut RegistrationSink,
) -> Result<(), String> {
    let register = declaration
        .register
        .ok_or("Plugin declaration missing register callback")?;
    let handle = HostRegistrarHandle(registration as *mut RegistrationSink as *mut c_void);
    let register_result = catch_unwind(AssertUnwindSafe(|| {
        // SAFETY: The declaration passed validation. The opaque handle points to
        // `registration`, which outlives this call, and HOST_CALLBACKS is static.
        unsafe { register(handle, &HOST_CALLBACKS) }
    }));
    match register_result {
        Ok(status) => status_to_result("register", status),
        Err(panic_info) => Err(format!(
            "Plugin panicked during register: {}",
            panic_payload_to_string(&panic_info)
        )),
    }
}

fn install_registration_assets(
    executor: &mut CommandExecutor,
    registration: &mut RegistrationSink,
    library: Arc<LibraryHandle>,
    gpu_cache: SharedGpuPluginCache,
) {
    for command in std::mem::take(&mut registration.commands) {
        executor
            .registry_mut()
            .register_boxed(Box::new(AbiCommandProxy::new(
                command,
                library.clone(),
                gpu_cache.clone(),
            )));
    }

    for script in std::mem::take(&mut registration.script_handlers) {
        executor.register_script_handler(
            script.extension.clone(),
            AbiScriptHandlerProxy::into_handler(script, library.clone()),
        );
    }

    for format in std::mem::take(&mut registration.format_handlers) {
        executor
            .register_format_handler(AbiFormatHandlerProxy::into_handler(format, library.clone()));
    }

    for (descriptor, store) in std::mem::take(&mut registration.settings) {
        if let Err(e) = executor.dynamic_settings_mut().register(descriptor, store) {
            log::warn!("Failed to register plugin setting: {}", e);
        }
    }
}

pub fn validate_declaration_versions(abi_version: u32, sdk_version: &str) -> Result<(), String> {
    if abi_version != ABI_VERSION {
        return Err(format!(
            "ABI version mismatch: plugin has {}, host expects {}",
            abi_version, ABI_VERSION
        ));
    }
    if sdk_version != SDK_VERSION {
        return Err(format!(
            "SDK version mismatch: plugin has {}, host expects {}",
            sdk_version, SDK_VERSION
        ));
    }
    Ok(())
}

pub(crate) fn validate_declaration(declaration: &PluginDeclaration) -> Result<(), String> {
    validate_declaration_versions(declaration.abi_version, &plugin_sdk_version(declaration)?)?;
    if declaration.register.is_none() {
        return Err("Plugin declaration missing register callback".to_string());
    }
    if declaration.capabilities & CAPABILITY_REGISTRATION == 0 {
        return Err("Plugin declaration missing registration capability".to_string());
    }
    let unknown = declaration.capabilities & !KNOWN_CAPABILITIES;
    if unknown != 0 {
        return Err(format!(
            "Plugin declaration has unknown capabilities: 0x{unknown:x}"
        ));
    }
    Ok(())
}

pub(crate) fn plugin_sdk_version(declaration: &PluginDeclaration) -> Result<String, String> {
    copy_str(declaration.sdk_version, "sdk_version")
}

fn status_to_result(phase: &str, status: AbiStatus) -> Result<(), String> {
    if status.is_ok() {
        return Ok(());
    }
    let message = match status.code {
        ABI_STATUS_PANIC => "panic",
        ABI_STATUS_INVALID => "invalid ABI input",
        ABI_STATUS_UNSUPPORTED => "unsupported feature",
        ABI_STATUS_HOST_ERROR => "host callback error",
        other => {
            return Err(format!(
                "Plugin {phase} returned unknown ABI status {other}"
            ))
        }
    };
    Err(format!("Plugin {phase} failed: {message}"))
}

#[derive(Default)]
struct RegistrationSink {
    metadata: Option<PluginMetadata>,
    commands: Vec<RegisteredCommand>,
    panels: Vec<RegisteredPanel>,
    message_handler: Option<RegisteredMessageHandler>,
    script_handlers: Vec<RegisteredScriptHandler>,
    format_handlers: Vec<RegisteredFormatHandler>,
    hotkeys: Vec<RegisteredHotkey>,
    settings: Vec<(DynamicSettingDescriptor, SharedSettingStore)>,
    diagnostics: Vec<String>,
}

impl RegistrationSink {
    fn new() -> Self {
        Self::default()
    }
}

struct RegisteredPanel {
    descriptor: PanelDescriptor,
    runtime_requirements: PanelRuntimeRequirements,
    handle: PluginPanelHandle,
    vtable: AbiPanelVTable,
}

unsafe extern "C" fn host_register_metadata(
    handle: HostRegistrarHandle,
    name: AbiStr,
    version: AbiStr,
    description: AbiStr,
) -> AbiStatus {
    callback_status(|| {
        let sink = sink_mut(handle)?;
        sink.metadata = Some(PluginMetadata::new(
            copy_str(name, "metadata.name")?,
            copy_str(version, "metadata.version")?,
            copy_str(description, "metadata.description")?,
        ));
        Ok(())
    })
}

unsafe extern "C" fn host_register_command(
    handle: HostRegistrarHandle,
    descriptor: *const AbiCommandDescriptor,
) -> AbiStatus {
    callback_status(|| {
        let sink = sink_mut(handle)?;
        if descriptor.is_null() {
            return Err("command descriptor pointer was null".to_string());
        }
        // SAFETY: The plugin provided a non-null descriptor pointer for the
        // duration of this callback. The host copies all borrowed fields here.
        let descriptor = unsafe { &*descriptor };
        sink.commands
            .push(RegisteredCommand::from_descriptor(descriptor)?);
        Ok(())
    })
}

unsafe extern "C" fn host_register_panel(
    handle: HostRegistrarHandle,
    descriptor: *const AbiPanelDescriptor,
) -> AbiStatus {
    callback_status(|| {
        let sink = sink_mut(handle)?;
        if descriptor.is_null() {
            return Err("panel descriptor pointer was null".to_string());
        }
        // SAFETY: The plugin provided a non-null descriptor pointer for the
        // duration of this callback. The host copies all borrowed fields here.
        let descriptor = unsafe { &*descriptor };
        validate_panel_runtime(descriptor)?;
        let panel_descriptor = panel_descriptor_from_abi(descriptor)?;
        sink.panels.push(RegisteredPanel {
            descriptor: panel_descriptor,
            runtime_requirements: PanelRuntimeRequirements::from_bits(
                descriptor.runtime_requirements,
            ),
            handle: descriptor.handle,
            vtable: descriptor.vtable,
        });
        Ok(())
    })
}

unsafe extern "C" fn host_register_setting(
    handle: HostRegistrarHandle,
    descriptor: *const AbiSettingDescriptor,
) -> AbiStatus {
    callback_status(|| {
        let sink = sink_mut(handle)?;
        if descriptor.is_null() {
            return Err("setting descriptor pointer was null".to_string());
        }
        // SAFETY: The plugin provided a non-null descriptor pointer for the
        // duration of this callback. The host copies all borrowed fields here.
        let descriptor = unsafe { &*descriptor };
        sink.settings.push(setting_from_abi(descriptor)?);
        Ok(())
    })
}

unsafe extern "C" fn host_register_message_handler(
    handle: HostRegistrarHandle,
    descriptor: *const AbiMessageHandlerDescriptor,
) -> AbiStatus {
    callback_status(|| {
        let sink = sink_mut(handle)?;
        if descriptor.is_null() {
            return Err("message handler descriptor pointer was null".to_string());
        }
        // SAFETY: The plugin provided a non-null descriptor pointer for the
        // duration of this callback. The host copies handle and vtable here.
        let descriptor = unsafe { &*descriptor };
        validate_message_handler_runtime(descriptor)?;
        sink.message_handler = Some(RegisteredMessageHandler {
            handle: descriptor.handle,
            vtable: descriptor.vtable,
        });
        Ok(())
    })
}

unsafe extern "C" fn host_register_script_handler(
    handle: HostRegistrarHandle,
    descriptor: *const AbiScriptHandlerDescriptor,
) -> AbiStatus {
    callback_status(|| {
        let sink = sink_mut(handle)?;
        if descriptor.is_null() {
            return Err("script handler descriptor pointer was null".to_string());
        }
        // SAFETY: The plugin provided a non-null descriptor pointer for the
        // duration of this callback. The host copies all borrowed fields here.
        let descriptor = unsafe { &*descriptor };
        validate_script_handler_runtime(descriptor)?;
        sink.script_handlers.push(RegisteredScriptHandler {
            extension: copy_str(descriptor.extension, "script.extension")?,
            handle: descriptor.handle,
            vtable: descriptor.vtable,
        });
        Ok(())
    })
}

unsafe extern "C" fn host_register_format_handler(
    handle: HostRegistrarHandle,
    descriptor: *const AbiFormatHandlerDescriptor,
) -> AbiStatus {
    callback_status(|| {
        let sink = sink_mut(handle)?;
        if descriptor.is_null() {
            return Err("format handler descriptor pointer was null".to_string());
        }
        // SAFETY: The plugin provided a non-null descriptor pointer for the
        // duration of this callback. The host copies all borrowed fields here.
        let descriptor = unsafe { &*descriptor };
        validate_format_handler_runtime(descriptor)?;
        sink.format_handlers.push(RegisteredFormatHandler {
            name: copy_str(descriptor.name, "format.name")?,
            extensions: copy_str_slice(descriptor.extensions, "format.extensions")?,
            handle: descriptor.handle,
            vtable: descriptor.vtable,
        });
        Ok(())
    })
}

unsafe extern "C" fn host_register_hotkey(
    handle: HostRegistrarHandle,
    descriptor: *const AbiHotkeyDescriptor,
) -> AbiStatus {
    callback_status(|| {
        let sink = sink_mut(handle)?;
        if descriptor.is_null() {
            return Err("hotkey descriptor pointer was null".to_string());
        }
        // SAFETY: The plugin provided a non-null descriptor pointer for the
        // duration of this callback. The host copies all borrowed fields here.
        let descriptor = unsafe { &*descriptor };
        sink.hotkeys
            .push(RegisteredHotkey::from_descriptor(descriptor)?);
        Ok(())
    })
}

unsafe extern "C" fn host_report_unsupported(
    handle: HostRegistrarHandle,
    feature: AbiStr,
    detail: AbiStr,
) -> AbiStatus {
    callback_status(|| {
        let sink = sink_mut(handle)?;
        let feature = copy_str(feature, "unsupported.feature")?;
        let detail = copy_str(detail, "unsupported.detail")?;
        sink.diagnostics.push(format!(
            "{feature} is unsupported by ABI v2 runtime: {detail}"
        ));
        Ok(())
    })
}

fn callback_status(f: impl FnOnce() -> Result<(), String>) -> AbiStatus {
    match catch_unwind(AssertUnwindSafe(f)) {
        Ok(Ok(())) => AbiStatus::OK,
        Ok(Err(error)) => {
            log::warn!("Plugin registration callback failed: {error}");
            AbiStatus::HOST_ERROR
        }
        Err(panic_info) => {
            log::error!(
                "Plugin registration callback panicked: {}",
                panic_payload_to_string(&panic_info)
            );
            AbiStatus::PANIC
        }
    }
}

fn sink_mut<'a>(handle: HostRegistrarHandle) -> Result<&'a mut RegistrationSink, String> {
    if handle.0.is_null() {
        return Err("host registrar handle was null".to_string());
    }
    // SAFETY: The handle is created from a live `RegistrationSink` in
    // `register_plugin` and is only used during that synchronous registration.
    Ok(unsafe { &mut *(handle.0.cast::<RegistrationSink>()) })
}

fn copy_str(value: AbiStr, field: &str) -> Result<String, String> {
    // SAFETY: The ABI contract requires string views to point to initialized
    // UTF-8 bytes. This call validates null and length before reading.
    let bytes = unsafe { value.as_bytes_checked(MAX_ABI_STRING_LEN) }
        .map_err(|error| abi_read_error(field, error))?;
    let value =
        std::str::from_utf8(bytes).map_err(|_| abi_read_error(field, AbiReadError::InvalidUtf8))?;
    Ok(value.to_string())
}

fn copy_str_slice(value: AbiStrSlice, field: &str) -> Result<Vec<String>, String> {
    if value.len > MAX_ABI_SLICE_LEN {
        return Err(format!("{field} slice length exceeds ABI limit"));
    }
    if value.len == 0 {
        return Ok(Vec::new());
    }
    if value.ptr.is_null() {
        return Err(format!("{field} slice pointer was null"));
    }
    // SAFETY: The ABI contract requires non-null slices to point to `len`
    // initialized AbiStr entries for the duration of the callback.
    let values = unsafe { std::slice::from_raw_parts(value.ptr, value.len) };
    values
        .iter()
        .enumerate()
        .map(|(index, value)| copy_str(*value, &format!("{field}[{index}]")))
        .collect()
}

fn copy_u8_slice(value: AbiU8Slice, field: &str) -> Result<Vec<u8>, String> {
    if value.len > MAX_ABI_SLICE_LEN {
        return Err(format!("{field} slice length exceeds ABI limit"));
    }
    if value.len == 0 {
        return Ok(Vec::new());
    }
    if value.ptr.is_null() {
        return Err(format!("{field} slice pointer was null"));
    }
    // SAFETY: The ABI contract requires non-null slices to point to `len`
    // initialized bytes for the duration of the callback.
    Ok(unsafe { std::slice::from_raw_parts(value.ptr, value.len) }.to_vec())
}

fn copy_runtime_u8_slice(value: AbiU8Slice, field: &str) -> Result<Vec<u8>, String> {
    if value.len > wire::MAX_WIRE_PAYLOAD_LEN {
        return Err(format!("{field} runtime payload exceeds ABI limit"));
    }
    if value.len == 0 {
        return Ok(Vec::new());
    }
    if value.ptr.is_null() {
        return Err(format!("{field} runtime payload pointer was null"));
    }
    // SAFETY: The ABI contract requires non-null runtime byte slices to point
    // to `len` initialized bytes for the duration of the callback.
    Ok(unsafe { std::slice::from_raw_parts(value.ptr, value.len) }.to_vec())
}

fn copy_hint_slice(
    value: patinae_plugin::ffi::AbiSettingValueHintSlice,
    field: &str,
) -> Result<Vec<(String, SettingValue)>, String> {
    if value.len > MAX_ABI_SLICE_LEN {
        return Err(format!("{field} slice length exceeds ABI limit"));
    }
    if value.len == 0 {
        return Ok(Vec::new());
    }
    if value.ptr.is_null() {
        return Err(format!("{field} slice pointer was null"));
    }
    // SAFETY: The ABI contract requires non-null slices to point to `len`
    // initialized hint entries for the duration of the callback.
    let values = unsafe { std::slice::from_raw_parts(value.ptr, value.len) };
    values
        .iter()
        .enumerate()
        .map(|(index, value)| setting_hint_from_abi(value, &format!("{field}[{index}]")))
        .collect()
}

fn abi_read_error(field: &str, error: AbiReadError) -> String {
    match error {
        AbiReadError::NullPointer => format!("{field} pointer was null"),
        AbiReadError::TooLong => format!("{field} length exceeds ABI limit"),
        AbiReadError::InvalidUtf8 => format!("{field} was not valid UTF-8"),
    }
}

#[derive(Default)]
struct RuntimeBytes {
    bytes: Vec<u8>,
}

unsafe extern "C" fn host_runtime_bytes_sink(
    user_data: *mut c_void,
    bytes: AbiU8Slice,
) -> AbiStatus {
    if user_data.is_null() {
        return AbiStatus::INVALID;
    }
    let copied = match copy_runtime_u8_slice(bytes, "runtime.output") {
        Ok(bytes) => bytes,
        Err(error) => {
            log::warn!("Plugin runtime byte sink rejected output: {}", error);
            return AbiStatus::INVALID;
        }
    };
    // SAFETY: The caller passes a pointer to a live `RuntimeBytes` for the
    // duration of this sink call. Null was checked above.
    unsafe { &mut *(user_data.cast::<RuntimeBytes>()) }.bytes = copied;
    AbiStatus::OK
}

fn invoke_runtime<I, O>(
    phase: &str,
    input: &I,
    call: impl FnOnce(AbiU8Slice, AbiBytesSinkFn, *mut c_void) -> AbiStatus,
) -> Result<O, String>
where
    I: Serialize,
    O: DeserializeOwned,
{
    let input_bytes = wire::encode(input)?;
    ensure_runtime_payload_within_limit(phase, input_bytes.len(), "runtime input")?;
    let input_slice = AbiU8Slice {
        ptr: input_bytes.as_ptr(),
        len: input_bytes.len(),
    };
    let mut output = RuntimeBytes::default();
    let output_ptr = (&mut output as *mut RuntimeBytes).cast::<c_void>();
    let status = catch_unwind(AssertUnwindSafe(|| {
        call(input_slice, host_runtime_bytes_sink, output_ptr)
    }))
    .map_err(|panic_info| {
        format!(
            "plugin {phase} panicked across host boundary: {}",
            panic_payload_to_string(&panic_info)
        )
    })?;
    status_to_result(phase, status)?;
    wire::decode(&output.bytes)
        .map_err(|error| format!("plugin {phase} returned malformed MessagePack: {error}"))
}

unsafe extern "C" fn host_command_runtime_request(
    handle: HostCommandRuntimeHandle,
    input: AbiU8Slice,
    sink: AbiBytesSinkFn,
    sink_user_data: *mut c_void,
) -> AbiStatus {
    if handle.0.is_null() {
        return AbiStatus::INVALID;
    }

    let result = catch_unwind(AssertUnwindSafe(|| {
        let request: WireCommandRuntimeRequest = decode_runtime_input(input)?;
        // SAFETY: `AbiCommandProxy::execute` passes a pointer to a live
        // `HostCommandRuntimeState` for the duration of the plugin command
        // call. Null was checked above.
        let state = unsafe { &mut *(handle.0.cast::<HostCommandRuntimeState<'static>>()) };
        let response = state.handle_request(request);
        send_runtime_output(sink, sink_user_data, &response)
    }));

    match result {
        Ok(Ok(())) => AbiStatus::OK,
        Ok(Err(error)) => {
            log::warn!("Plugin command runtime request failed: {}", error);
            AbiStatus::HOST_ERROR
        }
        Err(panic_info) => {
            log::warn!(
                "Plugin command runtime request panicked: {}",
                panic_payload_to_string(&panic_info)
            );
            AbiStatus::PANIC
        }
    }
}

fn decode_runtime_input<T: DeserializeOwned>(input: AbiU8Slice) -> Result<T, String> {
    if input.len == 0 {
        return wire::decode(&[]);
    }
    if input.ptr.is_null() {
        return Err("runtime request pointer was null".to_string());
    }
    if input.len > wire::MAX_WIRE_PAYLOAD_LEN {
        return Err("runtime request exceeds ABI limit".to_string());
    }
    // SAFETY: The plugin ABI contract requires non-null byte slices to remain
    // initialized and readable for the duration of this callback.
    let bytes = unsafe { std::slice::from_raw_parts(input.ptr, input.len) };
    wire::decode(bytes)
}

fn send_runtime_output<T: Serialize>(
    sink: AbiBytesSinkFn,
    sink_user_data: *mut c_void,
    output: &T,
) -> Result<(), String> {
    let bytes = wire::encode(output)?;
    ensure_runtime_payload_within_limit("command runtime request", bytes.len(), "runtime output")?;
    let slice = AbiU8Slice {
        ptr: bytes.as_ptr(),
        len: bytes.len(),
    };
    // SAFETY: `slice` points to `bytes`, which lives until the sink returns.
    // The plugin sink copies the bytes synchronously.
    let status = unsafe { sink(sink_user_data, slice) };
    status_to_result("command runtime request", status)
}

fn validate_runtime_wire_version(version: u32) -> Result<(), String> {
    if version == RUNTIME_WIRE_VERSION {
        Ok(())
    } else {
        Err(format!(
            "runtime wire version mismatch: plugin has {}, host expects {}",
            version, RUNTIME_WIRE_VERSION
        ))
    }
}

fn validate_command_runtime(descriptor: &AbiCommandDescriptor) -> Result<(), String> {
    if descriptor.handle.0.is_null() {
        return Err("command handle was null".to_string());
    }
    if descriptor.vtable.execute.is_none() {
        return Err("command execute callback was null".to_string());
    }
    if descriptor.vtable.destroy.is_none() {
        return Err("command destroy callback was null".to_string());
    }
    if descriptor.runtime_requirements & !KNOWN_COMMAND_RUNTIME_REQUIREMENTS != 0 {
        return Err("command runtime requirements contain unknown bits".to_string());
    }
    Ok(())
}

fn validate_panel_runtime(descriptor: &AbiPanelDescriptor) -> Result<(), String> {
    if descriptor.handle.0.is_null() {
        return Err("panel handle was null".to_string());
    }
    if descriptor.vtable.snapshot.is_none() {
        return Err("panel snapshot callback was null".to_string());
    }
    if descriptor.vtable.handle_event.is_none() {
        return Err("panel event callback was null".to_string());
    }
    if descriptor.vtable.destroy.is_none() {
        return Err("panel destroy callback was null".to_string());
    }
    if descriptor.runtime_requirements & !KNOWN_PANEL_RUNTIME_REQUIREMENTS != 0 {
        return Err("panel runtime requirements contain unknown bits".to_string());
    }
    Ok(())
}

fn validate_message_handler_runtime(
    descriptor: &AbiMessageHandlerDescriptor,
) -> Result<(), String> {
    if descriptor.handle.0.is_null() {
        return Err("message handler handle was null".to_string());
    }
    if descriptor.vtable.on_message.is_none() {
        return Err("message handler callback was null".to_string());
    }
    if descriptor.vtable.poll.is_none() {
        return Err("message poll callback was null".to_string());
    }
    if descriptor.vtable.destroy.is_none() {
        return Err("message handler destroy callback was null".to_string());
    }
    Ok(())
}

fn validate_script_handler_runtime(descriptor: &AbiScriptHandlerDescriptor) -> Result<(), String> {
    if descriptor.handle.0.is_null() {
        return Err("script handler handle was null".to_string());
    }
    if descriptor.vtable.run.is_none() {
        return Err("script handler callback was null".to_string());
    }
    if descriptor.vtable.destroy.is_none() {
        return Err("script handler destroy callback was null".to_string());
    }
    Ok(())
}

fn validate_format_handler_runtime(descriptor: &AbiFormatHandlerDescriptor) -> Result<(), String> {
    if descriptor.handle.0.is_null() {
        return Err("format handler handle was null".to_string());
    }
    if descriptor.vtable.read.is_none() && descriptor.vtable.write.is_none() {
        return Err("format handler has no reader or writer callback".to_string());
    }
    if descriptor.vtable.destroy.is_none() {
        return Err("format handler destroy callback was null".to_string());
    }
    Ok(())
}

fn validate_hotkey_callback_runtime(descriptor: &AbiHotkeyDescriptor) -> Result<(), String> {
    if descriptor.handle.0.is_null() {
        return Err("hotkey action and callback handle were both empty".to_string());
    }
    if descriptor.vtable.invoke.is_none() {
        return Err("hotkey callback was null".to_string());
    }
    if descriptor.vtable.destroy.is_none() {
        return Err("hotkey destroy callback was null".to_string());
    }
    Ok(())
}

struct RegisteredCommand {
    name: String,
    description: String,
    usage: String,
    arguments: String,
    help: String,
    aliases: &'static [&'static str],
    arg_hints: &'static [ArgHint],
    runtime_requirements: CommandRuntimeRequirements,
    handle: PluginCommandHandle,
    vtable: AbiCommandVTable,
}

impl RegisteredCommand {
    fn from_descriptor(descriptor: &AbiCommandDescriptor) -> Result<Self, String> {
        validate_command_runtime(descriptor)?;
        Ok(Self {
            name: copy_str(descriptor.name, "command.name")?,
            description: copy_str(descriptor.description, "command.description")?,
            usage: copy_str(descriptor.usage, "command.usage")?,
            arguments: copy_str(descriptor.arguments, "command.arguments")?,
            help: copy_str(descriptor.help, "command.help")?,
            aliases: leak_aliases(copy_str_slice(descriptor.aliases, "command.aliases")?),
            arg_hints: leak_arg_hints(
                copy_u8_slice(descriptor.arg_hints, "command.arg_hints")?
                    .into_iter()
                    .map(arg_hint_from_abi)
                    .collect(),
            ),
            runtime_requirements: CommandRuntimeRequirements::from_bits(
                descriptor.runtime_requirements,
            ),
            handle: descriptor.handle,
            vtable: descriptor.vtable,
        })
    }
}

struct RegisteredMessageHandler {
    handle: PluginMessageHandlerHandle,
    vtable: AbiMessageHandlerVTable,
}

struct RegisteredScriptHandler {
    extension: String,
    handle: PluginScriptHandlerHandle,
    vtable: AbiScriptHandlerVTable,
}

struct RegisteredFormatHandler {
    name: String,
    extensions: Vec<String>,
    handle: PluginFormatHandlerHandle,
    vtable: AbiFormatHandlerVTable,
}

struct RegisteredHotkey {
    key: String,
    action: Option<WireHotkeyAction>,
    handle: PluginHotkeyHandle,
    vtable: AbiHotkeyVTable,
}

impl RegisteredHotkey {
    fn from_descriptor(descriptor: &AbiHotkeyDescriptor) -> Result<Self, String> {
        let action_bytes = copy_runtime_u8_slice(descriptor.action, "hotkey.action")?;
        let action = if action_bytes.is_empty() {
            validate_hotkey_callback_runtime(descriptor)?;
            None
        } else {
            Some(
                wire::decode(&action_bytes)
                    .map_err(|error| format!("hotkey.action was not valid MessagePack: {error}"))?,
            )
        };
        Ok(Self {
            key: copy_str(descriptor.key, "hotkey.key")?,
            action,
            handle: descriptor.handle,
            vtable: descriptor.vtable,
        })
    }

    fn into_action(self, library: Arc<LibraryHandle>) -> PluginKeyAction {
        match self.action {
            Some(action) => plugin_key_action_from_wire(action),
            None => {
                let proxy = AbiHotkeyProxy::new(self.handle, self.vtable, library);
                PluginKeyAction::Callback(Box::new(move |ctx| {
                    if let Err(error) = proxy.invoke(ctx) {
                        ctx.set_notification(format!("plugin hotkey callback failed: {error}"));
                    }
                }))
            }
        }
    }
}

struct AbiCommandProxy {
    command: RegisteredCommand,
    _library: Arc<LibraryHandle>,
    gpu_cache: SharedGpuPluginCache,
}

impl AbiCommandProxy {
    fn new(
        command: RegisteredCommand,
        library: Arc<LibraryHandle>,
        gpu_cache: SharedGpuPluginCache,
    ) -> Self {
        Self {
            command,
            _library: library,
            gpu_cache,
        }
    }
}

// SAFETY: The proxy never dereferences the raw plugin handle locally. All
// access goes through the plugin vtable while the library is kept alive. Plugin
// commands implement `Command: Send + Sync` before the SDK exports a handle.
unsafe impl Send for AbiCommandProxy {}
// SAFETY: See the `Send` rationale; command callbacks are immutable command
// invocations and plugin command objects are required to be `Sync`.
unsafe impl Sync for AbiCommandProxy {}

impl Command for AbiCommandProxy {
    fn name(&self) -> &str {
        &self.command.name
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        let input =
            command_runtime_input_from_context(ctx, args, self.command.runtime_requirements)
                .map_err(|error| CmdError::execution(format!("plugin command input: {error}")))?;
        let execute = self
            .command
            .vtable
            .execute
            .ok_or_else(|| CmdError::execution("plugin command execute callback was null"))?;
        let output: WireCommandOutput = {
            let mut runtime_state = HostCommandRuntimeState::with_cache(
                ctx.viewer,
                self.command.runtime_requirements,
                self.gpu_cache.clone(),
            );
            invoke_runtime("command execute", &input.input, |slice, sink, user_data| {
                // SAFETY: The descriptor was validated during registration and the
                // library is kept alive by this proxy while the callback runs.
                unsafe {
                    execute(
                        self.command.handle,
                        slice,
                        &HOST_COMMAND_RUNTIME_CALLBACKS,
                        runtime_state.handle(),
                        sink,
                        user_data,
                    )
                }
            })
        }
        .map_err(CmdError::execution)?;
        apply_command_output(ctx, output, self.command.runtime_requirements)
    }

    fn description(&self) -> &str {
        &self.command.description
    }

    fn usage(&self) -> &str {
        &self.command.usage
    }

    fn arguments(&self) -> &str {
        &self.command.arguments
    }

    fn help(&self) -> &str {
        &self.command.help
    }

    fn aliases(&self) -> &[&str] {
        self.command.aliases
    }

    fn arg_hints(&self) -> &[ArgHint] {
        self.command.arg_hints
    }

    fn runtime_requirements(&self) -> CommandRuntimeRequirements {
        self.command.runtime_requirements
    }
}

impl Drop for AbiCommandProxy {
    fn drop(&mut self) {
        destroy_command_handle(self.command.handle, self.command.vtable);
    }
}

struct AbiPanelProxy {
    panel: RegisteredPanel,
    _library: Arc<LibraryHandle>,
}

impl AbiPanelProxy {
    fn new(panel: RegisteredPanel, library: Arc<LibraryHandle>) -> Self {
        Self {
            panel,
            _library: library,
        }
    }
}

// SAFETY: The proxy keeps the dynamic library alive and only calls validated
// vtable entries with ABI byte payloads. Panel methods require mutable access.
unsafe impl Send for AbiPanelProxy {}

impl PluginPanel for AbiPanelProxy {
    fn descriptor(&self) -> PanelDescriptor {
        self.panel.descriptor.clone()
    }

    fn snapshot(&mut self, ctx: &SharedContext<'_>) -> PanelSnapshot {
        match self.snapshot_inner(ctx) {
            Ok(snapshot) => snapshot,
            Err(error) => PanelSnapshot::new(vec![PanelControl::Text {
                id: "abi_v2_error".to_string(),
                text: format!("Plugin panel failed: {error}"),
            }]),
        }
    }

    fn handle_event(
        &mut self,
        event: PanelEvent,
        ctx: &SharedContext<'_>,
        _bus: &mut MessageBus,
    ) -> Vec<patinae_framework::plugin_ui::PanelAction> {
        match self.handle_event_inner(event, ctx) {
            Ok(actions) => actions,
            Err(error) => {
                log::warn!(
                    "Plugin panel '{}' event failed: {}",
                    self.panel.descriptor.id,
                    error
                );
                Vec::new()
            }
        }
    }
}

impl AbiPanelProxy {
    fn snapshot_inner(&mut self, ctx: &SharedContext<'_>) -> Result<PanelSnapshot, String> {
        let snapshot = shared_input_from_context(ctx, self.panel.runtime_requirements)?;
        let callback = self
            .panel
            .vtable
            .snapshot
            .ok_or("plugin panel snapshot callback was null")?;
        let output: WirePanelSnapshotOutput =
            invoke_runtime("panel snapshot", &snapshot, |slice, sink, user_data| {
                // SAFETY: The descriptor was validated during registration and
                // this proxy keeps the library alive while the callback runs.
                unsafe { callback(self.panel.handle, slice, sink, user_data) }
            })?;
        validate_runtime_wire_version(output.wire_version)?;
        Ok(output.snapshot)
    }

    fn handle_event_inner(
        &mut self,
        event: PanelEvent,
        ctx: &SharedContext<'_>,
    ) -> Result<Vec<patinae_framework::plugin_ui::PanelAction>, String> {
        let input = WirePanelEventInput {
            shared: shared_input_from_context(ctx, self.panel.runtime_requirements)?,
            event,
        };
        let callback = self
            .panel
            .vtable
            .handle_event
            .ok_or("plugin panel event callback was null")?;
        let output: WirePanelEventOutput =
            invoke_runtime("panel event", &input, |slice, sink, user_data| {
                // SAFETY: The descriptor was validated during registration and
                // this proxy keeps the library alive while the callback runs.
                unsafe { callback(self.panel.handle, slice, sink, user_data) }
            })?;
        validate_runtime_wire_version(output.wire_version)?;
        Ok(output.actions)
    }
}

impl Drop for AbiPanelProxy {
    fn drop(&mut self) {
        destroy_panel_handle(self.panel.handle, self.panel.vtable);
    }
}

struct AbiMessageHandlerProxy {
    handle: PluginMessageHandlerHandle,
    vtable: AbiMessageHandlerVTable,
    _library: Arc<LibraryHandle>,
}

impl AbiMessageHandlerProxy {
    fn new(handler: RegisteredMessageHandler, library: Arc<LibraryHandle>) -> Self {
        Self {
            handle: handler.handle,
            vtable: handler.vtable,
            _library: library,
        }
    }
}

// SAFETY: The proxy keeps the dynamic library alive and only calls validated
// vtable entries with ABI byte payloads. `MessageHandler` only requires `Send`.
unsafe impl Send for AbiMessageHandlerProxy {}

impl MessageHandler for AbiMessageHandlerProxy {
    fn on_message(&mut self, msg: &AppMessage, bus: &mut MessageBus) {
        if let Err(error) = self.on_message_inner(msg, bus) {
            log::warn!("Plugin message handler failed: {}", error);
        }
    }

    fn needs_poll(&self) -> bool {
        self.vtable.needs_poll != 0
    }

    fn poll(&mut self, ctx: &mut PollContext<'_>) {
        if let Err(error) = self.poll_inner(ctx) {
            ctx.set_notification(format!("plugin poll failed: {error}"));
        }
    }
}

impl AbiMessageHandlerProxy {
    fn on_message_inner(&mut self, msg: &AppMessage, bus: &mut MessageBus) -> Result<(), String> {
        let input = WireMessageInput {
            wire_version: RUNTIME_WIRE_VERSION,
            message: msg.clone(),
        };
        let callback = self
            .vtable
            .on_message
            .ok_or("plugin message callback was null")?;
        let output: WireMessageOutput =
            invoke_runtime("message handler", &input, |slice, sink, user_data| {
                // SAFETY: The descriptor was validated during registration and
                // this proxy keeps the library alive while the callback runs.
                unsafe { callback(self.handle, slice, sink, user_data) }
            })?;
        validate_runtime_wire_version(output.wire_version)?;
        send_bus_messages(output.messages, bus);
        Ok(())
    }

    fn poll_inner(&mut self, ctx: &mut PollContext<'_>) -> Result<(), String> {
        let input = poll_input_from_context(ctx)?;
        let callback = self.vtable.poll.ok_or("plugin poll callback was null")?;
        let output: WirePollOutput =
            invoke_runtime("message poll", &input, |slice, sink, user_data| {
                // SAFETY: The descriptor was validated during registration and
                // this proxy keeps the library alive while the callback runs.
                unsafe { callback(self.handle, slice, sink, user_data) }
            })?;
        apply_poll_output(ctx, output)
    }
}

impl Drop for AbiMessageHandlerProxy {
    fn drop(&mut self) {
        destroy_message_handler_handle(self.handle, self.vtable);
    }
}

struct AbiScriptHandlerProxy {
    extension: String,
    handle: PluginScriptHandlerHandle,
    vtable: AbiScriptHandlerVTable,
    _library: Arc<LibraryHandle>,
}

impl AbiScriptHandlerProxy {
    fn into_handler(script: RegisteredScriptHandler, library: Arc<LibraryHandle>) -> ScriptHandler {
        let proxy = Arc::new(Self {
            extension: script.extension,
            handle: script.handle,
            vtable: script.vtable,
            _library: library,
        });
        Arc::new(move |path: &str| proxy.run(path))
    }

    fn run(&self, path: &str) -> Result<(), String> {
        let input = WireScriptInput {
            wire_version: RUNTIME_WIRE_VERSION,
            path: path.to_string(),
        };
        let callback = self.vtable.run.ok_or("plugin script callback was null")?;
        let output: WireScriptOutput =
            invoke_runtime("script handler", &input, |slice, sink, user_data| {
                // SAFETY: The descriptor was validated during registration and
                // this proxy keeps the library alive while the callback runs.
                unsafe { callback(self.handle, slice, sink, user_data) }
            })?;
        validate_runtime_wire_version(output.wire_version)?;
        output
            .result
            .map_err(|error| format!("plugin script .{} failed: {error}", self.extension))
    }
}

// SAFETY: The proxy never dereferences the raw plugin handle locally and the
// SDK only exports script handlers from `Send + Sync` closures.
unsafe impl Send for AbiScriptHandlerProxy {}
// SAFETY: See the `Send` rationale; calls go through validated ABI callbacks.
unsafe impl Sync for AbiScriptHandlerProxy {}

impl Drop for AbiScriptHandlerProxy {
    fn drop(&mut self) {
        destroy_script_handler_handle(self.handle, self.vtable);
    }
}

struct AbiFormatHandlerProxy {
    name: String,
    handle: PluginFormatHandlerHandle,
    vtable: AbiFormatHandlerVTable,
    _library: Arc<LibraryHandle>,
}

impl AbiFormatHandlerProxy {
    fn into_handler(format: RegisteredFormatHandler, library: Arc<LibraryHandle>) -> FormatHandler {
        let name = format.name;
        let extensions = format.extensions;
        let vtable = format.vtable;
        let proxy = Arc::new(Self {
            name: name.clone(),
            handle: format.handle,
            vtable,
            _library: library,
        });
        let reader: Option<PluginReaderFn> = if vtable.read.is_some() {
            let proxy = proxy.clone();
            Some(Arc::new(move |reader: Box<dyn Read>| proxy.read(reader)))
        } else {
            None
        };
        let writer: Option<PluginWriterFn> = if vtable.write.is_some() {
            let proxy = proxy.clone();
            Some(Arc::new(
                move |writer: Box<dyn Write>, molecules: &[ObjectMolecule]| {
                    proxy.write(writer, molecules)
                },
            ))
        } else {
            None
        };
        FormatHandler {
            name,
            extensions,
            reader,
            writer,
        }
    }

    fn read(&self, mut reader: Box<dyn Read>) -> Result<Vec<ObjectMolecule>, String> {
        let mut bytes = Vec::new();
        reader
            .read_to_end(&mut bytes)
            .map_err(|error| format!("failed to read input bytes: {error}"))?;
        let input = WireFormatReadInput {
            wire_version: RUNTIME_WIRE_VERSION,
            extension: String::new(),
            bytes,
        };
        let callback = self
            .vtable
            .read
            .ok_or("plugin format reader callback was null")?;
        let output: WireFormatReadOutput =
            invoke_runtime("format reader", &input, |slice, sink, user_data| {
                // SAFETY: The descriptor was validated during registration and
                // this proxy keeps the library alive while the callback runs.
                unsafe { callback(self.handle, slice, sink, user_data) }
            })?;
        validate_runtime_wire_version(output.wire_version)?;
        output
            .result
            .map_err(|error| format!("plugin format '{}' read failed: {error}", self.name))
    }

    fn write(
        &self,
        mut writer: Box<dyn Write>,
        molecules: &[ObjectMolecule],
    ) -> Result<(), String> {
        let input = WireFormatWriteInput {
            wire_version: RUNTIME_WIRE_VERSION,
            extension: String::new(),
            molecules: molecules.to_vec(),
        };
        let callback = self
            .vtable
            .write
            .ok_or("plugin format writer callback was null")?;
        let output: WireFormatWriteOutput =
            invoke_runtime("format writer", &input, |slice, sink, user_data| {
                // SAFETY: The descriptor was validated during registration and
                // this proxy keeps the library alive while the callback runs.
                unsafe { callback(self.handle, slice, sink, user_data) }
            })?;
        validate_runtime_wire_version(output.wire_version)?;
        let bytes = output
            .result
            .map_err(|error| format!("plugin format '{}' write failed: {error}", self.name))?;
        writer
            .write_all(&bytes)
            .map_err(|error| format!("failed to write plugin output bytes: {error}"))
    }
}

// SAFETY: The proxy never dereferences the raw plugin handle locally and the
// SDK only exports format handlers from `Send + Sync` reader/writer closures.
unsafe impl Send for AbiFormatHandlerProxy {}
// SAFETY: See the `Send` rationale; calls go through validated ABI callbacks.
unsafe impl Sync for AbiFormatHandlerProxy {}

impl Drop for AbiFormatHandlerProxy {
    fn drop(&mut self) {
        destroy_format_handler_handle(self.handle, self.vtable);
    }
}

struct AbiHotkeyProxy {
    handle: PluginHotkeyHandle,
    vtable: AbiHotkeyVTable,
    _library: Arc<LibraryHandle>,
}

impl AbiHotkeyProxy {
    fn new(
        handle: PluginHotkeyHandle,
        vtable: AbiHotkeyVTable,
        library: Arc<LibraryHandle>,
    ) -> Self {
        Self {
            handle,
            vtable,
            _library: library,
        }
    }

    fn invoke(&self, ctx: &mut PollContext<'_>) -> Result<(), String> {
        let input = poll_input_from_context(ctx)?;
        let callback = self
            .vtable
            .invoke
            .ok_or("plugin hotkey callback was null")?;
        let output: WirePollOutput =
            invoke_runtime("hotkey callback", &input, |slice, sink, user_data| {
                // SAFETY: The descriptor was validated during registration and
                // this proxy keeps the library alive while the callback runs.
                unsafe { callback(self.handle, slice, sink, user_data) }
            })?;
        apply_poll_output(ctx, output)
    }
}

// SAFETY: The proxy is moved into a host-owned callback and only calls through
// validated ABI callbacks while keeping the dynamic library alive.
unsafe impl Send for AbiHotkeyProxy {}

impl Drop for AbiHotkeyProxy {
    fn drop(&mut self) {
        destroy_hotkey_handle(self.handle, self.vtable);
    }
}

struct CommandRuntimeInput {
    input: WireCommandInput,
    _spools: Vec<RuntimePayloadSpool>,
}

struct HostCommandRuntimeState<'a> {
    viewer: &'a mut dyn ViewerLike,
    runtime_requirements: CommandRuntimeRequirements,
    streams: Vec<HostTraceGeometryStream>,
    gpu_handles: GpuHandleRegistry,
    gpu_cache: SharedGpuPluginCache,
    snapshots: Vec<HostRenderArtifactSnapshot>,
    next_stream_id: u64,
    next_snapshot_id: u64,
}

impl<'a> HostCommandRuntimeState<'a> {
    #[cfg(test)]
    fn new(
        viewer: &'a mut dyn ViewerLike,
        runtime_requirements: CommandRuntimeRequirements,
    ) -> Self {
        Self::with_cache(
            viewer,
            runtime_requirements,
            Arc::new(Mutex::new(GpuPluginCache::new("test".to_string()))),
        )
    }

    fn with_cache(
        viewer: &'a mut dyn ViewerLike,
        runtime_requirements: CommandRuntimeRequirements,
        gpu_cache: SharedGpuPluginCache,
    ) -> Self {
        Self {
            viewer,
            runtime_requirements,
            streams: Vec::new(),
            gpu_handles: GpuHandleRegistry::new(),
            gpu_cache,
            snapshots: Vec::new(),
            next_stream_id: 1,
            next_snapshot_id: 1,
        }
    }

    fn handle(&mut self) -> HostCommandRuntimeHandle {
        HostCommandRuntimeHandle((self as *mut Self).cast::<c_void>())
    }

    fn handle_request(&mut self, request: WireCommandRuntimeRequest) -> WireCommandRuntimeResponse {
        match request {
            WireCommandRuntimeRequest::OpenTraceGeometryStream { id } => {
                let result = self.open_trace_geometry_stream();
                Self::response(id, result)
            }
            WireCommandRuntimeRequest::ReadTraceGeometryStream { id, stream_id } => {
                let result = self.read_trace_geometry_stream(stream_id);
                Self::response(id, result)
            }
            WireCommandRuntimeRequest::CloseTraceGeometryStream { id, stream_id } => {
                let result = self.close_trace_geometry_stream(stream_id);
                Self::response(id, result)
            }
            WireCommandRuntimeRequest::OpenRenderArtifactSnapshot { id } => {
                let result = self.open_render_artifact_snapshot();
                Self::response(id, result)
            }
            WireCommandRuntimeRequest::CloseRenderArtifactSnapshot { id, snapshot_id } => {
                let result = self.close_render_artifact_snapshot(snapshot_id);
                Self::response(id, result)
            }
            WireCommandRuntimeRequest::GpuDeviceLimits { id } => {
                let result = self.gpu_device_limits();
                Self::response(id, result)
            }
            WireCommandRuntimeRequest::GpuCreateBuffer {
                id,
                descriptor,
                initial_data,
            } => {
                let result = self.gpu_create_buffer(descriptor, initial_data);
                Self::response(id, result)
            }
            WireCommandRuntimeRequest::GpuWriteBuffer {
                id,
                buffer,
                offset,
                data,
            } => {
                let result = self.gpu_write_buffer(buffer, offset, &data);
                Self::response(id, result)
            }
            WireCommandRuntimeRequest::GpuCopyBufferToBuffer {
                id,
                source,
                source_offset,
                destination,
                destination_offset,
                size,
            } => {
                let result = self.gpu_copy_buffer_to_buffer(
                    source,
                    source_offset,
                    destination,
                    destination_offset,
                    size,
                );
                Self::response(id, result)
            }
            WireCommandRuntimeRequest::GpuReadBuffer {
                id,
                buffer,
                offset,
                size,
            } => {
                let result = self.gpu_read_buffer(buffer, offset, size);
                Self::response(id, result)
            }
            WireCommandRuntimeRequest::GpuCreateShaderModule { id, descriptor } => {
                let result = self.gpu_create_shader_module(descriptor);
                Self::response(id, result)
            }
            WireCommandRuntimeRequest::GpuCreateBindGroupLayout { id, descriptor } => {
                let result = self.gpu_create_bind_group_layout(descriptor);
                Self::response(id, result)
            }
            WireCommandRuntimeRequest::GpuCreatePipelineLayout { id, descriptor } => {
                let result = self.gpu_create_pipeline_layout(descriptor);
                Self::response(id, result)
            }
            WireCommandRuntimeRequest::GpuCreateComputePipeline { id, descriptor } => {
                let result = self.gpu_create_compute_pipeline(descriptor);
                Self::response(id, result)
            }
            WireCommandRuntimeRequest::GpuCreateCachedShaderModule { id, descriptor } => {
                let result = self.gpu_create_cached_shader_module(descriptor);
                Self::response(id, result)
            }
            WireCommandRuntimeRequest::GpuCreateCachedBindGroupLayout { id, descriptor } => {
                let result = self.gpu_create_cached_bind_group_layout(descriptor);
                Self::response(id, result)
            }
            WireCommandRuntimeRequest::GpuCreateCachedPipelineLayout { id, descriptor } => {
                let result = self.gpu_create_cached_pipeline_layout(descriptor);
                Self::response(id, result)
            }
            WireCommandRuntimeRequest::GpuCreateCachedComputePipeline { id, descriptor } => {
                let result = self.gpu_create_cached_compute_pipeline(descriptor);
                Self::response(id, result)
            }
            WireCommandRuntimeRequest::GpuCacheStats { id } => {
                let result = self.gpu_cache_stats();
                Self::response(id, result)
            }
            WireCommandRuntimeRequest::GpuDropPluginCache { id } => {
                let result = self.gpu_drop_plugin_cache();
                Self::response(id, result)
            }
            WireCommandRuntimeRequest::GpuCreateBindGroup { id, descriptor } => {
                let result = self.gpu_create_bind_group(descriptor);
                Self::response(id, result)
            }
            WireCommandRuntimeRequest::GpuDispatchCompute {
                id,
                pipeline,
                bind_groups,
                workgroups,
            } => {
                let result = self.gpu_dispatch_compute(pipeline, &bind_groups, workgroups);
                Self::response(id, result)
            }
            WireCommandRuntimeRequest::GpuSubmitBatch { id, batch } => {
                let result = self.gpu_submit_batch(batch);
                Self::response(id, result)
            }
            WireCommandRuntimeRequest::GpuDropHandles { id, handles } => {
                let result = self.gpu_drop_handles(&handles);
                Self::response(id, result)
            }
        }
    }

    fn response(
        id: u64,
        result: Result<WireCommandRuntimeValue, String>,
    ) -> WireCommandRuntimeResponse {
        WireCommandRuntimeResponse {
            wire_version: RUNTIME_WIRE_VERSION,
            id,
            result,
        }
    }

    fn ensure_requirement(&self, requirement: CommandRuntimeRequirements) -> Result<(), String> {
        if self.runtime_requirements.contains(requirement) {
            Ok(())
        } else {
            Err("plugin command did not declare the required runtime input".to_string())
        }
    }

    fn open_trace_geometry_stream(&mut self) -> Result<WireCommandRuntimeValue, String> {
        self.ensure_requirement(CommandRuntimeRequirements::TRACE_GEOMETRY_STREAM)?;
        let stream_id = self.next_stream_id;
        self.next_stream_id += 1;
        let start = std::time::Instant::now();
        let mut chunks = VecDeque::new();
        let mut chunk_count = 0usize;
        let mut primitive_count = 0usize;

        self.viewer.for_each_trace_geometry_chunk(
            &patinae_render::GeometryExportOptions::default(),
            &mut |chunk| {
                primitive_count += chunk.primitive_count();
                let wire_chunk = wire::trace_geometry_chunk_to_wire(chunk);
                let encoded = wire::encode(&wire_chunk)?;
                ensure_runtime_payload_within_limit(
                    "command trace geometry stream",
                    encoded.len(),
                    "trace geometry chunk",
                )?;
                chunks.push_back(encoded);
                chunk_count += 1;
                Ok(())
            },
        )?;

        if rt_profile_enabled() {
            log::info!(
                "patinae.rt_profile host.trace_export_ms={} chunks={} primitives={}",
                start.elapsed().as_millis(),
                chunk_count,
                primitive_count
            );
        }

        self.streams.push(HostTraceGeometryStream {
            id: stream_id,
            chunks,
        });
        Ok(WireCommandRuntimeValue::TraceGeometryOpened(
            WireTraceGeometryOpened { stream_id },
        ))
    }

    fn read_trace_geometry_stream(
        &mut self,
        stream_id: u64,
    ) -> Result<WireCommandRuntimeValue, String> {
        self.ensure_requirement(CommandRuntimeRequirements::TRACE_GEOMETRY_STREAM)?;
        let stream = self
            .streams
            .iter_mut()
            .find(|stream| stream.id == stream_id)
            .ok_or_else(|| format!("trace geometry stream {stream_id} was not open"))?;
        Ok(WireCommandRuntimeValue::TraceGeometryChunkBytes(
            stream.chunks.pop_front(),
        ))
    }

    fn close_trace_geometry_stream(
        &mut self,
        stream_id: u64,
    ) -> Result<WireCommandRuntimeValue, String> {
        self.ensure_requirement(CommandRuntimeRequirements::TRACE_GEOMETRY_STREAM)?;
        let index = self
            .streams
            .iter()
            .position(|stream| stream.id == stream_id)
            .ok_or_else(|| format!("trace geometry stream {stream_id} was not open"))?;
        self.streams.swap_remove(index);
        Ok(WireCommandRuntimeValue::TraceGeometryClosed)
    }

    fn open_render_artifact_snapshot(&mut self) -> Result<WireCommandRuntimeValue, String> {
        self.ensure_requirement(CommandRuntimeRequirements::RENDER_ARTIFACTS)?;
        let snapshot_id = self.next_snapshot_id;
        self.next_snapshot_id += 1;

        let limits = gpu_device_limits_from_viewer(self.viewer)?;
        let scene_generation = self
            .viewer
            .objects()
            .generation()
            .wrapping_add(self.viewer.selections().generation().rotate_left(1))
            .wrapping_add(self.viewer.movie().generation().rotate_left(2));
        let mut captured = None;
        self.viewer.visit_render_artifacts(&mut |snapshot| {
            let buffers = snapshot
                .buffers
                .iter()
                .map(|buffer| HostRenderArtifactBuffer {
                    role: map_render_artifact_buffer_role(buffer.role),
                    buffer: buffer.buffer.clone(),
                    size: buffer.size,
                    stride: buffer.stride,
                    element_count: buffer.element_count,
                })
                .collect();
            let reps = snapshot
                .reps
                .iter()
                .copied()
                .map(|rep| HostRenderArtifactRep {
                    object_id: rep.object_id,
                    rep_kind: map_render_artifact_rep_kind(rep.rep_kind),
                    topology: map_render_artifact_topology(rep.topology),
                    geometry_buffer_index: rep.geometry_buffer_index,
                    count_buffer_index: rep.count_buffer_index,
                    indirect_buffer_index: rep.indirect_buffer_index,
                    element_count: rep.element_count,
                    max_element_count: rep.max_element_count,
                    atom_offset: rep.atom_offset,
                    atom_count: rep.atom_count,
                    material_rgba: rep.material_rgba,
                    transparency: rep.transparency,
                })
                .collect();
            captured = Some(HostRenderArtifactSnapshot {
                id: snapshot_id,
                layout_version: snapshot.layout_version,
                scene_generation,
                scene_bounds_min: snapshot.scene_bounds_min,
                scene_bounds_max: snapshot.scene_bounds_max,
                cull_pass_initialized: snapshot.cull_pass_initialized,
                buffers,
                reps,
            });
            Ok(())
        })?;

        let snapshot = captured
            .ok_or_else(|| "renderer did not provide a render artifact snapshot".to_string())?;
        let mut buffer_descriptors = Vec::with_capacity(snapshot.buffers.len());
        for buffer in &snapshot.buffers {
            let handle = self
                .gpu_handles
                .insert_buffer(buffer.buffer.clone(), buffer.size);
            buffer_descriptors.push(RenderArtifactBufferDescriptor {
                handle,
                role: buffer.role,
                size: buffer.size,
                stride: buffer.stride,
                element_count: buffer.element_count,
            });
        }

        let mut rep_descriptors = Vec::with_capacity(snapshot.reps.len());
        for rep in &snapshot.reps {
            let geometry = buffer_descriptors
                .get(rep.geometry_buffer_index)
                .ok_or_else(|| {
                    "render artifact rep referenced missing geometry buffer".to_string()
                })?
                .handle;
            let count = rep
                .count_buffer_index
                .map(|index| {
                    buffer_descriptors
                        .get(index)
                        .map(|descriptor| descriptor.handle)
                        .ok_or_else(|| {
                            "render artifact rep referenced missing count buffer".to_string()
                        })
                })
                .transpose()?;
            let indirect = rep
                .indirect_buffer_index
                .map(|index| {
                    buffer_descriptors
                        .get(index)
                        .map(|descriptor| descriptor.handle)
                        .ok_or_else(|| {
                            "render artifact rep referenced missing indirect buffer".to_string()
                        })
                })
                .transpose()?;
            rep_descriptors.push(RenderArtifactRepDescriptor {
                object_id: rep.object_id,
                rep_kind: rep.rep_kind,
                topology: rep.topology,
                geometry,
                count,
                indirect,
                element_count: rep.element_count,
                max_element_count: rep.max_element_count,
                atom_offset: rep.atom_offset,
                atom_count: rep.atom_count,
                material_rgba: rep.material_rgba,
                transparency: rep.transparency,
            });
        }

        let descriptor = RenderArtifactSnapshotDescriptor {
            snapshot_id,
            layout_version: snapshot.layout_version,
            scene_generation: snapshot.scene_generation,
            scene_bounds_min: snapshot.scene_bounds_min,
            scene_bounds_max: snapshot.scene_bounds_max,
            cull_pass_initialized: snapshot.cull_pass_initialized,
            device_limits: limits,
            buffers: buffer_descriptors,
            reps: rep_descriptors,
        };
        self.snapshots.push(snapshot);
        Ok(WireCommandRuntimeValue::RenderArtifactSnapshotOpened(
            descriptor,
        ))
    }

    fn close_render_artifact_snapshot(
        &mut self,
        snapshot_id: u64,
    ) -> Result<WireCommandRuntimeValue, String> {
        self.ensure_requirement(CommandRuntimeRequirements::RENDER_ARTIFACTS)?;
        let index = self
            .snapshots
            .iter()
            .position(|snapshot| snapshot.id == snapshot_id)
            .ok_or_else(|| format!("render artifact snapshot {snapshot_id} was not open"))?;
        self.snapshots.swap_remove(index);
        Ok(WireCommandRuntimeValue::RenderArtifactSnapshotClosed)
    }

    fn gpu_device_limits(&mut self) -> Result<WireCommandRuntimeValue, String> {
        self.ensure_requirement(CommandRuntimeRequirements::GPU_COMMANDS)?;
        Ok(WireCommandRuntimeValue::GpuDeviceLimits(
            gpu_device_limits_from_viewer(self.viewer)?,
        ))
    }

    fn gpu_create_buffer(
        &mut self,
        descriptor: GpuBufferDescriptor,
        initial_data: Option<Vec<u8>>,
    ) -> Result<WireCommandRuntimeValue, String> {
        self.ensure_requirement(CommandRuntimeRequirements::GPU_COMMANDS)?;
        let device = self
            .viewer
            .gpu_device()
            .ok_or_else(|| "GPU command runtime requires an active device".to_string())?;
        let usage = wgpu_buffer_usage(descriptor.usage)?;
        let buffer = if let Some(data) = initial_data {
            if data.len() as u64 > descriptor.size {
                return Err("initial GPU buffer data exceeds descriptor size".to_string());
            }
            device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: descriptor.label.as_deref(),
                contents: &data,
                usage,
            })
        } else {
            device.create_buffer(&wgpu::BufferDescriptor {
                label: descriptor.label.as_deref(),
                size: descriptor.size,
                usage,
                mapped_at_creation: false,
            })
        };
        let handle = self.gpu_handles.insert_buffer(buffer, descriptor.size);
        Ok(WireCommandRuntimeValue::GpuHandle(handle))
    }

    fn gpu_write_buffer(
        &mut self,
        buffer: GpuHandle,
        offset: u64,
        data: &[u8],
    ) -> Result<WireCommandRuntimeValue, String> {
        self.ensure_requirement(CommandRuntimeRequirements::GPU_COMMANDS)?;
        let queue = self
            .viewer
            .gpu_queue()
            .ok_or_else(|| "GPU command runtime requires an active queue".to_string())?;
        let target = self.gpu_handles.buffer(buffer)?;
        validate_buffer_range(target.size, offset, data.len() as u64)?;
        queue.write_buffer(&target.buffer, offset, data);
        Ok(WireCommandRuntimeValue::GpuOk)
    }

    fn gpu_copy_buffer_to_buffer(
        &mut self,
        source: GpuHandle,
        source_offset: u64,
        destination: GpuHandle,
        destination_offset: u64,
        size: u64,
    ) -> Result<WireCommandRuntimeValue, String> {
        self.ensure_requirement(CommandRuntimeRequirements::GPU_COMMANDS)?;
        let device = self
            .viewer
            .gpu_device()
            .ok_or_else(|| "GPU command runtime requires an active device".to_string())?;
        let queue = self
            .viewer
            .gpu_queue()
            .ok_or_else(|| "GPU command runtime requires an active queue".to_string())?;
        let source = self.gpu_handles.buffer(source)?;
        let destination = self.gpu_handles.buffer(destination)?;
        validate_buffer_range(source.size, source_offset, size)?;
        validate_buffer_range(destination.size, destination_offset, size)?;
        let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("patinae.plugin.gpu.copy_buffer_to_buffer"),
        });
        encoder.copy_buffer_to_buffer(
            &source.buffer,
            source_offset,
            &destination.buffer,
            destination_offset,
            size,
        );
        queue.submit(std::iter::once(encoder.finish()));
        Ok(WireCommandRuntimeValue::GpuOk)
    }

    fn gpu_read_buffer(
        &mut self,
        buffer: GpuHandle,
        offset: u64,
        size: u64,
    ) -> Result<WireCommandRuntimeValue, String> {
        self.ensure_requirement(CommandRuntimeRequirements::GPU_COMMANDS)?;
        let device = self
            .viewer
            .gpu_device()
            .ok_or_else(|| "GPU command runtime requires an active device".to_string())?;
        let queue = self
            .viewer
            .gpu_queue()
            .ok_or_else(|| "GPU command runtime requires an active queue".to_string())?;
        let source = self.gpu_handles.buffer(buffer)?;
        validate_buffer_range(source.size, offset, size)?;
        let staging = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("patinae.plugin.gpu.readback"),
            size: size.max(4),
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });
        let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("patinae.plugin.gpu.readback.copy"),
        });
        encoder.copy_buffer_to_buffer(&source.buffer, offset, &staging, 0, size);
        queue.submit(std::iter::once(encoder.finish()));

        let bytes = map_readback_buffer(device, staging, size)?;
        Ok(WireCommandRuntimeValue::GpuBytes(bytes))
    }

    fn gpu_create_shader_module(
        &mut self,
        descriptor: GpuShaderModuleDescriptor,
    ) -> Result<WireCommandRuntimeValue, String> {
        self.ensure_requirement(CommandRuntimeRequirements::GPU_COMMANDS)?;
        let fingerprint = shader_module_fingerprint(&descriptor);
        let device = self
            .viewer
            .gpu_device()
            .ok_or_else(|| "GPU command runtime requires an active device".to_string())?;
        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: descriptor.label.as_deref(),
            source: wgpu::ShaderSource::Wgsl(descriptor.wgsl.into()),
        });
        let handle = self.gpu_handles.insert_shader_module(module, fingerprint);
        Ok(WireCommandRuntimeValue::GpuHandle(handle))
    }

    fn gpu_create_bind_group_layout(
        &mut self,
        descriptor: GpuBindGroupLayoutDescriptor,
    ) -> Result<WireCommandRuntimeValue, String> {
        self.ensure_requirement(CommandRuntimeRequirements::GPU_COMMANDS)?;
        let fingerprint = bind_group_layout_fingerprint(&descriptor);
        let device = self
            .viewer
            .gpu_device()
            .ok_or_else(|| "GPU command runtime requires an active device".to_string())?;
        let entries = descriptor
            .entries
            .iter()
            .copied()
            .map(wgpu_bind_group_layout_entry)
            .collect::<Result<Vec<_>, _>>()?;
        let layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: descriptor.label.as_deref(),
            entries: &entries,
        });
        let handle = self
            .gpu_handles
            .insert_bind_group_layout(layout, fingerprint);
        Ok(WireCommandRuntimeValue::GpuHandle(handle))
    }

    fn gpu_create_pipeline_layout(
        &mut self,
        descriptor: GpuPipelineLayoutDescriptor,
    ) -> Result<WireCommandRuntimeValue, String> {
        self.ensure_requirement(CommandRuntimeRequirements::GPU_COMMANDS)?;
        let fingerprint = self.pipeline_layout_fingerprint(&descriptor)?;
        let device = self
            .viewer
            .gpu_device()
            .ok_or_else(|| "GPU command runtime requires an active device".to_string())?;
        let layout = {
            let layouts = descriptor
                .bind_group_layouts
                .iter()
                .map(|handle| self.gpu_handles.bind_group_layout(*handle))
                .collect::<Result<Vec<_>, _>>()?;
            device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                label: descriptor.label.as_deref(),
                bind_group_layouts: &layouts,
                immediate_size: 0,
            })
        };
        let handle = self.gpu_handles.insert_pipeline_layout(layout, fingerprint);
        Ok(WireCommandRuntimeValue::GpuHandle(handle))
    }

    fn gpu_create_compute_pipeline(
        &mut self,
        descriptor: GpuComputePipelineDescriptor,
    ) -> Result<WireCommandRuntimeValue, String> {
        self.ensure_requirement(CommandRuntimeRequirements::GPU_COMMANDS)?;
        let fingerprint = self.compute_pipeline_fingerprint(&descriptor)?;
        let device = self
            .viewer
            .gpu_device()
            .ok_or_else(|| "GPU command runtime requires an active device".to_string())?;
        let pipeline = {
            let layout = self.gpu_handles.pipeline_layout(descriptor.layout)?;
            let module = self.gpu_handles.shader_module(descriptor.module)?;
            device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: descriptor.label.as_deref(),
                layout: Some(layout),
                module,
                entry_point: Some(&descriptor.entry_point),
                compilation_options: wgpu::PipelineCompilationOptions::default(),
                cache: None,
            })
        };
        let handle = self
            .gpu_handles
            .insert_compute_pipeline(pipeline, fingerprint);
        Ok(WireCommandRuntimeValue::GpuHandle(handle))
    }

    fn gpu_create_cached_shader_module(
        &mut self,
        descriptor: GpuShaderModuleDescriptor,
    ) -> Result<WireCommandRuntimeValue, String> {
        self.ensure_requirement(CommandRuntimeRequirements::GPU_COMMANDS)?;
        let fingerprint = shader_module_fingerprint(&descriptor);
        let device = self
            .viewer
            .gpu_device()
            .ok_or_else(|| "GPU command runtime requires an active device".to_string())?;
        let start = std::time::Instant::now();
        let mut cache = self.lock_gpu_cache()?;
        let key = cache.resource_key(device, fingerprint);
        let (module, status) = if let Some(object) = cache.cached_object(&key) {
            let module = object.into_shader_module()?;
            cache.record_hit();
            (module, GpuCacheStatus::Hit)
        } else {
            let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
                label: descriptor.label.as_deref(),
                source: wgpu::ShaderSource::Wgsl(descriptor.wgsl.clone().into()),
            });
            cache.insert(key, GpuCachedObject::ShaderModule(module.clone()));
            (module, GpuCacheStatus::Miss)
        };
        drop(cache);
        let handle = self.gpu_handles.insert_shader_module(module, fingerprint);
        log_gpu_cache_event("shader_module", status, start.elapsed().as_millis());
        Ok(WireCommandRuntimeValue::GpuCachedHandle(GpuCachedHandle {
            handle,
            status,
        }))
    }

    fn gpu_create_cached_bind_group_layout(
        &mut self,
        descriptor: GpuBindGroupLayoutDescriptor,
    ) -> Result<WireCommandRuntimeValue, String> {
        self.ensure_requirement(CommandRuntimeRequirements::GPU_COMMANDS)?;
        let fingerprint = bind_group_layout_fingerprint(&descriptor);
        let entries = descriptor
            .entries
            .iter()
            .copied()
            .map(wgpu_bind_group_layout_entry)
            .collect::<Result<Vec<_>, _>>()?;
        let device = self
            .viewer
            .gpu_device()
            .ok_or_else(|| "GPU command runtime requires an active device".to_string())?;
        let start = std::time::Instant::now();
        let mut cache = self.lock_gpu_cache()?;
        let key = cache.resource_key(device, fingerprint);
        let (layout, status) = if let Some(object) = cache.cached_object(&key) {
            let layout = object.into_bind_group_layout()?;
            cache.record_hit();
            (layout, GpuCacheStatus::Hit)
        } else {
            let layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: descriptor.label.as_deref(),
                entries: &entries,
            });
            cache.insert(key, GpuCachedObject::BindGroupLayout(layout.clone()));
            (layout, GpuCacheStatus::Miss)
        };
        drop(cache);
        let handle = self
            .gpu_handles
            .insert_bind_group_layout(layout, fingerprint);
        log_gpu_cache_event("bind_group_layout", status, start.elapsed().as_millis());
        Ok(WireCommandRuntimeValue::GpuCachedHandle(GpuCachedHandle {
            handle,
            status,
        }))
    }

    fn gpu_create_cached_pipeline_layout(
        &mut self,
        descriptor: GpuPipelineLayoutDescriptor,
    ) -> Result<WireCommandRuntimeValue, String> {
        self.ensure_requirement(CommandRuntimeRequirements::GPU_COMMANDS)?;
        let fingerprint = self.pipeline_layout_fingerprint(&descriptor)?;
        let device = self
            .viewer
            .gpu_device()
            .ok_or_else(|| "GPU command runtime requires an active device".to_string())?;
        let start = std::time::Instant::now();
        let mut cache = self.lock_gpu_cache()?;
        let key = cache.resource_key(device, fingerprint);
        let (layout, status) = if let Some(object) = cache.cached_object(&key) {
            let layout = object.into_pipeline_layout()?;
            cache.record_hit();
            (layout, GpuCacheStatus::Hit)
        } else {
            let layout = {
                let layouts = descriptor
                    .bind_group_layouts
                    .iter()
                    .map(|handle| self.gpu_handles.bind_group_layout(*handle))
                    .collect::<Result<Vec<_>, _>>()?;
                device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                    label: descriptor.label.as_deref(),
                    bind_group_layouts: &layouts,
                    immediate_size: 0,
                })
            };
            cache.insert(key, GpuCachedObject::PipelineLayout(layout.clone()));
            (layout, GpuCacheStatus::Miss)
        };
        drop(cache);
        let handle = self.gpu_handles.insert_pipeline_layout(layout, fingerprint);
        log_gpu_cache_event("pipeline_layout", status, start.elapsed().as_millis());
        Ok(WireCommandRuntimeValue::GpuCachedHandle(GpuCachedHandle {
            handle,
            status,
        }))
    }

    fn gpu_create_cached_compute_pipeline(
        &mut self,
        descriptor: GpuComputePipelineDescriptor,
    ) -> Result<WireCommandRuntimeValue, String> {
        self.ensure_requirement(CommandRuntimeRequirements::GPU_COMMANDS)?;
        let fingerprint = self.compute_pipeline_fingerprint(&descriptor)?;
        let device = self
            .viewer
            .gpu_device()
            .ok_or_else(|| "GPU command runtime requires an active device".to_string())?;
        let start = std::time::Instant::now();
        let mut cache = self.lock_gpu_cache()?;
        let key = cache.resource_key(device, fingerprint);
        let (pipeline, status) = if let Some(object) = cache.cached_object(&key) {
            let pipeline = object.into_compute_pipeline()?;
            cache.record_hit();
            (pipeline, GpuCacheStatus::Hit)
        } else {
            let pipeline = {
                let layout = self.gpu_handles.pipeline_layout(descriptor.layout)?;
                let module = self.gpu_handles.shader_module(descriptor.module)?;
                device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                    label: descriptor.label.as_deref(),
                    layout: Some(layout),
                    module,
                    entry_point: Some(&descriptor.entry_point),
                    compilation_options: wgpu::PipelineCompilationOptions::default(),
                    cache: None,
                })
            };
            cache.insert(key, GpuCachedObject::ComputePipeline(pipeline.clone()));
            (pipeline, GpuCacheStatus::Miss)
        };
        drop(cache);
        let handle = self
            .gpu_handles
            .insert_compute_pipeline(pipeline, fingerprint);
        log_gpu_cache_event("compute_pipeline", status, start.elapsed().as_millis());
        Ok(WireCommandRuntimeValue::GpuCachedHandle(GpuCachedHandle {
            handle,
            status,
        }))
    }

    fn gpu_cache_stats(&mut self) -> Result<WireCommandRuntimeValue, String> {
        self.ensure_requirement(CommandRuntimeRequirements::GPU_COMMANDS)?;
        let cache = self.lock_gpu_cache()?;
        Ok(WireCommandRuntimeValue::GpuCacheStats(cache.stats()))
    }

    fn gpu_drop_plugin_cache(&mut self) -> Result<WireCommandRuntimeValue, String> {
        self.ensure_requirement(CommandRuntimeRequirements::GPU_COMMANDS)?;
        self.lock_gpu_cache()?.clear();
        Ok(WireCommandRuntimeValue::GpuOk)
    }

    fn gpu_create_bind_group(
        &mut self,
        descriptor: GpuBindGroupDescriptor,
    ) -> Result<WireCommandRuntimeValue, String> {
        self.ensure_requirement(CommandRuntimeRequirements::GPU_COMMANDS)?;
        let device = self
            .viewer
            .gpu_device()
            .ok_or_else(|| "GPU command runtime requires an active device".to_string())?;
        let bind_group = {
            let layout = self.gpu_handles.bind_group_layout(descriptor.layout)?;
            let mut entries = Vec::with_capacity(descriptor.entries.len());
            for entry in descriptor.entries {
                entries.push(wgpu_bind_group_entry(&self.gpu_handles, entry)?);
            }
            device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: descriptor.label.as_deref(),
                layout,
                entries: &entries,
            })
        };
        let handle = self.gpu_handles.insert_bind_group(bind_group);
        Ok(WireCommandRuntimeValue::GpuHandle(handle))
    }

    fn gpu_dispatch_compute(
        &mut self,
        pipeline: GpuHandle,
        bind_groups: &[GpuHandle],
        workgroups: [u32; 3],
    ) -> Result<WireCommandRuntimeValue, String> {
        self.ensure_requirement(CommandRuntimeRequirements::GPU_COMMANDS)?;
        let device = self
            .viewer
            .gpu_device()
            .ok_or_else(|| "GPU command runtime requires an active device".to_string())?;
        let queue = self
            .viewer
            .gpu_queue()
            .ok_or_else(|| "GPU command runtime requires an active queue".to_string())?;
        validate_dispatch_workgroups(
            workgroups,
            device.limits().max_compute_workgroups_per_dimension,
        )?;
        let pipeline = self.gpu_handles.compute_pipeline(pipeline)?;
        let bind_groups = bind_groups
            .iter()
            .map(|handle| self.gpu_handles.bind_group(*handle))
            .collect::<Result<Vec<_>, _>>()?;
        let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("patinae.plugin.gpu.dispatch_compute"),
        });
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("patinae.plugin.gpu.dispatch_compute.pass"),
                timestamp_writes: None,
            });
            pass.set_pipeline(pipeline);
            for (index, bind_group) in bind_groups.iter().enumerate() {
                pass.set_bind_group(index as u32, Some(*bind_group), &[]);
            }
            pass.dispatch_workgroups(workgroups[0], workgroups[1], workgroups[2]);
        }
        queue.submit(std::iter::once(encoder.finish()));
        Ok(WireCommandRuntimeValue::GpuOk)
    }

    fn gpu_submit_batch(
        &mut self,
        batch: GpuSubmitBatch,
    ) -> Result<WireCommandRuntimeValue, String> {
        self.ensure_requirement(CommandRuntimeRequirements::GPU_COMMANDS)?;
        let device = self
            .viewer
            .gpu_device()
            .ok_or_else(|| "GPU command runtime requires an active device".to_string())?;
        let queue = self
            .viewer
            .gpu_queue()
            .ok_or_else(|| "GPU command runtime requires an active queue".to_string())?;
        let record_start = std::time::Instant::now();
        let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: batch.label.as_deref().or(Some("patinae.plugin.gpu.batch")),
        });
        let mut readbacks = Vec::new();
        let mut upload_staging = Vec::new();

        for command in batch.commands {
            match command {
                GpuBatchCommand::WriteBuffer {
                    buffer,
                    offset,
                    data,
                } => {
                    let target = self.gpu_handles.buffer(buffer)?;
                    validate_buffer_range(target.size, offset, data.len() as u64)?;
                    let staging = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                        label: Some("patinae.plugin.gpu.batch.upload"),
                        contents: &data,
                        usage: wgpu::BufferUsages::COPY_SRC,
                    });
                    encoder.copy_buffer_to_buffer(
                        &staging,
                        0,
                        &target.buffer,
                        offset,
                        data.len() as u64,
                    );
                    upload_staging.push(staging);
                }
                GpuBatchCommand::DispatchCompute {
                    pipeline,
                    bind_groups,
                    workgroups,
                } => {
                    validate_dispatch_workgroups(
                        workgroups,
                        device.limits().max_compute_workgroups_per_dimension,
                    )?;
                    let pipeline = self.gpu_handles.compute_pipeline(pipeline)?;
                    let bind_groups = bind_groups
                        .iter()
                        .map(|handle| self.gpu_handles.bind_group(*handle))
                        .collect::<Result<Vec<_>, _>>()?;
                    let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                        label: Some("patinae.plugin.gpu.batch.compute"),
                        timestamp_writes: None,
                    });
                    pass.set_pipeline(pipeline);
                    for (index, bind_group) in bind_groups.iter().enumerate() {
                        pass.set_bind_group(index as u32, Some(*bind_group), &[]);
                    }
                    pass.dispatch_workgroups(workgroups[0], workgroups[1], workgroups[2]);
                }
                GpuBatchCommand::CopyBufferToBuffer {
                    source,
                    source_offset,
                    destination,
                    destination_offset,
                    size,
                } => {
                    let source = self.gpu_handles.buffer(source)?;
                    let destination = self.gpu_handles.buffer(destination)?;
                    validate_buffer_range(source.size, source_offset, size)?;
                    validate_buffer_range(destination.size, destination_offset, size)?;
                    encoder.copy_buffer_to_buffer(
                        &source.buffer,
                        source_offset,
                        &destination.buffer,
                        destination_offset,
                        size,
                    );
                }
                GpuBatchCommand::ReadBuffer {
                    buffer,
                    offset,
                    size,
                } => {
                    let source = self.gpu_handles.buffer(buffer)?;
                    validate_buffer_range(source.size, offset, size)?;
                    let staging = device.create_buffer(&wgpu::BufferDescriptor {
                        label: Some("patinae.plugin.gpu.batch.readback"),
                        size: size.max(4),
                        usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
                        mapped_at_creation: false,
                    });
                    encoder.copy_buffer_to_buffer(&source.buffer, offset, &staging, 0, size);
                    readbacks.push((staging, size));
                }
            }
        }

        let record_ms = record_start.elapsed().as_millis();
        let submit_start = std::time::Instant::now();
        queue.submit(std::iter::once(encoder.finish()));
        let submit_wait_ms = submit_start.elapsed().as_millis();
        let read_start = std::time::Instant::now();
        let mut bytes = Vec::with_capacity(readbacks.len());
        for (buffer, size) in readbacks {
            bytes.push(map_readback_buffer(device, buffer, size)?);
        }
        if rt_profile_enabled() {
            log::info!(
                "patinae.rt_profile host.batch_record_ms={} host.batch_submit_wait_ms={} host.batch_readback_ms={} readbacks={}",
                record_ms,
                submit_wait_ms,
                read_start.elapsed().as_millis(),
                bytes.len()
            );
        }
        Ok(WireCommandRuntimeValue::GpuBatchResult(GpuBatchResult {
            readbacks: bytes,
        }))
    }

    fn gpu_drop_handles(
        &mut self,
        handles: &[GpuHandle],
    ) -> Result<WireCommandRuntimeValue, String> {
        self.ensure_requirement(CommandRuntimeRequirements::GPU_COMMANDS)?;
        for handle in handles {
            self.gpu_handles.drop_handle(*handle)?;
        }
        Ok(WireCommandRuntimeValue::GpuOk)
    }

    fn lock_gpu_cache(&self) -> Result<std::sync::MutexGuard<'_, GpuPluginCache>, String> {
        self.gpu_cache
            .lock()
            .map_err(|_| "GPU plugin cache was poisoned".to_string())
    }

    fn pipeline_layout_fingerprint(
        &self,
        descriptor: &GpuPipelineLayoutDescriptor,
    ) -> Result<GpuResourceFingerprint, String> {
        let layouts = descriptor
            .bind_group_layouts
            .iter()
            .map(|handle| {
                self.gpu_handles
                    .fingerprint(*handle, GpuHandleKind::BindGroupLayout)
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(pipeline_layout_fingerprint(&layouts))
    }

    fn compute_pipeline_fingerprint(
        &self,
        descriptor: &GpuComputePipelineDescriptor,
    ) -> Result<GpuResourceFingerprint, String> {
        let layout = self
            .gpu_handles
            .fingerprint(descriptor.layout, GpuHandleKind::PipelineLayout)?;
        let module = self
            .gpu_handles
            .fingerprint(descriptor.module, GpuHandleKind::ShaderModule)?;
        Ok(compute_pipeline_fingerprint(
            layout,
            module,
            &descriptor.entry_point,
        ))
    }
}

struct HostTraceGeometryStream {
    id: u64,
    chunks: VecDeque<Vec<u8>>,
}

struct HostRenderArtifactSnapshot {
    id: u64,
    layout_version: u32,
    scene_generation: u64,
    scene_bounds_min: [f32; 3],
    scene_bounds_max: [f32; 3],
    cull_pass_initialized: bool,
    buffers: Vec<HostRenderArtifactBuffer>,
    reps: Vec<HostRenderArtifactRep>,
}

struct HostRenderArtifactBuffer {
    role: RenderArtifactBufferRole,
    buffer: wgpu::Buffer,
    size: u64,
    stride: u64,
    element_count: u64,
}

struct HostRenderArtifactRep {
    object_id: u32,
    rep_kind: RenderArtifactRepKind,
    topology: RenderArtifactPrimitiveTopology,
    geometry_buffer_index: usize,
    count_buffer_index: Option<usize>,
    indirect_buffer_index: Option<usize>,
    element_count: u64,
    max_element_count: u64,
    atom_offset: u32,
    atom_count: u32,
    material_rgba: [f32; 4],
    transparency: f32,
}

#[derive(Clone)]
struct HostGpuBuffer {
    buffer: wgpu::Buffer,
    size: u64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct GpuResourceFingerprint {
    kind: GpuHandleKind,
    hash: u64,
}

#[derive(Clone)]
enum GpuHandleObject {
    Buffer(HostGpuBuffer),
    ShaderModule(wgpu::ShaderModule),
    BindGroupLayout(wgpu::BindGroupLayout),
    PipelineLayout(wgpu::PipelineLayout),
    BindGroup(wgpu::BindGroup),
    ComputePipeline(wgpu::ComputePipeline),
}

struct GpuHandleEntry {
    generation: u64,
    kind: GpuHandleKind,
    fingerprint: Option<GpuResourceFingerprint>,
    object: GpuHandleObject,
}

struct GpuHandleRegistry {
    next_id: u64,
    generation: u64,
    entries: HashMap<u64, GpuHandleEntry>,
}

impl GpuHandleRegistry {
    fn new() -> Self {
        Self {
            next_id: 1,
            generation: 1,
            entries: HashMap::new(),
        }
    }

    fn insert_buffer(&mut self, buffer: wgpu::Buffer, size: u64) -> GpuHandle {
        self.insert(
            GpuHandleKind::Buffer,
            GpuHandleObject::Buffer(HostGpuBuffer { buffer, size }),
            None,
        )
    }

    fn insert_shader_module(
        &mut self,
        module: wgpu::ShaderModule,
        fingerprint: GpuResourceFingerprint,
    ) -> GpuHandle {
        self.insert(
            GpuHandleKind::ShaderModule,
            GpuHandleObject::ShaderModule(module),
            Some(fingerprint),
        )
    }

    fn insert_bind_group_layout(
        &mut self,
        layout: wgpu::BindGroupLayout,
        fingerprint: GpuResourceFingerprint,
    ) -> GpuHandle {
        self.insert(
            GpuHandleKind::BindGroupLayout,
            GpuHandleObject::BindGroupLayout(layout),
            Some(fingerprint),
        )
    }

    fn insert_pipeline_layout(
        &mut self,
        layout: wgpu::PipelineLayout,
        fingerprint: GpuResourceFingerprint,
    ) -> GpuHandle {
        self.insert(
            GpuHandleKind::PipelineLayout,
            GpuHandleObject::PipelineLayout(layout),
            Some(fingerprint),
        )
    }

    fn insert_bind_group(&mut self, bind_group: wgpu::BindGroup) -> GpuHandle {
        self.insert(
            GpuHandleKind::BindGroup,
            GpuHandleObject::BindGroup(bind_group),
            None,
        )
    }

    fn insert_compute_pipeline(
        &mut self,
        pipeline: wgpu::ComputePipeline,
        fingerprint: GpuResourceFingerprint,
    ) -> GpuHandle {
        self.insert(
            GpuHandleKind::ComputePipeline,
            GpuHandleObject::ComputePipeline(pipeline),
            Some(fingerprint),
        )
    }

    fn insert(
        &mut self,
        kind: GpuHandleKind,
        object: GpuHandleObject,
        fingerprint: Option<GpuResourceFingerprint>,
    ) -> GpuHandle {
        let id = self.next_id;
        self.next_id += 1;
        self.entries.insert(
            id,
            GpuHandleEntry {
                generation: self.generation,
                kind,
                fingerprint,
                object,
            },
        );
        GpuHandle {
            id,
            kind,
            generation: self.generation,
        }
    }

    fn drop_handle(&mut self, handle: GpuHandle) -> Result<(), String> {
        self.validate(handle)?;
        self.entries.remove(&handle.id);
        Ok(())
    }

    fn buffer(&self, handle: GpuHandle) -> Result<&HostGpuBuffer, String> {
        match &self.entry(handle, GpuHandleKind::Buffer)?.object {
            GpuHandleObject::Buffer(buffer) => Ok(buffer),
            _ => Err("GPU handle kind mismatch for buffer".to_string()),
        }
    }

    fn shader_module(&self, handle: GpuHandle) -> Result<&wgpu::ShaderModule, String> {
        match &self.entry(handle, GpuHandleKind::ShaderModule)?.object {
            GpuHandleObject::ShaderModule(module) => Ok(module),
            _ => Err("GPU handle kind mismatch for shader module".to_string()),
        }
    }

    fn bind_group_layout(&self, handle: GpuHandle) -> Result<&wgpu::BindGroupLayout, String> {
        match &self.entry(handle, GpuHandleKind::BindGroupLayout)?.object {
            GpuHandleObject::BindGroupLayout(layout) => Ok(layout),
            _ => Err("GPU handle kind mismatch for bind-group layout".to_string()),
        }
    }

    fn pipeline_layout(&self, handle: GpuHandle) -> Result<&wgpu::PipelineLayout, String> {
        match &self.entry(handle, GpuHandleKind::PipelineLayout)?.object {
            GpuHandleObject::PipelineLayout(layout) => Ok(layout),
            _ => Err("GPU handle kind mismatch for pipeline layout".to_string()),
        }
    }

    fn bind_group(&self, handle: GpuHandle) -> Result<&wgpu::BindGroup, String> {
        match &self.entry(handle, GpuHandleKind::BindGroup)?.object {
            GpuHandleObject::BindGroup(bind_group) => Ok(bind_group),
            _ => Err("GPU handle kind mismatch for bind group".to_string()),
        }
    }

    fn compute_pipeline(&self, handle: GpuHandle) -> Result<&wgpu::ComputePipeline, String> {
        match &self.entry(handle, GpuHandleKind::ComputePipeline)?.object {
            GpuHandleObject::ComputePipeline(pipeline) => Ok(pipeline),
            _ => Err("GPU handle kind mismatch for compute pipeline".to_string()),
        }
    }

    fn fingerprint(
        &self,
        handle: GpuHandle,
        expected: GpuHandleKind,
    ) -> Result<GpuResourceFingerprint, String> {
        let entry = self.entry(handle, expected)?;
        entry.fingerprint.ok_or_else(|| {
            format!(
                "GPU handle {} of kind {:?} cannot be used in cached descriptor",
                handle.id, expected
            )
        })
    }

    fn entry(&self, handle: GpuHandle, expected: GpuHandleKind) -> Result<&GpuHandleEntry, String> {
        self.validate(handle)?;
        let entry = self
            .entries
            .get(&handle.id)
            .ok_or_else(|| format!("GPU handle {} was not found", handle.id))?;
        if entry.kind != expected {
            return Err(format!(
                "GPU handle {} has kind {:?}, expected {:?}",
                handle.id, entry.kind, expected
            ));
        }
        Ok(entry)
    }

    fn validate(&self, handle: GpuHandle) -> Result<(), String> {
        let entry = self
            .entries
            .get(&handle.id)
            .ok_or_else(|| format!("GPU handle {} was not found", handle.id))?;
        if entry.generation != handle.generation {
            return Err(format!(
                "GPU handle {} generation mismatch: got {}, expected {}",
                handle.id, handle.generation, entry.generation
            ));
        }
        if entry.kind != handle.kind {
            return Err(format!(
                "GPU handle {} kind mismatch: got {:?}, expected {:?}",
                handle.id, handle.kind, entry.kind
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct GpuCacheKey {
    plugin_id_hash: u64,
    device_generation: u64,
    layout_version: u32,
    fingerprint: GpuResourceFingerprint,
}

#[derive(Clone)]
enum GpuCachedObject {
    ShaderModule(wgpu::ShaderModule),
    BindGroupLayout(wgpu::BindGroupLayout),
    PipelineLayout(wgpu::PipelineLayout),
    ComputePipeline(wgpu::ComputePipeline),
}

impl GpuCachedObject {
    fn into_shader_module(self) -> Result<wgpu::ShaderModule, String> {
        match self {
            Self::ShaderModule(module) => Ok(module),
            _ => Err("GPU cache entry kind mismatch for shader module".to_string()),
        }
    }

    fn into_bind_group_layout(self) -> Result<wgpu::BindGroupLayout, String> {
        match self {
            Self::BindGroupLayout(layout) => Ok(layout),
            _ => Err("GPU cache entry kind mismatch for bind-group layout".to_string()),
        }
    }

    fn into_pipeline_layout(self) -> Result<wgpu::PipelineLayout, String> {
        match self {
            Self::PipelineLayout(layout) => Ok(layout),
            _ => Err("GPU cache entry kind mismatch for pipeline layout".to_string()),
        }
    }

    fn into_compute_pipeline(self) -> Result<wgpu::ComputePipeline, String> {
        match self {
            Self::ComputePipeline(pipeline) => Ok(pipeline),
            _ => Err("GPU cache entry kind mismatch for compute pipeline".to_string()),
        }
    }
}

struct GpuCacheEntry {
    object: GpuCachedObject,
}

struct GpuPluginCache {
    plugin_id_hash: u64,
    device_identity: Option<usize>,
    device_generation: u64,
    hits: u64,
    misses: u64,
    entries: HashMap<GpuCacheKey, GpuCacheEntry>,
}

impl GpuPluginCache {
    fn new(plugin_id: String) -> Self {
        Self {
            plugin_id_hash: semantic_hash(&plugin_id),
            device_identity: None,
            device_generation: 0,
            hits: 0,
            misses: 0,
            entries: HashMap::new(),
        }
    }

    fn resource_key(
        &mut self,
        device: &wgpu::Device,
        fingerprint: GpuResourceFingerprint,
    ) -> GpuCacheKey {
        self.resource_key_for_identity(device as *const wgpu::Device as usize, fingerprint)
    }

    fn resource_key_for_identity(
        &mut self,
        identity: usize,
        fingerprint: GpuResourceFingerprint,
    ) -> GpuCacheKey {
        let device_generation = self.ensure_device_identity(identity);
        GpuCacheKey {
            plugin_id_hash: self.plugin_id_hash,
            device_generation,
            layout_version: GPU_CACHE_LAYOUT_VERSION,
            fingerprint,
        }
    }

    fn cached_object(&self, key: &GpuCacheKey) -> Option<GpuCachedObject> {
        self.entries.get(key).map(|entry| entry.object.clone())
    }

    fn insert(&mut self, key: GpuCacheKey, object: GpuCachedObject) {
        self.entries.insert(key, GpuCacheEntry { object });
        self.record_miss();
    }

    fn record_hit(&mut self) {
        self.hits = self.hits.saturating_add(1);
    }

    fn record_miss(&mut self) {
        self.misses = self.misses.saturating_add(1);
    }

    fn stats(&self) -> GpuCacheStats {
        GpuCacheStats {
            hits: self.hits,
            misses: self.misses,
            entries: self.entries.len() as u64,
        }
    }

    fn clear(&mut self) {
        self.entries.clear();
    }

    fn ensure_device_identity(&mut self, identity: usize) -> u64 {
        if self.device_identity != Some(identity) {
            self.device_identity = Some(identity);
            self.device_generation = self.device_generation.saturating_add(1);
            self.entries.clear();
        }
        self.device_generation
    }
}

fn semantic_hash<T: Hash>(value: &T) -> u64 {
    let mut hasher = DefaultHasher::new();
    value.hash(&mut hasher);
    hasher.finish()
}

fn fingerprint_for(kind: GpuHandleKind, value: impl Hash) -> GpuResourceFingerprint {
    let mut hasher = DefaultHasher::new();
    GPU_CACHE_LAYOUT_VERSION.hash(&mut hasher);
    kind.hash(&mut hasher);
    value.hash(&mut hasher);
    GpuResourceFingerprint {
        kind,
        hash: hasher.finish(),
    }
}

fn shader_module_fingerprint(descriptor: &GpuShaderModuleDescriptor) -> GpuResourceFingerprint {
    fingerprint_for(GpuHandleKind::ShaderModule, &descriptor.wgsl)
}

fn bind_group_layout_fingerprint(
    descriptor: &GpuBindGroupLayoutDescriptor,
) -> GpuResourceFingerprint {
    fingerprint_for(GpuHandleKind::BindGroupLayout, &descriptor.entries)
}

fn pipeline_layout_fingerprint(
    bind_group_layouts: &[GpuResourceFingerprint],
) -> GpuResourceFingerprint {
    fingerprint_for(GpuHandleKind::PipelineLayout, bind_group_layouts)
}

fn compute_pipeline_fingerprint(
    layout: GpuResourceFingerprint,
    module: GpuResourceFingerprint,
    entry_point: &str,
) -> GpuResourceFingerprint {
    fingerprint_for(
        GpuHandleKind::ComputePipeline,
        (layout, module, entry_point),
    )
}

fn log_gpu_cache_event(resource: &str, status: GpuCacheStatus, elapsed_ms: u128) {
    if rt_profile_enabled() {
        log::info!(
            "patinae.rt_profile host.gpu_cache resource={} status={:?} lazy_create_ms={}",
            resource,
            status,
            elapsed_ms
        );
    }
}

fn gpu_device_limits_from_viewer(viewer: &dyn ViewerLike) -> Result<GpuDeviceLimits, String> {
    let device = viewer
        .gpu_device()
        .ok_or_else(|| "GPU command runtime requires an active device".to_string())?;
    let limits = device.limits();
    Ok(GpuDeviceLimits {
        max_buffer_size: limits.max_buffer_size,
        max_storage_buffer_binding_size: u64::from(limits.max_storage_buffer_binding_size),
        max_compute_workgroups_per_dimension: limits.max_compute_workgroups_per_dimension,
        max_compute_invocations_per_workgroup: limits.max_compute_invocations_per_workgroup,
        max_compute_workgroup_size_x: limits.max_compute_workgroup_size_x,
        max_compute_workgroup_size_y: limits.max_compute_workgroup_size_y,
        max_compute_workgroup_size_z: limits.max_compute_workgroup_size_z,
    })
}

fn validate_buffer_range(size: u64, offset: u64, len: u64) -> Result<(), String> {
    let end = offset
        .checked_add(len)
        .ok_or_else(|| "GPU buffer range overflow".to_string())?;
    if end > size {
        return Err(format!(
            "GPU buffer range {}..{} exceeds buffer size {}",
            offset, end, size
        ));
    }
    Ok(())
}

fn validate_dispatch_workgroups(
    workgroups: [u32; 3],
    max_workgroups_per_dimension: u32,
) -> Result<(), String> {
    for (axis, count) in ["x", "y", "z"].into_iter().zip(workgroups) {
        if count > max_workgroups_per_dimension {
            return Err(format!(
                "compute dispatch workgroups.{axis}={} exceeds GPU limit {}",
                count, max_workgroups_per_dimension
            ));
        }
    }
    Ok(())
}

fn map_readback_buffer(
    device: &wgpu::Device,
    buffer: wgpu::Buffer,
    size: u64,
) -> Result<Vec<u8>, String> {
    let slice = buffer.slice(0..size);
    let (tx, rx) = std::sync::mpsc::channel();
    slice.map_async(wgpu::MapMode::Read, move |result| {
        let _ = tx.send(result);
    });
    device
        .poll(wgpu::PollType::Wait {
            submission_index: None,
            timeout: None,
        })
        .map_err(|error| error.to_string())?;
    rx.recv()
        .map_err(|error| error.to_string())?
        .map_err(|error| error.to_string())?;
    let mapped = slice.get_mapped_range();
    let bytes = mapped.to_vec();
    drop(mapped);
    buffer.unmap();
    Ok(bytes)
}

fn wgpu_buffer_usage(usage: GpuBufferUsage) -> Result<wgpu::BufferUsages, String> {
    let known = GpuBufferUsage::MAP_READ
        .union(GpuBufferUsage::COPY_SRC)
        .union(GpuBufferUsage::COPY_DST)
        .union(GpuBufferUsage::INDEX)
        .union(GpuBufferUsage::VERTEX)
        .union(GpuBufferUsage::UNIFORM)
        .union(GpuBufferUsage::STORAGE)
        .union(GpuBufferUsage::INDIRECT)
        .union(GpuBufferUsage::QUERY_RESOLVE);
    if usage.bits & !known.bits != 0 {
        return Err(format!("unknown GPU buffer usage bits: {}", usage.bits));
    }
    let mut out = wgpu::BufferUsages::empty();
    if usage.contains(GpuBufferUsage::MAP_READ) {
        out |= wgpu::BufferUsages::MAP_READ;
    }
    if usage.contains(GpuBufferUsage::COPY_SRC) {
        out |= wgpu::BufferUsages::COPY_SRC;
    }
    if usage.contains(GpuBufferUsage::COPY_DST) {
        out |= wgpu::BufferUsages::COPY_DST;
    }
    if usage.contains(GpuBufferUsage::INDEX) {
        out |= wgpu::BufferUsages::INDEX;
    }
    if usage.contains(GpuBufferUsage::VERTEX) {
        out |= wgpu::BufferUsages::VERTEX;
    }
    if usage.contains(GpuBufferUsage::UNIFORM) {
        out |= wgpu::BufferUsages::UNIFORM;
    }
    if usage.contains(GpuBufferUsage::STORAGE) {
        out |= wgpu::BufferUsages::STORAGE;
    }
    if usage.contains(GpuBufferUsage::INDIRECT) {
        out |= wgpu::BufferUsages::INDIRECT;
    }
    if usage.contains(GpuBufferUsage::QUERY_RESOLVE) {
        out |= wgpu::BufferUsages::QUERY_RESOLVE;
    }
    Ok(out)
}

fn wgpu_shader_stages(stages: GpuShaderStages) -> Result<wgpu::ShaderStages, String> {
    let known = GpuShaderStages::VERTEX
        .union(GpuShaderStages::FRAGMENT)
        .union(GpuShaderStages::COMPUTE);
    if stages.bits & !known.bits != 0 {
        return Err(format!("unknown GPU shader stage bits: {}", stages.bits));
    }
    let mut out = wgpu::ShaderStages::empty();
    if stages.bits & GpuShaderStages::VERTEX.bits != 0 {
        out |= wgpu::ShaderStages::VERTEX;
    }
    if stages.bits & GpuShaderStages::FRAGMENT.bits != 0 {
        out |= wgpu::ShaderStages::FRAGMENT;
    }
    if stages.bits & GpuShaderStages::COMPUTE.bits != 0 {
        out |= wgpu::ShaderStages::COMPUTE;
    }
    Ok(out)
}

fn wgpu_bind_group_layout_entry(
    entry: patinae_scene::GpuBindGroupLayoutEntry,
) -> Result<wgpu::BindGroupLayoutEntry, String> {
    Ok(wgpu::BindGroupLayoutEntry {
        binding: entry.binding,
        visibility: wgpu_shader_stages(entry.visibility)?,
        ty: wgpu_binding_type(entry.ty)?,
        count: None,
    })
}

fn wgpu_binding_type(ty: GpuBindingType) -> Result<wgpu::BindingType, String> {
    match ty {
        GpuBindingType::Buffer {
            ty,
            has_dynamic_offset,
            min_binding_size,
        } => Ok(wgpu::BindingType::Buffer {
            ty: wgpu_buffer_binding_type(ty),
            has_dynamic_offset,
            min_binding_size: min_binding_size.and_then(wgpu::BufferSize::new),
        }),
    }
}

fn wgpu_buffer_binding_type(ty: GpuBufferBindingType) -> wgpu::BufferBindingType {
    match ty {
        GpuBufferBindingType::Uniform => wgpu::BufferBindingType::Uniform,
        GpuBufferBindingType::StorageReadOnly => {
            wgpu::BufferBindingType::Storage { read_only: true }
        }
        GpuBufferBindingType::StorageReadWrite => {
            wgpu::BufferBindingType::Storage { read_only: false }
        }
    }
}

fn wgpu_bind_group_entry<'a>(
    registry: &'a GpuHandleRegistry,
    entry: GpuBindGroupEntry,
) -> Result<wgpu::BindGroupEntry<'a>, String> {
    match entry.resource {
        GpuBindingResource::Buffer(binding) => {
            let buffer = registry.buffer(binding.buffer)?;
            let size = binding.size.and_then(wgpu::BufferSize::new);
            let len = binding
                .size
                .unwrap_or(buffer.size.saturating_sub(binding.offset));
            validate_buffer_range(buffer.size, binding.offset, len)?;
            Ok(wgpu::BindGroupEntry {
                binding: entry.binding,
                resource: wgpu::BindingResource::Buffer(wgpu::BufferBinding {
                    buffer: &buffer.buffer,
                    offset: binding.offset,
                    size,
                }),
            })
        }
    }
}

fn map_render_artifact_buffer_role(
    role: patinae_render::RenderArtifactBufferRole,
) -> RenderArtifactBufferRole {
    match role {
        patinae_render::RenderArtifactBufferRole::FrameUniforms => {
            RenderArtifactBufferRole::FrameUniforms
        }
        patinae_render::RenderArtifactBufferRole::SceneAtoms => {
            RenderArtifactBufferRole::SceneAtoms
        }
        patinae_render::RenderArtifactBufferRole::SceneCoords => {
            RenderArtifactBufferRole::SceneCoords
        }
        patinae_render::RenderArtifactBufferRole::SceneBonds => {
            RenderArtifactBufferRole::SceneBonds
        }
        patinae_render::RenderArtifactBufferRole::SceneColorLut => {
            RenderArtifactBufferRole::SceneColorLut
        }
        patinae_render::RenderArtifactBufferRole::SceneMaskLut => {
            RenderArtifactBufferRole::SceneMaskLut
        }
        patinae_render::RenderArtifactBufferRole::SceneMarkerLut => {
            RenderArtifactBufferRole::SceneMarkerLut
        }
        patinae_render::RenderArtifactBufferRole::SceneCsrOffsets => {
            RenderArtifactBufferRole::SceneCsrOffsets
        }
        patinae_render::RenderArtifactBufferRole::SceneCsrIndices => {
            RenderArtifactBufferRole::SceneCsrIndices
        }
        patinae_render::RenderArtifactBufferRole::SceneObjectTable => {
            RenderArtifactBufferRole::SceneObjectTable
        }
        patinae_render::RenderArtifactBufferRole::SphereInstances => {
            RenderArtifactBufferRole::SphereInstances
        }
        patinae_render::RenderArtifactBufferRole::StickInstances => {
            RenderArtifactBufferRole::StickInstances
        }
        patinae_render::RenderArtifactBufferRole::LineInstances => {
            RenderArtifactBufferRole::LineInstances
        }
        patinae_render::RenderArtifactBufferRole::StdVertices => {
            RenderArtifactBufferRole::StdVertices
        }
        patinae_render::RenderArtifactBufferRole::InstanceCount => {
            RenderArtifactBufferRole::InstanceCount
        }
        patinae_render::RenderArtifactBufferRole::IndirectDraw => {
            RenderArtifactBufferRole::IndirectDraw
        }
    }
}

fn map_render_artifact_topology(
    topology: patinae_render::RenderArtifactPrimitiveTopology,
) -> RenderArtifactPrimitiveTopology {
    match topology {
        patinae_render::RenderArtifactPrimitiveTopology::SphereInstances => {
            RenderArtifactPrimitiveTopology::SphereInstances
        }
        patinae_render::RenderArtifactPrimitiveTopology::CylinderInstances => {
            RenderArtifactPrimitiveTopology::CylinderInstances
        }
        patinae_render::RenderArtifactPrimitiveTopology::LineInstances => {
            RenderArtifactPrimitiveTopology::LineInstances
        }
        patinae_render::RenderArtifactPrimitiveTopology::TriangleList => {
            RenderArtifactPrimitiveTopology::TriangleList
        }
        patinae_render::RenderArtifactPrimitiveTopology::LineList => {
            RenderArtifactPrimitiveTopology::LineList
        }
    }
}

fn map_render_artifact_rep_kind(kind: patinae_render::RepKind) -> RenderArtifactRepKind {
    match kind {
        patinae_render::RepKind::Sphere => RenderArtifactRepKind::Sphere,
        patinae_render::RepKind::Stick => RenderArtifactRepKind::Stick,
        patinae_render::RepKind::Line => RenderArtifactRepKind::Line,
        patinae_render::RepKind::Cartoon => RenderArtifactRepKind::Cartoon,
        patinae_render::RepKind::Ribbon => RenderArtifactRepKind::Ribbon,
        patinae_render::RepKind::Surface => RenderArtifactRepKind::Surface,
        patinae_render::RepKind::Mesh => RenderArtifactRepKind::Mesh,
        patinae_render::RepKind::Dot => RenderArtifactRepKind::Dot,
        patinae_render::RepKind::Ellipsoid => RenderArtifactRepKind::Ellipsoid,
        _ => RenderArtifactRepKind::Other,
    }
}

struct RuntimePayloadSpool {
    path: PathBuf,
}

impl Drop for RuntimePayloadSpool {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.path);
    }
}

fn command_runtime_input_from_context<V: ViewerLike + ?Sized>(
    ctx: &mut CommandContext<'_, '_, V>,
    args: &ParsedCommand,
    runtime_requirements: CommandRuntimeRequirements,
) -> Result<CommandRuntimeInput, String> {
    let mut spools = Vec::new();
    let input = build_command_input_from_context(
        ctx,
        args,
        runtime_requirements,
        DisplayedGeometryInputMode::Spool(&mut spools),
    )?;
    Ok(CommandRuntimeInput {
        input,
        _spools: spools,
    })
}

#[cfg(test)]
pub(crate) fn command_input_from_context<V: ViewerLike + ?Sized>(
    ctx: &mut CommandContext<'_, '_, V>,
    args: &ParsedCommand,
    runtime_requirements: CommandRuntimeRequirements,
) -> Result<WireCommandInput, String> {
    build_command_input_from_context(
        ctx,
        args,
        runtime_requirements,
        DisplayedGeometryInputMode::Inline,
    )
}

enum DisplayedGeometryInputMode<'a> {
    #[cfg(test)]
    Inline,
    Spool(&'a mut Vec<RuntimePayloadSpool>),
}

fn build_command_input_from_context<V: ViewerLike + ?Sized>(
    ctx: &mut CommandContext<'_, '_, V>,
    args: &ParsedCommand,
    runtime_requirements: CommandRuntimeRequirements,
    geometry_mode: DisplayedGeometryInputMode<'_>,
) -> Result<WireCommandInput, String> {
    let (viewport_width, viewport_height) = ctx.viewer.viewport_size();
    let (displayed_geometry, displayed_geometry_spool) =
        command_displayed_geometry_input(ctx, runtime_requirements, geometry_mode)?;
    Ok(WireCommandInput {
        wire_version: RUNTIME_WIRE_VERSION,
        session: command_session_bytes(ctx.viewer, runtime_requirements)?,
        viewport_image: ctx.viewer.viewport_image_ref().cloned(),
        parsed: args.clone(),
        quiet: ctx.quiet,
        viewport_width,
        viewport_height,
        dynamic_settings: dynamic_settings_to_wire(ctx.dynamic_settings()),
        displayed_geometry,
        displayed_geometry_spool,
    })
}

fn command_displayed_geometry_input<V: ViewerLike + ?Sized>(
    ctx: &mut CommandContext<'_, '_, V>,
    runtime_requirements: CommandRuntimeRequirements,
    geometry_mode: DisplayedGeometryInputMode<'_>,
) -> Result<
    (
        Option<wire::WireDisplayedGeometry>,
        Option<wire::WireDisplayedGeometrySpool>,
    ),
    String,
> {
    if !runtime_requirements.contains(CommandRuntimeRequirements::DISPLAYED_GEOMETRY) {
        return Ok((None, None));
    }

    let displayed = ctx
        .viewer
        .export_displayed_geometry(&patinae_render::GeometryExportOptions::default())?;
    match geometry_mode {
        #[cfg(test)]
        DisplayedGeometryInputMode::Inline => {
            Ok((Some(wire::displayed_geometry_to_wire(&displayed)), None))
        }
        DisplayedGeometryInputMode::Spool(spools) => {
            let spool = write_displayed_geometry_spool(&displayed)?;
            let descriptor = wire::WireDisplayedGeometrySpool {
                path: spool.path.to_string_lossy().into_owned(),
                chunk_count: spool.chunk_count,
            };
            spools.push(RuntimePayloadSpool { path: spool.path });
            Ok((None, Some(descriptor)))
        }
    }
}

struct DisplayedGeometrySpool {
    path: PathBuf,
    chunk_count: usize,
}

fn write_displayed_geometry_spool(
    displayed: &patinae_render::DisplayedGeometry,
) -> Result<DisplayedGeometrySpool, String> {
    let (path, mut file) = create_displayed_geometry_spool_file()?;
    let cleanup = RuntimePayloadSpool { path: path.clone() };
    file.write_all(wire::DISPLAYED_GEOMETRY_SPOOL_MAGIC)
        .map_err(|error| format!("write displayed geometry spool header: {error}"))?;

    let mut chunk_count = 0usize;
    for object in &displayed.objects {
        write_displayed_object_spool_chunks(object, &mut file, &mut chunk_count)?;
    }
    file.flush()
        .map_err(|error| format!("flush displayed geometry spool: {error}"))?;

    std::mem::forget(cleanup);
    Ok(DisplayedGeometrySpool { path, chunk_count })
}

fn create_displayed_geometry_spool_file() -> Result<(PathBuf, std::fs::File), String> {
    for _ in 0..32 {
        let path = next_displayed_geometry_spool_path();
        match OpenOptions::new().write(true).create_new(true).open(&path) {
            Ok(file) => return Ok((path, file)),
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(error) => {
                return Err(format!("create displayed geometry spool: {error}"));
            }
        }
    }
    Err("create displayed geometry spool: exhausted unique path attempts".to_string())
}

fn next_displayed_geometry_spool_path() -> PathBuf {
    let id = DISPLAYED_GEOMETRY_SPOOL_COUNTER.fetch_add(1, Ordering::Relaxed);
    std::env::temp_dir().join(format!(
        "patinae-displayed-geometry-{}-{id}.mpack",
        std::process::id()
    ))
}

fn write_displayed_object_spool_chunks(
    object: &patinae_render::DisplayedObjectGeometry,
    file: &mut std::fs::File,
    chunk_count: &mut usize,
) -> Result<(), String> {
    let object_id = object.object_id.0;
    let mut primitives = Vec::new();
    for primitive in &object.primitives {
        match primitive {
            patinae_render::DisplayedPrimitive::Mesh { rep, mesh } => {
                write_displayed_spool_primitive_chunk(
                    object_id,
                    &mut primitives,
                    file,
                    chunk_count,
                )?;
                for vertices in mesh.vertices.chunks(DISPLAYED_GEOMETRY_SPOOL_MESH_VERTICES) {
                    let chunk_primitive = patinae_render::DisplayedPrimitive::Mesh {
                        rep: *rep,
                        mesh: patinae_render::DisplayedMesh {
                            vertices: vertices.to_vec(),
                        },
                    };
                    primitives.push(wire::displayed_primitive_to_wire(&chunk_primitive));
                    write_displayed_spool_primitive_chunk(
                        object_id,
                        &mut primitives,
                        file,
                        chunk_count,
                    )?;
                }
            }
            _ => {
                primitives.push(wire::displayed_primitive_to_wire(primitive));
                if primitives.len() >= DISPLAYED_GEOMETRY_SPOOL_PRIMITIVES {
                    write_displayed_spool_primitive_chunk(
                        object_id,
                        &mut primitives,
                        file,
                        chunk_count,
                    )?;
                }
            }
        }
    }

    write_displayed_spool_primitive_chunk(object_id, &mut primitives, file, chunk_count)
}

fn write_displayed_spool_primitive_chunk(
    object_id: u32,
    primitives: &mut Vec<wire::WireDisplayedPrimitive>,
    file: &mut std::fs::File,
    chunk_count: &mut usize,
) -> Result<(), String> {
    if primitives.is_empty() {
        return Ok(());
    }

    let chunk = wire::WireDisplayedGeometryChunk {
        objects: vec![wire::WireDisplayedObjectGeometry {
            object_id,
            primitives: std::mem::take(primitives),
        }],
    };
    write_displayed_spool_chunk(file, &chunk)?;
    *chunk_count += 1;
    Ok(())
}

fn write_displayed_spool_chunk(
    file: &mut std::fs::File,
    chunk: &wire::WireDisplayedGeometryChunk,
) -> Result<(), String> {
    let bytes = wire::encode(chunk)?;
    ensure_runtime_payload_within_limit(
        "command execute",
        bytes.len(),
        "displayed geometry chunk",
    )?;
    file.write_all(&(bytes.len() as u64).to_le_bytes())
        .map_err(|error| format!("write displayed geometry chunk length: {error}"))?;
    file.write_all(&bytes)
        .map_err(|error| format!("write displayed geometry chunk: {error}"))
}

pub(crate) fn apply_command_output<V: ViewerLike + ?Sized>(
    ctx: &mut CommandContext<'_, '_, V>,
    output: WireCommandOutput,
    runtime_requirements: CommandRuntimeRequirements,
) -> CmdResult {
    validate_runtime_wire_version(output.wire_version).map_err(CmdError::execution)?;
    if runtime_requirements.contains(CommandRuntimeRequirements::FULL_SESSION) {
        let session = wire::decode_session(&output.session).map_err(CmdError::execution)?;
        ctx.viewer.replace_session(session);
    }
    ctx.viewer.set_viewport_image(output.viewport_image);
    for message in output.output {
        apply_output_message(ctx, message);
    }
    for action in output.actions {
        apply_command_action(ctx, action);
    }
    output.result.map_err(CmdError::execution)
}

fn apply_output_message<V: ViewerLike + ?Sized>(
    ctx: &mut CommandContext<'_, '_, V>,
    message: patinae_cmd::OutputMessage,
) {
    match message.kind {
        MessageKind::Info => ctx.print(&message.text),
        MessageKind::Warning => ctx.print_warning(&message.text),
        MessageKind::Error => ctx.print_error(&message.text),
    }
}

fn apply_command_action<V: ViewerLike + ?Sized>(
    ctx: &mut CommandContext<'_, '_, V>,
    action: CommandAction,
) {
    match action {
        CommandAction::ShowPanel(panel) => ctx.show_panel(panel),
        CommandAction::HidePanel(panel) => ctx.hide_panel(panel),
        CommandAction::ClearOutput => ctx.clear_output(),
        CommandAction::Quit => ctx.quit(),
        CommandAction::RecordRecentFile { path, command } => ctx.record_recent_file(path, command),
    }
}

pub(crate) fn shared_input_from_context(
    ctx: &SharedContext<'_>,
    runtime_requirements: PanelRuntimeRequirements,
) -> Result<WireSharedInput, String> {
    let session = session_from_shared_context(ctx, runtime_requirements)?;
    Ok(WireSharedInput {
        wire_version: RUNTIME_WIRE_VERSION,
        session: wire::encode_session(&session)?,
        scene_generation: ctx.scene_generation,
        viewport_image: ctx.viewport_image.cloned(),
        command_names: ctx.command_names.to_vec(),
        setting_names: ctx
            .setting_names
            .iter()
            .map(|name| (*name).to_string())
            .collect(),
        dynamic_settings: dynamic_settings_to_wire(ctx.dynamic_settings),
    })
}

fn session_from_shared_context(
    ctx: &SharedContext<'_>,
    runtime_requirements: PanelRuntimeRequirements,
) -> Result<Session, String> {
    if !runtime_requirements.contains(PanelRuntimeRequirements::FULL_SESSION) {
        return lightweight_session_from_shared_context(ctx);
    }

    let mut session = Session::new();
    session.registry = patinae_scene::ObjectRegistry::from_snapshot(ctx.registry.to_snapshot());
    session.camera = ctx.camera.clone();
    session.selections = ctx.selections.clone();
    session.named_palette = clone_wire(ctx.named_palette)?;
    session.movie = clone_wire(ctx.movie)?;
    session.settings = ctx.settings.clone();
    session.clear_color = ctx.clear_color;
    session.viewport_image = ctx.viewport_image.cloned();
    ensure_full_session_payload_within_limit("panel snapshot", &session)?;
    Ok(session)
}

fn command_session_bytes<V: ViewerLike + ?Sized>(
    viewer: &V,
    runtime_requirements: CommandRuntimeRequirements,
) -> Result<Vec<u8>, String> {
    if runtime_requirements.contains(CommandRuntimeRequirements::FULL_SESSION) {
        let bytes = wire::encode_session(viewer.session())?;
        ensure_runtime_payload_within_limit("command execute", bytes.len(), "full scene state")?;
        return Ok(bytes);
    }

    let session = lightweight_session_from_viewer(viewer)?;
    wire::encode_session(&session)
}

fn lightweight_session_from_viewer<V: ViewerLike + ?Sized>(viewer: &V) -> Result<Session, String> {
    let source = viewer.session();
    let mut session = Session::new();
    session.camera = viewer.camera().clone();
    session.settings = viewer.settings().clone();
    session.named_palette = clone_wire(&source.named_palette)?;
    session.movie = clone_wire(&source.movie)?;
    session.clear_color = source.clear_color;
    session.clear_color_set = source.clear_color_set;
    session.viewport_image = viewer.viewport_image_ref().cloned();
    Ok(session)
}

fn lightweight_session_from_shared_context(ctx: &SharedContext<'_>) -> Result<Session, String> {
    let mut session = Session::new();
    session.camera = ctx.camera.clone();
    session.settings = ctx.settings.clone();
    session.named_palette = clone_wire(ctx.named_palette)?;
    session.movie = clone_wire(ctx.movie)?;
    session.clear_color = ctx.clear_color;
    session.viewport_image = ctx.viewport_image.cloned();
    Ok(session)
}

fn ensure_full_session_payload_within_limit(phase: &str, session: &Session) -> Result<(), String> {
    let bytes = wire::encode_session(session)?;
    ensure_runtime_payload_within_limit(phase, bytes.len(), "full scene state")
}

pub(crate) fn ensure_runtime_payload_within_limit(
    phase: &str,
    len: usize,
    requirement: &str,
) -> Result<(), String> {
    if len <= wire::MAX_WIRE_PAYLOAD_LEN {
        return Ok(());
    }
    let limit_mib = wire::MAX_WIRE_PAYLOAD_LEN / (1024 * 1024);
    Err(format!(
        "Plugin {phase} input requires {requirement}, but serialized payload is {len} bytes and exceeds the {limit_mib} MiB plugin ABI limit"
    ))
}

fn clone_wire<T>(value: &T) -> Result<T, String>
where
    T: Serialize + DeserializeOwned,
{
    let bytes = wire::encode(value)?;
    wire::decode(&bytes)
}

fn dynamic_settings_to_wire(
    dynamic_settings: Option<&DynamicSettingRegistry>,
) -> Vec<wire::WireDynamicSetting> {
    let Some(dynamic_settings) = dynamic_settings else {
        return Vec::new();
    };
    dynamic_settings
        .names()
        .iter()
        .filter_map(|name| {
            let entry = dynamic_settings.lookup(name)?;
            let value = entry
                .store
                .read()
                .ok()
                .and_then(|store| store.get(name).cloned());
            Some(wire::WireDynamicSetting {
                descriptor: (&entry.descriptor).into(),
                value,
            })
        })
        .collect()
}

fn poll_input_from_context(ctx: &PollContext<'_>) -> Result<WirePollInput, String> {
    Ok(WirePollInput {
        shared: ctx.poll_shared.clone(),
        command_results: ctx
            .command_results
            .iter()
            .map(|result| wire::WireCommandResult {
                id: result.id,
                result: result.result.clone(),
            })
            .collect(),
        host_query_results: ctx.host_query_results.to_vec(),
        dynamic_invocations: ctx.dynamic_invocations.to_vec(),
        plugin_dirs: ctx
            .plugin_dirs
            .iter()
            .map(|path| path.to_string_lossy().into_owned())
            .collect(),
    })
}

fn apply_poll_output(ctx: &mut PollContext<'_>, output: WirePollOutput) -> Result<(), String> {
    validate_runtime_wire_version(output.wire_version)?;
    send_bus_messages(output.messages, ctx.bus);
    for request in output.command_exec {
        ctx.execute_command(request.id, &request.command, request.silent);
    }
    for registration in output.dynamic_registrations {
        ctx.register_dynamic_command(
            registration.name,
            registration.description,
            registration.usage,
            registration.arguments,
        );
    }
    for name in output.dynamic_unregistrations {
        ctx.unregister_dynamic_command(&name);
    }
    for notification in output.notifications {
        ctx.set_notification(notification);
    }
    for registration in output.hotkey_registrations {
        ctx.register_hotkey(
            registration.key,
            plugin_key_action_from_wire(registration.action),
        );
    }
    for key in output.hotkey_unregistrations {
        ctx.unregister_hotkey(key);
    }
    for query in output.host_queries {
        ctx.query_host(query);
    }
    for action in output.viewer_actions {
        apply_viewer_action(ctx, action);
    }
    Ok(())
}

fn apply_viewer_action(ctx: &mut PollContext<'_>, action: WireViewerAction) {
    match action {
        WireViewerAction::ApplyAtomPropertyChanges(changes) => {
            ctx.queue_viewer_mutation(move |viewer| {
                let mut changed = false;
                for change in &changes {
                    let Some(mol_obj) = viewer.objects_mut().get_molecule_mut(&change.object)
                    else {
                        continue;
                    };
                    let Some(atom) = mol_obj
                        .molecule_mut()
                        .get_atom_mut(AtomIndex(change.atom_index))
                    else {
                        continue;
                    };
                    changed |= apply_atom_property_changes(atom, &change.changes);
                }
                if changed {
                    viewer.request_redraw();
                }
            });
        }
        WireViewerAction::SetViewportImage(image) => {
            ctx.queue_viewer_mutation(move |viewer| {
                viewer.set_viewport_image(Some(image));
            });
        }
        WireViewerAction::ClearViewportImage => {
            ctx.queue_viewer_mutation(|viewer| {
                viewer.set_viewport_image(None);
            });
        }
        WireViewerAction::RequestRedraw => {
            ctx.queue_viewer_mutation(|viewer| viewer.request_redraw());
        }
        WireViewerAction::RequestPanelUpdate => ctx.request_panel_update(),
    }
}

fn apply_atom_property_changes(
    atom: &mut Atom,
    changes: &[(String, wire::WireAtomPropertyValue)],
) -> bool {
    let mut changed = false;
    for (key, value) in changes {
        changed |= apply_atom_property_change(atom, key, value);
    }
    changed
}

fn apply_atom_property_change(
    atom: &mut Atom,
    key: &str,
    value: &wire::WireAtomPropertyValue,
) -> bool {
    match (key, value) {
        ("name", wire::WireAtomPropertyValue::Str(value)) => {
            atom.name = Arc::from(value.as_str());
        }
        ("b", wire::WireAtomPropertyValue::F32(value)) => atom.b_factor = *value,
        ("q", wire::WireAtomPropertyValue::F32(value)) => atom.occupancy = *value,
        ("vdw", wire::WireAtomPropertyValue::F32(value)) => atom.vdw = *value,
        ("partial_charge", wire::WireAtomPropertyValue::F32(value)) => {
            atom.partial_charge = *value;
        }
        ("formal_charge", wire::WireAtomPropertyValue::I8(value)) => {
            atom.formal_charge = *value;
        }
        ("color", wire::WireAtomPropertyValue::I32(value)) => atom.repr.colors.base = *value,
        ("elem", wire::WireAtomPropertyValue::Str(value)) => {
            let Some(element) = Element::from_symbol(value) else {
                return false;
            };
            atom.element = element;
        }
        ("ss", wire::WireAtomPropertyValue::Str(value)) => {
            atom.ss_type = match value.as_str() {
                "H" => SecondaryStructure::Helix,
                "S" => SecondaryStructure::Sheet,
                _ => SecondaryStructure::Loop,
            };
        }
        ("type", wire::WireAtomPropertyValue::Str(value)) => atom.state.hetatm = value == "HETATM",
        ("alt", wire::WireAtomPropertyValue::Str(value)) => {
            atom.alt = value.chars().next().unwrap_or(' ');
        }
        ("chain", wire::WireAtomPropertyValue::Str(value)) => {
            let mut residue = (*atom.residue).clone();
            residue.key.chain = value.clone();
            atom.residue = Arc::new(residue);
        }
        ("resn", wire::WireAtomPropertyValue::Str(value)) => {
            let mut residue = (*atom.residue).clone();
            residue.key.resn = value.clone();
            atom.residue = Arc::new(residue);
        }
        ("resv", wire::WireAtomPropertyValue::I32(value)) => {
            let mut residue = (*atom.residue).clone();
            residue.key.resv = *value;
            atom.residue = Arc::new(residue);
        }
        ("segi", wire::WireAtomPropertyValue::Str(value)) => {
            let mut residue = (*atom.residue).clone();
            residue.segi = value.clone();
            atom.residue = Arc::new(residue);
        }
        _ => return false,
    }
    true
}

fn send_bus_messages(messages: Vec<AppMessage>, bus: &mut MessageBus) {
    for message in messages {
        bus.send(message);
    }
}

fn plugin_key_action_from_wire(action: WireHotkeyAction) -> PluginKeyAction {
    match action {
        WireHotkeyAction::Command(command) => PluginKeyAction::Command(command),
        WireHotkeyAction::DynamicCommand { name, args } => {
            PluginKeyAction::DynamicCommand { name, args }
        }
        WireHotkeyAction::Custom { topic, payload } => PluginKeyAction::Custom { topic, payload },
    }
}

fn destroy_command_handle(handle: PluginCommandHandle, vtable: AbiCommandVTable) {
    if let Some(destroy) = vtable.destroy {
        destroy_runtime_handle("command destroy", || {
            // SAFETY: The handle and destroy callback were validated during
            // registration, and the proxy keeps the library alive during drop.
            unsafe { destroy(handle) }
        });
    }
}

fn destroy_panel_handle(handle: PluginPanelHandle, vtable: AbiPanelVTable) {
    if let Some(destroy) = vtable.destroy {
        destroy_runtime_handle("panel destroy", || {
            // SAFETY: The handle and destroy callback were validated during
            // registration, and the proxy keeps the library alive during drop.
            unsafe { destroy(handle) }
        });
    }
}

fn destroy_message_handler_handle(
    handle: PluginMessageHandlerHandle,
    vtable: AbiMessageHandlerVTable,
) {
    if let Some(destroy) = vtable.destroy {
        destroy_runtime_handle("message handler destroy", || {
            // SAFETY: The handle and destroy callback were validated during
            // registration, and the proxy keeps the library alive during drop.
            unsafe { destroy(handle) }
        });
    }
}

fn destroy_script_handler_handle(
    handle: PluginScriptHandlerHandle,
    vtable: AbiScriptHandlerVTable,
) {
    if let Some(destroy) = vtable.destroy {
        destroy_runtime_handle("script handler destroy", || {
            // SAFETY: The handle and destroy callback were validated during
            // registration, and the proxy keeps the library alive during drop.
            unsafe { destroy(handle) }
        });
    }
}

fn destroy_format_handler_handle(
    handle: PluginFormatHandlerHandle,
    vtable: AbiFormatHandlerVTable,
) {
    if let Some(destroy) = vtable.destroy {
        destroy_runtime_handle("format handler destroy", || {
            // SAFETY: The handle and destroy callback were validated during
            // registration, and the proxy keeps the library alive during drop.
            unsafe { destroy(handle) }
        });
    }
}

fn destroy_hotkey_handle(handle: PluginHotkeyHandle, vtable: AbiHotkeyVTable) {
    if let Some(destroy) = vtable.destroy {
        destroy_runtime_handle("hotkey destroy", || {
            // SAFETY: The handle and destroy callback were validated during
            // registration, and the proxy keeps the library alive during drop.
            unsafe { destroy(handle) }
        });
    }
}

fn destroy_runtime_handle(phase: &str, destroy: impl FnOnce() -> AbiStatus) {
    match catch_unwind(AssertUnwindSafe(destroy)) {
        Ok(status) if status.is_ok() => {}
        Ok(status) => {
            if let Err(error) = status_to_result(phase, status) {
                log::warn!("{error}");
            }
        }
        Err(panic_info) => log::error!(
            "Plugin {phase} panicked: {}",
            panic_payload_to_string(&panic_info)
        ),
    }
}

fn leak_aliases(values: Vec<String>) -> &'static [&'static str] {
    if values.is_empty() {
        return &[];
    }
    let aliases: Vec<&'static str> = values
        .into_iter()
        .map(|value| {
            let value: &'static mut str = Box::leak(value.into_boxed_str());
            &*value
        })
        .collect();
    Box::leak(aliases.into_boxed_slice())
}

fn leak_arg_hints(values: Vec<ArgHint>) -> &'static [ArgHint] {
    if values.is_empty() {
        return &[];
    }
    Box::leak(values.into_boxed_slice())
}

fn arg_hint_from_abi(value: u8) -> ArgHint {
    match value {
        ARG_HINT_PATH => ArgHint::Path,
        ARG_HINT_SELECTION => ArgHint::Selection,
        ARG_HINT_OBJECT => ArgHint::Object,
        ARG_HINT_REPRESENTATION => ArgHint::Representation,
        ARG_HINT_COLOR => ArgHint::Color,
        ARG_HINT_SETTING => ArgHint::Setting,
        ARG_HINT_SETTING_VALUE => ArgHint::SettingValue,
        ARG_HINT_NAMED_SELECTION => ArgHint::NamedSelection,
        ARG_HINT_LABEL_PROPERTY => ArgHint::LabelProperty,
        ARG_HINT_COMMAND => ArgHint::Command,
        _ => ArgHint::None,
    }
}

fn panel_descriptor_from_abi(descriptor: &AbiPanelDescriptor) -> Result<PanelDescriptor, String> {
    let placement = match descriptor.placement {
        PANEL_PLACEMENT_RIGHT => PanelPlacement::Right,
        PANEL_PLACEMENT_BOTTOM => PanelPlacement::Bottom,
        other => return Err(format!("unknown panel placement code {other}")),
    };
    Ok(PanelDescriptor {
        id: copy_str(descriptor.id, "panel.id")?,
        title: copy_str(descriptor.title, "panel.title")?,
        icon: copy_str(descriptor.icon, "panel.icon")?,
        placement,
        default_visible: descriptor.default_visible != 0,
    })
}

fn setting_from_abi(
    descriptor: &AbiSettingDescriptor,
) -> Result<(DynamicSettingDescriptor, SharedSettingStore), String> {
    let name = copy_str(descriptor.name, "setting.name")?;
    let default = setting_value_from_abi(descriptor.default_value, "setting.default")?;
    let side_effects = copy_u8_slice(descriptor.side_effects, "setting.side_effects")?
        .into_iter()
        .map(side_effect_from_abi)
        .collect::<Result<Vec<_>, _>>()?;
    let value_hints = copy_hint_slice(descriptor.value_hints, "setting.value_hints")?;
    let setting = DynamicSettingDescriptor {
        name: name.clone(),
        setting_type: setting_type_from_abi(descriptor.setting_type)?,
        default: default.clone(),
        min: (descriptor.has_min != 0).then_some(descriptor.min),
        max: (descriptor.has_max != 0).then_some(descriptor.max),
        value_hints,
        side_effects,
        object_overridable: descriptor.object_overridable != 0,
    };
    let mut store = DynamicSettingStore::new();
    store.set(&name, default);
    Ok((setting, Arc::new(RwLock::new(store))))
}

fn setting_hint_from_abi(
    hint: &AbiSettingValueHint,
    field: &str,
) -> Result<(String, SettingValue), String> {
    Ok((
        copy_str(hint.name, &format!("{field}.name"))?,
        setting_value_from_abi(hint.value, &format!("{field}.value"))?,
    ))
}

fn setting_type_from_abi(value: u8) -> Result<SettingType, String> {
    match value {
        SETTING_TYPE_BLANK => Ok(SettingType::Blank),
        SETTING_TYPE_BOOL => Ok(SettingType::Bool),
        SETTING_TYPE_INT => Ok(SettingType::Int),
        SETTING_TYPE_FLOAT => Ok(SettingType::Float),
        SETTING_TYPE_FLOAT3 => Ok(SettingType::Float3),
        SETTING_TYPE_COLOR => Ok(SettingType::Color),
        SETTING_TYPE_STRING => Ok(SettingType::String),
        other => Err(format!("unknown setting type code {other}")),
    }
}

fn setting_value_from_abi(value: AbiSettingValue, field: &str) -> Result<SettingValue, String> {
    match value.tag {
        SETTING_VALUE_NONE => Ok(SettingValue::String(String::new())),
        SETTING_VALUE_BOOL => Ok(SettingValue::Bool(value.bool_value != 0)),
        SETTING_VALUE_INT => Ok(SettingValue::Int(value.int_value)),
        SETTING_VALUE_FLOAT => Ok(SettingValue::Float(value.float_values[0])),
        SETTING_VALUE_FLOAT3 => Ok(SettingValue::Float3(value.float_values)),
        SETTING_VALUE_COLOR => Ok(SettingValue::Color(value.int_value)),
        SETTING_VALUE_STRING => Ok(SettingValue::String(copy_str(
            value.string_value,
            &format!("{field}.string"),
        )?)),
        other => Err(format!("{field} has unknown setting value code {other}")),
    }
}

fn side_effect_from_abi(value: u8) -> Result<SideEffectCategory, String> {
    match value {
        SIDE_EFFECT_SCENE_INVALIDATE => Ok(SideEffectCategory::SceneInvalidate),
        SIDE_EFFECT_SCENE_CHANGED => Ok(SideEffectCategory::SceneChanged),
        SIDE_EFFECT_SHADER_RELOAD => Ok(SideEffectCategory::ShaderReload),
        SIDE_EFFECT_SHADER_COMPUTE_LIGHTING => Ok(SideEffectCategory::ShaderComputeLighting),
        SIDE_EFFECT_ORTHO_DIRTY => Ok(SideEffectCategory::OrthoDirty),
        SIDE_EFFECT_SEQ_CHANGED => Ok(SideEffectCategory::SeqChanged),
        SIDE_EFFECT_STEREO_UPDATE => Ok(SideEffectCategory::StereoUpdate),
        SIDE_EFFECT_REPRESENTATION_REBUILD => Ok(SideEffectCategory::RepresentationRebuild),
        SIDE_EFFECT_COLOR_REBUILD => Ok(SideEffectCategory::ColorRebuild),
        SIDE_EFFECT_SURFACE_TRANSPARENCY => Ok(SideEffectCategory::SurfaceTransparency),
        SIDE_EFFECT_FULL_REBUILD => Ok(SideEffectCategory::FullRebuild),
        SIDE_EFFECT_VIEWPORT_UPDATE => Ok(SideEffectCategory::ViewportUpdate),
        other => Err(format!("unknown side-effect code {other}")),
    }
}

#[cfg(test)]
mod loader_tests {
    use super::*;

    struct TraceFixtureViewer {
        session: Session,
        chunks: Vec<patinae_render::TraceGeometryChunk>,
    }

    impl ViewerLike for TraceFixtureViewer {
        fn objects(&self) -> &patinae_scene::ObjectRegistry {
            &self.session.registry
        }

        fn objects_mut(&mut self) -> &mut patinae_scene::ObjectRegistry {
            &mut self.session.registry
        }

        fn camera(&self) -> &patinae_scene::Camera {
            &self.session.camera
        }

        fn camera_mut(&mut self) -> &mut patinae_scene::Camera {
            &mut self.session.camera
        }

        fn settings(&self) -> &patinae_settings::Settings {
            &self.session.settings
        }

        fn settings_mut(&mut self) -> &mut patinae_settings::Settings {
            &mut self.session.settings
        }

        fn request_redraw(&mut self) {}

        fn session(&self) -> &Session {
            &self.session
        }

        fn session_mut(&mut self) -> &mut Session {
            &mut self.session
        }

        fn replace_session(&mut self, session: Session) {
            self.session = session;
        }

        fn movie(&self) -> &patinae_scene::Movie {
            &self.session.movie
        }

        fn movie_mut(&mut self) -> &mut patinae_scene::Movie {
            &mut self.session.movie
        }

        fn scenes(&self) -> &patinae_scene::SceneManager {
            &self.session.scenes
        }

        fn scenes_mut(&mut self) -> &mut patinae_scene::SceneManager {
            &mut self.session.scenes
        }

        fn scene_store(&mut self, key: &str, storemask: u32) {
            let mask = patinae_scene::SceneStoreMask::from_bits_truncate(storemask);
            self.session
                .scenes
                .store(key, mask, &self.session.camera, &self.session.registry);
        }

        fn scene_recall(&mut self, key: &str, animate: bool, duration: f32) -> Result<(), String> {
            self.session
                .scenes
                .recall(
                    key,
                    &mut self.session.camera,
                    &mut self.session.registry,
                    animate,
                    duration,
                )
                .map_err(|error| error.to_string())
        }

        fn views(&self) -> &patinae_scene::ViewManager {
            &self.session.views
        }

        fn views_mut(&mut self) -> &mut patinae_scene::ViewManager {
            &mut self.session.views
        }

        fn view_recall(&mut self, key: &str, animate: f32) -> Result<(), String> {
            self.session
                .views
                .recall(key, &mut self.session.camera, animate)
                .map_err(|error| error.to_string())
        }

        fn selections(&self) -> &patinae_scene::SelectionManager {
            &self.session.selections
        }

        fn selections_mut(&mut self) -> &mut patinae_scene::SelectionManager {
            &mut self.session.selections
        }

        fn named_palette(&self) -> &patinae_scene::NamedPalette {
            &self.session.named_palette
        }

        fn named_palette_mut(&mut self) -> &mut patinae_scene::NamedPalette {
            &mut self.session.named_palette
        }

        fn clear_color(&self) -> [f32; 3] {
            self.session.clear_color
        }

        fn set_clear_color(&mut self, color: [f32; 3]) {
            self.session.clear_color = color;
        }

        fn viewport_image_ref(&self) -> Option<&patinae_scene::ViewportImage> {
            self.session.viewport_image.as_ref()
        }

        fn set_viewport_image_internal(&mut self, image: Option<patinae_scene::ViewportImage>) {
            self.session.viewport_image = image;
        }

        fn for_each_trace_geometry_chunk(
            &mut self,
            _options: &patinae_render::GeometryExportOptions,
            visitor: &mut dyn FnMut(patinae_render::TraceGeometryChunk) -> Result<(), String>,
        ) -> Result<(), String> {
            for chunk in self.chunks.clone() {
                visitor(chunk)?;
            }
            Ok(())
        }
    }

    fn trace_material() -> patinae_render::TraceMaterial {
        patinae_render::TraceMaterial {
            rgba: [0.2, 0.4, 0.8, 1.0],
            transparency: 0.0,
        }
    }

    fn trace_chunk(
        sphere_count: usize,
        cylinder_count: usize,
        triangle_count: usize,
    ) -> patinae_render::TraceGeometryChunk {
        let material = trace_material();
        patinae_render::TraceGeometryChunk {
            spheres: (0..sphere_count)
                .map(|i| patinae_render::TraceSphere {
                    center: [i as f32, 0.0, 0.0],
                    radius: 1.0,
                    material,
                })
                .collect(),
            cylinders: (0..cylinder_count)
                .map(|i| patinae_render::TraceCylinder {
                    start: [i as f32, 0.0, 0.0],
                    end: [i as f32, 1.0, 0.0],
                    radius: 0.2,
                    material_start: material,
                    material_end: material,
                })
                .collect(),
            triangles: (0..triangle_count)
                .map(|i| patinae_render::TraceTriangle {
                    positions: [
                        [i as f32, 0.0, 0.0],
                        [i as f32, 1.0, 0.0],
                        [i as f32, 0.0, 1.0],
                    ],
                    normals: [[0.0, 0.0, 1.0]; 3],
                    material,
                })
                .collect(),
            line_segments: Vec::new(),
            point_samples: Vec::new(),
        }
    }

    #[test]
    fn displayed_geometry_spool_splits_and_reads_mesh_chunks() {
        let material = patinae_render::DisplayedMaterial::from_rgba([0.2, 0.4, 0.8, 1.0]);
        let vertex_count = DISPLAYED_GEOMETRY_SPOOL_MESH_VERTICES + 3;
        let vertices = (0..vertex_count)
            .map(|i| patinae_render::DisplayedMeshVertex {
                position: [i as f32, 0.0, 0.0],
                normal: [0.0, 0.0, 1.0],
                owner_atom_id: i as u32,
                material,
                flags: 0,
            })
            .collect::<Vec<_>>();
        let displayed = patinae_render::DisplayedGeometry {
            objects: vec![patinae_render::DisplayedObjectGeometry {
                object_id: patinae_render::ObjectId(7),
                primitives: vec![patinae_render::DisplayedPrimitive::Mesh {
                    rep: patinae_render::RepKind::Cartoon,
                    mesh: patinae_render::DisplayedMesh {
                        vertices: vertices.clone(),
                    },
                }],
            }],
        };

        let spool = write_displayed_geometry_spool(&displayed).unwrap();
        let _cleanup = RuntimePayloadSpool {
            path: spool.path.clone(),
        };
        assert!(spool.chunk_count >= 2);

        let descriptor = wire::WireDisplayedGeometrySpool {
            path: spool.path.to_string_lossy().into_owned(),
            chunk_count: spool.chunk_count,
        };
        let mut chunk_vertex_counts = Vec::new();
        wire::for_each_displayed_geometry_spool_chunk(&descriptor, |chunk| {
            for object in chunk.objects {
                for primitive in object.primitives {
                    if let patinae_render::DisplayedPrimitive::Mesh { mesh, .. } = primitive {
                        chunk_vertex_counts.push(mesh.vertices.len());
                    }
                }
            }
            Ok(())
        })
        .unwrap();

        assert_eq!(chunk_vertex_counts.iter().sum::<usize>(), vertices.len());
        assert!(chunk_vertex_counts
            .iter()
            .all(|&count| count <= DISPLAYED_GEOMETRY_SPOOL_MESH_VERTICES));
    }

    #[test]
    fn trace_geometry_chunk_stays_below_abi_limit() {
        let chunk = trace_chunk(0, 0, 65_536);
        let wire_chunk = wire::trace_geometry_chunk_to_wire(chunk);
        let encoded = wire::encode(&wire_chunk).unwrap();

        ensure_runtime_payload_within_limit(
            "test trace geometry stream",
            encoded.len(),
            "trace geometry chunk",
        )
        .unwrap();
    }

    #[test]
    fn trace_geometry_stream_reads_multiple_chunks() {
        let mut viewer = TraceFixtureViewer {
            session: Session::new(),
            chunks: vec![trace_chunk(1, 2, 3), trace_chunk(4, 5, 6)],
        };
        let mut state = HostCommandRuntimeState::new(
            &mut viewer,
            CommandRuntimeRequirements::TRACE_GEOMETRY_STREAM,
        );

        let opened =
            state.handle_request(WireCommandRuntimeRequest::OpenTraceGeometryStream { id: 1 });
        let stream_id = match opened.result.unwrap() {
            WireCommandRuntimeValue::TraceGeometryOpened(opened) => opened.stream_id,
            _ => panic!("unexpected open response"),
        };

        let mut counts = (0usize, 0usize, 0usize);
        for id in [2, 3] {
            let response =
                state.handle_request(WireCommandRuntimeRequest::ReadTraceGeometryStream {
                    id,
                    stream_id,
                });
            let chunk = match response.result.unwrap() {
                WireCommandRuntimeValue::TraceGeometryChunk(Some(chunk)) => chunk,
                WireCommandRuntimeValue::TraceGeometryChunkBytes(Some(bytes)) => {
                    wire::decode(&bytes).unwrap()
                }
                _ => panic!("unexpected read response"),
            };
            counts.0 += chunk.spheres.len();
            counts.1 += chunk.cylinders.len();
            counts.2 += chunk.triangles.len();
        }

        let end = state.handle_request(WireCommandRuntimeRequest::ReadTraceGeometryStream {
            id: 4,
            stream_id,
        });
        assert!(matches!(
            end.result.unwrap(),
            WireCommandRuntimeValue::TraceGeometryChunkBytes(None)
        ));
        let closed = state.handle_request(WireCommandRuntimeRequest::CloseTraceGeometryStream {
            id: 5,
            stream_id,
        });
        assert!(matches!(
            closed.result.unwrap(),
            WireCommandRuntimeValue::TraceGeometryClosed
        ));
        assert_eq!(counts, (5, 7, 9));
    }

    #[test]
    fn gpu_runtime_request_requires_declared_requirement() {
        let mut viewer = TraceFixtureViewer {
            session: Session::new(),
            chunks: Vec::new(),
        };
        let mut state = HostCommandRuntimeState::new(&mut viewer, CommandRuntimeRequirements::NONE);

        let response = state.handle_request(WireCommandRuntimeRequest::GpuDeviceLimits { id: 9 });

        assert_eq!(response.id, 9);
        assert!(response.result.is_err());
    }

    #[test]
    fn gpu_batch_request_requires_declared_requirement() {
        let mut viewer = TraceFixtureViewer {
            session: Session::new(),
            chunks: Vec::new(),
        };
        let mut state = HostCommandRuntimeState::new(&mut viewer, CommandRuntimeRequirements::NONE);

        let response = state.handle_request(WireCommandRuntimeRequest::GpuSubmitBatch {
            id: 10,
            batch: GpuSubmitBatch {
                label: Some("test.batch".to_string()),
                commands: vec![GpuBatchCommand::ReadBuffer {
                    buffer: GpuHandle {
                        id: 1,
                        kind: GpuHandleKind::Buffer,
                        generation: 1,
                    },
                    offset: 0,
                    size: 4,
                }],
            },
        });

        assert_eq!(response.id, 10);
        assert!(response.result.is_err());
    }

    #[test]
    fn gpu_cache_fingerprints_ignore_debug_labels() {
        let first = shader_module_fingerprint(&GpuShaderModuleDescriptor {
            label: Some("first".to_string()),
            wgsl: "@compute @workgroup_size(1) fn main() {}".to_string(),
        });
        let second = shader_module_fingerprint(&GpuShaderModuleDescriptor {
            label: Some("second".to_string()),
            wgsl: "@compute @workgroup_size(1) fn main() {}".to_string(),
        });
        let different_source = shader_module_fingerprint(&GpuShaderModuleDescriptor {
            label: Some("first".to_string()),
            wgsl: "@compute @workgroup_size(2) fn main() {}".to_string(),
        });

        assert_eq!(first, second);
        assert_ne!(first, different_source);
    }

    #[test]
    fn gpu_cache_keys_are_plugin_scoped_and_device_scoped() {
        let fingerprint = shader_module_fingerprint(&GpuShaderModuleDescriptor {
            label: Some("shader".to_string()),
            wgsl: "@compute @workgroup_size(1) fn main() {}".to_string(),
        });
        let mut plugin_a = GpuPluginCache::new("plugin-a@1".to_string());
        let mut plugin_b = GpuPluginCache::new("plugin-b@1".to_string());

        let first = plugin_a.resource_key_for_identity(10, fingerprint);
        let repeated = plugin_a.resource_key_for_identity(10, fingerprint);
        let other_plugin = plugin_b.resource_key_for_identity(10, fingerprint);
        let other_device = plugin_a.resource_key_for_identity(11, fingerprint);

        assert_eq!(first, repeated);
        assert_ne!(first, other_plugin);
        assert_ne!(first, other_device);
        assert_eq!(plugin_a.stats().entries, 0);

        plugin_a.record_hit();
        plugin_a.record_miss();
        let stats = plugin_a.stats();
        assert_eq!(stats.hits, 1);
        assert_eq!(stats.misses, 1);
    }

    #[test]
    fn gpu_cache_stats_and_drop_require_gpu_runtime() {
        let mut viewer = TraceFixtureViewer {
            session: Session::new(),
            chunks: Vec::new(),
        };
        let mut state = HostCommandRuntimeState::new(&mut viewer, CommandRuntimeRequirements::NONE);

        let stats = state.handle_request(WireCommandRuntimeRequest::GpuCacheStats { id: 11 });
        let drop = state.handle_request(WireCommandRuntimeRequest::GpuDropPluginCache { id: 12 });

        assert_eq!(stats.id, 11);
        assert_eq!(drop.id, 12);
        assert!(stats.result.is_err());
        assert!(drop.result.is_err());
    }

    #[test]
    fn gpu_cache_stats_and_drop_roundtrip_without_device() {
        let mut viewer = TraceFixtureViewer {
            session: Session::new(),
            chunks: Vec::new(),
        };
        let mut state =
            HostCommandRuntimeState::new(&mut viewer, CommandRuntimeRequirements::GPU_COMMANDS);

        let stats = state.handle_request(WireCommandRuntimeRequest::GpuCacheStats { id: 21 });
        match stats.result.unwrap() {
            WireCommandRuntimeValue::GpuCacheStats(stats) => {
                assert_eq!(stats.hits, 0);
                assert_eq!(stats.misses, 0);
                assert_eq!(stats.entries, 0);
            }
            _ => panic!("unexpected cache stats response"),
        }

        let dropped =
            state.handle_request(WireCommandRuntimeRequest::GpuDropPluginCache { id: 22 });
        assert!(matches!(
            dropped.result.unwrap(),
            WireCommandRuntimeValue::GpuOk
        ));
    }

    #[test]
    fn gpu_dispatch_workgroups_are_validated_before_wgpu() {
        let err = validate_dispatch_workgroups([65_536, 1, 1], 65_535).unwrap_err();

        assert!(err.contains("workgroups.x=65536"));
        assert!(validate_dispatch_workgroups([65_535, 2, 1], 65_535).is_ok());
    }

    #[test]
    fn gpu_batch_buffer_ranges_are_validated_before_wgpu() {
        let err = validate_buffer_range(16, 12, 8).unwrap_err();

        assert!(err.contains("exceeds buffer size"));
        assert!(validate_buffer_range(16, 12, 4).is_ok());
    }
}
