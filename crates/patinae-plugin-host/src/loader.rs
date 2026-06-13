use std::ffi::c_void;
use std::io::{Read, Write};
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::Path;
use std::sync::{Arc, RwLock};

use libloading::Library;

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
    HostRegistrarHandle, PluginCommandHandle, PluginDeclaration, PluginFormatHandlerHandle,
    PluginHotkeyHandle, PluginMessageHandlerHandle, PluginPanelHandle, PluginScriptHandlerHandle,
    ABI_STATUS_HOST_ERROR, ABI_STATUS_INVALID, ABI_STATUS_PANIC, ABI_STATUS_UNSUPPORTED,
    ABI_VERSION, ARG_HINT_COLOR, ARG_HINT_COMMAND, ARG_HINT_LABEL_PROPERTY,
    ARG_HINT_NAMED_SELECTION, ARG_HINT_OBJECT, ARG_HINT_PATH, ARG_HINT_REPRESENTATION,
    ARG_HINT_SELECTION, ARG_HINT_SETTING, ARG_HINT_SETTING_VALUE, CAPABILITY_REGISTRATION,
    HOST_CALLBACKS_VERSION, KNOWN_CAPABILITIES, KNOWN_COMMAND_RUNTIME_REQUIREMENTS,
    KNOWN_PANEL_RUNTIME_REQUIREMENTS, MAX_ABI_SLICE_LEN, MAX_ABI_STRING_LEN,
    PANEL_PLACEMENT_BOTTOM, PANEL_PLACEMENT_RIGHT, SDK_VERSION, SETTING_TYPE_BLANK,
    SETTING_TYPE_BOOL, SETTING_TYPE_COLOR, SETTING_TYPE_FLOAT, SETTING_TYPE_FLOAT3,
    SETTING_TYPE_INT, SETTING_TYPE_STRING, SETTING_VALUE_BOOL, SETTING_VALUE_COLOR,
    SETTING_VALUE_FLOAT, SETTING_VALUE_FLOAT3, SETTING_VALUE_INT, SETTING_VALUE_NONE,
    SETTING_VALUE_STRING, SIDE_EFFECT_COLOR_REBUILD, SIDE_EFFECT_FULL_REBUILD,
    SIDE_EFFECT_ORTHO_DIRTY, SIDE_EFFECT_REPRESENTATION_REBUILD, SIDE_EFFECT_SCENE_CHANGED,
    SIDE_EFFECT_SCENE_INVALIDATE, SIDE_EFFECT_SEQ_CHANGED, SIDE_EFFECT_SHADER_COMPUTE_LIGHTING,
    SIDE_EFFECT_SHADER_RELOAD, SIDE_EFFECT_STEREO_UPDATE, SIDE_EFFECT_SURFACE_TRANSPARENCY,
    SIDE_EFFECT_VIEWPORT_UPDATE,
};
use patinae_plugin::registrar::{MessageHandler, PluginKeyAction, PluginMetadata, PollContext};
use patinae_plugin::wire::{
    self, WireCommandInput, WireCommandOutput, WireFormatReadInput, WireFormatReadOutput,
    WireFormatWriteInput, WireFormatWriteOutput, WireHotkeyAction, WireMessageInput,
    WireMessageOutput, WirePanelEventInput, WirePanelEventOutput, WirePanelSnapshotOutput,
    WirePollInput, WirePollOutput, WireScriptInput, WireScriptOutput, WireSharedInput,
    WireViewerAction, RUNTIME_WIRE_VERSION,
};
use patinae_scene::{parse_key_string, KeyBindings, Session, ViewerLike};
use patinae_settings::{
    DynamicSettingDescriptor, DynamicSettingStore, SettingType, SettingValue, SharedSettingStore,
    SideEffectCategory,
};
use serde::{de::DeserializeOwned, Serialize};

use crate::host::PluginHost;
use crate::panic::panic_payload_to_string;
use crate::paths::{is_plugin_library_path, PluginDiscovery};
use crate::plugin::{LibraryHandle, LoadedPanel, LoadedPlugin};

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
    let library = Arc::new(library_handle);

    install_registration_assets(executor, &mut registration, library.clone());

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
) {
    for command in std::mem::take(&mut registration.commands) {
        executor
            .registry_mut()
            .register_boxed(Box::new(AbiCommandProxy::new(command, library.clone())));
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
}

impl AbiCommandProxy {
    fn new(command: RegisteredCommand, library: Arc<LibraryHandle>) -> Self {
        Self {
            command,
            _library: library,
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
        let input = command_input_from_context(ctx, args, self.command.runtime_requirements)
            .map_err(|error| CmdError::execution(format!("plugin command input: {error}")))?;
        let execute = self
            .command
            .vtable
            .execute
            .ok_or_else(|| CmdError::execution("plugin command execute callback was null"))?;
        let output: WireCommandOutput =
            invoke_runtime("command execute", &input, |slice, sink, user_data| {
                // SAFETY: The descriptor was validated during registration and the
                // library is kept alive by this proxy while the callback runs.
                unsafe { execute(self.command.handle, slice, sink, user_data) }
            })
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

pub(crate) fn command_input_from_context<V: ViewerLike + ?Sized>(
    ctx: &mut CommandContext<'_, '_, V>,
    args: &ParsedCommand,
    runtime_requirements: CommandRuntimeRequirements,
) -> Result<WireCommandInput, String> {
    let (viewport_width, viewport_height) = ctx.viewer.viewport_size();
    let displayed_geometry =
        if runtime_requirements.contains(CommandRuntimeRequirements::DISPLAYED_GEOMETRY) {
            let displayed = ctx
                .viewer
                .export_displayed_geometry(&patinae_render::GeometryExportOptions::default())?;
            Some(wire::displayed_geometry_to_wire(&displayed))
        } else {
            None
        };
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
    })
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
