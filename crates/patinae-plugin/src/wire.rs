//! Versioned runtime wire payloads.
//!
//! Dynamic plugins exchange rich runtime data with the host as MessagePack
//! bytes. This module owns those DTOs so the FFI boundary stays limited to
//! pointer-length byte views and opaque handles.

use std::sync::{Arc, RwLock};

use patinae_cmd::{
    CommandAction, DynamicCommandInvocation, DynamicSettingRegistry, OutputMessage, ParsedCommand,
};
use patinae_framework::atom_stream::{AtomChunk, AtomStreamRequest};
use patinae_framework::message::AppMessage;
use patinae_framework::plugin_ui::{PanelAction, PanelEvent, PanelSnapshot};
use patinae_mol::ObjectMolecule;
use patinae_render::{
    DisplayedGeometry, DisplayedMaterial, DisplayedMesh, DisplayedMeshVertex,
    DisplayedObjectGeometry, DisplayedPrimitive, ObjectId, RepKind,
};
use patinae_scene::{Camera, MovieStateSnapshot, SceneView, Session, ViewportImage};
use patinae_settings::{
    DynamicSettingDescriptor, DynamicSettingStore, SettingType, SettingValue, SideEffectCategory,
};
use serde::{de::DeserializeOwned, Deserialize, Serialize};

/// Runtime wire version for MessagePack DTOs.
pub const RUNTIME_WIRE_VERSION: u32 = 3;

/// Maximum MessagePack payload copied across the runtime ABI.
pub const MAX_WIRE_PAYLOAD_LEN: usize = 64 * 1024 * 1024;

/// Serialized dynamic setting state.
#[derive(Clone, Serialize, Deserialize)]
pub struct WireDynamicSetting {
    /// Setting descriptor metadata.
    pub descriptor: WireDynamicSettingDescriptor,
    /// Explicit global value, if one has been set.
    pub value: Option<SettingValue>,
}

/// Serialized dynamic setting descriptor.
#[derive(Clone, Serialize, Deserialize)]
pub struct WireDynamicSettingDescriptor {
    /// Setting name.
    pub name: String,
    /// Value type.
    pub setting_type: SettingType,
    /// Default value.
    pub default: SettingValue,
    /// Optional minimum value.
    pub min: Option<f32>,
    /// Optional maximum value.
    pub max: Option<f32>,
    /// Named value hints.
    pub value_hints: Vec<(String, SettingValue)>,
    /// Side effects for host invalidation.
    pub side_effects: Vec<SideEffectCategory>,
    /// Whether object-specific values are accepted.
    pub object_overridable: bool,
}

/// Lightweight identity for a viewport image.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct WireViewportImageSummary {
    /// Image width in pixels.
    pub width: u32,
    /// Image height in pixels.
    pub height: u32,
    /// Image byte length.
    pub len: usize,
    /// Host-computed identity for change detection.
    pub signature: u64,
}

/// Lightweight host state used by per-frame poll callbacks.
#[derive(Clone, Serialize, Deserialize)]
pub struct WirePollSharedInput {
    /// Must equal [`RUNTIME_WIRE_VERSION`].
    pub wire_version: u32,
    /// Scene generation counter.
    pub scene_generation: u64,
    /// Loaded object names.
    pub object_names: Vec<String>,
    /// Current camera state.
    pub camera: Camera,
    /// Current movie state.
    pub movie: MovieStateSnapshot,
    /// Current global settings.
    pub settings: patinae_settings::Settings,
    /// Current background color.
    pub clear_color: [f32; 3],
    /// Viewport image identity without pixel bytes.
    pub viewport_image: Option<WireViewportImageSummary>,
    /// Known command names.
    pub command_names: Vec<String>,
    /// Known setting names.
    pub setting_names: Vec<String>,
    /// Host-owned dynamic setting state.
    pub dynamic_settings: Vec<WireDynamicSetting>,
}

/// Renderer-neutral displayed geometry for runtime command input.
#[derive(Clone, Serialize, Deserialize)]
pub struct WireDisplayedGeometry {
    /// Display objects in render order.
    pub objects: Vec<WireDisplayedObjectGeometry>,
}

/// Displayed geometry for one renderer object.
#[derive(Clone, Serialize, Deserialize)]
pub struct WireDisplayedObjectGeometry {
    /// Host renderer object id.
    pub object_id: u32,
    /// Displayed primitives for this object.
    pub primitives: Vec<WireDisplayedPrimitive>,
}

/// A displayed primitive in wire form.
#[derive(Clone, Serialize, Deserialize)]
pub enum WireDisplayedPrimitive {
    /// Non-indexed triangle-list mesh.
    Mesh {
        /// Raw [`RepKind`] tag.
        rep: u8,
        /// Mesh vertices.
        vertices: Vec<WireDisplayedMeshVertex>,
    },
    /// Analytic sphere.
    Sphere {
        /// Raw [`RepKind`] tag.
        rep: u8,
        /// Owning atom id.
        owner_atom_id: u32,
        /// Sphere center.
        center: [f32; 3],
        /// Sphere radius.
        radius: f32,
        /// Resolved material.
        material: WireDisplayedMaterial,
    },
    /// Analytic cylinder.
    Cylinder {
        /// Raw [`RepKind`] tag.
        rep: u8,
        /// Owning atom ids.
        owner_atom_ids: [u32; 2],
        /// Cylinder start.
        start: [f32; 3],
        /// Cylinder end.
        end: [f32; 3],
        /// Cylinder radius.
        radius: f32,
        /// Material at the start.
        material_start: WireDisplayedMaterial,
        /// Material at the end.
        material_end: WireDisplayedMaterial,
    },
    /// Semantic screen-space line.
    LineSegment {
        /// Raw [`RepKind`] tag.
        rep: u8,
        /// Owning atom ids.
        owner_atom_ids: [u32; 2],
        /// Line start.
        start: [f32; 3],
        /// Line end.
        end: [f32; 3],
        /// Line width in pixels.
        width_px: f32,
        /// Material at the start.
        material_start: WireDisplayedMaterial,
        /// Material at the end.
        material_end: WireDisplayedMaterial,
    },
    /// Semantic screen-space point.
    PointSample {
        /// Raw [`RepKind`] tag.
        rep: u8,
        /// Owning atom id.
        owner_atom_id: u32,
        /// Point position.
        position: [f32; 3],
        /// Point radius in pixels.
        radius_px: f32,
        /// Resolved material.
        material: WireDisplayedMaterial,
    },
}

/// Displayed mesh vertex in wire form.
#[derive(Clone, Serialize, Deserialize)]
pub struct WireDisplayedMeshVertex {
    /// World position.
    pub position: [f32; 3],
    /// World normal.
    pub normal: [f32; 3],
    /// Owning atom id.
    pub owner_atom_id: u32,
    /// Resolved material.
    pub material: WireDisplayedMaterial,
    /// Renderer flags.
    pub flags: u32,
}

/// Display material in wire form.
#[derive(Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct WireDisplayedMaterial {
    /// Host-resolved base color.
    pub base_rgba: [f32; 4],
    /// Representation color.
    pub rep_rgba: [f32; 4],
    /// Final display color.
    pub rgba: [f32; 4],
    /// PyMOL-style transparency.
    pub transparency: f32,
}

/// Portable atom property value.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum WireAtomPropertyValue {
    /// String value.
    Str(String),
    /// 32-bit float value.
    F32(f32),
    /// 32-bit integer value.
    I32(i32),
    /// 8-bit integer value.
    I8(i8),
    /// Boolean value.
    Bool(bool),
}

/// Portable atom property change.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WireAtomPropertyChange {
    /// Object name.
    pub object: String,
    /// Zero-based atom index.
    pub atom_index: u32,
    /// Property changes.
    pub changes: Vec<(String, WireAtomPropertyValue)>,
}

/// Portable viewer-side actions emitted by dynamic plugin polling.
#[derive(Clone, Serialize, Deserialize)]
pub enum WireViewerAction {
    /// Apply atom property changes.
    ApplyAtomPropertyChanges(Vec<WireAtomPropertyChange>),
    /// Store a viewport overlay image.
    SetViewportImage(ViewportImage),
    /// Clear the viewport overlay image.
    ClearViewportImage,
    /// Request a viewport redraw.
    RequestRedraw,
    /// Invalidate plugin panel snapshots.
    RequestPanelUpdate,
}

/// Host-side query requested by a dynamic plugin.
#[derive(Clone, Serialize, Deserialize)]
pub enum WireHostQuery {
    /// Return loaded object names.
    ObjectNames { id: u64 },
    /// Return the current camera view.
    View { id: u64 },
    /// Count atoms matching a selection expression.
    CountAtoms { id: u64, selection: String },
    /// Return full viewport image bytes.
    ViewportImage { id: u64 },
    /// Open a host-owned atom stream.
    OpenAtomStream {
        /// Correlation id.
        id: u64,
        /// Stream request.
        request: AtomStreamRequest,
    },
    /// Read rows from a host-owned atom stream.
    ReadAtomStream {
        /// Correlation id.
        id: u64,
        /// Host-owned stream id.
        stream_id: u64,
        /// Maximum rows to return.
        max_rows: usize,
    },
    /// Close a host-owned atom stream.
    CloseAtomStream {
        /// Correlation id.
        id: u64,
        /// Host-owned stream id.
        stream_id: u64,
    },
}

/// Host-side query result delivered on a later poll.
#[derive(Clone, Serialize, Deserialize)]
pub struct WireHostQueryResult {
    /// Correlation id.
    pub id: u64,
    /// Query value or portable error.
    pub result: Result<WireHostQueryValue, String>,
}

/// Host-side query value.
#[derive(Clone, Serialize, Deserialize)]
pub enum WireHostQueryValue {
    /// Loaded object names.
    ObjectNames(Vec<String>),
    /// Current camera view.
    View(SceneView),
    /// Atom count.
    CountAtoms(usize),
    /// Full viewport image bytes.
    ViewportImage(Option<ViewportImage>),
    /// Open atom stream metadata.
    AtomStreamOpened(WireAtomStreamOpened),
    /// Atom stream rows.
    AtomStreamChunk(AtomChunk),
    /// Atom stream closed.
    AtomStreamClosed,
}

/// Open atom stream metadata.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct WireAtomStreamOpened {
    /// Host-owned stream id.
    pub stream_id: u64,
    /// Total rows in the stream.
    pub total_count: usize,
}

impl From<&DynamicSettingDescriptor> for WireDynamicSettingDescriptor {
    fn from(value: &DynamicSettingDescriptor) -> Self {
        Self {
            name: value.name.clone(),
            setting_type: value.setting_type,
            default: value.default.clone(),
            min: value.min,
            max: value.max,
            value_hints: value.value_hints.clone(),
            side_effects: value.side_effects.clone(),
            object_overridable: value.object_overridable,
        }
    }
}

impl From<WireDynamicSettingDescriptor> for DynamicSettingDescriptor {
    fn from(value: WireDynamicSettingDescriptor) -> Self {
        Self {
            name: value.name,
            setting_type: value.setting_type,
            default: value.default,
            min: value.min,
            max: value.max,
            value_hints: value.value_hints,
            side_effects: value.side_effects,
            object_overridable: value.object_overridable,
        }
    }
}

/// Command execution input.
#[derive(Clone, Serialize, Deserialize)]
pub struct WireCommandInput {
    /// Must equal [`RUNTIME_WIRE_VERSION`].
    pub wire_version: u32,
    /// MessagePack-encoded [`Session`].
    pub session: Vec<u8>,
    /// Current viewport image snapshot.
    pub viewport_image: Option<ViewportImage>,
    /// Parsed command invocation.
    pub parsed: ParsedCommand,
    /// Whether user-visible output should be suppressed.
    pub quiet: bool,
    /// Current viewport width.
    pub viewport_width: u32,
    /// Current viewport height.
    pub viewport_height: u32,
    /// Host-owned dynamic setting state.
    pub dynamic_settings: Vec<WireDynamicSetting>,
    /// Host-resolved displayed geometry, when requested by the command.
    pub displayed_geometry: Option<WireDisplayedGeometry>,
}

/// Command execution output.
#[derive(Clone, Serialize, Deserialize)]
pub struct WireCommandOutput {
    /// Must equal [`RUNTIME_WIRE_VERSION`].
    pub wire_version: u32,
    /// Command result, using a string error for portability.
    pub result: Result<(), String>,
    /// Messages printed by the command.
    pub output: Vec<OutputMessage>,
    /// Host actions requested by the command.
    pub actions: Vec<CommandAction>,
    /// MessagePack-encoded updated [`Session`].
    pub session: Vec<u8>,
    /// Viewport image set by the plugin, if any.
    pub viewport_image: Option<ViewportImage>,
}

/// Shared host state input for panel and poll callbacks.
#[derive(Clone, Serialize, Deserialize)]
pub struct WireSharedInput {
    /// Must equal [`RUNTIME_WIRE_VERSION`].
    pub wire_version: u32,
    /// MessagePack-encoded [`Session`].
    pub session: Vec<u8>,
    /// Scene generation counter.
    pub scene_generation: u64,
    /// Viewport image snapshot.
    pub viewport_image: Option<ViewportImage>,
    /// Known command names.
    pub command_names: Vec<String>,
    /// Known setting names.
    pub setting_names: Vec<String>,
    /// Host-owned dynamic setting state.
    pub dynamic_settings: Vec<WireDynamicSetting>,
}

/// Panel event input.
#[derive(Clone, Serialize, Deserialize)]
pub struct WirePanelEventInput {
    /// Shared host state snapshot.
    pub shared: WireSharedInput,
    /// Frontend event.
    pub event: PanelEvent,
}

/// Panel event output.
#[derive(Clone, Serialize, Deserialize)]
pub struct WirePanelEventOutput {
    /// Must equal [`RUNTIME_WIRE_VERSION`].
    pub wire_version: u32,
    /// Panel actions emitted by the plugin.
    pub actions: Vec<PanelAction>,
}

/// Panel snapshot output.
#[derive(Clone, Serialize, Deserialize)]
pub struct WirePanelSnapshotOutput {
    /// Must equal [`RUNTIME_WIRE_VERSION`].
    pub wire_version: u32,
    /// Snapshot to render.
    pub snapshot: PanelSnapshot,
}

/// Command execution result delivered during polling.
#[derive(Clone, Serialize, Deserialize)]
pub struct WireCommandResult {
    /// Correlation ID.
    pub id: u64,
    /// Result of the command.
    pub result: Result<(), String>,
}

/// Deferred command execution request.
#[derive(Clone, Serialize, Deserialize)]
pub struct WireCommandExecRequest {
    /// Correlation ID.
    pub id: u64,
    /// Command text.
    pub command: String,
    /// Whether execution should suppress echo/output.
    pub silent: bool,
}

/// Dynamic command registration request.
#[derive(Clone, Serialize, Deserialize)]
pub struct WireDynCmdRegistration {
    /// Command name.
    pub name: String,
    /// Short description.
    pub description: String,
    /// Usage text.
    pub usage: String,
    /// Argument help text.
    pub arguments: String,
}

/// Portable hotkey action.
#[derive(Clone, Serialize, Deserialize)]
pub enum WireHotkeyAction {
    /// Execute a command string.
    Command(String),
    /// Invoke a dynamic command.
    DynamicCommand { name: String, args: Vec<String> },
    /// Publish a custom message.
    Custom { topic: String, payload: Vec<u8> },
}

/// Portable hotkey registration.
#[derive(Clone, Serialize, Deserialize)]
pub struct WireHotkeyRegistration {
    /// Key string to parse on the host side.
    pub key: String,
    /// Host-owned action.
    pub action: WireHotkeyAction,
}

/// Message handler poll input.
#[derive(Clone, Serialize, Deserialize)]
pub struct WirePollInput {
    /// Lightweight shared host state.
    pub shared: WirePollSharedInput,
    /// Results from previously queued commands.
    pub command_results: Vec<WireCommandResult>,
    /// Results from previously queued host queries.
    pub host_query_results: Vec<WireHostQueryResult>,
    /// Dynamic command invocations since the last poll.
    pub dynamic_invocations: Vec<DynamicCommandInvocation>,
    /// Plugin directories as UTF-8 paths.
    pub plugin_dirs: Vec<String>,
}

/// Message handler poll output.
#[derive(Clone, Serialize, Deserialize)]
pub struct WirePollOutput {
    /// Must equal [`RUNTIME_WIRE_VERSION`].
    pub wire_version: u32,
    /// Bus messages emitted by the plugin.
    pub messages: Vec<AppMessage>,
    /// Deferred command executions.
    pub command_exec: Vec<WireCommandExecRequest>,
    /// Dynamic command registrations.
    pub dynamic_registrations: Vec<WireDynCmdRegistration>,
    /// Dynamic command unregistrations.
    pub dynamic_unregistrations: Vec<String>,
    /// User notifications.
    pub notifications: Vec<String>,
    /// Hotkey registrations.
    pub hotkey_registrations: Vec<WireHotkeyRegistration>,
    /// Hotkey unregistrations.
    pub hotkey_unregistrations: Vec<String>,
    /// Portable host queries requested by the plugin.
    pub host_queries: Vec<WireHostQuery>,
    /// Portable viewer actions requested by the plugin.
    pub viewer_actions: Vec<WireViewerAction>,
}

/// Message delivery input.
#[derive(Clone, Serialize, Deserialize)]
pub struct WireMessageInput {
    /// Must equal [`RUNTIME_WIRE_VERSION`].
    pub wire_version: u32,
    /// Broadcast message.
    pub message: AppMessage,
}

/// Message delivery output.
#[derive(Clone, Serialize, Deserialize)]
pub struct WireMessageOutput {
    /// Must equal [`RUNTIME_WIRE_VERSION`].
    pub wire_version: u32,
    /// Bus messages emitted by the plugin.
    pub messages: Vec<AppMessage>,
}

/// Script handler input.
#[derive(Clone, Serialize, Deserialize)]
pub struct WireScriptInput {
    /// Must equal [`RUNTIME_WIRE_VERSION`].
    pub wire_version: u32,
    /// Script path.
    pub path: String,
}

/// Script handler output.
#[derive(Clone, Serialize, Deserialize)]
pub struct WireScriptOutput {
    /// Must equal [`RUNTIME_WIRE_VERSION`].
    pub wire_version: u32,
    /// Script dispatch result.
    pub result: Result<(), String>,
}

/// File format reader input.
#[derive(Clone, Serialize, Deserialize)]
pub struct WireFormatReadInput {
    /// Must equal [`RUNTIME_WIRE_VERSION`].
    pub wire_version: u32,
    /// File extension without a dot.
    pub extension: String,
    /// Full file bytes read by the host.
    pub bytes: Vec<u8>,
}

/// File format reader output.
#[derive(Clone, Serialize, Deserialize)]
pub struct WireFormatReadOutput {
    /// Must equal [`RUNTIME_WIRE_VERSION`].
    pub wire_version: u32,
    /// Parsed molecules or a portable error.
    pub result: Result<Vec<ObjectMolecule>, String>,
}

/// File format writer input.
#[derive(Clone, Serialize, Deserialize)]
pub struct WireFormatWriteInput {
    /// Must equal [`RUNTIME_WIRE_VERSION`].
    pub wire_version: u32,
    /// File extension without a dot.
    pub extension: String,
    /// Molecules selected by the host.
    pub molecules: Vec<ObjectMolecule>,
}

/// File format writer output.
#[derive(Clone, Serialize, Deserialize)]
pub struct WireFormatWriteOutput {
    /// Must equal [`RUNTIME_WIRE_VERSION`].
    pub wire_version: u32,
    /// File bytes to write or a portable error.
    pub result: Result<Vec<u8>, String>,
}

/// Converts displayed geometry into wire DTOs.
pub fn displayed_geometry_to_wire(displayed: &DisplayedGeometry) -> WireDisplayedGeometry {
    WireDisplayedGeometry {
        objects: displayed
            .objects
            .iter()
            .map(|object| WireDisplayedObjectGeometry {
                object_id: object.object_id.0,
                primitives: object
                    .primitives
                    .iter()
                    .map(displayed_primitive_to_wire)
                    .collect(),
            })
            .collect(),
    }
}

/// Converts wire DTOs into displayed geometry.
pub fn displayed_geometry_from_wire(displayed: WireDisplayedGeometry) -> DisplayedGeometry {
    DisplayedGeometry {
        objects: displayed
            .objects
            .into_iter()
            .map(|object| DisplayedObjectGeometry {
                object_id: ObjectId(object.object_id),
                primitives: object
                    .primitives
                    .into_iter()
                    .map(displayed_primitive_from_wire)
                    .collect(),
            })
            .collect(),
    }
}

fn displayed_primitive_to_wire(primitive: &DisplayedPrimitive) -> WireDisplayedPrimitive {
    match primitive {
        DisplayedPrimitive::Mesh { rep, mesh } => WireDisplayedPrimitive::Mesh {
            rep: rep.as_raw(),
            vertices: mesh
                .vertices
                .iter()
                .map(displayed_mesh_vertex_to_wire)
                .collect(),
        },
        DisplayedPrimitive::Sphere {
            rep,
            owner_atom_id,
            center,
            radius,
            material,
        } => WireDisplayedPrimitive::Sphere {
            rep: rep.as_raw(),
            owner_atom_id: *owner_atom_id,
            center: *center,
            radius: *radius,
            material: displayed_material_to_wire(*material),
        },
        DisplayedPrimitive::Cylinder {
            rep,
            owner_atom_ids,
            start,
            end,
            radius,
            material_start,
            material_end,
        } => WireDisplayedPrimitive::Cylinder {
            rep: rep.as_raw(),
            owner_atom_ids: *owner_atom_ids,
            start: *start,
            end: *end,
            radius: *radius,
            material_start: displayed_material_to_wire(*material_start),
            material_end: displayed_material_to_wire(*material_end),
        },
        DisplayedPrimitive::LineSegment {
            rep,
            owner_atom_ids,
            start,
            end,
            width_px,
            material_start,
            material_end,
        } => WireDisplayedPrimitive::LineSegment {
            rep: rep.as_raw(),
            owner_atom_ids: *owner_atom_ids,
            start: *start,
            end: *end,
            width_px: *width_px,
            material_start: displayed_material_to_wire(*material_start),
            material_end: displayed_material_to_wire(*material_end),
        },
        DisplayedPrimitive::PointSample {
            rep,
            owner_atom_id,
            position,
            radius_px,
            material,
        } => WireDisplayedPrimitive::PointSample {
            rep: rep.as_raw(),
            owner_atom_id: *owner_atom_id,
            position: *position,
            radius_px: *radius_px,
            material: displayed_material_to_wire(*material),
        },
    }
}

fn displayed_primitive_from_wire(primitive: WireDisplayedPrimitive) -> DisplayedPrimitive {
    match primitive {
        WireDisplayedPrimitive::Mesh { rep, vertices } => DisplayedPrimitive::Mesh {
            rep: RepKind::from_raw(rep),
            mesh: DisplayedMesh {
                vertices: vertices
                    .into_iter()
                    .map(displayed_mesh_vertex_from_wire)
                    .collect(),
            },
        },
        WireDisplayedPrimitive::Sphere {
            rep,
            owner_atom_id,
            center,
            radius,
            material,
        } => DisplayedPrimitive::Sphere {
            rep: RepKind::from_raw(rep),
            owner_atom_id,
            center,
            radius,
            material: displayed_material_from_wire(material),
        },
        WireDisplayedPrimitive::Cylinder {
            rep,
            owner_atom_ids,
            start,
            end,
            radius,
            material_start,
            material_end,
        } => DisplayedPrimitive::Cylinder {
            rep: RepKind::from_raw(rep),
            owner_atom_ids,
            start,
            end,
            radius,
            material_start: displayed_material_from_wire(material_start),
            material_end: displayed_material_from_wire(material_end),
        },
        WireDisplayedPrimitive::LineSegment {
            rep,
            owner_atom_ids,
            start,
            end,
            width_px,
            material_start,
            material_end,
        } => DisplayedPrimitive::LineSegment {
            rep: RepKind::from_raw(rep),
            owner_atom_ids,
            start,
            end,
            width_px,
            material_start: displayed_material_from_wire(material_start),
            material_end: displayed_material_from_wire(material_end),
        },
        WireDisplayedPrimitive::PointSample {
            rep,
            owner_atom_id,
            position,
            radius_px,
            material,
        } => DisplayedPrimitive::PointSample {
            rep: RepKind::from_raw(rep),
            owner_atom_id,
            position,
            radius_px,
            material: displayed_material_from_wire(material),
        },
    }
}

fn displayed_mesh_vertex_to_wire(vertex: &DisplayedMeshVertex) -> WireDisplayedMeshVertex {
    WireDisplayedMeshVertex {
        position: vertex.position,
        normal: vertex.normal,
        owner_atom_id: vertex.owner_atom_id,
        material: displayed_material_to_wire(vertex.material),
        flags: vertex.flags,
    }
}

fn displayed_mesh_vertex_from_wire(vertex: WireDisplayedMeshVertex) -> DisplayedMeshVertex {
    DisplayedMeshVertex {
        position: vertex.position,
        normal: vertex.normal,
        owner_atom_id: vertex.owner_atom_id,
        material: displayed_material_from_wire(vertex.material),
        flags: vertex.flags,
    }
}

fn displayed_material_to_wire(material: DisplayedMaterial) -> WireDisplayedMaterial {
    WireDisplayedMaterial {
        base_rgba: material.base_rgba,
        rep_rgba: material.rep_rgba,
        rgba: material.rgba,
        transparency: material.transparency,
    }
}

fn displayed_material_from_wire(material: WireDisplayedMaterial) -> DisplayedMaterial {
    DisplayedMaterial {
        base_rgba: material.base_rgba,
        rep_rgba: material.rep_rgba,
        rgba: material.rgba,
        transparency: material.transparency,
    }
}

/// Serialize a wire value to MessagePack bytes.
///
/// # Errors
///
/// Returns a string if MessagePack serialization fails.
pub fn encode<T: Serialize>(value: &T) -> Result<Vec<u8>, String> {
    rmp_serde::to_vec(value).map_err(|error| error.to_string())
}

/// Deserialize a wire value from MessagePack bytes.
///
/// # Errors
///
/// Returns a string when the payload is too large or malformed.
pub fn decode<T: DeserializeOwned>(bytes: &[u8]) -> Result<T, String> {
    if bytes.len() > MAX_WIRE_PAYLOAD_LEN {
        return Err("wire payload exceeds ABI limit".to_string());
    }
    rmp_serde::from_slice(bytes).map_err(|error| error.to_string())
}

/// Encode a session to MessagePack bytes.
///
/// # Errors
///
/// Returns a string if session serialization fails.
pub fn encode_session(session: &Session) -> Result<Vec<u8>, String> {
    encode(session)
}

/// Decode a session from MessagePack bytes.
///
/// # Errors
///
/// Returns a string if the payload is too large or malformed.
pub fn decode_session(bytes: &[u8]) -> Result<Session, String> {
    decode(bytes)
}

/// Build a dynamic setting registry from wire entries.
///
/// # Errors
///
/// Returns a string if a descriptor collides with an existing setting.
pub fn dynamic_registry_from_wire(
    entries: &[WireDynamicSetting],
) -> Result<DynamicSettingRegistry, String> {
    let mut registry = DynamicSettingRegistry::new();
    for entry in entries {
        let descriptor: DynamicSettingDescriptor = entry.descriptor.clone().into();
        let mut store = DynamicSettingStore::new();
        if let Some(value) = &entry.value {
            store.set(&descriptor.name, value.clone());
        }
        registry.register(descriptor, Arc::new(RwLock::new(store)))?;
    }
    Ok(registry)
}
