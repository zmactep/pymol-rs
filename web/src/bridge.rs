//! WebViewer — the main WASM-exported type.
//!
//! Mirrors the Python `StandaloneBackend` pattern: owns a `Session` +
//! `CommandExecutor`, adds GPU rendering via `RenderContext`.

use serde::Serialize;
use wasm_bindgen::prelude::*;

use pymol_cmd::{CommandExecutor, MessageKind};
use pymol_render::ShadingManager;
use pymol_scene::{CameraDelta, InputState, MoleculeObject, Object, Session, SessionAdapter};
use pymol_select::select;

use crate::event;
use crate::gpu::GpuState;
use crate::render_loop;

// ---------------------------------------------------------------------------
// Serializable types returned to JS
// ---------------------------------------------------------------------------

#[derive(Serialize)]
struct CmdOutput {
    messages: Vec<OutputMsg>,
}

#[derive(Serialize)]
struct OutputMsg {
    level: &'static str,
    text: String,
}

#[derive(Serialize)]
struct ObjectInfo {
    name: String,
    object_type: &'static str,
    atom_count: usize,
    enabled: bool,
}

#[derive(Serialize)]
struct SequenceChain {
    object_name: String,
    chain_id: String,
    residues: Vec<SequenceResidue>,
}

#[derive(Serialize)]
struct SequenceResidue {
    resn: String,
    resv: i32,
    one_letter: String,
}

#[derive(Serialize)]
struct SelectionInfo {
    name: String,
    expression: String,
    visible: bool,
}

#[derive(Serialize)]
struct MovieState {
    frame_count: usize,
    current_frame: usize,
    is_playing: bool,
}

#[derive(Serialize)]
struct LabelInfo {
    x: f32,
    y: f32,
    text: String,
    kind: &'static str,
}

// ---------------------------------------------------------------------------
// WebViewer
// ---------------------------------------------------------------------------

/// The main web viewer — owns scene state, command executor, and GPU resources.
#[wasm_bindgen]
pub struct WebViewer {
    session: Session,
    executor: CommandExecutor,
    gpu: Option<GpuState>,
    shading: ShadingManager,
    input: InputState,
    needs_redraw: bool,
    width: u32,
    height: u32,
}

#[wasm_bindgen]
impl WebViewer {
    /// Create a new WebViewer bound to a `<canvas>` element.
    ///
    /// This is async because WebGPU initialization requires awaiting the
    /// adapter and device.
    #[wasm_bindgen]
    pub async fn create(canvas_id: &str) -> Result<WebViewer, JsValue> {
        let gpu = GpuState::from_canvas(canvas_id)
            .await
            .map_err(|e| JsValue::from_str(&e))?;

        let width = gpu.surface_config.width;
        let height = gpu.surface_config.height;

        let mut session = Session::new();
        session.apply_default_settings();
        session.camera.set_aspect(width as f32 / height as f32);

        let shading = ShadingManager::new();

        Ok(WebViewer {
            session,
            executor: CommandExecutor::new(),
            gpu: Some(gpu),
            shading,
            input: InputState::new(),
            needs_redraw: true,
            width,
            height,
        })
    }

    // =======================================================================
    // Rendering
    // =======================================================================

    /// Render one frame to the canvas.
    #[wasm_bindgen]
    pub fn render_frame(&mut self) {
        let gpu = match &mut self.gpu {
            Some(g) => g,
            None => return,
        };

        if let Err(e) = render_loop::render_frame(gpu, &mut self.session, &mut self.shading) {
            log::warn!("Render error: {:?}", e);
        }
        self.needs_redraw = false;
    }

    /// Returns true when the scene has changed and needs a re-render.
    #[wasm_bindgen]
    pub fn needs_redraw(&self) -> bool {
        self.needs_redraw
    }

    /// Handle canvas resize.
    #[wasm_bindgen]
    pub fn resize(&mut self, width: u32, height: u32) {
        self.width = width.max(1);
        self.height = height.max(1);
        if let Some(gpu) = &mut self.gpu {
            gpu.resize(self.width, self.height);
        }
        self.session
            .camera
            .set_aspect(self.width as f32 / self.height as f32);
        self.needs_redraw = true;
    }

    // =======================================================================
    // Input events
    // =======================================================================

    #[wasm_bindgen]
    pub fn on_mouse_down(&mut self, x: f32, y: f32, button: u32, modifiers: u32) {
        event::handle_mouse_down(&mut self.input, x, y, button, modifiers);
        self.needs_redraw = true;
    }

    #[wasm_bindgen]
    pub fn on_mouse_move(&mut self, x: f32, y: f32, modifiers: u32) {
        event::handle_mouse_move(&mut self.input, x, y, modifiers);
        if self.input.any_button_pressed() {
            self.needs_redraw = true;
        }
    }

    #[wasm_bindgen]
    pub fn on_mouse_up(&mut self, x: f32, y: f32, button: u32) {
        event::handle_mouse_up(&mut self.input, x, y, button);
        self.needs_redraw = true;
    }

    #[wasm_bindgen]
    pub fn on_wheel(&mut self, delta_y: f32, modifiers: u32) {
        event::handle_wheel(&mut self.input, delta_y, modifiers);
        self.needs_redraw = true;
    }

    /// Process accumulated input deltas and update the camera.
    ///
    /// Call this once per frame before `render_frame()`.
    #[wasm_bindgen]
    pub fn process_input(&mut self) {
        let deltas = self.input.take_camera_deltas();
        let screen_vertex_scale = self.session.camera.screen_vertex_scale(self.height as f32);
        for delta in deltas {
            match delta {
                CameraDelta::Rotate { x, y } => {
                    self.session.camera.rotate_x(x);
                    self.session.camera.rotate_y(y);
                }
                CameraDelta::Translate(v) => {
                    let scaled = lin_alg::f32::Vec3::new(
                        v.x * screen_vertex_scale * self.input.pan_sensitivity,
                        v.y * screen_vertex_scale * self.input.pan_sensitivity,
                        v.z * screen_vertex_scale * self.input.pan_sensitivity,
                    );
                    self.session.camera.translate(scaled);
                }
                CameraDelta::Zoom(z) => {
                    self.session.camera.zoom(z);
                }
                CameraDelta::Clip { front, back } => {
                    let view = self.session.camera.view_mut();
                    view.clip_front = (view.clip_front + front).max(0.01);
                    view.clip_back = (view.clip_back + back).max(view.clip_front + 0.01);
                }
                CameraDelta::SlabScale(raw_delta) => {
                    let mws = self
                        .session
                        .settings
                        .get_float(pymol_settings::id::mouse_wheel_scale);
                    let scale = 1.0 + 0.04 * mws * raw_delta;

                    let view = self.session.camera.view_mut();
                    let avg = (view.clip_front + view.clip_back) * 0.5;
                    let half_width = (view.clip_back - avg).max(0.1);
                    let new_half = (half_width * scale).max(0.1);

                    view.clip_front = (avg - new_half).max(0.01);
                    view.clip_back = (avg + new_half).max(view.clip_front + 0.1);
                }
            }
            self.needs_redraw = true;
        }
    }

    // =======================================================================
    // Commands
    // =======================================================================

    /// Execute a PyMOL command string. Returns JSON with output messages.
    #[wasm_bindgen]
    pub fn execute(&mut self, command: &str) -> JsValue {
        let render_context = self.gpu.as_ref().map(|g| &g.render_context);
        let mut adapter = SessionAdapter {
            session: &mut self.session,
            render_context,
            default_size: (self.width, self.height),
            needs_redraw: &mut self.needs_redraw,
            async_fetch_fn: None,
        };

        let result = self.executor.do_with_options(&mut adapter, command, false);

        let messages = match result {
            Ok(output) => output
                .messages
                .iter()
                .map(|m| OutputMsg {
                    level: match m.kind {
                        MessageKind::Info => "info",
                        MessageKind::Warning => "warning",
                        MessageKind::Error => "error",
                    },
                    text: m.text.clone(),
                })
                .collect(),
            Err(e) => vec![OutputMsg {
                level: "error",
                text: e.to_string(),
            }],
        };

        let output = CmdOutput { messages };
        serde_wasm_bindgen::to_value(&output).unwrap_or(JsValue::NULL)
    }

    /// Load molecular or map data from bytes.
    ///
    /// `format` should be one of: "pdb", "xyz", "cif", "mmcif", "bcif",
    /// "ccp4", "map", "mrc"
    ///
    /// Gzip-compressed data is automatically detected and decompressed.
    #[wasm_bindgen]
    pub fn load_data(&mut self, data: &[u8], name: &str, format: &str) -> Result<(), JsValue> {
        // Decompress gzip if detected (magic bytes 0x1f 0x8b)
        let decompressed;
        let data = if data.len() >= 2 && data[0] == 0x1f && data[1] == 0x8b {
            use std::io::Read;
            let mut decoder = pymol_io::compress::gzip_reader(data);
            let mut buf = Vec::new();
            decoder.read_to_end(&mut buf).map_err(|e| {
                JsValue::from_str(&format!("Gzip decompression failed: {}", e))
            })?;
            decompressed = buf;
            decompressed.as_slice()
        } else {
            data
        };

        let fmt = format.to_lowercase();

        // CCP4/MRC density maps — binary format, returns MapObject
        if matches!(fmt.as_str(), "ccp4" | "map" | "mrc") {
            let ccp4 = pymol_io::ccp4::read_ccp4_from(std::io::Cursor::new(data))
                .map_err(|e| JsValue::from_str(&format!("CCP4 parse error: {}", e)))?;
            let grid = pymol_render::Grid3D::from_dims(
                ccp4.origin, ccp4.spacing, ccp4.dims, ccp4.values,
            );
            let map_data = pymol_scene::MapData::new(grid);
            let map_obj = pymol_scene::MapObject::from_map_data(name, map_data);
            self.session.registry.add(map_obj);
            if let Some((min, max)) = self.session.registry.extent() {
                self.session.camera.zoom_to(min, max, 0.0);
            }
            self.needs_redraw = true;
            return Ok(());
        }

        let mol = match fmt.as_str() {
            // Binary formats — parse directly from bytes
            "bcif" => pymol_io::bcif::read_bcif_bytes(data),
            // Text formats — require UTF-8
            _ => {
                let data_str = std::str::from_utf8(data)
                    .map_err(|_| JsValue::from_str("Data is not valid UTF-8"))?;
                match fmt.as_str() {
                    "pdb" => pymol_io::pdb::read_pdb_str(data_str),
                    "xyz" => pymol_io::xyz::read_xyz_str(data_str),
                    "cif" | "mmcif" => pymol_io::cif::read_cif_str(data_str),
                    _ => {
                        return Err(JsValue::from_str(&format!(
                            "Direct loading not yet supported for: {}. Use execute() instead.",
                            fmt
                        )))
                    }
                }
            }
        };

        match mol {
            Ok(molecule) => {
                let mol_obj = MoleculeObject::with_name(molecule, name);
                self.session.registry.add(mol_obj);
                // Zoom to fit
                if let Some((min, max)) = self.session.registry.extent() {
                    self.session.camera.zoom_to(min, max, 0.0);
                }
                self.needs_redraw = true;
                Ok(())
            }
            Err(e) => Err(JsValue::from_str(&format!("Parse error: {}", e))),
        }
    }

    // =======================================================================
    // State queries
    // =======================================================================

    /// Get loaded object names as a JSON array.
    #[wasm_bindgen]
    pub fn get_object_names(&self) -> JsValue {
        let names: Vec<String> = self
            .session
            .registry
            .names()
            .map(|s| s.to_string())
            .collect();
        serde_wasm_bindgen::to_value(&names).unwrap_or(JsValue::NULL)
    }

    /// Get info about a specific object.
    #[wasm_bindgen]
    pub fn get_object_info(&self, name: &str) -> JsValue {
        if let Some(mol_obj) = self.session.registry.get_molecule(name) {
            let info = ObjectInfo {
                name: name.to_string(),
                object_type: "molecule",
                atom_count: mol_obj.molecule().atom_count(),
                enabled: mol_obj.is_enabled(),
            };
            return serde_wasm_bindgen::to_value(&info).unwrap_or(JsValue::NULL);
        }
        if let Some(map_obj) = self.session.registry.get_map(name) {
            let info = ObjectInfo {
                name: name.to_string(),
                object_type: "map",
                atom_count: 0,
                enabled: map_obj.is_enabled(),
            };
            return serde_wasm_bindgen::to_value(&info).unwrap_or(JsValue::NULL);
        }
        JsValue::NULL
    }

    /// Get sequence data for all loaded molecules as JSON.
    #[wasm_bindgen]
    pub fn get_sequence_data(&self) -> JsValue {
        let mut chains: Vec<SequenceChain> = Vec::new();

        for name in self.session.registry.names() {
            if let Some(mol_obj) = self.session.registry.get_molecule(name) {
                let mol = mol_obj.molecule();
                let mut chain_map: std::collections::BTreeMap<String, Vec<SequenceResidue>> =
                    std::collections::BTreeMap::new();

                let mut last_resv: Option<(String, i32)> = None;
                for atom in mol.atoms() {
                    let chain = atom.residue.key.chain.clone();
                    let resv = atom.residue.key.resv;
                    if last_resv.as_ref() == Some(&(chain.clone(), resv)) {
                        continue;
                    }
                    last_resv = Some((chain.clone(), resv));

                    // Skip water and ions (like the desktop sequence builder)
                    if pymol_mol::is_water(&atom.residue.resn)
                        || pymol_mol::is_ion(&atom.residue.resn)
                    {
                        continue;
                    }

                    let one_letter = pymol_mol::residue_to_char(&atom.residue.resn).to_string();

                    chain_map
                        .entry(chain.clone())
                        .or_default()
                        .push(SequenceResidue {
                            resn: atom.residue.resn.clone(),
                            resv,
                            one_letter,
                        });
                }

                for (chain_id, residues) in chain_map {
                    chains.push(SequenceChain {
                        object_name: name.to_string(),
                        chain_id,
                        residues,
                    });
                }
            }
        }

        serde_wasm_bindgen::to_value(&chains).unwrap_or(JsValue::NULL)
    }

    /// Get current movie state as JSON.
    #[wasm_bindgen]
    pub fn get_movie_state(&self) -> JsValue {
        let state = MovieState {
            frame_count: self.session.movie.frame_count(),
            current_frame: self.session.movie.current_frame(),
            is_playing: self.session.movie.is_playing(),
        };
        serde_wasm_bindgen::to_value(&state).unwrap_or(JsValue::NULL)
    }

    /// Get all named selections as JSON array.
    #[wasm_bindgen]
    pub fn get_selection_list(&self) -> JsValue {
        let mut list: Vec<SelectionInfo> = self
            .session
            .selections
            .iter()
            .map(|(name, entry)| SelectionInfo {
                name: name.clone(),
                expression: entry.expression.clone(),
                visible: entry.visible,
            })
            .collect();
        list.sort_by(|a, b| a.name.cmp(&b.name));
        serde_wasm_bindgen::to_value(&list).unwrap_or(JsValue::NULL)
    }

    /// Count atoms matching a selection expression.
    #[wasm_bindgen]
    pub fn count_atoms(&self, selection: &str) -> Result<usize, JsValue> {
        let mut total = 0;
        for name in self.session.registry.names() {
            if let Some(mol_obj) = self.session.registry.get_molecule(name) {
                let mol = mol_obj.molecule();
                match select(mol, selection) {
                    Ok(mask) => total += mask.count(),
                    Err(e) => return Err(JsValue::from_str(&format!("Selection error: {}", e))),
                }
            }
        }
        Ok(total)
    }

    /// Get projected screen-space labels for overlay rendering.
    ///
    /// Returns a JSON array of `{ x, y, text, kind }` where coordinates
    /// are in physical pixels (divide by `devicePixelRatio` for CSS pixels).
    #[wasm_bindgen]
    pub fn get_labels(&self) -> JsValue {
        let viewport = (0.0, 0.0, self.width as f32, self.height as f32);
        let mut labels = Vec::new();

        for name in self.session.registry.names() {
            if let Some(mol_obj) = self.session.registry.get_molecule(name) {
                if !mol_obj.is_enabled() {
                    continue;
                }
                for (pos, text) in mol_obj.collect_labels() {
                    if let Some((sx, sy)) = self.session.camera.project_to_screen(pos, viewport) {
                        if sx >= 0.0
                            && sx <= self.width as f32
                            && sy >= 0.0
                            && sy <= self.height as f32
                        {
                            labels.push(LabelInfo {
                                x: sx,
                                y: sy,
                                text: text.to_string(),
                                kind: "atom",
                            });
                        }
                    }
                }
            }
        }

        for name in self.session.registry.names() {
            if let Some(meas_obj) = self.session.registry.get_measurement(name) {
                if !meas_obj.is_enabled() {
                    continue;
                }
                for (pos, text) in meas_obj.collect_labels() {
                    if let Some((sx, sy)) = self.session.camera.project_to_screen(pos, viewport) {
                        if sx >= 0.0
                            && sx <= self.width as f32
                            && sy >= 0.0
                            && sy <= self.height as f32
                        {
                            labels.push(LabelInfo {
                                x: sx,
                                y: sy,
                                text: text.to_string(),
                                kind: "measurement",
                            });
                        }
                    }
                }
            }
        }

        serde_wasm_bindgen::to_value(&labels).unwrap_or(JsValue::NULL)
    }
}
