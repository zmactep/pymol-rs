//! PyMOL-RS Web Viewer
//!
//! WASM module providing molecular visualization via WebGPU.
//! Exposes a `WebViewer` struct to JavaScript through wasm-bindgen.

mod bridge;
mod event;
mod gpu;
mod picking;
mod render_loop;

pub use bridge::WebViewer;

use wasm_bindgen::prelude::*;

/// Initialize the WASM module: panic hook + console logger.
///
/// Called automatically before `WebViewer::create()`, but can also
/// be called explicitly from JS if early logging is needed.
#[wasm_bindgen(start)]
pub fn init() {
    console_error_panic_hook::set_once();
    console_log::init_with_level(log::Level::Info).ok();
}
