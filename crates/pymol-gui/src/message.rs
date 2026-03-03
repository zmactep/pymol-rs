//! Unified message bus for inter-component communication
//!
//! All GUI components (command system, object list, viewport, sequence viewer,
//! movie player, IPC, plugins) communicate through a single [`AppMessage`] type
//! and a per-frame [`MessageBus`] queue.
//!
//! The centralized dispatch in [`App::dispatch`] uses exhaustive matching,
//! so the compiler ensures every message variant is handled.

use std::collections::HashSet;

use pymol_mol::RepMask;

use crate::model::sequence::ResidueRef;

/// Unified application message.
///
/// Produced by UI components, IPC handlers, and async tasks.
/// Consumed by the central dispatcher after each egui frame.
#[derive(Debug, Clone)]
pub enum AppMessage {
    // =====================================================================
    // Command Execution
    // =====================================================================
    /// Execute a PyMOL command string (from command line, sequence viewer,
    /// viewport click, IPC, or any other source).
    ExecuteCommand { command: String, silent: bool },

    // =====================================================================
    // Object Management
    // =====================================================================
    /// Toggle an object's enabled/disabled state.
    ToggleObject(String),
    /// Enable all objects.
    EnableAllObjects,
    /// Disable all objects.
    DisableAllObjects,
    /// Delete an object by name.
    DeleteObject(String),

    // =====================================================================
    // Representations
    // =====================================================================
    /// Show a representation on an object.
    ShowRepresentation { object: String, rep: RepMask },
    /// Hide a representation on an object.
    HideRepresentation { object: String, rep: RepMask },
    /// Show all representations on an object.
    ShowAllRepresentations(String),
    /// Hide all representations on an object.
    HideAllRepresentations(String),

    // =====================================================================
    // Coloring
    // =====================================================================
    /// Set the color of an object by color name.
    SetColor { object: String, color: String },

    // =====================================================================
    // Camera / Navigation
    // =====================================================================
    /// Zoom to an object.
    ZoomTo(String),
    /// Center on an object.
    CenterOn(String),

    // =====================================================================
    // Selections
    // =====================================================================
    /// Delete a named selection.
    DeleteSelection(String),
    /// Toggle a selection's visibility indicator.
    ToggleSelectionVisibility(String),
    /// A selection was modified (triggers cross-component sync).
    SelectionChanged { name: String },

    // =====================================================================
    // Output / Notifications
    // =====================================================================
    /// Print an informational message to the output log.
    PrintInfo(String),
    /// Print a warning message to the output log.
    PrintWarning(String),
    /// Print an error message to the output log.
    PrintError(String),
    /// Echo a command to the output log (green "PyMOL> ..." line).
    PrintCommand(String),

    // =====================================================================
    // Viewport Interaction
    // =====================================================================
    /// Hover state changed in the sequence viewer.
    HoverResidue(Option<ResidueRef>),

    // =====================================================================
    // Movie / Animation
    // =====================================================================
    /// Play the movie.
    MoviePlay,
    /// Pause the movie.
    MoviePause,
    /// Stop the movie (reset to frame 0).
    MovieStop,
    /// Jump to a specific frame.
    MovieGotoFrame(usize),
    /// Set playback FPS.
    MovieSetFps(f32),
    /// Set loop mode (0=once, 1=loop, 2=swing).
    MovieSetLoopMode(pymol_scene::LoopMode),

    // =====================================================================
    // System
    // =====================================================================
    /// Request a window redraw.
    RequestRedraw,
    /// Request application quit.
    Quit,
    /// Focus a panel by id.
    FocusPanel(String),

    // =====================================================================
    // Component System
    // =====================================================================
    /// Update sequence highlights (from selection sync).
    UpdateHighlights(HashSet<ResidueRef>),
    /// Toggle a panel's expanded/collapsed state by component ID.
    TogglePanel(String),
    /// Detach a docked panel into a floating window.
    FloatPanel(String),
    /// Dock a floating panel back to its original slot.
    DockPanel(String),
    /// Switch the active tab in a tabbed panel group.
    ActivateTab(String),

    // =====================================================================
    // Extensibility (Plugins)
    // =====================================================================
    /// Custom event for plugins.
    Custom { topic: String, payload: Vec<u8> },
}

/// Per-frame message queue with two-phase buffering.
///
/// Components write to the outbox during the egui frame. After the frame,
/// the app drains the outbox and dispatches messages. On the next frame,
/// drained messages become the inbox (readable by components if needed).
///
/// # Lifecycle per frame
///
/// 1. `begin_frame()` — moves outbox → inbox
/// 2. Components read inbox, write to outbox via `send()`
/// 3. `drain_outbox()` — app processes new messages
pub struct MessageBus {
    inbox: Vec<AppMessage>,
    outbox: Vec<AppMessage>,
}

impl MessageBus {
    /// Create an empty bus.
    pub fn new() -> Self {
        Self {
            inbox: Vec::with_capacity(8),
            outbox: Vec::with_capacity(8),
        }
    }

    /// Start a new frame: moves outbox → inbox for component reading.
    pub fn begin_frame(&mut self) {
        self.inbox.clear();
        std::mem::swap(&mut self.inbox, &mut self.outbox);
    }

    /// Send a message (writes to outbox).
    #[inline]
    pub fn send(&mut self, msg: AppMessage) {
        self.outbox.push(msg);
    }

    /// Read the inbox (messages from the previous frame).
    pub fn inbox(&self) -> &[AppMessage] {
        &self.inbox
    }

    /// Drain all outbox messages for dispatch.
    pub fn drain_outbox(&mut self) -> Vec<AppMessage> {
        std::mem::take(&mut self.outbox)
    }

    /// Inject a system event directly into the outbox.
    pub fn inject(&mut self, msg: AppMessage) {
        self.outbox.push(msg);
    }

    /// Whether the outbox has pending messages.
    pub fn has_pending(&self) -> bool {
        !self.outbox.is_empty()
    }

    // =================================================================
    // Convenience methods
    // =================================================================

    /// Send an `ExecuteCommand` message.
    pub fn execute_command(&mut self, cmd: impl Into<String>) {
        self.send(AppMessage::ExecuteCommand {
            command: cmd.into(),
            silent: false,
        });
    }

    /// Send a `PrintInfo` message.
    pub fn print_info(&mut self, msg: impl Into<String>) {
        self.send(AppMessage::PrintInfo(msg.into()));
    }

    /// Send a `PrintWarning` message.
    pub fn print_warning(&mut self, msg: impl Into<String>) {
        self.send(AppMessage::PrintWarning(msg.into()));
    }

    /// Send a `PrintError` message.
    pub fn print_error(&mut self, msg: impl Into<String>) {
        self.send(AppMessage::PrintError(msg.into()));
    }

    /// Send a `RequestRedraw` message.
    pub fn request_redraw(&mut self) {
        self.send(AppMessage::RequestRedraw);
    }
}

impl Default for MessageBus {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_two_phase_buffering() {
        let mut bus = MessageBus::new();

        // Frame 1: send a message
        bus.send(AppMessage::RequestRedraw);
        assert!(bus.has_pending());
        assert!(bus.inbox().is_empty());

        // Drain outbox (app dispatch)
        let messages = bus.drain_outbox();
        assert_eq!(messages.len(), 1);
        assert!(!bus.has_pending());

        // Frame 2: begin_frame moves nothing (outbox was drained)
        bus.begin_frame();
        assert!(bus.inbox().is_empty());
    }

    #[test]
    fn test_inbox_carries_over() {
        let mut bus = MessageBus::new();

        // Send messages but don't drain
        bus.send(AppMessage::PrintInfo("hello".into()));
        bus.send(AppMessage::RequestRedraw);

        // Begin next frame: outbox → inbox
        bus.begin_frame();
        assert_eq!(bus.inbox().len(), 2);
        assert!(!bus.has_pending());
    }

    #[test]
    fn test_convenience_methods() {
        let mut bus = MessageBus::new();

        bus.execute_command("load test.pdb");
        bus.print_info("loaded");
        bus.request_redraw();

        let messages = bus.drain_outbox();
        assert_eq!(messages.len(), 3);

        assert!(matches!(&messages[0], AppMessage::ExecuteCommand { command, silent: false } if command == "load test.pdb"));
        assert!(matches!(&messages[1], AppMessage::PrintInfo(msg) if msg == "loaded"));
        assert!(matches!(&messages[2], AppMessage::RequestRedraw));
    }
}
