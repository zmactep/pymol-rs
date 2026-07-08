//! CoreVideo render heartbeat for macOS.

use std::ffi::{c_int, c_void};
use std::ptr;

use slint::ComponentHandle;

use crate::AppWindow;

type CVReturn = i32;
type CVOptionFlags = u64;
type CVDisplayLinkRef = *mut c_void;
type CGDirectDisplayID = u32;
type CGDisplayChangeSummaryFlags = u32;

const SUCCESS: i32 = 0;

type CVDisplayLinkOutputCallback = extern "C" fn(
    display_link: CVDisplayLinkRef,
    now: *const c_void,
    output_time: *const c_void,
    flags_in: CVOptionFlags,
    flags_out: *mut CVOptionFlags,
    display_link_context: *mut c_void,
) -> CVReturn;

type CGDisplayReconfigurationCallBack = extern "C" fn(
    display: CGDirectDisplayID,
    flags: CGDisplayChangeSummaryFlags,
    user_info: *mut c_void,
);

#[link(name = "CoreVideo", kind = "framework")]
extern "C" {
    fn CVDisplayLinkCreateWithActiveCGDisplays(display_link_out: *mut CVDisplayLinkRef)
        -> CVReturn;
    fn CVDisplayLinkSetOutputCallback(
        display_link: CVDisplayLinkRef,
        callback: CVDisplayLinkOutputCallback,
        user_info: *mut c_void,
    ) -> CVReturn;
    fn CVDisplayLinkStart(display_link: CVDisplayLinkRef) -> CVReturn;
    fn CVDisplayLinkStop(display_link: CVDisplayLinkRef) -> CVReturn;
}

#[link(name = "CoreGraphics", kind = "framework")]
extern "C" {
    fn CGDisplayRegisterReconfigurationCallback(
        callback: CGDisplayReconfigurationCallBack,
        user_info: *mut c_void,
    ) -> c_int;
    fn CGDisplayRemoveReconfigurationCallback(
        callback: CGDisplayReconfigurationCallBack,
        user_info: *mut c_void,
    ) -> c_int;
}

#[link(name = "CoreFoundation", kind = "framework")]
extern "C" {
    fn CFRelease(cf: *const c_void);
}

/// Requests redraws from a CoreVideo display-link callback.
///
/// macOS can stop delivering Slint/winit redraws in clamshell external-display
/// transitions. This heartbeat stays vsync-aligned and queues redraw requests
/// onto Slint's event loop so Patinae still renders through Slint's normal
/// render notifier.
pub struct DisplayLinkHeartbeat {
    display_link: CVDisplayLinkRef,
    state: *mut HeartbeatState,
    reconfiguration_callback_registered: bool,
}

struct HeartbeatState {
    window: slint::Weak<AppWindow>,
}

impl DisplayLinkHeartbeat {
    /// Starts a display-link heartbeat for a Slint window.
    pub fn start(window: slint::Weak<AppWindow>) -> Result<Self, String> {
        let mut display_link = ptr::null_mut();
        let create_status = unsafe {
            // SAFETY: `display_link` points to writable storage for the
            // CoreVideo-created display-link reference.
            CVDisplayLinkCreateWithActiveCGDisplays(&mut display_link)
        };
        if create_status != SUCCESS {
            return Err(format!(
                "CVDisplayLinkCreateWithActiveCGDisplays failed with status {create_status}"
            ));
        }
        if display_link.is_null() {
            return Err("CVDisplayLinkCreateWithActiveCGDisplays returned null".to_string());
        }

        let state = Box::into_raw(Box::new(HeartbeatState { window }));
        let callback_status = unsafe {
            // SAFETY: `display_link` is a valid CoreVideo object from the
            // create call above, and `state` remains owned by this heartbeat
            // until the display link is stopped in `Drop`.
            CVDisplayLinkSetOutputCallback(display_link, display_link_callback, state.cast())
        };
        if callback_status != SUCCESS {
            unsafe {
                // SAFETY: `state` has not been given to a running display link,
                // and `display_link` is the owned object from CoreVideo.
                drop(Box::from_raw(state));
                CFRelease(display_link.cast_const());
            }
            return Err(format!(
                "CVDisplayLinkSetOutputCallback failed with status {callback_status}"
            ));
        }

        let register_status = unsafe {
            // SAFETY: The callback only dereferences the `state` pointer while
            // this heartbeat is alive; `Drop` unregisters it before freeing.
            CGDisplayRegisterReconfigurationCallback(display_reconfiguration_callback, state.cast())
        };
        let reconfiguration_callback_registered = register_status == SUCCESS;
        if !reconfiguration_callback_registered {
            log::warn!(
                "CGDisplayRegisterReconfigurationCallback failed with status {register_status}"
            );
        }

        let start_status = unsafe {
            // SAFETY: `display_link` is a valid CoreVideo object with a valid
            // output callback and context.
            CVDisplayLinkStart(display_link)
        };
        if start_status != SUCCESS {
            unsafe {
                // SAFETY: The display link never started; unregister any
                // registered display callback, release CoreVideo ownership,
                // and reclaim the boxed callback state.
                if reconfiguration_callback_registered {
                    CGDisplayRemoveReconfigurationCallback(
                        display_reconfiguration_callback,
                        state.cast(),
                    );
                }
                drop(Box::from_raw(state));
                CFRelease(display_link.cast_const());
            }
            return Err(format!(
                "CVDisplayLinkStart failed with status {start_status}"
            ));
        }

        Ok(Self {
            display_link,
            state,
            reconfiguration_callback_registered,
        })
    }
}

impl Drop for DisplayLinkHeartbeat {
    fn drop(&mut self) {
        unsafe {
            // SAFETY: `display_link` and `state` are owned by this heartbeat.
            // Stopping and unregistering callbacks before freeing the state
            // prevents CoreVideo/CoreGraphics from observing a dangling pointer.
            CVDisplayLinkStop(self.display_link);
            if self.reconfiguration_callback_registered {
                CGDisplayRemoveReconfigurationCallback(
                    display_reconfiguration_callback,
                    self.state.cast(),
                );
            }
            CFRelease(self.display_link.cast_const());
            drop(Box::from_raw(self.state));
        }
    }
}

extern "C" fn display_link_callback(
    _display_link: CVDisplayLinkRef,
    _now: *const c_void,
    _output_time: *const c_void,
    _flags_in: CVOptionFlags,
    _flags_out: *mut CVOptionFlags,
    display_link_context: *mut c_void,
) -> CVReturn {
    request_redraw(display_link_context);
    SUCCESS
}

extern "C" fn display_reconfiguration_callback(
    _display: CGDirectDisplayID,
    _flags: CGDisplayChangeSummaryFlags,
    user_info: *mut c_void,
) {
    request_redraw(user_info);
}

fn request_redraw(context: *mut c_void) {
    if context.is_null() {
        return;
    }

    let state = unsafe {
        // SAFETY: CoreVideo/CoreGraphics pass the context pointer registered by
        // `DisplayLinkHeartbeat::start`; it stays valid until callbacks are
        // stopped/unregistered in `Drop`.
        &*(context.cast::<HeartbeatState>())
    };
    let window = state.window.clone();
    let _ = window.upgrade_in_event_loop(|app| {
        app.window().request_redraw();
    });
}
