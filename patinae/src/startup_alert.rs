use slint::ComponentHandle;

use crate::StartupErrorWindow;

pub(crate) fn show_error(title: &str, message: &str) {
    if let Err(err) = show_error_with_slint(title, message) {
        log::warn!("Failed to show startup alert dialog with Slint software renderer: {err}");
        eprintln!("{title}: {message}");
    }
}

fn show_error_with_slint(title: &str, message: &str) -> Result<(), Box<dyn std::error::Error>> {
    slint::BackendSelector::new()
        .backend_name("winit".into())
        .renderer_name("software".into())
        .select()?;

    let window = StartupErrorWindow::new()?;
    window.set_dialog_title(title.into());
    window.set_message_text(message.into());
    let weak_window = window.as_weak();
    window.on_acknowledged(move || {
        if let Some(window) = weak_window.upgrade() {
            let _ = window.hide();
        }
        let _ = slint::quit_event_loop();
    });
    window.run()?;
    Ok(())
}
