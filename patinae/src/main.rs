mod app;
mod bridges;
mod clipboard;
mod components;
mod display_recovery;
mod fetch;
#[cfg(target_os = "macos")]
mod macos;
mod native_file_actions;
mod native_menu;
mod recent_files;
mod recent_thumbnails;

slint::include_modules!();

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info"))
        .filter_module("wgpu_core", log::LevelFilter::Warn)
        .filter_module("wgpu_hal", log::LevelFilter::Warn)
        .filter_module("naga", log::LevelFilter::Warn)
        .init();

    log::info!("Starting Patinae");
    app::run()
}
