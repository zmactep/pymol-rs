//! Native Menu Bar
//!
//! Creates a platform-native menu bar using `muda`. On macOS this appears
//! in the system menu bar; on Windows/Linux it attaches to the window.

use muda::{
    AboutMetadata, CheckMenuItem, Icon, Menu, MenuId, MenuItemKind, MenuItem,
    PredefinedMenuItem, Submenu,
    accelerator::{Accelerator, Code, Modifiers},
};

/// IDs for custom menu items that trigger application commands.
pub struct MenuIds {
    // File
    pub run: MenuId,
    pub open: MenuId,
    pub fetch: MenuId,
    pub save: MenuId,
    pub export_png: MenuId,
    pub export_movie: MenuId,

    // Edit
    pub select_all: MenuId,
    pub deselect_all: MenuId,

    // View
    pub reset_view: MenuId,
    pub zoom_all: MenuId,
    pub orient: MenuId,
    pub center: MenuId,
    pub opaque_background: MenuId,
    pub bg_white: MenuId,
    pub bg_black: MenuId,
    pub transparent_panels: MenuId,
    pub panels_submenu: MenuId,

    // Help
    pub help_commands: MenuId,
}

/// Application state that the menu needs to reflect.
/// Set values before calling `AppMenu::sync()`.
pub struct AppMenuState {
    pub transparent_panels: bool,
    pub opaque_background: bool,
    pub bg_color: [f32; 3],
    /// (component_id, title, visible) for each layout panel.
    pub panels: Vec<(String, String, bool)>,
}

/// Built menu and the IDs needed to dispatch events.
pub struct AppMenu {
    pub menu: Menu,
    pub ids: MenuIds,
    pub state: AppMenuState,
}

/// Load the application icon from the embedded PNG.
fn load_app_icon() -> Option<Icon> {
    let png_bytes = include_bytes!("../../../images/pymol-rs.png");
    let img = image::load_from_memory(png_bytes).ok()?.into_rgba8();
    let (width, height) = img.dimensions();
    Icon::from_rgba(img.into_raw(), width, height).ok()
}

/// Build the full application menu bar.
pub fn build_menu() -> AppMenu {
    let menu = Menu::new();

    // ── App menu (macOS only — first submenu becomes the app menu) ───────
    #[cfg(target_os = "macos")]
    {
        let app_menu = Submenu::new("PyMOL-RS", true);
        let _ = app_menu.append_items(&[
            &PredefinedMenuItem::about(
                None,
                Some(AboutMetadata {
                    icon: load_app_icon(),
                    name: Some("PyMOL-RS".into()),
                    version: Some(env!("CARGO_PKG_VERSION").into()),
                    copyright: Some("2026, Pavel Yakovlev".into()),
                    license: Some("BSD-3-Clause".into()),
                    website: Some("https://zmactep.github.io/pymol-rs".into()),
                    ..Default::default()
                }),
            ),
            &PredefinedMenuItem::separator(),
            &PredefinedMenuItem::services(None),
            &PredefinedMenuItem::separator(),
            &PredefinedMenuItem::hide(None),
            &PredefinedMenuItem::hide_others(None),
            &PredefinedMenuItem::show_all(None),
            &PredefinedMenuItem::separator(),
            &PredefinedMenuItem::quit(None),
        ]);
        let _ = menu.append(&app_menu);
    }

    // ── File ─────────────────────────────────────────────────────────────
    let run = MenuItem::new(
        "Run Script…",
        true,
        Some(Accelerator::new(Some(Modifiers::META), Code::KeyR)),
    );
    let open = MenuItem::new(
        "Open…",
        true,
        Some(Accelerator::new(Some(Modifiers::META), Code::KeyO)),
    );
    let fetch = MenuItem::new(
        "Fetch PDB…",
        true,
        Some(Accelerator::new(Some(Modifiers::META), Code::KeyF)),
    );
    let save = MenuItem::new(
        "Save…",
        true,
        Some(Accelerator::new(Some(Modifiers::META), Code::KeyS)),
    );
    let export_png = MenuItem::new(
        "Export PNG…",
        true,
        Some(Accelerator::new(
            Some(Modifiers::META | Modifiers::SHIFT),
            Code::KeyE,
        )),
    );
    let export_movie = MenuItem::new("Export Movie…", true, None);

    let file_menu = Submenu::new("File", true);
    let _ = file_menu.append_items(&[
        &run,
        &PredefinedMenuItem::separator(),
        &open,
        &fetch,
        &PredefinedMenuItem::separator(),
        &save,
        &export_png,
        &export_movie,
        &PredefinedMenuItem::separator(),
        #[cfg(not(target_os = "macos"))]
        &PredefinedMenuItem::quit(None),
    ]);
    let _ = menu.append(&file_menu);

    // ── Edit ─────────────────────────────────────────────────────────────
    let select_all = MenuItem::new(
        "Select All Atoms",
        true,
        Some(Accelerator::new(
            Some(Modifiers::META | Modifiers::SHIFT),
            Code::KeyA,
        )),
    );
    let deselect_all = MenuItem::new("Deselect All", true, None);

    let edit_menu = Submenu::new("Edit", true);
    let _ = edit_menu.append_items(&[
        &PredefinedMenuItem::undo(None),
        &PredefinedMenuItem::redo(None),
        &PredefinedMenuItem::separator(),
        &select_all,
        &deselect_all,
    ]);
    let _ = menu.append(&edit_menu);

    // ── View ─────────────────────────────────────────────────────────────
    let reset_view = MenuItem::new("Reset", true, None);
    let zoom_all = MenuItem::new(
        "Zoom All",
        true,
        Some(Accelerator::new(Some(Modifiers::META), Code::Digit0)),
    );
    let orient = MenuItem::new("Orient", true, None);
    let center = MenuItem::new("Center", true, None);

    let opaque_background = CheckMenuItem::new("Opaque", true, false, None);
    let bg_white = CheckMenuItem::new("White", true, false, None);
    let bg_black = CheckMenuItem::new("Black", true, false, None);

    let bg_submenu = Submenu::new("Background", true);
    let _ = bg_submenu.append_items(&[
        &opaque_background,
        &PredefinedMenuItem::separator(),
        &bg_white,
        &bg_black,
    ]);

    let transparent_panels = CheckMenuItem::new("Transparent Panels", true, false, None);

    let panels_submenu = Submenu::with_id("panels_submenu", "Panels", true);

    let view_menu = Submenu::new("View", true);
    let _ = view_menu.append_items(&[
        &reset_view,
        &zoom_all,
        &orient,
        &center,
        &PredefinedMenuItem::separator(),
        &bg_submenu,
        &PredefinedMenuItem::separator(),
        &transparent_panels,
        &panels_submenu,
        &PredefinedMenuItem::separator(),
        &PredefinedMenuItem::fullscreen(None)
    ]);
    let _ = menu.append(&view_menu);

    // ── Window (macOS standard) ──────────────────────────────────────────
    #[cfg(target_os = "macos")]
    {
        let window_menu = Submenu::new("Window", true);
        let _ = window_menu.append_items(&[
            &PredefinedMenuItem::minimize(None),
            &PredefinedMenuItem::maximize(None),
            &PredefinedMenuItem::separator(),
            &PredefinedMenuItem::close_window(None),
            &PredefinedMenuItem::separator(),
            &PredefinedMenuItem::bring_all_to_front(None),
        ]);
        let _ = menu.append(&window_menu);
    }

    // ── Help ─────────────────────────────────────────────────────────────
    let help_commands = MenuItem::new("Command Reference", true, None);

    let help_menu = Submenu::new("Help", true);
    let _ = help_menu.append_items(&[&help_commands]);
    let _ = menu.append(&help_menu);

    // ── Collect IDs ──────────────────────────────────────────────────────
    let ids = MenuIds {
        run: run.id().clone(),
        open: open.id().clone(),
        fetch: fetch.id().clone(),
        save: save.id().clone(),
        export_png: export_png.id().clone(),
        export_movie: export_movie.id().clone(),
        select_all: select_all.id().clone(),
        deselect_all: deselect_all.id().clone(),
        reset_view: reset_view.id().clone(),
        zoom_all: zoom_all.id().clone(),
        orient: orient.id().clone(),
        center: center.id().clone(),
        opaque_background: opaque_background.id().clone(),
        bg_white: bg_white.id().clone(),
        bg_black: bg_black.id().clone(),
        transparent_panels: transparent_panels.id().clone(),
        panels_submenu: panels_submenu.id().clone(),
        help_commands: help_commands.id().clone(),
    };

    let state = AppMenuState {
        transparent_panels: false,
        opaque_background: false,
        bg_color: [0.0, 0.0, 0.0],
        panels: Vec::new(),
    };

    AppMenu {
        menu,
        ids,
        state
    }
}

/// Find a `CheckMenuItem` by ID, recursing into submenus.
fn find_check_item(items: &[MenuItemKind], id: &MenuId) -> Option<CheckMenuItem> {
    for item in items {
        if item.id() == id {
            return item.as_check_menuitem().cloned();
        }
        if let Some(sub) = item.as_submenu() {
            if let Some(found) = find_check_item(&sub.items(), id) {
                return Some(found);
            }
        }
    }
    None
}

/// Find a `Submenu` by ID, recursing into submenus.
fn find_submenu(items: &[MenuItemKind], id: &MenuId) -> Option<Submenu> {
    for item in items {
        if let Some(sub) = item.as_submenu() {
            if sub.id() == id {
                return Some(sub.clone());
            }
            if let Some(found) = find_submenu(&sub.items(), id) {
                return Some(found);
            }
        }
    }
    None
}

impl AppMenu {
    /// Sync menu check marks with `self.state`.
    /// Call after updating `state` fields.
    pub fn sync(&mut self) {
        let items = self.menu.items();

        // Transparent panels
        if let Some(item) = find_check_item(&items, &self.ids.transparent_panels) {
            if item.is_checked() != self.state.transparent_panels {
                item.set_checked(self.state.transparent_panels);
            }
        }

        // Opaque background
        if let Some(item) = find_check_item(&items, &self.ids.opaque_background) {
            if item.is_checked() != self.state.opaque_background {
                item.set_checked(self.state.opaque_background);
            }
        }

        // Background color
        let cc = self.state.bg_color;
        let is_white = cc[0] >= 0.99 && cc[1] >= 0.99 && cc[2] >= 0.99;
        let is_black = cc[0] <= 0.01 && cc[1] <= 0.01 && cc[2] <= 0.01;
        if let Some(item) = find_check_item(&items, &self.ids.bg_white) {
            if item.is_checked() != is_white {
                item.set_checked(is_white);
            }
        }
        if let Some(item) = find_check_item(&items, &self.ids.bg_black) {
            if item.is_checked() != is_black {
                item.set_checked(is_black);
            }
        }

        // Panels submenu — rebuild when panel list changes
        if let Some(submenu) = find_submenu(&items, &self.ids.panels_submenu) {
            let current_items = submenu.items();
            let current_ids: Vec<_> = current_items
                .iter()
                .filter_map(|item| {
                    let check = item.as_check_menuitem()?;
                    Some(check.id().clone())
                })
                .collect();

            let desired_ids: Vec<_> = self
                .state
                .panels
                .iter()
                .map(|(id, _, _)| MenuId::new(id.as_str()))
                .collect();

            let needs_rebuild = current_ids != desired_ids;

            if needs_rebuild {
                // Remove all existing items (backwards to keep indices valid)
                let count = current_items.len();
                for i in (0..count).rev() {
                    submenu.remove_at(i);
                }
                // Add fresh items
                for (id, title, visible) in &self.state.panels {
                    let item = CheckMenuItem::with_id(
                        MenuId::new(id),
                        title,
                        true,
                        *visible,
                        None,
                    );
                    let _ = submenu.append(&item);
                }
            } else {
                // Just sync check states
                for item_kind in submenu.items() {
                    if let Some(check) = item_kind.as_check_menuitem() {
                        let check_id = check.id().clone();
                        if let Some((_, _, vis)) = self
                            .state
                            .panels
                            .iter()
                            .find(|(id, _, _)| MenuId::new(id.as_str()) == check_id)
                        {
                            if check.is_checked() != *vis {
                                check.set_checked(*vis);
                            }
                        }
                    }
                }
            }
        }
    }

    /// If a menu event ID matches a panel toggle item, return the component ID.
    pub fn panel_toggle_id(&self, menu_id: &MenuId) -> Option<String> {
        self.state
            .panels
            .iter()
            .find(|(id, _, _)| MenuId::new(id.as_str()) == *menu_id)
            .map(|(id, _, _)| id.clone())
    }
}
