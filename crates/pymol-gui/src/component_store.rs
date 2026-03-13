//! Component Store
//!
//! Dynamic registry of UI components, keyed by component ID.

use pymol_framework::component::{Component, EguiComponent, SharedContext};
use pymol_framework::message::{AppMessage, MessageBus};

/// Dynamic registry of UI components.
///
/// Components are stored in registration order and accessed either
/// generically (iteration, layout dispatch) or by concrete type via
/// [`get`](Self::get) / [`get_mut`](Self::get_mut) with downcasting.
pub struct ComponentStore {
    components: Vec<Box<dyn EguiComponent>>,
}

impl ComponentStore {
    pub fn new() -> Self {
        Self {
            components: Vec::new(),
        }
    }

    /// Register a component. Panics on duplicate IDs.
    pub fn add(&mut self, component: impl EguiComponent + 'static) {
        self.add_boxed(Box::new(component));
    }

    /// Register a boxed component. Panics on duplicate IDs.
    pub fn add_boxed(&mut self, component: Box<dyn EguiComponent>) {
        let id = component.id();
        assert!(
            !self.components.iter().any(|c| c.id() == id),
            "duplicate component ID: {id}"
        );
        self.components.push(component);
    }

    /// Downcast to a concrete component type (shared ref).
    pub fn get<T: Component + 'static>(&self) -> Option<&T> {
        self.components.iter().find_map(|c| {
            let any: &dyn std::any::Any = c.as_ref();
            any.downcast_ref::<T>()
        })
    }

    /// Downcast to a concrete component type (mutable ref).
    pub fn get_mut<T: Component + 'static>(&mut self) -> Option<&mut T> {
        self.components.iter_mut().find_map(|c| {
            let any: &mut dyn std::any::Any = c.as_mut();
            any.downcast_mut::<T>()
        })
    }

    /// Broadcast a message to all components.
    pub fn broadcast(&mut self, msg: &AppMessage) {
        for component in &mut self.components {
            component.on_message(msg);
        }
    }

    /// Get a component's display title by ID.
    pub fn get_title(&self, id: &str) -> Option<String> {
        self.components
            .iter()
            .find(|c| c.id() == id)
            .map(|c| c.title().to_owned())
    }

    /// Look up a component by ID and call `f` with it.
    ///
    /// Used by the layout engine to render a component into a panel.
    pub fn show_in_panel(
        &mut self,
        id: &str,
        f: impl FnOnce(&mut dyn EguiComponent, &SharedContext, &mut MessageBus),
        shared: &SharedContext,
        bus: &mut MessageBus,
    ) {
        if let Some(component) = self.components.iter_mut().find(|c| c.id() == id) {
            f(component.as_mut(), shared, bus);
        } else {
            log::warn!("Unknown component ID in layout: {id}");
        }
    }
}

impl Default for ComponentStore {
    fn default() -> Self {
        Self::new()
    }
}
