//! Keyboard binding system
//!
//! This module provides types for defining keyboard shortcuts with optional
//! modifier keys (Ctrl, Shift, Alt).
//!
//! The main types are:
//! - [`KeyBinding`] - Represents a key with optional modifiers (Ctrl, Shift, Alt)
//! - [`KeyBindings`] - Storage for key-to-action mappings

use std::collections::HashMap;

pub use winit::keyboard::KeyCode;

/// Represents a keyboard shortcut with optional modifier keys
///
/// A key binding consists of a main key code and optional modifier keys
/// (Ctrl, Shift, Alt). Use the builder-style methods to add modifiers.
///
/// # Examples
///
/// ```ignore
/// use pymol_scene::{KeyBinding, KeyCode};
///
/// // Simple key
/// let key = KeyBinding::new(KeyCode::KeyR);
///
/// // Key with modifiers
/// let ctrl_r = KeyBinding::new(KeyCode::KeyR).ctrl();
/// let ctrl_shift_s = KeyBinding::new(KeyCode::KeyS).ctrl().shift();
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct KeyBinding {
    /// The main key code
    pub key: KeyCode,
    /// Whether Ctrl (or Cmd on macOS) must be held
    pub ctrl: bool,
    /// Whether Shift must be held
    pub shift: bool,
    /// Whether Alt (or Option on macOS) must be held
    pub alt: bool,
}

impl KeyBinding {
    /// Create a new key binding for a single key without modifiers
    pub const fn new(key: KeyCode) -> Self {
        Self {
            key,
            ctrl: false,
            shift: false,
            alt: false,
        }
    }

    /// Require Ctrl (or Cmd on macOS) to be held
    pub const fn ctrl(mut self) -> Self {
        self.ctrl = true;
        self
    }

    /// Require Shift to be held
    pub const fn shift(mut self) -> Self {
        self.shift = true;
        self
    }

    /// Require Alt (or Option on macOS) to be held
    pub const fn alt(mut self) -> Self {
        self.alt = true;
        self
    }
}

/// Allow creating a KeyBinding from just a KeyCode for convenience
impl From<KeyCode> for KeyBinding {
    fn from(key: KeyCode) -> Self {
        Self::new(key)
    }
}

/// Storage for keyboard bindings
///
/// A generic container that maps [`KeyBinding`]s to actions of type `A`.
/// This is typically used with closure types for handling key events.
///
/// # Type Parameter
///
/// - `A`: The action type, typically `Arc<dyn Fn(&mut T) + Send + Sync>` for some target type `T`
///
/// # Example
///
/// ```ignore
/// use pymol_scene::{KeyBindings, KeyBinding, KeyCode};
/// use std::sync::Arc;
///
/// let mut bindings: KeyBindings<Arc<dyn Fn() + Send + Sync>> = KeyBindings::new();
/// bindings.bind(KeyCode::KeyR, Arc::new(|| println!("R pressed")));
/// bindings.bind(KeyBinding::new(KeyCode::KeyS).ctrl(), Arc::new(|| println!("Ctrl+S")));
/// ```
#[derive(Debug)]
pub struct KeyBindings<A> {
    bindings: HashMap<KeyBinding, A>,
}

impl<A> Default for KeyBindings<A> {
    fn default() -> Self {
        Self::new()
    }
}

impl<A> KeyBindings<A> {
    /// Create an empty key bindings container
    pub fn new() -> Self {
        Self {
            bindings: HashMap::new(),
        }
    }

    /// Bind a key or key combination to an action
    ///
    /// If the key was already bound, the previous binding is replaced.
    pub fn bind<K: Into<KeyBinding>>(&mut self, key: K, action: A) {
        self.bindings.insert(key.into(), action);
    }

    /// Remove a key binding
    ///
    /// Returns `true` if a binding was removed.
    pub fn unbind<K: Into<KeyBinding>>(&mut self, key: K) -> bool {
        self.bindings.remove(&key.into()).is_some()
    }

    /// Check if a key is bound to an action
    pub fn is_bound<K: Into<KeyBinding>>(&self, key: K) -> bool {
        self.bindings.contains_key(&key.into())
    }

    /// Get the action bound to a key binding
    pub fn get(&self, binding: &KeyBinding) -> Option<&A> {
        self.bindings.get(binding)
    }

    /// Get the number of bindings
    pub fn len(&self) -> usize {
        self.bindings.len()
    }

    /// Check if there are no bindings
    pub fn is_empty(&self) -> bool {
        self.bindings.is_empty()
    }

    /// Remove all bindings
    pub fn clear(&mut self) {
        self.bindings.clear();
    }
}

impl<A: Clone> KeyBindings<A> {
    /// Get a clone of the action bound to a key binding
    ///
    /// This is useful when the action needs to be executed while
    /// the bindings container is borrowed.
    pub fn get_cloned(&self, binding: &KeyBinding) -> Option<A> {
        self.bindings.get(binding).cloned()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_key_binding_new() {
        let binding = KeyBinding::new(KeyCode::KeyR);
        assert_eq!(binding.key, KeyCode::KeyR);
        assert!(!binding.ctrl);
        assert!(!binding.shift);
        assert!(!binding.alt);
    }

    #[test]
    fn test_key_binding_modifiers() {
        let binding = KeyBinding::new(KeyCode::KeyS).ctrl().shift();
        assert_eq!(binding.key, KeyCode::KeyS);
        assert!(binding.ctrl);
        assert!(binding.shift);
        assert!(!binding.alt);
    }

    #[test]
    fn test_key_binding_from_keycode() {
        let binding: KeyBinding = KeyCode::KeyA.into();
        assert_eq!(binding.key, KeyCode::KeyA);
        assert!(!binding.ctrl);
        assert!(!binding.shift);
        assert!(!binding.alt);
    }

    #[test]
    fn test_key_binding_equality() {
        let a = KeyBinding::new(KeyCode::KeyR).ctrl();
        let b = KeyBinding::new(KeyCode::KeyR).ctrl();
        let c = KeyBinding::new(KeyCode::KeyR);

        assert_eq!(a, b);
        assert_ne!(a, c);
    }

    #[test]
    fn test_key_bindings_new() {
        let bindings: KeyBindings<i32> = KeyBindings::new();
        assert!(bindings.is_empty());
        assert_eq!(bindings.len(), 0);
    }

    #[test]
    fn test_key_bindings_bind_and_get() {
        let mut bindings: KeyBindings<&str> = KeyBindings::new();
        bindings.bind(KeyCode::KeyR, "reset");
        bindings.bind(KeyBinding::new(KeyCode::KeyS).ctrl(), "save");

        assert_eq!(bindings.len(), 2);
        assert_eq!(bindings.get(&KeyBinding::new(KeyCode::KeyR)), Some(&"reset"));
        assert_eq!(
            bindings.get(&KeyBinding::new(KeyCode::KeyS).ctrl()),
            Some(&"save")
        );
        assert_eq!(bindings.get(&KeyBinding::new(KeyCode::KeyA)), None);
    }

    #[test]
    fn test_key_bindings_unbind() {
        let mut bindings: KeyBindings<i32> = KeyBindings::new();
        bindings.bind(KeyCode::KeyR, 1);
        bindings.bind(KeyCode::KeyS, 2);

        assert!(bindings.unbind(KeyCode::KeyR));
        assert!(!bindings.unbind(KeyCode::KeyR)); // Already removed
        assert_eq!(bindings.len(), 1);
    }

    #[test]
    fn test_key_bindings_is_bound() {
        let mut bindings: KeyBindings<i32> = KeyBindings::new();
        bindings.bind(KeyCode::KeyR, 1);

        assert!(bindings.is_bound(KeyCode::KeyR));
        assert!(!bindings.is_bound(KeyCode::KeyS));
    }

    #[test]
    fn test_key_bindings_clear() {
        let mut bindings: KeyBindings<i32> = KeyBindings::new();
        bindings.bind(KeyCode::KeyR, 1);
        bindings.bind(KeyCode::KeyS, 2);

        bindings.clear();
        assert!(bindings.is_empty());
    }

    #[test]
    fn test_key_bindings_modifier_distinction() {
        let mut bindings: KeyBindings<&str> = KeyBindings::new();
        bindings.bind(KeyCode::KeyR, "plain R");
        bindings.bind(KeyBinding::new(KeyCode::KeyR).ctrl(), "Ctrl+R");
        bindings.bind(KeyBinding::new(KeyCode::KeyR).ctrl().shift(), "Ctrl+Shift+R");

        assert_eq!(bindings.len(), 3);
        assert_eq!(bindings.get(&KeyBinding::new(KeyCode::KeyR)), Some(&"plain R"));
        assert_eq!(
            bindings.get(&KeyBinding::new(KeyCode::KeyR).ctrl()),
            Some(&"Ctrl+R")
        );
        assert_eq!(
            bindings.get(&KeyBinding::new(KeyCode::KeyR).ctrl().shift()),
            Some(&"Ctrl+Shift+R")
        );
    }
}
