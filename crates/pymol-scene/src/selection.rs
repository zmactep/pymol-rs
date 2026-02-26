//! Named selection management
//!
//! This module provides storage and evaluation of named selections,
//! similar to PyMOL's selection system.

use ahash::AHashMap;
use pymol_select::{EvalContext, SelectionOptions, SelectionResult};
use serde::{Deserialize, Serialize};

use crate::object::ObjectRegistry;

/// Entry for a named selection
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SelectionEntry {
    /// The selection expression
    pub expression: String,
    /// Whether the selection indicators are visible
    pub visible: bool,
}

impl SelectionEntry {
    /// Create a new selection entry with visibility enabled by default
    pub fn new(expression: String) -> Self {
        Self {
            expression,
            visible: true,
        }
    }

    /// Create a new selection entry with specified visibility
    pub fn with_visibility(expression: String, visible: bool) -> Self {
        Self { expression, visible }
    }
}

/// Manager for named selections
///
/// Handles storage, retrieval, and evaluation of named selection expressions.
#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub struct SelectionManager {
    /// Named selections (name -> entry with expression and visibility)
    selections: AHashMap<String, SelectionEntry>,
}

impl SelectionManager {
    /// Create a new empty selection manager
    pub fn new() -> Self {
        Self::default()
    }

    /// Get a selection entry by name
    pub fn get(&self, name: &str) -> Option<&SelectionEntry> {
        self.selections.get(name)
    }

    /// Get a mutable selection entry by name
    pub fn get_mut(&mut self, name: &str) -> Option<&mut SelectionEntry> {
        self.selections.get_mut(name)
    }

    /// Get a selection expression by name
    pub fn get_expression(&self, name: &str) -> Option<&str> {
        self.selections.get(name).map(|e| e.expression.as_str())
    }

    /// Define (store) a named selection expression
    ///
    /// If a selection with this name already exists, it will be replaced.
    pub fn define(&mut self, name: &str, expression: &str) {
        self.selections.insert(
            name.to_string(),
            SelectionEntry::new(expression.to_string()),
        );
    }

    /// Remove a named selection
    ///
    /// Returns true if the selection existed and was removed.
    pub fn remove(&mut self, name: &str) -> bool {
        self.selections.remove(name).is_some()
    }

    /// Get all selection names
    pub fn names(&self) -> Vec<String> {
        self.selections.keys().cloned().collect()
    }

    /// Rename a named selection
    ///
    /// Returns true if the selection existed and was renamed.
    pub fn rename(&mut self, old_name: &str, new_name: &str) -> bool {
        if let Some(entry) = self.selections.remove(old_name) {
            self.selections.insert(new_name.to_string(), entry);
            true
        } else {
            false
        }
    }

    /// Check if a selection exists
    pub fn contains(&self, name: &str) -> bool {
        self.selections.contains_key(name)
    }

    /// Get the number of selections
    pub fn len(&self) -> usize {
        self.selections.len()
    }

    /// Check if there are no selections
    pub fn is_empty(&self) -> bool {
        self.selections.is_empty()
    }

    /// Set the visibility of a named selection's indicators
    pub fn set_visible(&mut self, name: &str, visible: bool) {
        if let Some(entry) = self.selections.get_mut(name) {
            entry.visible = visible;
        }
    }

    /// Check if a named selection's indicators are visible
    pub fn is_visible(&self, name: &str) -> bool {
        self.selections
            .get(name)
            .map(|e| e.visible)
            .unwrap_or(false)
    }

    /// Set the "indicate" selection (shows pink indicators in the 3D view)
    ///
    /// This creates or updates the special "indicate" selection.
    pub fn indicate(&mut self, selection: &str) {
        self.selections.insert(
            "indicate".to_string(),
            SelectionEntry::new(selection.to_string()),
        );
    }

    /// Clear all selection indicators (hide all visible selections)
    pub fn clear_indication(&mut self) {
        for entry in self.selections.values_mut() {
            entry.visible = false;
        }
    }

    /// Get the currently indicated selection expression
    ///
    /// Returns the first visible selection expression for backwards compatibility.
    pub fn indicated_selection(&self) -> Option<&str> {
        self.selections
            .values()
            .find(|e| e.visible)
            .map(|e| e.expression.as_str())
    }

    /// Get an iterator over all selections
    pub fn iter(&self) -> impl Iterator<Item = (&String, &SelectionEntry)> {
        self.selections.iter()
    }

    /// Get an iterator over visible selections
    pub fn visible_selections(&self) -> impl Iterator<Item = (&str, &str)> {
        self.selections
            .iter()
            .filter(|(_, entry)| entry.visible)
            .map(|(name, entry)| (name.as_str(), entry.expression.as_str()))
    }

    /// Clear all selections
    pub fn clear(&mut self) {
        self.selections.clear();
    }

    /// Evaluate all visible selections for all molecules.
    ///
    /// Builds a full evaluation context per molecule (implicit object-name
    /// selections + pre-evaluated named selections) and combines all visible
    /// selection results with OR.
    ///
    /// Returns a vector of (object_name, SelectionResult) pairs where the
    /// SelectionResult is the union of all visible selections.
    pub fn evaluate_visible(
        &self,
        registry: &ObjectRegistry,
        options: SelectionOptions,
    ) -> Vec<(String, SelectionResult)> {
        let visible: Vec<(&str, &str)> = self.visible_selections().collect();
        if visible.is_empty() {
            return Vec::new();
        }

        let object_names: Vec<String> = registry.names().map(|s| s.to_string()).collect();
        let mut results: Vec<(String, SelectionResult)> = Vec::new();

        for obj_name in &object_names {
            let Some(mol_obj) = registry.get_molecule(obj_name) else {
                continue;
            };
            let ctx = self.build_eval_context(mol_obj.molecule(), obj_name, &object_names, options);

            // Evaluate and OR together all visible selections
            let mut combined: Option<SelectionResult> = None;
            for (sel_name, sel_expr) in &visible {
                let ast = match pymol_select::parse(sel_expr) {
                    Ok(e) => e,
                    Err(e) => {
                        log::debug!("Failed to parse selection '{}': {:?}", sel_name, e);
                        continue;
                    }
                };
                match pymol_select::evaluate(&ast, &ctx) {
                    Ok(result) => {
                        combined = Some(match combined {
                            Some(existing) => existing.union(&result),
                            None => result,
                        });
                    }
                    Err(e) => {
                        log::debug!(
                            "Failed to evaluate selection '{}' for '{}': {:?}",
                            sel_name, obj_name, e
                        );
                    }
                }
            }

            if let Some(result) = combined {
                if result.any() {
                    results.push((obj_name.clone(), result));
                }
            }
        }

        results
    }

    /// Build an `EvalContext` for a single molecule within the registry.
    ///
    /// The context includes:
    /// - Implicit object-name selections (all for the target, none for others)
    /// - Pre-evaluated named selections for reference resolution
    /// - The provided `SelectionOptions` for case-sensitivity control
    ///
    /// This is the shared context-building logic used by both command-level
    /// selection evaluation and rendering-level visible selection evaluation.
    pub fn build_eval_context<'a>(
        &self,
        mol: &'a pymol_mol::ObjectMolecule,
        obj_name: &str,
        all_object_names: &[String],
        options: SelectionOptions,
    ) -> EvalContext<'a> {
        let mut ctx = EvalContext::single(mol);
        ctx.options = options;

        // Implicit object-name selections
        let atom_count = mol.atom_count();
        for other in all_object_names {
            if other == obj_name {
                ctx.add_selection(other.clone(), SelectionResult::all(atom_count));
            } else {
                ctx.add_selection(other.clone(), SelectionResult::none(atom_count));
            }
        }

        // Pre-evaluate all named selections into the context
        for (sel_name, entry) in &self.selections {
            if let Ok(ast) = pymol_select::parse(&entry.expression) {
                if let Ok(result) = pymol_select::evaluate(&ast, &ctx) {
                    ctx.add_selection(sel_name.clone(), result);
                }
            }
        }

        ctx
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_selection_entry() {
        let entry = SelectionEntry::new("chain A".to_string());
        assert_eq!(entry.expression, "chain A");
        assert!(entry.visible);
    }

    #[test]
    fn test_selection_manager_basic() {
        let mut manager = SelectionManager::new();

        manager.define("sel1", "chain A");
        assert!(manager.contains("sel1"));
        assert_eq!(manager.get_expression("sel1"), Some("chain A"));

        manager.define("sel2", "resn ALA");
        assert_eq!(manager.len(), 2);

        assert!(manager.remove("sel1"));
        assert!(!manager.contains("sel1"));
        assert_eq!(manager.len(), 1);
    }

    #[test]
    fn test_selection_visibility() {
        let mut manager = SelectionManager::new();

        manager.define("sel1", "chain A");
        assert!(manager.is_visible("sel1"));

        manager.set_visible("sel1", false);
        assert!(!manager.is_visible("sel1"));

        manager.set_visible("sel1", true);
        assert!(manager.is_visible("sel1"));
    }

    #[test]
    fn test_indicate_selection() {
        let mut manager = SelectionManager::new();

        manager.indicate("chain B");
        assert!(manager.contains("indicate"));
        assert_eq!(manager.get_expression("indicate"), Some("chain B"));
        assert_eq!(manager.indicated_selection(), Some("chain B"));
    }

    #[test]
    fn test_clear_indication() {
        let mut manager = SelectionManager::new();

        manager.define("sel1", "chain A");
        manager.define("sel2", "chain B");

        manager.clear_indication();

        assert!(!manager.is_visible("sel1"));
        assert!(!manager.is_visible("sel2"));
        assert!(manager.indicated_selection().is_none());
    }
}
