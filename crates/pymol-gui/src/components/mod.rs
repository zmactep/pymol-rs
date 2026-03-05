//! Concrete Component Implementations
//!
//! Each module bundles a domain model, egui-specific UI state, and view logic
//! into a single struct that implements [`Component`](crate::component::Component).

pub mod movie;
pub mod object_list;
pub mod repl;
pub mod sequence;

pub use movie::MovieComponent;
pub use object_list::ObjectListComponent;
pub use repl::ReplComponent;
pub use sequence::SequenceComponent;

use pymol_framework::component::Component;

/// Create the default set of UI components.
///
/// Add new components here rather than in `App::new()`.
pub fn default_components() -> Vec<Box<dyn Component>> {
    vec![
        Box::new(ReplComponent::new()),
        Box::new(SequenceComponent::new()),
        Box::new(ObjectListComponent::new()),
        Box::new(MovieComponent::new()),
    ]
}
