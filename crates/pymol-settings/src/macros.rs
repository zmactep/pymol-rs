//! Declarative macros for generating typed settings groups.
//!
//! The `define_settings_group!` macro is the single source of truth for each
//! settings group. It generates:
//! 1. The settings struct with concrete field types
//! 2. `Default` impl with specified default values
//! 3. (For `group`) An `*Overrides` struct with `Option<T>` fields
//! 4. (For `group`) A `Merge` impl
//! 5. Descriptor entries for the string-based `set`/`get` command registry

use crate::error::SettingError;
use crate::setting::{SettingType, SettingValue};
use crate::side_effects::SideEffectCategory;

/// Trait that maps Rust types to `SettingType`/`SettingValue` conversions.
///
/// Implemented for `bool`, `i32`, `f32`, `[f32; 3]`, and any type implementing
/// `SettingEnum` (via blanket impl below).
pub trait SettingKind: Clone {
    fn setting_type() -> SettingType;
    fn wrap(&self) -> SettingValue;
    fn unwrap_value(v: SettingValue) -> Result<Self, SettingError>;
}

impl SettingKind for bool {
    fn setting_type() -> SettingType { SettingType::Bool }
    fn wrap(&self) -> SettingValue { SettingValue::Bool(*self) }
    fn unwrap_value(v: SettingValue) -> Result<Self, SettingError> {
        v.as_bool().ok_or(SettingError::type_mismatch("bool"))
    }
}

impl SettingKind for i32 {
    fn setting_type() -> SettingType { SettingType::Int }
    fn wrap(&self) -> SettingValue { SettingValue::Int(*self) }
    fn unwrap_value(v: SettingValue) -> Result<Self, SettingError> {
        v.as_int().ok_or(SettingError::type_mismatch("int"))
    }
}

impl SettingKind for f32 {
    fn setting_type() -> SettingType { SettingType::Float }
    fn wrap(&self) -> SettingValue { SettingValue::Float(*self) }
    fn unwrap_value(v: SettingValue) -> Result<Self, SettingError> {
        v.as_float().ok_or(SettingError::type_mismatch("float"))
    }
}

impl SettingKind for [f32; 3] {
    fn setting_type() -> SettingType { SettingType::Float3 }
    fn wrap(&self) -> SettingValue { SettingValue::Float3(*self) }
    fn unwrap_value(v: SettingValue) -> Result<Self, SettingError> {
        v.as_float3().ok_or(SettingError::type_mismatch("float3"))
    }
}

/// Blanket impl for any enum type implementing `SettingEnum` — stored as Int on the wire.
impl<T: crate::SettingEnum + Clone> SettingKind for T {
    fn setting_type() -> SettingType { SettingType::Int }
    fn wrap(&self) -> SettingValue { SettingValue::Int(self.to_i32()) }
    fn unwrap_value(v: SettingValue) -> Result<Self, SettingError> {
        v.as_int().map(T::from_i32).ok_or(SettingError::type_mismatch("int (enum)"))
    }
}

/// Metadata for a single setting field, used by the descriptor registry.
pub struct SettingDescriptor {
    /// String name for the `set`/`get` commands.
    pub name: &'static str,
    /// Value type.
    pub setting_type: SettingType,
    /// Read the current value from global `Settings`.
    pub get: fn(&crate::groups::Settings) -> SettingValue,
    /// Write a value to global `Settings`.
    pub set: fn(&mut crate::groups::Settings, SettingValue) -> Result<(), SettingError>,
    /// Read override from `ObjectOverrides`. `None` getter → global-only setting.
    pub get_override: Option<fn(&crate::overrides::ObjectOverrides) -> Option<SettingValue>>,
    /// Write override to `ObjectOverrides`.
    pub set_override:
        Option<fn(&mut crate::overrides::ObjectOverrides, SettingValue) -> Result<(), SettingError>>,
    /// Clear override (revert to global).
    pub unset_override: Option<fn(&mut crate::overrides::ObjectOverrides)>,
    /// Min constraint (numeric).
    pub min: Option<f32>,
    /// Max constraint (numeric).
    pub max: Option<f32>,
    /// Named value variants (for enum settings).
    pub value_hints: &'static [(&'static str, SettingValue)],
    /// Side effects triggered when this setting changes.
    pub side_effects: &'static [SideEffectCategory],
}

impl SettingDescriptor {
    /// Whether this setting has named value variants.
    pub fn has_value_hints(&self) -> bool {
        !self.value_hints.is_empty()
    }

    /// Look up a setting value from a string alias (case-insensitive).
    pub fn resolve_hint(&self, name: &str) -> Option<&'static SettingValue> {
        let name_lower = name.to_lowercase();
        self.value_hints
            .iter()
            .find(|(n, _)| n.to_lowercase() == name_lower)
            .map(|(_, v)| v)
    }

    /// Look up a display name from a setting value.
    pub fn hint_name(&self, value: &SettingValue) -> Option<&'static str> {
        self.value_hints
            .iter()
            .find(|(_, v)| v == value)
            .map(|(n, _)| *n)
    }

    /// Get all hint names (for autocomplete).
    pub fn hint_names(&self) -> impl Iterator<Item = &'static str> + '_ {
        self.value_hints.iter().map(|(n, _)| *n)
    }

    /// Whether this setting supports per-object overrides.
    pub fn is_object_overridable(&self) -> bool {
        self.set_override.is_some()
    }
}

/// Trait for merging a base settings struct with an overrides struct.
pub trait Merge: Clone {
    type Overrides: Default;
    /// Return a copy of `self` with non-None overrides applied.
    fn with_overrides(&self, overrides: &Self::Overrides) -> Self;
}

// =============================================================================
// Core macro: define_settings_group!
// =============================================================================

/// Define a typed settings group.
///
/// # Forms
///
/// **`group` — object-overridable** (generates struct + Overrides + Merge):
/// ```ignore
/// define_settings_group! {
///     group StickSettings / StickOverrides {
///         radius: f32 = 0.25, name = "stick_radius",
///             min = 0.01, max = 5.0,
///             side_effects = [RepresentationRebuild];
///         color: i32 = -1, name = "stick_color",
///             side_effects = [ColorRebuild];
///     }
/// }
/// ```
///
/// **`group_global` — global-only** (generates struct only, no Overrides):
/// ```ignore
/// define_settings_group! {
///     group_global RayTraceSettings {
///         max_passes: i32 = 25, name = "ray_max_passes";
///     }
/// }
/// ```
#[macro_export]
macro_rules! define_settings_group {
    // =========================================================================
    // group (object-overridable): generates struct + Overrides + Merge
    // =========================================================================
    (
        $(#[$meta:meta])*
        group $Name:ident / $Overrides:ident {
            $(
                $(#[$field_meta:meta])*
                $field:ident : $ty:ty = $default:expr,
                    name = $sname:expr
                    $(, min = $min:expr, max = $max:expr)?
                    $(, hints = $hints_ty:ty)?
                    $(, side_effects = [$($se:ident),* $(,)?])?
                    ;
            )*
        }
    ) => {
        $(#[$meta])*
        #[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
        pub struct $Name {
            $(
                $(#[$field_meta])*
                pub $field: $ty,
            )*
        }

        impl Default for $Name {
            fn default() -> Self {
                Self {
                    $( $field: $default, )*
                }
            }
        }

        #[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
        pub struct $Overrides {
            $(
                pub $field: Option<$ty>,
            )*
        }

        impl $crate::macros::Merge for $Name {
            type Overrides = $Overrides;

            fn with_overrides(&self, o: &$Overrides) -> Self {
                Self {
                    $( $field: o.$field.unwrap_or(self.$field.clone()), )*
                }
            }
        }

        impl $Name {
            /// Return descriptor entries for all fields in this group.
            #[allow(unused_imports)]
            pub fn descriptors() -> Vec<$crate::macros::SettingDescriptor> {
                use $crate::SettingEnum as _;
                vec![
                    $(
                        $crate::macros::SettingDescriptor {
                            name: $sname,
                            setting_type: <$ty as $crate::macros::SettingKind>::setting_type(),
                            get: |_s| { $crate::SettingValue::Int(0) }, // placeholder
                            set: |_s, _v| { Ok(()) }, // placeholder
                            get_override: None, // placeholder
                            set_override: None, // placeholder
                            unset_override: None, // placeholder
                            min: $crate::define_settings_group!(@opt_f32 $($min)?),
                            max: $crate::define_settings_group!(@opt_f32 $($max)?),
                            value_hints: $crate::define_settings_group!(@hints $($hints_ty)?),
                            side_effects: &[$($($crate::SideEffectCategory::$se),*)?],
                        },
                    )*
                ]
            }
        }
    };

    // =========================================================================
    // group_global: generates struct only, no Overrides
    // =========================================================================
    (
        $(#[$meta:meta])*
        group_global $Name:ident {
            $(
                $(#[$field_meta:meta])*
                $field:ident : $ty:ty = $default:expr,
                    name = $sname:expr
                    $(, min = $min:expr, max = $max:expr)?
                    $(, hints = $hints_ty:ty)?
                    $(, side_effects = [$($se:ident),* $(,)?])?
                    ;
            )*
        }
    ) => {
        $(#[$meta])*
        #[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
        pub struct $Name {
            $(
                $(#[$field_meta])*
                pub $field: $ty,
            )*
        }

        impl Default for $Name {
            fn default() -> Self {
                Self {
                    $( $field: $default, )*
                }
            }
        }

        impl $Name {
            /// Return descriptor entries for all fields in this group.
            #[allow(unused_imports)]
            pub fn descriptors() -> Vec<$crate::macros::SettingDescriptor> {
                use $crate::SettingEnum as _;
                vec![
                    $(
                        $crate::macros::SettingDescriptor {
                            name: $sname,
                            setting_type: <$ty as $crate::macros::SettingKind>::setting_type(),
                            get: |_s| { $crate::SettingValue::Int(0) }, // placeholder
                            set: |_s, _v| { Ok(()) }, // placeholder
                            get_override: None,
                            set_override: None,
                            unset_override: None,
                            min: $crate::define_settings_group!(@opt_f32 $($min)?),
                            max: $crate::define_settings_group!(@opt_f32 $($max)?),
                            value_hints: $crate::define_settings_group!(@hints $($hints_ty)?),
                            side_effects: &[$($($crate::SideEffectCategory::$se),*)?],
                        },
                    )*
                ]
            }
        }
    };

    // =========================================================================
    // Internal helpers
    // =========================================================================
    (@opt_f32) => { None };
    (@opt_f32 $val:expr) => { Some($val as f32) };

    (@hints) => { &[] };
    (@hints $ty:ty) => { <$ty as $crate::SettingEnum>::value_hints() };
}

#[cfg(test)]
mod tests {
    use super::*;

    // Test the `group` form (object-overridable): generates struct + Overrides + Merge
    crate::define_settings_group! {
        /// Test stick settings
        group TestStick / TestStickOverrides {
            radius: f32 = 0.25,
                name = "stick_radius",
                min = 0.01, max = 5.0,
                side_effects = [RepresentationRebuild];
            color: i32 = -1,
                name = "stick_color",
                side_effects = [ColorRebuild];
        }
    }

    // Test the `group_global` form (global-only): generates struct only
    crate::define_settings_group! {
        group_global TestGlobal {
            max_passes: i32 = 25,
                name = "ray_max_passes";
            enabled: bool = true,
                name = "test_enabled",
                side_effects = [SceneInvalidate];
        }
    }

    // Test with enum type + hints
    crate::define_settings_group! {
        group_global TestWithEnum {
            mode: crate::ShadingMode = crate::ShadingMode::Classic,
                name = "shading_mode",
                hints = crate::ShadingMode,
                side_effects = [ShaderReload];
        }
    }

    #[test]
    fn test_group_defaults() {
        let s = TestStick::default();
        assert_eq!(s.radius, 0.25);
        assert_eq!(s.color, -1);
    }

    #[test]
    fn test_overrides_default_is_none() {
        let o = TestStickOverrides::default();
        assert!(o.radius.is_none());
        assert!(o.color.is_none());
    }

    #[test]
    fn test_merge_no_overrides() {
        let base = TestStick::default();
        let o = TestStickOverrides::default();
        let merged = base.with_overrides(&o);
        assert_eq!(merged.radius, 0.25);
        assert_eq!(merged.color, -1);
    }

    #[test]
    fn test_merge_with_override() {
        let base = TestStick::default();
        let o = TestStickOverrides {
            radius: Some(0.5),
            color: None,
        };
        let merged = base.with_overrides(&o);
        assert_eq!(merged.radius, 0.5);
        assert_eq!(merged.color, -1);
    }

    #[test]
    fn test_global_defaults() {
        let s = TestGlobal::default();
        assert_eq!(s.max_passes, 25);
        assert!(s.enabled);
    }

    #[test]
    fn test_descriptors_group() {
        let descs = TestStick::descriptors();
        assert_eq!(descs.len(), 2);
        assert_eq!(descs[0].name, "stick_radius");
        assert_eq!(descs[0].setting_type, SettingType::Float);
        assert_eq!(descs[0].min, Some(0.01));
        assert_eq!(descs[0].max, Some(5.0));
        assert_eq!(descs[0].side_effects.len(), 1);

        assert_eq!(descs[1].name, "stick_color");
        assert_eq!(descs[1].setting_type, SettingType::Int);
        assert!(descs[1].min.is_none());
    }

    #[test]
    fn test_descriptors_global() {
        let descs = TestGlobal::descriptors();
        assert_eq!(descs.len(), 2);
        assert_eq!(descs[0].name, "ray_max_passes");
        assert!(descs[0].get_override.is_none());
    }

    #[test]
    fn test_enum_hints() {
        let descs = TestWithEnum::descriptors();
        assert_eq!(descs.len(), 1);
        assert_eq!(descs[0].name, "shading_mode");
        assert_eq!(descs[0].setting_type, SettingType::Int); // enum stored as Int
        assert_eq!(descs[0].value_hints.len(), 3); // classic, skripkin, full
        assert_eq!(descs[0].value_hints[0].0, "classic");
    }
}
