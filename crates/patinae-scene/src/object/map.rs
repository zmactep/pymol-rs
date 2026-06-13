//! Map object for electron density maps
//!
//! Provides storage and rendering of volumetric density data
//! with isomesh and isosurface visualization.

use std::sync::Arc;

use lin_alg::f32::Vec3;
use patinae_algos::surface::Grid3D;
use patinae_mol::Symmetry;
use patinae_settings::ObjectOverrides;

use serde::{Deserialize, Serialize};

use super::{Object, ObjectState, ObjectType};

/// Data for a single map state
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MapData {
    /// 3D grid storing density values
    pub grid: Arc<Grid3D>,
    /// Symmetry information (if crystallographic)
    pub symmetry: Option<Symmetry>,
}

impl MapData {
    /// Create a new map data from a grid
    pub fn new(grid: Grid3D) -> Self {
        Self {
            grid: Arc::new(grid),
            symmetry: None,
        }
    }

    /// Create map data with symmetry
    pub fn with_symmetry(grid: Grid3D, symmetry: Symmetry) -> Self {
        Self {
            grid: Arc::new(grid),
            symmetry: Some(symmetry),
        }
    }

    /// Create map data from shared grid storage.
    pub fn from_shared_grid(grid: Arc<Grid3D>) -> Self {
        Self {
            grid,
            symmetry: None,
        }
    }

    /// Get the bounding box of this map
    pub fn bounds(&self) -> ([f32; 3], [f32; 3]) {
        let vd = self.grid.vertex_dims();
        let min = self.grid.origin;
        let max = [
            min[0] + (vd[0] - 1) as f32 * self.grid.spacing[0],
            min[1] + (vd[1] - 1) as f32 * self.grid.spacing[1],
            min[2] + (vd[2] - 1) as f32 * self.grid.spacing[2],
        ];
        (min, max)
    }
}

/// Display mode for map objects
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum MapDisplayMode {
    /// Loaded map data without an active contour visual.
    #[default]
    None,
    /// Wireframe mesh at isovalue
    Isomesh,
    /// Solid surface at isovalue
    Isosurface,
    /// Dot cloud at isovalue crossings
    Isodot,
    /// Volume rendering (not yet implemented)
    Volume,
}

/// An electron density map object
///
/// Stores 3D volumetric data and renders it as isomesh or isosurface.
/// Supports multiple states (e.g., for time-series maps).
pub struct MapObject {
    /// Map name
    name: String,
    /// Visual state (enabled, color, etc.)
    state: ObjectState,
    /// Map data for each state
    states: Vec<MapData>,
    /// Current display state index
    current_state: usize,
    /// Isocontour level (sigma or absolute)
    level: f32,
    /// Display mode (mesh, surface, volume)
    display_mode: MapDisplayMode,
    /// Mesh color (RGBA)
    mesh_color: [f32; 4],
    /// Carve around selection radius (0 = no carve)
    carve_radius: f32,
    /// Positions to carve around (atoms from selection)
    carve_positions: Option<Vec<[f32; 3]>>,
    /// Per-object settings overrides
    overrides: Option<ObjectOverrides>,
    /// Whether the cached data needs rebuilding (consumed by patinae-render
    /// when its map mesh path is implemented).
    dirty: bool,
    /// Revision bumped when contour geometry inputs change.
    geometry_revision: u64,
    /// Revision bumped when material-only inputs change.
    material_revision: u64,
}

impl MapObject {
    /// Create a new map object from grid data
    pub fn new(name: &str, grid: Grid3D) -> Self {
        Self {
            name: name.to_string(),
            state: ObjectState::default(),
            states: vec![MapData::new(grid)],
            current_state: 0,
            level: 1.0,
            display_mode: MapDisplayMode::default(),
            mesh_color: [0.0, 0.5, 1.0, 1.0], // Blue
            carve_radius: 0.0,
            carve_positions: None,
            overrides: None,
            dirty: true,
            geometry_revision: 1,
            material_revision: 1,
        }
    }

    /// Create a map object with map data
    pub fn from_map_data(name: &str, data: MapData) -> Self {
        Self {
            name: name.to_string(),
            state: ObjectState::default(),
            states: vec![data],
            current_state: 0,
            level: 1.0,
            display_mode: MapDisplayMode::default(),
            mesh_color: [0.0, 0.5, 1.0, 1.0],
            carve_radius: 0.0,
            carve_positions: None,
            overrides: None,
            dirty: true,
            geometry_revision: 1,
            material_revision: 1,
        }
    }

    /// Create an empty map object
    pub fn empty(name: &str) -> Self {
        Self {
            name: name.to_string(),
            state: ObjectState::default(),
            states: Vec::new(),
            current_state: 0,
            level: 1.0,
            display_mode: MapDisplayMode::default(),
            mesh_color: [0.0, 0.5, 1.0, 1.0],
            carve_radius: 0.0,
            carve_positions: None,
            overrides: None,
            dirty: true,
            geometry_revision: 1,
            material_revision: 1,
        }
    }

    /// Add a state to the map
    pub fn add_state(&mut self, data: MapData) {
        self.states.push(data);
        self.invalidate_geometry();
    }

    /// Get the current map data
    pub fn current_data(&self) -> Option<&MapData> {
        self.states.get(self.current_state)
    }

    /// Get mutable access to current map data
    pub fn current_data_mut(&mut self) -> Option<&mut MapData> {
        self.invalidate_geometry();
        self.states.get_mut(self.current_state)
    }

    /// Get the grid for the current state
    pub fn grid(&self) -> Option<&Grid3D> {
        self.states
            .get(self.current_state)
            .map(|data| data.grid.as_ref())
    }

    /// Get shared access to the current state's grid storage.
    pub fn grid_shared(&self) -> Option<Arc<Grid3D>> {
        self.states
            .get(self.current_state)
            .map(|data| Arc::clone(&data.grid))
    }

    /// Get mutable access to the grid (marks dirty)
    pub fn grid_mut(&mut self) -> Option<&mut Grid3D> {
        self.invalidate_geometry();
        self.states
            .get_mut(self.current_state)
            .map(|data| Arc::make_mut(&mut data.grid))
    }

    /// Get the isocontour level
    pub fn level(&self) -> f32 {
        self.level
    }

    /// Set the isocontour level
    pub fn set_level(&mut self, level: f32) {
        if (self.level - level).abs() > 1e-6 {
            self.level = level;
            self.invalidate_geometry();
        }
    }

    /// Get the display mode
    pub fn display_mode(&self) -> MapDisplayMode {
        self.display_mode
    }

    /// Set the display mode
    pub fn set_display_mode(&mut self, mode: MapDisplayMode) {
        if self.display_mode != mode {
            self.display_mode = mode;
            self.invalidate_geometry();
        }
    }

    /// Get the mesh color
    pub fn mesh_color(&self) -> [f32; 4] {
        self.mesh_color
    }

    /// Set the mesh color
    pub fn set_mesh_color(&mut self, color: [f32; 4]) {
        if self.mesh_color != color {
            self.mesh_color = color;
            self.invalidate_material();
        }
    }

    /// Get the carve radius
    pub fn carve_radius(&self) -> f32 {
        self.carve_radius
    }

    /// Set the carve radius
    pub fn set_carve_radius(&mut self, radius: f32) {
        if (self.carve_radius - radius).abs() > 1e-6 {
            self.carve_radius = radius;
            self.invalidate_geometry();
        }
    }

    /// Set positions to carve around
    pub fn set_carve_positions(&mut self, positions: Vec<[f32; 3]>) {
        self.carve_positions = Some(positions);
        self.invalidate_geometry();
    }

    /// Clear carve positions
    pub fn clear_carve_positions(&mut self) {
        self.carve_positions = None;
        self.invalidate_geometry();
    }

    /// Get all map states
    pub fn states(&self) -> &[MapData] {
        &self.states
    }

    /// Get the current state index
    pub fn current_state(&self) -> usize {
        self.current_state
    }

    /// Set the current state index
    pub fn set_current_state(&mut self, state: usize) {
        if state < self.states.len() && self.current_state != state {
            self.current_state = state;
            self.invalidate_geometry();
        }
    }

    /// Get carve positions
    pub fn carve_positions(&self) -> Option<&Vec<[f32; 3]>> {
        self.carve_positions.as_ref()
    }

    /// Get per-object overrides
    pub fn overrides(&self) -> Option<&ObjectOverrides> {
        self.overrides.as_ref()
    }

    /// Set per-object overrides
    pub fn set_overrides(&mut self, overrides: ObjectOverrides) {
        self.overrides = Some(overrides);
        self.invalidate_material();
    }

    /// Check if dirty
    pub fn is_dirty(&self) -> bool {
        self.dirty
    }

    /// Mark as dirty
    pub fn invalidate(&mut self) {
        self.invalidate_geometry();
        self.invalidate_material();
    }

    /// Mark contour geometry as dirty.
    pub fn invalidate_geometry(&mut self) {
        self.geometry_revision = self.geometry_revision.wrapping_add(1).max(1);
        self.dirty = true;
    }

    /// Mark material-only map data as dirty.
    pub fn invalidate_material(&mut self) {
        self.material_revision = self.material_revision.wrapping_add(1).max(1);
        self.dirty = true;
    }

    /// Clear short-lived dirty state after render sync.
    pub fn clear_dirty(&mut self) {
        self.dirty = false;
    }

    /// Get the current geometry revision.
    pub fn geometry_revision(&self) -> u64 {
        self.geometry_revision
    }

    /// Get the current material revision.
    pub fn material_revision(&self) -> u64 {
        self.material_revision
    }

    /// Whether this map mode can produce renderable geometry.
    pub fn is_renderable(&self) -> bool {
        matches!(
            self.display_mode,
            MapDisplayMode::Isomesh | MapDisplayMode::Isosurface
        )
    }

    /// Check if the map has any data
    pub fn is_empty(&self) -> bool {
        self.states.is_empty()
    }

    /// Compute the bounding box
    fn compute_extent(&self) -> Option<(Vec3, Vec3)> {
        let data = self.states.get(self.current_state)?;
        let (min, max) = data.bounds();
        Some((
            Vec3::new(min[0], min[1], min[2]),
            Vec3::new(max[0], max[1], max[2]),
        ))
    }
}

impl Object for MapObject {
    fn name(&self) -> &str {
        &self.name
    }

    fn object_type(&self) -> ObjectType {
        ObjectType::Map
    }

    fn state(&self) -> &ObjectState {
        &self.state
    }

    fn state_mut(&mut self) -> &mut ObjectState {
        &mut self.state
    }

    fn extent(&self) -> Option<(Vec3, Vec3)> {
        self.compute_extent()
    }

    fn n_states(&self) -> usize {
        self.states.len().max(1)
    }

    fn current_state(&self) -> usize {
        self.current_state
    }

    fn set_current_state(&mut self, state: usize) -> bool {
        if state < self.states.len() {
            if self.current_state != state {
                self.current_state = state;
                self.invalidate_geometry();
            }
            true
        } else {
            false
        }
    }

    fn overrides(&self) -> Option<&ObjectOverrides> {
        self.overrides.as_ref()
    }

    fn overrides_mut(&mut self) -> Option<&mut ObjectOverrides> {
        self.overrides.as_mut()
    }

    fn get_or_create_overrides(&mut self) -> &mut ObjectOverrides {
        if self.overrides.is_none() {
            self.overrides = Some(ObjectOverrides::default());
        }
        self.overrides.as_mut().unwrap()
    }

    fn set_name(&mut self, name: String) {
        self.name = name;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_map_creation() {
        let grid = Grid3D::from_bounds([0.0, 0.0, 0.0], [10.0, 10.0, 10.0], 1.0, 0.0);
        let map = MapObject::new("test_map", grid);

        assert_eq!(map.name(), "test_map");
        assert_eq!(map.object_type(), ObjectType::Map);
        assert!(map.is_enabled());
        assert!(!map.is_empty());
    }

    #[test]
    fn test_map_level() {
        let grid = Grid3D::from_bounds([0.0, 0.0, 0.0], [10.0, 10.0, 10.0], 1.0, 0.0);
        let mut map = MapObject::new("test", grid);

        assert_eq!(map.level(), 1.0);
        map.set_level(2.5);
        assert_eq!(map.level(), 2.5);
        assert!(map.is_dirty());
    }

    #[test]
    fn test_map_display_mode() {
        let grid = Grid3D::from_bounds([0.0, 0.0, 0.0], [10.0, 10.0, 10.0], 1.0, 0.0);
        let mut map = MapObject::new("test", grid);

        assert_eq!(map.display_mode(), MapDisplayMode::None);
        map.set_display_mode(MapDisplayMode::Isosurface);
        assert_eq!(map.display_mode(), MapDisplayMode::Isosurface);
        assert!(map.is_renderable());
    }

    #[test]
    fn test_map_extent() {
        let grid = Grid3D::from_bounds([0.0, 0.0, 0.0], [10.0, 10.0, 10.0], 1.0, 0.0);
        let map = MapObject::new("test", grid);

        let extent = map.extent().expect("Should have extent");
        assert!(extent.0.x <= 0.0);
        assert!(extent.1.x >= 10.0);
    }

    #[test]
    fn test_empty_map() {
        let map = MapObject::empty("test");
        assert!(map.is_empty());
        assert!(map.extent().is_none());
    }

    #[test]
    fn test_map_dirty_revisions_split_geometry_and_material() {
        let grid = Grid3D::from_bounds([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 1.0, 0.0);
        let mut map = MapObject::new("test", grid);
        let geometry_revision = map.geometry_revision();
        let material_revision = map.material_revision();

        map.set_level(2.0);
        assert!(map.geometry_revision() > geometry_revision);
        assert_eq!(map.material_revision(), material_revision);

        let geometry_revision = map.geometry_revision();
        map.set_mesh_color([1.0, 0.0, 0.0, 0.5]);
        assert_eq!(map.geometry_revision(), geometry_revision);
        assert!(map.material_revision() > material_revision);

        map.clear_dirty();
        assert!(!map.is_dirty());
    }
}
