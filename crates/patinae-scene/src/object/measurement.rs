//! Measurement object for distance, angle, and dihedral visualization
//!
//! Measurement objects display dashed lines between atoms with labels showing
//! the measured value (distance in Å, angles in degrees).

use lin_alg::f32::Vec3;
use serde::{Deserialize, Serialize};

use super::{Object, ObjectState, ObjectType};

/// Arc radius as fraction of the shorter leg length
const ANGLE_SIZE: f32 = 0.4;
/// Scale factor to push label outside the arc
const LABEL_OFFSET: f32 = 1.3;

// ── Arc geometry ─────────────────────────────────────────────────────────────

/// Precomputed arc frame for angle/dihedral measurements.
///
/// Stores the center point and two radius-scaled perpendicular axes that
/// define the arc plane, plus the sweep angle. Computed once per measurement
/// and reused by vertex generators, solid-line generators, and label placement.
#[derive(Debug, Clone)]
struct ArcGeometry {
    /// Arc center point (vertex for angles, bond midpoint for dihedrals)
    center: Vec3,
    /// Radius-scaled axis toward the first reference direction
    x_axis: Vec3,
    /// Radius-scaled axis perpendicular to x_axis in the arc plane
    y_axis: Vec3,
    /// Sweep angle in radians (signed for dihedrals)
    arc_angle: f32,
}

impl ArcGeometry {
    /// Compute the arc frame for an angle measurement (vertex at p2).
    fn for_angle(p1: Vec3, p2: Vec3, p3: Vec3, value_degrees: f64) -> Option<Self> {
        let d1 = p1 - p2;
        let d2 = p3 - p2;
        let len1 = d1.magnitude();
        let len2 = d2.magnitude();
        if len1 < 1e-6 || len2 < 1e-6 {
            return None;
        }

        let radius = len1.min(len2) * ANGLE_SIZE;
        let (x_axis, y_axis) = perpendicular_frame(d1 / len1, d2 / len2, radius)?;

        Some(Self {
            center: p2,
            x_axis,
            y_axis,
            arc_angle: (value_degrees as f32).to_radians(),
        })
    }

    /// Compute the arc frame for a dihedral measurement.
    ///
    /// The arc is centered at the midpoint of the p2–p3 bond. The flanking
    /// bond vectors are projected onto the plane perpendicular to the central
    /// bond to define the arc axes.
    fn for_dihedral(p1: Vec3, p2: Vec3, p3: Vec3, p4: Vec3, value_degrees: f64) -> Option<Self> {
        let bond_dir = p3 - p2;
        let bond_len = bond_dir.magnitude();
        if bond_len < 1e-6 {
            return None;
        }
        let bd = bond_dir / bond_len;

        let proj1 = project_perpendicular(p1 - p2, bd);
        let proj2 = project_perpendicular(p4 - p3, bd);
        let proj1_len = proj1.magnitude();
        let proj2_len = proj2.magnitude();
        if proj1_len < 1e-6 || proj2_len < 1e-6 {
            return None;
        }

        let radius = proj1_len.min(proj2_len) * ANGLE_SIZE;
        let (x_axis, y_axis) = perpendicular_frame(proj1 / proj1_len, proj2 / proj2_len, radius)?;

        Some(Self {
            center: (p2 + p3) * 0.5,
            x_axis,
            y_axis,
            arc_angle: (value_degrees as f32).to_radians(),
        })
    }

    /// Label position: at the half-angle, pushed out beyond the arc radius.
    fn label_position(&self) -> Vec3 {
        let half = self.arc_angle * 0.5;
        self.center + (self.x_axis * half.cos() + self.y_axis * half.sin()) * LABEL_OFFSET
    }
}

// ── Vector helpers ───────────────────────────────────────────────────────────

/// Project `v` onto the plane perpendicular to `axis`: `v - (v·axis)*axis`.
fn project_perpendicular(v: Vec3, axis: Vec3) -> Vec3 {
    v - axis * v.dot(axis)
}

/// Build a radius-scaled perpendicular frame from two unit directions.
///
/// Returns `(x_axis, y_axis)` where `x_axis = n1 * radius` and `y_axis`
/// is the component of `n2` perpendicular to `n1`, scaled to `radius`.
/// Returns `None` if the directions are collinear.
fn perpendicular_frame(n1: Vec3, n2: Vec3, radius: f32) -> Option<(Vec3, Vec3)> {
    let perp = n2 - n1 * n2.dot(n1);
    let perp_len = perp.magnitude();
    if perp_len < 1e-6 {
        return None;
    }
    Some((n1 * radius, perp / perp_len * radius))
}

// ── Measurement types ────────────────────────────────────────────────────────

/// Type of measurement
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum MeasurementType {
    /// Distance between two atoms (Å)
    Distance,
    /// Angle formed by three atoms (degrees)
    Angle,
    /// Dihedral angle formed by four atoms (degrees)
    Dihedral,
}

/// A single measurement entry
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Measurement {
    /// Type of measurement
    pub kind: MeasurementType,
    /// Atom positions (2 for distance, 3 for angle, 4 for dihedral)
    #[serde(with = "vec_vec3_serde")]
    pub points: Vec<Vec3>,
    /// Computed value (Å for distance, degrees for angle/dihedral)
    pub value: f64,
    /// Label text (formatted value)
    pub label: String,
    /// Label position in world space
    #[serde(with = "crate::serde_helpers::vec3_serde")]
    pub label_position: Vec3,
    /// Color for the dashed lines (RGBA)
    pub color: [f32; 4],
}

impl Measurement {
    /// Create a distance measurement between two points
    pub fn distance(p1: Vec3, p2: Vec3, color: [f32; 4]) -> Self {
        let value = (p2 - p1).magnitude() as f64;
        Self {
            kind: MeasurementType::Distance,
            points: vec![p1, p2],
            value,
            label: format!("{:.1}", value),
            label_position: (p1 + p2) * 0.5,
            color,
        }
    }

    /// Create an angle measurement between three points (vertex at p2)
    pub fn angle(p1: Vec3, p2: Vec3, p3: Vec3, color: [f32; 4]) -> Self {
        let v1 = p1 - p2;
        let v2 = p3 - p2;
        let len1 = v1.magnitude();
        let len2 = v2.magnitude();
        let value = if len1 > 1e-6 && len2 > 1e-6 {
            let cos_angle = (v1.dot(v2) / (len1 * len2)).clamp(-1.0, 1.0);
            (cos_angle.acos() as f64).to_degrees()
        } else {
            0.0
        };

        let arc = ArcGeometry::for_angle(p1, p2, p3, value);
        let label_position = arc.as_ref().map(|a| a.label_position()).unwrap_or_else(|| {
            // Degenerate: fall back to bisector
            let n1 = if len1 > 1e-6 {
                v1 / len1
            } else {
                Vec3::new(0.0, 0.0, 0.0)
            };
            let n2 = if len2 > 1e-6 {
                v2 / len2
            } else {
                Vec3::new(0.0, 0.0, 0.0)
            };
            p2 + (n1 + n2) * 0.75
        });

        Self {
            kind: MeasurementType::Angle,
            points: vec![p1, p2, p3],
            value,
            label: format!("{:.1}", value),
            label_position,
            color,
        }
    }

    /// Create a dihedral measurement between four points
    pub fn dihedral(p1: Vec3, p2: Vec3, p3: Vec3, p4: Vec3, color: [f32; 4]) -> Self {
        let d21 = p2 - p1;
        let d32 = p3 - p2;
        let d43 = p4 - p3;

        let d32_len = d32.magnitude();
        let value = if d32_len > 1e-6 {
            let n1 = d21.cross(d32);
            let n2 = d32.cross(d43);
            let n1_len = n1.magnitude();
            let n2_len = n2.magnitude();
            if n1_len > 1e-6 && n2_len > 1e-6 {
                let cos_angle = (n1.dot(n2) / (n1_len * n2_len)).clamp(-1.0, 1.0);
                let mut angle = (cos_angle.acos() as f64).to_degrees();
                if n2.dot(d32.cross(n1)) < 0.0 {
                    angle = -angle;
                }
                angle
            } else {
                0.0
            }
        } else {
            0.0
        };

        let arc = ArcGeometry::for_dihedral(p1, p2, p3, p4, value);
        let label_position = arc
            .as_ref()
            .map(|a| a.label_position())
            .unwrap_or_else(|| (p2 + p3) * 0.5);

        Self {
            kind: MeasurementType::Dihedral,
            points: vec![p1, p2, p3, p4],
            value,
            label: format!("{:.1}", value),
            label_position,
            color,
        }
    }
}

// ── Measurement object ───────────────────────────────────────────────────────

/// A collection of measurements as a scene object
#[derive(Serialize, Deserialize)]
pub struct MeasurementObject {
    name: String,
    state: ObjectState,
    measurements: Vec<Measurement>,
    #[serde(skip, default = "default_dirty_true")]
    dirty: bool,
}

fn default_dirty_true() -> bool {
    true
}

impl MeasurementObject {
    /// Create a new empty measurement object
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            state: ObjectState::default(),
            measurements: Vec::new(),
            dirty: true,
        }
    }

    /// Get the measurements
    pub fn measurements(&self) -> &[Measurement] {
        &self.measurements
    }

    /// Add a measurement
    pub fn add_measurement(&mut self, measurement: Measurement) {
        self.measurements.push(measurement);
        self.dirty = true;
    }

    /// Get the number of measurements
    pub fn len(&self) -> usize {
        self.measurements.len()
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.measurements.is_empty()
    }

    /// Collect labels for screen-space rendering
    pub fn collect_labels(&self) -> Vec<(Vec3, &str)> {
        self.measurements
            .iter()
            .map(|m| (m.label_position, m.label.as_str()))
            .collect()
    }

    /// Compute bounding box of all measurement positions
    fn compute_extent(&self) -> Option<(Vec3, Vec3)> {
        if self.measurements.is_empty() {
            return None;
        }

        let mut min = Vec3::new(f32::MAX, f32::MAX, f32::MAX);
        let mut max = Vec3::new(f32::MIN, f32::MIN, f32::MIN);

        for measurement in &self.measurements {
            for p in &measurement.points {
                min.x = min.x.min(p.x);
                min.y = min.y.min(p.y);
                min.z = min.z.min(p.z);
                max.x = max.x.max(p.x);
                max.y = max.y.max(p.y);
                max.z = max.z.max(p.z);
            }
        }

        Some((min, max))
    }
}

impl Object for MeasurementObject {
    fn name(&self) -> &str {
        &self.name
    }

    fn object_type(&self) -> ObjectType {
        ObjectType::Measurement
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

    fn set_name(&mut self, name: String) {
        self.name = name;
    }
}

/// Serde helpers for `Vec<Vec3>`.
mod vec_vec3_serde {
    use lin_alg::f32::Vec3;
    use serde::{Deserialize, Deserializer, Serialize, Serializer};

    #[derive(Serialize, Deserialize)]
    struct Vec3Proxy {
        x: f32,
        y: f32,
        z: f32,
    }

    pub fn serialize<S>(vec: &[Vec3], serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let proxies: Vec<Vec3Proxy> = vec
            .iter()
            .map(|v| Vec3Proxy {
                x: v.x,
                y: v.y,
                z: v.z,
            })
            .collect();
        proxies.serialize(serializer)
    }

    pub fn deserialize<'de, D>(deserializer: D) -> Result<Vec<Vec3>, D::Error>
    where
        D: Deserializer<'de>,
    {
        let proxies: Vec<Vec3Proxy> = Vec::deserialize(deserializer)?;
        Ok(proxies
            .into_iter()
            .map(|p| Vec3::new(p.x, p.y, p.z))
            .collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const DEFAULT_COLOR: [f32; 4] = [1.0, 1.0, 0.0, 1.0];
    const DASH_LENGTH: f32 = 0.15;
    const GAP_LENGTH: f32 = 0.10;

    #[derive(Debug, Clone, Copy)]
    struct DashVertex {
        position: [f32; 3],
    }

    fn dash_vertex(position: Vec3, _color: [f32; 4]) -> DashVertex {
        DashVertex {
            position: [position.x, position.y, position.z],
        }
    }

    fn generate_dash_vertices(p1: Vec3, p2: Vec3, color: [f32; 4]) -> Vec<DashVertex> {
        let delta = p2 - p1;
        let len = delta.magnitude();
        if len < 1e-6 {
            return Vec::new();
        }

        let direction = delta / len;
        let mut verts = Vec::new();
        let mut offset = 0.0;
        while offset < len {
            let dash_end = (offset + DASH_LENGTH).min(len);
            if dash_end > offset {
                verts.push(dash_vertex(p1 + direction * offset, color));
                verts.push(dash_vertex(p1 + direction * dash_end, color));
            }
            offset += DASH_LENGTH + GAP_LENGTH;
        }
        verts
    }

    fn generate_arc_dash_vertices(
        center: Vec3,
        x_axis: Vec3,
        y_axis: Vec3,
        arc_angle: f32,
        color: [f32; 4],
    ) -> Vec<DashVertex> {
        let radius = x_axis.magnitude();
        if radius < 1e-6 || arc_angle.abs() < 1e-6 {
            return Vec::new();
        }

        let sign = arc_angle.signum();
        let abs_angle = arc_angle.abs();
        let dash_angle = DASH_LENGTH / radius;
        let gap_angle = GAP_LENGTH / radius;
        let mut verts = Vec::new();
        let mut angle = 0.0;

        while angle < abs_angle {
            let end = (angle + dash_angle).min(abs_angle);
            let a0 = angle * sign;
            let a1 = end * sign;
            let p0 = center + x_axis * a0.cos() + y_axis * a0.sin();
            let p1 = center + x_axis * a1.cos() + y_axis * a1.sin();
            verts.push(dash_vertex(p0, color));
            verts.push(dash_vertex(p1, color));
            angle += dash_angle + gap_angle;
        }
        verts
    }

    trait MeasurementVertexTestExt {
        fn generate_vertices(&self) -> Vec<DashVertex>;
    }

    impl MeasurementVertexTestExt for Measurement {
        fn generate_vertices(&self) -> Vec<DashVertex> {
            match self.kind {
                MeasurementType::Distance => {
                    generate_dash_vertices(self.points[0], self.points[1], self.color)
                }
                MeasurementType::Angle => {
                    let mut verts =
                        generate_dash_vertices(self.points[0], self.points[1], self.color);
                    verts.extend(generate_dash_vertices(
                        self.points[1],
                        self.points[2],
                        self.color,
                    ));
                    if let Some(arc) = ArcGeometry::for_angle(
                        self.points[0],
                        self.points[1],
                        self.points[2],
                        self.value,
                    ) {
                        verts.extend(generate_arc_dash_vertices(
                            arc.center,
                            arc.x_axis,
                            arc.y_axis,
                            arc.arc_angle,
                            self.color,
                        ));
                    }
                    verts
                }
                MeasurementType::Dihedral => {
                    let mut verts =
                        generate_dash_vertices(self.points[0], self.points[1], self.color);
                    verts.extend(generate_dash_vertices(
                        self.points[1],
                        self.points[2],
                        self.color,
                    ));
                    verts.extend(generate_dash_vertices(
                        self.points[2],
                        self.points[3],
                        self.color,
                    ));
                    if let Some(arc) = ArcGeometry::for_dihedral(
                        self.points[0],
                        self.points[1],
                        self.points[2],
                        self.points[3],
                        self.value,
                    ) {
                        verts.extend(generate_arc_dash_vertices(
                            arc.center,
                            arc.x_axis,
                            arc.y_axis,
                            arc.arc_angle,
                            self.color,
                        ));
                    }
                    verts
                }
            }
        }
    }

    #[test]
    fn test_distance_measurement() {
        let p1 = Vec3::new(0.0, 0.0, 0.0);
        let p2 = Vec3::new(3.0, 4.0, 0.0);
        let m = Measurement::distance(p1, p2, DEFAULT_COLOR);
        assert!((m.value - 5.0).abs() < 1e-6);
        assert_eq!(m.label, "5.0");
        assert_eq!(m.points.len(), 2);
    }

    #[test]
    fn test_angle_measurement() {
        let p1 = Vec3::new(1.0, 0.0, 0.0);
        let p2 = Vec3::new(0.0, 0.0, 0.0);
        let p3 = Vec3::new(0.0, 1.0, 0.0);
        let m = Measurement::angle(p1, p2, p3, DEFAULT_COLOR);
        assert!((m.value - 90.0).abs() < 0.1);
        assert_eq!(m.points.len(), 3);
    }

    #[test]
    fn test_dihedral_measurement() {
        // Simple 90-degree dihedral
        let p1 = Vec3::new(1.0, 0.0, 0.0);
        let p2 = Vec3::new(0.0, 0.0, 0.0);
        let p3 = Vec3::new(0.0, 1.0, 0.0);
        let p4 = Vec3::new(0.0, 1.0, 1.0);
        let m = Measurement::dihedral(p1, p2, p3, p4, DEFAULT_COLOR);
        assert!((m.value.abs() - 90.0).abs() < 0.1);
        assert_eq!(m.points.len(), 4);
    }

    #[test]
    fn test_measurement_object() {
        let mut obj = MeasurementObject::new("dist01");
        assert!(obj.is_empty());

        let m = Measurement::distance(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(5.0, 0.0, 0.0),
            DEFAULT_COLOR,
        );
        obj.add_measurement(m);
        assert_eq!(obj.len(), 1);
        assert!(!obj.is_empty());
    }

    #[test]
    fn test_measurement_extent() {
        let mut obj = MeasurementObject::new("test");
        assert!(obj.extent().is_none());

        obj.add_measurement(Measurement::distance(
            Vec3::new(1.0, 2.0, 3.0),
            Vec3::new(4.0, 5.0, 6.0),
            DEFAULT_COLOR,
        ));

        let (min, max) = obj.extent().unwrap();
        assert_eq!(min, Vec3::new(1.0, 2.0, 3.0));
        assert_eq!(max, Vec3::new(4.0, 5.0, 6.0));
    }

    #[test]
    fn test_collect_labels() {
        let mut obj = MeasurementObject::new("test");
        obj.add_measurement(Measurement::distance(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(10.0, 0.0, 0.0),
            DEFAULT_COLOR,
        ));

        let labels = obj.collect_labels();
        assert_eq!(labels.len(), 1);
        assert_eq!(labels[0].0, Vec3::new(5.0, 0.0, 0.0)); // midpoint
        assert_eq!(labels[0].1, "10.0"); // distance
    }

    #[test]
    fn test_dash_generation() {
        let p1 = Vec3::new(0.0, 0.0, 0.0);
        let p2 = Vec3::new(1.0, 0.0, 0.0);
        let verts = generate_dash_vertices(p1, p2, DEFAULT_COLOR);
        // With dash=0.15, gap=0.10, stride=0.25, over 1.0 Å: 4 dashes
        assert!(verts.len() >= 8); // at least 4 dash segments * 2 vertices
        assert_eq!(verts.len() % 2, 0); // always pairs
    }

    #[test]
    fn test_generate_vertices() {
        // Distance: straight dashes only
        let m = Measurement::distance(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            DEFAULT_COLOR,
        );
        let verts = m.generate_vertices();
        assert!(!verts.is_empty());
        assert_eq!(verts.len() % 2, 0);

        // Angle: two bond dashes + arc dashes
        let m = Measurement::angle(
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            DEFAULT_COLOR,
        );
        let verts = m.generate_vertices();
        assert!(!verts.is_empty());
        assert_eq!(verts.len() % 2, 0);
        // Should have more vertices than just two straight dashes (arc adds extra)
        let straight_only = generate_dash_vertices(
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            DEFAULT_COLOR,
        )
        .len()
            + generate_dash_vertices(
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                DEFAULT_COLOR,
            )
            .len();
        assert!(verts.len() > straight_only);

        // Dihedral: three bond dashes + arc dashes
        let m = Measurement::dihedral(
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 1.0, 1.0),
            DEFAULT_COLOR,
        );
        let verts = m.generate_vertices();
        assert!(!verts.is_empty());
        assert_eq!(verts.len() % 2, 0);
    }

    #[test]
    fn test_arc_dash_vertices() {
        let center = Vec3::new(0.0, 0.0, 0.0);
        let x_axis = Vec3::new(1.0, 0.0, 0.0);
        let y_axis = Vec3::new(0.0, 1.0, 0.0);

        // 90-degree arc
        let verts = generate_arc_dash_vertices(
            center,
            x_axis,
            y_axis,
            std::f32::consts::FRAC_PI_2,
            DEFAULT_COLOR,
        );
        assert!(!verts.is_empty());
        assert_eq!(verts.len() % 2, 0);

        // First vertex should be near (1,0,0), last near (0,1,0)
        assert!((verts[0].position[0] - 1.0).abs() < 0.01);
        assert!(verts[0].position[1].abs() < 0.01);

        // Zero angle should produce no vertices
        let verts = generate_arc_dash_vertices(center, x_axis, y_axis, 0.0, DEFAULT_COLOR);
        assert!(verts.is_empty());
    }
}
