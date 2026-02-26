//! Measurement object for distance, angle, and dihedral visualization
//!
//! Measurement objects display dashed lines between atoms with labels showing
//! the measured value (distance in Å, angles in degrees).

use lin_alg::f32::Vec3;
use pymol_render::{LineRep, LineVertex, MeshRep, RenderContext, Representation};
use serde::{Deserialize, Serialize};

use super::{Object, ObjectState, ObjectType};

/// Default dash length in Angstroms
const DASH_LENGTH: f32 = 0.15;
/// Default gap length in Angstroms
const DASH_GAP: f32 = 0.10;
/// Arc radius as fraction of the shorter leg length
const ANGLE_SIZE: f32 = 0.4;
/// How much longer the solid wing lines extend beyond the arc radius
const WING_EXTENSION: f32 = 1.2;
/// Radius of solid indicator lines (Angstroms)
const SOLID_LINE_RADIUS: f32 = 0.05;
/// Number of sides for the solid line tube
const TUBE_SIDES: usize = 12;
/// Number of latitude rings for sphere caps
const SPHERE_RINGS: usize = 6;
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
    fn for_dihedral(
        p1: Vec3,
        p2: Vec3,
        p3: Vec3,
        p4: Vec3,
        value_degrees: f64,
    ) -> Option<Self> {
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
        let (x_axis, y_axis) =
            perpendicular_frame(proj1 / proj1_len, proj2 / proj2_len, radius)?;

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

/// Build a perpendicular frame (u, v) for a given normalized direction vector.
fn perp_frame(d: Vec3) -> (Vec3, Vec3) {
    let ref_vec = if d.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let u = d.cross(ref_vec).to_normalized();
    let v = d.cross(u);
    (u, v)
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
    /// Precomputed arc frame (None for distances; recomputed after deserialization)
    #[serde(skip)]
    arc: Option<ArcGeometry>,
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
            arc: None,
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
        let label_position = arc
            .as_ref()
            .map(|a| a.label_position())
            .unwrap_or_else(|| {
                // Degenerate: fall back to bisector
                let n1 = if len1 > 1e-6 { v1 / len1 } else { Vec3::new(0.0, 0.0, 0.0) };
                let n2 = if len2 > 1e-6 { v2 / len2 } else { Vec3::new(0.0, 0.0, 0.0) };
                p2 + (n1 + n2) * 0.75
            });

        Self {
            kind: MeasurementType::Angle,
            points: vec![p1, p2, p3],
            value,
            label: format!("{:.1}", value),
            label_position,
            color,
            arc,
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
            arc,
        }
    }

    /// Recompute arc geometry from stored points and value.
    /// Called after deserialization when the `#[serde(skip)]` arc field is None.
    fn rebuild_arc(&mut self) {
        self.arc = match self.kind {
            MeasurementType::Distance => None,
            MeasurementType::Angle => {
                ArcGeometry::for_angle(self.points[0], self.points[1], self.points[2], self.value)
            }
            MeasurementType::Dihedral => ArcGeometry::for_dihedral(
                self.points[0],
                self.points[1],
                self.points[2],
                self.points[3],
                self.value,
            ),
        };
    }

    /// Generate all dashed line vertices for this measurement, including arcs.
    fn generate_vertices(&self) -> Vec<LineVertex> {
        match self.kind {
            MeasurementType::Distance => {
                generate_dash_vertices(self.points[0], self.points[1], self.color)
            }
            MeasurementType::Angle => {
                let (p1, p2, p3) = (self.points[0], self.points[1], self.points[2]);
                let mut verts = Vec::new();
                verts.extend(generate_dash_vertices(p1, p2, self.color));
                verts.extend(generate_dash_vertices(p2, p3, self.color));
                if let Some(ref arc) = self.arc {
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
                let (p1, p2, p3, p4) =
                    (self.points[0], self.points[1], self.points[2], self.points[3]);
                let mut verts = Vec::new();
                verts.extend(generate_dash_vertices(p1, p2, self.color));
                verts.extend(generate_dash_vertices(p2, p3, self.color));
                verts.extend(generate_dash_vertices(p3, p4, self.color));
                if let Some(ref arc) = self.arc {
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

    /// Generate solid (thick) wing lines for dihedral measurements.
    fn generate_solid_lines(&self, mesh: &mut MeshRep) {
        if self.kind != MeasurementType::Dihedral {
            return;
        }
        let arc = match self.arc {
            Some(ref a) => a,
            None => return,
        };

        // Wing from center to arc start (extended beyond arc radius)
        let arc_start = arc.center + arc.x_axis * WING_EXTENSION;
        generate_solid_line_tube(mesh, arc.center, arc_start, self.color);

        // Wing from center to arc end (extended beyond arc radius)
        let end_offset =
            arc.x_axis * arc.arc_angle.cos() + arc.y_axis * arc.arc_angle.sin();
        let arc_end = arc.center + end_offset * WING_EXTENSION;
        generate_solid_line_tube(mesh, arc.center, arc_end, self.color);
    }
}

// ── Dashed line generation ───────────────────────────────────────────────────

/// Generate dashed line vertices from a line segment.
fn generate_dash_vertices(p1: Vec3, p2: Vec3, color: [f32; 4]) -> Vec<LineVertex> {
    let mut vertices = Vec::new();
    let dir = p2 - p1;
    let total_length = dir.magnitude();
    if total_length < 1e-6 {
        return vertices;
    }

    let n = dir / total_length;
    let stride = DASH_LENGTH + DASH_GAP;
    let mut t = 0.0f32;

    while t < total_length {
        let t_end = (t + DASH_LENGTH).min(total_length);
        let start = p1 + n * t;
        let end = p1 + n * t_end;
        vertices.push(LineVertex {
            position: [start.x, start.y, start.z],
            color,
        });
        vertices.push(LineVertex {
            position: [end.x, end.y, end.z],
            color,
        });
        t += stride;
    }

    vertices
}

/// Generate dashed arc vertices along a circular arc.
fn generate_arc_dash_vertices(
    center: Vec3,
    x_axis: Vec3,
    y_axis: Vec3,
    arc_angle: f32,
    color: [f32; 4],
) -> Vec<LineVertex> {
    let mut vertices = Vec::new();
    let radius = x_axis.magnitude();
    if radius < 1e-6 || arc_angle.abs() < 1e-6 {
        return vertices;
    }

    let arc_length = arc_angle.abs() * radius;
    let stride = DASH_LENGTH + DASH_GAP;
    let sign = arc_angle.signum();
    let mut s = 0.0f32;

    while s < arc_length {
        let s_end = (s + DASH_LENGTH).min(arc_length);
        let t0 = sign * s / radius;
        let t1 = sign * s_end / radius;

        let p0 = center + x_axis * t0.cos() + y_axis * t0.sin();
        let p1 = center + x_axis * t1.cos() + y_axis * t1.sin();

        vertices.push(LineVertex {
            position: [p0.x, p0.y, p0.z],
            color,
        });
        vertices.push(LineVertex {
            position: [p1.x, p1.y, p1.z],
            color,
        });

        s += stride;
    }

    vertices
}

// ── Solid tube generation ────────────────────────────────────────────────────

/// Generate a capped tube between two points as mesh triangles with smooth normals.
fn generate_solid_line_tube(mesh: &mut MeshRep, p1: Vec3, p2: Vec3, color: [f32; 4]) {
    let dir = p2 - p1;
    let dir_len = dir.magnitude();
    if dir_len < 1e-6 {
        return;
    }
    let d = dir / dir_len;
    let (u, v) = perp_frame(d);
    let r = SOLID_LINE_RADIUS;

    // Compute ring positions and radial normals
    let ring_pos_norm = |center: Vec3| -> Vec<([f32; 3], [f32; 3])> {
        (0..TUBE_SIDES)
            .map(|i| {
                let angle = std::f32::consts::TAU * i as f32 / TUBE_SIDES as f32;
                let normal = u * angle.cos() + v * angle.sin();
                let pos = center + normal * r;
                ([pos.x, pos.y, pos.z], [normal.x, normal.y, normal.z])
            })
            .collect()
    };

    let ring1 = ring_pos_norm(p1);
    let ring2 = ring_pos_norm(p2);

    // Side faces with smooth radial normals
    for i in 0..TUBE_SIDES {
        let j = (i + 1) % TUBE_SIDES;
        mesh.add_triangle_smooth(
            [ring1[i].0, ring2[i].0, ring1[j].0],
            [ring1[i].1, ring2[i].1, ring1[j].1],
            color,
        );
        mesh.add_triangle_smooth(
            [ring1[j].0, ring2[i].0, ring2[j].0],
            [ring1[j].1, ring2[i].1, ring2[j].1],
            color,
        );
    }

    // Hemisphere caps
    generate_hemisphere_cap(mesh, p2, d, u, v, r, color);
    generate_hemisphere_cap(mesh, p1, -d, u, v, r, color);
}

/// Generate a hemisphere cap as mesh triangles with smooth normals.
fn generate_hemisphere_cap(
    mesh: &mut MeshRep,
    center: Vec3,
    pole: Vec3,
    u: Vec3,
    v: Vec3,
    r: f32,
    color: [f32; 4],
) {
    let half_pi = std::f32::consts::FRAC_PI_2;

    let ring_at = |lat_idx: usize| -> Vec<([f32; 3], [f32; 3])> {
        let lat = half_pi * lat_idx as f32 / SPHERE_RINGS as f32;
        let cos_lat = lat.cos();
        let sin_lat = lat.sin();
        (0..TUBE_SIDES)
            .map(|i| {
                let lon = std::f32::consts::TAU * i as f32 / TUBE_SIDES as f32;
                let radial = u * lon.cos() + v * lon.sin();
                let normal = radial * cos_lat + pole * sin_lat;
                let pos = center + normal * r;
                ([pos.x, pos.y, pos.z], [normal.x, normal.y, normal.z])
            })
            .collect()
    };

    let mut prev_ring = ring_at(0);
    for lat in 1..SPHERE_RINGS {
        let curr_ring = ring_at(lat);
        for i in 0..TUBE_SIDES {
            let j = (i + 1) % TUBE_SIDES;
            mesh.add_triangle_smooth(
                [prev_ring[i].0, curr_ring[i].0, prev_ring[j].0],
                [prev_ring[i].1, curr_ring[i].1, prev_ring[j].1],
                color,
            );
            mesh.add_triangle_smooth(
                [prev_ring[j].0, curr_ring[i].0, curr_ring[j].0],
                [prev_ring[j].1, curr_ring[i].1, curr_ring[j].1],
                color,
            );
        }
        prev_ring = curr_ring;
    }

    // Triangle fan to pole
    let pole_pos = center + pole * r;
    let pole_pt = [pole_pos.x, pole_pos.y, pole_pos.z];
    let pole_normal = [pole.x, pole.y, pole.z];
    for i in 0..TUBE_SIDES {
        let j = (i + 1) % TUBE_SIDES;
        mesh.add_triangle_smooth(
            [prev_ring[i].0, pole_pt, prev_ring[j].0],
            [prev_ring[i].1, pole_normal, prev_ring[j].1],
            color,
        );
    }
}

// ── Measurement object ───────────────────────────────────────────────────────

/// A collection of measurements as a scene object
#[derive(Serialize, Deserialize)]
pub struct MeasurementObject {
    name: String,
    state: ObjectState,
    measurements: Vec<Measurement>,
    #[serde(skip)]
    dash_lines: Option<LineRep>,
    #[serde(skip)]
    solid_lines: Option<MeshRep>,
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
            dash_lines: None,
            solid_lines: None,
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

    /// Build dashed line geometry and solid wing line mesh from measurements
    fn build_geometry(&mut self) {
        // Ensure arcs are computed (needed after deserialization)
        for measurement in &mut self.measurements {
            if measurement.arc.is_none() && measurement.kind != MeasurementType::Distance {
                measurement.rebuild_arc();
            }
        }

        let lines = self.dash_lines.get_or_insert_with(LineRep::new);
        let mut vertices = Vec::new();
        for measurement in &self.measurements {
            vertices.extend(measurement.generate_vertices());
        }
        lines.set_line_data(vertices);

        let mesh = self.solid_lines.get_or_insert_with(MeshRep::new);
        mesh.clear();
        for measurement in &self.measurements {
            measurement.generate_solid_lines(mesh);
        }
    }

    /// Prepare for rendering (rebuild + upload if dirty)
    pub fn prepare_render(&mut self, context: &RenderContext) {
        if !self.dirty {
            return;
        }

        self.build_geometry();

        if let Some(ref mut lines) = self.dash_lines {
            lines.upload(context.device(), context.queue());
        }
        if let Some(ref mut mesh) = self.solid_lines {
            mesh.upload(context.device(), context.queue());
        }

        self.dirty = false;
    }

    /// Render the measurement dashed lines and solid wing lines
    pub fn render<'a>(
        &'a self,
        render_pass: &mut wgpu::RenderPass<'a>,
        context: &'a RenderContext,
    ) {
        if !self.state.enabled || self.is_empty() {
            return;
        }

        if let Some(ref lines) = self.dash_lines {
            if !lines.is_empty() {
                let pipeline =
                    context.line_pipeline(pymol_render::pipeline::BlendMode::Opaque);
                render_pass.set_pipeline(&pipeline);
                render_pass.set_bind_group(0, context.uniform_bind_group(), &[]);
                render_pass.set_bind_group(1, context.shadow_bind_group(), &[]);
                lines.render(render_pass);
            }
        }

        if let Some(ref mesh) = self.solid_lines {
            if !mesh.is_empty() {
                let pipeline =
                    context.mesh_pipeline(pymol_render::pipeline::BlendMode::Opaque);
                render_pass.set_pipeline(&pipeline);
                render_pass.set_bind_group(0, context.uniform_bind_group(), &[]);
                render_pass.set_bind_group(1, context.shadow_bind_group(), &[]);
                mesh.render(render_pass);
            }
        }
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

/// Serde helpers for Vec<Vec3>
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
