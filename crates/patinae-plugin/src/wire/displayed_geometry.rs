use super::*;

pub fn displayed_geometry_to_wire(displayed: &DisplayedGeometry) -> WireDisplayedGeometry {
    WireDisplayedGeometry {
        objects: displayed
            .objects
            .iter()
            .map(displayed_object_geometry_to_wire)
            .collect(),
    }
}

/// Converts wire DTOs into displayed geometry.
pub fn displayed_geometry_from_wire(displayed: WireDisplayedGeometry) -> DisplayedGeometry {
    DisplayedGeometry {
        objects: displayed
            .objects
            .into_iter()
            .map(displayed_object_geometry_from_wire)
            .collect(),
    }
}

/// Converts one displayed-geometry chunk into host geometry.
pub fn displayed_geometry_chunk_from_wire(chunk: WireDisplayedGeometryChunk) -> DisplayedGeometry {
    DisplayedGeometry {
        objects: chunk
            .objects
            .into_iter()
            .map(displayed_object_geometry_from_wire)
            .collect(),
    }
}

/// Converts one displayed object into its wire DTO.
pub fn displayed_object_geometry_to_wire(
    object: &DisplayedObjectGeometry,
) -> WireDisplayedObjectGeometry {
    WireDisplayedObjectGeometry {
        object_id: object.object_id.0,
        primitives: object
            .primitives
            .iter()
            .map(displayed_primitive_to_wire)
            .collect(),
    }
}

fn displayed_object_geometry_from_wire(
    object: WireDisplayedObjectGeometry,
) -> DisplayedObjectGeometry {
    DisplayedObjectGeometry {
        object_id: ObjectId(object.object_id),
        primitives: object
            .primitives
            .into_iter()
            .map(displayed_primitive_from_wire)
            .collect(),
    }
}

/// Converts one displayed primitive into its wire DTO.
pub fn displayed_primitive_to_wire(primitive: &DisplayedPrimitive) -> WireDisplayedPrimitive {
    match primitive {
        DisplayedPrimitive::Mesh { rep, mesh } => WireDisplayedPrimitive::Mesh {
            rep: rep.as_raw(),
            vertices: mesh
                .vertices
                .iter()
                .map(displayed_mesh_vertex_to_wire)
                .collect(),
        },
        DisplayedPrimitive::Sphere {
            rep,
            owner_atom_id,
            center,
            radius,
            material,
        } => WireDisplayedPrimitive::Sphere {
            rep: rep.as_raw(),
            owner_atom_id: *owner_atom_id,
            center: *center,
            radius: *radius,
            material: displayed_material_to_wire(*material),
        },
        DisplayedPrimitive::Cylinder {
            rep,
            owner_atom_ids,
            start,
            end,
            radius,
            material_start,
            material_end,
        } => WireDisplayedPrimitive::Cylinder {
            rep: rep.as_raw(),
            owner_atom_ids: *owner_atom_ids,
            start: *start,
            end: *end,
            radius: *radius,
            material_start: displayed_material_to_wire(*material_start),
            material_end: displayed_material_to_wire(*material_end),
        },
        DisplayedPrimitive::LineSegment {
            rep,
            owner_atom_ids,
            start,
            end,
            width_px,
            material_start,
            material_end,
        } => WireDisplayedPrimitive::LineSegment {
            rep: rep.as_raw(),
            owner_atom_ids: *owner_atom_ids,
            start: *start,
            end: *end,
            width_px: *width_px,
            material_start: displayed_material_to_wire(*material_start),
            material_end: displayed_material_to_wire(*material_end),
        },
        DisplayedPrimitive::PointSample {
            rep,
            owner_atom_id,
            position,
            radius_px,
            material,
        } => WireDisplayedPrimitive::PointSample {
            rep: rep.as_raw(),
            owner_atom_id: *owner_atom_id,
            position: *position,
            radius_px: *radius_px,
            material: displayed_material_to_wire(*material),
        },
    }
}

fn displayed_primitive_from_wire(primitive: WireDisplayedPrimitive) -> DisplayedPrimitive {
    match primitive {
        WireDisplayedPrimitive::Mesh { rep, vertices } => DisplayedPrimitive::Mesh {
            rep: RepKind::from_raw(rep),
            mesh: DisplayedMesh {
                vertices: vertices
                    .into_iter()
                    .map(displayed_mesh_vertex_from_wire)
                    .collect(),
            },
        },
        WireDisplayedPrimitive::Sphere {
            rep,
            owner_atom_id,
            center,
            radius,
            material,
        } => DisplayedPrimitive::Sphere {
            rep: RepKind::from_raw(rep),
            owner_atom_id,
            center,
            radius,
            material: displayed_material_from_wire(material),
        },
        WireDisplayedPrimitive::Cylinder {
            rep,
            owner_atom_ids,
            start,
            end,
            radius,
            material_start,
            material_end,
        } => DisplayedPrimitive::Cylinder {
            rep: RepKind::from_raw(rep),
            owner_atom_ids,
            start,
            end,
            radius,
            material_start: displayed_material_from_wire(material_start),
            material_end: displayed_material_from_wire(material_end),
        },
        WireDisplayedPrimitive::LineSegment {
            rep,
            owner_atom_ids,
            start,
            end,
            width_px,
            material_start,
            material_end,
        } => DisplayedPrimitive::LineSegment {
            rep: RepKind::from_raw(rep),
            owner_atom_ids,
            start,
            end,
            width_px,
            material_start: displayed_material_from_wire(material_start),
            material_end: displayed_material_from_wire(material_end),
        },
        WireDisplayedPrimitive::PointSample {
            rep,
            owner_atom_id,
            position,
            radius_px,
            material,
        } => DisplayedPrimitive::PointSample {
            rep: RepKind::from_raw(rep),
            owner_atom_id,
            position,
            radius_px,
            material: displayed_material_from_wire(material),
        },
    }
}

fn displayed_mesh_vertex_to_wire(vertex: &DisplayedMeshVertex) -> WireDisplayedMeshVertex {
    WireDisplayedMeshVertex {
        position: vertex.position,
        normal: vertex.normal,
        owner_atom_id: vertex.owner_atom_id,
        material: displayed_material_to_wire(vertex.material),
        flags: vertex.flags,
    }
}

fn displayed_mesh_vertex_from_wire(vertex: WireDisplayedMeshVertex) -> DisplayedMeshVertex {
    DisplayedMeshVertex {
        position: vertex.position,
        normal: vertex.normal,
        owner_atom_id: vertex.owner_atom_id,
        material: displayed_material_from_wire(vertex.material),
        flags: vertex.flags,
    }
}

fn displayed_material_to_wire(material: DisplayedMaterial) -> WireDisplayedMaterial {
    WireDisplayedMaterial {
        base_rgba: material.base_rgba,
        rep_rgba: material.rep_rgba,
        rgba: material.rgba,
        transparency: material.transparency,
    }
}

fn displayed_material_from_wire(material: WireDisplayedMaterial) -> DisplayedMaterial {
    DisplayedMaterial {
        base_rgba: material.base_rgba,
        rep_rgba: material.rep_rgba,
        rgba: material.rgba,
        transparency: material.transparency,
    }
}
