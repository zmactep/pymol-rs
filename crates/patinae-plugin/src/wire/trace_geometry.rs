use super::*;

/// Converts compact trace geometry into wire DTOs.
pub fn trace_geometry_chunk_to_wire(chunk: TraceGeometryChunk) -> WireTraceGeometryChunk {
    WireTraceGeometryChunk {
        spheres: chunk
            .spheres
            .into_iter()
            .map(trace_sphere_to_wire)
            .collect(),
        cylinders: chunk
            .cylinders
            .into_iter()
            .map(trace_cylinder_to_wire)
            .collect(),
        triangles: chunk
            .triangles
            .into_iter()
            .map(trace_triangle_to_wire)
            .collect(),
        line_segments: chunk
            .line_segments
            .into_iter()
            .map(trace_line_segment_to_wire)
            .collect(),
        point_samples: chunk
            .point_samples
            .into_iter()
            .map(trace_point_sample_to_wire)
            .collect(),
    }
}

/// Converts compact trace wire DTOs into host trace geometry.
pub fn trace_geometry_chunk_from_wire(chunk: WireTraceGeometryChunk) -> TraceGeometryChunk {
    TraceGeometryChunk {
        spheres: chunk
            .spheres
            .into_iter()
            .map(trace_sphere_from_wire)
            .collect(),
        cylinders: chunk
            .cylinders
            .into_iter()
            .map(trace_cylinder_from_wire)
            .collect(),
        triangles: chunk
            .triangles
            .into_iter()
            .map(trace_triangle_from_wire)
            .collect(),
        line_segments: chunk
            .line_segments
            .into_iter()
            .map(trace_line_segment_from_wire)
            .collect(),
        point_samples: chunk
            .point_samples
            .into_iter()
            .map(trace_point_sample_from_wire)
            .collect(),
    }
}

fn trace_material_to_wire(material: TraceMaterial) -> WireTraceMaterial {
    WireTraceMaterial {
        rgba: material.rgba,
        transparency: material.transparency,
    }
}

fn trace_material_from_wire(material: WireTraceMaterial) -> TraceMaterial {
    TraceMaterial {
        rgba: material.rgba,
        transparency: material.transparency,
    }
}

fn trace_sphere_to_wire(sphere: TraceSphere) -> WireTraceSphere {
    WireTraceSphere {
        center: sphere.center,
        radius: sphere.radius,
        material: trace_material_to_wire(sphere.material),
    }
}

fn trace_sphere_from_wire(sphere: WireTraceSphere) -> TraceSphere {
    TraceSphere {
        center: sphere.center,
        radius: sphere.radius,
        material: trace_material_from_wire(sphere.material),
    }
}

fn trace_cylinder_to_wire(cylinder: TraceCylinder) -> WireTraceCylinder {
    WireTraceCylinder {
        start: cylinder.start,
        end: cylinder.end,
        radius: cylinder.radius,
        material_start: trace_material_to_wire(cylinder.material_start),
        material_end: trace_material_to_wire(cylinder.material_end),
    }
}

fn trace_cylinder_from_wire(cylinder: WireTraceCylinder) -> TraceCylinder {
    TraceCylinder {
        start: cylinder.start,
        end: cylinder.end,
        radius: cylinder.radius,
        material_start: trace_material_from_wire(cylinder.material_start),
        material_end: trace_material_from_wire(cylinder.material_end),
    }
}

fn trace_triangle_to_wire(triangle: TraceTriangle) -> WireTraceTriangle {
    WireTraceTriangle {
        positions: triangle.positions,
        normals: triangle.normals,
        material: trace_material_to_wire(triangle.material),
    }
}

fn trace_triangle_from_wire(triangle: WireTraceTriangle) -> TraceTriangle {
    TraceTriangle {
        positions: triangle.positions,
        normals: triangle.normals,
        material: trace_material_from_wire(triangle.material),
    }
}

fn trace_line_segment_to_wire(line: TraceLineSegment) -> WireTraceLineSegment {
    WireTraceLineSegment {
        start: line.start,
        end: line.end,
        width_px: line.width_px,
        material_start: trace_material_to_wire(line.material_start),
        material_end: trace_material_to_wire(line.material_end),
    }
}

fn trace_line_segment_from_wire(line: WireTraceLineSegment) -> TraceLineSegment {
    TraceLineSegment {
        start: line.start,
        end: line.end,
        width_px: line.width_px,
        material_start: trace_material_from_wire(line.material_start),
        material_end: trace_material_from_wire(line.material_end),
    }
}

fn trace_point_sample_to_wire(point: TracePointSample) -> WireTracePointSample {
    WireTracePointSample {
        position: point.position,
        radius_px: point.radius_px,
        material: trace_material_to_wire(point.material),
    }
}

fn trace_point_sample_from_wire(point: WireTracePointSample) -> TracePointSample {
    TracePointSample {
        position: point.position,
        radius_px: point.radius_px,
        material: trace_material_from_wire(point.material),
    }
}
