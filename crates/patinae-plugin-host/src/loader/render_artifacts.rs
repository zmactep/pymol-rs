use super::*;

pub(super) fn render_artifact_buffer_usage(role: RenderArtifactBufferRole) -> GpuBufferUsage {
    match role {
        RenderArtifactBufferRole::FrameUniforms => GpuBufferUsage::UNIFORM,
        RenderArtifactBufferRole::SphereInstances
        | RenderArtifactBufferRole::StickInstances
        | RenderArtifactBufferRole::LineInstances
        | RenderArtifactBufferRole::StdVertices => {
            GpuBufferUsage::VERTEX.union(GpuBufferUsage::STORAGE)
        }
        RenderArtifactBufferRole::IndirectDraw => {
            GpuBufferUsage::INDIRECT.union(GpuBufferUsage::STORAGE)
        }
        RenderArtifactBufferRole::InstanceCount => GpuBufferUsage::STORAGE,
        RenderArtifactBufferRole::SceneAtoms
        | RenderArtifactBufferRole::SceneCoords
        | RenderArtifactBufferRole::SceneBonds
        | RenderArtifactBufferRole::SceneColorLut
        | RenderArtifactBufferRole::SceneMaskLut
        | RenderArtifactBufferRole::SceneMarkerLut
        | RenderArtifactBufferRole::SceneCsrOffsets
        | RenderArtifactBufferRole::SceneCsrIndices
        | RenderArtifactBufferRole::SceneObjectTable => GpuBufferUsage::STORAGE,
    }
}

pub(super) fn map_render_artifact_buffer_role(
    role: patinae_render::RenderArtifactBufferRole,
) -> RenderArtifactBufferRole {
    match role {
        patinae_render::RenderArtifactBufferRole::FrameUniforms => {
            RenderArtifactBufferRole::FrameUniforms
        }
        patinae_render::RenderArtifactBufferRole::SceneAtoms => {
            RenderArtifactBufferRole::SceneAtoms
        }
        patinae_render::RenderArtifactBufferRole::SceneCoords => {
            RenderArtifactBufferRole::SceneCoords
        }
        patinae_render::RenderArtifactBufferRole::SceneBonds => {
            RenderArtifactBufferRole::SceneBonds
        }
        patinae_render::RenderArtifactBufferRole::SceneColorLut => {
            RenderArtifactBufferRole::SceneColorLut
        }
        patinae_render::RenderArtifactBufferRole::SceneMaskLut => {
            RenderArtifactBufferRole::SceneMaskLut
        }
        patinae_render::RenderArtifactBufferRole::SceneMarkerLut => {
            RenderArtifactBufferRole::SceneMarkerLut
        }
        patinae_render::RenderArtifactBufferRole::SceneCsrOffsets => {
            RenderArtifactBufferRole::SceneCsrOffsets
        }
        patinae_render::RenderArtifactBufferRole::SceneCsrIndices => {
            RenderArtifactBufferRole::SceneCsrIndices
        }
        patinae_render::RenderArtifactBufferRole::SceneObjectTable => {
            RenderArtifactBufferRole::SceneObjectTable
        }
        patinae_render::RenderArtifactBufferRole::SphereInstances => {
            RenderArtifactBufferRole::SphereInstances
        }
        patinae_render::RenderArtifactBufferRole::StickInstances => {
            RenderArtifactBufferRole::StickInstances
        }
        patinae_render::RenderArtifactBufferRole::LineInstances => {
            RenderArtifactBufferRole::LineInstances
        }
        patinae_render::RenderArtifactBufferRole::StdVertices => {
            RenderArtifactBufferRole::StdVertices
        }
        patinae_render::RenderArtifactBufferRole::InstanceCount => {
            RenderArtifactBufferRole::InstanceCount
        }
        patinae_render::RenderArtifactBufferRole::IndirectDraw => {
            RenderArtifactBufferRole::IndirectDraw
        }
    }
}

pub(super) fn map_render_artifact_topology(
    topology: patinae_render::RenderArtifactPrimitiveTopology,
) -> RenderArtifactPrimitiveTopology {
    match topology {
        patinae_render::RenderArtifactPrimitiveTopology::SphereInstances => {
            RenderArtifactPrimitiveTopology::SphereInstances
        }
        patinae_render::RenderArtifactPrimitiveTopology::CylinderInstances => {
            RenderArtifactPrimitiveTopology::CylinderInstances
        }
        patinae_render::RenderArtifactPrimitiveTopology::LineInstances => {
            RenderArtifactPrimitiveTopology::LineInstances
        }
        patinae_render::RenderArtifactPrimitiveTopology::TriangleList => {
            RenderArtifactPrimitiveTopology::TriangleList
        }
        patinae_render::RenderArtifactPrimitiveTopology::LineList => {
            RenderArtifactPrimitiveTopology::LineList
        }
    }
}

pub(super) fn map_render_artifact_rep_kind(kind: patinae_render::RepKind) -> RenderArtifactRepKind {
    match kind {
        patinae_render::RepKind::Sphere => RenderArtifactRepKind::Sphere,
        patinae_render::RepKind::Stick => RenderArtifactRepKind::Stick,
        patinae_render::RepKind::Line => RenderArtifactRepKind::Line,
        patinae_render::RepKind::Cartoon => RenderArtifactRepKind::Cartoon,
        patinae_render::RepKind::Ribbon => RenderArtifactRepKind::Ribbon,
        patinae_render::RepKind::Surface => RenderArtifactRepKind::Surface,
        patinae_render::RepKind::Mesh => RenderArtifactRepKind::Mesh,
        patinae_render::RepKind::Dot => RenderArtifactRepKind::Dot,
        patinae_render::RepKind::Ellipsoid => RenderArtifactRepKind::Ellipsoid,
        _ => RenderArtifactRepKind::Other,
    }
}
