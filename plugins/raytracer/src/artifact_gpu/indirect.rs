use patinae_plugin::prelude::ViewerLike;
use patinae_scene::{GpuHandle, RenderArtifactRepKind};

use super::layout::DRAW_INDIRECT_SIZE;

pub(super) fn read_indirect_draw_args(
    viewer: &mut dyn ViewerLike,
    indirect: GpuHandle,
    rep_kind: RenderArtifactRepKind,
) -> Result<[u32; 4], String> {
    let bytes = viewer
        .gpu_read_buffer(indirect, 0, DRAW_INDIRECT_SIZE)
        .map_err(|err| format!("failed to read {rep_kind:?} triangle indirect draw args: {err}"))?;
    decode_draw_indirect_args(&bytes, rep_kind)
}

pub(super) fn decode_draw_indirect_args(
    bytes: &[u8],
    rep_kind: RenderArtifactRepKind,
) -> Result<[u32; 4], String> {
    if bytes.len() != DRAW_INDIRECT_SIZE as usize {
        return Err(format!(
            "{rep_kind:?} triangle indirect draw args read returned {} bytes, expected {DRAW_INDIRECT_SIZE}",
            bytes.len()
        ));
    }

    let mut args = [0_u32; 4];
    for (slot, chunk) in args.iter_mut().zip(bytes.chunks_exact(4)) {
        let bytes: [u8; 4] = chunk
            .try_into()
            .map_err(|_| "draw indirect chunk size mismatch".to_string())?;
        *slot = u32::from_le_bytes(bytes);
    }
    Ok(args)
}

pub(super) fn active_triangle_vertex_count(draw_vertex_count: u32, vertex_capacity: u32) -> u32 {
    draw_vertex_count.min(vertex_capacity) / 3 * 3
}
