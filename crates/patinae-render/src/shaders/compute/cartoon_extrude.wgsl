// Run-aware cartoon extrusion.
//
// Run-based dispatch — each thread looks up its `RunDescriptor` via a
// binary-search-like lookup on the `vertex_offset` table, then branches
// on `car_type` to apply per-CartoonType emission rules:
//
//   - Loop/Oval (`extrude_tube`): perimeter-quad strip + triangle-fan
//     caps at run start/end.
//   - Rect (sheet body, no arrow): 4 face strips per pair + 4-vert
//     back/forward caps.
//   - Arrow (sheet body + arrow head): body up to body_end + arrow with
//     FIXED base frame + back cap + shoulder caps.
//   - NucleicRect: same as Rect.
//
// Frames (T, N, B) are computed inline from ExtrudePoint positions:
//   T = averaged unit chord
//   B = normalize(T × orientation)
//   N = B × T

// {{INCLUDE_OCTAHEDRAL}}

const SS_LOOP: u32 = 0u;
const SS_OVAL: u32 = 1u;        // CartoonType::Oval
const SS_RECT: u32 = 2u;        // CartoonType::Rect
const SS_ARROW: u32 = 3u;       // CartoonType::Arrow
const SS_NUCLEIC_RECT: u32 = 4u; // CartoonType::NucleicRect (matches Rust enum)

const TWO_SIDED_FLAG: u32 = 1u << 0u;

// Profile constants for cartoon cross-sections.
const HELIX_WIDTH: f32 = 0.25;
const HELIX_HEIGHT: f32 = 1.35;
const SHEET_WIDTH: f32 = 0.4;     // full thickness; half-thick = 0.2
const SHEET_HEIGHT: f32 = 1.4;    // wide half-extent
const LOOP_RADIUS: f32 = 0.2;
const ARROW_TIP_SCALE: f32 = 1.5;

struct ExtrudeParams {
    n_runs: u32,
    n_samples: u32,
    n_atoms: u32,
    // Tube profile segment count. MUST equal the `geom.quality` used by
    // `compute_run_vertex_layout` on the CPU side; otherwise the
    // `local_vidx → (pair, k)` decode below drifts vs. the vertex-slot
    // allocation and the tube develops gaps.
    quality: u32,
};

struct ExtrudePoint {
    position: vec3<f32>,
    atom_idx: u32,
    orientation: vec3<f32>,
    _pad: u32,
};

struct RunDescriptor {
    car_type: u32,
    sample_start: u32,
    sample_end: u32,        // INCLUSIVE
    vertex_offset: u32,
    vertex_count: u32,
    body_end: u32,          // body sample count (within run); 0 if no arrow
    flags: u32,
    _pad: u32,
};

struct StdVertex {
    px: f32,
    py: f32,
    pz: f32,
    normal_oct: u32,
    group_id: u32,
    flags: u32,
};

@group(0) @binding(0) var<uniform> params: ExtrudeParams;
@group(0) @binding(1) var<storage, read> extrude_points: array<ExtrudePoint>;
@group(0) @binding(2) var<storage, read> runs: array<RunDescriptor>;
@group(0) @binding(3) var<storage, read_write> vertices: array<StdVertex>;

fn make_vertex(p: vec3<f32>, n: vec3<f32>, group_id: u32, flags: u32) -> StdVertex {
    var v: StdVertex;
    v.px = p.x;
    v.py = p.y;
    v.pz = p.z;
    v.normal_oct = oct_encode(n);
    v.group_id = group_id;
    v.flags = flags;
    return v;
}

fn make_degen() -> StdVertex {
    return make_vertex(vec3<f32>(0.0, 0.0, 0.0), vec3<f32>(0.0, 0.0, 1.0), 0u, 0u);
}

fn normalize_safe(v: vec3<f32>) -> vec3<f32> {
    let l = length(v);
    if (l > 1e-6) {
        return v / l;
    }
    return vec3<f32>(0.0, 0.0, 1.0);
}

// Compute averaged-tangent at sample `s` within a run [run_start..=run_end].
fn tangent_at(s: u32, run_start: u32, run_end: u32) -> vec3<f32> {
    let p_self = extrude_points[s].position;
    if (s == run_start) {
        let p_next = extrude_points[s + 1u].position;
        return normalize_safe(p_next - p_self);
    }
    if (s == run_end) {
        let p_prev = extrude_points[s - 1u].position;
        return normalize_safe(p_self - p_prev);
    }
    let p_next = extrude_points[s + 1u].position;
    let p_prev = extrude_points[s - 1u].position;
    let d_next = normalize_safe(p_next - p_self);
    let d_prev = normalize_safe(p_self - p_prev);
    return normalize_safe(d_next + d_prev);
}

struct Frame {
    position: vec3<f32>,
    tangent: vec3<f32>,
    normal: vec3<f32>,
    binormal: vec3<f32>,
};

fn frame_at(s: u32, run_start: u32, run_end: u32) -> Frame {
    var f: Frame;
    f.position = extrude_points[s].position;
    f.tangent = tangent_at(s, run_start, run_end);
    let orient = extrude_points[s].orientation;
    f.binormal = normalize_safe(cross(f.tangent, orient));
    f.normal = cross(f.binormal, f.tangent);
    return f;
}

// Loop tubes get a tiny local centerline pass in the extrusion shader. It is
// endpoint-preserving and deliberately not shared with Oval/sheet/arrow, whose
// silhouettes come from different cartoon semantics.
fn loop_center_at(run: RunDescriptor, s: u32) -> vec3<f32> {
    let p = extrude_points[s].position;
    if (s == run.sample_start || s == run.sample_end) {
        return p;
    }
    return 0.25 * extrude_points[s - 1u].position + 0.5 * p + 0.25 * extrude_points[s + 1u].position;
}

fn loop_tangent_at(run: RunDescriptor, s: u32) -> vec3<f32> {
    let p_self = loop_center_at(run, s);
    if (s == run.sample_start) {
        return normalize_safe(loop_center_at(run, s + 1u) - p_self);
    }
    if (s == run.sample_end) {
        return normalize_safe(p_self - loop_center_at(run, s - 1u));
    }
    let d_next = normalize_safe(loop_center_at(run, s + 1u) - p_self);
    let d_prev = normalize_safe(p_self - loop_center_at(run, s - 1u));
    return normalize_safe(d_next + d_prev);
}

fn loop_frame_at(run: RunDescriptor, s: u32) -> Frame {
    var f: Frame;
    f.position = loop_center_at(run, s);
    f.tangent = loop_tangent_at(run, s);
    let orient = extrude_points[s].orientation;
    f.binormal = normalize_safe(cross(f.tangent, orient));
    f.normal = cross(f.binormal, f.tangent);
    return f;
}

fn tube_frame_at(run: RunDescriptor, s: u32, is_oval: bool) -> Frame {
    if (is_oval) {
        return frame_at(s, run.sample_start, run.sample_end);
    }
    return loop_frame_at(run, s);
}

// ============================================================================
// Profile evaluators — return (offset_2d, normal_2d) for vertex k of profile.
// ============================================================================

struct ProfileVert {
    offset: vec2<f32>,  // (x along normal, y along binormal)
    normal: vec2<f32>,  // surface normal in same basis
};

fn profile_circle(k: u32, segments: u32, radius: f32) -> ProfileVert {
    let theta = 6.2831853 * f32(k) / f32(segments);
    let c = cos(theta);
    let s = sin(theta);
    var v: ProfileVert;
    v.offset = vec2<f32>(radius * c, radius * s);
    v.normal = vec2<f32>(c, s);
    return v;
}

fn profile_ellipse(k: u32, segments: u32, half_normal: f32, half_binormal: f32) -> ProfileVert {
    let theta = 6.2831853 * f32(k) / f32(segments);
    let c = cos(theta);
    let s = sin(theta);
    var v: ProfileVert;
    v.offset = vec2<f32>(half_normal * s, half_binormal * c);
    let nx = half_binormal * s;
    let ny = half_normal * c;
    let nl = max(sqrt(nx * nx + ny * ny), 1e-6);
    v.normal = vec2<f32>(nx, ny) / nl;
    return v;
}

fn tube_profile_normal_at(run: RunDescriptor, s: u32, k: u32, is_oval: bool) -> vec3<f32> {
    let f = tube_frame_at(run, s, is_oval);
    let pv = profile_for_run(k, params.quality, is_oval);
    return pv.normal.x * f.normal + pv.normal.y * f.binormal;
}

fn smoothed_loop_profile_normal_at(run: RunDescriptor, s: u32, k: u32) -> vec3<f32> {
    var n = 2.0 * tube_profile_normal_at(run, s, k, false);
    if (s > run.sample_start) {
        n = n + tube_profile_normal_at(run, s - 1u, k, false);
    }
    if (s < run.sample_end) {
        n = n + tube_profile_normal_at(run, s + 1u, k, false);
    }
    return normalize_safe(n);
}

// ============================================================================
// Per-CartoonType vertex emission. Each function takes a global vertex index
// within the run [0..vertex_count) and emits one StdVertex.
// ============================================================================

fn emit_tube_vertex(
    run: RunDescriptor,
    local_vidx: u32,
    is_oval: bool,
    out_idx: u32,
) {
    let q = params.quality;
    let sample_count = run.sample_end - run.sample_start + 1u;
    let pair_count = sample_count - 1u;
    let body_total = pair_count * q * 6u;
    let cap_total = q * 3u;

    if (local_vidx < body_total) {
        // Body: one perimeter quad per (pair, k). 6 verts per quad.
        let body_idx = local_vidx;
        let pair_idx = body_idx / (q * 6u);
        let in_pair = body_idx - pair_idx * (q * 6u);
        let k = in_pair / 6u;
        let tri_v = in_pair - k * 6u;
        let s0 = run.sample_start + pair_idx;
        let s1 = s0 + 1u;
        let f0 = tube_frame_at(run, s0, is_oval);
        let f1 = tube_frame_at(run, s1, is_oval);
        let pv0_a = profile_for_run(k,        q, is_oval);
        let pv0_b = profile_for_run(k + 1u,   q, is_oval);
        let pv1_a = profile_for_run(k,        q, is_oval);
        let pv1_b = profile_for_run(k + 1u,   q, is_oval);
        // Position lifted via each frame's (N, B).
        let p00 = f0.position + pv0_a.offset.x * f0.normal + pv0_a.offset.y * f0.binormal;
        let p01 = f0.position + pv0_b.offset.x * f0.normal + pv0_b.offset.y * f0.binormal;
        let p10 = f1.position + pv1_a.offset.x * f1.normal + pv1_a.offset.y * f1.binormal;
        let p11 = f1.position + pv1_b.offset.x * f1.normal + pv1_b.offset.y * f1.binormal;
        let id0 = extrude_points[s0].atom_idx;
        // Triangulation: (p00, p10, p11), (p00, p11, p01).
        var pos: vec3<f32>;
        var nrm: vec3<f32>;
        var atom: u32;
        var normal_s: u32;
        var normal_k: u32;
        switch tri_v {
            case 0u, 3u: { pos = p00; normal_s = s0; normal_k = k; atom = id0; }
            case 1u: { pos = p10; normal_s = s1; normal_k = k; atom = id0; }
            case 2u, 4u: { pos = p11; normal_s = s1; normal_k = k + 1u; atom = id0; }
            case 5u: { pos = p01; normal_s = s0; normal_k = k + 1u; atom = id0; }
            default: { pos = p00; normal_s = s0; normal_k = k; atom = id0; }
        }
        if (is_oval) {
            nrm = tube_profile_normal_at(run, normal_s, normal_k, true);
        } else {
            nrm = smoothed_loop_profile_normal_at(run, normal_s, normal_k);
        }
        // Closed tube: leave TWO_SIDED clear (single-sided Lambert).
        vertices[out_idx] = make_vertex(pos, nrm, atom, 0u);
        return;
    }

    let cap_local = local_vidx - body_total;
    if (cap_local < cap_total) {
        // Start cap (triangle fan from center to perimeter at first sample).
        let s = run.sample_start;
        let f = tube_frame_at(run, s, is_oval);
        let k = cap_local / 3u;
        let tri_v = cap_local - k * 3u;
        let pv_a = profile_for_run(k,        q, is_oval);
        let pv_b = profile_for_run(k + 1u,   q, is_oval);
        let cap_n = -f.tangent;     // backward (closing the start)
        var pos: vec3<f32>;
        switch tri_v {
            case 0u: { pos = f.position; }
            case 1u: { pos = f.position + pv_b.offset.x * f.normal + pv_b.offset.y * f.binormal; }
            case 2u: { pos = f.position + pv_a.offset.x * f.normal + pv_a.offset.y * f.binormal; }
            default: { pos = f.position; }
        }
        vertices[out_idx] = make_vertex(pos, cap_n, extrude_points[s].atom_idx, 0u);
        return;
    }

    // End cap.
    let end_local = cap_local - cap_total;
    let s = run.sample_end;
    let f = tube_frame_at(run, s, is_oval);
    let k = end_local / 3u;
    let tri_v = end_local - k * 3u;
    let pv_a = profile_for_run(k,        q, is_oval);
    let pv_b = profile_for_run(k + 1u,   q, is_oval);
    let cap_n = f.tangent;          // forward (closing the end)
    var pos: vec3<f32>;
    switch tri_v {
        case 0u: { pos = f.position; }
        case 1u: { pos = f.position + pv_a.offset.x * f.normal + pv_a.offset.y * f.binormal; }
        case 2u: { pos = f.position + pv_b.offset.x * f.normal + pv_b.offset.y * f.binormal; }
        default: { pos = f.position; }
    }
    vertices[out_idx] = make_vertex(pos, cap_n, extrude_points[s].atom_idx, 0u);
}

fn profile_for_run(k: u32, segments: u32, is_oval: bool) -> ProfileVert {
    let kk = k % segments;
    if (is_oval) {
        return profile_ellipse(kk, segments, HELIX_WIDTH, HELIX_HEIGHT);
    }
    return profile_circle(kk, segments, LOOP_RADIUS);
}

// ============================================================================
// Sheet body / arrow / caps emission. Mirror generate_explicit_sheet.
// ============================================================================

struct SheetCorners {
    tl: vec3<f32>,
    tr: vec3<f32>,
    bl: vec3<f32>,
    br: vec3<f32>,
};

fn sheet_corners(p: vec3<f32>, n: vec3<f32>, b: vec3<f32>, hw: f32, half_thick: f32) -> SheetCorners {
    var c: SheetCorners;
    c.tl = p + n * hw + b * half_thick;
    c.tr = p - n * hw + b * half_thick;
    c.bl = p + n * hw - b * half_thick;
    c.br = p - n * hw - b * half_thick;
    return c;
}

// Emit one vertex of a face strip (4 face strips × 6 verts/quad = 24 verts/pair).
fn emit_sheet_body_vertex(
    run: RunDescriptor,
    pair_idx: u32,
    in_pair: u32,         // 0..23
    out_idx: u32,
) {
    let s0 = run.sample_start + pair_idx;
    let s1 = s0 + 1u;
    let f0 = frame_at(s0, run.sample_start, run.sample_end);
    let f1 = frame_at(s1, run.sample_start, run.sample_end);
    let half_width = SHEET_HEIGHT;        // wide half-extent along NORMAL
    let half_thick = SHEET_WIDTH * 0.5;   // thin half-extent along BINORMAL
    let c0 = sheet_corners(f0.position, f0.normal, f0.binormal, half_width, half_thick);
    let c1 = sheet_corners(f1.position, f1.normal, f1.binormal, half_width, half_thick);
    let id0 = extrude_points[s0].atom_idx;
    let id1 = extrude_points[s1].atom_idx;

    // 4 face strips × 6 verts:
    //   0..5  TOP    (face normal = +binormal): tl0, tr0, tr1, tl0, tr1, tl1
    //   6..11 BOTTOM (face normal = -binormal): br0, bl0, bl1, br0, bl1, br1
    //   12..17 LEFT  (face normal = +normal):   tl0, tl1, bl1, tl0, bl1, bl0
    //   18..23 RIGHT (face normal = -normal):   tr1, tr0, br0, tr1, br0, br1
    let face = in_pair / 6u;
    let tri_v = in_pair - face * 6u;
    var pos: vec3<f32>;
    var nrm: vec3<f32>;
    var atom: u32;
    switch face {
        case 0u: { // TOP +binormal
            switch tri_v {
                case 0u, 3u: { pos = c0.tl; nrm = f0.binormal; atom = id0; }
                case 1u: { pos = c0.tr; nrm = f0.binormal; atom = id0; }
                case 2u, 4u: { pos = c1.tr; nrm = f1.binormal; atom = id1; }
                case 5u: { pos = c1.tl; nrm = f1.binormal; atom = id1; }
                default: { pos = c0.tl; nrm = f0.binormal; atom = id0; }
            }
        }
        case 1u: { // BOTTOM -binormal
            switch tri_v {
                case 0u, 3u: { pos = c0.br; nrm = -f0.binormal; atom = id0; }
                case 1u: { pos = c0.bl; nrm = -f0.binormal; atom = id0; }
                case 2u, 4u: { pos = c1.bl; nrm = -f1.binormal; atom = id1; }
                case 5u: { pos = c1.br; nrm = -f1.binormal; atom = id1; }
                default: { pos = c0.br; nrm = -f0.binormal; atom = id0; }
            }
        }
        case 2u: { // LEFT +normal
            switch tri_v {
                case 0u, 3u: { pos = c0.tl; nrm = f0.normal; atom = id0; }
                case 1u: { pos = c1.tl; nrm = f1.normal; atom = id1; }
                case 2u, 4u: { pos = c1.bl; nrm = f1.normal; atom = id1; }
                case 5u: { pos = c0.bl; nrm = f0.normal; atom = id0; }
                default: { pos = c0.tl; nrm = f0.normal; atom = id0; }
            }
        }
        default: { // RIGHT -normal
            switch tri_v {
                case 0u, 3u: { pos = c1.tr; nrm = -f1.normal; atom = id1; }
                case 1u: { pos = c0.tr; nrm = -f0.normal; atom = id0; }
                case 2u, 4u: { pos = c0.br; nrm = -f0.normal; atom = id0; }
                case 5u: { pos = c1.br; nrm = -f1.normal; atom = id1; }
                default: { pos = c1.tr; nrm = -f1.normal; atom = id1; }
            }
        }
    }
    // Sheet faces are open ribbons — two-sided lighting flips back-face normals.
    vertices[out_idx] = make_vertex(pos, nrm, atom, TWO_SIDED_FLAG);
}

// Sheet back / forward cap (4 vertices = 1 quad). face_normal: +/- tangent.
fn emit_sheet_cap_vertex(
    run: RunDescriptor,
    is_back: bool,
    tri_v: u32,
    out_idx: u32,
) {
    let s = select(run.sample_end, run.sample_start, is_back);
    let f = frame_at(s, run.sample_start, run.sample_end);
    let half_width = SHEET_HEIGHT;
    let half_thick = SHEET_WIDTH * 0.5;
    let c = sheet_corners(f.position, f.normal, f.binormal, half_width, half_thick);
    let cap_n = select(f.tangent, -f.tangent, is_back);
    let id = extrude_points[s].atom_idx;
    // Quad order for back cap (face normal -tangent): tl → bl → br → tl → br → tr.
    // For forward cap (face normal +tangent):         tl → tr → br → tl → br → bl.
    var pos: vec3<f32>;
    if (is_back) {
        switch tri_v {
            case 0u, 3u: { pos = c.tl; }
            case 1u: { pos = c.bl; }
            case 2u, 4u: { pos = c.br; }
            case 5u: { pos = c.tr; }
            default: { pos = c.tl; }
        }
    } else {
        switch tri_v {
            case 0u, 3u: { pos = c.tl; }
            case 1u: { pos = c.tr; }
            case 2u, 4u: { pos = c.br; }
            case 5u: { pos = c.bl; }
            default: { pos = c.tl; }
        }
    }
    vertices[out_idx] = make_vertex(pos, cap_n, id, TWO_SIDED_FLAG);
}

// ============================================================================
// Main dispatcher: linear scan over runs[] to find which run owns this vertex.
// For typical structures n_runs is < 100, so linear is fine.
// ============================================================================

fn find_run_for_vertex(vidx: u32) -> i32 {
    for (var r = 0u; r < params.n_runs; r = r + 1u) {
        let run = runs[r];
        if (vidx >= run.vertex_offset && vidx < run.vertex_offset + run.vertex_count) {
            return i32(r);
        }
    }
    return -1;
}

@compute @workgroup_size(64)
fn cs_main(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(num_workgroups) nwg: vec3<u32>,
) {
    // 2D dispatch handles total_verts > 65535*64.
    let out_idx = gid.x + gid.y * nwg.x * 64u;
    let total_verts = arrayLength(&vertices);
    if (out_idx >= total_verts) {
        return;
    }
    let r = find_run_for_vertex(out_idx);
    if (r < 0) {
        // Beyond emitted runs — emit a degenerate vertex.
        vertices[out_idx] = make_degen();
        return;
    }
    let run = runs[u32(r)];
    let local_vidx = out_idx - run.vertex_offset;

    switch run.car_type {
        case 0u: { emit_tube_vertex(run, local_vidx, false, out_idx); } // Loop
        case 1u: { emit_tube_vertex(run, local_vidx, true, out_idx); }  // Oval
        case 2u, 4u: {
            // Rect / NucleicRect (CartoonType::Rect=2 / NucleicRect=4):
            // body + back cap + forward cap.
            let pair_count = run.sample_end - run.sample_start; // (n-1) pairs
            let body_total = pair_count * 4u * 6u;
            if (local_vidx < body_total) {
                let pair_idx = local_vidx / 24u;
                let in_pair = local_vidx - pair_idx * 24u;
                emit_sheet_body_vertex(run, pair_idx, in_pair, out_idx);
            } else if (local_vidx < body_total + 6u) {
                // Back cap (at sample_start).
                emit_sheet_cap_vertex(run, true, local_vidx - body_total, out_idx);
            } else {
                // Forward cap (at sample_end).
                emit_sheet_cap_vertex(run, false, local_vidx - body_total - 6u, out_idx);
            }
        }
        case 3u: {
            // Arrow: body + arrow + back cap + shoulder caps.
            // body_pairs = run.body_end (in pairs); sample_count - 1 - body_pairs = arrow_pairs.
            let pair_count = run.sample_end - run.sample_start;
            let body_pairs = run.body_end;
            let arrow_pairs = pair_count - body_pairs;
            let body_total = body_pairs * 4u * 6u;
            let arrow_total = arrow_pairs * 4u * 6u;
            let n_cap_total = 6u;
            let shoulder_total = 12u;
            if (local_vidx < body_total) {
                let pair_idx = local_vidx / 24u;
                let in_pair = local_vidx - pair_idx * 24u;
                emit_sheet_body_vertex(run, pair_idx, in_pair, out_idx);
            } else if (local_vidx < body_total + arrow_total) {
                // Arrow segments — emit with FIXED base frame.
                let arrow_local = local_vidx - body_total;
                let seg_idx = arrow_local / 24u;
                let in_seg = arrow_local - seg_idx * 24u;
                emit_arrow_segment_vertex(run, seg_idx, arrow_pairs, in_seg, out_idx);
            } else if (local_vidx < body_total + arrow_total + n_cap_total) {
                emit_sheet_cap_vertex(run, true, local_vidx - body_total - arrow_total, out_idx);
            } else {
                // Shoulder caps (2 quads = 12 verts).
                let sh_local = local_vidx - body_total - arrow_total - n_cap_total;
                emit_shoulder_cap_vertex(run, sh_local, out_idx);
            }
        }
        default: {
            // Unknown CartoonType — fall back to tube emission.
            emit_tube_vertex(run, local_vidx, true, out_idx);
        }
    }
}

fn emit_arrow_segment_vertex(
    run: RunDescriptor,
    seg_idx: u32,
    arrow_pairs: u32,
    in_seg: u32,        // 0..23
    out_idx: u32,
) {
    // Base = first arrow sample (= run.sample_start + run.body_end).
    let base_s = run.sample_start + run.body_end;
    let tip_s = run.sample_end;
    let base_f = frame_at(base_s, run.sample_start, run.sample_end);
    let n = base_f.normal;
    let b = base_f.binormal;
    let half_thick = SHEET_WIDTH * 0.5;
    let arrow_width = SHEET_HEIGHT * ARROW_TIP_SCALE;
    let start_pos = base_f.position;
    let end_pos = extrude_points[tip_s].position;
    let t0 = f32(seg_idx) / f32(arrow_pairs);
    let t1 = f32(seg_idx + 1u) / f32(arrow_pairs);
    let w0 = arrow_width * (1.0 - t0);
    let w1 = arrow_width * (1.0 - t1);
    let pos0 = start_pos + (end_pos - start_pos) * t0;
    let pos1 = start_pos + (end_pos - start_pos) * t1;
    let tl0 = pos0 + n * w0 + b * half_thick;
    let tr0 = pos0 - n * w0 + b * half_thick;
    let bl0 = pos0 + n * w0 - b * half_thick;
    let br0 = pos0 - n * w0 - b * half_thick;
    let tl1 = pos1 + n * w1 + b * half_thick;
    let tr1 = pos1 - n * w1 + b * half_thick;
    let bl1 = pos1 + n * w1 - b * half_thick;
    let br1 = pos1 - n * w1 - b * half_thick;
    let id0 = extrude_points[base_s + seg_idx].atom_idx;
    let id1 = extrude_points[base_s + seg_idx + 1u].atom_idx;
    let face = in_seg / 6u;
    let tri_v = in_seg - face * 6u;
    var pos: vec3<f32>;
    var nrm: vec3<f32>;
    var atom: u32;
    switch face {
        case 0u: { // TOP +binormal
            switch tri_v {
                case 0u, 3u: { pos = tl0; nrm = b; atom = id0; }
                case 1u: { pos = tr0; nrm = b; atom = id0; }
                case 2u, 4u: { pos = tr1; nrm = b; atom = id1; }
                case 5u: { pos = tl1; nrm = b; atom = id1; }
                default: { pos = tl0; nrm = b; atom = id0; }
            }
        }
        case 1u: { // BOTTOM -binormal
            switch tri_v {
                case 0u, 3u: { pos = br0; nrm = -b; atom = id0; }
                case 1u: { pos = bl0; nrm = -b; atom = id0; }
                case 2u, 4u: { pos = bl1; nrm = -b; atom = id1; }
                case 5u: { pos = br1; nrm = -b; atom = id1; }
                default: { pos = br0; nrm = -b; atom = id0; }
            }
        }
        case 2u: { // LEFT +normal
            switch tri_v {
                case 0u, 3u: { pos = tl0; nrm = n; atom = id0; }
                case 1u: { pos = tl1; nrm = n; atom = id1; }
                case 2u, 4u: { pos = bl1; nrm = n; atom = id1; }
                case 5u: { pos = bl0; nrm = n; atom = id0; }
                default: { pos = tl0; nrm = n; atom = id0; }
            }
        }
        default: { // RIGHT -normal
            switch tri_v {
                case 0u, 3u: { pos = tr1; nrm = -n; atom = id1; }
                case 1u: { pos = tr0; nrm = -n; atom = id0; }
                case 2u, 4u: { pos = br0; nrm = -n; atom = id0; }
                case 5u: { pos = br1; nrm = -n; atom = id1; }
                default: { pos = tr1; nrm = -n; atom = id1; }
            }
        }
    }
    vertices[out_idx] = make_vertex(pos, nrm, atom, TWO_SIDED_FLAG);
}

fn emit_shoulder_cap_vertex(
    run: RunDescriptor,
    sh_local: u32,       // 0..11
    out_idx: u32,
) {
    let base_s = run.sample_start + run.body_end;
    let f = frame_at(base_s, run.sample_start, run.sample_end);
    let n = f.normal;
    let b = f.binormal;
    let half_width = SHEET_HEIGHT;
    let half_thick = SHEET_WIDTH * 0.5;
    let arrow_width = half_width * ARROW_TIP_SCALE;
    let p = f.position;
    let body_tl = p + n * half_width + b * half_thick;
    let body_tr = p - n * half_width + b * half_thick;
    let body_bl = p + n * half_width - b * half_thick;
    let body_br = p - n * half_width - b * half_thick;
    let arrow_tl = p + n * arrow_width + b * half_thick;
    let arrow_tr = p - n * arrow_width + b * half_thick;
    let arrow_bl = p + n * arrow_width - b * half_thick;
    let arrow_br = p - n * arrow_width - b * half_thick;
    let cap_n = -f.tangent;
    let id = extrude_points[base_s].atom_idx;
    // Left shoulder (sh_local 0..5): body_tl → arrow_tl → arrow_bl → body_tl → arrow_bl → body_bl
    // Right shoulder (sh_local 6..11): arrow_tr → body_tr → body_br → arrow_tr → body_br → arrow_br
    var pos: vec3<f32>;
    if (sh_local < 6u) {
        switch sh_local {
            case 0u, 3u: { pos = body_tl; }
            case 1u: { pos = arrow_tl; }
            case 2u, 4u: { pos = arrow_bl; }
            case 5u: { pos = body_bl; }
            default: { pos = body_tl; }
        }
    } else {
        let local = sh_local - 6u;
        switch local {
            case 0u, 3u: { pos = arrow_tr; }
            case 1u: { pos = body_tr; }
            case 2u, 4u: { pos = body_br; }
            case 5u: { pos = arrow_br; }
            default: { pos = arrow_tr; }
        }
    }
    vertices[out_idx] = make_vertex(pos, cap_n, id, TWO_SIDED_FLAG);
}
