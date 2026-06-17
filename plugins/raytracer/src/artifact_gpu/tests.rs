use super::indirect::{active_triangle_vertex_count, decode_draw_indirect_args};
use super::layout::{
    ArtifactTriangle, ArtifactVisibleTriangleParams, ATOM_STRIDE, COLOR_LUT_STRIDE,
    EMPTY_STORAGE_BYTES, LINE_INSTANCE_STRIDE, RAY_LINE_RADIUS, SPHERE_INSTANCE_STRIDE,
    STD_VERTEX_STRIDE, STICK_INSTANCE_STRIDE, WORKGROUP_SIZE,
};
use super::resources::{checked_storage_buffer_size, storage_bytes_for_device};
use super::*;
use crate::primitive::GpuTriangle;
use crate::shaders;
use patinae_scene::{
    GpuHandle, GpuHandleKind, RenderArtifactBufferDescriptor, RenderArtifactBufferRole,
    RenderArtifactPrimitiveTopology, RenderArtifactRepDescriptor, RenderArtifactRepKind,
    RenderArtifactSnapshotDescriptor,
};

fn handle(id: u64) -> GpuHandle {
    GpuHandle {
        id,
        kind: GpuHandleKind::Buffer,
        generation: 1,
    }
}

fn scene_atoms_descriptor(id: u64, element_count: u64) -> RenderArtifactBufferDescriptor {
    RenderArtifactBufferDescriptor {
        handle: handle(id),
        role: RenderArtifactBufferRole::SceneAtoms,
        size: element_count * ATOM_STRIDE,
        stride: ATOM_STRIDE,
        element_count,
    }
}

fn assert_shader_has(shader: &str, needle: &str) {
    assert!(
        shader.contains(needle),
        "expected shader to contain `{needle}`"
    );
}

fn assert_shader_lacks(shader: &str, needle: &str) {
    assert!(
        !shader.contains(needle),
        "expected shader not to contain `{needle}`"
    );
}

// Dispatch, layout, and storage sizing.

#[test]
fn bvh_shape_pads_leaves_to_power_of_two() {
    let shape = bvh::bvh_shape_for(9).expect("shape");

    assert_eq!(shape.leaf_slots, 4);
    assert_eq!(shape.leaf_start, 3);
    assert_eq!(shape.node_count, 7);
}

#[test]
fn dispatch_grid_tiles_past_single_dimension_limit() {
    let items = 90_473 * WORKGROUP_SIZE;
    let grid = dispatch::dispatch_grid_for(items, 65_535).expect("dispatch grid");

    assert_eq!(grid.workgroups, [65_535, 2, 1]);
    assert_eq!(grid.invocation_width, 65_535 * WORKGROUP_SIZE);
}

#[test]
fn typed_storage_buffers_keep_one_element_for_empty_arrays() {
    assert_eq!(
        resources::storage_bytes_for::<GpuSphere>(0, "spheres").expect("sphere bytes"),
        std::mem::size_of::<GpuSphere>() as u64
    );
    assert_eq!(
        resources::storage_bytes_for::<GpuCylinder>(0, "cylinders").expect("cylinder bytes"),
        std::mem::size_of::<GpuCylinder>() as u64
    );
    assert_eq!(
        resources::storage_bytes_for::<GpuCapsule>(0, "capsules").expect("capsule bytes"),
        std::mem::size_of::<GpuCapsule>() as u64
    );
    assert_eq!(
        resources::storage_bytes_for::<GpuTriangle>(0, "triangles").expect("triangle bytes"),
        std::mem::size_of::<GpuTriangle>() as u64
    );
    assert_eq!(
        resources::storage_bytes_for::<ArtifactTriangle>(0, "artifact triangles")
            .expect("artifact triangle bytes"),
        std::mem::size_of::<ArtifactTriangle>() as u64
    );
}

#[test]
fn artifact_triangle_layout_is_compact() {
    assert_eq!(std::mem::size_of::<ArtifactTriangle>(), 80);
}

#[test]
fn artifact_triangle_storage_fits_reported_3j3q_visible_count() {
    let limits = patinae_scene::GpuDeviceLimits {
        max_buffer_size: 4_294_967_296,
        max_storage_buffer_binding_size: 2_147_483_647,
        max_compute_workgroups_per_dimension: 65_535,
        max_compute_invocations_per_workgroup: 256,
        max_compute_workgroup_size_x: 256,
        max_compute_workgroup_size_y: 1,
        max_compute_workgroup_size_z: 1,
    };
    let visible_triangles = 25_174_294;

    let compact_bytes = storage_bytes_for_device::<ArtifactTriangle>(
        visible_triangles,
        &limits,
        "ray artifact triangles",
    )
    .expect("compact artifact triangle storage fits");

    assert_eq!(compact_bytes, 2_013_943_520);
    assert!(storage_bytes_for_device::<GpuTriangle>(
        visible_triangles,
        &limits,
        "ray artifact triangles",
    )
    .unwrap_err()
    .contains(
        "ray artifact triangles buffer size 3222309632 exceeds GPU storage buffer limit 2147483647"
    ));
}

#[test]
fn indirect_triangle_active_vertex_count_uses_draw_count_not_capacity() {
    assert_eq!(active_triangle_vertex_count(300, 4_000_000), 300);
    assert_eq!(active_triangle_vertex_count(302, 4_000_000), 300);
    assert_eq!(
        active_triangle_vertex_count(4_000_000, 4_000_000),
        3_999_999
    );
    assert_eq!(active_triangle_vertex_count(9, 5), 3);
}

#[test]
fn draw_indirect_args_decode_as_little_endian_u32_lanes() {
    let bytes = [0x78, 0x56, 0x34, 0x12, 2, 0, 0, 0, 3, 0, 0, 0, 4, 0, 0, 0];

    let args =
        decode_draw_indirect_args(&bytes, RenderArtifactRepKind::Surface).expect("draw args");

    assert_eq!(args, [0x1234_5678, 2, 3, 4]);
}

#[test]
fn storage_buffer_size_check_uses_storage_binding_limit() {
    let limits = patinae_scene::GpuDeviceLimits {
        max_buffer_size: 1024,
        max_storage_buffer_binding_size: 512,
        max_compute_workgroups_per_dimension: 65_535,
        max_compute_invocations_per_workgroup: 256,
        max_compute_workgroup_size_x: 256,
        max_compute_workgroup_size_y: 1,
        max_compute_workgroup_size_z: 1,
    };

    assert!(checked_storage_buffer_size(512, &limits, "test storage").is_ok());
    let err = checked_storage_buffer_size(513, &limits, "test storage").unwrap_err();

    assert!(err.contains("test storage buffer size 513 exceeds GPU storage buffer limit 512"));
}

// Visibility and indirect draw helpers.

#[test]
fn visible_triangle_params_uniform_layout_matches_wgsl() {
    assert_eq!(std::mem::size_of::<ArtifactVisibleTriangleParams>(), 160);
}

#[test]
fn visible_triangle_counts_decode_as_little_endian_u32_lanes() {
    let bytes = [7, 0, 0, 0, 0x78, 0x56, 0x34, 0x12];

    let counts = visibility::decode_visible_triangle_counts(&bytes, 2).expect("counts");

    assert_eq!(counts, vec![7, 0x1234_5678]);
    assert!(visibility::decode_visible_triangle_counts(&bytes[..4], 2)
        .unwrap_err()
        .contains("returned 4 bytes, expected 8"));
}

// WGSL expansion and validation.

#[test]
fn raytrace_output_buffer_shader_uses_output_buffer_and_empty_node_guard() {
    let shader = shaders::artifact_raytrace_output_buffer();

    assert_shader_has(&shader, "output_pixels");
    assert_shader_has(
        &shader,
        "global_id.y * u32(uniforms.viewport.x) + global_id.x",
    );
    assert_shader_has(&shader, "node.left_or_first == 0xffffffffu");
    assert_shader_lacks(&shader, "texture_storage_2d");
    assert_shader_lacks(&shader, "textureStore(");
}

#[test]
fn artifact_triangle_shader_uses_resolved_vertex_count_without_indirect_draw_args() {
    let shader = shaders::artifact_triangles();

    assert_shader_lacks(&shader, "draw_args");
    assert_shader_has(&shader, "params.vertex_capacity / 3u");
    assert_shader_has(&shader, "scene_atoms");
    assert_shader_has(&shader, "atom.alpha_pack_b");
    assert_shader_has(&shader, "normal_oct: vec4<u32>");
}

#[test]
fn artifact_raytrace_shader_decodes_compact_triangle_normals() {
    let shader = shaders::artifact_raytrace_output_buffer();

    assert_shader_has(&shader, "normal_oct: vec4<u32>");
    assert_shader_has(&shader, "oct_decode(tri.normal_oct.x)");
    assert_shader_has(&shader, "hit.transparency = 1.0 - tri.color.a");
    assert_shader_has(&shader, "trace_ray(current_ray, opaque_only_trace)");
    assert_shader_has(&shader, "opaque_only_trace = true");
}

#[test]
fn artifact_wgsl_modules_parse_and_validate() {
    let modules: [(&str, String); 10] = [
        ("ray.main", shaders::RAYTRACE.to_string()),
        ("ray.artifact.spheres", shaders::artifact_spheres()),
        (
            "ray.artifact.stick_capsules",
            shaders::artifact_stick_capsules(),
        ),
        (
            "ray.artifact.line_cylinders",
            shaders::artifact_line_cylinders(),
        ),
        ("ray.artifact.triangles", shaders::artifact_triangles()),
        (
            "ray.artifact.count_visible_triangles",
            shaders::artifact_count_visible_triangles(),
        ),
        (
            "ray.artifact.visible_triangles",
            shaders::artifact_visible_triangles(),
        ),
        ("ray.artifact.bvh", shaders::artifact_bvh()),
        (
            "ray.artifact.raytrace",
            shaders::artifact_raytrace_output_buffer(),
        ),
        (
            "ray.artifact.downsample",
            shaders::ARTIFACT_DOWNSAMPLE.to_string(),
        ),
    ];

    let mut failures = Vec::new();

    for (label, source) in modules {
        let module = match naga::front::wgsl::parse_str(&source) {
            Ok(module) => module,
            Err(err) => {
                failures.push(format!(
                    "{label}: WGSL parse failed\n{}",
                    err.emit_to_string_with_path(&source, label)
                ));
                continue;
            }
        };

        let mut validator = naga::valid::Validator::new(
            naga::valid::ValidationFlags::all(),
            naga::valid::Capabilities::all(),
        );
        if let Err(err) = validator.validate(&module) {
            failures.push(format!("{label}: WGSL validation failed\n{err:#?}"));
        }
    }

    assert!(
        failures.is_empty(),
        "invalid ray artifact WGSL modules:\n{}",
        failures.join("\n\n")
    );
}

// Artifact planning.

#[test]
fn artifact_plan_uses_scene_color_lut_and_cartoon_slot() {
    let snapshot = RenderArtifactSnapshotDescriptor {
        snapshot_id: 1,
        layout_version: 2,
        scene_generation: 7,
        scene_bounds_min: [0.0, 0.0, 0.0],
        scene_bounds_max: [1.0, 1.0, 1.0],
        cull_pass_initialized: true,
        device_limits: patinae_scene::GpuDeviceLimits {
            max_buffer_size: 1024,
            max_storage_buffer_binding_size: 1024,
            max_compute_workgroups_per_dimension: 65_535,
            max_compute_invocations_per_workgroup: 256,
            max_compute_workgroup_size_x: 256,
            max_compute_workgroup_size_y: 1,
            max_compute_workgroup_size_z: 1,
        },
        buffers: vec![
            RenderArtifactBufferDescriptor {
                handle: handle(10),
                role: RenderArtifactBufferRole::SceneColorLut,
                size: 64,
                stride: COLOR_LUT_STRIDE,
                element_count: 1,
            },
            scene_atoms_descriptor(90, 8),
            RenderArtifactBufferDescriptor {
                handle: handle(11),
                role: RenderArtifactBufferRole::StdVertices,
                size: 72,
                stride: STD_VERTEX_STRIDE,
                element_count: 3,
            },
        ],
        reps: vec![RenderArtifactRepDescriptor {
            object_id: 1,
            rep_kind: RenderArtifactRepKind::Cartoon,
            topology: RenderArtifactPrimitiveTopology::TriangleList,
            geometry: handle(11),
            count: None,
            indirect: None,
            element_count: 3,
            max_element_count: 3,
            atom_offset: 5,
            atom_count: 1,
            material_rgba: [1.0, 0.0, 0.0, 1.0],
            transparency: 0.0,
        }],
    };

    let plan = plan::plan_artifact_primitives(&snapshot).expect("plan");

    assert_eq!(plan.color_lut, handle(10));
    assert_eq!(plan.scene_atoms, Some(handle(90)));
    assert_eq!(plan.triangle_count, 1);
    assert_eq!(plan.triangle_reps[0].rep_slot, 4);
    assert_eq!(plan.primitive_count().expect("primitive count"), 1);
}

#[test]
fn artifact_plan_accepts_native_instance_artifacts() {
    let snapshot = RenderArtifactSnapshotDescriptor {
        snapshot_id: 1,
        layout_version: 2,
        scene_generation: 7,
        scene_bounds_min: [0.0, 0.0, 0.0],
        scene_bounds_max: [1.0, 1.0, 1.0],
        cull_pass_initialized: true,
        device_limits: patinae_scene::GpuDeviceLimits {
            max_buffer_size: 4096,
            max_storage_buffer_binding_size: 4096,
            max_compute_workgroups_per_dimension: 65_535,
            max_compute_invocations_per_workgroup: 256,
            max_compute_workgroup_size_x: 256,
            max_compute_workgroup_size_y: 1,
            max_compute_workgroup_size_z: 1,
        },
        buffers: vec![
            RenderArtifactBufferDescriptor {
                handle: handle(10),
                role: RenderArtifactBufferRole::SceneColorLut,
                size: 64,
                stride: COLOR_LUT_STRIDE,
                element_count: 1,
            },
            RenderArtifactBufferDescriptor {
                handle: handle(11),
                role: RenderArtifactBufferRole::SphereInstances,
                size: 5 * SPHERE_INSTANCE_STRIDE,
                stride: SPHERE_INSTANCE_STRIDE,
                element_count: 5,
            },
            RenderArtifactBufferDescriptor {
                handle: handle(21),
                role: RenderArtifactBufferRole::StickInstances,
                size: 7 * STICK_INSTANCE_STRIDE,
                stride: STICK_INSTANCE_STRIDE,
                element_count: 7,
            },
            RenderArtifactBufferDescriptor {
                handle: handle(31),
                role: RenderArtifactBufferRole::LineInstances,
                size: 11 * LINE_INSTANCE_STRIDE,
                stride: LINE_INSTANCE_STRIDE,
                element_count: 11,
            },
        ],
        reps: vec![
            RenderArtifactRepDescriptor {
                object_id: 1,
                rep_kind: RenderArtifactRepKind::Sphere,
                topology: RenderArtifactPrimitiveTopology::SphereInstances,
                geometry: handle(11),
                count: Some(handle(12)),
                indirect: Some(handle(13)),
                element_count: 0,
                max_element_count: 5,
                atom_offset: 0,
                atom_count: 5,
                material_rgba: [1.0, 0.0, 0.0, 1.0],
                transparency: 0.0,
            },
            RenderArtifactRepDescriptor {
                object_id: 1,
                rep_kind: RenderArtifactRepKind::Stick,
                topology: RenderArtifactPrimitiveTopology::CylinderInstances,
                geometry: handle(21),
                count: Some(handle(22)),
                indirect: Some(handle(23)),
                element_count: 0,
                max_element_count: 7,
                atom_offset: 0,
                atom_count: 5,
                material_rgba: [0.0, 1.0, 0.0, 1.0],
                transparency: 0.25,
            },
            RenderArtifactRepDescriptor {
                object_id: 1,
                rep_kind: RenderArtifactRepKind::Line,
                topology: RenderArtifactPrimitiveTopology::LineInstances,
                geometry: handle(31),
                count: Some(handle(32)),
                indirect: Some(handle(33)),
                element_count: 0,
                max_element_count: 11,
                atom_offset: 0,
                atom_count: 5,
                material_rgba: [0.0, 0.0, 1.0, 1.0],
                transparency: 0.0,
            },
        ],
    };

    let plan = plan::plan_artifact_primitives(&snapshot).expect("plan");

    assert_eq!(plan.sphere_count, 5);
    assert_eq!(plan.cylinder_count, 11);
    assert_eq!(plan.capsule_count, 7);
    assert_eq!(plan.triangle_count, 0);
    assert_eq!(plan.sphere_reps[0].rep_slot, 0);
    assert_eq!(
        plan.sphere_reps[0].geometry_binding_size,
        5 * SPHERE_INSTANCE_STRIDE
    );
    assert_eq!(plan.capsule_reps[0].rep_slot, 1);
    assert_eq!(
        plan.capsule_reps[0].geometry_binding_size,
        7 * STICK_INSTANCE_STRIDE
    );
    assert_eq!(plan.cylinder_reps[0].rep_slot, 2);
    assert_eq!(plan.cylinder_reps[0].radius, RAY_LINE_RADIUS);
    assert_eq!(
        plan.cylinder_reps[0].geometry_binding_size,
        11 * LINE_INSTANCE_STRIDE
    );
    assert_eq!(plan.primitive_count().expect("primitive count"), 23);
}

#[test]
fn artifact_plan_skips_undersized_instance_geometry() {
    let snapshot = RenderArtifactSnapshotDescriptor {
        snapshot_id: 1,
        layout_version: 2,
        scene_generation: 7,
        scene_bounds_min: [0.0, 0.0, 0.0],
        scene_bounds_max: [1.0, 1.0, 1.0],
        cull_pass_initialized: true,
        device_limits: patinae_scene::GpuDeviceLimits {
            max_buffer_size: 4096,
            max_storage_buffer_binding_size: 4096,
            max_compute_workgroups_per_dimension: 65_535,
            max_compute_invocations_per_workgroup: 256,
            max_compute_workgroup_size_x: 256,
            max_compute_workgroup_size_y: 1,
            max_compute_workgroup_size_z: 1,
        },
        buffers: vec![
            RenderArtifactBufferDescriptor {
                handle: handle(10),
                role: RenderArtifactBufferRole::SceneColorLut,
                size: 64,
                stride: COLOR_LUT_STRIDE,
                element_count: 1,
            },
            RenderArtifactBufferDescriptor {
                handle: handle(21),
                role: RenderArtifactBufferRole::StickInstances,
                size: EMPTY_STORAGE_BYTES,
                stride: STICK_INSTANCE_STRIDE,
                element_count: 1,
            },
        ],
        reps: vec![RenderArtifactRepDescriptor {
            object_id: 1,
            rep_kind: RenderArtifactRepKind::Stick,
            topology: RenderArtifactPrimitiveTopology::CylinderInstances,
            geometry: handle(21),
            count: Some(handle(22)),
            indirect: Some(handle(23)),
            element_count: 0,
            max_element_count: 1,
            atom_offset: 0,
            atom_count: 2,
            material_rgba: [0.0, 1.0, 0.0, 1.0],
            transparency: 0.0,
        }],
    };

    let plan = plan::plan_artifact_primitives(&snapshot).expect("plan");

    assert_eq!(plan.cylinder_count, 0);
    assert_eq!(plan.capsule_count, 0);
    assert!(plan.cylinder_reps.is_empty());
    assert!(plan.capsule_reps.is_empty());
    assert_eq!(plan.primitive_count().expect("primitive count"), 0);
}

#[test]
fn artifact_plan_accepts_indirect_surface_capacity() {
    let surface_capacity = 4_000_000;
    let snapshot = RenderArtifactSnapshotDescriptor {
        snapshot_id: 1,
        layout_version: 2,
        scene_generation: 7,
        scene_bounds_min: [0.0, 0.0, 0.0],
        scene_bounds_max: [1.0, 1.0, 1.0],
        cull_pass_initialized: true,
        device_limits: patinae_scene::GpuDeviceLimits {
            max_buffer_size: 128 * 1024 * 1024,
            max_storage_buffer_binding_size: 128 * 1024 * 1024,
            max_compute_workgroups_per_dimension: 65_535,
            max_compute_invocations_per_workgroup: 256,
            max_compute_workgroup_size_x: 256,
            max_compute_workgroup_size_y: 1,
            max_compute_workgroup_size_z: 1,
        },
        buffers: vec![
            RenderArtifactBufferDescriptor {
                handle: handle(10),
                role: RenderArtifactBufferRole::SceneColorLut,
                size: 64,
                stride: COLOR_LUT_STRIDE,
                element_count: 1,
            },
            scene_atoms_descriptor(90, 16),
            RenderArtifactBufferDescriptor {
                handle: handle(11),
                role: RenderArtifactBufferRole::StdVertices,
                size: surface_capacity * STD_VERTEX_STRIDE,
                stride: STD_VERTEX_STRIDE,
                element_count: surface_capacity,
            },
        ],
        reps: vec![RenderArtifactRepDescriptor {
            object_id: 1,
            rep_kind: RenderArtifactRepKind::Surface,
            topology: RenderArtifactPrimitiveTopology::TriangleList,
            geometry: handle(11),
            count: None,
            indirect: Some(handle(12)),
            element_count: 0,
            max_element_count: surface_capacity,
            atom_offset: 5,
            atom_count: 8,
            material_rgba: [0.5, 0.5, 0.5, 1.0],
            transparency: 0.5,
        }],
    };

    let plan = plan::plan_artifact_primitives(&snapshot).expect("plan");

    assert_eq!(plan.triangle_count, 1_333_333);
    assert_eq!(plan.triangle_reps[0].rep_slot, 6);
    assert_eq!(plan.triangle_reps[0].vertex_count, 3_999_999);
    assert_eq!(plan.triangle_reps[0].visibility_counter_index, None);
    assert_eq!(
        plan.triangle_reps[0].geometry_binding_size,
        3_999_999 * STD_VERTEX_STRIDE
    );
}

// Surface visibility planning.

#[test]
fn surface_visibility_counts_compact_triangle_offsets() {
    let snapshot = RenderArtifactSnapshotDescriptor {
        snapshot_id: 1,
        layout_version: 2,
        scene_generation: 7,
        scene_bounds_min: [0.0, 0.0, 0.0],
        scene_bounds_max: [1.0, 1.0, 1.0],
        cull_pass_initialized: true,
        device_limits: patinae_scene::GpuDeviceLimits {
            max_buffer_size: 4096,
            max_storage_buffer_binding_size: 4096,
            max_compute_workgroups_per_dimension: 65_535,
            max_compute_invocations_per_workgroup: 256,
            max_compute_workgroup_size_x: 256,
            max_compute_workgroup_size_y: 1,
            max_compute_workgroup_size_z: 1,
        },
        buffers: vec![
            RenderArtifactBufferDescriptor {
                handle: handle(10),
                role: RenderArtifactBufferRole::SceneColorLut,
                size: 64,
                stride: COLOR_LUT_STRIDE,
                element_count: 1,
            },
            scene_atoms_descriptor(90, 16),
            RenderArtifactBufferDescriptor {
                handle: handle(11),
                role: RenderArtifactBufferRole::StdVertices,
                size: 6 * STD_VERTEX_STRIDE,
                stride: STD_VERTEX_STRIDE,
                element_count: 6,
            },
            RenderArtifactBufferDescriptor {
                handle: handle(12),
                role: RenderArtifactBufferRole::StdVertices,
                size: 9 * STD_VERTEX_STRIDE,
                stride: STD_VERTEX_STRIDE,
                element_count: 9,
            },
        ],
        reps: vec![
            RenderArtifactRepDescriptor {
                object_id: 1,
                rep_kind: RenderArtifactRepKind::Cartoon,
                topology: RenderArtifactPrimitiveTopology::TriangleList,
                geometry: handle(11),
                count: None,
                indirect: None,
                element_count: 6,
                max_element_count: 6,
                atom_offset: 5,
                atom_count: 1,
                material_rgba: [1.0, 0.0, 0.0, 1.0],
                transparency: 0.0,
            },
            RenderArtifactRepDescriptor {
                object_id: 1,
                rep_kind: RenderArtifactRepKind::Surface,
                topology: RenderArtifactPrimitiveTopology::TriangleList,
                geometry: handle(12),
                count: None,
                indirect: None,
                element_count: 9,
                max_element_count: 9,
                atom_offset: 6,
                atom_count: 1,
                material_rgba: [0.5, 0.5, 0.5, 1.0],
                transparency: 0.25,
            },
        ],
    };
    let mut plan = plan::plan_artifact_primitives(&snapshot).expect("plan");

    apply_visible_surface_triangle_counts(&mut plan, &[1]).expect("visibility counts");

    assert_eq!(plan.triangle_count, 3);
    assert_eq!(plan.triangle_reps.len(), 2);
    assert_eq!(plan.triangle_reps[0].triangle_offset, 0);
    assert_eq!(plan.triangle_reps[0].triangle_count, 2);
    assert_eq!(plan.triangle_reps[0].visibility_counter_index, None);
    assert_eq!(plan.triangle_reps[1].triangle_offset, 2);
    assert_eq!(plan.triangle_reps[1].triangle_count, 1);
    assert_eq!(plan.triangle_reps[1].source_triangle_count(), 3);
    assert_eq!(plan.triangle_reps[1].visibility_counter_index, Some(0));
    assert_eq!(
        plan.triangle_reps[1].geometry_binding_size,
        9 * STD_VERTEX_STRIDE
    );
}

// Planning rejection paths.

#[test]
fn artifact_plan_rejects_direct_triangle_count_not_divisible_by_three() {
    let snapshot = RenderArtifactSnapshotDescriptor {
        snapshot_id: 1,
        layout_version: 2,
        scene_generation: 7,
        scene_bounds_min: [0.0, 0.0, 0.0],
        scene_bounds_max: [1.0, 1.0, 1.0],
        cull_pass_initialized: true,
        device_limits: patinae_scene::GpuDeviceLimits {
            max_buffer_size: 4096,
            max_storage_buffer_binding_size: 4096,
            max_compute_workgroups_per_dimension: 65_535,
            max_compute_invocations_per_workgroup: 256,
            max_compute_workgroup_size_x: 256,
            max_compute_workgroup_size_y: 1,
            max_compute_workgroup_size_z: 1,
        },
        buffers: vec![
            RenderArtifactBufferDescriptor {
                handle: handle(10),
                role: RenderArtifactBufferRole::SceneColorLut,
                size: 64,
                stride: COLOR_LUT_STRIDE,
                element_count: 1,
            },
            RenderArtifactBufferDescriptor {
                handle: handle(11),
                role: RenderArtifactBufferRole::StdVertices,
                size: 4 * STD_VERTEX_STRIDE,
                stride: STD_VERTEX_STRIDE,
                element_count: 4,
            },
        ],
        reps: vec![RenderArtifactRepDescriptor {
            object_id: 1,
            rep_kind: RenderArtifactRepKind::Cartoon,
            topology: RenderArtifactPrimitiveTopology::TriangleList,
            geometry: handle(11),
            count: None,
            indirect: None,
            element_count: 4,
            max_element_count: 4,
            atom_offset: 5,
            atom_count: 1,
            material_rgba: [1.0, 0.0, 0.0, 1.0],
            transparency: 0.0,
        }],
    };

    let err = match plan::plan_artifact_primitives(&snapshot) {
        Ok(_) => panic!("plan should fail"),
        Err(err) => err,
    };

    assert!(err.contains("Cartoon artifact vertex count 4 is not divisible by 3"));
}

#[test]
fn artifact_plan_rejects_uninitialized_cull_counts() {
    let snapshot = RenderArtifactSnapshotDescriptor {
        snapshot_id: 1,
        layout_version: 2,
        scene_generation: 7,
        scene_bounds_min: [0.0, 0.0, 0.0],
        scene_bounds_max: [1.0, 1.0, 1.0],
        cull_pass_initialized: false,
        device_limits: patinae_scene::GpuDeviceLimits {
            max_buffer_size: 4096,
            max_storage_buffer_binding_size: 4096,
            max_compute_workgroups_per_dimension: 65_535,
            max_compute_invocations_per_workgroup: 256,
            max_compute_workgroup_size_x: 256,
            max_compute_workgroup_size_y: 1,
            max_compute_workgroup_size_z: 1,
        },
        buffers: vec![RenderArtifactBufferDescriptor {
            handle: handle(10),
            role: RenderArtifactBufferRole::SceneColorLut,
            size: 64,
            stride: COLOR_LUT_STRIDE,
            element_count: 1,
        }],
        reps: vec![RenderArtifactRepDescriptor {
            object_id: 1,
            rep_kind: RenderArtifactRepKind::Sphere,
            topology: RenderArtifactPrimitiveTopology::SphereInstances,
            geometry: handle(11),
            count: Some(handle(12)),
            indirect: Some(handle(13)),
            element_count: 0,
            max_element_count: 5,
            atom_offset: 0,
            atom_count: 5,
            material_rgba: [1.0, 0.0, 0.0, 1.0],
            transparency: 0.0,
        }],
    };

    let err = match plan::plan_artifact_primitives(&snapshot) {
        Ok(_) => panic!("plan should fail"),
        Err(err) => err,
    };

    assert!(err.contains("cull pass is not initialized"));
}
