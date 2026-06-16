use super::*;
use crate::shaders;
use patinae_scene::{
    GpuHandleKind, RenderArtifactBufferDescriptor, RenderArtifactBufferRole,
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
}

#[test]
fn raytrace_output_buffer_shader_uses_output_buffer_and_empty_node_guard() {
    let shader = shaders::ARTIFACT_RAYTRACE_OUTPUT_BUFFER;

    assert!(shader.contains("output_pixels"));
    assert!(shader.contains("global_id.y * u32(uniforms.viewport.x) + global_id.x"));
    assert!(shader.contains("node.left_or_first == 0xffffffffu"));
    assert!(!shader.contains("texture_storage_2d"));
    assert!(!shader.contains("textureStore("));
}

#[test]
fn artifact_wgsl_modules_parse_and_validate() {
    let modules: [(&str, String); 8] = [
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
        ("ray.artifact.bvh", shaders::ARTIFACT_BVH.to_string()),
        (
            "ray.artifact.raytrace",
            shaders::ARTIFACT_RAYTRACE_OUTPUT_BUFFER.to_string(),
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
                size: 96 * STD_VERTEX_STRIDE,
                stride: STD_VERTEX_STRIDE,
                element_count: 96,
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
            max_element_count: 96,
            atom_offset: 5,
            atom_count: 8,
            material_rgba: [0.5, 0.5, 0.5, 1.0],
            transparency: 0.5,
        }],
    };

    let plan = plan::plan_artifact_primitives(&snapshot).expect("plan");

    assert_eq!(plan.triangle_count, 32);
    assert_eq!(plan.triangle_reps[0].rep_slot, 6);
    assert_eq!(plan.triangle_reps[0].vertex_count, 96);
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
