//! Performance benchmarks for pymol-rs using a large biological assembly.
//!
//! Run with: `BENCH_FILE=path/to/structure.cif.gz cargo bench -p pymol-bench`

mod common;

use std::io::Read;
use std::time::Duration;

use criterion::{black_box, criterion_group, criterion_main, BatchSize, Criterion};

use pymol_algos::{Dssp, PyMolDss};
use pymol_mol::dss::assign_secondary_structure;
use pymol_mol::RepMask;
use pymol_render::picking::Ray;
use pymol_render::{
    CartoonRep, LineRep, Representation, SphereRep, StickRep, SurfaceRep,
};
use pymol_scene::{MoleculeObject, ObjectRegistry, Picker};
use pymol_select::{evaluate, parse, EvalContext};

use common::{ColorResolverOwned, MOLECULE, MOLECULE_WITH_DSS, RAW_TEXT};

// ---------------------------------------------------------------------------
// Group 1: CIF Loading
// ---------------------------------------------------------------------------

fn cif_loading(c: &mut Criterion) {
    let mut group = c.benchmark_group("cif_loading");
    group.sample_size(10);
    group.warm_up_time(Duration::from_secs(3));
    group.measurement_time(Duration::from_secs(15));

    let path = common::test_file_path();

    group.bench_function("decompress_only", |b| {
        b.iter(|| {
            let file = std::fs::File::open(&path).unwrap();
            let mut decoder = flate2::read::GzDecoder::new(file);
            let mut text = String::new();
            decoder.read_to_string(&mut text).unwrap();
            black_box(text.len());
        });
    });

    group.bench_function("tokenize_only", |b| {
        let text = &*RAW_TEXT;
        b.iter(|| {
            let tokens = pymol_io::cif::lexer::tokenize(text);
            black_box(tokens.len());
        });
    });

    group.bench_function("parse_full", |b| {
        let text = &*RAW_TEXT;
        b.iter(|| {
            let mol = pymol_io::cif::read_cif_str(text).unwrap();
            black_box(mol.atom_count());
        });
    });

    group.bench_function("load_from_disk", |b| {
        b.iter(|| {
            let mol = pymol_io::cif::read_cif(&path).unwrap();
            black_box(mol.atom_count());
        });
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Group 2: DSS (Secondary Structure Assignment)
// ---------------------------------------------------------------------------

fn dss_bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("dss");
    group.sample_size(20);
    group.warm_up_time(Duration::from_secs(2));
    group.measurement_time(Duration::from_secs(10));

    // Force lazy init before benchmarking
    let _ = &*MOLECULE;

    let assigner = PyMolDss::default();
    group.bench_function("assign_secondary_structure", |b| {
        b.iter_batched(
            || MOLECULE.clone(),
            |mut mol| {
                assign_secondary_structure(&mut mol, 0, &assigner);
                black_box(mol.atom_count());
            },
            BatchSize::LargeInput,
        );
    });

    let dssp_assigner = Dssp::default();
    group.bench_function("assign_secondary_structure_dssp", |b| {
        b.iter_batched(
            || MOLECULE.clone(),
            |mut mol| {
                assign_secondary_structure(&mut mol, 0, &dssp_assigner);
                black_box(mol.atom_count());
            },
            BatchSize::LargeInput,
        );
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Group 3: Selection Evaluation
// ---------------------------------------------------------------------------

fn selection_bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("selection");
    group.sample_size(20);
    group.warm_up_time(Duration::from_secs(2));
    group.measurement_time(Duration::from_secs(10));

    let mol = &*MOLECULE;

    let cases = [
        ("select_all", "all"),
        ("select_chain_A", "chain A"),
        ("select_protein", "protein"),
        ("select_resn_ala", "resn ALA"),
        ("select_within_5_chain_A", "within 5 of chain A"),
        (
            "select_complex",
            "protein and not elem H",
        ),
    ];

    for (name, sele_str) in &cases {
        let expr = parse(sele_str).unwrap_or_else(|e| {
            panic!("failed to parse selection '{}': {}", sele_str, e);
        });
        let ctx = EvalContext::single(mol);

        group.bench_function(*name, |b| {
            b.iter(|| {
                let result = evaluate(&expr, &ctx).unwrap();
                black_box(result);
            });
        });
    }

    group.finish();
}

// ---------------------------------------------------------------------------
// Group 4: Representation Build (CPU geometry generation)
// ---------------------------------------------------------------------------

fn representation_bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("representation_build");
    group.sample_size(10);
    group.warm_up_time(Duration::from_secs(2));
    group.measurement_time(Duration::from_secs(20));

    let mol = &*MOLECULE_WITH_DSS;
    let color_owner = ColorResolverOwned::new();
    let colors = color_owner.resolver();
    let settings = common::default_settings();

    // Prepare a molecule clone with all reps visible
    let mut mol_all_reps = mol.clone();
    for atom in mol_all_reps.atoms_mut() {
        atom.repr.visible_reps = RepMask::ALL;
    }
    let coord_set_all = mol_all_reps
        .get_coord_set(0)
        .expect("molecule has no coord set");

    // -- Surface benchmarks --

    group.bench_function("surface_marching_cubes_q-2", |b| {
        b.iter_batched(
            || {
                let mut rep = SurfaceRep::new();
                rep.set_quality(-2);
                rep
            },
            |mut rep| {
                rep.build(&mol_all_reps, coord_set_all, &colors, &settings);
                black_box(rep.primitive_count());
            },
            BatchSize::LargeInput,
        );
    });

    group.bench_function("surface_marching_cubes_q0", |b| {
        b.iter_batched(
            || {
                let mut rep = SurfaceRep::new();
                rep.set_quality(0);
                rep
            },
            |mut rep| {
                rep.build(&mol_all_reps, coord_set_all, &colors, &settings);
                black_box(rep.primitive_count());
            },
            BatchSize::LargeInput,
        );
    });

    // -- Other representations --

    group.bench_function("cartoon_build", |b| {
        b.iter_batched(
            CartoonRep::new,
            |mut rep| {
                rep.build(&mol_all_reps, coord_set_all, &colors, &settings);
                black_box(rep.primitive_count());
            },
            BatchSize::LargeInput,
        );
    });

    group.bench_function("stick_build", |b| {
        b.iter_batched(
            StickRep::new,
            |mut rep| {
                rep.build(&mol_all_reps, coord_set_all, &colors, &settings);
                black_box(rep.primitive_count());
            },
            BatchSize::LargeInput,
        );
    });

    group.bench_function("sphere_build", |b| {
        b.iter_batched(
            SphereRep::new,
            |mut rep| {
                rep.build(&mol_all_reps, coord_set_all, &colors, &settings);
                black_box(rep.primitive_count());
            },
            BatchSize::LargeInput,
        );
    });

    group.bench_function("line_build", |b| {
        b.iter_batched(
            LineRep::new,
            |mut rep| {
                rep.build(&mol_all_reps, coord_set_all, &colors, &settings);
                black_box(rep.primitive_count());
            },
            BatchSize::LargeInput,
        );
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Group 5: Color Resolution
// ---------------------------------------------------------------------------

fn color_bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("color_resolution");
    group.sample_size(50);
    group.warm_up_time(Duration::from_secs(1));
    group.measurement_time(Duration::from_secs(5));

    let mol = &*MOLECULE;
    let color_owner = ColorResolverOwned::new();
    let colors = color_owner.resolver();

    group.bench_function("resolve_all_atoms", |b| {
        b.iter(|| {
            let mut count = 0u32;
            for atom in mol.atoms() {
                black_box(colors.resolve_atom(atom));
                count += 1;
            }
            black_box(count);
        });
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Group 6: Picking
// ---------------------------------------------------------------------------

fn picking_bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("picking");
    group.sample_size(20);
    group.warm_up_time(Duration::from_secs(2));
    group.measurement_time(Duration::from_secs(10));

    let mol = &*MOLECULE;

    // Set up an ObjectRegistry with the molecule
    let mut registry = ObjectRegistry::new();
    registry.add(MoleculeObject::new(mol.clone()));

    // Create a ray pointing into the center of the molecule (simulates a click)
    // Use the molecule's bounding box center as target
    let cs = mol.get_coord_set(0).expect("molecule has no coord set");
    let (min, max) = cs.bounding_box().expect("empty coord set");
    let center = (min + max) * 0.5;

    // Ray from in front of the molecule, pointing toward center
    let origin = [center.x, center.y, center.z + 200.0];
    let direction = [0.0, 0.0, -1.0];
    let ray_hit = Ray::new(origin, direction);

    // Ray that misses entirely (pointing away)
    let ray_miss = Ray::new([1e6, 1e6, 1e6], [0.0, 0.0, -1.0]);

    group.bench_function("pick_ray_hit", |b| {
        let mut picker = Picker::new();
        b.iter(|| {
            let result = picker.pick_ray(&ray_hit, &mut registry);
            black_box(result);
        });
    });

    group.bench_function("pick_ray_miss", |b| {
        let mut picker = Picker::new();
        b.iter(|| {
            let result = picker.pick_ray(&ray_miss, &mut registry);
            black_box(result);
        });
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Criterion harness
// ---------------------------------------------------------------------------

criterion_group! {
    name = io_benches;
    config = Criterion::default();
    targets = cif_loading
}

criterion_group! {
    name = compute_benches;
    config = Criterion::default();
    targets = dss_bench, selection_bench, color_bench, picking_bench
}

criterion_group! {
    name = render_benches;
    config = Criterion::default();
    targets = representation_bench
}

criterion_main!(io_benches, compute_benches, render_benches);
