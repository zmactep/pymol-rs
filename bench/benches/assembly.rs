//! Performance benchmarks for patinae using a large biological assembly.
//!
//! Run with: `BENCH_FILE=path/to/structure.cif.gz cargo bench -p patinae-bench`

mod common;

use std::io::Read;
use std::time::Duration;

use criterion::{black_box, criterion_group, criterion_main, BatchSize, Criterion};

use patinae_algos::{Dssp, PyMolDss};
use patinae_mol::dss::assign_secondary_structure;
use patinae_mol::DEFAULT_BOND_TOLERANCE;
use patinae_select::{evaluate, parse, EvalContext};

use common::{ColorResolverOwned, MOLECULE, RAW_TEXT};

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
            let tokens = patinae_io::cif::lexer::tokenize(text);
            black_box(tokens.len());
        });
    });

    group.bench_function("parse_full", |b| {
        let text = &*RAW_TEXT;
        b.iter(|| {
            let mol = patinae_io::cif::read_cif_str(text).unwrap();
            black_box(mol.atom_count());
        });
    });

    group.bench_function("load_from_disk", |b| {
        b.iter(|| {
            let mol = patinae_io::cif::read_cif(&path).unwrap();
            black_box(mol.atom_count());
        });
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Group 2: Bond Generation
// ---------------------------------------------------------------------------

fn bond_generation(c: &mut Criterion) {
    let mut group = c.benchmark_group("bond_generation");
    group.sample_size(10);
    group.warm_up_time(Duration::from_secs(2));
    group.measurement_time(Duration::from_secs(15));

    let _ = &*MOLECULE;

    group.bench_function("generate_bonds", |b| {
        b.iter_batched(
            || {
                let mut mol = MOLECULE.clone();
                mol.clear_bonds();
                mol
            },
            |mut mol| {
                mol.generate_bonds(DEFAULT_BOND_TOLERANCE);
                black_box(mol.bond_count());
            },
            BatchSize::LargeInput,
        );
    });

    group.finish();
}

// ---------------------------------------------------------------------------
// Group 3: DSS (Secondary Structure Assignment)
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
// Group 4: Selection Evaluation
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
        ("select_complex", "protein and not elem H"),
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
    targets = bond_generation, dss_bench, selection_bench, color_bench
}

criterion_main!(io_benches, compute_benches);
