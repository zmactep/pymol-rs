# Patinae Benchmarks

The `patinae-bench` crate keeps the performance-sensitive paths visible without
turning Phase 8 into premature optimization work. Benchmarks are Criterion-based
and run under Cargo's `bench` profile, which keeps debug symbols enabled for
profiler-friendly builds.

## Hot Paths

- Parsing and loading: gzip decompression, CIF tokenization, full CIF parsing,
  object construction, and disk-backed `load_from_disk`.
- Molecular computation: secondary-structure assignment, selection evaluation,
  and per-atom color resolution.
- Scene/render sync: object-to-render input assembly, dirty flag handling,
  representation generation, and scene-store updates.
- Renderer passes: compute build, culling, shadows, atlas AO, picking,
  opaque/translucent rendering, marking, silhouette, SSAO, FXAA, and composite.

## CPU Benchmarks

`assembly` measures allocation-heavy and CPU-heavy molecular workloads. It
requires `BENCH_FILE` because results are only meaningful on a representative
structure.

```sh
BENCH_FILE="$TEST_STRUCTURES_DIR/7KP3-assembly1.cif.gz" cargo bench -p patinae-bench --bench assembly
```

Useful groups:

- `cif_loading`: decompression, tokenization, full parse, and disk load.
- `bond_generation`: distance-based bond generation on a parsed molecule.
- `dss`: PyMOL-style and DSSP-style secondary-structure assignment.
- `selection`: common selection expressions against the loaded molecule.
- `color_resolution`: per-atom color lookup throughput.

Use this benchmark first for allocator comparisons because parsing, molecule
cloning, and selection exercise allocation-heavy code paths.

## Render Loop Benchmarks

`render_loop` measures sync, render recording, queue submission, and resolved
GPU pass timings against a real molecule.

```sh
BENCH_FILE="$TEST_STRUCTURES_DIR/1fsd.cif" BENCH_SCENARIO=cartoon_only cargo bench -p patinae-bench --bench render_loop
BENCH_FILE="$TEST_STRUCTURES_DIR/3J3Q.cif" BENCH_SCENARIO=cartoon_only cargo bench -p patinae-bench --bench render_loop
```

Environment knobs:

- `BENCH_FILE`: absolute path or path relative to the repository root.
- `TEST_STRUCTURES_DIR`: directory containing local-only structure fixtures;
  fixture-gated tests skip when it is unset, and `render_loop` uses
  `TEST_STRUCTURES_DIR/1fsd.cif` when `BENCH_FILE` is unset.
- `BENCH_SCENARIO`: one of `cartoon_only`, `classic_small`, `classic_large`,
  `full_shadows`, `skripkin_cached`, `skripkin_forced_rebuild`, `sticks_heavy`,
  `marking`, or `surface_cartoon_large`. Use `cartoon_only` for web-vs-desktop
  picking parity; `classic_large` remains the `cartoon,sticks,spheres` stress
  scenario.
- `BENCH_REPS`: comma-separated representation override such as
  `cartoon`, `cartoon,sticks,spheres`, `sticks`, or `surface,cartoon`.
- `BENCH_MARKERS`: `none`, `selected`, `hover`, or `both`.
- `BENCH_W` / `BENCH_H`: viewport dimensions for fill-rate scaling.
- `BENCH_SHADOWS`, `BENCH_SHADOW_MAP`, `BENCH_SHADOW_BIAS`,
  `BENCH_SHADOW_INTENSITY`, `BENCH_SHADOW_PCF`: shadow-pass overrides.
- `BENCH_SKRIPKIN`, `BENCH_SKRIPKIN_DIRS`, `BENCH_SKRIPKIN_MAP`,
  `BENCH_SKRIPKIN_BIAS`, `BENCH_SKRIPKIN_INTENSITY`: atlas-AO overrides.

The Criterion score reports end-to-end frame-loop latency. The stderr summary
prints `FrameStats`: CPU sync/record/total timings plus per-pass GPU timings.
Negative GPU pass values mean the pass did not run in that frame, and missing
histograms usually mean the adapter does not expose timestamp queries.

## Raytrace Procedure Benchmarks

`raytrace_procedure` measures the raytracer's standalone GPU procedure with a
stage-level profile: CPU BVH build, main raytrace pipeline creation, primitive
and BVH uploads, output texture allocation, command recording, queue submit,
readback wait, and CPU downsample for antialiasing.

```sh
BENCH_RT_SCENE=surface_2k BENCH_RT_W=512 BENCH_RT_H=512 cargo bench -p patinae-bench --bench raytrace_procedure
BENCH_RT_SCENE=surface_8k BENCH_RT_W=1024 BENCH_RT_H=1024 BENCH_RT_AA=2 cargo bench -p patinae-bench --bench raytrace_procedure
```

Environment knobs:

- `BENCH_RT_SCENE`: `surface_512`, `surface_2k`, `surface_8k`,
  `surface_32k`, or `mixed_4k`.
- `BENCH_RT_W` / `BENCH_RT_H`: final output dimensions.
- `BENCH_RT_AA`: antialiasing factor passed to `RaytraceParams`.
- `BENCH_RT_SHADOWS`: set to `1` to include shadow rays.
- `BENCH_RT_MODE`: standalone edge/composite mode; keep `0` when comparing
  with the renderer-artifact command path, which currently supports normal mode.
- `BENCH_RT_TRANSPARENCY_MODE`: transparency mode uniform for shader sweeps.

For the live `ray` command's renderer-artifact path, set `PATINAE_RT_PROFILE=1`.
The plugin logs planning, indirect-count readback, surface-visibility counting,
temporary command-resource setup, BVH command setup, batch command count,
dispatch grid, batch/readback time, and total artifact-path time.

## Web Perf Harness

The local web harness lives at `web/examples/perf.html`. Keep local copies of
the acceptance fixtures in `web/assets/1fsd.cif` and `web/assets/3J3Q.cif`;
those large files are ignored, while `web/assets/.gitkeep` preserves the
directory.

```sh
cd web
npm run dev
```

Open `/examples/perf.html?fixture=3J3Q.cif&preset=cartoon_only` for the default
web-vs-desktop picking baseline. Use `preset=classic_large` only as a stress
smoke for `cartoon,sticks,spheres`; it is GPU-bound on large structures and is
not the parity gate.

## Allocator Policy

Do not add a global allocator to library crates. An app-only `mimalloc` change
belongs in the binary crate only, and only after benchmark comparison shows a
clear win: roughly 10 percent or better on allocation-heavy scenarios, with no
meaningful regression on the default render scenarios. If the comparison is
noisy or below that bar, keep the default allocator and record the decision in
the phase plan.

## Yield and Cancellation Policy

Do not add ad-hoc yield points to synchronous file load, DSS, scene sync, or
render paths. Moving synchronous load/DSS work into a cancellable background
task changes command execution, data ownership, and UI completion flow, so it is
separate architecture work. Existing async fetch paths are task-based with a
timeout, and the Python worker already has cancellation-style plumbing.
