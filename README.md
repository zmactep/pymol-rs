<p align="center">
  <img src="images/pymol-rs.png" alt="PyMOL-RS" width="200">
</p>

<h1 align="center">PyMOL-RS</h1>

<p align="center">
  <strong>PyMOL, reimagined in Rust.</strong><br>
  Same power. Modern core. Zero legacy baggage.
</p>

<p align="center">
  <img src="https://img.shields.io/badge/license-BSD--3--Clause-blue.svg" alt="License">
  <img src="https://img.shields.io/badge/rust-1.70%2B-orange.svg" alt="Rust">
  <img src="https://img.shields.io/badge/status-beta-blue.svg" alt="Status">
</p>

---

## Why?

PyMOL is the gold standard for molecular visualization — but it's 25 years of C/C++/Python accumulated into a monolith. PyMOL-RS is a clean-room rewrite that keeps the familiar command language and workflow while replacing the engine:

| | PyMOL (classic) | PyMOL-RS |
|---|---|---|
| **Language** | C / C++ / Python | Rust + wgpu |
| **Rendering** | OpenGL 2.x fixed pipeline | WebGPU (wgpu), GPU impostors |
| **Architecture** | Monolithic | Independent crates |
| **Python** | Embedded CPython | PyO3 bindings (optional) |
| **Memory safety** | Manual | Guaranteed at compile time |
| **Cross-platform** | Build scripts per OS | Single `cargo build` |

> **Beta status.** Core visualization, selection language, structural alignment, sessions, and Python API are fully functional. Electron density maps and movie export are on the roadmap.

## Quick Start

### Pre-built binaries

Grab the latest release for your platform from [Releases](https://github.com/zmactep/pymol-rs/releases/latest).

**Standalone executable** — no Python needed:
```bash
pymol-rs protein.pdb
```

**Python wheel** — includes both the CLI and `pymol_rs` module:
```bash
pip install pymol_rs-<version>-<platform>.whl
```

### Build from source

```bash
git clone https://github.com/zmactep/pymol-rs
cd pymol-rs
make release && make run
```

| Target | Command |
|--------|---------|
| Release build | `make release` |
| Debug build | `make debug` |
| Python wheel | `make python` (requires [maturin](https://github.com/PyO3/maturin)) |
| Tests | `make test` |

## What works

**Formats:** PDB · mmCIF · bCIF · MOL2 · SDF/MOL · XYZ · GRO (+ gzip)

**Representations:** Spheres (GPU impostors) · Sticks · Lines · Cartoon · Ribbon · Surface (SAS/SES/VdW) · Mesh · Dots · Labels

**Selection language** — full PyMOL-compatible syntax:
```
chain A and name CA
byres around 5 ligand
polymer and not solvent
```

**Commands** — familiar PyMOL verbs with tab completion: `load`, `show`, `hide`, `color`, `select`, `zoom`, `center`, `orient`, `png`, `ray`, …

**Structural analysis:**
- Alignment — Kabsch superposition & CE structural alignment
- Measurements — distance, angle, dihedral with visual feedback
- Crystallographic symmetry — `symexp` with all 230 space groups
- Secondary structure — `dss` assignment from geometry

**Sessions** — save and load your sessions with high-efficiency `.prs` file format or use your old PyMOL sessions with `.pse` parser.

**Ray tracing** — offline GPU ray tracing with BVH acceleration, shadows, and transparency.

**GUI** — egui-based interface with command line, object panel, sequence viewer, mouse picking, and viewport:

<p align="center">
  <img src="images/interface.png" alt="Interface" width="800">
</p>

## Python API

```python
from pymol_rs import cmd

cmd.load("protein.pdb")
cmd.show("cartoon")
cmd.color("green", "chain A")
cmd.select("site", "byres around 5 ligand")
cmd.png("output.png", width=1920, height=1080)

# Extend with custom commands
def highlight(selection):
    cmd.color("yellow", selection)

cmd.extend("highlight", highlight)
cmd.show_gui()
```

## Architecture

```
pymol-rs/
├── pymol-mol         Core data: Atom, Bond, Molecule
├── pymol-io          Format parsers & writers
├── pymol-select      Selection language (parser + evaluator)
├── pymol-color       Colors, schemes, ramps
├── pymol-settings    Configuration system
├── pymol-algos       Molecular algorithms
├── pymol-render      wgpu rendering engine
├── pymol-raytracer   Offline ray tracing
├── pymol-scene       Viewer, camera, scene graph
├── pymol-session     Sessions save and load (`.prs` and `.pse`)
├── pymol-cmd         Command parser & executor
├── pymol-gui         GUI (egui)
└── pymol-python      Python bindings (PyO3 / maturin)
```

Each crate is independently usable. Want just the selection parser? `pymol-select`. Need to read PDB files in your pipeline? `pymol-io` + `pymol-mol`. No GUI tax.

## Roadmap

- [ ] Electron density maps (isomesh / isosurface)
- [ ] Movie export (video formats)
- [ ] Object groups

## License

[BSD 3-Clause](LICENSE)

## Acknowledgments

Inspired by [PyMOL](https://pymol.org/), created by Warren Lyford DeLano. This is an independent reimplementation, not affiliated with Schrödinger, Inc.
