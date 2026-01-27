<p align="center">
  <img src="images/pymol-rs.png" alt="PyMOL-RS Logo" width="200">
</p>

<h1 align="center">PyMOL-RS</h1>

<p align="center">
  <strong>A Rust reimplementation of PyMOL molecular visualization</strong>
</p>

<p align="center">
  <a href="#features">Features</a> •
  <a href="#quick-start">Quick Start</a> •
  <a href="#usage">Usage</a> •
  <a href="#architecture">Architecture</a> •
  <a href="#roadmap">Roadmap</a> •
  <a href="#license">License</a>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/license-BSD--3--Clause-blue.svg" alt="License">
  <img src="https://img.shields.io/badge/rust-1.70%2B-orange.svg" alt="Rust">
  <img src="https://img.shields.io/badge/status-alpha-yellow.svg" alt="Status">
</p>

---

## Overview

PyMOL-RS is a modern molecular visualization system written in Rust. It aims to provide:

- **GPU-accelerated rendering** via wgpu (WebGPU) for cross-platform graphics
- **PyMOL-compatible command syntax** for familiar workflows
- **Full selection language** support for precise atom selection
- **High performance** through Rust's zero-cost abstractions

## Features

### File Formats

- **PDB** - Protein Data Bank format
- **mmCIF** - Macromolecular Crystallographic Information File
- **MOL2** - TRIPOS MOL2 format
- **SDF/MOL** - MDL Structure Data File
- **XYZ** - Simple coordinate format
- Automatic compression support (gzip)

### Molecular Representations

| Representation | Description |
|----------------|-------------|
| **Spheres** | Van der Waals spheres (GPU impostor rendering) |
| **Sticks** | Cylinder-based bonds |
| **Lines** | Wireframe bonds |
| **Cartoon** | Secondary structure (helices, sheets, loops) |
| **Ribbon** | Backbone ribbon trace |
| **Surface** | Molecular surface (SAS, SES, VdW) |
| **Mesh** | Triangle mesh surface |
| **Dots** | Dot surface representation |

### Selection Language

Full PyMOL-compatible selection syntax:

```
# Property selectors
name CA              # Atoms named CA
resn ALA             # Alanine residues
chain A              # Chain A
elem C               # Carbon atoms
ss H                 # Alpha helices

# Logical operators
chain A and name CA
not solvent
polymer or ligand

# Distance operators
around 5 of ligand   # Within 5Å of ligand
byres near_to 4 of site

# Grouping
bychain resn HEM     # Entire chains containing HEM
```

### Interactive Viewer

- Real-time 3D navigation (rotate, pan, zoom)
- Scene management and snapshots
- Screenshot export (PNG)
- Customizable key bindings

### Graphical Interface

PyMOL-RS includes a full graphical interface built with egui:

<p align="center">
  <img src="images/interface.png" alt="PyMOL-RS Interface" width="800">
</p>

The interface consists of three main areas:

| Component | Description |
|-----------|-------------|
| **Command Console** | Top area with command input, history, and colored output messages |
| **3D Viewer** | Central GPU-accelerated molecular rendering viewport |
| **Objects Panel** | Right sidebar listing loaded objects with quick action buttons |

The Objects panel provides quick access buttons for each object:

| Button | Action |
|--------|--------|
| **A** | Actions menu (delete, duplicate, rename) |
| **S** | Show representations |
| **H** | Hide representations |
| **L** | Label options |
| **C** | Color options |

## Quick Start

### Requirements

- Rust 1.70 or later
- GPU with Vulkan, Metal, or DirectX 12 support

### Build

```bash
git clone https://github.com/pymol-rs/pymol-rs
cd pymol-rs

# Build the GUI application
cargo build --release

# Or use make for convenience
make release
```

#### Build Targets

| Command | Description |
|---------|-------------|
| `make release` | Build Rust workspace (release) |
| `make debug` | Build Rust workspace (debug) |
| `make python` | Build Python wheel (requires [maturin](https://github.com/PyO3/maturin)) |
| `make all` | Build both Rust and Python |
| `make test` | Run tests |
| `make clean` | Clean all build artifacts |

### Run the GUI

```bash
# Run the GUI application
cargo run -p pymol-gui --release

# Or after building
./target/release/pymol-rs

# Load a file directly
cargo run -p pymol-gui --release -- protein.pdb
```

## Usage

### Programmatic API

```rust
use std::path::Path;
use pymol_io::read_file;
use pymol_mol::RepMask;
use pymol_scene::{run, MoleculeObject, Viewer};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Create viewer
    let mut viewer = Viewer::new();
    
    // Load molecule from file
    let mol = read_file(Path::new("protein.pdb"))?;
    
    // Create molecule object and configure representation
    let mut obj = MoleculeObject::new(mol);
    obj.show(RepMask::CARTOON);
    obj.show(RepMask::STICKS);
    
    // Add to viewer
    viewer.objects_mut().add(obj);
    viewer.center_all();
    
    // Run the viewer
    run(viewer)?;
    Ok(())
}
```

### Command System

```rust
use pymol_cmd::CommandExecutor;
use pymol_scene::Viewer;

let mut viewer = Viewer::new();
let mut executor = CommandExecutor::new();

// Execute PyMOL-style commands
executor.do_(&mut viewer, "load protein.pdb")?;
executor.do_(&mut viewer, "show cartoon")?;
executor.do_(&mut viewer, "hide lines")?;
executor.do_(&mut viewer, "color green, chain A")?;
executor.do_(&mut viewer, "zoom")?;
```

### Interactive Command Line

Run the interactive viewer and use PyMOL-style commands:

```bash
cargo run --example interactive -- 1IGT.cif
```

```
PyMOL> as cartoon
PyMOL> color by_chain
PyMOL> bg_color black
PyMOL> show sticks, chain A and resi 1-100
PyMOL> show surface, chain C
PyMOL> color atomic, rep sticks
PyMOL> hide cartoon, rep sticks
PyMOL> png output.png, 1920, 1080
PyMOL> quit
```

<p align="center">
  <img src="images/output.png" alt="PyMOL-RS Example Output" width="800">
</p>

### Python API

PyMOL-RS provides Python bindings with a familiar PyMOL-like API:

```bash
# Build and install the Python package
make python-dev
```

```python
from pymol import cmd

# Load a structure
cmd.load("protein.pdb")

# Show representations
cmd.show("cartoon")
cmd.show("sticks", "chain A and resi 1-100")

# Color by various schemes
cmd.color("green", "chain A")
cmd.color("cyan", "chain B")
cmd.color("atomic", "rep sticks")a

# Selection operations
cmd.select("active_site", "byres around 5 ligand")

# Export image
cmd.png("output.png", width=1920, height=1080)

# Access atom data
for atom in cmd.get_model("chain A and name CA").atom:
    print(f"{atom.resn} {atom.resi}: {atom.coord}")
```

## Controls

### Mouse

| Action | Control |
|--------|---------|
| Rotate | Left drag |
| Pan | Middle drag |
| Zoom | Right drag / Scroll |



## Architecture

PyMOL-RS is organized as a Rust workspace with modular crates:

```
pymol-rs/
├── pymol-mol       # Core molecular data structures (Atom, Bond, Molecule)
├── pymol-io        # File format parsers and writers
├── pymol-select    # Selection language parser and evaluator
├── pymol-color     # Color system (named colors, schemes, ramps)
├── pymol-render    # wgpu-based GPU rendering engine
├── pymol-scene     # Viewer, camera, and scene management
├── pymol-cmd       # Command parser and executor
├── pymol-settings  # Configuration and settings system
├── pymol-gui       # Graphical user interface (egui)
└── pymol-python    # Python bindings (PyO3, built separately with maturin)
```

### Layer Overview

```
┌─────────────────────────────────────────┐
│             Application Layer           │
│    pymol-gui, pymol-cmd, pymol-python   │
├─────────────────────────────────────────┤
│              Scene Layer                │
│             pymol-scene                 │
├─────────────────────────────────────────┤
│             Rendering Layer             │
│             pymol-render                │
├─────────────────────────────────────────┤
│              Domain Layer               │
│  pymol-mol, pymol-select, pymol-color   │
├─────────────────────────────────────────┤
│               I/O Layer                 │
│         pymol-io, pymol-settings        │
└─────────────────────────────────────────┘
```

## Roadmap

PyMOL-RS is in active development. The following features are planned but not yet available:

- **Labels** - Text labels for atoms and residues
- **Measurements** - Distance, angle, and dihedral measurements
- **Symmetry mates** - Crystallographic symmetry display
- **Electron density maps** - Map visualization and contouring
- **Movie export** - Video rendering
- **Session files** - Save/load PyMOL sessions (.pse)
- **Plugins** - Python plugin system
- **Ray tracing** - High-quality offline rendering

## License

This project is licensed under the [BSD 3-Clause License](LICENSE).

## Acknowledgments

PyMOL-RS is inspired by [PyMOL](https://pymol.org/), the original molecular visualization system created by Warren Lyford DeLano. This project is an independent Rust implementation and is not affiliated with Schrodinger, Inc.
