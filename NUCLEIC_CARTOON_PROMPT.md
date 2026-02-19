# Task: Nucleic Acid Cartoon Representation

## Goal

Implement cartoon rendering for nucleic acids (DNA/RNA) in `pymol-rs`, similar to PyMOL's default view.
Test structure: **6YOV** (nucleosome with histone octamer + DNA). Currently, only protein chains render a cartoon — DNA/RNA are invisible in cartoon mode.

## Reference

PyMOL's nucleic acid cartoon (see attached image) shows:
- **Smooth backbone tube** (orange): a constant-radius cylinder spline through the phosphate/sugar backbone
- **Base sticks** (blue): thin flat rods or rectangles from the sugar C1' atom out to the base ring center

For this implementation, start with **backbone tube only** (no base geometry for now — that's a follow-up). If time permits, add base sticks as a second pass.

---

## Codebase orientation

```
crates/
  pymol-mol/
    src/
      molecule.rs       — ObjectMolecule, classify_atoms() (C4' gets GUIDE flag)
      residue.rs        — is_nucleotide(), nucleotide_to_char(), NUCLEOTIDES list
      flags.rs          — AtomFlags: NUCLEIC, GUIDE, PROTEIN, POLYMER
  pymol-render/
    src/representation/cartoon/
      backbone.rs       — extract_backbone_segments() → protein only (line 146: if !residue.is_protein() { continue; })
      spline.rs         — Catmull-Rom / displacement-based interpolation
      frame.rs          — Frenet-Serret frames along spline
      geometry.rs       — Extrusion shapes (helix ribbon, sheet arrow, loop tube)
      pipeline.rs       — generate_segment_cartoon() — main cartoon per-segment builder
      mod.rs            — CartoonRep::build() entry point
```

### Key insight

`classify_atoms()` already sets `AtomFlags::GUIDE` on **C4'** atoms (or C4*) for every nucleotide residue.
`extract_backbone_segments()` skips any non-protein residue at line 146. The fix is to add an analogous extraction path for nucleic acids.

---

## Implementation plan

### Step 1 — Extract nucleic backbone segments

In `backbone.rs`, add a new public function:

```rust
pub fn extract_nucleic_segments(
    molecule: &ObjectMolecule,
    coord_set: &CoordSet,
    colors: &ColorResolver,
    gap_cutoff: i32,
    rep_mask: RepMask,
) -> Vec<BackboneSegment>
```

Logic (mirrors `extract_backbone_segments` for protein, but for nucleic):

For each chain → for each residue where `residue.is_nucleic()`:
1. Find **C4'** (or C4*) — the guide atom, like CA for proteins
2. Check `ca_atom.repr.visible_reps.is_visible(rep_mask)` (same as protein path)
3. Get C4' coordinates from `coord_set`
4. **Gap detection**: if `|resv - prev_resv| > gap_cutoff` → start new segment
5. **Orientation vector**: C4' → C1' direction (C1' connects the sugar to the base).
   - Try `residue.find_by_name("C1'")`, fallback to `"C1*"`, then `Vec3::new(0.0, 1.0, 0.0)`
6. **Color**: use `colors.resolve_cartoon(c4_atom, molecule)`
7. **SecondaryStructure**: nucleic acids have no helix/sheet — use `SecondaryStructure::Loop` for all residues
8. Create `GuidePoint` and push to current segment; push segment on gap or chain end

> Note: `ResidueView` methods available: `.is_nucleic()`, `.find_by_name(name)`, `.resv()`.
> Same structs used for protein: `GuidePoint`, `BackboneSegment` — reuse as-is.

---

### Step 2 — Integrate into CartoonRep::build()

In `mod.rs`, inside `CartoonRep::build()` (after the existing protein loop):

```rust
// --- Nucleic acid backbone tube ---
let nucleic_segments = extract_nucleic_segments(
    molecule, coord_set, colors, gap_cutoff, RepMask::CARTOON,
);

for mut segment in nucleic_segments {
    if segment.len() < 2 {
        continue;
    }
    // No helix/sheet smoothing for nucleic — just orientation smoothing
    smooth_orientations(&mut segment, smooth_settings.smooth_cycles);

    // Generate tube geometry (loop-style for all residues)
    // Use generate_segment_cartoon with all-Loop SS — this will produce a tube
    let (seg_verts, seg_idxs) = generate_segment_cartoon(
        &segment,
        &interp_settings,
        &pipeline_settings,
        index_offset,
    );
    index_offset += seg_verts.len() as u32;
    vertices.extend(seg_verts);
    indices.extend(seg_idxs);
}
```

> The existing cartoon pipeline (`generate_segment_cartoon`) already produces a **round tube** for `SecondaryStructure::Loop` residues. Since we're setting all nucleotide residues to `Loop`, the backbone will render as a smooth tube automatically — no new geometry code needed.

---

### Step 3 — Also handle RibbonRep

In `build_cartoon_geometry()` (used by `RibbonRep`), the same `extract_backbone_segments` is called. Add an analogous call to `extract_nucleic_segments` and append those segments to the list:

```rust
let mut segments = extract_backbone_segments(molecule, coord_set, colors, gap_cutoff, rep_mask);
let nucleic = extract_nucleic_segments(molecule, coord_set, colors, gap_cutoff, rep_mask);
segments.extend(nucleic);
```

---

### Step 4 — Base sticks (optional, do if time allows)

If backbone tube works, add base geometry. For each nucleotide residue in a segment:

1. Find `C1'` position (already used for orientation)
2. Find base attachment atom: `N9` for purines (A, G, DA, DG), `N1` for pyrimidines (C, T, U, DC, DT)
3. Find base ring center: average positions of the 5- or 6-membered ring atoms
4. Draw a flat filled rectangle from C1' → base ring center:
   - Width: ~2.5 Å (PyMOL uses `cartoon_ladder_color` for bases)
   - Normal: perpendicular to the C4'→C1' axis

This can be deferred — the backbone tube alone is a big improvement over nothing.

---

### Step 5 — BACKLOG.md update

Update `BACKLOG.md` in the repo root:
- Mark the GRO chain segmentation / nucleic cartoon features appropriately
- Add a new section for "Nucleic Acid Cartoon" describing what's done and what remains (base geometry)

---

## Testing

1. `cargo build` must succeed with no warnings in new code
2. Add a unit test in `backbone.rs`:

```rust
#[test]
fn test_extract_nucleic_segments_empty_on_protein() {
    // A pure protein molecule should yield zero nucleic segments
    // Construct a minimal ObjectMolecule with one ALA residue, verify extract_nucleic_segments returns []
}
```

3. Manual test: load `_tests/complex_fail.pdb` (has protein + water) and `fetch 6yov` — nucleic chains must render a visible backbone tube in cartoon mode.

---

## Constraints / style notes

- Do **not** touch the protein cartoon pipeline — additive changes only
- `SecondaryStructure::Loop` is the correct fallback for nucleotide residues (no DSS for nucleic acids)
- Atom name lookup: prefer `"C4'"` then `"C4*"` (some old PDB files use asterisk notation)
- Gap cutoff for nucleic acids: use the same `gap_cutoff` parameter as protein (default 10 in the codebase)
- No new crate dependencies
- Rust edition 2021, `#[allow(dead_code)]` where needed for new fields

---

## Files to modify

| File | Change |
|---|---|
| `crates/pymol-render/src/representation/cartoon/backbone.rs` | Add `extract_nucleic_segments()` |
| `crates/pymol-render/src/representation/cartoon/mod.rs` | Call nucleic extraction in `CartoonRep::build()` |
| `crates/pymol-render/src/representation/cartoon/mod.rs` | Update `build_cartoon_geometry()` for RibbonRep |
| `BACKLOG.md` | Update feature status |

---

## Success criteria

- `cargo test -p pymol-render` passes
- Loading `6yov.cif` shows a smooth orange-ish tube for the DNA chains alongside protein cartoon
- No regressions on existing protein cartoon tests
