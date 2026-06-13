#!/usr/bin/env python3
"""
Build a binary CCD (Chemical Component Dictionary) bond cache from
components-pub.sdf.gz for patinae runtime lookup.

Usage:
    python scripts/build_ccd_cache.py components-pub.sdf.gz [output.bin]

Output format: MessagePack dict  { resn: { "a": [atom_names], "b": [[a1, a2, order], ...] } }
Bond order encoding: 1=Single, 2=Double, 3=Triple, 4=Aromatic

The output file defaults to ~/.patinae/resources/components.bin
"""

import gzip
import os
import sys
from pathlib import Path

try:
    import msgpack
except ImportError:
    print("Error: msgpack not installed. Run: pip install msgpack", file=sys.stderr)
    sys.exit(1)


def parse_sdf_stream(fh):
    """Yield (component_id, atom_names, bonds) from an SDF stream.

    Each bond is (atom1_name, atom2_name, bond_order).
    Bond order: 1=single, 2=double, 3=triple, 4=aromatic.
    """
    while True:
        # Line 1: molecule name
        name_line = fh.readline()
        if not name_line:
            break
        comp_id = name_line.strip()
        if not comp_id:
            # skip blank lines between records
            continue

        # Lines 2-3: program info, comment
        fh.readline()
        fh.readline()

        # Line 4: counts line
        counts_line = fh.readline().strip()
        if not counts_line:
            continue

        parts = counts_line.split()
        try:
            n_atoms = int(parts[0])
            n_bonds = int(parts[1])
        except (IndexError, ValueError):
            # skip to next $$$$
            for line in fh:
                if line.strip() == "$$$$":
                    break
            continue

        # Atom block: n_atoms lines
        atom_names = []
        for _ in range(n_atoms):
            atom_line = fh.readline()
            # SDF atom line: x y z symbol ...
            # Columns (V2000): 0-9 x, 10-19 y, 20-29 z, 31-33 symbol
            # But CCD SDF may have atom names in different positions
            # Standard V2000: first 3 fields are coords, 4th is symbol
            aparts = atom_line.split()
            if len(aparts) >= 4:
                atom_names.append(aparts[3])  # element symbol as placeholder
            else:
                atom_names.append("?")

        # Bond block: n_bonds lines
        bonds = []
        for _ in range(n_bonds):
            bond_line = fh.readline()
            bparts = bond_line.split()
            if len(bparts) >= 3:
                try:
                    a1_idx = int(bparts[0]) - 1  # 1-based to 0-based
                    a2_idx = int(bparts[1]) - 1
                    order = int(bparts[2])
                    # Map SDF bond types: 1=single, 2=double, 3=triple, 4=aromatic
                    if order not in (1, 2, 3, 4):
                        order = 1
                    if 0 <= a1_idx < len(atom_names) and 0 <= a2_idx < len(atom_names):
                        bonds.append((atom_names[a1_idx], atom_names[a2_idx], order))
                except (ValueError, IndexError):
                    pass

        # Read until $$$$ (skip properties block)
        for line in fh:
            if line.strip() == "$$$$":
                break

        if comp_id and bonds:
            yield comp_id, atom_names, bonds


def parse_ccd_sdf_with_names(fh):
    """Parse CCD SDF where atom names are in the properties block (M  ...).

    CCD SDF files store PDB atom names in the properties section, not in the
    atom block (which only has element symbols). We need to look for atom name
    properties or use a different strategy.

    Actually, CCD SDF V2000 stores atom names/IDs via the _chem_comp_atom
    mapping. In practice, the SDF atom block has element symbols, and the
    actual PDB atom names are in a separate property or can be derived.

    For CCD specifically: the molecule name IS the component ID, and atom
    names in PDB convention are provided via data items. However, the
    standard SDF format doesn't natively support this.

    REVISED APPROACH: Parse the CCD CIF format instead, which has explicit
    _chem_comp_bond with atom names.
    """
    pass


def parse_ccd_cif(fh):
    """Parse components.cif (CCD in mmCIF format).

    Extracts _chem_comp_bond loops which contain:
    - _chem_comp_bond.comp_id
    - _chem_comp_bond.atom_id_1
    - _chem_comp_bond.atom_id_2
    - _chem_comp_bond.value_order  (SING, DOUB, TRIP, AROM)

    And _chem_comp_atom loops for atom name lists:
    - _chem_comp_atom.comp_id
    - _chem_comp_atom.atom_id
    """
    ORDER_MAP = {
        "SING": 1,
        "DOUB": 2,
        "TRIP": 3,
        "AROM": 4,
    }

    components = {}  # comp_id -> {"a": set(atom_names), "b": [(a1, a2, order)]}

    in_loop = False
    loop_tags = []
    current_category = None

    for raw_line in fh:
        line = (
            raw_line.strip()
            if isinstance(raw_line, str)
            else raw_line.decode("utf-8", errors="replace").strip()
        )

        if not line or line.startswith("#"):
            if in_loop and loop_tags:
                in_loop = False
                loop_tags = []
                current_category = None
            continue

        if line.startswith("data_"):
            in_loop = False
            loop_tags = []
            current_category = None
            continue

        if line == "loop_":
            in_loop = True
            loop_tags = []
            current_category = None
            continue

        if in_loop and line.startswith("_"):
            loop_tags.append(line)
            # Detect category
            if "_chem_comp_bond." in line:
                current_category = "bond"
            elif "_chem_comp_atom." in line:
                current_category = "atom"
            continue

        if in_loop and loop_tags and not line.startswith("_"):
            if current_category == "bond":
                _process_bond_line(line, loop_tags, components, ORDER_MAP)
            elif current_category == "atom":
                _process_atom_line(line, loop_tags, components)
            continue

        # Non-loop line starting with _ (single-value data item) — skip
        if line.startswith("_"):
            in_loop = False
            loop_tags = []
            current_category = None
            continue

        # If we hit a non-tag, non-data line outside a recognized loop, reset
        if not in_loop:
            continue

    return components


def _tokenize_cif_line(line):
    """Split a CIF data line into tokens, respecting quoted strings."""
    tokens = []
    i = 0
    n = len(line)
    while i < n:
        if line[i] in (" ", "\t"):
            i += 1
            continue
        if line[i] in ("'", '"'):
            quote = line[i]
            i += 1
            start = i
            while i < n and line[i] != quote:
                i += 1
            tokens.append(line[start:i])
            i += 1  # skip closing quote
        else:
            start = i
            while i < n and line[i] not in (" ", "\t"):
                i += 1
            tokens.append(line[start:i])
    return tokens


def _process_bond_line(line, loop_tags, components, order_map):
    """Process a single data line from a _chem_comp_bond loop."""
    tokens = _tokenize_cif_line(line)
    if len(tokens) < len(loop_tags):
        return

    tag_map = {tag: idx for idx, tag in enumerate(loop_tags)}
    comp_idx = tag_map.get("_chem_comp_bond.comp_id")
    a1_idx = tag_map.get("_chem_comp_bond.atom_id_1")
    a2_idx = tag_map.get("_chem_comp_bond.atom_id_2")
    order_idx = tag_map.get("_chem_comp_bond.value_order")

    if any(x is None for x in (comp_idx, a1_idx, a2_idx)):
        return

    comp_id = tokens[comp_idx]
    atom1 = tokens[a1_idx]
    atom2 = tokens[a2_idx]
    order_str = (
        tokens[order_idx]
        if order_idx is not None and order_idx < len(tokens)
        else "SING"
    )
    order = order_map.get(order_str, 1)

    if comp_id not in components:
        components[comp_id] = {"a": set(), "b": []}
    components[comp_id]["a"].add(atom1)
    components[comp_id]["a"].add(atom2)
    components[comp_id]["b"].append((atom1, atom2, order))


def _process_atom_line(line, loop_tags, components):
    """Process a single data line from a _chem_comp_atom loop."""
    tokens = _tokenize_cif_line(line)
    if len(tokens) < len(loop_tags):
        return

    tag_map = {tag: idx for idx, tag in enumerate(loop_tags)}
    comp_idx = tag_map.get("_chem_comp_atom.comp_id")
    atom_idx = tag_map.get("_chem_comp_atom.atom_id")

    if comp_idx is None or atom_idx is None:
        return

    comp_id = tokens[comp_idx]
    atom_name = tokens[atom_idx]

    if comp_id not in components:
        components[comp_id] = {"a": set(), "b": []}
    components[comp_id]["a"].add(atom_name)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    input_path = sys.argv[1]

    if len(sys.argv) >= 3:
        output_path = Path(sys.argv[2])
    else:
        output_path = Path.home() / ".patinae" / "resources" / "components.bin"

    print(f"Input:  {input_path}")
    print(f"Output: {output_path}")

    # Detect format from filename
    is_sdf = "sdf" in input_path.lower()
    is_gzipped = input_path.endswith(".gz")

    if is_sdf:
        print("Error: SDF format does not preserve PDB atom names.", file=sys.stderr)
        print("Please use components.cif.gz from wwPDB instead:", file=sys.stderr)
        print(
            "  https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz",
            file=sys.stderr,
        )
        sys.exit(1)

    print("Parsing CCD CIF (this may take a few minutes)...")

    if is_gzipped:
        fh = gzip.open(input_path, "rt", encoding="utf-8", errors="replace")
    else:
        fh = open(input_path, "r", encoding="utf-8", errors="replace")

    with fh:
        components = parse_ccd_cif(fh)

    # Convert sets to sorted lists for deterministic output
    result = {}
    for comp_id, data in components.items():
        if not data["b"]:
            continue  # skip components with no bonds
        result[comp_id] = {
            "a": sorted(data["a"]),
            "b": data["b"],  # list of [atom1, atom2, order]
        }

    print(f"Parsed {len(result)} components with bonds")

    # Serialize to msgpack
    output_path.parent.mkdir(parents=True, exist_ok=True)
    packed = msgpack.packb(result, use_bin_type=True)

    with open(output_path, "wb") as f:
        f.write(packed)

    size_mb = len(packed) / (1024 * 1024)
    print(f"Written {output_path} ({size_mb:.1f} MB)")

    # Print some stats
    hem = result.get("HEM")
    if hem:
        print(f"  HEM: {len(hem['a'])} atoms, {len(hem['b'])} bonds")
    atp = result.get("ATP")
    if atp:
        print(f"  ATP: {len(atp['a'])} atoms, {len(atp['b'])} bonds")


if __name__ == "__main__":
    main()
