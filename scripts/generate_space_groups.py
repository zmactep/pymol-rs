#!/usr/bin/env python3
"""Generate space_groups_data.rs from CCP4's syminfo.lib.

Source: International Tables for Crystallography, Vol. A (IUCr),
distributed as syminfo.lib in the CCP4 suite (derived from sgtbx).

Usage:
    python3 scripts/generate_space_groups.py scripts/syminfo.lib \
        > crates/pymol-algos/src/symmetry/space_groups_data.rs

    python3 scripts/generate_space_groups.py scripts/syminfo.lib --verify \
        crates/pymol-algos/src/symmetry/space_groups_data.rs
"""

from __future__ import annotations

import re
import sys
from fractions import Fraction
from typing import TypedDict


class SpaceGroup(TypedDict):
    xhm: str
    symops: list[str]
    cenops: list[str]
    old_names: list[str]


# ---------------------------------------------------------------------------
# Parsing syminfo.lib
# ---------------------------------------------------------------------------


def parse_syminfo(path: str) -> list[SpaceGroup]:
    """Parse syminfo.lib into a list of space group records."""
    groups: list[SpaceGroup] = []
    symops: list[str] = []
    cenops: list[str] = []
    old_names: list[str] = []
    xhm: str = ""
    in_block = False

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            if line == "begin_spacegroup":
                symops = []
                cenops = []
                old_names = []
                xhm = ""
                in_block = True
                continue
            if line == "end_spacegroup":
                if in_block and xhm:
                    groups.append(
                        SpaceGroup(
                            xhm=xhm,
                            symops=symops,
                            cenops=cenops,
                            old_names=old_names,
                        )
                    )
                in_block = False
                continue
            if not in_block:
                continue

            if line.startswith("symbol xHM"):
                m = re.search(r"'([^']*)'", line)
                if m:
                    xhm = m.group(1).strip()
            elif line.startswith("symbol old"):
                for m in re.finditer(r"'([^']*)'", line):
                    name = m.group(1).strip()
                    if name:
                        old_names.append(name)
            elif line.startswith("symop "):
                symops.append(line[6:].strip())
            elif line.startswith("cenop "):
                cenops.append(line[6:].strip())

    return groups


# ---------------------------------------------------------------------------
# Symmetry operation arithmetic (exact fractions)
# ---------------------------------------------------------------------------


def parse_symop_expr(expr: str) -> tuple[list[Fraction], Fraction]:
    """Parse one component like '-x+1/2' into ([cx, cy, cz], translation)."""
    coeffs = [Fraction(0)] * 3
    trans = Fraction(0)
    var_map = {"x": 0, "y": 1, "z": 2}

    chars = expr.replace(" ", "")
    n = len(chars)
    i = 0

    while i < n:
        sign = Fraction(1)
        if chars[i] == "-":
            sign = Fraction(-1)
            i += 1
        elif chars[i] == "+":
            i += 1

        if i >= n:
            break

        c = chars[i].lower()
        if c in var_map:
            coeffs[var_map[c]] = sign
            i += 1
        elif c.isdigit():
            j = i
            while j < n and chars[j].isdigit():
                j += 1
            num = int(chars[i:j])
            if j < n and chars[j] == "/":
                j += 1
                k = j
                while k < n and chars[k].isdigit():
                    k += 1
                den = int(chars[j:k])
                trans += sign * Fraction(num, den)
                i = k
            else:
                trans += sign * Fraction(num)
                i = j
        else:
            i += 1

    return coeffs, trans


def parse_cenop_translation(cenop: str) -> list[Fraction]:
    """Extract the translation vector from a cenop string."""
    parts = cenop.split(",")
    return [parse_symop_expr(p.strip())[1] for p in parts]


def parse_full_symop(op: str) -> list[tuple[list[Fraction], Fraction]]:
    """Parse a full symop string into components for x,y,z."""
    parts = op.split(",")
    return [parse_symop_expr(p.strip()) for p in parts]


def format_fraction(f: Fraction) -> str:
    """Format a fraction: '1/2', '2/3', etc."""
    if f.denominator == 1:
        return str(f.numerator)
    return f"{f.numerator}/{f.denominator}"


def format_component(coeffs: list[Fraction], trans: Fraction) -> str:
    """Format one component of a symmetry operation back to string."""
    var_names = ["x", "y", "z"]
    parts = []

    for var_idx, coeff in enumerate(coeffs):
        if coeff == 0:
            continue
        if coeff == 1:
            if parts:
                parts.append(f"+{var_names[var_idx]}")
            else:
                parts.append(var_names[var_idx])
        elif coeff == -1:
            parts.append(f"-{var_names[var_idx]}")
        else:
            s = format_fraction(coeff)
            if parts and coeff > 0:
                parts.append(f"+{s}*{var_names[var_idx]}")
            else:
                parts.append(f"{s}*{var_names[var_idx]}")

    if trans != 0:
        s = format_fraction(abs(trans))
        if trans > 0:
            parts.append(f"+{s}")
        else:
            parts.append(f"-{s}")

    return "".join(parts) if parts else "0"


def format_symop(components: list[tuple[list[Fraction], Fraction]]) -> str:
    """Format a full symmetry operation as 'expr,expr,expr'."""
    return ",".join(format_component(c, t) for c, t in components)


def combine_symop_cenop(symop: str, cenop_trans: list[Fraction]) -> str:
    """Apply centering translation to a symmetry operation."""
    components = parse_full_symop(symop)
    result = []
    for i, (coeffs, trans) in enumerate(components):
        new_trans = (trans + cenop_trans[i]) % 1
        result.append((coeffs, new_trans))
    return format_symop(result)


def expand_operations(symops: list[str], cenops: list[str]) -> list[str]:
    """Compute the full set of operations as cenops x symops."""
    cenop_translations = [parse_cenop_translation(c) for c in cenops]
    result = []
    for cenop_trans in cenop_translations:
        for symop in symops:
            result.append(combine_symop_cenop(symop, cenop_trans))
    return result


# ---------------------------------------------------------------------------
# Name handling
# ---------------------------------------------------------------------------


def upper_name(name: str) -> str:
    """Uppercase a space group name: 'R -3 c :H' -> 'R -3 C :H'."""
    return name.upper()


def compact(name: str) -> str:
    """Strip all spaces: 'P 21 21 21' -> 'P212121'."""
    return name.replace(" ", "")


def is_valid_old_name(name: str) -> bool:
    """Filter out invalid/garbage old names from syminfo.lib."""
    if not name:
        return False
    if "(" in name or ")" in name:
        return False
    # Skip full HM symbols with multiple '/' like "P 2/m 2/m 2/m"
    if name.count("/") > 1:
        return False
    return True


def monoclinic_short_forms(name: str) -> list[str]:
    """Generate short monoclinic aliases by dropping '1' positions.

    'P 1 2 1' -> ['P 2']
    'P 1 21 1' -> ['P 21']
    'C 1 2 1' -> ['C 2']
    'P 1 2/M 1' -> ['P 2/M']
    """
    # Match pattern: "X 1 Y 1" where X is a letter, Y is the content
    m = re.match(r"^([A-Z])\s+1\s+(.+?)\s+1$", name)
    if m:
        short = f"{m.group(1)} {m.group(2)}"
        return [short]
    return []


def cubic_short_forms(name: str) -> list[str]:
    """Generate aliases for cubic groups where numbers get merged.

    'I 2 3' -> ['I 23']
    'P 2 3' -> ['P 23']
    'F 2 3' -> ['F 23']
    """
    m = re.match(r"^([A-Z])\s+(\d+)\s+(\d+)$", name)
    if m:
        return [f"{m.group(1)} {m.group(2)}{m.group(3)}"]
    return []


# ---------------------------------------------------------------------------
# Main generation
# ---------------------------------------------------------------------------


def build_data(groups: list[SpaceGroup]) -> tuple[dict[str, list[str]], dict[str, str]]:
    """Build SYM_OPS and ALIASES from parsed space group data."""
    sym_ops: dict[str, list[str]] = {}
    aliases: dict[str, str] = {}
    # Track xHM -> ops for H/R aliasing
    xhm_ops: dict[str, list[str]] = {}

    for group in groups:
        xhm = group["xhm"]
        if not xhm:
            continue
        name = upper_name(xhm)
        ops = expand_operations(group["symops"], group["cenops"])

        # Store under canonical uppercase xHM name
        if name not in sym_ops:
            sym_ops[name] = ops
            xhm_ops[name] = ops

        # For names with setting suffixes (:1, :2, :H, :R),
        # also store variant without space before ':'
        no_space_colon = re.sub(r"\s+:", ":", name)
        if no_space_colon != name and no_space_colon not in sym_ops:
            sym_ops[no_space_colon] = ops

        # Compact form as alias
        comp = compact(name)
        if comp != name and comp not in sym_ops and comp not in aliases:
            aliases[comp] = name

        # Compact of no-space-colon variant
        if no_space_colon != name:
            comp_ns = compact(no_space_colon)
            if (
                comp_ns != no_space_colon
                and comp_ns not in sym_ops
                and comp_ns not in aliases
            ):
                aliases[comp_ns] = name

        # Old names as aliases
        for old_name in group["old_names"]:
            if not is_valid_old_name(old_name):
                continue
            old_upper = upper_name(old_name)
            if old_upper == name:
                continue
            if old_upper not in sym_ops and old_upper not in aliases:
                aliases[old_upper] = name
            # Compact form of old name
            old_comp = compact(old_upper)
            if (
                old_comp != old_upper
                and old_comp not in sym_ops
                and old_comp not in aliases
            ):
                aliases[old_comp] = name

        # Monoclinic short forms (P 2, P 21, C 2, C 21, etc.)
        for short in monoclinic_short_forms(name):
            if short not in sym_ops and short not in aliases:
                aliases[short] = name
            short_comp = compact(short)
            if (
                short_comp != short
                and short_comp not in sym_ops
                and short_comp not in aliases
            ):
                aliases[short_comp] = name

        # Cubic short forms (I 23, P 23, etc.)
        for short in cubic_short_forms(name):
            if short not in sym_ops and short not in aliases:
                aliases[short] = name
            short_comp = compact(short)
            if (
                short_comp != short
                and short_comp not in sym_ops
                and short_comp not in aliases
            ):
                aliases[short_comp] = name

    # --- Rhombohedral/Hexagonal aliases ---

    # For R groups in hexagonal setting (:H), generate H ... aliases
    for name in list(sym_ops.keys()):
        if not name.endswith(":H") and not name.endswith(" :H"):
            continue
        if not name.startswith("R "):
            continue

        # Extract core: "R -3 C :H" -> "-3 C", "R 3 :H" -> "3"
        m = re.match(r"R\s+(.+?)\s*:H$", name)
        if not m:
            continue
        core = m.group(1)

        # "H {core}" alias
        h_name = f"H {core}"
        if h_name not in sym_ops and h_name not in aliases:
            aliases[h_name] = name
        h_comp = compact(h_name)
        if h_comp != h_name and h_comp not in sym_ops and h_comp not in aliases:
            aliases[h_comp] = name

        # Lowercase glide plane variants: "H -3 C" -> "H -3 c"
        for letter in "ABCDEMNR":
            if letter in core:
                lower_core = core.replace(letter, letter.lower())
                h_lower = f"H {lower_core}"
                if h_lower not in sym_ops and h_lower not in aliases:
                    aliases[h_lower] = name

        # "H -3 2/C" form for "R -3 C :H"
        m2 = re.match(r"(-?\d+)\s+([A-Z])$", core)
        if m2:
            num, letter = m2.group(1), m2.group(2)
            for variant in [f"H {num} 2/{letter}", f"H {num} 2/{letter.lower()}"]:
                if variant not in sym_ops and variant not in aliases:
                    aliases[variant] = name

    # For unqualified R names (without :H/:R suffix), map to :H (PDB convention).
    # PDB files use "R 3", "R 32" etc. to mean the hexagonal setting.
    # Override any existing aliases from syminfo.lib old names (which map to :R).
    for name in list(xhm_ops.keys()):
        if not name.endswith(" :H"):
            continue
        if not name.startswith("R "):
            continue
        # "R 3 :H" -> unqualified "R 3"
        unqualified = name.rsplit(" :H", 1)[0]
        if unqualified not in sym_ops:
            aliases[unqualified] = name  # Override existing alias if any
        # Also compact form: "R3"
        uq_comp = compact(unqualified)
        if uq_comp != unqualified and uq_comp not in sym_ops:
            aliases[uq_comp] = name  # Override existing alias if any
        # Also number-merged forms: "R 3 2" -> "R 32"
        for short in cubic_short_forms(unqualified):
            if short not in sym_ops:
                aliases[short] = name

    # Clean up: remove aliases pointing to non-existent SYM_OPS keys
    aliases = {k: v for k, v in aliases.items() if v in sym_ops}
    # Remove aliases whose key is already a SYM_OPS key
    aliases = {k: v for k, v in aliases.items() if k not in sym_ops}

    return sym_ops, aliases


def generate_rust(sym_ops: dict[str, list[str]], aliases: dict[str, str]) -> str:
    """Generate the Rust source file content."""
    lines = []
    lines.append("//! Space group symmetry operation data")
    lines.append("//!")
    lines.append(
        "//! Generated from CCP4's syminfo.lib (International Tables for Crystallography, Vol. A)."
    )
    lines.append(
        "//! Do not edit manually. Regenerate with: scripts/generate_space_groups.py"
    )
    lines.append("")
    lines.append("use phf::phf_map;")
    lines.append("")
    lines.append("/// Mapping from space group name to symmetry operation strings.")
    lines.append(
        "pub static SYM_OPS: phf::Map<&'static str, &'static [&'static str]> = phf_map! {"
    )

    for name in sorted(sym_ops.keys()):
        ops = sym_ops[name]
        ops_str = ", ".join(f'"{op}"' for op in ops)
        lines.append(f'    "{name}" => &[{ops_str}],')

    lines.append("};")
    lines.append("")
    lines.append("/// Alias map for alternative space group name spellings.")
    lines.append(
        "pub static ALIASES: phf::Map<&'static str, &'static str> = phf_map! {"
    )

    for alias in sorted(aliases.keys()):
        target = aliases[alias]
        lines.append(f'    "{alias}" => "{target}",')

    lines.append("};")
    lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Verification mode
# ---------------------------------------------------------------------------


def parse_existing_rust(path: str) -> tuple[dict[str, list[str]], dict[str, str]]:
    """Parse an existing space_groups_data.rs to extract SYM_OPS and ALIASES."""
    with open(path) as f:
        content = f.read()

    sym_ops = {}
    aliases = {}

    sym_match = re.search(
        r"pub static SYM_OPS.*?phf_map!\s*\{(.*?)\};", content, re.DOTALL
    )
    if sym_match:
        for m in re.finditer(r'"([^"]+)"\s*=>\s*&\[(.*?)\]', sym_match.group(1)):
            name = m.group(1)
            ops_str = m.group(2)
            ops = [om.group(1) for om in re.finditer(r'"([^"]+)"', ops_str)]
            sym_ops[name] = ops

    alias_match = re.search(
        r"pub static ALIASES.*?phf_map!\s*\{(.*?)\};", content, re.DOTALL
    )
    if alias_match:
        for m in re.finditer(r'"([^"]+)"\s*=>\s*"([^"]+)"', alias_match.group(1)):
            aliases[m.group(1)] = m.group(2)

    return sym_ops, aliases


def verify(
    new_ops: dict[str, list[str]], new_aliases: dict[str, str], existing_path: str
) -> bool:
    """Compare generated data against existing file, report differences."""
    old_ops, old_aliases = parse_existing_rust(existing_path)

    ok = True
    old_keys = set(old_ops.keys())
    new_keys = set(new_ops.keys())

    missing = old_keys - new_keys
    extra = new_keys - old_keys

    if missing:
        print(f"SYM_OPS: {len(missing)} keys in old but not new:", file=sys.stderr)
        for k in sorted(missing):
            # Check if accessible via alias
            alias_target = new_aliases.get(k)
            suffix = f" (alias -> '{alias_target}')" if alias_target else ""
            print(f"  - '{k}' ({len(old_ops[k])} ops){suffix}", file=sys.stderr)
        ok = False

    if extra:
        print(
            f"SYM_OPS: {len(extra)} keys in new but not old (first 20):",
            file=sys.stderr,
        )
        for k in sorted(extra)[:20]:
            print(f"  + '{k}' ({len(new_ops[k])} ops)", file=sys.stderr)
        if len(extra) > 20:
            print(f"  ... and {len(extra) - 20} more", file=sys.stderr)

    # Check operations for common keys (compare as sets since order may differ)
    common = old_keys & new_keys
    order_mismatches = 0
    content_mismatches = 0
    for k in sorted(common):
        if old_ops[k] == new_ops[k]:
            continue
        old_set = set(old_ops[k])
        new_set = set(new_ops[k])
        if old_set == new_set:
            order_mismatches += 1
        else:
            content_mismatches += 1
            if content_mismatches <= 5:
                print(f"SYM_OPS content mismatch for '{k}':", file=sys.stderr)
                only_old = old_set - new_set
                only_new = new_set - old_set
                if only_old:
                    print(f"  only in old: {sorted(only_old)[:3]}", file=sys.stderr)
                if only_new:
                    print(f"  only in new: {sorted(only_new)[:3]}", file=sys.stderr)
                print(
                    f"  old count: {len(old_ops[k])}, new count: {len(new_ops[k])}",
                    file=sys.stderr,
                )
            ok = False

    if content_mismatches > 5:
        print(
            f"  ... and {content_mismatches - 5} more content mismatches",
            file=sys.stderr,
        )

    print(
        f"SYM_OPS: {len(common)} common, {order_mismatches} order-only diffs, "
        f"{content_mismatches} content diffs",
        file=sys.stderr,
    )

    # Check ALIASES
    old_akeys = set(old_aliases.keys())
    new_akeys = set(new_aliases.keys())

    a_missing = old_akeys - new_akeys
    a_extra = new_akeys - old_akeys

    if a_missing:
        print(f"ALIASES: {len(a_missing)} in old but not new:", file=sys.stderr)
        for k in sorted(a_missing):
            print(f"  - '{k}' => '{old_aliases[k]}'", file=sys.stderr)

    if a_extra:
        print(
            f"ALIASES: {len(a_extra)} in new but not old (first 20):", file=sys.stderr
        )
        for k in sorted(a_extra)[:20]:
            print(f"  + '{k}' => '{new_aliases[k]}'", file=sys.stderr)
        if len(a_extra) > 20:
            print(f"  ... and {len(a_extra) - 20} more", file=sys.stderr)

    a_value_mismatches = 0
    for k in sorted(old_akeys & new_akeys):
        if old_aliases[k] != new_aliases[k]:
            a_value_mismatches += 1
            if a_value_mismatches <= 10:
                print(
                    f"ALIASES value mismatch for '{k}': "
                    f"old='{old_aliases[k]}' new='{new_aliases[k]}'",
                    file=sys.stderr,
                )

    if a_value_mismatches > 10:
        print(
            f"  ... and {a_value_mismatches - 10} more alias mismatches",
            file=sys.stderr,
        )

    summary = "OK" if ok else "DIFFERENCES FOUND"
    print(
        f"\n{summary}: {len(new_keys)} SYM_OPS, {len(new_akeys)} ALIASES generated "
        f"(old had {len(old_keys)} SYM_OPS, {len(old_akeys)} ALIASES)",
        file=sys.stderr,
    )

    return ok


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main():
    if len(sys.argv) < 2:
        print(
            f"Usage: {sys.argv[0]} <syminfo.lib> [--verify <existing.rs>]",
            file=sys.stderr,
        )
        sys.exit(1)

    syminfo_path = sys.argv[1]
    verify_mode = "--verify" in sys.argv
    verify_path = None
    if verify_mode:
        idx = sys.argv.index("--verify")
        if idx + 1 < len(sys.argv):
            verify_path = sys.argv[idx + 1]
        else:
            print("--verify requires path to existing .rs file", file=sys.stderr)
            sys.exit(1)

    groups = parse_syminfo(syminfo_path)
    sym_ops, aliases = build_data(groups)

    if verify_mode:
        assert verify_path is not None
        ok = verify(sym_ops, aliases, verify_path)
        sys.exit(0 if ok else 1)
    else:
        print(generate_rust(sym_ops, aliases), end="")
        print(
            f"Generated {len(sym_ops)} SYM_OPS entries and {len(aliases)} ALIASES.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
