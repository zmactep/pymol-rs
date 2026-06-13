#!/usr/bin/env python3
"""Pack a macOS .iconset directory into an .icns file.

This uses the PNG-backed ICNS chunks that iconutil emits for modern app icons.
It intentionally avoids image libraries; the Makefile already creates the
correct PNG sizes with sips.
"""

from __future__ import annotations

import struct
import sys
from pathlib import Path


CHUNKS = (
    ("icp4", "icon_16x16.png", (16, 16)),
    ("icp5", "icon_32x32.png", (32, 32)),
    ("ic11", "icon_16x16@2x.png", (32, 32)),
    ("ic12", "icon_32x32@2x.png", (64, 64)),
    ("ic07", "icon_128x128.png", (128, 128)),
    ("ic13", "icon_128x128@2x.png", (256, 256)),
    ("ic08", "icon_256x256.png", (256, 256)),
    ("ic14", "icon_256x256@2x.png", (512, 512)),
    ("ic09", "icon_512x512.png", (512, 512)),
    ("ic10", "icon_512x512@2x.png", (1024, 1024)),
)


def png_size(data: bytes) -> tuple[int, int]:
    if len(data) < 24 or not data.startswith(b"\x89PNG\r\n\x1a\n"):
        raise ValueError("not a PNG")
    return struct.unpack(">II", data[16:24])


def pack_chunk(kind: str, path: Path, expected_size: tuple[int, int]) -> bytes:
    data = path.read_bytes()
    actual_size = png_size(data)
    if actual_size != expected_size:
        actual = f"{actual_size[0]}x{actual_size[1]}"
        expected = f"{expected_size[0]}x{expected_size[1]}"
        raise ValueError(f"{path} is {actual}, expected {expected}")
    return kind.encode("ascii") + struct.pack(">I", len(data) + 8) + data


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: iconset_to_icns.py <AppIcon.iconset> <AppIcon.icns>", file=sys.stderr)
        return 2

    iconset = Path(sys.argv[1])
    output = Path(sys.argv[2])

    if not iconset.is_dir():
        print(f"{iconset}: not an iconset directory", file=sys.stderr)
        return 2

    try:
        chunks = [pack_chunk(kind, iconset / name, size) for kind, name, size in CHUNKS]
    except (FileNotFoundError, ValueError) as exc:
        print(f"{iconset}: invalid iconset: {exc}", file=sys.stderr)
        return 1

    payload = b"".join(chunks)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(b"icns" + struct.pack(">I", len(payload) + 8) + payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
