"""
Command-line interface for PyMOL-RS.

This module provides the entry point for the `pymol-rs` command.
"""

import argparse
import sys
from pathlib import Path


def main() -> int:
    """Main entry point for the pymol-rs command."""
    parser = argparse.ArgumentParser(
        prog="pymol-rs",
        description="PyMOL-RS molecular visualization",
    )
    parser.add_argument(
        "--ipc",
        metavar="SOCKET_PATH",
        type=Path,
        help="Enable IPC server for external control (e.g., from Python)",
    )
    parser.add_argument(
        "--headless",
        action="store_true",
        help="Run in headless mode (no window displayed)",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress log output (used when spawned by Python code)",
    )
    parser.add_argument(
        "files",
        nargs="*",
        metavar="FILE",
        help="Files to load at startup",
    )

    args = parser.parse_args()

    # Import here to avoid circular imports and slow startup for --help
    from pymol_rs._pymol_rs import run_gui

    try:
        run_gui(
            ipc_socket=args.ipc,
            headless=args.headless,
            files=args.files if args.files else None,
            quiet=args.quiet,
        )
        return 0
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
