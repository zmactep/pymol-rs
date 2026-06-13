"""CLI entry point for patinae."""

import argparse
import sys


def main():
    parser = argparse.ArgumentParser(description="Patinae molecular visualization")
    parser.add_argument("files", nargs="*", help="Files to load at startup")
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="Suppress log output"
    )
    args = parser.parse_args()

    from patinae import cmd

    for f in args.files:
        cmd.load(f, quiet=args.quiet)

    # Interactive mode
    try:
        while True:
            try:
                line = input("Patinae> ")
            except EOFError:
                break
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line in ("quit", "exit"):
                break
            try:
                cmd.do(line)
            except Exception as e:
                print(f"Error: {e}", file=sys.stderr)
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
