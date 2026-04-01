#!/usr/bin/env bash
set -euo pipefail

if [ $# -ne 1 ]; then
    echo "Usage: $0 <version>"
    echo "Example: $0 0.3.0"
    exit 1
fi

NEW="$1"

if ! [[ "$NEW" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
    echo "Error: version must be in semver format X.Y.Z (e.g., 0.3.0)"
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$ROOT"

# Detect current version from workspace Cargo.toml
OLD=$(sed -n '/\[workspace\.package\]/,/\[/{s/^version = "\(.*\)"/\1/p;}' Cargo.toml)

if [ -z "$OLD" ]; then
    echo "Error: could not detect current version from Cargo.toml"
    exit 1
fi

if [ "$OLD" = "$NEW" ]; then
    echo "Version is already $NEW — nothing to do."
    exit 0
fi

echo "Updating version: $OLD -> $NEW"
echo ""

update_file() {
    local file="$1"
    local pattern="$2"
    local replacement="$3"
    if [ ! -f "$file" ]; then
        echo "  SKIP  $file (not found)"
        return
    fi
    local count
    count=$(grep -c "$pattern" "$file" || true)
    sed -i '' "s|$pattern|$replacement|g" "$file"
    echo "  OK    $file ($count replacement(s))"
}

# Workspace root Cargo.toml (workspace.package version + dependency versions)
update_file "Cargo.toml" \
    "version = \"$OLD\"" \
    "version = \"$NEW\""

# Python crate (excluded from workspace)
update_file "python/Cargo.toml" \
    "version = \"$OLD\"" \
    "version = \"$NEW\""

# Python pyproject.toml
update_file "python/pyproject.toml" \
    "version = \"$OLD\"" \
    "version = \"$NEW\""

# Web crate (excluded from workspace)
update_file "web/Cargo.toml" \
    "version = \"$OLD\"" \
    "version = \"$NEW\""

# Web package.json
update_file "web/package.json" \
    "\"version\": \"$OLD\"" \
    "\"version\": \"$NEW\""

# Makefile (match any existing version after VERSION :=)
if [ -f "Makefile" ]; then
    sed -i '' "s/^VERSION[[:space:]]*:= .*/VERSION        := $NEW/" Makefile
    echo "  OK    Makefile"
else
    echo "  SKIP  Makefile (not found)"
fi

echo ""
echo "Done. Version updated to $NEW."
