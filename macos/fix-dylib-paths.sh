#!/bin/bash
# Fix libpython dylib references in the .app bundle so the bundled
# Python is found at runtime via @executable_path.
#
# Usage: bash macos/fix-dylib-paths.sh <path-to-app-dir>

set -euo pipefail

APP_DIR="$1"

# Find the bundled libpython dylib
PYTHON_DYLIB=$(find "$APP_DIR/Contents/Resources/python/lib" -maxdepth 1 -name 'libpython*.dylib' | head -1)

if [ -z "$PYTHON_DYLIB" ]; then
    echo "No libpython dylib found in bundle — skipping dylib path fixups"
    exit 0
fi

PYTHON_DYLIB_NAME=$(basename "$PYTHON_DYLIB")
NEW_PATH="@executable_path/../Resources/python/lib/$PYTHON_DYLIB_NAME"

echo "Fixing dylib references: $PYTHON_DYLIB_NAME -> $NEW_PATH"

# Fix the dylib's own install name
install_name_tool -id "$NEW_PATH" "$PYTHON_DYLIB"
echo "  Fixed install name: $PYTHON_DYLIB"

# Fix references in the main binary and plugin dylibs
for binary in "$APP_DIR/Contents/MacOS/"* "$APP_DIR/Contents/PlugIns/"*.dylib; do
    [ -f "$binary" ] || continue
    OLD_PATH=$(otool -L "$binary" 2>/dev/null | grep 'libpython[0-9]' | awk '{print $1}' | head -1 || true)
    if [ -n "$OLD_PATH" ] && [ "$OLD_PATH" != "$NEW_PATH" ]; then
        install_name_tool -change "$OLD_PATH" "$NEW_PATH" "$binary"
        echo "  Fixed $(basename "$binary"): $OLD_PATH -> $NEW_PATH"
    fi
done

echo "Done"
